/******************************************************************************
 * Copyright 2019 The Apollo Authors. All Rights Reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *****************************************************************************/

#include "modules/planning/math/piecewise_jerk/piecewise_jerk_path_problem.h"

#include "cyber/common/log.h"
#include "modules/planning/common/planning_gflags.h"

namespace apollo {
namespace planning {

PiecewiseJerkPathProblem::PiecewiseJerkPathProblem(
    const size_t num_of_knots, const double delta_s,
    const std::array<double, 3>& x_init)
    : PiecewiseJerkProblem(num_of_knots, delta_s, x_init) {}



// 1、计算 Hessian 矩阵
// 2、转变成 OSQP求解器需要的素材
// 其中涉及到 csc_matrix 的相关参数计算，可以查看这个链接了解https://blog.csdn.net/weixin_34945803/article/details/106576629
void PiecewiseJerkPathProblem::CalculateKernel(std::vector<c_float>* P_data,
                                               std::vector<c_int>* P_indices,
                                               std::vector<c_int>* P_indptr) {
  const int n = static_cast<int>(num_of_knots_);  // 待平滑的点数
  const int num_of_variables = 3 * n;             // 决策变量的数目
  const int num_of_nonzeros = num_of_variables + (n - 1);
  std::vector<std::vector<std::pair<c_int, c_float>>> columns(num_of_variables);
  int value_index = 0;

  // Hessian矩阵 每个元素的计算参考 https://mp.weixin.qq.com/s?__biz=MzI1NjkxOTMyNQ==&mid=2247486492&idx=1&sn=2bfb421186364cbd6b00fb278a4a70fe&chksm=ea1e1e6edd699778417d0b1dc6e009b445c8c03218d555f5c872ecc08f49f4b1767fcabdd029&mpshare=1&scene=1&srcid=#rd
  
  // x(i)^2 * (w_x + w_x_ref[i]), w_x_ref might be a uniform value for all x(i) or piecewise values for different x(i)
  for (int i = 0; i < n - 1; ++i) {
    columns[i].emplace_back(i, (weight_x_ + weight_x_ref_vec_[i]) /
                               (scale_factor_[0] * scale_factor_[0]));
    ++value_index;
  }
  // x(n-1)^2 * (w_x + w_x_ref[n-1] + w_end_x)
  columns[n - 1].emplace_back(
      n - 1, (weight_x_ + weight_x_ref_vec_[n - 1] + weight_end_state_[0]) /
                 (scale_factor_[0] * scale_factor_[0]));
  ++value_index;


  // x(i)' ^2 * w_dx
  for (int i = 0; i < n - 1; ++i) {
    columns[n + i].emplace_back(
        n + i, weight_dx_ / (scale_factor_[1] * scale_factor_[1]));
    ++value_index;
  } 
  // x(n-1)' ^2 * (w_dx + w_end_dx)
  columns[2 * n - 1].emplace_back(2 * n - 1,
                                  (weight_dx_ + weight_end_state_[1]) /
                                      (scale_factor_[1] * scale_factor_[1]));
  ++value_index;

  // 在计算 l'' ^2 的权重时， 考虑了  l''' ^2 的权重
  // 其中 l''' = ( l''(i)-l''(i-1) ) / delta_s
  // l'''^2 = l''(i)^2 /delta_s^2   +  l''(i-1)^2 /delta_s^2  -  2 * l''(i)*l''(i-1) /delta_s^2
  auto delta_s_square = delta_s_ * delta_s_;
  // x(i)'' ^2  * (w_ddx + 2 * w_dddx / delta_s^2)
  // 所以每个 l''(i)^2 的权重有两部分组成，一个是w_ddx， 一个是w_dddx
  columns[2 * n].emplace_back(2 * n,
                              (weight_ddx_ + weight_dddx_ / delta_s_square) /
                                  (scale_factor_[2] * scale_factor_[2]));
  ++value_index;
  for (int i = 1; i < n - 1; ++i) {
    columns[2 * n + i].emplace_back(
        2 * n + i, (weight_ddx_ + 2.0 * weight_dddx_ / delta_s_square) /
                       (scale_factor_[2] * scale_factor_[2]));
    ++value_index;
  }
  columns[3 * n - 1].emplace_back(
      3 * n - 1,
      (weight_ddx_ + weight_dddx_ / delta_s_square + weight_end_state_[2]) /
          (scale_factor_[2] * scale_factor_[2]));
  ++value_index;


  // -2 * w_dddx / delta_s^2 * x(i)'' * x(i + 1)''
  // hession矩阵的 右下角这个 n*n的矩阵，除了对角线元素，它右上角还有一排元素
  /***       |    o                  | 
   *         |         o             |
   *         |              o  o     |
   *         |                 o  o  |
   * ***/    |                    0  |
  for (int i = 0; i < n - 1; ++i) {
    columns[2 * n + i].emplace_back(2 * n + i + 1,
                                    (-2.0 * weight_dddx_ / delta_s_square) /
                                        (scale_factor_[2] * scale_factor_[2]));
    ++value_index;
  }

  CHECK_EQ(value_index, num_of_nonzeros);

  // 这里比较特殊，目标函数的 hession 对角矩阵，每个元素对应的是权重，所以P_indptr直接累加就好，
  // 这里ind_p 便是记录到目前为止  非零元素的个数
  int ind_p = 0;
  for (int i = 0; i < num_of_variables; ++i) {
    P_indptr->push_back(ind_p);
    for (const auto& row_data_pair : columns[i]) {
      P_data->push_back(row_data_pair.second * 2.0); // P_data来记录 hession矩阵的元素        <------------------
      P_indices->push_back(row_data_pair.first);     // P_indices来记录各个元素所在列的行号                //    |
      ++ind_p;                                                                                           //    |
    }                                                                                                    //    |
  }                                                                                                      //    |
  P_indptr->push_back(ind_p);                                                                            //    |
}                                                                                                        //    |
                                                                                                         //    |
                                                                                                         //    |
                                                                                                         //    |
// 最后，我们将进行offset 补偿，这里主要指的是 最后的参考线要考虑referline这一因素，即初始解。                  //    |
// 保证尽可能不要有太大偏差，这样有可能给车辆带来不稳定因素，这里主要是给目标函数进行补偿，目标函数的 ref一项。    //   |
                                                                                                         //    |
// 其实目标函数在横向位移上有两项： l^2+(l-ref)^2,因此可以看到为什么在目标函数里，l^2的系数乘以2，--------------------|
// 在这里将第二项进行了拆解，于是有了offset。 即 -2ref*i，这个就对应了
// 至于为什么不考虑ref^2，因为它是个非负实数，并不包含任何变量，因此不影响梯度下降-->不影响求解-->此处省略
void PiecewiseJerkPathProblem::CalculateOffset(std::vector<c_float>* q) {
  CHECK_NOTNULL(q);
  const int n = static_cast<int>(num_of_knots_);
  const int kNumParam = 3 * n;
  q->resize(kNumParam, 0.0);

  if (has_x_ref_) {
    // l^2+(l-ref)^2 拆开项中的 -2ref*i
    for (int i = 0; i < n; ++i) {
      q->at(i) += -2.0 * weight_x_ref_vec_.at(i) * x_ref_[i] / scale_factor_[0];
    }
  }

  if (has_end_state_ref_) {
    q->at(n - 1) +=
        -2.0 * weight_end_state_[0] * end_state_ref_[0] / scale_factor_[0];
    q->at(2 * n - 1) +=
        -2.0 * weight_end_state_[1] * end_state_ref_[1] / scale_factor_[1];
    q->at(3 * n - 1) +=
        -2.0 * weight_end_state_[2] * end_state_ref_[2] / scale_factor_[2];
  }
}

}  // namespace planning
}  // namespace apollo

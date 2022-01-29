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

/*
 * @file
 */

#include "modules/planning/open_space/coarse_trajectory_generator/grid_search.h"

namespace apollo {
namespace planning {

GridSearch::GridSearch(const PlannerOpenSpaceConfig& open_space_conf) {
  xy_grid_resolution_ =
      open_space_conf.warm_start_config().grid_a_star_xy_resolution();
  node_radius_ = open_space_conf.warm_start_config().node_radius();
}


// 计算两点欧式距离
double GridSearch::EuclidDistance(const double x1, const double y1,
                                  const double x2, const double y2) {
  return std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}


//检测是否越界 or 存在碰撞
bool GridSearch::CheckConstraints(std::shared_ptr<Node2d> node) {
  const double node_grid_x = node->GetGridX();
  const double node_grid_y = node->GetGridY();
  // 超出边界，fasle
  if (node_grid_x > max_grid_x_ || node_grid_x < 0 ||
      node_grid_y > max_grid_y_ || node_grid_y < 0) {
    return false;
  }
  // 没有障碍物，true
  if (obstacles_linesegments_vec_.empty()) {
    return true;
  }
  // 与obstacles有重叠，false
  for (const auto& obstacle_linesegments : obstacles_linesegments_vec_) {
    for (const common::math::LineSegment2d& linesegment :
         obstacle_linesegments) {
      if (linesegment.DistanceTo({node->GetGridX(), node->GetGridY()}) <
          node_radius_) {
        return false;
      }
    }
  }
  return true;
}



// 扩展出周围8个节点，此处不考虑碰撞、越界，先拿到再做下一步的处理
std::vector<std::shared_ptr<Node2d>> 
GridSearch::GenerateNextNodes(std::shared_ptr<Node2d> current_node) {
    double current_node_x = current_node->GetGridX();
    double current_node_y = current_node->GetGridY();
    double current_node_path_cost = current_node->GetPathCost();
    double diagonal_distance = std::sqrt(2.0);

    // 依次扩展周围 8 个节点，并初始化path_cost
    // 但是父节点的设置是在GenerateAStarPath()进行
    std::vector<std::shared_ptr<Node2d>> next_nodes;
    std::shared_ptr<Node2d> up = std::make_shared<Node2d>(current_node_x, current_node_y + 1.0, XYbounds_);
    up->SetPathCost(current_node_path_cost + 1.0);

    std::shared_ptr<Node2d> up_right = std::make_shared<Node2d>(current_node_x + 1.0, current_node_y + 1.0, XYbounds_);
    up_right->SetPathCost(current_node_path_cost + diagonal_distance);

    std::shared_ptr<Node2d> right = std::make_shared<Node2d>(current_node_x + 1.0, current_node_y, XYbounds_);
    right->SetPathCost(current_node_path_cost + 1.0);

    std::shared_ptr<Node2d> down_right = std::make_shared<Node2d>(current_node_x + 1.0, current_node_y - 1.0, XYbounds_);
    down_right->SetPathCost(current_node_path_cost + diagonal_distance);

    std::shared_ptr<Node2d> down = std::make_shared<Node2d>(current_node_x, current_node_y - 1.0, XYbounds_);
    down->SetPathCost(current_node_path_cost + 1.0);

    std::shared_ptr<Node2d> down_left = std::make_shared<Node2d>(current_node_x - 1.0, current_node_y - 1.0, XYbounds_);
    down_left->SetPathCost(current_node_path_cost + diagonal_distance);

    std::shared_ptr<Node2d> left = std::make_shared<Node2d>(current_node_x - 1.0, current_node_y, XYbounds_);
    left->SetPathCost(current_node_path_cost + 1.0);

    std::shared_ptr<Node2d> up_left = std::make_shared<Node2d>(current_node_x - 1.0, current_node_y + 1.0, XYbounds_);
    up_left->SetPathCost(current_node_path_cost + diagonal_distance);

    // 以up为起点，顺时针转一圈
    next_nodes.emplace_back(up);
    next_nodes.emplace_back(up_right);
    next_nodes.emplace_back(right);
    next_nodes.emplace_back(down_right);
    next_nodes.emplace_back(down);
    next_nodes.emplace_back(down_left);
    next_nodes.emplace_back(left);
    next_nodes.emplace_back(up_left);
    return next_nodes;
}





// node数据格式
/*  class Node2d{
      ...
    private:
      int grid_x_ = 0;
      int grid_y_ = 0;
      double path_cost_ = 0.0;    g值
      double heuristic_ = 0.0;    h值
      double cost_ = 0.0;         f = g + h;值
      std::string index_;
      std::shared_ptr<Node2d> pre_node_ = nullptr;
    }
*/
// A* 搜索算法
bool GridSearch::GenerateAStarPath(
    const double sx, const double sy, const double ex, const double ey,
    const std::vector<double>& XYbounds,
    const std::vector<std::vector<common::math::LineSegment2d>>&
        obstacles_linesegments_vec,
    GridAStartResult* result) {

  // 初始化
  std::priority_queue<std::pair<std::string, double>,
                      std::vector<std::pair<std::string, double>>, cmp>
      open_pq;
  std::unordered_map<std::string, std::shared_ptr<Node2d>> open_set;
  std::unordered_map<std::string, std::shared_ptr<Node2d>> close_set;
  // 边界、起点、终点
  XYbounds_ = XYbounds;
  std::shared_ptr<Node2d> start_node = std::make_shared<Node2d>(sx, sy, xy_grid_resolution_, XYbounds_);
  std::shared_ptr<Node2d> end_node   = std::make_shared<Node2d>(ex, ey, xy_grid_resolution_, XYbounds_);
  std::shared_ptr<Node2d> final_node_ = nullptr;
  obstacles_linesegments_vec_ = obstacles_linesegments_vec;
  open_set.emplace(start_node->GetIndex(), start_node);
  open_pq.emplace(start_node->GetIndex(), start_node->GetCost());

  // Grid a star begins
  size_t explored_node_num = 0;
  while (!open_pq.empty()) {
    std::string current_id = open_pq.top().first;
    open_pq.pop();
    std::shared_ptr<Node2d> current_node = open_set[current_id];

    // 退出条件
    if (*(current_node) == *(end_node)) {
      final_node_ = current_node;
      break;
    }

    close_set.emplace(current_node->GetIndex(), current_node);
    // 得到周围节点
    std::vector<std::shared_ptr<Node2d>> next_nodes = std::move(GenerateNextNodes(current_node));
    for (auto& next_node : next_nodes) {
      // 超出边界、存在碰撞，跳过
      if (!CheckConstraints(next_node)) {
        continue;
      }
      // 在 close list中，跳过
      if (close_set.find(next_node->GetIndex()) != close_set.end()) {
        continue;
      }

      // 不在open list，初始化并压入
      if (open_set.find(next_node->GetIndex()) == open_set.end()) {
        ++explored_node_num;  // 记录已经扩展的node数目
        // 存入 h 值
        next_node->SetHeuristic(
            EuclidDistance(next_node->GetGridX(), next_node->GetGridY(),
                           end_node->GetGridX(), end_node->GetGridY()));
        // 存入父节点
        next_node->SetPreNode(current_node);
        // 压入open list
        open_set.emplace(next_node->GetIndex(), next_node);
        open_pq.emplace(next_node->GetIndex(), next_node->GetCost());
      } else {
        /*这里应该缺少了已经存在于open list,是否需要修改cost和父节点的情况*/
      }
    }
  }

  if (final_node_ == nullptr) {
    AERROR << "Grid A searching return null ptr(open_set ran out)";
    return false;
  }
  LoadGridAStarResult(result);
  ADEBUG << "explored node num is " << explored_node_num;
  return true;
}




// 基于  广度优先/Dijkstra  方法初始化heuristic，可以直接查表得到h，降低在线计算的耗时
// 注意：这个是给 HybridAStar 用的
bool GridSearch::GenerateDpMap(
    const double ex, const double ey, const std::vector<double>& XYbounds,
    const std::vector<std::vector<common::math::LineSegment2d>>&
        obstacles_linesegments_vec) {
  std::priority_queue<std::pair<std::string, double>,
                      std::vector<std::pair<std::string, double>>, 
                      cmp>
      open_pq;
  std::unordered_map<std::string, std::shared_ptr<Node2d>> open_set;
  dp_map_ = decltype(dp_map_)();
  XYbounds_ = XYbounds;
  // XYbounds with xmin, xmax, ymin, ymax
  max_grid_y_ = std::round((XYbounds_[3] - XYbounds_[2]) / xy_grid_resolution_);
  max_grid_x_ = std::round((XYbounds_[1] - XYbounds_[0]) / xy_grid_resolution_);
  std::shared_ptr<Node2d> end_node = std::make_shared<Node2d>(ex, ey, xy_grid_resolution_, XYbounds_);
  obstacles_linesegments_vec_ = obstacles_linesegments_vec;
  // 压入终点
  open_set.emplace(end_node->GetIndex(), end_node);
  open_pq.emplace(end_node->GetIndex(), end_node->GetCost());

  // Grid a star begins
  size_t explored_node_num = 0;
  while (!open_pq.empty()) {
    const std::string current_id = open_pq.top().first;
    open_pq.pop();
    std::shared_ptr<Node2d> current_node = open_set[current_id];
    dp_map_.emplace(current_node->GetIndex(), current_node);
    std::vector<std::shared_ptr<Node2d>> next_nodes = std::move(GenerateNextNodes(current_node));
    for (auto& next_node : next_nodes) {
      if (!CheckConstraints(next_node)) {
        continue;
      }
      // 在 dp_map_ 中，此处dp_mp_相当于close list
      if (dp_map_.find(next_node->GetIndex()) != dp_map_.end()) {
        continue;
      }
      if (open_set.find(next_node->GetIndex()) == open_set.end()) {
        // 未被扩展
        ++explored_node_num;
        next_node->SetPreNode(current_node);
        open_set.emplace(next_node->GetIndex(), next_node);
        open_pq.emplace(next_node->GetIndex(), next_node->GetCost());
      } else {
        // 在open lsit中，看要不要更新cost
        if (open_set[next_node->GetIndex()]->GetCost() > next_node->GetCost()) {
          open_set[next_node->GetIndex()]->SetCost(next_node->GetCost());
          open_set[next_node->GetIndex()]->SetPreNode(current_node);
        }
      }
    }
  }
  ADEBUG << "explored node num is " << explored_node_num;
  return true;
}


// 根据index 得到当前node到终点的  heuristic cost
double GridSearch::CheckDpMap(const double sx, const double sy) {
  std::string index = Node2d::CalcIndex(sx, sy, xy_grid_resolution_, XYbounds_);
  if (dp_map_.find(index) != dp_map_.end()) {
    return dp_map_[index]->GetCost() * xy_grid_resolution_;
  } else {
    return std::numeric_limits<double>::infinity();
  }
}



// 回溯找到父节点，得到路径
void GridSearch::LoadGridAStarResult(GridAStartResult* result) {
  (*result).path_cost = final_node_->GetPathCost() * xy_grid_resolution_;
  std::shared_ptr<Node2d> current_node = final_node_;
  std::vector<double> grid_a_x;
  std::vector<double> grid_a_y;
  while (current_node->GetPreNode() != nullptr) {
    grid_a_x.push_back(current_node->GetGridX() * xy_grid_resolution_ +
                       XYbounds_[0]);
    grid_a_y.push_back(current_node->GetGridY() * xy_grid_resolution_ +
                       XYbounds_[2]);
    current_node = current_node->GetPreNode();
  }
  std::reverse(grid_a_x.begin(), grid_a_x.end());
  std::reverse(grid_a_y.begin(), grid_a_y.end());
  (*result).x = std::move(grid_a_x);
  (*result).y = std::move(grid_a_y);
}
}  // namespace planning
}  // namespace apollo

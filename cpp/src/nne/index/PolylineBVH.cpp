#include "nne/index/PolylineBVH.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <utility>

namespace fastnsa::nne {
namespace {

constexpr std::size_t leaf_size = 8;

double center_x(const Envelope& bounds) noexcept {
    return bounds.min_x + (bounds.max_x - bounds.min_x) * 0.5;
}

double center_y(const Envelope& bounds) noexcept {
    return bounds.min_y + (bounds.max_y - bounds.min_y) * 0.5;
}

}  // namespace

PolylineBVH::PolylineBVH(const std::vector<Polyline>& polylines) : polylines_(&polylines) {
    if (polylines.empty()) {
        throw std::invalid_argument("Cannot build a nearest-edge index for an empty network.");
    }
    edge_order_.resize(polylines.size());
    std::iota(edge_order_.begin(), edge_order_.end(), 0);
    nodes_.reserve(polylines.size() * 2);
    build(0, edge_order_.size());
}

int PolylineBVH::build(std::size_t begin, std::size_t end) {
    const int node_index = static_cast<int>(nodes_.size());
    nodes_.push_back({});
    auto bounds = (*polylines_)[edge_order_[begin]].envelope();
    for (std::size_t index = begin + 1; index < end; ++index) {
        bounds.expand((*polylines_)[edge_order_[index]].envelope());
    }

    nodes_[node_index].bounds = bounds;
    nodes_[node_index].begin = begin;
    nodes_[node_index].end = end;
    if (end - begin <= leaf_size) {
        return node_index;
    }

    const bool split_x = (bounds.max_x - bounds.min_x) >= (bounds.max_y - bounds.min_y);
    const std::size_t middle = begin + (end - begin) / 2;
    std::nth_element(
        edge_order_.begin() + static_cast<std::ptrdiff_t>(begin),
        edge_order_.begin() + static_cast<std::ptrdiff_t>(middle),
        edge_order_.begin() + static_cast<std::ptrdiff_t>(end),
        [&](int first, int second) {
            const auto& first_bounds = (*polylines_)[first].envelope();
            const auto& second_bounds = (*polylines_)[second].envelope();
            const double first_center = split_x ? center_x(first_bounds) : center_y(first_bounds);
            const double second_center = split_x ? center_x(second_bounds) : center_y(second_bounds);
            return first_center == second_center ? first < second : first_center < second_center;
        }
    );
    const int left = build(begin, middle);
    const int right = build(middle, end);
    nodes_[node_index].left = left;
    nodes_[node_index].right = right;
    return node_index;
}

EdgeMatch PolylineBVH::nearest(const Point2D& point) const {
    using QueueEntry = std::pair<double, int>;
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> queue;
    queue.push({nodes_[0].bounds.distance_squared(point), 0});
    EdgeMatch best{-1, std::numeric_limits<double>::infinity(), 0.0};

    while (!queue.empty()) {
        const auto [lower_bound, node_index] = queue.top();
        queue.pop();
        if (lower_bound > best.distance_squared) {
            break;
        }
        const auto& node = nodes_[node_index];
        if (node.leaf()) {
            for (std::size_t index = node.begin; index < node.end; ++index) {
                const int edge_id = edge_order_[index];
                const auto projection = (*polylines_)[edge_id].project(point);
                if (projection.distance_squared < best.distance_squared ||
                    (projection.distance_squared == best.distance_squared && edge_id < best.edge_id)) {
                    best = {edge_id, projection.distance_squared, projection.offset};
                }
            }
        } else {
            for (const int child : {node.left, node.right}) {
                const double child_distance = nodes_[child].bounds.distance_squared(point);
                if (child_distance <= best.distance_squared) {
                    queue.push({child_distance, child});
                }
            }
        }
    }
    return best;
}

std::vector<int> PolylineBVH::candidates_within(
    const Point2D& point,
    double maximum_distance_squared
) const {
    std::vector<int> candidates;
    std::vector<int> stack{0};
    while (!stack.empty()) {
        const int node_index = stack.back();
        stack.pop_back();
        const auto& node = nodes_[node_index];
        if (node.bounds.distance_squared(point) > maximum_distance_squared) {
            continue;
        }
        if (node.leaf()) {
            for (std::size_t index = node.begin; index < node.end; ++index) {
                const int edge_id = edge_order_[index];
                if ((*polylines_)[edge_id].project(point).distance_squared <= maximum_distance_squared) {
                    candidates.push_back(edge_id);
                }
            }
        } else {
            stack.push_back(node.left);
            stack.push_back(node.right);
        }
    }
    return candidates;
}

std::vector<int> PolylineBVH::within_nearest_distance(
    const Point2D& point,
    double extra_distance,
    std::size_t maximum_count
) const {
    const auto closest = nearest(point);
    const double maximum_distance = std::sqrt(closest.distance_squared) + extra_distance;
    auto candidates = candidates_within(point, maximum_distance * maximum_distance);
    candidates.erase(
        std::remove_if(candidates.begin(), candidates.end(), [&](int edge_id) {
            return (*polylines_)[edge_id].project(point).distance_squared >=
                maximum_distance * maximum_distance;
        }),
        candidates.end()
    );
    std::sort(candidates.begin(), candidates.end(), [&](int first, int second) {
        const double first_distance = (*polylines_)[first].project(point).distance_squared;
        const double second_distance = (*polylines_)[second].project(point).distance_squared;
        return first_distance == second_distance ? first < second : first_distance < second_distance;
    });
    if (candidates.size() > maximum_count) {
        candidates.resize(maximum_count);
    }
    return candidates;
}

}  // namespace fastnsa::nne

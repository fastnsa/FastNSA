#include "nne/core/NNE.h"

#include "nne/core/SpatialGrid.h"
#include "nne/core/Threading.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <stdexcept>
#include <thread>

namespace fastnsa::nne {

NNE::NNE(std::vector<Polyline> polylines)
    : polylines_(std::move(polylines)), index_(polylines_) {}

std::vector<int> NNE::query_exact_mt(
    const std::vector<QueryPoint>& query_points,
    int num_threads
) const {
    std::vector<int> results(query_points.size(), -1);
    const std::size_t count = query_points.size();
    const int worker_count = effective_thread_count(num_threads, count, 256);
    auto query_range = [&](std::size_t begin, std::size_t end) {
        for (std::size_t index = begin; index < end; ++index) {
            results[index] = index_.nearest(query_points[index]).edge_id;
        }
    };
    if (worker_count == 1) {
        query_range(0, count);
        return results;
    }

    std::vector<std::thread> workers;
    workers.reserve(static_cast<std::size_t>(worker_count));
    const std::size_t chunk_size = (count + worker_count - 1) /
        static_cast<std::size_t>(worker_count);
    try {
        for (int worker = 0; worker < worker_count; ++worker) {
            const std::size_t begin = static_cast<std::size_t>(worker) * chunk_size;
            const std::size_t end = std::min(count, begin + chunk_size);
            workers.emplace_back(query_range, begin, end);
        }
    } catch (...) {
        for (auto& worker : workers) {
            worker.join();
        }
        throw;
    }
    for (auto& worker : workers) {
        worker.join();
    }
    return results;
}

std::vector<int> NNE::query_nearest_mt(
    const std::vector<QueryPoint>& query_points,
    int method,
    double grid_size,
    int num_threads
) const {
    if (method == 0) {
        return query_exact_mt(query_points, num_threads);
    }
    if (method != 1 && method != 3) {
        throw std::invalid_argument("Unsupported nearest-edge mapping method.");
    }

    SpatialGrid grid(grid_size);
    grid.set_points(&query_points);
    grid.build_grid_mt(num_threads);
    const auto centers = grid.centers();
    std::vector<int> results(grid.point_count(), -1);
    const std::size_t center_count = centers.size();
    const int worker_count = effective_thread_count(num_threads, center_count, 64);
    std::atomic<std::size_t> next_center{0};

    auto process_centers = [&]() {
        while (true) {
            const std::size_t center_index = next_center.fetch_add(1, std::memory_order_relaxed);
            if (center_index >= center_count) {
                return;
            }
            const auto& center = centers[center_index];
            const auto& point_indices = grid.cell(center.x, center.y);
            if (point_indices.empty()) {
                continue;
            }

            const auto nearest = index_.nearest(center);
            if (method == 1) {
                for (const int point_index : point_indices) {
                    results[static_cast<std::size_t>(point_index)] = nearest.edge_id;
                }
                continue;
            }

            const auto candidates = index_.within_nearest_distance(
                center, grid.epsilon() * 2.0, 200
            );
            const double initial_limit = std::sqrt(nearest.distance_squared) + grid.epsilon();
            for (const int point_index : point_indices) {
                const auto& point = query_points[static_cast<std::size_t>(point_index)];
                int best_edge = nearest.edge_id;
                double best_distance = initial_limit;
                for (const int edge_id : candidates) {
                    const double distance = std::sqrt(
                        polylines_[static_cast<std::size_t>(edge_id)]
                            .project(point).distance_squared
                    );
                    if (distance < best_distance ||
                        (distance == best_distance && edge_id < best_edge)) {
                        best_distance = distance;
                        best_edge = edge_id;
                    }
                }
                results[static_cast<std::size_t>(point_index)] = best_edge;
            }
        }
    };

    if (worker_count == 1) {
        process_centers();
    } else {
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(worker_count));
        try {
            for (int worker = 0; worker < worker_count; ++worker) {
                workers.emplace_back(process_centers);
            }
        } catch (...) {
            for (auto& worker : workers) {
                worker.join();
            }
            throw;
        }
        for (auto& worker : workers) {
            worker.join();
        }
    }
    return results;
}

double NNE::project_point_onto_edge(int edge_id, double x, double y) const {
    return edge(edge_id).project({x, y}).offset;
}

const Polyline& NNE::edge(int edge_id) const {
    if (edge_id < 0 || static_cast<std::size_t>(edge_id) >= polylines_.size()) {
        throw std::out_of_range("Edge ID is out of range.");
    }
    return polylines_[static_cast<std::size_t>(edge_id)];
}

}  // namespace fastnsa::nne

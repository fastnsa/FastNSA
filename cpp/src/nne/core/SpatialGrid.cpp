#include "nne/core/SpatialGrid.h"
#include "nne/core/Threading.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <thread>

namespace fastnsa::nne {
SpatialGrid::SpatialGrid(double grid_size) : width_(grid_size), height_(grid_size) {
    if (!std::isfinite(grid_size) || grid_size <= 0.0) {
        throw std::invalid_argument("Grid size must be finite and greater than zero.");
    }
}

std::size_t SpatialGrid::PairHash::operator()(const GridIndex& key) const noexcept {
    return std::hash<int>{}(key.first) ^ (std::hash<int>{}(key.second) << 1);
}

SpatialGrid::GridIndex SpatialGrid::index_for(double x, double y) const {
    return {static_cast<int>(std::floor(x / width_)),
            static_cast<int>(std::floor(y / height_))};
}

QueryPoint SpatialGrid::center_for(const GridIndex& index) const {
    return {index.first * width_ + width_ * 0.5,
            index.second * height_ + height_ * 0.5};
}

void SpatialGrid::set_points(const std::vector<QueryPoint>* points) {
    points_ = points;
}

void SpatialGrid::build_grid_mt(int num_threads) {
    if (!points_) {
        throw std::logic_error("Spatial grid points have not been initialized.");
    }
    const std::size_t count = points_->size();
    const int worker_count = effective_thread_count(num_threads, count, 512);
    std::vector<Grid> local_grids(static_cast<std::size_t>(worker_count));
    auto build_range = [&](int worker, std::size_t begin, std::size_t end) {
        auto& local_grid = local_grids[static_cast<std::size_t>(worker)];
        for (std::size_t index = begin; index < end; ++index) {
            const auto& point = (*points_)[index];
            local_grid[index_for(point.x, point.y)].push_back(static_cast<int>(index));
        }
    };

    if (worker_count == 1) {
        build_range(0, 0, count);
    } else {
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(worker_count));
        const std::size_t chunk_size = (count + worker_count - 1) /
            static_cast<std::size_t>(worker_count);
        try {
            for (int worker = 0; worker < worker_count; ++worker) {
                const std::size_t begin = static_cast<std::size_t>(worker) * chunk_size;
                const std::size_t end = std::min(count, begin + chunk_size);
                workers.emplace_back(build_range, worker, begin, end);
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

    grid_.clear();
    for (auto& local_grid : local_grids) {
        for (auto& entry : local_grid) {
            auto& destination = grid_[entry.first];
            destination.insert(destination.end(), entry.second.begin(), entry.second.end());
        }
    }
    point_count_ = points_->size();
}

std::vector<QueryPoint> SpatialGrid::centers() const {
    std::vector<QueryPoint> result;
    result.reserve(grid_.size());
    for (const auto& entry : grid_) {
        result.push_back(center_for(entry.first));
    }
    return result;
}

const std::vector<int>& SpatialGrid::cell(double x, double y) const {
    static const std::vector<int> empty;
    const auto iterator = grid_.find(index_for(x, y));
    return iterator == grid_.end() ? empty : iterator->second;
}

double SpatialGrid::epsilon() const noexcept {
    return std::hypot(width_, height_) * 0.5;
}

}  // namespace fastnsa::nne

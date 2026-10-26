#pragma once

#include "nne/geometry/QueryPoint.h"

#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fastnsa::nne {

class SpatialGrid {
public:
    using GridIndex = std::pair<int, int>;

    explicit SpatialGrid(double grid_size);

    void set_points(const std::vector<QueryPoint>* points);
    void build_grid_mt(int num_threads);
    std::vector<QueryPoint> centers() const;
    const std::vector<int>& cell(double x, double y) const;
    double epsilon() const noexcept;
    std::size_t point_count() const noexcept { return point_count_; }

private:
    struct PairHash {
        std::size_t operator()(const GridIndex& key) const noexcept;
    };

    using Grid = std::unordered_map<GridIndex, std::vector<int>, PairHash>;

    GridIndex index_for(double x, double y) const;
    QueryPoint center_for(const GridIndex& index) const;

    double width_;
    double height_;
    const std::vector<QueryPoint>* points_ = nullptr;
    Grid grid_;
    std::size_t point_count_ = 0;
};

}  // namespace fastnsa::nne

#pragma once

#include "nne/geometry/Polyline.h"

#include <cstddef>
#include <vector>

namespace fastnsa::nne {

struct EdgeMatch {
    int edge_id = -1;
    double distance_squared = 0.0;
    double offset = 0.0;
};

class PolylineBVH {
public:
    explicit PolylineBVH(const std::vector<Polyline>& polylines);

    EdgeMatch nearest(const Point2D& point) const;
    std::vector<int> within_nearest_distance(
        const Point2D& point,
        double extra_distance,
        std::size_t maximum_count
    ) const;

private:
    struct Node {
        Envelope bounds;
        int left = -1;
        int right = -1;
        std::size_t begin = 0;
        std::size_t end = 0;

        bool leaf() const noexcept { return left < 0; }
    };

    int build(std::size_t begin, std::size_t end);
    std::vector<int> candidates_within(
        const Point2D& point,
        double maximum_distance_squared
    ) const;

    const std::vector<Polyline>* polylines_;
    std::vector<int> edge_order_;
    std::vector<Node> nodes_;
};

}  // namespace fastnsa::nne

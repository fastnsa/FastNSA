#pragma once

#include "nne/geometry/Polyline.h"
#include "nne/geometry/QueryPoint.h"
#include "nne/index/PolylineBVH.h"

#include <vector>

namespace fastnsa::nne {

class NNE {
public:
    explicit NNE(std::vector<Polyline> polylines);

    std::vector<int> query_nearest_mt(
        const std::vector<QueryPoint>& query_points,
        int method,
        double grid_size,
        int num_threads
    ) const;

    double project_point_onto_edge(int edge_id, double x, double y) const;
    const Polyline& edge(int edge_id) const;

private:
    std::vector<int> query_exact_mt(
        const std::vector<QueryPoint>& query_points,
        int num_threads
    ) const;

    std::vector<Polyline> polylines_;
    PolylineBVH index_;
};

}  // namespace fastnsa::nne

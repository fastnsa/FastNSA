//
// Created by yehongwei on 2025/7/23.
//

#ifndef NNE_H
#define NNE_H
#include "nne/index/MySTRtree.h"
#include "nne/core/SpatialGrid.h"
#include "nne/geometry/QueryPoint.h"
#include "nne/geometry/distance.h"
// #include "Algorithm.h"


class NNE {
public:
explicit NNE(std::vector<std::unique_ptr<geos::geom::Geometry>> geoms);
explicit NNE(const std::string & path_to_geoms);

    // vector<int> query_nearest_str(std::string path_to_points, double grid_size);

    vector<int> query_nearest(const vector<QueryPoint> &query_points,int method, double grid_size, int max_k);

    vector<int> query_nearest_mt(
        const vector<QueryPoint> &query_points,
        int method,
        double grid_size,
        int num_threads); // deprecate gknn
    // vector<int> query_nearest_grid_enn(const vector<shared_ptr<geos::geom::Geometry>> &query_points, double grid_size, int max_k);

    vector<int> query_nearest_grid_enn(const vector<QueryPoint> &query_points, double grid_size); // default k = 200 //todo
    vector<int> query_nearest_grid(const vector<QueryPoint> &query_points, double grid_size);

private:
    std::vector<std::unique_ptr<geos::geom::Geometry>> edges;
    MySTRtree rtree;
    void build_rtree();
    static void build_grid(SpatialGrid &grid, const vector<QueryPoint>& query_points);
    static void build_grid_mt(SpatialGrid &grid, const vector<QueryPoint>& query_points);

    std::vector<int> query_nearest0(const std::vector<QueryPoint>& query_points);//TODO
    std::vector<int> query_nearest0_mt(const std::vector<QueryPoint>& query_points);

    std::vector<std::unique_ptr<geos::geom::Geometry>> read_geoms_wkt(const char *filename);
    void process_grid_center(const QueryPoint& query_point, int method,int max_k,SpatialGrid& grid,vector<int>& results);// todo
};



#endif //NNE_H

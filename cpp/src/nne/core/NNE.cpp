//
// Created by yehongwei on 2025/7/23.
//

#include "nne/core/NNE.h"
#include <geos/io/WKTReader.h>
#include <geos/linearref/LengthIndexedLine.h>
using namespace geos::geom;
using namespace geos::io;

NNE::NNE(std::vector<std::unique_ptr<geos::geom::Geometry>> geoms) {
    edges = std::move(geoms);
    build_rtree();
}

NNE::NNE(const std::string & path_to_geoms) {
    std::cout << "Loading geometries from: " << path_to_geoms << std::endl;
    edges = read_geoms_wkt(path_to_geoms.c_str());
    
    build_rtree(); 
}

// vector<int> NNE::query_nearest_str(std::string path_to_points, double grid_size)
// {
//     auto points = read_geoms_wkt(path_to_points.c_str());
//     return query_nearest_grid_enn(points, grid_size);
// }

void NNE::build_rtree() {
    if (edges.empty()) {
        throw std::invalid_argument("build_rtree: empty input geometries");
    }

    int idx = 0;
    for(const auto& geometry : edges) {
        assert(geometry);
        geometry->setSRID(idx);
        rtree.insert(geometry->getEnvelopeInternal(), geometry.get()); //test 改为存index
        ++idx;
    }
    rtree.build();
}

// test : optimize build_grid
void NNE::build_grid(SpatialGrid &grid, const vector<QueryPoint>& query_points) {
    for (const auto& point : query_points) {
        grid.insertPoint(point.x, point.y, grid.getPointsNum()); 
    }
    grid.set_points(&query_points);

    std::cerr<<"  LOG: cells number: "<<grid.getGridNums()<<std::endl;
}

void NNE::build_grid_mt(SpatialGrid &grid, const vector<QueryPoint> &query_points)
{
    grid.set_points(&query_points);
    grid.build_grid_mt();
}

vector<int> NNE::query_nearest0( const vector<QueryPoint> &query_points) {
    vector<int> results;
    results.reserve(query_points.size());

    for (const auto& qp : query_points) {

        const auto* result = static_cast<const Geometry*>(rtree.nearestNeighbour(&qp.env, &qp));

        results.emplace_back(result->getSRID());
    }
    return results;
}
vector<int> NNE::query_nearest0_mt(const vector<QueryPoint>& query_points) {
    const int n = static_cast<int>(query_points.size());
    vector<int> results(n);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        const auto* result = static_cast<const Geometry*>(
            rtree.nearestNeighbour(&query_points[i].env,&query_points[i])
        );
        results[i] = result->getSRID();
    }

    return results;
}


// TODO 优化grid查询部分
void NNE::process_grid_center(
    const QueryPoint& query_point,
    int method,
    int max_k,
    SpatialGrid& grid,
    vector<int>& results)
{
    if (method == 1) {
        // 单点最近邻（grid:center → edge）
        const auto* result = static_cast<const Geometry*>(
            rtree.nearestNeighbour(&query_point.env, &query_point)
        );
        auto idxs = grid.getCell(query_point.x, query_point.y);

        for (int idx : idxs) {
            results[idx] = result->getSRID();
        }
    }
    else {
        // KNN 或 epsilon NN
        vector<const void*> knn_raw;
       
        if (method == 2)
            knn_raw = rtree.kNearestNeighbours(&query_point.env, &query_point, max_k);
        else
            knn_raw = rtree.epsilonNearestNeighbours(&query_point.env, &query_point, max_k, grid.getGridEpsilon());

        vector<const LineString*> knn_edges;
        knn_edges.reserve(knn_raw.size());
        for (auto raw : knn_raw) {
                const auto* geom = static_cast<const Geometry*>(raw);
                knn_edges.push_back(static_cast<const LineString*>(geom));
        }

        auto idxs = grid.getCell(query_point.x, query_point.y); // MARK grid.getcell

        const auto* nearest_edge = knn_edges.front();
        // double lowerbound = query_point->distance(nearest_edge) + grid.getGridEpsilon(); // MARK grid.getgridepsilon
        

        double lowerbound = lineDistanceToQueryPoint(nearest_edge, query_point) + grid.getGridEpsilon(); // MARK grid.getgridepsilon


        for (int idx : idxs) {
            double minDist = lowerbound;
            auto* bestEdge = nearest_edge;

            for (auto* edge : knn_edges) {
                double d = grid.pointDistance(idx, edge); //Mark grid.pointDistance [需要原始点集，从而根据idx计算相应的点线距离]
                if (d < minDist) {
                    minDist = d;
                    bestEdge = edge;
                }
            }

            results[idx] = bestEdge->getSRID();
        }
    }
}

vector<QueryPoint> NNE::transfrom_coords_to_querypoints(double *xs, double *ys, int n_points)
{
    std::vector<QueryPoint> points;
    points.reserve(n_points);

    for (size_t i = 0; i < n_points; ++i) {
        points.push_back(QueryPoint{
            xs[i],
            ys[i]
        });
    }

    return points;
}

// test: query_nearest_mt
 vector<int> NNE::query_nearest_mt(const vector<QueryPoint> &query_points,
     const int method, const double grid_size, const int num_threads) {
        if(method == 0) return query_nearest0_mt(query_points);

        // cerr << "  LOG: start querying " << std::endl;

        // auto t_start = std::chrono::high_resolution_clock::now();

        auto grid = SpatialGrid(grid_size);
        build_grid_mt(grid,query_points);

        // auto t_end = std::chrono::high_resolution_clock::now();
        // auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
        // std::cout << "build grid elapsed: " << elapsed_ms << " ms" << std::endl;

        // cerr << "points num:" << grid.getPointsNum()<<endl;
        // cerr << "grid num:" << grid.getGridNums()<<endl;

        auto centers  = grid.getAllGridCenters();

        //这种情况下query_points是所有网格的中心点，需要用到SpatialGrid
        vector<int> results(grid.getPointsNum());

        int max_k = 100;

        // cout << "start processing  centers size():" << centers.size() <<endl;

        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (int i = 0; i < centers.size(); ++i) {
            process_grid_center(centers[i], method, max_k, grid, results);
        }

        return results;

    }



//test: query_nearest
 vector<int> NNE::query_nearest(const vector<QueryPoint> &query_points,
     const int method, const double grid_size, const int max_k) {

    if(method == 0) return query_nearest0(query_points);

    // cerr << "  LOG: start querying " << std::endl;

    // auto t_start = std::chrono::high_resolution_clock::now();

    auto grid = SpatialGrid(grid_size);
    build_grid(grid,query_points);

    // auto t_end = std::chrono::high_resolution_clock::now();
    // auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
    // std::cout << "build grid elapsed: " << elapsed_ms << " ms" << std::endl;

    // cerr << "points num:" << grid.getPointsNum()<<endl;
    // cerr << "grid num:" << grid.getGridNums()<<endl;

    auto query_centers = grid.getAllGridCenters();

    //这种情况下query_points是所有网格的中心点，需要用到SpatialGrid
    vector<int> results;
    results.resize(grid.getPointsNum());

    for(const auto& query_point : query_centers) {
        process_grid_center(query_point, method, max_k, grid, results);
    }

    return results;
}

// vector<int> NNE::query_nearest_grid_enn(const vector<shared_ptr<geos::geom::Geometry>> &query_points, const double grid_size,
//     const int max_k) {
//     return query_nearest(query_points, 3,grid_size,max_k);
// }

vector<int> NNE::query_nearest_grid_enn(const  vector<QueryPoint> &query_points, const double grid_size) {
    return query_nearest(query_points, 3,grid_size,200);
}

vector<int> NNE::query_nearest_grid(const  vector<QueryPoint> &query_points, const double grid_size) {

    return query_nearest(query_points, 1, grid_size,0);
}

// todo test
double NNE::project_point_onto_edge(int edge_id, double x, double y) {
    if (edge_id < 0 || edge_id >= edges.size()) {
        throw std::out_of_range("Edge ID is out of range.");
    }

    const auto& geom = edges[edge_id];
    geos::linearref::LengthIndexedLine lil(geom.get());
    double offset = lil.project(geos::geom::Coordinate(x, y));

    offset = std::clamp(offset, 0.0, geom->getLength());

    return offset;
}


std::vector<std::unique_ptr<geos::geom::Geometry>> NNE::read_geoms_wkt(const char *filename) {
    std::vector<std::unique_ptr<Geometry>> geoms;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件！" << std::endl;
        return geoms;
    }

    std::string line;

    while (std::getline(file, line)) {
        // 去除行首尾的空白字符
        line.erase(std::remove(line.begin(), line.end(), '"'), line.end()); // 删除双引号
        line.erase(0, line.find_first_not_of(" \t\n\r\f\v"));
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);

        if (line.empty()) {continue;}

        try {
            WKTReader reader;
            geoms.emplace_back(reader.read(line));

        } catch (const geos::util::GEOSException& e) {
            std::cerr << "解析错误: " << e.what() << "\n行内容: " << line << std::endl;
        }
    }
    file.close();

    return geoms;
}
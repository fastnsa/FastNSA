//
// Created by yehongwei on 2025/3/26.
//

#ifndef SPATIALGRID_H
#define SPATIALGRID_H

#include <geos/geom/LineString.h>
#include <geos/geom/Point.h>

#include <geos/geom/GeometryFactory.h>

#include <omp.h> // multi thread
#include "nne/geometry/QueryPoint.h"
#include "nne/geometry/distance.h"

#include "Library.h"


// using namespace geos::geom;

class SpatialGrid {
    using Geometry = geos::geom::Geometry;
    using LineString = geos::geom::LineString;
    using GridIndex = std::pair<int,int>;

private:
    double grid_w, grid_h;
    const std::vector<QueryPoint>* points = nullptr;

    struct PairHash {
        std::size_t operator()(const std::pair<int, int>& key) const {
            // 混合哈希值以减少冲突
            return std::hash<int>{}(key.first) ^ (std::hash<int>{}(key.second) << 1);
        }
    };

    using LocalGrid = std::unordered_map<GridIndex, std::vector<int>, PairHash>;

    unordered_map<pair<int, int>, vector<int>, PairHash> grid; // key: grid index, value: point indices in this grid
    size_t grid_point_count = 0;


public:
    // void insertPoint(unique_ptr<Geometry> p); 

    explicit SpatialGrid(const double g_size):grid_w(g_size),grid_h(g_size){}
    SpatialGrid(const double grid_w,const double grid_h):grid_w(grid_w),grid_h(grid_h){}

    void insertPoint(double coord_x, double coord_y, int point_index); 

    void set_points(const std::vector<QueryPoint>* pts);

    void build_grid_mt();

    /**
     * 根据网格的key获取某网格的中心点坐标
     * @param key 网格的key，由一个pair组成，key.first表示网格的第一个编号，key.second表示网格的第二个编号
     * @return 该网格的中心点坐标
     */
    pair<double, double> getGridCenter(const pair<int, int> &key) const;

    /**
     *获取输入位置所处的网格内的所有点。
     * @param coord
     * @return 输入点所处的cell内所有点index的vector
     */
    vector<int> getCell(double x, double y);


    vector<QueryPoint> getAllGridCenters() const;
    size_t getGridNums()const;
    size_t getPointsNum()const; 
    // const vector<unique_ptr<Geometry>>& getPoints();


    double getGridWidth()const;
    double getGridHeight()const;

    double getGridEpsilon()const;
    // double getGridSize()const;

    vector<int>& operator[](const pair<int,int> &index);

    // Point* getPointAt(int i) const;
    double pointDistance(int i, const LineString *&edge) const;


};


#endif //SPATIALGRID_H

//
// Created by yehongwei on 2025/3/26.
//

#include "nne/core/SpatialGrid.h"

using namespace geos::geom;


void SpatialGrid::insertPoint(double coord_x, double coord_y, int point_index) 
{
    int x = floor(coord_x/grid_w);
    int y = floor(coord_y/grid_h);

    pair index = {x, y};
    grid[index].emplace_back(point_index); //record point index // test this
    grid_point_count++;
}

void SpatialGrid::set_points(const std::vector<QueryPoint> *pts)
{
    points = pts;
}

void SpatialGrid::build_grid_mt()
{
    int n = (*points).size();

    int num_threads = omp_get_max_threads();  
    std::vector<LocalGrid> local_grids(num_threads);

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& mygrid = local_grids[tid];

        #pragma omp for schedule(static)
        for (int i = 0; i < n; i++) {
            auto& geom = (*points)[i];

            int gx = static_cast<int>(floor(geom.x / grid_w));
            int gy = static_cast<int>(floor(geom.y / grid_h));

            mygrid[{gx, gy}].push_back(i);
        }
    }

    // merge local grids
    for (auto& lg : local_grids) {
        for (auto& kv : lg) {
            auto& cell = grid[kv.first];
            cell.insert(cell.end(), kv.second.begin(), kv.second.end());
        }
    }

    grid_point_count = n;
}

size_t SpatialGrid::getGridNums() const {
    return grid.size();
}


pair<double, double> SpatialGrid::getGridCenter(const pair<int, int> &key) const {
    double x = key.first * grid_w + grid_w/2;
    double y = key.second * grid_h + grid_h/2;

    return {x,y};
}


vector<int> SpatialGrid::getCell(double x, double y){


    int col = floor(x/grid_w);
    int row = floor(y/grid_h);

    return grid[{col,row}];
}

vector<QueryPoint> SpatialGrid::getAllGridCenters() const { //todo getAllGridCenters optimize
    vector<QueryPoint> centers;
    for(const auto& entry:grid) {
        auto center = getGridCenter(entry.first);
        centers.emplace_back(QueryPoint(center.first,center.second));
    }
    return centers;
}

double SpatialGrid::getGridWidth() const {
    return grid_w;

}

double SpatialGrid::getGridHeight() const {
    return grid_h;
}

double SpatialGrid::getGridEpsilon() const {
    return sqrt(grid_w*grid_w+grid_h*grid_h)/2;
}

vector<int> & SpatialGrid::operator[](const pair<int, int> &index) {
    return grid[index];
}

size_t SpatialGrid::getPointsNum() const { //done getPointsNum [调整内部逻辑]
    return grid_point_count;
}


// TODO pointsDistance [参数、逻辑]
double SpatialGrid::pointDistance(int i,const LineString*& edge) const {
    if(i>=grid_point_count || i<0) {
        cerr<<i <<" is out of bound";
        return 0;
    }
    return lineDistanceToQueryPoint(edge, (*points)[i]);
}





//
// Created by yehongwei on 2025/4/3.
//

#ifndef POINTLINEDISTANCE_H
#define POINTLINEDISTANCE_H
#include <geos/geom/GeometryFactory.h>
#include <geos/index/strtree/STRtree.h>
#include "nne/geometry/QueryPoint.h"
#include "nne/geometry/distance.h"



class PointLineDistance : public geos::index::strtree::ItemDistance{
  public:
    double distance(const geos::index::strtree::ItemBoundable* item1, const geos::index::strtree::ItemBoundable* item2) override ;
};



#endif //POINTLINEDISTANCE_H

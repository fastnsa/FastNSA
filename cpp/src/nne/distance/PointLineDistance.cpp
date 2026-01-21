//
// Created by yehongwei on 2025/4/3.
//

#include "nne/distance/PointLineDistance.h"



using namespace geos::index::strtree;
using namespace geos::geom;




double PointLineDistance::distance( //test 重构distance
    const ItemBoundable* a,
    const ItemBoundable* b
) {
    const void* o1 = a->getItem();
    const void* o2 = b->getItem();

    auto line = static_cast<const LineString*>(o1); // test this
    auto qp = static_cast<const QueryPoint*>(o2);

    return lineDistanceToQueryPoint(
        line, *qp
    );

    // if (auto qp = static_cast<const QueryPoint*>(o1)) {
    //     if (auto geom = static_cast<const Geometry*>(o2)) {
    //         return lineDistanceToQueryPoint(
    //             dynamic_cast<const LineString*>(geom), *qp
    //         );
    //     }
    // }

    // if (auto qp = static_cast<const QueryPoint*>(o2)) {
    //     if (auto geom = static_cast<const Geometry*>(o1)) {
    //         return lineDistanceToQueryPoint(
    //             dynamic_cast<const LineString*>(geom), *qp
    //         );
    //     }
    // }

    // 兜底：Node–Node / Point–Point / 剪枝

    // const auto* e1 =
    //     static_cast<const geos::geom::Envelope*>(a->getBounds());
    // const auto* e2 =
    //     static_cast<const geos::geom::Envelope*>(b->getBounds());

    // return e1->distance(*e2);
}

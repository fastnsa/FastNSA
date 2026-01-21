//
// Created by yehongwei on 2025/4/3.
//

#ifndef MYSTRTREE_H
#define MYSTRTREE_H

#include <geos/index/strtree/STRtree.h>
#include <geos/index/strtree/AbstractNode.h>
#include <geos/index/strtree/Boundable.h>


// #include <geos/index/strtree/BoundablePair.h>
#include "nne/core/BoundablePair.h"
#include "nne/distance/PointLineDistance.h"


class MySTRtree : public geos::index::strtree::STRtree {
public:
  // using STRtree::STRtree;
  using Envelope = geos::geom::Envelope;

  std::vector<const void *> kNearestNeighbours(const Envelope *env, const void *item, int k);
  const void* nearestNeighbour(const Envelope *env, const void *item) ;

  std::vector<const void *> epsilonNearestNeighbours(const Envelope *env, const void *item, int k, double epsilon);

  std::vector<const void *> epsilonNearestNeighbours(const Envelope *env, const void *item, double epsilon);

  std::vector<const void *> kNearestNeighbours(const Envelope *env, const void *item, int k,double maxDistance);

private:
  PointLineDistance dist_ = PointLineDistance();

};



#endif //MYSTRTREE_H

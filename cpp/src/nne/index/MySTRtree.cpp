//
// Created by yehongwei on 2025/4/3.
//

#include "nne/index/MySTRtree.h"

using namespace geos::geom;
using namespace geos::index::strtree;

std::vector<const void *> MySTRtree::kNearestNeighbours(const Envelope *env, const void *item, int k) {

    return kNearestNeighbours(env,item,k,geos::DoubleInfinity);
}


std::vector<const void *> MySTRtree::kNearestNeighbours(const Envelope *env, const void *item, int k, double maxDistance) {//TODO
    build();
    ItemBoundable bnd = ItemBoundable(env,const_cast<void*>(item));
    nne::BoundablePair bp(getRoot(),&bnd,&dist_);
    nne::BoundablePair* initBndPair = &bp;

    double distanceLowerBound = maxDistance;

    nne::BoundablePair::BoundablePairQueue priQ;
    priQ.push(initBndPair);

    std::vector <const void*> results;


   while(!priQ.empty() && distanceLowerBound > 0.0) {
        nne::BoundablePair* bndPair = priQ.top();
        double currentDistance = bndPair->getDistance();

       if(results.size()==k) {
           break;
       }

        priQ.pop();

        if(bndPair->isLeaves()) {
            // distanceLowerBound = currentDistance;
            results.push_back(dynamic_cast<const ItemBoundable*>(bndPair->getBoundable(0))->getItem());
        }
        else {
            bndPair->expandToQueue(priQ, distanceLowerBound);
        }

        if(bndPair != initBndPair) {
            delete bndPair;
        }
    }

    /* Free any remaining BoundablePairs in the queue */
    while(!priQ.empty()) {
        nne::BoundablePair* bndPair = priQ.top();
        priQ.pop();
        if(bndPair != initBndPair) {
            delete bndPair;
        }
    }

    if(results.size() != k) {
         throw geos::util::GEOSException("Error computing nearest neighbor");
    }

    return results;
}

const void * MySTRtree::nearestNeighbour(const Envelope *env, const void *item) {
    return STRtree::nearestNeighbour(env,item,&dist_);

}

std::vector<const void *> MySTRtree::epsilonNearestNeighbours(const Envelope *env, const void *item, const int k, const double epsilon) {
    build();
    ItemBoundable bnd = ItemBoundable(env,const_cast<void*>(item));
    nne::BoundablePair bp(getRoot(),&bnd,&dist_);
    nne::BoundablePair* initBndPair = &bp;

    double distanceLowerBound = geos::DoubleInfinity;

    nne::BoundablePair::BoundablePairQueue priQ;
    priQ.push(initBndPair);

    std::vector <const void*> results;

    while(!priQ.empty() && distanceLowerBound > 0.0) {
        nne::BoundablePair* bndPair = priQ.top();
        double currentDistance = bndPair->getDistance();
        if(currentDistance >= distanceLowerBound || results.size() == k) {
            break;
        }

        priQ.pop();

        if(bndPair->isLeaves()) {
            // distanceLowerBound = currentDistance;
            //todo update lowerBound
            //首次找到的叶节点为nearest neighbor, 基于此距离设置2epsilon 阈值，过滤其它节点

            if(results.empty()) {
                distanceLowerBound = currentDistance + 2 * epsilon;
            }
            results.push_back(dynamic_cast<const ItemBoundable*>(bndPair->getBoundable(0))->getItem());
        }
        else {
            bndPair->expandToQueue(priQ, distanceLowerBound);
        }
        if(bndPair != initBndPair) {
            delete bndPair;
        }
    }

    /* Free any remaining BoundablePairs in the queue */
    while(!priQ.empty()) {
        nne::BoundablePair* bndPair = priQ.top();
        priQ.pop();
        if(bndPair != initBndPair) {
            delete bndPair;
        }
    }

    if(results.empty()) {
        throw geos::util::GEOSException("Error computing nearest neighbor");
    }

    return results;
}

std::vector<const void *> MySTRtree::epsilonNearestNeighbours(const Envelope *env, const void *item, const double epsilon) {
    return epsilonNearestNeighbours(env,item,200,epsilon);
}

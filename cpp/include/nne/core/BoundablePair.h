
#ifndef BOUNDABLEPAIR_H
#define BOUNDABLEPAIR_H

#include <geos/index/strtree/Boundable.h>
#include <geos/index/strtree/ItemDistance.h>
#include <queue>

namespace nne
{
    
class BoundablePair {
private:
    const geos::index::strtree::Boundable* boundable1;
    const geos::index::strtree::Boundable* boundable2;
    geos::index::strtree::ItemDistance* itemDistance;
    double mDistance;

public:
    struct BoundablePairQueueCompare {
        bool operator()(const BoundablePair* a, const BoundablePair* b) {
            return a->getDistance() > b->getDistance();
        }
    };

    using BoundablePairQueue = std::priority_queue<BoundablePair*, std::vector<BoundablePair*>, BoundablePairQueueCompare>;

    BoundablePair(const geos::index::strtree::Boundable* b1,
                  const geos::index::strtree::Boundable* b2,
                  geos::index::strtree::ItemDistance* dist);

    double getDistance() const;
    bool isLeaves() const;
    const geos::index::strtree::Boundable* getBoundable(int i) const;
    void expandToQueue(BoundablePairQueue& queue, double minDistance);
};

} // namespace nne

#endif
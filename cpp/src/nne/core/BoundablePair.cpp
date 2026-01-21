#include "nne/core/BoundablePair.h"

#include <geos/index/strtree/AbstractNode.h>
#include <geos/index/strtree/ItemBoundable.h>
#include <geos/geom/Envelope.h>
#include <geos/index/strtree/EnvelopeUtil.h>
#include <geos/util.h>
#include <geos/util/IllegalArgumentException.h>
#include <memory>

using namespace geos::index::strtree;

namespace nne{

BoundablePair::BoundablePair(const Boundable* b1, const Boundable* b2, ItemDistance* dist)
    : boundable1(b1), boundable2(b2), itemDistance(dist) {
    if (isLeaves()) {
        mDistance = itemDistance->distance(
            geos::detail::down_cast<const ItemBoundable*>(boundable1),
            geos::detail::down_cast<const ItemBoundable*>(boundable2));
    } else {
        const auto* e1 = static_cast<const geos::geom::Envelope*>(boundable1->getBounds());
        const auto* e2 = static_cast<const geos::geom::Envelope*>(boundable2->getBounds());
        if (!e1 || !e2) {
            throw geos::util::GEOSException("Cannot compute envelope of BoundablePair");
        }
        mDistance = e1->distance(*e2);
    }
}

double BoundablePair::getDistance() const {
    return mDistance;
}

bool BoundablePair::isLeaves() const {
    return boundable1->isLeaf() && boundable2->isLeaf();
}

const Boundable* BoundablePair::getBoundable(int i) const {
    return (i == 0) ? boundable1 : boundable2;
}

void BoundablePair::expandToQueue(BoundablePairQueue& queue, double minDistance) {
    bool isComp1 = !boundable1->isLeaf();
    bool isComp2 = !boundable2->isLeaf();

    auto expand = [&](const Boundable* composite, const Boundable* other, bool flip) {
        const auto& children = geos::detail::down_cast<const AbstractNode*>(composite)->getChildBoundables();
        for (Boundable* child : *children) {
            std::unique_ptr<BoundablePair> bp;
            if (flip) {
                bp.reset(new BoundablePair(other, child, itemDistance));
            } else {
                bp.reset(new BoundablePair(child, other, itemDistance));
            }
            if (minDistance == geos::DoubleInfinity || bp->getDistance() < minDistance) {
                queue.push(bp.release());
            }
        }
    };

    if (isComp1 && isComp2) {
        expand(boundable1, boundable2, false);
        expand(boundable2, boundable1, true);
    } else if (isComp1) {
        expand(boundable1, boundable2, false);
    } else if (isComp2) {
        expand(boundable2, boundable1, true);
    } else {
        throw geos::util::IllegalArgumentException("Neither boundable is composite");
    }
}


}
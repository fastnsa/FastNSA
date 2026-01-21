#include "nne/geometry/distance.h"

static inline double pointSegmentDistanceSquared(
    double px, double py,
    double x1, double y1,
    double x2, double y2
) {
    double dx = x2 - x1;
    double dy = y2 - y1;

    if (dx == 0.0 && dy == 0.0) {
        dx = px - x1;
        dy = py - y1;
        return std::sqrt(dx * dx + dy * dy);
    }

    double t = ((px - x1) * dx + (py - y1) * dy)
             / (dx * dx + dy * dy);

    t = std::max(0.0, std::min(1.0, t));

    double projx = x1 + t * dx;
    double projy = y1 + t * dy;

    dx = px - projx;
    dy = py - projy;
    return std::sqrt(dx * dx + dy * dy);
}


double lineDistanceToQueryPoint(
    const geos::geom::LineString* line,
    const QueryPoint& qp
) {
    if (!line) return std::numeric_limits<double>::infinity();

    const auto* coords = line->getCoordinatesRO();
    if (!coords) return std::numeric_limits<double>::infinity();


    if (!line) {
        std::cerr << "line == nullptr" << std::endl;
        std::abort();
    }

    const std::size_t n = coords->size();
    if (n < 2) return std::numeric_limits<double>::infinity();

    double minDist2 = std::numeric_limits<double>::infinity();

    for (std::size_t i = 1; i < n; ++i) {
        const auto& a = coords->getAt(i - 1);
        const auto& b = coords->getAt(i);

        const double d2 = pointSegmentDistanceSquared(
            qp.x, qp.y,
            a.x, a.y,
            b.x, b.y
        );

        minDist2 = std::min(minDist2, d2);
    }

    return minDist2;
}

// double lineDistanceToQueryPoint(
//     const geos::geom::LineString* line,
//     const QueryPoint& qp
// ) {
//     const auto* coords = line->getCoordinatesRO();
//     const std::size_t n = coords->size();

//     double minDist2 = std::numeric_limits<double>::infinity();

//     for (std::size_t i = 1; i < n; ++i) {
//         const auto& a = coords->getAt(i - 1);
//         const auto& b = coords->getAt(i);

//         const double d2 = pointSegmentDistanceSquared(
//             qp.x, qp.y,
//             a.x, a.y,
//             b.x, b.y
//         );

//         minDist2 = std::min(minDist2, d2);
//     }

//     return std::sqrt(minDist2);
// }

// double lineDistanceToQueryPoint(
//     const geos::geom::LineString* line,
//     const QueryPoint& qp
// ) {

//     geos::geom::Coordinate c(qp.x,qp.y);

//     return line->distance(&c);
// }

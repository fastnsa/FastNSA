#pragma once
#include <geos/geom/Envelope.h>

struct QueryPoint {
    double x;
    double y;

    // test queryPoint
    geos::geom::Envelope env;

    QueryPoint(double x_, double y_)
        : x(x_), y(y_), env(x_, x_, y_, y_) {}
};

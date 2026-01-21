#pragma once
#include <geos/geom/LineString.h>
#include "nne/geometry/QueryPoint.h"

double lineDistanceToQueryPoint(
    const geos::geom::LineString* line,
    const QueryPoint& qp
);

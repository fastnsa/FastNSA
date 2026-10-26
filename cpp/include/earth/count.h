#pragma once
#ifndef FASTNSA_EARTH_COUNT_H
#define FASTNSA_EARTH_COUNT_H

#include "earth/Network.h"

namespace fastnsa::earth{

double point_edge_count_same_edge(model& our_model, int p_index, int e_index);
double point_edge_count(model& our_model, int e_index, int node_u, int node_v);
void point_edge_count_same_edge_MTS(model& our_model, int e_index, int node_u, int node_v);
void point_edge_count_MTS(model& our_model, int e_index, int node_u, int node_v);

}
#endif

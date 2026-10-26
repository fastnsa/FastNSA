
#pragma once
#include "stars/Network.h"

namespace fastnsa::stars{

void sam_edge_edge_count_same_edge(model& our_model, int e_index);

void sam_edge_edge_count2(model& our_model, int wide_e_index, int e_index); //ns+sampling
void sam_edge_edge_count3(model& our_model, int wide_e_index, int e_index);

}

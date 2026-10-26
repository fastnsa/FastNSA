#pragma once


#include "kdv/Network.h"

namespace fastnsa::kdv {

double kernel_value(model& our_model, double dist);
double edge_KAF(model& our_model, int edge_index);
void correct_same_edge_shortcuts(model& our_model, Lixel& lixel, double to_n1, double to_n2);
void NKDV_basic(model& our_model);

//Used in aggregate distance augmentation (ADA)
void augment_preprocess(model& our_model);
//Used in interval augmentation (IA)
void augment_interval_preprocess(model& our_model);

} // namespace fastnsa::kdv


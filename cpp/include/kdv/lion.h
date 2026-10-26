#pragma once

#include "kdv/shortest_path.h"
#include "kdv/KAF.h"

namespace fastnsa::kdv {

void init_LION(model& our_model);
void lixel_augmentation(model& our_model);
void lixel_aggregation(model& our_model);
void init_one_D_KDV(model& our_model);
void one_D_KDV(model& our_model);


}


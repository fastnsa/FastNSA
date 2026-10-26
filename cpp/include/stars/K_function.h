#pragma once

#include "stars/Network.h"
#include "stars/SP.h"
#include "stars/count.h"
#include "stars/sam_count.h"

namespace fastnsa::stars{

void call_RQS(model& our_model);
void call_SPS(model& our_model);
void call_ASPS(model& our_model);
void call_STARS(model& our_model);
void compute_K_function(model& our_model);

}

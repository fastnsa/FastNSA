#pragma once
#ifndef FASTNSA_EARTH_K_FUNCTION_H
#define FASTNSA_EARTH_K_FUNCTION_H

#include "earth/Network.h"
#include "earth/SP.h"
#include "earth/count.h"

namespace fastnsa::earth{
void compute_K_function(model& our_model);
}
#endif

#pragma once

#include "kdv/Network.h"

#include <string>
#include <vector>
#include <functional>

namespace fastnsa::kdv {

enum class Method : int {
    BASIC = 1,
    SHARED_SPD = 2,
    ADA = 3,
    IA = 4,
    ADA_IA = 5,
    LION = 6,
};

enum class Kernel : int {
    GAUSSIAN = 0,
    TRIANGULAR = 1,
    EPANECHNIKOV = 2,
    QUARTIC = 3,
};

struct Result {
    std::vector<int> edge_ids;
    std::vector<double> offsets_n1;
    std::vector<double> offsets_n2;
    std::vector<double> density;
};

using ProgressCallback = std::function<bool(const char*, int, int)>;

int count_lixels(const std::vector<double>& edge_lengths, double lixel_length);

bool compute_from_arrays(
    int n_nodes,
    const std::vector<int>& edge_u,
    const std::vector<int>& edge_v,
    const std::vector<double>& edge_lengths,
    const std::vector<int>& point_edge_ids,
    const std::vector<double>& point_offsets,
    const std::vector<double>& point_weights,
    Method method,
    double lixel_length,
    Kernel kernel,
    double bandwidth,
    Result& out,
    std::string* error = nullptr,
    ProgressCallback progress = {}
);

}  // namespace fastnsa::kdv

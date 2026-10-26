#include "test_support.h"

#include <algorithm>

namespace {

void test_network_k_function_contract() {
    auto network = fastnsa::test::make_two_edge_network();
    FASTNSA_CHECK(network.get() != nullptr);
    auto points = fastnsa::test::PointSetHandle(fastnsa_pointset_create(network.get()));
    FASTNSA_CHECK(points.get() != nullptr);

    const std::vector<int> edge_ids{0, 0, 1, 1};
    const std::vector<double> offsets{2.0, 12.0, 2.0, 8.0};
    FASTNSA_CHECK(fastnsa_pointset_set(
        points.get(), edge_ids.data(), offsets.data(), static_cast<int>(edge_ids.size())
    ) == 0);

    const std::vector<double> radii{0.0, 5.0, 15.0, 30.0};
    std::vector<double> observed(radii.size(), -1.0);
    FASTNSA_CHECK(fastnsa_compute_k_function_checked(
        network.get(), points.get(), radii.data(), static_cast<int>(radii.size()),
        FASTNSA_K_STARS, 0, observed.data(), nullptr, 0.1
    ) == 0);
    FASTNSA_CHECK(std::all_of(observed.begin(), observed.end(), [](double value) {
        return std::isfinite(value) && value >= 0.0;
    }));
    FASTNSA_CHECK(std::is_sorted(observed.begin(), observed.end()));

    const std::vector<double> invalid_radii{5.0, 4.0};
    std::vector<double> invalid_output(invalid_radii.size());
    FASTNSA_CHECK(fastnsa_compute_k_function_checked(
        network.get(), points.get(), invalid_radii.data(),
        static_cast<int>(invalid_radii.size()), FASTNSA_K_STARS, 0,
        invalid_output.data(), nullptr, 0.1
    ) != 0);
    FASTNSA_CHECK(std::string(fastnsa_last_error()).find("threshold") != std::string::npos);
}

}  // namespace

int main() {
    return fastnsa::test::run("C API network K-function", [] {
        test_network_k_function_contract();
    });
}

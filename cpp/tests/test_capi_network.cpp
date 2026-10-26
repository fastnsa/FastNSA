#include "test_support.h"

#include <limits>

namespace {

void test_network_and_pointset_contracts() {
    FASTNSA_CHECK(fastnsa_abi_version() == 2);
    FASTNSA_CHECK(fastnsa_network_create(nullptr, 1, nullptr, nullptr, 1) == nullptr);
    FASTNSA_CHECK(std::string(fastnsa_last_error()).find("Invalid network") != std::string::npos);

    auto network = fastnsa::test::make_two_edge_network();
    FASTNSA_CHECK(network.get() != nullptr);

    const std::vector<double> coordinates{
        0.0, 0.0, 10.0, 0.0, 10.0, 10.0,
        10.0, 10.0, 20.0, 10.0,
    };
    const std::vector<int> geometry_offsets{0, 3, 5};
    FASTNSA_CHECK(fastnsa_network_set_polylines(
        network.get(), coordinates.data(), geometry_offsets.data(), 2
    ) == 0);

    auto pointset = fastnsa::test::PointSetHandle(fastnsa_pointset_create(network.get()));
    FASTNSA_CHECK(pointset.get() != nullptr);
    const std::vector<int> edge_ids{0, 1};
    const std::vector<double> offsets{7.5, 4.0};
    FASTNSA_CHECK(fastnsa_pointset_set(
        pointset.get(), edge_ids.data(), offsets.data(), 2
    ) == 0);
    FASTNSA_CHECK(fastnsa_pointset_size(pointset.get()) == 2);

    std::vector<int> copied_edges(2, -1);
    std::vector<double> copied_offsets(2, -1.0);
    std::vector<double> copied_times(2, -1.0);
    FASTNSA_CHECK(fastnsa_pointset_copy_edge_ids(
        pointset.get(), copied_edges.data(), 2
    ) == 0);
    FASTNSA_CHECK(fastnsa_pointset_copy_offsets(
        pointset.get(), copied_offsets.data(), 2
    ) == 0);
    FASTNSA_CHECK(copied_edges == edge_ids);
    FASTNSA_CHECK(fastnsa::test::close(copied_offsets[0], offsets[0]));
    FASTNSA_CHECK(fastnsa::test::close(copied_offsets[1], offsets[1]));

    const std::vector<double> times{100.0, 250.0};
    FASTNSA_CHECK(fastnsa_pointset_set_times(
        pointset.get(), times.data(), 2
    ) == 0);
    FASTNSA_CHECK(fastnsa_pointset_copy_times(
        pointset.get(), copied_times.data(), 2
    ) == 0);
    FASTNSA_CHECK(copied_times == times);

    const std::vector<double> invalid_offsets{21.0, 4.0};
    FASTNSA_CHECK(fastnsa_pointset_set(
        pointset.get(), edge_ids.data(), invalid_offsets.data(), 2
    ) == -2);
    const std::vector<double> invalid_times{
        1.0, std::numeric_limits<double>::quiet_NaN()
    };
    FASTNSA_CHECK(fastnsa_pointset_set_times(
        pointset.get(), invalid_times.data(), 2
    ) == -2);
}

}  // namespace

int main() {
    return fastnsa::test::run("C API network and point-set contracts", [] {
        test_network_and_pointset_contracts();
    });
}

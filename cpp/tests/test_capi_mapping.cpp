#include "test_support.h"

#include <cstdint>
#include <cstring>
#include <limits>

namespace {

template <typename T>
void append_scalar(std::vector<std::uint8_t>& destination, const T& value) {
    const auto* bytes = reinterpret_cast<const std::uint8_t*>(&value);
    destination.insert(destination.end(), bytes, bytes + sizeof(T));
}

std::vector<std::uint8_t> make_linestring_wkb(
    const std::vector<double>& interleaved_xy
) {
    std::vector<std::uint8_t> result;
    result.push_back(1);
    const std::uint32_t type = 2;
    const std::uint32_t point_count = static_cast<std::uint32_t>(interleaved_xy.size() / 2);
    append_scalar(result, type);
    append_scalar(result, point_count);
    for (const double coordinate : interleaved_xy) {
        append_scalar(result, coordinate);
    }
    return result;
}

struct MappingResult {
    std::vector<int> edge_ids;
    std::vector<double> offsets;
};

MappingResult map(void* network, int threads) {
    const std::vector<double> x{9.0, 16.0, 10.0};
    const std::vector<double> y{7.0, 12.0, 10.0};
    fastnsa::test::PointSetHandle points(fastnsa_map_points(
        network, x.data(), y.data(), static_cast<int>(x.size()),
        FASTNSA_PLNNS_RTREE, 10.0, threads
    ));
    FASTNSA_CHECK(points.get() != nullptr);
    MappingResult result{std::vector<int>(x.size()), std::vector<double>(x.size())};
    FASTNSA_CHECK(fastnsa_pointset_copy_edge_ids(
        points.get(), result.edge_ids.data(), static_cast<int>(x.size())
    ) == 0);
    FASTNSA_CHECK(fastnsa_pointset_copy_offsets(
        points.get(), result.offsets.data(), static_cast<int>(x.size())
    ) == 0);
    return result;
}

void check_expected(const MappingResult& result) {
    FASTNSA_CHECK(result.edge_ids == std::vector<int>({0, 1, 0}));
    FASTNSA_CHECK(fastnsa::test::close(result.offsets[0], 17.0));
    FASTNSA_CHECK(fastnsa::test::close(result.offsets[1], 6.0));
    FASTNSA_CHECK(fastnsa::test::close(result.offsets[2], 20.0));
}

void test_packed_wkb_and_threads() {
    auto packed_network = fastnsa::test::make_two_edge_network();
    FASTNSA_CHECK(packed_network.get() != nullptr);
    const std::vector<double> coordinates{
        0.0, 0.0, 10.0, 0.0, 10.0, 10.0,
        10.0, 10.0, 20.0, 10.0,
    };
    const std::vector<int> geometry_offsets{0, 3, 5};
    FASTNSA_CHECK(fastnsa_network_set_polylines(
        packed_network.get(), coordinates.data(), geometry_offsets.data(), 2
    ) == 0);

    const auto baseline = map(packed_network.get(), 1);
    check_expected(baseline);
    for (const int threads : {0, 4, 32, 512}) {
        const auto result = map(packed_network.get(), threads);
        FASTNSA_CHECK(result.edge_ids == baseline.edge_ids);
        for (std::size_t index = 0; index < result.offsets.size(); ++index) {
            FASTNSA_CHECK(fastnsa::test::close(
                result.offsets[index], baseline.offsets[index]
            ));
        }
    }

    auto wkb_network = fastnsa::test::make_two_edge_network();
    FASTNSA_CHECK(wkb_network.get() != nullptr);
    const auto first = make_linestring_wkb({0.0, 0.0, 10.0, 0.0, 10.0, 10.0});
    const auto second = make_linestring_wkb({10.0, 10.0, 20.0, 10.0});
    std::vector<std::uint8_t> wkbs = first;
    wkbs.insert(wkbs.end(), second.begin(), second.end());
    const std::vector<int> wkb_offsets{
        0, static_cast<int>(first.size()), static_cast<int>(wkbs.size())
    };
    FASTNSA_CHECK(fastnsa_network_set_wkbs_checked(
        wkb_network.get(), wkbs.data(), wkb_offsets.data(), 2
    ) == 0);
    const auto wkb_result = map(wkb_network.get(), 0);
    FASTNSA_CHECK(wkb_result.edge_ids == baseline.edge_ids);
    for (std::size_t index = 0; index < wkb_result.offsets.size(); ++index) {
        FASTNSA_CHECK(fastnsa::test::close(
            wkb_result.offsets[index], baseline.offsets[index]
        ));
    }

    const double invalid_x = std::numeric_limits<double>::quiet_NaN();
    const double valid_y = 0.0;
    FASTNSA_CHECK(fastnsa_map_points(
        packed_network.get(), &invalid_x, &valid_y, 1,
        FASTNSA_PLNNS_RTREE, 10.0, 0
    ) == nullptr);
}

void test_endpoint_roundoff_and_atomic_setters() {
    const double length = std::nextafter(100.0, 0.0);
    double edges[]{0.0, 1.0, length};
    int adjacency_offsets[]{0, 1, 2};
    int indices[]{0, 0};
    fastnsa::test::NetworkHandle network(fastnsa_network_create(
        edges, 1, adjacency_offsets, indices, 2
    ));
    FASTNSA_CHECK(network.get() != nullptr);
    const double xy[]{0.0, 0.0, 100.0, 0.0};
    const int geometry_offsets[]{0, 2};
    FASTNSA_CHECK(fastnsa_network_set_polylines(
        network.get(), xy, geometry_offsets, 1
    ) == 0);
    const double x[]{100.0, 0.0}, y[]{0.0, 0.0};
    fastnsa::test::PointSetHandle points(fastnsa_map_points(
        network.get(), x, y, 2, FASTNSA_PLNNS_RTREE, 10.0, 1
    ));
    FASTNSA_CHECK(points.get() != nullptr);
    double copied[2];
    FASTNSA_CHECK(fastnsa_pointset_copy_offsets(points.get(), copied, 2) == 0);
    FASTNSA_CHECK(copied[0] == length && copied[1] == 0.0);

    const int ids[]{0, 0};
    const double times[]{0.0, 1.0};
    const double roundoff[]{-1e-11, length + 1e-11};
    for (const bool temporal : {false, true}) {
        auto set = [&](const double* offsets) {
            return temporal
                ? fastnsa_pointset_set_st(points.get(), ids, offsets, times, 2)
                : fastnsa_pointset_set(points.get(), ids, offsets, 2);
        };
        FASTNSA_CHECK(set(roundoff) == 0);
        FASTNSA_CHECK(roundoff[0] == -1e-11 && roundoff[1] > length);
        FASTNSA_CHECK(fastnsa_pointset_copy_offsets(points.get(), copied, 2) == 0);
        FASTNSA_CHECK(copied[0] == 0.0 && copied[1] == length);
        const std::vector<std::vector<double>> invalid{
            {-1e-5, length}, {0.0, length + 1e-5},
            {std::numeric_limits<double>::quiet_NaN(), length},
            {0.0, std::numeric_limits<double>::infinity()},
        };
        for (const auto& offsets : invalid) {
            FASTNSA_CHECK(set(offsets.data()) == -2);
            FASTNSA_CHECK(fastnsa_pointset_copy_offsets(points.get(), copied, 2) == 0);
            FASTNSA_CHECK(copied[0] == 0.0 && copied[1] == length);
        }
    }
}

}  // namespace

int main() {
    return fastnsa::test::run("C API nearest-edge mapping", [] {
        test_packed_wkb_and_threads();
        test_endpoint_roundoff_and_atomic_setters();
    });
}

#include "test_support.h"

#include <algorithm>

namespace {

fastnsa::test::NetworkHandle make_chain(int nodes, double length) {
    std::vector<double> edges;
    std::vector<int> offsets(nodes + 1, 0);
    std::vector<int> indices;
    for (int edge = 0; edge < nodes - 1; ++edge) {
        edges.insert(edges.end(), {double(edge), double(edge + 1), length});
    }
    for (int node = 0; node < nodes; ++node) {
        if (node > 0) indices.push_back(node - 1);
        if (node + 1 < nodes) indices.push_back(node);
        offsets[node + 1] = static_cast<int>(indices.size());
    }
    return fastnsa::test::NetworkHandle(fastnsa_network_create(
        edges.data(), nodes - 1, offsets.data(), indices.data(), nodes
    ));
}

struct Output {
    std::vector<int> edges;
    std::vector<double> first, second, density;

    explicit Output(int count)
        : edges(count, -1), first(count, -1.0), second(count, -1.0),
          density(count, -1.0) {}
};

struct Fixture {
    fastnsa::test::NetworkHandle network;
    fastnsa::test::PointSetHandle points;
    std::vector<double> weights{1.0, 2.0, 0.5};
    double resolution, bandwidth;
    int count;

    Fixture(int nodes, double length)
        : network(make_chain(nodes, length)),
          points(fastnsa_pointset_create(network.get())),
          resolution(length * 0.6), bandwidth(nodes * length),
          count(fastnsa_nkdv_count_lixels(network.get(), resolution)) {
        FASTNSA_CHECK(network.get() != nullptr);
        FASTNSA_CHECK(points.get() != nullptr);
        FASTNSA_CHECK(count > 1);
        const int edges[]{0, nodes / 2, nodes - 2};
        const double offsets[]{length * 0.2, length * 0.7, length};
        FASTNSA_CHECK(fastnsa_pointset_set(points.get(), edges, offsets, 3) == 0);
    }

    int compute(fastnsa_nkdv_method method, Output& out,
                fastnsa_progress_callback callback = nullptr,
                void* data = nullptr, bool weighted = false) const {
        if (weighted) {
            return fastnsa_compute_nkdv_weighted(
                network.get(), points.get(), weights.data(), int(weights.size()),
                method, resolution, FASTNSA_NKDV_KERNEL_EPANECHNIKOV, bandwidth,
                count, out.edges.data(), out.first.data(), out.second.data(),
                out.density.data(), callback, data
            );
        }
        return fastnsa_compute_nkdv_with_progress(
            network.get(), points.get(), method, resolution,
            FASTNSA_NKDV_KERNEL_EPANECHNIKOV, bandwidth,
            count, out.edges.data(), out.first.data(), out.second.data(),
            out.density.data(), callback, data
        );
    }
};

void same_output(const Output& actual, const Output& expected) {
    FASTNSA_CHECK(actual.edges == expected.edges);
    FASTNSA_CHECK(actual.first == expected.first);
    FASTNSA_CHECK(actual.second == expected.second);
    FASTNSA_CHECK(actual.density.size() == expected.density.size());
    for (std::size_t i = 0; i < actual.density.size(); ++i) {
        FASTNSA_CHECK(std::isfinite(actual.density[i]));
        FASTNSA_CHECK(fastnsa::test::close(actual.density[i], expected.density[i], 1e-9));
    }
}

struct Cancellation {
    int polls = 0;
    int cancel_at;

    static int callback(const char*, int completed, int total, void* data) {
        auto& state = *static_cast<Cancellation*>(data);
        // All phase-level totals in this fixture exceed one. This targets
        // Dijkstra itself, while its queue still has unfinished work.
        if (completed == 0 && total == 1) {
            return ++state.polls >= state.cancel_at;
        }
        return 0;
    }
};

void test_cancel_inside_search_and_recover() {
    Fixture fixture(600, 10.0);
    for (const auto method : {FASTNSA_NKDV_BASIC, FASTNSA_NKDV_SHARED_SPD,
                              FASTNSA_NKDV_ADA, FASTNSA_NKDV_LION}) {
        for (const bool weighted : {false, true}) {
            Output baseline(fixture.count);
            FASTNSA_CHECK(fixture.compute(method, baseline, nullptr, nullptr, weighted) == 0);
            // Each full search polls three times. Poll five exercises the
            // second endpoint search for shared-SPD, ADA and LION.
            for (const int cancel_at : {2, 5}) {
                Output canceled(fixture.count);
                Cancellation state{0, cancel_at};
                FASTNSA_CHECK(fixture.compute(
                    method, canceled, Cancellation::callback, &state, weighted
                ) != 0);
                FASTNSA_CHECK(std::string(fastnsa_last_error()).find("canceled") != std::string::npos);
                FASTNSA_CHECK(state.polls == cancel_at);
                FASTNSA_CHECK(std::all_of(canceled.edges.begin(), canceled.edges.end(),
                    [](int value) { return value == -1; }));
                for (const auto* buffer : {&canceled.first, &canceled.second, &canceled.density}) {
                    FASTNSA_CHECK(std::all_of(buffer->begin(), buffer->end(),
                        [](double value) { return value == -1.0; }));
                }
                Output recovered(fixture.count);
                FASTNSA_CHECK(fixture.compute(method, recovered, nullptr, nullptr, weighted) == 0);
                same_output(recovered, baseline);
            }
        }
    }
}

void test_concurrent_ada_lion() {
    Fixture first(65, 11.0), second(91, 7.0);
    Output baseline_first(first.count), baseline_second(second.count);
    FASTNSA_CHECK(first.compute(FASTNSA_NKDV_ADA, baseline_first) == 0);
    FASTNSA_CHECK(second.compute(FASTNSA_NKDV_LION, baseline_second, nullptr, nullptr, true) == 0);
    fastnsa::test::run_concurrently(4, [&](std::size_t task) {
        for (int iteration = 0; iteration < 8; ++iteration) {
            if (task % 2 == 0) {
                Output out(first.count);
                FASTNSA_CHECK(first.compute(FASTNSA_NKDV_ADA, out) == 0);
                same_output(out, baseline_first);
            } else {
                Output out(second.count);
                FASTNSA_CHECK(second.compute(FASTNSA_NKDV_LION, out, nullptr, nullptr, true) == 0);
                same_output(out, baseline_second);
            }
        }
    });
}

}  // namespace

int main() {
    return fastnsa::test::run("C API NKDV isolation and cancellation", [] {
        test_cancel_inside_search_and_recover();
        test_concurrent_ada_lion();
    });
}

#include "test_support.h"

#include <algorithm>
#include <atomic>

namespace {

fastnsa::test::NetworkHandle make_chain(int nodes, double length) {
    std::vector<double> edges;
    std::vector<int> offsets(nodes + 1, 0), indices;
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
    std::vector<double> observed, simulations;
    Output(bool temporal, int count)
        : observed(temporal ? 4 : 2, -1.0),
          simulations(count * observed.size(), -1.0) {}
};

struct Fixture {
    fastnsa::test::NetworkHandle network;
    fastnsa::test::PointSetHandle points;
    std::vector<double> radii, times{0.0, 10.0};

    Fixture(int nodes, double length)
        : network(make_chain(nodes, length)),
          points(fastnsa_pointset_create(network.get())),
          radii{nodes * length, (nodes + 1) * length} {
        FASTNSA_CHECK(network.get() != nullptr);
        FASTNSA_CHECK(points.get() != nullptr);
        const int edges[]{0, 0, nodes / 2, nodes - 2};
        const double offsets[]{length * 0.2, length * 0.8, length * 0.3, length * 0.9};
        const double event_times[]{0.0, 3.0, 7.0, 10.0};
        FASTNSA_CHECK(fastnsa_pointset_set_st(
            points.get(), edges, offsets, event_times, 4
        ) == 0);
    }

    int compute(bool temporal, fastnsa_k_method method, int count, Output& out,
                fastnsa_progress_callback callback = nullptr, void* data = nullptr) const {
        double* simulated = out.simulations.empty() ? nullptr : out.simulations.data();
        if (temporal) {
            return fastnsa_compute_spatiotemporal_k_function_seeded(
                network.get(), points.get(), radii.data(), int(radii.size()),
                times.data(), int(times.size()), method, count,
                out.observed.data(), simulated, 1234u, callback, data
            );
        }
        return fastnsa_compute_k_function_seeded(
            network.get(), points.get(), radii.data(), int(radii.size()),
            method, count, out.observed.data(), simulated, 0.1, 1234u, callback, data
        );
    }
};

void same_output(const Output& actual, const Output& expected) {
    FASTNSA_CHECK(actual.observed == expected.observed);
    FASTNSA_CHECK(actual.simulations == expected.simulations);
    for (const auto* values : {&actual.observed, &actual.simulations}) {
        FASTNSA_CHECK(std::all_of(values->begin(), values->end(), [](double value) {
            return std::isfinite(value) && value >= 0.0;
        }));
    }
}

struct Cancellation {
    const char* stage;
    int cancel_after;
    int previous_completed = -1;
    int previous_total = -1;
    int repeated = 0;
    bool canceled = false;

    static int callback(const char* stage, int completed, int total, void* data) {
        auto& state = *static_cast<Cancellation*>(data);
        if (std::string(stage) != state.stage || completed <= 0) return 0;
        if (completed != state.previous_completed || total != state.previous_total) {
            state.previous_completed = completed;
            state.previous_total = total;
            state.repeated = 0;
        }
        // Outer edge callbacks advance progress; repeated positive progress
        // here comes from polling within a still-running shortest-path search.
        state.canceled = ++state.repeated >= state.cancel_after;
        return state.canceled;
    }
};

void test_cancellation_and_recovery(bool temporal, fastnsa_k_method method) {
    Fixture fixture(600, 1.0);
    for (const int count : {0, 2}) {
        Output baseline(temporal, count);
        FASTNSA_CHECK(fixture.compute(temporal, method, count, baseline) == 0);
        for (const int cancel_after : {3, 6}) {
            // RQS has one search per event; simulated events usually occupy
            // separate source edges. Its observed first edge has two events.
            if (count > 0 && temporal && method == FASTNSA_K_RQS && cancel_after == 6) continue;
            Output canceled(temporal, count);
            Cancellation state{count > 0 ? "Monte Carlo batch" : "backend compute", cancel_after};
            FASTNSA_CHECK(fixture.compute(
                temporal, method, count, canceled, Cancellation::callback, &state
            ) != 0);
            FASTNSA_CHECK(state.canceled);
            FASTNSA_CHECK(std::string(fastnsa_last_error()).find("canceled") != std::string::npos);
            for (const auto* values : {&canceled.observed, &canceled.simulations}) {
                FASTNSA_CHECK(std::all_of(values->begin(), values->end(),
                    [](double value) { return value == -1.0; }));
            }
            Output recovered(temporal, count);
            FASTNSA_CHECK(fixture.compute(temporal, method, count, recovered) == 0);
            same_output(recovered, baseline);
        }
    }
}

struct Job {
    bool temporal;
    fastnsa_k_method method;
};

const std::vector<Job> jobs{
    {false, FASTNSA_K_NS}, {false, FASTNSA_K_STARS}, {false, FASTNSA_K_SPS},
    {true, FASTNSA_K_RQS}, {true, FASTNSA_K_SPS}, {true, FASTNSA_K_EAR},
    {true, FASTNSA_K_MTS}, {true, FASTNSA_K_EARTH},
};

void test_concurrent_methods() {
    Fixture first(45, 7.0), second(63, 3.0);
    std::vector<Output> baseline;
    for (std::size_t i = 0; i < jobs.size(); ++i) {
        baseline.emplace_back(jobs[i].temporal, 2);
        const auto& fixture = i % 2 == 0 ? first : second;
        FASTNSA_CHECK(fixture.compute(jobs[i].temporal, jobs[i].method, 2, baseline.back()) == 0);
    }
    std::cout << "Serial baselines complete; starting independent worker threads" << std::endl;
    fastnsa::test::run_concurrently(jobs.size(), [&](std::size_t i) {
        const auto& fixture = i % 2 == 0 ? first : second;
        for (int repeat = 0; repeat < 4; ++repeat) {
            Output actual(jobs[i].temporal, 2);
            FASTNSA_CHECK(fixture.compute(jobs[i].temporal, jobs[i].method, 2, actual) == 0);
            same_output(actual, baseline[i]);
        }
    });
}

void test_concurrent_worker_failure_is_propagated() {
    std::atomic<int> completed{0};
    bool caught = false;
    try {
        fastnsa::test::run_concurrently(8, [&](std::size_t i) {
            ++completed;
            if (i == 0) throw std::runtime_error("expected worker failure");
        });
    } catch (const std::runtime_error& error) {
        caught = std::string(error.what()) == "expected worker failure";
    }
    FASTNSA_CHECK(caught);
    FASTNSA_CHECK(completed.load() == 8);
}

}  // namespace

int main() {
    return fastnsa::test::run("K-function queue isolation and cancellation", [] {
        std::cout << "Checking worker exception propagation" << std::endl;
        test_concurrent_worker_failure_is_propagated();
        for (const auto& job : jobs) {
            std::cout << "Checking " << (job.temporal ? "ST K" : "Network K")
                      << " method " << int(job.method) << " cancellation/recovery" << std::endl;
            test_cancellation_and_recovery(job.temporal, job.method);
        }
        std::cout << "Checking eight concurrent K-function calls" << std::endl;
        test_concurrent_methods();
    });
}

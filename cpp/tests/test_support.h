#pragma once

#include "fastnsa_capi.h"

#include <cmath>
#include <exception>
#include <future>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace fastnsa::test {

inline void require(bool condition, const char* expression, const char* file, int line) {
    if (condition) {
        return;
    }
    std::ostringstream message;
    message << file << ':' << line << ": check failed: " << expression;
    const char* native_error = fastnsa_last_error();
    if (native_error && *native_error) {
        message << " (native error: " << native_error << ')';
    }
    throw std::runtime_error(message.str());
}

inline bool close(double first, double second, double tolerance = 1e-10) {
    return std::abs(first - second) <= tolerance;
}

class NetworkHandle {
public:
    explicit NetworkHandle(void* pointer = nullptr) : pointer_(pointer) {}
    ~NetworkHandle() { fastnsa_network_destroy(pointer_); }
    NetworkHandle(const NetworkHandle&) = delete;
    NetworkHandle& operator=(const NetworkHandle&) = delete;
    void* get() const noexcept { return pointer_; }

private:
    void* pointer_;
};

class PointSetHandle {
public:
    explicit PointSetHandle(void* pointer = nullptr) : pointer_(pointer) {}
    ~PointSetHandle() { fastnsa_pointset_destroy(pointer_); }
    PointSetHandle(const PointSetHandle&) = delete;
    PointSetHandle& operator=(const PointSetHandle&) = delete;
    void* get() const noexcept { return pointer_; }

private:
    void* pointer_;
};

inline NetworkHandle make_two_edge_network() {
    std::vector<double> edges{
        0.0, 1.0, 20.0,
        1.0, 2.0, 10.0,
    };
    std::vector<int> offsets{0, 1, 3, 4};
    std::vector<int> indices{0, 0, 1, 1};
    return NetworkHandle(fastnsa_network_create(
        edges.data(), 2, offsets.data(), indices.data(), 3
    ));
}

template <typename Function>
void run_concurrently(std::size_t count, Function function) {
    std::promise<void> start;
    const auto gate = start.get_future().share();
    std::vector<std::thread> workers;
    std::vector<std::future<void>> results;
    workers.reserve(count);
    results.reserve(count);
    try {
        for (std::size_t i = 0; i < count; ++i) {
            // Do not use std::async here: MSVC schedules it through a task
            // runtime. Waiting for other tasks can starve that runtime's pool.
            std::packaged_task<void()> task([gate, &function, i] {
                gate.wait();
                function(i);
            });
            results.push_back(task.get_future());
            workers.emplace_back(std::move(task));
        }
    } catch (...) {
        // Thread creation can fail partway through. Unblock and join every
        // worker already started before allowing the exception to escape.
        start.set_value();
        for (auto& worker : workers) worker.join();
        throw;
    }
    start.set_value();
    for (auto& worker : workers) worker.join();
    // Report worker assertions only after all threads have finished safely.
    for (auto& result : results) result.get();
}

template <typename Function>
int run(const char* name, Function function) {
    try {
        function();
        std::cout << "PASS: " << name << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << name << ": " << error.what() << '\n';
        return 1;
    }
}

}  // namespace fastnsa::test

#define FASTNSA_CHECK(expression) \
    ::fastnsa::test::require(static_cast<bool>(expression), #expression, __FILE__, __LINE__)

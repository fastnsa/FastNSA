#include "nne/core/Threading.h"

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <limits>
#include <thread>

namespace fastnsa::nne {
namespace {

constexpr int automatic_thread_limit = 8;
constexpr int default_hard_thread_limit = 32;

int hardware_threads() noexcept {
    const unsigned int detected = std::thread::hardware_concurrency();
    if (detected == 0) {
        return 1;
    }
    return static_cast<int>(std::min<unsigned int>(
        detected,
        static_cast<unsigned int>(std::numeric_limits<int>::max())
    ));
}

int configured_thread_limit() noexcept {
    const char* value = std::getenv("FASTNSA_MAX_THREADS");
    if (!value || *value == '\0') {
        return default_hard_thread_limit;
    }
    errno = 0;
    char* end = nullptr;
    const long parsed = std::strtol(value, &end, 10);
    if (errno != 0 || end == value || *end != '\0' || parsed <= 0) {
        return default_hard_thread_limit;
    }
    return static_cast<int>(std::min<long>(parsed, std::numeric_limits<int>::max()));
}

}  // namespace

int effective_thread_count(
    int requested_threads,
    std::size_t work_items,
    std::size_t minimum_items_per_thread
) noexcept {
    if (work_items == 0) {
        return 1;
    }
    const int available = std::max(
        1,
        std::min(hardware_threads(), configured_thread_limit())
    );
    const int requested = requested_threads == 0
        ? std::min(available, automatic_thread_limit)
        : std::clamp(requested_threads, 1, available);
    const std::size_t granularity = std::max<std::size_t>(1, minimum_items_per_thread);
    const std::size_t useful_workers = 1 + (work_items - 1) / granularity;
    return std::max(1, std::min<int>(
        requested,
        static_cast<int>(std::min<std::size_t>(
            useful_workers,
            static_cast<std::size_t>(std::numeric_limits<int>::max())
        ))
    ));
}

}  // namespace fastnsa::nne

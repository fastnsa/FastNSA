#pragma once

#include <cstddef>

namespace fastnsa::nne {

// Resolve 0 (automatic) and explicit requests to a portable worker count.
// FASTNSA_MAX_THREADS can lower or raise the default hard limit of 32.
int effective_thread_count(
    int requested_threads,
    std::size_t work_items,
    std::size_t minimum_items_per_thread
) noexcept;

}  // namespace fastnsa::nne

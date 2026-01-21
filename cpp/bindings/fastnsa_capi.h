#pragma once
#include <cstdint>

#ifdef _WIN32
#define FASTNSA_API __declspec(dllexport)
#else
#define FASTNSA_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ---------- Network ---------- */

FASTNSA_API void* fastnsa_network_create(double* edges, int n_edges, int* offsets, int* indices, int n_nodes);
FASTNSA_API void  fastnsa_network_destroy(void* net);

FASTNSA_API void fastnsa_network_set_wkbs(
    void* net_ptr,
    const uint8_t* wkb,
    const int* offsets,
    int n_edges
);


/* ---------- PointSet on Network ---------- */

FASTNSA_API void* fastnsa_pointset_create(void* net);
FASTNSA_API void  fastnsa_pointset_destroy(void* pts);

FASTNSA_API int fastnsa_pointset_set(
    void* pts,
    const int* edge_ids,
    const double* offsets,
    int n_points
);


// Map points from (x, y) coordinates to (edge_id, offset) on the network
FASTNSA_API void* fastnsa_map_points(
    void* net,
    const double* xs,
    const double* ys,
    int n_points
);



/* ---------- Network K-function ---------- */
// TODO: 修改对应值
typedef enum {
    FASTNSA_K_NS = 0,
    FASTNSA_K_SPS = 1,
    FASTNSA_K_STARS = 2
} fastnsa_k_method;


FASTNSA_API void fastnsa_compute_k_function(
    void* net,
    void* pts,
    const double* r_values,
    int n_r,
    fastnsa_k_method method,
    int simulations,
    double* out_k,          // (n_r)
    double* out_k_simulations  // (simulations, n_r) or NULL if not requested
);

/* ---------- Network Kernel Density Visualization ---------- */
// TODO: implement this function in the future
FASTNSA_API int fastnsa_compute_nkdv(
    void* net,
    void* pts,
    const double* distances,
    int n,
    double r_max,
    int n_bins,
    double* out_k
);



#ifdef __cplusplus
}
#endif

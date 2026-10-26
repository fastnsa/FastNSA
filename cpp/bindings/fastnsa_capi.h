#pragma once
#include <stdint.h>

#if defined(_WIN32) || defined(__CYGWIN__)
#  if defined(FASTNSA_BUILDING_DLL)
#    define FASTNSA_API __declspec(dllexport)
#  else
#    define FASTNSA_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define FASTNSA_API __attribute__((visibility("default")))
#else
#  define FASTNSA_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

FASTNSA_API int fastnsa_abi_version(void);
FASTNSA_API const char* fastnsa_last_error(void);

/* Return non-zero from the callback to request cooperative cancellation. */
typedef int (*fastnsa_progress_callback)(
    const char* stage,
    int completed,
    int total,
    void* user_data
);

/* ---------- Network ---------- */

FASTNSA_API void* fastnsa_network_create(double* edges, int n_edges, int* offsets, int* indices, int n_nodes);
FASTNSA_API void  fastnsa_network_destroy(void* net);

FASTNSA_API void fastnsa_network_set_wkbs(
    void* net_ptr,
    const uint8_t* wkb,
    const int* offsets,
    int n_edges
);
FASTNSA_API int fastnsa_network_set_wkbs_checked(
    void* net_ptr,
    const uint8_t* wkb,
    const int* offsets,
    int n_edges
);
/*
 * Set edge geometry without WKB serialization. ``xy`` contains interleaved
 * x/y coordinates and ``vertex_offsets`` contains n_edges + 1 vertex offsets.
 * Each edge must contain at least two vertices. The arrays are copied before
 * this function returns.
 */
FASTNSA_API int fastnsa_network_set_polylines(
    void* net_ptr,
    const double* xy,
    const int* vertex_offsets,
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
FASTNSA_API int fastnsa_pointset_set_st(
    void* pts,
    const int* edge_ids,
    const double* offsets,
    const double* times,
    int n_points
);
FASTNSA_API int fastnsa_pointset_set_times(
    void* pts,
    const double* times,
    int n_points
);
FASTNSA_API int fastnsa_pointset_size(void* pts);
FASTNSA_API int fastnsa_pointset_copy_edge_ids(void* pts, int* out_edge_ids, int n_points);
FASTNSA_API int fastnsa_pointset_copy_offsets(void* pts, double* out_offsets, int n_points);
FASTNSA_API int fastnsa_pointset_copy_times(void* pts, double* out_times, int n_points);



typedef enum {
    FASTNSA_PLNNS_RTREE = 0,
    FASTNSA_PLNNS_COMP_APPROX = 1,
    
    FASTNSA_PLNNS_COMP_EXACT = 3
} fastnsa_plnns_method;


// Map points from (x, y) coordinates to (edge_id, offset) on the network
FASTNSA_API void* fastnsa_map_points(
    void* net,
    const double* xs,
    const double* ys,
    int n_points,
    fastnsa_plnns_method method,
    double comp_eps, // only used if method == FASTNSA_METHOD_COMP_*
    int n_threads /* 0 selects the bounded automatic policy */
);


FASTNSA_API int fastnsa_pointset_dump_text (void* pts, const char* filename);
/* ST text times are elapsed from the earliest event, in the input time unit.
 * No seconds-to-days conversion is performed. */
FASTNSA_API int fastnsa_pointset_dump_text_st(void* pts, const char* filename);

/* ---------- Network K-function ---------- */
// TODO: 修改对应值
typedef enum {
    FASTNSA_K_RQS = 1,
    FASTNSA_K_SPS = 2,
    FASTNSA_K_EAR = 3,
    FASTNSA_K_MTS = 4,
    FASTNSA_K_EARTH = 5,
    FASTNSA_K_NS = 6,
    FASTNSA_K_STARS = 7
} fastnsa_k_method;


FASTNSA_API void fastnsa_compute_k_function(
    void* net,
    void* pts,
    const double* r_values,
    int n_r,
    fastnsa_k_method method,
    int simulations,
    double* out_k,          // (n_r)
    double* out_k_simulations,  // (simulations, n_r) or NULL if not requested
    const double  epsilon   // tolerance for approximate methods
);

FASTNSA_API int fastnsa_compute_k_function_checked(
    void* net,
    void* pts,
    const double* r_values,
    int n_r,
    fastnsa_k_method method,
    int simulations,
    double* out_k,
    double* out_k_simulations,
    const double epsilon
);
FASTNSA_API int fastnsa_compute_k_function_with_progress(
    void* net, void* pts, const double* r_values, int n_r,
    fastnsa_k_method method, int simulations, double* out_k,
    double* out_k_simulations, double epsilon,
    fastnsa_progress_callback callback, void* user_data
);
FASTNSA_API int fastnsa_compute_k_function_seeded(
    void* net, void* pts, const double* r_values, int n_r,
    fastnsa_k_method method, int simulations, double* out_k,
    double* out_k_simulations, double epsilon, uint32_t random_seed,
    fastnsa_progress_callback callback, void* user_data
);

FASTNSA_API void fastnsa_compute_spatiotemporal_k_function(
    void* net,
    void* pts,
    const double* s_values,
    int n_s,
    const double* t_values,
    int n_t,
    fastnsa_k_method method,
    int simulations,
    double* out_k,              // (n_s, n_t), row-major
    double* out_k_simulations   // (simulations, n_s, n_t), row-major or NULL
);

FASTNSA_API int fastnsa_compute_spatiotemporal_k_function_checked(
    void* net,
    void* pts,
    const double* s_values,
    int n_s,
    const double* t_values,
    int n_t,
    fastnsa_k_method method,
    int simulations,
    double* out_k,
    double* out_k_simulations
);
FASTNSA_API int fastnsa_compute_spatiotemporal_k_function_with_progress(
    void* net, void* pts, const double* s_values, int n_s,
    const double* t_values, int n_t, fastnsa_k_method method,
    int simulations, double* out_k, double* out_k_simulations,
    fastnsa_progress_callback callback, void* user_data
);
FASTNSA_API int fastnsa_compute_spatiotemporal_k_function_seeded(
    void* net, void* pts, const double* s_values, int n_s,
    const double* t_values, int n_t, fastnsa_k_method method,
    int simulations, double* out_k, double* out_k_simulations,
    uint32_t random_seed, fastnsa_progress_callback callback, void* user_data
);

/* ---------- Network Kernel Density Visualization ---------- */
typedef enum {
    FASTNSA_NKDV_BASIC = 1,
    FASTNSA_NKDV_SHARED_SPD = 2,
    FASTNSA_NKDV_ADA = 3,
    FASTNSA_NKDV_IA = 4,
    FASTNSA_NKDV_ADA_IA = 5,
    FASTNSA_NKDV_LION = 6
} fastnsa_nkdv_method;

typedef enum {
    FASTNSA_NKDV_KERNEL_GAUSSIAN = 0,
    FASTNSA_NKDV_KERNEL_TRIANGULAR = 1,
    FASTNSA_NKDV_KERNEL_EPANECHNIKOV = 2,
    FASTNSA_NKDV_KERNEL_QUARTIC = 3
} fastnsa_nkdv_kernel;

FASTNSA_API int fastnsa_nkdv_count_lixels(
    void* net,
    double lixel_length
);

FASTNSA_API int fastnsa_compute_nkdv(
    void* net,
    void* pts,
    fastnsa_nkdv_method method,
    double lixel_length,
    fastnsa_nkdv_kernel kernel,
    double bandwidth,
    int n_lixels,
    int* out_edge_ids,
    double* out_offsets_n1,
    double* out_offsets_n2,
    double* out_density
);
FASTNSA_API int fastnsa_compute_nkdv_with_progress(
    void* net, void* pts, fastnsa_nkdv_method method,
    double lixel_length, fastnsa_nkdv_kernel kernel, double bandwidth,
    int n_lixels, int* out_edge_ids, double* out_offsets_n1,
    double* out_offsets_n2, double* out_density,
    fastnsa_progress_callback callback, void* user_data
);
FASTNSA_API int fastnsa_compute_nkdv_weighted(
    void* net, void* pts, const double* point_weights, int n_weights,
    fastnsa_nkdv_method method, double lixel_length,
    fastnsa_nkdv_kernel kernel, double bandwidth, int n_lixels,
    int* out_edge_ids, double* out_offsets_n1, double* out_offsets_n2,
    double* out_density, fastnsa_progress_callback callback, void* user_data
);



#ifdef __cplusplus
}
#endif

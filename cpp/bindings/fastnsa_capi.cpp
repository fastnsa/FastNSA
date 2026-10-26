#include "fastnsa_capi.h"
#include "edge_offset.h"
#include "kdv/compute.h"
#include "stars/K_function.h"
#include "earth/K_function.h"
#include "nne/core/NNE.h"
#include "nne/geometry/Polyline.h"

#include <algorithm>
#include <vector>
#include <memory>
#include <new>
#include <cmath>
#include <exception>
#include <limits>

#include <assert.h>
#include <cstring>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <string>

struct Network {
    int n_nodes;
    int n_edges;

    std::vector<int> edge_u, edge_v;
    std::vector<double> edge_len;

    std::vector<int> csr_offsets;
    std::vector<int> csr_indices;

    std::unique_ptr<fastnsa::nne::NNE> nearest_edge_index;
};

struct PointSet {
    const Network* network; // const Network* network_ptr;
    std::vector<int> edge_ids;
    std::vector<double> offsets;
    std::vector<double> times;
};

namespace {
thread_local std::string last_error_message;
thread_local fastnsa_progress_callback active_progress_callback = nullptr;
thread_local void* active_progress_user_data = nullptr;
thread_local uint32_t active_random_seed = 5489u;
thread_local const double* active_point_weights = nullptr;
thread_local int active_point_weight_count = 0;

void set_last_error(const std::string& message) {
    last_error_message = message;
}

std::unique_ptr<fastnsa::nne::NNE> make_nearest_edge_index(
    std::vector<fastnsa::nne::Polyline> polylines,
    int expected_edges
) {
    if (static_cast<int>(polylines.size()) != expected_edges) {
        throw std::invalid_argument("Network geometry count does not match the edge count.");
    }
    return std::make_unique<fastnsa::nne::NNE>(std::move(polylines));
}

bool validate_thresholds(const double* values, int count, bool strictly_increasing) {
    for (int index = 0; index < count; ++index) {
        if (!std::isfinite(values[index]) || values[index] < 0.0) {
            return false;
        }
        if (index > 0 && (strictly_increasing
                ? values[index] <= values[index - 1]
                : values[index] < values[index - 1])) {
            return false;
        }
    }
    return true;
}

bool report_progress(const char* stage, int completed, int total) {
    return active_progress_callback &&
        active_progress_callback(stage, completed, total, active_progress_user_data) != 0;
}

struct ProgressScope {
    ProgressScope(fastnsa_progress_callback callback, void* user_data) {
        active_progress_callback = callback;
        active_progress_user_data = user_data;
    }
    ~ProgressScope() {
        active_progress_callback = nullptr;
        active_progress_user_data = nullptr;
    }
};

struct RandomSeedScope {
    explicit RandomSeedScope(uint32_t seed) : previous(active_random_seed) {
        active_random_seed = seed;
    }
    ~RandomSeedScope() { active_random_seed = previous; }
    uint32_t previous;
};

struct PointWeightScope {
    PointWeightScope(const double* weights, int count)
        : previous_weights(active_point_weights), previous_count(active_point_weight_count) {
        active_point_weights = weights;
        active_point_weight_count = count;
    }
    ~PointWeightScope() {
        active_point_weights = previous_weights;
        active_point_weight_count = previous_count;
    }
    const double* previous_weights;
    int previous_count;
};

struct StarsCleanupScope {
    explicit StarsCleanupScope(fastnsa::stars::model& value) : model(value) {}
    ~StarsCleanupScope() { fastnsa::stars::kfun_cleanup(model); }
    fastnsa::stars::model& model;
};

struct EarthCleanupScope {
    explicit EarthCleanupScope(fastnsa::earth::model& value) : model(value) {}
    ~EarthCleanupScope() { fastnsa::earth::kfun_cleanup(model); }
    fastnsa::earth::model& model;
};
}

extern "C" {

int fastnsa_abi_version() {
    return 2;
}

const char* fastnsa_last_error() {
    return last_error_message.c_str();
}

void* fastnsa_network_create(double* edges, int n_edges, int* offsets, int* indices, int n_nodes) {
    last_error_message.clear();
    if (!edges || !offsets || !indices || n_edges <= 0 || n_nodes <= 0) {
        set_last_error("Invalid network arrays or dimensions.");
        return nullptr;
    }
    if (offsets[0] != 0 || offsets[n_nodes] < 0) {
        set_last_error("Invalid CSR offsets.");
        return nullptr;
    }
    for (int node = 0; node < n_nodes; ++node) {
        if (offsets[node] > offsets[node + 1]) {
            set_last_error("CSR offsets must be monotonic.");
            return nullptr;
        }
    }

    try {
    auto net = std::make_unique<Network>();
    net->n_nodes = n_nodes;
    net->n_edges = n_edges;

        // --------- 1. edges: u, v, length ----------
    net->edge_u.resize(n_edges);
    net->edge_v.resize(n_edges);
    net->edge_len.resize(n_edges);

    for (int e = 0; e < n_edges; ++e) {
        // edges layout: [u, v, len]
        if (!std::isfinite(edges[3 * e]) || !std::isfinite(edges[3 * e + 1]) ||
            !std::isfinite(edges[3 * e + 2]) ||
            std::floor(edges[3 * e]) != edges[3 * e] ||
            std::floor(edges[3 * e + 1]) != edges[3 * e + 1]) {
            set_last_error("Network edges contain invalid node IDs or lengths.");
            return nullptr;
        }
        int u = static_cast<int>(edges[3 * e + 0]);
        int v = static_cast<int>(edges[3 * e + 1]);
        double len = edges[3 * e + 2];

        if (u < 0 || u >= n_nodes || v < 0 || v >= n_nodes || len < 0.0) {
            set_last_error("Network edge values are outside valid bounds.");
            return nullptr;
        }

        net->edge_u[e] = u;
        net->edge_v[e] = v;
        net->edge_len[e] = len;
    }

    // --------- 2. CSR adjacency ----------
    // offsets: size n_nodes + 1
    // indices: size offsets[n_nodes]

    if (offsets != nullptr && indices != nullptr) {
        net->csr_offsets.assign(offsets, offsets + n_nodes + 1);

        int nnz = offsets[n_nodes];
        for (int i = 0; i < nnz; ++i) {
            if (indices[i] < 0 || indices[i] >= n_edges) {
                set_last_error("CSR indices contain an invalid edge ID.");
                return nullptr;
            }
        }
        net->csr_indices.assign(indices, indices + nnz);
    }


    return net.release();
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return nullptr;
    } catch (...) {
        set_last_error("Unknown network construction failure.");
        return nullptr;
    }
}

void fastnsa_network_destroy(void* net) {
    delete static_cast<Network*>(net);
}

int fastnsa_network_set_wkbs_checked(
    void* net_ptr,
    const uint8_t* wkb,
    const int* offsets,
    int n_edges
){
    last_error_message.clear();
    auto* net = static_cast<Network*>(net_ptr);
    if (!net || !wkb || !offsets || n_edges != net->n_edges || offsets[0] != 0) {
        set_last_error("Invalid network WKB buffer.");
        return -1;
    }
    for (int edge = 0; edge < n_edges; ++edge) {
        if (offsets[edge] > offsets[edge + 1]) {
            set_last_error("Network WKB offsets must be monotonic.");
            return -2;
        }
    }

    try {
        std::vector<fastnsa::nne::Polyline> polylines;
        polylines.reserve(static_cast<std::size_t>(n_edges));
        for (int edge = 0; edge < n_edges; ++edge) {
            const int begin = offsets[edge];
            const int end = offsets[edge + 1];
            if (begin < 0 || end <= begin) {
                set_last_error("Network WKB contains an empty edge geometry.");
                return -2;
            }
            polylines.push_back(fastnsa::nne::parse_linestring_wkb(
                wkb + begin,
                static_cast<std::size_t>(end - begin)
            ));
        }
        net->nearest_edge_index = make_nearest_edge_index(std::move(polylines), n_edges);
        return 0;
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return -3;
    } catch (...) {
        set_last_error("Unknown network geometry initialization failure.");
        return -3;
    }

}

int fastnsa_network_set_polylines(
    void* net_ptr,
    const double* xy,
    const int* vertex_offsets,
    int n_edges
) {
    last_error_message.clear();
    auto* net = static_cast<Network*>(net_ptr);
    if (!net || !xy || !vertex_offsets || n_edges != net->n_edges ||
        vertex_offsets[0] != 0) {
        set_last_error("Invalid packed network polyline arrays.");
        return -1;
    }
    try {
        std::vector<fastnsa::nne::Polyline> polylines;
        polylines.reserve(static_cast<std::size_t>(n_edges));
        for (int edge = 0; edge < n_edges; ++edge) {
            const int begin = vertex_offsets[edge];
            const int end = vertex_offsets[edge + 1];
            if (begin < 0 || end - begin < 2) {
                set_last_error("Each network edge must contain at least two vertices.");
                return -2;
            }
            std::vector<fastnsa::nne::Point2D> points;
            points.reserve(static_cast<std::size_t>(end - begin));
            for (int vertex = begin; vertex < end; ++vertex) {
                points.push_back({xy[2 * vertex], xy[2 * vertex + 1]});
            }
            polylines.emplace_back(std::move(points));
        }
        net->nearest_edge_index = make_nearest_edge_index(std::move(polylines), n_edges);
        return 0;
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return -3;
    } catch (...) {
        set_last_error("Unknown packed network geometry initialization failure.");
        return -3;
    }
}

void fastnsa_network_set_wkbs(
    void* net_ptr,
    const uint8_t* wkb,
    const int* offsets,
    int n_edges
) {
    fastnsa_network_set_wkbs_checked(net_ptr, wkb, offsets, n_edges);
}


void* fastnsa_pointset_create(void* net) {
    last_error_message.clear();
    if (!net) {
        set_last_error("Cannot create a point set without a network.");
        return nullptr;
    }
    try {
        auto p = std::make_unique<PointSet>();
        p->network = static_cast<Network*>(net);
        return p.release();
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return nullptr;
    } catch (...) {
        set_last_error("Unknown point set construction failure.");
        return nullptr;
    }
}

void fastnsa_pointset_destroy(void* pts) {
    delete static_cast<PointSet*>(pts);
}

int fastnsa_pointset_set(
    void* pts,
    const int* edge_ids,
    const double* offsets,
    int n_points
) {
    last_error_message.clear();
    auto* p = static_cast<PointSet*>(pts);
    if (!p || !p->network || !edge_ids || !offsets || n_points < 0) {
        set_last_error("Invalid point set inputs.");
        return -3;
    }

    for (int i = 0; i < n_points; i++) {
        if (edge_ids[i] < 0 || edge_ids[i] >= p->network->n_edges) {
            set_last_error("Point set contains an invalid edge ID.");
            return -1;
        }
        double len = p->network->edge_len[edge_ids[i]];
        double canonical;
        if (!fastnsa::canonical_edge_offset(offsets[i], len, canonical)) {
            set_last_error("Point set contains an invalid edge offset.");
            return -2;
        }
    }

    try {
        std::vector<int> new_edge_ids(edge_ids, edge_ids + n_points);
        std::vector<double> new_offsets(offsets, offsets + n_points);
        for (int i = 0; i < n_points; ++i)
            fastnsa::canonical_edge_offset(offsets[i], p->network->edge_len[edge_ids[i]], new_offsets[i]);
        std::vector<double> new_times(n_points);
        for (int i = 0; i < n_points; ++i)
            new_times[i] = static_cast<double>(i);
        p->edge_ids.swap(new_edge_ids);
        p->offsets.swap(new_offsets);
        p->times.swap(new_times);
        return 0;
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return -5;
    } catch (...) {
        set_last_error("Unknown point set initialization failure.");
        return -5;
    }
}

int fastnsa_pointset_set_st(
    void* pts,
    const int* edge_ids,
    const double* offsets,
    const double* times,
    int n_points
) {
    last_error_message.clear();
    auto* p = static_cast<PointSet*>(pts);
    if (!p || !p->network || !edge_ids || !offsets || !times || n_points < 0) {
        set_last_error("Invalid spatiotemporal point set inputs.");
        return -3;
    }

    for (int i = 0; i < n_points; i++) {
        if (edge_ids[i] < 0 || edge_ids[i] >= p->network->n_edges) {
            set_last_error("Point set contains an invalid edge ID.");
            return -1;
        }
        double len = p->network->edge_len[edge_ids[i]];
        double canonical;
        if (!fastnsa::canonical_edge_offset(offsets[i], len, canonical)) {
            set_last_error("Point set contains an invalid edge offset.");
            return -2;
        }
        if (!std::isfinite(times[i])) {
            set_last_error("Point set contains an invalid time value.");
            return -4;
        }
    }

    try {
        std::vector<int> new_edge_ids(edge_ids, edge_ids + n_points);
        std::vector<double> new_offsets(offsets, offsets + n_points);
        for (int i = 0; i < n_points; ++i)
            fastnsa::canonical_edge_offset(offsets[i], p->network->edge_len[edge_ids[i]], new_offsets[i]);
        std::vector<double> new_times(times, times + n_points);
        p->edge_ids.swap(new_edge_ids);
        p->offsets.swap(new_offsets);
        p->times.swap(new_times);
        return 0;
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return -5;
    } catch (...) {
        set_last_error("Unknown spatiotemporal point set initialization failure.");
        return -5;
    }
}

int fastnsa_pointset_set_times(
    void* pts,
    const double* times,
    int n_points
) {
    last_error_message.clear();
    auto* p = static_cast<PointSet*>(pts);
    if (!p || !times || n_points < 0) {
        set_last_error("Invalid point time inputs.");
        return -4;
    }
    if (static_cast<int>(p->edge_ids.size()) != n_points ||
        static_cast<int>(p->offsets.size()) != n_points) {
        set_last_error("Point time count does not match the point set size.");
        return -3;
    }
    for (int i = 0; i < n_points; ++i) {
        if (!std::isfinite(times[i])) {
            set_last_error("Point set contains an invalid time value.");
            return -2;
        }
    }
    try {
        std::vector<double> new_times(times, times + n_points);
        p->times.swap(new_times);
        return 0;
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return -5;
    } catch (...) {
        set_last_error("Unknown point time initialization failure.");
        return -5;
    }
}

int fastnsa_pointset_size(void* pts) {
    if (!pts) return -1;
    auto* p = static_cast<PointSet*>(pts);
    return static_cast<int>(p->edge_ids.size());
}

int fastnsa_pointset_copy_edge_ids(void* pts, int* out_edge_ids, int n_points) {
    if (!pts || !out_edge_ids || n_points < 0) return -1;
    auto* p = static_cast<PointSet*>(pts);
    if (static_cast<int>(p->edge_ids.size()) != n_points) return -2;
    std::memcpy(out_edge_ids, p->edge_ids.data(), sizeof(int) * n_points);
    return 0;
}

int fastnsa_pointset_copy_offsets(void* pts, double* out_offsets, int n_points) {
    if (!pts || !out_offsets || n_points < 0) return -1;
    auto* p = static_cast<PointSet*>(pts);
    if (static_cast<int>(p->offsets.size()) != n_points) return -2;
    std::memcpy(out_offsets, p->offsets.data(), sizeof(double) * n_points);
    return 0;
}

int fastnsa_pointset_copy_times(void* pts, double* out_times, int n_points) {
    if (!pts || !out_times || n_points < 0) return -1;
    auto* p = static_cast<PointSet*>(pts);
    if (static_cast<int>(p->times.size()) != n_points) return -2;
    std::memcpy(out_times, p->times.data(), sizeof(double) * n_points);
    return 0;
}

// todo 迁移逻辑到src/内部
void* fastnsa_map_points(
    void* net,
    const double* xs,
    const double* ys,
    int n_points,
    fastnsa_plnns_method method,
    double comp_eps, // only used if method == FASTNSA_METHOD_COMP_*
    int n_threads 
) {
    last_error_message.clear();
    if (!net || !xs || !ys || n_points <= 0 || !std::isfinite(comp_eps) ||
        comp_eps <= 0.0 || n_threads < 0) {
        set_last_error("Invalid point mapping inputs.");
        return nullptr;
    }
    if (method != FASTNSA_PLNNS_RTREE &&
        method != FASTNSA_PLNNS_COMP_APPROX &&
        method != FASTNSA_PLNNS_COMP_EXACT) {
        set_last_error("Unsupported point mapping method.");
        return nullptr;
    }
    try {
    auto pts = std::make_unique<PointSet>();
    auto* network = static_cast<Network*>(net);
    pts->network = network;

    std::vector<fastnsa::nne::QueryPoint> points;
    points.reserve(n_points);

    for (int i = 0; i < n_points; ++i) {
        if (!std::isfinite(xs[i]) || !std::isfinite(ys[i])) {
            set_last_error("Point coordinates contain NaN or infinite values.");
            return nullptr;
        }
        points.push_back(fastnsa::nne::QueryPoint{
            xs[i],
            ys[i]
        });
    }

    if (!network->nearest_edge_index) {
        set_last_error("Network geometries have not been initialized.");
        return nullptr;
    }

    double grid_size = comp_eps/sqrt(2.0);


    auto results = network->nearest_edge_index->query_nearest_mt(
        points, static_cast<int>(method), grid_size, n_threads
    );
    if (results.size() != static_cast<size_t>(n_points)) {
        set_last_error("Nearest-edge matching returned an unexpected result count.");
        return nullptr;
    }


    pts->edge_ids.resize(n_points);
    pts->offsets.resize(n_points);
    pts->times.resize(n_points);

    for (int i = 0; i < n_points; ++i) {
        int edge_id = results[i];
        if (edge_id < 0 || edge_id >= network->n_edges) {
            set_last_error("Nearest-edge matching returned an invalid edge ID.");
            return nullptr;
        }

        double offset = network->nearest_edge_index->project_point_onto_edge(
            edge_id, xs[i], ys[i]
        );
        double canonical;
        if (!fastnsa::canonical_edge_offset(offset, network->edge_len[edge_id], canonical)) {
            set_last_error("Nearest-edge offset exceeds its edge length; check geometry and network length units.");
            return nullptr;
        }

        pts->edge_ids[i] = edge_id;
        pts->offsets[i] = canonical;
        pts->times[i] = static_cast<double>(i);
    }

    // Here we would implement the actual mapping logic.
    return pts.release();
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return nullptr;
    } catch (...) {
        set_last_error("Unknown point mapping failure.");
        return nullptr;
    }
}

// todo 迁移逻辑到src/内部
int fastnsa_pointset_dump_text(void* pts, const char* filename){
    if (!pts) return -1;

    auto* ps = static_cast<PointSet*>(pts);
    auto* net = ps->network;
    if (!net) return -1;

    std::ofstream out(filename);
    if (!out.is_open())
        return -2;

    out << std::fixed << std::setprecision(4);


    std::vector<std::vector<double>> edge_offsets(net->n_edges);

    for (size_t i = 0; i < ps->edge_ids.size(); i++) {
        int eid = ps->edge_ids[i];
        if (eid < 0 || eid >= net->n_edges) continue;
        edge_offsets[eid].push_back(ps->offsets[i]);
    }

    for (auto& v : edge_offsets) {
        std::sort(v.begin(), v.end());
    }

    out << net->n_nodes << " " << net->n_edges << "\n";

    for (int e = 0; e < net->n_edges; e++) {
        out << net->edge_u[e] << " "
            << net->edge_v[e] << " "
            << net->edge_len[e] << " "
            << edge_offsets[e].size();

        for (double off : edge_offsets[e]) {
            out << " " << off;
        }
        out << "\n";
    }
    
    out.close();
    return 0;
}

// todo 迁移逻辑到src/内部
int fastnsa_pointset_dump_text_st(void* pts, const char* filename){
    if (!pts) return -1;

    auto* ps = static_cast<PointSet*>(pts);
    auto* net = ps->network;
    if (!net) return -1;
    if (ps->edge_ids.size() != ps->offsets.size()) return -3;
    if (!ps->times.empty() && ps->times.size() != ps->edge_ids.size()) return -3;

    std::ofstream out(filename);
    if (!out.is_open())
        return -2;

    out << std::setprecision(std::numeric_limits<double>::max_digits10);

    const bool has_times = !ps->times.empty();
    double min_time = 0.0;
    if (has_times) {
        min_time = ps->times[0];
        for (size_t i = 1; i < ps->times.size(); ++i)
            if (ps->times[i] < min_time)
                min_time = ps->times[i];
    }

    std::vector<std::vector<std::pair<double, double>>> edge_pairs(net->n_edges);

    for (size_t i = 0; i < ps->edge_ids.size(); i++) {
        int eid = ps->edge_ids[i];
        if (eid < 0 || eid >= net->n_edges) continue;
        double t = has_times ? ps->times[i] - min_time : 0.0;
        edge_pairs[eid].emplace_back(ps->offsets[i], t);
    }

    for (auto& v : edge_pairs) {
        std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });
    }

    out << net->n_nodes << " " << net->n_edges << "\n";

    for (int e = 0; e < net->n_edges; e++) {
        out << net->edge_u[e] << " "
            << net->edge_v[e] << " "
            << net->edge_len[e] << " "
            << edge_pairs[e].size();

        for (const auto& p : edge_pairs[e]) {
            out << " " << p.first << " " << p.second;
        }
        out << "\n";
    }

    out.close();
    return 0;
}


int fastnsa_compute_k_function_checked(
    void* net,
    void* pts,
    const double* r_values,
    int n_r,
    fastnsa_k_method method,
    int simulations,
    double* out_k,          // (n_r)
    double* out_k_simulations,  // (simulations, n_r) or NULL if not requested
    const double  epsilon
) {
    last_error_message.clear();
    if (!net || !pts || !r_values || !out_k || n_r <= 0 || simulations < 0) {
        set_last_error("Invalid network K-function inputs.");
        return -1;
    }
    Network* network = static_cast<Network*>(net);
    PointSet* pointset = static_cast<PointSet*>(pts);
    if (pointset->network != network || pointset->edge_ids.empty() ||
        pointset->edge_ids.size() != pointset->offsets.size()) {
        set_last_error("The point set is empty, invalid, or belongs to another network.");
        return -1;
    }
    if (pointset->edge_ids.size() < 2) {
        set_last_error("Network K-function requires at least two event points.");
        return -1;
    }
    const double network_length = std::accumulate(
        network->edge_len.begin(), network->edge_len.end(), 0.0
    );
    if (!std::isfinite(network_length) || network_length <= 0.0) {
        set_last_error("Network K-function requires a positive finite network length.");
        return -1;
    }
    if (method != FASTNSA_K_SPS && method != FASTNSA_K_NS && method != FASTNSA_K_STARS) {
        set_last_error("Unsupported network K-function method.");
        return -1;
    }
    if (!std::isfinite(epsilon) || epsilon <= 0.0 ||
        !validate_thresholds(r_values, n_r, false)) {
        set_last_error("K-function thresholds and epsilon must be finite and valid.");
        return -1;
    }

    fastnsa::stars::model m{};
    StarsCleanupScope cleanup(m);

    try {
	m.progress = [](const char* stage, int completed, int total) {
	    return report_progress(stage, completed, total);
	};
    fastnsa::stars::kfun_init_parameters(
        m, r_values, n_r, method, simulations, epsilon, active_random_seed
    );
    fastnsa::stars::kfun_build_model(m, 
        network->n_nodes,  network->n_edges, 
        network->edge_u,
        network->edge_v,
        network->edge_len,
        pointset->edge_ids,
        pointset->offsets
    );
    fastnsa::stars::compute_K_function(m);

    const double point_count = static_cast<double>(pointset->edge_ids.size());
    const double k_scale = network_length / (point_count * (point_count - 1.0));

    for (int t = 0; t < n_r; ++t)
        out_k[t] = m.multiple_K_values[0][t] * k_scale;

    // simulations（如果需要）
    if (out_k_simulations) {
        for (int l = 1; l < m.L; ++l) {
            for (int t = 0; t < n_r; ++t) {
                out_k_simulations[(l - 1) * n_r + t]
                    = m.multiple_K_values[l][t] * k_scale;
            }
        }
    }

    return 0;
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return -2;
    } catch (...) {
        set_last_error("Unknown network K-function failure.");
        return -2;
    }
}

int fastnsa_compute_k_function_with_progress(
    void* net, void* pts, const double* r_values, int n_r,
    fastnsa_k_method method, int simulations, double* out_k,
    double* out_k_simulations, double epsilon,
    fastnsa_progress_callback callback, void* user_data
) {
    ProgressScope scope(callback, user_data);
    return fastnsa_compute_k_function_checked(
        net, pts, r_values, n_r, method, simulations,
        out_k, out_k_simulations, epsilon
    );
}

int fastnsa_compute_k_function_seeded(
    void* net, void* pts, const double* r_values, int n_r,
    fastnsa_k_method method, int simulations, double* out_k,
    double* out_k_simulations, double epsilon, uint32_t random_seed,
    fastnsa_progress_callback callback, void* user_data
) {
    ProgressScope progress_scope(callback, user_data);
    RandomSeedScope seed_scope(random_seed);
    return fastnsa_compute_k_function_checked(
        net, pts, r_values, n_r, method, simulations,
        out_k, out_k_simulations, epsilon
    );
}

void fastnsa_compute_k_function(
    void* net,
    void* pts,
    const double* r_values,
    int n_r,
    fastnsa_k_method method,
    int simulations,
    double* out_k,
    double* out_k_simulations,
    const double epsilon
) {
    fastnsa_compute_k_function_checked(
        net, pts, r_values, n_r, method, simulations,
        out_k, out_k_simulations, epsilon
    );
}

int fastnsa_compute_spatiotemporal_k_function_checked(
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
) {
    last_error_message.clear();
    if (!net || !pts || !s_values || !t_values || !out_k ||
        n_s <= 0 || n_t <= 0 || simulations < 0) {
        set_last_error("Invalid spatiotemporal K-function inputs.");
        return -1;
    }
    Network* network = static_cast<Network*>(net);
    PointSet* pointset = static_cast<PointSet*>(pts);
    if (pointset->network != network || pointset->edge_ids.empty() ||
        pointset->edge_ids.size() != pointset->offsets.size() ||
        pointset->edge_ids.size() != pointset->times.size()) {
        set_last_error("The spatiotemporal point set is invalid or belongs to another network.");
        return -1;
    }
    if (pointset->edge_ids.size() < 2) {
        set_last_error("Spatiotemporal K-function requires at least two event points.");
        return -1;
    }
    if (method < FASTNSA_K_RQS || method > FASTNSA_K_EARTH) {
        set_last_error("Unsupported spatiotemporal K-function method.");
        return -1;
    }
    const bool regular_grid_method = method == FASTNSA_K_MTS || method == FASTNSA_K_EARTH;
    if ((regular_grid_method && (n_s < 2 || n_t < 2)) ||
        !validate_thresholds(s_values, n_s, regular_grid_method) ||
        !validate_thresholds(t_values, n_t, regular_grid_method)) {
        set_last_error("Spatiotemporal thresholds are invalid for the selected method.");
        return -1;
    }
    if (regular_grid_method) {
        const double spatial_step = s_values[1] - s_values[0];
        const double temporal_step = t_values[1] - t_values[0];
        for (int index = 2; index < n_s; ++index) {
            if (std::abs((s_values[index] - s_values[index - 1]) - spatial_step) >
                1e-9 * std::max(1.0, std::abs(spatial_step))) {
                set_last_error("The selected method requires evenly spaced spatial thresholds.");
                return -1;
            }
        }
        for (int index = 2; index < n_t; ++index) {
            if (std::abs((t_values[index] - t_values[index - 1]) - temporal_step) >
                1e-9 * std::max(1.0, std::abs(temporal_step))) {
                set_last_error("The selected method requires evenly spaced temporal thresholds.");
                return -1;
            }
        }
    }

    const double network_length = std::accumulate(
        network->edge_len.begin(), network->edge_len.end(), 0.0
    );
    const auto time_bounds = std::minmax_element(
        pointset->times.begin(), pointset->times.end()
    );
    const double time_span = *time_bounds.second - *time_bounds.first;
    const double point_count = static_cast<double>(pointset->edge_ids.size());
    if (!std::isfinite(network_length) || network_length <= 0.0 ||
        !std::isfinite(time_span) || time_span <= 0.0) {
        set_last_error("Spatiotemporal K-function requires positive network length and time span.");
        return -1;
    }
    const double k_scale = network_length * time_span /
        (point_count * (point_count - 1.0));

    fastnsa::earth::model m{};
    EarthCleanupScope cleanup(m);

    try {
	m.progress = [](const char* stage, int completed, int total) {
	    return report_progress(stage, completed, total);
	};
    fastnsa::earth::kfun_init_parameters(
        m, s_values, n_s, t_values, n_t, method, simulations, active_random_seed
    );
    fastnsa::earth::kfun_build_model(
        m,
        network->n_nodes,
        network->n_edges,
        network->edge_u,
        network->edge_v,
        network->edge_len,
        pointset->edge_ids,
        pointset->offsets,
        pointset->times
    );
    fastnsa::earth::compute_K_function(m);

    for (int s = 0; s < n_s; ++s) {
        for (int t = 0; t < n_t; ++t) {
            out_k[s * n_t + t] = m.multiple_K_values[0][s][t] * k_scale;
        }
    }

    if (out_k_simulations) {
        for (int l = 1; l < m.L; ++l) {
            for (int s = 0; s < n_s; ++s) {
                for (int t = 0; t < n_t; ++t) {
                    out_k_simulations[(l - 1) * n_s * n_t + s * n_t + t] =
                        m.multiple_K_values[l][s][t] * k_scale;
                }
            }
        }
    }

    return 0;
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return -2;
    } catch (...) {
        set_last_error("Unknown spatiotemporal K-function failure.");
        return -2;
    }
}

int fastnsa_compute_spatiotemporal_k_function_with_progress(
    void* net, void* pts, const double* s_values, int n_s,
    const double* t_values, int n_t, fastnsa_k_method method,
    int simulations, double* out_k, double* out_k_simulations,
    fastnsa_progress_callback callback, void* user_data
) {
    ProgressScope scope(callback, user_data);
    return fastnsa_compute_spatiotemporal_k_function_checked(
        net, pts, s_values, n_s, t_values, n_t, method,
        simulations, out_k, out_k_simulations
    );
}

int fastnsa_compute_spatiotemporal_k_function_seeded(
    void* net, void* pts, const double* s_values, int n_s,
    const double* t_values, int n_t, fastnsa_k_method method,
    int simulations, double* out_k, double* out_k_simulations,
    uint32_t random_seed, fastnsa_progress_callback callback, void* user_data
) {
    ProgressScope progress_scope(callback, user_data);
    RandomSeedScope seed_scope(random_seed);
    return fastnsa_compute_spatiotemporal_k_function_checked(
        net, pts, s_values, n_s, t_values, n_t, method,
        simulations, out_k, out_k_simulations
    );
}

void fastnsa_compute_spatiotemporal_k_function(
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
) {
    fastnsa_compute_spatiotemporal_k_function_checked(
        net, pts, s_values, n_s, t_values, n_t, method,
        simulations, out_k, out_k_simulations
    );
}

int fastnsa_nkdv_count_lixels(
    void* net,
    double lixel_length
) {
    last_error_message.clear();
    if (!net) {
        set_last_error("Cannot count lixels without a network.");
        return -1;
    }
    if (!std::isfinite(lixel_length) || lixel_length <= 0.0) {
        set_last_error("Lixel length must be finite and positive.");
        return -2;
    }

    try {
        Network* network = static_cast<Network*>(net);
        const int count = fastnsa::kdv::count_lixels(network->edge_len, lixel_length);
        if (count < 0) {
            set_last_error("The requested lixel resolution produces an invalid or excessive output size.");
        }
        return count;
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return -3;
    } catch (...) {
        set_last_error("Unknown lixel count failure.");
        return -3;
    }
}

int fastnsa_compute_nkdv(
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
) {
    last_error_message.clear();
    if (!net || !pts) {
        set_last_error("NKDE requires a network and point set.");
        return -1;
    }
    if (!std::isfinite(lixel_length) || !std::isfinite(bandwidth) ||
        lixel_length <= 0.0 || bandwidth <= 0.0) {
        set_last_error("Lixel length and bandwidth must be finite and positive.");
        return -2;
    }
    if (!out_edge_ids || !out_offsets_n1 || !out_offsets_n2 || !out_density) {
        set_last_error("NKDE output buffers are invalid.");
        return -3;
    }

    Network* network = static_cast<Network*>(net);
    PointSet* pointset = static_cast<PointSet*>(pts);
    if (pointset->network != network || pointset->edge_ids.empty()) {
        set_last_error("The point set is empty or belongs to another network.");
        return -1;
    }

    const int expected_lixels = fastnsa::kdv::count_lixels(network->edge_len, lixel_length);
    if (expected_lixels < 0) {
        set_last_error("Could not determine the number of lixels.");
        return -4;
    }
    if (n_lixels != expected_lixels) {
        set_last_error("The NKDE output buffer size does not match the lixel count.");
        return -5;
    }

    try {
    fastnsa::kdv::Result result;
    std::string error;
    std::vector<double> point_weights;
    if (active_point_weights != nullptr) {
        if (active_point_weight_count != static_cast<int>(pointset->edge_ids.size())) {
            set_last_error("NKDE weight count does not match the point count.");
            return -5;
        }
        point_weights.assign(
            active_point_weights, active_point_weights + active_point_weight_count
        );
    }
    const bool ok = fastnsa::kdv::compute_from_arrays(
        network->n_nodes,
        network->edge_u,
        network->edge_v,
        network->edge_len,
        pointset->edge_ids,
        pointset->offsets,
        point_weights,
        static_cast<fastnsa::kdv::Method>(method),
        lixel_length,
        static_cast<fastnsa::kdv::Kernel>(kernel),
        bandwidth,
        result,
        &error,
        [](const char* stage, int completed, int total) {
            return report_progress(stage, completed, total);
        }
    );
    if (!ok) {
        set_last_error(error.empty() ? "NKDE computation failed." : error);
        return -6;
    }
    if (static_cast<int>(result.edge_ids.size()) != n_lixels ||
        static_cast<int>(result.offsets_n1.size()) != n_lixels ||
        static_cast<int>(result.offsets_n2.size()) != n_lixels ||
        static_cast<int>(result.density.size()) != n_lixels) {
        set_last_error("NKDE returned inconsistent lixel output arrays.");
        return -7;
    }

    std::memcpy(out_edge_ids, result.edge_ids.data(), sizeof(int) * n_lixels);
    std::memcpy(out_offsets_n1, result.offsets_n1.data(), sizeof(double) * n_lixels);
    for (int index = 0; index < n_lixels; ++index) {
        const double value = result.density[index];
        if (!std::isfinite(value)) {
            set_last_error("NKDE returned a non-finite density value.");
            return -7;
        }
        out_density[index] = std::max(0.0, value);
    }
    std::memcpy(out_offsets_n2, result.offsets_n2.data(), sizeof(double) * n_lixels);

    return 0;
    } catch (const std::exception& exc) {
        set_last_error(exc.what());
        return -8;
    } catch (...) {
        set_last_error("Unknown NKDE failure.");
        return -8;
    }
}

int fastnsa_compute_nkdv_with_progress(
    void* net, void* pts, fastnsa_nkdv_method method,
    double lixel_length, fastnsa_nkdv_kernel kernel, double bandwidth,
    int n_lixels, int* out_edge_ids, double* out_offsets_n1,
    double* out_offsets_n2, double* out_density,
    fastnsa_progress_callback callback, void* user_data
) {
    ProgressScope scope(callback, user_data);
    return fastnsa_compute_nkdv(
        net, pts, method, lixel_length, kernel, bandwidth, n_lixels,
        out_edge_ids, out_offsets_n1, out_offsets_n2, out_density
    );
}

int fastnsa_compute_nkdv_weighted(
    void* net, void* pts, const double* point_weights, int n_weights,
    fastnsa_nkdv_method method, double lixel_length,
    fastnsa_nkdv_kernel kernel, double bandwidth, int n_lixels,
    int* out_edge_ids, double* out_offsets_n1, double* out_offsets_n2,
    double* out_density, fastnsa_progress_callback callback, void* user_data
) {
    if (!point_weights || n_weights <= 0) {
        set_last_error("Weighted NKDE requires a non-empty weight array.");
        return -1;
    }
    ProgressScope progress_scope(callback, user_data);
    PointWeightScope weight_scope(point_weights, n_weights);
    return fastnsa_compute_nkdv(
        net, pts, method, lixel_length, kernel, bandwidth, n_lixels,
        out_edge_ids, out_offsets_n1, out_offsets_n2, out_density
    );
}

}

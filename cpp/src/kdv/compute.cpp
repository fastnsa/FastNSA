#include "kdv/compute.h"
#include "edge_offset.h"

#include "kdv/alg_NKDV.h"
#include "kdv/KAF.h"
#include "kdv/lion.h"
#include "kdv/shortest_path.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace fastnsa::kdv {
namespace {

bool set_error(std::string* error, const std::string& message) {
    if (error) {
        *error = message;
    }
    return false;
}

bool validate_parameters(Method method, double lixel_length, Kernel kernel, double bandwidth, std::string* error) {
    const int method_value = static_cast<int>(method);
    const int kernel_value = static_cast<int>(kernel);

    if (method_value < static_cast<int>(Method::BASIC) || method_value > static_cast<int>(Method::LION)) {
        return set_error(error, "Unsupported NKDV method.");
    }
    if (kernel_value < static_cast<int>(Kernel::GAUSSIAN) || kernel_value > static_cast<int>(Kernel::QUARTIC)) {
        return set_error(error, "Unsupported NKDV kernel.");
    }
    if (!std::isfinite(lixel_length) || lixel_length <= 0.0) {
        return set_error(error, "Lixel length must be positive.");
    }
    if (!std::isfinite(bandwidth) || bandwidth <= 0.0) {
        return set_error(error, "Bandwidth must be positive.");
    }
    if (kernel != Kernel::EPANECHNIKOV) {
        return set_error(error, "The current NKDV backend only supports the Epanechnikov kernel.");
    }
    return true;
}

void init_parameters(model& m, Method method, double lixel_length, Kernel kernel, double bandwidth) {
    m.method = static_cast<int>(method);
    m.lixel_reg_length = lixel_length;
    m.k_type = static_cast<int>(kernel);
    m.bandwidth = bandwidth;
    m.gamma = 1.0;
    m.network_fileName = nullptr;
    m.out_NKDV_fileName = nullptr;
    m.smallest_interval_size = 0.0;
    m.cur_edge_index = -1;
    m.num_deg = 0;
    m.lixel_total_number = 0;
    m.canceled = false;
    m.ell.a_L = nullptr;
    m.ell.a_R = nullptr;

    if (m.k_type == static_cast<int>(Kernel::GAUSSIAN)) {
        m.gamma = 9.0 / (2.0 * bandwidth * bandwidth);
    } else if (m.k_type == static_cast<int>(Kernel::TRIANGULAR)) {
        m.gamma = 1.0 / bandwidth;
    } else if (m.k_type == static_cast<int>(Kernel::EPANECHNIKOV) ||
               m.k_type == static_cast<int>(Kernel::QUARTIC)) {
        m.gamma = 1.0 / (bandwidth * bandwidth);
    }
}

bool load_network_from_arrays(
    model& m,
    int n_nodes,
    const std::vector<int>& edge_u,
    const std::vector<int>& edge_v,
    const std::vector<double>& edge_lengths,
    const std::vector<int>& point_edge_ids,
    const std::vector<double>& point_offsets,
    const std::vector<double>& point_weights,
    std::string* error
) {
    if (n_nodes <= 0) {
        return set_error(error, "Network must contain at least one node.");
    }
    if (edge_u.size() != edge_v.size() || edge_u.size() != edge_lengths.size()) {
        return set_error(error, "Edge arrays must have the same length.");
    }
    if (edge_u.empty() || edge_u.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        return set_error(error, "Network must contain a supported number of edges.");
    }
    if (point_edge_ids.size() != point_offsets.size()) {
        return set_error(error, "Point mapping arrays must have the same length.");
    }
    if (!point_weights.empty() && point_weights.size() != point_edge_ids.size()) {
        return set_error(error, "Point weights must be empty or match the point count.");
    }

    m.n = n_nodes;
    m.m = static_cast<int>(edge_u.size());
    m.edge_set = new Edge[static_cast<std::size_t>(m.m)];

    for (int e = 0; e < m.m; ++e) {
        if (edge_u[e] < 0 || edge_u[e] >= m.n || edge_v[e] < 0 || edge_v[e] >= m.n) {
            return set_error(error, "Edge endpoint is out of range.");
        }
        if (!std::isfinite(edge_lengths[e]) || edge_lengths[e] < 0.0) {
            return set_error(error, "Edge length cannot be negative.");
        }

        Edge& edge = m.edge_set[e];
        edge.n1 = edge_u[e];
        edge.n2 = edge_v[e];
        edge.length = edge_lengths[e];
        edge.is_IA = false;
        edge.num_intervals = 0;
        edge.min_dist_diff = 0.0;
        edge.last_interval_size = 0.0;
    }

    for (size_t i = 0; i < point_edge_ids.size(); ++i) {
        const int edge_index = point_edge_ids[i];
        if (edge_index < 0 || edge_index >= m.m) {
            return set_error(error, "Point edge_id is out of range.");
        }

        Edge& edge = m.edge_set[edge_index];
        double offset;
        if (!fastnsa::canonical_edge_offset(point_offsets[i], edge.length, offset)) {
            return set_error(error, "Point offset is out of range for its edge.");
        }

        Point p{};
        p.edge_index = edge_index;
        p.dist_n1 = offset;
        p.dist_n2 = edge.length - offset;
        p.weight = point_weights.empty() ? 1.0 : point_weights[i];
        if (!std::isfinite(p.weight) || p.weight < 0.0) {
            return set_error(error, "Point weights must be finite and non-negative.");
        }
        edge.PS.push_back(p);
        edge.aug_dist_diff_vec.push_back(-inf);
        edge.dist_n1_vec.push_back(-1.0);
        edge.dist_n2_vec.push_back(-1.0);
    }

    // ADA's prefix/suffix arrays and LION's same-edge sweep both require
    // events to be ordered along their hosting edge. Mapped input normally
    // arrives in source feature order and need not satisfy that invariant.
    for (int e = 0; e < m.m; ++e) {
        auto& points = m.edge_set[e].PS;
        std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
            return a.dist_n1 < b.dist_n1;
        });
    }

    m.Network.assign(m.n, {});
    for (int e = 0; e < m.m; ++e) {
        m.Network[m.edge_set[e].n1].push_back(e);
        m.Network[m.edge_set[e].n2].push_back(e);
    }

    m.sp_node_vec.resize(m.n);
    if (m.method >= static_cast<int>(Method::SHARED_SPD) && m.method <= static_cast<int>(Method::LION)) {
        m.sp_node_vec_node_a.resize(m.n);
        m.sp_node_vec_node_b.resize(m.n);
    }

    m.is_scan_edge_list = new bool[m.m];
    for (int e = 0; e < m.m; ++e) {
        m.is_scan_edge_list[e] = false;
    }

    return true;
}

bool report_progress(const ProgressCallback& progress, const char* stage, int completed, int total) {
    return progress && progress(stage, completed, total);
}

bool run_algorithm(model& m, const ProgressCallback& progress) {
    int prev_edge_index = -1;

    if (m.method == static_cast<int>(Method::BASIC)) {
        for (int l = 0; l < static_cast<int>(m.lixel_set.size()); ++l) {
            if (report_progress(progress, "backend compute", l, static_cast<int>(m.lixel_set.size()))) return false;
            clear_scan_edge_list(m);
            m.cur_l = m.lixel_set[l];
            dijkstra(m);
            if (m.canceled) return false;
            NKDV_basic(m);
            m.lixel_set[l] = m.cur_l;
        }
        return true;
    }

    if (m.method >= static_cast<int>(Method::SHARED_SPD) && m.method <= static_cast<int>(Method::ADA_IA)) {
        if (m.method == static_cast<int>(Method::ADA) || m.method == static_cast<int>(Method::ADA_IA)) {
            augment_preprocess(m);
        }
        if (m.method == static_cast<int>(Method::IA) || m.method == static_cast<int>(Method::ADA_IA)) {
            augment_interval_preprocess(m);
        }

        for (int l = 0; l < static_cast<int>(m.lixel_set.size()); ++l) {
            if (report_progress(progress, "backend compute", l, static_cast<int>(m.lixel_set.size()))) return false;
            m.cur_l = m.lixel_set[l];
            const int edge_index = m.cur_l.edge_index;

            if (edge_index != prev_edge_index) {
                clear_scan_edge_list(m);
                prev_edge_index = edge_index;
                m.sel_node_index = m.edge_set[edge_index].n1;
                m.node_index_pi_a = m.sel_node_index;
                dijkstra(m);
                if (m.canceled) return false;
                copy_sp_info(m, true);

                m.sel_node_index = m.edge_set[edge_index].n2;
                m.node_index_pi_b = m.sel_node_index;
                dijkstra(m);
                if (m.canceled) return false;
                copy_sp_info(m, false);
            }

            NKDV_basic(m);
            m.lixel_set[l] = m.cur_l;
        }
        return true;
    }

    if (m.method == static_cast<int>(Method::LION)) {
        if (report_progress(progress, "backend compute", 0, 3)) return false;
        init_LION(m);
        if (m.canceled) return false;
        if (report_progress(progress, "backend compute", 1, 3)) return false;
        lixel_augmentation(m);
        if (m.canceled) return false;
        if (report_progress(progress, "backend compute", 2, 3)) return false;
        lixel_aggregation(m);
        if (m.canceled) return false;
        one_D_KDV(m);
        if (m.canceled) return false;
    }
    return !report_progress(progress, "backend compute", 3, 3);
}

void collect_result(model& m, Result& out) {
    out.edge_ids.clear();
    out.offsets_n1.clear();
    out.offsets_n2.clear();
    out.density.clear();

    if (m.method >= static_cast<int>(Method::BASIC) && m.method <= static_cast<int>(Method::ADA_IA)) {
        out.edge_ids.reserve(m.lixel_set.size());
        out.offsets_n1.reserve(m.lixel_set.size());
        out.offsets_n2.reserve(m.lixel_set.size());
        out.density.reserve(m.lixel_set.size());
        for (const auto& lixel : m.lixel_set) {
            out.edge_ids.push_back(lixel.edge_index);
            out.offsets_n1.push_back(lixel.dist_n1);
            out.offsets_n2.push_back(lixel.dist_n2);
            out.density.push_back(lixel.KDE_value);
        }
        return;
    }

    if (m.method == static_cast<int>(Method::LION)) {
        out.edge_ids.reserve(m.lixel_total_number);
        out.offsets_n1.reserve(m.lixel_total_number);
        out.offsets_n2.reserve(m.lixel_total_number);
        out.density.reserve(m.lixel_total_number);
        for (int edge_index = 0; edge_index < m.m; ++edge_index) {
            for (const auto& lixel : m.edge_set[edge_index].lixel_vec) {
                out.edge_ids.push_back(edge_index);
                out.offsets_n1.push_back(lixel.dist_n1);
                out.offsets_n2.push_back(lixel.dist_n2);
                out.density.push_back(lixel.KDE_value);
            }
        }
    }
}

void cleanup_model(model& m) {
    m.lixel_set.clear();

    if (m.edge_set) {
        for (int e = 0; e < m.m; ++e) {
            if (m.method == static_cast<int>(Method::ADA) || m.method == static_cast<int>(Method::ADA_IA)) {
                m.edge_set[e].dist_n1_vec.clear();
                m.edge_set[e].dist_n2_vec.clear();
                m.edge_set[e].aug_dist_diff_vec.clear();
            }

            if (m.method == static_cast<int>(Method::LION)) {
                for (auto& lixel : m.edge_set[e].lixel_vec) {
                    delete[] lixel.alpha_B_x_deg_vec;
                    delete[] lixel.alpha_B_y_deg_vec;
                    delete[] lixel.alpha_R_x_deg_vec;
                    delete[] lixel.alpha_R_y_deg_vec;
                    lixel.alpha_B_x_deg_vec = nullptr;
                    lixel.alpha_B_y_deg_vec = nullptr;
                    lixel.alpha_R_x_deg_vec = nullptr;
                    lixel.alpha_R_y_deg_vec = nullptr;
                }
                m.edge_set[e].lixel_vec.clear();
            }

            m.edge_set[e].PS.clear();
            m.edge_set[e].aug_sum_dist_c_I_deg_1_vec.clear();
            m.edge_set[e].aug_sum_dist_c_I_deg_2_vec.clear();
            m.edge_set[e].aug_sum_dist_d_I_deg_1_vec.clear();
            m.edge_set[e].aug_sum_dist_d_I_deg_2_vec.clear();
            m.edge_set[e].interval_point_vec.clear();
            m.edge_set[e].weight_vec.clear();
            m.edge_set[e].agg_weight_c_vec.clear();
            m.edge_set[e].agg_weight_d_vec.clear();
        }
        delete[] m.edge_set;
        m.edge_set = nullptr;
    }

    if (m.is_scan_edge_list) {
        delete[] m.is_scan_edge_list;
        m.is_scan_edge_list = nullptr;
    }

    if (m.method == static_cast<int>(Method::LION)) {
        delete[] m.ell.a_L;
        delete[] m.ell.a_R;
        m.ell.a_L = nullptr;
        m.ell.a_R = nullptr;
    }

    m.Network.clear();
    m.sp_node_vec.clear();
    m.sp_node_vec_node_a.clear();
    m.sp_node_vec_node_b.clear();
    m.access_edge_list.clear();
}

}  // namespace

int count_lixels(const std::vector<double>& edge_lengths, double lixel_length) {
    if (!std::isfinite(lixel_length) || lixel_length <= 0.0) {
        return -1;
    }

    int total = 0;
    for (double length : edge_lengths) {
        if (!std::isfinite(length) || length < 0.0) {
            return -1;
        }
        if (length == 0.0) {
            continue;
        }
        const double edge_lixels = std::ceil(length / lixel_length);
        if (!std::isfinite(edge_lixels) ||
            edge_lixels > static_cast<double>(std::numeric_limits<int>::max() - total)) {
            return -1;
        }
        total += static_cast<int>(edge_lixels);
    }
    return total;
}

bool compute_from_arrays(
    int n_nodes,
    const std::vector<int>& edge_u,
    const std::vector<int>& edge_v,
    const std::vector<double>& edge_lengths,
    const std::vector<int>& point_edge_ids,
    const std::vector<double>& point_offsets,
    const std::vector<double>& point_weights,
    Method method,
    double lixel_length,
    Kernel kernel,
    double bandwidth,
    Result& out,
    std::string* error,
    ProgressCallback progress
) {
    if (!validate_parameters(method, lixel_length, kernel, bandwidth, error)) {
        return false;
    }

    model m{};
    struct CleanupGuard {
        model& value;
        ~CleanupGuard() { cleanup_model(value); }
    } cleanup{m};

    init_parameters(m, method, lixel_length, kernel, bandwidth);
    m.progress = progress;
    if (!load_network_from_arrays(
            m, n_nodes, edge_u, edge_v, edge_lengths,
            point_edge_ids, point_offsets, point_weights, error)) {
        return false;
    }

    obtain_lixel_set(m);
    if (!run_algorithm(m, progress)) {
        return set_error(error, "FastNSA operation canceled.");
    }
    collect_result(m, out);
    return true;
}

}  // namespace fastnsa::kdv

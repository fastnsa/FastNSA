#include "fastnsa_capi.h"
#include "stars/K_function.h"
#include "nne/core/NNE.h"

#include <vector>
#include <memory>

#include <geos/io/WKBReader.h>
#include <geos/linearref/LengthIndexedLine.h>
#include <assert.h>
#include <cstring>
#include <cstdint>

struct Network {
    int n_nodes;
    int n_edges;

    std::vector<int> edge_u, edge_v;
    std::vector<double> edge_len;

    std::vector<int> csr_offsets;
    std::vector<int> csr_indices;

    // std::vector<std::vector<unsigned char>> edge_wkbs;
    // edge wkbs 
    std::vector<uint8_t> wkb_data;
    std::vector<int> wkb_offsets;
};

struct PointSet {
    const Network* network; // const Network* network_ptr;
    std::vector<int> edge_ids;
    std::vector<double> offsets;
};

extern "C" {

void* fastnsa_network_create(double* edges, int n_edges, int* offsets, int* indices, int n_nodes) {
    auto* net = new Network();
    net->n_nodes = n_nodes;
    net->n_edges = n_edges;

        // --------- 1. edges: u, v, length ----------
    net->edge_u.resize(n_edges);
    net->edge_v.resize(n_edges);
    net->edge_len.resize(n_edges);

    for (int e = 0; e < n_edges; ++e) {
        // edges layout: [u, v, len]
        assert(std::floor(edges[3*e]) == edges[3*e]); // 未来考虑更改为 int 类型输入
        int u = static_cast<int>(edges[3 * e + 0]);
        int v = static_cast<int>(edges[3 * e + 1]);
        double len = edges[3 * e + 2];

        net->edge_u[e] = u;
        net->edge_v[e] = v;
        net->edge_len[e] = len;
    }

    // --------- 2. CSR adjacency ----------
    // offsets: size n_nodes + 1
    // indices: size offsets[n_nodes]

    net->csr_offsets.assign(offsets, offsets + n_nodes + 1);

    int nnz = offsets[n_nodes];
    net->csr_indices.assign(indices, indices + nnz);

    // return ;

    return net;
}

void fastnsa_network_destroy(void* net) {
    delete static_cast<Network*>(net);
}

 void fastnsa_network_set_wkbs(
    void* net_ptr,
    const uint8_t* wkb,
    const int* offsets,
    int n_edges
){
    auto* net = static_cast<Network*>(net_ptr);

    net->n_edges = n_edges;
    net->wkb_offsets.assign(offsets, offsets + n_edges + 1);

    int total = net->wkb_offsets.back();
    net->wkb_data.assign(wkb, wkb + total);

}


void* fastnsa_pointset_create(void* net) {
    auto* p = new PointSet();
    p->network = static_cast<Network*>(net);
    return p;
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
    auto* p = static_cast<PointSet*>(pts);
    assert(p->edge_ids.size() == p->offsets.size());

    for (int i = 0; i < n_points; i++) {
        if (edge_ids[i] < 0 || edge_ids[i] >= p->network->n_edges) // (edge_id越界检查)
            return -1;
        double len = p->network->edge_len[edge_ids[i]];
        if (offsets[i] < 0 || offsets[i] > len) //（offset合法性检查）
            return -2;
    }

    p->edge_ids.assign(edge_ids, edge_ids + n_points);
    p->offsets.assign(offsets, offsets + n_points);
    return 0;
}


void* fastnsa_map_points(
    void* net,
    const double* xs,
    const double* ys,
    int n_points
) {
    auto* pts = new PointSet();
    auto* network = static_cast<Network*>(net);
    pts->network = network;

    std::vector<QueryPoint> points;
    points.reserve(n_points);

    for (size_t i = 0; i < n_points; ++i) {
        points.push_back(QueryPoint{
            xs[i],
            ys[i]
        });
    }

    geos::io::WKBReader reader;
    std::vector<std::unique_ptr<geos::geom::Geometry>> geoms;
    for (int e = 0; e < network->n_edges; e++) {

        const uint8_t* wkb_ptr = network->wkb_data.data() + network->wkb_offsets[e];
        size_t wkb_len = network->wkb_offsets[e + 1] - network->wkb_offsets[e];
        geoms.emplace_back(reader.read(wkb_ptr, wkb_len));

        cout << "Edge " << e << " WKB: " << geoms.back()->toString()<< endl; //test

        // 用 geom 做 point-line projection
        // geom 在作用域结束自动释放
    }

    cout << "Building NNE index..." << endl; // test

    auto nne_index = NNE(std::move(geoms));

    cout << "Building NNE index done" << endl; // test

    int method = 0;  //single 
    double grid_size = 100.0;
    int num_threads = 4;

    // todo : 这里的 method, grid_size, num_threads 可以通过参数传入
    auto results = nne_index.query_nearest_mt(points, method, grid_size, num_threads); 

    cout << "done querying, results size:" << results.size() << endl; // test


    pts->edge_ids.resize(n_points);
    pts->offsets.resize(n_points);

    // vector<geos::linearref::LengthIndexedLine> indexed_lines;
    // indexed_lines.reserve(geoms.size());
    // for (const auto& geom : geoms) {
    //     indexed_lines.emplace_back(geom.get());
    // }

    cout << "done build lines" << endl; // test


    for (size_t i = 0; i < n_points; ++i) {
        int edge_id = results[i];
        cout << edge_id << endl; // test

        double offset = nne_index.project_point_onto_edge(edge_id, xs[i], ys[i]);

        cout << "mapped offset: " << offset << endl; // test

        pts->edge_ids[i] = edge_id;
        pts->offsets[i] = offset;
    }

    // Here we would implement the actual mapping logic.
    return pts;
}

void fastnsa_compute_k_function(
    void* net,
    void* pts,
    const double* r_values,
    int n_r,
    fastnsa_k_method method,
    int simulations,
    double* out_k,          // (n_r)
    double* out_k_simulations  // (simulations, n_r) or NULL if not requested
) {
    Network* network = static_cast<Network*>(net);
    PointSet* pointset = static_cast<PointSet*>(pts);

    model m;
    memset(&m, 0, sizeof(model));

    kfun_init_parameters(m, r_values, n_r, method, simulations);
    kfun_build_model(m, 
    network->n_nodes,  network->n_edges, 
    network->edge_u,
    network->edge_v,
    network->edge_len,
    pointset->edge_ids,
    pointset->offsets
    );
    compute_K_function(m);

    for (int t = 0; t < n_r; ++t)
        out_k[t] = m.multiple_K_values[0][t];

    // simulations（如果需要）
    if (out_k_simulations) {
        for (int l = 1; l < m.L; ++l) {
            for (int t = 0; t < n_r; ++t) {
                out_k_simulations[(l - 1) * n_r + t] 
                    = m.multiple_K_values[l][t];
            }
        }
    }

    kfun_cleanup(m);
}

}

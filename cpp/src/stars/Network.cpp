#include "stars/Network.h"
#include "fastnsa.h"

void init_PS_memory_edges(model& our_model, int num_edges)
{
	vector<Point> temp_vec{};

	for (int e = 0; e < num_edges; e++)
		for (int s = 0; s < our_model.L; s++){
			our_model.edge_set[e].PS_vec.push_back(temp_vec);
			// our_model.edge_set[e].sam_PS_vec.push_back(temp_vec);
		}
}
void init_sampling_memory(model& our_model)
{
	int sam_size=ceil(1.0/our_model.epsilon);
	int size=our_model.D*sam_size;

	our_model.minheap.resize(2);

	our_model.sam_edge=vector<Point>{};
}

void generate_points(model& our_model)
{
	default_random_engine generator_int;
	default_random_engine generator_double;
	int edge_index;
	Point pt;

	uniform_int_distribution<int> distribution(0, our_model.num_edges - 1);

	for (int l = 1; l < our_model.L; l++)
	{
		for (int i = 0; i < our_model.num_points; i++)
		{
			edge_index = distribution(generator_int);
			vector<Point>& PS_cur = our_model.edge_set[edge_index].PS_vec[l];
			uniform_real_distribution<double> distribution_double(0.0, our_model.edge_set[edge_index].length);
			pt.dist_n1 = distribution_double(generator_double);
			pt.dist_n2 = our_model.edge_set[edge_index].length - pt.dist_n1;
			pt.edge_index = edge_index;

			PS_cur.push_back(pt);
		}
	}
}

void sort_data_points(model& our_model)
{
	//Sorting the data points
	if (our_model.method >= 3 && our_model.method <= 8)
	{
		for (int e = 0; e < our_model.num_edges; e++)
			for(int l=0; l < our_model.L; l++){
				// if (our_model.edge_set[e].PS_vec[l].size()==0)
				// 	continue;
				sort(our_model.edge_set[e].PS_vec[l].begin(), our_model.edge_set[e].PS_vec[l].end());
			}
	}
}


void sort_selection(model& our_model, int e_index) {
	vector<Point>& PS=our_model.edge_set[e_index].PS_vec[our_model.cur_data_index];
	our_model.sam_edge=vector<Point>{};
	int sam_size=ceil(1.0/our_model.epsilon);
	if(sam_size>=PS.size()){
		our_model.sam_edge=PS;
		return;
	}
	for(int i=1;i<sam_size+1;i++){
		int sam_idx=ceil((i-1.0/2)*1.0/sam_size*PS.size()) - 1;
		our_model.sam_edge.push_back(PS[sam_idx]);
	}
}

void kfun_cleanup(model& m) {
    delete[] m.tau_vec;

    for (int l = 0; l < m.L; ++l)
        delete[] m.multiple_K_values[l];
    delete[] m.multiple_K_values;

    delete[] m.edge_set;

    delete[] m.is_scan_edge_list_first;
    delete[] m.is_scan_edge_list_second;
    delete[] m.detect_DoubleScan_edge_list;
}

void kfun_init_parameters(
    model& m,
    const double* r_values,
    int n_r,
    int method,
    int simulations
) {
    m.method = method;
    m.D = n_r;
    m.tau_vec = new double[m.D];
    for (int t = 0; t < m.D; ++t)
        m.tau_vec[t] = r_values[t];

    // m.incr_tau = m.tau_vec[1] - m.tau_vec[0];  //testing
    // m.start_tau = m.tau_vec[0];

    m.L = simulations + 1;   // 第 0 次是真实点集

    m.multiple_K_values = new double*[m.L];
    for (int l = 0; l < m.L; ++l) {
        m.multiple_K_values[l] = new double[m.D];
        std::fill(
            m.multiple_K_values[l],
            m.multiple_K_values[l] + m.D,
            0.0
        );
    }

    // ---- optional flags ----
    m.is_sam = false;  //TODO 添加 sampling flag
    m.epsilon = 0.1;   //TODO 添加 sampling epsilon
}

void kfun_build_model(model& m, 
    int n_nodes, int n_edges, 
    const std::vector<int>& edge_u,
    const std::vector<int>& edge_v,
    const std::vector<double>& edge_len,
    vector<int>& point_edge_ids,
    vector<double>& point_offsets
) {
    // This function would convert the Network and PointSet into the model structure used by the K-function computation.
    m.num_nodes = n_nodes;
    m.num_edges = n_edges;

    m.edge_set = new Edge[m.num_edges];
    init_PS_memory_edges(m, m.num_edges);

    if(m.method==7) //todo: 7 is STARS method, need to init sampling memory
		init_sampling_memory(m);

    for (int e = 0; e < n_edges; ++e) {
        m.edge_set[e].n1 = edge_u[e];
        m.edge_set[e].n2 = edge_v[e];
        m.edge_set[e].length = edge_len[e];
    }

    // Map points from the PointSet to the model's edge set
    for (size_t i = 0; i < point_edge_ids.size(); ++i) {
        int e = point_edge_ids[i];
        double d = point_offsets[i];

        Point pt;
        pt.edge_index = e;
        pt.dist_n1 = d;
        pt.dist_n2 = m.edge_set[e].length - d;

        m.edge_set[e].PS_vec[0].push_back(pt);
    }

    m.num_points += point_edge_ids.size();

    generate_points(m);

    // Build the adjacency list for the network
    m.Network.clear();
    m.Network.resize(n_nodes);
    for (int e = 0; e < n_edges; ++e) {
        m.Network[m.edge_set[e].n1].push_back(e);
        m.Network[m.edge_set[e].n2].push_back(e);
    }


}
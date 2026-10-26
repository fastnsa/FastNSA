#include "earth/Network.h"

namespace fastnsa::earth
{


void init_PS_memory_edges(model& our_model, int num_edges)
{
	vector<Point> temp_vec;

	for (int e = 0; e < num_edges; e++)
		for (int s = 0; s < our_model.L; s++)
			our_model.edge_set[e].PS_vec.push_back(temp_vec);
}

void generate_points(model& our_model)
{
	mt19937 generator_int(our_model.random_seed);
	mt19937 generator_double(our_model.random_seed + 1u);
	mt19937 generator_double_time(our_model.random_seed + 2u);
	int edge_index;
	Point pt;

	vector<double> cumulative_lengths(our_model.num_edges);
	double network_length = 0.0;
	for (int e = 0; e < our_model.num_edges; ++e)
	{
		network_length += our_model.edge_set[e].length;
		cumulative_lengths[e] = network_length;
	}
	if (network_length <= 0.0)
		return;

	uniform_real_distribution<double> distribution(0.0, network_length);

	for (int l = 1; l < our_model.L; l++)
	{
		for (int i = 0; i < our_model.num_points; i++)
		{
			const double network_offset = distribution(generator_int);
			auto edge_it = upper_bound(
				cumulative_lengths.begin(), cumulative_lengths.end(), network_offset
			);
			edge_index = static_cast<int>(edge_it - cumulative_lengths.begin());
			if (edge_index >= our_model.num_edges)
				edge_index = our_model.num_edges - 1;
			vector<Point>& PS_cur = our_model.edge_set[edge_index].PS_vec[l];
			uniform_real_distribution<double> distribution_double(0.0, our_model.edge_set[edge_index].length);
			uniform_real_distribution<double> distribution_double_time(our_model.min_time, our_model.max_time);
			pt.dist_n1 = distribution_double(generator_double);
			pt.dist_n2 = our_model.edge_set[edge_index].length - pt.dist_n1;
			pt.time = distribution_double_time(generator_double_time);
			pt.edge_index = edge_index;

			PS_cur.push_back(pt);
		}
	}
}

void sort_data_points(model& our_model)
{
	//Sorting the data points
	for (int e : our_model.active_source_edges)
		for (int l = 0; l < our_model.L; l++)
			sort(our_model.edge_set[e].PS_vec[l].begin(), our_model.edge_set[e].PS_vec[l].end());
}

void init_tree(Edge& e)
{
	int num_points;
	int L = (int)e.PS_vec.size();
	Range_Tree t;

	for (int l = 0; l < L; l++)
	{
		num_points = (int)e.PS_vec[l].size();
		e.t_vec.push_back(t);
		e.t_vec[l].num_points = num_points;
		e.t_vec[l].Point_Matrix_distance = new double*[num_points];
		for (int i = 0; i < num_points; i++)
		{
			e.t_vec[l].Point_Matrix_distance[i] = new double[2];
			e.t_vec[l].Point_Matrix_distance[i][0] = e.PS_vec[l][i].time;
			e.t_vec[l].Point_Matrix_distance[i][1] = e.PS_vec[l][i].dist_n1;
		}
	}
}

void init_prefix_count_array(model& our_model)
{
	our_model.prefix_count_one_D = new double[2 * our_model.T];
	our_model.prefix_count = new double*[2 * our_model.T];

	for (int t = 0; t < 2 * our_model.T; t++)
		our_model.prefix_count[t] = new double[our_model.M];

	clear_prefix_count_array(our_model, our_model.prefix_count, our_model.M);
	clear_prefix_count_one_D_array(our_model);
}

void construct_prefix_count_array(model& our_model, double**prefix_count, int d_max_index)
{
	double current_value;
	for (int t = 0; t < 2 * our_model.T; t++)
	{
		for (int d = 0; d < d_max_index; d++)
		{
			current_value = prefix_count[t][d];

			if (t != 0 && d != 0)
				prefix_count[t][d] = prefix_count[t - 1][d] + prefix_count[t][d - 1] - prefix_count[t - 1][d - 1] + current_value;
			if (t == 0 && d != 0)
				prefix_count[t][d] = prefix_count[t][d - 1] + current_value;
			if (t != 0 && d == 0)
				prefix_count[t][d] = prefix_count[t - 1][d] + current_value;
		}
	}
}

void clear_prefix_count_array(model& our_model, double**prefix_count, int d_max_index)
{
	for (int t = 0; t < 2 * our_model.T; t++)
		for (int d = 0; d < d_max_index; d++)
			prefix_count[t][d] = 0;
}

void construct_prefix_count_one_D_array(model& our_model)
{
	double current_value;
	for (int t = 0; t < 2 * our_model.T; t++)
	{
		current_value = our_model.prefix_count_one_D[t];
		if (t != 0)
			our_model.prefix_count_one_D[t] = our_model.prefix_count_one_D[t - 1] + current_value;
	}
}

void clear_prefix_count_one_D_array(model& our_model)
{
	for (int t = 0; t < 2 * our_model.T; t++)
		our_model.prefix_count_one_D[t] = 0;
}

namespace {

inline void free_matrix(double** matrix, int rows)
{
	if (matrix == nullptr)
		return;
	for (int i = 0; i < rows; ++i)
		delete[] matrix[i];
	delete[] matrix;
}

void free_time_tree_nodes(Range_Node* node)
{
	if (node == nullptr)
		return;
	if (node->index_u > node->index_l)
	{
		free_time_tree_nodes(node->left_child);
		free_time_tree_nodes(node->right_child);
	}
	delete node;
}

void free_base_tree(Range_Node* node)
{
	if (node == nullptr)
		return;
	if (node->index_u > node->index_l)
	{
		free_base_tree(node->left_child);
		free_base_tree(node->right_child);
	}

	Range_Tree& tt = node->time_tree;
	free_matrix(tt.Point_Matrix_time, tt.num_points);
	free_time_tree_nodes(tt.rootNode);

	delete node;
}

}

void kfun_cleanup(model& m) {
	delete[] m.s_tau_vec;
	m.s_tau_vec = nullptr;
	delete[] m.t_tau_vec;
	m.t_tau_vec = nullptr;

	if (m.multiple_K_values != nullptr)
	{
		for (int l = 0; l < m.L; ++l)
		{
			if (m.multiple_K_values[l] == nullptr)
				continue;
			for (int s = 0; s < m.M; ++s)
				delete[] m.multiple_K_values[l][s];
			delete[] m.multiple_K_values[l];
		}
		delete[] m.multiple_K_values;
		m.multiple_K_values = nullptr;
	}

	if (m.edge_set != nullptr)
	{
		for (int e = 0; e < m.num_edges; ++e)
		{
			Edge& edge = m.edge_set[e];
			for (Range_Tree& tree : edge.t_vec)
			{
				free_matrix(tree.Point_Matrix_distance, tree.num_points);
				if (tree.num_points > 0)
					free_base_tree(tree.rootNode);
				else
					delete tree.rootNode;
			}
		}
		delete[] m.edge_set;
		m.edge_set = nullptr;
	}

	delete[] m.is_scan_edge_list_first;
	m.is_scan_edge_list_first = nullptr;
	delete[] m.is_scan_edge_list_second;
	m.is_scan_edge_list_second = nullptr;
	delete[] m.detect_DoubleScan_edge_list;
	m.detect_DoubleScan_edge_list = nullptr;

	if (m.prefix_count != nullptr)
	{
		for (int t = 0; t < 2 * m.T; ++t)
			delete[] m.prefix_count[t];
		delete[] m.prefix_count;
		m.prefix_count = nullptr;
	}

	delete[] m.prefix_count_one_D;
	m.prefix_count_one_D = nullptr;

}

void kfun_init_parameters(
	model& m,
	const double* s_values,
	int n_s,
	const double* t_values,
	int n_t,
	int method,
		int simulations,
		unsigned int random_seed
) {
			m.random_seed = random_seed;
	m.method = method;
	m.M = n_s;
	m.T = n_t;
	m.L = simulations + 1;
	m.num_points = 0;
	m.progress_completed = 0;
	m.progress_total = 0;
	m.min_time = inf;
	m.max_time = -inf;

	m.start_s_tau = (n_s > 0) ? s_values[0] : 0.0;
	m.incr_s_tau = (n_s > 1) ? (s_values[1] - s_values[0]) : 0.0;
	m.start_t_tau = (n_t > 0) ? t_values[0] : 0.0;
	m.incr_t_tau = (n_t > 1) ? (t_values[1] - t_values[0]) : 0.0;

	m.s_tau_vec = new double[m.M];
	for (int s = 0; s < m.M; ++s)
		m.s_tau_vec[s] = s_values[s];

	m.t_tau_vec = new double[m.T];
	for (int t = 0; t < m.T; ++t)
		m.t_tau_vec[t] = t_values[t];

		m.multiple_K_values = new double**[m.L]();
	for (int l = 0; l < m.L; ++l)
	{
			m.multiple_K_values[l] = new double*[m.M]();
		for (int s = 0; s < m.M; ++s)
		{
			m.multiple_K_values[l][s] = new double[m.T];
			std::fill(m.multiple_K_values[l][s], m.multiple_K_values[l][s] + m.T, 0.0);
		}
	}
}

void kfun_build_model(
	model& m,
	int n_nodes,
	int n_edges,
	const std::vector<int>& edge_u,
	const std::vector<int>& edge_v,
	const std::vector<double>& edge_len,
	std::vector<int>& point_edge_ids,
	std::vector<double>& point_offsets,
	std::vector<double>& point_times
) {
	if (n_nodes <= 0 || n_edges <= 0) {
		throw std::invalid_argument("Network must contain nodes and edges.");
	}
	m.num_nodes = n_nodes;
	m.num_edges = n_edges;

	m.edge_set = new Edge[static_cast<std::size_t>(m.num_edges)];
	init_PS_memory_edges(m, m.num_edges);

	for (int e = 0; e < n_edges; ++e)
	{
		m.edge_set[e].n1 = edge_u[e];
		m.edge_set[e].n2 = edge_v[e];
		m.edge_set[e].length = edge_len[e];
	}

	const int n_points = static_cast<int>(point_edge_ids.size());
	const bool has_time = static_cast<int>(point_times.size()) == n_points;
	for (int i = 0; i < n_points; ++i)
	{
		const int e = point_edge_ids[i];
		const double d = point_offsets[i];
		const double tm = has_time ? point_times[i] : static_cast<double>(i);

		Point pt;
		pt.edge_index = e;
		pt.dist_n1 = d;
		pt.dist_n2 = m.edge_set[e].length - d;
		pt.time = tm;
		m.edge_set[e].PS_vec[0].push_back(pt);

		if (tm < m.min_time)
			m.min_time = tm;
		if (tm > m.max_time)
			m.max_time = tm;
	}

	m.num_points += n_points;
	if (n_points == 0)
	{
		m.min_time = 0.0;
		m.max_time = 0.0;
	}

	generate_points(m);

	m.active_source_edges.clear();
	m.active_source_edges.reserve(std::min(
		static_cast<size_t>(m.num_edges),
		static_cast<size_t>(m.num_points) * static_cast<size_t>(m.L)
	));
	for (int e = 0; e < m.num_edges; ++e)
	{
		bool has_points = false;
		for (int l = 0; l < m.L; ++l)
		{
			if (!m.edge_set[e].PS_vec[l].empty())
			{
				has_points = true;
				break;
			}
		}
		if (has_points)
			m.active_source_edges.push_back(e);
	}

	m.Network.clear();
	m.Network.resize(n_nodes);
	for (int e = 0; e < n_edges; ++e)
	{
		m.Network[m.edge_set[e].n1].push_back(e);
		m.Network[m.edge_set[e].n2].push_back(e);
	}
}





}

#pragma once
#ifndef FASTNSA_EARTH_NETWORK_H
#define FASTNSA_EARTH_NETWORK_H

#include "Library.h"
#include "earth/range_tree.h"

namespace fastnsa::earth{

const double inf = 999999999999999;
const double eps = 0.00000000001;

struct Point
{
	double dist_n1;
	double dist_n2;
	double time;
	int edge_index;

	bool operator<(const Point& p) const {
		return dist_n1 < p.dist_n1;
	}
};

struct Edge
{
	int n1;
	int n2;
	double length;
	vector<vector<Point>> PS_vec;

	//Used in the EAR method
	vector<Range_Tree> t_vec;

	//Used in the EARTH method
	double COST_EAR;
	double COST_MTS;
};

struct sp_node
{
	int node_index;
	double cur_sp_value;
	bool is_opt;
};

//This is the minimum heap for Dijkstra's algorithm
struct comparePriority
{
	bool operator()(sp_node& n1, sp_node& n2)
	{
		return n1.cur_sp_value > n2.cur_sp_value;
	}
};

struct model
{
	int num_nodes;
	int num_edges;
	char*network_fileName; //file Name
	int method; //method

	int L;
	int M;
	int T;
	int num_points;
	double start_s_tau;
	double incr_s_tau;
	double start_t_tau;
	double incr_t_tau;
	double*s_tau_vec;
	double*t_tau_vec;
	double***multiple_K_values; //L x M x T K-function values
	char*out_fileName;
	int cur_data_index;
	int cur_s_tau_index;
	int cur_t_tau_index;
	double min_time;
	double max_time;
		unsigned int random_seed;

	vector<vector<int>> Network; //road network
	Edge*edge_set; //edge set
	vector<int> active_source_edges; //edges containing observed or simulated points

	//Used in the Shortest path algorithm
	vector<sp_node> sp_node_vec_first;
	vector<sp_node> sp_node_vec_second;
	vector<int> access_edge_list_first; //store the edges that are within the distance s from the data point/node 
	vector<int> access_edge_list_second; //store the edges that are within the distance s from the node
	bool*is_scan_edge_list_first; //avoid to push the same edge in access_edge_list_first
	bool*is_scan_edge_list_second; //avoid to push the same edge in access_edge_list_second

	int cur_edge_index;
	int cur_node_index;
	int cur_point_index;

	//Used in SPS
	//********************************//
	bool is_first;
	bool*detect_DoubleScan_edge_list;
	//********************************//

	//Used in MTS
	double delta_s_threshold;
	double delta_t_threshold;
	double**prefix_count;
	double*prefix_count_one_D;
		function<bool(const char*, int, int)> progress;
	int progress_completed;
	int progress_total;
	//double**count_matrix;
};

void init_PS_memory_edges(model& our_model, int num_edges);
void generate_points(model& our_model);
void sort_data_points(model& our_model);
void init_tree(Edge& e);

//Prefix-count array
void init_prefix_count_array(model& our_model);
void construct_prefix_count_array(model& our_model, double**prefix_count, int d_max_index);
void clear_prefix_count_array(model& our_model, double**prefix_count, int d_max_index);
void construct_prefix_count_one_D_array(model& our_model);
void clear_prefix_count_one_D_array(model& our_model);
void kfun_cleanup(model& m);
void kfun_init_parameters(
	model& m,
	const double* s_values,
	int n_s,
	const double* t_values,
	int n_t,
	int method,
		int simulations,
		unsigned int random_seed = 5489u
);
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
);

}

#endif

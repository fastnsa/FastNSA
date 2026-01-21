
#ifndef NETWORK_H
#define NETWORK_H

#include "Library.h"

const double inf = 999999999999999;
const double eps = 0.00000000001;

struct Point
{
	double dist_n1;
	double dist_n2;
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
	// vector<vector<Point>> sam_PS_vec;//method 7 samping data
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

//Used in the multi-way merge sort, store the value of tau_u (or tau_v) and its corresponding position in the original array.
struct heapValue
{
    double val;
    
    //position
    int i;
    int j;
    heapValue(double val, int idx_i, int idx_j) : val(val), i(idx_i), j(idx_j) {}
};

//This is the minimum heap for k_merge
struct mergeCompare
{
	bool operator()(heapValue& v1, heapValue& v2)
	{
        return v1.val>v2.val;
	}
};


struct model
{
	int num_nodes;
	int num_edges;
	char*network_fileName; //file Name
	int method; //method

	int L;
	int D;
	int num_points;
	int num_sam_points;
	double start_tau;
	double incr_tau;
	double*tau_vec;
	double**multiple_K_values;
	char*out_fileName;
	int cur_data_index;
	int cur_tau_index;

	double runtime;//Record the running time of the algorithm

	//used in STARS method，
	double epsilon;
	bool is_sam;

    vector<priority_queue<heapValue,vector<heapValue>,mergeCompare>>minheap;//minheap[0] handle u,
    vector<Point> sam_edge; //store the sampling points of the edge

	//double tau; //tau
	//double K_value; //K-function value

	vector<vector<int>> Network; //road network
	Edge*edge_set; //edge set

	//Used in the Shortest path algorithm
	vector<sp_node> sp_node_vec_first;
	vector<sp_node> sp_node_vec_second;
	vector<int> access_edge_list_first; //store the edges that are within the distance tau from the data point/node 
	vector<int> access_edge_list_second; //store the edges that are within the distance tau from the node
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
};
void init_sampling_memory(model& our_model);
void init_PS_memory_edges(model& our_model, int num_edges);
void generate_points(model& our_model);
void sort_data_points(model& our_model);
// void sort_selection_sampling(model& our_model);
void sort_selection(model& our_model,int e_index);
void kfun_cleanup(model& m);
void kfun_init_parameters(model& m, const double* r_values, int n_r, int method, int simulations);
void kfun_build_model(model& m, 
    int n_nodes, int n_edges, 
    const std::vector<int>& edge_u,
    const std::vector<int>& edge_v,
    const std::vector<double>& edge_len,
    vector<int>& point_edge_ids,
    vector<double>& point_offsets
);
#endif

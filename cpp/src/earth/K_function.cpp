#include "earth/K_function.h"

namespace fastnsa::earth
{

void call_RQS(model& our_model)
{
	//Preprocessing
	int e_index;
	int node_u, node_v;
	our_model.multiple_K_values[our_model.cur_data_index][our_model.cur_s_tau_index][our_model.cur_t_tau_index] = 0;
	Edge*edge_set = our_model.edge_set;

	for (int wide_e_index : our_model.active_source_edges)
	{
		const char* stage = our_model.cur_data_index > 0 ? "Monte Carlo batch" : "backend compute";
		if (our_model.progress && our_model.progress(
			stage, our_model.progress_completed++, our_model.progress_total
		))
			throw runtime_error("FastNSA operation canceled.");
		our_model.cur_edge_index = wide_e_index;
		vector<Point>& PS = edge_set[wide_e_index].PS_vec[our_model.cur_data_index];
		for (int p = 0; p < (int)PS.size(); p++)
		{
			our_model.cur_point_index = p;
			Dijkstra(our_model, our_model.sp_node_vec_first, our_model.access_edge_list_first, our_model.is_scan_edge_list_first);
			for (int e = 0; e < (int)our_model.access_edge_list_first.size(); e++)
			{
				e_index = our_model.access_edge_list_first[e];

				if (wide_e_index != e_index)
				{
					node_u = edge_set[e_index].n1;
					node_v = edge_set[e_index].n2;
					if (our_model.sp_node_vec_first[node_u].cur_sp_value > our_model.s_tau_vec[our_model.cur_s_tau_index] &&
						our_model.sp_node_vec_first[node_v].cur_sp_value > our_model.s_tau_vec[our_model.cur_s_tau_index])
						continue;

					//call point-edge_count function
					our_model.multiple_K_values[our_model.cur_data_index][our_model.cur_s_tau_index][our_model.cur_t_tau_index] += point_edge_count(our_model, e_index, node_u, node_v);
				}
				else
					our_model.multiple_K_values[our_model.cur_data_index][our_model.cur_s_tau_index][our_model.cur_t_tau_index] += point_edge_count_same_edge(our_model, p, e_index);
			}

			clear_Dijkstra(edge_set, our_model.sp_node_vec_first, our_model.access_edge_list_first, our_model.is_scan_edge_list_first);
		}
	}
}

void call_SPS(model& our_model)
{
	//Preprocessing
	int e_index;
	int node_u, node_v;
	our_model.multiple_K_values[our_model.cur_data_index][our_model.cur_s_tau_index][our_model.cur_t_tau_index] = 0;
	Edge*edge_set = our_model.edge_set;
	//int edges_totalNum;
	int num_edges_two_sets;

	for (int wide_e_index : our_model.active_source_edges)
	{
		const char* stage = our_model.cur_data_index > 0 ? "Monte Carlo batch" : "backend compute";
		if (our_model.progress && our_model.progress(
			stage, our_model.progress_completed++, our_model.progress_total
		))
			throw runtime_error("FastNSA operation canceled.");
		our_model.cur_edge_index = wide_e_index;
		vector<Point>& PS = edge_set[wide_e_index].PS_vec[our_model.cur_data_index];
		if (PS.size() == 0)
			continue;

		our_model.is_first = true;
		Dijkstra(our_model, our_model.sp_node_vec_first, our_model.access_edge_list_first, our_model.is_scan_edge_list_first);
		our_model.is_first = false;
		Dijkstra(our_model, our_model.sp_node_vec_second, our_model.access_edge_list_second, our_model.is_scan_edge_list_second);
		num_edges_two_sets = (int)our_model.access_edge_list_first.size() + (int)our_model.access_edge_list_second.size();

		for (int e = 0; e < num_edges_two_sets; e++)
		{
			if (e < (int)our_model.access_edge_list_first.size())
			{
				e_index = our_model.access_edge_list_first[e];
				our_model.detect_DoubleScan_edge_list[e_index] = true;
			}
			else
			{
				e_index = our_model.access_edge_list_second[e - our_model.access_edge_list_first.size()];
				if (our_model.detect_DoubleScan_edge_list[e_index] == true)
					continue;
				else
					our_model.detect_DoubleScan_edge_list[e_index] = true;
			}

			if (wide_e_index != e_index)
			{
				node_u = edge_set[e_index].n1;
				node_v = edge_set[e_index].n2;
				if (our_model.sp_node_vec_first[node_u].cur_sp_value > our_model.s_tau_vec[our_model.cur_s_tau_index] &&
					our_model.sp_node_vec_first[node_v].cur_sp_value > our_model.s_tau_vec[our_model.cur_s_tau_index] &&
					our_model.sp_node_vec_second[node_u].cur_sp_value > our_model.s_tau_vec[our_model.cur_s_tau_index] &&
					our_model.sp_node_vec_second[node_v].cur_sp_value > our_model.s_tau_vec[our_model.cur_s_tau_index])
					continue;

				for (int p = 0; p < (int)PS.size(); p++)
				{
					our_model.cur_point_index = p;
					our_model.multiple_K_values[our_model.cur_data_index][our_model.cur_s_tau_index][our_model.cur_t_tau_index] += point_edge_count(our_model, e_index, node_u, node_v);
				}
			}
			else
			{
				for (int p = 0; p < (int)PS.size(); p++)
					our_model.multiple_K_values[our_model.cur_data_index][our_model.cur_s_tau_index][our_model.cur_t_tau_index] += point_edge_count_same_edge(our_model, p, e_index);
			}
		}

		clear_double_scan_edgeList(our_model.access_edge_list_first, our_model.detect_DoubleScan_edge_list);
		clear_Dijkstra(edge_set, our_model.sp_node_vec_first, our_model.access_edge_list_first, our_model.is_scan_edge_list_first);
		clear_double_scan_edgeList(our_model.access_edge_list_second, our_model.detect_DoubleScan_edge_list);
		clear_Dijkstra(edge_set, our_model.sp_node_vec_second, our_model.access_edge_list_second, our_model.is_scan_edge_list_second);
	}
}

void call_EAR_edge(model& our_model, int wide_e_index, int e_index, int node_u, int node_v)
{
	Edge*edge_set = our_model.edge_set;

	for (int l = 0; l < our_model.L; l++)
	{
		our_model.cur_data_index = l;
		vector<Point>& PS = edge_set[wide_e_index].PS_vec[l];

		if (PS.size() == 0)
			continue;

		for (int s = our_model.M - 1; s >= 0; s--)
		{
			our_model.cur_s_tau_index = s;

			if (wide_e_index != e_index)
			{
				if (our_model.sp_node_vec_first[node_u].cur_sp_value > our_model.s_tau_vec[s] &&
					our_model.sp_node_vec_first[node_v].cur_sp_value > our_model.s_tau_vec[s] &&
					our_model.sp_node_vec_second[node_u].cur_sp_value > our_model.s_tau_vec[s] &&
					our_model.sp_node_vec_second[node_v].cur_sp_value > our_model.s_tau_vec[s])
					break;
			}

			for (int t = our_model.T - 1; t >= 0; t--)
			{
				our_model.cur_t_tau_index = t;
				double edge_count_value;
				if (wide_e_index != e_index)
				{
					edge_count_value = 0;
					for (int p = 0; p < (int)PS.size(); p++)
					{
						our_model.cur_point_index = p;
						edge_count_value += point_edge_count(our_model, e_index, node_u, node_v);
					}
				}
				else
				{
					edge_count_value = 0;
					for (int p = 0; p < (int)PS.size(); p++)
						edge_count_value += point_edge_count_same_edge(our_model, p, e_index);
				}

				our_model.multiple_K_values[l][s][t] += edge_count_value;

				if (edge_count_value < 0.1)
					break;
			}
		}
	}
}

void call_MTS_edge(model& our_model, int wide_e_index, int e_index, int node_u, int node_v)
{
	Edge*edge_set = our_model.edge_set;

	for (int l = 0; l < our_model.L; l++)
	{
		our_model.cur_data_index = l;
		vector<Point>& PS = edge_set[wide_e_index].PS_vec[l];

		if (PS.size() == 0)
			continue;

		if (wide_e_index != e_index)
		{
			for (int p = 0; p < (int)PS.size(); p++)
			{
				our_model.cur_point_index = p;
				point_edge_count_MTS(our_model, e_index, node_u, node_v);
			}
		}
		else
		{
			for (int p = 0; p < (int)PS.size(); p++)
			{
				our_model.cur_point_index = p;
				point_edge_count_same_edge_MTS(our_model, e_index, node_u, node_v);
			}
		}
	}
}

void call_ASPS(model& our_model)
{
	int e_index;
	int node_u, node_v;
	Edge*edge_set = our_model.edge_set;
	int num_edges_two_sets;
	our_model.progress_completed = 0;
	our_model.progress_total = (int)our_model.active_source_edges.size();

	if (our_model.progress && our_model.progress(
		"backend compute", our_model.progress_completed, our_model.progress_total
	))
		throw runtime_error("FastNSA operation canceled.");

	if (our_model.method == 3) //EAR method (preprocessing step)
	{
		for (int e : our_model.active_source_edges)
		{
			init_tree(our_model.edge_set[e]);
			for (int l = 0; l < our_model.L; l++)
				our_model.edge_set[e].t_vec[l].build_range_tree();
		}
	}

	if (our_model.method == 4) //MTS method (preprocessing step)
	{
		our_model.delta_s_threshold = our_model.s_tau_vec[1] - our_model.s_tau_vec[0];
		our_model.delta_t_threshold = our_model.t_tau_vec[1] - our_model.t_tau_vec[0];
		init_prefix_count_array(our_model);
	}

	if (our_model.method == 5) //EARTH method (preprocessing step)
	{
		int edge_set_size;
		our_model.delta_s_threshold = our_model.s_tau_vec[1] - our_model.s_tau_vec[0];
		our_model.delta_t_threshold = our_model.t_tau_vec[1] - our_model.t_tau_vec[0];
		init_prefix_count_array(our_model);

		for (int e : our_model.active_source_edges)
		{
			edge_set[e].COST_EAR = 0;
			edge_set[e].COST_MTS = 0;
			for (int l = 0; l < our_model.L; l++)
			{
				edge_set_size = (int)edge_set[e].PS_vec[l].size();
				if (edge_set_size == 0)
					continue;
				else
				{
					edge_set[e].COST_EAR += edge_set_size * log(edge_set_size) + our_model.num_points * our_model.T * our_model.M * log(edge_set_size);
					edge_set[e].COST_MTS += our_model.num_points * (our_model.T*our_model.M + edge_set_size);
				}
			}

			if (edge_set[e].COST_EAR <= edge_set[e].COST_MTS || our_model.M == 1 || our_model.T == 1)
			{
				init_tree(our_model.edge_set[e]);
				for (int l = 0; l < our_model.L; l++)
					our_model.edge_set[e].t_vec[l].build_range_tree();
			}
		}
	}

	for (int wide_e_index : our_model.active_source_edges)
	{
		const char* stage = our_model.L > 1 ? "Monte Carlo batch" : "backend compute";
		if (our_model.progress && our_model.progress(
			stage, our_model.progress_completed++, our_model.progress_total
		))
			throw runtime_error("FastNSA operation canceled.");
		our_model.cur_edge_index = wide_e_index;
		our_model.cur_s_tau_index = our_model.M - 1;

		our_model.is_first = true;
		Dijkstra(our_model, our_model.sp_node_vec_first, our_model.access_edge_list_first, our_model.is_scan_edge_list_first);
		our_model.is_first = false;
		Dijkstra(our_model, our_model.sp_node_vec_second, our_model.access_edge_list_second, our_model.is_scan_edge_list_second);
		num_edges_two_sets = (int)our_model.access_edge_list_first.size() + (int)our_model.access_edge_list_second.size();

		for (int e = 0; e < num_edges_two_sets; e++)
		{
			if (e < (int)our_model.access_edge_list_first.size())
			{
				e_index = our_model.access_edge_list_first[e];
				our_model.detect_DoubleScan_edge_list[e_index] = true;
			}
			else
			{
				e_index = our_model.access_edge_list_second[e - our_model.access_edge_list_first.size()];
				if (our_model.detect_DoubleScan_edge_list[e_index] == true)
					continue;
				else
					our_model.detect_DoubleScan_edge_list[e_index] = true;
			}

			node_u = edge_set[e_index].n1;
			node_v = edge_set[e_index].n2;

			if (our_model.method == 3)
				call_EAR_edge(our_model, wide_e_index, e_index, node_u, node_v);
			if (our_model.method == 4)
				call_MTS_edge(our_model, wide_e_index, e_index, node_u, node_v);
			if (our_model.method == 5)
			{
				if (edge_set[e_index].COST_EAR <= edge_set[e_index].COST_MTS || our_model.M == 1 || our_model.T == 1)
					call_EAR_edge(our_model, wide_e_index, e_index, node_u, node_v);
				else
					call_MTS_edge(our_model, wide_e_index, e_index, node_u, node_v);
			}
		}

		clear_double_scan_edgeList(our_model.access_edge_list_first, our_model.detect_DoubleScan_edge_list);
		clear_Dijkstra(edge_set, our_model.sp_node_vec_first, our_model.access_edge_list_first, our_model.is_scan_edge_list_first);
		clear_double_scan_edgeList(our_model.access_edge_list_second, our_model.detect_DoubleScan_edge_list);
		clear_Dijkstra(edge_set, our_model.sp_node_vec_second, our_model.access_edge_list_second, our_model.is_scan_edge_list_second);
	}
}

void compute_K_function(model& our_model)
{
	double run_time;

	sort_data_points(our_model);

	auto start_s = chrono::high_resolution_clock::now();
	// Share one workspace across all thresholds and simulated point sets.
	Dijkstra_preprocess(our_model);

	if (our_model.method >= 1 && our_model.method <= 2)
	{
		our_model.progress_completed = 0;
		our_model.progress_total = our_model.L * our_model.M * our_model.T
			* (int)our_model.active_source_edges.size();
		for (int l = 0; l < our_model.L; l++)
		{
			const char* stage = l > 0 ? "Monte Carlo batch" : "backend compute";
			if (our_model.progress && our_model.progress(
				stage, our_model.progress_completed, our_model.progress_total
			))
				throw runtime_error("FastNSA operation canceled.");
			our_model.cur_data_index = l;
			for (int s = 0; s < our_model.M; s++)
			{
				our_model.cur_s_tau_index = s;
				for (int t = 0; t < our_model.T; t++)
				{
					our_model.cur_t_tau_index = t;

					if (our_model.method == 1) //RQS
						call_RQS(our_model);
					if (our_model.method == 2) //SPS
						call_SPS(our_model);
				}
			}
		}
	}

	if (our_model.method == 3) //EAR (+ ASPS)
		call_ASPS(our_model);
	if (our_model.method == 4) //MTS (+ ASPS)
		call_ASPS(our_model);
	if (our_model.method == 5) //EARTH (+ ASPS)
		call_ASPS(our_model);

	auto end_s = chrono::high_resolution_clock::now();
	run_time = (chrono::duration_cast<chrono::nanoseconds>(end_s - start_s).count()) / 1000000000.0;
	(void)run_time;
}

}

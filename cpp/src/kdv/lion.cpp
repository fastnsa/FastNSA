#include "kdv/lion.h"

namespace fastnsa::kdv {

void init_LION(model& our_model)
{
	int num_deg;

	if (our_model.k_type == 1) //Triangular kernel
		our_model.num_deg = 2;
	if (our_model.k_type == 2) //Epanechnikov kernel
		our_model.num_deg = 3;
	if (our_model.k_type == 3) //Quartic kernel
		our_model.num_deg = 5;

	num_deg = our_model.num_deg;
	for (int e = 0; e < our_model.m; e++)
	{
		if (our_model.progress && our_model.progress("backend compute", e, our_model.m * 4)) {
			our_model.canceled = true;
			return;
		}
		Edge& edge = our_model.edge_set[e];

		for (int q = 0; q < (int)edge.lixel_vec.size(); q++)
		{
			//degree = 3 (0, 1, and 2)
			edge.lixel_vec[q].alpha_B_x_deg_vec = new double[num_deg];
			edge.lixel_vec[q].alpha_B_y_deg_vec = new double[num_deg];
			edge.lixel_vec[q].alpha_R_x_deg_vec = new double[num_deg];
			edge.lixel_vec[q].alpha_R_y_deg_vec = new double[num_deg];
		}
	}

	our_model.ell.a_L = new double[our_model.num_deg];
	our_model.ell.a_R = new double[our_model.num_deg];
}

void lixel_augmentation(model& our_model)
{
	int num_deg = our_model.num_deg;
	int e_index;
	int node_index_pi_a; int node_index_pi_b;
	int node_index_pi_u; int node_index_pi_v;
	double temp_value;
	double dist_a_p;
	double dist_b_p;

	our_model.lixel_total_number = 0;
	for (int wide_e_index = 0; wide_e_index < our_model.m; wide_e_index++)
	{
		if (our_model.progress && our_model.progress("backend compute", our_model.m + wide_e_index, our_model.m * 4)) {
			our_model.canceled = true;
			return;
		}
		Edge& wide_edge = our_model.edge_set[wide_e_index];
		our_model.lixel_total_number += (int)wide_edge.lixel_vec.size();

		for (int q = 0; q < (int)wide_edge.lixel_vec.size(); q++)
		{
			for (int d = 0; d < num_deg; d++)
			{
				wide_edge.lixel_vec[q].alpha_B_x_deg_vec[d] = 0;
				wide_edge.lixel_vec[q].alpha_B_y_deg_vec[d] = 0;
				wide_edge.lixel_vec[q].alpha_R_x_deg_vec[d] = 0;
				wide_edge.lixel_vec[q].alpha_R_y_deg_vec[d] = 0;
			}
		}

		our_model.cur_edge_index = wide_e_index;

		//Obtain SPD from node x (node a) to all other nodes.
		clear_scan_edge_list(our_model);
		our_model.sel_node_index = wide_edge.n1;
		our_model.node_index_pi_a = our_model.sel_node_index;
		dijkstra(our_model);
		if (our_model.canceled) return;
		copy_sp_info(our_model, true);

		//Obtain SPD from node y (node b) to all other nodes.
		our_model.sel_node_index = wide_edge.n2;
		our_model.node_index_pi_b = our_model.sel_node_index;
		dijkstra(our_model);
		if (our_model.canceled) return;
		copy_sp_info(our_model, false);
		wide_edge.endpoint_shortest = our_model.sp_node_vec_node_a[wide_edge.n2].cur_sp_value;

		node_index_pi_a = our_model.node_index_pi_a;
		node_index_pi_b = our_model.node_index_pi_b;
		for (int e = 0; e < (int)our_model.access_edge_list.size(); e++)
		{
			e_index = our_model.access_edge_list[e];

			if (e_index == wide_e_index)
				continue;
			
			Edge& edge = our_model.edge_set[e_index];
			node_index_pi_u = edge.n1;
			node_index_pi_v = edge.n2;

			for (int p = 0; p < (int)edge.PS.size(); p++)
			{
				dist_a_p = min(our_model.sp_node_vec_node_a[node_index_pi_u].cur_sp_value + edge.PS[p].dist_n1,
					our_model.sp_node_vec_node_a[node_index_pi_v].cur_sp_value + edge.PS[p].dist_n2);
				dist_b_p = min(our_model.sp_node_vec_node_b[node_index_pi_u].cur_sp_value + edge.PS[p].dist_n1,
					our_model.sp_node_vec_node_b[node_index_pi_v].cur_sp_value + edge.PS[p].dist_n2);

				// Use actual lixel centres. The former regular-grid index arithmetic
				// was off by one near a short final lixel and could over-subtract.
				const int lixel_count = static_cast<int>(wide_edge.lixel_vec.size());
				const double switch_offset = (wide_edge.length - dist_a_p + dist_b_p) / 2.0;
				int last_left = -1;
				if (dist_a_p <= our_model.bandwidth) {
					const double left_limit = min(our_model.bandwidth - dist_a_p, switch_offset);
					last_left = static_cast<int>(upper_bound(
						wide_edge.lixel_vec.begin(), wide_edge.lixel_vec.end(), left_limit + eps,
						[](double value, const Lixel& lixel) { return value < lixel.dist_n1; })
						- wide_edge.lixel_vec.begin()) - 1;
				}

				int first_right = lixel_count;
				if (dist_b_p <= our_model.bandwidth) {
					const double support_start = wide_edge.length - (our_model.bandwidth - dist_b_p);
					const int by_support = static_cast<int>(lower_bound(
						wide_edge.lixel_vec.begin(), wide_edge.lixel_vec.end(), support_start - eps,
						[](const Lixel& lixel, double value) { return lixel.dist_n1 < value; })
						- wide_edge.lixel_vec.begin());
					const int after_switch = static_cast<int>(upper_bound(
						wide_edge.lixel_vec.begin(), wide_edge.lixel_vec.end(), switch_offset + eps,
						[](double value, const Lixel& lixel) { return value < lixel.dist_n1; })
						- wide_edge.lixel_vec.begin());
					first_right = max(by_support, after_switch);
				}

				if (last_left >= 0) {
					temp_value = edge.PS[p].weight;
					for (int d = 0; d < our_model.num_deg; d++) {
						wide_edge.lixel_vec[last_left].alpha_B_x_deg_vec[d] += temp_value;
						temp_value *= dist_a_p;
					}
				}
				if (first_right < lixel_count) {
					temp_value = edge.PS[p].weight;
					for (int d = 0; d < our_model.num_deg; d++) {
						wide_edge.lixel_vec[first_right].alpha_B_y_deg_vec[d] += temp_value;
						temp_value *= dist_b_p;
					}
				}

			}
		}
	}
}

void lixel_aggregation(model& our_model)
{
	double*alpha_left_vec = new double[our_model.num_deg];
	double*alpha_right_vec = new double[our_model.num_deg];

	for (int wide_e_index = 0; wide_e_index < our_model.m; wide_e_index++)
	{
		if (our_model.progress && our_model.progress("backend compute", our_model.m * 2 + wide_e_index, our_model.m * 4)) {
			our_model.canceled = true;
			delete[] alpha_left_vec;
			delete[] alpha_right_vec;
			return;
		}
		Edge& wide_edge = our_model.edge_set[wide_e_index];
		for (int d = 0; d < our_model.num_deg; d++)
		{
			alpha_left_vec[d] = 0;
			alpha_right_vec[d] = 0;
		}

		for (int q = wide_edge.lixel_vec.size() - 1; q >= 0; q--)
		{
			for (int d = 0; d < our_model.num_deg; d++)
			{
				alpha_left_vec[d] += wide_edge.lixel_vec[q].alpha_B_x_deg_vec[d];
				wide_edge.lixel_vec[q].alpha_R_x_deg_vec[d] = alpha_left_vec[d];
			}
		}

		for (int q = 0; q < (int)wide_edge.lixel_vec.size(); q++)
		{
			for (int d = 0; d < our_model.num_deg; d++)
			{
				alpha_right_vec[d] += wide_edge.lixel_vec[q].alpha_B_y_deg_vec[d];
				wide_edge.lixel_vec[q].alpha_R_y_deg_vec[d] = alpha_right_vec[d];
			}
		}

		for (int q = 0; q < (int)wide_edge.lixel_vec.size(); q++)
		{
			if (our_model.k_type == 2) //Epanechnikov kernel
			{
				wide_edge.lixel_vec[q].KDE_value = (1 - (wide_edge.lixel_vec[q].dist_n1 * wide_edge.lixel_vec[q].dist_n1) / (our_model.bandwidth * our_model.bandwidth)) * wide_edge.lixel_vec[q].alpha_R_x_deg_vec[0]
					- (2.0 * wide_edge.lixel_vec[q].dist_n1 * wide_edge.lixel_vec[q].alpha_R_x_deg_vec[1]) / (our_model.bandwidth*our_model.bandwidth)
					- (wide_edge.lixel_vec[q].alpha_R_x_deg_vec[2] / (our_model.bandwidth * our_model.bandwidth))
					+ (1 - (wide_edge.lixel_vec[q].dist_n2 * wide_edge.lixel_vec[q].dist_n2) / (our_model.bandwidth * our_model.bandwidth)) * wide_edge.lixel_vec[q].alpha_R_y_deg_vec[0]
					- (2.0 * wide_edge.lixel_vec[q].dist_n2 * wide_edge.lixel_vec[q].alpha_R_y_deg_vec[1]) / (our_model.bandwidth*our_model.bandwidth)
					- (wide_edge.lixel_vec[q].alpha_R_y_deg_vec[2] / (our_model.bandwidth * our_model.bandwidth));
			}

			//Other kernels (Code here)
		}
	}

	delete[] alpha_left_vec;
	delete[] alpha_right_vec;
}

void init_one_D_KDV(model& our_model)
{
	for (int d = 0; d < our_model.num_deg; d++)
	{
		our_model.ell.a_L[d] = 0;
		our_model.ell.a_R[d] = 0;
	}
}

void obtain_sorted_end_vec(model& our_model, Edge& edge, vector<end_point>& sorted_end_vec)
{
	vector<end_point> end_vec_left;
	vector<end_point> end_vec_right;
	end_point temp_e_point;
	int end_point_num = 2 * edge.PS.size();
	int left_num = 0;
	int right_num = 0;

	for (int p = 0; p < (int)edge.PS.size(); p++)
	{
		temp_e_point.index = p;
		temp_e_point.isLeft = true;
		temp_e_point.value = edge.PS[p].dist_n1 - our_model.bandwidth;
		
		end_vec_left.push_back(temp_e_point);

		temp_e_point.isLeft = false;
		temp_e_point.value = edge.PS[p].dist_n1 + our_model.bandwidth;
		
		end_vec_right.push_back(temp_e_point);
	}

	for (int i = 0; i < end_point_num; i++)
	{
		if (left_num == (int)edge.PS.size())
		{
			sorted_end_vec.push_back(end_vec_right[right_num]);
			right_num++;
			continue;
		}
		
		if (right_num == (int)edge.PS.size())
		{
			sorted_end_vec.push_back(end_vec_left[left_num]);
			left_num++;
			continue;
		}

		if (end_vec_left[left_num].value < end_vec_right[right_num].value)
		{
			sorted_end_vec.push_back(end_vec_left[left_num]);
			left_num++;
		}
		else
		{
			sorted_end_vec.push_back(end_vec_right[right_num]);
			right_num++;
		}
	}
}

void sweep_line_algorithm(model& our_model, Edge& edge, vector<end_point>& sorted_end_vec)
{
	int lixel_size = edge.lixel_vec.size();
	int end_point_num = 2 * edge.PS.size();
	int lixel_index = 0;
	int end_point_index = 0;
	double temp_value;

	for (int i = 0; i < lixel_size + end_point_num; i++)
	{
		if (lixel_index == lixel_size || end_point_index == end_point_num)
			break;

		//if (edge.lixel_vec[lixel_index].dist_n1 < edge.PS[sorted_end_vec[end_point_index].index].dist_n1)
		if (edge.lixel_vec[lixel_index].dist_n1 < sorted_end_vec[end_point_index].value)
		{
			if (our_model.k_type == 2) //Epanechnikov kernel
			{
				edge.lixel_vec[lixel_index].KDE_value += (1.0 - (edge.lixel_vec[lixel_index].dist_n1 * edge.lixel_vec[lixel_index].dist_n1) / (our_model.bandwidth * our_model.bandwidth)) * our_model.ell.a_L[0]
					+ (2.0 * edge.lixel_vec[lixel_index].dist_n1 * our_model.ell.a_L[1]) / (our_model.bandwidth * our_model.bandwidth)
					- our_model.ell.a_L[2] / (our_model.bandwidth * our_model.bandwidth)
					- (1.0 - (edge.lixel_vec[lixel_index].dist_n1 * edge.lixel_vec[lixel_index].dist_n1) / (our_model.bandwidth * our_model.bandwidth)) * our_model.ell.a_R[0]
					- (2.0 * edge.lixel_vec[lixel_index].dist_n1 * our_model.ell.a_R[1]) / (our_model.bandwidth * our_model.bandwidth)
					+ our_model.ell.a_R[2] / (our_model.bandwidth * our_model.bandwidth);
			}

			//Other kernels (Code here)

			lixel_index++;
		}
		else
		{
			if (sorted_end_vec[end_point_index].isLeft == true)
			{
				temp_value = edge.PS[sorted_end_vec[end_point_index].index].weight;
				for (int d = 0; d < our_model.num_deg; d++)
				{
					our_model.ell.a_L[d] += temp_value;
					temp_value *= edge.PS[sorted_end_vec[end_point_index].index].dist_n1;
				}
			}
			else
			{
				temp_value = edge.PS[sorted_end_vec[end_point_index].index].weight;
				for (int d = 0; d < our_model.num_deg; d++)
				{
					our_model.ell.a_R[d] += temp_value;
					temp_value *= edge.PS[sorted_end_vec[end_point_index].index].dist_n1;
				}
			}

			end_point_index++;
		}
	}
}

void one_D_KDV(model& our_model) //Section 3.1
{
	vector<end_point> sorted_end_vec;

	for (int wide_e = 0; wide_e < our_model.m; wide_e++)
	{
		if (our_model.progress && our_model.progress("backend compute", our_model.m * 3 + wide_e, our_model.m * 4)) {
			our_model.canceled = true;
			return;
		}
		init_one_D_KDV(our_model);
		Edge& edge = our_model.edge_set[wide_e];

		if (edge.PS.size() == 0)
			continue;

		obtain_sorted_end_vec(our_model, edge, sorted_end_vec);
		sweep_line_algorithm(our_model, edge, sorted_end_vec);
		sorted_end_vec.clear();
		if (edge.endpoint_shortest < edge.length - eps && edge.endpoint_shortest <= our_model.bandwidth)
		{
			for (Lixel& lixel : edge.lixel_vec)
				correct_same_edge_shortcuts(our_model, lixel,
					min(lixel.dist_n1, lixel.dist_n2 + edge.endpoint_shortest),
					min(lixel.dist_n2, lixel.dist_n1 + edge.endpoint_shortest));
		}
	}
}


}

#include "kdv/shortest_path.h"
namespace fastnsa::kdv {
inline void clearHeap(PQ& pq)
{
	int heapSize = (int)pq.size();
	for (int h = 0; h < heapSize; h++)
		pq.pop();
}

void init_dijkstra(model& our_model, PQ& pq)
{
	int n1, n2;
	Lixel& cur_l = our_model.cur_l;

	for (int v = 0; v < our_model.n; v++)
	{
		our_model.sp_node_vec[v].node_index = v;
		our_model.sp_node_vec[v].cur_sp_value = inf;
		our_model.sp_node_vec[v].is_opt = false;
	}

	if (our_model.method >= 1 && our_model.method <= 5)
	{
		if (our_model.is_scan_edge_list[cur_l.edge_index] == false)
		{
			our_model.is_scan_edge_list[cur_l.edge_index] = true;
			our_model.access_edge_list.push_back(cur_l.edge_index);
		}
	}
	
	if (our_model.method == 6)
	{
		if (our_model.is_scan_edge_list[our_model.cur_edge_index] == false)
		{
			our_model.is_scan_edge_list[our_model.cur_edge_index] = true;
			our_model.access_edge_list.push_back(our_model.cur_edge_index);
		}
	}

	if (our_model.method == 1) //basic method
	{
		n1 = our_model.edge_set[cur_l.edge_index].n1;
		n2 = our_model.edge_set[cur_l.edge_index].n2;
		our_model.sp_node_vec[n1].cur_sp_value = cur_l.dist_n1;
		pq.push(our_model.sp_node_vec[n1]);
		if (cur_l.dist_n2 < our_model.sp_node_vec[n2].cur_sp_value)
		{
			our_model.sp_node_vec[n2].cur_sp_value = cur_l.dist_n2;
			pq.push(our_model.sp_node_vec[n2]);
		}
	}
	if (our_model.method >= 2 && our_model.method <= 6) //share SPD value
	{
		our_model.sp_node_vec[our_model.sel_node_index].cur_sp_value = 0;
		pq.push(our_model.sp_node_vec[our_model.sel_node_index]);
	}
}

void dijkstra(model& our_model)
{
	// Each task owns its queue, including any entries left on cancellation.
	PQ pq;
	int node_index, node_index_neighbor;
	int edge_index;
	double length;
	sp_node SP_node;

	init_dijkstra(our_model, pq);

	std::size_t visited = 0;
	while (pq.size() > 0)
	{
		if ((visited++ & 255u) == 0 && our_model.progress &&
			our_model.progress("backend compute", 0, 1)) {
			our_model.canceled = true;
			return;
		}
		SP_node = pq.top();
		pq.pop();

		//discard outdated node
		if (our_model.sp_node_vec[SP_node.node_index].cur_sp_value < SP_node.cur_sp_value - eps)
			continue;

		//discard node with value more than the bandwidth
		if (our_model.sp_node_vec[SP_node.node_index].cur_sp_value > our_model.bandwidth)
			continue;

		our_model.sp_node_vec[SP_node.node_index].is_opt = true;
		node_index = SP_node.node_index;

		for (int e = 0; e < (int)our_model.Network[node_index].size(); e++)
		{
			//Relax (code here)
			edge_index = our_model.Network[node_index][e];
			if (our_model.is_scan_edge_list[edge_index] == false)
			{
				our_model.is_scan_edge_list[edge_index] = true;
				our_model.access_edge_list.push_back(edge_index);
			}

			if (node_index != our_model.edge_set[edge_index].n1)
				node_index_neighbor = our_model.edge_set[edge_index].n1;
			else
				node_index_neighbor = our_model.edge_set[edge_index].n2;

			if (our_model.sp_node_vec[node_index_neighbor].is_opt == true)
				continue;

			length = our_model.edge_set[edge_index].length;
			if (our_model.sp_node_vec[node_index_neighbor].cur_sp_value > our_model.sp_node_vec[node_index].cur_sp_value + length)
			{
				our_model.sp_node_vec[node_index_neighbor].cur_sp_value = our_model.sp_node_vec[node_index].cur_sp_value + length;
				pq.push(our_model.sp_node_vec[node_index_neighbor]);
			}
		}
	}
}

void copy_sp_info(model& our_model, bool isFirst)
{
	for (int v = 0; v < our_model.n; v++)
	{
		if (isFirst == true)
		{
			our_model.sp_node_vec_node_a[v].node_index = our_model.sp_node_vec[v].node_index;
			our_model.sp_node_vec_node_a[v].is_opt = our_model.sp_node_vec[v].is_opt;
			our_model.sp_node_vec_node_a[v].cur_sp_value = our_model.sp_node_vec[v].cur_sp_value;
		}
		else
		{
			our_model.sp_node_vec_node_b[v].node_index = our_model.sp_node_vec[v].node_index;
			our_model.sp_node_vec_node_b[v].is_opt = our_model.sp_node_vec[v].is_opt;
			our_model.sp_node_vec_node_b[v].cur_sp_value = our_model.sp_node_vec[v].cur_sp_value;
		}
	}
}

void clear_scan_edge_list(model& our_model)
{
	int e_index;
	for (int e = 0; e < (int)our_model.access_edge_list.size(); e++)
	{
		e_index = our_model.access_edge_list[e];
		our_model.is_scan_edge_list[e_index] = false;
	}
	
	our_model.access_edge_list.clear();
}


}

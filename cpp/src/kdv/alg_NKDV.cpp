#include "kdv/alg_NKDV.h"
namespace fastnsa::kdv {
void add_lixels_for_edge(int e, model& our_model)
{
    double cur_dist = 0.0;
    const double length = our_model.edge_set[e].length;

    while (cur_dist < length)
    {
        double next_dist = cur_dist + our_model.lixel_reg_length;
        if (next_dist > length)
            next_dist = length;

        Lixel l;
        const double middle_dist = (cur_dist + next_dist) / 2.0;
        l.dist_n1 = middle_dist;
        l.dist_n2 = length - middle_dist;
        l.edge_index = e;
        l.KDE_value = -100;

        if (our_model.method >= 1 && our_model.method <= 5)
            our_model.lixel_set.push_back(l);
        else
            our_model.edge_set[e].lixel_vec.push_back(l);

        cur_dist += our_model.lixel_reg_length;
    }
}

void obtain_lixel_set(model& our_model)
{
    for (int e = 0; e < our_model.m; e++)
        add_lixels_for_edge(e, our_model);
}

}

#pragma once

#include "core/problem.hpp"
#include "tntp.hpp"

namespace transport
{

using namespace t_opt;

struct SDM : public Problem
{
    SDM(const String& problem_name, const String& data_path, const String& data_name);

    void
    apply_flow(double k, bool use_max = true);

    void
    apply_total_flow(double k);

protected:
    inline double&
    T(DVector& v, uint32_t i, uint32_t j)
    {
        // i and j are come from edges which has 1 as start value, but we (C++) uses 0 as start index
        return v[(i - 1) * data.max_node_index + j - 1];
    }

    inline const double&
    T(const DVector& v, uint32_t i, uint32_t j)
    {
        // i and j are come from edges which has 1 as start value, but we (C++) uses 0 as start index
        return v[(i - 1) * data.max_node_index + j - 1];
    }

    std::pair<double, double>
    calc_exp_sum(const DVector& x, const tntp::Edge& edge, double mu);

    tntp::Data data;
    DVector work;
};

}

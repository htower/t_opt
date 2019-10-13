#include "sdm.hpp"

#include "core/chrono.hpp"

namespace transport
{

SDM::SDM(const String& problem_name, const String& data_path, const String& data_name)
    : Problem(problem_name, 0, ProblemProperty::Gradient) // pass 0 as size, real size will be calculated later
{
    auto time_0 = chrono::now();
    tntp::load_tntp_data(data_path, data_name, data);
    auto time_i = chrono::s(time_0);
    printf("Load time = %.4f s.\n", time_i);

    // FIXME we need sparse problem instead
    m_size = data.max_node_index * data.max_node_index;

    // FIXME

    printf("total flow = %e\n", data.total_flow);

    // FIXME for test

    work.resize(data.sources.size());
}

void
SDM::apply_flow(double k, bool use_max)
{
    if ((size_t)data.flow.size() == (size_t)data.edges.size())
    {
        for (int i = 0; i < data.flow.size(); ++i)
        {
            auto capacity = data.flow[i] * k;
            if (use_max)
            {
                capacity = std::max(capacity, data.edges[i].capacity);
            }

            data.edges[i].capacity = capacity;
        }

        printf("flow applied with k = %g\n", k);
    }
}

void
SDM::apply_total_flow(double k)
{
    for (auto& e : data.edges)
    {
        e.capacity = data.total_flow * k;
    }

    printf("total flow applied with k = %g\n", k);
}

std::pair<double, double>
SDM::calc_exp_sum(const DVector& x, const tntp::Edge& edge, double mu)
{
    auto u_max = 0.0; // we should fix only big POSITIVE value from exp(value)

    size_t i = 0;
    for (auto s : data.sources)
    {
        work[i] = (T(x, s, edge.target) - T(x, s, edge.source) - edge.free_flow_time) / (edge.free_flow_time * mu);
        u_max = std::max(u_max, work[i]);
        ++i;
    }

    auto exp_sum = 0.0;
    for (size_t i = 0; i < data.sources.size(); ++i)
    {
        work[i] = std::exp(work[i] - u_max);
        exp_sum += work[i];
    }
    exp_sum += std::exp(-u_max);

    return std::make_pair(u_max, exp_sum);
}

}

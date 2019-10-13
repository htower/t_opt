#include "smvsdm2.hpp"

#include "core/blas.hpp"

namespace transport
{

SmVSDM2::SmVSDM2(const String & path, const String & name)
    : SDM("SmVSDM2", path, name)
{
    m_properties |= ProblemProperty::LipschitzConstant;
    set_mu(1.0);

//     m_dual_size = m_size * m_size; // FIXME add sparsity
    m_dual_size = data.edges.size();

    double d_min = std::numeric_limits<double>::max();
    double d_max = std::numeric_limits<double>::lowest();
    double d_med = 0.0;

    for (const auto& d : data.trips)
    {
        d_min = std::min(d_min, d.flow);
        d_max = std::max(d_max, d.flow);
        d_med += d.flow;
    }

    d_med /= data.trips.size();
    printf("d: min = %e; max = %e; med = %e\n", d_min, d_max, d_med);
}

void
SmVSDM2::set_mu(double mu)
{
    this->mu = mu;

    auto max_f = limits<double>::lowest();
    auto min_tf = limits<double>::max();

    for (const auto& e : data.edges)
    {
        max_f = std::max(max_f, e.capacity);
        min_tf = std::min(min_tf, e.capacity * e.free_flow_time);
    }
    max_f *= max_f;

    m_l = max_f / (mu * min_tf);

    printf("mu = %e L = %e\n", mu, m_l);
}

void
SmVSDM2::restore_flow(Point& point)
{
//     const auto& x = point.x;
    auto& x = point.x;

    // FIXME Yura's idea - really not needed ?
//     for (size_t i = 1; i <= data.max_node_index; ++i)
//     {
//         for (size_t j = 1; j <= data.max_node_index; ++j)
//         {
//             if (j != i)
//             {
//                 T(x, i, j) -= T(x, i, i);
//             }
//         }
//         T(x, i, i) = 0.0;
//     }
//     printf("T_ii\n\n");

    DVector flow(data.edges.size());
    for (size_t i = 0; i < data.edges.size(); ++i)
    {
        auto pair = calc_exp_sum(x, data.edges[i], mu);
        auto u_max = pair.first;
        auto exp_sum = pair.second;

        flow[i] = data.edges[i].capacity * (1.0 - std::exp(-u_max) / exp_sum);
    }

    auto t_max = limits<double>::lowest();
    auto t_min = limits<double>::max();
    for (size_t i = 0; i < (size_t)x.size(); ++i)
    {
        t_max = std::max(t_max, x[i]);
        t_min = std::min(t_min, x[i]);
    }

    printf("t_max = %e t_min = %e\n", t_max, t_min);

//     if self.data.flow.len() == self.data.edges.len() {
//             let mut d_nrm1 = 0.0;
//             let mut d_nrm2 = 0.0;
//
//             for (_e, (f_o, f_c)) in self.data.edges.iter().zip(self.data.flow.iter().zip(flow.iter())) {
// //                println!("{:6} -> {:6} : {:.10e} {:.10e} {:.10e}", _e.source, _e.target, f_o, f_c, f_o - f_c);
//                 let d = f_o - f_c;
//                 d_nrm1 += d.abs();
//                 d_nrm2 += d * d;
//             }
//
//             println!("d_nrm1 = {:12.4e}", d_nrm1);
//             println!("d_nrm2 = {:12.4e}", d_nrm2);
//
//             return d_nrm2;
//         } else {
//             for (e, f_c) in self.data.edges.iter().zip(flow.iter()) {
//                 println!("{:6} -> {:6} : {:.6} {:.6}", e.source, e.target, f_c, e.capacity);
//             }
//         }

//     for (auto a : boost::combine(edges, trips))
//     {
//         auto& z = a.get<0>();
//         z.b = 0.0;
//     }

    for (size_t i = 0; i < data.edges.size(); ++i)
    {
        if (std::abs(flow[i]) > 1e-3)
        {
            const auto& edge = data.edges[i];
            printf("%d -> %d : %8.2f / %8.2f (%7.3f %%) {% e}\n",
//             printf("%d -> %d : %e / %e (%6.3f %%)\n",
                   edge.source, edge.target,
                   flow[i], edge.capacity,
                   flow[i] / edge.capacity * 100.0,
                   edge.capacity - flow[i]);
        }
    }
}

void
SmVSDM2::f(Point& p)
{
    // FIXME move into class
//     static const double l_log = std::log(1.0 / (1.0 + data.sources.size()));
//     static const double l_log1 = std::log(1.0 + data.sources.size());

    const auto& x = p.x;

    p.f = 0.0;
    if (mu > 0.0)
    {
        for (const auto& edge : data.edges)
        {
            auto pair = calc_exp_sum(x, edge, mu);
            auto u_max = pair.first;
            auto exp_sum = pair.second;

            p.f += edge.free_flow_time * edge.capacity * (u_max + std::log(exp_sum));
//             p.f += edge.free_flow_time * edge.capacity * (u_max + std::log(exp_sum) + l_log);
//             p.f += edge.free_flow_time * edge.capacity * (u_max + std::log(exp_sum) - l_log1);
        }
        p.f *= mu;
    }

    for (const auto& trip : data.trips)
    {
        const auto s = trip.source;
        const auto k = trip.target;

        p.f -= trip.flow * (T(x, s, k) - T(x, s, s));
    }
}

void
SmVSDM2::df(Point& p)
{
    // FIXME move into class
//     static const double l_log = std::log(1.0 / (1.0 + data.sources.size()));

    const auto& x = p.x;
    auto& g = p.g;

    for (const auto& edge : data.edges)
    {
        const auto i = edge.source;
        const auto j = edge.target;
        const auto f = edge.capacity;

        const auto pair = calc_exp_sum(x, edge, mu);
        const auto exp_sum = pair.second;
        const auto exp_sum_inv = 1.0 / exp_sum;

        size_t k = 0;
        for (auto s : data.sources)
        {
            const auto g_sji = work[k] * f * exp_sum_inv;
            T(g, s, j) += g_sji;
            T(g, s, i) -= g_sji;

            ++k;
        }
    }

    for (const auto& trip : data.trips)
    {
        const auto s = trip.source;
        const auto k = trip.target;

        T(g, s, k) -= trip.flow;
        T(g, s, s) += trip.flow;
    }
}

void
SmVSDM2::dual_x(Point& p, Point& dual_p)
{
    auto& flow = dual_p.x;
    for (size_t i = 0; i < data.edges.size(); ++i)
    {
        auto pair = calc_exp_sum(p.x, data.edges[i], mu);
        auto u_max = pair.first;
        auto exp_sum = pair.second;

        flow[i] = data.edges[i].capacity * (1.0 - std::exp(-u_max) / exp_sum);
    }
}

void
SmVSDM2::dual_f(Point& dual_p)
{
    auto z = [this](double flow, double free_flow_time, double capacity) -> double
    {
        double res = 0.0;

        res += mu * (flow - capacity);
        if (std::fabs(res) < 1e-8)
        {
            res = 0.0;
        }
        else
        {
            res *= log(capacity / (capacity - flow));
        }

        res += (mu + 1.0) * flow;
        res *= free_flow_time;

        return res;
    };

//     auto z = [this](double flow, double free_flow_time, double capacity) -> double
//     {
//         if (flow <= capacity)
//         {
//             return free_flow_time * flow;
//         }
//
//         return 1e300 * flow;
//     };

    // FIXME move into class
    static const double l_log = std::log(1.0 / (1.0 + data.sources.size()));

    dual_p.f = 0.0;
    for (size_t i = 0; i < data.edges.size(); ++i)
    {
        dual_p.f += z(dual_p.x[i], data.edges[i].free_flow_time, data.edges[i].capacity);
        dual_p.f -= z(0.0,         data.edges[i].free_flow_time, data.edges[i].capacity);

//         dual_p.f += std::min(dual_p.x[i], data.edges[i].capacity) * data.edges[i].free_flow_time;
    }
}

}

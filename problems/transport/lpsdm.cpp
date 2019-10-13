#include "lpsdm.hpp"

namespace transport
{

LPSDM::LPSDM(const String& path, const String& name)
    : SDM("LPSDM", path, name)
{
    T_size = m_size;
    m_size += data.edges.size();
}

void
LPSDM::f(Point& p)
{
    const auto& x = p.x;
    uint32_t i;

    p.f = 0.0;

    double phi1 = 0.0;
    for (const auto& trip : data.trips)
    {
        const auto i = trip.source;
        const auto j = trip.target;
        const auto d = trip.flow;

        phi1 -= d * T(x, i, j);
    }
    p.f += phi1;

    double phi2 = 0.0;
    i = 0;
    for (const auto& edge : data.edges)
    {
        const auto f = edge.capacity;
        const auto t = edge.free_flow_time;

        phi2 += f * (x[T_size + i] - t);
        ++i;
    }
    p.f += phi2;

    double phi3 = 0.0;
    for (auto s : data.sources)
    {
        // T_ss = 0
        const auto z = T(x, s, s);
        phi3 += K * z * z;
    }
    p.f += phi3;

    double phi4 = 0.0;
    i = 0;
    for (const auto& edge : data.edges)
    {
        const auto t = edge.free_flow_time;
        const auto z = t - x[T_size + i];

        // t_ >= t
        if (z > 0)
        {
            phi4 += K * z * z;
        }
        ++i;
    }
    p.f += phi4;

    double phi5 = 0.0;
    uint32_t k = 0;
    for (const auto& edge : data.edges)
    {
        const auto i = edge.source;
        const auto j = edge.target;

        for (auto s : data.sources)
        {
            const auto z = T(x, s, j) - T(x, s, i) - x[T_size + k];

            if (z > 0)
            {
                phi5 += K * z * z;
            }
        }
        ++k;
    }
    p.f += phi5;
}

void
LPSDM::df(Point& p)
{
    const auto& x = p.x;
    auto& g = p.g;
    uint32_t i;

    for (const auto& trip : data.trips)
    {
        const auto i = trip.source;
        const auto j = trip.target;
        const auto d = trip.flow;

        T(g, i, j) -= d;
    }

    i = 0;
    for (const auto& edge : data.edges)
    {
        const auto f = edge.capacity;

        g[T_size + i] += f;
        ++i;
    }

    for (auto s : data.sources)
    {
        T(g, s, s) += 2.0 * K * T(x, s, s);
    }

    i = 0;
    for (const auto& edge : data.edges)
    {
        const auto t = edge.free_flow_time;
        const auto z = t - x[T_size + i];

        if (z > 0)
        {
            g[T_size + i] -= 2.0 * K * z;
        }
        ++i;
    }

    uint32_t k = 0;
    for (const auto& edge : data.edges)
    {
        const auto i = edge.source;
        const auto j = edge.target;

        for (auto s : data.sources)
        {
            const auto z = T(x, s, j) - T(x, s, i) - x[T_size + k];

            if (z > 0)
            {
                T(g, s, j) += 2.0 * K * z;
                T(g, s, i) -= 2.0 * K * z;
                g[T_size + k] -= 2.0 * K * z;
            }
        }
        ++k;
    }
}

void
LPSDM::emoe(Point& point)
{
    const auto& x = point.x;
    uint32_t i;

    double phi1 = 0.0;
    for (const auto& trip : data.trips)
    {
        const auto i = trip.source;
        const auto j = trip.target;
        const auto d = trip.flow;

        phi1 -= d * T(x, i, j);
    }
    printf("phi1 = %e\n", phi1);

    double phi2 = 0.0;
    i = 0;
    for (const auto& edge : data.edges)
    {
        const auto f = edge.capacity;
        const auto t = edge.free_flow_time;

        phi2 += f * (x[T_size + i] - t);
        ++i;
    }
    printf("phi2 = %e\n", phi2);

    double phi3 = 0.0;
    for (auto s : data.sources)
    {
        // T_ss = 0
        const auto z = T(x, s, s);
        phi3 += K * z * z;
    }
    printf("phi3 = %e\n", phi3 / K);

    double phi4 = 0.0;
    i = 0;
    for (const auto& edge : data.edges)
    {
        const auto t = edge.free_flow_time;
        const auto z = t - x[T_size + i];

        // t_ >= t
        if (z > 0)
        {
            phi4 += K * z * z;
        }
        ++i;
    }
    printf("phi4 = %e\n", phi4 / K);

    double phi5 = 0.0;
    uint32_t k = 0;
    for (const auto& edge : data.edges)
    {
        const auto i = edge.source;
        const auto j = edge.target;

        for (auto s : data.sources)
        {
            const auto z = T(x, s, j) - T(x, s, i) - x[T_size + k];

            if (z > 0)
            {
                phi5 += K * z * z;
            }
        }
        ++k;
    }
    printf("phi5 = %e\n", phi5 / K);
}

void
LPSDM::restore_flow(Point& point)
{
    const auto& x = point.x;

    for (const auto& edge : data.edges)
    {
        const auto i = edge.source;
        const auto j = edge.target;

        const auto f = edge.capacity;
//         const auto t = edge.free_flow_time;

        auto pair = calc_exp_sum(x, edge, 1e-8);
        auto u_max = pair.first;
        auto exp_sum = pair.second;

        auto flow = f * (1.0 - std::exp(-u_max) / exp_sum);

//         double max = 0.0;
//         for (auto s : sources)
//         {
//             const auto value = T(x, s, j) - T(x, s, i) - t;
//             max = std::max(max, value);
//         }
//         auto flow = f * (1.0 - 1.0 / max);

        if (std::abs(flow) > 1e-3)
        {
            printf("%d -> %d : %e / %e (%6.3f %%)\n",
                    i, j,
                    flow, f,
                    flow / f * 100.0);
        }
    }


}

}

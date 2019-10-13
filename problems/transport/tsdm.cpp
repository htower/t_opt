#include "tsdm.hpp"

#include "core/blas.hpp"

namespace transport
{

#define BETTER_GRADIENT

TSDM::TSDM(const String& path, const String& name)
    : SDM("TSDM", path, name)
{
}

void
TSDM::restore_flow(Point& point)
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

void
TSDM::f(Point& p)
{
    const auto& x = p.x;

    double phi1 = 0.0;
    for (const auto& trip : data.trips)
    {
        const auto s = trip.source;
        const auto k = trip.target;
        const auto d = trip.flow;

        phi1 -= d * (T(x, s, k) - T(x, s, s));
    }

    double phi2 = 0.0;
    for (const auto& edge : data.edges)
    {
        const auto i = edge.source;
        const auto j = edge.target;

        const auto f = edge.capacity;
        const auto t = edge.free_flow_time;

        double max = 0.0;
        for (auto s : data.sources)
        {
            const auto value = T(x, s, j) - T(x, s, i) - t;
            max = std::max(max, value);
        }

        phi2 += f * max;
    }

    p.f = phi1 + phi2;
}

#ifdef BETTER_GRADIENT

// this version is better by result function value (tested on SiouxFalls)
// but slightly slower than upper ("plain") - about 12% slowdown

void
TSDM::df(Point& p)
{
    const auto& x = p.x;
    auto& g = p.g;

    for (const auto& trip : data.trips)
    {
        const auto s = trip.source;
        const auto k = trip.target;
        const auto d = trip.flow;

        T(g, s, k) -= d;
        T(g, s, s) += d;
    }

    for (const auto& edge : data.edges)
    {
        const auto i = edge.source;
        const auto j = edge.target;

        const auto f = edge.capacity;
        const auto t = edge.free_flow_time;

        double max = 0.0;
        bool has_max = false;
        for (uint32_t k = 0; k < data.sources.size(); ++k)
        {
            const auto s = data.sources[k];
            work[k] = T(x, s, j) - T(x, s, i) - t;
            if (work[k] > max)
            {
                max = work[k];
                has_max = true;
            }
        }

        if (has_max == false)
        {
            continue;
        }

        for (uint32_t k = 0; k < data.sources.size(); ++k)
        {
            if (work[k] == max)
            {
                const auto s = data.sources[k];
                T(g, s, j) += f;
                T(g, s, i) -= f;
            }
        }
    }
}

#else

void
TSDM::df(Point& p)
{
    const auto& x = p.x;
    auto& g = p.g;

    for (const auto& trip : data.trips)
    {
        const auto s = trip.source;
        const auto k = trip.target;
        const auto d = trip.flow;

        T(g, s, k) -= d;
        T(g, s, s) += d;
    }

    for (const auto& edge : data.edges)
    {
        const auto i = edge.source;
        const auto j = edge.target;

        const auto f = edge.capacity;
        const auto t = edge.free_flow_time;

        double max = 0.0;
        int64_t s_max = -1; // since source index is uint32_t
        for (auto s : data.sources)
        {
            const auto value = T(x, s, j) - T(x, s, i) - t;
            if (value > max)
            {
                max = value;
                s_max = s;
            }
        }

        if (s_max > 0)
        {
            T(g, s_max, j) += f;
            T(g, s_max, i) -= f;
        }
    }
}

#endif

}


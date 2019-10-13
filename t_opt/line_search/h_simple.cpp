#include "h_simple.hpp"
#include "core/blas.hpp"
#include "core/problem.hpp"

namespace t_opt
{

namespace line_search
{

HSimple::HSimple()
    : HSimple(0.5, 1.5)
{
}

HSimple::HSimple(double step_minus_k, double step_plus_k)
    : step_minus_k(step_minus_k)
    , step_plus_k(step_plus_k)
    , step_min(limits<double>::epsilon())
    , fixed_step(limits<double>::quiet_NaN())
{
}

HSimple::HSimple(double fixed_step)
    : step_minus_k(limits<double>::quiet_NaN())
    , step_plus_k(limits<double>::quiet_NaN())
    , step_min(limits<double>::epsilon())
    , fixed_step(fixed_step)
{
}

void
HSimple::setup(Problem& problem)
{
    probe_point.resize(problem);
}

LineSearchProbe
HSimple::probe(Problem& problem, const Point& point, const DVector& dir, double step)
{
    blas::axpyz(step, dir, point.x, probe_point.x);
    problem.f(probe_point);
//     printf("step = %g f = %g\n", step, probe_point.f);
    return { step, probe_point.f };
}

LineSearchProbe
HSimple::search(Problem& problem, const Point& point, const DVector& dir, bool dir_is_gradient, double start_step)
{
    if (std::isnan(fixed_step) == false)
    {
        fixed_step = fix_step(fixed_step, dir_is_gradient);
        return probe(problem, point, dir, fixed_step);
    }

//     printf("\n");
    start_step = fix_step(start_step, dir_is_gradient);
    auto probe_0 = probe(problem, point, dir, start_step);

    if (probe_0.f < point.f)
    {
        if (step_plus_k != 0.0)
        {
            auto probe_1 = probe(problem, point, dir, start_step * step_plus_k);
            if (probe_1.f < probe_0.f)
            {
                return probe_1;
            }
        }

        return probe_0;
    }

    while (probe_0.f >= point.f)
    {
        if (std::abs(probe_0.step) < step_min)
        {
            return { 0.0, point.f };
        };

        probe_0 = probe(problem, point, dir, probe_0.step * step_minus_k);
    }

    return probe_0;
}

}

}

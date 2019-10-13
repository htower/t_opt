#include "parabolic.hpp"

#include "core/blas.hpp"
#include "core/problem.hpp"

namespace t_opt
{

namespace line_search
{

// FIXME made this line searcher "plain parabolic"
// 1) pass single bool use_gradient parameter and set probes to 2 when true and 3 when false
// 2) add new h_parabolic searcher which combines 2 approaches

Parabolic::Parabolic(uint8_t probes, bool use_gradient)
    : probes(probes)
    , max_probes(20) // FIXME
    , use_gradient(use_gradient)
{
}

void
Parabolic::setup(Problem& problem)
{
    probe_point.resize(problem);

//         if self.use_gradient {
//             if !problem.has_gradient() {
//                 assert!(
//                     false,
//                     "Parabolic line search method with enabled gradient, can't be applied to the problem {}",
//                     problem.name());
//             }
//
//             if self.probes < 2 {
//                 assert!(
//                     false,
//                     "Parabolic line search method with enabled gradient requires at least 2 probes per search");
//             }
//         } else {
//             if self.probes < 3 {
//                 assert!(
//                     false,
//                     "Parabolic line search method with disabled gradient requires at least 3 probes per search");
//             }
//         }
}

void
Parabolic::probe(Problem& problem, const Point& point, const DVector& dir, uint8_t i)
{
    blas::axpyz(p[i].step, dir, point.x, probe_point.x);
    problem.f(probe_point);
    p[i].f = probe_point.f;
//     printf("step = %g f = %g\n", p[i].step, p[i].f);
}

void
Parabolic::sort_probes(uint8_t count)
{
    // Simple bubble sort

    for (uint8_t j = 0; j < (count - 1); ++j)
    {
        for (uint8_t i = 0; i < (count - j - 1); ++i)
        {
            if (p[i].f > p[i + 1].f)
            {
                std::swap(p[i], p[i + 1]);
            }
        }
    }
}

// FIXME place result inside vector instead return?
double
solve_3(const LineSearchProbe* p)
{
    const double det_a =
        p[0].f * (p[1].step - p[2].step) -
        p[0].step * (p[1].f - p[2].f) +
        (p[1].f * p[2].step - p[2].f * p[1].step);

    const double det_b =
        p[0].step * p[0].step * (p[1].f - p[2].f) -
        p[0].f * (p[1].step * p[1].step - p[2].step * p[2].step) +
        (p[1].step * p[1].step * p[2].f - p[2].step * p[2].step * p[1].f);

    return -0.5 * det_b / det_a;
}

LineSearchProbe
Parabolic::search(Problem& problem, const Point& point, const DVector& dir, bool dir_is_gradient, double start_step)
{
//     printf("\n");
    // FIXME check start_step with isfinite()
    start_step = fix_step(start_step, dir_is_gradient);

    p[0].step = 0.0;
    p[0].f = point.f;

    bool was_sorted = false;
    for (uint8_t i = 1; i <= probes; ++i)
    {
        if (i == 1)
        {
            if (use_gradient)
            {
                p[1].step = start_step;
            }
            else
            {
                p[1].step = -0.5 * start_step;
            }

            probe(problem, point, dir, 1);
        }
        else if (i == 2)
        {
            if (use_gradient)
            {
                // projection of gradient onto direction
                auto gdp = dir_is_gradient ? point.g_nrm2_2 : blas::dot(point.g, dir);

                p[2].step = -0.5 *
                            (gdp * p[1].step * p[1].step) /
                            (p[1].f - p[0].f - gdp * p[1].step);

                if (std::isfinite(p[2].step))
                {
                    probe(problem, point, dir, 2);
                }
                else
                {
                    // FIXME try use h_simple strategy
                    p[2].step = 0.0;
                    p[2].f = point.f;
                }
            }
            else
            {
                p[2].step = -p[1].step;
                probe(problem, point, dir, 2);
            };
        }
        else
        {
            p[3].step = solve_3(p);

            // FIXME move this check into probe() ?
            if (std::isfinite(p[2].step))
            {
                probe(problem, point, dir, 3);
            }
            else
            {
                p[3].step = 0.0;
                p[3].f = point.f;
            }

            sort_probes(4);
            was_sorted = true;
        }
    }

    if (was_sorted == false)
    {
        sort_probes(std::min(probes + 1, 4));
    }

    if (p[0].f < point.f)
    {
        return p[0];
    }

    for (uint8_t i = probes; i < max_probes; ++i)
    {
        p[3].step = solve_3(p);
        probe(problem, point, dir, 3);
        sort_probes(4);

        if (p[0].f < point.f)
        {
            return p[0];
        }
    }

    return {0.0, point.f};
}

}

}

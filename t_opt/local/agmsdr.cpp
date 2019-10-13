#include "agmsdr.hpp"

#include "core/blas.hpp"
#include "core/line_search.hpp"
#include "core/problem.hpp"

namespace t_opt
{

namespace local
{

AGMsDR::AGMsDR(LineSearchMethod& ls, double epsilon)
    : Method("AGMsDR", ProblemProperty::Gradient)
    , ls(ls)
    , ls_start_step(1.0)
    , ls_step(1.0)
    , epsilon(epsilon)
    , A(0.0)
{
}

void
AGMsDR::before(Problem& problem, Point& point)
{
    ls.setup(problem);

    // FIXME test this
//     ls_start_step = -1.0 / point.g_nrm2;
//     ls_start_step = -1.0 / point.g_nrm2_2;
    ls_start_step = 1.0;

    ls_step = ls_start_step;

    A = 0.0;
    y = Point(point); // FIXME add operator=() ?

    v.resize(problem.size());
    blas::copy(point.x, v);

    vx.resize(problem.size());
}

bool
AGMsDR::iteration(Problem& problem, Point& point, size_t iter)
{
    // At first iteration (== 0) y.x will be equal to point.x, so we set it to
    // start point inside the before() method.

    if (iter > 0)
    {
        // Original paper starts search from v to x_k (== point.x here). But this will
        // require calculation of f()/df() at v point. So we reverse search start point
        // and direction to avoid extra calculations.

        // vx = v - point.x
        blas::xmyz(v, point.x, vx);

        // here we need some other line searcher ?
        const auto probe = ls.search(problem, point, vx, false, 1.0);

        // beta(step) should be in [0, 1] interval
        // TODO use std::clamp() when C++17 will be used
        auto beta = std::min(probe.step, 1.0);
        beta = std::max(beta, 0.0);

        // y.x = point.x + beta * vx
        blas::axpyz(beta, vx, point.x, y.x);

        if (beta == probe.step)
        {
            y.f = probe.f;
        }
        else
        {
            problem.f(y);
        }
        problem.df(y);
    }

    const auto probe = ls.search(problem, y, y.g, true, ls_step);
    if (probe.step == 0.0)
    {
        printf("\nLS = 0; ERROR! %e\n", ls_step);
        return false;
    }

    ls_step = probe.step;

    // p.x = y.x + ls_step * y.g
    blas::axpyz(ls_step, y.g, y.x, point.x);

    point.f = probe.f;
    problem.df(point);

    // Here we solve internal quadratic equation
    double a;
    {
        const auto eq_a = -y.g_nrm2_2;
        const auto eq_b = 2.0 * (y.f - point.f);
        const auto eq_c = A * eq_b;

        const auto eq_d = eq_b * eq_b - 4.0 * eq_a * eq_c;
        if (eq_d < 0.0)
        {
            printf("\nD = %g; ERROR!", eq_d);
            return false;
        }

        const auto eq_x_1 = 0.5 * (-eq_b + std::sqrt(eq_d)) / eq_a;
        const auto eq_x_2 = 0.5 * (-eq_b - std::sqrt(eq_d)) / eq_a;

        a = std::max(eq_x_1, eq_x_2);
        A += a;
    }

    // v -= a * y.g
    blas::axpy(-a, y.g, v);

    return true;
}

}

}

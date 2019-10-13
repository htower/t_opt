#include "lbfgs.hpp"

#include "core/blas.hpp"
#include "core/line_search.hpp"
#include "core/logger.hpp"
#include "core/problem.hpp"

#include <fmt/core.h>

namespace t_opt
{

namespace local
{

LBFGS::LBFGS(uint32_t m, LineSearchMethod& ls)
    : Method(fmt::format("LBFGS_{}", m), ProblemProperty::Gradient)
    , m(m)
    , ls(ls)
{
    l_s.resize(m);
    l_y.resize(m);

    l_ys.resize(m);
    l_alpha.resize(m);
}

void
LBFGS::before(Problem & problem, Point & point)
{
    ls.setup(problem);
    ls_step = ls_start_step;

    for (auto& s : l_s)
    {
        s.resize(problem.size());
    }

    for (auto& y : l_y)
    {
        y.resize(problem.size());
    }


    x_p.resize(problem.size());
    g_p.resize(problem.size());
    dir.resize(problem.size());

    l_end = 0;

    // Reset direction to (normalized) antigradient
    blas::scal_copy(-1.0 / point.g_nrm2, point.g, dir);
}

void
LBFGS::log(Logger& logger, LoggerMode mode) const
{
    logger.put_ls_step(mode, ls_step);
}

bool
LBFGS::iteration(Problem& problem, Point& point, size_t iter)
{
    auto bound = std::min(m, (uint32_t)iter + 1); // +1 since iter == 0 at start

    blas::copy(point.x, x_p);
    blas::copy(point.g, g_p);

    auto probe = ls.search(problem, point, dir, false, ls_step);
    if (probe.step == 0.0)
    {
        if (iter == 0)
        {
            // Direction is already antigradient
            return false;
        }

        // Reset direction to (normalized) antigradient
        blas::scal_copy(-1.0 / point.g_nrm2, point.g, dir);

        probe = ls.search(problem, point, dir, false, ls_start_step);
        if (probe.step == 0.0)
        {
            return false;
        }

        bound = 1;
    }

    ls_step = probe.step;

    blas::axpy(ls_step, dir, point.x);
    point.f = probe.f;
    problem.df(point);

    blas::xmyz(point.x, x_p, l_s[l_end]);
    blas::xmyz(point.g, g_p, l_y[l_end]);

    auto ys_end = blas::dot(l_y[l_end], l_s[l_end]);
    auto yy_end = blas::dot(l_y[l_end], l_y[l_end]);

    l_ys[l_end] = ys_end;

    blas::scal_copy(-1.0, point.g, dir);

    l_end = (l_end + 1) % m;
    auto j = l_end;

    for (uint32_t i = 0; i < bound; ++i)
    {
        j = (j + m - 1) % m;

        l_alpha[j] = blas::dot(l_s[j], dir) / l_ys[j];

        blas::axpy(-l_alpha[j], l_y[j], dir);
    }

    blas::scal(ys_end / yy_end, dir);

    for (uint32_t i = 0; i < bound; ++i)
    {
        auto l_beta = blas::dot(l_y[j], dir) / l_ys[j];

        blas::axpy(l_alpha[j] - l_beta, l_s[j], dir);

        j = (j + 1) % m;
    }

    return true;
}

}

}

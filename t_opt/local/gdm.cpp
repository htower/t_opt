#include "gdm.hpp"

#include "core/blas.hpp"
#include "core/logger.hpp"
#include "core/problem.hpp"
#include "core/line_search.hpp"

namespace t_opt
{

namespace local
{

GDM::GDM(LineSearchMethod& ls)
    : Method("GDM", ProblemProperty::Gradient)
    , ls(ls)
{
}

void
GDM::before(Problem & problem, Point & point)
{
    ls.setup(problem);

    // FIXME test this
    ls_step = ls_start_step = -1.0 / point.g_nrm2;
//     ls_step = ls_start_step = -1.0 / point.g_nrm2_2; // too small start step ?
}

bool
GDM::iteration(Problem& problem, Point& point, size_t /*iter*/)
{
    auto probe = ls.search(problem, point, point.g, true, ls_step);
    if (probe.step == 0.0)
    {
        return false;
    }

    ls_step = probe.step;
    blas::axpy(ls_step, point.g, point.x);
    point.f = probe.f;

    problem.df(point);

    return true;
}

void
GDM::log(Logger& logger, LoggerMode mode) const
{
    logger.put_ls_step(mode, -ls_step);
}

}

}

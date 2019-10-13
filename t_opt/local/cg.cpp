#include "cg.hpp"

#include "core/blas.hpp"
#include "core/line_search.hpp"
#include "core/logger.hpp"
#include "core/problem.hpp"

#include <fmt/format.h>

namespace t_opt
{

namespace local
{

// TODO for CG:
// * code cleanup
// * use single-calculated bools instead calling cg_use_* everytime - seems no speedup
// * normal code for restart
// * when restart do only necessary actions (check before_step() - x/g/copy and so on)
// * unify restart tecnhiques (dir,step reset, etc) for CG/LBFGS/etc

// FIXME PRP(+)/CD/LS is broken for Quadratic(3e4) ??? - check with parabolic line search
// PRP with non-normalized dir is better than normalized one - WTF???
// when reset dir to -g ls_step should be resetted to ls_start_step !!!

String
to_string(CgVariant variant)
{
    switch (variant)
    {
        case CgVariant::HS:
        {
            static const String s("HS");
            return s;
        }

        case CgVariant::FR:
        {
            static const String s("FR");
            return s;
        }

        case CgVariant::PRP:
        {
            static const String s("PRP");
            return s;
        }

        case CgVariant::PRPplus:
        {
            static const String s("PRP+");
            return s;
        }

        case CgVariant::CD:
        {
            static const String s("CD");
            return s;
        }

        case CgVariant::LS:
        {
            static const String s("LS");
            return s;
        }

        case CgVariant::DY:
        {
            static const String s("DY");
            return s;
        }
    }

    // -Wreturn-type warning fix
    return {};
}

bool
cg_use_y(CgVariant variant)
{
    switch (variant)
    {
        case CgVariant::HS:
        case CgVariant::PRP:
        case CgVariant::PRPplus:
        case CgVariant::LS:
        case CgVariant::DY:
            return true;

        case CgVariant::FR:
        case CgVariant::CD:
            return false;
    }

    // FIXME Wreturn-type warning fix
    return false;
}

bool
cg_use_s(CgVariant variant)
{
    switch (variant)
    {
        case CgVariant::HS:
        case CgVariant::DY:
            return true;

        case CgVariant::FR:
        case CgVariant::PRP:
        case CgVariant::PRPplus:
        case CgVariant::CD:
        case CgVariant::LS:
            return false;
    }

    // FIXME Wreturn-type warning fix
    return false;
}

bool
cg_use_g_prev(CgVariant variant)
{
    switch (variant)
    {
        case CgVariant::CD:
        case CgVariant::LS:
            return true;

        case CgVariant::HS:
        case CgVariant::FR:
        case CgVariant::PRP:
        case CgVariant::PRPplus:
        case CgVariant::DY:
            return false;
    }

    // FIXME Wreturn-type warning fix
    return false;
}

bool
cg_dir_by_s(CgVariant variant)
{
    switch (variant)
    {
        case CgVariant::HS:
        case CgVariant::DY:
            return true;

        case CgVariant::FR:
        case CgVariant::PRP:
        case CgVariant::PRPplus:
        case CgVariant::CD:
        case CgVariant::LS:
            return false;
    }

    // FIXME Wreturn-type warning fix
    return false;
}

CG::CG(CgVariant variant, LineSearchMethod& ls)
    : Method(fmt::format("CG_{}", to_string(variant)), ProblemProperty::Gradient)
    , variant(variant)
    , ls(ls)
    , ls_start_step(1.0) // FIXME
{
//     use_s = cg_use_s(variant);
//     use_y = cg_use_y(variant);
//     use_g_prev = cg_use_g_prev(variant);
//     dir_by_s = cg_dir_by_s(variant);
}

void
CG::before(Problem& problem, Point& point)
{
//     ls_start_step = 1.0 / point.g_nrm2;

    ls.setup(problem);
    ls_step = ls_start_step;

    // FIXME remove ?
    g_nrm2 = g_nrm2_prev = limits<double>::quiet_NaN();

    dir.resize(problem.size());
    dir_prev.resize(problem.size());
    dir_normalized.resize(problem.size());

//     if (use_s)
    if (cg_use_s(variant))
    {
        s.resize(problem.size());
    }

//     if (use_y)
    if (cg_use_y(variant))
    {
        y.resize(problem.size());
    }

//     if (use_g_prev)
    if (cg_use_g_prev(variant))
    {
        g_prev.resize(problem.size());
    }
}

void
CG::log(Logger& logger, LoggerMode mode) const
{
    logger.put_ls_step(mode, ls_step);
}

void
CG::calc_dir(const DVector& g)
{
    // dir = -g
    blas::scal_copy(-1.0, g, dir);
    if (dir_reset)
    {
        return;
    }

    double beta = 1.0;
    switch (variant)
    {
        case CgVariant::HS:
            beta = blas::dot(y, g) / blas::dot(y, s);
            break;

        case CgVariant::FR:
            beta = g_nrm2 * g_nrm2 / (g_nrm2_prev * g_nrm2_prev);
            break;

        case CgVariant::PRP:
            beta = blas::dot(y, g) / (g_nrm2_prev * g_nrm2_prev);
            break;

        case CgVariant::PRPplus:
            beta = blas::dot(y, g) / (g_nrm2_prev * g_nrm2_prev);
            beta = std::max(0.0, beta);
            break;

        case CgVariant::CD:
            beta = -g_nrm2 * g_nrm2 / blas::dot(g_prev, dir_prev);
            break;

        case CgVariant::LS:
            beta = -blas::dot(y, g) / blas::dot(g_prev, dir_prev);
            break;

        case CgVariant::DY:
            beta = g_nrm2 * g_nrm2 / blas::dot(y, s);
            break;
    }

//     if (dir_by_s)
    if (cg_dir_by_s(variant))
    {
        // dir += beta * s
        blas::axpy(beta, s, dir);
    }
    else
    {
        // dir += beta * dir_prev
        blas::axpy(beta, dir_prev, dir);
    }
}

void
CG::before_step(const Point& point)
{
    // dir_prev = dir
    blas::copy(dir, dir_prev);
    g_nrm2_prev = point.g_nrm2;

//     if (use_s)
    if (cg_use_s(variant))
    {
        // s = point.x
        blas::copy(point.x, s);
    }

//     if (use_y)
    if (cg_use_y(variant))
    {
        // y = point.g
        blas::copy(point.g, y);
    }

//     if (use_g_prev)
    if (cg_use_g_prev(variant))
    {
        // g_prev = point.g
        blas::copy(point.g, g_prev);
    }
}

void
CG::after_step(const Point& point)
{
    g_nrm2 = point.g_nrm2;

//     if (use_s)
    if (cg_use_s(variant))
    {
        // s = x - s;
        blas::xmy(point.x, s);
    }

//     if (use_y)
    if (cg_use_y(variant))
    {
        // y = g - y;
        blas::xmy(point.g, y);
    }
}

bool
CG::iteration(Problem& problem, Point& point, size_t iter)
{
    dir_reset = ((iter % 100) == 0); // FIXME add reset parameter

    auto do_step = [this, &problem, &point]()
    {
        calc_dir(point.g);
        before_step(point);

        // FIXME add parameter for normalized/plain mode
        if (true)
        {
            // dir_normalized = dir / blas::nrm2(dir)
            blas::scal_copy(1.0 / blas::nrm2(dir), dir, dir_normalized);
        }
        else
        {
            // dir_normalized = dir
            blas::copy(dir, dir_normalized);
        }

        return ls.search(problem, point, dir_normalized, false, ls_step);
    };

    auto probe = do_step();
    if (probe.step == 0.0)
    {
        if (dir_reset)
        {
            return false;
        }

        // FIXME test and fix
        dir_reset = true;
        ls_step = ls_start_step;
        probe = do_step();
        if (probe.step == 0.0)
        {
            return false;
        }
    }

    ls_step = probe.step;
    point.f = probe.f;

    blas::axpy(ls_step, dir_normalized, point.x);
    problem.df(point);

    after_step(point);

    return true;
}

}

}

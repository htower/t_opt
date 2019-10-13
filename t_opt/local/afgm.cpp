#include "afgm.hpp"

#include "core/blas.hpp"
#include "core/line_search.hpp"
#include "core/logger.hpp"
#include "core/problem.hpp"

namespace t_opt
{

namespace local
{

AFGM::AFGM(LineSearchMethod& ls)
    : Method("AFGM", ProblemProperty::Gradient)
    , ls(ls)
    , ls_start_step(1.0)
    , ls_step(1.0)
{
}

void
AFGM::before(Problem& problem, Point& point)
{
    ls.setup(problem);

    // FIXME test this
    ls_start_step = -1.0 / point.g_nrm2;
//     ls_start_step = -1.0 / point.g_nrm2_2;
//     ls_start_step = 1.0;

    ls_step = ls_start_step;

    alpha = alpha_sum = 0.0;

    // FIXME add operator=() ?
    x = Point(point);
    z = Point(point);

//     x.x.setZero();
//     z.x.setZero();

    zy.resize(problem.size());

    // FIXME add method to problem class? or add dual size method ?
    dual_p.x.resize(problem.dual_size());
    dual_p_sum.x.resize(problem.dual_size());
    dual_p_sum.x.setZero();

    pd_delta = 0.0;
}

bool
AFGM::iteration(Problem& problem, Point& y, size_t iter)
{
    if (iter > 0)
    {
        // zy = z - y
        blas::xmyz(z.x, y.x, zy);

        // here we need some other line searcher ?
        const auto probe = ls.search(problem, y, zy, false, 1.0);
//         const auto probe = ls.search(problem, y, zy, false, ls_step);

        // beta(step) should be in [0, 1] interval
        // TODO use std::clamp() when C++17 will be used
        auto tau = std::min(probe.step, 1.0);
        tau = std::max(tau, 0.0);

        // x^{k+1} = \tau z^k + (1 - \tau) * y^k = y^k + \tau (z^k - y^k)
//         blas::axpbyz(tau, z.x, 1.0 - tau, y.x, x.x);
        blas::axpbyz(1.0, y.x, tau, zy, x.x);

        if (tau == probe.step)
        {
            x.f = probe.f;
        }
        else
        {
            // FIXME add check for tau == 0.0
            problem.f(x);

        }
//         if (tau == probe.step)
//         {
//             printf("xf = %e pf = %e dd = %e\n", x.f, probe.f, x.f - probe.f);
//         }
        problem.df(x);
    }

    const auto probe = ls.search(problem, x, x.g, true, ls_step);
    if (probe.step == 0.0)
    {
//         printf("\nLS = 0; ERROR! %e\n", ls_step);
//         return false;
    }
    else
    {
        ls_step = probe.step;
    }

    // y^{k+1} = x^{k+1} - h^{k+1} \nabla f(x^{k+1})
    blas::axpbyz(1.0, x.x, ls_step, x.g, y.x);
    y.f = probe.f;

    problem.df(y); // for logging + line_search

    const double L = x.g_nrm2_2 / (2.0 * (x.f - y.f));
    alpha = 0.5 / L + sqrt(0.25 / (L * L) + alpha * alpha);


    blas::axpy(-alpha, x.g, z.x);

    if (false)
    {
        problem.dual_x(y, dual_p);

        alpha_sum += alpha;
        dual_p_sum.x += alpha * dual_p.x;
        dual_p.x = (1.0 / alpha_sum) * dual_p_sum.x;

        problem.dual_f(dual_p);
        pd_delta = y.f + dual_p.f;
    }

    return true;
}

void
AFGM::log(Logger& logger, LoggerMode mode) const
{
    return;
    // FIXME add method for prime-dual f into logger ?
    char line[1024];

    static const char* F_NAME = "PD_delta";

    static const int F_GN_PRECISION = 16;
    static const int F_GN_WIDTH = F_GN_PRECISION + 7;

    logger.put_delimiter();

//     auto pd_delta = std::fabs(this->pd_delta); // FIXME

    switch (mode)
    {
        case LoggerMode::Header :
            switch (logger.device())
            {
                case LoggerDevice::Stdout :
                    sprintf(line, " %-*s", F_GN_WIDTH - 1, F_NAME);
                    logger.put_value(line);
                    break;

                case LoggerDevice::CsvFile :
                    logger.put_value(F_NAME);
                    break;
            }
            break;

        case LoggerMode::Value :
            switch (logger.device())
            {
                case LoggerDevice::Stdout :
                    sprintf(line, "% -*.*e", F_GN_WIDTH, F_GN_PRECISION, pd_delta);
                    logger.put_value(line);
                    break;

                case LoggerDevice::CsvFile :
                    sprintf(line, "%.*e", 6, pd_delta);
                    logger.put_value(line);
                    break;
            }
            break;
    };
}

}

}

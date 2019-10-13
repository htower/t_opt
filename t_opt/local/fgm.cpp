#include "fgm.hpp"

#include "core/blas.hpp"
#include "core/logger.hpp"
#include "core/problem.hpp"

namespace t_opt
{

namespace local
{

FGM::FGM() : Method("FGM", ProblemProperty::Gradient | ProblemProperty::LipschitzConstant)
{
}

void
FGM::before(Problem& problem, Point& point)
{
    x.resize(problem.size());
    y.resize(problem.size());

    // y = point.x
    blas::copy(point.x, y);

    step = 1.0/problem.L();

    // FIXME add method to problem class? or add dual size method ?
    dual_p.x.resize(problem.dual_size());
    dual_p_sum.x.resize(problem.dual_size());
    dual_p_sum.x.setZero();

    pd_delta = 0.0;
    alpha_sum = 0.0;
}

bool
FGM::iteration(Problem& problem, Point& point, size_t iter)
{
    // x = point.x
    blas::copy(point.x, x);

    // point.x = y - 1/L * point.g
    blas::axpyz(-step, point.g, y, point.x);
    problem.f(point);
    problem.df(point);

    // iter + 1, since iter starts from 0
    const double kk = (iter + 1.0) / (iter + 4.0);

    // y = point.x + k / (k + 3) * (point.x - x) =
    //   = point.x + kk * (point.x - x) =
    //   = point.x + kk * point.x - kk * x =
    //   = (kk + 1) * point.x - kk * x
    blas::axpbyz(kk + 1.0, point.x, -kk, x, y);

    if (false)
    {
        constexpr double alpha = 1.0;

        problem.dual_x(point, dual_p);

        alpha_sum += alpha;
        dual_p_sum.x += alpha * dual_p.x;
        dual_p.x = (1.0 / alpha_sum) * dual_p_sum.x;

        problem.dual_f(dual_p);
        pd_delta = point.f + dual_p.f;
    }

    return true;
}

void
FGM::log(Logger& logger, LoggerMode mode) const
{
    return;
    logger.put_ls_step(mode, step);
    // =============================================================================================


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

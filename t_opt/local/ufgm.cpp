#include "ufgm.hpp"

#include "core/blas.hpp"
#include "core/logger.hpp"
#include "core/problem.hpp"

namespace t_opt
{

namespace local
{

UFGM::UFGM(double epsilon)
    : Method("UFGM", ProblemProperty::Gradient)
    , epsilon(epsilon)
    , alpha_k(limits<double>::quiet_NaN())
    , alpha_kp1(limits<double>::quiet_NaN())
    , l_k(limits<double>::quiet_NaN())
    , l_kp1(limits<double>::quiet_NaN())
{
}

void
UFGM::before(Problem& problem, Point& point)
{
    v_k.resize(problem.size());
    blas::copy(point.x, v_k);

    dyx.resize(problem.size());

    x_kp1.resize(problem);
    y_kp1.resize(problem);
    z_kp1.resize(problem.size());

    alpha_k = alpha_kp1 = 0.0;
    l_k = l_kp1 = 1.0;


    // FIXME add method to problem class? or add dual size method ?
    dual_p.x.resize(problem.dual_size());
    dual_p_sum.x.resize(problem.dual_size());
    dual_p_sum.x.setZero();

    pd_delta = 0.0;
    alpha_sum = 0.0;
}

bool
UFGM::iteration(Problem& problem, Point& point, size_t iter)
{
    l_kp1 = 0.5 * l_k;

    while (true)
    {
        alpha_kp1  = 0.25 / (l_kp1 * l_kp1);
        alpha_kp1 += alpha_k * alpha_k * l_k / l_kp1;
        alpha_kp1  = std::sqrt(alpha_kp1) + 0.5 / l_kp1;

        auto tau_k = 1.0 / (alpha_kp1 * l_kp1);

        // We return y_kp1 as a result, therefore here y_k == point

        // x_kp1.x = tau_k * v_k + (1.0 - tau_k) * y_k.x
        blas::axpbyz(tau_k, v_k, 1.0 - tau_k, point.x, x_kp1.x);

        problem.f(x_kp1);
        problem.df(x_kp1);

        // z_kp1 = v_k - alpha_kp1 * x_kp1.g
        blas::axpyz(-alpha_kp1, x_kp1.g, v_k, z_kp1);

        // y_kp1.x = tau_k * z_kp1 + (1.0 - tau_k) * y_k.x
        blas::axpbyz(tau_k, z_kp1, 1.0 - tau_k, point.x, y_kp1.x);

        problem.f(y_kp1);

        blas::xmyz(y_kp1.x, x_kp1.x, dyx);

        auto cond = x_kp1.f - y_kp1.f;
        cond += blas::dot(x_kp1.g, dyx);
        cond += 0.5 * l_kp1 * blas::dot(dyx, dyx);
        cond += 0.5 * (tau_k * epsilon);

        if (cond >= 0.0)
        {
            break;
        }

        l_kp1 *= 2.0;
    }

    l_k     = l_kp1;
    alpha_k = alpha_kp1;

    // update v_k
    blas::axpy(-alpha_kp1, x_kp1.g, v_k);

    // return y_kp1 as a result
    problem.df(y_kp1);
    point.swap(y_kp1);

    if (false)
    {
        problem.dual_x(point, dual_p);

        const double alpha = alpha_k;
//         const double alpha = iter + 1.0;

        alpha_sum += alpha;
        dual_p_sum.x += alpha * dual_p.x;
        dual_p.x = (1.0 / alpha_sum) * dual_p_sum.x;

        problem.dual_f(dual_p);
        pd_delta = point.f + dual_p.f;
    }

    return true;
}

void
UFGM::log(Logger& logger, LoggerMode mode) const
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

#include "method.hpp"

#include "blas.hpp"
#include "chrono.hpp"
#include "logger.hpp"
#include "problem.hpp"

namespace t_opt
{

struct WrappedProblem : public Problem
{
    WrappedProblem(Problem& problem, State& state)
        : Problem(problem)
        , m_problem(problem)
        , m_state(state)
    {
    }

    inline void
    f(Point& p) override
    {
        m_problem.f(p);
        if (std::isfinite(p.f) == false)
        {
//             printf("\nF IS NOT FINITE\n");
            p.f = limits<double>::max();
        }
//         printf("\n");
//         m_problem.emoe(p); // FIXME

        m_state.f_count += 1;
    }

    inline void
    df(Point& p) override
    {
        blas::set_zero(p.g);
        m_problem.df(p);

//         p.g_nrm2_2 = blas::dot(p.g, p.g);
//         p.g_nrm2 = sqrt(p.g_nrm2_2);

        p.g_nrm_1 = p.g_nrm2_2 = p.g_nrm_inf = 0.0;
        for (size_t i = 0; i < size(); ++i)
        {
            auto value = p.g[i];
            p.g_nrm2_2 += value * value;

            value = std::fabs(value);
            p.g_nrm_1 += value;
            p.g_nrm_inf = std::max(p.g_nrm_inf, value);
        }
        p.g_nrm2 = sqrt(p.g_nrm2_2);

//         printf("norms = %e %e %e\n", p.g_nrm_1, p.g_nrm2, p.g_nrm_inf);
        m_state.g_count += 1;
    }

    void
    dual_x(Point & p, Point & dual_p) override
    {
        m_problem.dual_x(p, dual_p);
    }

    void
    dual_f(Point & dual_p) override
    {
        m_problem.dual_f(dual_p);
    }

private:
    Problem& m_problem;
    State& m_state;
};

// FIXME make public, add to_string(), add reason + state tuple as optimize result
enum class ExitReason
{
    NoRelaxation,
    Iterations,
    Time,
    FunctionValue,
    GradientNormValue,
};

String
to_string(ExitReason reason)
{
    // TODO fix text
    switch (reason)
    {
        case ExitReason::NoRelaxation:;
        {
            static const String s("no relaxation");
            return s;
        }

        case ExitReason::Iterations:
        {
            static const String s("max iteration");
            return s;
        }

        case ExitReason::Time:
        {
            static const String s("max time");
            return s;
        }

        case ExitReason::FunctionValue:
        {
            static const String s("function value");
            return s;
        }

        case ExitReason::GradientNormValue:
        {
            static const String s("gradient norm value");
            return s;
        }
    }

    // -Wreturn-type warning fix
    return {};
}

Method::Method(String name, ProblemPropertyFlags properties)
    : m_name(name)
    , m_properties(properties)
{
}

void
Method::before(Problem& /*problem*/, Point& /*point*/)
{
}

void
Method::after(Problem& /*problem*/, Point& /*point*/)
{
}

void
Method::log(Logger& /*logger*/, LoggerMode /*mode*/) const
{
}

void
Method::optimize(Problem& original_problem, Point& point, const MethodSettings& settings, State& state)
{
    auto print_error = [this, &original_problem](ProblemProperty p)
    {
        printf("Selected method '%s' requires %s, but the problem '%s' doesn't provide it\n",
                m_name.c_str(),
                to_string(p).c_str(), // TODO bright yellow ?
                original_problem.name().c_str());
    };

    for (auto p : m_properties)
    {
        if (original_problem.has(p) == false)
        {
            print_error(p);
            return;
        }
    }

    if (use(ProblemProperty::LipschitzConstant) && std::isnan(original_problem.L()))
    {
        print_error(ProblemProperty::LipschitzConstant);
        return;
    }

    if (settings.resume == false)
    {
        // FIXME add reset() method into State class ?
        state.f_count = state.g_count = 0;
        state.t_total = 0.0;
        state.iter_total = 0;
    }

    auto stdout_logger = Logger::stdout(*this);
    auto csv_logger = Logger::csv(*this, settings.resume);

    auto problem = WrappedProblem(original_problem, state);

    auto t_0 = chrono::now();
    auto t_p = t_0;
    size_t iter = state.iter_total; // FIXME remove variable ?

    problem.f(point);
    if (problem.has(ProblemProperty::Gradient))
    {
        problem.df(point);
    }

    before(problem, point);

    auto t_i = chrono::s(t_0) + state.t_total;

    log_all(stdout_logger, LoggerMode::Header, point, state, t_i, iter);
    stdout_logger.put_new_line();

    log_all(stdout_logger, LoggerMode::Value, point, state, t_i, iter);
    stdout_logger.put_new_line();

    if (settings.resume == false)
    {
        log_all(csv_logger, LoggerMode::Header, point, state, t_i, iter);
        log_all(csv_logger, LoggerMode::Value, point, state, t_i, iter);
    }

    auto exit_reason = ExitReason::NoRelaxation;

    double Z = 1e-1;
    bool zb = false;

    while (true)
    {
        // check start value of f and g_nrm2 before running the method
        if (point.f < settings.f_min)
        {
//             printf("\nF_MIN %e %e\n", point.f, settings.f_min);
            exit_reason = ExitReason::FunctionValue;
            break;
        }

        if (use(ProblemProperty::Gradient) && (point.g_nrm2 < settings.g_nrm2_min))
        {
            exit_reason = ExitReason::GradientNormValue;
            break;
        }

        auto iteration_result = iteration(problem, point, iter);
        iter += 1;

        auto t = chrono::now();
        auto t_i = chrono::s(t_0, t) + state.t_total;
        auto t_p_i = chrono::s(t_p, t);

        // TODO add settings.log_iterval_time ?
        log_all(csv_logger, LoggerMode::Value, point, state, t_i, iter);

        if (point.g_nrm_inf <= Z)
        {
            zb = true;
            Z *= 0.1;
        }

        if (zb || (t_p_i > settings.print_iterval_time))
        {
            t_p = t;

            log_all(stdout_logger, LoggerMode::Value, point, state, t_i, iter);

            // FIXME add log flushing interval parameter ?
            // FIXME flush stdout only when print_iterval_time is not too small ?
            stdout_logger.flush();
            csv_logger.flush();

            // FIXME add logging level
//            stdout_logger.put_new_line(); // FIXME

            if (zb)
            {
                stdout_logger.put_new_line();
                zb = false;
            }
        }

        if (iteration_result == false)
        {
            break;
        }

        if (iter >= settings.iter_max + state.iter_total) // FIXME integer overflow (see default value of iter_max)
        {
            exit_reason = ExitReason::Iterations;
            break;
        }

        if (t_i >= settings.time_max + state.t_total)
        {
            exit_reason = ExitReason::Time;
            break;
        }
    }

    after(problem, point);

//     t_i = chrono::s(t_0);
//     log_all(stdout_logger, LoggerMode::Value, point, state, t_i + state.t_total, iter);
    state.t_total += chrono::s(t_0);
    state.iter_total = iter;
    log_all(stdout_logger, LoggerMode::Value, point, state, state.t_total, state.iter_total);
    stdout_logger.put_new_line();

    printf("Exit by: %s\n", to_string(exit_reason).c_str());

    // FIXME add this into logger's destructor
    csv_logger.flush();
}

void
Method::log_all(Logger& logger, LoggerMode mode, const Point& point, const State& state, double time, size_t iter)
{
    if (logger.device() == LoggerDevice::Stdout)
    {
        logger.put_carriage_return();
    }
//     else
//     {
//         return;
//     }

    logger.put_method_name(mode);
    logger.put_time(mode, time);
    logger.put_iter(mode, iter);
    logger.put_point(mode, point);
    logger.put_state(mode, state);

    log(logger, mode);

    if (logger.device() == LoggerDevice::CsvFile)
    {
        logger.put_new_line();
    }
}

}

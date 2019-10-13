#pragma once

#include "types.hpp"

namespace t_opt
{

class Problem;

class Logger;
enum class LoggerMode : uint8_t;

class LineSearchMethod;

struct MethodSettings
{
    size_t iter_max = limits<size_t>::max();
    double time_max = 3600.0;

    double f_min = limits<double>::lowest();
    double g_nrm2_min = 1e-6;

    double print_iterval_time = 0.1;

    bool resume = false;
};

// FIXME rename
struct State
{
    size_t f_count = 0;
    size_t g_count = 0;

    double t_total = 0.0;
    size_t iter_total = 0;

    State&
    operator+=(const State& other)
    {
        f_count += other.f_count;
        g_count += other.g_count;

        return *this;
    }
};

class Method
{
public:
    inline const String&
    name() const
    {
        return m_name;
    }

    inline bool
    use(ProblemProperty p) const
    {
        return bool(m_properties & p);
    }

    void
    optimize(Problem& problem, Point& point, const MethodSettings& settings, State& state);

protected:
    Method(String name, ProblemPropertyFlags properties = ProblemPropertyFlags());

    virtual void
    before(Problem& problem, Point& point);

    virtual void
    after(Problem& problem, Point& point);

    virtual bool
    iteration(Problem& problem, Point& point, size_t iter) = 0;

    virtual void
    log(Logger& logger, LoggerMode mode) const;

    void
    log_all(Logger& logger, LoggerMode mode, const Point& point, const State& state, double time, size_t iter);

    String m_name;
    ProblemPropertyFlags m_properties;
};

}

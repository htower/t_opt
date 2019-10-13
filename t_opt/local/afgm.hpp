#pragma once

#include "core/method.hpp"

namespace t_opt
{

namespace local
{

class AFGM : public Method
{
public:
    explicit
    AFGM(LineSearchMethod& ls);

protected:
    void
    before(Problem& problem, Point& point) override;

    bool
    iteration(Problem& problem, Point& point, size_t iter) override;

    void
    log(Logger& logger, LoggerMode mode) const override;

private:
    LineSearchMethod& ls;
    double ls_start_step = 1.0;
    double ls_step = 1.0;

    double alpha;
    double alpha_sum;

    Point x;
    Point z;
    DVector zy;

    Point dual_p;
    Point dual_p_sum;
    double pd_delta;
};

}

}

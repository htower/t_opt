#pragma once

#include "core/method.hpp"

namespace t_opt
{

namespace local
{

class FGM : public Method
{
public:
    FGM();

protected:
    void
    before(Problem& problem, Point& point) override;

    void
    log(Logger& logger, LoggerMode mode) const override;

    bool
    iteration(Problem& problem, Point& point, size_t iter) override;

private:
    DVector x;
    DVector y;
    double step;

    Point dual_p;
    Point dual_p_sum;
    double pd_delta;
    double alpha_sum;
};

}

}


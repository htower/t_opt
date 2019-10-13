#pragma once

#include "core/method.hpp"

namespace t_opt
{

namespace local
{

class LBFGS : public Method
{
public:
    explicit
    LBFGS(uint32_t m, LineSearchMethod& ls);

protected:
    void
    before(Problem& problem, Point& point) override;

    void
    log(Logger& logger, LoggerMode mode) const override;

    bool
    iteration(Problem& problem, Point& point, size_t iter) override;

private:
    uint32_t m;

    LineSearchMethod& ls;
    double ls_step;
    double ls_start_step = 1.0; // FIXME make public ?

    std::vector<DVector> l_s;
    std::vector<DVector> l_y;

    DVector l_ys;
    DVector l_alpha;

    DVector x_p;
    DVector g_p;
    DVector dir;

    uint32_t l_end;
};

}

}

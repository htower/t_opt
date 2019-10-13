#pragma once

#include "core/method.hpp"

namespace t_opt
{

namespace local
{

class GDM : public Method
{
public:
    explicit
    GDM(LineSearchMethod& ls);

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
};

}

}

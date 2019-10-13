#pragma once

#include "core/method.hpp"

namespace t_opt
{

namespace local
{

class UFGM : public Method
{
public:
    UFGM(double epsilon);

protected:
    void
    before(Problem& problem, Point& point) override;

    bool
    iteration(Problem& problem, Point& point, size_t iter) override;

    void
    log(Logger& logger, LoggerMode mode) const override;

    double epsilon;

private:
    DVector v_k;
    DVector dyx;

    Point x_kp1;
    Point y_kp1;
    DVector z_kp1;

    double alpha_k;
    double alpha_kp1;
    double l_k;
    double l_kp1;


    Point dual_p;
    Point dual_p_sum;
    double pd_delta;
    double alpha_sum;
};

}

}

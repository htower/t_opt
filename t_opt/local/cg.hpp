#pragma once

#include "core/method.hpp"

namespace t_opt
{

namespace local
{

enum class CgVariant : uint8_t
{
    // Hestenes-Stiefel
    // (M.R. Hestenes and E.L. Stiefel, Methods of conjugate gradients for solving linear systems,
    // J.Research Nat. Bur. Standards, 49 (1952), pp.409-436)
    HS,

    // Fletcher - Reeves
    // (R. Fletcher and C. Reeves, Function minimization by conjugate gradients,
    // Comput. J., 7 (1964), pp.149-154)
    FR,

    // Polak - Ribière - Polyak
    // (E. Polak and G. Ribière, Note sur la convergence de directions conjuguée,
    // Rev. Francaise Informat Recherche Operationelle, 3e Année 16 (1969), pp.35-43)
    // (B.T. Polyak, The conjugate gradient method in extreme problems. USSR Comp. Math.
    // Math. Phys., 9 (1969), pp.94-112)
    PRP,

    // Polak - Ribière - Polyak plus
    // (M.J.D. Powell, Nonconvex minimization calculations and the conjugate gradient method.
    // Numerical Analysis (Dundee, 1983) Lecture Notes in mathematics, vol.1066, Springer-Verlag,
    // Berlin, 1984, pp.122-141)
    PRPplus,

    // Conjugate Descent - Fletcher
    // (R. Fletcher, Practical Methods of Optimization, vol. 1: Unconstrained Optimization,
    // John Wiley & Sons, New York, 1987)
    CD,

    // Liu - Storey
    // (Y. Liu, and C. Storey, Efficient generalized conjugate gradient algorithms,
    // Part 1: Theory. JOTA, 69 (1991), pp.129-137)
    LS,

    // Dai - Yuan
    // (Y.H. Dai and Y. Yuan, An efficient hybrid conjugate gradient method for unconstrained
    // optimization, Ann. Oper. Res., 103 (2001), pp. 33-47)
    DY,
};

String
to_string(CgVariant variant);

class CG : public Method
{
public:
    CG(CgVariant variant, LineSearchMethod& ls);

protected:
    void
    before(Problem& problem, Point& point) override;

    void
    log(Logger& logger, LoggerMode mode) const override;

    bool
    iteration(Problem& problem, Point& point, size_t iter) override;

private:
    void
    calc_dir(const DVector& g);

    void
    before_step(const Point& point);

    void
    after_step(const Point& point);

    CgVariant variant;

    LineSearchMethod& ls;
    double ls_step;
    double ls_start_step = 1.0; // FIXME make public ?

    double g_nrm2;
    double g_nrm2_prev;

    DVector dir;
    DVector dir_normalized;
    DVector dir_prev;
    bool dir_reset;

    DVector y;
    DVector s;
    DVector g_prev;

//     bool use_y;
//     bool use_s;
//     bool use_g_prev;
//     bool dir_by_s;
};

}

}

#include "core/blas.hpp"

#include "line_search/h_simple.hpp"
#include "line_search/parabolic.hpp"

#include "local/afgm.hpp"
#include "local/agmsdr.hpp"
#include "local/cg.hpp"
#include "local/gdm.hpp"
#include "local/lbfgs.hpp"
#include "local/fgm.hpp"
#include "local/ufgm.hpp"

#include "problems/quadratic.hpp"

#include "problems/transport/smvsdm2.hpp"
#include "problems/transport/tsdm.hpp"
#include "problems/transport/lpsdm.hpp"

#include "utils.hpp"

void
run()
{
    using namespace t_opt;

    auto path = "/home/htower/science/data/TNTP/TransportationNetworks";
    auto name =
//         "test0";
//         "test1"; // OK, works well
//         "test_s";
        "SiouxFalls"; // 24 / 24
//         "Anaheim"; // 416 / 38

//         "Chicago-Sketch"; // 933 / 387 // PARSE ERROR
//         "Winnipeg"; // 1040 / 147 // very slow convergence
//         "Barcelona";  // 1020 / 110 // very slow convergence

//         "Berlin-Tiergarten"; // 361 / 26 // FAIL at start with mu = 1e-2 (nan) - WTF????
//         "Berlin-Prenzlauerberg-Center"; // 352 / 38 // FAIL at start with mu = 1e-2 (nan) - WTF????
//         "Berlin-Mitte-Center"; // 398 / 36 // FAIL at start with mu = 1e-2 (nan) - WTF????
//         "Austin";

    auto problem = transport::SmVSDM2(path, name);
    double mu = 1e1;
    problem.set_mu(mu);

    std::string asd;
    asd.shrink_to_fit();

//     auto problem = transport::TSDM(path, name);
//     auto problem = transport::LPSDM(path, name);
//     problem.K = 1e0;

//     problem.apply_flow(1.1, true);
//     problem.apply_total_flow(1.0);

    auto point = Point(problem);
    auto state = State();
    auto settings = MethodSettings();
//     settings.iter_max = 1000;
    settings.g_nrm2_min = 1e-5;
//     settings.g_nrm2_min = 0.0;
//     settings.time_max = 5.0;
//     settings.time_max = 120.0;
//     settings.iter_max = 5e6;
//     settings.iter_max = 1;
//     settings.iter_max = 1e2;
//     settings.iter_max = 100;
//     settings.iter_max = 300;
//     settings.iter_max = 80;
//     settings.f_min = 1e-5;
//     settings.print_iterval_time = 0.0;
//     settings.g_nrm2_min = 1e-3;

    auto ls = line_search::HSimple();
//     auto ls = line_search::Parabolic(2, true);
//     auto ls = line_search::Parabolic(3, false);

    auto method =
//         local::FGM();
//         local::UFGM(1e-4);
        local::AFGM(ls);
//         local::UFGM(1e-8); // much faster on SiouxFalls - WHY ????

//     auto method = local::AGMsDR(ls, 1e-4);
//     auto method = local::GDM(ls);
//     auto method = local::CG(local::CgVariant::PRP, ls);
//     auto method = local::LBFGS(3, ls);

//     speed_test(problem, point, 1e4); return 0;

//     for (int i = 0; i < 2; ++i)
//     {
//         point.x.setRandom();
//         g_test(problem, point);
//     }
//     return;

//     Point dual_p(problem);
//     dual_p.x.resize(1e6); // FIXME
//     {
//         problem.f(point);
//         problem.dual_x(point, dual_p);
//         problem.dual_f(dual_p);
//         printf("\n\nF = %e; DUAL_F = %e; DELTA = %e\n\n", point.f, dual_p.f, point.f + dual_p.f);
//     }

//     problem.restore_flow(point);
    method.optimize(problem, point, settings, state);
//     problem.emoe(point);

    return;
    problem.restore_flow(point);

//     {
//         problem.dual_x(point, dual_p);
//         problem.dual_f(dual_p);
//         printf("\n\nF = %e; DUAL_F = %e; DELTA = %e\n\n", point.f, dual_p.f, point.f + dual_p.f);
//
// //         dual_p.x[0] = 100.0;
// //         dual_p.x[1] = 200.0;
// //         dual_p.x[2] = 200.0;
//
// //         problem.set_mu(0.0);
// //         problem.f(point);
// //         problem.dual_f(dual_p);
// //         printf("\n\nF = %e; DUAL_F = %e; DELTA = %e\n\n", point.f, dual_p.f, point.f + dual_p.f);
//     }

    settings.resume = true;
    for (int i = 0; i < 0; ++i)
    {
        mu *= 0.1;
        problem.set_mu(mu);
//         printf("K = %g\n", problem.K);
        method.optimize(problem, point, settings, state);
//         problem.emoe(point);
// //         problem.restore_flow(point);
    }

//     local::UFGM(1e-4).optimize(problem, point, settings);
}

int
main(int /*argc*/, char** /*argv*/)
{
    try
    {
        run();
        return 0;
    }
    catch (std::exception& e)
    {
        printf("%s\n", e.what());
        return -1;
    }
}

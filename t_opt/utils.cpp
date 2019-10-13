#include "utils.hpp"

#include "core/blas.hpp"
#include "core/chrono.hpp"
#include "core/types.hpp"
#include "core/problem.hpp"

namespace t_opt
{

void
speed_test(Problem& problem, Point& point, long count)
{
    {
        double f_sum = 0.0;
        auto t_0 = chrono::now();
        for (long i = 0; i < count; ++i)
        {
            problem.f(point);
            f_sum += point.f;
        }
        auto f_t = chrono::ms(t_0);
        printf("f time = %.3f ms. f_sum = %e\n", f_t, f_sum);
    }

    {
        double g_nrm2_sum = 0.0;
        auto t_0 = chrono::now();
        for (long i = 0; i < count; ++i)
        {
            problem.df(point);
            g_nrm2_sum += point.g_nrm2;
        }
        auto g_t = chrono::ms(t_0);
        printf("g time = %.3f ms. g_nrm2_sum = %e\n", g_t, g_nrm2_sum);
    }
}

void
g_test(Problem& problem, Point& point)
{
    const double h_max = 1e-1;
    const double h_min = limits<double>::epsilon();

    problem.f(point);
    blas::set_zero(point.g); // FIXME use wrapped problem here
    problem.df(point);

    const auto f0 = point.f;
    DVector g0(problem.size());
    blas::copy(point.g, g0);

    DVector g1(problem.size());
    DVector g2(problem.size());
    DVector g3(problem.size());

    auto delta = [](const DVector& x, const DVector& y)
    {
        double res = 0.0;
        for (int i = 0; i < x.size(); ++i)
        {
            const auto d = x[i] - y[i];
            res += d * d;
        }

        return std::sqrt(res);
    };

    printf("h       g(+) - g     g(-) - g     g(+-) - g\n");
    for (double h = h_max; h >= h_min; h *= 0.1)
    {
        for (size_t i = 0; i < problem.size(); ++i)
        {
            const auto x_save = point.x[i];

            // g+ =============================================================
            point.x[i] = x_save + h;
            problem.f(point);
            g1[i] = (point.f - f0) / h;

            // g-
            point.x[i] = x_save - h;
            problem.f(point);
            g2[i] = (f0 - point.f) / h;

            // g+- ============================================================
            point.x[i] = x_save + h;
            problem.f(point);
            g3[i] = point.f;

            point.x[i] = x_save - h;
            problem.f(point);
            g3[i] -= point.f;
            g3[i] /= 2.0 * h;

            // ================================================================
            point.x[i] = x_save;
        }

        printf("%5.0e : %e %e %e\n",
            (double)h,
            (double)delta(g0, g1),
            (double)delta(g0, g2),
            (double)delta(g0, g3));
    }
}

}

#pragma once

#include "src/sdm.hpp"

namespace transport
{

class SmVSDM2 : public SDM
{
public:
    SmVSDM2(const String& path, const String& name);

    void
    f(Point& p) override;

    void
    df(Point& p) override;

    void
    set_mu(double mu);

    // FIXME remove this, replaced with dual_x()
    void
    restore_flow(Point& point);

    void
    dual_x(Point & p, Point & dual_p) override;

    void
    dual_f(Point & dual_p) override;

private:
    double mu;
};

}

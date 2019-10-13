#pragma once

#include "src/sdm.hpp"

namespace transport
{

class LPSDM : public SDM
{
public:
    LPSDM(const String& path, const String& name);

    void
    f(Point& p) override;

    void
    df(Point& p) override;

    void
    restore_flow(Point& point);

    void
    emoe(Point& point) override;

    double K = 1.0;

private:
    size_t T_size;
};

}

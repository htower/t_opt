#pragma once

#include "src/sdm.hpp"

namespace transport
{

struct TSDM : public SDM
{
    TSDM(const String& path, const String& name);

    void
    f(Point& p) override;

    void
    df(Point& p) override;

    void
    restore_flow(Point& point);
};

}

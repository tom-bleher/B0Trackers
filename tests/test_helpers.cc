#include "B0TrackersHelpers.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

int main() {
    using namespace b0trk;

    {
        const auto [p, ok] = momentumFromQOverP(0.25);
        assert(ok);
        assert(std::abs(p - 4.0) < 1e-12);
    }
    {
        const auto [p, ok] = momentumFromQOverP(0.0);
        assert(!ok);
        assert(std::isnan(p));
    }
    {
        const auto [p, ok] = momentumFromQOverP(std::numeric_limits<double>::quiet_NaN());
        assert(!ok);
        assert(std::isnan(p));
    }

    assert(stationFromLayerId(1, true) == 1);
    assert(stationFromLayerId(2, true) == 1);
    assert(stationFromLayerId(3, true) == 2);
    assert(stationFromLayerId(8, true) == 4);
    assert(stationFromLayerId(1, false) == 1);
    assert(stationFromLayerId(4, false) == 4);
    assert(stationFromLayerId(0, true) == -1);

    {
        const auto [st, side] = parseLayerName("B0Tracker_layer1_front_P");
        assert(st == 1);
        assert(side == 1);
    }
    {
        const auto [st, side] = parseLayerName("B0Tracker_layer2_back_P");
        assert(st == 2);
        assert(side == 0);
    }
    {
        const auto [st, side] = parseLayerName("B0Tracker_layer3_P");
        assert(st == 3);
        assert(side == -1);
    }
    {
        const auto [st, side] = parseLayerName("not_a_layer");
        assert(st == -1);
        assert(side == -1);
    }

    {
        const double wrapped = wrapPi(3.2);
        assert(wrapped < 3.14159265358979323846 && wrapped > -3.14159265358979323846);
        assert(std::abs(std::cos(wrapped) - std::cos(3.2)) < 1e-12);
    }
    assert(std::abs(wrapPi(0.1) - 0.1) < 1e-12);

    std::cout << "B0TrackersHelpers tests passed\n";
    return 0;
}

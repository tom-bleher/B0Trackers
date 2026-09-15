// Checks must survive -DNDEBUG: a Release build silently compiles out C assert,
// and a test that asserts nothing still reports success.

#include "B0TrackersHelpers.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace {

int g_failures = 0;

bool g_quiet = false;

void check(bool condition, const char* expr, const char* file, int line) {
    if (!condition) {
        ++g_failures;
        if (!g_quiet) {
            std::cerr << file << ":" << line << ": CHECK failed: " << expr << "\n";
        }
    }
}

bool closeTo(double a, double b, double tol = 1e-12) { return std::abs(a - b) <= tol; }

} // namespace

#define CHECK(expr) check(static_cast<bool>(expr), #expr, __FILE__, __LINE__)

int main() {
    using namespace b0trk;

    {
        const auto [p, ok] = momentumFromQOverP(0.25);
        CHECK(ok);
        CHECK(closeTo(p, 4.0));
    }
    {
        const auto [p, ok] = momentumFromQOverP(-0.25);
        CHECK(ok);
        CHECK(closeTo(p, 4.0));
    }
    {
        const auto [p, ok] = momentumFromQOverP(0.0);
        CHECK(!ok);
        CHECK(std::isnan(p));
    }
    {
        const auto [p, ok] = momentumFromQOverP(std::numeric_limits<double>::quiet_NaN());
        CHECK(!ok);
        CHECK(std::isnan(p));
    }

    CHECK(stationFromLayerId(1, true) == 1);
    CHECK(stationFromLayerId(2, true) == 1);
    CHECK(stationFromLayerId(3, true) == 2);
    CHECK(stationFromLayerId(8, true) == 4);
    CHECK(stationFromLayerId(1, false) == 1);
    CHECK(stationFromLayerId(4, false) == 4);
    CHECK(stationFromLayerId(0, true) == -1);

    {
        const auto [st, side] = parseLayerName("B0Tracker_layer1_front_P");
        CHECK(st == 1);
        CHECK(side == 1);
    }
    {
        const auto [st, side] = parseLayerName("B0Tracker_layer2_back_P");
        CHECK(st == 2);
        CHECK(side == 0);
    }
    {
        const auto [st, side] = parseLayerName("B0Tracker_layer3_P");
        CHECK(st == 3);
        CHECK(side == -1);
    }
    {
        const auto [st, side] = parseLayerName("not_a_layer");
        CHECK(st == -1);
        CHECK(side == -1);
    }

    {
        const double wrapped = wrapPi(3.2);
        CHECK(wrapped < 3.14159265358979323846 && wrapped > -3.14159265358979323846);
        CHECK(closeTo(std::cos(wrapped), std::cos(3.2)));
    }
    CHECK(closeTo(wrapPi(0.1), 0.1));

    CHECK(kSchemaVersion == 3);
    CHECK(preferredEstimateKind(true, false, false) == kEstimatePredicted);
    CHECK(preferredEstimateKind(true, true, false) == kEstimateFiltered);
    CHECK(preferredEstimateKind(true, true, true) == kEstimateSmoothed);
    CHECK(preferredEstimateKind(false, false, false) == kEstimateNone);
    CHECK(closeTo(nativeTimeToNs(299.792458, 299.792458), 1.0));
    CHECK(std::isnan(nativeTimeToNs(1.0, 0.0)));

    // The test harness itself must be able to fail; a check that never trips is
    // indistinguishable from one that was compiled away.
    {
        const int before = g_failures;
        g_quiet            = true;
        CHECK(1 == 2);
        g_quiet                 = false;
        const bool harnessWorks = (g_failures == before + 1);
        g_failures              = before;
        CHECK(harnessWorks);
    }

    if (g_failures != 0) {
        std::cerr << "B0TrackersHelpers tests FAILED (" << g_failures << " check(s))\n";
        return EXIT_FAILURE;
    }
    std::cout << "B0TrackersHelpers tests passed\n";
    return EXIT_SUCCESS;
}

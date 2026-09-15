#pragma once

#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

// Pure helpers (no JANA/DD4hep). Used by the plugin and by standalone tests.

namespace b0trk {

inline constexpr int kSchemaVersion = 3;

inline constexpr int kMapUnresolved = 0;
inline constexpr int kMapExact      = 1;
inline constexpr int kMapFallback   = 2;

enum StateEstimateKind {
    kEstimateNone      = 0,
    kEstimatePredicted = 1,
    kEstimateFiltered  = 2,
    kEstimateSmoothed  = 3,
};

inline double quietNaN() { return std::numeric_limits<double>::quiet_NaN(); }

inline int preferredEstimateKind(bool hasPredicted, bool hasFiltered, bool hasSmoothed) {
    if (hasSmoothed) {
        return kEstimateSmoothed;
    }
    if (hasFiltered) {
        return kEstimateFiltered;
    }
    if (hasPredicted) {
        return kEstimatePredicted;
    }
    return kEstimateNone;
}

inline double nativeTimeToNs(double nativeTime, double actsNsUnit) {
    if (!std::isfinite(nativeTime) || !std::isfinite(actsNsUnit) || actsNsUnit == 0.0) {
        return quietNaN();
    }
    return nativeTime / actsNsUnit;
}

// q/p == 0 is unbounded momentum, not p = 0.
inline std::pair<double, bool> momentumFromQOverP(double qOverP) {
    if (!std::isfinite(qOverP) || qOverP == 0.0) {
        return {quietNaN(), false};
    }
    return {std::abs(1.0 / qOverP), true};
}

inline double wrapPi(double angle) {
    constexpr double pi = 3.14159265358979323846;
    constexpr double twoPi = 2.0 * pi;
    while (angle > pi) {
        angle -= twoPi;
    }
    while (angle < -pi) {
        angle += twoPi;
    }
    return angle;
}

// Official 4-layer: station = layer id.
// Per-plane 8-layer: layer = 2*(station-1)+(1=back|2=front) => station = (layer+1)/2.
inline int stationFromLayerId(int layer, bool perPlaneFrontBack) {
    if (layer <= 0) {
        return -1;
    }
    return perPlaneFrontBack ? (layer + 1) / 2 : layer;
}

// Layer DetElement names:
//   B0Tracker_layer{station}_{front|back}_P  (realistic)
//   B0Tracker_layer{N}_P                     (official; N is the station)
// Returns {station, side} with side 1=front, 0=back, -1=unknown.
inline std::pair<int, int> parseLayerName(const std::string& name) {
    const auto pos = name.find("layer");
    if (pos == std::string::npos || pos + 5 >= name.size()) {
        return {-1, -1};
    }
    std::size_t i = pos + 5;
    if (!std::isdigit(static_cast<unsigned char>(name[i]))) {
        return {-1, -1};
    }
    int n = 0;
    while (i < name.size() && std::isdigit(static_cast<unsigned char>(name[i]))) {
        n = n * 10 + (name[i++] - '0');
    }
    int side = -1;
    if (name.find("front", i) != std::string::npos) {
        side = 1;
    } else if (name.find("back", i) != std::string::npos) {
        side = 0;
    }
    return {n, side};
}

} // namespace b0trk

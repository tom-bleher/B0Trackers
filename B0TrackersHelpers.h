#pragma once

#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

// Pure helpers (no JANA/DD4hep). Used by the plugin and by standalone tests.

namespace b0trk {

inline constexpr int kSchemaVersion = 4;

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

// Per-seed CKF failure diagnostics (schema 4).
//
// The status codes mirror eicrecon::b0counters::ckfdiag
// (EICrecon/src/algorithms/tracking/B0ReconstructionCounters.h), which owns
// the canonical definition written into the unfiltered ACTS containers by
// CKFTracking. The values are duplicated here -- rather than included --
// because the standalone helper test builds without EICrecon headers, so a
// mismatch would surface as wrong diagnostics, not a build error. Keep them
// in sync by hand.
namespace ckfdiag {

inline constexpr unsigned kAccepted = 0;
inline constexpr unsigned kNoValidMeasurement = 1;
inline constexpr unsigned kTooFewHits = 2;
inline constexpr unsigned kTooFewStations = 3;
inline constexpr unsigned kSmoothingFailed = 4;
inline constexpr unsigned kExtrapolationFailed = 5;
inline constexpr unsigned kFindFailed = 6;
inline constexpr unsigned kNoCandidates = 7;

inline constexpr unsigned kFindErrNone = 0;
inline constexpr unsigned kFindErrCkf = 1;
inline constexpr unsigned kFindErrPropagation = 2;

// Per-seed CKF outcome stage: -1 means the unfiltered ACTS containers were
// unavailable, so the stage is unknowable rather than failed.
inline constexpr int kStageUnknown = -1;
inline constexpr int kStageAccepted = 0;
inline constexpr int kStageAllRejected = 1;
inline constexpr int kStageFindFailed = 2;
inline constexpr int kStageFindEmpty = 3;

inline unsigned packFindError(unsigned errorClass, unsigned errorValue) {
    return errorClass * 1000u + (errorValue % 1000u);
}

inline unsigned findErrorClass(unsigned packed) { return packed / 1000u; }

inline unsigned findErrorValue(unsigned packed) { return packed % 1000u; }

// Best-candidate selection: an accepted candidate beats any rejected one,
// then most measurements, fewest holes, lowest status. With haveBest == false
// the candidate is selected unconditionally.
inline bool isBetterCandidate(bool accepted, int nMeas, int nHoles, unsigned status,
                              bool bestAccepted, int bestMeas, int bestHoles,
                              unsigned bestStatus, bool haveBest) {
    if (!haveBest) {
        return true;
    }
    if (accepted != bestAccepted) {
        return accepted;
    }
    if (nMeas != bestMeas) {
        return nMeas > bestMeas;
    }
    if (nHoles != bestHoles) {
        return nHoles < bestHoles;
    }
    return status < bestStatus;
}

} // namespace ckfdiag

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

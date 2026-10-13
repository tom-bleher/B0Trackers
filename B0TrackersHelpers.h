#pragma once

#include <cctype>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <string>
#include <utility>

// Pure helpers (no JANA/DD4hep). Used by the plugin and by standalone tests.

namespace b0trk {

inline constexpr int kSchemaVersion = 10;

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

// Innovation compatibility for one track state: r = m - H x_pred,
// S = V + H C_pred H^T with (v00, v01, v11) the measurement covariance
// and (c00, c01, c11) the projected predicted covariance. Mirrors the
// per-state computation in B0Trackers.cc so the filtered and unfiltered
// dumps share one definition. measDim follows ACTS calibratedSize (1 or 2
// for B0 pixels); anything else yields NaNs. Pulls use the S diagonal even
// when S is singular (then chi2 is NaN), matching what the CKF measurement
// selector computed from the same quantities.
struct InnovCompat {
    double pull0 = quietNaN();
    double pull1 = quietNaN();
    double chi2 = quietNaN();
};

inline InnovCompat innovCompat(int measDim, double r0, double r1, double v00, double v01,
                               double v11, double c00, double c01, double c11) {
    InnovCompat out;
    if (measDim == 1) {
        const double s00 = v00 + c00;
        if (s00 > 0.0 && std::isfinite(s00)) {
            out.pull0 = r0 / std::sqrt(s00);
            out.chi2  = r0 * r0 / s00;
        }
        return out;
    }
    if (measDim == 2) {
        const double s00 = v00 + c00;
        const double s01 = v01 + c01;
        const double s11 = v11 + c11;
        if (!std::isfinite(s00) || !std::isfinite(s01) || !std::isfinite(s11)) {
            return out;
        }
        if (s00 <= 0.0 || s11 <= 0.0) {
            return out;
        }
        out.pull0 = r0 / std::sqrt(s00);
        out.pull1 = r1 / std::sqrt(s11);
        const double det = s00 * s11 - s01 * s01;
        if (det > 0.0 && std::isfinite(det)) {
            out.chi2 = (s11 * r0 * r0 - 2.0 * s01 * r0 * r1 + s00 * r1 * r1) / det;
        }
        return out;
    }
    return out;
}

// Local single-propagation cells (schema 10): a bounded straight-ray
// intersection (BoundaryTolerance::None) crossed with the SurfaceReached
// propagation from a few mm upstream of a SimTrackerHit.
namespace tsloc {

inline constexpr int kLocalOk = 0;             // propagation reached
inline constexpr int kLocalTargetBug = 1;      // ray hits finite surface, propagation misses
inline constexpr int kLocalBoundsMapping = 2;  // ray misses finite surface, propagation misses
inline constexpr int kLocalInconsistent = 3;   // ray misses, propagation reaches

inline int classifyLocal(bool ixReachable, bool propReached) {
    if (propReached) {
        return ixReachable ? kLocalOk : kLocalInconsistent;
    }
    return ixReachable ? kLocalTargetBug : kLocalBoundsMapping;
}

// Charge numerator for common PDGs; 0 = neutral or unknown.
inline int chargeFromPdg(int pdg) {
    const int a = (pdg < 0) ? -pdg : pdg;
    int q = 0;
    if (a == 11 || a == 13) {
        q = -1;
    } else if (a == 211 || a == 321 || a == 2212) {
        q = +1;
    } else if (a == 2112 || a == 22) {
        q = 0;
    } else {
        return 0;
    }
    if (pdg < 0 && q != 0) {
        q = -q;
    }
    return q;
}

} // namespace tsloc

// Layer DetElement names:
//   B0Tracker_station{N}_{front|back}_P  (realistic, eic/epic b0-realistic-geometry)
//   B0Tracker_layer{N}_{front|back}_P    (earlier realistic prototypes)
//   B0Tracker_layer{N}_P                 (official; N is the station)
// Returns {station, side} with side 1=front, 0=back, -1=unknown. The side
// comes from the name, never from the cellID layer parity, because which face
// gets the odd layer id is a geometry convention that has changed before.
inline std::pair<int, int> parseLayerName(const std::string& name) {
    std::size_t i = std::string::npos;
    for (const char* token : {"station", "layer"}) {
        const auto pos = name.find(token);
        const std::size_t len = std::char_traits<char>::length(token);
        if (pos != std::string::npos && pos + len < name.size() &&
            std::isdigit(static_cast<unsigned char>(name[pos + len]))) {
            i = pos + len;
            break;
        }
    }
    if (i == std::string::npos) {
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

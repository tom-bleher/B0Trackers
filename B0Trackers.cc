#include "B0Trackers.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <initializer_list>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <Eigen/Cholesky>

#include <JANA/JException.h>
#include <JANA/Services/JGlobalRootLock.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TGeoMatrix.h>
#include <TTree.h>

#include <Acts/ActsVersion.hpp>
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/ProxyAccessor.hpp>
#include <Acts/EventData/TrackContainer.hpp>
#include <Acts/EventData/TrackStateType.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>

#include <edm4eic/MCRecoTrackParticleAssociation.h>
#include <edm4eic/MCRecoTrackerHitAssociation.h>
#include <edm4eic/Measurement2D.h>
#include <edm4eic/RawTrackerHit.h>
#include <edm4eic/Track.h>
#include <edm4eic/TrackParameters.h>
#include <edm4eic/TrackSeed.h>
#include <edm4eic/TrackerHit.h>
#include <edm4eic/Trajectory.h>

#include <DD4hep/Objects.h>

#include <services/rootfile/RootFile_service.h>
#include <algorithms/tracking/ActsGeometryProvider.h>
#include <services/geometry/acts/ACTSGeo_service.h>
#include <services/log/Log_service.h>

// B0Trackers/hits branch map (schema 3):
//   trk_*                  = B0TrackerCKFTruthSeeded* (filtered), IP-perigee
//   ckf_trk_*              = B0TrackerCKF* (filtered), IP-perigee
//   *_oracle_best_*        = min |p_reco-p_sel| (cheat; also aliased as best_trk_*)
//   *_truth_matched_*      = assoc -> selected MC, ranked by weight
//   *_reco_best_*          = nMeasurements / nHoles / chi2, no truth
//   station                = from compact layer name/id (not a hardcoded layer/2)
//   matchesPrimarySelector = every MC matching primary_pdg/status (alias: isPrimary)
//   isSelPrimary           = selected highest-p primary only
//   rec_* / raw_*          = B0TrackerRecHits / B0TrackerRawHits
//   firstHit/lastHit       = aliases of Entry/Exit (min/max-time SimHits)
//   Output tree is in -Phistsfile (default eicrecon.root), directory B0Trackers.

extern "C" {
    void InitPlugin(JApplication* app) {
        InitJANAPlugin(app);
        app->Add(new B0Trackers);
    }
}

namespace {

// Collection access is deliberately split into required and optional.
//
// A blanket catch(...) makes six different situations indistinguishable in the
// ntuple -- genuinely empty, factory not registered, collection renamed,
// upstream factory threw, plugin/EICrecon mismatch, real reconstruction bug --
// and writes a zero count for all of them. That is exactly how "B0 tracks are
// zero" hides its own cause. Required products therefore fail loudly, and
// optional products are reported through has_* branches rather than by
// conflating absence with emptiness.

// Required: any failure propagates and aborts the job.
template <typename T>
std::vector<const T*> getRequired(const std::shared_ptr<const JEvent>& event, const char* name) {
    return event->Get<T>(name);
}

// Optional: returns false only when no factory is registered for (T, name).
// An exception thrown *by* an existing factory still propagates -- that is a
// reconstruction failure, not an absent collection.
template <typename T>
bool getOpt(const std::shared_ptr<const JEvent>& event, const char* name,
            std::vector<const T*>& out) {
    out.clear();
    if (event->GetFactory<T>(name, false) == nullptr) {
        return false;
    }
    out = event->Get<T>(name);
    return true;
}

// ACTS 45 replaced the TrackStateFlag enumerators with named accessors, so a
// direct flags.test(Acts::TrackStateFlag::...) call stops compiling there.
//
// The mask below is built only from accessors that are pure single-bit tests.
// isMeasurement() and isMaterial() are compound predicates in ACTS 45
// (measurement-and-not-outlier, material-and-not-measurement) and using them
// would silently change what this diagnostic column means between versions.
// The >= 45 branch follows the documented ACTS 45 API; the image in use here
// is ACTS 44.4.0, so only the < 45 branch is exercised by our own builds.
enum StateTypeBit {
    kStateMeasurement   = 1 << 0,
    kStateParameter     = 1 << 1,
    kStateOutlier       = 1 << 2,
    kStateHole          = 1 << 3,
    kStateMaterial      = 1 << 4,
    kStateSharedHit     = 1 << 5,
    kStateSplitHit      = 1 << 6,
    kStateNoExpectedHit = 1 << 7,
};

template <typename FlagsT>
int trackStateTypeMask(const FlagsT& flags) {
#if Acts_VERSION_MAJOR >= 45
    return (flags.hasMeasurement()   ? kStateMeasurement   : 0) |
           (flags.hasParameters()    ? kStateParameter     : 0) |
           (flags.isOutlier()        ? kStateOutlier       : 0) |
           (flags.isHole()           ? kStateHole          : 0) |
           (flags.hasMaterial()      ? kStateMaterial      : 0) |
           (flags.isSharedHit()      ? kStateSharedHit     : 0) |
           (flags.isSplitHit()       ? kStateSplitHit      : 0) |
           (flags.hasNoExpectedHit() ? kStateNoExpectedHit : 0);
#else
    return (flags.test(Acts::TrackStateFlag::MeasurementFlag)   ? kStateMeasurement   : 0) |
           (flags.test(Acts::TrackStateFlag::ParameterFlag)     ? kStateParameter     : 0) |
           (flags.test(Acts::TrackStateFlag::OutlierFlag)       ? kStateOutlier       : 0) |
           (flags.test(Acts::TrackStateFlag::HoleFlag)          ? kStateHole          : 0) |
           (flags.test(Acts::TrackStateFlag::MaterialFlag)      ? kStateMaterial      : 0) |
           (flags.test(Acts::TrackStateFlag::SharedHitFlag)     ? kStateSharedHit     : 0) |
           (flags.test(Acts::TrackStateFlag::SplitHitFlag)      ? kStateSplitHit      : 0) |
           (flags.test(Acts::TrackStateFlag::NoExpectedHitFlag) ? kStateNoExpectedHit : 0);
#endif
}

bool decoderHasField(const dd4hep::DDSegmentation::BitFieldCoder* decoder, const char* field) {
    if (decoder == nullptr) {
        return false;
    }
    try {
        dd4hep::CellID tmp = 0;
        decoder->set(tmp, field, 0);
        return true;
    } catch (...) {
        return false;
    }
}

} // namespace

int B0Trackers::stationOf(int layer) const {
    const auto it = m_layerToStation.find(layer);
    if (it != m_layerToStation.end()) {
        return it->second;
    }
    return b0trk::stationFromLayerId(layer, m_perPlaneFrontBack);
}

void B0Trackers::Init() {
    auto* app = GetApplication();
    m_log = app->GetService<Log_service>()->logger("B0Trackers");

    app->SetDefaultParameter("B0Trackers:primary_pdg", m_primaryPdg,
                             "PDG code used to tag the selected primary particle");
    app->SetDefaultParameter("B0Trackers:primary_status", m_primaryStatus,
                             "Generator status used to tag the selected primary particle");
    app->SetDefaultParameter("B0Trackers:min_measurement_stations", m_minMeasurementStations,
                             "Minimum distinct selected-primary measurement stations for reconstructability");
    app->SetDefaultParameter("B0Trackers:fallback_max_normal_mm", m_fallbackMaxNormalMm,
                             "Max |localZ| (mm) allowed for nearest-sensor fallback");
    app->SetDefaultParameter("B0Trackers:fail_on_empty_sensor_map", m_failOnEmptySensorMap,
                             "Throw in Init if the B0 sensor map is empty");
    app->SetDefaultParameter("B0Trackers:fail_on_incomplete_surface_map",
                             m_failOnIncompleteSurfaceMap,
                             "Throw in Init unless every B0 sensor maps to an ACTS surface");
    app->SetDefaultParameter("B0Trackers:enable_truth_seeded_chain", m_enableTruthSeededChain,
                             "Request and analyze the truth-seeded B0 CKF chain");
    app->SetDefaultParameter("B0Trackers:enable_stub_seeded_chain", m_enableStubSeededChain,
                             "Request and analyze the stub-seeded B0 CKF chain");
    app->SetDefaultParameter("B0Trackers:write_track_states", m_writeTrackStates,
                             "Request ACTS track containers and write per-state diagnostics");

    if (const char* cfg = std::getenv("DETECTOR_CONFIG")) {
        m_geometryName = cfg;
    }
    if (const char* path = std::getenv("DETECTOR_PATH")) {
        m_detectorPath = path;
    }

    // RAII: several geometry and service calls below can throw, and a manual
    // release on every path is one edit away from leaking the global lock and
    // deadlocking the job.
    auto rootLock = app->GetService<JGlobalRootLock>();
    struct RootWriteLock {
        std::shared_ptr<JGlobalRootLock> lock;
        TDirectory* prevDir = nullptr;
        explicit RootWriteLock(std::shared_ptr<JGlobalRootLock> l) : lock(std::move(l)) {
            lock->acquire_write_lock();
            prevDir = gDirectory;
        }
        ~RootWriteLock() {
            if (prevDir != nullptr) {
                prevDir->cd();
            }
            lock->release_lock();
        }
        RootWriteLock(const RootWriteLock&) = delete;
        RootWriteLock& operator=(const RootWriteLock&) = delete;
    } rootWriteLock{rootLock};
    auto rf_svc = app->GetService<RootFile_service>();
    TFile* outfile = rf_svc->GetHistFile()->GetFile();
    TDirectory* dir = outfile->GetDirectory("B0Trackers");
    if (dir == nullptr) {
        dir = outfile->mkdir("B0Trackers");
    }
    if (dir == nullptr) {
        throw JException("B0Trackers: failed to create output directory");
    }
    dir->cd();

    m_tree = new TTree("hits", "B0 truth, hits, and tracks");

    m_tree->Branch("schema_version", &m_schemaVersion);
    m_tree->Branch("geometry_name", &m_geometryName);
    m_tree->Branch("detector_path", &m_detectorPath);
    m_tree->Branch("eventNumber", &m_eventNumber);
    m_tree->Branch("config_enable_truth_seeded_chain", &m_enableTruthSeededChain);
    m_tree->Branch("config_enable_stub_seeded_chain", &m_enableStubSeededChain);
    m_tree->Branch("config_write_track_states", &m_writeTrackStates);

    m_tree->Branch("xR",        &vm_xR);
    m_tree->Branch("yR",        &vm_yR);
    m_tree->Branch("zR",        &vm_zR);
    m_tree->Branch("xT",        &vm_xT);
    m_tree->Branch("yT",        &vm_yT);
    m_tree->Branch("zT",        &vm_zT);
    m_tree->Branch("aclgad_xPixT", &vm_aclgad_xPixT);
    m_tree->Branch("aclgad_yPixT", &vm_aclgad_yPixT);
    m_tree->Branch("aclgad_zPixT", &vm_aclgad_zPixT);
    m_tree->Branch("aclgad_dxT", &vm_aclgad_dxT);
    m_tree->Branch("aclgad_dyT", &vm_aclgad_dyT);
    m_tree->Branch("aclgad_dzT", &vm_aclgad_dzT);
    m_tree->Branch("aclgad_xPixR", &vm_aclgad_xPixR);
    m_tree->Branch("aclgad_yPixR", &vm_aclgad_yPixR);
    m_tree->Branch("aclgad_zPixR", &vm_aclgad_zPixR);
    m_tree->Branch("aclgad_dxR", &vm_aclgad_dxR);
    m_tree->Branch("aclgad_dyR", &vm_aclgad_dyR);
    m_tree->Branch("aclgad_dzR", &vm_aclgad_dzR);
    m_tree->Branch("plane",     &vm_plane);
    m_tree->Branch("station",   &vm_station);
    m_tree->Branch("module",    &vm_module);
    m_tree->Branch("side",      &vm_side);
    m_tree->Branch("sensor",    &vm_sensor);
    m_tree->Branch("pixX",      &vm_pixX);
    m_tree->Branch("pixY",      &vm_pixY);
    m_tree->Branch("pixZ",      &vm_pixZ);
    m_tree->Branch("aclgad_pixXT", &vm_aclgad_pixXT);
    m_tree->Branch("aclgad_pixYT", &vm_aclgad_pixYT);
    m_tree->Branch("aclgad_pixXR", &vm_aclgad_pixXR);
    m_tree->Branch("aclgad_pixYR", &vm_aclgad_pixYR);
    m_tree->Branch("cellID",    &vm_cellID);
    m_tree->Branch("detX",      &vm_detX);
    m_tree->Branch("detY",      &vm_detY);
    m_tree->Branch("detZ",      &vm_detZ);
    m_tree->Branch("eDep",      &vm_eDep);
    m_tree->Branch("time",      &vm_time);
    m_tree->Branch("path",      &vm_path);
    m_tree->Branch("pdg",       &vm_pdg);
    m_tree->Branch("mcIndex",   &vm_mcIndex);
    m_tree->Branch("mcCollectionID", &vm_mcCollectionID);
    m_tree->Branch("px",        &vm_px);
    m_tree->Branch("py",        &vm_py);
    m_tree->Branch("pz",        &vm_pz);
    m_tree->Branch("p",         &vm_p);
    m_tree->Branch("pT",        &vm_pT);
    m_tree->Branch("status",    &vm_status);
    m_tree->Branch("isPrimary", &vm_isPrimary);
    m_tree->Branch("matchesPrimarySelector", &vm_isPrimary);
    m_tree->Branch("isSelPrimary", &vm_isSelPrimary);
    m_tree->Branch("cell_fired", &vm_cellFired);
    m_tree->Branch("xP",        &vm_xP);
    m_tree->Branch("yP",        &vm_yP);
    m_tree->Branch("zP",        &vm_zP);
    m_tree->Branch("pathP",     &vm_pathP);
    m_tree->Branch("timeP",     &vm_timeP);
    m_tree->Branch("planeP",    &vm_planeP);
    m_tree->Branch("stationP",  &vm_stationP);
    m_tree->Branch("moduleP",   &vm_moduleP);
    m_tree->Branch("sideP",     &vm_sideP);
    m_tree->Branch("sensorP",   &vm_sensorP);
    m_tree->Branch("pdgP",      &vm_pdgP);
    m_tree->Branch("statusP",   &vm_statusP);
    m_tree->Branch("isPrimaryP", &vm_isPrimaryP);
    m_tree->Branch("matchesPrimarySelectorP", &vm_isPrimaryP);
    m_tree->Branch("isSelPrimaryP", &vm_isSelPrimaryP);
    m_tree->Branch("mcIndexP",  &vm_mcIndexP);
    m_tree->Branch("mcCollectionIDP", &vm_mcCollectionIDP);
    m_tree->Branch("pxP",       &vm_pxP);
    m_tree->Branch("pyP",       &vm_pyP);
    m_tree->Branch("pzP",       &vm_pzP);
    m_tree->Branch("pP",        &vm_pP);
    m_tree->Branch("pTP",       &vm_pTP);
    m_tree->Branch("beampx",    &beam_px);
    m_tree->Branch("beampy",    &beam_py);
    m_tree->Branch("beampz",    &beam_pz);
    m_tree->Branch("beamp",     &beam_p);
    m_tree->Branch("beampT",    &beam_pT);
    m_tree->Branch("beam_pdg",  &beam_pdg);
    m_tree->Branch("primary_pdg", &m_primaryPdgOut);
    m_tree->Branch("primary_status", &m_primaryStatusOut);
    m_tree->Branch("primary_mcIndex", &m_primaryMcIndex);
    m_tree->Branch("primary_mcCollectionID", &m_primaryMcCollectionID);
    m_tree->Branch("primary_px", &m_primaryPx);
    m_tree->Branch("primary_py", &m_primaryPy);
    m_tree->Branch("primary_pz", &m_primaryPz);
    m_tree->Branch("primary_p",  &m_primaryP);
    m_tree->Branch("primary_pT", &m_primaryPT);
    m_tree->Branch("genPpx",    &m_genPpx);
    m_tree->Branch("genPpy",    &m_genPpy);
    m_tree->Branch("genPpz",    &m_genPpz);
    m_tree->Branch("genPp",     &m_genPp);
    m_tree->Branch("genPpT",    &m_genPpT);
    m_tree->Branch("genBeamP",  &m_genBeamP);
    m_tree->Branch("genBeamPT", &m_genBeamPT);
    m_tree->Branch("genBeamPx", &m_genBeamPx);
    m_tree->Branch("genBeamPy", &m_genBeamPy);
    m_tree->Branch("genBeamPz", &m_genBeamPz);
    m_tree->Branch("genBeamPP", &m_genBeamPP);
    m_tree->Branch("genBeamPPT",&m_genBeamPPT);
    m_tree->Branch("genBeamPPx",&m_genBeamPPx);
    m_tree->Branch("genBeamPPy",&m_genBeamPPy);
    m_tree->Branch("genBeamPPz",&m_genBeamPPz);
    m_tree->Branch("beam_proton_px", &m_genPpx);
    m_tree->Branch("beam_proton_py", &m_genPpy);
    m_tree->Branch("beam_proton_pz", &m_genPpz);
    m_tree->Branch("beam_proton_p",  &m_genPp);
    m_tree->Branch("beam_proton_pT", &m_genPpT);
    m_tree->Branch("scattered_e_px", &m_genBeamPx);
    m_tree->Branch("scattered_e_py", &m_genBeamPy);
    m_tree->Branch("scattered_e_pz", &m_genBeamPz);
    m_tree->Branch("scattered_e_p",  &m_genBeamP);
    m_tree->Branch("scattered_e_pT", &m_genBeamPT);
    m_tree->Branch("scattered_p_px", &m_genBeamPPx);
    m_tree->Branch("scattered_p_py", &m_genBeamPPy);
    m_tree->Branch("scattered_p_pz", &m_genBeamPPz);
    m_tree->Branch("scattered_p_p",  &m_genBeamPP);
    m_tree->Branch("scattered_p_pT", &m_genBeamPPT);
    m_tree->Branch("sel_primary_mcIndex", &m_selPrimaryMcIndex);
    m_tree->Branch("sel_primary_mcCollectionID", &m_selPrimaryMcCollectionID);
    m_tree->Branch("sel_primary_px", &m_selPrimaryPx);
    m_tree->Branch("sel_primary_py", &m_selPrimaryPy);
    m_tree->Branch("sel_primary_pz", &m_selPrimaryPz);
    m_tree->Branch("sel_primary_p",  &m_selPrimaryP);
    m_tree->Branch("sel_primary_pT", &m_selPrimaryPT);
    m_tree->Branch("sel_primary_thscat_mrad", &m_selPrimaryThscatMrad);
    m_tree->Branch("sel_primary_charge", &m_selPrimaryCharge);
    m_tree->Branch("n_stations_primary", &m_nStationsPrimary);
    m_tree->Branch("min_measurement_stations_required", &m_minMeasurementStations);
    m_tree->Branch("n_measurements_selected_primary", &m_nSelectedPrimaryMeasurements);
    m_tree->Branch("n_measurement_stations_selected_primary",
                   &m_nMeasurementStationsSelectedPrimary);
    m_tree->Branch("sel_primary_measurement_reconstructable",
                   &m_selPrimaryMeasurementReconstructable);
    m_tree->Branch("sel_primary_has_seed", &m_selPrimaryHasSeed);
    m_tree->Branch("sel_primary_has_unfiltered_track", &m_selPrimaryHasUnfilteredTrack);
    m_tree->Branch("sel_primary_has_filtered_track", &m_selPrimaryHasFilteredTrack);
    m_tree->Branch("sel_primary_has_truth_matched_track", &m_selPrimaryHasTruthMatchedTrack);

    bindTrackChain("trk_", m_ts);
    bindTrackChain("ckf_trk_", m_ckf);

    m_tree->Branch("xEntry",     &vm_xEntry);
    m_tree->Branch("yEntry",     &vm_yEntry);
    m_tree->Branch("zEntry",     &vm_zEntry);
    m_tree->Branch("timeEntry",  &vm_timeEntry);
    m_tree->Branch("pxEntry",    &vm_pxEntry);
    m_tree->Branch("pyEntry",    &vm_pyEntry);
    m_tree->Branch("pzEntry",    &vm_pzEntry);
    m_tree->Branch("pEntry",     &vm_pEntry);
    m_tree->Branch("pTEntry",    &vm_pTEntry);
    m_tree->Branch("xExit",      &vm_xExit);
    m_tree->Branch("yExit",      &vm_yExit);
    m_tree->Branch("zExit",      &vm_zExit);
    m_tree->Branch("timeExit",   &vm_timeExit);
    m_tree->Branch("pxExit",     &vm_pxExit);
    m_tree->Branch("pyExit",     &vm_pyExit);
    m_tree->Branch("pzExit",     &vm_pzExit);
    m_tree->Branch("pExit",      &vm_pExit);
    m_tree->Branch("pTExit",     &vm_pTExit);
    m_tree->Branch("xFirstHit",  &vm_xEntry);
    m_tree->Branch("yFirstHit",  &vm_yEntry);
    m_tree->Branch("zFirstHit",  &vm_zEntry);
    m_tree->Branch("timeFirstHit", &vm_timeEntry);
    m_tree->Branch("xLastHit",   &vm_xExit);
    m_tree->Branch("yLastHit",   &vm_yExit);
    m_tree->Branch("zLastHit",   &vm_zExit);
    m_tree->Branch("timeLastHit", &vm_timeExit);
    m_tree->Branch("sensorEntry", &vm_sensorEntry);
    m_tree->Branch("sensorExit",  &vm_sensorExit);
    m_tree->Branch("cellIDEntry", &vm_cellIDEntry);
    m_tree->Branch("cellIDExit",  &vm_cellIDExit);
    m_tree->Branch("planeEE",    &vm_planeEE);
    m_tree->Branch("stationEE",  &vm_stationEE);
    m_tree->Branch("moduleEE",   &vm_moduleEE);
    m_tree->Branch("sideEE",     &vm_sideEE);
    m_tree->Branch("pdgEE",      &vm_pdgEE);
    m_tree->Branch("mcIndexEE",  &vm_mcIndexEE);
    m_tree->Branch("mcCollectionIDEE", &vm_mcCollectionIDEE);
    m_tree->Branch("nStepsEE",   &vm_nStepsEE);
    m_tree->Branch("nSimHitsEE", &vm_nStepsEE);
    m_tree->Branch("statusEE",   &vm_statusEE);
    m_tree->Branch("isPrimaryEE", &vm_isPrimaryEE);
    m_tree->Branch("matchesPrimarySelectorEE", &vm_isPrimaryEE);
    m_tree->Branch("isSelPrimaryEE", &vm_isSelPrimaryEE);

    m_tree->Branch("raw_cellID", &vm_raw_cellID);
    m_tree->Branch("raw_charge", &vm_raw_charge);
    m_tree->Branch("raw_timeStamp", &vm_raw_timeStamp);
    m_tree->Branch("raw_plane", &vm_raw_plane);
    m_tree->Branch("raw_station", &vm_raw_station);
    m_tree->Branch("raw_module", &vm_raw_module);
    m_tree->Branch("raw_side", &vm_raw_side);
    m_tree->Branch("raw_sensor", &vm_raw_sensor);
    m_tree->Branch("raw_mcIndex", &vm_raw_mcIndex);
    m_tree->Branch("raw_mcCollectionID", &vm_raw_mcCollectionID);
    m_tree->Branch("raw_nContribSim", &vm_raw_nContribSim);
    m_tree->Branch("raw_nContribMc", &vm_raw_nContribMc);
    m_tree->Branch("raw_dominantFrac", &vm_raw_dominantFrac);
    m_tree->Branch("raw_totalEdep", &vm_raw_totalEdep);
    m_tree->Branch("raw_mixedCell", &vm_raw_mixedCell);

    m_tree->Branch("rec_x", &vm_rec_x);
    m_tree->Branch("rec_y", &vm_rec_y);
    m_tree->Branch("rec_z", &vm_rec_z);
    m_tree->Branch("rec_covxx", &vm_rec_covxx);
    m_tree->Branch("rec_covyy", &vm_rec_covyy);
    m_tree->Branch("rec_covzz", &vm_rec_covzz);
    m_tree->Branch("rec_time", &vm_rec_time);
    m_tree->Branch("rec_time_err", &vm_rec_time_err);
    m_tree->Branch("rec_edep", &vm_rec_edep);
    m_tree->Branch("rec_edep_err", &vm_rec_edep_err);
    m_tree->Branch("rec_cellID", &vm_rec_cellID);
    m_tree->Branch("rec_plane", &vm_rec_plane);
    m_tree->Branch("rec_station", &vm_rec_station);
    m_tree->Branch("rec_module", &vm_rec_module);
    m_tree->Branch("rec_side", &vm_rec_side);
    m_tree->Branch("rec_sensor", &vm_rec_sensor);
    m_tree->Branch("rec_pixX", &vm_rec_pixX);
    m_tree->Branch("rec_pixY", &vm_rec_pixY);
    m_tree->Branch("rec_pixZ", &vm_rec_pixZ);
    m_tree->Branch("rec_mcIndex", &vm_rec_mcIndex);
    m_tree->Branch("rec_mcCollectionID", &vm_rec_mcCollectionID);
    m_tree->Branch("rec_nContribSim", &vm_rec_nContribSim);
    m_tree->Branch("rec_nContribMc", &vm_rec_nContribMc);
    m_tree->Branch("rec_dominantFrac", &vm_rec_dominantFrac);
    m_tree->Branch("rec_totalEdep", &vm_rec_totalEdep);
    m_tree->Branch("rec_mixedCell", &vm_rec_mixedCell);

    m_tree->Branch("seed_quality", &vm_seed_quality);
    m_tree->Branch("seed_p", &vm_seed_p);
    m_tree->Branch("seed_qOverP", &vm_seed_qOverP);
    m_tree->Branch("seed_theta", &vm_seed_theta);
    m_tree->Branch("seed_phi", &vm_seed_phi);
    m_tree->Branch("seed_loc0", &vm_seed_loc0);
    m_tree->Branch("seed_loc1", &vm_seed_loc1);
    m_tree->Branch("seed_sigma_qOverP", &vm_seed_sigma_qOverP);
    m_tree->Branch("seed_sigma_theta", &vm_seed_sigma_theta);
    m_tree->Branch("seed_sigma_phi", &vm_seed_sigma_phi);
    m_tree->Branch("seed_nHits", &vm_seed_nHits);
    m_tree->Branch("seed_charge", &vm_seed_charge);
    m_tree->Branch("seed_momentum_resolved", &vm_seed_momentum_resolved);
    m_tree->Branch("seed_became_track", &vm_seed_became_track);
    m_tree->Branch("seed_made_unfiltered_track", &vm_seed_made_unfiltered_track);
    m_tree->Branch("seed_survived_ambiguity", &vm_seed_survived_ambiguity);
    m_tree->Branch("seed_n_unfiltered_tracks", &vm_seed_n_unfiltered_tracks);
    m_tree->Branch("seed_n_filtered_tracks", &vm_seed_n_filtered_tracks);
    m_tree->Branch("seed_assoc_mcIndex", &vm_seed_assoc_mcIndex);
    m_tree->Branch("seed_assoc_mcCollectionID", &vm_seed_assoc_mcCollectionID);
    m_tree->Branch("seed_assoc_weight", &vm_seed_assoc_weight);
    m_tree->Branch("truth_seed_quality", &vm_truth_seed_quality);
    m_tree->Branch("truth_seed_p", &vm_truth_seed_p);
    m_tree->Branch("truth_seed_qOverP", &vm_truth_seed_qOverP);
    m_tree->Branch("truth_seed_theta", &vm_truth_seed_theta);
    m_tree->Branch("truth_seed_phi", &vm_truth_seed_phi);
    m_tree->Branch("truth_seed_loc0", &vm_truth_seed_loc0);
    m_tree->Branch("truth_seed_loc1", &vm_truth_seed_loc1);
    m_tree->Branch("truth_seed_sigma_qOverP", &vm_truth_seed_sigma_qOverP);
    m_tree->Branch("truth_seed_sigma_theta", &vm_truth_seed_sigma_theta);
    m_tree->Branch("truth_seed_sigma_phi", &vm_truth_seed_sigma_phi);
    m_tree->Branch("truth_seed_nHits", &vm_truth_seed_nHits);
    m_tree->Branch("truth_seed_charge", &vm_truth_seed_charge);
    m_tree->Branch("truth_seed_momentum_resolved", &vm_truth_seed_momentum_resolved);
    m_tree->Branch("truth_seed_became_track", &vm_truth_seed_became_track);
    m_tree->Branch("truth_seed_made_unfiltered_track", &vm_truth_seed_made_unfiltered_track);
    m_tree->Branch("truth_seed_survived_ambiguity", &vm_truth_seed_survived_ambiguity);
    m_tree->Branch("truth_seed_n_unfiltered_tracks", &vm_truth_seed_n_unfiltered_tracks);
    m_tree->Branch("truth_seed_n_filtered_tracks", &vm_truth_seed_n_filtered_tracks);

    m_tree->Branch("n_simhits", &m_nSimHits);
    m_tree->Branch("n_rawhits", &m_nRawHits);
    m_tree->Branch("n_rechits", &m_nRecHits);
    m_tree->Branch("n_measurements", &m_nMeasurements);
    m_tree->Branch("n_truth_seeds", &m_nTruthSeeds);
    m_tree->Branch("n_stub_seeds", &m_nStubSeeds);
    m_tree->Branch("n_ts_unfiltered", &m_nTsUnfiltered);
    m_tree->Branch("n_ts_filtered", &m_nTsFiltered);
    m_tree->Branch("n_ckf_unfiltered", &m_nCkfUnfiltered);
    m_tree->Branch("n_ckf_filtered", &m_nCkfFiltered);
    m_tree->Branch("n_simhits_unresolved_cellid", &m_nSimHitsUnresolvedCellID);
    m_tree->Branch("n_missing_mc_relation", &m_nMissingMcRelation);
    m_tree->Branch("n_sensor_map_exact", &m_nSensorMapExact);
    m_tree->Branch("n_sensor_map_fallback", &m_nSensorMapFallback);
    m_tree->Branch("n_sensor_map_failed", &m_nSensorMapFailed);
    m_tree->Branch("n_pixel_snap_failed", &m_nPixelSnapFailed);

    // Input availability. A false here means the factory is not registered at
    // all -- distinct from a registered factory that produced nothing.
    m_tree->Branch("has_raw_assocs", &m_hasRawAssocs);
    m_tree->Branch("has_stub_seeds", &m_hasStubSeeds);
    m_tree->Branch("has_truth_seeds", &m_hasTruthSeeds);
    m_tree->Branch("has_ts_track_params", &m_hasTsTrackParams);
    m_tree->Branch("has_ts_trajectories", &m_hasTsTrajectories);
    m_tree->Branch("has_ts_trajectories_unfiltered", &m_hasTsTrajectoriesUnfiltered);
    m_tree->Branch("has_ts_tracks", &m_hasTsTracks);
    m_tree->Branch("has_ts_assocs", &m_hasTsAssocs);
    m_tree->Branch("has_ts_acts_states", &m_hasTsActsStates);
    m_tree->Branch("has_ts_acts_tracks", &m_hasTsActsTracks);
    m_tree->Branch("has_ts_tracks_unfiltered", &m_hasTsTracksUnfiltered);
    m_tree->Branch("has_ckf_track_params", &m_hasCkfTrackParams);
    m_tree->Branch("has_ckf_trajectories", &m_hasCkfTrajectories);
    m_tree->Branch("has_ckf_trajectories_unfiltered", &m_hasCkfTrajectoriesUnfiltered);
    m_tree->Branch("has_ckf_tracks", &m_hasCkfTracks);
    m_tree->Branch("has_ckf_assocs", &m_hasCkfAssocs);
    m_tree->Branch("has_ckf_acts_states", &m_hasCkfActsStates);
    m_tree->Branch("has_ckf_acts_tracks", &m_hasCkfActsTracks);
    m_tree->Branch("has_ckf_tracks_unfiltered", &m_hasCkfTracksUnfiltered);
    m_tree->Branch("has_ckf_assocs_unfiltered", &m_hasCkfAssocsUnfiltered);

    m_geoSvc = app->GetService<DD4hep_service>();
    m_actsGeoSvc = app->GetService<ACTSGeo_service>();
    m_actsGeoProvider = m_actsGeoSvc ? m_actsGeoSvc->actsGeoProvider() : nullptr;
    auto ro  = m_geoSvc->detector()->readout("B0TrackerHits");
    m_segmentation = ro.segmentation();
    m_decoder = ro.idSpec().decoder();
    m_volman  = m_geoSvc->detector()->volumeManager();
    m_hasPixX = decoderHasField(m_decoder, "x");
    m_hasPixY = decoderHasField(m_decoder, "y");
    m_hasPixZ = decoderHasField(m_decoder, "z");

    const auto sensorHalfExtentMm = [this](std::initializer_list<const char*> names,
                                           double fallbackMm) -> double {
        for (const char* name : names) {
            try {
                return 5.0 * m_geoSvc->detector()->constant<double>(name);
            } catch (const std::exception&) {
            }
        }
        return fallbackMm;
    };
    m_sensorHalfX = sensorHalfExtentMm({"B0TrackerSensorWidth", "SensorWidth"}, 8.0);
    m_sensorHalfY = sensorHalfExtentMm({"B0TrackerSensorLength", "SensorLength"}, 8.0);

    auto b0Det = m_geoSvc->detector()->detector("B0Tracker");
    if (!b0Det.isValid()) {
        throw JException("B0Trackers: detector B0Tracker is not valid");
    }

    const auto findVolID = [](const dd4hep::PlacedVolume& pv, const std::string& name) -> int {
        for (const auto& id : pv.volIDs()) {
            if (id.first == name) return id.second;
        }
        return -1;
    };
    const auto packVolumeID = [this, &b0Det](int layer_id, int module_id, int sensor_id,
                                             bool& ok) -> std::uint64_t {
        ok = false;
        dd4hep::CellID refCellID = 0;
        try {
            m_decoder->set(refCellID, "system", b0Det.id());
            m_decoder->set(refCellID, "layer", layer_id);
            m_decoder->set(refCellID, "module", module_id);
            m_decoder->set(refCellID, "sensor", sensor_id);
        } catch (const std::exception&) {
            return 0;
        }
        for (const char* field : {"x", "y", "z"}) {
            if ((field[0] == 'x' && !m_hasPixX) || (field[0] == 'y' && !m_hasPixY) ||
                (field[0] == 'z' && !m_hasPixZ)) {
                continue;
            }
            try {
                m_decoder->set(refCellID, field, 0);
            } catch (const std::exception&) {
            }
        }
        ok = true;
        return static_cast<std::uint64_t>(refCellID);
    };

    for (const auto& [layerName, layerDE] : b0Det.children()) {
        const auto layer_pv = layerDE.placement();
        if (!layer_pv.isValid()) continue;
        const int layer_id = findVolID(layer_pv, "layer");
        if (layer_id < 0) continue;

        const auto [nameStation, nameSide] = b0trk::parseLayerName(layerName);
        if (nameSide >= 0) {
            m_perPlaneFrontBack = true;
        }
        const int station = (nameStation > 0)
            ? nameStation
            : b0trk::stationFromLayerId(layer_id, nameSide >= 0);
        m_layerToStation[layer_id] = station;

        int layerSide = nameSide;
        for (const auto& [modName, modDE] : layerDE.children()) {
            const auto pv = modDE.placement();
            if (!pv.isValid()) continue;
            const int module_id = findVolID(pv, "module");
            if (module_id < 0) continue;
            int side = layerSide;
            if (side < 0) {
                const TGeoMatrix& mat = pv.matrix();
                const double* tr = mat.GetTranslation();
                side = (tr[2] > 0.0) ? 1 : 0;
            }
            m_moduleToSide[{layer_id, module_id}] = side;

            for (const auto& [sensorName, sensorDE] : modDE.children()) {
                const auto sensorPV = sensorDE.placement();
                if (!sensorPV.isValid()) continue;
                const int sensor_id = findVolID(sensorPV, "sensor");
                if (sensor_id < 0) continue;
                bool packed = false;
                const auto refCellID = packVolumeID(layer_id, module_id, sensor_id, packed);
                if (!packed) continue;
                m_sensorRefs.push_back({refCellID, layer_id, module_id, side, sensor_id, station, sensorDE});
            }
        }
    }

    if (m_actsGeoProvider) {
        std::unordered_map<std::uint64_t, std::size_t> cellToIdx;
        for (std::size_t i = 0; i < m_sensorRefs.size(); ++i) {
            cellToIdx.emplace(m_sensorRefs[i].cellID, i);
        }
        for (const auto& [volumeID, surface] : m_actsGeoProvider->surfaceMap()) {
            if (surface == nullptr) continue;
            const auto it = cellToIdx.find(volumeID);
            if (it == cellToIdx.end()) continue;
            m_surfaceToSensorIdx.emplace(surface->geometryId().value(), it->second);
        }
    }

    m_log->info("schema {}  geometry='{}'  sensors={}  surfaceMap={}  "
                "pixFields x={} y={} z={}  perPlaneFrontBack={}",
                m_schemaVersion, m_geometryName, m_sensorRefs.size(), m_surfaceToSensorIdx.size(),
                m_hasPixX, m_hasPixY, m_hasPixZ, m_perPlaneFrontBack);

    if (m_sensorRefs.empty() && m_failOnEmptySensorMap) {
        throw JException("B0Trackers: empty sensor map (check B0TrackerHits id fields / geometry)");
    }
    // An incomplete ACTS->DD4hep surface map is silently absorbed by the
    // nearest-sensor fallback, which is a debugging facility rather than a
    // physics-grade mapping. Production runs should require exact coverage.
    if (m_failOnIncompleteSurfaceMap && !m_sensorRefs.empty() &&
        m_surfaceToSensorIdx.size() != m_sensorRefs.size()) {
        throw JException(
            "B0Trackers: incomplete ACTS surface map (" +
            std::to_string(m_surfaceToSensorIdx.size()) + " of " +
            std::to_string(m_sensorRefs.size()) +
            " sensors mapped exactly); set B0Trackers:fail_on_incomplete_surface_map=0 "
            "to fall back to nearest-sensor matching");
    }
}

void B0Trackers::Process(const std::shared_ptr<const JEvent>& event) {
    // Required -- a missing or throwing factory here aborts the job.
    auto mcparticles  = getRequired<edm4hep::MCParticle>(event, "MCParticles");
    auto simHits      = getRequired<edm4hep::SimTrackerHit>(event, "B0TrackerHits");
    auto rawHits      = getRequired<edm4eic::RawTrackerHit>(event, "B0TrackerRawHits");
    auto recHits      = getRequired<edm4eic::TrackerHit>(event, "B0TrackerRecHits");
    auto measurements = getRequired<edm4eic::Measurement2D>(event, "B0TrackerMeasurements");

    // Optional -- a configuration may legitimately not run these stages.
    std::vector<const edm4eic::MCRecoTrackerHitAssociation*> rawAssocs;
    std::vector<const edm4eic::TrackSeed*> stubSeeds;
    std::vector<const edm4eic::TrackSeed*> truthSeeds;
    const bool hasRawAssocs =
        getOpt(event, "B0TrackerRawHitAssociations", rawAssocs);
    bool hasStubSeeds = false;
    bool hasTruthSeeds = false;
    if (m_enableStubSeededChain) hasStubSeeds = getOpt(event, "B0TrackerSeeds", stubSeeds);
    if (m_enableTruthSeededChain) hasTruthSeeds = getOpt(event, "B0TrackerTruthSeeds", truthSeeds);

    std::vector<const edm4eic::TrackParameters*> tsTracks;
    std::vector<const edm4eic::Trajectory*> tsTrajectories;
    std::vector<const edm4eic::Trajectory*> tsTrajectoriesUnfiltered;
    std::vector<const edm4eic::Track*> tsEdmTracks;
    std::vector<const edm4eic::MCRecoTrackParticleAssociation*> tsAssocs;
    std::vector<const Acts::ConstVectorMultiTrajectory*> tsActsTrackStates;
    std::vector<const Acts::ConstVectorTrackContainer*> tsActsTracks;
    std::vector<const edm4eic::Track*> tsUnfiltered;
    bool hasTsTrackParams = false;
    bool hasTsTrajectories = false;
    bool hasTsTrajectoriesUnfiltered = false;
    bool hasTsTracks = false;
    bool hasTsAssocs = false;
    bool hasTsActsStates = false;
    bool hasTsActsTracks = false;
    bool hasTsTracksUnfiltered = false;
    if (m_enableTruthSeededChain) {
        hasTsTrackParams = getOpt(event, "B0TrackerCKFTruthSeededTrackParameters", tsTracks);
        hasTsTrajectories = getOpt(event, "B0TrackerCKFTruthSeededTrajectories", tsTrajectories);
        hasTsTrajectoriesUnfiltered = getOpt(
            event, "B0TrackerCKFTruthSeededTrajectoriesUnfiltered", tsTrajectoriesUnfiltered);
        hasTsTracks = getOpt(event, "B0TrackerCKFTruthSeededTracks", tsEdmTracks);
        hasTsAssocs = getOpt(event, "B0TrackerCKFTruthSeededTrackAssociations", tsAssocs);
        hasTsTracksUnfiltered = getOpt(event, "B0TrackerCKFTruthSeededTracksUnfiltered", tsUnfiltered);
        if (m_writeTrackStates) {
            hasTsActsStates = getOpt(event, "B0TrackerCKFTruthSeededActsTrackStates", tsActsTrackStates);
            hasTsActsTracks = getOpt(event, "B0TrackerCKFTruthSeededActsTracks", tsActsTracks);
        }
    }

    std::vector<const edm4eic::TrackParameters*> ckfTracks;
    std::vector<const edm4eic::Trajectory*> ckfTrajectories;
    std::vector<const edm4eic::Trajectory*> ckfTrajectoriesUnfiltered;
    std::vector<const edm4eic::Track*> ckfEdmTracks;
    std::vector<const edm4eic::MCRecoTrackParticleAssociation*> ckfAssocs;
    std::vector<const Acts::ConstVectorMultiTrajectory*> ckfActsTrackStates;
    std::vector<const Acts::ConstVectorTrackContainer*> ckfActsTracks;
    std::vector<const edm4eic::Track*> ckfUnfiltered;
    std::vector<const edm4eic::MCRecoTrackParticleAssociation*> ckfUnfilteredAssocs;
    bool hasCkfTrackParams = false;
    bool hasCkfTrajectories = false;
    bool hasCkfTrajectoriesUnfiltered = false;
    bool hasCkfTracks = false;
    bool hasCkfAssocs = false;
    bool hasCkfActsStates = false;
    bool hasCkfActsTracks = false;
    bool hasCkfTracksUnfiltered = false;
    bool hasCkfAssocsUnfiltered = false;
    if (m_enableStubSeededChain) {
        hasCkfTrackParams = getOpt(event, "B0TrackerCKFTrackParameters", ckfTracks);
        hasCkfTrajectories = getOpt(event, "B0TrackerCKFTrajectories", ckfTrajectories);
        hasCkfTrajectoriesUnfiltered = getOpt(
            event, "B0TrackerCKFTrajectoriesUnfiltered", ckfTrajectoriesUnfiltered);
        hasCkfTracks = getOpt(event, "B0TrackerCKFTracks", ckfEdmTracks);
        hasCkfAssocs = getOpt(event, "B0TrackerCKFTrackAssociations", ckfAssocs);
        hasCkfTracksUnfiltered = getOpt(event, "B0TrackerCKFTracksUnfiltered", ckfUnfiltered);
        hasCkfAssocsUnfiltered = getOpt(
            event, "B0TrackerCKFTrackUnfilteredAssociations", ckfUnfilteredAssocs);
        if (m_writeTrackStates) {
            hasCkfActsStates = getOpt(event, "B0TrackerCKFActsTrackStates", ckfActsTrackStates);
            hasCkfActsTracks = getOpt(event, "B0TrackerCKFActsTracks", ckfActsTracks);
        }
    }

    // Reusable payload owned by this JANA worker thread. All event analysis
    // happens here; only branch-buffer swapping and TTree::Fill() are serialized.
    thread_local EventBuffers local;

    local.m_hasRawAssocs               = hasRawAssocs;
    local.m_hasStubSeeds               = hasStubSeeds;
    local.m_hasTruthSeeds              = hasTruthSeeds;
    local.m_hasTsTrackParams           = hasTsTrackParams;
    local.m_hasTsTrajectories          = hasTsTrajectories;
    local.m_hasTsTrajectoriesUnfiltered = hasTsTrajectoriesUnfiltered;
    local.m_hasTsTracks                = hasTsTracks;
    local.m_hasTsAssocs                = hasTsAssocs;
    local.m_hasTsActsStates            = hasTsActsStates;
    local.m_hasTsActsTracks            = hasTsActsTracks;
    local.m_hasTsTracksUnfiltered      = hasTsTracksUnfiltered;
    local.m_hasCkfTrackParams          = hasCkfTrackParams;
    local.m_hasCkfTrajectories         = hasCkfTrajectories;
    local.m_hasCkfTrajectoriesUnfiltered = hasCkfTrajectoriesUnfiltered;
    local.m_hasCkfTracks               = hasCkfTracks;
    local.m_hasCkfAssocs               = hasCkfAssocs;
    local.m_hasCkfActsStates           = hasCkfActsStates;
    local.m_hasCkfActsTracks           = hasCkfActsTracks;
    local.m_hasCkfTracksUnfiltered     = hasCkfTracksUnfiltered;
    local.m_hasCkfAssocsUnfiltered     = hasCkfAssocsUnfiltered;

    local.m_eventNumber = event->GetEventNumber();

    local.vm_xR.clear();    local.vm_yR.clear();    local.vm_zR.clear();
    local.vm_xT.clear();    local.vm_yT.clear();    local.vm_zT.clear();
    local.vm_aclgad_xPixT.clear(); local.vm_aclgad_yPixT.clear(); local.vm_aclgad_zPixT.clear();
    local.vm_aclgad_dxT.clear(); local.vm_aclgad_dyT.clear(); local.vm_aclgad_dzT.clear();
    local.vm_aclgad_xPixR.clear(); local.vm_aclgad_yPixR.clear(); local.vm_aclgad_zPixR.clear();
    local.vm_aclgad_dxR.clear(); local.vm_aclgad_dyR.clear(); local.vm_aclgad_dzR.clear();
    local.vm_detX.clear();  local.vm_detY.clear();  local.vm_detZ.clear();
    local.vm_plane.clear(); local.vm_station.clear(); local.vm_module.clear(); local.vm_side.clear(); local.vm_sensor.clear();
    local.vm_pixX.clear();  local.vm_pixY.clear();   local.vm_pixZ.clear();
    local.vm_aclgad_pixXT.clear(); local.vm_aclgad_pixYT.clear();
    local.vm_aclgad_pixXR.clear(); local.vm_aclgad_pixYR.clear();
    local.vm_cellID.clear(); local.vm_mcIndex.clear(); local.vm_mcCollectionID.clear();
    local.vm_eDep.clear();  local.vm_time.clear();   local.vm_path.clear();
    local.vm_pdg.clear();   local.vm_status.clear();
    local.vm_isPrimary.clear(); local.vm_isSelPrimary.clear(); local.vm_cellFired.clear();
    local.vm_px.clear();    local.vm_py.clear();    local.vm_pz.clear();   local.vm_p.clear();    local.vm_pT.clear();

    local.vm_xP.clear();      local.vm_yP.clear();      local.vm_zP.clear();      local.vm_pathP.clear();   local.vm_timeP.clear();
    local.vm_planeP.clear();  local.vm_stationP.clear(); local.vm_moduleP.clear(); local.vm_sideP.clear(); local.vm_sensorP.clear();
    local.vm_pdgP.clear();    local.vm_statusP.clear(); local.vm_isPrimaryP.clear(); local.vm_isSelPrimaryP.clear();
    local.vm_mcIndexP.clear(); local.vm_mcCollectionIDP.clear();
    local.vm_pxP.clear();     local.vm_pyP.clear();     local.vm_pzP.clear();     local.vm_pP.clear();     local.vm_pTP.clear();

    local.vm_xEntry.clear();   local.vm_yEntry.clear();  local.vm_zEntry.clear();  local.vm_timeEntry.clear();
    local.vm_pxEntry.clear();  local.vm_pyEntry.clear(); local.vm_pzEntry.clear(); local.vm_pEntry.clear(); local.vm_pTEntry.clear();
    local.vm_xExit.clear();    local.vm_yExit.clear();   local.vm_zExit.clear();   local.vm_timeExit.clear();
    local.vm_pxExit.clear();   local.vm_pyExit.clear();  local.vm_pzExit.clear();  local.vm_pExit.clear();  local.vm_pTExit.clear();
    local.vm_sensorEntry.clear(); local.vm_sensorExit.clear();
    local.vm_cellIDEntry.clear(); local.vm_cellIDExit.clear();
    local.vm_planeEE.clear();  local.vm_stationEE.clear(); local.vm_moduleEE.clear(); local.vm_sideEE.clear();  local.vm_pdgEE.clear();
    local.vm_isPrimaryEE.clear(); local.vm_isSelPrimaryEE.clear();
    local.vm_mcIndexEE.clear(); local.vm_mcCollectionIDEE.clear(); local.vm_nStepsEE.clear();
    local.vm_statusEE.clear();

    local.vm_raw_cellID.clear(); local.vm_raw_charge.clear(); local.vm_raw_timeStamp.clear();
    local.vm_raw_plane.clear(); local.vm_raw_station.clear(); local.vm_raw_module.clear();
    local.vm_raw_side.clear(); local.vm_raw_sensor.clear();
    local.vm_raw_mcIndex.clear(); local.vm_raw_mcCollectionID.clear();
    local.vm_raw_nContribSim.clear(); local.vm_raw_nContribMc.clear(); local.vm_raw_mixedCell.clear();
    local.vm_raw_dominantFrac.clear(); local.vm_raw_totalEdep.clear();
    local.vm_rec_x.clear(); local.vm_rec_y.clear(); local.vm_rec_z.clear();
    local.vm_rec_covxx.clear(); local.vm_rec_covyy.clear(); local.vm_rec_covzz.clear();
    local.vm_rec_time.clear(); local.vm_rec_time_err.clear(); local.vm_rec_edep.clear(); local.vm_rec_edep_err.clear();
    local.vm_rec_cellID.clear();
    local.vm_rec_plane.clear(); local.vm_rec_station.clear(); local.vm_rec_module.clear();
    local.vm_rec_side.clear(); local.vm_rec_sensor.clear();
    local.vm_rec_pixX.clear(); local.vm_rec_pixY.clear(); local.vm_rec_pixZ.clear();
    local.vm_rec_mcIndex.clear(); local.vm_rec_mcCollectionID.clear();
    local.vm_rec_nContribSim.clear(); local.vm_rec_nContribMc.clear(); local.vm_rec_mixedCell.clear();
    local.vm_rec_dominantFrac.clear(); local.vm_rec_totalEdep.clear();

    local.vm_seed_quality.clear(); local.vm_seed_p.clear(); local.vm_seed_qOverP.clear();
    local.vm_seed_theta.clear(); local.vm_seed_phi.clear(); local.vm_seed_loc0.clear(); local.vm_seed_loc1.clear();
    local.vm_seed_sigma_qOverP.clear(); local.vm_seed_sigma_theta.clear(); local.vm_seed_sigma_phi.clear();
    local.vm_seed_nHits.clear(); local.vm_seed_charge.clear();
    local.vm_seed_momentum_resolved.clear(); local.vm_seed_became_track.clear();
    local.vm_seed_made_unfiltered_track.clear(); local.vm_seed_survived_ambiguity.clear();
    local.vm_seed_n_unfiltered_tracks.clear(); local.vm_seed_n_filtered_tracks.clear();
    local.vm_seed_assoc_mcIndex.clear(); local.vm_seed_assoc_mcCollectionID.clear(); local.vm_seed_assoc_weight.clear();
    local.vm_truth_seed_quality.clear(); local.vm_truth_seed_p.clear(); local.vm_truth_seed_qOverP.clear();
    local.vm_truth_seed_theta.clear(); local.vm_truth_seed_phi.clear();
    local.vm_truth_seed_loc0.clear(); local.vm_truth_seed_loc1.clear();
    local.vm_truth_seed_sigma_qOverP.clear(); local.vm_truth_seed_sigma_theta.clear(); local.vm_truth_seed_sigma_phi.clear();
    local.vm_truth_seed_nHits.clear(); local.vm_truth_seed_charge.clear();
    local.vm_truth_seed_momentum_resolved.clear(); local.vm_truth_seed_became_track.clear();
    local.vm_truth_seed_made_unfiltered_track.clear(); local.vm_truth_seed_survived_ambiguity.clear();
    local.vm_truth_seed_n_unfiltered_tracks.clear(); local.vm_truth_seed_n_filtered_tracks.clear();

    const double nan = b0trk::quietNaN();
    local.m_ts.clear(nan);
    local.m_ckf.clear(nan);
    local.m_selPrimaryMcIndex = -1;
    local.m_selPrimaryMcCollectionID = 0;
    local.m_selPrimaryPx = nan;
    local.m_selPrimaryPy = nan;
    local.m_selPrimaryPz = nan;
    local.m_selPrimaryP = nan;
    local.m_selPrimaryPT = nan;
    local.m_selPrimaryThscatMrad = nan;
    local.m_selPrimaryCharge = nan;
    local.m_nStationsPrimary = 0;
    local.m_nSelectedPrimaryMeasurements = hasRawAssocs ? 0 : -1;
    local.m_nMeasurementStationsSelectedPrimary = hasRawAssocs ? 0 : -1;
    local.m_selPrimaryMeasurementReconstructable = hasRawAssocs ? 0 : -1;
    local.m_selPrimaryHasSeed = 0;
    local.m_selPrimaryHasUnfilteredTrack = 0;
    local.m_selPrimaryHasFilteredTrack = 0;
    local.m_selPrimaryHasTruthMatchedTrack = 0;

    local.beam_px.clear();  local.beam_py.clear();  local.beam_pz.clear(); local.beam_p.clear(); local.beam_pT.clear();
    local.beam_pdg.clear();
    local.m_primaryPx.clear(); local.m_primaryPy.clear(); local.m_primaryPz.clear(); local.m_primaryP.clear(); local.m_primaryPT.clear();
    local.m_primaryCharge.clear();
    local.m_primaryPdgOut.clear(); local.m_primaryStatusOut.clear(); local.m_primaryMcIndex.clear();
    local.m_primaryMcCollectionID.clear();
    local.m_genPpx.clear();    local.m_genPpy.clear();    local.m_genPpz.clear();    local.m_genPp.clear();    local.m_genPpT.clear();
    local.m_genBeamPx.clear(); local.m_genBeamPy.clear(); local.m_genBeamPz.clear(); local.m_genBeamP.clear(); local.m_genBeamPT.clear();
    local.m_genBeamPPx.clear();local.m_genBeamPPy.clear();local.m_genBeamPPz.clear();local.m_genBeamPP.clear();local.m_genBeamPPT.clear();

    local.m_nSimHits = static_cast<int>(simHits.size());
    local.m_nRawHits = static_cast<int>(rawHits.size());
    local.m_nRecHits = static_cast<int>(recHits.size());
    local.m_nMeasurements = static_cast<int>(measurements.size());
    local.m_nTruthSeeds = static_cast<int>(truthSeeds.size());
    local.m_nStubSeeds = static_cast<int>(stubSeeds.size());
    local.m_nTsUnfiltered = static_cast<int>(tsUnfiltered.size());
    local.m_nTsFiltered = static_cast<int>(tsEdmTracks.size());
    local.m_nCkfUnfiltered = static_cast<int>(ckfUnfiltered.size());
    local.m_nCkfFiltered = static_cast<int>(ckfEdmTracks.size());
    local.m_nSimHitsUnresolvedCellID = 0;
    local.m_nMissingMcRelation = 0;
    local.m_nSensorMapExact = 0;
    local.m_nSensorMapFallback = 0;
    local.m_nSensorMapFailed = 0;
    local.m_nPixelSnapFailed = 0;

    std::map<std::tuple<uint32_t, int, int, int, int, int>, std::size_t> penetrationIndex;
    std::map<std::tuple<uint32_t, int, int, int, int>, std::size_t> entryExitIndex;

    const auto getFieldOr = [this](std::uint64_t cellID, const char* field, int fallback) -> int {
        try {
            return static_cast<int>(m_decoder->get(cellID, field));
        } catch (...) {
            return fallback;
        }
    };
    const auto decodeIds = [&](std::uint64_t cid, int& plane, int& module, int& sensor, int& side) {
        plane = getFieldOr(cid, "layer", -1);
        module = getFieldOr(cid, "module", -1);
        sensor = getFieldOr(cid, "sensor", -1);
        const auto sideIt = m_moduleToSide.find({plane, module});
        side = (sideIt != m_moduleToSide.end()) ? sideIt->second : -1;
    };

    struct PixelSnap {
        double x = b0trk::quietNaN();
        double y = b0trk::quietNaN();
        double z = b0trk::quietNaN();
        double dx = b0trk::quietNaN();
        double dy = b0trk::quietNaN();
        double dz = b0trk::quietNaN();
        int pixX = -1;
        int pixY = -1;
        int pixZ = -1;
        std::uint64_t cellID = 0;
        bool ok = false;
    };
    const auto snapToAclgadPixel = [this, &getFieldOr](double xMm, double yMm, double zMm,
                                                             std::uint64_t referenceCellID,
                                                             dd4hep::DetElement detElementHint = {}) -> PixelSnap {
        PixelSnap snap;
        if (!std::isfinite(xMm) || !std::isfinite(yMm) || !std::isfinite(zMm) ||
            referenceCellID == 0 || m_decoder == nullptr) {
            return snap;
        }
        try {
            const dd4hep::Position globalPosition(0.1 * xMm, 0.1 * yMm, 0.1 * zMm);
            const auto detElement = detElementHint.isValid()
                ? detElementHint
                : m_volman.lookupDetElement(referenceCellID);
            const auto localPosition = detElement.nominal().worldToLocal(globalPosition);
            const auto volumeID = m_segmentation.volumeID(static_cast<dd4hep::CellID>(referenceCellID));
            const auto snappedCellID =
                m_segmentation.cellID(localPosition, globalPosition, volumeID);
            const auto localCenter = m_segmentation.position(snappedCellID);
            const auto globalCenter = detElement.nominal().localToWorld(localCenter);

            snap.x = 10.0 * globalCenter.x();
            snap.y = 10.0 * globalCenter.y();
            snap.z = 10.0 * globalCenter.z();
            snap.dx = xMm - snap.x;
            snap.dy = yMm - snap.y;
            snap.dz = zMm - snap.z;
            snap.pixX = getFieldOr(static_cast<std::uint64_t>(snappedCellID), "x", -1);
            snap.pixY = getFieldOr(static_cast<std::uint64_t>(snappedCellID), "y", -1);
            snap.pixZ = getFieldOr(static_cast<std::uint64_t>(snappedCellID), "z", -1);
            snap.cellID = static_cast<std::uint64_t>(snappedCellID);
            snap.ok = true;
        } catch (...) {
            snap = PixelSnap{};
        }
        return snap;
    };

    struct ClosestHit {
        const SensorRef* ref = nullptr;
        double localZ = 0.0;
        double outsideX = 0.0;
        double outsideY = 0.0;
        double score = std::numeric_limits<double>::infinity();
    };
    const auto closestSensor = [this](double xMm, double yMm, double zMm) -> ClosestHit {
        ClosestHit best;
        if (!std::isfinite(xMm) || !std::isfinite(yMm) || !std::isfinite(zMm)) {
            return best;
        }
        const dd4hep::Position globalPosition(0.1 * xMm, 0.1 * yMm, 0.1 * zMm);
        for (const auto& sensorRef : m_sensorRefs) {
            try {
                const auto detElement = sensorRef.detElement.isValid()
                    ? sensorRef.detElement
                    : m_volman.lookupDetElement(sensorRef.cellID);
                const auto localPosition = detElement.nominal().worldToLocal(globalPosition);
                const double localX = 10.0 * localPosition.x();
                const double localY = 10.0 * localPosition.y();
                const double localZ = 10.0 * localPosition.z();
                const double outsideX = std::max(0.0, std::abs(localX) - m_sensorHalfX);
                const double outsideY = std::max(0.0, std::abs(localY) - m_sensorHalfY);
                const double score = localZ * localZ + outsideX * outsideX + outsideY * outsideY;
                if (score < best.score) {
                    best.score = score;
                    best.ref = &sensorRef;
                    best.localZ = localZ;
                    best.outsideX = outsideX;
                    best.outsideY = outsideY;
                }
            } catch (...) {
                continue;
            }
        }
        return best;
    };

    std::unordered_set<std::uint64_t> rawCellIDs;
    for (const auto* raw : rawHits) {
        if (raw != nullptr) {
            rawCellIDs.insert(raw->getCellID());
        }
    }

    // Truth attribution of a cell must rank particles by their *total*
    // contribution, not by their largest single Geant4 step. Two 6 keV steps
    // from one particle outweigh a single 10 keV step from another, and a
    // largest-single-deposit rule gets that backwards. The dominant fraction is
    // kept alongside the label so a mixed cell is visible rather than implied.
    struct SimLink {
        int mcIndex = -1;
        std::uint32_t mcCollectionID = 0;
        double eDep = 0.0;         // summed deposit of the dominant particle
        double totalEDep = 0.0;    // summed deposit of the whole cell
        int nContribSim = 0;
        int nContribMc = 0;
    };
    std::unordered_map<std::uint64_t, SimLink> simLinkByCell;
    std::unordered_map<std::uint64_t, std::map<std::pair<std::uint32_t, int>, double>> cellParticleEDepByCell;
    {
        // (cellID, MC ObjectID) -> summed deposit.
        std::map<std::pair<std::uint64_t, std::pair<std::uint32_t, int>>, double> cellParticleEDep;
        for (const auto* assoc : rawAssocs) {
            if (assoc == nullptr) continue;
            const auto raw = assoc->getRawHit();
            const auto sim = assoc->getSimHit();
            if (!raw.isAvailable() || !sim.isAvailable()) continue;
            const auto particle = sim.getParticle();
            if (!particle.isAvailable()) continue;
            const auto id  = particle.id();
            const auto cid = raw.getCellID();
            cellParticleEDep[{cid, {id.collectionID, id.index}}] += sim.getEDep();
            cellParticleEDepByCell[cid][{id.collectionID, id.index}] += sim.getEDep();
            auto& link = simLinkByCell[cid];
            link.totalEDep += sim.getEDep();
            ++link.nContribSim;
        }
        for (const auto& [key, summed] : cellParticleEDep) {
            const auto cid = key.first;
            auto& link     = simLinkByCell[cid];
            ++link.nContribMc;
            if (summed > link.eDep) {
                link.eDep           = summed;
                link.mcCollectionID = key.second.first;
                link.mcIndex        = key.second.second;
            }
        }
    }
    const auto dominantFraction = [nan](const SimLink& link) {
        return link.totalEDep > 0.0 ? link.eDep / link.totalEDep : nan;
    };

    std::vector<std::pair<std::uint32_t, int>> primaryIds;
    for (const auto* part : mcparticles) {
        const int pdg = part->getPDG();
        const int status = part->getGeneratorStatus();
        if (pdg != m_primaryPdg || status != m_primaryStatus) continue;

        const auto id = part->id();
        const auto p = part->getMomentum();
        const double pmag = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        const double pT = std::hypot(p.x, p.y);

        primaryIds.emplace_back(id.collectionID, id.index);
        local.m_primaryPdgOut.push_back(pdg);
        local.m_primaryStatusOut.push_back(status);
        local.m_primaryMcIndex.push_back(id.index);
        local.m_primaryMcCollectionID.push_back(id.collectionID);
        local.m_primaryPx.push_back(p.x);
        local.m_primaryPy.push_back(p.y);
        local.m_primaryPz.push_back(p.z);
        local.m_primaryP.push_back(pmag);
        local.m_primaryPT.push_back(pT);
        // Charge comes from the simulation record, not from a hand-maintained
        // PDG table: primary_pdg is configurable and a negative species would
        // otherwise get a silently wrong-signed q/p pull.
        local.m_primaryCharge.push_back(part->getCharge());
    }

    const auto isSelectedPrimary = [&primaryIds](std::uint32_t collectionID, int index) -> int {
        for (const auto& selected : primaryIds) {
            if (selected.first == collectionID && selected.second == index) return 1;
        }
        return 0;
    };

    std::size_t primaryRef = 0;
    for (std::size_t i = 1; i < local.m_primaryP.size(); ++i) {
        if (local.m_primaryP[i] > local.m_primaryP[primaryRef]) primaryRef = i;
    }
    if (!local.m_primaryP.empty()) {
        local.m_selPrimaryMcIndex = local.m_primaryMcIndex[primaryRef];
        local.m_selPrimaryMcCollectionID = local.m_primaryMcCollectionID[primaryRef];
        local.m_selPrimaryPx = local.m_primaryPx[primaryRef];
        local.m_selPrimaryPy = local.m_primaryPy[primaryRef];
        local.m_selPrimaryPz = local.m_primaryPz[primaryRef];
        local.m_selPrimaryP = local.m_primaryP[primaryRef];
        local.m_selPrimaryPT = local.m_primaryPT[primaryRef];
        local.m_selPrimaryCharge = local.m_primaryCharge[primaryRef];
    }
    const auto isSelPrimary = [this](std::uint32_t collectionID, int index) -> int {
        return (local.m_selPrimaryMcIndex >= 0 &&
                collectionID == local.m_selPrimaryMcCollectionID &&
                index == local.m_selPrimaryMcIndex) ? 1 : 0;
    };

    // A truth crossing is only geometrical acceptance. Reconstructability also
    // requires digitized/reconstructed measurements attributable to the selected
    // primary in enough distinct physical B0 stations. Build the measurement
    // truth label from the summed truth composition of its constituent raw cells.
    local.m_nSelectedPrimaryMeasurements = hasRawAssocs ? 0 : -1;
    local.m_nMeasurementStationsSelectedPrimary = hasRawAssocs ? 0 : -1;
    local.m_selPrimaryMeasurementReconstructable = hasRawAssocs ? 0 : -1;
    if (hasRawAssocs && local.m_selPrimaryMcIndex >= 0) {
        std::set<int> selectedMeasurementStations;
        for (const auto* measurement : measurements) {
            if (measurement == nullptr) continue;
            std::map<std::pair<std::uint32_t, int>, double> measurementTruth;
            std::set<std::uint64_t> seenCells;
            std::set<int> measurementStations;
            for (const auto& hit : measurement->getHits()) {
                const auto raw = hit.getRawHit();
                if (!raw.isAvailable()) continue;
                const auto cid = static_cast<std::uint64_t>(raw.getCellID());
                if (!seenCells.insert(cid).second) continue;
                int plane = -1, module = -1, sensor = -1, side = -1;
                decodeIds(cid, plane, module, sensor, side);
                const int station = stationOf(plane);
                if (station > 0) measurementStations.insert(station);
                const auto truthIt = cellParticleEDepByCell.find(cid);
                if (truthIt == cellParticleEDepByCell.end()) continue;
                for (const auto& [mcId, edep] : truthIt->second) {
                    measurementTruth[mcId] += edep;
                }
            }
            if (measurementTruth.empty()) continue;
            const auto dominant = std::max_element(
                measurementTruth.begin(), measurementTruth.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });
            if (dominant == measurementTruth.end()) continue;
            if (dominant->first.first != local.m_selPrimaryMcCollectionID ||
                dominant->first.second != local.m_selPrimaryMcIndex) {
                continue;
            }
            ++local.m_nSelectedPrimaryMeasurements;
            selectedMeasurementStations.insert(measurementStations.begin(), measurementStations.end());
        }
        local.m_nMeasurementStationsSelectedPrimary =
            static_cast<int>(selectedMeasurementStations.size());
        local.m_selPrimaryMeasurementReconstructable =
            local.m_nMeasurementStationsSelectedPrimary >= m_minMeasurementStations ? 1 : 0;
    }

    for (const auto* h : simHits) {
        const auto mc = h->getParticle();
        if (!mc.isAvailable()) {
            ++local.m_nMissingMcRelation;
            continue;
        }

        const auto mom = h->getMomentum();
        const double pmag = std::sqrt(mom.x*mom.x + mom.y*mom.y + mom.z*mom.z);
        const double pT = std::hypot(mom.x, mom.y);

        const uint64_t cid = h->getCellID();
        dd4hep::Position gpos;
        dd4hep::Position lpos;
        try {
            gpos = m_geoSvc->converter()->position(cid);
            lpos = m_volman.lookupDetElement(cid).nominal()
                       .worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z()));
        } catch (const std::exception&) {
            ++local.m_nSimHitsUnresolvedCellID;
            continue;
        }
        const auto truthPos = h->getPosition();
        int plane = -1, module = -1, sensor = -1, side = -1;
        decodeIds(cid, plane, module, sensor, side);
        const int pixX   = getFieldOr(cid, "x", -1);
        const int pixY   = getFieldOr(cid, "y", -1);
        const int pixZ   = getFieldOr(cid, "z", -1);
        const double path = h->getPathLength();
        const auto id = mc.id();
        const int primaryFlag = isSelectedPrimary(id.collectionID, id.index);
        const int selFlag = isSelPrimary(id.collectionID, id.index);
        const auto truthPixel = snapToAclgadPixel(truthPos.x, truthPos.y, truthPos.z, cid);
        const auto readoutPixel =
            snapToAclgadPixel(10. * gpos.x(), 10. * gpos.y(), 10. * gpos.z(), cid);
        if (!truthPixel.ok || !readoutPixel.ok) {
            ++local.m_nPixelSnapFailed;
        }

        local.vm_xR.push_back(10. * gpos.x());
        local.vm_yR.push_back(10. * gpos.y());
        local.vm_zR.push_back(10. * gpos.z());
        local.vm_xT.push_back(truthPos.x);
        local.vm_yT.push_back(truthPos.y);
        local.vm_zT.push_back(truthPos.z);
        local.vm_aclgad_xPixT.push_back(truthPixel.x);
        local.vm_aclgad_yPixT.push_back(truthPixel.y);
        local.vm_aclgad_zPixT.push_back(truthPixel.z);
        local.vm_aclgad_dxT.push_back(truthPixel.dx);
        local.vm_aclgad_dyT.push_back(truthPixel.dy);
        local.vm_aclgad_dzT.push_back(truthPixel.dz);
        local.vm_aclgad_xPixR.push_back(readoutPixel.x);
        local.vm_aclgad_yPixR.push_back(readoutPixel.y);
        local.vm_aclgad_zPixR.push_back(readoutPixel.z);
        local.vm_aclgad_dxR.push_back(readoutPixel.dx);
        local.vm_aclgad_dyR.push_back(readoutPixel.dy);
        local.vm_aclgad_dzR.push_back(readoutPixel.dz);
        local.vm_detX.push_back(10. * lpos.x());
        local.vm_detY.push_back(10. * lpos.y());
        local.vm_detZ.push_back(10. * lpos.z());
        local.vm_plane .push_back(plane);
        local.vm_station.push_back(stationOf(plane));
        local.vm_module.push_back(module);
        local.vm_side  .push_back(side);
        local.vm_sensor.push_back(sensor);
        local.vm_pixX  .push_back(pixX);
        local.vm_pixY  .push_back(pixY);
        local.vm_pixZ  .push_back(pixZ);
        local.vm_aclgad_pixXT.push_back(truthPixel.pixX);
        local.vm_aclgad_pixYT.push_back(truthPixel.pixY);
        local.vm_aclgad_pixXR.push_back(readoutPixel.pixX);
        local.vm_aclgad_pixYR.push_back(readoutPixel.pixY);
        local.vm_cellID.push_back(cid);
        local.vm_eDep  .push_back(h->getEDep());
        local.vm_time  .push_back(h->getTime());
        local.vm_path  .push_back(path);
        local.vm_pdg   .push_back(mc.getPDG());
        local.vm_mcIndex.push_back(id.index);
        local.vm_mcCollectionID.push_back(id.collectionID);
        local.vm_px    .push_back(mom.x);
        local.vm_py    .push_back(mom.y);
        local.vm_pz    .push_back(mom.z);
        local.vm_p     .push_back(pmag);
        local.vm_pT    .push_back(pT);
        local.vm_status.push_back(mc.getGeneratorStatus());
        local.vm_isPrimary.push_back(primaryFlag);
        local.vm_isSelPrimary.push_back(selFlag);
        // Cell-level, not SimHit-level: the generic digitizer associates every
        // SimHit sharing a fired cellID, including subthreshold ones, so this
        // says the cell produced a RawHit -- not that this deposit did.
        local.vm_cellFired.push_back(rawCellIDs.count(cid) ? 1 : 0);

        const auto key = std::make_tuple(id.collectionID, id.index, plane, side, module, sensor);
        const auto existing = penetrationIndex.find(key);
        if (existing == penetrationIndex.end()) {
            penetrationIndex.emplace(key, local.vm_xP.size());
            local.vm_xP     .push_back(truthPos.x);
            local.vm_yP     .push_back(truthPos.y);
            local.vm_zP     .push_back(truthPos.z);
            local.vm_pathP  .push_back(path);
            local.vm_timeP  .push_back(h->getTime());
            local.vm_planeP .push_back(plane);
            local.vm_stationP.push_back(stationOf(plane));
            local.vm_moduleP.push_back(module);
            local.vm_sideP  .push_back(side);
            local.vm_sensorP.push_back(sensor);
            local.vm_pdgP   .push_back(mc.getPDG());
            local.vm_statusP.push_back(mc.getGeneratorStatus());
            local.vm_isPrimaryP.push_back(primaryFlag);
            local.vm_isSelPrimaryP.push_back(selFlag);
            local.vm_mcIndexP.push_back(id.index);
            local.vm_mcCollectionIDP.push_back(id.collectionID);
            local.vm_pxP    .push_back(mom.x);
            local.vm_pyP    .push_back(mom.y);
            local.vm_pzP    .push_back(mom.z);
            local.vm_pP     .push_back(pmag);
            local.vm_pTP    .push_back(pT);
        } else if (h->getTime() < local.vm_timeP[existing->second]) {
            const std::size_t idx = existing->second;
            local.vm_xP[idx]    = truthPos.x;
            local.vm_yP[idx]    = truthPos.y;
            local.vm_zP[idx]    = truthPos.z;
            local.vm_pathP[idx] = path;
            local.vm_timeP[idx] = h->getTime();
            local.vm_pxP[idx]   = mom.x;
            local.vm_pyP[idx]   = mom.y;
            local.vm_pzP[idx]   = mom.z;
            local.vm_pP[idx]    = pmag;
            local.vm_pTP[idx]   = pT;
        }

        if (side >= 0) {
            const auto eeKey = std::make_tuple(id.collectionID, id.index, plane, side, module);
            const double thisTime = h->getTime();
            const auto eeIt = entryExitIndex.find(eeKey);
            if (eeIt == entryExitIndex.end()) {
                entryExitIndex.emplace(eeKey, local.vm_xEntry.size());
                local.vm_xEntry .push_back(truthPos.x);   local.vm_yEntry .push_back(truthPos.y);
                local.vm_zEntry .push_back(truthPos.z);   local.vm_timeEntry.push_back(thisTime);
                local.vm_pxEntry.push_back(mom.x);        local.vm_pyEntry.push_back(mom.y);
                local.vm_pzEntry.push_back(mom.z);        local.vm_pEntry .push_back(pmag);
                local.vm_pTEntry.push_back(pT);
                local.vm_xExit  .push_back(truthPos.x);   local.vm_yExit  .push_back(truthPos.y);
                local.vm_zExit  .push_back(truthPos.z);   local.vm_timeExit.push_back(thisTime);
                local.vm_pxExit .push_back(mom.x);        local.vm_pyExit .push_back(mom.y);
                local.vm_pzExit .push_back(mom.z);        local.vm_pExit  .push_back(pmag);
                local.vm_pTExit .push_back(pT);
                local.vm_sensorEntry.push_back(sensor);   local.vm_sensorExit.push_back(sensor);
                local.vm_cellIDEntry.push_back(cid);      local.vm_cellIDExit.push_back(cid);
                local.vm_planeEE.push_back(plane);
                local.vm_stationEE.push_back(stationOf(plane));
                local.vm_moduleEE.push_back(module);
                local.vm_sideEE.push_back(side);
                local.vm_pdgEE  .push_back(mc.getPDG());
                local.vm_isPrimaryEE.push_back(primaryFlag);
                local.vm_isSelPrimaryEE.push_back(selFlag);
                local.vm_mcIndexEE.push_back(id.index);
                local.vm_mcCollectionIDEE.push_back(id.collectionID);
                local.vm_nStepsEE.push_back(1);
                local.vm_statusEE.push_back(mc.getGeneratorStatus());
            } else {
                const std::size_t idx = eeIt->second;
                ++local.vm_nStepsEE[idx];
                if (thisTime < local.vm_timeEntry[idx]) {
                    local.vm_xEntry[idx]    = truthPos.x;
                    local.vm_yEntry[idx]    = truthPos.y;
                    local.vm_zEntry[idx]    = truthPos.z;
                    local.vm_timeEntry[idx] = thisTime;
                    local.vm_pxEntry[idx]   = mom.x;
                    local.vm_pyEntry[idx]   = mom.y;
                    local.vm_pzEntry[idx]   = mom.z;
                    local.vm_pEntry[idx]    = pmag;
                    local.vm_pTEntry[idx]   = pT;
                    local.vm_sensorEntry[idx] = sensor;
                    local.vm_cellIDEntry[idx] = cid;
                }
                if (thisTime > local.vm_timeExit[idx]) {
                    local.vm_xExit[idx]    = truthPos.x;
                    local.vm_yExit[idx]    = truthPos.y;
                    local.vm_zExit[idx]    = truthPos.z;
                    local.vm_timeExit[idx] = thisTime;
                    local.vm_pxExit[idx]   = mom.x;
                    local.vm_pyExit[idx]   = mom.y;
                    local.vm_pzExit[idx]   = mom.z;
                    local.vm_pExit[idx]    = pmag;
                    local.vm_pTExit[idx]   = pT;
                    local.vm_sensorExit[idx] = sensor;
                    local.vm_cellIDExit[idx] = cid;
                }
            }
        }
    }

    for (const auto* part : mcparticles) {
        const auto p = part->getMomentum();
        const int  pdg    = part->getPDG();
        const int  status = part->getGeneratorStatus();
        const double pmag = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
        const double pT = std::hypot(p.x, p.y);

        if (status == 4 &&
            (pdg == 22 || pdg == 11 || pdg == -11 || pdg == 2212)) {
            local.beam_px .push_back(p.x);
            local.beam_py .push_back(p.y);
            local.beam_pz .push_back(p.z);
            local.beam_p  .push_back(pmag);
            local.beam_pT .push_back(pT);
            local.beam_pdg.push_back(pdg);
        }
        if (pdg == 2212 && status == 4) {
            local.m_genPpx.push_back(p.x);
            local.m_genPpy.push_back(p.y);
            local.m_genPpz.push_back(p.z);
            local.m_genPp .push_back(pmag);
            local.m_genPpT.push_back(pT);
        }
        if (pdg == 2212 && status == 1) {
            local.m_genBeamPPx.push_back(p.x);
            local.m_genBeamPPy.push_back(p.y);
            local.m_genBeamPPz.push_back(p.z);
            local.m_genBeamPP .push_back(pmag);
            local.m_genBeamPPT.push_back(pT);
        }
        if (pdg == 11 && status == 1) {
            local.m_genBeamPx.push_back(p.x);
            local.m_genBeamPy.push_back(p.y);
            local.m_genBeamPz.push_back(p.z);
            local.m_genBeamP .push_back(pmag);
            local.m_genBeamPT.push_back(pT);
        }
    }

    if (std::isfinite(local.m_selPrimaryP) && !local.m_genPp.empty()) {
        const double bx = local.m_genPpx[0];
        const double by = local.m_genPpy[0];
        const double bz = local.m_genPpz[0];
        const double bn = std::sqrt(bx * bx + by * by + bz * bz);
        const double pn = local.m_selPrimaryP;
        if (bn > 0.0 && pn > 0.0) {
            const double cosang = std::clamp(
                (local.m_selPrimaryPx * bx + local.m_selPrimaryPy * by + local.m_selPrimaryPz * bz) / (pn * bn),
                -1.0, 1.0);
            local.m_selPrimaryThscatMrad = 1.0e3 * std::acos(cosang);
        }
    }
    {
        std::set<int> stations;
        for (std::size_t i = 0; i < local.vm_stationP.size(); ++i) {
            if (local.vm_isSelPrimaryP[i] == 1 && local.vm_stationP[i] > 0) {
                stations.insert(local.vm_stationP[i]);
            }
        }
        local.m_nStationsPrimary = static_cast<int>(stations.size());
    }

    for (const auto* raw : rawHits) {
        if (raw == nullptr) continue;
        const auto cid = raw->getCellID();
        int plane = -1, module = -1, sensor = -1, side = -1;
        decodeIds(cid, plane, module, sensor, side);
        const auto linkIt = simLinkByCell.find(cid);
        local.vm_raw_cellID.push_back(cid);
        local.vm_raw_charge.push_back(raw->getCharge());
        local.vm_raw_timeStamp.push_back(raw->getTimeStamp());
        local.vm_raw_plane.push_back(plane);
        local.vm_raw_station.push_back(stationOf(plane));
        local.vm_raw_module.push_back(module);
        local.vm_raw_side.push_back(side);
        local.vm_raw_sensor.push_back(sensor);
        local.vm_raw_mcIndex.push_back(linkIt == simLinkByCell.end() ? -1 : linkIt->second.mcIndex);
        local.vm_raw_mcCollectionID.push_back(linkIt == simLinkByCell.end() ? 0 : linkIt->second.mcCollectionID);
        local.vm_raw_nContribSim.push_back(linkIt == simLinkByCell.end() ? 0 : linkIt->second.nContribSim);
        local.vm_raw_nContribMc.push_back(linkIt == simLinkByCell.end() ? 0 : linkIt->second.nContribMc);
        local.vm_raw_dominantFrac.push_back(linkIt == simLinkByCell.end() ? nan : dominantFraction(linkIt->second));
        local.vm_raw_totalEdep.push_back(linkIt == simLinkByCell.end() ? nan : linkIt->second.totalEDep);
        local.vm_raw_mixedCell.push_back(linkIt == simLinkByCell.end() ? -1 : (linkIt->second.nContribMc > 1 ? 1 : 0));
    }

    for (const auto* rec : recHits) {
        if (rec == nullptr) continue;
        const auto cid = rec->getCellID();
        int plane = -1, module = -1, sensor = -1, side = -1;
        decodeIds(cid, plane, module, sensor, side);
        const auto pos = rec->getPosition();
        const auto err = rec->getPositionError();
        const auto linkIt = simLinkByCell.find(cid);
        local.vm_rec_x.push_back(pos.x);
        local.vm_rec_y.push_back(pos.y);
        local.vm_rec_z.push_back(pos.z);
        local.vm_rec_covxx.push_back(err.xx);
        local.vm_rec_covyy.push_back(err.yy);
        local.vm_rec_covzz.push_back(err.zz);
        local.vm_rec_time.push_back(rec->getTime());
        local.vm_rec_time_err.push_back(rec->getTimeError());
        local.vm_rec_edep.push_back(rec->getEdep());
        local.vm_rec_edep_err.push_back(rec->getEdepError());
        local.vm_rec_cellID.push_back(cid);
        local.vm_rec_plane.push_back(plane);
        local.vm_rec_station.push_back(stationOf(plane));
        local.vm_rec_module.push_back(module);
        local.vm_rec_side.push_back(side);
        local.vm_rec_sensor.push_back(sensor);
        local.vm_rec_pixX.push_back(getFieldOr(cid, "x", -1));
        local.vm_rec_pixY.push_back(getFieldOr(cid, "y", -1));
        local.vm_rec_pixZ.push_back(getFieldOr(cid, "z", -1));
        local.vm_rec_mcIndex.push_back(linkIt == simLinkByCell.end() ? -1 : linkIt->second.mcIndex);
        local.vm_rec_mcCollectionID.push_back(linkIt == simLinkByCell.end() ? 0 : linkIt->second.mcCollectionID);
        local.vm_rec_nContribSim.push_back(linkIt == simLinkByCell.end() ? 0 : linkIt->second.nContribSim);
        local.vm_rec_nContribMc.push_back(linkIt == simLinkByCell.end() ? 0 : linkIt->second.nContribMc);
        local.vm_rec_dominantFrac.push_back(linkIt == simLinkByCell.end() ? nan : dominantFraction(linkIt->second));
        local.vm_rec_totalEdep.push_back(linkIt == simLinkByCell.end() ? nan : linkIt->second.totalEDep);
        local.vm_rec_mixedCell.push_back(linkIt == simLinkByCell.end() ? -1 : (linkIt->second.nContribMc > 1 ? 1 : 0));
    }

    // Seed survival is two distinct stages: the CKF either builds a trajectory
    // from the seed or it does not, and the ambiguity solver then either keeps
    // that trajectory or drops it. Counting only the filtered collection merges
    // a CKF failure with an ambiguity rejection, which are different problems
    // with different fixes.
    const auto countSeedUse = [](const auto& trajectories,
                                 std::unordered_map<int, int>& counts) {
        for (const auto* traj : trajectories) {
            if (traj == nullptr) continue;
            const auto seed = traj->getSeed();
            if (seed.isAvailable()) {
                ++counts[seed.id().index];
            }
        }
    };

    const auto fillSeeds = [&](const auto& seeds, const auto& trajectories,
                               const auto& trajectoriesUnfiltered, bool haveUnfiltered,
                               auto& quality, auto& p, auto& qOverP, auto& theta, auto& phi,
                               auto& loc0, auto& loc1, auto& sigQ, auto& sigTh, auto& sigPh,
                               auto& nHits, auto& charge, auto& resolved, auto& became,
                               auto& madeUnfiltered, auto& survivedAmbiguity,
                               auto& nUnfiltered, auto& nFiltered) {
        std::unordered_map<int, int> filteredCounts;
        std::unordered_map<int, int> unfilteredCounts;
        countSeedUse(trajectories, filteredCounts);
        countSeedUse(trajectoriesUnfiltered, unfilteredCounts);
        for (std::size_t i = 0; i < seeds.size(); ++i) {
            const auto* seed = seeds[i];
            if (seed == nullptr) continue;
            // Key on the PODIO ObjectID rather than the loop position: the two
            // coincide for an unfiltered collection but the invariant is not
            // guaranteed and must not be assumed.
            const int seedIndex = seed->id().index;
            const auto filteredIt   = filteredCounts.find(seedIndex);
            const auto unfilteredIt = unfilteredCounts.find(seedIndex);
            const int nFilt   = filteredIt == filteredCounts.end() ? 0 : filteredIt->second;
            const int nUnfilt = unfilteredIt == unfilteredCounts.end() ? 0 : unfilteredIt->second;
            auto tp = seed->getParams();
            double qop = nan, th = nan, ph = nan, l0 = nan, l1 = nan;
            double sq = nan, st = nan, sp = nan;
            int ch = 0;
            bool momOk = false;
            double pm = nan;
            if (tp.isAvailable()) {
                qop = tp.getQOverP();
                th = tp.getTheta();
                ph = tp.getPhi();
                l0 = tp.getLoc().a;
                l1 = tp.getLoc().b;
                const auto& cov = tp.getCovariance();
                const auto sig = [&](unsigned k) {
                    const double var = cov(k, k);
                    return (std::isfinite(var) && var >= 0.0) ? std::sqrt(var) : nan;
                };
                sq = sig(4); st = sig(3); sp = sig(2);
                const auto [pp, ok] = b0trk::momentumFromQOverP(qop);
                pm = pp;
                momOk = ok;
                ch = (qop > 0.0) ? 1 : ((qop < 0.0) ? -1 : 0);
            }
            quality.push_back(seed->getQuality());
            p.push_back(pm);
            qOverP.push_back(qop);
            theta.push_back(th);
            phi.push_back(ph);
            loc0.push_back(l0);
            loc1.push_back(l1);
            sigQ.push_back(sq);
            sigTh.push_back(st);
            sigPh.push_back(sp);
            nHits.push_back(static_cast<int>(seed->hits_size()));
            charge.push_back(ch);
            resolved.push_back(momOk ? 1 : 0);
            became.push_back(nFilt > 0 ? 1 : 0);
            // -1 marks "unknowable": without the unfiltered collection the two
            // stages cannot be separated, and reporting 0 would be a claim that
            // the CKF failed.
            madeUnfiltered.push_back(haveUnfiltered ? (nUnfilt > 0 ? 1 : 0) : -1);
            survivedAmbiguity.push_back(
                haveUnfiltered ? (nUnfilt > 0 ? (nFilt > 0 ? 1 : 0) : -1) : -1);
            nUnfiltered.push_back(haveUnfiltered ? nUnfilt : -1);
            nFiltered.push_back(nFilt);
        }
    };
    fillSeeds(stubSeeds, ckfTrajectories, ckfTrajectoriesUnfiltered, hasCkfTrajectoriesUnfiltered,
              local.vm_seed_quality, local.vm_seed_p, local.vm_seed_qOverP,
              local.vm_seed_theta, local.vm_seed_phi, local.vm_seed_loc0, local.vm_seed_loc1,
              local.vm_seed_sigma_qOverP, local.vm_seed_sigma_theta, local.vm_seed_sigma_phi,
              local.vm_seed_nHits, local.vm_seed_charge, local.vm_seed_momentum_resolved, local.vm_seed_became_track,
              local.vm_seed_made_unfiltered_track, local.vm_seed_survived_ambiguity,
              local.vm_seed_n_unfiltered_tracks, local.vm_seed_n_filtered_tracks);
    fillSeeds(truthSeeds, tsTrajectories, tsTrajectoriesUnfiltered, hasTsTrajectoriesUnfiltered,
              local.vm_truth_seed_quality, local.vm_truth_seed_p, local.vm_truth_seed_qOverP,
              local.vm_truth_seed_theta, local.vm_truth_seed_phi, local.vm_truth_seed_loc0, local.vm_truth_seed_loc1,
              local.vm_truth_seed_sigma_qOverP, local.vm_truth_seed_sigma_theta, local.vm_truth_seed_sigma_phi,
              local.vm_truth_seed_nHits, local.vm_truth_seed_charge, local.vm_truth_seed_momentum_resolved,
              local.vm_truth_seed_became_track,
              local.vm_truth_seed_made_unfiltered_track, local.vm_truth_seed_survived_ambiguity,
              local.vm_truth_seed_n_unfiltered_tracks, local.vm_truth_seed_n_filtered_tracks);

    // Attribute each stub seed to MC truth from its constituent TrackerHits.
    // Each hit/cell contributes one unit distributed among the cell's MC energy fractions,
    // matching the per-measurement convention used by ActsToTracks.
    for (const auto* seed : stubSeeds) {
        std::map<std::pair<std::uint32_t, int>, double> weights;
        if (seed != nullptr) {
            for (std::size_t ih = 0; ih < seed->hits_size(); ++ih) {
                const auto hit = seed->getHits(ih);
                if (!hit.isAvailable()) continue;
                const auto found = cellParticleEDepByCell.find(hit.getCellID());
                if (found == cellParticleEDepByCell.end()) continue;
                double total = 0.0;
                for (const auto& [id, edep] : found->second) total += edep;
                if (!(total > 0.0)) continue;
                for (const auto& [id, edep] : found->second) weights[id] += edep / total;
            }
        }
        int bestIndex = -1;
        std::uint32_t bestCollection = 0;
        double bestWeight = -1.0;
        double totalWeight = 0.0;
        for (const auto& [id, weight] : weights) {
            totalWeight += weight;
            if (weight > bestWeight) {
                bestWeight = weight;
                bestCollection = id.first;
                bestIndex = id.second;
            }
        }
        const double normalized = totalWeight > 0.0 ? bestWeight / totalWeight : nan;
        local.vm_seed_assoc_mcIndex.push_back(bestIndex);
        local.vm_seed_assoc_mcCollectionID.push_back(bestCollection);
        local.vm_seed_assoc_weight.push_back(normalized);
        if (bestIndex == local.m_selPrimaryMcIndex && bestCollection == local.m_selPrimaryMcCollectionID &&
            local.m_selPrimaryMcIndex >= 0) {
            local.m_selPrimaryHasSeed = 1;
        }
    }

    const auto collectionHasSelectedPrimary = [this](const auto& assocs) {
        std::map<std::pair<std::uint32_t, int>, std::tuple<double, std::uint32_t, int>> best;
        for (const auto* assoc : assocs) {
            if (assoc == nullptr) continue;
            const auto rec = assoc->getRec();
            const auto sim = assoc->getSim();
            if (!rec.isAvailable() || !sim.isAvailable()) continue;
            const auto rid = rec.id();
            const auto sid = sim.id();
            const auto key = std::make_pair(rid.collectionID, rid.index);
            const double weight = assoc->getWeight();
            const auto it = best.find(key);
            if (it == best.end() || weight > std::get<0>(it->second)) {
                best[key] = {weight, sid.collectionID, sid.index};
            }
        }
        for (const auto& [key, value] : best) {
            if (std::get<1>(value) == local.m_selPrimaryMcCollectionID &&
                std::get<2>(value) == local.m_selPrimaryMcIndex && local.m_selPrimaryMcIndex >= 0) return true;
        }
        return false;
    };
    local.m_selPrimaryHasUnfilteredTrack =
        hasCkfAssocsUnfiltered && collectionHasSelectedPrimary(ckfUnfilteredAssocs) ? 1 : 0;
    local.m_selPrimaryHasFilteredTrack =
        hasCkfAssocs && collectionHasSelectedPrimary(ckfAssocs) ? 1 : 0;

    auto fillChain = [&](TrackChain& out,
                         const auto& trajectories,
                         const auto& tracks,
                         const auto& edmTracks,
                         const auto& assocs,
                         const auto& actsTracks,
                         const auto& actsStates,
                         const auto& chainSeeds) {
        struct AssocHit {
            int mcIndex = -1;
            std::uint32_t mcCollectionID = 0;
            double weight = 0.0;
        };
        std::map<std::pair<std::uint32_t, int>, AssocHit> bestAssoc;
        for (const auto* assoc : assocs) {
            if (assoc == nullptr) continue;
            const auto rec = assoc->getRec();
            const auto sim = assoc->getSim();
            if (!rec.isAvailable() || !sim.isAvailable()) continue;
            const auto recId = rec.id();
            const double weight = assoc->getWeight();
            const auto recKey = std::make_pair(recId.collectionID, recId.index);
            auto it = bestAssoc.find(recKey);
            if (it == bestAssoc.end() || weight > it->second.weight) {
                bestAssoc[recKey] = {sim.id().index, sim.id().collectionID, weight};
            }
        }

        std::map<std::pair<std::uint32_t, int>, const edm4eic::Track*> edmTrackByTrajectory;
        std::map<std::pair<std::uint32_t, int>, std::vector<std::pair<std::uint32_t, int>>> trackObjectsBySeed;
        std::map<std::pair<std::uint32_t, int>, int> pdgByTrackObject;
        for (const auto* edmTrack : edmTracks) {
            if (edmTrack == nullptr) continue;
            const auto objectId = edmTrack->id();
            pdgByTrackObject[{objectId.collectionID, objectId.index}] = edmTrack->getPdg();
            const auto trajectory = edmTrack->getTrajectory();
            if (!trajectory.isAvailable()) continue;
            const auto trajectoryId = trajectory.id();
            edmTrackByTrajectory[{trajectoryId.collectionID, trajectoryId.index}] = edmTrack;
            const auto seed = trajectory.getSeed();
            if (seed.isAvailable()) {
                const auto seedId = seed.id();
                trackObjectsBySeed[{seedId.collectionID, seedId.index}].push_back(
                    {objectId.collectionID, objectId.index});
            }
        }

        const double primaryP = local.m_selPrimaryP;
        const double primaryPT = local.m_selPrimaryPT;
        const double truthQOverP =
            (std::isfinite(primaryP) && primaryP > 0.0 && std::isfinite(local.m_selPrimaryCharge))
                ? (local.m_selPrimaryCharge / primaryP) : nan;
        const double truthTheta = (std::isfinite(primaryP) && primaryP > 0.0)
            ? std::acos(std::clamp(local.m_selPrimaryPz / primaryP, -1.0, 1.0)) : nan;
        const double truthPhi = (std::isfinite(local.m_selPrimaryPx) && std::isfinite(local.m_selPrimaryPy))
            ? std::atan2(local.m_selPrimaryPy, local.m_selPrimaryPx) : nan;

        double bestAbsDeltaP = std::numeric_limits<double>::infinity();
        double bestTruthWeight = -1.0;
        int bestTruthMeas = -1;
        double bestTruthChi2 = std::numeric_limits<double>::infinity();
        int bestRecoMeas = -1;
        int bestRecoHoles = std::numeric_limits<int>::max();
        double bestRecoChi2Ndf = std::numeric_limits<double>::infinity();

        const auto assignBest = [&](BestSel& dest, int idx, double p, double pT, double dP, double dPT,
                                    double theta, double phi, int nStates, int nMeas, int nOut,
                                    int nHoles, int nShared, double chi2, int ndf, int assocMc,
                                    double assocW, int pdg, int charge, int momOk,
                                    double pullQ, double pullTh, double pullPh) {
            dest.index = idx;
            dest.p = p;
            dest.pT = pT;
            dest.deltaP = dP;
            dest.deltaPT = dPT;
            dest.theta = theta;
            dest.phi = phi;
            dest.nStates = nStates;
            dest.nMeasurements = nMeas;
            dest.nOutliers = nOut;
            dest.nHoles = nHoles;
            dest.nSharedHits = nShared;
            dest.chi2 = chi2;
            dest.ndf = ndf;
            dest.assocMcIndex = assocMc;
            dest.assocWeight = assocW;
            dest.pdg = pdg;
            dest.charge = charge;
            dest.momentumResolved = momOk;
            dest.pullQOverP = pullQ;
            dest.pullTheta = pullTh;
            dest.pullPhi = pullPh;
        };

        for (std::size_t trajIndex = 0; trajIndex < trajectories.size(); ++trajIndex) {
            const auto* trajectory = trajectories[trajIndex];
            if (trajectory == nullptr) continue;
            auto tp = (trajectory->trackParameters_size() > 0)
                ? trajectory->getTrackParameters(0)
                : edm4eic::TrackParameters::makeEmpty();
            // Schema 3 never guesses trajectory identity from parallel collection
            // positions. If the trajectory does not carry fitted parameters, leave
            // this object unresolved rather than borrowing tracks[trajIndex].
            if (!tp.isAvailable()) continue;

            const float theta  = tp.getTheta();
            const float phi    = tp.getPhi();
            const float qOverP = tp.getQOverP();
            const auto [p, momOk] = b0trk::momentumFromQOverP(qOverP);
            const double pT = (momOk && std::isfinite(theta)) ? std::abs(p * std::sin(theta)) : nan;
            const int charge = (qOverP > 0.f) ? 1 : ((qOverP < 0.f) ? -1 : 0);
            const double deltaP = (momOk && std::isfinite(primaryP)) ? p - primaryP : nan;
            const double deltaPT = (std::isfinite(pT) && std::isfinite(primaryPT)) ? pT - primaryPT : nan;
            const int currentTrackIndex = static_cast<int>(trajIndex); // legacy positional index
            int objectIndex = -1;
            std::uint32_t objectCollectionID = 0;
            int seedIndex = -1;
            std::uint32_t seedCollectionID = 0;
            int identityValid = 0;
            const auto trajectoryId = trajectory->id();
            const auto stableTrackIt = edmTrackByTrajectory.find({trajectoryId.collectionID, trajectoryId.index});
            const edm4eic::Track* stableEdmTrack =
                stableTrackIt == edmTrackByTrajectory.end() ? nullptr : stableTrackIt->second;
            if (stableEdmTrack != nullptr) {
                const auto objectId = stableEdmTrack->id();
                objectIndex = objectId.index;
                objectCollectionID = objectId.collectionID;
                identityValid = 1;
            }
            const auto trajectorySeed = trajectory->getSeed();
            if (trajectorySeed.isAvailable()) {
                const auto seedId = trajectorySeed.id();
                seedIndex = seedId.index;
                seedCollectionID = seedId.collectionID;
            }
            const int nStates = static_cast<int>(trajectory->getNStates());
            const int nMeasurements = static_cast<int>(trajectory->getNMeasurements());
            const int nOutliers = static_cast<int>(trajectory->getNOutliers());
            const int nHoles = static_cast<int>(trajectory->getNHoles());
            const int nSharedHits = static_cast<int>(trajectory->getNSharedHits());

            double chi2 = nan;
            int ndf = -1;
            int pdg = 0;
            if (stableEdmTrack != nullptr) {
                chi2 = stableEdmTrack->getChi2();
                ndf = static_cast<int>(stableEdmTrack->getNdf());
                pdg = stableEdmTrack->getPdg();
            }
            int assocMc = -1;
            std::uint32_t assocCol = 0;
            double assocW = 0.0;
            if (objectIndex >= 0) {
                if (const auto it = bestAssoc.find({objectCollectionID, objectIndex}); it != bestAssoc.end()) {
                    assocMc = it->second.mcIndex;
                    assocCol = it->second.mcCollectionID;
                    assocW = it->second.weight;
                }
            }

            const auto& loc = tp.getLoc();
            const auto& cov = tp.getCovariance();
            const auto diagSigma = [&cov, nan](unsigned k) -> double {
                const double var = cov(k, k);
                return (std::isfinite(var) && var >= 0.0) ? std::sqrt(var) : nan;
            };
            const double sigQ = diagSigma(4);
            const double sigTh = diagSigma(3);
            const double sigPh = diagSigma(2);
            const bool truthAssoc = (assocMc == local.m_selPrimaryMcIndex &&
                                     assocCol == local.m_selPrimaryMcCollectionID &&
                                     local.m_selPrimaryMcIndex >= 0);
            const double pullQ = (truthAssoc && std::isfinite(qOverP) && std::isfinite(truthQOverP) &&
                                  std::isfinite(sigQ) && sigQ > 0.0)
                ? (qOverP - truthQOverP) / sigQ : nan;
            const double pullTh = (truthAssoc && std::isfinite(theta) && std::isfinite(truthTheta) &&
                                   std::isfinite(sigTh) && sigTh > 0.0)
                ? (theta - truthTheta) / sigTh : nan;
            const double pullPh = (truthAssoc && std::isfinite(phi) && std::isfinite(truthPhi) &&
                                   std::isfinite(sigPh) && sigPh > 0.0)
                ? b0trk::wrapPi(phi - truthPhi) / sigPh : nan;

            out.p.push_back(p);
            out.pT.push_back(pT);
            out.delta_p.push_back(deltaP);
            out.delta_pT.push_back(deltaPT);
            out.theta.push_back(theta);
            out.phi.push_back(phi);
            out.px.push_back(std::isfinite(pT) ? pT * std::cos(phi) : nan);
            out.py.push_back(std::isfinite(pT) ? pT * std::sin(phi) : nan);
            out.pz.push_back(momOk ? p * std::cos(theta) : nan);
            out.qOverP.push_back(qOverP);
            out.momentum_resolved.push_back(momOk ? 1 : 0);
            out.loc0.push_back(loc.a);
            out.loc1.push_back(loc.b);
            out.sigma_loc0.push_back(diagSigma(0));
            out.sigma_loc1.push_back(diagSigma(1));
            out.sigma_phi.push_back(sigPh);
            out.sigma_theta.push_back(sigTh);
            out.sigma_qOverP.push_back(sigQ);
            out.sigma_time.push_back(diagSigma(5));
            out.pull_qOverP.push_back(pullQ);
            out.pull_theta.push_back(pullTh);
            out.pull_phi.push_back(pullPh);
            out.chi2.push_back(chi2);
            out.ndf.push_back(ndf);
            out.charge.push_back(charge);
            out.index.push_back(currentTrackIndex);
            out.object_index.push_back(objectIndex);
            out.object_collectionID.push_back(objectCollectionID);
            out.seed_index.push_back(seedIndex);
            out.seed_collectionID.push_back(seedCollectionID);
            out.identity_valid.push_back(identityValid);
            {
                int k = 0;
                for (unsigned i = 0; i < 6; ++i) {
                    for (unsigned j = i; j < 6; ++j) out.cov_upper[k++].push_back(cov(i, j));
                }
            }
            out.type.push_back(tp.getType());
            out.time.push_back(tp.getTime());
            out.pdg.push_back(pdg);
            out.nStates.push_back(nStates);
            out.nMeasurements.push_back(nMeasurements);
            out.nOutliers.push_back(nOutliers);
            out.nHoles.push_back(nHoles);
            out.nSharedHits.push_back(nSharedHits);
            out.assoc_mcIndex.push_back(assocMc);
            out.assoc_mcCollectionID.push_back(assocCol);
            out.assoc_weight.push_back(assocW);

            if (std::isfinite(deltaP) && std::abs(deltaP) < bestAbsDeltaP) {
                bestAbsDeltaP = std::abs(deltaP);
                assignBest(out.oracle, currentTrackIndex, p, pT, deltaP, deltaPT, theta, phi,
                           nStates, nMeasurements, nOutliers, nHoles, nSharedHits, chi2, ndf,
                           assocMc, assocW, pdg, charge, momOk ? 1 : 0, pullQ, pullTh, pullPh);
                out.oracle.objectIndex = objectIndex; out.oracle.objectCollectionID = objectCollectionID;
                out.oracle.seedIndex = seedIndex; out.oracle.seedCollectionID = seedCollectionID;
                out.oracle.identityValid = identityValid;
                out.oracle.assocMcCollectionID = assocCol;
            }
            if (truthAssoc) {
                const bool better =
                    (assocW > bestTruthWeight) ||
                    (assocW == bestTruthWeight && nMeasurements > bestTruthMeas) ||
                    (assocW == bestTruthWeight && nMeasurements == bestTruthMeas &&
                     std::isfinite(chi2) && chi2 < bestTruthChi2);
                if (better) {
                    bestTruthWeight = assocW;
                    bestTruthMeas = nMeasurements;
                    bestTruthChi2 = chi2;
                    assignBest(out.truthMatched, currentTrackIndex, p, pT, deltaP, deltaPT, theta, phi,
                               nStates, nMeasurements, nOutliers, nHoles, nSharedHits, chi2, ndf,
                               assocMc, assocW, pdg, charge, momOk ? 1 : 0, pullQ, pullTh, pullPh);
                    out.truthMatched.objectIndex = objectIndex;
                    out.truthMatched.objectCollectionID = objectCollectionID;
                    out.truthMatched.seedIndex = seedIndex;
                    out.truthMatched.seedCollectionID = seedCollectionID;
                    out.truthMatched.identityValid = identityValid;
                    out.truthMatched.assocMcCollectionID = assocCol;
                }
            }
            if (nMeasurements > 0) {
                const double chi2ndf = (ndf > 0 && std::isfinite(chi2)) ? chi2 / ndf : chi2;
                const bool better =
                    (nMeasurements > bestRecoMeas) ||
                    (nMeasurements == bestRecoMeas && nHoles < bestRecoHoles) ||
                    (nMeasurements == bestRecoMeas && nHoles == bestRecoHoles &&
                     std::isfinite(chi2ndf) && chi2ndf < bestRecoChi2Ndf);
                if (better) {
                    bestRecoMeas = nMeasurements;
                    bestRecoHoles = nHoles;
                    bestRecoChi2Ndf = chi2ndf;
                    assignBest(out.recoBest, currentTrackIndex, p, pT, deltaP, deltaPT, theta, phi,
                               nStates, nMeasurements, nOutliers, nHoles, nSharedHits, chi2, ndf,
                               assocMc, assocW, pdg, charge, momOk ? 1 : 0, pullQ, pullTh, pullPh);
                    out.recoBest.objectIndex = objectIndex;
                    out.recoBest.objectCollectionID = objectCollectionID;
                    out.recoBest.seedIndex = seedIndex;
                    out.recoBest.seedCollectionID = seedCollectionID;
                    out.recoBest.identityValid = identityValid;
                    out.recoBest.assocMcCollectionID = assocCol;
                }
            }
        }
        out.hasTrack = out.p.empty() ? 0 : 1;

        auto actsGeoProvider = m_actsGeoProvider ? m_actsGeoProvider
                                                 : (m_actsGeoSvc ? m_actsGeoSvc->actsGeoProvider() : nullptr);
        const auto* actsGeoCtx = actsGeoProvider ? &actsGeoProvider->getActsGeometryContext() : nullptr;
        if (actsTracks.empty() || actsTracks.front() == nullptr ||
            actsStates.empty() || actsStates.front() == nullptr ||
            actsGeoCtx == nullptr) {
            return;
        }
        Acts::TrackContainer<Acts::ConstVectorTrackContainer,
                             Acts::ConstVectorMultiTrajectory,
                             Acts::detail::ConstRefHolder>
            trackContainer(*actsTracks.front(), *actsStates.front());
        const auto nActsTracks = static_cast<int>(trackContainer.size());
        Acts::ConstProxyAccessor<unsigned int> seedNumber("seed");
        for (int actsTrackIndex = 0; actsTrackIndex < nActsTracks; ++actsTrackIndex) {
            const auto track = trackContainer.getTrack(actsTrackIndex);
            int parentSeedIndex = -1;
            std::uint32_t parentSeedCollectionID = 0;
            try {
                const auto seedPosition = static_cast<std::size_t>(seedNumber(track));
                if (seedPosition < chainSeeds.size() && chainSeeds[seedPosition] != nullptr) {
                    const auto seedId = chainSeeds[seedPosition]->id();
                    parentSeedIndex = seedId.index;
                    parentSeedCollectionID = seedId.collectionID;
                }
            } catch (...) {}
            int parentTrackIndex = -1;
            std::uint32_t parentTrackCollectionID = 0;
            int parentIdentityValid = 0;
            if (const auto found = trackObjectsBySeed.find({parentSeedCollectionID, parentSeedIndex});
                parentSeedIndex >= 0 && found != trackObjectsBySeed.end() &&
                found->second.size() == 1) {
                parentTrackCollectionID = found->second.front().first;
                parentTrackIndex = found->second.front().second;
                parentIdentityValid = 1;
            }
            const std::size_t stateBegin = out.state_index.size();
            int stateIndex = 0;
            for (const auto& state : track.trackStatesReversed()) {
                if (!state.hasReferenceSurface()) {
                    continue;
                }

                const int estimateKind = b0trk::preferredEstimateKind(
                    state.hasPredicted(), state.hasFiltered(), state.hasSmoothed());
                const auto selectParams = [&state](auto&& fill) {
                    if (state.hasSmoothed()) {
                        fill(state.smoothed());
                    } else if (state.hasFiltered()) {
                        fill(state.filtered());
                    } else if (state.hasPredicted()) {
                        fill(state.predicted());
                    }
                };

                double loc0 = nan;
                double loc1 = nan;
                double stateTheta = nan;
                double statePhi = nan;
                double stateQOverP = nan;
                double stateTime = nan;
                selectParams([&](const auto& params) {
                    loc0 = params[Acts::eBoundLoc0];
                    loc1 = params[Acts::eBoundLoc1];
                    stateTheta = params[Acts::eBoundTheta];
                    statePhi = params[Acts::eBoundPhi];
                    stateQOverP = params[Acts::eBoundQOverP];
                    stateTime = b0trk::nativeTimeToNs(params[Acts::eBoundTime], Acts::UnitConstants::ns);
                });

                double predLoc0 = nan, predLoc1 = nan, predTheta = nan, predPhi = nan;
                double predQOverP = nan, predTime = nan;
                std::array<double, 21> predCovUpper{};
                predCovUpper.fill(nan);
                int measDim = 0, proj0 = -1, proj1 = -1;
                double measCov00 = nan, measCov01 = nan, measCov11 = nan;
                double innov0 = nan, innov1 = nan, innovPull0 = nan, innovPull1 = nan, innovChi2 = nan;
                if (state.hasPredicted()) {
                    const auto pred = state.predicted();
                    const auto predCov = state.predictedCovariance();
                    predLoc0 = pred[Acts::eBoundLoc0]; predLoc1 = pred[Acts::eBoundLoc1];
                    predPhi = pred[Acts::eBoundPhi]; predTheta = pred[Acts::eBoundTheta];
                    predQOverP = pred[Acts::eBoundQOverP];
                    predTime = b0trk::nativeTimeToNs(pred[Acts::eBoundTime], Acts::UnitConstants::ns);
                    const std::array<double, 6> scale{
                        1.0 / Acts::UnitConstants::mm, 1.0 / Acts::UnitConstants::mm,
                        1.0 / Acts::UnitConstants::rad, 1.0 / Acts::UnitConstants::rad,
                        Acts::UnitConstants::GeV, 1.0 / Acts::UnitConstants::ns};
                    int k = 0;
                    for (unsigned i = 0; i < 6; ++i) {
                        for (unsigned j = i; j < 6; ++j) predCovUpper[k++] = predCov(i, j) * scale[i] * scale[j];
                    }
                    if (state.hasCalibrated()) {
                        try {
                            const auto meas = state.effectiveCalibrated();
                            const auto measCov = state.effectiveCalibratedCovariance();
                            const auto projector = state.projectorSubspaceIndices();
                            measDim = static_cast<int>(state.calibratedSize());
                            if (measDim > 0) {
                                proj0 = static_cast<int>(projector[0]);
                                measCov00 = measCov(0, 0);
                                innov0 = meas[0] - pred[projector[0]];
                            }
                            if (measDim > 1) {
                                proj1 = static_cast<int>(projector[1]);
                                measCov01 = measCov(0, 1);
                                measCov11 = measCov(1, 1);
                                innov1 = meas[1] - pred[projector[1]];
                            }
                            if (measDim == 1) {
                                const double s00 = measCov(0, 0) + predCov(projector[0], projector[0]);
                                if (s00 > 0.0 && std::isfinite(s00)) {
                                    innovPull0 = innov0 / std::sqrt(s00);
                                    innovChi2 = innov0 * innov0 / s00;
                                }
                            } else if (measDim == 2) {
                                Eigen::Matrix2d S;
                                S(0, 0) = measCov(0, 0) + predCov(projector[0], projector[0]);
                                S(0, 1) = measCov(0, 1) + predCov(projector[0], projector[1]);
                                S(1, 0) = measCov(1, 0) + predCov(projector[1], projector[0]);
                                S(1, 1) = measCov(1, 1) + predCov(projector[1], projector[1]);
                                if (S.allFinite() && S(0, 0) > 0.0 && S(1, 1) > 0.0) {
                                    innovPull0 = innov0 / std::sqrt(S(0, 0));
                                    innovPull1 = innov1 / std::sqrt(S(1, 1));
                                    const Eigen::Vector2d r(innov0, innov1);
                                    const auto ldlt = S.ldlt();
                                    if (ldlt.info() == Eigen::Success) innovChi2 = r.dot(ldlt.solve(r));
                                }
                            }
                        } catch (...) {}
                    }
                }

                const auto flags   = state.typeFlags();
                const int typeMask = trackStateTypeMask(flags);
                const bool physicsState =
                    (typeMask & (kStateMeasurement | kStateOutlier | kStateHole)) != 0;

                double globalX = nan;
                double globalY = nan;
                double globalZ = nan;
                if (std::isfinite(loc0) && std::isfinite(loc1) &&
                    std::isfinite(stateTheta) && std::isfinite(statePhi)) {
                    const auto& surface = state.referenceSurface();
                    const Acts::Vector2 local(loc0, loc1);
                    const Acts::Vector3 direction(std::sin(stateTheta) * std::cos(statePhi),
                                                  std::sin(stateTheta) * std::sin(statePhi),
                                                  std::cos(stateTheta));
                    const auto global = surface.localToGlobal(*actsGeoCtx, local, direction);
                    globalX = global.x();
                    globalY = global.y();
                    globalZ = global.z();
                }

                const SensorRef* sensorRef = nullptr;
                int mapping = b0trk::kMapUnresolved;
                const auto surfIt =
                    m_surfaceToSensorIdx.find(state.referenceSurface().geometryId().value());
                if (surfIt != m_surfaceToSensorIdx.end()) {
                    sensorRef = &m_sensorRefs[surfIt->second];
                    mapping = b0trk::kMapExact;
                    ++out.nMapExact;
                } else if (physicsState) {
                    const auto hit = closestSensor(globalX, globalY, globalZ);
                    if (hit.ref != nullptr && hit.outsideX == 0.0 && hit.outsideY == 0.0 &&
                        std::abs(hit.localZ) <= m_fallbackMaxNormalMm) {
                        sensorRef = hit.ref;
                        mapping = b0trk::kMapFallback;
                        ++out.nMapFallback;
                    } else {
                        ++out.nMapFailedPhysics;
                    }
                } else {
                    ++out.nMapUnmappedNonPhysics;
                }

                PixelSnap trackPixel{};
                if (sensorRef != nullptr) {
                    trackPixel = snapToAclgadPixel(globalX, globalY, globalZ, sensorRef->cellID,
                                                   sensorRef->detElement);
                    if (!trackPixel.ok) {
                        ++local.m_nPixelSnapFailed;
                    }
                }

                double measLoc0 = nan;
                double measLoc1 = nan;
                double resid0 = nan;
                double resid1 = nan;
                if (state.hasCalibrated()) {
                    try {
                        const auto meas = state.effectiveCalibrated();
                        if (meas.size() >= 2) {
                            measLoc0 = meas[0];
                            measLoc1 = meas[1];
                            if (std::isfinite(loc0) && std::isfinite(measLoc0)) {
                                resid0 = loc0 - measLoc0;
                            }
                            if (std::isfinite(loc1) && std::isfinite(measLoc1)) {
                                resid1 = loc1 - measLoc1;
                            }
                        }
                    } catch (...) {
                    }
                }

                out.state_track_index.push_back(actsTrackIndex);
                out.state_parent_seed_index.push_back(parentSeedIndex);
                out.state_parent_seed_collectionID.push_back(parentSeedCollectionID);
                out.state_parent_track_index.push_back(parentTrackIndex);
                out.state_parent_track_collectionID.push_back(parentTrackCollectionID);
                out.state_parent_identity_valid.push_back(parentIdentityValid);
                out.state_estimate_kind.push_back(estimateKind);
                out.state_index.push_back(stateIndex++);
                out.state_acts_index.push_back(static_cast<int>(state.index()));
                out.state_type.push_back(typeMask);
                out.state_mapping_method.push_back(mapping);
                out.state_surface.push_back(state.referenceSurface().geometryId().value());
                out.state_loc0.push_back(loc0);
                out.state_loc1.push_back(loc1);
                out.state_meas_loc0.push_back(measLoc0);
                out.state_meas_loc1.push_back(measLoc1);
                out.state_resid_loc0.push_back(resid0);
                out.state_resid_loc1.push_back(resid1);
                out.state_meas_dim.push_back(measDim);
                out.state_proj_index0.push_back(proj0); out.state_proj_index1.push_back(proj1);
                out.state_meas_cov00.push_back(measCov00); out.state_meas_cov01.push_back(measCov01);
                out.state_meas_cov11.push_back(measCov11);
                out.state_pred_loc0.push_back(predLoc0); out.state_pred_loc1.push_back(predLoc1);
                out.state_pred_theta.push_back(predTheta); out.state_pred_phi.push_back(predPhi);
                out.state_pred_qOverP.push_back(predQOverP); out.state_pred_time.push_back(predTime);
                for (std::size_t k = 0; k < predCovUpper.size(); ++k) out.state_pred_cov_upper[k].push_back(predCovUpper[k]);
                out.state_innov0.push_back(innov0); out.state_innov1.push_back(innov1);
                out.state_innov_pull0.push_back(innovPull0); out.state_innov_pull1.push_back(innovPull1);
                out.state_innov_chi2.push_back(innovChi2);
                out.x_on_plane.push_back(globalX);
                out.y_on_plane.push_back(globalY);
                out.z_on_plane.push_back(globalZ);
                out.aclgad_xPix.push_back(trackPixel.x);
                out.aclgad_yPix.push_back(trackPixel.y);
                out.aclgad_zPix.push_back(trackPixel.z);
                out.aclgad_dx.push_back(trackPixel.dx);
                out.aclgad_dy.push_back(trackPixel.dy);
                out.aclgad_dz.push_back(trackPixel.dz);
                out.aclgad_pixX.push_back(trackPixel.pixX);
                out.aclgad_pixY.push_back(trackPixel.pixY);
                out.aclgad_pixZ.push_back(trackPixel.pixZ);
                out.aclgad_plane.push_back(sensorRef ? sensorRef->plane : -1);
                out.aclgad_station.push_back(sensorRef ? sensorRef->station : -1);
                out.aclgad_module.push_back(sensorRef ? sensorRef->module : -1);
                out.aclgad_side.push_back(sensorRef ? sensorRef->side : -1);
                out.aclgad_sensor.push_back(sensorRef ? sensorRef->sensor : -1);
                out.aclgad_cellID.push_back(trackPixel.cellID);
                out.state_theta.push_back(stateTheta);
                out.state_phi.push_back(statePhi);
                out.state_qOverP.push_back(stateQOverP);
                out.state_time.push_back(stateTime);
                int statePdg = 0;
                if (parentIdentityValid) {
                    if (const auto found = pdgByTrackObject.find({parentTrackCollectionID, parentTrackIndex});
                        found != pdgByTrackObject.end()) statePdg = found->second;
                }
                out.state_pdg.push_back(statePdg);
            }
            const int nPushed = stateIndex;
            for (int i = 0; i < nPushed; ++i) {
                out.state_index[stateBegin + static_cast<std::size_t>(i)] = nPushed - 1 - i;
            }
        }
    };

    fillChain(local.m_ts, tsTrajectories, tsTracks, tsEdmTracks, tsAssocs, tsActsTracks,
              tsActsTrackStates, truthSeeds);
    fillChain(local.m_ckf, ckfTrajectories, ckfTracks, ckfEdmTracks, ckfAssocs, ckfActsTracks,
              ckfActsTrackStates, stubSeeds);
    local.m_selPrimaryHasTruthMatchedTrack = local.m_ckf.truthMatched.index >= 0 ? 1 : 0;
    local.m_nSensorMapExact = local.m_ts.nMapExact + local.m_ckf.nMapExact;
    local.m_nSensorMapFallback = local.m_ts.nMapFallback + local.m_ckf.nMapFallback;
    local.m_nSensorMapFailed = local.m_ts.nMapFailedPhysics + local.m_ckf.nMapFailedPhysics;

    {
        std::lock_guard<std::mutex> lock(m_fillMutex);
        using std::swap;
        swap(this->m_selPrimaryMeasurementReconstructable, local.m_selPrimaryMeasurementReconstructable);
        swap(this->m_nMeasurementStationsSelectedPrimary, local.m_nMeasurementStationsSelectedPrimary);
        swap(this->vm_truth_seed_made_unfiltered_track, local.vm_truth_seed_made_unfiltered_track);
        swap(this->vm_truth_seed_n_unfiltered_tracks, local.vm_truth_seed_n_unfiltered_tracks);
        swap(this->m_selPrimaryHasTruthMatchedTrack, local.m_selPrimaryHasTruthMatchedTrack);
        swap(this->vm_truth_seed_survived_ambiguity, local.vm_truth_seed_survived_ambiguity);
        swap(this->vm_truth_seed_momentum_resolved, local.vm_truth_seed_momentum_resolved);
        swap(this->vm_truth_seed_n_filtered_tracks, local.vm_truth_seed_n_filtered_tracks);
        swap(this->m_hasCkfTrajectoriesUnfiltered, local.m_hasCkfTrajectoriesUnfiltered);
        swap(this->m_nSelectedPrimaryMeasurements, local.m_nSelectedPrimaryMeasurements);
        swap(this->m_selPrimaryHasUnfilteredTrack, local.m_selPrimaryHasUnfilteredTrack);
        swap(this->m_hasTsTrajectoriesUnfiltered, local.m_hasTsTrajectoriesUnfiltered);
        swap(this->vm_seed_made_unfiltered_track, local.vm_seed_made_unfiltered_track);
        swap(this->vm_seed_assoc_mcCollectionID, local.vm_seed_assoc_mcCollectionID);
        swap(this->m_selPrimaryHasFilteredTrack, local.m_selPrimaryHasFilteredTrack);
        swap(this->vm_seed_n_unfiltered_tracks, local.vm_seed_n_unfiltered_tracks);
        swap(this->m_selPrimaryMcCollectionID, local.m_selPrimaryMcCollectionID);
        swap(this->vm_truth_seed_sigma_qOverP, local.vm_truth_seed_sigma_qOverP);
        swap(this->vm_truth_seed_became_track, local.vm_truth_seed_became_track);
        swap(this->vm_seed_survived_ambiguity, local.vm_seed_survived_ambiguity);
        swap(this->m_nSimHitsUnresolvedCellID, local.m_nSimHitsUnresolvedCellID);
        swap(this->vm_truth_seed_sigma_theta, local.vm_truth_seed_sigma_theta);
        swap(this->vm_seed_momentum_resolved, local.vm_seed_momentum_resolved);
        swap(this->vm_seed_n_filtered_tracks, local.vm_seed_n_filtered_tracks);
        swap(this->m_hasCkfTracksUnfiltered, local.m_hasCkfTracksUnfiltered);
        swap(this->m_hasCkfAssocsUnfiltered, local.m_hasCkfAssocsUnfiltered);
        swap(this->m_primaryMcCollectionID, local.m_primaryMcCollectionID);
        swap(this->vm_truth_seed_sigma_phi, local.vm_truth_seed_sigma_phi);
        swap(this->m_hasTsTracksUnfiltered, local.m_hasTsTracksUnfiltered);
        swap(this->m_selPrimaryThscatMrad, local.m_selPrimaryThscatMrad);
        swap(this->vm_truth_seed_quality, local.vm_truth_seed_quality);
        swap(this->vm_raw_mcCollectionID, local.vm_raw_mcCollectionID);
        swap(this->vm_rec_mcCollectionID, local.vm_rec_mcCollectionID);
        swap(this->vm_seed_assoc_mcIndex, local.vm_seed_assoc_mcIndex);
        swap(this->m_nMissingMcRelation, local.m_nMissingMcRelation);
        swap(this->vm_seed_assoc_weight, local.vm_seed_assoc_weight);
        swap(this->vm_seed_became_track, local.vm_seed_became_track);
        swap(this->m_nSensorMapFallback, local.m_nSensorMapFallback);
        swap(this->vm_truth_seed_qOverP, local.vm_truth_seed_qOverP);
        swap(this->vm_truth_seed_charge, local.vm_truth_seed_charge);
        swap(this->vm_seed_sigma_qOverP, local.vm_seed_sigma_qOverP);
        swap(this->m_hasCkfTrajectories, local.m_hasCkfTrajectories);
        swap(this->m_hasCkfTrackParams, local.m_hasCkfTrackParams);
        swap(this->vm_truth_seed_theta, local.vm_truth_seed_theta);
        swap(this->m_hasTsTrajectories, local.m_hasTsTrajectories);
        swap(this->m_selPrimaryHasSeed, local.m_selPrimaryHasSeed);
        swap(this->vm_raw_dominantFrac, local.vm_raw_dominantFrac);
        swap(this->vm_truth_seed_nHits, local.vm_truth_seed_nHits);
        swap(this->vm_rec_dominantFrac, local.vm_rec_dominantFrac);
        swap(this->vm_seed_sigma_theta, local.vm_seed_sigma_theta);
        swap(this->vm_mcCollectionIDEE, local.vm_mcCollectionIDEE);
        swap(this->m_selPrimaryMcIndex, local.m_selPrimaryMcIndex);
        swap(this->m_hasCkfActsTracks, local.m_hasCkfActsTracks);
        swap(this->vm_rec_nContribSim, local.vm_rec_nContribSim);
        swap(this->m_nSensorMapFailed, local.m_nSensorMapFailed);
        swap(this->vm_truth_seed_loc1, local.vm_truth_seed_loc1);
        swap(this->m_hasTsTrackParams, local.m_hasTsTrackParams);
        swap(this->vm_truth_seed_loc0, local.vm_truth_seed_loc0);
        swap(this->vm_mcCollectionIDP, local.vm_mcCollectionIDP);
        swap(this->m_primaryStatusOut, local.m_primaryStatusOut);
        swap(this->m_nStationsPrimary, local.m_nStationsPrimary);
        swap(this->m_selPrimaryCharge, local.m_selPrimaryCharge);
        swap(this->m_hasCkfActsStates, local.m_hasCkfActsStates);
        swap(this->vm_raw_nContribSim, local.vm_raw_nContribSim);
        swap(this->m_nPixelSnapFailed, local.m_nPixelSnapFailed);
        swap(this->vm_rec_nContribMc, local.vm_rec_nContribMc);
        swap(this->vm_isSelPrimaryEE, local.vm_isSelPrimaryEE);
        swap(this->vm_mcCollectionID, local.vm_mcCollectionID);
        swap(this->vm_raw_nContribMc, local.vm_raw_nContribMc);
        swap(this->m_hasTsActsStates, local.m_hasTsActsStates);
        swap(this->m_hasTsActsTracks, local.m_hasTsActsTracks);
        swap(this->vm_truth_seed_phi, local.vm_truth_seed_phi);
        swap(this->m_nSensorMapExact, local.m_nSensorMapExact);
        swap(this->vm_seed_sigma_phi, local.vm_seed_sigma_phi);
        swap(this->vm_isSelPrimaryP, local.vm_isSelPrimaryP);
        swap(this->m_primaryMcIndex, local.m_primaryMcIndex);
        swap(this->vm_rec_mixedCell, local.vm_rec_mixedCell);
        swap(this->vm_rec_totalEdep, local.vm_rec_totalEdep);
        swap(this->vm_raw_mixedCell, local.vm_raw_mixedCell);
        swap(this->vm_raw_totalEdep, local.vm_raw_totalEdep);
        swap(this->vm_raw_timeStamp, local.vm_raw_timeStamp);
        swap(this->m_nCkfUnfiltered, local.m_nCkfUnfiltered);
        swap(this->vm_aclgad_yPixT, local.vm_aclgad_yPixT);
        swap(this->m_nTsUnfiltered, local.m_nTsUnfiltered);
        swap(this->vm_aclgad_pixXR, local.vm_aclgad_pixXR);
        swap(this->vm_truth_seed_p, local.vm_truth_seed_p);
        swap(this->vm_aclgad_pixXT, local.vm_aclgad_pixXT);
        swap(this->vm_aclgad_zPixR, local.vm_aclgad_zPixR);
        swap(this->vm_aclgad_zPixT, local.vm_aclgad_zPixT);
        swap(this->vm_isSelPrimary, local.vm_isSelPrimary);
        swap(this->vm_aclgad_pixYR, local.vm_aclgad_pixYR);
        swap(this->vm_seed_quality, local.vm_seed_quality);
        swap(this->m_nMeasurements, local.m_nMeasurements);
        swap(this->vm_rec_time_err, local.vm_rec_time_err);
        swap(this->vm_rec_edep_err, local.vm_rec_edep_err);
        swap(this->m_hasTruthSeeds, local.m_hasTruthSeeds);
        swap(this->vm_aclgad_yPixR, local.vm_aclgad_yPixR);
        swap(this->m_primaryCharge, local.m_primaryCharge);
        swap(this->m_primaryPdgOut, local.m_primaryPdgOut);
        swap(this->vm_aclgad_pixYT, local.vm_aclgad_pixYT);
        swap(this->vm_aclgad_xPixR, local.vm_aclgad_xPixR);
        swap(this->vm_aclgad_xPixT, local.vm_aclgad_xPixT);
        swap(this->vm_isPrimaryEE, local.vm_isPrimaryEE);
        swap(this->vm_raw_station, local.vm_raw_station);
        swap(this->vm_cellIDEntry, local.vm_cellIDEntry);
        swap(this->vm_seed_charge, local.vm_seed_charge);
        swap(this->m_hasStubSeeds, local.m_hasStubSeeds);
        swap(this->vm_seed_qOverP, local.vm_seed_qOverP);
        swap(this->vm_rec_station, local.vm_rec_station);
        swap(this->vm_rec_mcIndex, local.vm_rec_mcIndex);
        swap(this->m_selPrimaryPT, local.m_selPrimaryPT);
        swap(this->m_selPrimaryPz, local.m_selPrimaryPz);
        swap(this->m_hasCkfAssocs, local.m_hasCkfAssocs);
        swap(this->vm_raw_mcIndex, local.vm_raw_mcIndex);
        swap(this->m_selPrimaryPx, local.m_selPrimaryPx);
        swap(this->m_hasRawAssocs, local.m_hasRawAssocs);
        swap(this->m_selPrimaryPy, local.m_selPrimaryPy);
        swap(this->m_nCkfFiltered, local.m_nCkfFiltered);
        swap(this->vm_sensorEntry, local.vm_sensorEntry);
        swap(this->m_hasCkfTracks, local.m_hasCkfTracks);
        swap(this->m_selPrimaryP, local.m_selPrimaryP);
        swap(this->vm_sensorExit, local.vm_sensorExit);
        swap(this->vm_raw_charge, local.vm_raw_charge);
        swap(this->vm_aclgad_dyT, local.vm_aclgad_dyT);
        swap(this->vm_isPrimaryP, local.vm_isPrimaryP);
        swap(this->vm_rec_cellID, local.vm_rec_cellID);
        swap(this->m_hasTsTracks, local.m_hasTsTracks);
        swap(this->vm_aclgad_dxR, local.vm_aclgad_dxR);
        swap(this->vm_raw_module, local.vm_raw_module);
        swap(this->vm_rec_module, local.vm_rec_module);
        swap(this->m_nTruthSeeds, local.m_nTruthSeeds);
        swap(this->vm_raw_cellID, local.vm_raw_cellID);
        swap(this->vm_cellIDExit, local.vm_cellIDExit);
        swap(this->vm_seed_theta, local.vm_seed_theta);
        swap(this->vm_aclgad_dxT, local.vm_aclgad_dxT);
        swap(this->m_nTsFiltered, local.m_nTsFiltered);
        swap(this->vm_rec_sensor, local.vm_rec_sensor);
        swap(this->m_eventNumber, local.m_eventNumber);
        swap(this->vm_aclgad_dzR, local.vm_aclgad_dzR);
        swap(this->vm_aclgad_dzT, local.vm_aclgad_dzT);
        swap(this->vm_aclgad_dyR, local.vm_aclgad_dyR);
        swap(this->vm_raw_sensor, local.vm_raw_sensor);
        swap(this->vm_seed_nHits, local.vm_seed_nHits);
        swap(this->m_hasTsAssocs, local.m_hasTsAssocs);
        swap(this->vm_seed_loc0, local.vm_seed_loc0);
        swap(this->vm_timeEntry, local.vm_timeEntry);
        swap(this->vm_rec_plane, local.vm_rec_plane);
        swap(this->m_genBeamPPy, local.m_genBeamPPy);
        swap(this->vm_stationEE, local.vm_stationEE);
        swap(this->vm_rec_covxx, local.vm_rec_covxx);
        swap(this->m_nStubSeeds, local.m_nStubSeeds);
        swap(this->vm_cellFired, local.vm_cellFired);
        swap(this->vm_seed_loc1, local.vm_seed_loc1);
        swap(this->vm_mcIndexEE, local.vm_mcIndexEE);
        swap(this->vm_rec_covzz, local.vm_rec_covzz);
        swap(this->m_genBeamPPz, local.m_genBeamPPz);
        swap(this->vm_raw_plane, local.vm_raw_plane);
        swap(this->m_genBeamPPx, local.m_genBeamPPx);
        swap(this->vm_isPrimary, local.vm_isPrimary);
        swap(this->m_genBeamPPT, local.m_genBeamPPT);
        swap(this->vm_rec_covyy, local.vm_rec_covyy);
        swap(this->m_genBeamPz, local.m_genBeamPz);
        swap(this->vm_seed_phi, local.vm_seed_phi);
        swap(this->vm_rec_side, local.vm_rec_side);
        swap(this->vm_moduleEE, local.vm_moduleEE);
        swap(this->vm_mcIndexP, local.vm_mcIndexP);
        swap(this->m_primaryPy, local.m_primaryPy);
        swap(this->m_primaryPz, local.m_primaryPz);
        swap(this->m_primaryPT, local.m_primaryPT);
        swap(this->m_genBeamPx, local.m_genBeamPx);
        swap(this->vm_rec_pixY, local.vm_rec_pixY);
        swap(this->m_primaryPx, local.m_primaryPx);
        swap(this->vm_rec_edep, local.vm_rec_edep);
        swap(this->m_genBeamPT, local.m_genBeamPT);
        swap(this->vm_rec_pixZ, local.vm_rec_pixZ);
        swap(this->m_genBeamPy, local.m_genBeamPy);
        swap(this->vm_rec_time, local.vm_rec_time);
        swap(this->m_genBeamPP, local.m_genBeamPP);
        swap(this->vm_raw_side, local.vm_raw_side);
        swap(this->vm_nStepsEE, local.vm_nStepsEE);
        swap(this->vm_timeExit, local.vm_timeExit);
        swap(this->vm_stationP, local.vm_stationP);
        swap(this->vm_rec_pixX, local.vm_rec_pixX);
        swap(this->vm_statusEE, local.vm_statusEE);
        swap(this->vm_mcIndex, local.vm_mcIndex);
        swap(this->vm_pyEntry, local.vm_pyEntry);
        swap(this->m_nRawHits, local.m_nRawHits);
        swap(this->m_nSimHits, local.m_nSimHits);
        swap(this->vm_pTEntry, local.vm_pTEntry);
        swap(this->m_genBeamP, local.m_genBeamP);
        swap(this->m_nRecHits, local.m_nRecHits);
        swap(this->vm_pxEntry, local.vm_pxEntry);
        swap(this->vm_moduleP, local.vm_moduleP);
        swap(this->vm_pzEntry, local.vm_pzEntry);
        swap(this->m_primaryP, local.m_primaryP);
        swap(this->vm_station, local.vm_station);
        swap(this->vm_statusP, local.vm_statusP);
        swap(this->vm_planeEE, local.vm_planeEE);
        swap(this->vm_sensorP, local.vm_sensorP);
        swap(this->vm_status, local.vm_status);
        swap(this->vm_pTExit, local.vm_pTExit);
        swap(this->vm_sideEE, local.vm_sideEE);
        swap(this->vm_seed_p, local.vm_seed_p);
        swap(this->vm_xEntry, local.vm_xEntry);
        swap(this->vm_zEntry, local.vm_zEntry);
        swap(this->vm_yEntry, local.vm_yEntry);
        swap(this->vm_planeP, local.vm_planeP);
        swap(this->vm_pzExit, local.vm_pzExit);
        swap(this->vm_pxExit, local.vm_pxExit);
        swap(this->vm_sensor, local.vm_sensor);
        swap(this->vm_module, local.vm_module);
        swap(this->vm_pEntry, local.vm_pEntry);
        swap(this->vm_pyExit, local.vm_pyExit);
        swap(this->vm_cellID, local.vm_cellID);
        swap(this->vm_yExit, local.vm_yExit);
        swap(this->vm_rec_x, local.vm_rec_x);
        swap(this->vm_pdgEE, local.vm_pdgEE);
        swap(this->vm_pathP, local.vm_pathP);
        swap(this->vm_pExit, local.vm_pExit);
        swap(this->vm_rec_z, local.vm_rec_z);
        swap(this->vm_timeP, local.vm_timeP);
        swap(this->m_genPpy, local.m_genPpy);
        swap(this->m_genPpx, local.m_genPpx);
        swap(this->vm_rec_y, local.vm_rec_y);
        swap(this->vm_sideP, local.vm_sideP);
        swap(this->vm_plane, local.vm_plane);
        swap(this->m_genPpz, local.m_genPpz);
        swap(this->m_genPpT, local.m_genPpT);
        swap(this->vm_xExit, local.vm_xExit);
        swap(this->vm_zExit, local.vm_zExit);
        swap(this->beam_pdg, local.beam_pdg);
        swap(this->vm_pixX, local.vm_pixX);
        swap(this->vm_path, local.vm_path);
        swap(this->vm_pdgP, local.vm_pdgP);
        swap(this->vm_detZ, local.vm_detZ);
        swap(this->beam_pz, local.beam_pz);
        swap(this->vm_time, local.vm_time);
        swap(this->vm_pixZ, local.vm_pixZ);
        swap(this->beam_pT, local.beam_pT);
        swap(this->m_genPp, local.m_genPp);
        swap(this->beam_py, local.beam_py);
        swap(this->vm_pixY, local.vm_pixY);
        swap(this->beam_px, local.beam_px);
        swap(this->vm_detY, local.vm_detY);
        swap(this->vm_side, local.vm_side);
        swap(this->vm_eDep, local.vm_eDep);
        swap(this->vm_detX, local.vm_detX);
        swap(this->beam_p, local.beam_p);
        swap(this->vm_pyP, local.vm_pyP);
        swap(this->vm_pTP, local.vm_pTP);
        swap(this->vm_pzP, local.vm_pzP);
        swap(this->vm_pdg, local.vm_pdg);
        swap(this->vm_pxP, local.vm_pxP);
        swap(this->m_ckf, local.m_ckf);
        swap(this->vm_zR, local.vm_zR);
        swap(this->vm_zT, local.vm_zT);
        swap(this->vm_py, local.vm_py);
        swap(this->vm_px, local.vm_px);
        swap(this->vm_pz, local.vm_pz);
        swap(this->vm_yT, local.vm_yT);
        swap(this->vm_yP, local.vm_yP);
        swap(this->vm_pP, local.vm_pP);
        swap(this->vm_xP, local.vm_xP);
        swap(this->vm_pT, local.vm_pT);
        swap(this->vm_yR, local.vm_yR);
        swap(this->vm_xT, local.vm_xT);
        swap(this->vm_xR, local.vm_xR);
        swap(this->vm_zP, local.vm_zP);
        swap(this->vm_p, local.vm_p);
        swap(this->m_ts, local.m_ts);
        m_tree->Fill();
    }
}

void B0Trackers::BestSel::reset(double nan) {
    index = -1;
    nStates = -1;
    nMeasurements = -1;
    nOutliers = -1;
    nHoles = -1;
    nSharedHits = -1;
    ndf = -1;
    assocMcIndex = -1;
    assocMcCollectionID = 0;
    objectIndex = -1; objectCollectionID = 0; seedIndex = -1; seedCollectionID = 0; identityValid = 0;
    pdg = 0;
    charge = 0;
    momentumResolved = 0;
    p = nan;
    pT = nan;
    deltaP = nan;
    deltaPT = nan;
    theta = nan;
    phi = nan;
    chi2 = nan;
    assocWeight = nan;
    pullQOverP = nan;
    pullTheta = nan;
    pullPhi = nan;
}

void B0Trackers::TrackChain::clear(double nan) {
    p.clear(); pT.clear(); delta_p.clear(); delta_pT.clear();
    px.clear(); py.clear(); pz.clear();
    theta.clear(); phi.clear();
    qOverP.clear(); time.clear(); momentum_resolved.clear();
    loc0.clear(); loc1.clear();
    sigma_loc0.clear(); sigma_loc1.clear(); sigma_phi.clear();
    sigma_theta.clear(); sigma_qOverP.clear(); sigma_time.clear();
    pull_qOverP.clear(); pull_theta.clear(); pull_phi.clear();
    chi2.clear(); ndf.clear();
    index.clear(); charge.clear(); type.clear(); pdg.clear();
    object_index.clear(); object_collectionID.clear(); seed_index.clear(); seed_collectionID.clear();
    identity_valid.clear();
    for (auto& v : cov_upper) v.clear();
    nStates.clear(); nMeasurements.clear(); nOutliers.clear(); nHoles.clear(); nSharedHits.clear();
    assoc_mcIndex.clear(); assoc_mcCollectionID.clear(); assoc_weight.clear();
    state_track_index.clear(); state_index.clear(); state_acts_index.clear();
    state_parent_seed_index.clear(); state_parent_seed_collectionID.clear();
    state_parent_track_index.clear();
    state_parent_track_collectionID.clear(); state_parent_identity_valid.clear();
    state_estimate_kind.clear();
    state_type.clear(); state_pdg.clear(); state_mapping_method.clear();
    state_surface.clear();
    state_loc0.clear(); state_loc1.clear();
    state_meas_loc0.clear(); state_meas_loc1.clear();
    state_resid_loc0.clear(); state_resid_loc1.clear();
    state_meas_dim.clear(); state_proj_index0.clear(); state_proj_index1.clear();
    state_meas_cov00.clear(); state_meas_cov01.clear(); state_meas_cov11.clear();
    state_pred_loc0.clear(); state_pred_loc1.clear(); state_pred_theta.clear(); state_pred_phi.clear();
    state_pred_qOverP.clear(); state_pred_time.clear();
    for (auto& v : state_pred_cov_upper) v.clear();
    state_innov0.clear(); state_innov1.clear(); state_innov_pull0.clear(); state_innov_pull1.clear();
    state_innov_chi2.clear();
    x_on_plane.clear(); y_on_plane.clear(); z_on_plane.clear();
    aclgad_xPix.clear(); aclgad_yPix.clear(); aclgad_zPix.clear();
    aclgad_dx.clear(); aclgad_dy.clear(); aclgad_dz.clear();
    aclgad_pixX.clear(); aclgad_pixY.clear(); aclgad_pixZ.clear();
    aclgad_plane.clear(); aclgad_module.clear(); aclgad_side.clear();
    aclgad_sensor.clear(); aclgad_station.clear();
    aclgad_cellID.clear();
    state_theta.clear(); state_phi.clear(); state_qOverP.clear(); state_time.clear();
    oracle.reset(nan);
    truthMatched.reset(nan);
    recoBest.reset(nan);
    nMapExact = 0; nMapFallback = 0; nMapFailedPhysics = 0; nMapUnmappedNonPhysics = 0;
    hasTrack = 0;
}

void B0Trackers::bindBestSel(const std::string& prefix, BestSel& b) {
    auto br = [this](const std::string& name, auto* ptr) {
        m_tree->Branch(name.c_str(), ptr);
    };
    br(prefix + "index", &b.index);
    br(prefix + "p", &b.p);
    br(prefix + "pT", &b.pT);
    br(prefix + "delta_p", &b.deltaP);
    br(prefix + "delta_pT", &b.deltaPT);
    br(prefix + "theta", &b.theta);
    br(prefix + "phi", &b.phi);
    br(prefix + "nStates", &b.nStates);
    br(prefix + "nMeasurements", &b.nMeasurements);
    br(prefix + "nOutliers", &b.nOutliers);
    br(prefix + "nHoles", &b.nHoles);
    br(prefix + "nSharedHits", &b.nSharedHits);
    br(prefix + "chi2", &b.chi2);
    br(prefix + "ndf", &b.ndf);
    br(prefix + "assoc_mcIndex", &b.assocMcIndex);
    br(prefix + "assoc_mcCollectionID", &b.assocMcCollectionID);
    br(prefix + "object_index", &b.objectIndex);
    br(prefix + "object_collectionID", &b.objectCollectionID);
    br(prefix + "seed_index", &b.seedIndex);
    br(prefix + "seed_collectionID", &b.seedCollectionID);
    br(prefix + "identity_valid", &b.identityValid);
    br(prefix + "assoc_weight", &b.assocWeight);
    br(prefix + "pdg", &b.pdg);
    br(prefix + "charge", &b.charge);
    br(prefix + "momentum_resolved", &b.momentumResolved);
    br(prefix + "pull_qOverP", &b.pullQOverP);
    br(prefix + "pull_theta", &b.pullTheta);
    br(prefix + "pull_phi", &b.pullPhi);
}

void B0Trackers::bindTrackChain(const std::string& trkPrefix, TrackChain& c) {
    auto br = [this](const std::string& name, auto* ptr) {
        m_tree->Branch(name.c_str(), ptr);
    };
    br(trkPrefix + "p", &c.p);
    br(trkPrefix + "pT", &c.pT);
    br(trkPrefix + "delta_p", &c.delta_p);
    br(trkPrefix + "delta_pT", &c.delta_pT);
    br(trkPrefix + "theta", &c.theta);
    br(trkPrefix + "phi", &c.phi);
    br(trkPrefix + "px", &c.px);
    br(trkPrefix + "py", &c.py);
    br(trkPrefix + "pz", &c.pz);
    br(trkPrefix + "qOverP", &c.qOverP);
    br(trkPrefix + "momentum_resolved", &c.momentum_resolved);
    br(trkPrefix + "loc0", &c.loc0);
    br(trkPrefix + "loc1", &c.loc1);
    br(trkPrefix + "sigma_loc0", &c.sigma_loc0);
    br(trkPrefix + "sigma_loc1", &c.sigma_loc1);
    br(trkPrefix + "sigma_phi", &c.sigma_phi);
    br(trkPrefix + "sigma_theta", &c.sigma_theta);
    br(trkPrefix + "sigma_qOverP", &c.sigma_qOverP);
    br(trkPrefix + "sigma_time", &c.sigma_time);
    br(trkPrefix + "pull_qOverP", &c.pull_qOverP);
    br(trkPrefix + "pull_theta", &c.pull_theta);
    br(trkPrefix + "pull_phi", &c.pull_phi);
    br(trkPrefix + "chi2", &c.chi2);
    br(trkPrefix + "ndf", &c.ndf);
    br(trkPrefix + "charge", &c.charge);
    br(trkPrefix + "index", &c.index);
    br(trkPrefix + "object_index", &c.object_index);
    br(trkPrefix + "object_collectionID", &c.object_collectionID);
    br(trkPrefix + "seed_index", &c.seed_index);
    br(trkPrefix + "seed_collectionID", &c.seed_collectionID);
    br(trkPrefix + "identity_valid", &c.identity_valid);
    {
        static const std::array<std::pair<int,int>,21> ij{{
            {0,0},{0,1},{0,2},{0,3},{0,4},{0,5},{1,1},{1,2},{1,3},{1,4},{1,5},
            {2,2},{2,3},{2,4},{2,5},{3,3},{3,4},{3,5},{4,4},{4,5},{5,5}}};
        for (std::size_t k = 0; k < ij.size(); ++k)
            br(trkPrefix + "cov_" + std::to_string(ij[k].first) + std::to_string(ij[k].second), &c.cov_upper[k]);
    }
    br(trkPrefix + "type", &c.type);
    br(trkPrefix + "time", &c.time);
    br(trkPrefix + "pdg", &c.pdg);
    br(trkPrefix + "nStates", &c.nStates);
    br(trkPrefix + "nMeasurements", &c.nMeasurements);
    br(trkPrefix + "nOutliers", &c.nOutliers);
    br(trkPrefix + "nHoles", &c.nHoles);
    br(trkPrefix + "nSharedHits", &c.nSharedHits);
    br(trkPrefix + "assoc_mcIndex", &c.assoc_mcIndex);
    br(trkPrefix + "assoc_mcCollectionID", &c.assoc_mcCollectionID);
    br(trkPrefix + "assoc_weight", &c.assoc_weight);
    br(trkPrefix + "state_track_index", &c.state_track_index);
    br(trkPrefix + "state_parent_seed_index", &c.state_parent_seed_index);
    br(trkPrefix + "state_parent_seed_collectionID", &c.state_parent_seed_collectionID);
    br(trkPrefix + "state_parent_track_index", &c.state_parent_track_index);
    br(trkPrefix + "state_parent_track_collectionID", &c.state_parent_track_collectionID);
    br(trkPrefix + "state_parent_identity_valid", &c.state_parent_identity_valid);
    br(trkPrefix + "state_estimate_kind", &c.state_estimate_kind);
    br(trkPrefix + "state_index", &c.state_index);
    br(trkPrefix + "state_acts_index", &c.state_acts_index);
    br(trkPrefix + "state_type", &c.state_type);
    br(trkPrefix + "state_mapping_method", &c.state_mapping_method);
    br(trkPrefix + "state_surface", &c.state_surface);
    br(trkPrefix + "state_loc0", &c.state_loc0);
    br(trkPrefix + "state_loc1", &c.state_loc1);
    br(trkPrefix + "state_meas_loc0", &c.state_meas_loc0);
    br(trkPrefix + "state_meas_loc1", &c.state_meas_loc1);
    br(trkPrefix + "state_resid_loc0", &c.state_resid_loc0);
    br(trkPrefix + "state_resid_loc1", &c.state_resid_loc1);
    br(trkPrefix + "state_meas_dim", &c.state_meas_dim);
    br(trkPrefix + "state_proj_index0", &c.state_proj_index0); br(trkPrefix + "state_proj_index1", &c.state_proj_index1);
    br(trkPrefix + "state_meas_cov00", &c.state_meas_cov00); br(trkPrefix + "state_meas_cov01", &c.state_meas_cov01);
    br(trkPrefix + "state_meas_cov11", &c.state_meas_cov11);
    br(trkPrefix + "state_pred_loc0", &c.state_pred_loc0); br(trkPrefix + "state_pred_loc1", &c.state_pred_loc1);
    br(trkPrefix + "state_pred_theta", &c.state_pred_theta); br(trkPrefix + "state_pred_phi", &c.state_pred_phi);
    br(trkPrefix + "state_pred_qOverP", &c.state_pred_qOverP); br(trkPrefix + "state_pred_time", &c.state_pred_time);
    {
        static const std::array<std::pair<int,int>,21> ij{{
            {0,0},{0,1},{0,2},{0,3},{0,4},{0,5},{1,1},{1,2},{1,3},{1,4},{1,5},
            {2,2},{2,3},{2,4},{2,5},{3,3},{3,4},{3,5},{4,4},{4,5},{5,5}}};
        for (std::size_t k = 0; k < ij.size(); ++k)
            br(trkPrefix + "state_pred_cov_" + std::to_string(ij[k].first) + std::to_string(ij[k].second), &c.state_pred_cov_upper[k]);
    }
    br(trkPrefix + "state_innov0", &c.state_innov0); br(trkPrefix + "state_innov1", &c.state_innov1);
    br(trkPrefix + "state_innov_pull0", &c.state_innov_pull0); br(trkPrefix + "state_innov_pull1", &c.state_innov_pull1);
    br(trkPrefix + "state_innov_chi2", &c.state_innov_chi2);
    br(trkPrefix + "x_on_plane", &c.x_on_plane);
    br(trkPrefix + "y_on_plane", &c.y_on_plane);
    br(trkPrefix + "z_on_plane", &c.z_on_plane);
    br(trkPrefix + "aclgad_xPix", &c.aclgad_xPix);
    br(trkPrefix + "aclgad_yPix", &c.aclgad_yPix);
    br(trkPrefix + "aclgad_zPix", &c.aclgad_zPix);
    br(trkPrefix + "aclgad_dx", &c.aclgad_dx);
    br(trkPrefix + "aclgad_dy", &c.aclgad_dy);
    br(trkPrefix + "aclgad_dz", &c.aclgad_dz);
    br(trkPrefix + "aclgad_pixX", &c.aclgad_pixX);
    br(trkPrefix + "aclgad_pixY", &c.aclgad_pixY);
    br(trkPrefix + "aclgad_pixZ", &c.aclgad_pixZ);
    br(trkPrefix + "aclgad_plane", &c.aclgad_plane);
    br(trkPrefix + "aclgad_station", &c.aclgad_station);
    br(trkPrefix + "aclgad_module", &c.aclgad_module);
    br(trkPrefix + "aclgad_side", &c.aclgad_side);
    br(trkPrefix + "aclgad_sensor", &c.aclgad_sensor);
    br(trkPrefix + "aclgad_cellID", &c.aclgad_cellID);
    br(trkPrefix + "state_theta", &c.state_theta);
    br(trkPrefix + "state_phi", &c.state_phi);
    br(trkPrefix + "state_qOverP", &c.state_qOverP);
    br(trkPrefix + "state_time", &c.state_time);
    br(trkPrefix + "state_pdg", &c.state_pdg);
    br(trkPrefix + "n_sensor_map_exact", &c.nMapExact);
    br(trkPrefix + "n_sensor_map_fallback", &c.nMapFallback);
    br(trkPrefix + "n_sensor_map_failed_physics", &c.nMapFailedPhysics);
    br(trkPrefix + "n_sensor_map_unmapped_nonphysics", &c.nMapUnmappedNonPhysics);
    br(trkPrefix + "has_track", &c.hasTrack);
    bindBestSel(trkPrefix == "trk_" ? "oracle_best_trk_" : "ckf_oracle_best_trk_", c.oracle);
    bindBestSel(trkPrefix == "trk_" ? "best_trk_" : "ckf_best_trk_", c.oracle);
    bindBestSel(trkPrefix == "trk_" ? "truth_matched_trk_" : "ckf_truth_matched_trk_", c.truthMatched);
    bindBestSel(trkPrefix == "trk_" ? "reco_best_trk_" : "ckf_reco_best_trk_", c.recoBest);
}

void B0Trackers::Finish() {
}

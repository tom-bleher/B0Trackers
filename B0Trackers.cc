#include "B0Trackers.h"

#include <algorithm>
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

#include <JANA/JException.h>
#include <JANA/Services/JGlobalRootLock.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TGeoMatrix.h>
#include <TTree.h>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
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

// B0Trackers/hits branch map (schema 2):
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

template <typename T>
std::vector<const T*> getOpt(const std::shared_ptr<const JEvent>& event, const char* name) {
    try {
        return event->Get<T>(name);
    } catch (...) {
        return {};
    }
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
    app->SetDefaultParameter("B0Trackers:fallback_max_normal_mm", m_fallbackMaxNormalMm,
                             "Max |localZ| (mm) allowed for nearest-sensor fallback");
    app->SetDefaultParameter("B0Trackers:fail_on_empty_sensor_map", m_failOnEmptySensorMap,
                             "Throw in Init if the B0 sensor map is empty");

    if (const char* cfg = std::getenv("DETECTOR_CONFIG")) {
        m_geometryName = cfg;
    }
    if (const char* path = std::getenv("DETECTOR_PATH")) {
        m_detectorPath = path;
    }

    auto rootLock = app->GetService<JGlobalRootLock>();
    rootLock->acquire_write_lock();

    TDirectory* prevDir = gDirectory;
    auto rf_svc = app->GetService<RootFile_service>();
    TFile* outfile = rf_svc->GetHistFile()->GetFile();
    TDirectory* dir = outfile->GetDirectory("B0Trackers");
    if (dir == nullptr) {
        dir = outfile->mkdir("B0Trackers");
    }
    if (dir == nullptr) {
        rootLock->release_lock();
        throw JException("B0Trackers: failed to create output directory");
    }
    dir->cd();

    m_tree = new TTree("hits", "B0 truth, hits, and tracks");

    m_tree->Branch("schema_version", &m_schemaVersion);
    m_tree->Branch("geometry_name", &m_geometryName);
    m_tree->Branch("detector_path", &m_detectorPath);
    m_tree->Branch("eventNumber", &m_eventNumber);

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
    m_tree->Branch("survivedDigi", &vm_survivedDigi);
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
    m_tree->Branch("n_stations_primary", &m_nStationsPrimary);

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
        if (prevDir != nullptr) {
            prevDir->cd();
        }
        rootLock->release_lock();
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

    if (prevDir != nullptr) {
        prevDir->cd();
    }
    rootLock->release_lock();

    if (m_sensorRefs.empty() && m_failOnEmptySensorMap) {
        throw JException("B0Trackers: empty sensor map (check B0TrackerHits id fields / geometry)");
    }
}

void B0Trackers::Process(const std::shared_ptr<const JEvent>& event) {
    auto mcparticles = getOpt<edm4hep::MCParticle>(event, "MCParticles");
    auto simHits     = getOpt<edm4hep::SimTrackerHit>(event, "B0TrackerHits");
    auto rawHits     = getOpt<edm4eic::RawTrackerHit>(event, "B0TrackerRawHits");
    auto recHits     = getOpt<edm4eic::TrackerHit>(event, "B0TrackerRecHits");
    auto rawAssocs   = getOpt<edm4eic::MCRecoTrackerHitAssociation>(event, "B0TrackerRawHitAssociations");
    auto measurements = getOpt<edm4eic::Measurement2D>(event, "B0TrackerMeasurements");
    auto stubSeeds   = getOpt<edm4eic::TrackSeed>(event, "B0TrackerSeeds");
    auto truthSeeds  = getOpt<edm4eic::TrackSeed>(event, "B0TrackerTruthSeeds");

    auto tsTracks = getOpt<edm4eic::TrackParameters>(event, "B0TrackerCKFTruthSeededTrackParameters");
    auto tsTrajectories = getOpt<edm4eic::Trajectory>(event, "B0TrackerCKFTruthSeededTrajectories");
    auto tsEdmTracks = getOpt<edm4eic::Track>(event, "B0TrackerCKFTruthSeededTracks");
    auto tsAssocs = getOpt<edm4eic::MCRecoTrackParticleAssociation>(
        event, "B0TrackerCKFTruthSeededTrackAssociations");
    auto tsActsTrackStates =
        getOpt<Acts::ConstVectorMultiTrajectory>(event, "B0TrackerCKFTruthSeededActsTrackStates");
    auto tsActsTracks = getOpt<Acts::ConstVectorTrackContainer>(event, "B0TrackerCKFTruthSeededActsTracks");
    auto tsUnfiltered = getOpt<edm4eic::Track>(event, "B0TrackerCKFTruthSeededTracksUnfiltered");

    auto ckfTracks = getOpt<edm4eic::TrackParameters>(event, "B0TrackerCKFTrackParameters");
    auto ckfTrajectories = getOpt<edm4eic::Trajectory>(event, "B0TrackerCKFTrajectories");
    auto ckfEdmTracks = getOpt<edm4eic::Track>(event, "B0TrackerCKFTracks");
    auto ckfAssocs = getOpt<edm4eic::MCRecoTrackParticleAssociation>(event, "B0TrackerCKFTrackAssociations");
    auto ckfActsTrackStates =
        getOpt<Acts::ConstVectorMultiTrajectory>(event, "B0TrackerCKFActsTrackStates");
    auto ckfActsTracks = getOpt<Acts::ConstVectorTrackContainer>(event, "B0TrackerCKFActsTracks");
    auto ckfUnfiltered = getOpt<edm4eic::Track>(event, "B0TrackerCKFTracksUnfiltered");

    std::lock_guard<std::mutex> lock(m_fillMutex);

    m_eventNumber = event->GetEventNumber();

    vm_xR.clear();    vm_yR.clear();    vm_zR.clear();
    vm_xT.clear();    vm_yT.clear();    vm_zT.clear();
    vm_aclgad_xPixT.clear(); vm_aclgad_yPixT.clear(); vm_aclgad_zPixT.clear();
    vm_aclgad_dxT.clear(); vm_aclgad_dyT.clear(); vm_aclgad_dzT.clear();
    vm_aclgad_xPixR.clear(); vm_aclgad_yPixR.clear(); vm_aclgad_zPixR.clear();
    vm_aclgad_dxR.clear(); vm_aclgad_dyR.clear(); vm_aclgad_dzR.clear();
    vm_detX.clear();  vm_detY.clear();  vm_detZ.clear();
    vm_plane.clear(); vm_station.clear(); vm_module.clear(); vm_side.clear(); vm_sensor.clear();
    vm_pixX.clear();  vm_pixY.clear();   vm_pixZ.clear();
    vm_aclgad_pixXT.clear(); vm_aclgad_pixYT.clear();
    vm_aclgad_pixXR.clear(); vm_aclgad_pixYR.clear();
    vm_cellID.clear(); vm_mcIndex.clear(); vm_mcCollectionID.clear();
    vm_eDep.clear();  vm_time.clear();   vm_path.clear();
    vm_pdg.clear();   vm_status.clear();
    vm_isPrimary.clear(); vm_isSelPrimary.clear(); vm_survivedDigi.clear();
    vm_px.clear();    vm_py.clear();    vm_pz.clear();   vm_p.clear();    vm_pT.clear();

    vm_xP.clear();      vm_yP.clear();      vm_zP.clear();      vm_pathP.clear();   vm_timeP.clear();
    vm_planeP.clear();  vm_stationP.clear(); vm_moduleP.clear(); vm_sideP.clear(); vm_sensorP.clear();
    vm_pdgP.clear();    vm_statusP.clear(); vm_isPrimaryP.clear(); vm_isSelPrimaryP.clear();
    vm_mcIndexP.clear(); vm_mcCollectionIDP.clear();
    vm_pxP.clear();     vm_pyP.clear();     vm_pzP.clear();     vm_pP.clear();     vm_pTP.clear();

    vm_xEntry.clear();   vm_yEntry.clear();  vm_zEntry.clear();  vm_timeEntry.clear();
    vm_pxEntry.clear();  vm_pyEntry.clear(); vm_pzEntry.clear(); vm_pEntry.clear(); vm_pTEntry.clear();
    vm_xExit.clear();    vm_yExit.clear();   vm_zExit.clear();   vm_timeExit.clear();
    vm_pxExit.clear();   vm_pyExit.clear();  vm_pzExit.clear();  vm_pExit.clear();  vm_pTExit.clear();
    vm_sensorEntry.clear(); vm_sensorExit.clear();
    vm_cellIDEntry.clear(); vm_cellIDExit.clear();
    vm_planeEE.clear();  vm_stationEE.clear(); vm_moduleEE.clear(); vm_sideEE.clear();  vm_pdgEE.clear();
    vm_isPrimaryEE.clear(); vm_isSelPrimaryEE.clear();
    vm_mcIndexEE.clear(); vm_mcCollectionIDEE.clear(); vm_nStepsEE.clear();
    vm_statusEE.clear();

    vm_raw_cellID.clear(); vm_raw_charge.clear(); vm_raw_timeStamp.clear();
    vm_raw_plane.clear(); vm_raw_station.clear(); vm_raw_module.clear();
    vm_raw_side.clear(); vm_raw_sensor.clear();
    vm_raw_mcIndex.clear(); vm_raw_mcCollectionID.clear();
    vm_rec_x.clear(); vm_rec_y.clear(); vm_rec_z.clear();
    vm_rec_covxx.clear(); vm_rec_covyy.clear(); vm_rec_covzz.clear();
    vm_rec_time.clear(); vm_rec_time_err.clear(); vm_rec_edep.clear(); vm_rec_edep_err.clear();
    vm_rec_cellID.clear();
    vm_rec_plane.clear(); vm_rec_station.clear(); vm_rec_module.clear();
    vm_rec_side.clear(); vm_rec_sensor.clear();
    vm_rec_pixX.clear(); vm_rec_pixY.clear(); vm_rec_pixZ.clear();
    vm_rec_mcIndex.clear(); vm_rec_mcCollectionID.clear();

    vm_seed_quality.clear(); vm_seed_p.clear(); vm_seed_qOverP.clear();
    vm_seed_theta.clear(); vm_seed_phi.clear(); vm_seed_loc0.clear(); vm_seed_loc1.clear();
    vm_seed_sigma_qOverP.clear(); vm_seed_sigma_theta.clear(); vm_seed_sigma_phi.clear();
    vm_seed_nHits.clear(); vm_seed_charge.clear();
    vm_seed_momentum_resolved.clear(); vm_seed_became_track.clear();
    vm_truth_seed_quality.clear(); vm_truth_seed_p.clear(); vm_truth_seed_qOverP.clear();
    vm_truth_seed_theta.clear(); vm_truth_seed_phi.clear();
    vm_truth_seed_loc0.clear(); vm_truth_seed_loc1.clear();
    vm_truth_seed_sigma_qOverP.clear(); vm_truth_seed_sigma_theta.clear(); vm_truth_seed_sigma_phi.clear();
    vm_truth_seed_nHits.clear(); vm_truth_seed_charge.clear();
    vm_truth_seed_momentum_resolved.clear(); vm_truth_seed_became_track.clear();

    const double nan = b0trk::quietNaN();
    m_ts.clear(nan);
    m_ckf.clear(nan);
    m_selPrimaryMcIndex = -1;
    m_selPrimaryMcCollectionID = 0;
    m_selPrimaryPx = nan;
    m_selPrimaryPy = nan;
    m_selPrimaryPz = nan;
    m_selPrimaryP = nan;
    m_selPrimaryPT = nan;
    m_selPrimaryThscatMrad = nan;
    m_nStationsPrimary = 0;

    beam_px.clear();  beam_py.clear();  beam_pz.clear(); beam_p.clear(); beam_pT.clear();
    beam_pdg.clear();
    m_primaryPx.clear(); m_primaryPy.clear(); m_primaryPz.clear(); m_primaryP.clear(); m_primaryPT.clear();
    m_primaryPdgOut.clear(); m_primaryStatusOut.clear(); m_primaryMcIndex.clear();
    m_primaryMcCollectionID.clear();
    m_genPpx.clear();    m_genPpy.clear();    m_genPpz.clear();    m_genPp.clear();    m_genPpT.clear();
    m_genBeamPx.clear(); m_genBeamPy.clear(); m_genBeamPz.clear(); m_genBeamP.clear(); m_genBeamPT.clear();
    m_genBeamPPx.clear();m_genBeamPPy.clear();m_genBeamPPz.clear();m_genBeamPP.clear();m_genBeamPPT.clear();

    m_nSimHits = static_cast<int>(simHits.size());
    m_nRawHits = static_cast<int>(rawHits.size());
    m_nRecHits = static_cast<int>(recHits.size());
    m_nMeasurements = static_cast<int>(measurements.size());
    m_nTruthSeeds = static_cast<int>(truthSeeds.size());
    m_nStubSeeds = static_cast<int>(stubSeeds.size());
    m_nTsUnfiltered = static_cast<int>(tsUnfiltered.size());
    m_nTsFiltered = static_cast<int>(tsEdmTracks.size());
    m_nCkfUnfiltered = static_cast<int>(ckfUnfiltered.size());
    m_nCkfFiltered = static_cast<int>(ckfEdmTracks.size());
    m_nSimHitsUnresolvedCellID = 0;
    m_nMissingMcRelation = 0;
    m_nSensorMapExact = 0;
    m_nSensorMapFallback = 0;
    m_nSensorMapFailed = 0;
    m_nPixelSnapFailed = 0;

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

    struct SimLink {
        int mcIndex = -1;
        std::uint32_t mcCollectionID = 0;
        double eDep = -1.0;
    };
    std::unordered_map<std::uint64_t, SimLink> simLinkByCell;
    for (const auto* assoc : rawAssocs) {
        if (assoc == nullptr) continue;
        const auto raw = assoc->getRawHit();
        const auto sim = assoc->getSimHit();
        if (!raw.isAvailable() || !sim.isAvailable()) continue;
        const auto particle = sim.getParticle();
        if (!particle.isAvailable()) continue;
        const auto id = particle.id();
        auto& link = simLinkByCell[raw.getCellID()];
        if (sim.getEDep() > link.eDep) {
            link.eDep = sim.getEDep();
            link.mcIndex = id.index;
            link.mcCollectionID = id.collectionID;
        }
    }

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
        m_primaryPdgOut.push_back(pdg);
        m_primaryStatusOut.push_back(status);
        m_primaryMcIndex.push_back(id.index);
        m_primaryMcCollectionID.push_back(id.collectionID);
        m_primaryPx.push_back(p.x);
        m_primaryPy.push_back(p.y);
        m_primaryPz.push_back(p.z);
        m_primaryP.push_back(pmag);
        m_primaryPT.push_back(pT);
    }

    const auto isSelectedPrimary = [&primaryIds](std::uint32_t collectionID, int index) -> int {
        for (const auto& selected : primaryIds) {
            if (selected.first == collectionID && selected.second == index) return 1;
        }
        return 0;
    };

    std::size_t primaryRef = 0;
    for (std::size_t i = 1; i < m_primaryP.size(); ++i) {
        if (m_primaryP[i] > m_primaryP[primaryRef]) primaryRef = i;
    }
    if (!m_primaryP.empty()) {
        m_selPrimaryMcIndex = m_primaryMcIndex[primaryRef];
        m_selPrimaryMcCollectionID = m_primaryMcCollectionID[primaryRef];
        m_selPrimaryPx = m_primaryPx[primaryRef];
        m_selPrimaryPy = m_primaryPy[primaryRef];
        m_selPrimaryPz = m_primaryPz[primaryRef];
        m_selPrimaryP = m_primaryP[primaryRef];
        m_selPrimaryPT = m_primaryPT[primaryRef];
    }
    const auto isSelPrimary = [this](std::uint32_t collectionID, int index) -> int {
        return (m_selPrimaryMcIndex >= 0 &&
                collectionID == m_selPrimaryMcCollectionID &&
                index == m_selPrimaryMcIndex) ? 1 : 0;
    };

    for (const auto* h : simHits) {
        const auto mc = h->getParticle();
        if (!mc.isAvailable()) {
            ++m_nMissingMcRelation;
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
            ++m_nSimHitsUnresolvedCellID;
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
            ++m_nPixelSnapFailed;
        }

        vm_xR.push_back(10. * gpos.x());
        vm_yR.push_back(10. * gpos.y());
        vm_zR.push_back(10. * gpos.z());
        vm_xT.push_back(truthPos.x);
        vm_yT.push_back(truthPos.y);
        vm_zT.push_back(truthPos.z);
        vm_aclgad_xPixT.push_back(truthPixel.x);
        vm_aclgad_yPixT.push_back(truthPixel.y);
        vm_aclgad_zPixT.push_back(truthPixel.z);
        vm_aclgad_dxT.push_back(truthPixel.dx);
        vm_aclgad_dyT.push_back(truthPixel.dy);
        vm_aclgad_dzT.push_back(truthPixel.dz);
        vm_aclgad_xPixR.push_back(readoutPixel.x);
        vm_aclgad_yPixR.push_back(readoutPixel.y);
        vm_aclgad_zPixR.push_back(readoutPixel.z);
        vm_aclgad_dxR.push_back(readoutPixel.dx);
        vm_aclgad_dyR.push_back(readoutPixel.dy);
        vm_aclgad_dzR.push_back(readoutPixel.dz);
        vm_detX.push_back(10. * lpos.x());
        vm_detY.push_back(10. * lpos.y());
        vm_detZ.push_back(10. * lpos.z());
        vm_plane .push_back(plane);
        vm_station.push_back(stationOf(plane));
        vm_module.push_back(module);
        vm_side  .push_back(side);
        vm_sensor.push_back(sensor);
        vm_pixX  .push_back(pixX);
        vm_pixY  .push_back(pixY);
        vm_pixZ  .push_back(pixZ);
        vm_aclgad_pixXT.push_back(truthPixel.pixX);
        vm_aclgad_pixYT.push_back(truthPixel.pixY);
        vm_aclgad_pixXR.push_back(readoutPixel.pixX);
        vm_aclgad_pixYR.push_back(readoutPixel.pixY);
        vm_cellID.push_back(cid);
        vm_eDep  .push_back(h->getEDep());
        vm_time  .push_back(h->getTime());
        vm_path  .push_back(path);
        vm_pdg   .push_back(mc.getPDG());
        vm_mcIndex.push_back(id.index);
        vm_mcCollectionID.push_back(id.collectionID);
        vm_px    .push_back(mom.x);
        vm_py    .push_back(mom.y);
        vm_pz    .push_back(mom.z);
        vm_p     .push_back(pmag);
        vm_pT    .push_back(pT);
        vm_status.push_back(mc.getGeneratorStatus());
        vm_isPrimary.push_back(primaryFlag);
        vm_isSelPrimary.push_back(selFlag);
        vm_survivedDigi.push_back(rawCellIDs.count(cid) ? 1 : 0);

        const auto key = std::make_tuple(id.collectionID, id.index, plane, side, module, sensor);
        const auto existing = penetrationIndex.find(key);
        if (existing == penetrationIndex.end()) {
            penetrationIndex.emplace(key, vm_xP.size());
            vm_xP     .push_back(truthPos.x);
            vm_yP     .push_back(truthPos.y);
            vm_zP     .push_back(truthPos.z);
            vm_pathP  .push_back(path);
            vm_timeP  .push_back(h->getTime());
            vm_planeP .push_back(plane);
            vm_stationP.push_back(stationOf(plane));
            vm_moduleP.push_back(module);
            vm_sideP  .push_back(side);
            vm_sensorP.push_back(sensor);
            vm_pdgP   .push_back(mc.getPDG());
            vm_statusP.push_back(mc.getGeneratorStatus());
            vm_isPrimaryP.push_back(primaryFlag);
            vm_isSelPrimaryP.push_back(selFlag);
            vm_mcIndexP.push_back(id.index);
            vm_mcCollectionIDP.push_back(id.collectionID);
            vm_pxP    .push_back(mom.x);
            vm_pyP    .push_back(mom.y);
            vm_pzP    .push_back(mom.z);
            vm_pP     .push_back(pmag);
            vm_pTP    .push_back(pT);
        } else if (h->getTime() < vm_timeP[existing->second]) {
            const std::size_t idx = existing->second;
            vm_xP[idx]    = truthPos.x;
            vm_yP[idx]    = truthPos.y;
            vm_zP[idx]    = truthPos.z;
            vm_pathP[idx] = path;
            vm_timeP[idx] = h->getTime();
            vm_pxP[idx]   = mom.x;
            vm_pyP[idx]   = mom.y;
            vm_pzP[idx]   = mom.z;
            vm_pP[idx]    = pmag;
            vm_pTP[idx]   = pT;
        }

        if (side >= 0) {
            const auto eeKey = std::make_tuple(id.collectionID, id.index, plane, side, module);
            const double thisTime = h->getTime();
            const auto eeIt = entryExitIndex.find(eeKey);
            if (eeIt == entryExitIndex.end()) {
                entryExitIndex.emplace(eeKey, vm_xEntry.size());
                vm_xEntry .push_back(truthPos.x);   vm_yEntry .push_back(truthPos.y);
                vm_zEntry .push_back(truthPos.z);   vm_timeEntry.push_back(thisTime);
                vm_pxEntry.push_back(mom.x);        vm_pyEntry.push_back(mom.y);
                vm_pzEntry.push_back(mom.z);        vm_pEntry .push_back(pmag);
                vm_pTEntry.push_back(pT);
                vm_xExit  .push_back(truthPos.x);   vm_yExit  .push_back(truthPos.y);
                vm_zExit  .push_back(truthPos.z);   vm_timeExit.push_back(thisTime);
                vm_pxExit .push_back(mom.x);        vm_pyExit .push_back(mom.y);
                vm_pzExit .push_back(mom.z);        vm_pExit  .push_back(pmag);
                vm_pTExit .push_back(pT);
                vm_sensorEntry.push_back(sensor);   vm_sensorExit.push_back(sensor);
                vm_cellIDEntry.push_back(cid);      vm_cellIDExit.push_back(cid);
                vm_planeEE.push_back(plane);
                vm_stationEE.push_back(stationOf(plane));
                vm_moduleEE.push_back(module);
                vm_sideEE.push_back(side);
                vm_pdgEE  .push_back(mc.getPDG());
                vm_isPrimaryEE.push_back(primaryFlag);
                vm_isSelPrimaryEE.push_back(selFlag);
                vm_mcIndexEE.push_back(id.index);
                vm_mcCollectionIDEE.push_back(id.collectionID);
                vm_nStepsEE.push_back(1);
                vm_statusEE.push_back(mc.getGeneratorStatus());
            } else {
                const std::size_t idx = eeIt->second;
                ++vm_nStepsEE[idx];
                if (thisTime < vm_timeEntry[idx]) {
                    vm_xEntry[idx]    = truthPos.x;
                    vm_yEntry[idx]    = truthPos.y;
                    vm_zEntry[idx]    = truthPos.z;
                    vm_timeEntry[idx] = thisTime;
                    vm_pxEntry[idx]   = mom.x;
                    vm_pyEntry[idx]   = mom.y;
                    vm_pzEntry[idx]   = mom.z;
                    vm_pEntry[idx]    = pmag;
                    vm_pTEntry[idx]   = pT;
                    vm_sensorEntry[idx] = sensor;
                    vm_cellIDEntry[idx] = cid;
                }
                if (thisTime > vm_timeExit[idx]) {
                    vm_xExit[idx]    = truthPos.x;
                    vm_yExit[idx]    = truthPos.y;
                    vm_zExit[idx]    = truthPos.z;
                    vm_timeExit[idx] = thisTime;
                    vm_pxExit[idx]   = mom.x;
                    vm_pyExit[idx]   = mom.y;
                    vm_pzExit[idx]   = mom.z;
                    vm_pExit[idx]    = pmag;
                    vm_pTExit[idx]   = pT;
                    vm_sensorExit[idx] = sensor;
                    vm_cellIDExit[idx] = cid;
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
            beam_px .push_back(p.x);
            beam_py .push_back(p.y);
            beam_pz .push_back(p.z);
            beam_p  .push_back(pmag);
            beam_pT .push_back(pT);
            beam_pdg.push_back(pdg);
        }
        if (pdg == 2212 && status == 4) {
            m_genPpx.push_back(p.x);
            m_genPpy.push_back(p.y);
            m_genPpz.push_back(p.z);
            m_genPp .push_back(pmag);
            m_genPpT.push_back(pT);
        }
        if (pdg == 2212 && status == 1) {
            m_genBeamPPx.push_back(p.x);
            m_genBeamPPy.push_back(p.y);
            m_genBeamPPz.push_back(p.z);
            m_genBeamPP .push_back(pmag);
            m_genBeamPPT.push_back(pT);
        }
        if (pdg == 11 && status == 1) {
            m_genBeamPx.push_back(p.x);
            m_genBeamPy.push_back(p.y);
            m_genBeamPz.push_back(p.z);
            m_genBeamP .push_back(pmag);
            m_genBeamPT.push_back(pT);
        }
    }

    if (std::isfinite(m_selPrimaryP) && !m_genPp.empty()) {
        const double bx = m_genPpx[0];
        const double by = m_genPpy[0];
        const double bz = m_genPpz[0];
        const double bn = std::sqrt(bx * bx + by * by + bz * bz);
        const double pn = m_selPrimaryP;
        if (bn > 0.0 && pn > 0.0) {
            const double cosang = std::clamp(
                (m_selPrimaryPx * bx + m_selPrimaryPy * by + m_selPrimaryPz * bz) / (pn * bn),
                -1.0, 1.0);
            m_selPrimaryThscatMrad = 1.0e3 * std::acos(cosang);
        }
    }
    {
        std::set<int> stations;
        for (std::size_t i = 0; i < vm_stationP.size(); ++i) {
            if (vm_isSelPrimaryP[i] == 1 && vm_stationP[i] > 0) {
                stations.insert(vm_stationP[i]);
            }
        }
        m_nStationsPrimary = static_cast<int>(stations.size());
    }

    for (const auto* raw : rawHits) {
        if (raw == nullptr) continue;
        const auto cid = raw->getCellID();
        int plane = -1, module = -1, sensor = -1, side = -1;
        decodeIds(cid, plane, module, sensor, side);
        const auto linkIt = simLinkByCell.find(cid);
        vm_raw_cellID.push_back(cid);
        vm_raw_charge.push_back(raw->getCharge());
        vm_raw_timeStamp.push_back(raw->getTimeStamp());
        vm_raw_plane.push_back(plane);
        vm_raw_station.push_back(stationOf(plane));
        vm_raw_module.push_back(module);
        vm_raw_side.push_back(side);
        vm_raw_sensor.push_back(sensor);
        vm_raw_mcIndex.push_back(linkIt == simLinkByCell.end() ? -1 : linkIt->second.mcIndex);
        vm_raw_mcCollectionID.push_back(linkIt == simLinkByCell.end() ? 0 : linkIt->second.mcCollectionID);
    }

    for (const auto* rec : recHits) {
        if (rec == nullptr) continue;
        const auto cid = rec->getCellID();
        int plane = -1, module = -1, sensor = -1, side = -1;
        decodeIds(cid, plane, module, sensor, side);
        const auto pos = rec->getPosition();
        const auto err = rec->getPositionError();
        const auto linkIt = simLinkByCell.find(cid);
        vm_rec_x.push_back(pos.x);
        vm_rec_y.push_back(pos.y);
        vm_rec_z.push_back(pos.z);
        vm_rec_covxx.push_back(err.xx);
        vm_rec_covyy.push_back(err.yy);
        vm_rec_covzz.push_back(err.zz);
        vm_rec_time.push_back(rec->getTime());
        vm_rec_time_err.push_back(rec->getTimeError());
        vm_rec_edep.push_back(rec->getEdep());
        vm_rec_edep_err.push_back(rec->getEdepError());
        vm_rec_cellID.push_back(cid);
        vm_rec_plane.push_back(plane);
        vm_rec_station.push_back(stationOf(plane));
        vm_rec_module.push_back(module);
        vm_rec_side.push_back(side);
        vm_rec_sensor.push_back(sensor);
        vm_rec_pixX.push_back(getFieldOr(cid, "x", -1));
        vm_rec_pixY.push_back(getFieldOr(cid, "y", -1));
        vm_rec_pixZ.push_back(getFieldOr(cid, "z", -1));
        vm_rec_mcIndex.push_back(linkIt == simLinkByCell.end() ? -1 : linkIt->second.mcIndex);
        vm_rec_mcCollectionID.push_back(linkIt == simLinkByCell.end() ? 0 : linkIt->second.mcCollectionID);
    }

    const auto fillSeeds = [&](const auto& seeds, const auto& trajectories,
                               auto& quality, auto& p, auto& qOverP, auto& theta, auto& phi,
                               auto& loc0, auto& loc1, auto& sigQ, auto& sigTh, auto& sigPh,
                               auto& nHits, auto& charge, auto& resolved, auto& became) {
        std::unordered_set<int> usedSeedIndex;
        for (const auto* traj : trajectories) {
            if (traj == nullptr) continue;
            const auto seed = traj->getSeed();
            if (seed.isAvailable()) {
                usedSeedIndex.insert(seed.id().index);
            }
        }
        for (std::size_t i = 0; i < seeds.size(); ++i) {
            const auto* seed = seeds[i];
            if (seed == nullptr) continue;
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
            became.push_back(usedSeedIndex.count(static_cast<int>(i)) ? 1 : 0);
        }
    };
    fillSeeds(stubSeeds, ckfTrajectories, vm_seed_quality, vm_seed_p, vm_seed_qOverP,
              vm_seed_theta, vm_seed_phi, vm_seed_loc0, vm_seed_loc1,
              vm_seed_sigma_qOverP, vm_seed_sigma_theta, vm_seed_sigma_phi,
              vm_seed_nHits, vm_seed_charge, vm_seed_momentum_resolved, vm_seed_became_track);
    fillSeeds(truthSeeds, tsTrajectories, vm_truth_seed_quality, vm_truth_seed_p, vm_truth_seed_qOverP,
              vm_truth_seed_theta, vm_truth_seed_phi, vm_truth_seed_loc0, vm_truth_seed_loc1,
              vm_truth_seed_sigma_qOverP, vm_truth_seed_sigma_theta, vm_truth_seed_sigma_phi,
              vm_truth_seed_nHits, vm_truth_seed_charge, vm_truth_seed_momentum_resolved,
              vm_truth_seed_became_track);

    auto fillChain = [&](TrackChain& out,
                         const auto& trajectories,
                         const auto& tracks,
                         const auto& edmTracks,
                         const auto& assocs,
                         const auto& actsTracks,
                         const auto& actsStates) {
        struct AssocHit {
            int mcIndex = -1;
            std::uint32_t mcCollectionID = 0;
            double weight = 0.0;
        };
        std::unordered_map<int, AssocHit> bestAssoc;
        for (const auto* assoc : assocs) {
            if (assoc == nullptr) continue;
            const auto rec = assoc->getRec();
            const auto sim = assoc->getSim();
            if (!rec.isAvailable() || !sim.isAvailable()) continue;
            const int recIndex = rec.id().index;
            const double weight = assoc->getWeight();
            auto it = bestAssoc.find(recIndex);
            if (it == bestAssoc.end() || weight > it->second.weight) {
                bestAssoc[recIndex] = {sim.id().index, sim.id().collectionID, weight};
            }
        }

        const double primaryP = m_selPrimaryP;
        const double primaryPT = m_selPrimaryPT;
        const double truthQOverP = (std::isfinite(primaryP) && primaryP > 0.0)
            ? (1.0 / primaryP) : nan;
        const double truthTheta = (std::isfinite(primaryP) && primaryP > 0.0)
            ? std::acos(std::clamp(m_selPrimaryPz / primaryP, -1.0, 1.0)) : nan;
        const double truthPhi = (std::isfinite(m_selPrimaryPx) && std::isfinite(m_selPrimaryPy))
            ? std::atan2(m_selPrimaryPy, m_selPrimaryPx) : nan;

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
            if (!tp.isAvailable() && trajIndex < tracks.size() && tracks[trajIndex] != nullptr) {
                tp = *tracks[trajIndex];
            }
            if (!tp.isAvailable()) continue;

            const float theta  = tp.getTheta();
            const float phi    = tp.getPhi();
            const float qOverP = tp.getQOverP();
            const auto [p, momOk] = b0trk::momentumFromQOverP(qOverP);
            const double pT = (momOk && std::isfinite(theta)) ? std::abs(p * std::sin(theta)) : nan;
            const int charge = (qOverP > 0.f) ? 1 : ((qOverP < 0.f) ? -1 : 0);
            const double deltaP = (momOk && std::isfinite(primaryP)) ? p - primaryP : nan;
            const double deltaPT = (std::isfinite(pT) && std::isfinite(primaryPT)) ? pT - primaryPT : nan;
            const int currentTrackIndex = static_cast<int>(trajIndex);
            const int nStates = static_cast<int>(trajectory->getNStates());
            const int nMeasurements = static_cast<int>(trajectory->getNMeasurements());
            const int nOutliers = static_cast<int>(trajectory->getNOutliers());
            const int nHoles = static_cast<int>(trajectory->getNHoles());
            const int nSharedHits = static_cast<int>(trajectory->getNSharedHits());

            double chi2 = nan;
            int ndf = -1;
            int pdg = 0;
            if (trajIndex < edmTracks.size() && edmTracks[trajIndex] != nullptr) {
                chi2 = edmTracks[trajIndex]->getChi2();
                ndf = static_cast<int>(edmTracks[trajIndex]->getNdf());
                pdg = edmTracks[trajIndex]->getPdg();
            }
            int assocMc = -1;
            std::uint32_t assocCol = 0;
            double assocW = 0.0;
            if (const auto it = bestAssoc.find(currentTrackIndex); it != bestAssoc.end()) {
                assocMc = it->second.mcIndex;
                assocCol = it->second.mcCollectionID;
                assocW = it->second.weight;
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
            const bool truthAssoc = (assocMc == m_selPrimaryMcIndex &&
                                     assocCol == m_selPrimaryMcCollectionID &&
                                     m_selPrimaryMcIndex >= 0);
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
            out.type.push_back(tp.getType());
            out.surface.push_back(tp.getSurface());
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
        for (int actsTrackIndex = 0; actsTrackIndex < nActsTracks; ++actsTrackIndex) {
            const auto track = trackContainer.getTrack(actsTrackIndex);
            const std::size_t stateBegin = out.state_index.size();
            int stateIndex = 0;
            for (const auto& state : track.trackStatesReversed()) {
                if (!state.hasReferenceSurface()) {
                    continue;
                }

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
                    stateTime = params[Acts::eBoundTime];
                });

                int typeMask = 0;
                const auto flags = state.typeFlags();
                if (flags.test(Acts::TrackStateFlag::MeasurementFlag)) typeMask |= (1 << 0);
                if (flags.test(Acts::TrackStateFlag::ParameterFlag))   typeMask |= (1 << 1);
                if (flags.test(Acts::TrackStateFlag::OutlierFlag))     typeMask |= (1 << 2);
                if (flags.test(Acts::TrackStateFlag::HoleFlag))        typeMask |= (1 << 3);
                if (flags.test(Acts::TrackStateFlag::MaterialFlag))    typeMask |= (1 << 4);
                if (flags.test(Acts::TrackStateFlag::SharedHitFlag))   typeMask |= (1 << 5);
                if (flags.test(Acts::TrackStateFlag::SplitHitFlag))    typeMask |= (1 << 6);
                if (flags.test(Acts::TrackStateFlag::NoExpectedHitFlag)) typeMask |= (1 << 7);
                const bool physicsState =
                    flags.test(Acts::TrackStateFlag::MeasurementFlag) ||
                    flags.test(Acts::TrackStateFlag::OutlierFlag) ||
                    flags.test(Acts::TrackStateFlag::HoleFlag);

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
                    ++m_nSensorMapExact;
                } else if (physicsState) {
                    const auto hit = closestSensor(globalX, globalY, globalZ);
                    if (hit.ref != nullptr && hit.outsideX == 0.0 && hit.outsideY == 0.0 &&
                        std::abs(hit.localZ) <= m_fallbackMaxNormalMm) {
                        sensorRef = hit.ref;
                        mapping = b0trk::kMapFallback;
                        ++m_nSensorMapFallback;
                    } else {
                        ++m_nSensorMapFailed;
                    }
                } else {
                    ++m_nSensorMapFailed;
                }

                PixelSnap trackPixel{};
                if (sensorRef != nullptr) {
                    trackPixel = snapToAclgadPixel(globalX, globalY, globalZ, sensorRef->cellID,
                                                   sensorRef->detElement);
                    if (!trackPixel.ok) {
                        ++m_nPixelSnapFailed;
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
                out.state_pdg.push_back(
                    (actsTrackIndex >= 0 && static_cast<std::size_t>(actsTrackIndex) < out.pdg.size())
                        ? out.pdg[static_cast<std::size_t>(actsTrackIndex)]
                        : 0);
            }
            const int nPushed = stateIndex;
            for (int i = 0; i < nPushed; ++i) {
                out.state_index[stateBegin + static_cast<std::size_t>(i)] = nPushed - 1 - i;
            }
        }
    };

    fillChain(m_ts, tsTrajectories, tsTracks, tsEdmTracks, tsAssocs, tsActsTracks, tsActsTrackStates);
    fillChain(m_ckf, ckfTrajectories, ckfTracks, ckfEdmTracks, ckfAssocs, ckfActsTracks, ckfActsTrackStates);

    m_tree->Fill();
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
    index.clear(); charge.clear(); type.clear(); pdg.clear(); surface.clear();
    nStates.clear(); nMeasurements.clear(); nOutliers.clear(); nHoles.clear(); nSharedHits.clear();
    assoc_mcIndex.clear(); assoc_mcCollectionID.clear(); assoc_weight.clear();
    state_track_index.clear(); state_index.clear(); state_acts_index.clear();
    state_type.clear(); state_pdg.clear(); state_mapping_method.clear();
    state_surface.clear();
    state_loc0.clear(); state_loc1.clear();
    state_meas_loc0.clear(); state_meas_loc1.clear();
    state_resid_loc0.clear(); state_resid_loc1.clear();
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
    br(trkPrefix + "type", &c.type);
    br(trkPrefix + "surface", &c.surface);
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
    br(trkPrefix + "has_track", &c.hasTrack);
    bindBestSel(trkPrefix == "trk_" ? "oracle_best_trk_" : "ckf_oracle_best_trk_", c.oracle);
    bindBestSel(trkPrefix == "trk_" ? "best_trk_" : "ckf_best_trk_", c.oracle);
    bindBestSel(trkPrefix == "trk_" ? "truth_matched_trk_" : "ckf_truth_matched_trk_", c.truthMatched);
    bindBestSel(trkPrefix == "trk_" ? "reco_best_trk_" : "ckf_reco_best_trk_", c.recoBest);
}

void B0Trackers::Finish() {
}

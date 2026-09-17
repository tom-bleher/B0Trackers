#pragma once

#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <JANA/JEventProcessor.h>

#include <DD4hep/DetElement.h>
#include <DD4hep/Segmentations.h>
#include <DD4hep/VolumeManager.h>
#include <DDSegmentation/BitFieldCoder.h>

#include <services/geometry/dd4hep/DD4hep_service.h>

#include "B0TrackersHelpers.h"

class TTree;
class ACTSGeo_service;
class ActsGeometryProvider;

namespace spdlog {
class logger;
}

class B0Trackers : public JEventProcessor {
public:
    void Init() override;
    void Process(const std::shared_ptr<const JEvent>& event) override;
    void Finish() override;

private:
    // Serializes only the final branch-buffer swap + TTree::Fill(). Event
    // analysis is performed in thread-local EventBuffers outside this lock.
    std::mutex m_fillMutex;

    std::shared_ptr<spdlog::logger> m_log;

    std::shared_ptr<DD4hep_service> m_geoSvc;
    std::shared_ptr<ACTSGeo_service> m_actsGeoSvc;
    std::shared_ptr<const ActsGeometryProvider> m_actsGeoProvider;
    dd4hep::Segmentation                         m_segmentation;
    const dd4hep::DDSegmentation::BitFieldCoder*  m_decoder = nullptr;
    dd4hep::VolumeManager                         m_volman;
    double m_sensorHalfX = 8.0; // mm
    double m_sensorHalfY = 8.0; // mm
    bool m_perPlaneFrontBack = false;
    bool m_hasPixX = false;
    bool m_hasPixY = false;
    bool m_hasPixZ = false;
    double m_fallbackMaxNormalMm = 2.0;
    bool m_failOnEmptySensorMap = true;
    bool m_failOnIncompleteSurfaceMap = false;
    bool m_enableTruthSeededChain = true;
    bool m_enableStubSeededChain = true;
    bool m_writeTrackStates = true;
    std::string m_dumpSurfaceMap;
    void dumpSurfaceMap(const std::string& path);

    struct SensorRef {
        std::uint64_t cellID = 0;
        int plane = -1;
        int module = -1;
        int side = -1;
        int sensor = -1;
        int station = -1;
        dd4hep::DetElement detElement;
    };
    std::vector<SensorRef> m_sensorRefs;
    std::unordered_map<std::uint64_t, std::size_t> m_surfaceToSensorIdx;
    std::map<std::pair<int, int>, int> m_moduleToSide;
    std::map<int, int> m_layerToStation;

    int stationOf(int layer) const;

    struct BestSel {
        int index = -1;
        int nStates = -1;
        int nMeasurements = -1;
        int nOutliers = -1;
        int nHoles = -1;
        int nSharedHits = -1;
        int ndf = -1;
        int assocMcIndex = -1;
        std::uint32_t assocMcCollectionID = 0;
        int objectIndex = -1;
        std::uint32_t objectCollectionID = 0;
        int seedIndex = -1;
        std::uint32_t seedCollectionID = 0;
        int identityValid = 0;
        int pdg = 0;
        int charge = 0;
        int momentumResolved = 0;
        double p = 0.0;
        double pT = 0.0;
        double deltaP = 0.0;
        double deltaPT = 0.0;
        double theta = 0.0;
        double phi = 0.0;
        double chi2 = 0.0;
        double assocWeight = 0.0;
        double pullQOverP = 0.0;
        double pullTheta = 0.0;
        double pullPhi = 0.0;
        void reset(double nan);
    };

    // trk_* = B0TrackerCKFTruthSeeded*; ckf_trk_* = B0TrackerCKF*.
    struct TrackChain {
        std::vector<double> p, pT, delta_p, delta_pT, px, py, pz, theta, phi;
        std::vector<double> qOverP, time;
        std::vector<int>    momentum_resolved;
        std::vector<double> loc0, loc1;
        std::vector<double> sigma_loc0, sigma_loc1, sigma_phi, sigma_theta, sigma_qOverP, sigma_time;
        std::vector<double> pull_qOverP, pull_theta, pull_phi;
        std::vector<double> chi2;
        std::vector<int>    ndf;
        std::vector<int>    index, charge, type, pdg;
        std::vector<int>    object_index, seed_index, identity_valid;
        std::vector<std::uint32_t> object_collectionID, seed_collectionID;
        std::array<std::vector<double>, 21> cov_upper;
        std::vector<int>    nStates, nMeasurements, nOutliers, nHoles, nSharedHits;
        std::vector<int>    assoc_mcIndex;
        std::vector<std::uint32_t> assoc_mcCollectionID;
        std::vector<double> assoc_weight;

        std::vector<int>    state_track_index, state_index, state_acts_index, state_type, state_pdg;
        std::vector<int>    state_parent_seed_index, state_parent_track_index, state_parent_identity_valid;
        std::vector<std::uint32_t> state_parent_seed_collectionID, state_parent_track_collectionID;
        std::vector<int>    state_estimate_kind;
        std::vector<int>    state_mapping_method;
        std::vector<std::uint64_t> state_surface;
        std::vector<double> state_loc0, state_loc1;
        std::vector<double> state_meas_loc0, state_meas_loc1, state_resid_loc0, state_resid_loc1;
        std::vector<int>    state_meas_dim, state_proj_index0, state_proj_index1;
        std::vector<double> state_meas_cov00, state_meas_cov01, state_meas_cov11;
        std::vector<double> state_pred_loc0, state_pred_loc1, state_pred_theta, state_pred_phi;
        std::vector<double> state_pred_qOverP, state_pred_time;
        std::array<std::vector<double>, 21> state_pred_cov_upper;
        std::vector<double> state_innov0, state_innov1, state_innov_pull0, state_innov_pull1;
        std::vector<double> state_innov_chi2;
        std::vector<double> x_on_plane, y_on_plane, z_on_plane;
        std::vector<double> aclgad_xPix, aclgad_yPix, aclgad_zPix;
        std::vector<double> aclgad_dx, aclgad_dy, aclgad_dz;
        std::vector<int>    aclgad_pixX, aclgad_pixY, aclgad_pixZ;
        std::vector<int>    aclgad_plane, aclgad_module, aclgad_side, aclgad_sensor, aclgad_station;
        std::vector<std::uint64_t> aclgad_cellID;
        std::vector<double> state_theta, state_phi, state_qOverP, state_time;

        int nMapExact = 0;
        int nMapFallback = 0;
        int nMapFailedPhysics = 0;
        int nMapUnmappedNonPhysics = 0;

        BestSel oracle;
        BestSel truthMatched;
        BestSel recoBest;
        int hasTrack = 0;

        void clear(double nan);
    };

    void bindTrackChain(const std::string& trkPrefix, TrackChain& chain);
    void bindBestSel(const std::string& prefix, BestSel& b);

    // Per-event branch payload. Process() fills this outside the ROOT writer
    // critical section, then swaps it into the stable branch-backed members.
    struct EventBuffers {
        std::uint64_t m_eventNumber = 0;


        // SimTrackerHit-level (cell centers of truth hits, not RecHits).
        std::vector<double> vm_xR, vm_yR, vm_zR;
        std::vector<double> vm_xT, vm_yT, vm_zT;
        std::vector<double> vm_aclgad_xPixT, vm_aclgad_yPixT, vm_aclgad_zPixT;
        std::vector<double> vm_aclgad_dxT, vm_aclgad_dyT, vm_aclgad_dzT;
        std::vector<double> vm_aclgad_xPixR, vm_aclgad_yPixR, vm_aclgad_zPixR;
        std::vector<double> vm_aclgad_dxR, vm_aclgad_dyR, vm_aclgad_dzR;
        std::vector<double> vm_detX, vm_detY, vm_detZ;
        std::vector<int>    vm_plane, vm_station, vm_module, vm_side, vm_sensor, vm_pixX, vm_pixY, vm_pixZ;
        std::vector<int>    vm_aclgad_pixXT, vm_aclgad_pixYT, vm_aclgad_pixXR, vm_aclgad_pixYR;
        std::vector<int>    vm_pdg, vm_status;
        std::vector<int>    vm_isPrimary;              // alias: matches primary_pdg/status
        std::vector<int>    vm_isSelPrimary;           // selected (highest-p) primary only
        std::vector<int>    vm_cellFired;   // cell produced a RawHit (not: this SimHit did)
        std::vector<double> vm_px, vm_py, vm_pz, vm_p, vm_pT;
        std::vector<std::uint64_t> vm_cellID;
        std::vector<int>    vm_mcIndex;
        std::vector<std::uint32_t> vm_mcCollectionID;
        std::vector<double> vm_eDep, vm_time, vm_path;

        std::vector<double> vm_xP, vm_yP, vm_zP, vm_pathP, vm_timeP;
        std::vector<int>    vm_planeP, vm_stationP, vm_moduleP, vm_sideP, vm_sensorP;
        std::vector<int>    vm_pdgP, vm_statusP, vm_isPrimaryP, vm_isSelPrimaryP, vm_mcIndexP;
        std::vector<std::uint32_t> vm_mcCollectionIDP;
        std::vector<double> vm_pxP, vm_pyP, vm_pzP, vm_pP, vm_pTP;

        // first/last SimTrackerHit in time per (mc, layer, side, module).
        std::vector<double> vm_xEntry, vm_yEntry, vm_zEntry, vm_timeEntry;
        std::vector<double> vm_pxEntry, vm_pyEntry, vm_pzEntry, vm_pEntry, vm_pTEntry;
        std::vector<double> vm_xExit,  vm_yExit,  vm_zExit,  vm_timeExit;
        std::vector<double> vm_pxExit, vm_pyExit, vm_pzExit, vm_pExit, vm_pTExit;
        std::vector<int>    vm_sensorEntry, vm_sensorExit;
        std::vector<std::uint64_t> vm_cellIDEntry, vm_cellIDExit;
        std::vector<int>    vm_planeEE, vm_stationEE, vm_moduleEE, vm_sideEE, vm_pdgEE;
        std::vector<int>    vm_isPrimaryEE, vm_isSelPrimaryEE, vm_mcIndexEE, vm_nStepsEE;
        std::vector<int>    vm_statusEE;
        std::vector<std::uint32_t> vm_mcCollectionIDEE;

        // Digitized RawHits / reconstructed RecHits.
        std::vector<std::uint64_t> vm_raw_cellID;
        std::vector<int>    vm_raw_charge, vm_raw_timeStamp;
        std::vector<int>    vm_raw_plane, vm_raw_station, vm_raw_module, vm_raw_side, vm_raw_sensor;
        std::vector<int>    vm_raw_mcIndex;
        std::vector<std::uint32_t> vm_raw_mcCollectionID;
        std::vector<int>    vm_raw_nContribSim, vm_raw_nContribMc, vm_raw_mixedCell;
        std::vector<double> vm_raw_dominantFrac, vm_raw_totalEdep;

        std::vector<double> vm_rec_x, vm_rec_y, vm_rec_z;
        std::vector<double> vm_rec_covxx, vm_rec_covyy, vm_rec_covzz;
        std::vector<double> vm_rec_time, vm_rec_time_err, vm_rec_edep, vm_rec_edep_err;
        std::vector<std::uint64_t> vm_rec_cellID;
        std::vector<int>    vm_rec_plane, vm_rec_station, vm_rec_module, vm_rec_side, vm_rec_sensor;
        std::vector<int>    vm_rec_pixX, vm_rec_pixY, vm_rec_pixZ;
        std::vector<int>    vm_rec_mcIndex;
        std::vector<std::uint32_t> vm_rec_mcCollectionID;
        std::vector<int>    vm_rec_nContribSim, vm_rec_nContribMc, vm_rec_mixedCell;
        std::vector<double> vm_rec_dominantFrac, vm_rec_totalEdep;

        // Seeds.
        std::vector<double> vm_seed_quality, vm_seed_p, vm_seed_qOverP, vm_seed_theta, vm_seed_phi;
        std::vector<double> vm_seed_loc0, vm_seed_loc1;
        std::vector<double> vm_seed_sigma_qOverP, vm_seed_sigma_theta, vm_seed_sigma_phi;
        std::vector<int>    vm_seed_nHits, vm_seed_charge, vm_seed_momentum_resolved, vm_seed_became_track;
        std::vector<int>    vm_seed_made_unfiltered_track, vm_seed_survived_ambiguity;
        std::vector<int>    vm_seed_n_unfiltered_tracks, vm_seed_n_filtered_tracks;
        std::vector<int>    vm_seed_assoc_mcIndex;
        std::vector<std::uint32_t> vm_seed_assoc_mcCollectionID;
        std::vector<double> vm_seed_assoc_weight;
        std::vector<double> vm_truth_seed_quality, vm_truth_seed_p, vm_truth_seed_qOverP;
        std::vector<double> vm_truth_seed_theta, vm_truth_seed_phi, vm_truth_seed_loc0, vm_truth_seed_loc1;
        std::vector<double> vm_truth_seed_sigma_qOverP, vm_truth_seed_sigma_theta, vm_truth_seed_sigma_phi;
        std::vector<int>    vm_truth_seed_nHits, vm_truth_seed_charge, vm_truth_seed_momentum_resolved;
        std::vector<int>    vm_truth_seed_became_track;
        std::vector<int>    vm_truth_seed_made_unfiltered_track, vm_truth_seed_survived_ambiguity;
        std::vector<int>    vm_truth_seed_n_unfiltered_tracks, vm_truth_seed_n_filtered_tracks;

        // Per-stub-seed CKF failure diagnostics (schema 4). Parallel to the
        // stub-seed input order (including null entries); use
        // ckfdiag_seed_index (PODIO ObjectID) for a robust join.
        std::vector<int>    vm_ckfdiag_seed_index, vm_ckfdiag_seed_n_stations;
        std::vector<std::uint32_t> vm_ckfdiag_seed_collectionID;
        std::vector<double> vm_ckfdiag_truth_p, vm_ckfdiag_truth_theta, vm_ckfdiag_truth_phi;
        std::vector<double> vm_ckfdiag_truth_thscat_mrad, vm_ckfdiag_seed_assoc_weight;
        std::vector<int>    vm_ckfdiag_n_candidates, vm_ckfdiag_n_accepted, vm_ckfdiag_stage;
        std::vector<int>    vm_ckfdiag_find_err_class, vm_ckfdiag_find_err_value;
        std::vector<int>    vm_ckfdiag_best_status;
        std::vector<int>    vm_ckfdiag_best_n_states, vm_ckfdiag_best_n_meas;
        std::vector<int>    vm_ckfdiag_best_n_holes, vm_ckfdiag_best_n_outliers;
        std::vector<int>    vm_ckfdiag_best_last_station, vm_ckfdiag_best_first_hole_station;
        std::vector<int>    vm_ckfdiag_best_station_mask;
        std::vector<std::vector<int>> vm_ckfdiag_cand_per_station;
        int m_ckfdiagMaxStation = -1;

        TrackChain m_ts;
        TrackChain m_ckf;

        std::vector<double> beam_px, beam_py, beam_pz, beam_p, beam_pT;
        std::vector<int>    beam_pdg;
        std::vector<double> m_primaryPx, m_primaryPy, m_primaryPz, m_primaryP, m_primaryPT;
        std::vector<double> m_primaryCharge;
        std::vector<int>    m_primaryPdgOut, m_primaryStatusOut, m_primaryMcIndex;
        std::vector<std::uint32_t> m_primaryMcCollectionID;
        std::vector<double> m_genPpx, m_genPpy, m_genPpz, m_genPp, m_genPpT;
        std::vector<double> m_genBeamPx, m_genBeamPy, m_genBeamPz, m_genBeamP, m_genBeamPT;
        std::vector<double> m_genBeamPPx, m_genBeamPPy, m_genBeamPPz, m_genBeamPP, m_genBeamPPT;

        int m_selPrimaryMcIndex = -1;
        std::uint32_t m_selPrimaryMcCollectionID = 0;
        double m_selPrimaryPx = 0.0;
        double m_selPrimaryPy = 0.0;
        double m_selPrimaryPz = 0.0;
        double m_selPrimaryP = 0.0;
        double m_selPrimaryPT = 0.0;
        double m_selPrimaryThscatMrad = 0.0;
        double m_selPrimaryCharge = 0.0;
        int m_nStationsPrimary = 0; // selected primary only
        int m_nSelectedPrimaryMeasurements = -1;
        int m_nMeasurementStationsSelectedPrimary = -1;
        int m_selPrimaryMeasurementReconstructable = -1;
        int m_selPrimaryHasSeed = 0;
        int m_selPrimaryHasUnfilteredTrack = 0;
        int m_selPrimaryHasFilteredTrack = 0;
        int m_selPrimaryHasTruthMatchedTrack = 0;

        int m_nSimHits = 0;
        int m_nRawHits = 0;
        int m_nRecHits = 0;
        int m_nMeasurements = 0;
        int m_nTruthSeeds = 0;
        int m_nStubSeeds = 0;
        int m_nTsUnfiltered = 0;
        int m_nTsFiltered = 0;
        int m_nCkfUnfiltered = 0;
        int m_nCkfFiltered = 0;
        int m_nSimHitsUnresolvedCellID = 0;
        int m_nMissingMcRelation = 0;
        int m_nSensorMapExact = 0;
        int m_nSensorMapFallback = 0;
        int m_nSensorMapFailed = 0;
        int m_nPixelSnapFailed = 0;

        // Optional-input presence. Distinguishes "the factory is not registered"
        // from "the collection is genuinely empty"; required inputs throw instead.
        bool m_hasRawAssocs = false;
        bool m_hasStubSeeds = false;
        bool m_hasTruthSeeds = false;
        bool m_hasTsTrackParams = false;
        bool m_hasTsTrajectories = false;
        bool m_hasTsTrajectoriesUnfiltered = false;
        bool m_hasTsTracks = false;
        bool m_hasTsAssocs = false;
        bool m_hasTsActsStates = false;
        bool m_hasTsActsTracks = false;
        bool m_hasTsTracksUnfiltered = false;
        bool m_hasCkfTrackParams = false;
        bool m_hasCkfTrajectories = false;
        bool m_hasCkfTrajectoriesUnfiltered = false;
        bool m_hasCkfTracks = false;
        bool m_hasCkfAssocs = false;
        bool m_hasCkfActsStates = false;
        bool m_hasCkfActsTracks = false;
        bool m_hasCkfTracksUnfiltered = false;
        bool m_hasCkfAssocsUnfiltered = false;
        bool m_hasCkfActsStatesUnfiltered = false;
        bool m_hasCkfActsTracksUnfiltered = false;
    };

    TTree* m_tree = nullptr;

    int m_schemaVersion = b0trk::kSchemaVersion;
    std::string m_geometryName;
    std::string m_detectorPath;

    std::uint64_t m_eventNumber = 0;

    int m_primaryPdg = 2212;
    int m_primaryStatus = 1;
    int m_minMeasurementStations = 3;

    // SimTrackerHit-level (cell centers of truth hits, not RecHits).
    std::vector<double> vm_xR, vm_yR, vm_zR;
    std::vector<double> vm_xT, vm_yT, vm_zT;
    std::vector<double> vm_aclgad_xPixT, vm_aclgad_yPixT, vm_aclgad_zPixT;
    std::vector<double> vm_aclgad_dxT, vm_aclgad_dyT, vm_aclgad_dzT;
    std::vector<double> vm_aclgad_xPixR, vm_aclgad_yPixR, vm_aclgad_zPixR;
    std::vector<double> vm_aclgad_dxR, vm_aclgad_dyR, vm_aclgad_dzR;
    std::vector<double> vm_detX, vm_detY, vm_detZ;
    std::vector<int>    vm_plane, vm_station, vm_module, vm_side, vm_sensor, vm_pixX, vm_pixY, vm_pixZ;
    std::vector<int>    vm_aclgad_pixXT, vm_aclgad_pixYT, vm_aclgad_pixXR, vm_aclgad_pixYR;
    std::vector<int>    vm_pdg, vm_status;
    std::vector<int>    vm_isPrimary;              // alias: matches primary_pdg/status
    std::vector<int>    vm_isSelPrimary;           // selected (highest-p) primary only
    std::vector<int>    vm_cellFired;   // cell produced a RawHit (not: this SimHit did)
    std::vector<double> vm_px, vm_py, vm_pz, vm_p, vm_pT;
    std::vector<std::uint64_t> vm_cellID;
    std::vector<int>    vm_mcIndex;
    std::vector<std::uint32_t> vm_mcCollectionID;
    std::vector<double> vm_eDep, vm_time, vm_path;

    std::vector<double> vm_xP, vm_yP, vm_zP, vm_pathP, vm_timeP;
    std::vector<int>    vm_planeP, vm_stationP, vm_moduleP, vm_sideP, vm_sensorP;
    std::vector<int>    vm_pdgP, vm_statusP, vm_isPrimaryP, vm_isSelPrimaryP, vm_mcIndexP;
    std::vector<std::uint32_t> vm_mcCollectionIDP;
    std::vector<double> vm_pxP, vm_pyP, vm_pzP, vm_pP, vm_pTP;

    // first/last SimTrackerHit in time per (mc, layer, side, module).
    std::vector<double> vm_xEntry, vm_yEntry, vm_zEntry, vm_timeEntry;
    std::vector<double> vm_pxEntry, vm_pyEntry, vm_pzEntry, vm_pEntry, vm_pTEntry;
    std::vector<double> vm_xExit,  vm_yExit,  vm_zExit,  vm_timeExit;
    std::vector<double> vm_pxExit, vm_pyExit, vm_pzExit, vm_pExit, vm_pTExit;
    std::vector<int>    vm_sensorEntry, vm_sensorExit;
    std::vector<std::uint64_t> vm_cellIDEntry, vm_cellIDExit;
    std::vector<int>    vm_planeEE, vm_stationEE, vm_moduleEE, vm_sideEE, vm_pdgEE;
    std::vector<int>    vm_isPrimaryEE, vm_isSelPrimaryEE, vm_mcIndexEE, vm_nStepsEE;
    std::vector<int>    vm_statusEE;
    std::vector<std::uint32_t> vm_mcCollectionIDEE;

    // Digitized RawHits / reconstructed RecHits.
    std::vector<std::uint64_t> vm_raw_cellID;
    std::vector<int>    vm_raw_charge, vm_raw_timeStamp;
    std::vector<int>    vm_raw_plane, vm_raw_station, vm_raw_module, vm_raw_side, vm_raw_sensor;
    std::vector<int>    vm_raw_mcIndex;
    std::vector<std::uint32_t> vm_raw_mcCollectionID;
    std::vector<int>    vm_raw_nContribSim, vm_raw_nContribMc, vm_raw_mixedCell;
    std::vector<double> vm_raw_dominantFrac, vm_raw_totalEdep;

    std::vector<double> vm_rec_x, vm_rec_y, vm_rec_z;
    std::vector<double> vm_rec_covxx, vm_rec_covyy, vm_rec_covzz;
    std::vector<double> vm_rec_time, vm_rec_time_err, vm_rec_edep, vm_rec_edep_err;
    std::vector<std::uint64_t> vm_rec_cellID;
    std::vector<int>    vm_rec_plane, vm_rec_station, vm_rec_module, vm_rec_side, vm_rec_sensor;
    std::vector<int>    vm_rec_pixX, vm_rec_pixY, vm_rec_pixZ;
    std::vector<int>    vm_rec_mcIndex;
    std::vector<std::uint32_t> vm_rec_mcCollectionID;
    std::vector<int>    vm_rec_nContribSim, vm_rec_nContribMc, vm_rec_mixedCell;
    std::vector<double> vm_rec_dominantFrac, vm_rec_totalEdep;

    // Seeds.
    std::vector<double> vm_seed_quality, vm_seed_p, vm_seed_qOverP, vm_seed_theta, vm_seed_phi;
    std::vector<double> vm_seed_loc0, vm_seed_loc1;
    std::vector<double> vm_seed_sigma_qOverP, vm_seed_sigma_theta, vm_seed_sigma_phi;
    std::vector<int>    vm_seed_nHits, vm_seed_charge, vm_seed_momentum_resolved, vm_seed_became_track;
    std::vector<int>    vm_seed_made_unfiltered_track, vm_seed_survived_ambiguity;
    std::vector<int>    vm_seed_n_unfiltered_tracks, vm_seed_n_filtered_tracks;
    std::vector<int>    vm_seed_assoc_mcIndex;
    std::vector<std::uint32_t> vm_seed_assoc_mcCollectionID;
    std::vector<double> vm_seed_assoc_weight;
    std::vector<double> vm_truth_seed_quality, vm_truth_seed_p, vm_truth_seed_qOverP;
    std::vector<double> vm_truth_seed_theta, vm_truth_seed_phi, vm_truth_seed_loc0, vm_truth_seed_loc1;
    std::vector<double> vm_truth_seed_sigma_qOverP, vm_truth_seed_sigma_theta, vm_truth_seed_sigma_phi;
    std::vector<int>    vm_truth_seed_nHits, vm_truth_seed_charge, vm_truth_seed_momentum_resolved;
    std::vector<int>    vm_truth_seed_became_track;
    std::vector<int>    vm_truth_seed_made_unfiltered_track, vm_truth_seed_survived_ambiguity;
    std::vector<int>    vm_truth_seed_n_unfiltered_tracks, vm_truth_seed_n_filtered_tracks;

    std::vector<int>    vm_ckfdiag_seed_index, vm_ckfdiag_seed_n_stations;
    std::vector<std::uint32_t> vm_ckfdiag_seed_collectionID;
    std::vector<double> vm_ckfdiag_truth_p, vm_ckfdiag_truth_theta, vm_ckfdiag_truth_phi;
    std::vector<double> vm_ckfdiag_truth_thscat_mrad, vm_ckfdiag_seed_assoc_weight;
    std::vector<int>    vm_ckfdiag_n_candidates, vm_ckfdiag_n_accepted, vm_ckfdiag_stage;
    std::vector<int>    vm_ckfdiag_find_err_class, vm_ckfdiag_find_err_value;
    std::vector<int>    vm_ckfdiag_best_status;
    std::vector<int>    vm_ckfdiag_best_n_states, vm_ckfdiag_best_n_meas;
    std::vector<int>    vm_ckfdiag_best_n_holes, vm_ckfdiag_best_n_outliers;
    std::vector<int>    vm_ckfdiag_best_last_station, vm_ckfdiag_best_first_hole_station;
    std::vector<int>    vm_ckfdiag_best_station_mask;
    std::vector<std::vector<int>> vm_ckfdiag_cand_per_station;
    int m_ckfdiagMaxStation = -1;

    TrackChain m_ts;
    TrackChain m_ckf;

    std::vector<double> beam_px, beam_py, beam_pz, beam_p, beam_pT;
    std::vector<int>    beam_pdg;
    std::vector<double> m_primaryPx, m_primaryPy, m_primaryPz, m_primaryP, m_primaryPT;
    std::vector<double> m_primaryCharge;
    std::vector<int>    m_primaryPdgOut, m_primaryStatusOut, m_primaryMcIndex;
    std::vector<std::uint32_t> m_primaryMcCollectionID;
    std::vector<double> m_genPpx, m_genPpy, m_genPpz, m_genPp, m_genPpT;
    std::vector<double> m_genBeamPx, m_genBeamPy, m_genBeamPz, m_genBeamP, m_genBeamPT;
    std::vector<double> m_genBeamPPx, m_genBeamPPy, m_genBeamPPz, m_genBeamPP, m_genBeamPPT;

    int m_selPrimaryMcIndex = -1;
    std::uint32_t m_selPrimaryMcCollectionID = 0;
    double m_selPrimaryPx = 0.0;
    double m_selPrimaryPy = 0.0;
    double m_selPrimaryPz = 0.0;
    double m_selPrimaryP = 0.0;
    double m_selPrimaryPT = 0.0;
    double m_selPrimaryThscatMrad = 0.0;
    double m_selPrimaryCharge = 0.0;
    int m_nStationsPrimary = 0; // selected primary only
    int m_nSelectedPrimaryMeasurements = -1;
    int m_nMeasurementStationsSelectedPrimary = -1;
    int m_selPrimaryMeasurementReconstructable = -1;
    int m_selPrimaryHasSeed = 0;
    int m_selPrimaryHasUnfilteredTrack = 0;
    int m_selPrimaryHasFilteredTrack = 0;
    int m_selPrimaryHasTruthMatchedTrack = 0;

    int m_nSimHits = 0;
    int m_nRawHits = 0;
    int m_nRecHits = 0;
    int m_nMeasurements = 0;
    int m_nTruthSeeds = 0;
    int m_nStubSeeds = 0;
    int m_nTsUnfiltered = 0;
    int m_nTsFiltered = 0;
    int m_nCkfUnfiltered = 0;
    int m_nCkfFiltered = 0;
    int m_nSimHitsUnresolvedCellID = 0;
    int m_nMissingMcRelation = 0;
    int m_nSensorMapExact = 0;
    int m_nSensorMapFallback = 0;
    int m_nSensorMapFailed = 0;
    int m_nPixelSnapFailed = 0;

    // Optional-input presence. Distinguishes "the factory is not registered"
    // from "the collection is genuinely empty"; required inputs throw instead.
    bool m_hasRawAssocs = false;
    bool m_hasStubSeeds = false;
    bool m_hasTruthSeeds = false;
    bool m_hasTsTrackParams = false;
    bool m_hasTsTrajectories = false;
    bool m_hasTsTrajectoriesUnfiltered = false;
    bool m_hasTsTracks = false;
    bool m_hasTsAssocs = false;
    bool m_hasTsActsStates = false;
    bool m_hasTsActsTracks = false;
    bool m_hasTsTracksUnfiltered = false;
    bool m_hasCkfTrackParams = false;
    bool m_hasCkfTrajectories = false;
    bool m_hasCkfTrajectoriesUnfiltered = false;
    bool m_hasCkfTracks = false;
    bool m_hasCkfAssocs = false;
    bool m_hasCkfActsStates = false;
    bool m_hasCkfActsTracks = false;
    bool m_hasCkfTracksUnfiltered = false;
    bool m_hasCkfAssocsUnfiltered = false;
    bool m_hasCkfActsStatesUnfiltered = false;
    bool m_hasCkfActsTracksUnfiltered = false;
};

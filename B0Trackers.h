#pragma once

#include <memory>
#include <mutex>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <JANA/JEventProcessor.h>

#include <DD4hep/DetElement.h>
#include <DD4hep/VolumeManager.h>
#include <DD4hep/Segmentations.h>
#include <DDSegmentation/BitFieldCoder.h>

#include <services/geometry/dd4hep/DD4hep_service.h>

class TTree;
class ACTSGeo_service;
class ActsGeometryProvider;

class B0Trackers : public JEventProcessor {
public:
    void Init() override;
    void Process(const std::shared_ptr<const JEvent>& event) override;
    void Finish() override;

private:
    // Serializes the entire Process() body across JANA worker threads:
    // protects both TTree::Fill() and the per-event member vectors below.
    std::mutex m_fillMutex;

    // Cached geometry handles (set once in Init).
    std::shared_ptr<DD4hep_service> m_geoSvc;
    std::shared_ptr<ACTSGeo_service> m_actsGeoSvc;
    std::shared_ptr<const ActsGeometryProvider> m_actsGeoProvider;
    dd4hep::Segmentation                         m_segmentation;
    const dd4hep::DDSegmentation::BitFieldCoder*  m_decoder = nullptr;
    dd4hep::VolumeManager                         m_volman;
    double m_sensorHalfX = 8.0; // mm; overwritten from geometry constants when available.
    double m_sensorHalfY = 8.0; // mm; overwritten from geometry constants when available.

    struct SensorRef {
        std::uint64_t cellID = 0;
        int plane = -1;
        int module = -1;
        int side = -1;
        int sensor = -1;
        // Cached from the Init tree walk so per-hit/per-state transforms skip
        // the VolumeManager lookup.
        dd4hep::DetElement detElement;
    };
    std::vector<SensorRef> m_sensorRefs;

    // Acts surface geometryId value -> index into m_sensorRefs. Built in Init by
    // inverting ActsGeometryProvider::surfaceMap(); gives the exact sensor for a
    // track state's reference surface (the nearest-sensor scan is only a fallback).
    std::unordered_map<std::uint64_t, std::size_t> m_surfaceToSensorIdx;

    // (layer, module-volID) -> side: 1=front, 0=back. Built once in Init by
    // walking the B0Tracker DetElement tree. Side comes from the layer
    // DetElement name ("..._front_P"/"..._back_P"; each plane is its own cellID
    // layer in the per-plane geometry). For older layouts whose layer names
    // don't carry a side, fall back to the sign of the module translation
    // (modules at +/-ModuleOffsetFromSupport in the station frame). The key
    // uses the physVolIDs that the cellID decoder reports so it works across
    // dev / official B0 geometries (their DetElement::id() conventions differ).
    std::map<std::pair<int, int>, int> m_moduleToSide;

    // Compact B0 layer id -> physical station. layer = 2*(station-1)+(1=back|2=front).
    static int stationFromLayer(int layer) {
        return (layer > 0) ? (layer + 1) / 2 : -1;
    }

    // One reconstructed B0 CKF chain. `trk_*` on the tree is truth-seeded
    // (B0TrackerCKFTruthSeeded*); `ckf_trk_*` is stub-seeded (B0TrackerCKF*).
    struct TrackChain {
        std::vector<double> p, pT, delta_p, delta_pT, px, py, pz, theta, phi;
        std::vector<double> qOverP, time;
        std::vector<double> loc0, loc1;
        std::vector<double> sigma_loc0, sigma_loc1, sigma_phi, sigma_theta, sigma_qOverP, sigma_time;
        std::vector<double> chi2;
        std::vector<int>    ndf;
        std::vector<int>    index, charge, type, pdg;
        std::vector<int>    nStates, nMeasurements, nOutliers, nHoles, nSharedHits;
        std::vector<int>    assoc_mcIndex;
        std::vector<double> assoc_weight;
        std::vector<std::uint64_t> surface;

        std::vector<int>    state_track_index, state_index, state_acts_index, state_type, state_pdg;
        std::vector<std::uint64_t> state_surface;
        std::vector<double> state_loc0, state_loc1;
        std::vector<double> x_on_plane, y_on_plane, z_on_plane;
        std::vector<double> aclgad_xPix, aclgad_yPix, aclgad_zPix;
        std::vector<double> aclgad_dx, aclgad_dy, aclgad_dz;
        std::vector<int>    aclgad_pixX, aclgad_pixY;
        std::vector<int>    aclgad_plane, aclgad_module, aclgad_side, aclgad_sensor, aclgad_station;
        std::vector<std::uint64_t> aclgad_cellID;
        std::vector<double> state_theta, state_phi, state_qOverP, state_time;

        int bestIndex = -1;
        int bestNStates = -1;
        int bestNMeasurements = -1;
        int bestNOutliers = -1;
        int bestNHoles = -1;
        int bestNSharedHits = -1;
        int bestNdf = -1;
        int bestAssocMcIndex = -1;
        double bestP = 0.0;
        double bestPT = 0.0;
        double bestDeltaP = 0.0;
        double bestDeltaPT = 0.0;
        double bestTheta = 0.0;
        double bestPhi = 0.0;
        double bestChi2 = 0.0;
        double bestAssocWeight = 0.0;
        int hasTrack = 0;

        void clear(double nan);
    };

    void bindTrackChain(const std::string& trkPrefix, const std::string& bestPrefix, TrackChain& chain);

    // ROOT output (owned by the file via TDirectory::cd() in Init).
    TTree* m_tree = nullptr;

    // Event identifier for this TTree entry.
    std::uint64_t m_eventNumber = 0;

    // Primary-particle selector. Defaults to the scattered final-state proton.
    int m_primaryPdg = 2212;
    int m_primaryStatus = 1;

    // Per-event hit vectors.
    std::vector<double> vm_xR, vm_yR, vm_zR;        // readout cell-center, mm, global
    std::vector<double> vm_xT, vm_yT, vm_zT;        // truth step position, mm, global
    std::vector<double> vm_aclgad_xPixT, vm_aclgad_yPixT, vm_aclgad_zPixT;
    std::vector<double> vm_aclgad_dxT, vm_aclgad_dyT, vm_aclgad_dzT;
    std::vector<double> vm_aclgad_xPixR, vm_aclgad_yPixR, vm_aclgad_zPixR;
    std::vector<double> vm_aclgad_dxR, vm_aclgad_dyR, vm_aclgad_dzR;
    std::vector<double> vm_detX, vm_detY, vm_detZ;  // local sensor frame, mm
    std::vector<int>    vm_plane, vm_station, vm_module, vm_side, vm_sensor, vm_pixX, vm_pixY, vm_pixZ;
    std::vector<int>    vm_aclgad_pixXT, vm_aclgad_pixYT, vm_aclgad_pixXR, vm_aclgad_pixYR;
    std::vector<int>    vm_pdg, vm_status;
    std::vector<int>    vm_isPrimary;
    std::vector<double> vm_px, vm_py, vm_pz, vm_p, vm_pT;
    std::vector<std::uint64_t> vm_cellID;
    std::vector<int>    vm_mcIndex;
    std::vector<std::uint32_t> vm_mcCollectionID;
    std::vector<double> vm_eDep, vm_time, vm_path;

    // One crossing per MC particle per sensitive volume (layer, side, module, sensor).
    // Each *P entry is the first-entry Geant4 step (smallest time) into that sensor.
    std::vector<double> vm_xP, vm_yP, vm_zP, vm_pathP, vm_timeP; // mm, mm, mm, step path length, ns
    std::vector<int>    vm_planeP, vm_stationP, vm_moduleP, vm_sideP, vm_sensorP;
    std::vector<int>    vm_pdgP, vm_statusP, vm_isPrimaryP, vm_mcIndexP;
    std::vector<std::uint32_t> vm_mcCollectionIDP;
    std::vector<double> vm_pxP, vm_pyP, vm_pzP, vm_pP, vm_pTP;

    // Entry/exit summary: one row per (mc particle, disk, side, module). Side is
    // {0=back, 1=front}. Entry = smallest-time hit in this group; exit = largest-time.
    // nStepsEE counts how many SimTrackerHits contributed (entry==exit if 1).
    std::vector<double> vm_xEntry, vm_yEntry, vm_zEntry, vm_timeEntry;
    std::vector<double> vm_pxEntry, vm_pyEntry, vm_pzEntry, vm_pEntry, vm_pTEntry;
    std::vector<double> vm_xExit,  vm_yExit,  vm_zExit,  vm_timeExit;
    std::vector<double> vm_pxExit, vm_pyExit, vm_pzExit, vm_pExit, vm_pTExit;
    std::vector<int>    vm_sensorEntry, vm_sensorExit;
    std::vector<std::uint64_t> vm_cellIDEntry, vm_cellIDExit;
    std::vector<int>    vm_planeEE, vm_stationEE, vm_moduleEE, vm_sideEE, vm_pdgEE, vm_isPrimaryEE, vm_mcIndexEE, vm_nStepsEE;
    std::vector<int>    vm_statusEE;   // linked MCParticle::generatorStatus
    std::vector<std::uint32_t> vm_mcCollectionIDEE;

    // trk_* / best_trk_* = B0TrackerCKFTruthSeeded (truth-seeded CKF, filtered).
    // ckf_trk_* / ckf_best_trk_* = B0TrackerCKF (stub-seeded CKF, filtered).
    TrackChain m_ts;
    TrackChain m_ckf;

    // Per-event MC truth.
    std::vector<double> beam_px, beam_py, beam_pz, beam_p, beam_pT;
    std::vector<int>    beam_pdg;
    std::vector<double> m_primaryPx, m_primaryPy, m_primaryPz, m_primaryP, m_primaryPT;
    std::vector<int>    m_primaryPdgOut, m_primaryStatusOut, m_primaryMcIndex;
    std::vector<std::uint32_t> m_primaryMcCollectionID;
    std::vector<double> m_genPpx, m_genPpy, m_genPpz, m_genPp, m_genPpT;          // beam proton (status 4)
    std::vector<double> m_genBeamPx, m_genBeamPy, m_genBeamPz, m_genBeamP, m_genBeamPT;     // scattered electron (status 1)
    std::vector<double> m_genBeamPPx, m_genBeamPPy, m_genBeamPPz, m_genBeamPP, m_genBeamPPT; // scattered proton (status 1)

    // Highest-p status-1 proton (the one used for best-track |Δp| matching).
    int m_selPrimaryMcIndex = -1;
    std::uint32_t m_selPrimaryMcCollectionID = 0;
    double m_selPrimaryPx = 0.0;
    double m_selPrimaryPy = 0.0;
    double m_selPrimaryPz = 0.0;
    double m_selPrimaryP = 0.0;
    double m_selPrimaryPT = 0.0;
    double m_selPrimaryThscatMrad = 0.0; // beam-frame scattering angle vs status-4 proton
    int m_nStationsPrimary = 0;          // unique physical stations hit by selected primary
};

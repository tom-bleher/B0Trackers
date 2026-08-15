#include "B0Trackers.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>

#include <JANA/Services/JGlobalRootLock.h>

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
#include <edm4eic/Track.h>
#include <edm4eic/TrackParameters.h>
#include <edm4eic/Trajectory.h>

#include <DD4hep/Objects.h>

#include <services/rootfile/RootFile_service.h>
#include <algorithms/tracking/ActsGeometryProvider.h>
#include <services/geometry/acts/ACTSGeo_service.h>

// B0Trackers/hits branch map (p13 names kept):
//   trk_* / best_trk_*     = B0TrackerCKFTruthSeeded* (truth-seeded CKF, filtered)
//   ckf_trk_* / ckf_best_* = B0TrackerCKF*            (stub-seeded CKF, filtered)
//   station*               = physical disk 1-4 from layer id (layer+1)/2
//   sel_primary_*          = highest-p status-1 proton (used for |Δp| matching)
//   beam_proton_*          = status-4 proton   (alias of genP*)
//   scattered_e_*          = status-1 electron (alias of genBeamP*; not the beam)
//   scattered_p_*          = status-1 proton   (alias of genBeamPP*)
//   n_stations_primary     = unique stations with a primary (status-1 proton) truth hit
extern "C" {
    void InitPlugin(JApplication* app) {
        InitJANAPlugin(app);
        app->Add(new B0Trackers);
    }
}

void B0Trackers::Init() {
    auto* app = GetApplication();

    app->SetDefaultParameter("B0Trackers:primary_pdg", m_primaryPdg,
                             "PDG code used to tag the selected primary particle");
    app->SetDefaultParameter("B0Trackers:primary_status", m_primaryStatus,
                             "Generator status used to tag the selected primary particle");

    auto rootLock = app->GetService<JGlobalRootLock>();
    rootLock->acquire_write_lock();

    auto rf_svc = app->GetService<RootFile_service>();
    TFile* outfile = rf_svc->GetHistFile()->GetFile();
    outfile->mkdir("B0Trackers")->cd();

    m_tree = new TTree("hits", "Truth vs read-out");

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
    m_tree->Branch("plane",     &vm_plane);          // ACTS layer id 1-8
    m_tree->Branch("station",   &vm_station);        // physical disk 1-4; (layer+1)/2
    m_tree->Branch("module",    &vm_module);
    m_tree->Branch("side",      &vm_side);           // -1=unknown, 0=back, 1=front
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
    m_tree->Branch("xP",        &vm_xP);
    m_tree->Branch("yP",        &vm_yP);
    m_tree->Branch("zP",        &vm_zP);
    m_tree->Branch("pathP",     &vm_pathP);
    m_tree->Branch("timeP",     &vm_timeP);
    m_tree->Branch("planeP",    &vm_planeP);
    m_tree->Branch("stationP",  &vm_stationP);
    m_tree->Branch("moduleP",   &vm_moduleP);
    m_tree->Branch("sideP",     &vm_sideP);          // -1=unknown, 0=back, 1=front
    m_tree->Branch("sensorP",   &vm_sensorP);
    m_tree->Branch("pdgP",      &vm_pdgP);
    m_tree->Branch("statusP",   &vm_statusP);
    m_tree->Branch("isPrimaryP", &vm_isPrimaryP);
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
    // Historical names kept. Prefer the explicit aliases below.
    m_tree->Branch("genPpx",    &m_genPpx);           // beam proton, status 4
    m_tree->Branch("genPpy",    &m_genPpy);
    m_tree->Branch("genPpz",    &m_genPpz);
    m_tree->Branch("genPp",     &m_genPp);
    m_tree->Branch("genPpT",    &m_genPpT);
    m_tree->Branch("genBeamP",  &m_genBeamP);         // scattered electron, status 1
    m_tree->Branch("genBeamPT", &m_genBeamPT);
    m_tree->Branch("genBeamPx", &m_genBeamPx);
    m_tree->Branch("genBeamPy", &m_genBeamPy);
    m_tree->Branch("genBeamPz", &m_genBeamPz);
    m_tree->Branch("genBeamPP", &m_genBeamPP);        // scattered proton, status 1
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

    // trk_* = B0TrackerCKFTruthSeeded (truth-seeded). ckf_trk_* = B0TrackerCKF (stub-seeded).
    bindTrackChain("trk_", "best_trk_", m_ts);
    bindTrackChain("ckf_trk_", "ckf_best_trk_", m_ckf);

    // Entry/exit summary branches (one row per mc-particle/disk/side/module).
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
    m_tree->Branch("sensorEntry", &vm_sensorEntry);
    m_tree->Branch("sensorExit",  &vm_sensorExit);
    m_tree->Branch("cellIDEntry", &vm_cellIDEntry);
    m_tree->Branch("cellIDExit",  &vm_cellIDExit);
    m_tree->Branch("planeEE",    &vm_planeEE);
    m_tree->Branch("stationEE",  &vm_stationEE);
    m_tree->Branch("moduleEE",   &vm_moduleEE);
    m_tree->Branch("sideEE",     &vm_sideEE);          // 0=back, 1=front
    m_tree->Branch("pdgEE",      &vm_pdgEE);
    m_tree->Branch("mcIndexEE",  &vm_mcIndexEE);
    m_tree->Branch("mcCollectionIDEE", &vm_mcCollectionIDEE);
    m_tree->Branch("nStepsEE",   &vm_nStepsEE);
    m_tree->Branch("statusEE",   &vm_statusEE);   // linked MCParticle::generatorStatus
    m_tree->Branch("isPrimaryEE", &vm_isPrimaryEE);

    // B0TrackerHits readout fields vary by geometry version; decode optional pixel fields safely.
    m_geoSvc = app->GetService<DD4hep_service>();
    m_actsGeoSvc = app->GetService<ACTSGeo_service>();
    m_actsGeoProvider = m_actsGeoSvc ? m_actsGeoSvc->actsGeoProvider() : nullptr;
    auto ro  = m_geoSvc->detector()->readout("B0TrackerHits");
    m_segmentation = ro.segmentation();
    m_decoder = ro.idSpec().decoder();
    m_volman  = m_geoSvc->detector()->volumeManager();
    // Half-extent in mm from the compact constant: dd4hep constants come back in
    // native units (cm), so mm = x10 and half = x0.5, i.e. a factor 5. The main
    // geometry names the constants B0TrackerSensor*, older match_* variants use
    // the bare Sensor* names.
    const auto sensorHalfExtentMm = [this](std::initializer_list<const char*> names,
                                           double fallbackMm) -> double {
        for (const char* name : names) {
            try {
                return 5.0 * m_geoSvc->detector()->constant<double>(name);
            } catch (const std::exception&) {
                // Constant not defined in this geometry; try the next name.
            }
        }
        return fallbackMm;
    };
    m_sensorHalfX = sensorHalfExtentMm({"B0TrackerSensorWidth", "SensorWidth"}, 8.0);
    m_sensorHalfY = sensorHalfExtentMm({"B0TrackerSensorLength", "SensorLength"}, 8.0);

    // Build (layer, module) -> side map by walking the B0Tracker DetElement tree.
    // In the per-plane geometry each front/back plane is its own layer DetElement
    // named "..._front_P"/"..._back_P", and modules sit at z ~ 0 inside it, so the
    // side must come from the layer name. Older layouts placed modules at
    // z = +/-ModuleOffsetFromSupport inside a station layer; keep the sign of the
    // module translation as a fallback for those. Key uses the cellID's physVolID
    // values (not the DetElement id, which has different semantics in dev vs
    // official geometries).
    auto b0Det = m_geoSvc->detector()->detector("B0Tracker");
    if (b0Det.isValid()) {
        const auto findVolID = [](const dd4hep::PlacedVolume& pv, const std::string& name) -> int {
            for (const auto& id : pv.volIDs()) {
                if (id.first == name) return id.second;
            }
            return -1;
        };
        for (const auto& [layerName, layerDE] : b0Det.children()) {
            const auto layer_pv = layerDE.placement();
            if (!layer_pv.isValid()) continue;
            const int layer_id = findVolID(layer_pv, "layer");
            if (layer_id < 0) continue;

            int layerSide = -1;
            if (layerName.find("front") != std::string::npos)     layerSide = 1;
            else if (layerName.find("back") != std::string::npos) layerSide = 0;

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

                    dd4hep::CellID refCellID = 0;
                    try {
                        m_decoder->set(refCellID, "system", b0Det.id());
                        m_decoder->set(refCellID, "layer", layer_id);
                        m_decoder->set(refCellID, "module", module_id);
                        m_decoder->set(refCellID, "sensor", sensor_id);
                        m_decoder->set(refCellID, "x", 0);
                        m_decoder->set(refCellID, "y", 0);
                    } catch (const std::exception&) {
                        continue;
                    }
                    m_sensorRefs.push_back({static_cast<std::uint64_t>(refCellID), layer_id, module_id, side, sensor_id, sensorDE});
                }
            }
        }
    }

    // Invert the ACTS surface map (sensor volumeID -> surface) so a track state's
    // reference surface resolves to its sensor exactly, without a proximity scan.
    // The volumeID keys carry the same physVolID fields we packed into the
    // SensorRef cellIDs, so a plain lookup filters out all non-B0 surfaces.
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

    rootLock->release_lock();
}

void B0Trackers::Process(const std::shared_ptr<const JEvent>& event) {
    // Fetch before locking: event->Get() may run expensive factories (CKF),
    // and JANA already serializes that work internally.
    auto mcparticles = event->Get<edm4hep::MCParticle>("MCParticles");
    auto simHits     = event->Get<edm4hep::SimTrackerHit>("B0TrackerHits");
    auto tsTracks = event->Get<edm4eic::TrackParameters>("B0TrackerCKFTruthSeededTrackParameters");
    auto tsTrajectories = event->Get<edm4eic::Trajectory>("B0TrackerCKFTruthSeededTrajectories");
    auto tsEdmTracks = event->Get<edm4eic::Track>("B0TrackerCKFTruthSeededTracks");
    auto tsAssocs = event->Get<edm4eic::MCRecoTrackParticleAssociation>(
        "B0TrackerCKFTruthSeededTrackAssociations");
    auto tsActsTrackStates =
        event->Get<Acts::ConstVectorMultiTrajectory>("B0TrackerCKFTruthSeededActsTrackStates");
    auto tsActsTracks = event->Get<Acts::ConstVectorTrackContainer>("B0TrackerCKFTruthSeededActsTracks");
    auto ckfTracks = event->Get<edm4eic::TrackParameters>("B0TrackerCKFTrackParameters");
    auto ckfTrajectories = event->Get<edm4eic::Trajectory>("B0TrackerCKFTrajectories");
    auto ckfEdmTracks = event->Get<edm4eic::Track>("B0TrackerCKFTracks");
    auto ckfAssocs = event->Get<edm4eic::MCRecoTrackParticleAssociation>("B0TrackerCKFTrackAssociations");
    auto ckfActsTrackStates =
        event->Get<Acts::ConstVectorMultiTrajectory>("B0TrackerCKFActsTrackStates");
    auto ckfActsTracks = event->Get<Acts::ConstVectorTrackContainer>("B0TrackerCKFActsTracks");

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
    vm_isPrimary.clear();
    vm_px.clear();    vm_py.clear();    vm_pz.clear();   vm_p.clear();    vm_pT.clear();

    vm_xP.clear();      vm_yP.clear();      vm_zP.clear();      vm_pathP.clear();   vm_timeP.clear();
    vm_planeP.clear();  vm_stationP.clear(); vm_moduleP.clear(); vm_sideP.clear(); vm_sensorP.clear();

    vm_xEntry.clear();   vm_yEntry.clear();  vm_zEntry.clear();  vm_timeEntry.clear();
    vm_pxEntry.clear();  vm_pyEntry.clear(); vm_pzEntry.clear(); vm_pEntry.clear(); vm_pTEntry.clear();
    vm_xExit.clear();    vm_yExit.clear();   vm_zExit.clear();   vm_timeExit.clear();
    vm_pxExit.clear();   vm_pyExit.clear();  vm_pzExit.clear();  vm_pExit.clear();  vm_pTExit.clear();
    vm_sensorEntry.clear(); vm_sensorExit.clear();
    vm_cellIDEntry.clear(); vm_cellIDExit.clear();
    vm_planeEE.clear();  vm_stationEE.clear(); vm_moduleEE.clear(); vm_sideEE.clear();  vm_pdgEE.clear();
    vm_isPrimaryEE.clear();
    vm_mcIndexEE.clear(); vm_mcCollectionIDEE.clear(); vm_nStepsEE.clear();
    vm_statusEE.clear();
    vm_pdgP.clear();    vm_statusP.clear(); vm_isPrimaryP.clear(); vm_mcIndexP.clear(); vm_mcCollectionIDP.clear();
    vm_pxP.clear();     vm_pyP.clear();     vm_pzP.clear();     vm_pP.clear();     vm_pTP.clear();

    const double nan = std::numeric_limits<double>::quiet_NaN();
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

    // Key: (mcCollectionID, mcIndex, layer, side, module, sensor) — one *P entry
    // per (particle, sensitive-volume) pair.
    std::map<std::tuple<uint32_t, int, int, int, int, int>, std::size_t> penetrationIndex;

    // Key: (mcCollectionID, mcIndex, layer, side, module) — one entry/exit row per
    // (particle, disk side, module). Tracks the smallest- and largest-time
    // SimTrackerHits contributing to that group so we can report where the particle
    // entered and exited the silicon in each module. Note: a particle that loops
    // back through the same module still gets a single row spanning both passes
    // (nStepsEE counts the steps of every pass).
    std::map<std::tuple<uint32_t, int, int, int, int>, std::size_t> entryExitIndex;
    const auto getFieldOr = [this](std::uint64_t cellID, const char* field, int fallback) -> int {
        try {
            return static_cast<int>(m_decoder->get(cellID, field));
        } catch (...) {
            return fallback;
        }
    };
    struct PixelSnap {
        double x = std::numeric_limits<double>::quiet_NaN();
        double y = std::numeric_limits<double>::quiet_NaN();
        double z = std::numeric_limits<double>::quiet_NaN();
        double dx = std::numeric_limits<double>::quiet_NaN();
        double dy = std::numeric_limits<double>::quiet_NaN();
        double dz = std::numeric_limits<double>::quiet_NaN();
        int pixX = -1;
        int pixY = -1;
        std::uint64_t cellID = 0;
    };
    const auto snapToAclgadPixel = [this, &nan, &getFieldOr](double xMm, double yMm, double zMm,
                                                             std::uint64_t referenceCellID,
                                                             dd4hep::DetElement detElementHint = {}) -> PixelSnap {
        PixelSnap snap;
        if (!std::isfinite(xMm) || !std::isfinite(yMm) || !std::isfinite(zMm) ||
            referenceCellID == 0) {
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
            snap.cellID = static_cast<std::uint64_t>(snappedCellID);
        } catch (...) {
            snap.x = snap.y = snap.z = nan;
            snap.dx = snap.dy = snap.dz = nan;
            snap.pixX = snap.pixY = -1;
            snap.cellID = 0;
        }
        return snap;
    };
    const auto closestSensorRef = [this](double xMm, double yMm, double zMm) -> const SensorRef* {
        if (!std::isfinite(xMm) || !std::isfinite(yMm) || !std::isfinite(zMm)) {
            return nullptr;
        }

        const dd4hep::Position globalPosition(0.1 * xMm, 0.1 * yMm, 0.1 * zMm);
        const SensorRef* best = nullptr;
        double bestScore = std::numeric_limits<double>::infinity();
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
                if (score < bestScore) {
                    bestScore = score;
                    best = &sensorRef;
                }
            } catch (...) {
                continue;
            }
        }
        return best;
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

    for (const auto* h : simHits) {
        const auto mc = h->getParticle();
        // Accessing an unavailable podio relation is UB; skip instead.
        if (!mc.isAvailable()) continue;

        const auto mom = h->getMomentum();
        const double pmag = std::sqrt(mom.x*mom.x + mom.y*mom.y + mom.z*mom.z);
        const double pT = std::hypot(mom.x, mom.y);

        const uint64_t cid = h->getCellID();
        dd4hep::Position gpos;
        dd4hep::Position lpos;
        try {
            gpos = m_geoSvc->converter()->position(cid);                     // cm, global
            lpos = m_volman.lookupDetElement(cid).nominal()
                       .worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z()));
        } catch (const std::exception&) {
            continue;  // DD4hep can't resolve this cellID; skip the hit instead of failing the event.
        }
        const auto truthPos = h->getPosition();
        const int plane  = m_decoder->get(cid, "layer");
        const int module = m_decoder->get(cid, "module");
        const int sensor = m_decoder->get(cid, "sensor");
        const auto sideIt = m_moduleToSide.find({plane, module});
        const int side = (sideIt != m_moduleToSide.end()) ? sideIt->second : -1;
        const int pixX   = getFieldOr(cid, "x", -1);
        const int pixY   = getFieldOr(cid, "y", -1);
        const int pixZ   = getFieldOr(cid, "z", -1);
        const double path = h->getPathLength();
        const auto id = mc.id();
        const int primaryFlag = isSelectedPrimary(id.collectionID, id.index);
        const auto truthPixel = snapToAclgadPixel(truthPos.x, truthPos.y, truthPos.z, cid);
        const auto readoutPixel =
            snapToAclgadPixel(10. * gpos.x(), 10. * gpos.y(), 10. * gpos.z(), cid);

        vm_xR.push_back(10. * gpos.x());     // cm -> mm
        vm_yR.push_back(10. * gpos.y());
        vm_zR.push_back(10. * gpos.z());
        vm_xT.push_back(truthPos.x);          // mm (EDM4hep convention)
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
        vm_detX.push_back(10. * lpos.x());   // cm -> mm
        vm_detY.push_back(10. * lpos.y());
        vm_detZ.push_back(10. * lpos.z());
        vm_plane .push_back(plane);
        vm_station.push_back(stationFromLayer(plane));
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

        const auto key = std::make_tuple(id.collectionID, id.index, plane, side, module, sensor);
        const auto existing = penetrationIndex.find(key);
        if (existing == penetrationIndex.end()) {
            const std::size_t idx = vm_xP.size();
            penetrationIndex.emplace(key, idx);
            vm_xP     .push_back(truthPos.x);
            vm_yP     .push_back(truthPos.y);
            vm_zP     .push_back(truthPos.z);
            vm_pathP  .push_back(path);
            vm_timeP  .push_back(h->getTime());
            vm_planeP .push_back(plane);
            vm_stationP.push_back(stationFromLayer(plane));
            vm_moduleP.push_back(module);
            vm_sideP  .push_back(side);
            vm_sensorP.push_back(sensor);
            vm_pdgP   .push_back(mc.getPDG());
            vm_statusP.push_back(mc.getGeneratorStatus());
            vm_isPrimaryP.push_back(primaryFlag);
            vm_mcIndexP.push_back(id.index);
            vm_mcCollectionIDP.push_back(id.collectionID);
            vm_pxP    .push_back(mom.x);
            vm_pyP    .push_back(mom.y);
            vm_pzP    .push_back(mom.z);
            vm_pP     .push_back(pmag);
            vm_pTP    .push_back(pT);
        } else if (h->getTime() < vm_timeP[existing->second]) {
            // Earlier step into the same sensor — overwrite position/momentum with entry values.
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

        // Entry/exit upsert: one row per (particle, disk, side, module).
        // Skip hits from modules we couldn't classify (shouldn't happen for
        // B0TrackerHits since Init walked the full tree, but be defensive).
        if (side >= 0) {
            const auto eeKey = std::make_tuple(id.collectionID, id.index, plane, side, module);
            const double thisTime = h->getTime();
            const auto eeIt = entryExitIndex.find(eeKey);
            if (eeIt == entryExitIndex.end()) {
                const std::size_t idx = vm_xEntry.size();
                entryExitIndex.emplace(eeKey, idx);
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
                vm_stationEE.push_back(stationFromLayer(plane));
                vm_moduleEE.push_back(module);
                vm_sideEE.push_back(side);
                vm_pdgEE  .push_back(mc.getPDG());
                vm_isPrimaryEE.push_back(primaryFlag);
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

    // Generator-status convention (HepMC/Pythia8): 1 = stable final-state, 4 = beam.
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

    std::size_t primaryRef = 0;
    for (std::size_t i = 1; i < m_primaryP.size(); ++i) {
        if (m_primaryP[i] > m_primaryP[primaryRef]) primaryRef = i;
    }
    const double primaryP = m_primaryP.empty() ? nan : m_primaryP[primaryRef];
    const double primaryPT = m_primaryPT.empty() ? nan : m_primaryPT[primaryRef];
    if (!m_primaryP.empty()) {
        m_selPrimaryMcIndex = m_primaryMcIndex[primaryRef];
        m_selPrimaryMcCollectionID = m_primaryMcCollectionID[primaryRef];
        m_selPrimaryPx = m_primaryPx[primaryRef];
        m_selPrimaryPy = m_primaryPy[primaryRef];
        m_selPrimaryPz = m_primaryPz[primaryRef];
        m_selPrimaryP = primaryP;
        m_selPrimaryPT = primaryPT;
    }
    if (std::isfinite(primaryP) && !m_genPp.empty()) {
        const double bx = m_genPpx[0];
        const double by = m_genPpy[0];
        const double bz = m_genPpz[0];
        const double bn = std::sqrt(bx * bx + by * by + bz * bz);
        const double pn = primaryP;
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
            if (vm_isPrimaryP[i] == 1 && vm_stationP[i] > 0) {
                stations.insert(vm_stationP[i]);
            }
        }
        m_nStationsPrimary = static_cast<int>(stations.size());
    }

    auto fillChain = [&](TrackChain& out,
                         const auto& trajectories,
                         const auto& tracks,
                         const auto& edmTracks,
                         const auto& assocs,
                         const auto& actsTracks,
                         const auto& actsStates) {
        std::unordered_map<int, std::pair<int, double>> bestAssoc;
        for (const auto* assoc : assocs) {
            if (assoc == nullptr) continue;
            const auto rec = assoc->getRec();
            const auto sim = assoc->getSim();
            if (!rec.isAvailable() || !sim.isAvailable()) continue;
            const int recIndex = rec.id().index;
            const double weight = assoc->getWeight();
            auto it = bestAssoc.find(recIndex);
            if (it == bestAssoc.end() || weight > it->second.second) {
                bestAssoc[recIndex] = {sim.id().index, weight};
            }
        }

        double bestAbsDeltaP = std::numeric_limits<double>::infinity();
        for (std::size_t trajIndex = 0; trajIndex < trajectories.size(); ++trajIndex) {
            const auto* trajectory = trajectories[trajIndex];
            auto tp = (trajectory->trackParameters_size() > 0)
                ? trajectory->getTrackParameters(0)
                : edm4eic::TrackParameters::makeEmpty();
            if (!tp.isAvailable() && trajIndex < tracks.size()) {
                tp = *tracks[trajIndex];
            }
            if (!tp.isAvailable()) continue;

            const float theta  = tp.getTheta();
            const float phi    = tp.getPhi();
            const float qOverP = tp.getQOverP();
            const double p = (qOverP != 0.f) ? std::abs(1.0 / qOverP) : 0.0;
            const double pT = std::abs(p * std::sin(theta));
            const int charge = (qOverP > 0.f) ? 1 : ((qOverP < 0.f) ? -1 : 0);
            const double deltaP = p - primaryP;
            const double deltaPT = pT - primaryPT;
            const int currentTrackIndex = static_cast<int>(trajIndex);
            const int nStates = static_cast<int>(trajectory->getNStates());
            const int nMeasurements = static_cast<int>(trajectory->getNMeasurements());
            const int nOutliers = static_cast<int>(trajectory->getNOutliers());
            const int nHoles = static_cast<int>(trajectory->getNHoles());
            const int nSharedHits = static_cast<int>(trajectory->getNSharedHits());

            double chi2 = nan;
            int ndf = -1;
            if (trajIndex < edmTracks.size() && edmTracks[trajIndex] != nullptr) {
                chi2 = edmTracks[trajIndex]->getChi2();
                ndf = static_cast<int>(edmTracks[trajIndex]->getNdf());
            }
            int assocMc = -1;
            double assocW = 0.0;
            if (const auto it = bestAssoc.find(currentTrackIndex); it != bestAssoc.end()) {
                assocMc = it->second.first;
                assocW = it->second.second;
            }

            out.p.push_back(p);
            out.pT.push_back(pT);
            out.delta_p.push_back(deltaP);
            out.delta_pT.push_back(deltaPT);
            out.theta.push_back(theta);
            out.phi.push_back(phi);
            out.px.push_back(pT * std::cos(phi));
            out.py.push_back(pT * std::sin(phi));
            out.pz.push_back(p * std::cos(theta));
            out.qOverP.push_back(qOverP);

            const auto& loc = tp.getLoc();
            const auto& cov = tp.getCovariance();
            const auto diagSigma = [&cov](unsigned k) -> double {
                const double var = cov(k, k);
                return (std::isfinite(var) && var >= 0.0) ? std::sqrt(var)
                                                          : std::numeric_limits<double>::quiet_NaN();
            };
            out.loc0.push_back(loc.a);
            out.loc1.push_back(loc.b);
            out.sigma_loc0.push_back(diagSigma(0));
            out.sigma_loc1.push_back(diagSigma(1));
            out.sigma_phi.push_back(diagSigma(2));
            out.sigma_theta.push_back(diagSigma(3));
            out.sigma_qOverP.push_back(diagSigma(4));
            out.sigma_time.push_back(diagSigma(5));
            out.chi2.push_back(chi2);
            out.ndf.push_back(ndf);
            out.charge.push_back(charge);
            out.index.push_back(currentTrackIndex);
            out.type.push_back(tp.getType());
            out.surface.push_back(tp.getSurface());
            out.time.push_back(tp.getTime());
            out.pdg.push_back(tp.getPdg());
            out.nStates.push_back(nStates);
            out.nMeasurements.push_back(nMeasurements);
            out.nOutliers.push_back(nOutliers);
            out.nHoles.push_back(nHoles);
            out.nSharedHits.push_back(nSharedHits);
            out.assoc_mcIndex.push_back(assocMc);
            out.assoc_weight.push_back(assocW);

            if (std::isfinite(deltaP) && std::abs(deltaP) < bestAbsDeltaP) {
                bestAbsDeltaP = std::abs(deltaP);
                out.bestIndex = currentTrackIndex;
                out.bestP = p;
                out.bestPT = pT;
                out.bestDeltaP = deltaP;
                out.bestDeltaPT = deltaPT;
                out.bestTheta = theta;
                out.bestPhi = phi;
                out.bestNStates = nStates;
                out.bestNMeasurements = nMeasurements;
                out.bestNOutliers = nOutliers;
                out.bestNHoles = nHoles;
                out.bestNSharedHits = nSharedHits;
                out.bestChi2 = chi2;
                out.bestNdf = ndf;
                out.bestAssocMcIndex = assocMc;
                out.bestAssocWeight = assocW;
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
                const auto surfIt =
                    m_surfaceToSensorIdx.find(state.referenceSurface().geometryId().value());
                if (surfIt != m_surfaceToSensorIdx.end()) {
                    sensorRef = &m_sensorRefs[surfIt->second];
                } else {
                    sensorRef = closestSensorRef(globalX, globalY, globalZ);
                }
                const auto trackPixel = sensorRef
                    ? snapToAclgadPixel(globalX, globalY, globalZ, sensorRef->cellID,
                                        sensorRef->detElement)
                    : PixelSnap{};

                out.state_track_index.push_back(actsTrackIndex);
                out.state_index.push_back(stateIndex++);
                out.state_acts_index.push_back(static_cast<int>(state.index()));
                out.state_type.push_back(typeMask);
                out.state_surface.push_back(state.referenceSurface().geometryId().value());
                out.state_loc0.push_back(loc0);
                out.state_loc1.push_back(loc1);
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
                out.aclgad_plane.push_back(sensorRef ? sensorRef->plane : -1);
                out.aclgad_station.push_back(sensorRef ? stationFromLayer(sensorRef->plane) : -1);
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
        }
    };

    fillChain(m_ts, tsTrajectories, tsTracks, tsEdmTracks, tsAssocs, tsActsTracks, tsActsTrackStates);
    fillChain(m_ckf, ckfTrajectories, ckfTracks, ckfEdmTracks, ckfAssocs, ckfActsTracks, ckfActsTrackStates);

    m_tree->Fill();
}

void B0Trackers::TrackChain::clear(double nan) {
    p.clear(); pT.clear(); delta_p.clear(); delta_pT.clear();
    px.clear(); py.clear(); pz.clear();
    theta.clear(); phi.clear();
    qOverP.clear(); time.clear();
    loc0.clear(); loc1.clear();
    sigma_loc0.clear(); sigma_loc1.clear(); sigma_phi.clear();
    sigma_theta.clear(); sigma_qOverP.clear(); sigma_time.clear();
    chi2.clear(); ndf.clear();
    index.clear(); charge.clear(); type.clear(); pdg.clear(); surface.clear();
    nStates.clear(); nMeasurements.clear(); nOutliers.clear(); nHoles.clear(); nSharedHits.clear();
    assoc_mcIndex.clear(); assoc_weight.clear();
    state_track_index.clear(); state_index.clear(); state_acts_index.clear();
    state_type.clear(); state_pdg.clear();
    state_surface.clear();
    state_loc0.clear(); state_loc1.clear();
    x_on_plane.clear(); y_on_plane.clear(); z_on_plane.clear();
    aclgad_xPix.clear(); aclgad_yPix.clear(); aclgad_zPix.clear();
    aclgad_dx.clear(); aclgad_dy.clear(); aclgad_dz.clear();
    aclgad_pixX.clear(); aclgad_pixY.clear();
    aclgad_plane.clear(); aclgad_module.clear(); aclgad_side.clear();
    aclgad_sensor.clear(); aclgad_station.clear();
    aclgad_cellID.clear();
    state_theta.clear(); state_phi.clear(); state_qOverP.clear(); state_time.clear();

    bestIndex = -1;
    bestNStates = -1;
    bestNMeasurements = -1;
    bestNOutliers = -1;
    bestNHoles = -1;
    bestNSharedHits = -1;
    bestNdf = -1;
    bestAssocMcIndex = -1;
    bestP = nan;
    bestPT = nan;
    bestDeltaP = nan;
    bestDeltaPT = nan;
    bestTheta = nan;
    bestPhi = nan;
    bestChi2 = nan;
    bestAssocWeight = nan;
    hasTrack = 0;
}

void B0Trackers::bindTrackChain(const std::string& trkPrefix, const std::string& bestPrefix, TrackChain& c) {
    auto br = [this](const std::string& name, auto* ptr) {
        m_tree->Branch(name.c_str(), ptr);
    };
    br(trkPrefix + "p", &c.p);
    br(trkPrefix + "pT", &c.pT);
    br(trkPrefix + "delta_p", &c.delta_p);
    br(trkPrefix + "delta_pT", &c.delta_pT);
    br(trkPrefix + "theta", &c.theta);   // lab polar angle, not beam-frame scatter
    br(trkPrefix + "phi", &c.phi);
    br(trkPrefix + "px", &c.px);
    br(trkPrefix + "py", &c.py);
    br(trkPrefix + "pz", &c.pz);
    br(trkPrefix + "qOverP", &c.qOverP);
    br(trkPrefix + "loc0", &c.loc0);
    br(trkPrefix + "loc1", &c.loc1);
    br(trkPrefix + "sigma_loc0", &c.sigma_loc0);
    br(trkPrefix + "sigma_loc1", &c.sigma_loc1);
    br(trkPrefix + "sigma_phi", &c.sigma_phi);
    br(trkPrefix + "sigma_theta", &c.sigma_theta);
    br(trkPrefix + "sigma_qOverP", &c.sigma_qOverP);
    br(trkPrefix + "sigma_time", &c.sigma_time);
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
    br(trkPrefix + "assoc_weight", &c.assoc_weight);
    br(trkPrefix + "state_track_index", &c.state_track_index);
    br(trkPrefix + "state_index", &c.state_index);
    br(trkPrefix + "state_acts_index", &c.state_acts_index);
    br(trkPrefix + "state_type", &c.state_type);
    br(trkPrefix + "state_surface", &c.state_surface);
    br(trkPrefix + "state_loc0", &c.state_loc0);
    br(trkPrefix + "state_loc1", &c.state_loc1);
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
    br(bestPrefix + "index", &c.bestIndex);
    br(bestPrefix + "p", &c.bestP);
    br(bestPrefix + "pT", &c.bestPT);
    br(bestPrefix + "delta_p", &c.bestDeltaP);
    br(bestPrefix + "delta_pT", &c.bestDeltaPT);
    br(bestPrefix + "theta", &c.bestTheta);
    br(bestPrefix + "phi", &c.bestPhi);
    br(bestPrefix + "nStates", &c.bestNStates);
    br(bestPrefix + "nMeasurements", &c.bestNMeasurements);
    br(bestPrefix + "nOutliers", &c.bestNOutliers);
    br(bestPrefix + "nHoles", &c.bestNHoles);
    br(bestPrefix + "nSharedHits", &c.bestNSharedHits);
    br(bestPrefix + "chi2", &c.bestChi2);
    br(bestPrefix + "ndf", &c.bestNdf);
    br(bestPrefix + "assoc_mcIndex", &c.bestAssocMcIndex);
    br(bestPrefix + "assoc_weight", &c.bestAssocWeight);
}

void B0Trackers::Finish() {
    // TTree is owned by the output TFile; RootFile_service writes/closes it.
}

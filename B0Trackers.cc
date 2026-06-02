#include "B0Trackers.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <mutex>
#include <tuple>

#include <JANA/Services/JGlobalRootLock.h>

#include <TFile.h>
#include <TGeoMatrix.h>
#include <TTree.h>

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>

#include <edm4eic/TrackParameters.h>
#include <edm4eic/Trajectory.h>

#include <DD4hep/Objects.h>

#include <services/rootfile/RootFile_service.h>

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
    m_tree->Branch("plane",     &vm_plane);
    m_tree->Branch("module",    &vm_module);
    m_tree->Branch("side",      &vm_side);           // -1=unknown, 0=back, 1=front
    m_tree->Branch("sensor",    &vm_sensor);
    m_tree->Branch("pixX",      &vm_pixX);
    m_tree->Branch("pixY",      &vm_pixY);
    m_tree->Branch("pixZ",      &vm_pixZ);
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
    m_tree->Branch("trk_p",     &trk_p);
    m_tree->Branch("trk_pT",    &trk_pT);
    m_tree->Branch("trk_delta_p", &trk_delta_p);
    m_tree->Branch("trk_delta_pT", &trk_delta_pT);
    m_tree->Branch("trk_theta", &trk_theta);
    m_tree->Branch("trk_phi",   &trk_phi);
    m_tree->Branch("trk_px",    &trk_px);
    m_tree->Branch("trk_py",    &trk_py);
    m_tree->Branch("trk_pz",    &trk_pz);
    m_tree->Branch("trk_qOverP",&trk_qOverP);
    m_tree->Branch("trk_charge",&trk_charge);
    m_tree->Branch("trk_index", &trk_index);
    m_tree->Branch("trk_type",  &trk_type);
    m_tree->Branch("trk_surface", &trk_surface);
    m_tree->Branch("trk_time",  &trk_time);
    m_tree->Branch("trk_pdg",   &trk_pdg);
    m_tree->Branch("trk_nStates", &trk_nStates);
    m_tree->Branch("trk_nMeasurements", &trk_nMeasurements);
    m_tree->Branch("trk_nOutliers", &trk_nOutliers);
    m_tree->Branch("trk_nHoles", &trk_nHoles);
    m_tree->Branch("trk_nSharedHits", &trk_nSharedHits);
    m_tree->Branch("best_trk_index", &m_bestTrkIndex);
    m_tree->Branch("best_trk_p", &m_bestTrkP);
    m_tree->Branch("best_trk_pT", &m_bestTrkPT);
    m_tree->Branch("best_trk_delta_p", &m_bestTrkDeltaP);
    m_tree->Branch("best_trk_delta_pT", &m_bestTrkDeltaPT);
    m_tree->Branch("best_trk_theta", &m_bestTrkTheta);
    m_tree->Branch("best_trk_phi", &m_bestTrkPhi);
    m_tree->Branch("best_trk_nStates", &m_bestTrkNStates);
    m_tree->Branch("best_trk_nMeasurements", &m_bestTrkNMeasurements);
    m_tree->Branch("best_trk_nOutliers", &m_bestTrkNOutliers);
    m_tree->Branch("best_trk_nHoles", &m_bestTrkNHoles);
    m_tree->Branch("best_trk_nSharedHits", &m_bestTrkNSharedHits);

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
    auto ro  = m_geoSvc->detector()->readout("B0TrackerHits");
    m_decoder = ro.idSpec().decoder();
    m_volman  = m_geoSvc->detector()->volumeManager();

    // Build (layer, module) -> side map by walking the B0Tracker DetElement tree.
    // Each module is placed in the layer Assembly at z = +ModuleOffsetFromSupport
    // (front) or -offset (back); the placement's translation z carries the sign.
    // Key uses the cellID's physVolID values (not the DetElement id, which has
    // different semantics in dev vs official geometries).
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
            for (const auto& [modName, modDE] : layerDE.children()) {
                const auto pv = modDE.placement();
                if (!pv.isValid()) continue;
                const int module_id = findVolID(pv, "module");
                if (module_id < 0) continue;
                const TGeoMatrix& mat = pv.matrix();
                const double* tr = mat.GetTranslation();
                m_moduleToSide[{layer_id, module_id}] = (tr[2] > 0.0) ? 1 : 0;
            }
        }
    }

    rootLock->release_lock();
}

void B0Trackers::Process(const std::shared_ptr<const JEvent>& event) {
    // Fetch before locking: event->Get() may run expensive factories (CKF),
    // and JANA already serializes that work internally.
    auto mcparticles = event->Get<edm4hep::MCParticle>("MCParticles");
    auto simHits     = event->Get<edm4hep::SimTrackerHit>("B0TrackerHits");
    auto tracks      = event->Get<edm4eic::TrackParameters>("B0TrackerCKFTruthSeededTrackParameters");
    auto trajectories = event->Get<edm4eic::Trajectory>("B0TrackerCKFTruthSeededTrajectories");

    std::lock_guard<std::mutex> lock(m_fillMutex);

    m_eventNumber = event->GetEventNumber();

    vm_xR.clear();    vm_yR.clear();    vm_zR.clear();
    vm_xT.clear();    vm_yT.clear();    vm_zT.clear();
    vm_detX.clear();  vm_detY.clear();  vm_detZ.clear();
    vm_plane.clear(); vm_module.clear(); vm_side.clear();   vm_sensor.clear();
    vm_pixX.clear();  vm_pixY.clear();   vm_pixZ.clear();
    vm_cellID.clear(); vm_mcIndex.clear(); vm_mcCollectionID.clear();
    vm_eDep.clear();  vm_time.clear();   vm_path.clear();
    vm_pdg.clear();   vm_status.clear();
    vm_isPrimary.clear();
    vm_px.clear();    vm_py.clear();    vm_pz.clear();   vm_p.clear();    vm_pT.clear();

    vm_xP.clear();      vm_yP.clear();      vm_zP.clear();      vm_pathP.clear();   vm_timeP.clear();
    vm_planeP.clear();  vm_moduleP.clear(); vm_sideP.clear(); vm_sensorP.clear();

    vm_xEntry.clear();   vm_yEntry.clear();  vm_zEntry.clear();  vm_timeEntry.clear();
    vm_pxEntry.clear();  vm_pyEntry.clear(); vm_pzEntry.clear(); vm_pEntry.clear(); vm_pTEntry.clear();
    vm_xExit.clear();    vm_yExit.clear();   vm_zExit.clear();   vm_timeExit.clear();
    vm_pxExit.clear();   vm_pyExit.clear();  vm_pzExit.clear();  vm_pExit.clear();  vm_pTExit.clear();
    vm_sensorEntry.clear(); vm_sensorExit.clear();
    vm_cellIDEntry.clear(); vm_cellIDExit.clear();
    vm_planeEE.clear();  vm_moduleEE.clear(); vm_sideEE.clear();  vm_pdgEE.clear();
    vm_isPrimaryEE.clear();
    vm_mcIndexEE.clear(); vm_mcCollectionIDEE.clear(); vm_nStepsEE.clear();
    vm_statusEE.clear();
    vm_pdgP.clear();    vm_statusP.clear(); vm_isPrimaryP.clear(); vm_mcIndexP.clear(); vm_mcCollectionIDP.clear();
    vm_pxP.clear();     vm_pyP.clear();     vm_pzP.clear();     vm_pP.clear();     vm_pTP.clear();

    trk_p.clear();    trk_pT.clear();   trk_delta_p.clear(); trk_delta_pT.clear();
    trk_px.clear();   trk_py.clear();   trk_pz.clear();
    trk_theta.clear();trk_phi.clear();
    trk_qOverP.clear(); trk_time.clear();
    trk_index.clear(); trk_charge.clear(); trk_type.clear(); trk_pdg.clear(); trk_surface.clear();
    trk_nStates.clear(); trk_nMeasurements.clear(); trk_nOutliers.clear(); trk_nHoles.clear(); trk_nSharedHits.clear();

    const double nan = std::numeric_limits<double>::quiet_NaN();
    m_bestTrkIndex = -1;
    m_bestTrkNStates = -1;
    m_bestTrkNMeasurements = -1;
    m_bestTrkNOutliers = -1;
    m_bestTrkNHoles = -1;
    m_bestTrkNSharedHits = -1;
    m_bestTrkP = nan;
    m_bestTrkPT = nan;
    m_bestTrkDeltaP = nan;
    m_bestTrkDeltaPT = nan;
    m_bestTrkTheta = nan;
    m_bestTrkPhi = nan;

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
    // entered and exited the silicon in each module.
    std::map<std::tuple<uint32_t, int, int, int, int>, std::size_t> entryExitIndex;
    const auto getFieldOr = [this](std::uint64_t cellID, const char* field, int fallback) -> int {
        try {
            return static_cast<int>(m_decoder->get(cellID, field));
        } catch (...) {
            return fallback;
        }
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

        vm_xR.push_back(10. * gpos.x());     // cm -> mm
        vm_yR.push_back(10. * gpos.y());
        vm_zR.push_back(10. * gpos.z());
        vm_xT.push_back(truthPos.x);          // mm (EDM4hep convention)
        vm_yT.push_back(truthPos.y);
        vm_zT.push_back(truthPos.z);
        vm_detX.push_back(10. * lpos.x());   // cm -> mm
        vm_detY.push_back(10. * lpos.y());
        vm_detZ.push_back(10. * lpos.z());
        vm_plane .push_back(plane);
        vm_module.push_back(module);
        vm_side  .push_back(side);
        vm_sensor.push_back(sensor);
        vm_pixX  .push_back(pixX);
        vm_pixY  .push_back(pixY);
        vm_pixZ  .push_back(pixZ);
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
                vm_planeEE.push_back(plane);        vm_moduleEE.push_back(module);
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

    int trackIndex = 0;
    std::size_t primaryRef = 0;
    for (std::size_t i = 1; i < m_primaryP.size(); ++i) {
        if (m_primaryP[i] > m_primaryP[primaryRef]) primaryRef = i;
    }
    const double primaryP = m_primaryP.empty() ? nan : m_primaryP[primaryRef];
    const double primaryPT = m_primaryPT.empty() ? nan : m_primaryPT[primaryRef];
    double bestAbsDeltaP = std::numeric_limits<double>::infinity();
    for (const auto* tp : tracks) {
        const float theta  = tp->getTheta();
        const float phi    = tp->getPhi();
        const float qOverP = tp->getQOverP();
        const double p = (qOverP != 0.f) ? std::abs(1.0 / qOverP) : 0.0;
        const double pT = std::abs(p * std::sin(theta));
        const int charge = (qOverP > 0.f) ? 1 : ((qOverP < 0.f) ? -1 : 0);

        const double deltaP = p - primaryP;
        const double deltaPT = pT - primaryPT;
        const int currentTrackIndex = trackIndex++;
        const bool hasTrajectory = currentTrackIndex >= 0 &&
                                   static_cast<std::size_t>(currentTrackIndex) < trajectories.size();
        const int nStates = hasTrajectory ? static_cast<int>(trajectories[currentTrackIndex]->getNStates()) : -1;
        const int nMeasurements = hasTrajectory ? static_cast<int>(trajectories[currentTrackIndex]->getNMeasurements()) : -1;
        const int nOutliers = hasTrajectory ? static_cast<int>(trajectories[currentTrackIndex]->getNOutliers()) : -1;
        const int nHoles = hasTrajectory ? static_cast<int>(trajectories[currentTrackIndex]->getNHoles()) : -1;
        const int nSharedHits = hasTrajectory ? static_cast<int>(trajectories[currentTrackIndex]->getNSharedHits()) : -1;

        trk_p    .push_back(p);
        trk_pT   .push_back(pT);
        trk_delta_p .push_back(deltaP);
        trk_delta_pT.push_back(deltaPT);
        trk_theta.push_back(theta);
        trk_phi  .push_back(phi);
        trk_px   .push_back(pT * std::cos(phi));
        trk_py   .push_back(pT * std::sin(phi));
        trk_pz   .push_back(p * std::cos(theta));
        trk_qOverP.push_back(qOverP);
        trk_charge.push_back(charge);
        trk_index .push_back(currentTrackIndex);
        trk_type  .push_back(tp->getType());
        trk_surface.push_back(tp->getSurface());
        trk_time  .push_back(tp->getTime());
        trk_pdg   .push_back(tp->getPdg());
        trk_nStates.push_back(nStates);
        trk_nMeasurements.push_back(nMeasurements);
        trk_nOutliers.push_back(nOutliers);
        trk_nHoles.push_back(nHoles);
        trk_nSharedHits.push_back(nSharedHits);

        if (std::isfinite(deltaP) && std::abs(deltaP) < bestAbsDeltaP) {
            bestAbsDeltaP = std::abs(deltaP);
            m_bestTrkIndex = currentTrackIndex;
            m_bestTrkP = p;
            m_bestTrkPT = pT;
            m_bestTrkDeltaP = deltaP;
            m_bestTrkDeltaPT = deltaPT;
            m_bestTrkTheta = theta;
            m_bestTrkPhi = phi;
            m_bestTrkNStates = nStates;
            m_bestTrkNMeasurements = nMeasurements;
            m_bestTrkNOutliers = nOutliers;
            m_bestTrkNHoles = nHoles;
            m_bestTrkNSharedHits = nSharedHits;
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

    m_tree->Fill();
}

void B0Trackers::Finish() {
    // TTree is owned by the output TFile; RootFile_service writes/closes it.
}

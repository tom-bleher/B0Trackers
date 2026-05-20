#include "B0Trackers.h"

#include <cmath>
#include <cstdint>
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
    m_tree->Branch("status",    &vm_status);
    m_tree->Branch("xP",        &vm_xP);
    m_tree->Branch("yP",        &vm_yP);
    m_tree->Branch("zP",        &vm_zP);
    m_tree->Branch("pathP",     &vm_pathP);
    m_tree->Branch("timeP",     &vm_timeP);
    m_tree->Branch("planeP",    &vm_planeP);
    m_tree->Branch("moduleP",   &vm_moduleP);
    m_tree->Branch("sensorP",   &vm_sensorP);
    m_tree->Branch("pdgP",      &vm_pdgP);
    m_tree->Branch("statusP",   &vm_statusP);
    m_tree->Branch("mcIndexP",  &vm_mcIndexP);
    m_tree->Branch("mcCollectionIDP", &vm_mcCollectionIDP);
    m_tree->Branch("pxP",       &vm_pxP);
    m_tree->Branch("pyP",       &vm_pyP);
    m_tree->Branch("pzP",       &vm_pzP);
    m_tree->Branch("pP",        &vm_pP);
    m_tree->Branch("beampx",    &beam_px);
    m_tree->Branch("beampy",    &beam_py);
    m_tree->Branch("beampz",    &beam_pz);
    m_tree->Branch("beamp",     &beam_p);
    m_tree->Branch("beam_pdg",  &beam_pdg);
    m_tree->Branch("genPpx",    &m_genPpx);
    m_tree->Branch("genPpy",    &m_genPpy);
    m_tree->Branch("genPpz",    &m_genPpz);
    m_tree->Branch("genPp",     &m_genPp);
    m_tree->Branch("genBeamP",  &m_genBeamP);
    m_tree->Branch("genBeamPx", &m_genBeamPx);
    m_tree->Branch("genBeamPy", &m_genBeamPy);
    m_tree->Branch("genBeamPz", &m_genBeamPz);
    m_tree->Branch("genBeamPP", &m_genBeamPP);
    m_tree->Branch("genBeamPPx",&m_genBeamPPx);
    m_tree->Branch("genBeamPPy",&m_genBeamPPy);
    m_tree->Branch("genBeamPPz",&m_genBeamPPz);
    m_tree->Branch("trk_p",     &trk_p);
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

    // Entry/exit summary branches (one row per mc-particle/disk/side).
    m_tree->Branch("xEntry",     &vm_xEntry);
    m_tree->Branch("yEntry",     &vm_yEntry);
    m_tree->Branch("zEntry",     &vm_zEntry);
    m_tree->Branch("timeEntry",  &vm_timeEntry);
    m_tree->Branch("pxEntry",    &vm_pxEntry);
    m_tree->Branch("pyEntry",    &vm_pyEntry);
    m_tree->Branch("pzEntry",    &vm_pzEntry);
    m_tree->Branch("pEntry",     &vm_pEntry);
    m_tree->Branch("xExit",      &vm_xExit);
    m_tree->Branch("yExit",      &vm_yExit);
    m_tree->Branch("zExit",      &vm_zExit);
    m_tree->Branch("timeExit",   &vm_timeExit);
    m_tree->Branch("pxExit",     &vm_pxExit);
    m_tree->Branch("pyExit",     &vm_pyExit);
    m_tree->Branch("pzExit",     &vm_pzExit);
    m_tree->Branch("pExit",      &vm_pExit);
    m_tree->Branch("planeEE",    &vm_planeEE);
    m_tree->Branch("sideEE",     &vm_sideEE);          // 0=back, 1=front
    m_tree->Branch("pdgEE",      &vm_pdgEE);
    m_tree->Branch("mcIndexEE",  &vm_mcIndexEE);
    m_tree->Branch("mcCollectionIDEE", &vm_mcCollectionIDEE);
    m_tree->Branch("nStepsEE",   &vm_nStepsEE);
    m_tree->Branch("statusEE",   &vm_statusEE);   // 1=primary final-state, 0=secondary, 4=beam

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

    std::lock_guard<std::mutex> lock(m_fillMutex);

    m_eventNumber = event->GetEventNumber();

    vm_xR.clear();    vm_yR.clear();    vm_zR.clear();
    vm_xT.clear();    vm_yT.clear();    vm_zT.clear();
    vm_detX.clear();  vm_detY.clear();  vm_detZ.clear();
    vm_plane.clear(); vm_module.clear(); vm_sensor.clear();
    vm_pixX.clear();  vm_pixY.clear();   vm_pixZ.clear();
    vm_cellID.clear(); vm_mcIndex.clear(); vm_mcCollectionID.clear();
    vm_eDep.clear();  vm_time.clear();   vm_path.clear();
    vm_pdg.clear();   vm_status.clear();
    vm_px.clear();    vm_py.clear();    vm_pz.clear();   vm_p.clear();

    vm_xP.clear();      vm_yP.clear();      vm_zP.clear();      vm_pathP.clear();   vm_timeP.clear();
    vm_planeP.clear();  vm_moduleP.clear(); vm_sensorP.clear();

    vm_xEntry.clear();   vm_yEntry.clear();  vm_zEntry.clear();  vm_timeEntry.clear();
    vm_pxEntry.clear();  vm_pyEntry.clear(); vm_pzEntry.clear(); vm_pEntry.clear();
    vm_xExit.clear();    vm_yExit.clear();   vm_zExit.clear();   vm_timeExit.clear();
    vm_pxExit.clear();   vm_pyExit.clear();  vm_pzExit.clear();  vm_pExit.clear();
    vm_planeEE.clear();  vm_sideEE.clear();  vm_pdgEE.clear();
    vm_mcIndexEE.clear(); vm_mcCollectionIDEE.clear(); vm_nStepsEE.clear();
    vm_statusEE.clear();
    vm_pdgP.clear();    vm_statusP.clear(); vm_mcIndexP.clear(); vm_mcCollectionIDP.clear();
    vm_pxP.clear();     vm_pyP.clear();     vm_pzP.clear();     vm_pP.clear();

    trk_p.clear();    trk_px.clear();   trk_py.clear();  trk_pz.clear();
    trk_theta.clear();trk_phi.clear();
    trk_qOverP.clear(); trk_time.clear();
    trk_index.clear(); trk_charge.clear(); trk_type.clear(); trk_pdg.clear(); trk_surface.clear();

    beam_px.clear();  beam_py.clear();  beam_pz.clear(); beam_p.clear();
    beam_pdg.clear();
    m_genPpx.clear();    m_genPpy.clear();    m_genPpz.clear();    m_genPp.clear();
    m_genBeamPx.clear(); m_genBeamPy.clear(); m_genBeamPz.clear(); m_genBeamP.clear();
    m_genBeamPPx.clear();m_genBeamPPy.clear();m_genBeamPPz.clear();m_genBeamPP.clear();

    // Key: (mcCollectionID, mcIndex, layer, module, sensor) — one *P entry per
    // (particle, sensitive-volume) pair. `module` is globally unique so front-
    // vs-back placements on the same disk land in distinct buckets.
    std::map<std::tuple<uint32_t, int, int, int, int>, std::size_t> penetrationIndex;

    // Key: (mcCollectionID, mcIndex, layer, side) — one entry/exit row per
    // (particle, disk side). Tracks the smallest- and largest-time SimTrackerHits
    // contributing to that group so we can report where the particle entered
    // and exited the silicon on each side of each disk.
    std::map<std::tuple<uint32_t, int, int, int>, std::size_t> entryExitIndex;
    const auto getFieldOr = [this](std::uint64_t cellID, const char* field, int fallback) -> int {
        try {
            return static_cast<int>(m_decoder->get(cellID, field));
        } catch (...) {
            return fallback;
        }
    };

    for (const auto* h : simHits) {
        const auto mc = h->getParticle();
        // Accessing an unavailable podio relation is UB; skip instead.
        if (!mc.isAvailable()) continue;

        const auto mom = mc.getMomentum();
        const double pmag = std::sqrt(mom.x*mom.x + mom.y*mom.y + mom.z*mom.z);

        const uint64_t cid = h->getCellID();
        const auto gpos = m_geoSvc->converter()->position(cid);              // cm, global
        const auto lpos = m_volman.lookupDetElement(cid).nominal()
                              .worldToLocal(dd4hep::Position(gpos.x(), gpos.y(), gpos.z())); // cm, local
        const auto truthPos = h->getPosition();
        const int plane  = m_decoder->get(cid, "layer");
        const int module = m_decoder->get(cid, "module");
        const int sensor = m_decoder->get(cid, "sensor");
        const int pixX   = getFieldOr(cid, "x", -1);
        const int pixY   = getFieldOr(cid, "y", -1);
        const int pixZ   = getFieldOr(cid, "z", -1);
        const double path = h->getPathLength();
        const auto id = mc.id();

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
        vm_status.push_back(mc.getGeneratorStatus());

        const auto key = std::make_tuple(id.collectionID, id.index, plane, module, sensor);
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
            vm_sensorP.push_back(sensor);
            vm_pdgP   .push_back(mc.getPDG());
            vm_statusP.push_back(mc.getGeneratorStatus());
            vm_mcIndexP.push_back(id.index);
            vm_mcCollectionIDP.push_back(id.collectionID);
            vm_pxP    .push_back(mom.x);
            vm_pyP    .push_back(mom.y);
            vm_pzP    .push_back(mom.z);
            vm_pP     .push_back(pmag);
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
        }

        // Entry/exit upsert: one row per (particle, disk, side).
        // Skip hits from modules we couldn't classify (shouldn't happen for
        // B0TrackerHits since Init walked the full tree, but be defensive).
        const auto sideIt = m_moduleToSide.find({plane, module});
        if (sideIt != m_moduleToSide.end()) {
            const int side = sideIt->second;
            const auto eeKey = std::make_tuple(id.collectionID, id.index, plane, side);
            const double thisTime = h->getTime();
            const auto eeIt = entryExitIndex.find(eeKey);
            if (eeIt == entryExitIndex.end()) {
                const std::size_t idx = vm_xEntry.size();
                entryExitIndex.emplace(eeKey, idx);
                vm_xEntry .push_back(truthPos.x);   vm_yEntry .push_back(truthPos.y);
                vm_zEntry .push_back(truthPos.z);   vm_timeEntry.push_back(thisTime);
                vm_pxEntry.push_back(mom.x);        vm_pyEntry.push_back(mom.y);
                vm_pzEntry.push_back(mom.z);        vm_pEntry .push_back(pmag);
                vm_xExit  .push_back(truthPos.x);   vm_yExit  .push_back(truthPos.y);
                vm_zExit  .push_back(truthPos.z);   vm_timeExit.push_back(thisTime);
                vm_pxExit .push_back(mom.x);        vm_pyExit .push_back(mom.y);
                vm_pzExit .push_back(mom.z);        vm_pExit  .push_back(pmag);
                vm_planeEE.push_back(plane);        vm_sideEE.push_back(side);
                vm_pdgEE  .push_back(mc.getPDG());
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
                }
            }
        }
    }

    int trackIndex = 0;
    for (const auto* tp : tracks) {
        const float theta  = tp->getTheta();
        const float phi    = tp->getPhi();
        const float qOverP = tp->getQOverP();
        const double p = (qOverP != 0.f) ? std::abs(1.0 / qOverP) : 0.0;
        const int charge = (qOverP > 0.f) ? 1 : ((qOverP < 0.f) ? -1 : 0);

        trk_p    .push_back(p);
        trk_theta.push_back(theta);
        trk_phi  .push_back(phi);
        trk_px   .push_back(p * std::sin(theta) * std::cos(phi));
        trk_py   .push_back(p * std::sin(theta) * std::sin(phi));
        trk_pz   .push_back(p * std::cos(theta));
        trk_qOverP.push_back(qOverP);
        trk_charge.push_back(charge);
        trk_index .push_back(trackIndex++);
        trk_type  .push_back(tp->getType());
        trk_surface.push_back(tp->getSurface());
        trk_time  .push_back(tp->getTime());
        trk_pdg   .push_back(tp->getPdg());
    }

    // Generator-status convention (HepMC/Pythia8): 1 = stable final-state, 4 = beam.
    for (const auto* part : mcparticles) {
        const auto p = part->getMomentum();
        const int  pdg    = part->getPDG();
        const int  status = part->getGeneratorStatus();
        const double pmag = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);

        if (status == 4 &&
            (pdg == 22 || pdg == 11 || pdg == -11 || pdg == 2212)) {
            beam_px .push_back(p.x);
            beam_py .push_back(p.y);
            beam_pz .push_back(p.z);
            beam_p  .push_back(pmag);
            beam_pdg.push_back(pdg);
        }
        if (pdg == 2212 && status == 4) {
            m_genPpx.push_back(p.x);
            m_genPpy.push_back(p.y);
            m_genPpz.push_back(p.z);
            m_genPp .push_back(pmag);
        }
        if (pdg == 2212 && status == 1) {
            m_genBeamPPx.push_back(p.x);
            m_genBeamPPy.push_back(p.y);
            m_genBeamPPz.push_back(p.z);
            m_genBeamPP .push_back(pmag);
        }
        if (pdg == 11 && status == 1) {
            m_genBeamPx.push_back(p.x);
            m_genBeamPy.push_back(p.y);
            m_genBeamPz.push_back(p.z);
            m_genBeamP .push_back(pmag);
        }
    }

    m_tree->Fill();
}

void B0Trackers::Finish() {
    // TTree is owned by the output TFile; RootFile_service writes/closes it.
}

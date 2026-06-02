#pragma once

#include <memory>
#include <mutex>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include <JANA/JEventProcessor.h>

#include <DD4hep/VolumeManager.h>
#include <DDSegmentation/BitFieldCoder.h>

#include <services/geometry/dd4hep/DD4hep_service.h>

class TTree;

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
    const dd4hep::DDSegmentation::BitFieldCoder*  m_decoder = nullptr;
    dd4hep::VolumeManager                         m_volman;

    // (layer, module-volID) -> side: 1=front (placed at +ModuleOffsetFromSupport
    // in layer frame), 0=back (placed at -offset). Built once in Init by walking
    // the B0Tracker DetElement tree. The key uses the physVolIDs that the cellID
    // decoder reports so it works across dev / official B0 geometries (their
    // DetElement::id() conventions differ).
    std::map<std::pair<int, int>, int> m_moduleToSide;

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
    std::vector<double> vm_detX, vm_detY, vm_detZ;  // local sensor frame, mm
    std::vector<int>    vm_plane, vm_module, vm_side, vm_sensor, vm_pixX, vm_pixY, vm_pixZ;
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
    std::vector<int>    vm_planeP, vm_moduleP, vm_sideP, vm_sensorP;
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
    std::vector<int>    vm_planeEE, vm_moduleEE, vm_sideEE, vm_pdgEE, vm_isPrimaryEE, vm_mcIndexEE, vm_nStepsEE;
    std::vector<int>    vm_statusEE;   // linked MCParticle::generatorStatus
    std::vector<std::uint32_t> vm_mcCollectionIDEE;

    // Per-event ACTS track vectors (B0 CKF, truth-seeded).
    std::vector<double> trk_p, trk_pT, trk_px, trk_py, trk_pz, trk_theta, trk_phi;
    std::vector<double> trk_qOverP, trk_time;
    std::vector<int>    trk_index, trk_charge, trk_type, trk_pdg;
    std::vector<std::uint64_t> trk_surface;

    // Per-event MC truth.
    std::vector<double> beam_px, beam_py, beam_pz, beam_p, beam_pT;
    std::vector<int>    beam_pdg;
    std::vector<double> m_primaryPx, m_primaryPy, m_primaryPz, m_primaryP, m_primaryPT;
    std::vector<int>    m_primaryPdgOut, m_primaryStatusOut, m_primaryMcIndex;
    std::vector<std::uint32_t> m_primaryMcCollectionID;
    std::vector<double> m_genPpx, m_genPpy, m_genPpz, m_genPp, m_genPpT;          // beam proton (status 4)
    std::vector<double> m_genBeamPx, m_genBeamPy, m_genBeamPz, m_genBeamP, m_genBeamPT;     // scattered electron (status 1)
    std::vector<double> m_genBeamPPx, m_genBeamPPy, m_genBeamPPz, m_genBeamPP, m_genBeamPPT; // scattered proton (status 1)
};

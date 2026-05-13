#include <algorithm>
#include <bitset>
#include <spdlog/spdlog.h>
#include <cstdlib>
#include <iostream>

//#include <JANA/JEventProcessorSequentialRoot.h>
#include <JANA/JEventProcessor.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <THashList.h>

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4hep/SimCalorimeterHit.h>
#include <edm4hep/RawCalorimeterHit.h>

#include <edm4eic/Cluster.h>
//#include <edm4eic/RawCalorimeterHit.h>
#include <edm4eic/ProtoCluster.h>

#include <services/geometry/dd4hep/DD4hep_service.h>

#include "variables.h"
#include <edm4eic/TrackParameters.h>
#include <edm4eic/TrackSeedCollection.h>
#include <edm4eic/TrackParametersCollection.h>
#include <edm4eic/Trajectory.h>

using namespace std;

//class B0Trackers : public JEventProcessorSequentialRoot {
   
class B0Trackers : public JEventProcessor {

private:
    // Data objects we will need from JANA e.g.
/*    PrefetchT<edm4hep::MCParticle> MCParticles          = {this, "MCParticles"};

    
    PrefetchT<edm4hep::SimTrackerHit> Tracker_hits      = {this, "B0TrackerHits"};


// NEW: ACTS output (names match tracking plugin defaults)
//PrefetchT<edm4eic::TrackParameters> TrackParameters = {this, "TrackParameters"};
//PrefetchT<edm4eic::Trajectory> CentralCKFTrajectories = {this, "CentralCKFTrajectories"};
 PrefetchT<edm4eic::TrackParameters> B0CKFTrackParams = {this, "B0TrackerCKFTrackParameters"};
*/
/*
auto MCParticles = event->Get<edm4hep::MCParticle>("MCParticles");
auto Tracker_hits    = event->Get<edm4hep::SimTrackerHit>("B0TrackerHits");
auto B0CKFTrackParams  = event->Get<edm4eic::TrackParameters>("B0TrackerCKFTrackParameters");
*/

    // Geometry decoder & segmentation
    const dd4hep::DDSegmentation::BitFieldCoder* m_decoder = nullptr;

    const dd4hep::DDSegmentation::Segmentation* m_seg = nullptr;


    std::shared_ptr<DD4hep_service> m_geoSvc = nullptr;
    dd4hep::VolumeManager m_volman;



    std::vector<double> trk_p;
    std::vector<double> trk_px;
    std::vector<double> trk_py;
    std::vector<double> trk_pz;
std::vector<double> trk_theta;
std::vector<double> trk_phi;

public:
//    B0Trackers()             { SetTypeName(NAME_OF_THIS); }

//    void InitWithGlobalRootLock() override;
//    void ProcessSequential(const std::shared_ptr<const JEvent>& event) override;
    
//    void FinishWithGlobalRootLock() override;

    void Init() override;
    void Process(const std::shared_ptr<const JEvent>& event) override;
    void Finish() override;




    void MCgenAnalysis(const std::vector<const edm4hep::MCParticle*>& mcparts);
  
    
    // ---- ROOT objects ----
    TTree* m_tree = nullptr;

    double m_xR, m_yR, m_zR;
    double m_xT, m_yT, m_zT;
    double m_detX, m_detY, m_detZ;
    double m_eDep;
    double m_time;
    double m_HitPath;

    double m_Truepx;
    double m_Truepy;
    double m_Truepz;
    double m_Truep;
    int m_TruePDG;
    
    int m_plane, m_pixX, m_pixZ;
    int m_module, m_sensor;
    int m_primary;

	vector<TrackClass> m_TrackerHits;
	vector<HitClass> m_HitsClas;
	VertexDataTrue m_VertexTrue;


	vector<double>	  m_xi;
	vector<double>    m_yi;
	vector<double>    m_zi;
	vector<double>    m_pxi;
	vector<double>    m_pyi;
	vector<double>    m_pzi;
	vector<double>    m_pathL;
	
	std::vector<double> vm_xR;
	std::vector<double> vm_yR;
	std::vector<double> vm_zR;
	std::vector<double> vm_xT;
	std::vector<double> vm_yT;
	std::vector<double> vm_zT;
	std::vector<int> vm_plane;
	std::vector<int> vm_module;
	std::vector<int> vm_sensor;
	std::vector<double> vm_detX;
	std::vector<double> vm_detY;
	std::vector<double> vm_detZ;
	std::vector<int> vm_pdg;
	std::vector<double> vm_px;
	std::vector<double> vm_py;
	std::vector<double> vm_pz;
	std::vector<double> vm_p;
	std::vector<int> vm_status;
	std::vector<double> m_genPpx;
	std::vector<double> m_genPpy;
	std::vector<double> m_genPpz;
	std::vector<double> m_genPp;

	std::vector<double> m_genBeamP;
        std::vector<double> m_genBeamPx;
        std::vector<double> m_genBeamPy;
        std::vector<double> m_genBeamPz;

        std::vector<double> m_genBeamPP;
        std::vector<double> m_genBeamPPx;
        std::vector<double> m_genBeamPPy;
        std::vector<double> m_genBeamPPz;

	std::vector<double> beam_px;
	std::vector<double> beam_py;
	std::vector<double> beam_pz;
	std::vector<double> beam_p;
	std::vector<int> beam_pdg;


};

#include <services/rootfile/RootFile_service.h>
#include "B0Trackers.h"
//#include "TrackingUtils.h"
#include <cmath>
#include <vector>
#include "variables.h"  // for HitClass and Track


// The following just makes this a JANA plugin
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app->Add(new B0Trackers);
  }
}

//-------------------------------------------
// InitWithGlobalRootLock
//-------------------------------------------
//void B0Trackers::InitWithGlobalRootLock() {
void B0Trackers::Init(){   
	
    auto rf_svc  = GetApplication()->GetService<RootFile_service>();
//    auto outfile = rf_svc->GetHistFile();
//    auto outfile = rf_svc->GetFile();
//    outfile->mkdir("B0Trackers")->cd();

    TDirectory* histdir =  rf_svc->GetHistFile();
    TFile* outfile =histdir->GetFile();

    outfile->mkdir("B0Trackers")->cd();

    m_tree = new TTree("hits", "Truth vs read-out");
/*    m_tree->Branch("xT", &m_xT); 
    m_tree->Branch("yT", &m_yT); 
    m_tree->Branch("zT", &m_zT);
    m_tree->Branch("xR", &m_xR); 
    m_tree->Branch("yR", &m_yR); 
    m_tree->Branch("zR", &m_zR);
    m_tree->Branch("detX", &m_detX);
    m_tree->Branch("detY", &m_detY);
    m_tree->Branch("detZ", &m_detZ);
    m_tree->Branch("plane", &m_plane);
    m_tree->Branch("module", &m_module);
    m_tree->Branch("sensor", &m_sensor);
    m_tree->Branch("primary",&m_primary);
    m_tree->Branch("pixX", &m_pixX);
    m_tree->Branch("pixY", &m_pixY);
    m_tree->Branch("pixZ", &m_pixZ);
    m_tree->Branch("eDep",&m_eDep);
    m_tree->Branch("time",&m_time);
    m_tree->Branch("hitPath",&m_HitPath);
    m_tree->Branch("truePx",&m_Truepx);
    m_tree->Branch("truePy",&m_Truepy);
    m_tree->Branch("truePz",&m_Truepz);
    m_tree->Branch("trueP",&m_Truep);
	m_tree->Branch("truePDG",&m_TruePDG);
	m_tree->Branch("truePDG",&m_VertexTrue.pdg);
	m_tree->Branch("trueE",&m_VertexTrue.E);
	m_tree->Branch("trueTheta",&m_VertexTrue.theta);
	m_tree->Branch("truePhi",&m_VertexTrue.phi);
	m_tree->Branch("trueVx",&m_VertexTrue.vx);
	m_tree->Branch("trueVy",&m_VertexTrue.vy);
	m_tree->Branch("trueVz",&m_VertexTrue.vz);
	m_tree->Branch("trueStatus",&m_VertexTrue.status);
*/

    m_tree->Branch("xR",&vm_xR);
    m_tree->Branch("yR",&vm_yR);
    m_tree->Branch("zR",&vm_zR);
    m_tree->Branch("xT",&vm_xT);
    m_tree->Branch("yT",&vm_yT);
    m_tree->Branch("zT",&vm_zT);
    m_tree->Branch("plane",&vm_plane);
    m_tree->Branch("module",&vm_module);
    m_tree->Branch("sensor",&vm_sensor);
    m_tree->Branch("detX",&vm_detX);
    m_tree->Branch("detY",&vm_detY);
    m_tree->Branch("detZ",&vm_detZ);
    m_tree->Branch("pdg",&vm_pdg);
    m_tree->Branch("px",&vm_px);
    m_tree->Branch("py",&vm_py);
    m_tree->Branch("pz",&vm_pz);
    m_tree->Branch("p",&vm_p);
    m_tree->Branch("status",&vm_status);
    m_tree->Branch("beampx",&beam_px);
    m_tree->Branch("beampy",&beam_py);
    m_tree->Branch("beampz",&beam_pz);
    m_tree->Branch("beamp",&beam_p);
    m_tree->Branch("beam_pdg",&beam_pdg);
    m_tree->Branch("genPpx",&m_genPpx);
    m_tree->Branch("genPpy",&m_genPpy);
    m_tree->Branch("genPpz",&m_genPpz);
    m_tree->Branch("genPp",&m_genPp);
    m_tree->Branch("genBeamP",&m_genBeamP);
    m_tree->Branch("genBeamPx",&m_genBeamPx);
    m_tree->Branch("genBeamPy",&m_genBeamPy);
    m_tree->Branch("genBeamPz",&m_genBeamPz);

    m_tree->Branch("genBeamPP",&m_genBeamPP);
    m_tree->Branch("genBeamPPx",&m_genBeamPPx);
    m_tree->Branch("genBeamPPy",&m_genBeamPPy);
    m_tree->Branch("genBeamPPz",&m_genBeamPPz);


    m_tree->Branch("trk_p",     &trk_p);
m_tree->Branch("trk_theta", &trk_theta);
m_tree->Branch("trk_phi",   &trk_phi);
m_tree->Branch("trk_px", &trk_px);
m_tree->Branch("trk_py", &trk_py);
m_tree->Branch("trk_pz", &trk_pz);

    // Geometry and segmentation
    m_geoSvc = GetApplication()->GetService<DD4hep_service>();
    auto ro = m_geoSvc->detector()->readout("B0TrackerHits");

    m_seg = ro.segmentation().segmentation();
    
    m_decoder = ro.idSpec().decoder();

    auto& detector = dd4hep::Detector::getInstance();
    const dd4hep::VolumeManager& volman = detector.volumeManager();
    m_volman = dd4hep::Detector::getInstance().volumeManager();


}
void B0Trackers::MCgenAnalysis(const std::vector<const edm4hep::MCParticle*>& mcparts){
	m_VertexTrue.pdg.clear();
	m_VertexTrue.status.clear();
	m_VertexTrue.E.clear();
	m_VertexTrue.theta.clear();
	m_VertexTrue.phi.clear();
	m_VertexTrue.vx.clear();
	m_VertexTrue.vy.clear();
	m_VertexTrue.vz.clear();
	m_VertexTrue.px.clear();
	m_VertexTrue.py.clear();
	m_VertexTrue.pz.clear();
	m_VertexTrue.p.clear();

	for(auto particle : mcparts){
		edm4hep::Vector3d p = particle->getMomentum();
		edm4hep::Vector3d v = particle->getVertex();

		if((particle->getPDG() == 22 ||
		   particle->getPDG() == 11 ||
		   particle->getPDG() == -11 ||
            particle->getPDG() == 2212) &&
		   (particle->getGeneratorStatus() == 1 ||  particle->getGeneratorStatus() == 4)){

			m_VertexTrue.pdg.push_back(particle->getPDG());
			m_VertexTrue.status.push_back(particle->getGeneratorStatus());
			m_VertexTrue.E.push_back(particle->getEnergy());
			double theta = atan2( sqrt(pow(p.x,2) + pow(p.y,2)), p.z );
			m_VertexTrue.theta.push_back(theta);
			m_VertexTrue.phi.push_back(atan2( p.y, p.x ));
			m_VertexTrue.vx.push_back(v.x);
			m_VertexTrue.vy.push_back(v.y);
			m_VertexTrue.vz.push_back(v.z);
			m_VertexTrue.px.push_back(p.x);
			m_VertexTrue.py.push_back(p.y);
			m_VertexTrue.pz.push_back(p.z);
			m_VertexTrue.p.push_back(sqrt(pow(p.x,2) + pow(p.y,2) + pow(p.z,2)));
		}
	}

}



//void B0Trackers::ProcessSequential(const std::shared_ptr<const JEvent>& event) {
void B0Trackers::Process(const std::shared_ptr<const JEvent>& event){
   


// to see the available factories ====
/*
    for (auto* fac : event->GetAllFactories()) {
        std::cout << "Factory: object=\"" << fac->GetObjectName()
                  << "\" tag=\"" << fac->GetTag() << "\""
                  << std::endl;
    }
*/

auto mcparticles = event->Get<edm4hep::MCParticle>("MCParticles");
auto tracker_hits    = event->Get<edm4hep::SimTrackerHit>("B0TrackerHits");
//auto tracks  = event->Get<edm4eic::TrackParameters>("B0TrackerCKFTrackParameters");
auto tracks  = event->Get<edm4eic::TrackParameters>("B0TrackerCKFTruthSeededTrackParameters");
auto seeds = event->Get<edm4eic::TrackSeed>("B0TrackerTruthSeeds");



    std::array<std::vector<const edm4hep::SimTrackerHit*>,4> hitBuckets;
    m_HitsClas.clear();
    beam_px.clear();
    beam_py.clear();
    beam_pz.clear();
    beam_p.clear();
    beam_pdg.clear();

    vm_xR.clear();
    vm_yR.clear();
    vm_zR.clear();
    vm_xT.clear();
    vm_yT.clear();
    vm_zT.clear();
    vm_plane.clear();
    vm_module.clear();
    vm_sensor.clear();
    vm_detX.clear();
    vm_detY.clear();
    vm_detZ.clear();
    vm_pdg.clear();
    vm_px.clear();
    vm_py.clear();
    vm_pz.clear();
    vm_p.clear();
    vm_status.clear(); 


    m_genPpx.clear();
    m_genPpy.clear();
    m_genPpz.clear();
    m_genPp.clear();

m_genBeamP.clear();
m_genBeamPx.clear();
m_genBeamPy.clear();
m_genBeamPz.clear();


m_genBeamPP.clear();
m_genBeamPPx.clear();
m_genBeamPPy.clear();
m_genBeamPPz.clear();


    trk_p.clear();
    trk_px.clear();
    trk_py.clear();
    trk_pz.clear();
trk_theta.clear();
trk_phi.clear();

    //
    // ----- clear vectors -----
    //
    m_xi.clear();
    m_yi.clear();
    m_zi.clear();
    m_pxi.clear();
    m_pyi.clear();
    m_pzi.clear();
    m_pathL.clear();
    //-------------
 
    auto app = GetApplication();
    m_geoSvc = app->template GetService<DD4hep_service>();
    
    
//    const auto& simHits = *event->GetCollection<edm4hep::SimTrackerHit>("B0TrackerHits");
    auto simHits = event->Get<edm4hep::SimTrackerHit>("B0TrackerHits");


    std::cout<<" ============ "<<std::endl;
    int tempCount = 0;
//    for (const auto& h : simHits) {
    for (const auto *h:simHits){

    	HitClass tempHit;
        auto  mc = h->getParticle();
        if (!mc.isAvailable()) {
            // SimTrackerHit has no associated MCParticle; skip to avoid
            // contaminating the TTree with pdg=0, p=0 rows.
            continue;
        }
        int genStat = mc.getGeneratorStatus();
        m_primary = genStat;
        auto  mom = mc.getMomentum();
        m_Truepx = mom.x;
        m_Truepy = mom.y;
        m_Truepz = mom.z;

        m_Truep = std::sqrt(mom.x*mom.x + mom.y*mom.y + mom.z*mom.z);
        m_TruePDG = mc.getPDG();

        uint64_t cid = h->getCellID();

        // readout centre ------------------------------
        auto cpos = m_seg->position(cid);
        //auto local = m_seg->position(cid);

        const auto gpos = m_geoSvc->converter()->position(cid); // cm

    	const auto volman 	= m_geoSvc->detector()->volumeManager();
    	const auto alignment = volman.lookupDetElement(cid).nominal();
    	const auto lpos 	= alignment.worldToLocal( dd4hep::Position( gpos.x(), gpos.y(), gpos.z() ) ); // cm



/*
	static int debug_hit_prints = 0;
	if (debug_hit_prints < 1000) { // only print for first 10 hits total
    		auto de = volman.lookupDetElement(cid);
    		std::cout << "[DBG] cid=" << cid
              	<< " system=" << m_decoder->get(cid, "system")
              	<< " layer="  << m_decoder->get(cid, "layer")
              //	<< " plane="  << m_decoder->get(cid, "plane")
              	<< " module=" << m_decoder->get(cid, "module")
              	<< " sensor=" << m_decoder->get(cid, "sensor")
              	<< "  DetElement path='" << de.path() << "'"
              	<< " volume='" << de.volume().name() << "'"
              	<< std::endl;
    		++debug_hit_prints;
	}

*/


        m_xR = 10.*gpos.x();
        m_yR = 10.*gpos.y();
        m_zR = 10.*gpos.z();

        m_eDep = h->getEDep();
        m_time = h->getTime();
        m_HitPath = h->getPathLength();

	std::cout<<tempCount++<<" "<<genStat<<" p = "<<m_Truep<<" PDG = "<<m_TruePDG<<std::endl;
        // truth hit position --------------------------
        m_xT = h->getPosition().x;   // same as h.x()
        m_yT = h->getPosition().y;   // etc.
        m_zT = h->getPosition().z;

        // bit-fields (readout is CartesianGridXZ: system,layer,module,sensor,x,z)
        m_plane  = m_decoder->get(cid,"layer");
        m_module = m_decoder->get(cid,"module");
        m_sensor = m_decoder->get(cid,"sensor");
        m_pixX   = m_decoder->get(cid,"x");
        m_pixZ   = m_decoder->get(cid,"z");


	vm_xR.push_back(m_xR);
	vm_yR.push_back(m_yR);
	vm_zR.push_back(m_zR);
        vm_xT.push_back(m_xT);
        vm_yT.push_back(m_yT);
        vm_zT.push_back(m_zT);
	vm_plane.push_back(m_plane);
	vm_module.push_back(m_module);
	vm_sensor.push_back(m_sensor);
	vm_detX.push_back(lpos.x());
	vm_detY.push_back(lpos.y());
	vm_detZ.push_back(lpos.z());
	vm_pdg.push_back(mc.getPDG());
	vm_px.push_back(m_Truepx);
        vm_py.push_back(m_Truepy);
        vm_pz.push_back(m_Truepz);
        vm_p.push_back(m_Truep);
	vm_status.push_back(m_primary);
	
	
        tempHit.xR = m_xR;
        tempHit.yR = m_yR;
        tempHit.zR = m_zR;
        tempHit.xT = m_xT;
        tempHit.yT = m_yT;
        tempHit.zT = m_zT;
        //tempHit.plane = m_plane;
        tempHit.mod = m_module;
        tempHit.prime = m_primary;
        tempHit.eDep = m_eDep;
        tempHit.time = m_time;
        tempHit.pathL = m_HitPath;
        tempHit.cellID = cid;
        tempHit.px = m_Truepx;
        tempHit.py = m_Truepy;
        tempHit.pz = m_Truepz;
        tempHit.p = m_Truep;
        tempHit.theta = atan2( sqrt(pow(m_Truepx,2) + pow(m_Truepy,2)), m_Truepz );
        tempHit.phi = atan2( m_Truepy, m_Truepx );
        tempHit.xDet = lpos.x();
        tempHit.yDet = lpos.y();
        tempHit.zDet = lpos.z();



	m_detX = lpos.x();
	m_detY = lpos.y();
	m_detZ = lpos.z();

        m_HitsClas.push_back(tempHit);

	//BuildTracks(m_HitsClas);
        // ... fill your TTree or do residual = (xT-xR) etc.
        // --- write one row ---
	//
	
	
        /*MCgenAnalysis();

	for(int ii=0;ii<m_VertexTrue.status.size();ii++){
		if(m_VertexTrue.status.at(ii) == 4){
		}
		std::cout<<ii<<" stat = "<<m_VertexTrue.status.at(ii)<<std::endl;
		
	}

        m_tree->Fill();*/

	


    }

        MCgenAnalysis(mcparticles);


        for(int ii=0;ii<m_VertexTrue.status.size();ii++){
                if(m_VertexTrue.status.at(ii) == 4){
			beam_px.push_back(m_VertexTrue.px.at(ii));
			beam_py.push_back(m_VertexTrue.py.at(ii));
			beam_pz.push_back(m_VertexTrue.pz.at(ii));
			beam_p.push_back(m_VertexTrue.p.at(ii));
			beam_pdg.push_back(m_VertexTrue.pdg.at(ii));
                }

        }


// Access ACTS track parameters (one set per fitted track)
for (const auto *tp : tracks){//B0CKFTrackParams()) {

	cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<endl<<endl<<endl<<endl;
//    if (!tp) continue;  // just in case


    cout<<" =================\t\t\t\t\t\t\t\t\t\t>>>>>>>>>>>>>>>>>>>>>> +_+_++_+_+_++_+__ "<<endl;
    // Extract angles and q/p
//    auto mom = tp->getMomentum();
//    cout<<" track momentum  = "<<std::sqrt(mom.x*mom.x + mom.y*mom.y + mom.z*mom.z)<<endl;
    float theta  = tp->getTheta();
    float phi    = tp->getPhi();
    float qOverP = tp->getQOverP();

    cout<< "momentum  = "<<std::abs(1.0 / qOverP)<<endl;
    // reconstruct |p| in GeV (ignore sign)
    double p = (qOverP != 0.f) ? std::abs(1.0 / qOverP) : 0.0;

    trk_p.push_back(p);
    trk_theta.push_back(theta);
    trk_phi.push_back(phi);

    // If you later want components:
    double px = p * std::sin(theta) * std::cos(phi);
    double py = p * std::sin(theta) * std::sin(phi);
    double pz = p * std::cos(theta);

    trk_px.push_back(px);
    trk_py.push_back(py);
    trk_pz.push_back(pz);
}


// generate generated:

	for(auto particle : mcparticles){
		edm4hep::Vector3d p = particle->getMomentum();
		edm4hep::Vector3d v = particle->getVertex();

		if(particle->getPDG() == 2212 && particle->getGeneratorStatus() == 4){

			m_genPpx.push_back(p.x);
			m_genPpy.push_back(p.y);
			m_genPpz.push_back(p.z);
			m_genPp.push_back(sqrt(pow(p.x,2)+pow(p.y,2)+pow(p.z,2)));
		}
		if(particle->getPDG() == 2212 && particle->getGeneratorStatus() == 1){
			m_genBeamPPx.push_back(p.x);
                        m_genBeamPPy.push_back(p.y);
                        m_genBeamPPz.push_back(p.z);
			m_genBeamPP.push_back(sqrt(pow(p.x,2)+pow(p.y,2)+pow(p.z,2)));
		}
		if(particle->getPDG() == 11 && particle->getGeneratorStatus() == 1){
			m_genBeamPx.push_back(p.x);
                        m_genBeamPy.push_back(p.y);
                        m_genBeamPz.push_back(p.z);
                        m_genBeamP.push_back(sqrt(pow(p.x,2)+pow(p.y,2)+pow(p.z,2)));
		}
	}	


        m_tree->Fill();

    
}
//-------------------------------------------
// FinishWithGlobalRootLock
//-------------------------------------------
//void B0Trackers::FinishWithGlobalRootLock() {
void B0Trackers::Finish(){

}

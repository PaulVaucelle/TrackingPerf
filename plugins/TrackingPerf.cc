// system include files
#include <memory>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <utility>
#include <TNtuple.h>
#include <bitset>
#include <memory>

// user include files
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "boost/functional/hash.hpp"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"

#include "CondFormats/GeometryObjects/interface/HcalParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//---------------------------------------------------------------------//

#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DetectorDescription/Core/interface/DDCompactView.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/RegexMatch.h"

#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"

#include "Geometry/HcalCommonData/interface/HcalParametersFromDD.h"
#include "Geometry/Records/interface/HcalParametersRcd.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "RecoLocalTracker/SiStripClusterizer/interface/SiStripClusterInfo.h"

#include "RecoTracker/FinalTrackSelectors/plugins/getBestVertex.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//$$
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
//$$
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimTracker/TrackAssociation/interface/TrackingParticleIP.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"

#include "TrackingPerf/TrackingPerf/interface/Proto.h"
#include "TrackingPerf/TrackingPerf/interface/DeltaFunc.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;
using TrackingParticleRefKeyToIndex = std::unordered_map<reco::RecoToSimCollection::index_type, size_t>;


class TrackingPerf : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
    explicit TrackingPerf(const edm::ParameterSet&);
    ~TrackingPerf();
    
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    
    // ----------member data ---------------------------
    edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
    
    // ----------member data ---------------------------
    
    void clearVariables();
    
    void getHLTPathNames(const edm::TriggerNames& triggerNames);
    void fillTriggerBits(std::vector<std::string>& HLTPathsByName_, edm::Handle<edm::TriggerResults>& triggerBits, const edm::TriggerNames &names);
    
    edm::InputTag tkTraj_;
    
    //------------------------------------
    //// GETTOKEN DECLARATIONS
    //------------------------------------
    
    std::string weightFile_;
    
    //------------------------------------
    // track ( and event ) information
    //------------------------------------
    // const edm::EDGetTokenT<edm::View<reco::Track> > trackToken_;
    const edm::EDGetTokenT<edm::View<pat::PackedCandidate> > trackToken_;
    // const edm::EDGetTokenT<edm::View<reco::Track> > trackSrc_;
    const edm::EDGetTokenT<edm::View<pat::PackedCandidate> > trackSrc_;
    const edm::EDGetTokenT<std::vector<Trajectory> > trajSrc_;
    const edm::EDGetTokenT<TrajTrackAssociationCollection> trajTrackAssociationSrc_;
    const edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> trackAssociatorToken_;
    const edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;
    
    const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
    
    edm::EDGetTokenT<TrackingParticleCollection> trackingParticleToken_;
    edm::EDGetTokenT<TrackingParticleRefVector> trackingParticleRefToken_;
    
    //------------------------------------
    // jet information
    //------------------------------------
    // const edm::EDGetTokenT<edm::View<reco::Jet> > ak4slimmedJetToken_;
    // const edm::EDGetTokenT<edm::View<reco::Jet> > ak4PFJetToken_;
    // const edm::EDGetTokenT<edm::View<reco::Jet> > CaloJetToken_;
    // const edm::EDGetTokenT<edm::View<reco::Jet> > ak8jetToken_;
    //   const edm::EDGetTokenT<edm::View<reco::Jet> > ak8CaloJetToken_;

    const edm::EDGetTokenT<edm::View<pat::Jet> > ak4slimmedJetToken_;
    const edm::EDGetTokenT<edm::View<pat::Jet> > ak8jetToken_;
    
    //------------------------------------
    // gen information
    //------------------------------------
    const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    const edm::EDGetTokenT<edm::View<reco::GenJet> >     genJetToken_;
    const edm::EDGetTokenT<reco::GenJetCollection>      ak8GenJetToken_;
    const edm::EDGetTokenT<reco::GenJetCollection>      genTTXJetsToken_;
    const edm::EDGetTokenT<GenEventInfoProduct>         genEventInfoToken_;
    const edm::EDGetTokenT<LHEEventProduct>             LHEEventProductToken_;
    
    //pat information
    const edm::EDGetTokenT<edm::View<pat::PackedCandidate> > pfcandsToken_;
    
    std::string parametersDefinerName_;
    
    //------------------------------------
    // electrons
    //------------------------------------
    const edm::EDGetTokenT<pat::ElectronCollection> electronPATToken_;
    //------------------------------------
    // muons
    //------------------------------------
    const edm::EDGetTokenT<pat::MuonCollection> slimmedmuonToken_;
    //------------------------------------
    // MET
    //------------------------------------
    const edm::EDGetTokenT<pat::METCollection> metToken_;
    //------------------------------------
    
    //------------------------------------
    // Zmu skim bit
    //------------------------------------
    std::vector<std::string> filterTriggerNames_;
    std::vector<std::string> HLTPathsByName_;
    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    //reco dataformat//
    // const edm::EDGetTokenT<vector<reco::CompositeCandidate> > ZmumuCandidate_;
    
    //additional
    bool useCluster_;
    bool runOnData_;
    
    TTree *smalltree;
    
    edm::Service<TFileService> fs;
    
    std::string ttrhbuilder_;
    
    edm::ESHandle<MagneticField> bField;
    
    edm::ParameterSet kvfPSet;
       
    const edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
    
    /////// VARIABLES FOR FILLING THE TREE
    
    int runNumber, eventNumber, lumiBlock;
    int  tree_NbrOfZCand;
    bool tree_passesHTFilter;
    int  tree_nTracks; 
    int  nBC = 0, nFromC = 0, nFromB = 0; 
    
    //-----------------------
    //trigger variable
    std::vector<string > tree_trigger_names;
    std::vector<bool >   tree_trigger_bits;
    
    //-----------------------
    //fill the tree per track
    //std::vector<int> tree_track_nclusters;
    std::vector< float > tree_track_pt;
    // std::vector< float > tree_track_ptError;
    std::vector< float > tree_track_outerPt;
    std::vector< float > tree_track_eta;
    std::vector< float > tree_track_phi;
    std::vector<int>     tree_track_charge;
    std::vector<int>     tree_track_nhits;
    std::vector<float >  tree_track_NChi2;
    //std::vector<float >  tree_track_Quality;
    std::vector<bool >   tree_track_isHighQuality;
    std::vector<bool >   tree_track_isLoose;
    std::vector<bool >   tree_track_isTight;
    std::vector< float>  tree_track_dxy; // dxy with respect to beam spot position
    std::vector< float>  tree_track_dxyError;
    std::vector< float>  tree_track_dz;
    std::vector< float>  tree_track_dzError ;
    std::vector<int>     tree_track_numberOfLostHits;
    std::vector<int>     tree_track_numberOfValidHits;
    std::vector<unsigned int>    tree_track_originalAlgo; // definition as comments at the end of the file,
    //http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_1_3/doc/html/d8/df2/classreco_1_1TrackBase.html#aca7611bd1a33d535cefc72b6e497ece8
    std::vector<unsigned int>    tree_track_algo;
    std::vector<unsigned short>  tree_track_stopReason;
    
    std::vector<int>     tree_track_numberOfValidPixelHits;
    std::vector<int>     tree_track_numberOfValidStripHits;
    std::vector<int>     tree_track_numberOfValidStripTIBHits;
    std::vector<int>     tree_track_numberOfValidStripTIDHits;
    std::vector<int>     tree_track_numberOfValidStripTOBHits;
    std::vector<int>     tree_track_numberOfValidStripTECHits;
    std::vector<int>     tree_track_numberOfValidPixelBarrelHits;
    std::vector<int>     tree_track_numberOfValidPixelEndcapHits;
    std::vector<int>     tree_track_hasValidHitInPixelLayer;
    std::vector<int>     tree_track_trackerLayersWithMeasurement;
    std::vector<int>     tree_track_pixelLayersWithMeasurement;
    std::vector<int>     tree_track_stripTECLayersWithMeasurement ;
    std::vector<int>     tree_track_stripTIBLayersWithMeasurement;
    std::vector<int>     tree_track_stripTIDLayersWithMeasurement;
    std::vector<int>     tree_track_stripTOBLayersWithMeasurement;
    //  std::vector<int>      tree_track_nPixel;
    //  std::vector<int>      tree_track_nStrip;
    
    std::vector< float >  tree_track_vx;
    std::vector< float >  tree_track_vy;
    std::vector< float >  tree_track_vz;
    std::vector<float>    tree_track_firsthit_X;
    std::vector<float>    tree_track_firsthit_Y;
    std::vector<float>    tree_track_firsthit_Z;
    // std::vector<float>    tree_track_firsthit_phi;
    std::vector<float>    tree_track_ntrk10;
    std::vector<float>    tree_track_ntrk20;
    std::vector<float>    tree_track_ntrk30;
    
    std::vector<double>    tree_track_MVAval;
//$$
    std::vector<int>       tree_track_Hemi;
    std::vector<double>    tree_track_Hemi_dR;
    std::vector<double>    tree_track_Hemi_mva_NChi2;
    std::vector<int>       tree_track_Hemi_LLP;
//$$
    
    std::vector<int>    tree_track_recoVertex_idx;
    std::vector<int>    tree_track_recoAK4SlimmedJet_idx;
    std::vector<int>    tree_track_recoAK4PFJet_idx;
    std::vector<int>    tree_track_reco08Jet_idx;
    std::vector<int>    tree_track_recoCaloJet_idx;
    //   std::vector<int>      tree_track_reco08CaloJet_idx;
    
    std::vector<int>      tree_track_nSimHits;
    std::vector<bool>     tree_track_isSimMatched;
    
    std::vector< int >    tree_track_sim_charge;
    std::vector< float >  tree_track_sim_pt;
    std::vector< float >  tree_track_sim_eta  ;
    std::vector< float >  tree_track_sim_phi  ;
    std::vector<bool>     tree_track_sim_longLived      ;
    // std::vector<int>   tree_track_sim_matchedHit    ;
    std::vector<int>      tree_track_sim_pdgId;
    std::vector<int>      tree_track_sim_numberOfTrackerHits  ;
    std::vector<int>      tree_track_sim_numberOfTrackerLayers;
    std::vector<float>    tree_track_sim_mass  ;
    std::vector<int>      tree_track_sim_status;
    
    std::vector<float>    tree_track_sim_vx;
    std::vector<float>    tree_track_sim_vy;
    std::vector<float>    tree_track_sim_vz;
    std::vector<int>      tree_track_sim_isFromLLP;
//$$
    std::vector<int>      tree_track_sim_isFromBC;
    std::vector<int>      tree_track_sim_isFromBC_mother_pdgId;
    std::vector<int>      tree_track_sim_isFromBC_LLP;
    std::vector<float>    tree_track_sim_dV;
    std::vector<float>    tree_track_sim_dist;
//$$
    
    //--------------------------------
    // Tracking Particle infos -------
    //--------------------------------
    // std::vector< int >           tree_simtrack_charge;
    // std::vector< float >         tree_simtrack_pt;
    // std::vector< float >         tree_simtrack_eta;
    // std::vector< float >         tree_simtrack_phi;
    // std::vector<bool>            tree_simtrack_longLived;
    // std::vector<int>             tree_simtrack_pdgId;
    // std::vector<int>             tree_simtrack_numberOfTrackerHits;
    // std::vector<int>             tree_simtrack_numberOfTrackerLayers;
    // std::vector<float>           tree_simtrack_mass;
    // std::vector<int>             tree_simtrack_status;
    
    // std::vector<float>           tree_simtrack_vx;
    // std::vector<float>           tree_simtrack_vy;
    // std::vector<float>           tree_simtrack_vz;
    // std::vector<bool>            tree_simtrack_isRecoMatched;
    // std::vector<float>           tree_simtrack_pca_dxy;
    // std::vector<float>           tree_simtrack_pca_dz   ;
    // std::vector<std::vector<int> > tree_simtrack_trkIdx;
    // std::vector<float>           tree_simtrack_isRecoMatched_pt;
    
    //--------------------------------
    // Beam spot infos -------
    //--------------------------------
    float       tree_bs_PosX;
    float       tree_bs_PosY;
    float       tree_bs_PosZ;
    
    //--------------------------------
    // vertex infos -------
    //--------------------------------
    // for the main vertex
    std::vector<float> tree_vtx_PosX;
    std::vector<float> tree_vtx_PosY;
    std::vector<float> tree_vtx_PosZ;
    std::vector<float> tree_vtx_NChi2;
    std::vector<float> tree_vtx_PosXError;
    std::vector<float> tree_vtx_PosYError;
    std::vector<float> tree_vtx_PosZError;
    
    //--------------------------------
    // jet infos -------
    //--------------------------------
    
    std::vector<float> tree_AK4Slimmedjet_E;
    std::vector<float> tree_AK4Slimmedjet_pt;
    std::vector<float> tree_AK4Slimmedjet_eta;
    std::vector<float> tree_AK4Slimmedjet_phi;
//$$    std::vector<int>   tree_AK4Slimmedjet_idxTrack;
    
    std::vector<float> tree_AK4PFjet_E;
    std::vector<float> tree_AK4PFjet_pt;
    std::vector<float> tree_AK4PFjet_eta;
    std::vector<float> tree_AK4PFjet_phi;
//$$    std::vector<int>   tree_AK4PFjet_idxTrack;
    
    std::vector<float> tree_CaloJet_E;
    std::vector<float> tree_CaloJet_pt;
    std::vector<float> tree_CaloJet_eta;
    std::vector<float> tree_CaloJet_phi;
//$$    std::vector<int>   tree_CaloJet_idxTrack;
    
    std::vector<float> tree_jet08_E;
    std::vector<float> tree_jet08_pt;
    std::vector<float> tree_jet08_eta;
    std::vector<float> tree_jet08_phi;
//$$    std::vector<int>   tree_jet08_idxTrack;
    
    //--------------------------------
    // met infos -------
    //--------------------------------
    float tree_PFMet_et;
    float tree_PFMet_phi;
    float tree_PFMet_sig;
    
    //--------------------------------
    // gen infos -------
    //--------------------------------
    float tree_GenPVx;
    float tree_GenPVy;
    float tree_GenPVz;
    
    std::vector< float > tree_genParticle_pt;
    std::vector< float > tree_genParticle_eta;
    std::vector< float > tree_genParticle_phi;
    std::vector< float > tree_genParticle_charge;
    std::vector< int >   tree_genParticle_pdgId;
    std::vector< float > tree_genParticle_vx;
    std::vector< float > tree_genParticle_vy;
    std::vector< float > tree_genParticle_vz;
    std::vector< float > tree_genParticle_mass;
    std::vector< int >   tree_genParticle_statusCode;
    std::vector< int >   tree_genParticle_mother_pdgId;

//$$
    std::vector< float > tree_genFromC_pt;
    std::vector< float > tree_genFromC_eta;
    std::vector< float > tree_genFromC_phi;
    std::vector< float > tree_genFromC_charge;
    std::vector< int >   tree_genFromC_pdgId;
    std::vector< float > tree_genFromC_vx;
    std::vector< float > tree_genFromC_vy;
    std::vector< float > tree_genFromC_vz;
    std::vector< int >   tree_genFromC_mother_pdgId;
    std::vector< int >   tree_genFromC_generation;
    std::vector< int >   tree_genFromC_LLP;

    std::vector< float > tree_genFromB_pt;
    std::vector< float > tree_genFromB_eta;
    std::vector< float > tree_genFromB_phi;
    std::vector< float > tree_genFromB_charge;
    std::vector< int >   tree_genFromB_pdgId;
    std::vector< float > tree_genFromB_vx;
    std::vector< float > tree_genFromB_vy;
    std::vector< float > tree_genFromB_vz;
    std::vector< int >   tree_genFromB_mother_pdgId;
    std::vector< int >   tree_genFromB_generation;
    std::vector< int >   tree_genFromB_LLP;
//$$
   
    //--------------------------------
    // gen jet infos -------
    //--------------------------------
    std::vector<float> tree_genJet_pt;
    std::vector<float> tree_genJet_eta;
    std::vector<float> tree_genJet_phi;
    std::vector<float> tree_genJet_mass;
    std::vector<float> tree_genJet_energy;
    
    std::vector<float> tree_ak8GenJet_pt;
    std::vector<float> tree_ak8GenJet_eta;
    std::vector<float> tree_ak8GenJet_phi;
    std::vector<float> tree_ak8GenJet_mass;
    std::vector<float> tree_ak8GenJet_energy;
    
    //--------------------------------
    // gen event info -------
    //--------------------------------
    
    //--------------------------------
    // lhe event infos -------
    //--------------------------------
    
    //--------------------------------
    // PF infos -------
    //--------------------------------
    
    //--------------------------------
    // electrons infos -------
    //--------------------------------
    std::vector<float> tree_electron_pt;
    std::vector<float> tree_electron_eta;
    std::vector<float> tree_electron_phi;
    std::vector<float> tree_electron_vx;
    std::vector<float> tree_electron_vy;
    std::vector<float> tree_electron_vz;
    std::vector<float> tree_electron_energy;
    std::vector< int > tree_electron_charge;
    
    //--------------------------------
    // muons infos -------
    //--------------------------------
    // slimmed muons
    std::vector<float> tree_slimmedmuon_pt;
    std::vector<float> tree_slimmedmuon_eta;
    std::vector<float> tree_slimmedmuon_phi;
    std::vector<float> tree_slimmedmuon_vx;
    std::vector<float> tree_slimmedmuon_vy;
    std::vector<float> tree_slimmedmuon_vz;
    std::vector<float> tree_slimmedmuon_energy;
    std::vector<float> tree_slimmedmuon_dxy;
    std::vector<float> tree_slimmedmuon_dxyError;
    std::vector<float> tree_slimmedmuon_dz;
    std::vector<float> tree_slimmedmuon_dzError;
    std::vector< int > tree_slimmedmuon_charge;
    std::vector<bool>  tree_slimmedmuon_PFisoVeryTight;
    std::vector<bool>  tree_slimmedmuon_PFisoTight;
    std::vector<bool>  tree_slimmedmuon_PFisoMedium;
    std::vector<bool>  tree_slimmedmuon_PFisoLoose;
    std::vector<bool>  tree_slimmedmuon_MVAisoLoose;
    std::vector<bool>  tree_slimmedmuon_MVAisoMedium;
    std::vector<bool>  tree_slimmedmuon_MVAisoTight;
    std::vector<bool>  tree_slimmedmuon_isGlobalMuon;
    std::vector<bool>  tree_slimmedmuon_isStandAloneMuon;
    std::vector<bool>  tree_slimmedmuon_CutBasedIdLoose;
    std::vector<bool>  tree_slimmedmuon_CutBasedIdMedium;
    std::vector<bool>  tree_slimmedmuon_CutBasedIdMediumPrompt;
    std::vector<bool>  tree_slimmedmuon_CutBasedIdTight;
    
    int nEvent ;
    
//$$
    //-----------------------
    // generated LLPs 
    //-----------------------
    float LLP1_pt, LLP1_eta, LLP1_phi, LLP2_pt, LLP2_eta, LLP2_phi;
    float LLP1_x, LLP1_y, LLP1_z, LLP2_x, LLP2_y, LLP2_z;
    float LLP1_dist, LLP2_dist;
    int   LLP1_nTrks = 0, LLP2_nTrks = 0;
    int   tree_nLLP = -1;

    std::vector< int >   tree_LLP;
    std::vector< float > tree_LLP_pt;
    std::vector< float > tree_LLP_eta;
    std::vector< float > tree_LLP_phi;
    std::vector< float > tree_LLP_x;
    std::vector< float > tree_LLP_y;
    std::vector< float > tree_LLP_z;
    std::vector< int >   tree_LLP_nTrks;
    std::vector< int >   tree_LLP_Vtx_nTrks;
    std::vector< float > tree_LLP_Vtx_NChi2;
    std::vector< float > tree_LLP_Vtx_dx;
    std::vector< float > tree_LLP_Vtx_dy;
    std::vector< float > tree_LLP_Vtx_dz;
    
    //-----------------------
    //Analysis with the two hemispheres
    //-----------------------
    std::vector< int >   tree_Hemi;
    std::vector< int >   tree_Hemi_njet;
    std::vector< float > tree_Hemi_eta;
    std::vector< float > tree_Hemi_phi;
    std::vector< float > tree_Hemi_dR;
    std::vector< int >   tree_Hemi_nTrks;
    std::vector< int >   tree_Hemi_nTrks_sig;
    std::vector< int >   tree_Hemi_nTrks_bad;
    std::vector< int >   tree_Hemi_LLP;
    std::vector< float > tree_Hemi_LLP_pt;
    std::vector< float > tree_Hemi_LLP_eta;
    std::vector< float > tree_Hemi_LLP_phi;
    std::vector< float > tree_Hemi_LLP_dist;
    std::vector< float > tree_Hemi_LLP_x;
    std::vector< float > tree_Hemi_LLP_y;
    std::vector< float > tree_Hemi_LLP_z;
    std::vector< float > tree_Hemi_Vtx_NChi2;
    std::vector< int >   tree_Hemi_Vtx_nTrks;
    std::vector< int >   tree_Hemi_Vtx_nTrks_sig;
    std::vector< int >   tree_Hemi_Vtx_nTrks_bad;
    std::vector< float > tree_Hemi_Vtx_x;
    std::vector< float > tree_Hemi_Vtx_y;
    std::vector< float > tree_Hemi_Vtx_z;
    std::vector< float > tree_Hemi_Vtx_dx;
    std::vector< float > tree_Hemi_Vtx_dy;
    std::vector< float > tree_Hemi_Vtx_dz;
    std::vector< float > tree_Hemi_dR12;
    std::vector< float > tree_Hemi_LLP_dR12;
//$$

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//


//
// constructors and destructor
//
TrackingPerf::TrackingPerf(const edm::ParameterSet& iConfig):

weightFile_( iConfig.getUntrackedParameter<std::string>("weightFileMVA") ),

// trackToken_(       consumes<edm::View<reco::Track> >(           iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
trackToken_(       consumes<edm::View<pat::PackedCandidate> >(           iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
// trackSrc_(         consumes<edm::View<reco::Track> >(           iConfig.getParameter<edm::InputTag>("trackLabel") )),
trackSrc_(         consumes<edm::View<pat::PackedCandidate> >(           iConfig.getParameter<edm::InputTag>("trackLabel") )),
trackAssociatorToken_( consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssociator"))),
puInfoToken_(consumes<vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileup"))),
beamSpotToken_(     consumes<reco::BeamSpot>(               iConfig.getUntrackedParameter<edm::InputTag>("beamSpot"))),
vertexToken_(      consumes<reco::VertexCollection>(           iConfig.getUntrackedParameter<edm::InputTag>("vertices"))),
ak4slimmedJetToken_(   consumes<edm::View<pat::Jet> >(           iConfig.getParameter<edm::InputTag>("ak4slimmedJetInput"))),
// ak4PFJetToken_(     consumes<edm::View<pat::Jet> >(           iConfig.getParameter<edm::InputTag>("ak4PFJetInput"))),
// CaloJetToken_(     consumes<edm::View<pat::Jet> >(           iConfig.getParameter<edm::InputTag>("caloJetInput"))),
ak8jetToken_(      consumes<edm::View<pat::Jet> >(           iConfig.getParameter<edm::InputTag>("ak8jetInput"))),
//   ak8CaloJetToken_(     consumes<edm::View<pat::Jet> >(           iConfig.getParameter<edm::InputTag>("ak8CaloJetInput"))),
genParticlesToken_(    consumes<reco::GenParticleCollection>(            iConfig.getParameter<edm::InputTag>("genParticles"))),
genJetToken_(          consumes<edm::View<reco::GenJet>>(                 iConfig.getParameter<edm::InputTag>("genJetInput"))),
ak8GenJetToken_(          consumes<reco::GenJetCollection>(              iConfig.getParameter<edm::InputTag>("ak8GenJetInput"))),
genEventInfoToken_(    consumes<GenEventInfoProduct>(                    iConfig.getParameter<edm::InputTag>("genEventInfoInput"))),
LHEEventProductToken_( consumes<LHEEventProduct>(                        iConfig.getParameter<edm::InputTag>("LHEEventProductInput"))),
pfcandsToken_(         consumes<edm::View<pat::PackedCandidate>>(         iConfig.getParameter<edm::InputTag>("pfcands"))),
parametersDefinerName_(                           iConfig.getUntrackedParameter<std::string>("parametersDefiner")),
electronPATToken_(     consumes<pat::ElectronCollection>(                iConfig.getParameter<edm::InputTag>("electronInput"))),
slimmedmuonToken_(     consumes<pat::MuonCollection>(               iConfig.getParameter<edm::InputTag>("slimmedmuonInput"))),
//  recomuonToken_(        consumes<reco::MuonCollection>(               iConfig.getParameter<edm::InputTag>("recomuonInput"))),
metToken_(             consumes<pat::METCollection>(               iConfig.getParameter<edm::InputTag>("metInput"))),
filterTriggerNames_(                                                     iConfig.getUntrackedParameter<std::vector<std::string> >("filterTriggerNames")),
triggerBits_(          consumes<edm::TriggerResults>(                    edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("HLT"))) ),
// ZmumuCandidate_(       consumes<vector<reco::CompositeCandidate>>(       iConfig.getParameter<edm::InputTag>("Zmumucand"))),
useCluster_(                                         iConfig.getUntrackedParameter<bool>("useCluster")),
runOnData_ (                                           iConfig.getUntrackedParameter<bool>("runOnData")),
kvfPSet( iConfig.getParameter<edm::ParameterSet>("KVFParameters"))
{
    
    nEvent = 0;
    usesResource("TFileService");
    
    //trackSrc_ = consumes<reco::TrackCollection>(src_);
    ttrhbuilder_ = iConfig.getParameter<std::string>("TTRHBuilder");
    
    smalltree = fs->make<TTree>("ttree", "ttree");
    
    // event info
    smalltree->Branch("runNumber",  &runNumber,  "runNumber/I");
    smalltree->Branch("eventNumber",&eventNumber,"eventNumber/I");
    smalltree->Branch("lumiBlock"  ,&lumiBlock,  "lumiBlock/I");
    smalltree->Branch("tree_bs_PosX", &tree_bs_PosX,  "tree_bs_PosX/F"  );
    smalltree->Branch("tree_bs_PosY", &tree_bs_PosY,  "tree_bs_PosY/F" );
    smalltree->Branch("tree_bs_PosZ", &tree_bs_PosZ,  "tree_bs_PosZ/F"  );
    
    // primary vertex info
    smalltree->Branch("tree_vtx_PosX", &tree_vtx_PosX);
    smalltree->Branch("tree_vtx_PosY", &tree_vtx_PosY);
    smalltree->Branch("tree_vtx_PosZ", &tree_vtx_PosZ);
    smalltree->Branch("tree_vtx_NChi2", &tree_vtx_NChi2);
    smalltree->Branch("tree_vtx_PosXError", &tree_vtx_PosXError);
    smalltree->Branch("tree_vtx_PosYError", &tree_vtx_PosYError);
    smalltree->Branch("tree_vtx_PosZError", &tree_vtx_PosZError);
    
    // trigger info
    smalltree->Branch("tree_trigger_names", &tree_trigger_names);
    smalltree->Branch("tree_trigger_bits",  &tree_trigger_bits);
    
    smalltree->Branch("tree_NbrOfZCand",  &tree_NbrOfZCand,  "tree_NbrOfZCand/I");
    smalltree->Branch("tree_passesHTFilter", &tree_passesHTFilter);
    smalltree->Branch("tree_nTracks", &tree_nTracks, "tree_nTracks/I"); 
    
    // track
    smalltree->Branch("tree_track_pt",            &tree_track_pt);
    // smalltree->Branch("tree_track_ptError",&tree_track_ptError);
    smalltree->Branch("tree_track_outerPt",           &tree_track_outerPt );
    smalltree->Branch("tree_track_eta",                  &tree_track_eta );
    smalltree->Branch("tree_track_phi",                  &tree_track_phi );
    smalltree->Branch("tree_track_charge",            &tree_track_charge );
    smalltree->Branch("tree_track_nhits",                &tree_track_nhits);
    smalltree->Branch("tree_track_NChi2",                &tree_track_NChi2);
    smalltree->Branch("tree_track_isHighQuality",     &tree_track_isHighQuality);
    smalltree->Branch("tree_track_isLoose",        &tree_track_isLoose);
    smalltree->Branch("tree_track_isTight",        &tree_track_isTight);
    smalltree->Branch("tree_track_dxy",                  &tree_track_dxy );
    smalltree->Branch("tree_track_dxyError",        &tree_track_dxyError);
    smalltree->Branch("tree_track_dz",            &tree_track_dz);
    smalltree->Branch("tree_track_dzError",        &tree_track_dzError  );
    smalltree->Branch("tree_track_numberOfLostHits",  &tree_track_numberOfLostHits);
    smalltree->Branch("tree_track_numberOfValidHits", &tree_track_numberOfValidHits);
    smalltree->Branch("tree_track_originalAlgo",      &tree_track_originalAlgo);
    smalltree->Branch("tree_track_algo",             &tree_track_algo);
    smalltree->Branch("tree_track_stopReason",        &tree_track_stopReason);
    smalltree->Branch("tree_track_isSimMatched",      &tree_track_isSimMatched    );
    
    smalltree->Branch("tree_track_numberOfValidPixelHits",        &tree_track_numberOfValidPixelHits);
    smalltree->Branch("tree_track_numberOfValidStripHits",        &tree_track_numberOfValidStripHits);
    smalltree->Branch("tree_track_numberOfValidStripTIBHits",     &tree_track_numberOfValidStripTIBHits);
    smalltree->Branch("tree_track_numberOfValidStripTIDHits",     &tree_track_numberOfValidStripTIDHits);
    smalltree->Branch("tree_track_numberOfValidStripTOBHits",     &tree_track_numberOfValidStripTOBHits);
    smalltree->Branch("tree_track_numberOfValidStripTECHits",     &tree_track_numberOfValidStripTECHits);
    smalltree->Branch("tree_track_numberOfValidPixelBarrelHits",  &tree_track_numberOfValidPixelBarrelHits);
    smalltree->Branch("tree_track_numberOfValidPixelEndcapHits",  &tree_track_numberOfValidPixelEndcapHits);
    smalltree->Branch("tree_track_hasValidHitInPixelLayer",       &tree_track_hasValidHitInPixelLayer);
    smalltree->Branch("tree_track_trackerLayersWithMeasurement",  &tree_track_trackerLayersWithMeasurement);
    smalltree->Branch("tree_track_pixelLayersWithMeasurement",    &tree_track_pixelLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTECLayersWithMeasurement", &tree_track_stripTECLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTIBLayersWithMeasurement", &tree_track_stripTIBLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTIDLayersWithMeasurement", &tree_track_stripTIDLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTOBLayersWithMeasurement", &tree_track_stripTOBLayersWithMeasurement);
    
    smalltree->Branch("tree_track_vx",           &tree_track_vx );
    smalltree->Branch("tree_track_vy",           &tree_track_vy );
    smalltree->Branch("tree_track_vz",           &tree_track_vz );
    smalltree->Branch("tree_track_firsthit_X",   &tree_track_firsthit_X);
    smalltree->Branch("tree_track_firsthit_Y",   &tree_track_firsthit_Y);
    smalltree->Branch("tree_track_firsthit_Z",   &tree_track_firsthit_Z);
    // smalltree->Branch("tree_track_firsthit_phi", &tree_track_firsthit_phi);
    smalltree->Branch("tree_track_ntrk10",&tree_track_ntrk10);
    smalltree->Branch("tree_track_ntrk20",&tree_track_ntrk20);
    smalltree->Branch("tree_track_ntrk30",&tree_track_ntrk30);
    
    smalltree->Branch("tree_track_MVAval", &tree_track_MVAval);
//$$
    smalltree->Branch("tree_track_Hemi", &tree_track_Hemi);
    smalltree->Branch("tree_track_Hemi_dR", &tree_track_Hemi_dR);
    smalltree->Branch("tree_track_Hemi_mva_NChi2", &tree_track_Hemi_mva_NChi2);
    smalltree->Branch("tree_track_Hemi_LLP", &tree_track_Hemi_LLP);
//$$
    
    smalltree->Branch("tree_track_recoVertex_idx", &tree_track_recoVertex_idx);
    smalltree->Branch("tree_track_recoAK4SlimmedJet_idx", &tree_track_recoAK4SlimmedJet_idx);
    smalltree->Branch("tree_track_recoAK4PFJet_idx",      &tree_track_recoAK4PFJet_idx);
    smalltree->Branch("tree_track_reco08Jet_idx",         &tree_track_reco08Jet_idx);
    smalltree->Branch("tree_track_recoCaloJet_idx",       &tree_track_recoCaloJet_idx);
    //   smalltree->Branch("tree_track_reco08CaloJet_idx",     &tree_track_reco08CaloJet_idx);
    
    // info about the simulated track matched to the reco track
    smalltree->Branch("tree_track_sim_charge",                &tree_track_sim_charge );
    smalltree->Branch("tree_track_sim_pt",                    &tree_track_sim_pt );
    smalltree->Branch("tree_track_sim_eta",                   &tree_track_sim_eta  );
    smalltree->Branch("tree_track_sim_phi",                   &tree_track_sim_phi  );
    smalltree->Branch("tree_track_sim_longLived",             &tree_track_sim_longLived );
    smalltree->Branch("tree_track_sim_pdgId",                 &tree_track_sim_pdgId );
    smalltree->Branch("tree_track_sim_numberOfTrackerHits",   &tree_track_sim_numberOfTrackerHits   );
    smalltree->Branch("tree_track_sim_numberOfTrackerLayers", &tree_track_sim_numberOfTrackerLayers );
    smalltree->Branch("tree_track_sim_mass",                  &tree_track_sim_mass   );
    smalltree->Branch("tree_track_sim_status",                &tree_track_sim_status );
    smalltree->Branch("tree_track_sim_vx",                    &tree_track_sim_vx );
    smalltree->Branch("tree_track_sim_vy",                    &tree_track_sim_vy );
    smalltree->Branch("tree_track_sim_vz",                    &tree_track_sim_vz );
    smalltree->Branch("tree_track_sim_isFromLLP",             &tree_track_sim_isFromLLP );
//$$
    smalltree->Branch("tree_track_sim_isFromBC",              &tree_track_sim_isFromBC );
    smalltree->Branch("tree_track_sim_isFromBC_mother_pdgId", &tree_track_sim_isFromBC_mother_pdgId );
    smalltree->Branch("tree_track_sim_isFromBC_LLP",          &tree_track_sim_isFromBC_LLP );
    smalltree->Branch("tree_track_sim_dV",                    &tree_track_sim_dV );
    smalltree->Branch("tree_track_sim_dist",                  &tree_track_sim_dist );
//$$
    
    // Tracking particle info
    
    // smalltree->Branch("tree_simtrack_charge",    &tree_simtrack_charge );
    // smalltree->Branch("tree_simtrack_pt",    &tree_simtrack_pt );
    // smalltree->Branch("tree_simtrack_eta",    &tree_simtrack_eta  );
    // smalltree->Branch("tree_simtrack_phi",    &tree_simtrack_phi  );
    // smalltree->Branch("tree_simtrack_longLived",    &tree_simtrack_longLived );
    // smalltree->Branch("tree_simtrack_pdgId",    &tree_simtrack_pdgId );
    // smalltree->Branch("tree_simtrack_mass",    &tree_simtrack_mass   );
    // smalltree->Branch("tree_simtrack_status",    &tree_simtrack_status );
    // smalltree->Branch("tree_simtrack_vx",     &tree_simtrack_vx );
    // smalltree->Branch("tree_simtrack_vy",     &tree_simtrack_vy );
    // smalltree->Branch("tree_simtrack_vz",     &tree_simtrack_vz );
    // smalltree->Branch("tree_simtrack_isRecoMatched",      &tree_simtrack_isRecoMatched  );
    // smalltree->Branch("tree_simtrack_pca_dxy",            &tree_simtrack_pca_dxy);
    // smalltree->Branch("tree_simtrack_pca_dz",             &tree_simtrack_pca_dz);
    // smalltree->Branch("tree_simtrack_trkIdx",             &tree_simtrack_trkIdx);
    // smalltree->Branch("tree_simtrack_isRecoMatched_pt",   &tree_simtrack_isRecoMatched_pt);/*!*/
    
    smalltree->Branch("tree_AK4Slimmedjet_E"  ,       &tree_AK4Slimmedjet_E);
    smalltree->Branch("tree_AK4Slimmedjet_pt"  ,      &tree_AK4Slimmedjet_pt);
    smalltree->Branch("tree_AK4Slimmedjet_eta" ,      &tree_AK4Slimmedjet_eta);
    smalltree->Branch("tree_AK4Slimmedjet_phi" ,      &tree_AK4Slimmedjet_phi);
//$$    smalltree->Branch("tree_AK4Slimmedjet_idxTrack" , &tree_AK4Slimmedjet_idxTrack);
    
    smalltree->Branch("tree_AK4PFjet_E"  ,       &tree_AK4PFjet_E);
    smalltree->Branch("tree_AK4PFjet_pt"  ,      &tree_AK4PFjet_pt);
    smalltree->Branch("tree_AK4PFjet_eta" ,      &tree_AK4PFjet_eta);
    smalltree->Branch("tree_AK4PFjet_phi" ,      &tree_AK4PFjet_phi);
//$$    smalltree->Branch("tree_AK4PFjet_idxTrack" , &tree_AK4PFjet_idxTrack);
    
    smalltree->Branch("tree_CaloJet_E"  ,      &tree_CaloJet_E);
    smalltree->Branch("tree_CaloJet_pt"  ,     &tree_CaloJet_pt);
    smalltree->Branch("tree_CaloJet_eta" ,     &tree_CaloJet_eta);
    smalltree->Branch("tree_CaloJet_phi" ,     &tree_CaloJet_phi);
//$$    smalltree->Branch("tree_CaloJet_idxTrack", &tree_CaloJet_idxTrack);
    
    smalltree->Branch("tree_jet08_E"  ,      &tree_jet08_E);
    smalltree->Branch("tree_jet08_pt"  ,     &tree_jet08_pt);
    smalltree->Branch("tree_jet08_eta" ,     &tree_jet08_eta);
    smalltree->Branch("tree_jet08_phi" ,     &tree_jet08_phi);
//$$    smalltree->Branch("tree_jet08_idxTrack", &tree_jet08_idxTrack);
    
    // met info
    smalltree->Branch("tree_PFMet_et" ,  &tree_PFMet_et);
    smalltree->Branch("tree_PFMet_phi" , &tree_PFMet_phi);
    smalltree->Branch("tree_PFMet_sig" , &tree_PFMet_sig);
    
    // gen info
    smalltree->Branch("tree_GenPVx" ,  &tree_GenPVx);
    smalltree->Branch("tree_GenPVy" ,  &tree_GenPVy);
    smalltree->Branch("tree_GenPVz" ,  &tree_GenPVz);
    
    smalltree->Branch("tree_genParticle_pt"  ,          &tree_genParticle_pt);
    smalltree->Branch("tree_genParticle_eta" ,          &tree_genParticle_eta);
    smalltree->Branch("tree_genParticle_phi" ,          &tree_genParticle_phi);
    smalltree->Branch("tree_genParticle_charge" ,       &tree_genParticle_charge);
    smalltree->Branch("tree_genParticle_pdgId" ,        &tree_genParticle_pdgId);
    smalltree->Branch("tree_genParticle_vx"  ,          &tree_genParticle_vx);
    smalltree->Branch("tree_genParticle_vy" ,           &tree_genParticle_vy);
    smalltree->Branch("tree_genParticle_vz" ,           &tree_genParticle_vz);
    smalltree->Branch("tree_genParticle_mass" ,         &tree_genParticle_mass);
    smalltree->Branch("tree_genParticle_statusCode",    &tree_genParticle_statusCode);
    smalltree->Branch("tree_genParticle_mother_pdgId" , &tree_genParticle_mother_pdgId);

//$$
    smalltree->Branch("nFromC",  &nFromC,  "nFromC/I");
    smalltree->Branch("tree_genFromC_pt"  ,          &tree_genFromC_pt);
    smalltree->Branch("tree_genFromC_eta" ,          &tree_genFromC_eta);
    smalltree->Branch("tree_genFromC_phi" ,          &tree_genFromC_phi);
    smalltree->Branch("tree_genFromC_charge" ,       &tree_genFromC_charge);
    smalltree->Branch("tree_genFromC_pdgId" ,        &tree_genFromC_pdgId);
    smalltree->Branch("tree_genFromC_vx"  ,          &tree_genFromC_vx);
    smalltree->Branch("tree_genFromC_vy" ,           &tree_genFromC_vy);
    smalltree->Branch("tree_genFromC_vz" ,           &tree_genFromC_vz);
    smalltree->Branch("tree_genFromC_mother_pdgId" , &tree_genFromC_mother_pdgId);
    smalltree->Branch("tree_genFromC_generation" ,   &tree_genFromC_generation);
    smalltree->Branch("tree_genFromC_LLP" ,          &tree_genFromC_LLP);

    smalltree->Branch("nFromB",  &nFromB,  "nFromB/I");
    smalltree->Branch("tree_genFromB_pt"  ,	     &tree_genFromB_pt);
    smalltree->Branch("tree_genFromB_eta" ,	     &tree_genFromB_eta);
    smalltree->Branch("tree_genFromB_phi" ,	     &tree_genFromB_phi);
    smalltree->Branch("tree_genFromB_charge" ,	     &tree_genFromB_charge);
    smalltree->Branch("tree_genFromB_pdgId" ,	     &tree_genFromB_pdgId);
    smalltree->Branch("tree_genFromB_vx"  ,	     &tree_genFromB_vx);
    smalltree->Branch("tree_genFromB_vy" ,	     &tree_genFromB_vy);
    smalltree->Branch("tree_genFromB_vz" ,	     &tree_genFromB_vz);
    smalltree->Branch("tree_genFromB_mother_pdgId" , &tree_genFromB_mother_pdgId);
    smalltree->Branch("tree_genFromB_generation" ,   &tree_genFromB_generation);
    smalltree->Branch("tree_genFromB_LLP" ,          &tree_genFromB_LLP);
//$$
    
    // genJet info
    smalltree->Branch("tree_genJet_pt"  ,   &tree_genJet_pt);
    smalltree->Branch("tree_genJet_eta" ,   &tree_genJet_eta);
    smalltree->Branch("tree_genJet_phi" ,   &tree_genJet_phi);
    smalltree->Branch("tree_genJet_mass",   &tree_genJet_mass);
    smalltree->Branch("tree_genJet_energy", &tree_genJet_energy);
    
    smalltree->Branch("tree_ak8GenJet_pt"  ,   &tree_ak8GenJet_pt);
    smalltree->Branch("tree_ak8GenJet_eta" ,   &tree_ak8GenJet_eta);
    smalltree->Branch("tree_ak8GenJet_phi" ,   &tree_ak8GenJet_phi);
    smalltree->Branch("tree_ak8GenJet_mass",   &tree_ak8GenJet_mass);
    smalltree->Branch("tree_ak8GenJet_energy", &tree_ak8GenJet_energy);
    
    // electrons info
    smalltree->Branch("tree_electron_pt"  ,   &tree_electron_pt);
    smalltree->Branch("tree_electron_eta" ,   &tree_electron_eta);
    smalltree->Branch("tree_electron_phi" ,   &tree_electron_phi);
    smalltree->Branch("tree_electron_vx"  ,   &tree_electron_vx);
    smalltree->Branch("tree_electron_vy" ,    &tree_electron_vy);
    smalltree->Branch("tree_electron_vz" ,    &tree_electron_vz);
    smalltree->Branch("tree_electron_energy", &tree_electron_energy);
    smalltree->Branch("tree_electron_charge", &tree_electron_charge);
    
    // muons info
    smalltree->Branch("tree_slimmedmuon_pt"  ,                   &tree_slimmedmuon_pt);
    smalltree->Branch("tree_slimmedmuon_eta" ,                   &tree_slimmedmuon_eta);
    smalltree->Branch("tree_slimmedmuon_phi" ,                   &tree_slimmedmuon_phi);
    smalltree->Branch("tree_slimmedmuon_vx"  ,                   &tree_slimmedmuon_vx);
    smalltree->Branch("tree_slimmedmuon_vy" ,                    &tree_slimmedmuon_vy);
    smalltree->Branch("tree_slimmedmuon_vz" ,                    &tree_slimmedmuon_vz);
    smalltree->Branch("tree_slimmedmuon_energy",                 &tree_slimmedmuon_energy);
    smalltree->Branch("tree_slimmedmuon_dxy",                    &tree_slimmedmuon_dxy);
    smalltree->Branch("tree_slimmedmuon_dxyError",               &tree_slimmedmuon_dxyError);
    smalltree->Branch("tree_slimmedmuon_dz",                     &tree_slimmedmuon_dz);
    smalltree->Branch("tree_slimmedmuon_dzError",                &tree_slimmedmuon_dzError);
    smalltree->Branch("tree_slimmedmuon_charge",                 &tree_slimmedmuon_charge);
    smalltree->Branch("tree_slimmedmuon_PFisoVeryTight",         &tree_slimmedmuon_PFisoVeryTight);
    smalltree->Branch("tree_slimmedmuon_PFisoTight",             &tree_slimmedmuon_PFisoTight);
    smalltree->Branch("tree_slimmedmuon_PFisoMedium",            &tree_slimmedmuon_PFisoMedium);
    smalltree->Branch("tree_slimmedmuon_PFisoLoose",             &tree_slimmedmuon_PFisoLoose);
    smalltree->Branch("tree_slimmedmuon_MVAisoLoose",            &tree_slimmedmuon_MVAisoLoose);
    smalltree->Branch("tree_slimmedmuon_MVAisoMedium",           &tree_slimmedmuon_MVAisoMedium);
    smalltree->Branch("tree_slimmedmuon_MVAisoTight",            &tree_slimmedmuon_MVAisoTight);
    smalltree->Branch("tree_slimmedmuon_isGlobalMuon",           &tree_slimmedmuon_isGlobalMuon);
    smalltree->Branch("tree_slimmedmuon_isStandAloneMuon",       &tree_slimmedmuon_isStandAloneMuon);
    smalltree->Branch("tree_slimmedmuon_CutBasedIdLoose",        &tree_slimmedmuon_CutBasedIdLoose        );
    smalltree->Branch("tree_slimmedmuon_CutBasedIdMedium",       &tree_slimmedmuon_CutBasedIdMedium       );
    smalltree->Branch("tree_slimmedmuon_CutBasedIdMediumPrompt", &tree_slimmedmuon_CutBasedIdMediumPrompt );
    smalltree->Branch("tree_slimmedmuon_CutBasedIdTight",        &tree_slimmedmuon_CutBasedIdTight        );
    
//$$
    smalltree->Branch("tree_nLLP",&tree_nLLP);
    smalltree->Branch("tree_LLP",&tree_LLP);
    smalltree->Branch("tree_LLP_pt" ,&tree_LLP_pt);
    smalltree->Branch("tree_LLP_eta",&tree_LLP_eta);
    smalltree->Branch("tree_LLP_phi",&tree_LLP_phi);
    smalltree->Branch("tree_LLP_x",&tree_LLP_x);
    smalltree->Branch("tree_LLP_y",&tree_LLP_y);
    smalltree->Branch("tree_LLP_z",&tree_LLP_z);
    smalltree->Branch("tree_LLP_nTrks",&tree_LLP_nTrks);
    smalltree->Branch("tree_LLP_Vtx_nTrks",&tree_LLP_Vtx_nTrks);
    smalltree->Branch("tree_LLP_Vtx_NChi2",&tree_LLP_Vtx_NChi2);
    smalltree->Branch("tree_LLP_Vtx_dx",&tree_LLP_Vtx_dx);
    smalltree->Branch("tree_LLP_Vtx_dy",&tree_LLP_Vtx_dy);
    smalltree->Branch("tree_LLP_Vtx_dz",&tree_LLP_Vtx_dz);

    smalltree->Branch("tree_Hemi",       &tree_Hemi);
    smalltree->Branch("tree_Hemi_njet",  &tree_Hemi_njet);
    smalltree->Branch("tree_Hemi_eta",   &tree_Hemi_eta);
    smalltree->Branch("tree_Hemi_phi",   &tree_Hemi_phi);
    smalltree->Branch("tree_Hemi_dR",    &tree_Hemi_dR);
    smalltree->Branch("tree_Hemi_nTrks", &tree_Hemi_nTrks);
    smalltree->Branch("tree_Hemi_nTrks_sig", &tree_Hemi_nTrks_sig);
    smalltree->Branch("tree_Hemi_nTrks_bad", &tree_Hemi_nTrks_bad);
    smalltree->Branch("tree_Hemi_LLP",       &tree_Hemi_LLP);
    smalltree->Branch("tree_Hemi_LLP_pt",    &tree_Hemi_LLP_pt);
    smalltree->Branch("tree_Hemi_LLP_eta",   &tree_Hemi_LLP_eta);
    smalltree->Branch("tree_Hemi_LLP_phi",   &tree_Hemi_LLP_phi);
    smalltree->Branch("tree_Hemi_LLP_dist",  &tree_Hemi_LLP_dist);
    smalltree->Branch("tree_Hemi_LLP_x",     &tree_Hemi_LLP_x);
    smalltree->Branch("tree_Hemi_LLP_y",     &tree_Hemi_LLP_y);
    smalltree->Branch("tree_Hemi_LLP_z",     &tree_Hemi_LLP_z);
    smalltree->Branch("tree_Hemi_Vtx_NChi2", &tree_Hemi_Vtx_NChi2);
    smalltree->Branch("tree_Hemi_Vtx_nTrks", &tree_Hemi_Vtx_nTrks);
    smalltree->Branch("tree_Hemi_Vtx_nTrks_sig", &tree_Hemi_Vtx_nTrks_sig);
    smalltree->Branch("tree_Hemi_Vtx_nTrks_bad", &tree_Hemi_Vtx_nTrks_bad);
    smalltree->Branch("tree_Hemi_Vtx_x",     &tree_Hemi_Vtx_x);
    smalltree->Branch("tree_Hemi_Vtx_y",     &tree_Hemi_Vtx_y);
    smalltree->Branch("tree_Hemi_Vtx_z",     &tree_Hemi_Vtx_z);
    smalltree->Branch("tree_Hemi_Vtx_dx",    &tree_Hemi_Vtx_dx);
    smalltree->Branch("tree_Hemi_Vtx_dy",    &tree_Hemi_Vtx_dy);
    smalltree->Branch("tree_Hemi_Vtx_dz",    &tree_Hemi_Vtx_dz);
    smalltree->Branch("tree_Hemi_dR12",      &tree_Hemi_dR12);
    smalltree->Branch("tree_Hemi_LLP_dR12",  &tree_Hemi_LLP_dR12);
//$$

    tree_NbrOfZCand= 0;
    
    runNumber = 0;
    eventNumber = 0;
    lumiBlock = 0;
    
    tree_bs_PosX= 0;
    tree_bs_PosY= 0;
    tree_bs_PosZ= 0;
    
    const bool tpRef = iConfig.getUntrackedParameter<bool>("trackingParticlesRef");
    const auto tpTag = iConfig.getUntrackedParameter<edm::InputTag>("trackingParticles");
    if(tpRef && !runOnData_) {
        trackingParticleRefToken_ = consumes<TrackingParticleRefVector>(tpTag);
    }
    else {
        trackingParticleToken_ = consumes<TrackingParticleCollection>(tpTag);
    }
}


TrackingPerf::~TrackingPerf()
{
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackingPerf::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  clearVariables();
//$$
  bool showlog = true;
//$$
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  runNumber = iEvent.id().run();
  //std::cout << "runNumber = " << runNumber << std::endl;
  eventNumber = iEvent.id().event();
  std::cout << "eventNumber = "<< eventNumber <<std::endl;
  lumiBlock = iEvent.luminosityBlock();
  

  //// HANDLES /////
  // Triggers
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_,triggerBits);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  
  if ( nEvent == 0 ) getHLTPathNames(names);
  fillTriggerBits(HLTPathsByName_, triggerBits, names);
  
  // edm::Handle<edm::View<reco::Track> > TracksForRes;
  edm::Handle<edm::View<pat::PackedCandidate> > TracksForRes;
  iEvent.getByToken(trackSrc_, TracksForRes);

  // edm::Handle<edm::View<reco::Track> > tracksHandle;
  edm::Handle<edm::View<pat::PackedCandidate> > tracksHandle;
  iEvent.getByToken(trackToken_, tracksHandle);
  // const edm::View<reco::Track>& tracks = *tracksHandle;
  const edm::View<pat::PackedCandidate>& tracks = *tracksHandle;

  edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
  iSetup.get<TransientRecHitRecord>().get(ttrhbuilder_,theTrackerRecHitBuilder);
  
  Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(beamSpotToken_, recoBeamSpotHandle);
  
  // edm::Handle<reco::VertexCollection> vertices;
  edm::Handle<vector<reco::Vertex> > vertices;
  iEvent.getByToken(vertexToken_, vertices);
  
  // edm::Handle<edm::View<reco::Jet> > ak4slimmedJets;
  // edm::Handle<edm::View<reco::Jet> > ak4PFJets;
  // edm::Handle<edm::View<reco::Jet> > CaloJets;
  // edm::Handle<edm::View<reco::Jet> > ak8jets;
  
  // iEvent.getByToken(ak4slimmedJetToken_,ak4slimmedJets);
  // iEvent.getByToken(ak4PFJetToken_,ak4PFJets);
  // iEvent.getByToken(CaloJetToken_,CaloJets);
  // iEvent.getByToken(ak8jetToken_,ak8jets);
  edm::Handle<edm::View<pat::Jet> > ak4slimmedJets;
  iEvent.getByToken(ak4slimmedJetToken_,ak4slimmedJets);

  edm::Handle<edm::View<pat::Jet> > ak8jets;
  iEvent.getByToken(ak8jetToken_,ak8jets);

  edm::Handle<pat::METCollection> PFMETs;
  iEvent.getByToken(metToken_, PFMETs);
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  if(!runOnData_) iEvent.getByToken(genParticlesToken_,genParticles);
  
  edm::Handle<edm::View<reco::GenJet>>genJets;
  if(!runOnData_) iEvent.getByToken(genJetToken_,genJets);
  
  edm::Handle<reco::GenJetCollection> ak8GenJets;
  if(!runOnData_) iEvent.getByToken(ak8GenJetToken_,ak8GenJets);
  
  edm::Handle<GenEventInfoProduct> genEventInfo;
  if(!runOnData_) iEvent.getByToken(genEventInfoToken_,genEventInfo);
  
  edm::Handle<LHEEventProduct> lheEventProduct;
  if(!runOnData_) iEvent.getByToken(LHEEventProductToken_,lheEventProduct);
  
  edm::Handle<pat::PackedCandidate> pfcands;
  iEvent.getByToken(pfcandsToken_,pfcands);
  
  edm::Handle<pat::ElectronCollection> electronsPAT;
  iEvent.getByToken(electronPATToken_,electronsPAT);
  
  edm::Handle<pat::MuonCollection> slimmedmuons;
  iEvent.getByToken(slimmedmuonToken_,slimmedmuons);
  
    //Not compatible with miniaod
  
  // edm::Handle<vector<reco::CompositeCandidate> > dimuons;//RECO

  // if ( !runOnData_ )
  // {
  //     iEvent.getByToken(ZmumuCandidate_, dimuons);
  //     const vector<reco::CompositeCandidate> theDimuonToChecks = *dimuons;//RECO
      
  //         std::cout << "found " << theDimuonToChecks.size() << " z candidates " << endl;
  //     tree_NbrOfZCand = theDimuonToChecks.size();
  // }

  //geometry & magnetif field
  edm::ESHandle<TrackerGeometry> pDD;
  iSetup.get<TrackerDigiGeometryRecord>().get(pDD);
  
  iSetup.get<IdealMagneticFieldRecord>().get(bField);
  
  //	edm::ESHandle<TrackerTopology> tTopoHandle;
  //	iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  //	const TrackerTopology* const tTopo = tTopoHandle.product();
  
  //association between tracks
  edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
  
  if(!runOnData_) iEvent.getByToken(trackAssociatorToken_, theAssociator);
  // const reco::TrackToTrackingParticleAssociator& associatorByHits = *theAssociator;
  
  //tracking Particles info and matching to reco
  
  edm::ESHandle<ParametersDefinerForTP> parametersDefinerH;
  if(!runOnData_) {
      iSetup.get<TrackAssociatorRecord>().get(parametersDefinerName_, parametersDefinerH);
  }
  
  // const ParametersDefinerForTP *parametersDefiner = parametersDefinerH.product();
  
  TrackingParticleRefVector tmpTP;
  const TrackingParticleRefVector *tmpTPptr = nullptr;
  edm::Handle<TrackingParticleCollection>  TPCollectionH;
  edm::Handle<TrackingParticleRefVector>  TPCollectionHRefVector;
  if ( !runOnData_ ) {
      if( !trackingParticleToken_.isUninitialized() ) {
  	  iEvent.getByToken(trackingParticleToken_, TPCollectionH);
  	  for(size_t i=0, size=TPCollectionH->size(); i<size; ++i) {
  	      tmpTP.push_back(TrackingParticleRef(TPCollectionH, i));
  	  }
  	  tmpTPptr = &tmpTP;
      }
      else {
  	  iEvent.getByToken(trackingParticleRefToken_, TPCollectionHRefVector);
  	  tmpTPptr = TPCollectionHRefVector.product();
      }
  }
  
  const TrackingParticleRefVector& tpCollection = *tmpTPptr;
  TrackingParticleRefKeyToIndex tpKeyToIndex;
  
  if ( !runOnData_ ) {
      for (size_t i=0; i<tpCollection.size(); ++i) {
  	  tpKeyToIndex[tpCollection[i].key()] = i;
      }
  }
  
  //////////////////////////////////
  //////////////////////////////////
  //////////    BS     /////////////
  //////////////////////////////////
  //////////////////////////////////
  
  BeamSpot const & bs = *recoBeamSpotHandle;
  tree_bs_PosX = bs.x0();
  tree_bs_PosY = bs.y0();
  tree_bs_PosZ = bs.z0();
  
  //////////////////////////////////
  //////////////////////////////////
  //////   Primary Vertices  ///////
  //////////////////////////////////
  //////////////////////////////////
  
  int nPV = vertices->size();
  
  for (auto const & vertex : *vertices) {
      tree_vtx_PosX.push_back(vertex.x());
      tree_vtx_PosY.push_back(vertex.y());
      tree_vtx_PosZ.push_back(vertex.z());
      tree_vtx_NChi2.push_back(vertex.normalizedChi2());
      tree_vtx_PosXError.push_back(vertex.xError());
      tree_vtx_PosYError.push_back(vertex.yError());
      tree_vtx_PosZError.push_back(vertex.zError());
      nPV++;
  }
  float PVx = tree_vtx_PosX[0];//l'index 0 donne le PV!
  float PVy = tree_vtx_PosY[0];
  float PVz = tree_vtx_PosZ[0];

  //////////////////////////////////
  //////////////////////////////////
  ////////   MET   /////////////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_PFMet_et  = -10.;
  tree_PFMet_phi = -10.;
  tree_PFMet_sig = -10.;
  if ( PFMETs->size() > 0 ) {
      const pat::MET &themet = PFMETs->front();
      tree_PFMet_et  = themet.et();
      tree_PFMet_phi = themet.phi();
      tree_PFMet_sig = themet.significance();
  }
  
  //////////////////////////////////
  //////////////////////////////////
  ///////////	Jets   /////////////
  //////////////////////////////////
  //////////////////////////////////
  
  //cout << "number of jets "<<nJet<<endl;
  //for ak4
  for (int ij=0;ij<int(ak4slimmedJets->size());ij++)
  {
      const Jet& jet = ak4slimmedJets->at(ij);
      if ( jet.pt() < 20. ) continue;
      tree_AK4Slimmedjet_E.push_back(jet.energy());
      tree_AK4Slimmedjet_pt.push_back(jet.pt());
      tree_AK4Slimmedjet_eta.push_back(jet.eta());
      tree_AK4Slimmedjet_phi.push_back(jet.phi());
//$$
//       TLorentzVector thejet, track;
//       thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
//       for (unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++)
//       {
//   	  track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
//   	  if(thejet.DeltaR(track) < 0.4 ) tree_AK4Slimmedjet_idxTrack.push_back(iTrack);
//       }
//$$
  }
  
  // for (int ij=0;ij<int(ak4PFJets->size());ij++)
  // {
  //     const Jet& jet = ak4PFJets->at(ij);
  //     if ( jet.pt() < 20. ) continue;
  //     tree_AK4PFjet_E.push_back(jet.energy());
  //     tree_AK4PFjet_pt.push_back(jet.pt());
  //     tree_AK4PFjet_eta.push_back(jet.eta());
  //     tree_AK4PFjet_phi.push_back(jet.phi());
//$$
//       TLorentzVector thejet, track;
//       thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
//       for (unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++)
//       {
//   	  track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
//   	  if(thejet.DeltaR(track) < 0.4 ) tree_AK4PFjet_idxTrack.push_back(iTrack);
//       }
//$$
  // }
  
  // for (int ij=0;ij< int(CaloJets->size());ij++)//plus simpliste, utilisete just eles cellules calorimetriques
  // {
  //     const Jet& jet = CaloJets->at(ij);
  //     if ( jet.pt() < 20. ) continue;
  //     tree_CaloJet_E.push_back(jet.energy());
  //     tree_CaloJet_pt.push_back(jet.pt());
  //     tree_CaloJet_eta.push_back(jet.eta());
  //     tree_CaloJet_phi.push_back(jet.phi());
//$$
//       TLorentzVector thejet, track;
//       thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
//       for (unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) {
//   	  track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
//   	  if(thejet.DeltaR(track) < 0.4 ) tree_CaloJet_idxTrack.push_back(iTrack);
//       }
//$$
  // }
  
  for (int ij=0;ij< int(ak8jets->size() );ij++)
  {
      const Jet& jet = ak8jets->at(ij);
      if ( jet.pt() < 20. ) continue;
      tree_jet08_E.push_back(jet.energy());
      tree_jet08_pt.push_back(jet.pt());
      tree_jet08_eta.push_back(jet.eta());
      tree_jet08_phi.push_back(jet.phi());
//$$
//       TLorentzVector thejet, track;
//       thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
//       for (unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++)
//       {
//   	  track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
//   	  if(thejet.DeltaR(track) < 0.8 ) tree_jet08_idxTrack.push_back(iTrack);
//       }
//$$
  }

  //////////////////////////////////
  //////////////////////////////////
  ////////   Electrons   ///////////
  //////////////////////////////////
  //////////////////////////////////
  
  //int nElectrons = electronsPAT->size();
  
  for (auto const & electron : *electronsPAT)
  {
      tree_electron_pt.push_back(electron.pt());
      tree_electron_eta.push_back(electron.eta());
      tree_electron_phi.push_back(electron.phi());
      tree_electron_vx.push_back(electron.vx());
      tree_electron_vy.push_back(electron.vy());
      tree_electron_vz.push_back(electron.vz());
      tree_electron_energy.push_back(electron.energy());
      tree_electron_charge.push_back(electron.charge());
  }
  
  //////////////////////////////////
  //////////////////////////////////
  ////////////   Muons   ///////////
  //////////////////////////////////
  //////////////////////////////////
  
  int muon_idx=0;
  // bool ZMuRec = false;
  // int PFiso, MVAiso, CutBasedId;
  float dVr, dVz;
  float Mmumu = 0.;
  int nmurec = 0, imu1 = -1, imu2 = -1;
  float mupt1, mueta1, muphi1, mupt2, mueta2, muphi2;
  float mu_mass = 0.1057;
  TLorentzVector v1, v2, v;
  
  for (auto const & muon : *slimmedmuons)
  {
      if ( muon.pt() < 3. ) continue;
      //std::cout<< "midx: "<<muon_idx<<std::endl;
      muon_idx++;
      int nmu = slimmedmuons->size();
      tree_slimmedmuon_pt.push_back( muon.pt());
      tree_slimmedmuon_eta.push_back(muon.eta());
      tree_slimmedmuon_phi.push_back(muon.phi());
      tree_slimmedmuon_vx.push_back( muon.vx());
      tree_slimmedmuon_vy.push_back( muon.vy());
      tree_slimmedmuon_vz.push_back( muon.vz());
      tree_slimmedmuon_energy.push_back(   muon.energy());
      tree_slimmedmuon_dxy.push_back(	  muon.bestTrack()->dxy(bs));
      tree_slimmedmuon_dxyError.push_back( muon.bestTrack()->dxyError());
      tree_slimmedmuon_dz.push_back(	 muon.bestTrack()->dz(bs.position()));
      tree_slimmedmuon_dzError.push_back(  muon.bestTrack()->dzError());
      tree_slimmedmuon_charge.push_back(   muon.charge());
      
      tree_slimmedmuon_PFisoVeryTight.  push_back(   muon.passed(reco::Muon::PFIsoVeryTight));
      tree_slimmedmuon_PFisoTight.    push_back(   muon.passed(reco::Muon::PFIsoTight) );
      tree_slimmedmuon_PFisoMedium.    push_back(   muon.passed(reco::Muon::PFIsoMedium));
      tree_slimmedmuon_PFisoLoose.    push_back(   muon.passed(reco::Muon::PFIsoLoose ));
      tree_slimmedmuon_MVAisoLoose.    push_back(   muon.passed(reco::Muon::MvaLoose )  );
      tree_slimmedmuon_MVAisoMedium.	push_back(   muon.passed(reco::Muon::MvaMedium)  );
      tree_slimmedmuon_MVAisoTight.    push_back(   muon.passed(reco::Muon::MvaTight )  );
      tree_slimmedmuon_isGlobalMuon.	push_back(   muon.passed(muon.isGlobalMuon())  );
      tree_slimmedmuon_isStandAloneMuon.push_back(   muon.passed(muon.isStandAloneMuon()));
      tree_slimmedmuon_CutBasedIdLoose        .push_back(    muon.passed(reco::Muon::CutBasedIdLoose));
      tree_slimmedmuon_CutBasedIdMedium     .push_back(   muon.passed(reco::Muon::CutBasedIdMedium));
      tree_slimmedmuon_CutBasedIdMediumPrompt	  .push_back(	muon.passed(reco::Muon::CutBasedIdMediumPrompt));
      tree_slimmedmuon_CutBasedIdTight        .push_back(    muon.passed(reco::Muon::CutBasedIdTight));
      
      if ( !muon.passed(muon.isGlobalMuon()) ) continue; //Need global muons
      mupt1  = muon.pt();
      if ( mupt1 < 10. ) continue; //Zmu filter
      dVr = TMath::Sqrt( (muon.vx()-PVx)*(muon.vx()-PVx) + (muon.vy()-PVy)*(muon.vy()-PVy) );
      dVz = muon.vz()-PVz;
      if ( dVr > 0.1 || abs(dVz) > 0.2 ) continue;// on veut un bon fit pour nos PV d'ou un seuil maximum sur les distances
      mueta1 = muon.eta();
      muphi1 = muon.phi();
      
      v1.SetPtEtaPhiM(mupt1,mueta1,muphi1,mu_mass);
      nmurec++;
      //std::cout<<"nmu : "<<nmu<<std::endl;
      if ( muon_idx == nmu-1 ) continue;
      for ( int muon2_idx=muon_idx+1;muon2_idx<nmu;muon2_idx++) //recompte des éléments déjà comptés
      {       // Loop on other reco muons
      //std::cout<<"muon_idx:  "<< muon_idx<< "//"<<"muon2_idx : "<<muon2_idx<<std::endl;
  	  if ( !tree_slimmedmuon_isGlobalMuon[muon2_idx] ) continue;
  	  if ( tree_slimmedmuon_charge[muon_idx] == tree_slimmedmuon_charge[muon2_idx] ) continue;
  	  dVr = TMath::Sqrt( (tree_slimmedmuon_vx[muon2_idx]-PVx)*(tree_slimmedmuon_vx[muon2_idx]-PVx) + (tree_slimmedmuon_vy[muon2_idx]-PVy)*(tree_slimmedmuon_vy[muon2_idx]-PVy) );
  	  dVz = tree_slimmedmuon_vz[muon2_idx]-PVz;
  	  if ( dVr > 0.1 || abs(dVz) > 0.2 ) continue;
  	  mupt2  = tree_slimmedmuon_pt[muon2_idx];
  	  if ( mupt2 < 10. ) continue;
  	  if ( mupt1 < 28. && mupt2 < 28. ) continue; //Zmu FIlter
  	  mueta2 = tree_slimmedmuon_eta[muon2_idx];
  	  muphi2 = tree_slimmedmuon_phi[muon2_idx];
  	  muon2_idx++;
  	  v2.SetPtEtaPhiM(mupt2,mueta2,muphi2,mu_mass);
  	  v = v1 + v2;
  	  if ( v.Mag() > Mmumu )
  	  { //Mag pour masse invariante (magnitude)
  	      Mmumu = v.Mag();
  	      imu1 = muon_idx;
  	      imu2 = muon2_idx;
  	  }
      }
          
      // if ( Mmumu < 60. ) continue;
      // std::cout<<"imu1 "<<imu1<<std::endl;
      // std::cout<<"imu2 "<< imu2<< std::endl;
  }
  //  std::cout<< "number of muons first method : " << tree_slimmedmuon_pt.size() << std::endl;

  if ( tree_slimmedmuon_pt[imu2] > tree_slimmedmuon_pt[imu1] ) {
      int imu0 = imu2;
      imu2 = imu1; //muons reco with imu1 having the highest pt
      imu1 = imu0;
  }
  

  //////////////////////////////////
  //////////////////////////////////
  //////// gen Information /////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_nLLP = -1;
  tree_GenPVx = -1.;
  tree_GenPVy = -1.;
  tree_GenPVz = -20.;

  int nLLP = 0;
  int nllp = 0;
  nBC = 0; 
      
  // Gen Information  for event axis //
  float  Gen_neu1_eta=-10, Gen_neu1_phi=-10;
  float  Gen_neu2_eta=-10, Gen_neu2_phi=-10;
  //  float deltaPhi;
  int  neu[2],  nneu = 0;
  TLorentzVector vneu[2];
  
  float dRneuneu = 0.;
  
  for (int k=0; k<2; k++) {
    neu[k] = -1;
  }
  /////////////////////////////////////////////////////////
  
  if ( !runOnData_ ) {
    //      int nGenParticles = genParticles->size();
    //cout << "number of gen particles "<<nGenParticles<<endl;
    
    /// GEN PARTICLES

    int genParticle_idx=0;
    for (auto const & genParticle : *genParticles) 
    {
      genParticle_idx++;
      int pdgid     = genParticle.pdgId();
      float Gen_pt  = genParticle.pt();
      float Gen_eta = genParticle.eta();
      float Gen_phi = genParticle.phi();
      float Gen_m   = genParticle.mass();
      // smuon
      if ( pdgid == 1000013 ) {
	tree_GenPVx = genParticle.vx();
	tree_GenPVy = genParticle.vy();
	tree_GenPVz = genParticle.vz();
      }
      const Candidate * mom = genParticle.mother();
      
      // neutralino from smuon
      if ( pdgid == 1000023 && abs(mom->pdgId()) == 1000013 ) {
	nLLP++;
	if ( nLLP == 1 ) {
	  LLP1_pt  = Gen_pt;
	  LLP1_eta = Gen_eta;
	  LLP1_phi = Gen_phi;
	}
	if ( nLLP == 2 ) {
	  LLP2_pt  = Gen_pt;
	  LLP2_eta = Gen_eta;
	  LLP2_phi = Gen_phi;
	}
	if ( neu[0] < 0 ) {
	  neu[0] = genParticle_idx;
	  vneu[0].SetPtEtaPhiM( Gen_pt, Gen_eta, Gen_phi, Gen_m );
	  Gen_neu1_eta = Gen_eta;
	  Gen_neu1_phi = Gen_phi;
	}
	else if ( neu[1] < 0 ) {
	  neu[1] = genParticle_idx;
	  vneu[1].SetPtEtaPhiM( Gen_pt, Gen_eta, Gen_phi, Gen_m );
	  Gen_neu2_eta = Gen_eta;
	  Gen_neu2_phi = Gen_phi;
	}
	nneu++;
      }
      
      if ( nneu == 2 ) {
	dRneuneu = Deltar( Gen_neu1_eta, Gen_neu1_phi, Gen_neu2_eta, Gen_neu2_phi );
	// deltaPhi = Deltaphi( Gen_neu1_phi, Gen_neu2_phi );
      }
      
      // quarks from neutralino
      if ( abs(pdgid) >= 1 && abs(pdgid) <= 6 && abs(mom->pdgId()) == 1000023 ) {
	if ( nllp >= 2 ) {
	  float dV1 = (genParticle.vx() - LLP1_x)*(genParticle.vx() - LLP1_x)
	            + (genParticle.vy() - LLP1_y)*(genParticle.vy() - LLP1_y)
	            + (genParticle.vz() - LLP1_z)*(genParticle.vz() - LLP1_z); // dV1 is equal to dV from nllp==1
	  float dV2 = (genParticle.vx() - LLP2_x)*(genParticle.vx() - LLP2_x)
	            + (genParticle.vy() - LLP2_y)*(genParticle.vy() - LLP2_y)
	            + (genParticle.vz() - LLP2_z)*(genParticle.vz() - LLP2_z);
//             std::cout<<"dV1 for nllp>=2 : "<<dV1<<"  dV2 for nllp>=2 : "<<dV2<<std::endl;
	  if ( dV1 > 0.01 && dV2 > 0.01 ) nllp++; // should be == 2, so just to check : dV2 is always equal to 0 here
	}
	if ( nllp == 1 ) {
	  float dV = (genParticle.vx() - LLP1_x)*(genParticle.vx() - LLP1_x)
	           + (genParticle.vy() - LLP1_y)*(genParticle.vy() - LLP1_y)
	           + (genParticle.vz() - LLP1_z)*(genParticle.vz() - LLP1_z);
//             std::cout<<"dV for nllp==1 : "<<dV<<std::endl;
	  if ( dV > 0.01 ) {
	    nllp = 2;
	    LLP2_x = genParticle.vx();
	    LLP2_y = genParticle.vy();
	    LLP2_z = genParticle.vz();
//$$
	    LLP2_dist = TMath::Sqrt( (LLP2_x - tree_GenPVx)*(LLP2_x - tree_GenPVx) 
				   + (LLP2_y - tree_GenPVy)*(LLP2_y - tree_GenPVy) 
				   + (LLP2_z - tree_GenPVz)*(LLP2_z - tree_GenPVz) ); 
//$$
	  }
	}
	if ( nllp == 0 ) {
	  nllp = 1;
	  LLP1_x = genParticle.vx();
	  LLP1_y = genParticle.vy();
	  LLP1_z = genParticle.vz();
//$$
	  LLP1_dist = TMath::Sqrt( (LLP1_x - tree_GenPVx)*(LLP1_x - tree_GenPVx) 
				 + (LLP1_y - tree_GenPVy)*(LLP1_y - tree_GenPVy) 
				 + (LLP1_z - tree_GenPVz)*(LLP1_z - tree_GenPVz) ); 
//$$
	}
      }
      
    if (genParticle.pt() < 0.9 || fabs(genParticle.eta()) > 4.0) continue;
      
      tree_genParticle_pt.push_back(genParticle.pt());
      tree_genParticle_eta.push_back(genParticle.eta());
      tree_genParticle_phi.push_back(genParticle.phi());
      tree_genParticle_charge.push_back(genParticle.charge());
      tree_genParticle_pdgId.push_back(genParticle.pdgId());
      tree_genParticle_vx.push_back(genParticle.vx());
      tree_genParticle_vy.push_back(genParticle.vy());
      tree_genParticle_vz.push_back(genParticle.vz());
      tree_genParticle_mass.push_back(genParticle.mass());
      tree_genParticle_statusCode.push_back(genParticle.status());
      tree_genParticle_mother_pdgId.push_back( mom ? mom->pdgId() :  -1 );

    } // end loop on gen particles
    
    tree_nLLP = nllp;
    
//$$
    nFromC = 0; 
    nFromB = 0;
    for (size_t i = 0; i < genParticles->size(); ++i) { // loop on gen particles
      const GenParticle & genIt = (*genParticles)[i];
      int ID = abs(genIt.pdgId());
      unsigned int nDaughters = genIt.numberOfDaughters();

      int fromLLP = -1;
      if ( (ID/100)%10 == 4 || (ID/1000)%10 == 4 || (ID/100)%10 == 5 || (ID/1000)%10 == 5 ) {
	float dV0 = (genIt.vx() - tree_GenPVx)*(genIt.vx() - tree_GenPVx)
	    	  + (genIt.vy() - tree_GenPVy)*(genIt.vy() - tree_GenPVy)
	    	  + (genIt.vz() - tree_GenPVz)*(genIt.vz() - tree_GenPVz);
	float dV1 = (genIt.vx() - LLP1_x)*(genIt.vx() - LLP1_x)
	    	  + (genIt.vy() - LLP1_y)*(genIt.vy() - LLP1_y)
	    	  + (genIt.vz() - LLP1_z)*(genIt.vz() - LLP1_z);
	float dV2 = (genIt.vx() - LLP2_x)*(genIt.vx() - LLP2_x)
	    	  + (genIt.vy() - LLP2_y)*(genIt.vy() - LLP2_y)
	    	  + (genIt.vz() - LLP2_z)*(genIt.vz() - LLP2_z);
        if      ( dV1 < dV2 && dV1 < 0.01 ) fromLLP = 1;
        else if ( dV2 < dV1 && dV2 < 0.01 ) fromLLP = 2;
        else if ( dV0 < 0.01 )              fromLLP = 0;
      }

      // Final c Hadron and get all its final charged particles
      bool isFinalD = false;
      if ( (ID/100)%10 == 4 || (ID/1000)%10 == 4 ) {
        isFinalD = true;
        for (unsigned int d1=0; d1<nDaughters; d1++) {
          const Candidate* gen1 = genIt.daughter(d1);
          int ID1 = abs(gen1->pdgId());
          if ( (ID1/100)%10 == 4 || (ID1/1000)%10 == 4 ) isFinalD = false;
        }
      }
      if ( isFinalD && abs(genIt.eta()) < 4. ) {
        for (unsigned int d1=0; d1<nDaughters; d1++) {
          const Candidate* gen1 = genIt.daughter(d1);
          unsigned int nDaughters2 = gen1->numberOfDaughters();
          if ( nDaughters2 == 0 ) {
            if ( gen1->charge() != 0 && gen1->pt() > 0.9 ) {
              nFromC++;
              tree_genFromC_pt.push_back(    gen1->pt());
              tree_genFromC_eta.push_back(   gen1->eta());
              tree_genFromC_phi.push_back(   gen1->phi());
              tree_genFromC_charge.push_back(gen1->charge());
              tree_genFromC_pdgId.push_back( gen1->pdgId());
              tree_genFromC_vx.push_back(    gen1->vx());
              tree_genFromC_vy.push_back(    gen1->vy());
              tree_genFromC_vz.push_back(    gen1->vz());
              tree_genFromC_mother_pdgId.push_back(    genIt.pdgId());
              tree_genFromC_generation.push_back(1);
              tree_genFromC_LLP.push_back(fromLLP);
            }
          }
          else {
            for (unsigned int d2=0; d2<nDaughters2; d2++) {
              const Candidate* gen2 = gen1->daughter(d2);
              unsigned int nDaughters3 = gen2->numberOfDaughters();
              if ( nDaughters3 == 0 ) {
                if ( gen2->charge() != 0 && gen2->pt() > 0.9 ) {
                  nFromC++;
                  tree_genFromC_pt.push_back(    gen2->pt());
                  tree_genFromC_eta.push_back(   gen2->eta());
                  tree_genFromC_phi.push_back(   gen2->phi());
                  tree_genFromC_charge.push_back(gen2->charge());
                  tree_genFromC_pdgId.push_back( gen2->pdgId());
                  tree_genFromC_vx.push_back(    gen2->vx());
                  tree_genFromC_vy.push_back(    gen2->vy());
                  tree_genFromC_vz.push_back(    gen2->vz());
                  tree_genFromC_mother_pdgId.push_back( gen1 ? gen1->pdgId() :  -1 );
                  tree_genFromC_generation.push_back(2);
                  tree_genFromC_LLP.push_back(fromLLP);
                }
              }
              else {
        	for (unsigned int d3=0; d3<nDaughters3; d3++) {
        	  const Candidate* gen3 = gen2->daughter(d3);
                  unsigned int nDaughters4 = gen3->numberOfDaughters();
                  if ( nDaughters4 == 0 ) {
                    if ( gen3->charge() != 0 && gen3->pt() > 0.9 ) {
        	      nFromC++;
                      tree_genFromC_pt.push_back(    gen3->pt());
                      tree_genFromC_eta.push_back(   gen3->eta());
                      tree_genFromC_phi.push_back(   gen3->phi());
                      tree_genFromC_charge.push_back(gen3->charge());
                      tree_genFromC_pdgId.push_back( gen3->pdgId());
                      tree_genFromC_vx.push_back(    gen3->vx());
                      tree_genFromC_vy.push_back(    gen3->vy());
                      tree_genFromC_vz.push_back(    gen3->vz());
                      tree_genFromC_mother_pdgId.push_back( gen2 ? gen2->pdgId() :  -1 );
                      tree_genFromC_generation.push_back(3);
                      tree_genFromC_LLP.push_back(fromLLP);
        	    }
        	  }
		  else {
        	    for (unsigned int d4=0; d4<nDaughters4; d4++) {
        	      const Candidate* gen4 = gen3->daughter(d4);
                      unsigned int nDaughters5 = gen4->numberOfDaughters();
                      if ( nDaughters5 == 0 ) {
                        if ( gen4->charge() != 0 && gen4->pt() > 0.9 ) {
        	          nFromC++;
                          tree_genFromC_pt.push_back(	 gen4->pt());
                          tree_genFromC_eta.push_back(	 gen4->eta());
                          tree_genFromC_phi.push_back(	 gen4->phi());
                          tree_genFromC_charge.push_back(gen4->charge());
                          tree_genFromC_pdgId.push_back( gen4->pdgId());
                          tree_genFromC_vx.push_back(	 gen4->vx());
                          tree_genFromC_vy.push_back(	 gen4->vy());
                          tree_genFromC_vz.push_back(	 gen4->vz());
                          tree_genFromC_mother_pdgId.push_back( gen3 ? gen3->pdgId() :  -1 );
                          tree_genFromC_generation.push_back(4);
                          tree_genFromC_LLP.push_back(fromLLP);
        	        }
        	      }
		      else {
        	        for (unsigned int d5=0; d5<nDaughters5; d5++) {
        	          const Candidate* gen5 = gen4->daughter(d5);
                          unsigned int nDaughters6 = gen5->numberOfDaughters();
                          if ( nDaughters6 == 0 ) {
                            if ( gen5->charge() != 0 && gen5->pt() > 0.9 ) {
        	              nFromC++;
                              tree_genFromC_pt.push_back(    gen5->pt());
                              tree_genFromC_eta.push_back(   gen5->eta());
                              tree_genFromC_phi.push_back(   gen5->phi());
                              tree_genFromC_charge.push_back(gen5->charge());
                              tree_genFromC_pdgId.push_back( gen5->pdgId());
                              tree_genFromC_vx.push_back(    gen5->vx());
                              tree_genFromC_vy.push_back(    gen5->vy());
                              tree_genFromC_vz.push_back(    gen5->vz());
                              tree_genFromC_mother_pdgId.push_back( gen4 ? gen4->pdgId() :  -1 );
                              tree_genFromC_generation.push_back(5);
                              tree_genFromC_LLP.push_back(fromLLP);
        	            }
        	          }
		          else {
        	            for (unsigned int d6=0; d6<nDaughters6; d6++) {
        	              const Candidate* gen6 = gen5->daughter(d6);
                              unsigned int nDaughters7 = gen6->numberOfDaughters();
                              if ( nDaughters7 == 0 ) {
                                if ( gen6->charge() != 0 && gen6->pt() > 0.9 ) {
        	                  nFromC++;
                                  tree_genFromC_pt.push_back(    gen6->pt());
                                  tree_genFromC_eta.push_back(   gen6->eta());
                                  tree_genFromC_phi.push_back(   gen6->phi());
                                  tree_genFromC_charge.push_back(gen6->charge());
                                  tree_genFromC_pdgId.push_back( gen6->pdgId());
                                  tree_genFromC_vx.push_back(    gen6->vx());
                                  tree_genFromC_vy.push_back(    gen6->vy());
                                  tree_genFromC_vz.push_back(    gen6->vz());
                                  tree_genFromC_mother_pdgId.push_back( gen5 ? gen5->pdgId() :  -1 );
                                  tree_genFromC_generation.push_back(6);
                                  tree_genFromC_LLP.push_back(fromLLP);
        	                }
        	              }
		              else {
        	                for (unsigned int d7=0; d7<nDaughters7; d7++) {
        	                  const Candidate* gen7 = gen6->daughter(d7);
                                  unsigned int nDaughters8 = gen7->numberOfDaughters();
                                  if ( nDaughters8 == 0 ) {
                                    if ( gen7->charge() != 0 && gen7->pt() > 0.9 ) {
        	                      nFromC++;
                                      tree_genFromC_pt.push_back(    gen7->pt());
                                      tree_genFromC_eta.push_back(   gen7->eta());
                                      tree_genFromC_phi.push_back(   gen7->phi());
                                      tree_genFromC_charge.push_back(gen7->charge());
                                      tree_genFromC_pdgId.push_back( gen7->pdgId());
                                      tree_genFromC_vx.push_back(    gen7->vx());
                                      tree_genFromC_vy.push_back(    gen7->vy());
                                      tree_genFromC_vz.push_back(    gen7->vz());
                                      tree_genFromC_mother_pdgId.push_back( gen6 ? gen6->pdgId() :  -1 );
                                      tree_genFromC_generation.push_back(7);
                                      tree_genFromC_LLP.push_back(fromLLP);
        	                    }
        	                  }
        	                }
        	              }
        	            }
        	          }
        	        }
        	      }
        	    }
        	  }
        	}
              }
            }
          }
        }
      } // final D hadron

      // Final b Hadron and get all its final charged particles
      bool isFinalB = false;
      if ( (ID/100)%10 == 5 || (ID/1000)%10 == 5 ) {
        isFinalB = true;
        for (unsigned int d1=0; d1<nDaughters; d1++) {
          const Candidate* gen1 = genIt.daughter(d1);
          int ID1 = abs(gen1->pdgId());
          if ( (ID1/100)%10 == 5 || (ID1/1000)%10 == 5 ) isFinalB = false;
        }
      }
      if ( isFinalB && abs(genIt.eta()) < 4. ) {
        for (unsigned int d1=0; d1<nDaughters; d1++) {
          const Candidate* gen1 = genIt.daughter(d1);
          unsigned int nDaughters2 = gen1->numberOfDaughters();
          if ( nDaughters2 == 0 ) {
            if ( gen1->charge() != 0 && gen1->pt() > 0.9 ) {
              nFromB++;
              tree_genFromB_pt.push_back(    gen1->pt());
              tree_genFromB_eta.push_back(   gen1->eta());
              tree_genFromB_phi.push_back(   gen1->phi());
              tree_genFromB_charge.push_back(gen1->charge());
              tree_genFromB_pdgId.push_back( gen1->pdgId());
              tree_genFromB_vx.push_back(    gen1->vx());
              tree_genFromB_vy.push_back(    gen1->vy());
              tree_genFromB_vz.push_back(    gen1->vz());
              tree_genFromB_mother_pdgId.push_back(    genIt.pdgId());
              tree_genFromB_generation.push_back(1);
              tree_genFromB_LLP.push_back(fromLLP);
            }
          }
          else {
            for (unsigned int d2=0; d2<nDaughters2; d2++) {
              const Candidate* gen2 = gen1->daughter(d2);
              unsigned int nDaughters3 = gen2->numberOfDaughters();
              if ( nDaughters3 == 0 ) {
                if ( gen2->charge() != 0 && gen2->pt() > 0.9 ) {
                  nFromB++;
                  tree_genFromB_pt.push_back(    gen2->pt());
                  tree_genFromB_eta.push_back(   gen2->eta());
                  tree_genFromB_phi.push_back(   gen2->phi());
                  tree_genFromB_charge.push_back(gen2->charge());
                  tree_genFromB_pdgId.push_back( gen2->pdgId());
                  tree_genFromB_vx.push_back(    gen2->vx());
                  tree_genFromB_vy.push_back(    gen2->vy());
                  tree_genFromB_vz.push_back(    gen2->vz());
                  tree_genFromB_mother_pdgId.push_back( gen1 ? gen1->pdgId() :  -1 );
                  tree_genFromB_generation.push_back(2);
                  tree_genFromB_LLP.push_back(fromLLP);
                }
              }
              else {
        	for (unsigned int d3=0; d3<nDaughters3; d3++) {
        	  const Candidate* gen3 = gen2->daughter(d3);
                  unsigned int nDaughters4 = gen3->numberOfDaughters();
                  if ( nDaughters4 == 0 ) {
                    if ( gen3->charge() != 0 && gen3->pt() > 0.9 ) {
        	      nFromB++;
                      tree_genFromB_pt.push_back(    gen3->pt());
                      tree_genFromB_eta.push_back(   gen3->eta());
                      tree_genFromB_phi.push_back(   gen3->phi());
                      tree_genFromB_charge.push_back(gen3->charge());
                      tree_genFromB_pdgId.push_back( gen3->pdgId());
                      tree_genFromB_vx.push_back(    gen3->vx());
                      tree_genFromB_vy.push_back(    gen3->vy());
                      tree_genFromB_vz.push_back(    gen3->vz());
                      tree_genFromB_mother_pdgId.push_back( gen2 ? gen2->pdgId() :  -1 );
                      tree_genFromB_generation.push_back(3);
                      tree_genFromB_LLP.push_back(fromLLP);
        	    }
        	  }
		  else {
        	    for (unsigned int d4=0; d4<nDaughters4; d4++) {
        	      const Candidate* gen4 = gen3->daughter(d4);
                      unsigned int nDaughters5 = gen4->numberOfDaughters();
                      if ( nDaughters5 == 0 ) {
                        if ( gen4->charge() != 0 && gen4->pt() > 0.9 ) {
        	          nFromB++;
                          tree_genFromB_pt.push_back(    gen4->pt());
                          tree_genFromB_eta.push_back(   gen4->eta());
                          tree_genFromB_phi.push_back(   gen4->phi());
                          tree_genFromB_charge.push_back(gen4->charge());
                          tree_genFromB_pdgId.push_back( gen4->pdgId());
                          tree_genFromB_vx.push_back(    gen4->vx());
                          tree_genFromB_vy.push_back(    gen4->vy());
                          tree_genFromB_vz.push_back(    gen4->vz());
                          tree_genFromB_mother_pdgId.push_back( gen3 ? gen3->pdgId() :  -1 );
                          tree_genFromB_generation.push_back(4);
                          tree_genFromB_LLP.push_back(fromLLP);
        	        }
        	      }
		      else {
        	        for (unsigned int d5=0; d5<nDaughters5; d5++) {
        	          const Candidate* gen5 = gen4->daughter(d5);
                          unsigned int nDaughters6 = gen5->numberOfDaughters();
                          if ( nDaughters6 == 0 ) {
                            if ( gen5->charge() != 0 && gen5->pt() > 0.9 ) {
        	              nFromB++;
                              tree_genFromB_pt.push_back(    gen5->pt());
                              tree_genFromB_eta.push_back(   gen5->eta());
                              tree_genFromB_phi.push_back(   gen5->phi());
                              tree_genFromB_charge.push_back(gen5->charge());
                              tree_genFromB_pdgId.push_back( gen5->pdgId());
                              tree_genFromB_vx.push_back(    gen5->vx());
                              tree_genFromB_vy.push_back(    gen5->vy());
                              tree_genFromB_vz.push_back(    gen5->vz());
                              tree_genFromB_mother_pdgId.push_back( gen4 ? gen4->pdgId() :  -1 );
                              tree_genFromB_generation.push_back(5);
                              tree_genFromB_LLP.push_back(fromLLP);
        	            }
        	          }
		          else {
        	            for (unsigned int d6=0; d6<nDaughters6; d6++) {
        	              const Candidate* gen6 = gen5->daughter(d6);
                              unsigned int nDaughters7 = gen6->numberOfDaughters();
                              if ( nDaughters7 == 0 ) {
                                if ( gen6->charge() != 0 && gen6->pt() > 0.9 ) {
        	                  nFromB++;
                                  tree_genFromB_pt.push_back(    gen6->pt());
                                  tree_genFromB_eta.push_back(   gen6->eta());
                                  tree_genFromB_phi.push_back(   gen6->phi());
                                  tree_genFromB_charge.push_back(gen6->charge());
                                  tree_genFromB_pdgId.push_back( gen6->pdgId());
                                  tree_genFromB_vx.push_back(    gen6->vx());
                                  tree_genFromB_vy.push_back(    gen6->vy());
                                  tree_genFromB_vz.push_back(    gen6->vz());
                                  tree_genFromB_mother_pdgId.push_back( gen5 ? gen5->pdgId() :  -1 );
                                  tree_genFromB_generation.push_back(6);
                                  tree_genFromB_LLP.push_back(fromLLP);
                  	        }
        	              }
		              else {
        	                for (unsigned int d7=0; d7<nDaughters7; d7++) {
        	                  const Candidate* gen7 = gen6->daughter(d7);
                                  unsigned int nDaughters8 = gen7->numberOfDaughters();
                                  if ( nDaughters8 == 0 ) {
                                    if ( gen7->charge() != 0 && gen7->pt() > 0.9 ) {
        	                      nFromB++;
                                      tree_genFromB_pt.push_back(    gen7->pt());
                                      tree_genFromB_eta.push_back(   gen7->eta());
                                      tree_genFromB_phi.push_back(   gen7->phi());
                                      tree_genFromB_charge.push_back(gen7->charge());
                                      tree_genFromB_pdgId.push_back( gen7->pdgId());
                                      tree_genFromB_vx.push_back(    gen7->vx());
                                      tree_genFromB_vy.push_back(    gen7->vy());
                                      tree_genFromB_vz.push_back(    gen7->vz());
                                      tree_genFromB_mother_pdgId.push_back( gen6 ? gen6->pdgId() :  -1 );
                                      tree_genFromB_generation.push_back(7);
                                      tree_genFromB_LLP.push_back(fromLLP);
                  	            }
        	                  }
		                  else {
        	                    for (unsigned int d8=0; d8<nDaughters8; d8++) {
        	                      const Candidate* gen8 = gen7->daughter(d8);
                                      unsigned int nDaughters9 = gen8->numberOfDaughters();
                                      if ( nDaughters9 == 0 ) {
                                        if ( gen8->charge() != 0 && gen8->pt() > 0.9 ) {
        	                          nFromB++;
                                          tree_genFromB_pt.push_back(    gen8->pt());
                                          tree_genFromB_eta.push_back(   gen8->eta());
                                          tree_genFromB_phi.push_back(   gen8->phi());
                                          tree_genFromB_charge.push_back(gen8->charge());
                                          tree_genFromB_pdgId.push_back( gen8->pdgId());
                                          tree_genFromB_vx.push_back(    gen8->vx());
                                          tree_genFromB_vy.push_back(    gen8->vy());
                                          tree_genFromB_vz.push_back(    gen8->vz());
                                          tree_genFromB_mother_pdgId.push_back( gen7 ? gen7->pdgId() :  -1 );
                                          tree_genFromB_generation.push_back(8);
                                          tree_genFromB_LLP.push_back(fromLLP);
                  	                }
        	                      }
		                      else {
        	                        for (unsigned int d9=0; d9<nDaughters9; d9++) {
        	                          const Candidate* gen9 = gen8->daughter(d9);
                                          unsigned int nDaughters10 = gen9->numberOfDaughters();
                                          if ( nDaughters10 == 0 ) {
                                            if ( gen9->charge() != 0 && gen9->pt() > 0.9 ) {
        	                              nFromB++;
                                              tree_genFromB_pt.push_back(    gen9->pt());
                                              tree_genFromB_eta.push_back(   gen9->eta());
                                              tree_genFromB_phi.push_back(   gen9->phi());
                                              tree_genFromB_charge.push_back(gen9->charge());
                                              tree_genFromB_pdgId.push_back( gen9->pdgId());
                                              tree_genFromB_vx.push_back(    gen9->vx());
                                              tree_genFromB_vy.push_back(    gen9->vy());
                                              tree_genFromB_vz.push_back(    gen9->vz());
                                              tree_genFromB_mother_pdgId.push_back( gen8 ? gen8->pdgId() :  -1 );
                                              tree_genFromB_generation.push_back(9);
                                              tree_genFromB_LLP.push_back(fromLLP);
                  	                    }
        	                          }
        	                        }
		                      }
        	                    }
		                  }
        	                }
		              }
        	            }
		          }
        	        }
		      }
        	    }
		  }
        	}
              }
            }
          }
        }
      } // final B hadron

    } // end loop on gen particles
//$$
    

    // GEN JETS
    for (auto const & genJet : *genJets)
    {
    if ( genJet.pt() < 20. ) continue;
      tree_genJet_pt.push_back(genJet.pt());
      tree_genJet_eta.push_back(genJet.eta());
      tree_genJet_phi.push_back(genJet.phi());
      tree_genJet_mass.push_back(genJet.mass());
      tree_genJet_energy.push_back(genJet.energy());
    }
    
    for (auto const & genJet : *ak8GenJets)
    {
    if ( genJet.pt() < 20. ) continue;
      tree_ak8GenJet_pt.push_back(genJet.pt());
      tree_ak8GenJet_eta.push_back(genJet.eta());
      tree_ak8GenJet_phi.push_back(genJet.phi());
      tree_ak8GenJet_mass.push_back(genJet.mass());
      tree_ak8GenJet_energy.push_back(genJet.energy());
    }
  } // endif simulation
    

//$$
  //////////////////////////////////
  //////////////////////////////////
  //////// HT FILTER CHECK /////////
  //////////////////////////////////
  //////////////////////////////////
  
  double HT_val = 0;
  bool MuonSup10GeV = false;
  bool MuonSup28GeV = false;
  
  for (unsigned int ij=0; ij<ak4slimmedJets->size(); ij++){
      const Jet& jet = ak4slimmedJets->at(ij);
      if ( fabs(jet.eta()) < 2.4 && jet.pt() > 20. ) HT_val += jet.pt();  // used in HTFilter !
  }
  
  int NumberMuonsSup10GeV = 0;
  int NumberMuonsSup28GeV = 0;
  bool invMassGood = false;
  
  for (auto const & muon : *slimmedmuons) {
      if ( !muon.isGlobalMuon() ) continue;
      if ( muon.pt() < 10. ) continue;
      NumberMuonsSup10GeV++;
      if ( muon.pt() > 28. ) NumberMuonsSup28GeV++;
      v1.SetPtEtaPhiM(muon.pt(), muon.eta(), muon.phi(), mu_mass);
      
      for (auto const & muon2 : *slimmedmuons) {
  	  if ( !muon2.isGlobalMuon() ) continue;
  	  if ( muon2.pt() < 28. ) continue;
  	  if ( muon2.pt() ==  muon.pt() ) continue;
  	  if ( muon2.charge() ==  muon.charge() ) continue;
  	  v2.SetPtEtaPhiM(muon2.pt(), muon2.eta(), muon2.phi(), mu_mass);
  	  v = v1 + v2;
  	  if ( v.Mag() > 60. ) {
  	      invMassGood = true;
  	      break;
  	  }
      }
  }
  if ( NumberMuonsSup10GeV >= 2 ) MuonSup10GeV = true;
  if ( NumberMuonsSup28GeV >= 1 ) MuonSup28GeV = true;
  
//$$
  // if ( tree_NbrOfZCand > 0 && HT_val > 180. && MuonSup10GeV && MuonSup28GeV && invMassGood ) {
//$$
if ( HT_val > 180. && MuonSup10GeV && MuonSup28GeV && invMassGood ) {
    tree_passesHTFilter = true;
  }
  else tree_passesHTFilter = false;
  tree_nTracks = -1; 
    
  if ( tree_passesHTFilter ) {
//$$

    //////////////////////////////////
    //////////////////////////////////
    //////////   Tracks  /////////////
    //////////////////////////////////
    //////////////////////////////////
    ///Tracks and their association to other objects
    
    edm::RefToBaseVector<pat::PackedCandidate> trackRefs;
    // edm::RefToBaseVector<edm::View<pat::PackedCandidate> > trackRefs;//does not work 

    // edm::RefToBaseVector<reco::Track> trackRefs;
    size_t idxTrack = 0;
    
    std::map<size_t , int > trackToVextexMap;
    std::map<size_t , int > trackToAK4SlimmedJetMap;
    std::map<size_t , int > trackToAK4PFJetMap;
    std::map<size_t , int > trackTo08JetMap;
    std::map<size_t , int > trackToCaloJetMap;
    std::map<size_t , int > jetToTrackMap;
    std::map<size_t , int > CaloJetToTrackMap;
    std::map<size_t , int > ak8jetToTrackMap;
    std::map<size_t , int > trackToPFJetMap;
    
    ///// TRACK ASSOCIATION TO VERTICES AND JETS
    
    for (edm::View<Track>::size_type i=0; i<tracks.size(); ++i)
    {
      //---------------------------
      //minimum selection on tracks
    if ( tracks[i].pt() < 0.9 || fabs(tracks[i].eta()) > 2.5 ) continue;
      trackRefs.push_back(tracks.refAt(i));
      
      //---------------------------
      //loop on vertex to do track-vertex association
      //---------------------------

      int idxvtx = 0;
      int thematchidx = -1;
      float dzmin = std::numeric_limits<float>::max();
      for (auto const & vertex : *vertices) {
        float dz = std::abs(tracks[i].dz( vertex.position()));
        if ( dz < dzmin ) {
          dzmin = dz;
          thematchidx = idxvtx;
        }
        idxvtx++;
      }
      trackToVextexMap[idxTrack] = thematchidx;
      
      //---------------------------
      //loop on jets to do track-jet association
      //---------------------------
      
      int idxSlimmedJet = 0;
      bool found_match = false;
      for (unsigned int ij=0;ij<ak4slimmedJets->size();ij++) {
        const Jet& jet = ak4slimmedJets->at(ij);
      if ( jet.pt() < 20. ) continue;
        TLorentzVector jet4m, track4m;
        jet4m.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0);
        track4m.SetPtEtaPhiM(tracks[i].pt(), tracks[i].eta(), tracks[i].phi(), 0);
        if ( jet4m.DeltaR(track4m) < 0.4 ) {
          found_match = true;
          break;
        }
        else idxSlimmedJet++;
      }
      if (found_match) trackToAK4SlimmedJetMap[idxTrack] = idxSlimmedJet;
      else	       trackToAK4SlimmedJetMap[idxTrack] = -1;
      
      // int idxAK4PFJet=0;
      // bool found_match_ak4 = false;
      // for (unsigned int ij=0;ij<ak4PFJets->size();ij++) {
      //   const Jet& jet = ak4PFJets->at(ij);
      // if ( jet.pt() < 20. ) continue;
      //   TLorentzVector jet4m, track4m;
      //   jet4m.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0);
      //   track4m.SetPtEtaPhiM(tracks[i].pt(), tracks[i].eta(), tracks[i].phi(), 0);
      //   if ( jet4m.DeltaR(track4m) < 0.4 ) {
      //     found_match_ak4 = true;
      //     break;
      //   }
      //   else idxAK4PFJet++;
      // }
      // if (found_match_ak4) trackToAK4PFJetMap[idxTrack] = idxAK4PFJet;
      // else		   trackToAK4PFJetMap[idxTrack] = -1;
      
      int idxJet08=0;
      bool found_match08 = false;
      for (unsigned int ij=0;ij<ak8jets->size();ij++) {
        const Jet& jet = ak8jets->at(ij);
      if ( jet.pt() < 20. ) continue;
        TLorentzVector jet4m, track4m;
        jet4m.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0);
        track4m.SetPtEtaPhiM(tracks[i].pt(), tracks[i].eta(), tracks[i].phi(), 0);
        if ( jet4m.DeltaR(track4m) < 0.8 ) {
          found_match08 = true;
          break;
        }
        else idxJet08++;
      }
      if (found_match08) trackTo08JetMap[idxTrack] = idxJet08;
      else		 trackTo08JetMap[idxTrack] = -1;
      
      // int idxCaloJet=0;
      // bool found_match_calo = false;
      // for (unsigned int ij=0;ij<CaloJets->size();ij++) {
      //   const Jet& jet = CaloJets->at(ij);
      // if ( jet.pt() < 20. ) continue;
      //   TLorentzVector jet4m, track4m;
      //   jet4m.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0);
      //   track4m.SetPtEtaPhiM(tracks[i].pt(), tracks[i].eta(), tracks[i].phi(), 0);
      //   if ( jet4m.DeltaR(track4m) < 0.4 ) {
      //     found_match_calo = true;
      //     break;
      //   }
      //   else idxCaloJet++;
      // }
      // if (found_match_calo) trackToCaloJetMap[idxTrack] = idxCaloJet;
      // else		    trackToCaloJetMap[idxTrack] = -1;
      
      idxTrack++;
    }
    
    //reco::BeamSpot vertexBeamSpot= *recoBeamSpotHandle;
    edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);
    
    /////////////////
    //prepare association to tracks by hit
    // reco::RecoToSimCollection recSimColl ;
    // if ( !runOnData_ ) recSimColl= associatorByHits.associateRecoToSim(trackRefs, tpCollection);
    //fake rate determination : when a reco track has no matched simtrack
    
    //    int nUnmatchTrack_fromPU = 0;
    //    int nUnmatchTrack_fromPV = 0;
    //    int nUnmatchTrack= 0;
    //    int nPUTrack= 0;
    
    //---------------
    // loop on tracks
    //---------------
    
    tree_nTracks = 0; 
    for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack) {
      // for(View<pat::PackedCandidate>::const_iterator itTrack = tracks->begin();
      const auto& itTrack = trackRefs[iTrack];
      if (itTrack->hasTrackDetails())
        {
      tree_nTracks++; 
      //------------------------
      //general track properties
      //------------------------
      //std::cout << "  " << itTrack->pt()  << std::endl;
      tree_track_charge.push_back(itTrack->charge());
      tree_track_pt.push_back(itTrack->pt());
      // tree_track_ptError.push_back(itTrack->ptError());//reco
      tree_track_eta.push_back(itTrack->eta());
      tree_track_phi.push_back(itTrack->phi());
      // tree_track_nhits.push_back(itTrack->numberOfValidHits());//reco
      tree_track_nhits.push_back(itTrack->numberOfHits());
      // tree_track_NChi2.push_back(itTrack->normalizedChi2());//reco

          tree_track_NChi2.push_back(itTrack->pseudoTrack().normalizedChi2());
      // tree_track_outerPt.push_back(itTrack->outerPt());//reco

        tree_track_vx.push_back(itTrack->vx());
        tree_track_vy.push_back(itTrack->vy());
        tree_track_vz.push_back(itTrack->vz());

      tree_track_firsthit_X.push_back(itTrack->firstHit());//miniaod
      tree_track_firsthit_Y.push_back(itTrack->firstHit());
      tree_track_firsthit_Z.push_back(itTrack->firstHit());
      
      // tree_track_firsthit_X.push_back(itTrack->innerPosition().X());//reco
      // tree_track_firsthit_Y.push_back(itTrack->innerPosition().Y());
      // tree_track_firsthit_Z.push_back(itTrack->innerPosition().Z());
      // tree_track_firsthit_phi.push_back(itTrack->innerPosition().phi());
      
      //----------------reco--------------------------//
      // if( itTrack->quality(reco::TrackBase::highPurity) ){tree_track_isHighQuality.push_back(true);}
      // else {tree_track_isHighQuality.push_back(false);}
      // if( itTrack->quality(reco::TrackBase::loose) )	 {tree_track_isLoose.push_back(true);}
      // else {tree_track_isLoose.push_back(false);}
      // if( itTrack->quality(reco::TrackBase::tight))	 {tree_track_isTight.push_back(true);}
      // else {tree_track_isTight.push_back(false);}

      //----------------miniaod--------------------//
      if( itTrack->trackHighPurity() ){tree_track_isHighQuality.push_back(true);}
      else {tree_track_isHighQuality.push_back(false);}

      //--------------------------------------------//

      tree_track_dxy.push_back( 	itTrack->dxy(bs.position()));
      tree_track_dxyError.push_back(	 itTrack->dxyError());
      tree_track_dz.push_back(         itTrack->dz(bs.position()));
      tree_track_dzError.push_back(	itTrack->dzError());
      // tree_track_numberOfLostHits.push_back( itTrack->numberOfLostHits());//reco
      // tree_track_numberOfValidHits.push_back(itTrack->numberOfValidHits());//reco=> nhits
      
      // tree_track_originalAlgo.push_back(itTrack->originalAlgo());
      // tree_track_algo.push_back(itTrack->algo());
      // tree_track_stopReason.push_back(itTrack->stopReason());
      
      //--------------------------------
      //general hit properties of tracks
      //--------------------------------
      
      // const reco::HitPattern& hp = itTrack->hitPattern();
      //     tree_track_nPixel   .push_back(hp.numberOfValidPixelHits());
      //     tree_track_nStrip   .push_back(hp.numberOfValidStripHits());
      
      // tree_track_numberOfValidPixelHits.push_back(hp.numberOfValidPixelHits());
      // tree_track_numberOfValidStripHits.push_back(hp.numberOfValidStripHits());
      // tree_track_numberOfValidStripTIBHits.push_back(hp.numberOfValidStripTIBHits());
      // tree_track_numberOfValidStripTIDHits.push_back(hp.numberOfValidStripTIDHits());
      // tree_track_numberOfValidStripTOBHits.push_back(hp.numberOfValidStripTOBHits());
      // tree_track_numberOfValidStripTECHits.push_back(hp.numberOfValidStripTECHits());
      // tree_track_numberOfValidPixelBarrelHits.push_back(hp.numberOfValidPixelBarrelHits());
      // tree_track_numberOfValidPixelEndcapHits.push_back(hp.numberOfValidPixelEndcapHits());
      // tree_track_trackerLayersWithMeasurement.push_back(hp.trackerLayersWithMeasurement());
      // tree_track_pixelLayersWithMeasurement.push_back(hp.pixelLayersWithMeasurement());
      // tree_track_stripTECLayersWithMeasurement.push_back(hp.stripTECLayersWithMeasurement() );
      // tree_track_stripTIBLayersWithMeasurement.push_back(hp.stripTIBLayersWithMeasurement());
      // tree_track_stripTIDLayersWithMeasurement.push_back(hp.stripTIDLayersWithMeasurement());
      // tree_track_stripTOBLayersWithMeasurement.push_back(hp.stripTOBLayersWithMeasurement());
      
      // int hitPixelLayer = 0;
      // if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) )  hitPixelLayer += 1;
      // if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2) )  hitPixelLayer += 10;
      // if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3) )  hitPixelLayer += 100;
      // if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 4) )  hitPixelLayer += 1000;
      
      // if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1) )  hitPixelLayer += 2;
      // if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 2) )  hitPixelLayer += 20;
      // if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 3) )  hitPixelLayer += 200;
      
      // tree_track_hasValidHitInPixelLayer.push_back(hitPixelLayer);
      //----------MINIAOD CONSTRAINTS---------------------//
      //----------------------------
      //matching to simulated tracks
      //----------------------------
      //matching par hits
      
//       if ( !runOnData_ ) 
//       {
//         int nSimHits = 0;
//         bool isSimMatched = false;
//         std::vector<int> tpIdx;
//         std::vector<float> sharedFraction;
//         std::vector<float> tpChi2;
        
//         //initialized values for trackingParticle
//         float sim_vx = -1000;
//         float sim_vy = -1000;
//         float sim_vz = -1000;
        
//         int   sim_charge      =-1000;
//         float sim_pt	      =-1;
//         float sim_eta	      =-1000;
//         float sim_phi	      =-1000;
//         bool  sim_longLived   = false;
//         //int	sim_matchedHit	= 0;
//         int   sim_pdgId	      = -1000;
//         int   sim_numberOfTrackerHits   = 0;
//         int   sim_numberOfTrackerLayers = 0;
//         float sim_mass		= 0.;
//         int   sim_status	= -1000;
//         int   sim_isFromLLP	= 0;
// //$$
//         int   sim_isFromBC	        = 0;
//         int   sim_isFromBC_mother_pdgId = -1000;
//         int   sim_isFromBC_LLP	        = -1;
//         float sim_dV	                = 1000.;
//         float sim_dist	                = -10.;

        
//         auto foundTPs = recSimColl.find(itTrack->pseudoTrack());
//         if ( foundTPs != recSimColl.end() ) 
//         {
//           //if (!foundTPs->val.empty()) {
//           isSimMatched = true;
//           TrackingParticleRef tpr = foundTPs->val[0].first;
//           nSimHits = tpr->numberOfTrackerHits();
          
//           sim_charge	         = tpr->charge();
//           sim_pt		 = tpr->pt();
//           sim_eta  	         = tpr->eta();
//           sim_phi  	         = tpr->phi();
//           sim_longLived	         = tpr->longLived();
//           // sim_matchedHit	    = tpr->matchedHit();
//           sim_pdgId	         = tpr->pdgId();
//           sim_numberOfTrackerHits   = tpr->numberOfTrackerHits();
//           sim_numberOfTrackerLayers = tpr->numberOfTrackerLayers();
//           sim_mass 	         = tpr->mass();
//           sim_status	         = tpr->status();
          
//           //determine x,y,z position of the genVertex which produced the associated simtrack
//           sim_vx	 = tpr->vx();
//           sim_vy	 = tpr->vy();
//           sim_vz	 = tpr->vz();
          
//           // std::cout<< "nllp : "<< nllp<<std::endl;
// 	  float dV1 = 1000.;
//           if ( nllp >= 1 ) {
//             dV1 = (sim_vx - LLP1_x)*(sim_vx - LLP1_x)
//                 + (sim_vy - LLP1_y)*(sim_vy - LLP1_y)
//                 + (sim_vz - LLP1_z)*(sim_vz - LLP1_z);
//             // std::cout<< "genVertexPosX : "<<sim_vx<<" & LLP1_X : "<< LLP1_x<< " & dV1 : "<< dV1<<std::endl;
//             if ( dV1 < 0.01 ) sim_isFromLLP = 1; //ATTENTION 0.01
//           }
// 	  float dV2 = 1000.;
//           if ( nllp >= 2 && sim_isFromLLP != 1 ) {
//             dV2 = (sim_vx - LLP2_x)*(sim_vx - LLP2_x)
//             	+ (sim_vy - LLP2_y)*(sim_vy - LLP2_y)
//             	+ (sim_vz - LLP2_z)*(sim_vz - LLP2_z);
//             if ( dV2 < 0.01 ) sim_isFromLLP = nllp; //ATTENTION 0.01
//           }

// //$$
//           sim_dist = TMath::Sqrt( (sim_vx-tree_GenPVx)*(sim_vx-tree_GenPVx) 
// 	                        + (sim_vy-tree_GenPVy)*(sim_vy-tree_GenPVy)
// 	                        + (sim_vz-tree_GenPVz)*(sim_vz-tree_GenPVz) );
//           if ( dV1 < dV2 ) {
// 	    sim_dV = TMath::Sqrt(dV1);
// 	    if ( sim_dist < LLP1_dist ) sim_dV = -sim_dV;
// 	  }
//           else {
// 	    sim_dV = TMath::Sqrt(dV2);
// 	    if ( sim_dist < LLP2_dist ) sim_dV = -sim_dV;
// 	  }
	
//           // match simtrack to genparticle from B hadron
// 	  bool matchB = false;
//           if ( nFromB > 0 ) {
//             for (int k=0; k<nFromB; k++) { // loop on genFromB
// 	    if ( tree_genFromB_pdgId[k] != sim_pdgId ) continue;
// 	      float dpt  = sim_pt / tree_genFromB_pt[k] - 1.;
// 	      float deta = sim_eta - tree_genFromB_eta[k];
//   	      float dphi = Deltaphi( sim_phi , tree_genFromB_phi[k] );
//               if ( abs(deta) < 0.01 && abs(dphi) < 0.01 && abs(dpt) < 0.01 ) {
// 	        matchB = true;
// 		sim_isFromBC_mother_pdgId = tree_genFromB_mother_pdgId[k];
// 	        sim_isFromBC_LLP = tree_genFromB_LLP[k];
//               }
// 	    } // end loop on genFromB
// 	  }
	 
//           // match simtrack to genparticle from C hadron
// 	  bool matchC = false;
//           if ( nFromC > 0 ) {
//             for (int k=0; k<nFromC; k++) { // loop on genFromC
// 	    if ( tree_genFromC_pdgId[k] != sim_pdgId ) continue;
// 	      float dpt  = sim_pt / tree_genFromC_pt[k] - 1.;
// 	      float deta = sim_eta - tree_genFromC_eta[k];
//   	      float dphi = Deltaphi( sim_phi , tree_genFromC_phi[k] );
//               if ( abs(deta) < 0.01 && abs(dphi) < 0.01 && abs(dpt) < 0.01 ) {
// 	        matchC = true;
// 		if ( !matchB ) {
// 		  sim_isFromBC_mother_pdgId = tree_genFromC_mother_pdgId[k];
// 	          sim_isFromBC_LLP = tree_genFromC_LLP[k];
// 		}
//               }
// 	    } // end loop on genFromC
//           }
	
// 	  if	  ( matchB && matchC)  sim_isFromBC = 6;
// 	  else if ( matchB && !matchC) sim_isFromBC = 5;
// 	  else if ( matchC )	       sim_isFromBC = 4;
// //$$
//         }
        
//         tree_track_nSimHits	             .push_back(nSimHits);
//         tree_track_isSimMatched              .push_back(isSimMatched);
        
//         tree_track_sim_charge	             .push_back(sim_charge);
//         tree_track_sim_pt  	             .push_back(sim_pt);
//         tree_track_sim_eta 	             .push_back(sim_eta);
//         tree_track_sim_phi 	             .push_back(sim_phi);
//         tree_track_sim_longLived             .push_back(sim_longLived);
//         tree_track_sim_pdgId	             .push_back(sim_pdgId);
//         tree_track_sim_numberOfTrackerHits   .push_back(sim_numberOfTrackerHits);
//         tree_track_sim_numberOfTrackerLayers .push_back(sim_numberOfTrackerLayers);
//         tree_track_sim_mass	             .push_back(sim_mass);
//         tree_track_sim_status	             .push_back(sim_status);
        
//         tree_track_sim_vx	             .push_back(sim_vx);
//         tree_track_sim_vy	             .push_back(sim_vy);
//         tree_track_sim_vz	             .push_back(sim_vz);
//         tree_track_sim_isFromLLP	     .push_back(sim_isFromLLP);
// //$$
//         tree_track_sim_isFromBC	             .push_back(sim_isFromBC);
//         tree_track_sim_isFromBC_mother_pdgId .push_back(sim_isFromBC_mother_pdgId);
//         tree_track_sim_isFromBC_LLP          .push_back(sim_isFromBC_LLP);
//         tree_track_sim_dV	             .push_back(sim_dV);
//         tree_track_sim_dist                  .push_back(sim_dist);
// //$$
//       }
      
      /*if( !isSimMatched &&  trackToVextexMap[iTrack] != 0 ) nUnmatchTrack_fromPU++;
       if( !isSimMatched &&  trackToVextexMap[iTrack] == 0 ) nUnmatchTrack_fromPV++;
       if( !isSimMatched		    ) nUnmatchTrack++;
       if( trackToVextexMap[iTrack] != 0	    ) nPUTrack++;*/
      
      tree_track_recoVertex_idx.push_back(trackToVextexMap[iTrack]);
      tree_track_recoAK4SlimmedJet_idx.push_back(trackToAK4SlimmedJetMap[iTrack]);
      tree_track_recoAK4PFJet_idx.push_back(trackToAK4PFJetMap[iTrack]);
      tree_track_reco08Jet_idx.push_back(trackTo08JetMap[iTrack]);
      tree_track_recoCaloJet_idx.push_back(trackToCaloJetMap[iTrack]);
    }
    } //end loop on tracks
    
    /*std::cout << "---------------------"  << std::endl;
     std::cout << "nUnmatchTrack_fromPU " << nUnmatchTrack_fromPU << std::endl;
     std::cout << "nUnmatchTrack_fromPV " << nUnmatchTrack_fromPV << std::endl;
     std::cout << "nUnmatchTrack        " << nUnmatchTrack<< std::endl;
     std::cout << "nTrack               " << tree_track_charge.size() << std::endl;
     std::cout << "nPUTrack             " << nPUTrack<< std::endl;
     std::cout << "nNonPUTrack          " << tree_track_charge.size()- nPUTrack<< std::endl;*/
    
    //----------MINIAOD CONSTRAINTS---------------------//
    //------------------------------
    // loops on tracking particles
    //------------------------------
    
    // if ( !runOnData_ ) {
    //   reco::SimToRecoCollection simRecColl = associatorByHits.associateSimToReco(trackRefs, tpCollection);
    //   //Reco efficiecny estimation by looking at the number of simtracks matched with reco tracks
    
    //   for (const TrackingParticleRef& tp: tpCollection)
    //   {
    //     tree_simtrack_charge           .push_back( tp->charge());    //infos sur toutes les traces simul�es,
    //     //ex efficacit� reco, voir trace sim � trace reco, compar�e aux pt de otoutes les traces sim /*!*///finito
    //     tree_simtrack_pt            .push_back( tp->pt());        //avec un rapport (attention de bien prendre un pt simul�)/*!*/
    //     tree_simtrack_eta           .push_back( tp->eta());
    //     tree_simtrack_phi           .push_back( tp->phi());
    //     tree_simtrack_longLived      .push_back( tp->longLived());
    //     //tree_simtrack_matchedHit         .push_back( tp->matchedHit());
    //     tree_simtrack_pdgId           .push_back( tp->pdgId());
    //     tree_simtrack_numberOfTrackerHits   .push_back( tp->numberOfTrackerHits());
    //     tree_simtrack_numberOfTrackerLayers .push_back( tp->numberOfTrackerLayers());
    //     tree_simtrack_mass           .push_back( tp->mass());
    //     tree_simtrack_status           .push_back( tp->status());
    //     tree_simtrack_vx           .push_back( tp->vx());
    //     tree_simtrack_vy           .push_back( tp->vy());
    //     tree_simtrack_vz           .push_back( tp->vz());
    
    //     bool isRecoMatched = false;
    
    //     std::vector<int> tkIdx;
    //     auto foundTracks = simRecColl.find(tp);
    //     if ( foundTracks != simRecColl.end() )
    //     {
    //       isRecoMatched = true;
    //       for (const auto trackQuality: foundTracks->val)
    //             {
    //                tkIdx.push_back(trackQuality.first.key());
    //             }
    //     }
    
    //     tree_simtrack_isRecoMatched.push_back(isRecoMatched);//tree_simtrack_pt when isRecoMatched =true
    //      //Calcualte the impact parameters w.r.t. PCA
    //      tree_simtrack_trkIdx.push_back(tkIdx);
    
    //      TrackingParticle::Vector momentum = parametersDefiner->momentum(iEvent,iSetup,tp);
    //      TrackingParticle::Point vertex = parametersDefiner->vertex(iEvent,iSetup,tp);
    //      auto dxySim = TrackingParticleIP::dxy(vertex, momentum, bs.position());
    //      auto dzSim  = TrackingParticleIP::dz(vertex, momentum, bs.position());
    
    //      tree_simtrack_pca_dxy      .push_back(dxySim);
    //      tree_simtrack_pca_dz      .push_back(dzSim);
    
    //   } //end loop on tracking particles
    
    // }
    

    /////////////////////////////////////////////////////////
    //-----------------------------------------
    // Jets for event axes
    //-----------------------------------------
    /////////////////////////////////////////////////////////
    
    int njet = 0, njet1 = 0, njet2 = 0;
    bool isjet[99], isjet1[99], isjet2[99];
    TLorentzVector vaxis1, vaxis2, vjet[99];
    float PtMin = 20;   // (GeV) minimum jet pt is optimum
    float EtaMax = 10.; // no cut on eta is optimum
    
    int njetall = ak4slimmedJets->size();
    for ( int jetidx=0; jetidx<njetall; jetidx++)  // Loop on jet
    {
      const Jet& jet = ak4slimmedJets->at(jetidx);
      float jet_pt  = jet.pt();
      float jet_eta = jet.eta();
      float jet_phi = jet.phi();
      isjet[jetidx]  = false;
      isjet1[jetidx] = false; // first neutralino jets
      isjet2[jetidx] = false; // second neutralino jets
      v.SetPtEtaPhiM( jet_pt, jet_eta, jet_phi, 0. ); //set the axis
      
    if ( jet_pt < PtMin ) continue;
    if ( abs(jet_eta) > EtaMax ) continue;
      
      // look if prompt muon inside
      float deltaR1 = 1000., deltaR2 = 1000.;
      if ( imu1 >= 0 ) deltaR1 = Deltar( jet_eta, jet_phi, tree_slimmedmuon_eta[imu1], tree_slimmedmuon_phi[imu1] );
      if ( imu2 >= 0 ) deltaR2 = Deltar( jet_eta, jet_phi, tree_slimmedmuon_eta[imu2], tree_slimmedmuon_phi[imu2] );
      if ( deltaR1 < 0.4 || deltaR2 < 0.4 )
      {
        if ( deltaR1 < 0.4 )
        { //if muon is inside, we remove the muons infomation from the jet
          v1.SetPtEtaPhiM( tree_slimmedmuon_pt[imu1],
          		  tree_slimmedmuon_eta[imu1],
          		  tree_slimmedmuon_phi[imu1],
          		  0 );
          v -= v1; //v TLorentzFactor being just above, defined by jet data
        }
        if ( deltaR2 < 0.4 )
        {
          v2.SetPtEtaPhiM( tree_slimmedmuon_pt[imu2],
          		  tree_slimmedmuon_eta[imu2],
          		  tree_slimmedmuon_phi[imu2],
          		  0 );
          v -= v2;
        }
        jet_pt  = v.Pt(); //Update jet data by removing the muons information (muons that could be in the jet)
        jet_eta = v.Eta(); //+ we do not want muons data to build the two axis since they come from the PV
        jet_phi = v.Phi();
      }
      
      njet++;
      isjet[jetidx] = true;
      vjet[jetidx] = v; // Only jet data (with  possible muons being removed)
      if ( njet1 == 0 && jet_pt > PtMin && abs(jet_eta) < EtaMax )
      {
        njet1 = 1;
        isjet1[jetidx] = true;
        vaxis1 = v;
      }
    } // End Loop on jets
    

    //////////////////
    //------------------------------
    // Event Axes
    //------------------------------
    //////////////////
    
    float dR, dR1 = 10., dR2 = 10.;
//$$
    float dRcut_hemis  = 1.5; // subjective choice
    float dRcut_tracks = 10.; // no cut is better (could bias low track pT and high LLP ct) 
//$$
     
    for (int i=0; i<njetall; i++) // Loop on jet
    {
    if ( !isjet[i] ) continue;
      // float jet_pt  = vjet[i].Pt();
      float jet_eta = vjet[i].Eta();
      float jet_phi = vjet[i].Phi();
      if ( njet1 > 0 ) dR1 = Deltar( jet_eta, jet_phi, vaxis1.Eta(), vaxis1.Phi() );
      if ( njet2 > 0 ) dR2 = Deltar( jet_eta, jet_phi, vaxis2.Eta(), vaxis2.Phi() );
      // axis 1
      if ( njet1 > 0 && !isjet2[i]  && dR1 < dRcut_hemis) {	     //
        njet1++;
        vaxis1 += vjet[i];
        isjet1[i] = true;
      }
      // axis 2
      if ( njet2 == 0 && !isjet1[i] ) {
        njet2 = 1;
        vaxis2 = vjet[i];
        isjet2[i] = true;
      }
      else if ( njet2 > 0 && !isjet1[i] && !isjet2[i] && dR2 < dRcut_hemis ) {//
        njet2++;
        vaxis2 += vjet[i];
        isjet2[i] = true;
      }
    }       // end Loop on jet
    std::cout<<"axes has details1"<<std::endl;
//$$$$
//     // force the axes to the true LLP
//     vaxis1 = vneu[0];
//     vaxis2 = vneu[1];
//$$$$

    ///////////////////////////////
    // compare with neutralino axis
    ///////////////////////////////
    
    int iLLPrec1 = 1, iLLPrec2 = 2;
    float axis1_eta = vaxis1.Eta();
    float axis1_phi = vaxis1.Phi();
    if ( neu[0] >= 0 ) dR1 = Deltar( axis1_eta, axis1_phi, Gen_neu1_eta, Gen_neu1_phi ); //dR between reco axis of jets and gen neutralino
    if ( neu[1] >= 0 ) dR2 = Deltar( axis1_eta, axis1_phi, Gen_neu2_eta, Gen_neu2_phi );
    dR = dR1;
    if ( dR2 < dR1 )
    { // make sure that the reco axis defined matches well with the axis of the gen neutralino, if not it is swapped
      iLLPrec1 = 2;
      iLLPrec2 = 1;
      dR = dR2;
    }
    float axis1_dR = dR;
    float axis2_eta = vaxis2.Eta();
    float axis2_phi = vaxis2.Phi();
    if ( njet2 == 0 )
    {  // compute an axis 2 even without jet, by taking the opposite in phi to axis 1
      axis2_eta = axis1_eta;
      axis2_phi = axis1_phi - 3.14159;
      if ( axis1_phi < 0 ) axis2_phi = axis1_phi + 3.14159;
    }
    if ( iLLPrec2 == 1 ) dR = Deltar( axis2_eta, axis2_phi, Gen_neu1_eta, Gen_neu1_phi );
    else                 dR = Deltar( axis2_eta, axis2_phi, Gen_neu2_eta, Gen_neu2_phi );
    float axis2_dR = dR;

    float dR_axis12 = Deltar(axis1_eta,axis1_phi,axis2_eta,axis2_phi);

    if ( showlog ) {
        cout << " njet1 " << njet1 << " and njet2" << njet2 << endl;
        cout << " axis1_eta " << axis1_eta << " and axis2_eta" << axis2_eta << endl;
        cout << " axis1_phi " << axis1_phi << " and axis2_phi" << axis2_phi << endl;
        cout << " axis1_dR " << axis1_dR << " and axis2_dR" << axis2_dR << endl;
        cout << " dR_axis12 " << dR_axis12 << endl;
    }
    

    //////////////////
    //------------------------------
    //selection of displaced tracks
    //------------------------------
    //////////////////
    
    if (showlog) cout << "//////////////////////////"<< endl;
    if (showlog) cout << "start to select displaced tracks " << endl;

    vector<reco::TransientTrack> displacedTracks_llp1_mva, displacedTracks_llp2_mva; // Control Tracks
    vector<reco::TransientTrack> displacedTracks_Hemi1_mva, displacedTracks_Hemi2_mva; // Tracks selected wrt the hemisphere
    if (showlog) cout << "||||||||||||||||||||||||| " << endl;
    ///// MVA for track selection coming from displaced tops
    
    //ajoute par Paul /*!*/
    float drSig, isinjet;
    // float dptSig; /*!*/
    int jet; /*!*/
    float ntrk10=-1, ntrk20=-1, ntrk30=-1; /*!*/
    float firsthit_X, firsthit_Y, firsthit_Z, dxy, dxyError, pt, eta,phi, NChi2, nhits;
    float algo;
    // float  track_dR;
    // float track_dRmax ; /*!*/
    double bdtval = -100.;

//$$
    LLP1_nTrks = 0;
    LLP2_nTrks = 0;

    int nTrks_axis1 = 0, nTrks_axis1_sig=0, nTrks_axis1_bad=0;
    int nTrks_axis2 = 0, nTrks_axis2_sig=0, nTrks_axis2_bad=0;
    int nTrks_axis1_sig_mva=0, nTrks_axis1_bad_mva=0;
    int nTrks_axis2_sig_mva=0, nTrks_axis2_bad_mva=0;
//$$
    
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

    // reader->AddVariable( "mva_track_firstHit_x", &firsthit_X );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_y", &firsthit_Y );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_z", &firsthit_Z );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_dxy", &dxy );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_dxyError", &dxyError );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_dz", &dz );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_dzError", &dzError );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    reader->AddVariable( "mva_track_pt", &pt );
    reader->AddVariable( "mva_track_eta", &eta );
    reader->AddVariable( "mva_track_nchi2", &NChi2 );
    reader->AddVariable( "mva_track_nhits", &nhits );
    reader->AddVariable( "mva_track_algo", &algo);
    //add the variables from my BDT(Paul)
    reader->AddVariable( "mva_ntrk10", &ntrk10);
    reader->AddVariable( "mva_drSig", &drSig); /*!*/
    reader->AddVariable( "mva_track_isinjet", &isinjet); /*!*/
    //reader->AddVariable("mva_track_dR",&track_dR);
    // reader->AddVariable("mva_track_dRmax",&track_dRmax);

    reader->BookMVA( "BDTG", weightFile_ );
    
//$$
    float pt_Cut = 1.;
    float NChi2_Cut = 5.;
    float drSig_Cut = 5.;
    float bdtcut = -0.0401; // optimal value w/o track association to axis: -0.0401
//$$

    int counter_track = -1;
    std::cout<< "size of tracks "<< trackRefs.size() << "pls help m edear god of physics" << std::endl;

    for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack)  // Loop on all the tracks
    {
      std::cout<<"b4 details"<<std::endl;
      counter_track++;
      const auto& itTrack = trackRefs[iTrack];
       if (itTrack->hasTrackDetails())
        {
          std::cout<<"after details"<<std::endl;
      firsthit_X = tree_track_firsthit_X[counter_track];
      firsthit_Y = tree_track_firsthit_Y[counter_track];
      firsthit_Z = tree_track_firsthit_Z[counter_track];
      std::cout<<"after details1"<<std::endl;
      dxy	 = tree_track_dxy[counter_track];
      dxyError   = tree_track_dxyError[counter_track];
      pt	 = tree_track_pt[counter_track];
      // ptError    = tree_track_ptError[counter_track];
      eta	 = tree_track_eta[counter_track];
      phi	 = tree_track_phi[counter_track];
      std::cout<<"after details2bis"<<std::endl;
      NChi2	 = tree_track_NChi2[counter_track];
      nhits	 = tree_track_nhits[counter_track];
      std::cout<<"after details2"<<std::endl;
      // algo	 = tree_track_algo[counter_track];
      algo=-1;
      // float dptSig=-1;
      // if (pt>0) dptSig=ptError/pt;
      
      //Ajoute par Paul /*!*/
      drSig = -1.;
      if ( dxyError > 0 ) drSig = abs(dxy) / dxyError; /*!*/
//       tree_track_drSig.push_back(drSig);
      ntrk10 = 0;
      ntrk20 = 0;
      ntrk30 = 0;
      bdtval = -10.;
      dR = -1.;
      std::cout<<"after details3"<<std::endl;
      int tracks_axis = 0; // flag to check which axis is the closest from the track

//$$
      if ( pt > pt_Cut && NChi2 < NChi2_Cut && drSig > drSig_Cut ) // preselection : pt > 1. && NChi2 < 5. && drSig > 5.
//$$
      { // On regarde les track_selec[i] qui sont True donc de potentielles tracks secondaires
        
        jet = tree_track_recoAK4SlimmedJet_idx[counter_track]; /*!*/
        isinjet = 0.;
        if ( jet >= 0 ) isinjet = 1.; /*!*/
        std::cout<<"after details4"<<std::endl;
        // int isFromLLP    = tree_track_sim_isFromLLP[counter_track]; // sim_isFromLLP is done way above that
        // int isFromBC_LLP = tree_track_sim_isFromBC_LLP[counter_track];
        int isFromLLP = 1;
        int isFromBC_LLP = 1;
        // float bdtval=tree_track_MVAval->at(i);//bdt cut is done below
        //check the dR between the tracks and the second axis (without any selection on the tracks)
        float dR1  = Deltar( eta, phi, axis1_eta, axis1_phi ); // axis1_phi and axis1_eta for the first axis
        float dR2  = Deltar( eta, phi, axis2_eta, axis2_phi );

        tracks_axis = 1;
        dR = dR1;
        if ( dR2 < dR1 ) { // a restriction could be added on the value of dR to assign the value Tracks_axis  (avoid some background???)
          tracks_axis = 2;
          dR = dR2;
        }

        //Computation of the distances needed for the BDT
        int counter_othertrack = -1; ///*!*/
        for (size_t iTrack2=0; iTrack2<trackRefs.size(); ++iTrack2)    // Loop on all the other Tracks/*!*/
        {
          counter_othertrack++;
                 if (itTrack->hasTrackDetails())
        {
          if ( counter_othertrack == counter_track ) continue;
          float pt2  = tree_track_pt[counter_othertrack];
          float dr2 = abs(tree_track_dxy[counter_othertrack]);
          float drSig2 = -1.;
          if ( tree_track_dxyError[counter_othertrack] > 0 ) drSig2 = dr2 / tree_track_dxyError[counter_othertrack];
          NChi2 = tree_track_NChi2[counter_othertrack];
          std::cout<<"b4 slection cut on other tracks"<<std::endl;
        if ( !(pt2 > pt_Cut && NChi2 < NChi2_Cut && drSig2 > drSig_Cut ) ) continue; // On regarde les autres track_selec[i] qui sont True donc de potnetielles tracks secondaires
//$$
          std::cout<<"after slection cut on other tracks"<<std::endl;
          float x2 = tree_track_firsthit_X[counter_othertrack];
          float y2 = tree_track_firsthit_Y[counter_othertrack];
          float z2 = tree_track_firsthit_Z[counter_othertrack];
          float dist = TMath::Sqrt( (firsthit_X-x2)*(firsthit_X-x2) + (firsthit_Y-y2)*(firsthit_Y-y2) + (firsthit_Z-z2)*(firsthit_Z-z2) );//pour chaque reconstruite, on regarde les autres tracks,
          if ( dist < 10. )	     {ntrk10++;} // les sctocker les 3 , on teste sur une seule couche quand on regarde vers l'avant
          if ( dist < 20. )	     {ntrk20++;}
          if ( dist < 30. )	     {ntrk30++;}
        }
        }  // end Loop on other Tracks

        if ( dR < dRcut_tracks ) 
	{
	  if ( isFromLLP == 1 || isFromBC_LLP == 1 ) LLP1_nTrks++;
	  if ( isFromLLP == 2 || isFromBC_LLP == 2 ) LLP2_nTrks++;
	
          bdtval = reader->EvaluateMVA( "BDTG" ); //default value = -10 (no -10 observed and -999 comes from EvaluateMVA)
          //cout << "BDT VAL " << bdtval <<endl;

          if ( tracks_axis == 1 ) {
	    nTrks_axis1++;
            if ( isFromLLP == iLLPrec1 || isFromBC_LLP == iLLPrec1 ) nTrks_axis1_sig++;
            else if ( isFromLLP >= 1 || isFromBC_LLP >= 1 )          nTrks_axis1_bad++;
          }
        
          if ( tracks_axis == 2 ) {
	    nTrks_axis2++;
            if ( isFromLLP == iLLPrec2 || isFromBC_LLP == iLLPrec2 ) nTrks_axis2_sig++;
            else if ( isFromLLP >= 1 || isFromBC_LLP >= 1 )          nTrks_axis2_bad++;
          }
        
          // if (bdtval < 0.0674  ) { //optimal cut for 10 cm, trained on 19k events, <s>=22 and <b>=240 //J�r�my
          // if (bdtval < 0.07381 ) { //optimal cut for 50 cm, trained on 10k events, <s>=15 and <b>=234 //J�r�my
          if ( bdtval > bdtcut )      //optimal cut for 50 cm ,trained on 10k events //Paul, expected to change
          { 
            bool Detailed =  itTrack->hasTrackDetails() ;
            std::cout<<"track is detailed : "<<Detailed<<std::endl;
                  if (itTrack->hasTrackDetails())
        {
            ////--------------Control tracks-----------------////
            if ( isFromLLP == 1 || isFromBC_LLP == 1 )
            {
              // displacedTracks_llp1_mva.push_back(theTransientTrackBuilder->build(*itTrack));
              displacedTracks_llp1_mva.push_back(theTransientTrackBuilder->build(itTrack->pseudoTrack()));
            }
            if ( isFromLLP == 2 || isFromBC_LLP == 2 )
            {
              // displacedTracks_llp2_mva.push_back(theTransientTrackBuilder->build(*itTrack));
              displacedTracks_llp2_mva.push_back(theTransientTrackBuilder->build(itTrack->pseudoTrack()));
            }
          
            if ( tracks_axis == 1 )
            {
              // displacedTracks_Hemi1_mva.push_back(theTransientTrackBuilder->build(*itTrack));
              displacedTracks_Hemi1_mva.push_back(theTransientTrackBuilder->build(itTrack->pseudoTrack()));
              if ( isFromLLP == iLLPrec1 || isFromBC_LLP == iLLPrec1 ) nTrks_axis1_sig_mva++;
              else if ( isFromLLP >= 1 || isFromBC_LLP >= 1 )          nTrks_axis1_bad_mva++;
            }
 
            if ( tracks_axis == 2 )
            {
              // displacedTracks_Hemi2_mva.push_back(theTransientTrackBuilder->build(*itTrack));
              displacedTracks_Hemi2_mva.push_back(theTransientTrackBuilder->build(itTrack->pseudoTrack()));
              if ( isFromLLP == iLLPrec2 || isFromBC_LLP == iLLPrec2 ) nTrks_axis2_sig_mva++;
              else if ( isFromLLP >= 1 || isFromBC_LLP >= 1 )          nTrks_axis2_bad_mva++;
            }
        }
          }
        }
      }
      
      tree_track_ntrk10.push_back(ntrk10);
      tree_track_ntrk20.push_back(ntrk20);
      tree_track_ntrk30.push_back(ntrk30);
      tree_track_MVAval.push_back(bdtval);
//$$
      tree_track_Hemi.push_back(tracks_axis);
      tree_track_Hemi_dR.push_back(dR);
      if      ( tracks_axis == 1 ) tree_track_Hemi_LLP.push_back(iLLPrec1);
      else if ( tracks_axis == 2 ) tree_track_Hemi_LLP.push_back(iLLPrec2);
      else		           tree_track_Hemi_LLP.push_back(0);
//$$
        }
    } //End loop on all the tracks
    
    if ( showlog ) {
        cout << " displaced tracks LLP1 " << LLP1_nTrks << " and with mva" << displacedTracks_llp1_mva.size() << endl;
        cout << " displaced tracks LLP2 " << LLP2_nTrks << " and with mva" << displacedTracks_llp2_mva.size() << endl;
        cout << " displaced tracks Hemi1 " << nTrks_axis1 << " and with mva" << displacedTracks_Hemi1_mva.size() << endl;
        cout << " displaced tracks Hemi2 " << nTrks_axis2 << " and with mva" << displacedTracks_Hemi2_mva.size() << endl;
    }
    
    //---------------------------------------------------------------------------------------//
    
    
    //////////////////////////////////
    //--------------------------------
    // Vertex fitting 
    //--------------------------------
    //////////////////////////////////
    
    int   Vtx_ntk = 0;
    float Vtx_x = 0., Vtx_y = 0., Vtx_z= 0., Vtx_chi = 0.;
    
//$$
    // parameters for the Adaptive Vertex Fitter (AVF)
    double maxshift        = 0.0001;
    unsigned int maxstep   = 30;
    double maxlpshift      = 0.1;
    double weightThreshold = 0.001;
    double sigmacut        = 5.;
    double Tini            = 256.;
    double ratio           = 0.25;
//$$
    

    //------------------------------- FIRST LLP WITH MVA ----------------------------------//
    
//$$    KalmanVertexFitter theFitter_vertex_llp1_mva(kvfPSet); //One or less vertex: either Valid or NotValid /*!*/
//$$
    static AdaptiveVertexFitter 
    theFitter_vertex_llp1_mva(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
    theFitter_vertex_llp1_mva.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
//$$
    
    Vtx_ntk = displacedTracks_llp1_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -1.;
    
    if ( Vtx_ntk > 1 )
    {
      TransientVertex displacedVertex_llp1_mva = theFitter_vertex_llp1_mva.vertex(displacedTracks_llp1_mva); // fitted vertex
      
      // std::cout<< "displacedVertex_llp1_mva is built" << std::endl;
      
      if ( displacedVertex_llp1_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      {
        Vtx_ntk = displacedTracks_llp1_mva.size();
        Vtx_x = displacedVertex_llp1_mva.position().x();
        Vtx_y = displacedVertex_llp1_mva.position().y();
        Vtx_z = displacedVertex_llp1_mva.position().z();
        Vtx_chi = displacedVertex_llp1_mva.normalisedChiSquared();
      }
    }

//$$
    tree_LLP.push_back(1);
    tree_LLP_pt.push_back(   LLP1_pt);
    tree_LLP_eta.push_back(  LLP1_eta);
    tree_LLP_phi.push_back(  LLP1_phi);
    tree_LLP_x.push_back(    LLP1_x);
    tree_LLP_y.push_back(    LLP1_y);
    tree_LLP_z.push_back(    LLP1_z);
    tree_LLP_nTrks.push_back(LLP1_nTrks);
    tree_LLP_Vtx_nTrks.push_back(Vtx_ntk);
    tree_LLP_Vtx_NChi2.push_back(Vtx_chi);
    tree_LLP_Vtx_dx.push_back(Vtx_x - LLP1_x);
    tree_LLP_Vtx_dy.push_back(Vtx_y - LLP1_y);
    tree_LLP_Vtx_dz.push_back(Vtx_z - LLP1_z);
//$$


    //-------------------------- SECOND LLP WITH MVA -------------------------------------//
    
//$$    KalmanVertexFitter theFitter_vertex_llp2_mva(kvfPSet); //One or less vertex: either Valid or NotValid /*!*/
//$$
    static AdaptiveVertexFitter 
    theFitter_vertex_llp2_mva(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
    theFitter_vertex_llp2_mva.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
//$$
    
    Vtx_ntk = displacedTracks_llp2_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -1.;
    
    if ( Vtx_ntk > 1 )
    {
      TransientVertex displacedVertex_llp2_mva = theFitter_vertex_llp2_mva.vertex(displacedTracks_llp2_mva); // fitted vertex
      
      if ( displacedVertex_llp2_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      {
        Vtx_x = displacedVertex_llp2_mva.position().x();
        Vtx_y = displacedVertex_llp2_mva.position().y();
        Vtx_z = displacedVertex_llp2_mva.position().z();
        Vtx_chi = displacedVertex_llp2_mva.normalisedChiSquared();
      }
    }

//$$
    tree_LLP.push_back(2);
    tree_LLP_pt.push_back(   LLP2_pt);
    tree_LLP_eta.push_back(  LLP2_eta);
    tree_LLP_phi.push_back(  LLP2_phi);
    tree_LLP_x.push_back(    LLP2_x);
    tree_LLP_y.push_back(    LLP2_y);
    tree_LLP_z.push_back(    LLP2_z);
    tree_LLP_nTrks.push_back(LLP2_nTrks);
    tree_LLP_Vtx_nTrks.push_back(Vtx_ntk);
    tree_LLP_Vtx_NChi2.push_back(Vtx_chi);
    tree_LLP_Vtx_dx.push_back(Vtx_x - LLP2_x);
    tree_LLP_Vtx_dy.push_back(Vtx_y - LLP2_y);
    tree_LLP_Vtx_dz.push_back(Vtx_z - LLP2_z);
//$$


    //--------------------------- FIRST HEMISPHERE WITH MVA -------------------------------------//
    
//$$    KalmanVertexFitter theFitter_Vertex_Hemi1_mva(kvfPSet); //One or less vertex: either Valid or NotValid /*!*/
//$$
    static AdaptiveVertexFitter 
    theFitter_Vertex_Hemi1_mva(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
    theFitter_Vertex_Hemi1_mva.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
//$$
    
    Vtx_ntk = displacedTracks_Hemi1_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -1.;
	
    if ( Vtx_ntk > 1 )
    {
      TransientVertex displacedVertex_Hemi1_mva = theFitter_Vertex_Hemi1_mva.vertex(displacedTracks_Hemi1_mva); // fitted vertex
      if ( displacedVertex_Hemi1_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      { 
        Vtx_x = displacedVertex_Hemi1_mva.position().x();
        Vtx_y = displacedVertex_Hemi1_mva.position().y();
        Vtx_z = displacedVertex_Hemi1_mva.position().z();
        Vtx_chi = displacedVertex_Hemi1_mva.normalisedChiSquared();
        Vtx_ntk = displacedTracks_Hemi1_mva.size();
      }
    }
   
//$$
    float Vtx_chi1 = Vtx_chi;
    tree_Hemi.push_back(1);
    tree_Hemi_njet.push_back(njet1);
    tree_Hemi_eta.push_back(axis1_eta);
    tree_Hemi_phi.push_back(axis1_phi);
    tree_Hemi_dR.push_back(axis1_dR);
    tree_Hemi_nTrks.push_back(nTrks_axis1);
    tree_Hemi_nTrks_sig.push_back(nTrks_axis1_sig);
    tree_Hemi_nTrks_bad.push_back(nTrks_axis1_bad);
    tree_Hemi_Vtx_NChi2.push_back(Vtx_chi);
    tree_Hemi_Vtx_nTrks.push_back(Vtx_ntk);
    tree_Hemi_Vtx_nTrks_sig.push_back(nTrks_axis1_sig_mva);
    tree_Hemi_Vtx_nTrks_bad.push_back(nTrks_axis1_bad_mva);
    tree_Hemi_Vtx_x.push_back(Vtx_x);
    tree_Hemi_Vtx_y.push_back(Vtx_y);
    tree_Hemi_Vtx_z.push_back(Vtx_z);
    if ( iLLPrec1 == 1 ) {
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP1_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP1_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP1_z);
    }
    else {
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP2_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP2_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP2_z);
    }

    float dist = 0.;
    if ( iLLPrec1 == 1 ) {
      dist = (LLP1_x - tree_GenPVx)*(LLP1_x - tree_GenPVx)
           + (LLP1_y - tree_GenPVy)*(LLP1_y - tree_GenPVy)
           + (LLP1_z - tree_GenPVz)*(LLP1_z - tree_GenPVz);
      tree_Hemi_LLP_pt.push_back( LLP1_pt);
      tree_Hemi_LLP_eta.push_back(LLP1_eta);
      tree_Hemi_LLP_phi.push_back(LLP1_phi);
      tree_Hemi_LLP_x.push_back(LLP1_x);
      tree_Hemi_LLP_y.push_back(LLP1_y);
      tree_Hemi_LLP_z.push_back(LLP1_z);
    }
    else {
      dist = (LLP2_x - tree_GenPVx)*(LLP2_x - tree_GenPVx)
           + (LLP2_y - tree_GenPVy)*(LLP2_y - tree_GenPVy)
           + (LLP2_z - tree_GenPVz)*(LLP2_z - tree_GenPVz);
      tree_Hemi_LLP_pt.push_back( LLP2_pt);
      tree_Hemi_LLP_eta.push_back(LLP2_eta);
      tree_Hemi_LLP_phi.push_back(LLP2_phi);
      tree_Hemi_LLP_x.push_back(LLP2_x);
      tree_Hemi_LLP_y.push_back(LLP2_y);
      tree_Hemi_LLP_z.push_back(LLP2_z);
    }
    tree_Hemi_LLP.push_back(iLLPrec1);
    tree_Hemi_LLP_dist.push_back(TMath::Sqrt(dist));
//$$
    

    //--------------------------- SECOND HEMISPHERE WITH MVA -------------------------------------//
    
//$$    KalmanVertexFitter theFitter_Vertex_Hemi2_mva(kvfPSet); //One or less vertex: either Valid or NotValid /*!*/
//$$
    static AdaptiveVertexFitter 
    theFitter_Vertex_Hemi2_mva(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
    theFitter_Vertex_Hemi2_mva.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
//$$
    
    Vtx_ntk = displacedTracks_Hemi2_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -1.;
    
    if ( Vtx_ntk > 1 )
    {
      TransientVertex displacedVertex_Hemi2_mva = theFitter_Vertex_Hemi2_mva.vertex(displacedTracks_Hemi2_mva); // fitted vertex
      
      if ( displacedVertex_Hemi2_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      {
        Vtx_x = displacedVertex_Hemi2_mva.position().x();
        Vtx_y = displacedVertex_Hemi2_mva.position().y();
        Vtx_z = displacedVertex_Hemi2_mva.position().z();
        Vtx_chi = displacedVertex_Hemi2_mva.normalisedChiSquared();
      }
    }
    
//$$
    float Vtx_chi2 = Vtx_chi;
    tree_Hemi.push_back(2);
    tree_Hemi_njet.push_back(njet2);
    tree_Hemi_eta.push_back(axis2_eta);
    tree_Hemi_phi.push_back(axis2_phi);
    tree_Hemi_dR.push_back(axis2_dR);
    tree_Hemi_nTrks.push_back(nTrks_axis2);
    tree_Hemi_nTrks_sig.push_back(nTrks_axis2_sig);
    tree_Hemi_nTrks_bad.push_back(nTrks_axis2_bad);
    tree_Hemi_Vtx_NChi2.push_back(Vtx_chi);
    tree_Hemi_Vtx_nTrks.push_back(Vtx_ntk);
    tree_Hemi_Vtx_nTrks_sig.push_back(nTrks_axis2_sig_mva);
    tree_Hemi_Vtx_nTrks_bad.push_back(nTrks_axis2_bad_mva);
    tree_Hemi_Vtx_x.push_back(Vtx_x);
    tree_Hemi_Vtx_y.push_back(Vtx_y);
    tree_Hemi_Vtx_z.push_back(Vtx_z);
    if ( iLLPrec2 == 1 ) {
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP1_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP1_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP1_z);
    }
    else {
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP2_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP2_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP2_z);
    }

    dist = 0.;
    if ( iLLPrec2 == 1 ) {
      dist = (LLP1_x - tree_GenPVx)*(LLP1_x - tree_GenPVx)
           + (LLP1_y - tree_GenPVy)*(LLP1_y - tree_GenPVy)
           + (LLP1_z - tree_GenPVz)*(LLP1_z - tree_GenPVz);
      tree_Hemi_LLP_pt.push_back( LLP1_pt);
      tree_Hemi_LLP_eta.push_back(LLP1_eta);
      tree_Hemi_LLP_phi.push_back(LLP1_phi);
      tree_Hemi_LLP_x.push_back(LLP1_x);
      tree_Hemi_LLP_y.push_back(LLP1_y);
      tree_Hemi_LLP_z.push_back(LLP1_z);
    }
    else {
      dist = (LLP2_x - tree_GenPVx)*(LLP2_x - tree_GenPVx)
           + (LLP2_y - tree_GenPVy)*(LLP2_y - tree_GenPVy)
           + (LLP2_z - tree_GenPVz)*(LLP2_z - tree_GenPVz);
      tree_Hemi_LLP_pt.push_back( LLP2_pt);
      tree_Hemi_LLP_eta.push_back(LLP2_eta);
      tree_Hemi_LLP_phi.push_back(LLP2_phi);
      tree_Hemi_LLP_x.push_back(LLP2_x);
      tree_Hemi_LLP_y.push_back(LLP2_y);
      tree_Hemi_LLP_z.push_back(LLP2_z);
    }
    tree_Hemi_LLP.push_back(iLLPrec2);
    tree_Hemi_LLP_dist.push_back(TMath::Sqrt(dist));
    
    counter_track = -1;
    for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack) { // Loop on all the tracks
      counter_track++;
      const auto& itTrack = trackRefs[iTrack];
      int hemi      = tree_track_Hemi[counter_track];
      double MVAval = tree_track_MVAval[counter_track];
      Vtx_chi = -1000.;
      if      ( hemi == 1 && MVAval > bdtcut ) Vtx_chi = Vtx_chi1;
      else if ( hemi == 2 && MVAval > bdtcut ) Vtx_chi = Vtx_chi2;
      tree_track_Hemi_mva_NChi2.push_back(Vtx_chi);
    } //End loop on all the tracks
    
    tree_Hemi_dR12.push_back(dR_axis12);
    tree_Hemi_dR12.push_back(dR_axis12);
    tree_Hemi_LLP_dR12.push_back(dRneuneu);
    tree_Hemi_LLP_dR12.push_back(dRneuneu);
//$$

    
//$$
  } // endif ZMu and HT selections
//$$
  smalltree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
TrackingPerf::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TrackingPerf::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackingPerf::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
    
    //Specify that only 'tracks' is allowed
    //To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(TrackingPerf);

void TrackingPerf::clearVariables() {
    
    tree_trigger_names.clear();
    tree_trigger_bits.clear();
    
    //-----------------------
    //fill the tree per track
    //tree_track_nclusters.clear();
    tree_track_charge.clear();
    tree_track_pt.clear();
    tree_track_outerPt.clear();
    tree_track_eta.clear();
    tree_track_phi.clear();
    tree_track_charge.clear();
    tree_track_nhits.clear();
    tree_track_NChi2.clear();
    tree_track_isHighQuality.clear();
    tree_track_isLoose.clear();
    tree_track_isTight.clear();
    
    tree_track_dxy.clear();
    tree_track_dxyError.clear();
    tree_track_dz.clear();
    tree_track_dzError.clear();
    tree_track_numberOfLostHits.clear();
    tree_track_numberOfValidHits.clear();
    tree_track_originalAlgo.clear();
    tree_track_algo.clear();
    // tree_track_stopReason.clear();
    //   tree_track_nPixel.clear();
    //   tree_track_nStrip.clear();
    tree_track_numberOfValidPixelHits.clear();
    tree_track_numberOfValidStripHits.clear();
    tree_track_numberOfValidStripTIBHits.clear();
    tree_track_numberOfValidStripTIDHits.clear();
    tree_track_numberOfValidStripTOBHits.clear();
    tree_track_numberOfValidStripTECHits.clear();
    tree_track_numberOfValidPixelBarrelHits.clear();
    tree_track_numberOfValidPixelEndcapHits.clear();
    tree_track_hasValidHitInPixelLayer.clear();
    tree_track_trackerLayersWithMeasurement.clear();
    tree_track_pixelLayersWithMeasurement.clear();
    tree_track_stripTECLayersWithMeasurement .clear();
    tree_track_stripTIBLayersWithMeasurement.clear();
    tree_track_stripTIDLayersWithMeasurement.clear();
    tree_track_stripTOBLayersWithMeasurement.clear();
    
    tree_track_vx.clear();
    tree_track_vy.clear();
    tree_track_vz.clear();
    tree_track_firsthit_X.clear();
    tree_track_firsthit_Y.clear();
    tree_track_firsthit_Z.clear();
    // tree_track_firsthit_phi.clear();
    tree_track_ntrk10.clear();
    tree_track_ntrk20.clear();
    tree_track_ntrk30.clear();
    
    tree_track_MVAval.clear();
//$$
    tree_track_Hemi.clear();
    tree_track_Hemi_dR.clear();
    tree_track_Hemi_mva_NChi2.clear();
    tree_track_Hemi_LLP.clear();
//$$
    
    tree_track_recoVertex_idx.clear();
    tree_track_recoAK4SlimmedJet_idx.clear();
    tree_track_recoAK4PFJet_idx.clear();
    tree_track_reco08Jet_idx.clear();
    tree_track_recoCaloJet_idx.clear();
    //    tree_track_reco08CaloJet_idx.clear();
    
    tree_track_nSimHits.clear();
    tree_track_isSimMatched.clear();
    
    tree_track_sim_charge.clear();
    tree_track_sim_pt.clear();
    tree_track_sim_eta.clear();
    tree_track_sim_phi.clear();
    tree_track_sim_longLived.clear();
    //tree_track_sim_matchedHit .clear();
    tree_track_sim_pdgId.clear();
    tree_track_sim_numberOfTrackerHits.clear();
    tree_track_sim_numberOfTrackerLayers.clear();
    tree_track_sim_mass.clear();
    tree_track_sim_status.clear();
    
    tree_track_sim_vx.clear();
    tree_track_sim_vy.clear();
    tree_track_sim_vz.clear();
    tree_track_sim_isFromLLP.clear();
//$$
    tree_track_sim_isFromBC.clear();
    tree_track_sim_isFromBC_mother_pdgId.clear();
    tree_track_sim_isFromBC_LLP.clear();
    tree_track_sim_dV.clear();
    tree_track_sim_dist.clear();
//$$
    
    //  tree_simtrack_charge.clear();
    //  tree_simtrack_pt.clear();
    //  tree_simtrack_eta.clear();
    //  tree_simtrack_phi.clear();
    //  tree_simtrack_longLived.clear();
    
    //  tree_simtrack_pdgId.clear();
    //  tree_simtrack_numberOfTrackerHits.clear();
    //  tree_simtrack_numberOfTrackerLayers.clear();
    //  tree_simtrack_mass.clear();
    //  tree_simtrack_status.clear();
    
    //  tree_simtrack_vx.clear();
    //  tree_simtrack_vy.clear();
    //  tree_simtrack_vz.clear();
    
    //  tree_simtrack_isRecoMatched.clear();
    //  tree_simtrack_pca_dxy.clear();
    //  tree_simtrack_pca_dz.clear();
    //  tree_simtrack_trkIdx.clear();
    
    tree_vtx_PosX.clear();
    tree_vtx_PosY.clear();
    tree_vtx_PosZ.clear();
    tree_vtx_NChi2.clear();
    
    tree_vtx_PosXError.clear();
    tree_vtx_PosYError.clear();
    tree_vtx_PosZError.clear();
    
    tree_AK4Slimmedjet_E.clear();
    tree_AK4Slimmedjet_pt.clear();
    tree_AK4Slimmedjet_eta.clear();
    tree_AK4Slimmedjet_phi.clear();
//$$    tree_AK4Slimmedjet_idxTrack.clear();
    
    tree_AK4PFjet_E.clear();
    tree_AK4PFjet_pt.clear();
    tree_AK4PFjet_eta.clear();
    tree_AK4PFjet_phi.clear();
//$$    tree_AK4PFjet_idxTrack.clear();
    
    tree_CaloJet_E.clear();
    tree_CaloJet_pt.clear();
    tree_CaloJet_eta.clear();
    tree_CaloJet_phi.clear();
//$$    tree_CaloJet_idxTrack.clear();
    
    tree_jet08_E.clear();
    tree_jet08_pt.clear();
    tree_jet08_eta.clear();
    tree_jet08_phi.clear();
//$$    tree_jet08_idxTrack.clear();
    
    tree_genParticle_pt.clear();
    tree_genParticle_eta.clear();
    tree_genParticle_phi.clear();
    tree_genParticle_charge.clear();
    tree_genParticle_pdgId.clear();
    tree_genParticle_vx.clear();
    tree_genParticle_vy.clear();
    tree_genParticle_vz.clear();
    tree_genParticle_mass.clear();
    tree_genParticle_statusCode.clear();
    tree_genParticle_mother_pdgId.clear();

//$$
    tree_genFromC_pt.clear();
    tree_genFromC_eta.clear();
    tree_genFromC_phi.clear();
    tree_genFromC_charge.clear();
    tree_genFromC_pdgId.clear();
    tree_genFromC_vx.clear();
    tree_genFromC_vy.clear();
    tree_genFromC_vz.clear();
    tree_genFromC_mother_pdgId.clear();
    tree_genFromC_generation.clear();
    tree_genFromC_LLP.clear();

    tree_genFromB_pt.clear();
    tree_genFromB_eta.clear();
    tree_genFromB_phi.clear();
    tree_genFromB_charge.clear();
    tree_genFromB_pdgId.clear();
    tree_genFromB_vx.clear();
    tree_genFromB_vy.clear();
    tree_genFromB_vz.clear();
    tree_genFromB_mother_pdgId.clear();
    tree_genFromB_generation.clear();
    tree_genFromB_LLP.clear();
//$$

    tree_genJet_pt.clear();
    tree_genJet_eta.clear();
    tree_genJet_phi.clear();
    tree_genJet_mass.clear();
    tree_genJet_energy.clear();
    
    tree_ak8GenJet_pt.clear();
    tree_ak8GenJet_eta.clear();
    tree_ak8GenJet_phi.clear();
    tree_ak8GenJet_mass.clear();
    tree_ak8GenJet_energy.clear();
    
    tree_electron_pt.clear();
    tree_electron_eta.clear();
    tree_electron_phi.clear();
    tree_electron_vx.clear();
    tree_electron_vy.clear();
    tree_electron_vz.clear();
    tree_electron_energy.clear();
    tree_electron_charge.clear();
    
    tree_slimmedmuon_pt.clear();
    tree_slimmedmuon_eta.clear();
    tree_slimmedmuon_phi.clear();
    tree_slimmedmuon_vx.clear();
    tree_slimmedmuon_vy.clear();
    tree_slimmedmuon_vz.clear();
    tree_slimmedmuon_energy.clear();
    tree_slimmedmuon_dxy.clear();
    tree_slimmedmuon_dxyError.clear();
    tree_slimmedmuon_dz.clear();
    tree_slimmedmuon_dzError.clear();
    tree_slimmedmuon_charge.clear();
    tree_slimmedmuon_PFisoVeryTight.clear();
    tree_slimmedmuon_PFisoTight.clear();
    tree_slimmedmuon_PFisoMedium.clear();
    tree_slimmedmuon_PFisoLoose.clear();
    tree_slimmedmuon_MVAisoLoose.clear();
    tree_slimmedmuon_MVAisoMedium.clear();
    tree_slimmedmuon_MVAisoTight.clear();
    tree_slimmedmuon_isGlobalMuon.clear();
    tree_slimmedmuon_isStandAloneMuon.clear();
    tree_slimmedmuon_CutBasedIdLoose.clear();
    tree_slimmedmuon_CutBasedIdMedium.clear();
    tree_slimmedmuon_CutBasedIdMediumPrompt.clear();
    tree_slimmedmuon_CutBasedIdTight.clear();
    
//$$
    tree_LLP.clear();
    tree_LLP_pt.clear();
    tree_LLP_eta.clear();
    tree_LLP_phi.clear();
    tree_LLP_x.clear();
    tree_LLP_y.clear();
    tree_LLP_z.clear();
    tree_LLP_nTrks.clear();
    tree_LLP_Vtx_NChi2.clear();
    tree_LLP_Vtx_nTrks.clear();
    tree_LLP_Vtx_dx.clear();
    tree_LLP_Vtx_dy.clear();
    tree_LLP_Vtx_dz.clear();

    tree_Hemi.clear();
    tree_Hemi_njet.clear();
    tree_Hemi_eta.clear();
    tree_Hemi_phi.clear();
    tree_Hemi_dR.clear();
    tree_Hemi_nTrks.clear();
    tree_Hemi_nTrks_sig.clear();
    tree_Hemi_nTrks_bad.clear();
    tree_Hemi_LLP.clear();
    tree_Hemi_LLP_pt.clear();
    tree_Hemi_LLP_eta.clear();
    tree_Hemi_LLP_phi.clear();
    tree_Hemi_LLP_dist.clear();
    tree_Hemi_LLP_x.clear();
    tree_Hemi_LLP_y.clear();
    tree_Hemi_LLP_z.clear();
    tree_Hemi_Vtx_NChi2.clear();
    tree_Hemi_Vtx_nTrks.clear();
    tree_Hemi_Vtx_nTrks_sig.clear();
    tree_Hemi_Vtx_nTrks_bad.clear();
    tree_Hemi_Vtx_x.clear();
    tree_Hemi_Vtx_y.clear();
    tree_Hemi_Vtx_z.clear();
    tree_Hemi_Vtx_dx.clear();
    tree_Hemi_Vtx_dy.clear();
    tree_Hemi_Vtx_dz.clear();
    tree_Hemi_dR12.clear();
    tree_Hemi_LLP_dR12.clear();
//$$
}


void TrackingPerf::getHLTPathNames(const edm::TriggerNames& triggerNames ){
    for (auto const& pattern : filterTriggerNames_) {
        if (edm::is_glob(pattern)) {
            // found a glob pattern, expand it
            std::vector<std::vector<std::string>::const_iterator> matches =
            edm::regexMatch(triggerNames.triggerNames(), pattern);
            // store the matching patterns
            for (auto const& match : matches)  HLTPathsByName_.push_back(*match);
        } else {
            // found a trigger name, just copy it
            HLTPathsByName_.push_back(pattern);
        }
    }
}


void TrackingPerf::fillTriggerBits(std::vector<std::string>& HLTPathsByName_, edm::Handle<edm::TriggerResults>& triggerBits, const edm::TriggerNames &names ){
    for(unsigned int iName=0; iName < HLTPathsByName_.size(); iName++){
        for (unsigned int i = 0; i < triggerBits->size(); i++){
            
            std::string triggerName = names.triggerName(i);
            if(triggerName == HLTPathsByName_[iName]){
                tree_trigger_names.push_back(HLTPathsByName_[iName]);
                tree_trigger_bits.push_back(triggerBits->accept(i) ? true : false);
            }
        }
    }
}



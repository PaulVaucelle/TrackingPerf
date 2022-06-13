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

// user include files

#include <memory> 

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Utilities/interface/RegexMatch.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/InvalidTrackingRecHit.h"

#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"

#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "RecoLocalTracker/SiStripClusterizer/interface/SiStripClusterInfo.h"

#include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"

#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>
#include <TLorentzVector.h>

//includes for track-sim association
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimTracker/TrackAssociation/interface/TrackingParticleIP.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h" 
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoTracker/FinalTrackSelectors/plugins/getBestVertex.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h" 

#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include <DataFormats/VertexReco/interface/Vertex.h>

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Common/interface/TriggerNames.h>

#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include "DataFormats/VertexReco/interface/Vertex.h"

#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"

#include <DataFormats/Math/interface/LorentzVector.h>
#include <DataFormats/VertexReco/interface/Vertex.h>

#include <DataFormats/Common/interface/MergeableCounter.h>
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "Geometry/Records/interface/HcalParametersRcd.h"
#include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "CondFormats/GeometryObjects/interface/HcalParameters.h"
#include "DetectorDescription/Core/interface/DDCompactView.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HcalCommonData/interface/HcalParametersFromDD.h"

#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "boost/functional/hash.hpp"

#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/Jet.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include <DataFormats/HepMCCandidate/interface/GenStatusFlags.h>

#include "TrackingPerf/TrackingPerf/interface/Proto.h"

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
  
  //$$
  std::string weightFile_;
  //$$
  
  //------------------------------------
  // track ( and event ) information 
  //------------------------------------
  const edm::EDGetTokenT<edm::View<reco::Track> > trackToken_;
  const edm::EDGetTokenT<edm::View<reco::Track> > trackSrc_;
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
  const edm::EDGetTokenT<edm::View<reco::Jet> > ak4slimmedJetToken_;
  const edm::EDGetTokenT<edm::View<reco::Jet> > ak4PFJetToken_;
  const edm::EDGetTokenT<edm::View<reco::Jet> > CaloJetToken_;
  
  const edm::EDGetTokenT<edm::View<reco::Jet> > ak8jetToken_;
  //   const edm::EDGetTokenT<edm::View<reco::Jet> > ak8CaloJetToken_;
  
  //------------------------------------
  // gen information
  //------------------------------------
  const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  const edm::EDGetTokenT<reco::GenJetCollection>      genJetToken_;
  const edm::EDGetTokenT<reco::GenJetCollection>      ak8GenJetToken_;
  const edm::EDGetTokenT<reco::GenJetCollection>      genTTXJetsToken_;
  const edm::EDGetTokenT<GenEventInfoProduct>	      genEventInfoToken_;
  const edm::EDGetTokenT<LHEEventProduct>	     LHEEventProductToken_;
  
  //pat information
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;
  
  std::string parametersDefinerName_;
  
  //------------------------------------
  // electrons
  //------------------------------------
  const edm::EDGetTokenT<pat::ElectronCollection> electronPATToken_;
  //------------------------------------
  // muons 
  //------------------------------------
  const edm::EDGetTokenT<pat::MuonCollection> slimmedmuonToken_;
  //$$  const edm::EDGetTokenT<reco::MuonCollection> recomuonToken_;
  //------------------------------------
  // MET
  //------------------------------------
  const edm::EDGetTokenT<pat::METCollection> metToken_;
  
  
  //------------------------------------
  // Zmu skim bit
  //------------------------------------
  std::vector<std::string> filterTriggerNames_;
  std::vector<std::string> HLTPathsByName_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  const edm::EDGetTokenT<vector<reco::CompositeCandidate> > ZmumuCandidate_;
  
  //additional 
  bool useCluster_;
  bool runOnData_;
  
  TTree *smalltree;
  
  edm::Service<TFileService> fs;
  
  std::string ttrhbuilder_; 
  
  edm::ESHandle<MagneticField> bField;
  
  edm::ParameterSet kvfPSet;
  
  int tree_NbrOfZCand;
  
  
  /////// VARIABLES FOR FILLING THE TREE 
  
  //-----------------------
  //trigger variable
  std::vector<string > tree_trigger_names;
  std::vector<bool >   tree_trigger_bits;
  
  std::vector<bool>    tree_passesHTFilter; 
  
  //-----------------------
  //fill the tree per track
  //std::vector<int> tree_track_nclusters;
  std::vector< float > tree_track_pt;
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
  std::vector<float>	tree_track_firsthit_X; 
  std::vector<float>	tree_track_firsthit_Y; 
  std::vector<float>	tree_track_firsthit_Z;
  std::vector<float>	tree_track_firsthit_phi; 
  std::vector<float>	tree_track_ntrk10;
  std::vector<float>	tree_track_ntrk20;
  std::vector<float>	tree_track_ntrk30;
  
  std::vector<double>	tree_track_MVAval; 
  
  std::vector<int>	tree_track_recoVertex_idx; 
  std::vector<int>	tree_track_recoAK4SlimmedJet_idx; 
  std::vector<int>	tree_track_recoAK4PFJet_idx; 
  std::vector<int>	tree_track_reco08Jet_idx; 
  std::vector<int>	tree_track_recoCaloJet_idx; 
  //   std::vector<int>      tree_track_reco08CaloJet_idx; 
  
  std::vector<int>      tree_track_nSimHits; 
  std::vector<bool>     tree_track_isSimMatched;

  std::vector< int >	tree_track_simtrack_charge;		   
  std::vector< float >  tree_track_simtrack_pt; 	   
  std::vector< float >  tree_track_simtrack_eta  ;	   
  std::vector< float >  tree_track_simtrack_phi  ;	   
  std::vector<bool>	tree_track_simtrack_longLived	      ;
  // std::vector<int>	   tree_track_simtrack_matchedHit    ;
  std::vector<int>	tree_track_simtrack_pdgId;	     
  std::vector<int>	tree_track_simtrack_numberOfTrackerHits  ;
  std::vector<int>	tree_track_simtrack_numberOfTrackerLayers;
  std::vector<float>	tree_track_simtrack_mass  ;		   
  std::vector<int>	tree_track_simtrack_status;		   

  std::vector<float>	tree_track_genVertexPos_X;	      
  std::vector<float>	tree_track_genVertexPos_Y;	      
  std::vector<float>	tree_track_genVertexPos_Z;	      
  std::vector<float>	tree_track_simtrack_llp1_dV;
  std::vector<float>	tree_track_simtrack_llp2_dV;
  std::vector<bool>	tree_track_simtrack_isFromDispTop;
  std::vector<int>	tree_track_simtrack_isFromLLP;
  
  
  int runNumber, eventNumber, lumiBlock; 
  
  //--------------------------------
  // Tracking Particle infos ------- 
  //--------------------------------
  
  std::vector< int >	       tree_simtrack_simtrack_charge;			  
  std::vector< float >         tree_simtrack_simtrack_pt;		  
  std::vector< float >         tree_simtrack_simtrack_eta;	  
  std::vector< float >         tree_simtrack_simtrack_phi;	  
  std::vector<bool>	       tree_simtrack_simtrack_longLived;
  std::vector<int>	       tree_simtrack_simtrack_pdgId;		   
  std::vector<int>	       tree_simtrack_simtrack_numberOfTrackerHits;
  std::vector<int>	       tree_simtrack_simtrack_numberOfTrackerLayers;
  std::vector<float>	       tree_simtrack_simtrack_mass;		  
  std::vector<int>	       tree_simtrack_simtrack_status;			  
  
  std::vector<float>	       tree_simtrack_genVertexPos_X;		    
  std::vector<float>	       tree_simtrack_genVertexPos_Y;		    
  std::vector<float>	       tree_simtrack_genVertexPos_Z;	    
  std::vector<bool>	       tree_simtrack_isRecoMatched;
  std::vector<float>	       tree_simtrack_pca_dxy;
  std::vector<float>	       tree_simtrack_pca_dz   ; 
  std::vector<std::vector<int> > tree_simtrack_trkIdx;  
  std::vector<float>	       tree_simtrack_isRecoMatched_pt;

  //--------------------------------
  // Beam spot infos ------- 
  //--------------------------------
  float 	  tree_bs_PosX;        
  float 	  tree_bs_PosY;        
  float 	  tree_bs_PosZ; 
  
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
  
  //secondary vertices 
  
  std::vector<float> tree_secondaryVtx_X; 
  std::vector<float> tree_secondaryVtx_Y; 
  std::vector<float> tree_secondaryVtx_Z; 
  
  std::vector<float> tree_secondaryVtx_diff_X; 
  std::vector<float> tree_secondaryVtx_diff_Y;
  std::vector<float> tree_secondaryVtx_diff_Z;
  std::vector<int>   tree_secondaryVtx_nTracks;
  std::vector<int>   tree_secondaryVtx_isValid; 
  std::vector<float> tree_secondaryVtx_NChi2; 
  
//   std::vector<float> tree_secondaryVtx_iterative_X; 
//   std::vector<float> tree_secondaryVtx_iterative_Y; 
//   std::vector<float> tree_secondaryVtx_iterative_Z; 
//   std::vector<int>   tree_secondaryVtx_iterative_nTracks;
//   std::vector<float> tree_secondaryVtx_iterative_NChi2; 
//   std::vector<bool>  tree_secondaryVtx_iterative_match_top;
//   std::vector<bool>  tree_secondaryVtx_iterative_isSelected; 

  std::vector<float> tree_genTop_X; 
  std::vector<float> tree_genTop_Y; 
  std::vector<float> tree_genTop_Z; 
  std::vector<float> tree_genTop_charge; 

//   std::vector<float> tree_seedVtx_X; 
//   std::vector<float> tree_seedVtx_Y; 
//   std::vector<float> tree_seedVtx_Z; 
//   std::vector<float> tree_seedVtx_dd; 
//   std::vector<float> tree_seedVtx_dphi; 
//   std::vector<float> tree_seedVtx_distance2track; 
//   std::vector<float> tree_seedVtx_normChi2; 

  //--------------------------------
  // jet infos ------- 
  //--------------------------------
  
  std::vector<float> tree_AK4Slimmedjet_E;
  std::vector<float> tree_AK4Slimmedjet_pt;
  std::vector<float> tree_AK4Slimmedjet_eta;
  std::vector<float> tree_AK4Slimmedjet_phi;
  std::vector<int>   tree_AK4Slimmedjet_idxTrack;
  
  std::vector<float> tree_AK4PFjet_E;
  std::vector<float> tree_AK4PFjet_pt;
  std::vector<float> tree_AK4PFjet_eta;
  std::vector<float> tree_AK4PFjet_phi;
  std::vector<int>   tree_AK4PFjet_idxTrack;
  
  std::vector<float> tree_CaloJet_E;
  std::vector<float> tree_CaloJet_pt;
  std::vector<float> tree_CaloJet_eta;
  std::vector<float> tree_CaloJet_phi;
  std::vector<int>   tree_CaloJet_idxTrack;
  
  std::vector<float> tree_jet08_E;
  std::vector<float> tree_jet08_pt;
  std::vector<float> tree_jet08_eta;
  std::vector<float> tree_jet08_phi;
  std::vector<int>   tree_jet08_idxTrack;
  
  //  std::vector<float> tree_CaloJet08_E;
  //  std::vector<float> tree_CaloJet08_pt;
  //  std::vector<float> tree_CaloJet08_eta;
  //  std::vector<float> tree_CaloJet08_phi;
  //  std::vector<int> tree_CaloJet08_idxTrack; 
  
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
  
  std::vector< float > tree_genParticle_charge;
  std::vector< float > tree_genParticle_pt;
  std::vector< float > tree_genParticle_eta;
  std::vector< float > tree_genParticle_phi;
  std::vector< float > tree_genParticle_pdgId;
  std::vector< float > tree_genParticle_vx;
  std::vector< float > tree_genParticle_vy;
  std::vector< float > tree_genParticle_vz;
  std::vector< float > tree_genParticle_mass;
  std::vector< int >   tree_genParticle_statusCode; 
  std::vector< int >   tree_genParticle_mother_pdgId;
  std::vector< float > tree_genParticle_mother_pt;
  std::vector< float > tree_genParticle_mother_eta;
  std::vector< float > tree_genParticle_mother_phi;
  
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
  std::vector<float> tree_slimmedmuon_charge; 
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

  //added by Paul 

  // about generated LLPs (top is used as the physics process first studied when looking for tops coming from neutralinos)
  int   tree_nLLP = 0;
  float tree_LLP1_pt;
  float tree_LLP1_eta;
  float tree_LLP1_phi;
  float tree_LLP1_x;
  float tree_LLP1_y;
  float tree_LLP1_z;

  float tree_LLP2_pt;
  float tree_LLP2_eta;
  float tree_LLP2_phi;
  float tree_LLP2_x;
  float tree_LLP2_y;
  float tree_LLP2_z;

  std::vector<int> tree_NSim;
  std::vector<float> tree_RecoEfficiency;

  std::vector<float> tree_FakeRate;
  std::vector<float> tree_RecoVertex1_Eff;
  std::vector<float> tree_RecoVertex2_Eff;
  std::vector<float> tree_VtxReco_Eff;
  
  int	tree_Vtx_LLP1_nTrks = 0;
  float	tree_Vtx_LLP1_x;
  float tree_Vtx_LLP1_y;
  float tree_Vtx_LLP1_z;
  float	tree_Vtx_LLP1_NChi2;

  int	tree_Vtx_LLP2_nTrks = 0;
  float	tree_Vtx_LLP2_x;
  float tree_Vtx_LLP2_y;
  float tree_Vtx_LLP2_z;
  float	tree_Vtx_LLP2_NChi2;

  int	tree_Vtx_mva_LLP1_nTrks = 0;
  float	tree_Vtx_mva_LLP1_x;
  float tree_Vtx_mva_LLP1_y;
  float tree_Vtx_mva_LLP1_z;
  float	tree_Vtx_mva_LLP1_NChi2;

  int	tree_Vtx_mva_LLP2_nTrks = 0;
  float	tree_Vtx_mva_LLP2_x;
  float tree_Vtx_mva_LLP2_y;
  float tree_Vtx_mva_LLP2_z;
  float	tree_Vtx_mva_LLP2_NChi2;

//     //First top analysis
// 
//   std::vector<float>		    tree_seedVtx_X_top1;
//   std::vector<float>		    tree_seedVtx_Y_top1;
//   std::vector<float>		    tree_seedVtx_Z_top1;
//   std::vector<float>		    tree_seedVtx_dd_top1;
//   std::vector<float>		    tree_seedVtx_dphi_top1;
//   std::vector<float>		    tree_seedVtx_distance2track_top1; 
//   std::vector<float>		    tree_seedVtx_normChi2_top1;
//   std::vector<float>		    tree_VtxReco_Eff_top1;
// 
//      //Second top analysis
// 
//   std::vector<float>		    tree_seedVtx_X_top2;
//   std::vector<float>		    tree_seedVtx_Y_top2;
//   std::vector<float>		    tree_seedVtx_Z_top2;
//   std::vector<float>		    tree_seedVtx_dd_top2;
//   std::vector<float>		    tree_seedVtx_dphi_top2;
//   std::vector<float>		    tree_seedVtx_distance2track_top2; 
//   std::vector<float>		    tree_seedVtx_normChi2_top2;
//   std::vector<float>		    tree_VtxReco_Eff_top2;

//   std::vector<unsigned int>	    tree_DVertex_top1_nTrks;
//   std::vector<unsigned int>	    tree_DVertex_top2_nTrks;
//   std::vector<unsigned int>	    tree_DVertex_nTrks;

  std::vector<float>		    tree_seedVtx_dr_top1;
  std::vector<float>		    tree_seedVtx_dz_top1;
  std::vector<float>		    tree_seedVtx_dr_top2;
  std::vector<float>		    tree_seedVtx_dz_top2;
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
   //$$
  weightFile_( iConfig.getUntrackedParameter<std::string>("weightFileMVA") ),
  //$$
  trackToken_(  	 consumes<edm::View<reco::Track> >(		   iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
  trackSrc_(		 consumes<edm::View<reco::Track> >(		   iConfig.getParameter<edm::InputTag>("trackLabel") )),
  trackAssociatorToken_( consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssociator"))),
  beamSpotToken_(	 consumes<reco::BeamSpot>(			   iConfig.getUntrackedParameter<edm::InputTag>("beamSpot"))),
  vertexToken_( 	 consumes<reco::VertexCollection>(		   iConfig.getUntrackedParameter<edm::InputTag>("vertices"))),
  ak4slimmedJetToken_(   consumes<edm::View<reco::Jet> >(		   iConfig.getParameter<edm::InputTag>("ak4slimmedJetInput"))),
  ak4PFJetToken_(	 consumes<edm::View<reco::Jet> >(		   iConfig.getParameter<edm::InputTag>("ak4PFJetInput"))),
  CaloJetToken_(	 consumes<edm::View<reco::Jet> >(		   iConfig.getParameter<edm::InputTag>("caloJetInput"))),
  ak8jetToken_( 	 consumes<edm::View<reco::Jet> >(		   iConfig.getParameter<edm::InputTag>("ak8jetInput"))),
  //   ak8CaloJetToken_(	 consumes<edm::View<reco::Jet> >(		   iConfig.getParameter<edm::InputTag>("ak8CaloJetInput"))),
  genParticlesToken_(    consumes<reco::GenParticleCollection>(            iConfig.getParameter<edm::InputTag>("genParticles"))),
  genJetToken_(          consumes<reco::GenJetCollection>(                 iConfig.getParameter<edm::InputTag>("genJetInput"))),
  ak8GenJetToken_(          consumes<reco::GenJetCollection>(              iConfig.getParameter<edm::InputTag>("ak8GenJetInput"))),
  genEventInfoToken_(    consumes<GenEventInfoProduct>(                    iConfig.getParameter<edm::InputTag>("genEventInfoInput"))),
  LHEEventProductToken_( consumes<LHEEventProduct>(                        iConfig.getParameter<edm::InputTag>("LHEEventProductInput"))),
  pfcandsToken_(         consumes<pat::PackedCandidateCollection>(         iConfig.getParameter<edm::InputTag>("pfcands"))),
  parametersDefinerName_(						   iConfig.getUntrackedParameter<std::string>("parametersDefiner")),
  electronPATToken_(     consumes<pat::ElectronCollection>(                iConfig.getParameter<edm::InputTag>("electronInput"))),
  slimmedmuonToken_(     consumes<pat::MuonCollection>(			   iConfig.getParameter<edm::InputTag>("slimmedmuonInput"))), 
  //$$  recomuonToken_(        consumes<reco::MuonCollection>(	           iConfig.getParameter<edm::InputTag>("recomuonInput"))), 
  metToken_(             consumes<pat::METCollection>(			   iConfig.getParameter<edm::InputTag>("metInput"))),
  filterTriggerNames_(                                                     iConfig.getUntrackedParameter<std::vector<std::string> >("filterTriggerNames")),
  triggerBits_(          consumes<edm::TriggerResults>(                    edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("HLT"))) ),
  ZmumuCandidate_(       consumes<vector<reco::CompositeCandidate>>(	   iConfig.getParameter<edm::InputTag>("Zmumucand"))),
  useCluster_(          					           iConfig.getUntrackedParameter<bool>("useCluster")),
  runOnData_ (        						           iConfig.getUntrackedParameter<bool>("runOnData")),
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
  
  // track
  
  smalltree->Branch("tree_track_pt",		    &tree_track_pt);
  smalltree->Branch("tree_track_outerPt",           &tree_track_outerPt ); 		      
  smalltree->Branch("tree_track_eta",  	            &tree_track_eta ); 	      
  smalltree->Branch("tree_track_phi",  	            &tree_track_phi ); 
  smalltree->Branch("tree_track_charge",            &tree_track_charge );  	      
  smalltree->Branch("tree_track_nhits",	            &tree_track_nhits);	      
  smalltree->Branch("tree_track_NChi2",	            &tree_track_NChi2);	      
  smalltree->Branch("tree_track_isHighQuality",     &tree_track_isHighQuality);
  smalltree->Branch("tree_track_isLoose",	    &tree_track_isLoose);
  smalltree->Branch("tree_track_isTight",	    &tree_track_isTight);	 
  smalltree->Branch("tree_track_dxy",  	            &tree_track_dxy );        
  smalltree->Branch("tree_track_dxyError",	    &tree_track_dxyError);	   
  smalltree->Branch("tree_track_dz",		    &tree_track_dz);	       
  smalltree->Branch("tree_track_dzError",	    &tree_track_dzError  );	   
  smalltree->Branch("tree_track_numberOfLostHits",  &tree_track_numberOfLostHits); 
  smalltree->Branch("tree_track_numberOfValidHits", &tree_track_numberOfValidHits); 
  smalltree->Branch("tree_track_originalAlgo",      &tree_track_originalAlgo);      
  smalltree->Branch("tree_track_algo", 	            &tree_track_algo); 
  smalltree->Branch("tree_track_stopReason",	    &tree_track_stopReason);   
  smalltree->Branch("tree_track_isSimMatched",      &tree_track_isSimMatched	); 
   
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
   
  //  smalltree->Branch("tree_track_nPixel",	   &tree_track_nPixel );
  //  smalltree->Branch("tree_track_nStrip",	   &tree_track_nStrip );

  smalltree->Branch("tree_track_vx",           &tree_track_vx ); 
  smalltree->Branch("tree_track_vy",           &tree_track_vy ); 
  smalltree->Branch("tree_track_vz",           &tree_track_vz ); 
  smalltree->Branch("tree_track_firsthit_X",   &tree_track_firsthit_X); 
  smalltree->Branch("tree_track_firsthit_Y",   &tree_track_firsthit_Y);
  smalltree->Branch("tree_track_firsthit_Z",   &tree_track_firsthit_Z);
  smalltree->Branch("tree_track_firsthit_phi", &tree_track_firsthit_phi); 
  smalltree->Branch("tree_track_ntrk10",&tree_track_ntrk10);
  smalltree->Branch("tree_track_ntrk20",&tree_track_ntrk20);
  smalltree->Branch("tree_track_ntrk30",&tree_track_ntrk30);

  smalltree->Branch("tree_track_MVAval", &tree_track_MVAval); 

  smalltree->Branch("tree_track_recoVertex_idx", &tree_track_recoVertex_idx); 
  smalltree->Branch("tree_track_recoAK4SlimmedJet_idx", &tree_track_recoAK4SlimmedJet_idx); 
  smalltree->Branch("tree_track_recoAK4PFJet_idx",      &tree_track_recoAK4PFJet_idx); 
  smalltree->Branch("tree_track_reco08Jet_idx",         &tree_track_reco08Jet_idx); 
  smalltree->Branch("tree_track_recoCaloJet_idx",       &tree_track_recoCaloJet_idx); 
  //   smalltree->Branch("tree_track_reco08CaloJet_idx",     &tree_track_reco08CaloJet_idx); 
  
  // info about the simulated track matched to the reco track
   
  smalltree->Branch("tree_track_simtrack_charge",		 &tree_track_simtrack_charge );		      
  smalltree->Branch("tree_track_simtrack_pt",  		         &tree_track_simtrack_pt );	      
  smalltree->Branch("tree_track_simtrack_eta", 		         &tree_track_simtrack_eta  );  	      
  smalltree->Branch("tree_track_simtrack_phi", 		         &tree_track_simtrack_phi  );  	      
  smalltree->Branch("tree_track_simtrack_longLived",		 &tree_track_simtrack_longLived );
  smalltree->Branch("tree_track_simtrack_pdgId",		 &tree_track_simtrack_pdgId ); 
  smalltree->Branch("tree_track_simtrack_numberOfTrackerHits",   &tree_track_simtrack_numberOfTrackerHits   );
  smalltree->Branch("tree_track_simtrack_numberOfTrackerLayers", &tree_track_simtrack_numberOfTrackerLayers );
  smalltree->Branch("tree_track_simtrack_mass",		         &tree_track_simtrack_mass   );		      
  smalltree->Branch("tree_track_simtrack_status",		 &tree_track_simtrack_status );		      
 
  smalltree->Branch("tree_track_genVertexPos_X", 		 &tree_track_genVertexPos_X ); 		  
  smalltree->Branch("tree_track_genVertexPos_Y", 		 &tree_track_genVertexPos_Y ); 		  
  smalltree->Branch("tree_track_genVertexPos_Z", 		 &tree_track_genVertexPos_Z ); 		  
  smalltree->Branch("tree_track_simtrack_llp1_dV",               &tree_track_simtrack_llp1_dV );
  smalltree->Branch("tree_track_simtrack_llp2_dV",               &tree_track_simtrack_llp2_dV );
  smalltree->Branch("tree_track_simtrack_isFromLLP",	         &tree_track_simtrack_isFromLLP );
  smalltree->Branch("tree_track_simtrack_isFromDispTop",	 &tree_track_simtrack_isFromDispTop );
  
  smalltree->Branch("tree_genTop_X", &tree_genTop_X); 
  smalltree->Branch("tree_genTop_Y", &tree_genTop_Y); 
  smalltree->Branch("tree_genTop_Z", &tree_genTop_Z); 
  smalltree->Branch("tree_genTop_charge", &tree_genTop_charge); 

//   smalltree->Branch("tree_seedVtx_X", &tree_seedVtx_X); 
//   smalltree->Branch("tree_seedVtx_Y", &tree_seedVtx_Y); 
//   smalltree->Branch("tree_seedVtx_Z", &tree_seedVtx_Z); 
//   smalltree->Branch("tree_seedVtx_dd",             &tree_seedVtx_dd); 
//   smalltree->Branch("tree_seedVtx_dphi",           &tree_seedVtx_dphi); 
//   smalltree->Branch("tree_seedVtx_distance2track", &tree_seedVtx_distance2track); 
//   smalltree->Branch("tree_seedVtx_normChi2",       &tree_seedVtx_normChi2); 
      
  // secondary  vertices 
  
  smalltree->Branch("tree_secondaryVtx_X", &tree_secondaryVtx_X);
  smalltree->Branch("tree_secondaryVtx_Y", &tree_secondaryVtx_Y);
  smalltree->Branch("tree_secondaryVtx_Z", &tree_secondaryVtx_Z);
  
  smalltree->Branch("tree_secondaryVtx_diff_X",  &tree_secondaryVtx_diff_X);
  smalltree->Branch("tree_secondaryVtx_diff_Y",  &tree_secondaryVtx_diff_Y);
  smalltree->Branch("tree_secondaryVtx_diff_Z",  &tree_secondaryVtx_diff_Z); 
  smalltree->Branch("tree_secondaryVtx_nTracks", &tree_secondaryVtx_nTracks); 
  smalltree->Branch("tree_secondaryVtx_isValid", &tree_secondaryVtx_isValid); 
  smalltree->Branch("tree_secondaryVtx_NChi2",   &tree_secondaryVtx_NChi2); 
  
  //// iterative 
  
//   smalltree->Branch("tree_secondaryVtx_iterative_X"      ,    &tree_secondaryVtx_iterative_X); 
//   smalltree->Branch("tree_secondaryVtx_iterative_Y"      ,    &tree_secondaryVtx_iterative_Y); 
//   smalltree->Branch("tree_secondaryVtx_iterative_Z"      ,    &tree_secondaryVtx_iterative_Z); 
//   smalltree->Branch("tree_secondaryVtx_iterative_nTracks",    &tree_secondaryVtx_iterative_nTracks);
//   smalltree->Branch("tree_secondaryVtx_iterative_NChi2"  ,    &tree_secondaryVtx_iterative_NChi2); 
//   smalltree->Branch("tree_secondaryVtx_iterative_isSelected", &tree_secondaryVtx_iterative_isSelected);  
  
  // Tracking particle info  
  
  smalltree->Branch("tree_simtrack_simtrack_charge",	&tree_simtrack_simtrack_charge );			 
  smalltree->Branch("tree_simtrack_simtrack_pt",	&tree_simtrack_simtrack_pt );		 
  smalltree->Branch("tree_simtrack_simtrack_eta",	&tree_simtrack_simtrack_eta  ); 	 
  smalltree->Branch("tree_simtrack_simtrack_phi",	&tree_simtrack_simtrack_phi  ); 	 
  smalltree->Branch("tree_simtrack_simtrack_longLived",	&tree_simtrack_simtrack_longLived );
  smalltree->Branch("tree_simtrack_simtrack_pdgId",	&tree_simtrack_simtrack_pdgId );	   
  smalltree->Branch("tree_simtrack_simtrack_mass",	&tree_simtrack_simtrack_mass   );			 
  smalltree->Branch("tree_simtrack_simtrack_status",	&tree_simtrack_simtrack_status );			  
  smalltree->Branch("tree_simtrack_genVertexPos_X",     &tree_simtrack_genVertexPos_X );		  
  smalltree->Branch("tree_simtrack_genVertexPos_Y",     &tree_simtrack_genVertexPos_Y );		  
  smalltree->Branch("tree_simtrack_genVertexPos_Z",     &tree_simtrack_genVertexPos_Z );		  
  smalltree->Branch("tree_simtrack_isRecoMatched",      &tree_simtrack_isRecoMatched  );
  smalltree->Branch("tree_simtrack_pca_dxy",            &tree_simtrack_pca_dxy);
  smalltree->Branch("tree_simtrack_pca_dz",             &tree_simtrack_pca_dz);
  smalltree->Branch("tree_simtrack_trkIdx",             &tree_simtrack_trkIdx); 
  smalltree->Branch("tree_simtrack_isRecoMatched_pt",   &tree_simtrack_isRecoMatched_pt);/*!*/

  smalltree->Branch("tree_AK4Slimmedjet_E"  ,       &tree_AK4Slimmedjet_E);
  smalltree->Branch("tree_AK4Slimmedjet_pt"  ,      &tree_AK4Slimmedjet_pt);
  smalltree->Branch("tree_AK4Slimmedjet_eta" ,      &tree_AK4Slimmedjet_eta);
  smalltree->Branch("tree_AK4Slimmedjet_phi" ,      &tree_AK4Slimmedjet_phi);
  smalltree->Branch("tree_AK4Slimmedjet_idxTrack" , &tree_AK4Slimmedjet_idxTrack); 

  smalltree->Branch("tree_AK4PFjet_E"  ,       &tree_AK4PFjet_E);
  smalltree->Branch("tree_AK4PFjet_pt"  ,      &tree_AK4PFjet_pt);
  smalltree->Branch("tree_AK4PFjet_eta" ,      &tree_AK4PFjet_eta);
  smalltree->Branch("tree_AK4PFjet_phi" ,      &tree_AK4PFjet_phi);
  smalltree->Branch("tree_AK4PFjet_idxTrack" , &tree_AK4PFjet_idxTrack); 
  
  smalltree->Branch("tree_CaloJet_E"  ,      &tree_CaloJet_E);
  smalltree->Branch("tree_CaloJet_pt"  ,     &tree_CaloJet_pt);
  smalltree->Branch("tree_CaloJet_eta" ,     &tree_CaloJet_eta);
  smalltree->Branch("tree_CaloJet_phi" ,     &tree_CaloJet_phi);
  smalltree->Branch("tree_CaloJet_idxTrack", &tree_CaloJet_idxTrack); 

  smalltree->Branch("tree_jet08_E"  ,      &tree_jet08_E);
  smalltree->Branch("tree_jet08_pt"  ,     &tree_jet08_pt);
  smalltree->Branch("tree_jet08_eta" ,     &tree_jet08_eta);
  smalltree->Branch("tree_jet08_phi" ,     &tree_jet08_phi);
  smalltree->Branch("tree_jet08_idxTrack", &tree_jet08_idxTrack); 

  //   smalltree->Branch("tree_CaloJet08_E"   , &tree_CaloJet08_E);
  //   smalltree->Branch("tree_CaloJet08_pt"  , &tree_CaloJet08_pt);
  //   smalltree->Branch("tree_CaloJet08_eta" , &tree_CaloJet08_eta);
  //   smalltree->Branch("tree_CaloJet08_phi" , &tree_CaloJet08_phi);
  //   smalltree->Branch("tree_CaloJet08_idxTrack", &tree_CaloJet08_idxTrack);  
  
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
  smalltree->Branch("tree_genParticle_mother_pt" ,    &tree_genParticle_mother_pt);
  smalltree->Branch("tree_genParticle_mother_eta" ,   &tree_genParticle_mother_eta);
  smalltree->Branch("tree_genParticle_mother_phi" ,   &tree_genParticle_mother_phi);

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

// from Paul:
  smalltree->Branch("tree_nLLP",&tree_nLLP);
  smalltree->Branch("tree_LLP1_x",&tree_LLP1_x);
  smalltree->Branch("tree_LLP1_y",&tree_LLP1_y);
  smalltree->Branch("tree_LLP1_z",&tree_LLP1_z);    
  smalltree->Branch("tree_LLP1_pt",&tree_LLP1_pt);
  smalltree->Branch("tree_LLP1_eta",&tree_LLP1_eta);
  smalltree->Branch("tree_LLP1_phi",&tree_LLP1_phi);

  smalltree->Branch("tree_LLP2_x",&tree_LLP2_x);
  smalltree->Branch("tree_LLP2_y",&tree_LLP2_y);
  smalltree->Branch("tree_LLP2_z",&tree_LLP2_z);
  smalltree->Branch("tree_LLP2_pt",&tree_LLP2_pt);
  smalltree->Branch("tree_LLP2_eta",&tree_LLP2_eta);
  smalltree->Branch("tree_LLP2_phi",&tree_LLP2_phi);

  smalltree->Branch("tree_NSim",&tree_NSim);
  smalltree->Branch("tree_RecoEfficiency",&tree_RecoEfficiency);
  smalltree->Branch("tree_FakeRate",&tree_FakeRate);
  smalltree->Branch("tree_RecoVertex1_Eff",&tree_RecoVertex1_Eff);
  smalltree->Branch("tree_RecoVertex2_Eff",&tree_RecoVertex2_Eff);
  smalltree->Branch("tree_VtxReco_Eff",&tree_VtxReco_Eff);

  smalltree->Branch("tree_Vtx_LLP1_nTrks",&tree_Vtx_LLP1_nTrks);
  smalltree->Branch("tree_Vtx_LLP1_x",&tree_Vtx_LLP1_x);
  smalltree->Branch("tree_Vtx_LLP1_y",&tree_Vtx_LLP1_y);
  smalltree->Branch("tree_Vtx_LLP1_z",&tree_Vtx_LLP1_z);
  smalltree->Branch("tree_Vtx_LLP1_NChi2",&tree_Vtx_LLP1_NChi2);

  smalltree->Branch("tree_Vtx_LLP2_nTrks",&tree_Vtx_LLP2_nTrks);
  smalltree->Branch("tree_Vtx_LLP2_x",&tree_Vtx_LLP2_x);
  smalltree->Branch("tree_Vtx_LLP2_y",&tree_Vtx_LLP2_y);
  smalltree->Branch("tree_Vtx_LLP2_z",&tree_Vtx_LLP2_z);
  smalltree->Branch("tree_Vtx_LLP2_NChi2",&tree_Vtx_LLP2_NChi2);

  smalltree->Branch("tree_Vtx_mva_LLP1_nTrks",&tree_Vtx_mva_LLP1_nTrks);
  smalltree->Branch("tree_Vtx_mva_LLP1_x",&tree_Vtx_mva_LLP1_x);
  smalltree->Branch("tree_Vtx_mva_LLP1_y",&tree_Vtx_mva_LLP1_y);
  smalltree->Branch("tree_Vtx_mva_LLP1_z",&tree_Vtx_mva_LLP1_z);
  smalltree->Branch("tree_Vtx_mva_LLP1_NChi2",&tree_Vtx_mva_LLP1_NChi2);

  smalltree->Branch("tree_Vtx_mva_LLP2_nTrks",&tree_Vtx_mva_LLP2_nTrks);
  smalltree->Branch("tree_Vtx_mva_LLP2_x",&tree_Vtx_mva_LLP2_x);
  smalltree->Branch("tree_Vtx_mva_LLP2_y",&tree_Vtx_mva_LLP2_y);
  smalltree->Branch("tree_Vtx_mva_LLP2_z",&tree_Vtx_mva_LLP2_z);
  smalltree->Branch("tree_Vtx_mva_LLP2_NChi2",&tree_Vtx_mva_LLP2_NChi2);

//   //First top analysis
// 
//   smalltree->Branch("tree_seedVtx_X_top1",&tree_seedVtx_X_top1);
//   smalltree->Branch("tree_seedVtx_Y_top1",&tree_seedVtx_Y_top1);
//   smalltree->Branch("tree_seedVtx_Z_top1",&tree_seedVtx_Z_top1);
//   smalltree->Branch("tree_seedVtx_dd_top1",&tree_seedVtx_dd_top1);
//   smalltree->Branch("tree_seedVtx_dphi_top1",&tree_seedVtx_dphi_top1);
//   smalltree->Branch("tree_seedVtx_distance2track_top1",&tree_seedVtx_distance2track_top1); 
//   smalltree->Branch("tree_seedVtx_normChi2_top1",&tree_seedVtx_normChi2_top1);
//   smalltree->Branch("tree_VtxReco_Eff_top1",&tree_VtxReco_Eff_top1);
// 
//   //Second top analysis
// 
//   smalltree->Branch("tree_seedVtx_X_top2",&tree_seedVtx_X_top2);
//   smalltree->Branch("tree_seedVtx_Y_top2",&tree_seedVtx_Y_top2);
//   smalltree->Branch("tree_seedVtx_Z_top2",&tree_seedVtx_Z_top2);
//   smalltree->Branch("tree_seedVtx_dd_top2",&tree_seedVtx_dd_top2);
//   smalltree->Branch("tree_seedVtx_dphi_top2",&tree_seedVtx_dphi_top2);
//   smalltree->Branch("tree_seedVtx_distance2track_top2",&tree_seedVtx_distance2track_top2); 
//   smalltree->Branch("tree_seedVtx_normChi2_top2",&tree_seedVtx_normChi2_top2);
//   smalltree->Branch("tree_VtxReco_Eff_top2",&tree_VtxReco_Eff_top2);

//   smalltree->Branch("tree_DVertex_top1_nTrks",&tree_DVertex_top1_nTrks);
//   smalltree->Branch("tree_DVertex_top2_nTrks",&tree_DVertex_top2_nTrks);
//   smalltree->Branch("tree_DVertex_nTrks",&tree_DVertex_nTrks);


  smalltree->Branch("tree_seedVtx_dr_top1",&tree_seedVtx_dr_top1);
  smalltree->Branch("tree_seedVtx_dz_top1",&tree_seedVtx_dz_top1);
  smalltree->Branch("tree_seedVtx_dr_top2",&tree_seedVtx_dr_top2);
  smalltree->Branch("tree_seedVtx_dz_top2",&tree_seedVtx_dz_top2);

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
  bool showlog = true;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  runNumber = iEvent.id().run();
  //std::cout << "runNumber = " << runNumber << std::endl;
  eventNumber = iEvent.id().event();
  //std::cout << "eventNumber = "<< eventNumber <<std::endl; 
  lumiBlock = iEvent.luminosityBlock();
  
  //// HANDLES /////
  // Triggers
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_,triggerBits);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  
  if ( nEvent == 0 ) getHLTPathNames(names);
  fillTriggerBits(HLTPathsByName_, triggerBits, names);
  
  edm::Handle<  edm::View<reco::Track>  > TracksForRes;
  iEvent.getByToken(trackSrc_, TracksForRes);
  
  edm::Handle<edm::View<reco::Track> > tracksHandle;
  iEvent.getByToken(trackToken_, tracksHandle);
  const edm::View<reco::Track>& tracks = *tracksHandle;
  
  edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
  iSetup.get<TransientRecHitRecord>().get(ttrhbuilder_,theTrackerRecHitBuilder);
  
  Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(beamSpotToken_, recoBeamSpotHandle);
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexToken_, vertices);
  
  edm::Handle<edm::View<reco::Jet> > ak4slimmedJets;
  edm::Handle<edm::View<reco::Jet> > ak4PFJets;
  edm::Handle<edm::View<reco::Jet> > CaloJets;
  edm::Handle<edm::View<reco::Jet> > ak8jets;
  //   edm::Handle<edm::View<reco::Jet> > ak8CaloJets;
  
  iEvent.getByToken(ak4slimmedJetToken_,ak4slimmedJets); 
  iEvent.getByToken(ak4PFJetToken_,ak4PFJets); 
  iEvent.getByToken(CaloJetToken_,CaloJets);  
  iEvent.getByToken(ak8jetToken_,ak8jets);  
  //   iEvent.getByToken(ak8CaloJetToken_,ak8CaloJets);  
  
  edm::Handle<pat::METCollection> PFMETs;
  iEvent.getByToken(metToken_, PFMETs);
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  if(!runOnData_) iEvent.getByToken(genParticlesToken_,genParticles);

  edm::Handle<reco::GenJetCollection> genJets;
  if(!runOnData_) iEvent.getByToken(genJetToken_,genJets);
  
  edm::Handle<reco::GenJetCollection> ak8GenJets;
  if(!runOnData_) iEvent.getByToken(ak8GenJetToken_,ak8GenJets);

  edm::Handle<GenEventInfoProduct> genEventInfo;
  if(!runOnData_) iEvent.getByToken(genEventInfoToken_,genEventInfo);

  edm::Handle<LHEEventProduct> lheEventProduct;
  if(!runOnData_) iEvent.getByToken(LHEEventProductToken_,lheEventProduct);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfcandsToken_,pfcands);

  edm::Handle<pat::ElectronCollection> electronsPAT;
  iEvent.getByToken(electronPATToken_,electronsPAT);

  edm::Handle<pat::MuonCollection> slimmedmuons;
  iEvent.getByToken(slimmedmuonToken_,slimmedmuons);

  edm::Handle<vector<reco::CompositeCandidate> > dimuons;
  if ( !runOnData_ ) 
  { 
    iEvent.getByToken(ZmumuCandidate_, dimuons);
    const vector<reco::CompositeCandidate> theDimuonToChecks = *dimuons;
    //     std::cout << "found " << theDimuonToChecks.size() << " z candidates " << endl;
    tree_NbrOfZCand = theDimuonToChecks.size();
  }
  
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
  const reco::TrackToTrackingParticleAssociator& associatorByHits = *theAssociator;
  
  //tracking Particles info and matching to reco 
  
  edm::ESHandle<ParametersDefinerForTP> parametersDefinerH;
  if(!runOnData_) {
    iSetup.get<TrackAssociatorRecord>().get(parametersDefinerName_, parametersDefinerH);
  }
  
  const ParametersDefinerForTP *parametersDefiner = parametersDefinerH.product();
  
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
  ////////	BS     /////////////
  //////////////////////////////////
  //////////////////////////////////
  
  BeamSpot const & bs = *recoBeamSpotHandle;
  tree_bs_PosX = bs.x0();
  tree_bs_PosY = bs.y0();
  tree_bs_PosZ = bs.z0();
  
  //////////////////////////////////
  //////////////////////////////////
  ////////   Tracks  /////////////
  //////////////////////////////////
  //////////////////////////////////
  ///Tracks and their association to other objects 
  
  edm::RefToBaseVector<reco::Track> trackRefs;
  size_t idxTrack = 0;
  
  std::map<size_t , int > trackToVextexMap;
  std::map<size_t , int > trackToAK4SlimmedJetMap;
  std::map<size_t , int > trackToAK4PFJetMap;
  std::map<size_t , int > trackTo08JetMap; 
  std::map<size_t , int > trackToCaloJetMap; 
  std::map<size_t , int > trackTo08CaloJetMap; 
  std::map<size_t , int > jetToTrackMap;
  std::map<size_t , int > CaloJetToTrackMap;
  std::map<size_t , int > ak8jetToTrackMap;
  //	std::map<size_t , int > ak8CaloJetToTrackMap;
  std::map<size_t , int > trackToPFJetMap;
  
  ///// TRACK ASSOCIATION TO VERTICES AND JETS 
  
  //	int nRecoTracks = tracks.size(); 
  for (edm::View<Track>::size_type i=0; i<tracks.size(); ++i) 
  {
      //---------------------------
      //minimum selection on tracks
      if ( tracks[i].pt() < 0.9 || fabs(tracks[i].eta()) > 2.5 ) continue;
      //float trackpT = tracks[i].pt();
      trackRefs.push_back(tracks.refAt(i));
      
      //---------------------------
      //loop on vertex to do track-vertex association
      //---------------------------
      int idxvtx = 0;
      int thematchidx = -1;
      float dzmin = std::numeric_limits<float>::max();
      for (auto const & vertex : *vertices) {
  	//---------------------------
  	//loop on vertex to do track-vertex association;
  	float dz = std::abs(tracks[i].dz( vertex.position()));
  	if ( dz < dzmin ) {
  	  dzmin = dz;
  	  thematchidx = idxvtx;
  	 }
  	 //---------------------------
  	 idxvtx++;
      }
      
      trackToVextexMap[idxTrack] = thematchidx;
      
      //---------------------------
      //loop on jets to do track-jet association
      //---------------------------
      
      int idxSlimmedJet=0;
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
      
      int idxAK4PFJet=0;
      bool found_match_ak4 = false;
      for (unsigned int ij=0;ij<ak4PFJets->size();ij++) {
        const Jet& jet = ak4PFJets->at(ij);
        if ( jet.pt() < 20. ) continue; 
        TLorentzVector jet4m, track4m;
        jet4m.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0);
        track4m.SetPtEtaPhiM(tracks[i].pt(), tracks[i].eta(), tracks[i].phi(), 0);
        if ( jet4m.DeltaR(track4m) < 0.4 ) {
          found_match_ak4 = true;
          break;
        }
        else idxAK4PFJet++;
      }
      if (found_match_ak4) trackToAK4PFJetMap[idxTrack] = idxAK4PFJet;
      else		   trackToAK4PFJetMap[idxTrack] = -1;
      
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
      
      int idxCaloJet=0;
      bool found_match_calo = false;
      for (unsigned int ij=0;ij<CaloJets->size();ij++) {
        const Jet& jet = CaloJets->at(ij);
        if ( jet.pt() < 20. ) continue; 
        TLorentzVector jet4m, track4m;
        jet4m.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0);
        track4m.SetPtEtaPhiM(tracks[i].pt(), tracks[i].eta(), tracks[i].phi(), 0);
        if ( jet4m.DeltaR(track4m) < 0.4 ) {
          found_match_calo = true;
          break;
        }
        else idxCaloJet++;
      }
      if (found_match_calo) trackToCaloJetMap[idxTrack] = idxCaloJet;
      else		    trackToCaloJetMap[idxTrack] = -1;
      
      //       int idxCaloJet08=0;
      //       bool found_match_calo08 = false;
      //       for(unsigned int ij=0;ij<ak8CaloJets ->size();ij++){
      //	   const Jet& jet = ak8CaloJets->at(ij);
      //	   if(jet.pt() < 20) continue; 
      //	   TLorentzVector jet4m, track4m;
      //	   jet4m.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0);
      //	   track4m.SetPtEtaPhiM(tracks[i].pt(), tracks[i].eta(), tracks[i].phi(), 0);
      //	   if( jet4m.DeltaR(track4m) < 0.8){
      //	   found_match_calo08 = true;
      //	   break;
      //	 }
      //	 else idxCaloJet08++;
      //       }
      //       if(found_match_calo08) trackTo08CaloJetMap[idxTrack] = idxCaloJet08;
      //       else trackTo08CaloJetMap[idxTrack] = -1;
      
      idxTrack++;
    } 
  
  //////////////////////////////////
  //////////////////////////////////
  //////   primary vertices  ///////
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
  
  //////////////////////////////////
  //////////////////////////////////
  ////////   MET   /////////////////
  //////////////////////////////////
  //////////////////////////////////
  
  //const Met& met = PFMET->at(0); 
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
  //////// gen Information /////////
  //////////////////////////////////
  //////////////////////////////////
  
//$$$$
  int nLLP = 0;
  int nllp = 0;
  
  tree_GenPVx = 0.;
  tree_GenPVy = 0.;
  tree_GenPVz = 0.;

  if ( !runOnData_ ) {
    //      int nGenParticles = genParticles->size(); 
    //cout << "number of gen particles "<<nGenParticles<<endl; 
    
    ///GEN PARTICLES 
  
    for (auto const & genParticle : *genParticles) {
      int pdgid = abs(genParticle.pdgId());

      // smuon
      if ( pdgid == 1000013 ) {  
        tree_GenPVx = genParticle.vx();
        tree_GenPVy = genParticle.vy();
        tree_GenPVz = genParticle.vz();
      }

      // neutralino from smuon
      if ( pdgid == 1000023 ) {  
        const Candidate * mom = genParticle.mother();
        if ( abs(mom->pdgId()) == 1000013 ) {  
  	  nLLP++;
  	  if ( nLLP == 1 ) {
  	    tree_LLP1_pt  = genParticle.pt();
  	    tree_LLP1_eta = genParticle.eta();
  	    tree_LLP1_phi = genParticle.phi();	 
  	  } 
  	  if ( nLLP == 2 ) {
  	    tree_LLP2_pt  = genParticle.pt();
  	    tree_LLP2_eta = genParticle.eta();
  	    tree_LLP2_phi = genParticle.phi();	 
  	  } 
        }
      }
        
      // quarks from neutralino
      if ( pdgid >= 1 && pdgid <= 6 ) {
        const Candidate * mom = genParticle.mother();
        if ( abs(mom->pdgId()) == 1000023 ) { 
          if ( nllp >= 2 ) {
            float dV1 = (genParticle.vx() - tree_LLP1_x)*(genParticle.vx() - tree_LLP1_x)
        	      + (genParticle.vy() - tree_LLP1_y)*(genParticle.vy() - tree_LLP1_y)
        	      + (genParticle.vz() - tree_LLP1_z)*(genParticle.vz() - tree_LLP1_z);
            float dV2 = (genParticle.vx() - tree_LLP2_x)*(genParticle.vx() - tree_LLP2_x)
        	      + (genParticle.vy() - tree_LLP2_y)*(genParticle.vy() - tree_LLP2_y)
        	      + (genParticle.vz() - tree_LLP2_z)*(genParticle.vz() - tree_LLP2_z);
            if ( dV1 > 0.01 && dV2 > 0.01 ) nllp++; // should be == 2, so just to check
          } 
          if ( nllp == 1 ) {
            float dV = (genParticle.vx() - tree_LLP1_x)*(genParticle.vx() - tree_LLP1_x)
        	     + (genParticle.vy() - tree_LLP1_y)*(genParticle.vy() - tree_LLP1_y)
        	     + (genParticle.vz() - tree_LLP1_z)*(genParticle.vz() - tree_LLP1_z);
            if ( dV > 0.01 ) {
              nllp = 2;
              tree_LLP2_x = genParticle.vx();
              tree_LLP2_y = genParticle.vy();
              tree_LLP2_z = genParticle.vz();
            } 
          } 
          if ( nllp == 0 ) {
            nllp = 1;
            tree_LLP1_x = genParticle.vx();
            tree_LLP1_y = genParticle.vy();
            tree_LLP1_z = genParticle.vz();
          } 
        } 
      } 

//$$
    if (genParticle.pt() < 0.9 || fabs(genParticle.eta()) > 4.0) continue;
//$$
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
      const Candidate * mom = genParticle.mother();
      tree_genParticle_mother_pdgId.push_back( mom ? mom->pdgId() :  -1 );// llp donne neutralino, peut faire des masse sinvariantes si besoin
      tree_genParticle_mother_pt.push_back( mom ? mom->pt() :  -1 );
      tree_genParticle_mother_eta.push_back( mom ? mom->eta() :  -10 );
      tree_genParticle_mother_phi.push_back( mom ? mom->phi() :  -10 );
    }
    
    tree_nLLP = nllp;
    
    //	  cout << " nllp " << nllp << "   LLP1 " << tree_LLP1_x << " " << tree_LLP1_y << " " << tree_LLP1_z  
    //				   << "   LLP2 " << tree_LLP2_x << " " << tree_LLP2_y << " " << tree_LLP2_z << endl; 
//$$$$

    ///GEN JETS 
    
    //       int nGenJets = genJets->size();
    //cout << "number of gen Jets "<<nGenJets<<endl; 
    
    //FIX ME : genJet vx, vy and vz not filled in tree
    
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
  }
   
  /// Electrons 
  
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
    }
  
  //  ///  reco Muons

  //  int nmuon=0;
  //  float muon1_x = 0, muon1_y = 0, muon1_z = 0; //LLP1
  //  float muon2_x = 0, muon2_y = 0, muon2_z = 0; //LLP2 avoir ces informations l� dnas le NTuple
  //  float muon1_pt = -1, muon1_eta= -10, muon1_phi=-10;
  //  float muon2_pt = -1, muon2_eta= -10, muon2_phi=-10;
  //  float muon2_E = ;float muon2_dxy = ;float muon2_dxyError =  ;float muon2_dz =  ;float muon2_dzError =  ;float muon2_charge =  ;
  //  bool muon2_PFisoVeryTight = ;bool muon2_PFisoTight = ;bool muon2_PFisoMedium = ;bool muon2_PFisoLoose = ;
  //  bool muon2_MVAisoLoose = ;bool muon2_MVAisoMedium = ;bool muon2_MVAisoTight = ;bool muon2_isGlobalMuon = ;
  //  bool muon2_isStandAloneMuon = ;bool muon2_CutBasedIdLoose = ;bool muon2_CutBasedIdMedium = ;
  //  bool muon2_CutBasedIdMediumPrompt = ;bool muon2_CutBasedIdTight =; 
   
  for (auto const & muon : *slimmedmuons)
    {
      if ( muon.pt() < 3. ) continue;
      tree_slimmedmuon_pt.push_back( muon.pt());
      tree_slimmedmuon_eta.push_back(muon.eta());
      tree_slimmedmuon_phi.push_back(muon.phi());
      tree_slimmedmuon_vx.push_back( muon.vx());
      tree_slimmedmuon_vy.push_back( muon.vy());
      tree_slimmedmuon_vz.push_back( muon.vz());
      tree_slimmedmuon_energy.push_back(   muon.energy());
      tree_slimmedmuon_dxy.push_back(	 muon.bestTrack()->dxy(bs));
      tree_slimmedmuon_dxyError.push_back( muon.bestTrack()->dxyError());
      tree_slimmedmuon_dz.push_back(	 muon.bestTrack()->dz(bs.position()));
      tree_slimmedmuon_dzError.push_back(  muon.bestTrack()->dzError());
      tree_slimmedmuon_charge.push_back(   muon.charge()); 
      
      tree_slimmedmuon_PFisoVeryTight.  push_back(   muon.passed(reco::Muon::PFIsoVeryTight));
      tree_slimmedmuon_PFisoTight.	push_back(   muon.passed(reco::Muon::PFIsoTight) );
      tree_slimmedmuon_PFisoMedium.	push_back(   muon.passed(reco::Muon::PFIsoMedium));
      tree_slimmedmuon_PFisoLoose.	push_back(   muon.passed(reco::Muon::PFIsoLoose ));
      tree_slimmedmuon_MVAisoLoose.	push_back(   muon.passed(reco::Muon::MvaLoose )  );
      tree_slimmedmuon_MVAisoMedium.	push_back(   muon.passed(reco::Muon::MvaMedium)  );
      tree_slimmedmuon_MVAisoTight.	push_back(   muon.passed(reco::Muon::MvaTight )  );
      tree_slimmedmuon_isGlobalMuon.	push_back(   muon.passed(muon.isGlobalMuon())  );
      tree_slimmedmuon_isStandAloneMuon.push_back(   muon.passed(muon.isStandAloneMuon()));
      tree_slimmedmuon_CutBasedIdLoose  	  .push_back(	muon.passed(reco::Muon::CutBasedIdLoose));
      tree_slimmedmuon_CutBasedIdMedium     .push_back(   muon.passed(reco::Muon::CutBasedIdMedium));
      tree_slimmedmuon_CutBasedIdMediumPrompt	 .push_back(   muon.passed(reco::Muon::CutBasedIdMediumPrompt));
      tree_slimmedmuon_CutBasedIdTight  	  .push_back(	muon.passed(reco::Muon::CutBasedIdTight));
    }
     //  std::cout<< "number of muons first method : " << tree_slimmedmuon_pt.size() << std::endl;


  //reco::BeamSpot vertexBeamSpot= *recoBeamSpotHandle;
  edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);
  
  /////////////////
  //prepare association to tracks by hit
  reco::RecoToSimCollection recSimColl ;
  if(!runOnData_) recSimColl= associatorByHits.associateRecoToSim(trackRefs, tpCollection);
  //fake rate determination : when a reco track has no matched simtrack
  
  int nTracks = 0; 
  //	int nUnmatchTrack_fromPU = 0;
  //	int nUnmatchTrack_fromPV = 0;
  //	int nUnmatchTrack= 0;
  //	int nPUTrack= 0;
  
  //---------------
  //loops on tracks 
  //---------------

  for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack) {
    const auto& itTrack = trackRefs[iTrack];
    //------------------------
    //general track properties
    //------------------------
    //std::cout << "  " << itTrack->pt()  << std::endl;

    tree_track_charge.push_back(itTrack->charge());
    tree_track_pt.push_back(itTrack->pt());
    tree_track_eta.push_back(itTrack->eta());
    tree_track_phi.push_back(itTrack->phi());
    tree_track_nhits.push_back(itTrack->numberOfValidHits());
    tree_track_NChi2.push_back(itTrack->normalizedChi2());
    tree_track_outerPt.push_back(itTrack->outerPt()); 
    tree_track_vx.push_back(itTrack->vx()); 
    tree_track_vy.push_back(itTrack->vy()); 
    tree_track_vz.push_back(itTrack->vz()); 
    tree_track_firsthit_X.push_back(itTrack->innerPosition().X());  
    tree_track_firsthit_Y.push_back(itTrack->innerPosition().Y()); 
    tree_track_firsthit_Z.push_back(itTrack->innerPosition().Z()); 
    tree_track_firsthit_phi.push_back(itTrack->innerPosition().phi()); 
    
    if( itTrack->quality(reco::TrackBase::highPurity) ){tree_track_isHighQuality.push_back(true);} 
    else {tree_track_isHighQuality.push_back(false);}
    if( itTrack->quality(reco::TrackBase::loose) )     {tree_track_isLoose.push_back(true);} 
    else {tree_track_isLoose.push_back(false);}  
    if( itTrack->quality(reco::TrackBase::tight))      {tree_track_isTight.push_back(true);}
    else {tree_track_isTight.push_back(false);}        
    
    tree_track_dxy.push_back(		 itTrack->dxy(bs.position()));
    tree_track_dxyError.push_back(	 itTrack->dxyError());
    tree_track_dz.push_back(		 itTrack->dz(bs.position()));
    tree_track_dzError.push_back(	 itTrack->dzError());
    tree_track_numberOfLostHits.push_back( itTrack->numberOfLostHits());
    tree_track_numberOfValidHits.push_back(itTrack->numberOfValidHits());
    
    tree_track_originalAlgo.push_back(itTrack->originalAlgo());
    tree_track_algo.push_back(itTrack->algo());
    tree_track_stopReason.push_back(itTrack->stopReason());
    
    //--------------------------------
    //general hit properties of tracks
    //--------------------------------
    
    const reco::HitPattern& hp = itTrack->hitPattern();
    //     tree_track_nPixel   .push_back(hp.numberOfValidPixelHits());
    //     tree_track_nStrip   .push_back(hp.numberOfValidStripHits());
    
    tree_track_numberOfValidPixelHits.push_back(hp.numberOfValidPixelHits());
    tree_track_numberOfValidStripHits.push_back(hp.numberOfValidStripHits());
    tree_track_numberOfValidStripTIBHits.push_back(hp.numberOfValidStripTIBHits());
    tree_track_numberOfValidStripTIDHits.push_back(hp.numberOfValidStripTIDHits());
    tree_track_numberOfValidStripTOBHits.push_back(hp.numberOfValidStripTOBHits());
    tree_track_numberOfValidStripTECHits.push_back(hp.numberOfValidStripTECHits());
    tree_track_numberOfValidPixelBarrelHits.push_back(hp.numberOfValidPixelBarrelHits());
    tree_track_numberOfValidPixelEndcapHits.push_back(hp.numberOfValidPixelEndcapHits());
    tree_track_trackerLayersWithMeasurement.push_back(hp.trackerLayersWithMeasurement());
    tree_track_pixelLayersWithMeasurement.push_back(hp.pixelLayersWithMeasurement());
    tree_track_stripTECLayersWithMeasurement.push_back(hp.stripTECLayersWithMeasurement() );
    tree_track_stripTIBLayersWithMeasurement.push_back(hp.stripTIBLayersWithMeasurement());
    tree_track_stripTIDLayersWithMeasurement.push_back(hp.stripTIDLayersWithMeasurement());
    tree_track_stripTOBLayersWithMeasurement.push_back(hp.stripTOBLayersWithMeasurement());
    
    int hitPixelLayer = 0;
    if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) )  hitPixelLayer += 1;
    if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2) )  hitPixelLayer += 10;
    if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3) )  hitPixelLayer += 100;
    if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 4) )  hitPixelLayer += 1000;
    
    if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1) )  hitPixelLayer += 2;
    if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 2) )  hitPixelLayer += 20;
    if(hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 3) )  hitPixelLayer += 200;
    
    tree_track_hasValidHitInPixelLayer.push_back(hitPixelLayer);
    
    //----------------------------
    //matching to simulated tracks
    //----------------------------
    //matching par hits

    if ( !runOnData_ ) {
      float NReco=trackRefs.size();
      float NRecotoSim=0;
      float ratioReco=0;
  
      int nSimHits = 0;
      bool isSimMatched = false;
      std::vector<int> tpIdx;
      std::vector<float> sharedFraction;
      std::vector<float> tpChi2;
  
      //initialized values for trackingParticle
  
      float genVertexPos_X = -1000;
      float genVertexPos_Y = -1000;
      float genVertexPos_Z = -1000;
  
      int   simtrack_charge		 =-1000;
      float simtrack_pt 		 =-1;
      float simtrack_eta		 =-1000;
      float simtrack_phi		 =-1000;
      bool  simtrack_longLived  	 = false;
      //int   simtrack_matchedHit	 = 0;
      int   simtrack_pdgId		 = -1000;
      int   simtrack_numberOfTrackerHits   = 0;
      int   simtrack_numberOfTrackerLayers = 0;
      float simtrack_mass = 0.;
      int   simtrack_status= -1000; 
      bool  simtrack_isFromDispTop	   = false;
      int   simtrack_isFromLLP  	  = 0;
      float llp1_Diff_dV = -1.;
      float llp2_Diff_dV = -1.;
  
      auto foundTPs = recSimColl.find(itTrack);
      if ( foundTPs != recSimColl.end() ) {
        //if (!foundTPs->val.empty()) {
        isSimMatched = true;
        TrackingParticleRef tpr = foundTPs->val[0].first;
        NRecotoSim++;
        nSimHits = tpr->numberOfTrackerHits();
        
        simtrack_charge 	       = tpr->charge();     
        simtrack_pt		       = tpr->pt();	    
        simtrack_eta		       = tpr->eta();	       
        simtrack_phi		       = tpr->phi();	       
        simtrack_longLived	       = tpr->longLived();    
        //simtrack_matchedHit		 = tpr->matchedHit();	     
        simtrack_pdgId  	       = tpr->pdgId();    
        simtrack_numberOfTrackerHits   = tpr->numberOfTrackerHits(); 
        simtrack_numberOfTrackerLayers = tpr->numberOfTrackerLayers();
        simtrack_mass		       = tpr->mass();
        simtrack_status 	       = tpr->status();
        
        //determine x,y,z position of the genVertex which produced the associated simtrack
        genVertexPos_X         = tpr->vx();
        genVertexPos_Y         = tpr->vy();
        genVertexPos_Z         = tpr->vz();
        
//$$$$
        // std::cout<< "nllp : "<< nllp<<std::endl;
        if ( nllp >= 1 ) {
          float dV1 = (genVertexPos_X - tree_LLP1_x)*(genVertexPos_X - tree_LLP1_x)
                    + (genVertexPos_Y - tree_LLP1_y)*(genVertexPos_Y - tree_LLP1_y)
                    + (genVertexPos_Z - tree_LLP1_z)*(genVertexPos_Z - tree_LLP1_z);
          // std::cout<< "genVertexPosX : "<<genVertexPos_X<<" & LLP1_X : "<< LLP1_x<< " & dV1 : "<< dV1<<std::endl;
          llp1_Diff_dV = TMath::Sqrt( dV1 ); 
          if ( dV1 < 0.01 ) simtrack_isFromLLP = 1;
        }
        if ( nllp >= 2 && simtrack_isFromLLP != 1 ) {
          float dV2 = (genVertexPos_X - tree_LLP2_x)*(genVertexPos_X - tree_LLP2_x)
                    + (genVertexPos_Y - tree_LLP2_y)*(genVertexPos_Y - tree_LLP2_y)
                    + (genVertexPos_Z - tree_LLP2_z)*(genVertexPos_Z - tree_LLP2_z);
          llp2_Diff_dV = TMath::Sqrt( dV2 ); 
          if ( dV2 < 0.01 ) simtrack_isFromLLP = nllp;
        }
//$$$$
    	
        //ex: s�lectionner les traces venant du premeir top. Regarder � quelles poitns le  vertex se rapprotche du vertex g�n�r�
        //aussi pour le deuxieme top
        //Matcher les informations avec les x y z calcul�s pr�c�demment

        for (auto genIterator=genParticles->begin(); genIterator!=genParticles->end(); genIterator++)
        {
          if (abs(genIterator->pdgId())==6 && genIterator->status()==22) //if the gen particle is a displaced top 
          // if (abs(genIterator->pdgId())<=6 && abs(genIterator->pdgId())>=1)
          {  //Doublon de la m�thode de  Daniel ligne ~1500//
            double top_x=genIterator->vx(); 
            double top_y=genIterator->vy(); 
            double top_z=genIterator->vz(); 
            //cout << "genVertexPos_X - top_x " << genVertexPos_X - top_x << endl;
            //cout << "genVertexPos_Y - top_y " << genVertexPos_Y - top_y << endl;
            //cout << "genVertexPos_Z - top_z " << genVertexPos_Z - top_z << endl;
            if( fabs(genVertexPos_X - top_x)< 0.01 && fabs(genVertexPos_Y - top_y) < 0.01 && fabs(genVertexPos_Z - top_z) < 0.01) simtrack_isFromDispTop = true;
          }
        }
      }

      ratioReco = NRecotoSim/NReco;  
      tree_FakeRate.push_back(ratioReco);
      
      tree_track_nSimHits		       .push_back(nSimHits); 
      tree_track_isSimMatched		       .push_back(isSimMatched);

      tree_track_simtrack_charge	       .push_back(simtrack_charge);	 
      tree_track_simtrack_pt		       .push_back(simtrack_pt); 		 
      tree_track_simtrack_eta		       .push_back(simtrack_eta);	      
      tree_track_simtrack_phi		       .push_back(simtrack_phi);		      
      tree_track_simtrack_longLived	       .push_back(simtrack_longLived);  		  
      tree_track_simtrack_pdgId 	       .push_back(simtrack_pdgId);	
      tree_track_simtrack_numberOfTrackerHits  .push_back(simtrack_numberOfTrackerHits); 
      tree_track_simtrack_numberOfTrackerLayers.push_back(simtrack_numberOfTrackerLayers); 
      tree_track_simtrack_mass  	       .push_back(simtrack_mass        );   
      tree_track_simtrack_status	       .push_back(simtrack_status      );   
      
      tree_track_genVertexPos_X 	       .push_back(genVertexPos_X);		 
      tree_track_genVertexPos_Y 	       .push_back(genVertexPos_Y);		 
      tree_track_genVertexPos_Z 	       .push_back(genVertexPos_Z); 
      tree_track_simtrack_llp1_dV	       .push_back(llp1_Diff_dV);
      tree_track_simtrack_llp2_dV	       .push_back(llp2_Diff_dV);
      tree_track_simtrack_isFromLLP	       .push_back(simtrack_isFromLLP); 
      tree_track_simtrack_isFromDispTop        .push_back(simtrack_isFromDispTop      );    
    }

    nTracks++; 
     
    //if( itTrack->pt() < 1 ) std::cout << " there is a pt problem " << itTrack->pt() << std::endl;
    
    /*if( !isSimMatched &&  trackToVextexMap[iTrack] != 0 ) nUnmatchTrack_fromPU++;
    if( !isSimMatched &&  trackToVextexMap[iTrack] == 0 ) nUnmatchTrack_fromPV++;
    if( !isSimMatched					) nUnmatchTrack++;
    if( trackToVextexMap[iTrack] != 0			) nPUTrack++;*/
    
    tree_track_recoVertex_idx.push_back(trackToVextexMap[iTrack]);
    tree_track_recoAK4SlimmedJet_idx.push_back(trackToAK4SlimmedJetMap[iTrack]);
    tree_track_recoAK4PFJet_idx.push_back(trackToAK4PFJetMap[iTrack]);
    tree_track_reco08Jet_idx.push_back(trackTo08JetMap[iTrack]); 
    tree_track_recoCaloJet_idx.push_back(trackToCaloJetMap[iTrack]); 
    //     tree_track_reco08CaloJet_idx.push_back(trackTo08CaloJetMap[iTrack]); 
    
  } //end loop on tracks

  /*std::cout << "---------------------"  << std::endl;
  std::cout << "nUnmatchTrack_fromPU " << nUnmatchTrack_fromPU << std::endl;
  std::cout << "nUnmatchTrack_fromPV " << nUnmatchTrack_fromPV << std::endl;
  std::cout << "nUnmatchTrack        " << nUnmatchTrack<< std::endl;
  std::cout << "nTrack               " << tree_track_charge.size() << std::endl;
  std::cout << "nPUTrack             " << nPUTrack<< std::endl;
  std::cout << "nNonPUTrack          " << tree_track_charge.size()- nPUTrack<< std::endl;*/
  

  //------------------------------
  // loops on tracking particles 
  //------------------------------
  
  if ( !runOnData_ ) {
    reco::SimToRecoCollection simRecColl = associatorByHits.associateSimToReco(trackRefs, tpCollection);
    //Reco efficiecny estimation by looking at the number of simtracks matched with reco tracks
      float NSim=0;/*!*/
      float NSimtoReco=0;/*!*/
      float ratio=0;/*!*/

    for (const TrackingParticleRef& tp: tpCollection) 
    {
      tree_simtrack_simtrack_charge		   .push_back( tp->charge());	//infos sur toutes les traces simul�es, 
      //ex efficacit� reco, voir trace sim � trace reco, compar�e aux pt de otoutes les traces sim /*!*///finito
      tree_simtrack_simtrack_pt 		   .push_back( tp->pt());	    //avec un rapport (attention de bien prendre un pt simul�)/*!*/
      tree_simtrack_simtrack_eta		   .push_back( tp->eta());	    
      tree_simtrack_simtrack_phi		   .push_back( tp->phi());	    
      tree_simtrack_simtrack_longLived      .push_back( tp->longLived());   
      //tree_simtrack_simtrack_matchedHit	     .push_back( tp->matchedHit());    
      tree_simtrack_simtrack_pdgId		   .push_back( tp->pdgId());	
      tree_simtrack_simtrack_numberOfTrackerHits   .push_back( tp->numberOfTrackerHits()); 
      tree_simtrack_simtrack_numberOfTrackerLayers .push_back( tp->numberOfTrackerLayers());
      tree_simtrack_simtrack_mass		   .push_back( tp->mass());
      tree_simtrack_simtrack_status		   .push_back( tp->status());
      tree_simtrack_genVertexPos_X		   .push_back( tp->vx());
      tree_simtrack_genVertexPos_Y		   .push_back( tp->vy());
      tree_simtrack_genVertexPos_Z		   .push_back( tp->vz());
      
      bool isRecoMatched = false;
      NSim = NSim + 1;/*!*/
      std::vector<int> tkIdx;
      auto foundTracks = simRecColl.find(tp);
      if ( foundTracks != simRecColl.end() ) 
      {
    	isRecoMatched = true;
    	NSimtoReco=NSimtoReco+1;/*!*/
    	for (const auto trackQuality: foundTracks->val) 
              {
        	     tkIdx.push_back(trackQuality.first.key());
              }
      }
       
      tree_simtrack_isRecoMatched.push_back(isRecoMatched);//tree_simtrack_pt when isRecoMatched =true
       //Calcualte the impact parameters w.r.t. PCA
       tree_simtrack_trkIdx.push_back(tkIdx);

       TrackingParticle::Vector momentum = parametersDefiner->momentum(iEvent,iSetup,tp);
       TrackingParticle::Point vertex = parametersDefiner->vertex(iEvent,iSetup,tp);
       auto dxySim = TrackingParticleIP::dxy(vertex, momentum, bs.position());
       auto dzSim  = TrackingParticleIP::dz(vertex, momentum, bs.position());

       tree_simtrack_pca_dxy	  .push_back(dxySim);
       tree_simtrack_pca_dz	  .push_back(dzSim);

    } //end loop on tracking particles 

    ratio= NSimtoReco/NSim;/*!*/
    //  std::cout<< "ratio of NSimtoReco/Nsim = "<< ratio << std::endl;
    tree_NSim.push_back(NSim);
    tree_RecoEfficiency.push_back(ratio);//change en fonction du pt et �ta
  }

   
  //////////////////////////////////
  //////////////////////////////////
  ///////////   jets   /////////////
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
      TLorentzVector thejet, track;
      thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
      for (unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) 
        {
          track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
          if(thejet.DeltaR(track) < 0.4 ) tree_AK4Slimmedjet_idxTrack.push_back(iTrack);    
        }
    }
  
  for (int ij=0;ij<int(ak4PFJets->size());ij++)
    {
      const Jet& jet = ak4PFJets->at(ij); 
      if ( jet.pt() < 20. ) continue;
      tree_AK4PFjet_E.push_back(jet.energy());
      tree_AK4PFjet_pt.push_back(jet.pt());
      tree_AK4PFjet_eta.push_back(jet.eta());
      tree_AK4PFjet_phi.push_back(jet.phi());
      TLorentzVector thejet, track;
      thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
      for (unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) 
      {
  	track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
  	if(thejet.DeltaR(track) < 0.4 ) tree_AK4PFjet_idxTrack.push_back(iTrack);     
      }
    }
  
  for (int ij=0;ij< int(CaloJets->size());ij++)//plus simpliste, utilisete just eles cellules calorimetriques
    {
      const Jet& jet = CaloJets->at(ij);
      if ( jet.pt() < 20. ) continue; 
      tree_CaloJet_E.push_back(jet.energy());
      tree_CaloJet_pt.push_back(jet.pt());
      tree_CaloJet_eta.push_back(jet.eta());
      tree_CaloJet_phi.push_back(jet.phi());
      TLorentzVector thejet, track;
      thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
      for (unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) {
  	track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
  	if(thejet.DeltaR(track) < 0.4 ) tree_CaloJet_idxTrack.push_back(iTrack);    
       }
    }

    for (int ij=0;ij< int(ak8jets->size() );ij++)
    {
      const Jet& jet = ak8jets->at(ij);
      if ( jet.pt() < 20. ) continue;
      tree_jet08_E.push_back(jet.energy());
      tree_jet08_pt.push_back(jet.pt());
      tree_jet08_eta.push_back(jet.eta());
      tree_jet08_phi.push_back(jet.phi());
      TLorentzVector thejet, track;
      thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
      for (unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) 
      {
  	track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
  	if(thejet.DeltaR(track) < 0.8 ) tree_jet08_idxTrack.push_back(iTrack);    
      }
    }
  
  //	for(int ij=0;ij< int(ak8CaloJets->size()) ;ij++)
  //	 {
  //	  const Jet& jet = ak8CaloJets->at(ij);
  //	if ( jet.pt() < 20. ) continue;
  //	  tree_CaloJet08_E.push_back(jet.energy());
  //	  tree_CaloJet08_pt.push_back(jet.pt());
  //	  tree_CaloJet08_eta.push_back(jet.eta());
  //	  tree_CaloJet08_phi.push_back(jet.phi());
  //	  TLorentzVector thejet, track;
  //	  thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
  //	  for(unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) 
  //	   {
  //	    track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
  //	    if(thejet.DeltaR(track) < 0.8 ) tree_CaloJet08_idxTrack.push_back(iTrack);    
  //	   }
  //	 }
  
  //Nothing is done in Daniel's code to treat secondary. This is treated
  //in the next part by J�r�my, this is why there are "doublons"
  /////SECONDARY VERTEX RECONSTRUCTION USING DISPLACED TRACKS (J�r�my's method) here LLPs are only tops
   TransientVertex displacedVertex_top1_general;
   TransientVertex displacedVertex_top2_general;

   vector<reco::TransientTrack> displacedTracks_top1;
   vector<reco::TransientTrack> displacedTracks_top2;

  if ( !runOnData_ )
  {
   if(showlog) cout << "  --------- " <<endl;  
   if(showlog) cout << "SECONDARY VERTEX REFITTING	 "<<endl; 
   
   int n_displacedtop=0; 
   
   float secondaryVtx_diff_X = -1.; 
   float secondaryVtx_diff_Y = -1.;    
   float secondaryVtx_diff_Z = -1.; 
   int secondaryVtx_nTracks = -1; 
   int secondaryVtx_isValid = -3; 
   float secondaryVtx_NChi2=-1000; 

   float top1_x=0; 
   float top1_y=0; 
   float top1_z=0; 

   float top2_x=0; 
   float top2_y=0; 
   float top2_z=0;

   // genParticle.pdgId()
   for (auto genIterator=genParticles->begin(); genIterator!=genParticles->end(); genIterator++) //loop on gen particles
   {
     if (abs(genIterator->pdgId())==6 && genIterator->status()==22) //if the gen particle is a displaced top ,
     // if ( abs(genIterator->pdgId()) >= 1 && abs(genIterator->pdgId()) <= 6 )//daniel's condition, taking into account all quarks
             { 
               n_displacedtop++; 
               int num_track=-1; 
               if (n_displacedtop==1) 
  	 {
  	   top1_x=genIterator->vx(); //doublon ligne 1450 : loop on genParticles 
  	   top1_y=genIterator->vy(); // but here it is more restrective 
  	   top1_z=genIterator->vz(); // as all final quarks are not cosnidered and we are considering only top and not all LLPs
  	   if(showlog) cout << "FIRST TOP GEN VERTEX X "<<genIterator->vx()<<endl; 
  	   if(showlog) cout << "FIRST TOP GEN VERTEX Y "<<genIterator->vy()<<endl; 
  	   if(showlog) cout << "FIRST TOP GEN VERTEX Z "<<genIterator->vz()<<endl; 
  	   tree_genTop_X.push_back(top1_x);
  	   tree_genTop_Y.push_back(top1_y);
  	   tree_genTop_Z.push_back(top1_z);
  	   tree_genTop_charge.push_back(genIterator->charge());//only 0s??
               }
            
               if (n_displacedtop==2) 
  	 {
  	   top2_x=genIterator->vx();  //doublon ligne 1450 : loop on genParticles
  	   top2_y=genIterator->vy();  // but here it is more restrective 
  	   top2_z=genIterator->vz();  // as all final quarks are cosnidered and we are considering only top and not all LLPs
  	   if(showlog) cout << "------------------------" <<endl; 
  	   if(showlog) cout << "SECOND TOP GEN VERTEX X "<<genIterator->vx()<<endl; 
  	   if(showlog) cout << "SECOND TOP GEN VERTEX Y "<<genIterator->vy()<<endl; 
  	   if(showlog) cout << "SECOND TOP GEN VERTEX Z "<<genIterator->vz()<<endl; 
  	   tree_genTop_X.push_back(top2_x);
  	   tree_genTop_Y.push_back(top2_y);
  	   tree_genTop_Z.push_back(top2_z);
  	   tree_genTop_charge.push_back(genIterator->charge());//only 0s??
  	 }
  	 if (n_displacedtop>2)
  	   {
  	     std::cout<<"More than two tops for this event : "<< n_displacedtop <<std::endl;
  	   }
  		cout << "------------------------" <<endl; 
  	 
  	 //using gen information for vertex refitting
              for (auto trackIterator=tracks.begin(); trackIterator!=tracks.end(); trackIterator++) //loop on reco tracks
              {
        	     if(trackIterator->pt() < 1 || fabs(trackIterator->eta()) > 2.4) continue;
        	     num_track++; 
        	     if (tree_track_isSimMatched[num_track]==0) continue; //keeping only the reco tracks matched to sim tracks 
        	     if (abs(tree_track_genVertexPos_X[num_track]-genIterator->vx())<0.01 && abs(tree_track_genVertexPos_Y[num_track]-genIterator->vy())<0.01 && abs(tree_track_genVertexPos_Z[num_track]-genIterator->vz())<0.01)
  	 //same restriction as Daniel's method
  	 //if the position of the sim track associated to the reco track is the same as the gen top vertex
        	     {///genVertexPos X which produces de associated simtrack
        	       TransientTrack  transientTrack = theTransientTrackBuilder->build(*trackIterator); 
        	       if (n_displacedtop==1) displacedTracks_top1.push_back(transientTrack);
        	       if (n_displacedtop==2) displacedTracks_top2.push_back(transientTrack);
        	     }
              }
            
             }
      }//end loop on gen Particles of j�r�my's method
    
   if(showlog) cout << "------------------------" <<endl;
   if(showlog) cout << "Number of reco tracks for first top "<<displacedTracks_top1.size()<<endl; 
   if(showlog) cout << "Number of reco tracks for second top "<<displacedTracks_top2.size()<<endl; 
   int nReconstructibleVertex1=0;
   int nReconstructibleVertex2=0;
   int nRecoVertex1=0;
   int nRecoVertex2=0;
   if ( displacedTracks_top1.size() > 1 )
      {
       nReconstructibleVertex1++;
             KalmanVertexFitter theFitter_top1; 
             //AdaptiveVertexFitter  theFitter_top1;	
       
             TransientVertex displacedVertex_top1 = theFitter_top1.vertex(displacedTracks_top1);  // if you want the beam constraint
       displacedVertex_top1_general = displacedVertex_top1;//useless for the moment
             //TransientVertex myVertex = theFitter.vertex(mytracks);  // if you don't want the beam constraint
             // now you have a new vertex, can e.g. be compared with the original
             if (displacedVertex_top1.isValid()) // NotValid if the max number of steps has been exceedeor the fitted position is out of tracker bounds.
             {    
  	 nRecoVertex1++;
  	 auto lineariser = DefaultLinearizationPointFinder();
               auto linearization_point = lineariser.getLinearizationPoint(displacedTracks_top1);
               if(showlog) std::cout << "seed top1 is: " << linearization_point.x() << " - " << linearization_point.y() << " - " << linearization_point.z() << '\n';
        	     if(showlog) cout << "------------------------" <<endl;
        	     if(showlog) std::cout << " FIRST TOP DISPLACED VERTEX X POS " << displacedVertex_top1.position().x() << std::endl;//reco
        	     if(showlog) std::cout << " FIRST TOP DISPLACED VERTEX Y POS " << displacedVertex_top1.position().y() << std::endl;
        	     if(showlog) std::cout << " FIRST TOP DISPLACED VERTEX Z POS " << displacedVertex_top1.position().z() << std::endl; 
        	     if(showlog) std::cout << " FIRST TOP DISPLACED VERTEX CHI2/ndof  " << displacedVertex_top1.normalisedChiSquared()<< std::endl; 
        	     if(showlog) std::cout << " NUMBER OF ORIGINAL TRACKS      " << displacedVertex_top1.originalTracks().size() << std::endl; 
           
               secondaryVtx_diff_X = abs(displacedVertex_top1.position().x()-top1_x);//Kalman fitted vertex-original genVertex
               secondaryVtx_diff_Y = abs(displacedVertex_top1.position().y()-top1_y);
               secondaryVtx_diff_Z = abs(displacedVertex_top1.position().z()-top1_z);
               secondaryVtx_nTracks = displacedVertex_top1.originalTracks().size();
               secondaryVtx_isValid = 1;
               secondaryVtx_NChi2 = displacedVertex_top1.normalisedChiSquared();
            
               tree_secondaryVtx_X.push_back(displacedVertex_top1.position().x()); 
               tree_secondaryVtx_Y.push_back(displacedVertex_top1.position().y()); 
               tree_secondaryVtx_Z.push_back(displacedVertex_top1.position().z()); 
            
               tree_secondaryVtx_diff_X.push_back(secondaryVtx_diff_X);//Kalman fitted vertex
               tree_secondaryVtx_diff_Y.push_back(secondaryVtx_diff_Y);
               tree_secondaryVtx_diff_Z.push_back(secondaryVtx_diff_Z);
               tree_secondaryVtx_nTracks.push_back(secondaryVtx_nTracks); 
               tree_secondaryVtx_isValid.push_back(secondaryVtx_isValid);
               tree_secondaryVtx_NChi2.push_back(secondaryVtx_NChi2);//Nombre  entr�es �trange.
             } 
             else { 
               if(showlog) cout << "------------------------" <<endl;
  	 if(showlog) cout << "Refitted vertex for top 1 is not okay "<<endl;
  	 secondaryVtx_isValid=0; 
  	 tree_secondaryVtx_isValid.push_back(secondaryVtx_isValid);
       }
   }
   else{
           //	   cout << "------------------------" <<endl;
           //	   std::cout << "not enough tracks left for top 1" <<std::endl;
           secondaryVtx_isValid=-1;
           tree_secondaryVtx_isValid.push_back(secondaryVtx_isValid);
   }
    
   if ( displacedTracks_top2.size() > 1 )
   { 
       nReconstructibleVertex2++;
             auto lineariser = DefaultLinearizationPointFinder();
             auto linearization_point = lineariser.getLinearizationPoint(displacedTracks_top2);
             if(showlog) std::cout << "seed top 2 is: " << linearization_point.x() << " - " << linearization_point.y() << " - " << linearization_point.z() << '\n';
             KalmanVertexFitter theFitter_top2; 
             //AdaptiveVertexFitter  theFitter_top2;
        
       TransientVertex displacedVertex_top2 = theFitter_top2.vertex(displacedTracks_top2);  // if you want the beam constraint
       displacedVertex_top2_general = displacedVertex_top2;
       //TransientVertex myVertex = theFitter.vertex(mytracks);  // if you don't want the beam constraint
       // now you have a new vertex, can e.g. be compared with the original
       if (displacedVertex_top2.isValid()) //NotValid if the max number of steps has been exceedeor the fitted position is out of tracker bounds.
       { 
  	 nRecoVertex2++;   
  	 if(showlog)	    cout << "------------------------" <<endl;
  	 if(showlog) std::cout << " SECOND TOP DISPLACED VERTEX X POS " << displacedVertex_top2.position().x() << std::endl;
  	 if(showlog) std::cout << " SECOND TOP DISPLACED VERTEX Y POS " << displacedVertex_top2.position().y() << std::endl;
  	 if(showlog) std::cout << " SECOND TOP DISPLACED VERTEX Z POS " << displacedVertex_top2.position().z() << std::endl; 
  	 if(showlog) std::cout << " SECOND TOP DISPLACED VERTEX CHI2/ndof  " << displacedVertex_top2.normalisedChiSquared()<< std::endl; 
  	 if(showlog) std::cout << " NUMBER OF ORIGINAL TRACKS	     " << displacedVertex_top2.originalTracks().size() << std::endl;
  	 
  	 tree_secondaryVtx_X.push_back(displacedVertex_top2.position().x()); 
  	 tree_secondaryVtx_Y.push_back(displacedVertex_top2.position().y()); 
  	 tree_secondaryVtx_Z.push_back(displacedVertex_top2.position().z()); 
  	 
  	 secondaryVtx_diff_X = abs(displacedVertex_top2.position().x()-top2_x);//Kalman fitted vertex-original genVertex
  	 secondaryVtx_diff_Y = abs(displacedVertex_top2.position().y()-top2_y);
  	 secondaryVtx_diff_Z = abs(displacedVertex_top2.position().z()-top2_z);
  	 secondaryVtx_nTracks = displacedVertex_top2.originalTracks().size();
  	 secondaryVtx_isValid = 1;
  	 secondaryVtx_NChi2 = displacedVertex_top2.normalisedChiSquared();
  	 
  	 tree_secondaryVtx_diff_X.push_back(secondaryVtx_diff_X);
  	 tree_secondaryVtx_diff_Y.push_back(secondaryVtx_diff_Y);
  	 tree_secondaryVtx_diff_Z.push_back(secondaryVtx_diff_Z);
  	 tree_secondaryVtx_nTracks.push_back(secondaryVtx_nTracks); 
  	 tree_secondaryVtx_isValid.push_back(secondaryVtx_isValid);
  	 tree_secondaryVtx_NChi2.push_back(secondaryVtx_NChi2);//Nombre  entr�es �trange.
       } 
       else{
  	 secondaryVtx_isValid=0;
  	 tree_secondaryVtx_isValid.push_back(secondaryVtx_isValid);
  	 if(showlog) cout << "------------------------" <<endl;
  	 if(showlog) cout << "Refitted vertex for top 2 is not valid "<<endl;
       }
   }
      else
   {
  	secondaryVtx_isValid=-1;
  	tree_secondaryVtx_isValid.push_back(secondaryVtx_isValid);
  	//     cout << "------------------------" <<endl;
  	//     std::cout << "not enough tracks left for top 2" <<std::endl;
   }

   //Vertex Reconstruction efficiency : Reconstructed/Reconstructible vertices
   //Reconstructible: if more than one track is available for the recosntruction
   //Reconstructed : isValid()
   if (nReconstructibleVertex1!=0)
     {
       // std::cout<<"nreconstrible 1"<< std::endl;
       float RecoVertex1_Eff = nRecoVertex1/nReconstructibleVertex1;tree_RecoVertex1_Eff.push_back(RecoVertex1_Eff);
     }
   else  {tree_RecoVertex1_Eff.push_back(0);
     // std::cout<<"push back 0 pour top1"<< std::endl;
   }

   if(nReconstructibleVertex2!=0)
   {  //std::cout<<"nrecostrible 2"<< std::endl;
     float RecoVertex2_Eff = nRecoVertex2/nReconstructibleVertex2;
     tree_RecoVertex2_Eff.push_back(RecoVertex2_Eff);}
   else
   {  tree_RecoVertex2_Eff.push_back(0);
     //std::cout<<"push back 0 pour top2"<< std::endl;
   }
  }
  

  //////////////////
  //------------------------------
  //selection of displaced tracks
  //------------------------------
  //////////////////

  if(showlog) cout << "//////////////////////////"<< endl; 
  if(showlog) cout << "start to select displaced tracks " << endl;
  vector<reco::TransientTrack> displacedTracks_llp1, displacedTracks_llp2, displacedTracks_mva_llp1, displacedTracks_mva_llp2;
    
  ///// MVA for track selection coming from displaced tops 

  //ajout� par Paul /*!*/
  float drSig, isinjet; /*!*/
  int jet; /*!*/
  float ntrk10, ntrk20, ntrk30; /*!*/
  float firsthit_X, firsthit_Y, firsthit_Z, dxy, dxyError, pt, eta, NChi2, nhits, algo; /*!*/
  double bdtval = -100.;
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
  //$$
  reader->BookMVA( "BDTG", weightFile_ );
  //$$
   
  int counter_track = -1; 
  for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack)  // Loop on all the tracks
  {
    counter_track++; 
    const auto& itTrack = trackRefs[iTrack];
    firsthit_X = tree_track_firsthit_X[counter_track]; 
    firsthit_Y = tree_track_firsthit_Y[counter_track]; 
    firsthit_Z = tree_track_firsthit_Z[counter_track]; 
    dxy        = tree_track_dxy[counter_track]; 
    dxyError   = tree_track_dxyError[counter_track]; 
    pt         = tree_track_pt[counter_track]; 
    eta        = tree_track_eta[counter_track]; 
    NChi2      = tree_track_NChi2[counter_track]; 
    nhits      = tree_track_nhits[counter_track]; 
    algo       = tree_track_algo[counter_track]; 

    //Ajout� par Paul /*!*/
    drSig = -1.;
    if ( dxyError > 0 ) drSig = abs(dxy) / dxyError; /*!*/
    ntrk10 = 0;
    ntrk20 = 0;
    ntrk30 = 0;
    bdtval = -10.;

    if ( pt > 1. && NChi2 < 5. && drSig > 5. ) { // On regarde les track_selec[i] qui sont True donc de potentielles tracks secondaires

      jet = tree_track_recoAK4SlimmedJet_idx[counter_track]; /*!*/
      isinjet = 0.;
      if ( jet >= 0 ) isinjet = 1.; /*!*/

      //Computation of the distances needed for the BDT
      int counter_othertrack = -1; ///*!*/
      for (size_t iTrack2=0; iTrack2<trackRefs.size(); ++iTrack2)    // Loop on all the other Tracks/*!*/
      {
      	counter_othertrack++;
      if ( counter_othertrack == counter_track ) continue;

    	float pt2  = tree_track_pt[counter_othertrack];
    	float dr2 = abs(tree_track_dxy[counter_othertrack]);
    	float drSig2 = -1.;
    	if ( tree_track_dxyError[counter_othertrack] > 0 ) drSig2 = dr2 / tree_track_dxyError[counter_othertrack];
    	NChi2 = tree_track_NChi2[counter_othertrack];

      if ( !(pt2 > 1. && NChi2 < 5. && drSig2 > 5.) ) continue; // On regarde les autres track_selec[i] qui sont True donc de potnetielles tracks secondaires

    	float x2 = tree_track_firsthit_X[counter_othertrack];
    	float y2 = tree_track_firsthit_Y[counter_othertrack];
    	float z2 = tree_track_firsthit_Z[counter_othertrack];
    	float dist = TMath::Sqrt( (firsthit_X-x2)*(firsthit_X-x2) + (firsthit_Y-y2)*(firsthit_Y-y2) + (firsthit_Z-z2)*(firsthit_Z-z2) );//pour chaque reconstruite, on regarde les autres tracks, 
    	if ( dist < 10. )	      {ntrk10++;} // les sctocker les 3 , on teste sur une seule couche quand on regarde vers l'avant
    	if ( dist < 20. )	      {ntrk20++;} 
    	if ( dist < 30. )	      {ntrk30++;} 
      }  // end Loop on other Tracks

      if ( tree_track_simtrack_isFromLLP[counter_track] == 1 ) {
        displacedTracks_llp1.push_back(theTransientTrackBuilder->build(*itTrack));
      }
      if ( tree_track_simtrack_isFromLLP[counter_track] == 2 ) {
        displacedTracks_llp2.push_back(theTransientTrackBuilder->build(*itTrack));
      }

      bdtval = reader->EvaluateMVA( "BDTG" ); //default value = -10 (no -10 observed and -999 comes from EvaluateMVA)
      //cout << "BDT VAL " << bdtval <<endl; 
    
//$$
      // if (bdtval < 0.0674  ) { //optimal cut for 10 cm, trained on 19k events, <s>=22 and <b>=240 //J�r�my
      // if (bdtval < 0.07381 ) { //optimal cut for 50 cm, trained on 10k events, <s>=15 and <b>=234 //J�r�my
      if ( bdtval > -0.0401 ) { //optimal cut for 50 cm ,trained on 10k events //Paul, expected to change
//$$

        if ( tree_track_simtrack_isFromLLP[counter_track] == 1 ) {
          displacedTracks_mva_llp1.push_back(theTransientTrackBuilder->build(*itTrack));
        }
        if ( tree_track_simtrack_isFromLLP[counter_track] == 2 ) {
          displacedTracks_mva_llp2.push_back(theTransientTrackBuilder->build(*itTrack));
        }
        //cout << "kept track and it is from displaced top "<<tree_track_simtrack_isFromDispTop[counter_track]<<endl; 
      }
    }

    tree_track_ntrk10.push_back(ntrk10);
    tree_track_ntrk20.push_back(ntrk20);
    tree_track_ntrk30.push_back(ntrk30);
    tree_track_MVAval.push_back(bdtval); 

  } //End loop on all the tracks

  if ( showlog ) {
    cout << " displaced tracks LLP1 " << displacedTracks_llp1.size() << " " << displacedTracks_mva_llp1.size() << endl;
    cout << " displaced tracks LLP2 " << displacedTracks_llp2.size() << " " << displacedTracks_mva_llp2.size() << endl;
  }

  //---------------------------------------------------------------------------------------//


  //////////////////////////////////
  //--------------------------------
  // Vertex fitting separately for the two LLPs
  //--------------------------------
  //////////////////////////////////

  int Vtx_ntk;
  float Vtx_x, Vtx_y, Vtx_z, Vtx_chi;
  
  //--------------------------- FIRST LLP -------------------------------------//

  KalmanVertexFitter theFitter_vertex_llp1(kvfPSet); //One or less vertex: either Valid or NotValid /*!*/
  
  Vtx_ntk = displacedTracks_llp1.size();
  Vtx_x = -10.;
  Vtx_y = -10.;
  Vtx_z = -10.;
  Vtx_chi = -1.;

  if ( Vtx_ntk > 1 )
  {
    TransientVertex displacedVertex_llp1 = theFitter_vertex_llp1.vertex(displacedTracks_llp1); // fitted vertex

    if ( displacedVertex_llp1.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
    {
      Vtx_x = displacedVertex_llp1.position().x();
      Vtx_y = displacedVertex_llp1.position().y();
      Vtx_z = displacedVertex_llp1.position().z();
      Vtx_chi = displacedVertex_llp1.normalisedChiSquared();
    }
  }
  tree_Vtx_LLP1_nTrks = Vtx_ntk;	    
  tree_Vtx_LLP1_x = Vtx_x;
  tree_Vtx_LLP1_y = Vtx_y;
  tree_Vtx_LLP1_z = Vtx_z;
  tree_Vtx_LLP1_NChi2 = Vtx_chi;
    
 //--------------------------- SECOND LLP -------------------------------------//

  KalmanVertexFitter theFitter_vertex_llp2(kvfPSet); //One or less vertex: either Valid or NotValid /*!*/
  
  Vtx_ntk = displacedTracks_llp2.size();
  Vtx_x = -10.;
  Vtx_y = -10.;
  Vtx_z = -10.;
  Vtx_chi = -1.;

  if ( Vtx_ntk > 1 )
  {
    TransientVertex displacedVertex_llp2 = theFitter_vertex_llp2.vertex(displacedTracks_llp2); // fitted vertex

    if ( displacedVertex_llp2.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
    {
      Vtx_x = displacedVertex_llp2.position().x();
      Vtx_y = displacedVertex_llp2.position().y();
      Vtx_z = displacedVertex_llp2.position().z();
      Vtx_chi = displacedVertex_llp2.normalisedChiSquared();
    }
  }
  tree_Vtx_LLP2_nTrks = Vtx_ntk;	    
  tree_Vtx_LLP2_x = Vtx_x;
  tree_Vtx_LLP2_y = Vtx_y;
  tree_Vtx_LLP2_z = Vtx_z;
  tree_Vtx_LLP2_NChi2 = Vtx_chi;
    
  //------------------------------- FIRST LLP WITH MVA ----------------------------------//
  
  KalmanVertexFitter theFitter_vertex_mva_llp1(kvfPSet); //One or less vertex: either Valid or NotValid /*!*/
  
//   std::vector< std::pair< std::vector<reco::TransientTrack>, TransientVertex> > proto_vertex_mva_llp1; // container of the vertex and associated tracks
//   Proto PVtx_llp1(proto_vertex_mva_llp1); /*!*/

  Vtx_ntk = displacedTracks_mva_llp1.size();
  Vtx_x = -10.;
  Vtx_y = -10.;
  Vtx_z = -10.;
  Vtx_chi = -1.;

  if ( Vtx_ntk > 1 )
  {
    TransientVertex displacedVertex_mva_llp1 = theFitter_vertex_mva_llp1.vertex(displacedTracks_mva_llp1); // fitted vertex

    // std::cout<< "displacedVertex_mva_llp1 is built" << std::endl;

    if ( displacedVertex_mva_llp1.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
    {
      // std::cout<< "displacedVertex_mva_llp1 is valid" << std::endl;
//       for (unsigned int itracks = 0; itracks < displacedTracks_mva_llp1.size()-1; itracks++)
//       {	 
//     	for(unsigned int itracks2 = itracks+1; itracks2 < displacedTracks_mva_llp1.size(); itracks2++)
//         {
//     	  double r1 = sqrt(pow(displacedTracks_mva_llp1[itracks].innermostMeasurementState().globalPosition().x(),2)+pow(displacedTracks_mva_llp1[itracks].innermostMeasurementState().globalPosition().y(),2)); 
//     	  double r2 = sqrt(pow(displacedTracks_mva_llp1[itracks2].innermostMeasurementState().globalPosition().x(),2)+pow(displacedTracks_mva_llp1[itracks2].innermostMeasurementState().globalPosition().y(),2)); 
//     	  double phi1 = TMath::ATan2(displacedTracks_mva_llp1[itracks].innermostMeasurementState().globalPosition().y(),displacedTracks_mva_llp1[itracks].innermostMeasurementState().globalPosition().x());
//     	  double phi2 = TMath::ATan2(displacedTracks_mva_llp1[itracks2].innermostMeasurementState().globalPosition().y(),displacedTracks_mva_llp1[itracks2].innermostMeasurementState().globalPosition().x());
//     	  double dr = abs(r1-r2); 
//     	  double dz = abs(displacedTracks_mva_llp1[itracks].innermostMeasurementState().globalPosition().z()-displacedTracks_mva_llp1[itracks2].innermostMeasurementState().globalPosition().z());
//     	  double dd = sqrt(pow(dr,2)+pow(dz,2)); 
//     	  double dphi = 0.; 
//     	  if (abs(phi1 - phi2) < (TMath::Pi()) ) dphi = abs(phi1 - phi2); 
//     	  if (abs(phi1 - phi2) > (TMath::Pi()) ) dphi = 2*(TMath::Pi()) - abs(phi1 - phi2); 
//     	  double dist_to_tracks = pow( 
//     	    pow(displacedTracks_mva_llp1[itracks].innermostMeasurementState().globalPosition().x()- displacedTracks_mva_llp1[itracks2].innermostMeasurementState().globalPosition().x(), 2) +
//     	    pow(displacedTracks_mva_llp1[itracks].innermostMeasurementState().globalPosition().y()- displacedTracks_mva_llp1[itracks2].innermostMeasurementState().globalPosition().y(), 2) +
//     	    pow(displacedTracks_mva_llp1[itracks].innermostMeasurementState().globalPosition().z()- displacedTracks_mva_llp1[itracks2].innermostMeasurementState().globalPosition().z(), 2) 
//     	    , 0.5);
//     	  tree_seedVtx_dr_llp1.push_back(dr); 
//     	  tree_seedVtx_dz_llp1.push_back(dz); 
//     	  tree_seedVtx_dd_llp1.push_back(dd); 
//     	  tree_seedVtx_dphi_llp1.push_back(dphi); 
//     	  tree_seedVtx_distance2track_llp1.push_back(dist_to_tracks); 
//     	}
//       }		 

      Vtx_ntk = displacedTracks_mva_llp1.size();
      Vtx_x = displacedVertex_mva_llp1.position().x();
      Vtx_y = displacedVertex_mva_llp1.position().y();
      Vtx_z = displacedVertex_mva_llp1.position().z();
      Vtx_chi = displacedVertex_mva_llp1.normalisedChiSquared();
    }
  }
  tree_Vtx_mva_LLP1_nTrks = Vtx_ntk;	    
  tree_Vtx_mva_LLP1_x = Vtx_x;
  tree_Vtx_mva_LLP1_y = Vtx_y;
  tree_Vtx_mva_LLP1_z = Vtx_z;
  tree_Vtx_mva_LLP1_NChi2 = Vtx_chi;

  //-------------------------- SECOND LLP WITH MVA -------------------------------------//

  KalmanVertexFitter theFitter_vertex_mva_llp2(kvfPSet); //One or less vertex: either Valid or NotValid /*!*/
  
  Vtx_ntk = displacedTracks_mva_llp2.size();
  Vtx_x = -10.;
  Vtx_y = -10.;
  Vtx_z = -10.;
  Vtx_chi = -1.;

  if ( Vtx_ntk > 1 )
  {
    TransientVertex displacedVertex_mva_llp2 = theFitter_vertex_mva_llp2.vertex(displacedTracks_mva_llp2); // fitted vertex

    if ( displacedVertex_mva_llp2.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
    {
      Vtx_x = displacedVertex_mva_llp2.position().x();
      Vtx_y = displacedVertex_mva_llp2.position().y();
      Vtx_z = displacedVertex_mva_llp2.position().z();
      Vtx_chi = displacedVertex_mva_llp2.normalisedChiSquared();
    }
  }
  tree_Vtx_mva_LLP2_nTrks = Vtx_ntk;	    
  tree_Vtx_mva_LLP2_x = Vtx_x;
  tree_Vtx_mva_LLP2_y = Vtx_y;
  tree_Vtx_mva_LLP2_z = Vtx_z;
  tree_Vtx_mva_LLP2_NChi2 = Vtx_chi;
    

//Former version of the code (J�r�my)

//------------------------------------------------------------
   //reconstruct seed vertex from pair of tracks //Not used for now (Paul)
   //select those who have a chi2/ndf < 100
   //------------------------------------------------------------
    ///////////////////////////////////////////////////////////////////////Not USED ATM////////////////////////////////////////////////////
  //  KalmanVertexFitter theFitter_vertex(kvfPSet);//One or less vertex /*!*/ //Commenter �a
  //  std::vector< std::pair< std::vector<reco::TransientTrack>, TransientVertex> > proto_vertex; // container of the vertex and associated tracks
  //  Proto PVtx(proto_vertex);/*!*/
  //  for(unsigned int itracks = 0; itracks < displacedTTracks.size(); itracks++)
  //     {
  //     for(unsigned int itracks2 = itracks+1; itracks2 < displacedTTracks.size(); itracks2++)
  //       {
  //       double r1 = sqrt(pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().x(),2)+pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().y(),2)); 
  //       double r2 = sqrt(pow(displacedTTracks[itracks2].innermostMeasurementState().globalPosition().x(),2)+pow(displacedTTracks[itracks2].innermostMeasurementState().globalPosition().y(),2)); 
  //       double phi1 = TMath::ATan2(displacedTTracks[itracks].innermostMeasurementState().globalPosition().y(),displacedTTracks[itracks].innermostMeasurementState().globalPosition().x());
  //       double phi2 = TMath::ATan2(displacedTTracks[itracks2].innermostMeasurementState().globalPosition().y(),displacedTTracks[itracks2].innermostMeasurementState().globalPosition().x());
  //       double dr = abs(r1-r2); 
  //       double dz = abs(displacedTTracks[itracks].innermostMeasurementState().globalPosition().z()-displacedTTracks[itracks2].innermostMeasurementState().globalPosition().z());
  //       double dd = sqrt(pow(dr,2)+pow(dz,2)); 
  //       double dphi = 0.; 
       
  //       if (abs(phi1 - phi2) < (TMath::Pi()) ) dphi = abs(phi1 - phi2); 
  //       if (abs(phi1 - phi2) > (TMath::Pi()) ) dphi = 2*(TMath::Pi()) - abs(phi1 - phi2); 
  //       double dist_to_tracks = pow( 
  //       pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().x()- displacedTTracks[itracks2].innermostMeasurementState().globalPosition().x(), 2) +
  //       pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().y()- displacedTTracks[itracks2].innermostMeasurementState().globalPosition().y(), 2) +
  //       pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().z()- displacedTTracks[itracks2].innermostMeasurementState().globalPosition().z(), 2) 
  //       , 0.5
  //       );
  //       // !=dd
  //       //if (dd>40 || dphi>1.4) continue; 
  //       //if(dist_to_tracks > 30) continue;
  //       std::vector<reco::TransientTrack> temp_coll_tracks;//In former's version, was in the second loop and we were studying vertex generated by pairs of tracks
  //       temp_coll_tracks.push_back(displacedTTracks[itracks]);
  //       temp_coll_tracks.push_back(displacedTTracks[itracks2]); 
  //       TransientVertex displacedVertex_vertex = theFitter_vertex.vertex(temp_coll_tracks); 
  //       if(displacedVertex_vertex.isValid() ){//reconstructible, //peut prendre un temps fou
  //       //NotValid if the max number of steps has been exceedeor the fitted position is out of tracker bounds.
  //        std::pair< std::vector<reco::TransientTrack>, TransientVertex> pair_iterators_vertex(temp_coll_tracks,displacedVertex_vertex) ;
  //        tree_seedVtx_X.push_back(displacedVertex_vertex.position().x());
  //        tree_seedVtx_Y.push_back(displacedVertex_vertex.position().y());
  //        tree_seedVtx_Z.push_back(displacedVertex_vertex.position().z());
  //        tree_seedVtx_dd.push_back(dd); 
  //        tree_seedVtx_dphi.push_back(dphi); 
  //        tree_seedVtx_distance2track.push_back(dist_to_tracks); 
  //        tree_seedVtx_normChi2.push_back(displacedVertex_vertex.normalisedChiSquared()); 
  //        proto_vertex.push_back(pair_iterators_vertex);
  //        //  PVtx.GetProto().push_back(pair_iterators_vertex);/*!*//*This doesn't  work*/
  //        PVtx.PushBack(pair_iterators_vertex);/*!*//*it does*/ 
              
  //      }
  //    }
  //    tree_DVertex_nTrks.push_back(PVtx.SizeTTracks(itracks)); 
  //  }
    
  // if(showlog) cout << "-----------------" <<endl; 
  // //  if(showlog) cout << proto_vertex.size() << " proto-vertex selected " << endl;
  // if(showlog) cout << PVtx.Size() << " proto-vertex selected " << endl;/*!*/

//-----------------------------------next step-not started-------------------------------------//

   //displacedTTracks & displacedVertex_vertex
  /*!*/ // Diviser l'event en deux : tmva ou axe sprivil�gi�s dans l'event (voir root: vectuer o uaxe privil�gi�) pt eta phi?
  /*!*/ // regarder le vecteur des muons (corr�l�s) � celui des neutralinos
  // => kalman filter pour chaque hemisphere 
  //terminer le rapport l�-dessus : efficacit� et reconstruction des vertxex des neutralions
  // Daniel is looking at this atm


   //for testing
   /*
   int n_proto_match_top = 0;
   for (unsigned int i=0; i<proto_vertex.size(); i++)
     {
     bool proto_match_top = false; 
     if(   fabs(proto_vertex[i].second.position().x() - tree_gen_secondartVtxTop_X[0])  < 0.1 &&
     fabs(proto_vertex[i].second.position().y() - tree_gen_secondartVtxTop_Y[0])  < 0.1 &&
     fabs(proto_vertex[i].second.position().z() - tree_gen_secondartVtxTop_Z[0])  < 0.5  
     ) proto_match_top = true;
     if(   fabs(proto_vertex[i].second.position().x() - tree_gen_secondartVtxTop_X[1])  < 0.1 &&
     fabs(proto_vertex[i].second.position().y() - tree_gen_secondartVtxTop_Y[1])  < 0.1 &&
     fabs(proto_vertex[i].second.position().z() - tree_gen_secondartVtxTop_Z[1])  < 0.5  
     ) proto_match_top = true; 
     
     if (proto_match_top==true) n_proto_match_top++;
     
     tree_secondaryVtx_iterative_step1_NChi2.push_back(proto_vertex[i].second.normalisedChiSquared()); 
     tree_secondaryVtx_iterative_step1_nTracks.push_back(proto_vertex[i].second.originalTracks().size()); 
     tree_secondaryVtx_iterative_step1_match_top.push_back(proto_match_top); 
     tree_secondaryVtx_iterative_step1_x.push_back(proto_vertex[i].second.position().x());
     tree_secondaryVtx_iterative_step1_y.push_back(proto_vertex[i].second.position().y());
     tree_secondaryVtx_iterative_step1_z.push_back(proto_vertex[i].second.position().z());
     
     }
     cout << "-----------------"<<endl; 
     cout << n_proto_match_top << " proto vertex good "<<endl; 
     
     tree_secondaryVtx_iterative_step1_nGoodProtoVtx.push_back(n_proto_match_top); 
   */


   //------------------------------------------------------------
   //loop on all proto-vertex
   //  loop on tracks
   //      reconstruct new vertex
   //      filter new vertex based on chi2/ndf
   //      keep or reject tracks
   //------------------------------------------------------------
   //kalman fitter sur toutes les traces juste �a 
   

  //  for(unsigned int  idx_proto_vertex= 0; idx_proto_vertex < PVtx.Size(); idx_proto_vertex++){  /*!*/
  //   //  loop on track, reco vertex, test and reject or select track
  //    for(unsigned int itracks = 0; itracks < displacedTTracks.size(); itracks++){
  //      std::vector<reco::TransientTrack> temp_coll_tracks = PVtx.TTrack(idx_proto_vertex);/*!*/
  //     //  first hit distance to vertex
  //      double dist_to_vertex = pow( 
  //        pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().x()- PVtx.x(idx_proto_vertex), 2) +
  //        pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().y()- PVtx.y(idx_proto_vertex), 2) +
  //        pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().z()- PVtx.z(idx_proto_vertex), 2) 
  //        , 0.5);
  //     //  cout << "distance check to vertex " << dist_to_vertex << endl;
  //     // Selection criteria tochange 
  //      if( //dist_to_vertex < 30 && 
  //         fabs(displacedTTracks[itracks].innermostMeasurementState().globalPosition().x() - temp_coll_tracks[0].innermostMeasurementState().globalPosition().x()  ) > 0.001 &&  
  //         fabs(displacedTTracks[itracks].innermostMeasurementState().globalPosition().x() - temp_coll_tracks[1].innermostMeasurementState().globalPosition().x()  ) > 0.001 ){
  //           temp_coll_tracks.push_back(displacedTTracks[itracks]);
  //           TransientVertex displacedVertex_vertex = theFitter_vertex.vertex(temp_coll_tracks); 
  //           // /chi2 < 100 in original 
  //           if(displacedVertex_vertex.isValid() && displacedVertex_vertex.normalisedChiSquared() < 50){
  //             PVtx.TTrack(idx_proto_vertex)=temp_coll_tracks;/*!*/
  //             PVtx.TVertex(idx_proto_vertex)= displacedVertex_vertex;/*!*/
  //           } else temp_coll_tracks.pop_back();
  //       }
  //    }
  //  }
  // //reconstruit
  // int nVtxIsReco = PVtx.Size();
  // float VtxReco_Eff=0;
  // if(showlog)  cout << "-----------------"<<endl; 
  // //  if(showlog) cout << " number of vertex after adding tracks " << proto_vertex.size() << endl;
  // if(showlog) cout << " number of vertex after adding tracks " << PVtx.Size() << endl;/*!*/

  // if (nVtxIsValid!=0)
  // {VtxReco_Eff = nVtxIsReco/nVtxIsValid;
  // tree_VtxReco_Eff.push_back(VtxReco_Eff);}
  // else
  // {tree_VtxReco_Eff.push_back(0);}
  


  // //top1
  //  for(unsigned int  idx_proto_vertex= 0; idx_proto_vertex < PVtx_top1.Size(); idx_proto_vertex++){  /*!*/
  //   //  loop on track, reco vertex, test and reject or select track
  //    for(unsigned int itracks = 0; itracks < displacedTracks_top1.size(); itracks++){
  //      std::vector<reco::TransientTrack> temp_coll_tracks = PVtx_top1.TTrack(idx_proto_vertex);/*!*/
  //     //  first hit distance to vertex
  //      double dist_to_vertex = pow( 
  //        pow(displacedTracks_top1[itracks].innermostMeasurementState().globalPosition().x()- PVtx_top1.x(idx_proto_vertex), 2) +
  //        pow(displacedTracks_top1[itracks].innermostMeasurementState().globalPosition().y()- PVtx_top1.y(idx_proto_vertex), 2) +
  //        pow(displacedTracks_top1[itracks].innermostMeasurementState().globalPosition().z()- PVtx_top1.z(idx_proto_vertex), 2) 
  //        , 0.5);
  //     //  cout << "distance check to vertex " << dist_to_vertex << endl;
  //     // Selection criteria tochange 
  //      if( //dist_to_vertex < 30 && 
  //         fabs(displacedTracks_top1[itracks].innermostMeasurementState().globalPosition().x() - temp_coll_tracks[0].innermostMeasurementState().globalPosition().x()  ) > 0.001 &&  
  //         fabs(displacedTracks_top1[itracks].innermostMeasurementState().globalPosition().x() - temp_coll_tracks[1].innermostMeasurementState().globalPosition().x()  ) > 0.001 ){
  //           temp_coll_tracks.push_back(displacedTracks_top1[itracks]);
  //           TransientVertex displacedVertex_vertex = theFitter_vertex_top1.vertex(temp_coll_tracks); 
  //           // /chi2 < 100 in original 
  //           if(displacedVertex_vertex.isValid() && displacedVertex_vertex.normalisedChiSquared() < 50){
  //             PVtx_top1.TTrack(idx_proto_vertex) = temp_coll_tracks;/*!*/
  //             PVtx_top1.TVertex(idx_proto_vertex) = displacedVertex_vertex;/*!*/
  //           } else temp_coll_tracks.pop_back();
  //       }
  //    }
  //  }
  // //reconstruit
  // int nVtxIsReco_top1 = PVtx_top1.Size();
  // float VtxReco_Eff_top1=0;
  // if(showlog)  cout << "-----------------"<<endl; 
  // //  if(showlog) cout << " number of vertex after adding tracks " << proto_vertex.size() << endl;
  // if(showlog) cout << " number of vertex after adding tracks " << PVtx_top1.Size() << endl;/*!*/

  // if (nVtx_top1_IsValid!=0)
  // {VtxReco_Eff_top1 = nVtxIsReco_top1/nVtx_top1_IsValid;
  // tree_VtxReco_Eff_top1.push_back(VtxReco_Eff_top1);}
  // else
  // {tree_VtxReco_Eff_top1.push_back(0);}

  // //top2 
  //  for(unsigned int  idx_proto_vertex= 0; idx_proto_vertex < PVtx_top2.Size(); idx_proto_vertex++){  /*!*/
  //   //  loop on track, reco vertex, test and reject or select track
  //    for(unsigned int itracks = 0; itracks < displacedTracks_top2.size(); itracks++){
  //      std::vector<reco::TransientTrack> temp_coll_tracks = PVtx.TTrack(idx_proto_vertex);/*!*/
  //     //  first hit distance to vertex
  //      double dist_to_vertex = pow( 
  //        pow(displacedTracks_top2[itracks].innermostMeasurementState().globalPosition().x()- PVtx_top2.x(idx_proto_vertex), 2) +
  //        pow(displacedTracks_top2[itracks].innermostMeasurementState().globalPosition().y()- PVtx_top2.y(idx_proto_vertex), 2) +
  //        pow(displacedTracks_top2[itracks].innermostMeasurementState().globalPosition().z()- PVtx_top2.z(idx_proto_vertex), 2) 
  //        , 0.5);
  //     //  cout << "distance check to vertex " << dist_to_vertex << endl;
  //     // Selection criteria tochange 
  //      if( //dist_to_vertex < 30 && 
  //         fabs(displacedTracks_top2[itracks].innermostMeasurementState().globalPosition().x() - temp_coll_tracks[0].innermostMeasurementState().globalPosition().x()  ) > 0.001 &&  
  //         fabs(displacedTracks_top2[itracks].innermostMeasurementState().globalPosition().x() - temp_coll_tracks[1].innermostMeasurementState().globalPosition().x()  ) > 0.001 ){
  //           temp_coll_tracks.push_back(displacedTracks_top2[itracks]);
  //           TransientVertex displacedVertex_vertex = theFitter_vertex.vertex(temp_coll_tracks); 
  //           // /chi2 < 100 in original 
  //           if(displacedVertex_vertex.isValid() && displacedVertex_vertex.normalisedChiSquared() < 50){
  //             PVtx_top2.TTrack(idx_proto_vertex)=temp_coll_tracks;/*!*/
  //             PVtx_top2.TVertex(idx_proto_vertex)= displacedVertex_vertex;/*!*/
  //           } else temp_coll_tracks.pop_back();
  //       }
  //    }
  //  }
  // //reconstruit
  // int nVtxIsReco_top2 = PVtx_top2.Size();
  // float VtxReco_Eff_top2=0;
  // if(showlog)  cout << "-----------------"<<endl; 
  // //  if(showlog) cout << " number of vertex after adding tracks " << proto_vertex.size() << endl;
  // if(showlog) cout << " number of vertex after adding tracks " << PVtx_top2.Size() << endl;/*!*/

  // if (nVtx_top2_IsValid!=0)
  // {VtxReco_Eff_top2 = nVtxIsReco_top2/nVtx_top2_IsValid;
  // tree_VtxReco_Eff_top2.push_back(VtxReco_Eff_top2);}
  // else
  // {tree_VtxReco_Eff_top2.push_back(0);}

  //for testing
   /*
  int n_proto2_match_top = 0;
  for (unsigned int i=0; i<proto_vertex.size(); i++)
     {
       bool proto2_match_top = false; 
       if(   fabs(proto_vertex[i].second.position().x() - tree_gen_secondartVtxTop_X[0])  < 0.1 &&
       fabs(proto_vertex[i].second.position().y() - tree_gen_secondartVtxTop_Y[0])  < 0.1 &&
       fabs(proto_vertex[i].second.position().z() - tree_gen_secondartVtxTop_Z[0])  < 0.5  ) proto2_match_top = true;
       if(   fabs(proto_vertex[i].second.position().x() - tree_gen_secondartVtxTop_X[1])  < 0.1 &&
       fabs(proto_vertex[i].second.position().y() - tree_gen_secondartVtxTop_Y[1])  < 0.1 &&
       fabs(proto_vertex[i].second.position().z() - tree_gen_secondartVtxTop_Z[1])  < 0.5  ) proto2_match_top = true; 
     
       if (proto2_match_top==true) n_proto2_match_top++; 
     
       tree_secondaryVtx_iterative_step2_NChi2.push_back(proto_vertex[i].second.normalisedChiSquared()); 
       tree_secondaryVtx_iterative_step2_nTracks.push_back(proto_vertex[i].second.originalTracks().size()); 
       tree_secondaryVtx_iterative_step2_match_top.push_back(proto2_match_top); 
       tree_secondaryVtx_iterative_step2_x.push_back(proto_vertex[i].second.position().x());
       tree_secondaryVtx_iterative_step2_y.push_back(proto_vertex[i].second.position().y());
       tree_secondaryVtx_iterative_step2_z.push_back(proto_vertex[i].second.position().z());
       
     }
     
     cout << n_proto2_match_top << " good vertex after adding tracks "<<endl; 
   */

   






  //Nothing has been done by Paul after that except the clearing part and the replacement of objects 
  //using the new class Proto

   //------------------------------------------------------------
   //loop on all proto-vertex
      //remove duplicate
   //------------------------------------------------------------
   
  //  std::vector< std::pair< std::vector<reco::TransientTrack>, TransientVertex> > reco_vertex;  
  //   Proto RecoVtx(reco_vertex);/*!*/
    /***************
   for(unsigned int  idx_proto_vertex= 0; idx_proto_vertex < PVtx.Size(); idx_proto_vertex++){
     
     bool found_vertex_in_list = false;
     int idx_in_list = idx_proto_vertex;
     double bestchi2 = PVtx.Nchi2(idx_proto_vertex);
     //----------------------------------------
     //check if there is another close by vertex
     //if there is, take the best chi2
     for(unsigned int  idx_proto_vertex2= 0; idx_proto_vertex2 < PVtx.Size(); idx_proto_vertex2++){
       if(idx_proto_vertex == idx_proto_vertex2) continue;
       if( 
         fabs(PVtx.x(idx_proto_vertex) - PVtx.x(idx_proto_vertex2)) < 1 && 
         fabs(PVtx.y(idx_proto_vertex) - PVtx.y(idx_proto_vertex2)) < 1 && 
         fabs(PVtx.z(idx_proto_vertex) - PVtx.z(idx_proto_vertex2)) < 5 
        ) {
          if(PVtx.Nchi2(idx_proto_vertex2) < bestchi2 ){
            idx_in_list = idx_proto_vertex2;
            found_vertex_in_list = true;
            bestchi2 = PVtx.Nchi2(idx_proto_vertex2);
          }
        }
     }
     //---------------------------------------
     //check if already in final vertex vector 
     bool already_in_vect = false;
     for(unsigned int  idx_reco_vertex= 0; idx_reco_vertex < RecoVtx.Size(); idx_reco_vertex++){
       
       if( 
         fabs(PVtx.x(idx_in_list) - RecoVtx.x(idx_reco_vertex)) < 0.0001 && 
         fabs(PVtx.y(idx_in_list) - RecoVtx.y(idx_reco_vertex)) < 0.0001 && 
         fabs(PVtx.z(idx_in_list) - RecoVtx.z(idx_reco_vertex)) < 0.0001 
         ){
           already_in_vect = true;
       }
     }
     //cout << "found_vertex_in_list " << found_vertex_in_list << endl;
     //cout << "already_in_vect      " << already_in_vect      << endl;
     
     if( found_vertex_in_list  && !already_in_vect)  reco_vertex.push_back(proto_vertex[idx_in_list]);
     if(!found_vertex_in_list  && !already_in_vect)  reco_vertex.push_back(proto_vertex[idx_proto_vertex]);
    //  if( found_vertex_in_list  && !already_in_vect)  RecoVtx.PushBack(PVtx.Pair(idx_in_list));
    //  if(!found_vertex_in_list  && !already_in_vect)  RecoVtx.PushBack(PVtx.Pair(idx_proto_vertex));
   }*///////////////////
   
  //  if(showlog) cout << reco_vertex.size() << " number of vertex reco after duplicate removal " << endl;
  //  if(showlog) cout << RecoVtx.Size() << " number of vertex reco after duplicate removal " << endl;/*!*/
   ////for testing
   
   /*int n_reco_match_top = 0;
     for (unsigned int i=0; i<reco_vertex.size(); i++)
     {
     bool reco_match_top = false; 
     if(   fabs(reco_vertex[i].second.position().x() - tree_gen_secondartVtxTop_X[0])  < 0.1 &&
     fabs(reco_vertex[i].second.position().y() - tree_gen_secondartVtxTop_Y[0])  < 0.1 &&
     fabs(reco_vertex[i].second.position().z() - tree_gen_secondartVtxTop_Z[0])  < 0.5  
     ) reco_match_top = true;
     if(   fabs(reco_vertex[i].second.position().x() - tree_gen_secondartVtxTop_X[1])  < 0.1 &&
     fabs(reco_vertex[i].second.position().y() - tree_gen_secondartVtxTop_Y[1])  < 0.1 &&
     fabs(reco_vertex[i].second.position().z() - tree_gen_secondartVtxTop_Z[1])  < 0.5  
     ) reco_match_top = true; 
     
     if (reco_match_top==true) n_reco_match_top++; 
     
     tree_secondaryVtx_iterative_step3_NChi2.push_back(reco_vertex[i].second.normalisedChiSquared()); 
     tree_secondaryVtx_iterative_step3_nTracks.push_back(reco_vertex[i].second.originalTracks().size()); 
     tree_secondaryVtx_iterative_step3_match_top.push_back(reco_match_top); 
     tree_secondaryVtx_iterative_step3_x.push_back(reco_vertex[i].second.position().x());
     tree_secondaryVtx_iterative_step3_y.push_back(reco_vertex[i].second.position().y());
     tree_secondaryVtx_iterative_step3_z.push_back(reco_vertex[i].second.position().z());
     }
     
     //  cout << n_reco_match_top << " good vertex after duplicate removal "<<endl; */
   
   //------------------------------------------------------------
   //select final vertex, filtering #tracks and chi2/ndf
   //------------------------------------------------------------



  //  std::vector< std::pair< std::vector<reco::TransientTrack>, TransientVertex> > selected_vertex;
  // Proto SelVtx(selected_vertex);/*!*/  
   /*****
   for(unsigned int  idx_reco_vertex= 0; idx_reco_vertex < reco_vertex.size(); idx_reco_vertex++){
     if( reco_vertex[idx_reco_vertex].second.normalisedChiSquared() < 100 &&  reco_vertex[idx_reco_vertex].second.originalTracks().size() > 3) 
       {
         selected_vertex.push_back(reco_vertex[idx_reco_vertex]);
         tree_secondaryVtx_iterative_isSelected.push_back(true);  
       } 
       else 
       {
         tree_secondaryVtx_iterative_isSelected.push_back(false); 
       }
     
       tree_secondaryVtx_iterative_X.push_back(reco_vertex[idx_reco_vertex].second.position().x()); 
       tree_secondaryVtx_iterative_Y.push_back(reco_vertex[idx_reco_vertex].second.position().y()); 
       tree_secondaryVtx_iterative_Z.push_back(reco_vertex[idx_reco_vertex].second.position().z()); 
       tree_secondaryVtx_iterative_nTracks.push_back(reco_vertex[idx_reco_vertex].second.originalTracks().size());
       tree_secondaryVtx_iterative_NChi2.push_back(reco_vertex[idx_reco_vertex].second.normalisedChiSquared() ); 
       //------------------------------------------------------------
       //check if reco vertex matches "reconstructable" gen vertex 
       //------------------------------------------------------------
       bool match_top = false; 
       if(!runOnData_){
          if(displacedVertex_top1_general.isValid() ){
            if(  
              fabs(reco_vertex[idx_reco_vertex].second.position().x() - displacedVertex_top1_general.position().x())  < 0.1 &&
              fabs(reco_vertex[idx_reco_vertex].second.position().y() - displacedVertex_top1_general.position().y())  < 0.1 &&
              fabs(reco_vertex[idx_reco_vertex].second.position().z() - displacedVertex_top1_general.position().z())  < 0.5
              //fabs(reco_vertex[idx_reco_vertex].second.position().x() - tree_gen_secondartVtxTop_X[0])  < 0.1 &&
              //fabs(reco_vertex[idx_reco_vertex].second.position().y() - tree_gen_secondartVtxTop_Y[0])  < 0.1 &&
              //fabs(reco_vertex[idx_reco_vertex].second.position().z() - tree_gen_secondartVtxTop_Z[0])  < 0.5  
            ) match_top = true;
       }
       if(displacedVertex_top2_general.isValid() ){
       if(  
         fabs(reco_vertex[idx_reco_vertex].second.position().x() - displacedVertex_top2_general.position().x())  < 0.1 &&
         fabs(reco_vertex[idx_reco_vertex].second.position().y() - displacedVertex_top2_general.position().y())  < 0.1 &&
         fabs(reco_vertex[idx_reco_vertex].second.position().z() - displacedVertex_top2_general.position().z())  < 0.5 
         //fabs(reco_vertex[idx_reco_vertex].second.position().x() - tree_gen_secondartVtxTop_X[1])  < 0.1 &&
         //fabs(reco_vertex[idx_reco_vertex].second.position().y() - tree_gen_secondartVtxTop_Y[1])  < 0.1 &&
         //fabs(reco_vertex[idx_reco_vertex].second.position().z() - tree_gen_secondartVtxTop_Z[1])  < 0.5  
         ) match_top = true;
       }
    }
    tree_secondaryVtx_iterative_match_top.push_back(match_top);
   }
   ****/ 
   
   
   //  cout << "-----------------"<<endl; 
   //  cout << selected_vertex.size() <<" number of vertex reco after quality cuts  " << endl;
  //  if(showlog)   cout << "number of vertex reco after quality cuts  " << selected_vertex.size() << endl;
  //  if(showlog)   cout << "number of vertex reco after quality cuts  " << SelVtx.Size() << endl;/*!*/
   
   /****for(unsigned int  idx_proto_vertex= 0; idx_proto_vertex < selected_vertex.size(); idx_proto_vertex++){
    //dump information for validation
    if(showlog) cout << "------------------------" <<endl;
    if(showlog) std::cout << " Interative Reco  TOP DISPLACED VERTEX X POS      "  << selected_vertex[idx_proto_vertex].second.position().x() << std::endl;
    if(showlog) std::cout << " Interative Reco  TOP DISPLACED VERTEX Y POS      "  << selected_vertex[idx_proto_vertex].second.position().y() << std::endl;
    if(showlog) std::cout << " Interative Reco  TOP DISPLACED VERTEX Z POS      "  << selected_vertex[idx_proto_vertex].second.position().z() << std::endl; 
    if(showlog) std::cout << " Interative Reco  TOP DISPLACED VERTEX CHI2/ndof  "  << selected_vertex[idx_proto_vertex].second.normalisedChiSquared()<< std::endl; 
    if(showlog) std::cout << " Interative Reco NUMBER OF ORIGINAL TRACKS        "  << selected_vertex[idx_proto_vertex].second.originalTracks().size() << std::endl;
  }****/

   /*int n_selected_match_top = 0; 
     for (unsigned int i =0; i<selected_vertex.size(); i++)
     {
     bool selected_match_top=false; 
     if( fabs(selected_vertex[i].second.position().x() - tree_gen_secondartVtxTop_X[0])  < 0.1 &&
     fabs(selected_vertex[i].second.position().y() - tree_gen_secondartVtxTop_Y[0])  < 0.1 &&
     fabs(selected_vertex[i].second.position().z() - tree_gen_secondartVtxTop_Z[0])  < 0.5  
     ) selected_match_top = true;
     if(   fabs(selected_vertex[i].second.position().x() - tree_gen_secondartVtxTop_X[1])  < 0.1 &&
     fabs(selected_vertex[i].second.position().y() - tree_gen_secondartVtxTop_Y[1])  < 0.1 &&
     fabs(selected_vertex[i].second.position().z() - tree_gen_secondartVtxTop_Z[1])  < 0.5  
     ) selected_match_top = true; 
     if (selected_match_top==true) n_selected_match_top++; 
     
     tree_secondaryVtx_iterative_step4_NChi2.push_back(selected_vertex[i].second.normalisedChiSquared()); 
     tree_secondaryVtx_iterative_step4_nTracks.push_back(selected_vertex[i].second.originalTracks().size()); 
     tree_secondaryVtx_iterative_step4_match_top.push_back(selected_match_top); 
     tree_secondaryVtx_iterative_step4_x.push_back(selected_vertex[i].second.position().x());
     tree_secondaryVtx_iterative_step4_y.push_back(selected_vertex[i].second.position().y());
     tree_secondaryVtx_iterative_step4_z.push_back(selected_vertex[i].second.position().z());
     
     }
     
     cout << n_selected_match_top << " number of good vertex reco after quality cuts  " << endl;
   */  
   
   ////////HT FILTER CHECK////////

   double HT_val = 0;
   bool MuonSup10GeV = false; 
   bool MuonSup28GeV = false; 
   
   for (unsigned int ij=0; ij<ak4slimmedJets->size(); ij++){
     const Jet& jet = ak4slimmedJets->at(ij);
     if ( fabs(jet.eta()) < 2.4 && jet.pt() > 20. ) HT_val += jet.pt();
   }
   
   int NumberMuonsSup10GeV = 0; 
   int NumberMuonsSup28GeV = 0; 
   float mu_mass = 0.1057;
   bool invMassGood = false; 
   TLorentzVector v1, v2, v;
   
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
   
   if( HT_val > 180. && MuonSup10GeV && MuonSup28GeV && invMassGood ) tree_passesHTFilter.push_back(true); 
   else tree_passesHTFilter.push_back(false);

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

  tree_passesHTFilter.clear(); 

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
   tree_track_stopReason.clear();
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
   tree_track_firsthit_phi.clear(); 
   tree_track_ntrk10.clear();
   tree_track_ntrk20.clear();
   tree_track_ntrk30.clear();

   tree_track_MVAval.clear(); 

   tree_track_recoVertex_idx.clear();
   tree_track_recoAK4SlimmedJet_idx.clear();
   tree_track_recoAK4PFJet_idx.clear();
   tree_track_reco08Jet_idx.clear(); 
   tree_track_recoCaloJet_idx.clear(); 
//    tree_track_reco08CaloJet_idx.clear(); 
   
   tree_track_nSimHits.clear();			 
   tree_track_isSimMatched.clear();		 

   tree_track_simtrack_charge.clear();		 
   tree_track_simtrack_pt.clear();  		 
   tree_track_simtrack_eta.clear();		 
   tree_track_simtrack_phi.clear(); 		 
   tree_track_simtrack_longLived.clear();		 
   //tree_track_simtrack_matchedHit .clear(); 	 
   tree_track_simtrack_pdgId.clear();		 
   tree_track_simtrack_numberOfTrackerHits.clear();
   tree_track_simtrack_numberOfTrackerLayers.clear();
   tree_track_simtrack_mass.clear();		 
   tree_track_simtrack_status.clear();		

   tree_track_genVertexPos_X.clear();
   tree_track_genVertexPos_Y.clear();
   tree_track_genVertexPos_Z.clear();
   tree_track_simtrack_llp1_dV.clear();
   tree_track_simtrack_llp2_dV.clear();
   tree_track_simtrack_isFromLLP.clear();	 
   tree_track_simtrack_isFromDispTop.clear();

   tree_simtrack_simtrack_charge.clear();
   tree_simtrack_simtrack_pt.clear();
   tree_simtrack_simtrack_eta.clear();
   tree_simtrack_simtrack_phi.clear();
   tree_simtrack_simtrack_longLived.clear();

   tree_simtrack_simtrack_pdgId.clear();
   tree_simtrack_simtrack_numberOfTrackerHits.clear();
   tree_simtrack_simtrack_numberOfTrackerLayers.clear();
   tree_simtrack_simtrack_mass.clear();
   tree_simtrack_simtrack_status.clear();

   tree_simtrack_genVertexPos_X.clear();
   tree_simtrack_genVertexPos_Y.clear();
   tree_simtrack_genVertexPos_Z.clear();

   tree_simtrack_isRecoMatched.clear();
   tree_simtrack_pca_dxy.clear();
   tree_simtrack_pca_dz.clear();
   tree_simtrack_trkIdx.clear();
   
   tree_genTop_X.clear(); 
   tree_genTop_Y.clear(); 
   tree_genTop_Z.clear(); 
   tree_genTop_charge.clear(); 

   tree_vtx_PosX.clear();
   tree_vtx_PosY.clear();
   tree_vtx_PosZ.clear();
   tree_vtx_NChi2.clear(); 

   tree_vtx_PosXError.clear();
   tree_vtx_PosYError.clear();
   tree_vtx_PosZError.clear();

   tree_secondaryVtx_X.clear(); 
   tree_secondaryVtx_Y.clear(); 
   tree_secondaryVtx_Z.clear(); 
   tree_secondaryVtx_diff_X.clear(); 
   tree_secondaryVtx_diff_Y.clear(); 
   tree_secondaryVtx_diff_Z.clear(); 
   tree_secondaryVtx_nTracks.clear(); 
   tree_secondaryVtx_isValid.clear(); 
   
//    tree_seedVtx_X.clear(); 
//    tree_seedVtx_Y.clear(); 
//    tree_seedVtx_Z.clear(); 
//    tree_seedVtx_dd.clear(); 
//    tree_seedVtx_dphi.clear(); 
//    tree_seedVtx_distance2track.clear(); 
//    tree_seedVtx_normChi2.clear(); 

//    tree_secondaryVtx_iterative_X.clear(); 
//    tree_secondaryVtx_iterative_Y.clear(); 
//    tree_secondaryVtx_iterative_Z.clear(); 
//    tree_secondaryVtx_iterative_nTracks.clear();
//    tree_secondaryVtx_iterative_NChi2.clear(); 
//    tree_secondaryVtx_iterative_match_top.clear(); 
//    tree_secondaryVtx_iterative_isSelected.clear(); 

   tree_AK4Slimmedjet_E.clear(); 
   tree_AK4Slimmedjet_pt.clear();
   tree_AK4Slimmedjet_eta.clear();
   tree_AK4Slimmedjet_phi.clear();
   tree_AK4Slimmedjet_idxTrack.clear();

   tree_AK4PFjet_E.clear(); 
   tree_AK4PFjet_pt.clear();
   tree_AK4PFjet_eta.clear();
   tree_AK4PFjet_phi.clear();
   tree_AK4PFjet_idxTrack.clear();

   tree_CaloJet_E.clear(); 
   tree_CaloJet_pt.clear();
   tree_CaloJet_eta.clear();
   tree_CaloJet_phi.clear();
   tree_CaloJet_idxTrack.clear(); 

   tree_jet08_E.clear(); 
   tree_jet08_pt.clear();
   tree_jet08_eta.clear();
   tree_jet08_phi.clear();
   tree_jet08_idxTrack.clear(); 

//    tree_CaloJet08_E.clear(); 
//    tree_CaloJet08_pt.clear();
//    tree_CaloJet08_eta.clear();
//    tree_CaloJet08_phi.clear();
//    tree_CaloJet08_idxTrack.clear(); 

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
   tree_genParticle_mother_pt.clear();
   tree_genParticle_mother_eta.clear();
   tree_genParticle_mother_phi.clear();

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

  //added by Paul

  tree_NSim.clear();
  tree_RecoEfficiency.clear();
  tree_simtrack_isRecoMatched_pt.clear();
  tree_FakeRate.clear();

//   tree_DVertex_top1_nTrks.clear();
//   tree_DVertex_top2_nTrks.clear();
//   tree_DVertex_nTrks.clear();

  tree_seedVtx_dz_top1.clear();
  tree_seedVtx_dz_top2.clear();
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

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
      //// GETTOKENDECLARATIONS 
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
      //jet information 
      //------------------------------------
      const edm::EDGetTokenT<edm::View<reco::Jet> > ak4slimmedJetToken_;
      const edm::EDGetTokenT<edm::View<reco::Jet> > ak4PFJetToken_;
      const edm::EDGetTokenT<edm::View<reco::Jet> > CaloJetToken_;
      
      const edm::EDGetTokenT<edm::View<reco::Jet> > ak8jetToken_;
      //   const edm::EDGetTokenT<edm::View<reco::Jet> > ak8CaloJetToken_;
      
      //------------------------------------
      //gen information
      //------------------------------------
      const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
      const edm::EDGetTokenT<reco::GenJetCollection>      genJetToken_;
      const edm::EDGetTokenT<reco::GenJetCollection>      ak8GenJetToken_;
      const edm::EDGetTokenT<reco::GenJetCollection>      genTTXJetsToken_;
      const edm::EDGetTokenT<GenEventInfoProduct>         genEventInfoToken_;
      const edm::EDGetTokenT<LHEEventProduct>            LHEEventProductToken_;
      
      //pat information
      const edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;
      
      std::string parametersDefinerName_;
      
      //------------------------------------
      //electrons
      //------------------------------------
      const edm::EDGetTokenT<pat::ElectronCollection> electronPATToken_;
      //------------------------------------
      //muons 
      //------------------------------------
      const edm::EDGetTokenT<pat::MuonCollection> slimmedmuonToken_;
      //$$  const edm::EDGetTokenT<reco::MuonCollection> recomuonToken_;
      //------------------------------------
      //pet MET
      //------------------------------------
      const edm::EDGetTokenT<pat::METCollection> metToken_;
      
      
      //------------------------------------
      //Zmu skim bit
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
      std::vector<int>      tree_track_nSimHits   ; 
      std::vector<bool>     tree_track_isSimMatched;
      //$$  std::vector<int>      tree_track_nPixel;
      //$$  std::vector<int>      tree_track_nStrip;
      std::vector< float >  tree_track_vx;
      std::vector< float >  tree_track_vy; 
      std::vector< float >  tree_track_vz;
      std::vector<float>    tree_track_firsthit_X; 
      std::vector<float>    tree_track_firsthit_Y; 
      std::vector<float>    tree_track_firsthit_Z;
      std::vector<float>    tree_track_firsthit_phi; 
      
      std::vector< int >    tree_track_simtrack_charge;		       
      std::vector< float >  tree_track_simtrack_pt; 	       
      std::vector< float >  tree_track_simtrack_eta  ;	       
      std::vector< float >  tree_track_simtrack_phi  ;	       
      std::vector<bool>     tree_track_simtrack_longLived	  ;
      std::vector<bool>     tree_track_simtrack_isFromDispTop;
      std::vector<int>      tree_track_simtrack_isFromLLP;
      // std::vector<int>      tree_track_simtrack_matchedHit	 ;
      std::vector<int>      tree_track_simtrack_pdgId;		 
      std::vector<int>      tree_track_simtrack_numberOfTrackerHits  ;
      std::vector<int>      tree_track_simtrack_numberOfTrackerLayers;
      std::vector<float>    tree_track_simtrack_mass  ;		       
      std::vector<int>      tree_track_simtrack_status;		       
      
      std::vector<float>    tree_track_genVertexPos_X;		  
      std::vector<float>    tree_track_genVertexPos_Y;		  
      std::vector<float>    tree_track_genVertexPos_Z;		  
      
      std::vector<int>      tree_track_recoVertex_idx; 
      std::vector<int>      tree_track_recoAK4SlimmedJet_idx; 
      std::vector<int>      tree_track_recoAK4PFJet_idx; 
      std::vector<int>      tree_track_reco08Jet_idx; 
      std::vector<int>      tree_track_recoCaloJet_idx; 
      //   std::vector<int>      tree_track_reco08CaloJet_idx; 
      
      std::vector<double>   tree_track_MVAval_FromDispTop; 
      
      int runNumber,eventNumber,lumiBlock; 
      
      //--------------------------------
      // Tracking Particle infos ------- 
      //--------------------------------
      
      std::vector< int >           tree_simtrack_simtrack_charge;	              
      std::vector< float >         tree_simtrack_simtrack_pt;	              
      std::vector< float >         tree_simtrack_simtrack_eta;	      
      std::vector< float >         tree_simtrack_simtrack_phi;	      
      std::vector<bool>            tree_simtrack_simtrack_longLived;
      std::vector<int>             tree_simtrack_simtrack_pdgId;	       
      std::vector<int>             tree_simtrack_simtrack_numberOfTrackerHits;
      std::vector<int>             tree_simtrack_simtrack_numberOfTrackerLayers;
      std::vector<float>           tree_simtrack_simtrack_mass;	              
      std::vector<int>             tree_simtrack_simtrack_status;	              
      
      std::vector<float>           tree_simtrack_genVertexPos_X;                
      std::vector<float>           tree_simtrack_genVertexPos_Y;                
      std::vector<float>           tree_simtrack_genVertexPos_Z;        
      std::vector<bool>            tree_simtrack_isRecoMatched;
      std::vector<float>           tree_simtrack_pca_dxy;
      std::vector<float>           tree_simtrack_pca_dz   ; 
      std::vector<std::vector<int> > tree_simtrack_trkIdx;  
      std::vector<float>           tree_simtrack_isRecoMatched_pt;
      //--------------------------------
      // Beam spot infos ------- 
      //--------------------------------
      float	      tree_bs_PosX;	   
      float	      tree_bs_PosY;	   
      float	      tree_bs_PosZ; 
      
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
      
      std::vector<float> tree_secondaryVtx_iterative_X; 
      std::vector<float> tree_secondaryVtx_iterative_Y; 
      std::vector<float> tree_secondaryVtx_iterative_Z; 
      std::vector<int>   tree_secondaryVtx_iterative_nTracks;
      std::vector<float> tree_secondaryVtx_iterative_NChi2; 
      std::vector<bool>  tree_secondaryVtx_iterative_match_top;
      std::vector<bool>  tree_secondaryVtx_iterative_isSelected; 


      std::vector<float> tree_genTop_X; 
      std::vector<float> tree_genTop_Y; 
      std::vector<float> tree_genTop_Z; 
      std::vector<float> tree_genTop_charge; 


      std::vector<float> tree_seedVtx_X; 
      std::vector<float> tree_seedVtx_Y; 
      std::vector<float> tree_seedVtx_Z; 
      std::vector<float> tree_seedVtx_dd; 
      std::vector<float> tree_seedVtx_dphi; 
      std::vector<float> tree_seedVtx_distance2track; 
      std::vector<float> tree_seedVtx_normChi2; 
      

      //added by Paul 

      /* about generated LLPs (top is used as the physics process first studied is lookinf for tops coming from neutralinos)*/
      std::vector<float> tree_top1_x;
      std::vector<float> tree_top1_y;
      std::vector<float> tree_top1_z;
      std::vector<float> tree_top2_x;
      std::vector<float> tree_top2_y;
      std::vector<float> tree_top2_z;

      std::vector<float> tree_top1_pt;
      std::vector<float> tree_top1_eta;
      std::vector<float> tree_top1_phi;

      std::vector<float> tree_top2_pt;
      std::vector<float> tree_top2_eta;
      std::vector<float> tree_top2_phi;
      std::vector<int> tree_ntop;

      std::vector<float> tree_simtrack_genVertexPosDiff_top1_X;
      std::vector<float> tree_simtrack_genVertexPosDiff_top1_Y;
      std::vector<float> tree_simtrack_genVertexPosDiff_top1_Z;
      std::vector<float> tree_simtrack_genVertexPosDiff_top1_dV;
      std::vector<float> tree_simtrack_genVertexPosDiff_top2_X;
      std::vector<float> tree_simtrack_genVertexPosDiff_top2_Y;
      std::vector<float> tree_simtrack_genVertexPosDiff_top2_Z;
      std::vector<float> tree_simtrack_genVertexPosDiff_top2_dV;

      std::vector<int> tree_NSimtoReco;
      std::vector<int> tree_NSim;
      std::vector<float> tree_Ratio_NStRNS;
      std::vector<int> tree_nTracks_ADsel;
      std::vector<int> tree_nTracks_b4sel;
      std::vector<float> tree_Ratio_b4ADsel;
      std::vector<float> tree_displacedTTracks_pt;
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
      
      std::vector<float> tree_PFMet_et; 
      std::vector<float> tree_PFMet_phi; 
      std::vector<float> tree_PFMet_sig; 
      
      //--------------------------------
      // gen infos ------- 
      //-------------------------------- 
      
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
  
  //event info 
  
  smalltree->Branch("runNumber",  &runNumber,  "runNumber/I");
  smalltree->Branch("eventNumber",&eventNumber,"eventNumber/I");
  smalltree->Branch("lumiBlock"  ,&lumiBlock,  "lumiBlock/I");
  smalltree->Branch("tree_bs_PosX", &tree_bs_PosX,  "tree_bs_PosX/F"  );
  smalltree->Branch("tree_bs_PosY", &tree_bs_PosY,  "tree_bs_PosY/F" );
  smalltree->Branch("tree_bs_PosZ", &tree_bs_PosZ,  "tree_bs_PosZ/F"  );
  
  //primary vertex info
  
  smalltree->Branch("tree_vtx_PosX", &tree_vtx_PosX);        
  smalltree->Branch("tree_vtx_PosY", &tree_vtx_PosY);        
  smalltree->Branch("tree_vtx_PosZ", &tree_vtx_PosZ); 
  smalltree->Branch("tree_vtx_NChi2", &tree_vtx_NChi2); 
  
  smalltree->Branch("tree_vtx_PosXError", &tree_vtx_PosXError);        
  smalltree->Branch("tree_vtx_PosYError", &tree_vtx_PosYError);        
  smalltree->Branch("tree_vtx_PosZError", &tree_vtx_PosZError); 
  
  //trigger info
  smalltree->Branch("tree_trigger_names", &tree_trigger_names);
  smalltree->Branch("tree_trigger_bits",  &tree_trigger_bits);
  
  smalltree->Branch("tree_NbrOfZCand",  &tree_NbrOfZCand,  "tree_NbrOfZCand/I");
  smalltree->Branch("tree_passesHTFilter", &tree_passesHTFilter); 
  
  //track
  
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
   
  //$$  smalltree->Branch("tree_track_nPixel",       &tree_track_nPixel );
  //$$  smalltree->Branch("tree_track_nStrip",       &tree_track_nStrip );
  smalltree->Branch("tree_track_vx",           &tree_track_vx ); 
  smalltree->Branch("tree_track_vy",           &tree_track_vy ); 
  smalltree->Branch("tree_track_vz",           &tree_track_vz ); 
  smalltree->Branch("tree_track_firsthit_X",   &tree_track_firsthit_X); 
  smalltree->Branch("tree_track_firsthit_Y",   &tree_track_firsthit_Y);
  smalltree->Branch("tree_track_firsthit_Z",   &tree_track_firsthit_Z);
  smalltree->Branch("tree_track_firsthit_phi", &tree_track_firsthit_phi); 

  //info about the simulated track matched to the reco track
   
  smalltree->Branch("tree_track_simtrack_charge",		 &tree_track_simtrack_charge );		      
  smalltree->Branch("tree_track_simtrack_pt",  		         &tree_track_simtrack_pt );	      
  smalltree->Branch("tree_track_simtrack_eta", 		         &tree_track_simtrack_eta  );  	      
  smalltree->Branch("tree_track_simtrack_phi", 		         &tree_track_simtrack_phi  );  	      
  smalltree->Branch("tree_track_simtrack_longLived",		 &tree_track_simtrack_longLived );
  smalltree->Branch("tree_track_simtrack_pdgId",		 &tree_track_simtrack_pdgId ); 
  smalltree->Branch("tree_track_simtrack_isFromLLP",	         &tree_track_simtrack_isFromLLP );
  smalltree->Branch("tree_track_simtrack_isFromDispTop",	         &tree_track_simtrack_isFromDispTop );
  smalltree->Branch("tree_track_simtrack_numberOfTrackerHits",   &tree_track_simtrack_numberOfTrackerHits   );
  smalltree->Branch("tree_track_simtrack_numberOfTrackerLayers", &tree_track_simtrack_numberOfTrackerLayers );
  smalltree->Branch("tree_track_simtrack_mass",		         &tree_track_simtrack_mass   );		      
  smalltree->Branch("tree_track_simtrack_status",		 &tree_track_simtrack_status );		      
 
  smalltree->Branch("tree_track_genVertexPos_X", &tree_track_genVertexPos_X ); 		  
  smalltree->Branch("tree_track_genVertexPos_Y", &tree_track_genVertexPos_Y ); 		  
  smalltree->Branch("tree_track_genVertexPos_Z", &tree_track_genVertexPos_Z ); 		  
  smalltree->Branch("tree_track_recoVertex_idx", &tree_track_recoVertex_idx); 

  smalltree->Branch("tree_track_recoAK4SlimmedJet_idx", &tree_track_recoAK4SlimmedJet_idx); 
  smalltree->Branch("tree_track_recoAK4PFJet_idx",      &tree_track_recoAK4PFJet_idx); 
  smalltree->Branch("tree_track_reco08Jet_idx",         &tree_track_reco08Jet_idx); 
  smalltree->Branch("tree_track_recoCaloJet_idx",       &tree_track_recoCaloJet_idx); 
  //   smalltree->Branch("tree_track_reco08CaloJet_idx",     &tree_track_reco08CaloJet_idx); 
  
  smalltree->Branch("tree_track_MVAval_FromDispTop", &tree_track_MVAval_FromDispTop); 
  
  smalltree->Branch("tree_genTop_X", &tree_genTop_X); 
  smalltree->Branch("tree_genTop_Y", &tree_genTop_Y); 
  smalltree->Branch("tree_genTop_Z", &tree_genTop_Z); 
  smalltree->Branch("tree_genTop_charge", &tree_genTop_charge); 


  smalltree->Branch("tree_seedVtx_X", &tree_seedVtx_X); 
  smalltree->Branch("tree_seedVtx_Y", &tree_seedVtx_Y); 
  smalltree->Branch("tree_seedVtx_Z", &tree_seedVtx_Z); 
  smalltree->Branch("tree_seedVtx_dd",             &tree_seedVtx_dd); 
  smalltree->Branch("tree_seedVtx_dphi",           &tree_seedVtx_dphi); 
  smalltree->Branch("tree_seedVtx_distance2track", &tree_seedVtx_distance2track); 
  smalltree->Branch("tree_seedVtx_normChi2",       &tree_seedVtx_normChi2); 
      
      
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
  
  smalltree->Branch("tree_secondaryVtx_iterative_X"      ,    &tree_secondaryVtx_iterative_X); 
  smalltree->Branch("tree_secondaryVtx_iterative_Y"      ,    &tree_secondaryVtx_iterative_Y); 
  smalltree->Branch("tree_secondaryVtx_iterative_Z"      ,    &tree_secondaryVtx_iterative_Z); 
  smalltree->Branch("tree_secondaryVtx_iterative_nTracks",    &tree_secondaryVtx_iterative_nTracks);
  smalltree->Branch("tree_secondaryVtx_iterative_NChi2"  ,    &tree_secondaryVtx_iterative_NChi2); 
  smalltree->Branch("tree_secondaryVtx_iterative_isSelected", &tree_secondaryVtx_iterative_isSelected);  
  
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

  smalltree->Branch("tree_simtrack_isRecoMatched_pt",&tree_simtrack_isRecoMatched_pt);/*!*/

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
  
  ////met info 
  smalltree->Branch("tree_PFMet_et" ,  &tree_PFMet_et); 
  smalltree->Branch("tree_PFMet_phi" , &tree_PFMet_phi); 
  smalltree->Branch("tree_PFMet_sig" , &tree_PFMet_sig); 

  ////gen info 
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

  //genJet info 
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

   //electrons info 
  smalltree->Branch("tree_electron_pt"  ,   &tree_electron_pt);
  smalltree->Branch("tree_electron_eta" ,   &tree_electron_eta);
  smalltree->Branch("tree_electron_phi" ,   &tree_electron_phi);
  smalltree->Branch("tree_electron_vx"  ,   &tree_electron_vx);
  smalltree->Branch("tree_electron_vy" ,    &tree_electron_vy);
  smalltree->Branch("tree_electron_vz" ,    &tree_electron_vz); 
  smalltree->Branch("tree_electron_energy", &tree_electron_energy); 

   //muons info 
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

  smalltree->Branch("tree_top1_x",&tree_top1_x);
  smalltree->Branch("tree_top1_y",&tree_top1_y);
  smalltree->Branch("tree_top1_z",&tree_top1_z);    
  smalltree->Branch("tree_top2_x",&tree_top2_x);
  smalltree->Branch("tree_top2_y",&tree_top2_y);
  smalltree->Branch("tree_top2_z",&tree_top2_z);

  smalltree->Branch("tree_top1_pt",&tree_top1_pt);
  smalltree->Branch("tree_top1_eta",&tree_top1_eta);
  smalltree->Branch("tree_top1_phi",&tree_top1_phi);

  smalltree->Branch("tree_top2_pt",&tree_top2_pt);
  smalltree->Branch("tree_top2_eta",&tree_top2_eta);
  smalltree->Branch("tree_top2_phi",&tree_top2_phi);
  smalltree->Branch("tree_ntop",&tree_ntop);

  smalltree->Branch("tree_simtrack_genVertexPosDiff_top1_X",&tree_simtrack_genVertexPosDiff_top1_X);
  smalltree->Branch("tree_simtrack_genVertexPosDiff_top1_X",&tree_simtrack_genVertexPosDiff_top1_Y);
  smalltree->Branch("tree_simtrack_genVertexPosDiff_top1_X",&tree_simtrack_genVertexPosDiff_top1_Z);
  smalltree->Branch("tree_simtrack_genVertexPosDiff_top1_dV",&tree_simtrack_genVertexPosDiff_top1_dV);
  smalltree->Branch("tree_simtrack_genVertexPosDiff_top2_X",&tree_simtrack_genVertexPosDiff_top2_X);
  smalltree->Branch("tree_simtrack_genVertexPosDiff_top2_Y",&tree_simtrack_genVertexPosDiff_top2_Y);
  smalltree->Branch("tree_simtrack_genVertexPosDiff_top2_Z",&tree_simtrack_genVertexPosDiff_top2_Z);
  smalltree->Branch("tree_simtrack_genVertexPosDiff_top2_dV",&tree_simtrack_genVertexPosDiff_top2_dV);

  smalltree->Branch("tree_NSimtoReco",&tree_NSimtoReco);
  smalltree->Branch("tree_NSim",&tree_NSim);
  smalltree->Branch("tree_Ratio_NStRNS",&tree_Ratio_NStRNS);
  smalltree->Branch("tree_nTracks_ADsel",&tree_nTracks_ADsel);
  smalltree->Branch("tree_nTracks_b4sel",&tree_nTracks_b4sel);
  smalltree->Branch("tree_Ratio_b4ADsel",&tree_Ratio_b4ADsel);
  smalltree->Branch("tree_displacedTTracks_pt",&tree_displacedTTracks_pt);

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
   
   //    edm::ESHandle<TrackerTopology> tTopoHandle;
   //    iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
   //    const TrackerTopology* const tTopo = tTopoHandle.product();
   
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
   ////////      BS     /////////////
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
   //    std::map<size_t , int > ak8CaloJetToTrackMap;
   std::map<size_t , int > trackToPFJetMap;
   
   ///// TRACK ASSOCIATION TO VERTICES AND JETS 
   
   //    int nRecoTracks = tracks.size(); 
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
       else             trackToAK4SlimmedJetMap[idxTrack] = -1;
       
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
       else                 trackToAK4PFJetMap[idxTrack] = -1;
       
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
       else               trackTo08JetMap[idxTrack] = -1;
       
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
       else                  trackToCaloJetMap[idxTrack] = -1;
       
       //       int idxCaloJet08=0;
       //       bool found_match_calo08 = false;
       //       for(unsigned int ij=0;ij<ak8CaloJets ->size();ij++){
       //           const Jet& jet = ak8CaloJets->at(ij);
       //           if(jet.pt() < 20) continue; 
       //           TLorentzVector jet4m, track4m;
       //           jet4m.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0);
       //           track4m.SetPtEtaPhiM(tracks[i].pt(), tracks[i].eta(), tracks[i].phi(), 0);
       //           if( jet4m.DeltaR(track4m) < 0.8){
       //           found_match_calo08 = true;
       //           break;
       //         }
       //         else idxCaloJet08++;
       //       }
       //       if(found_match_calo08) trackTo08CaloJetMap[idxTrack] = idxCaloJet08;
       //       else trackTo08CaloJetMap[idxTrack] = -1;
       
       idxTrack++;
     } 
   
   //////////////////////////////////
   //////////////////////////////////
   ////////   vertices  /////////////
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
   if ( PFMETs->size() > 0 ) {
     const pat::MET &themet = PFMETs->front();
     tree_PFMet_et.push_back(themet.et());
     tree_PFMet_phi.push_back(themet.phi()); 
     tree_PFMet_sig.push_back(themet.significance()); 
   }
   
   
   //////////////////////////////////
   //////////////////////////////////
   //////// gen Information /////////
   //////////////////////////////////
   //////////////////////////////////
   
   ///GEN PARTICLES 
   
   //$$
   int ntop = 0;
   float top1_x = 0, top1_y = 0, top1_z = 0; //LLP1
   float top2_x = 0, top2_y = 0, top2_z = 0; //LLP2 avoir ces informations là dnas le NTuple
   float top1_pt = -1, top1_eta= -10, top1_phi=-10;
   float top2_pt = -1, top2_eta= -10, top2_phi=-10;
  //  //Avoir LLP1_pt, eta, phi same forLLP2, "mom", x, y,z
   //for mom see below the loops
   //ajouter ntop (nLLP)

   // everything is done here
   //$$
   
   if ( !runOnData_ ) {
     //      int nGenParticles = genParticles->size(); 
     //cout << "number of gen particles "<<nGenParticles<<endl; 
     
     for (auto const & genParticle : *genParticles)
       {
	 //$$
	 if ( abs(genParticle.pdgId()) >= 1 && abs(genParticle.pdgId()) <= 6 ) { //Version de Daniel, less restrective on quarks
	   const Candidate * mom = genParticle.mother();
	   if ( abs(mom->pdgId()) == 1000023 ) { // quark from neutralino 
	     if ( ntop >= 2 ) {
	       float dV1 = (genParticle.vx() - top1_x)*(genParticle.vx() - top1_x)
		 + (genParticle.vy() - top1_y)*(genParticle.vy() - top1_y)
		 + (genParticle.vz() - top1_z)*(genParticle.vz() - top1_z);
	       float dV2 = (genParticle.vx() - top2_x)*(genParticle.vx() - top2_x)
		 + (genParticle.vy() - top2_y)*(genParticle.vy() - top2_y)
		 + (genParticle.vz() - top2_z)*(genParticle.vz() - top2_z);
	       if ( dV1 > 0.01 && dV2 > 0.01 ) ntop++;
	     } 
	     if ( ntop == 1 ) {
	       float dV = (genParticle.vx() - top1_x)*(genParticle.vx() - top1_x)
		 + (genParticle.vy() - top1_y)*(genParticle.vy() - top1_y)
		 + (genParticle.vz() - top1_z)*(genParticle.vz() - top1_z);
	       if ( dV > 0.01 ) {
		 ntop = 2;
		 top2_x = genParticle.vx();
		 top2_y = genParticle.vy();
		 top2_z = genParticle.vz();
     top2_pt = genParticle.pt();
     top2_eta = genParticle.eta();
     top2_phi =genParticle.phi();
	       } 
	     } 
	     if ( ntop == 0 ) {
	       ntop = 1;
	       top1_x = genParticle.vx();
	       top1_y = genParticle.vy();
	       top1_z = genParticle.vz();
         top1_pt = genParticle.pt();
         top1_eta = genParticle.eta();
         top1_phi =genParticle.phi();        
	     } 
	   } 
	 } 
	 //$$
	tree_top1_x.push_back(top1_x);
  tree_top1_y.push_back(top1_y);
  tree_top1_z.push_back(top1_z);
  tree_top2_x.push_back(top2_x);
  tree_top2_y.push_back(top2_y);
  tree_top2_z.push_back(top2_z);

  tree_top1_pt.push_back(top1_pt);
  tree_top1_eta.push_back(top1_eta);
  tree_top1_phi.push_back(top1_phi);

  tree_top2_pt.push_back(top2_pt);
  tree_top2_eta.push_back(top2_eta);
  tree_top2_phi.push_back(top2_phi);
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
	 tree_genParticle_mother_pdgId.push_back( mom ? mom->pdgId() :  -1 );
	 tree_genParticle_mother_pt.push_back( mom ? mom->pt() :  -1 );
	 tree_genParticle_mother_eta.push_back( mom ? mom->eta() :  -10 );
	 tree_genParticle_mother_phi.push_back( mom ? mom->phi() :  -10 );
       }
     
     //      cout << " ntop " << ntop << "   top1 " << top1_x << " " << top1_y << " " << top1_z  
     //                               << "   top2 " << top2_x << " " << top2_y << " " << top2_z << endl; 

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
  //  float muon2_x = 0, muon2_y = 0, muon2_z = 0; //LLP2 avoir ces informations là dnas le NTuple
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
       tree_slimmedmuon_dxy.push_back(	  muon.bestTrack()->dxy(bs));
       tree_slimmedmuon_dxyError.push_back( muon.bestTrack()->dxyError());
       tree_slimmedmuon_dz.push_back(	  muon.bestTrack()->dz(bs.position()));
       tree_slimmedmuon_dzError.push_back(  muon.bestTrack()->dzError());
       tree_slimmedmuon_charge.push_back(   muon.charge()); 
       
       tree_slimmedmuon_PFisoVeryTight.  push_back(   muon.passed(reco::Muon::PFIsoVeryTight));
       tree_slimmedmuon_PFisoTight.      push_back(   muon.passed(reco::Muon::PFIsoTight) );
       tree_slimmedmuon_PFisoMedium.     push_back(   muon.passed(reco::Muon::PFIsoMedium));
       tree_slimmedmuon_PFisoLoose.      push_back(   muon.passed(reco::Muon::PFIsoLoose ));
       tree_slimmedmuon_MVAisoLoose.     push_back(   muon.passed(reco::Muon::MvaLoose )  );
       tree_slimmedmuon_MVAisoMedium.    push_back(   muon.passed(reco::Muon::MvaMedium)  );
       tree_slimmedmuon_MVAisoTight.     push_back(   muon.passed(reco::Muon::MvaTight )  );
       tree_slimmedmuon_isGlobalMuon.    push_back(   muon.passed(muon.isGlobalMuon())	);
       tree_slimmedmuon_isStandAloneMuon.push_back(   muon.passed(muon.isStandAloneMuon()));
       tree_slimmedmuon_CutBasedIdLoose		   .push_back(   muon.passed(reco::Muon::CutBasedIdLoose));
       tree_slimmedmuon_CutBasedIdMedium     .push_back(   muon.passed(reco::Muon::CutBasedIdMedium));
       tree_slimmedmuon_CutBasedIdMediumPrompt	  .push_back(	muon.passed(reco::Muon::CutBasedIdMediumPrompt));
       tree_slimmedmuon_CutBasedIdTight		   .push_back(   muon.passed(reco::Muon::CutBasedIdTight));
     }
       std::cout<< "number of muons first method : " << tree_slimmedmuon_pt.size() << std::endl;

  //  for (auto const & muon : *slimmedmuons)
  //    {
  //      if ( muon.pt() < 3. ) continue;
  //      if  (nmuon >=2)
  //       {
  //         float dv1 = (muon.vx() - muon1_x)*(muon.vx() - muon1_x)
	// 	      + (muon.vy() - muon1_y)*(muon.vy() - muon1_y)
	// 	      + (muon.vz() - muon1_z)*(muon.vz() - muon1_z);
	//         float dv2 = (muon.vx() - muon2_x)*(muon.vx() - muon2_x)
	// 	      + (muon.vy() - muon2_y)*(muon.vy() - muon2_y)
	// 	      + (muon.vz() - muon2_z)*(muon.vz() - muon2_z);
	//         if ( dv1 > 0.01 && dv2 > 0.01 ) {nmuon++;}
  //         else {std::cout<< "muon discarded from selection "<< std::endl;}
	//       } 
	//      if ( nmuon == 1 ) 
  //       {
	//         float dv = (muon.vx() - muon1_x)*(muon.vx() - muon1_x)
	// 	      + (muon.vy() - muon1_y)*(muon.vy() - muon1_y)
	// 	      + (muon.vz() - muon1_z)*(muon.vz() - muon1_z);
	//         if ( dv > 0.01 ) 
  //           {
	// 	          nmuon = 2;
	// 	          muon2_x = muon.vx();
	// 	          muon2_y = muon.vy();
	// 	          muon2_z = muon.vz();
  //             muon2_pt = muon.pt();
  //             muon2_eta = muon.eta();
  //             muon2_phi = muon.phi();
  //             muon2_E = muon.energy();
  //             muon2_dxy = muon.bestTrack()->dxy(bs);
  //             muon2_dxyError =  muon.bestTrack()->dxyError();
  //             muon2_dz = muon.bestTrack()->dz(bs.position());
  //             muon2_dzError = muon.bestTrack()->dzError();
  //             muon2_charge = muon.charge();
  //             muon2_PFisoVeryTight =muon.passed(reco::Muon::PFIsoVeryTight);
  //             muon2_PFisoTight =muon.passed(reco::Muon::PFIsoTight) ;
  //             muon2_PFisoMedium =muon.passed(reco::Muon::PFIsoMedium);
  //             muon2_PFisoLoose =muon.passed(reco::Muon::PFIsoLoose );
  //             muon2_MVAisoLoose =muon.passed(reco::Muon::MvaLoose );
  //             muon2_MVAisoMedium =muon.passed(reco::Muon::MvaMedium);
  //             muon2_MVAisoTight =muon.passed(reco::Muon::MvaTight );
  //             muon2_isGlobalMuon =muon.passed(muon.isGlobalMuon());
  //             muon2_isStandAloneMuon =muon.passed(muon.isStandAloneMuon());
  //             muon2_CutBasedIdLoose =muon.passed(reco::Muon::CutBasedIdLoose);
  //             muon2_CutBasedIdMedium =muon.passed(reco::Muon::CutBasedIdMedium);
  //             muon2_CutBasedIdMediumPrompt =muon.passed(reco::Muon::CutBasedIdMediumPrompt);
  //             muon2_CutBasedIdTight =muon.passed(reco::Muon::CutBasedIdTight);
	//           }
  //           else {std::cout<< "muon discarded from selection "<< std::endl;} 
	//       } 
	//      if ( nmuon == 0 ) 
  //       {
	//         nmuon = 1;
	//         muon1_x = muon.vx();
	//         muon1_y = muon.vy();
	//         muon1_z = muon.vz();
  //         muon1_pt = muon.pt();
  //         muon1_eta = muon.eta();
  //         muon1_phi = muon.phi();
  //         muon1_E = muon.energy();
  //         muon1_dxy = muon.bestTrack()->dxy(bs);
  //         muon1_dxyError =  muon.bestTrack()->dxyError();
  //         muon1_dz = muon.bestTrack()->dz(bs.position());
  //         muon1_dzError = muon.bestTrack()->dzError();
  //         muon1_charge = muon.charge();
  //         muon1_PFisoVeryTight =muon.passed(reco::Muon::PFIsoVeryTight);
  //         muon1_PFisoTight =muon.passed(reco::Muon::PFIsoTight) ;
  //         muon1_PFisoMedium =muon.passed(reco::Muon::PFIsoMedium);
  //         muon1_PFisoLoose =muon.passed(reco::Muon::PFIsoLoose );
  //         muon1_MVAisoLoose =muon.passed(reco::Muon::MvaLoose );
  //         muon1_MVAisoMedium =muon.passed(reco::Muon::MvaMedium);
  //         muon1_MVAisoTight =muon.passed(reco::Muon::MvaTight );
  //         muon1_isGlobalMuon =muon.passed(muon.isGlobalMuon());
  //         muon1_isStandAloneMuon =muon.passed(muon.isStandAloneMuon());
  //         muon1_CutBasedIdLoose =muon.passed(reco::Muon::CutBasedIdLoose);
  //         muon1_CutBasedIdMedium =muon.passed(reco::Muon::CutBasedIdMedium);
  //         muon1_CutBasedIdMediumPrompt =muon.passed(reco::Muon::CutBasedIdMediumPrompt);
  //         muon1_CutBasedIdTight =muon.passed(reco::Muon::CutBasedIdTight);
	//       } 
	//    }
    
  // std::cout<< "number of muons second method : " << nmuon << std::endl;
  // //$$
	// tree_slimmedmuon1_vx.push_back(muon1_x);
  // tree_slimmedmuon1_vy.push_back(muon1_y);
  // tree_slimmedmuon1_vz.push_back(muon1_z);
  // tree_slimmedmuon1_pt.push_back(muon1_pt);
  // tree_slimmedmuon1_eta.push_back(muon1_eta);
  // tree_slimmedmuon1_phi.push_back(muon1_phi);

  // tree_slimmedmuon1_energy.push_back(   muon1_E);
  // tree_slimmedmuon1_dxy.push_back(muon1_dxy);
  // tree_slimmedmuon1_dxyError.push_back( muon1_dxyError);
  // tree_slimmedmuon1_dz.push_back(	  muon1_dz);
  // tree_slimmedmuon1_dzError.push_back(  muon1_dzError);
  // tree_slimmedmuon1_charge.push_back(   muon1_charge); 
       
  // tree_slimmedmuon1_PFisoVeryTight.  push_back(   muon1_PFisoVeryTight);
  // tree_slimmedmuon1_PFisoTight.      push_back(   muon1_PFisoTight );
  // tree_slimmedmuon1_PFisoMedium.     push_back(   muon1_PFisoMedium);
  // tree_slimmedmuon1_PFisoLoose.      push_back(  muon1_PFisoLoose);
  // tree_slimmedmuon1_MVAisoLoose.     push_back(   muon1_MVAisoLoose  );
  // tree_slimmedmuon1_MVAisoMedium.    push_back(   muon1_MVAisoMedium  );
  // tree_slimmedmuon1_MVAisoTight.     push_back(   muon1_MVAisoTight  );
  // tree_slimmedmuon1_isGlobalMuon.    push_back(   muon1_isGlobalMuon	);
  // tree_slimmedmuon1_isStandAloneMuon.push_back(   muon1_isStandAloneMuon);
  // tree_slimmedmuon1_CutBasedIdLoose		   .push_back(   muon1_CutBasedIdLoose);
  // tree_slimmedmuon1_CutBasedIdMedium     .push_back(   muon1_CutBasedIdMedium);
  // tree_slimmedmuon1_CutBasedIdMediumPrompt	  .push_back(	muon1_CutBasedIdMediumPrompt);
  // tree_slimmedmuon1_CutBasedIdTight		   .push_back(   muon1_CutBasedIdTight);  


  // tree_slimmedmuon2_x.push_back(muon2_x);
  // tree_slimmedmuon2_y.push_back(muon2_y);
  // tree_slimmedmuon2_z.push_back(muon2_z);
  // tree_slimmedmuon2_pt.push_back(muon2_pt);
  // tree_slimmedmuon2_eta.push_back(muon2_eta);
  // tree_slimmedmuon2_phi.push_back(muon2_phi);
  // tree_slimmedmuon2_energy.push_back(   muon2_E);
  // tree_slimmedmuon2_dxy.push_back(muon2_dxy);
  // tree_slimmedmuon2_dxyError.push_back( muon2_dxyError);
  // tree_slimmedmuon2_dz.push_back(	  muon2_dz);
  // tree_slimmedmuon2_dzError.push_back(  muon2_dzError);
  // tree_slimmedmuon2_charge.push_back(   muon2_charge); 
       
  // tree_slimmedmuon2_PFisoVeryTight.  push_back(   muon2_PFisoVeryTight);
  // tree_slimmedmuon2_PFisoTight.      push_back(   muon2_PFisoTight );
  // tree_slimmedmuon2_PFisoMedium.     push_back(   muon2_PFisoMedium);
  // tree_slimmedmuon2_PFisoLoose.      push_back(  muon2_PFisoLoose);
  // tree_slimmedmuon2_MVAisoLoose.     push_back(   muon2_MVAisoLoose  );
  // tree_slimmedmuon2_MVAisoMedium.    push_back(   muon2_MVAisoMedium  );
  // tree_slimmedmuon2_MVAisoTight.     push_back(   muon2_MVAisoTight  );
  // tree_slimmedmuon2_isGlobalMuon.    push_back(   muon2_isGlobalMuon	);
  // tree_slimmedmuon2_isStandAloneMuon.push_back(   muon2_isStandAloneMuon);
  // tree_slimmedmuon2_CutBasedIdLoose		   .push_back(   muon2_CutBasedIdLoose);
  // tree_slimmedmuon2_CutBasedIdMedium     .push_back(   muon2_CutBasedIdMedium);
  // tree_slimmedmuon2_CutBasedIdMediumPrompt	  .push_back(	muon2_CutBasedIdMediumPrompt);
  // tree_slimmedmuon2_CutBasedIdTight		   .push_back(   muon2_CutBasedIdTight);
	 //$$


   //reco::BeamSpot vertexBeamSpot= *recoBeamSpotHandle;
   edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);
   
   /////////////////
   //prepare association to tracks by hit
   reco::RecoToSimCollection recSimColl ;
   if(!runOnData_) recSimColl= associatorByHits.associateRecoToSim(trackRefs, tpCollection);
   
   int nTracks = 0; 
   //    int nUnmatchTrack_fromPU = 0;
   //    int nUnmatchTrack_fromPV = 0;
   //    int nUnmatchTrack= 0;
   //    int nPUTrack= 0;
   
   //---------------
   //loops on tracks 
   //---------------
   
   for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack) 
     {
       const auto& itTrack = trackRefs[iTrack];
       //------------------------
       //general track properties
       //------------------------
       //std::cout << "  " << itTrack->pt()  << std::endl;
       if(itTrack->pt() < 1) continue;
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
       
       tree_track_dxy.push_back(  	    itTrack->dxy(bs.position()));
       tree_track_dxyError.push_back(	    itTrack->dxyError());
       tree_track_dz.push_back(		    itTrack->dz(bs.position()));
       tree_track_dzError.push_back(	    itTrack->dzError());
       tree_track_numberOfLostHits.push_back( itTrack->numberOfLostHits());
       tree_track_numberOfValidHits.push_back(itTrack->numberOfValidHits());
       
       tree_track_originalAlgo.push_back(itTrack->originalAlgo());
       tree_track_algo.push_back(itTrack->algo());
       tree_track_stopReason.push_back(itTrack->stopReason());
       
       //--------------------------------
       //general hit properties of tracks
       //--------------------------------
       
       const reco::HitPattern& hp = itTrack->hitPattern();
       //$$     tree_track_nPixel   .push_back(hp.numberOfValidPixelHits());
       //$$     tree_track_nStrip   .push_back(hp.numberOfValidStripHits());
       
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
         float NReco=0;
         float NRecotoSim;
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
         
         int   simtrack_charge		    =-1000;
         float simtrack_pt		    =-1;
         float simtrack_eta		    =-1000;
         float simtrack_phi		    =-1000;
         bool  simtrack_longLived 	    = false;
         //int   simtrack_matchedHit	    = 0;
         int   simtrack_pdgId		    = -1000;
         int   simtrack_numberOfTrackerHits   = 0;
         int   simtrack_numberOfTrackerLayers = 0;
         float simtrack_mass = 0.;
         int   simtrack_status= -1000; 
         bool  simtrack_isFromDispTop         = false;
         int   simtrack_isFromLLP            = 0;
         
         auto foundTPs = recSimColl.find(itTrack);
         if ( foundTPs != recSimColl.end() ) {
           //if (!foundTPs->val.empty()) {
            isSimMatched = true;
            TrackingParticleRef tpr = foundTPs->val[0].first;
            
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
            
            //$$

            float dV1=0;
            float dV2=0;
            float top1_Diff_X = 0;
            float top1_Diff_Y = 0;
            float top1_Diff_Z = 0;
            float top2_Diff_X = 0;
            float top2_Diff_Y = 0;
            float top2_Diff_Z = 0;


            if ( ntop >= 1 ) {
              float dV1 = (genVertexPos_X - top1_x)*(genVertexPos_X - top1_x)
              + (genVertexPos_Y - top1_y)*(genVertexPos_Y - top1_y)
              + (genVertexPos_Z - top1_z)*(genVertexPos_Z - top1_z);
              top1_Diff_X = genVertexPos_X - top1_x;
              top1_Diff_Y = genVertexPos_Y - top1_y;
              top1_Diff_Z = genVertexPos_Z - top1_z;
              if ( dV1 < 0.01 ) simtrack_isFromLLP = 1;
            }
            if ( ntop >= 2 && simtrack_isFromLLP != 1 ) {
              float dV2 = (genVertexPos_X - top2_x)*(genVertexPos_X - top2_x)
              + (genVertexPos_Y - top2_y)*(genVertexPos_Y - top2_y)
              + (genVertexPos_Z - top2_z)*(genVertexPos_Z - top2_z);
              top2_Diff_X = genVertexPos_X - top2_x;
              top2_Diff_Y = genVertexPos_Y - top2_y;
              top2_Diff_Z = genVertexPos_Z - top2_z;
              if ( dV2 < 0.01 ) simtrack_isFromLLP = ntop;
            }

            tree_simtrack_genVertexPosDiff_top1_X.push_back(top1_Diff_X);
            tree_simtrack_genVertexPosDiff_top1_Y.push_back(top1_Diff_Y);
            tree_simtrack_genVertexPosDiff_top1_Z.push_back(top1_Diff_Z);
            tree_simtrack_genVertexPosDiff_top1_dV.push_back(dV1);//if ( dV1 < 0.01 ) simtrack_isFromLLP = 1;

            tree_simtrack_genVertexPosDiff_top2_X.push_back(top2_Diff_X);
            tree_simtrack_genVertexPosDiff_top2_Y.push_back(top2_Diff_Y);
            tree_simtrack_genVertexPosDiff_top2_Z.push_back(top2_Diff_Z);
            tree_simtrack_genVertexPosDiff_top2_dV.push_back(dV2);//if ( dV2 < 0.01 ) simtrack_isFromLLP = ntop;
            //comparer les infos de simtrack avec gen (diff_X,...)
	   
            //ex: sélectionner les traces venant du premeir top. Regarder à quelles poitns le  vertex se rapprotche du vertex généré
            //aussi pour le deuxieme top
            //Matcher les informations avec les x y z calculés précédemment
	   
            for (auto genIterator=genParticles->begin(); genIterator!=genParticles->end(); genIterator++) //loop on gen particles (Jérémy's method)
            {
              if (abs(genIterator->pdgId())==6 && genIterator->status()==22) //if the gen particle is a displaced top 
              {  //Doublonde la méthode de  Daniel ligne ~1500//
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
	 
          tree_track_simtrack_charge 	            .push_back(simtrack_charge);	     
          tree_track_simtrack_pt		                .push_back(simtrack_pt);  		     
          tree_track_simtrack_eta		              .push_back(simtrack_eta); 	     
          tree_track_simtrack_phi		              .push_back(simtrack_phi); 		     
          tree_track_simtrack_longLived	          .push_back(simtrack_longLived);		      
          tree_track_simtrack_pdgId  	            .push_back(simtrack_pdgId);      
          tree_track_simtrack_isFromLLP            .push_back(simtrack_isFromLLP); 
          tree_track_simtrack_numberOfTrackerHits  .push_back(simtrack_numberOfTrackerHits); 
          tree_track_simtrack_numberOfTrackerLayers.push_back(simtrack_numberOfTrackerLayers); 
          tree_track_simtrack_mass		              .push_back(simtrack_mass        );	
          tree_track_simtrack_status 	            .push_back(simtrack_status      );	
          tree_track_simtrack_isFromDispTop        .push_back(simtrack_isFromDispTop      );	
          
          tree_track_genVertexPos_X  	      .push_back(genVertexPos_X);	     
          tree_track_genVertexPos_Y  	      .push_back(genVertexPos_Y);	     
          tree_track_genVertexPos_Z  	      .push_back(genVertexPos_Z); 
          
          tree_track_nSimHits		            .push_back(nSimHits); 
          tree_track_isSimMatched		        .push_back(isSimMatched);
        }
        nTracks++; 
        
        //if( itTrack->pt() < 1 ) std::cout << " there is a pt problem " << itTrack->pt() << std::endl;
        
        /*if( !isSimMatched &&  trackToVextexMap[iTrack] != 0 ) nUnmatchTrack_fromPU++;
        if( !isSimMatched &&  trackToVextexMap[iTrack] == 0 ) nUnmatchTrack_fromPV++;
        if( !isSimMatched                                   ) nUnmatchTrack++;
        if( trackToVextexMap[iTrack] != 0                   ) nPUTrack++;*/
        
        
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
  //loops on tracking particles 
  //------------------------------
  
  if ( !runOnData_ ) {
     reco::SimToRecoCollection simRecColl = associatorByHits.associateSimToReco(trackRefs, tpCollection);
       float NSim=0;/*!*/
       float NSimtoReco=0;/*!*/
       float ratio=0;/*!*/
     for(const TrackingParticleRef& tp: tpCollection) 
     {
       tree_simtrack_simtrack_charge  		     .push_back( tp->charge());	 //infos sur toutes les traces simulées, 
       //ex efficacité reco, voir trace sim à trace reco, comparée aux pt de otoutes les traces sim /*!*///finito
       tree_simtrack_simtrack_pt		           .push_back( tp->pt());		 //avec un rapport (attention de bien prendre un pt simulé)/*!*/
       tree_simtrack_simtrack_eta		         .push_back( tp->eta());	 
       tree_simtrack_simtrack_phi		         .push_back( tp->phi());	 
       tree_simtrack_simtrack_longLived	     .push_back( tp->longLived());   
       //tree_simtrack_simtrack_matchedHit	      .push_back( tp->matchedHit());	
       tree_simtrack_simtrack_pdgId		       .push_back( tp->pdgId());    
       tree_simtrack_simtrack_numberOfTrackerHits   .push_back( tp->numberOfTrackerHits()); 
       tree_simtrack_simtrack_numberOfTrackerLayers .push_back( tp->numberOfTrackerLayers());
       tree_simtrack_simtrack_mass			             .push_back( tp->mass());
       tree_simtrack_simtrack_status  		           .push_back( tp->status());
       
       tree_simtrack_genVertexPos_X			 .push_back( tp->vx());
       tree_simtrack_genVertexPos_Y			 .push_back( tp->vy());
       tree_simtrack_genVertexPos_Z			 .push_back( tp->vz());
       
       bool isRecoMatched = false;
       NSim = NSim + 1;/*!*/
       std::vector<int> tkIdx;
       auto foundTracks = simRecColl.find(tp);
       if(foundTracks != simRecColl.end()) 
       {
         isRecoMatched = true;
         NSimtoReco=NSimtoReco+1;/*!*/
         tree_simtrack_isRecoMatched_pt.push_back(tp->pt());/*!*/ //CHeck the ratio of histos: ptsimtoreco/ptsim(the one from)
         //simtrack_simtrack is close to 1
         for(const auto trackQuality: foundTracks->val) 
	       {
		      tkIdx.push_back(trackQuality.first.key());
	       }
        }
        
        tree_simtrack_isRecoMatched.push_back(isRecoMatched);
	      //Calcualte the impact parameters w.r.t. PCA
	      tree_simtrack_trkIdx.push_back(tkIdx);
	 
	      TrackingParticle::Vector momentum = parametersDefiner->momentum(iEvent,iSetup,tp);
	      TrackingParticle::Point vertex = parametersDefiner->vertex(iEvent,iSetup,tp);
	      auto dxySim = TrackingParticleIP::dxy(vertex, momentum, bs.position());
	      auto dzSim  = TrackingParticleIP::dz(vertex, momentum, bs.position());
	 
	      tree_simtrack_pca_dxy      .push_back(dxySim);
	      tree_simtrack_pca_dz	    .push_back(dzSim);

       } //end loop on tracking particles 
     ratio= NSimtoReco/NSim;/*!*/
     std::cout<< "ratio of NSimtoReco/Nsim = "<< ratio << std::endl;
     tree_NSimtoReco.push_back(NSimtoReco);
     tree_NSim.push_back(NSim);
     tree_Ratio_NStRNS.push_back(ratio);
   }

   
   //////////////////////////////////
   //////////////////////////////////
   ////////   jets ///  /////////////
   //////////////////////////////////
   //////////////////////////////////
   
   //cout << "number of jets "<<nJet<<endl;  
   //for ak4
   for(int ij=0;ij<int(ak4slimmedJets->size());ij++)
     {
       const Jet& jet = ak4slimmedJets->at(ij); 
       if ( jet.pt() < 20. ) continue;
       tree_AK4Slimmedjet_E.push_back(jet.energy());
       tree_AK4Slimmedjet_pt.push_back(jet.pt());
       tree_AK4Slimmedjet_eta.push_back(jet.eta());
       tree_AK4Slimmedjet_phi.push_back(jet.phi());
       TLorentzVector thejet, track;
       thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
       for(unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) 
	 {
	   track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
	   if(thejet.DeltaR(track) < 0.4 ) tree_AK4Slimmedjet_idxTrack.push_back(iTrack);    
	 }
     }
   
   for(int ij=0;ij<int(ak4PFJets->size());ij++)
     {
       const Jet& jet = ak4PFJets->at(ij); 
       if ( jet.pt() < 20. ) continue;
       tree_AK4PFjet_E.push_back(jet.energy());
       tree_AK4PFjet_pt.push_back(jet.pt());
       tree_AK4PFjet_eta.push_back(jet.eta());
       tree_AK4PFjet_phi.push_back(jet.phi());
       TLorentzVector thejet, track;
       thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
       for(unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) 
       {
         track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
         if(thejet.DeltaR(track) < 0.4 ) tree_AK4PFjet_idxTrack.push_back(iTrack);     
       }
     }
   
   for(int ij=0;ij< int(CaloJets->size());ij++)//plus simpliste, utilisete just eles cellules calorimetriques
     {
       const Jet& jet = CaloJets->at(ij);
       if ( jet.pt() < 20. ) continue; 
       tree_CaloJet_E.push_back(jet.energy());
       tree_CaloJet_pt.push_back(jet.pt());
       tree_CaloJet_eta.push_back(jet.eta());
       tree_CaloJet_phi.push_back(jet.phi());
       TLorentzVector thejet, track;
       thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
       for(unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) {
         track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
         if(thejet.DeltaR(track) < 0.4 ) tree_CaloJet_idxTrack.push_back(iTrack);    
        }
     }
     //for ak8
     for(int ij=0;ij< int(ak8jets->size() );ij++)
     {
       const Jet& jet = ak8jets->at(ij);
       if ( jet.pt() < 20. ) continue;
       tree_jet08_E.push_back(jet.energy());
       tree_jet08_pt.push_back(jet.pt());
       tree_jet08_eta.push_back(jet.eta());
       tree_jet08_phi.push_back(jet.phi());
       TLorentzVector thejet, track;
       thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
       for(unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) 
       {
         track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
         if(thejet.DeltaR(track) < 0.8 ) tree_jet08_idxTrack.push_back(iTrack);    
        }
     }
   
   //    for(int ij=0;ij< int(ak8CaloJets->size()) ;ij++)
   //     {
   //      const Jet& jet = ak8CaloJets->at(ij);
   //    if ( jet.pt() < 20. ) continue;
   //      tree_CaloJet08_E.push_back(jet.energy());
   //      tree_CaloJet08_pt.push_back(jet.pt());
   //      tree_CaloJet08_eta.push_back(jet.eta());
   //      tree_CaloJet08_phi.push_back(jet.phi());
   //      TLorentzVector thejet, track;
   //      thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
   //      for(unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) 
   //       {
   //        track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
   //        if(thejet.DeltaR(track) < 0.8 ) tree_CaloJet08_idxTrack.push_back(iTrack);    
   //       }
   //     }
   
   /////SECONDARY VERTEX RECONSTRUCTION USING DISPLACED TRACKS 
    TransientVertex displacedVertex_top1_general;
    TransientVertex displacedVertex_top2_general;
   if(!runOnData_)
   {
    if(showlog) cout << "  --------- " <<endl;  
    if(showlog) cout << "SECONDARY VERTEX REFITTING	  "<<endl; 
     
    vector<reco::TransientTrack> displacedTracks_top1;
    vector<reco::TransientTrack> displacedTracks_top2;
    
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

    //Jérémy's method 
    
    for (auto genIterator=genParticles->begin(); genIterator!=genParticles->end(); genIterator++) //loop on gen particles
    {
      if (abs(genIterator->pdgId())==6 && genIterator->status()==22) //if the gen particle is a displaced top 
	      { 
	        n_displacedtop++; 
	        int num_track=-1; 
	        if (n_displacedtop==1) 
          {
            top1_x=genIterator->vx(); 
            top1_y=genIterator->vy(); 
            top1_z=genIterator->vz(); 
            if(showlog) cout << "FIRST TOP GEN VERTEX X "<<genIterator->vx()<<endl; 
            if(showlog) cout << "FIRST TOP GEN VERTEX Y "<<genIterator->vy()<<endl; 
            if(showlog) cout << "FIRST TOP GEN VERTEX Z "<<genIterator->vz()<<endl; 
            tree_genTop_X.push_back(top1_x);
            tree_genTop_Y.push_back(top1_y);
            tree_genTop_Z.push_back(top1_z);
            tree_genTop_charge.push_back(genIterator->charge());
	       }
	     
	        if (n_displacedtop==2) {
            top2_x=genIterator->vx(); 
            top2_y=genIterator->vy(); 
            top2_z=genIterator->vz(); 
            if(showlog) cout << "------------------------" <<endl; 
            if(showlog) cout << "SECOND TOP GEN VERTEX X "<<genIterator->vx()<<endl; 
            if(showlog) cout << "SECOND TOP GEN VERTEX Y "<<genIterator->vy()<<endl; 
            if(showlog) cout << "SECOND TOP GEN VERTEX Z "<<genIterator->vz()<<endl; 
            tree_genTop_X.push_back(top2_x);
            tree_genTop_Y.push_back(top2_y);
            tree_genTop_Z.push_back(top2_z);
            tree_genTop_charge.push_back(genIterator->charge());
          }
	     
          //        cout << "------------------------" <<endl; 
          
          //using gen information for vertex refitting
	       for (auto trackIterator=tracks.begin(); trackIterator!=tracks.end(); trackIterator++) //loop on reco tracks
	       {
		      if(trackIterator->pt() < 1 || fabs(trackIterator->eta()) > 2.4) continue;
		      num_track++; 
		      if (tree_track_isSimMatched[num_track]==0) continue; //keeping only the reco tracks matched to sim tracks 
		      if (abs(tree_track_genVertexPos_X[num_track]-genIterator->vx())<0.01 && abs(tree_track_genVertexPos_Y[num_track]-genIterator->vy())<0.01 && abs(tree_track_genVertexPos_Z[num_track]-genIterator->vz())<0.01) //if the position of the sim track associated to the reco track is the same as the gen top vertex
		      {
		        TransientTrack  transientTrack = theTransientTrackBuilder->build(*trackIterator); 
		        if (n_displacedtop==1) displacedTracks_top1.push_back(transientTrack);
		        if (n_displacedtop==2) displacedTracks_top2.push_back(transientTrack);
		      }
	       }
	     
	      }
       }
     
    if(showlog) cout << "------------------------" <<endl;
    if(showlog) cout << "Number of reco tracks for first top "<<displacedTracks_top1.size()<<endl; 
    if(showlog) cout << "Number of reco tracks for second top "<<displacedTracks_top2.size()<<endl; 
     
    if ( displacedTracks_top1.size() > 1 )
       {
	 
	      KalmanVertexFitter theFitter_top1; 
	      //AdaptiveVertexFitter  theFitter_top1;    
	 
	      TransientVertex displacedVertex_top1 = theFitter_top1.vertex(displacedTracks_top1);  // if you want the beam constraint
        displacedVertex_top1_general = displacedVertex_top1;//useless for the moment
	      //TransientVertex myVertex = theFitter.vertex(mytracks);  // if you don't want the beam constraint
	      // now you have a new vertex, can e.g. be compared with the original
	      if (displacedVertex_top1.isValid())
	      {    
          auto lineariser = DefaultLinearizationPointFinder();
	        auto linearization_point = lineariser.getLinearizationPoint(displacedTracks_top1);
	        if(showlog) std::cout << "seed top1 is: " << linearization_point.x() << " - " << linearization_point.y() << " - " << linearization_point.z() << '\n';
		      if(showlog) cout << "------------------------" <<endl;
		      if(showlog) std::cout << " FIRST TOP DISPLACED VERTEX X POS " << displacedVertex_top1.position().x() << std::endl;//reco
		      if(showlog) std::cout << " FIRST TOP DISPLACED VERTEX Y POS " << displacedVertex_top1.position().y() << std::endl;
		      if(showlog) std::cout << " FIRST TOP DISPLACED VERTEX Z POS " << displacedVertex_top1.position().z() << std::endl; 
		      if(showlog) std::cout << " FIRST TOP DISPLACED VERTEX CHI2/ndof  " << displacedVertex_top1.normalisedChiSquared()<< std::endl; 
		      if(showlog) std::cout << " NUMBER OF ORIGINAL TRACKS	" << displacedVertex_top1.originalTracks().size() << std::endl; 
	    
	        secondaryVtx_diff_X = abs(displacedVertex_top1.position().x()-top1_x);//Kalman fitted vertex-original genVertex
	        secondaryVtx_diff_Y = abs(displacedVertex_top1.position().y()-top1_y);
	        secondaryVtx_diff_Z = abs(displacedVertex_top1.position().z()-top1_z);
	        secondaryVtx_nTracks = displacedVertex_top1.originalTracks().size();
	        secondaryVtx_isValid = 1;
	        secondaryVtx_NChi2 = displacedVertex_top1.normalisedChiSquared();
	     
	        tree_secondaryVtx_X.push_back(displacedVertex_top1.position().x()); 
	        tree_secondaryVtx_Y.push_back(displacedVertex_top1.position().y()); 
	        tree_secondaryVtx_Z.push_back(displacedVertex_top1.position().z()); 
	     
	        tree_secondaryVtx_diff_X.push_back(secondaryVtx_diff_X);
	        tree_secondaryVtx_diff_Y.push_back(secondaryVtx_diff_Y);
	        tree_secondaryVtx_diff_Z.push_back(secondaryVtx_diff_Z);
	        tree_secondaryVtx_nTracks.push_back(secondaryVtx_nTracks); 
	        tree_secondaryVtx_isValid.push_back(secondaryVtx_isValid);
	        tree_secondaryVtx_NChi2.push_back(secondaryVtx_NChi2);
	      } 
	      else { 
	        if(showlog) cout << "------------------------" <<endl;
          if(showlog) cout << "Refitted vertex for top 1 is not okay "<<endl;
          secondaryVtx_isValid=0; 
          tree_secondaryVtx_isValid.push_back(secondaryVtx_isValid);
        }
    }
    else{
	    //      cout << "------------------------" <<endl;
	    //      std::cout << "not enough tracks left for top 1" <<std::endl;
	    secondaryVtx_isValid=-1;
	    tree_secondaryVtx_isValid.push_back(secondaryVtx_isValid);
    }
     
    if ( displacedTracks_top2.size() > 1 )
    { 
	      auto lineariser = DefaultLinearizationPointFinder();
	      auto linearization_point = lineariser.getLinearizationPoint(displacedTracks_top2);
	      if(showlog) std::cout << "seed top 2 is: " << linearization_point.x() << " - " << linearization_point.y() << " - " << linearization_point.z() << '\n';
	      KalmanVertexFitter theFitter_top2; 
	      //AdaptiveVertexFitter  theFitter_top2;
	 
        TransientVertex displacedVertex_top2 = theFitter_top2.vertex(displacedTracks_top2);  // if you want the beam constraint
        displacedVertex_top2_general = displacedVertex_top2;
        //TransientVertex myVertex = theFitter.vertex(mytracks);  // if you don't want the beam constraint
        // now you have a new vertex, can e.g. be compared with the original
        if (displacedVertex_top2.isValid())
        {    
           if(showlog)        cout << "------------------------" <<endl;
          if(showlog) std::cout << " SECOND TOP DISPLACED VERTEX X POS " << displacedVertex_top2.position().x() << std::endl;
          if(showlog) std::cout << " SECOND TOP DISPLACED VERTEX Y POS " << displacedVertex_top2.position().y() << std::endl;
          if(showlog) std::cout << " SECOND TOP DISPLACED VERTEX Z POS " << displacedVertex_top2.position().z() << std::endl; 
          if(showlog) std::cout << " SECOND TOP DISPLACED VERTEX CHI2/ndof  " << displacedVertex_top2.normalisedChiSquared()<< std::endl; 
          if(showlog) std::cout << " NUMBER OF ORIGINAL TRACKS        " << displacedVertex_top2.originalTracks().size() << std::endl;
          
	     
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
          tree_secondaryVtx_NChi2.push_back(secondaryVtx_NChi2);//Nombre  entéres étrange.
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
         //	cout << "------------------------" <<endl;
         //	std::cout << "not enough tracks left for top 2" <<std::endl;
    }   
   }
   
  //////////////////
   
  //------------------------------
  //selection of displaced tracks
  //------------------------------
  if(showlog) cout << "//////////////////////////"<< endl; 
  if(showlog) cout << "start to select displaced tracks " << endl;
  vector<reco::TransientTrack> displacedTTracks;
   
  //   int n_track=-1; 
  //   int nn_track=0;   
  ///// MVA for track selection coming from displaced tops 

  //ajouté par Paul /*!*/
  float drSig,dzSig, isinjet,tmpjet;/*!*/
  int jet;/*!*/
  float ntrk10, ntrk1020 ,ntrk20,ntrk2030,ntrk30,ntrk3040,ntrk40,ntrk40XX;/*!*/
  float firsthit_X, firsthit_Y, firsthit_Z, dxy, dxyError, dz, dzError, pt, eta, NChi2, nhits, algo; /*!*/
  float nTracks_b4sel=trackRefs.size();
  float nTracks_ADsel=0;
  float Ratio_b4ADsel=0;
  double bdtval = -10;
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
  /*!*/
  //$$
  reader->BookMVA( "BDTG", weightFile_ );
  //reader->BookMVA("BDT",weightFile_)/*!*/
  //$$
   
  int counter_track=-1; 
  int counter_othertrack=-1;///*!*/
  for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack)//Loop on all the tracks
    {
      counter_track++; 
      const auto& itTrack = trackRefs[iTrack];
      //double track_sig =  abs(tree_track_dxy[counter_track]/tree_track_dxyError[counter_track]); 
      firsthit_X = tree_track_firsthit_X[counter_track]; 
      firsthit_Y = tree_track_firsthit_Y[counter_track]; 
      firsthit_Z = tree_track_firsthit_Z[counter_track]; 
      dxy        = tree_track_dxy[counter_track]; 
      dxyError   = tree_track_dxyError[counter_track]; 
      dz         = tree_track_dz[counter_track]; 
      dzError    = tree_track_dzError[counter_track]; 
      pt         = tree_track_pt[counter_track]; 
      eta        = tree_track_eta[counter_track]; 
      NChi2      = tree_track_NChi2[counter_track]; 
      nhits      = tree_track_nhits[counter_track]; 
      algo       = tree_track_algo[counter_track]; 
      //Ajouté par Paul /*!*/
      if ( dxyError > 0 ) drSig = dxy / dxyError;/*!*/
      if ( dzError > 0 ) dzSig = dz / dzError;/*!*/
      jet = tree_track_recoAK4SlimmedJet_idx[counter_track];/*!*/
      if ( jet >= 0 ) isinjet = 1;/*!*/
      tmpjet=isinjet+2;//bs
      if (tmpjet==3){dzSig=dzSig-1;}//bs
      dzSig=dzSig+1;//bs
      //Computation of the distances needed for the BDT
      for (size_t iTrack2=0; iTrack2<trackRefs.size(); ++iTrack2)    // Loop on all the other Tracks/*!*/
        {
          
          counter_othertrack++;
          if ( counter_othertrack == counter_track ) continue;
          float pt2  = tree_track_pt[counter_othertrack];
          float dr2 = abs(tree_track_dxy[counter_othertrack]);
          float drSig2 = -1.;
          if ( tree_track_dxyError[counter_othertrack] > 0 ) drSig2 = dr2 / tree_track_dxyError[counter_othertrack];
          float chi2 = tree_track_NChi2[counter_othertrack];
          if ( !(pt2 >= 1. && chi2 <= 5. && drSig2 >= 5.) ) continue;//On regarde les autres track_selec[i] qui sont True donc de potnetielles tracks secondaires
          float x2 = tree_track_firsthit_X[counter_othertrack];
          float y2 = tree_track_firsthit_Y[counter_othertrack];
          float z2 = tree_track_firsthit_Z[counter_othertrack];
          float dist = TMath::Sqrt( (firsthit_X-x2)*(firsthit_X-x2) + (firsthit_Y-y2)*(firsthit_Y-y2) + (firsthit_Z-z2)*(firsthit_Z-z2) );//pour chaque reconstruite, on regarde les autres tracks, 
          if ( dist < 10. ) {ntrk10++;}
          if (dist > 10 && dist < 20) {ntrk1020++;}
          if ( dist < 20. ) {ntrk20++;}
          if (dist > 20 && dist < 30) {ntrk2030++;} 
          if ( dist < 30. ) {ntrk30++;} 
          if (dist > 30 && dist < 40) {ntrk3040++;}
          if ( dist < 40. ) {ntrk40++;}//used in BDT
          if (dist > 40) {ntrk40XX++;}
        }  // end Loop on other Tracks

       bdtval = reader->EvaluateMVA( "BDTG" );//default value = -10 (no -10 observed and -999 comes from EvaluateMVA)
       tree_track_MVAval_FromDispTop.push_back(bdtval); 
       //cout << "BDT VAL " << bdtval <<endl; 
       
       if(tree_track_pt[counter_track] < 1 || fabs(tree_track_eta[counter_track]) > 2.4 || tree_track_recoAK4PFJet_idx[counter_track] == -1  || drSig<5 || tree_track_NChi2[counter_track]>5) continue;
       //if (bdtval < 0.0674  ) continue; //optimal cut for 10 cm, trained on 19k events, <s>=22 and <b>=240 //Jérémy
       //  if (bdtval < 0.07381 ) continue; //optimal cut for 50 cm, trained on 10k events, <s>=15 and <b>=234 //Jérémy
       if (bdtval< -0.0401) continue; //optimal cut for 50 cm ,trained on 10k events //Paul, expected to change
       nTracks_ADsel++;
        
       displacedTTracks.push_back(theTransientTrackBuilder->build(*itTrack));
       tree_displacedTTracks_pt.push_back(tree_track_pt[counter_track]);
       //cout << "kept track and it is from displaced top "<<tree_track_simtrack_isFromDispTop[counter_track]<<endl; 
    }  //End loop on all the tracks
    Ratio_b4ADsel = nTracks_ADsel/nTracks_b4sel;
    std::cout << " Ratio of tracks before and after sel with tmva = "<< Ratio_b4ADsel << std::endl;
    tree_nTracks_ADsel.push_back(nTracks_ADsel);
    tree_nTracks_b4sel.push_back(nTracks_b4sel);
    tree_Ratio_b4ADsel.push_back(Ratio_b4ADsel);

  //-----------------------------!!!!!!!!!!!!!!!!---------------------------------//
  //Add Block here to reduce the number of tracks/vertices to treat as it is taking
  //  a lot of computing time or add resctrictions in the following loop of reconstruction of
  // seed vertex from pair of tracks
  // add kalman filters.
  //-----------------------------!!!!!!!!!!!!!!!!---------------------------------//


  //   int n_tracks_test=0;
  //   int n_tracks_selected = 0; 
   
  if(showlog) cout << displacedTTracks.size() << " displaced tracks selected " << endl;
   
   //------------------------------------------------------------
   //reconstruct seed vertex from pair of tracks
   //select those who have a chi2/ndf < 100
   //------------------------------------------------------------
   /*!*/ //get all the track for 1 track (that pass the selection and tmva) for 1 top: comapre to Jérémy's method and check the efficiccy loss
   KalmanVertexFitter theFitter_vertex(kvfPSet);//One or less vertex /*!*/
   
   std::vector< std::pair< std::vector<reco::TransientTrack>, TransientVertex> > proto_vertex; // container of the vertex and associated tracks
      //Paul implémentation:/*!*/
   Proto PVtx(proto_vertex);/*!*/

   for(unsigned int itracks = 0; itracks < displacedTTracks.size(); itracks++){
     
     for(unsigned int itracks2 = itracks+1; itracks2 < displacedTTracks.size(); itracks2++){
       std::vector<reco::TransientTrack> temp_coll_tracks;
       
       double r1 = sqrt(pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().x(),2)+pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().y(),2)); 
       double r2 = sqrt(pow(displacedTTracks[itracks2].innermostMeasurementState().globalPosition().x(),2)+pow(displacedTTracks[itracks2].innermostMeasurementState().globalPosition().y(),2)); 
       double phi1 = TMath::ATan2(displacedTTracks[itracks].innermostMeasurementState().globalPosition().y(),displacedTTracks[itracks].innermostMeasurementState().globalPosition().x());
       double phi2 = TMath::ATan2(displacedTTracks[itracks2].innermostMeasurementState().globalPosition().y(),displacedTTracks[itracks2].innermostMeasurementState().globalPosition().x());
       double dr = abs(r1-r2); 
       double dz = abs(displacedTTracks[itracks].innermostMeasurementState().globalPosition().z()-displacedTTracks[itracks2].innermostMeasurementState().globalPosition().z());
       double dd = sqrt(pow(dr,2)+pow(dz,2)); 
       double dphi = 0.; 
       
       if (abs(phi1 - phi2) < (TMath::Pi()) ) dphi = abs(phi1 - phi2); 
       if (abs(phi1 - phi2) > (TMath::Pi()) ) dphi = 2*(TMath::Pi()) - abs(phi1 - phi2); 
       double dist_to_tracks = pow( 
        pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().x()- displacedTTracks[itracks2].innermostMeasurementState().globalPosition().x(), 2) +
        pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().y()- displacedTTracks[itracks2].innermostMeasurementState().globalPosition().y(), 2) +
        pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().z()- displacedTTracks[itracks2].innermostMeasurementState().globalPosition().z(), 2) 
        , 0.5
      );
      // !=dd

       //if (dd>40 || dphi>1.4) continue; 
       //if(dist_to_tracks > 30) continue;
       
       temp_coll_tracks.push_back(displacedTTracks[itracks]); 
       temp_coll_tracks.push_back(displacedTTracks[itracks2]); 
       TransientVertex displacedVertex_vertex = theFitter_vertex.vertex(temp_coll_tracks); 
       //if(displacedVertex_vertex.isValid() && displacedVertex_vertex.normalisedChiSquared() < 1000){
       if(displacedVertex_vertex.isValid() ){
         std::pair< std::vector<reco::TransientTrack>, TransientVertex> pair_iterators_vertex(temp_coll_tracks,displacedVertex_vertex) ;
         tree_seedVtx_X.push_back(displacedVertex_vertex.position().x());
         tree_seedVtx_Y.push_back(displacedVertex_vertex.position().y());
         tree_seedVtx_Z.push_back(displacedVertex_vertex.position().z());
         tree_seedVtx_dd.push_back(dd); 
         tree_seedVtx_dphi.push_back(dphi); 
         tree_seedVtx_distance2track.push_back(dist_to_tracks); 
         tree_seedVtx_normChi2.push_back(displacedVertex_vertex.normalisedChiSquared()); 
         

         proto_vertex.push_back(pair_iterators_vertex);
         //  PVtx.GetProto().push_back(pair_iterators_vertex);/*!*//*This doesn't  work*/
         PVtx.PushBack(pair_iterators_vertex);/*!*//*it does*/        
       }
     }
   }
   //displacedTTracks & displacedVertex_vertex
  /*!*/ // Diviser l'event en deux : tmva ou axe sprivilégiés dans l'event (voir root: vectuer o uaxe privilégié) pt eta phi?
  /*!*/ // regarder le vecteur des muons (corrélés) à celui des neutralinos
  // => kalman filter pour chaque hemisphere 
  //terminer le rapport là-dessus : efficacité et reconstruction des vertxex des neutralions
  // Daniel is looking at this atm

   if(showlog) cout << "-----------------" <<endl; 
  //  if(showlog) cout << proto_vertex.size() << " proto-vertex selected " << endl;
   if(showlog) cout << PVtx.Size() << " proto-vertex selected " << endl;/*!*/
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
   
   
  //  for(unsigned int  idx_proto_vertex= 0; idx_proto_vertex < proto_vertex.size(); idx_proto_vertex++){
   for(unsigned int  idx_proto_vertex= 0; idx_proto_vertex < PVtx.Size(); idx_proto_vertex++){  /*!*/
    //  loop on track, reco vertex, test and reject or select track
     for(unsigned int itracks = 0; itracks < displacedTTracks.size(); itracks++){
      //  std::vector<reco::TransientTrack> temp_coll_tracks = proto_vertex[idx_proto_vertex].first;
       std::vector<reco::TransientTrack> temp_coll_tracks = PVtx.TTrack(idx_proto_vertex);/*!*/
      //  first hit distance to vertex
      //  double dist_to_vertex = pow( 
      //    pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().x()- proto_vertex[idx_proto_vertex].second.position().x(), 2) +
      //    pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().y()- proto_vertex[idx_proto_vertex].second.position().y(), 2) +
      //    pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().z()- proto_vertex[idx_proto_vertex].second.position().z(), 2) 
      //    , 0.5);
       double dist_to_vertex = pow( 
         pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().x()- PVtx.x(idx_proto_vertex), 2) +
         pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().y()- PVtx.y(idx_proto_vertex), 2) +
         pow(displacedTTracks[itracks].innermostMeasurementState().globalPosition().z()- PVtx.z(idx_proto_vertex), 2) 
         , 0.5);
      //  cout << "distance check to vertex " << dist_to_vertex << endl;
      // Selection criteria tochange 
       if( //dist_to_vertex < 30 && 
          fabs(displacedTTracks[itracks].innermostMeasurementState().globalPosition().x() - temp_coll_tracks[0].innermostMeasurementState().globalPosition().x()  ) > 0.001 &&  
          fabs(displacedTTracks[itracks].innermostMeasurementState().globalPosition().x() - temp_coll_tracks[1].innermostMeasurementState().globalPosition().x()  ) > 0.001 ){
            temp_coll_tracks.push_back(displacedTTracks[itracks]);
            TransientVertex displacedVertex_vertex = theFitter_vertex.vertex(temp_coll_tracks); 
            // /chi2 < 100 in original 
            if(displacedVertex_vertex.isValid() && displacedVertex_vertex.normalisedChiSquared() < 50){
              // proto_vertex[idx_proto_vertex].first = temp_coll_tracks;
              // proto_vertex[idx_proto_vertex].second = displacedVertex_vertex;
              PVtx.TTrack(idx_proto_vertex)=temp_coll_tracks;/*!*/
              PVtx.TVertex(idx_proto_vertex)= displacedVertex_vertex;/*!*/
            } else temp_coll_tracks.pop_back();
        }
     }
   }
 

   if(showlog)  cout << "-----------------"<<endl; 
  //  if(showlog) cout << " number of vertex after adding tracks " << proto_vertex.size() << endl;
   if(showlog) cout << " number of vertex after adding tracks " << PVtx.Size() << endl;/*!*/
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
   /****
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
   
   //$$
   if( HT_val > 180. && MuonSup10GeV && MuonSup28GeV && invMassGood ) tree_passesHTFilter.push_back(true); 
   else tree_passesHTFilter.push_back(false);
   //$$
   *******/

   //BDT to make tracks clustering to get the best tracks for each vertex==> change TrackingPerf to get two MVA files 
  //  in input//

  //  //ajouté par Paul /*!*/
  // // --------- BDT Variables-------------//
  //  float drSig,dzSig, isinjet,tmpjet;/*!*/
  //  int jet;/*!*/
  //  float ntrk10, ntrk1020 ,ntrk20,ntrk2030,ntrk30,ntrk3040,ntrk40,ntrk40XX;/*!*/
  //  float firsthit_X, firsthit_Y, firsthit_Z, dxy, dxyError, dz, dzError, pt, eta, NChi2, nhits, algo; /*!*/
  // //-------------- New variables compared to last BDT ------------------------//
// float dphi, dR, dd,ddError,ddSig, dr,drError,drSig,dz,dzError,dzSig ;

  // //---------------------------------------------------------------------------// //
  //  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
   
  // // reader->AddVariable( "mva_track_firstHit_x", &firsthit_X );//to be exluded
  // // reader->AddVariable( "mva_track_firstHit_y", &firsthit_Y );
  // // reader->AddVariable( "mva_track_firstHit_z", &firsthit_Z );
  // // reader->AddVariable( "mva_track_firstHit_dxy", &dxy );
  // // reader->AddVariable( "mva_track_firstHit_dxyError", &dxyError );
  // // reader->AddVariable( "mva_track_firstHit_dz", &dz );
  // // reader->AddVariable( "mva_track_firstHit_dzError", &dzError );
  //  reader->AddVariable( "mva_track_pt", &pt );
  //  reader->AddVariable( "mva_track_eta", &eta );
  //  reader->AddVariable( "mva_track_nchi2", &NChi2 );
  //  reader->AddVariable( "mva_track_nhits", &nhits );
  //  reader->AddVariable( "mva_track_algo", &algo); 
  //    //add the variables from my BDT(Paul)
  // reader->AddVariable( "mva_ntrk10", &ntrk10);
  // reader->AddVariable( "mva_drSig", &drSig); /*!*/
  // reader->AddVariable( "mva_track_isinjet", &isinjet); /*!*/
  // /*!*/
  // //-------------- New variables compared to last BDT ------------------------//
  // reader->AddVariable("mva_track_dphi",&dphi);
  // reader->AddVariable("mva_track_dR",&dR);


  // //---------------------------------------------------------------------------// //
  //  //$$
  //  reader->BookMVA( "BDTG", weightFile_ );
  //  //$$
   
  //  int counter_track=-1; 
  //  int counter_othertrack=-1;//(Paul)
  //  for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack)
  //    {
  //      counter_track++; 
  //      const auto& itTrack = trackRefs[iTrack];
  //     //    double track_sig =  abs(tree_track_dxy[counter_track]/tree_track_dxyError[counter_track]);
       
  //      firsthit_X = tree_track_firsthit_X[counter_track]; 
  //      firsthit_Y = tree_track_firsthit_Y[counter_track]; 
  //      firsthit_Z = tree_track_firsthit_Z[counter_track]; 
  //      dxy        = tree_track_dxy[counter_track]; 
  //      dxyError   = tree_track_dxyError[counter_track]; 
  //      dz         = tree_track_dz[counter_track]; 
  //      dzError    = tree_track_dzError[counter_track]; 
  //      pt         = tree_track_pt[counter_track]; 
  //      eta        = tree_track_eta[counter_track]; 
  //      NChi2      = tree_track_NChi2[counter_track]; 
  //      nhits      = tree_track_nhits[counter_track]; 
  //      algo       = tree_track_algo[counter_track]; 
  //      //Ajouté par Paul
  //      if ( dxyError > 0 ) drSig = dxy / dxyError;
  //      if ( dzError > 0 ) dzSig = dz / dzError;
  //      jet = tree_track_recoAK4SlimmedJet_idx[counter_track];
  //      if ( jet >= 0 ) isinjet = 1;
  //      tmpjet=isinjet+2;
  //      if (tmpjet==3){dzSig=dzSig-1;}
  //     dzSig=dzSig+1;
  //   //Distances
  //   for (size_t iTrack2=0; iTrack2<trackRefs.size(); ++iTrack2)    // Loop on other Tracks
  //    {
  //      counter_othertrack++;
  //    if ( counter_othertrack == counter_track ) continue;
  //      float pt2  = tree_track_pt[counter_othertrack];
  //      float dr2 = abs(tree_track_dxy[counter_othertrack]);
  //      float drSig2 = -1.;
  //      if ( tree_track_dxyError[counter_othertrack] > 0 ) drSig2 = dr2 / tree_track_dxyError[counter_othertrack];
  //      float chi2 = tree_track_NChi2[counter_othertrack];
  //    if ( !(pt2 >= 1. && chi2 <= 5. && drSig2 >= 5.) ) continue;//On regarde les autres track_selec[i] qui sont True donc de potnetielles tracks secondaires
  //      float x2 = tree_track_firsthit_X[counter_othertrack];
  //      float y2 = tree_track_firsthit_Y[counter_othertrack];
  //      float z2 = tree_track_firsthit_Z[counter_othertrack];
  //      float dist = TMath::Sqrt( (firsthit_X-x2)*(firsthit_X-x2) + (firsthit_Y-y2)*(firsthit_Y-y2) + (firsthit_Z-z2)*(firsthit_Z-z2) );//pour chaque reconstruite, on regarde les autres tracks, 
  //      if ( dist < 10. ) {ntrk10++;}
  //      if (dist > 10 && dist < 20) {ntrk1020++;}
  //      if ( dist < 20. ) {ntrk20++;}
  //      if (dist > 20 && dist < 30) {ntrk2030++;} 
  //      if ( dist < 30. ) {ntrk30++;} 
  //      if (dist > 30 && dist < 40) {ntrk3040++;}
  //      if ( dist < 40. ) {ntrk40++;}//used in BDT
  //      if (dist > 40) {ntrk40XX++;}
  //    }  // end Loop on other Tracks
  //      double bdtval = reader->EvaluateMVA( "BDTG" );
  //      tree_track_MVAval_FromDispTop.push_back(bdtval); 
  //      //cout << "BDT VAL " << bdtval <<endl; 
       
  //      if(tree_track_pt[counter_track] < 1 || fabs(tree_track_eta[counter_track]) > 2.4 || tree_track_recoAK4PFJet_idx[counter_track] == -1  || drSig<5 || tree_track_NChi2[counter_track]>5) continue;
  //      //if (bdtval < 0.0674  ) continue; //optimal cut for 10 cm, trained on 19k events, <s>=22 and <b>=240
  //     //  if (bdtval < 0.07381 ) continue; //optimal cut for 50 cm, trained on 10k events, <s>=15 and <b>=234 
  //     if (bdtval< -0.0401) continue; //optimal cut for50 cm ,(Paul BDTG)
  //      displacedTTracks.push_back(theTransientTrackBuilder->build(*itTrack));
  //      //cout << "kept track and it is from displaced top "<<tree_track_simtrack_isFromDispTop[counter_track]<<endl; 
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
   tree_track_nSimHits.clear(); 
   tree_track_isSimMatched.clear();
//$$   tree_track_nPixel.clear();
//$$   tree_track_nStrip.clear();
   tree_track_vx.clear(); 
   tree_track_vy.clear(); 
   tree_track_vz.clear(); 
   tree_track_firsthit_X.clear(); 
   tree_track_firsthit_Y.clear(); 
   tree_track_firsthit_Z.clear(); 
   tree_track_firsthit_phi.clear(); 
   
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
   tree_track_simtrack_isFromLLP.clear();	 
   tree_track_simtrack_isFromDispTop.clear();

   tree_track_genVertexPos_X.clear();		 
   tree_track_genVertexPos_Y.clear();		 
   tree_track_genVertexPos_Z.clear();

   tree_track_recoVertex_idx.clear();
   tree_track_recoAK4SlimmedJet_idx.clear();
   tree_track_recoAK4PFJet_idx.clear();
   tree_track_reco08Jet_idx.clear(); 
   tree_track_recoCaloJet_idx.clear(); 
//    tree_track_reco08CaloJet_idx.clear(); 

   tree_track_MVAval_FromDispTop.clear(); 

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
   tree_simtrack_genVertexPos_Z .clear();

   tree_simtrack_isRecoMatched.clear();
   tree_simtrack_pca_dxy.clear();
   tree_simtrack_pca_dz .clear();
   tree_simtrack_trkIdx  .clear();
   

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
   
   tree_seedVtx_X.clear(); 
   tree_seedVtx_Y.clear(); 
   tree_seedVtx_Z.clear(); 
   tree_seedVtx_dd.clear(); 
   tree_seedVtx_dphi.clear(); 
   tree_seedVtx_distance2track.clear(); 
   tree_seedVtx_normChi2.clear(); 

   tree_secondaryVtx_iterative_X.clear(); 
   tree_secondaryVtx_iterative_Y.clear(); 
   tree_secondaryVtx_iterative_Z.clear(); 
   tree_secondaryVtx_iterative_nTracks.clear();
   tree_secondaryVtx_iterative_NChi2.clear(); 
   tree_secondaryVtx_iterative_match_top.clear(); 
   tree_secondaryVtx_iterative_isSelected.clear(); 

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

   tree_PFMet_et.clear(); 
   tree_PFMet_phi.clear(); 
   tree_PFMet_sig.clear(); 

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

  
  tree_top1_x.clear();
  tree_top1_y.clear();
  tree_top1_z.clear();    
  tree_top2_x.clear();
  tree_top2_y.clear();
  tree_top2_z.clear();

  tree_top1_pt.clear();
  tree_top1_eta.clear();
  tree_top1_phi.clear();

  tree_top2_pt.clear();
  tree_top2_eta.clear();
  tree_top2_phi.clear();
  tree_ntop.clear();

  tree_simtrack_genVertexPosDiff_top1_X.clear();
  tree_simtrack_genVertexPosDiff_top1_Y.clear();
  tree_simtrack_genVertexPosDiff_top1_Z.clear();
  tree_simtrack_genVertexPosDiff_top1_dV.clear();
  tree_simtrack_genVertexPosDiff_top2_X.clear();
  tree_simtrack_genVertexPosDiff_top2_Y.clear();
  tree_simtrack_genVertexPosDiff_top2_Z.clear();
  tree_simtrack_genVertexPosDiff_top2_dV.clear();

  tree_NSimtoReco.clear();
  tree_NSim.clear();
  tree_Ratio_NStRNS.clear();
  tree_nTracks_ADsel.clear();
  tree_nTracks_b4sel.clear();
  tree_Ratio_b4ADsel.clear();
  tree_displacedTTracks_pt.clear();
  tree_simtrack_isRecoMatched_pt.clear();


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

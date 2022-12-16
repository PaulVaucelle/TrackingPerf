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

//---------------------------------!!!!-----------------------------//
                      //-----------Transient Track/Vtx--------//
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
              //-------------Propagators------------//
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/GeomPropagators/interface/SmartPropagator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include  "TrackingTools/GeomPropagators/interface/StraightLineCylinderCrossing.h"
#include  "TrackingTools/GeomPropagators/interface/StraightLineBarrelCylinderCrossing.h"
#include "TrackingTools/GeomPropagators/interface/StraightLinePlaneCrossing.h"
              //-------------Surfaces---------------//
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/GeometryVector/interface/Point3DBase.h"
              //----------------?-----------------//
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
              //----------------BField--------------//
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"

//------------------------------!!!!------------------------//

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
    const edm::EDGetTokenT<GenEventInfoProduct>         genEventInfoToken_;
    const edm::EDGetTokenT<LHEEventProduct>             LHEEventProductToken_;
    
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
       
    
    /////// VARIABLES FOR FILLING THE TREE
    
    int  runNumber, eventNumber, lumiBlock;
    int  tree_NbrOfZCand;
    bool tree_passesHTFilter;
    int  tree_nTracks = 0, tree_nFromC = 0, tree_nFromB = 0; 
    
    //-----------------------
    //trigger variable
    std::vector<string > tree_trigger_names;
    std::vector<bool >   tree_trigger_bits;
    
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
    int   tree_nPV;
    float tree_PV_x;
    float tree_PV_y;
    float tree_PV_z;
    float tree_PV_ez;
    float tree_PV_NChi2;
    float tree_PV_ndf;
    
    std::vector<float> tree_vtx_PosX;
    std::vector<float> tree_vtx_PosY;
    std::vector<float> tree_vtx_PosZ;
    std::vector<float> tree_vtx_NChi2;
    std::vector<float> tree_vtx_PosXError;
    std::vector<float> tree_vtx_PosYError;
    std::vector<float> tree_vtx_PosZError;
    std::vector< int > tree_vtx_ndf;
    
    //--------------------------------
    // met infos -------
    //--------------------------------
    float tree_PFMet_et;
    float tree_PFMet_phi;
    float tree_PFMet_sig;
    
    //--------------------------------
    // jet infos -------
    //--------------------------------
    
    int tree_njet;
    std::vector<float> tree_jet_E;
    std::vector<float> tree_jet_pt;
    std::vector<float> tree_jet_eta;
    std::vector<float> tree_jet_phi;
//$$    std::vector<int>   tree_jet_idxTrack;
    
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
    // electrons infos -------
    //--------------------------------
    std::vector<float> tree_electron_pt;
    std::vector<float> tree_electron_eta;
    std::vector<float> tree_electron_phi;
    std::vector<float> tree_electron_x;
    std::vector<float> tree_electron_y;
    std::vector<float> tree_electron_z;
    std::vector<float> tree_electron_energy;
    std::vector< int > tree_electron_charge;
    
    //--------------------------------
    // muons infos -------
    //--------------------------------
    // slimmed muons
    std::vector<float> tree_muon_pt;
    std::vector<float> tree_muon_eta;
    std::vector<float> tree_muon_phi;
    std::vector<float> tree_muon_x;
    std::vector<float> tree_muon_y;
    std::vector<float> tree_muon_z;
    std::vector<float> tree_muon_energy;
    std::vector<float> tree_muon_dxy;
    std::vector<float> tree_muon_dxyError;
    std::vector<float> tree_muon_dz;
    std::vector<float> tree_muon_dzError;
    std::vector< int > tree_muon_charge;
    std::vector<bool>  tree_muon_isLoose;
    std::vector<bool>  tree_muon_isTight;
    std::vector<bool>  tree_muon_isGlobal;
    std::vector<bool>  tree_muon_PFisoVeryTight;
    std::vector<bool>  tree_muon_PFisoTight;
    std::vector<bool>  tree_muon_PFisoMedium;
    std::vector<bool>  tree_muon_PFisoLoose;
    std::vector<bool>  tree_muon_MVAisoLoose;
    std::vector<bool>  tree_muon_MVAisoMedium;
    std::vector<bool>  tree_muon_MVAisoTight;
    std::vector<bool>  tree_muon_isStandAloneMuon;
    std::vector<bool>  tree_muon_CutBasedIdMedium;
    std::vector<bool>  tree_muon_CutBasedIdMediumPrompt;
    
    //-----------------------
    //fill the tree per track
    std::vector< float > tree_track_pz;
    std::vector< float > tree_track_pt;
    std::vector< float > tree_track_eta;
    std::vector< float > tree_track_phi;
    std::vector< int >   tree_track_charge;
    std::vector<float >  tree_track_NChi2;
    std::vector<bool >   tree_track_isHighPurity;
    std::vector< float>  tree_track_dxy;
    std::vector< float>  tree_track_dxyError;
    std::vector< float>  tree_track_dz;
    std::vector< float>  tree_track_dzError ;
    std::vector<int>     tree_track_nHit;
    std::vector<int>     tree_track_nHitPixel;
    std::vector<int>     tree_track_nHitTIB;
    std::vector<int>     tree_track_nHitTID;
    std::vector<int>     tree_track_nHitTOB;
    std::vector<int>     tree_track_nHitTEC;
    std::vector<int>     tree_track_nHitPXB;
    std::vector<int>     tree_track_nHitPXF;
    std::vector<int>     tree_track_isHitPixel;
    std::vector<int>     tree_track_nLayers;
    std::vector<int>     tree_track_nLayersPixel;
    std::vector< float > tree_track_x;
    std::vector< float > tree_track_y;
    std::vector< float > tree_track_z;
    std::vector< int >   tree_track_firstHit;
    std::vector<float>   tree_track_firstHit_x;
    std::vector<float>   tree_track_firstHit_y;
    std::vector<float>   tree_track_firstHit_z;
    std::vector<float>   tree_track_firstHit_phi;
    //
    std::vector<float>   tree_track_ddxy;
    std::vector<float>   tree_track_ddxyError;
    std::vector<float>                    tree_track_dpt;
     std::vector<float>     tree_track_deta;
    std::vector<float>                tree_track_dphi;
     std::vector<float>     tree_track_dNChi2;
     std::vector<float>               tree_track_dvx;
     std::vector<float>     tree_track_dvy;
      std::vector<float>              tree_track_dvz;
      std::vector<float>    tree_track_ddz;
       std::vector<float>             tree_track_ddzError;
//$$$$
    std::vector<float>   tree_track_extraTrue_dx;
    std::vector<float>   tree_track_extraTrue_dy;
    std::vector<float>   tree_track_extraTrue_dz;
    std::vector<float>   tree_track_extraTrue_dphi;
    std::vector<float>   tree_track_fromLayer_dx;
    std::vector<float>   tree_track_fromLayer_dy;
    std::vector<float>   tree_track_fromLayer_dz;
    std::vector<float>   tree_track_fromLayer_dphi;
//$$$$
    std::vector<int>     tree_track_iJet;
    std::vector<float>   tree_track_ntrk10;
    std::vector<float>   tree_track_ntrk20;
    std::vector<float>   tree_track_ntrk30;
    std::vector<double>  tree_track_MVAval;
    std::vector<int>     tree_track_Hemi;
    std::vector<double>  tree_track_Hemi_dR;
    std::vector<double>  tree_track_Hemi_mva_NChi2;
    std::vector<int>     tree_track_Hemi_LLP;
    
    std::vector< float > tree_track_ptError;
    std::vector< float > tree_track_outerPt;
    std::vector<bool >   tree_track_isLoose;
    std::vector<bool >   tree_track_isTight;
    std::vector<int>     tree_track_numberOfLostHits;
    std::vector<unsigned int>    tree_track_originalAlgo; // definition as comments at the end of the file,
    //http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_1_3/doc/html/d8/df2/classreco_1_1TrackBase.html#aca7611bd1a33d535cefc72b6e497ece8
    std::vector<unsigned int>    tree_track_algo;
    std::vector<unsigned short>  tree_track_stopReason;
    std::vector<int>     tree_track_numberOfValidStripHits;
    std::vector<int>     tree_track_stripTECLayersWithMeasurement ;
    std::vector<int>     tree_track_stripTIBLayersWithMeasurement;
    std::vector<int>     tree_track_stripTIDLayersWithMeasurement;
    std::vector<int>     tree_track_stripTOBLayersWithMeasurement;
    std::vector<int>     tree_track_recoVertex_idx;
    std::vector<int>     tree_track_recoAK4PFJet_idx;
    std::vector<int>     tree_track_reco08Jet_idx;
    std::vector<int>     tree_track_recoCaloJet_idx;
    std::vector<float>   tree_track_drSig;
//!!!!

    std::vector<float>   tree_track_GeoBarrel_firsthit_x;
    std::vector<float>   tree_track_GeoBarrel_firsthit_y;
    std::vector<float>   tree_track_GeoBarrel_firsthit_z;
    std::vector<float>   tree_track_RECOvsMINI_GeoBarrel_firsthit_x;
    std::vector<float>   tree_track_RECOvsMINI_GeoBarrel_firsthit_y;
    std::vector<float>   tree_track_RECOvsMINI_GeoBarrel_firsthit_z;

    std::vector<float>                         tree_track_PropBarrel_firsthit_x;
    std::vector<float>                         tree_track_PropBarrel_firsthit_y;
    std::vector<float>                         tree_track_PropBarrel_firsthit_z;
    std::vector<float>                         tree_track_RECOvsMINI_PropBarrel_firsthit_x;
    std::vector<float>                         tree_track_RECOvsMINI_PropBarrel_firsthit_y;
    std::vector<float>                         tree_track_RECOvsMINI_PropBarrel_firsthit_z;

    std::vector<float>                         tree_track_PropBarrel_firsthit_opti_x;
    std::vector<float>                         tree_track_PropBarrel_firsthit_opti_y;
    std::vector<float>                         tree_track_PropBarrel_firsthit_opti_z;
    std::vector<float>                         tree_track_RECOvsMINI_PropBarrel_firsthit_opti_x;
    std::vector<float>                         tree_track_RECOvsMINI_PropBarrel_firsthit_opti_y;
    std::vector<float>                         tree_track_RECOvsMINI_PropBarrel_firsthit_opti_z;

    std::vector<float>                         tree_track_GeoDisk_firsthit_x;
    std::vector<float>                         tree_track_GeoDisk_firsthit_y;
    std::vector<float>                         tree_track_GeoDisk_firsthit_z;
    std::vector<float>                         tree_track_RECOvsMINI_GeoDisk_firsthit_x;
    std::vector<float>                         tree_track_RECOvsMINI_GeoDisk_firsthit_y;
    std::vector<float>                         tree_track_RECOvsMINI_GeoDisk_firsthit_z;
                            
    std::vector<float>                         tree_track_GeoDisk_firsthit_opti_x;
    std::vector<float>                         tree_track_GeoDisk_firsthit_opti_y;
    std::vector<float>                         tree_track_GeoDisk_firsthit_opti_z;
    std::vector<float>                         tree_track_RECOvsMINI_GeoDisk_firsthit_opti_x;
    std::vector<float>                         tree_track_RECOvsMINI_GeoDisk_firsthit_opti_y;
    std::vector<float>                         tree_track_RECOvsMINI_GeoDisk_firsthit_opti_z;

    std::vector<float>                         tree_track_PropDisk_SLCC_firsthit_x;
    std::vector<float>                         tree_track_PropDisk_SLCC_firsthit_y;
    std::vector<float>                         tree_track_PropDisk_SLCC_firsthit_z;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_x;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_y;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_z;
                                            
    std::vector<float>                         tree_track_PropDisk_SLCC_firsthit_opti_x;
    std::vector<float>                         tree_track_PropDisk_SLCC_firsthit_opti_y;
    std::vector<float>                         tree_track_PropDisk_SLCC_firsthit_opti_z;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_x;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_y;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_z;

                                
    std::vector<float>                         tree_track_PropDisk_SLPC_firsthit_x;
    std::vector<float>                         tree_track_PropDisk_SLPC_firsthit_y;
    std::vector<float>                         tree_track_PropDisk_SLPC_firsthit_z;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_x;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_y;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_z;

    std::vector<float>                         tree_track_PropDisk_SLPC_firsthit_opti_x;
    std::vector<float>                         tree_track_PropDisk_SLPC_firsthit_opti_y;
    std::vector<float>                         tree_track_PropDisk_SLPC_firsthit_opti_z;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_x;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_y;
    std::vector<float>                         tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_z;

    std::vector<float>   tree_track_theta;
    std::vector<int>     tree_track_surf;
    std::vector<int>     tree_track_hitpattern;

//!!!!

    //std::vector<int> tree_track_nclusters;
    std::vector<int>     tree_track_sim_LLP;
    std::vector<bool>    tree_track_sim_isFromB;
    std::vector<bool>    tree_track_sim_isFromC;
    std::vector< float > tree_track_sim_pt;
    std::vector< float > tree_track_sim_eta  ;
    std::vector< float > tree_track_sim_phi  ;
    std::vector< int >   tree_track_sim_charge;
    std::vector<int>     tree_track_sim_pdgId;
    std::vector<float>   tree_track_sim_mass  ;
    std::vector<float>   tree_track_sim_x;
    std::vector<float>   tree_track_sim_y;
    std::vector<float>   tree_track_sim_z;

    std::vector<int>     tree_track_nSimHits;
    std::vector<bool>    tree_track_isSimMatched;
    std::vector<bool>    tree_track_sim_longLived      ;
    // std::vector<int>  tree_track_sim_matchedHit    ;
    std::vector<int>     tree_track_sim_numberOfTrackerHits  ;
    std::vector<int>     tree_track_sim_numberOfTrackerLayers;
    std::vector<int>     tree_track_sim_status;
    std::vector<int>     tree_track_sim_isFromBC_mother_pdgId;
    std::vector<float>   tree_track_sim_dV;
    std::vector<float>   tree_track_sim_dist;
    
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
    std::vector< float > tree_genParticle_x;
    std::vector< float > tree_genParticle_y;
    std::vector< float > tree_genParticle_z;
    std::vector< float > tree_genParticle_mass;
    std::vector< int >   tree_genParticle_statusCode;
    std::vector< int >   tree_genParticle_mother_pdgId;
    std::vector< int >   tree_genParticle_LLP;

    std::vector< float > tree_genFromC_pt;
    std::vector< float > tree_genFromC_eta;
    std::vector< float > tree_genFromC_phi;
    std::vector< float > tree_genFromC_charge;
    std::vector< int >   tree_genFromC_pdgId;
    std::vector< float > tree_genFromC_x;
    std::vector< float > tree_genFromC_y;
    std::vector< float > tree_genFromC_z;
    std::vector< int >   tree_genFromC_mother_pdgId;
    std::vector< int >   tree_genFromC_generation;
    std::vector< int >   tree_genFromC_LLP;

    std::vector< float > tree_genFromB_pt;
    std::vector< float > tree_genFromB_eta;
    std::vector< float > tree_genFromB_phi;
    std::vector< float > tree_genFromB_charge;
    std::vector< int >   tree_genFromB_pdgId;
    std::vector< float > tree_genFromB_x;
    std::vector< float > tree_genFromB_y;
    std::vector< float > tree_genFromB_z;
    std::vector< int >   tree_genFromB_mother_pdgId;
    std::vector< int >   tree_genFromB_generation;
    std::vector< int >   tree_genFromB_LLP;
   
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
    
    int nEvent ;
    
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

trackToken_(       consumes<edm::View<reco::Track> >(           iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
trackSrc_(         consumes<edm::View<reco::Track> >(           iConfig.getParameter<edm::InputTag>("trackLabel") )),
trackAssociatorToken_( consumes<reco::TrackToTrackingParticleAssociator>(iConfig.getUntrackedParameter<edm::InputTag>("trackAssociator"))),
beamSpotToken_(     consumes<reco::BeamSpot>(               iConfig.getUntrackedParameter<edm::InputTag>("beamSpot"))),
vertexToken_(      consumes<reco::VertexCollection>(           iConfig.getUntrackedParameter<edm::InputTag>("vertices"))),
ak4slimmedJetToken_(   consumes<edm::View<reco::Jet> >(           iConfig.getParameter<edm::InputTag>("ak4slimmedJetInput"))),
ak4PFJetToken_(     consumes<edm::View<reco::Jet> >(           iConfig.getParameter<edm::InputTag>("ak4PFJetInput"))),
CaloJetToken_(     consumes<edm::View<reco::Jet> >(           iConfig.getParameter<edm::InputTag>("caloJetInput"))),
ak8jetToken_(      consumes<edm::View<reco::Jet> >(           iConfig.getParameter<edm::InputTag>("ak8jetInput"))),
//   ak8CaloJetToken_(     consumes<edm::View<reco::Jet> >(           iConfig.getParameter<edm::InputTag>("ak8CaloJetInput"))),
genParticlesToken_(    consumes<reco::GenParticleCollection>(            iConfig.getParameter<edm::InputTag>("genParticles"))),
genJetToken_(          consumes<reco::GenJetCollection>(                 iConfig.getParameter<edm::InputTag>("genJetInput"))),
ak8GenJetToken_(          consumes<reco::GenJetCollection>(              iConfig.getParameter<edm::InputTag>("ak8GenJetInput"))),
genEventInfoToken_(    consumes<GenEventInfoProduct>(                    iConfig.getParameter<edm::InputTag>("genEventInfoInput"))),
LHEEventProductToken_( consumes<LHEEventProduct>(                        iConfig.getParameter<edm::InputTag>("LHEEventProductInput"))),
pfcandsToken_(         consumes<pat::PackedCandidateCollection>(         iConfig.getParameter<edm::InputTag>("pfcands"))),
parametersDefinerName_(                           iConfig.getUntrackedParameter<std::string>("parametersDefiner")),
electronPATToken_(     consumes<pat::ElectronCollection>(                iConfig.getParameter<edm::InputTag>("electronInput"))),
slimmedmuonToken_(     consumes<pat::MuonCollection>(               iConfig.getParameter<edm::InputTag>("slimmedmuonInput"))),
//  recomuonToken_(        consumes<reco::MuonCollection>(               iConfig.getParameter<edm::InputTag>("recomuonInput"))),
metToken_(             consumes<pat::METCollection>(               iConfig.getParameter<edm::InputTag>("metInput"))),
filterTriggerNames_(                                                     iConfig.getUntrackedParameter<std::vector<std::string> >("filterTriggerNames")),
triggerBits_(          consumes<edm::TriggerResults>(                    edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("HLT"))) ),
ZmumuCandidate_(       consumes<vector<reco::CompositeCandidate>>(       iConfig.getParameter<edm::InputTag>("Zmumucand"))),
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
    smalltree->Branch("tree_nPV",      &tree_nPV);
    smalltree->Branch("tree_PV_x",     &tree_PV_x);
    smalltree->Branch("tree_PV_y",     &tree_PV_y);
    smalltree->Branch("tree_PV_z",     &tree_PV_z);    
    smalltree->Branch("tree_PV_ez",    &tree_PV_ez);    
    smalltree->Branch("tree_PV_NChi2", &tree_PV_NChi2);    
    smalltree->Branch("tree_PV_ndf",   &tree_PV_ndf);    
    
    smalltree->Branch("tree_vtx_PosX",      &tree_vtx_PosX);
    smalltree->Branch("tree_vtx_PosY",      &tree_vtx_PosY);
    smalltree->Branch("tree_vtx_PosZ",      &tree_vtx_PosZ);
    smalltree->Branch("tree_vtx_NChi2",     &tree_vtx_NChi2);
    smalltree->Branch("tree_vtx_PosXError", &tree_vtx_PosXError);
    smalltree->Branch("tree_vtx_PosYError", &tree_vtx_PosYError);
    smalltree->Branch("tree_vtx_PosZError", &tree_vtx_PosZError);
    smalltree->Branch("tree_vtx_ndf",       &tree_vtx_ndf);
    
    // trigger info
    smalltree->Branch("tree_trigger_names", &tree_trigger_names);
    smalltree->Branch("tree_trigger_bits",  &tree_trigger_bits);
    
    smalltree->Branch("tree_NbrOfZCand",     &tree_NbrOfZCand,  "tree_NbrOfZCand/I");
    smalltree->Branch("tree_passesHTFilter", &tree_passesHTFilter);
    
    // met info
    smalltree->Branch("tree_PFMet_et" ,  &tree_PFMet_et);
    smalltree->Branch("tree_PFMet_phi" , &tree_PFMet_phi);
    smalltree->Branch("tree_PFMet_sig" , &tree_PFMet_sig);
    
    smalltree->Branch("tree_njet"  ,        &tree_njet);
    smalltree->Branch("tree_jet_E"  ,       &tree_jet_E);
    smalltree->Branch("tree_jet_pt"  ,      &tree_jet_pt);
    smalltree->Branch("tree_jet_eta" ,      &tree_jet_eta);
    smalltree->Branch("tree_jet_phi" ,      &tree_jet_phi);
//$$    smalltree->Branch("tree_jet_idxTrack" , &tree_jet_idxTrack);
    
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
    
    // electrons info
    smalltree->Branch("tree_electron_pt"  ,   &tree_electron_pt);
    smalltree->Branch("tree_electron_eta" ,   &tree_electron_eta);
    smalltree->Branch("tree_electron_phi" ,   &tree_electron_phi);
    smalltree->Branch("tree_electron_x"  ,    &tree_electron_x);
    smalltree->Branch("tree_electron_y" ,     &tree_electron_y);
    smalltree->Branch("tree_electron_z" ,     &tree_electron_z);
    smalltree->Branch("tree_electron_energy", &tree_electron_energy);
    smalltree->Branch("tree_electron_charge", &tree_electron_charge);
    
    // muons info
    smalltree->Branch("tree_muon_pt"  ,       &tree_muon_pt);
    smalltree->Branch("tree_muon_eta" ,       &tree_muon_eta);
    smalltree->Branch("tree_muon_phi" ,       &tree_muon_phi);
    smalltree->Branch("tree_muon_x"  ,        &tree_muon_x);
    smalltree->Branch("tree_muon_y" ,         &tree_muon_y);
    smalltree->Branch("tree_muon_z" ,         &tree_muon_z);
    smalltree->Branch("tree_muon_energy",     &tree_muon_energy);
    smalltree->Branch("tree_muon_dxy",        &tree_muon_dxy);
    smalltree->Branch("tree_muon_dxyError",   &tree_muon_dxyError);
    smalltree->Branch("tree_muon_dz",         &tree_muon_dz);
    smalltree->Branch("tree_muon_dzError",    &tree_muon_dzError);
    smalltree->Branch("tree_muon_charge",     &tree_muon_charge);
    smalltree->Branch("tree_muon_isLoose",    &tree_muon_isLoose);
    smalltree->Branch("tree_muon_isTight",    &tree_muon_isTight);
    smalltree->Branch("tree_muon_isGlobal",   &tree_muon_isGlobal);
    smalltree->Branch("tree_muon_PFisoVeryTight",         &tree_muon_PFisoVeryTight);
    smalltree->Branch("tree_muon_PFisoTight",             &tree_muon_PFisoTight);
    smalltree->Branch("tree_muon_PFisoMedium",            &tree_muon_PFisoMedium);
    smalltree->Branch("tree_muon_PFisoLoose",             &tree_muon_PFisoLoose);
    smalltree->Branch("tree_muon_MVAisoLoose",            &tree_muon_MVAisoLoose);
    smalltree->Branch("tree_muon_MVAisoMedium",           &tree_muon_MVAisoMedium);
    smalltree->Branch("tree_muon_MVAisoTight",            &tree_muon_MVAisoTight);
    smalltree->Branch("tree_muon_isStandAloneMuon",       &tree_muon_isStandAloneMuon);
    smalltree->Branch("tree_muon_CutBasedIdMedium",       &tree_muon_CutBasedIdMedium       );
    smalltree->Branch("tree_muon_CutBasedIdMediumPrompt", &tree_muon_CutBasedIdMediumPrompt );
    
    // track
    smalltree->Branch("tree_nTracks",            &tree_nTracks, "tree_nTracks/I"); 
    smalltree->Branch("tree_track_pz",&tree_track_pz);
    smalltree->Branch("tree_track_pt",           &tree_track_pt );
    smalltree->Branch("tree_track_eta",          &tree_track_eta );
    smalltree->Branch("tree_track_phi",          &tree_track_phi );
    smalltree->Branch("tree_track_charge",       &tree_track_charge );
    smalltree->Branch("tree_track_NChi2",        &tree_track_NChi2 );
    smalltree->Branch("tree_track_isHighPurity", &tree_track_isHighPurity );
    smalltree->Branch("tree_track_dxy",          &tree_track_dxy );
    smalltree->Branch("tree_track_dxyError",     &tree_track_dxyError );
    smalltree->Branch("tree_track_dz",           &tree_track_dz );
    smalltree->Branch("tree_track_dzError",      &tree_track_dzError );
    smalltree->Branch("tree_track_nHit",         &tree_track_nHit );
    smalltree->Branch("tree_track_nHitPixel",    &tree_track_nHitPixel);
    smalltree->Branch("tree_track_nHitTIB",      &tree_track_nHitTIB);
    smalltree->Branch("tree_track_nHitTID",      &tree_track_nHitTID);
    smalltree->Branch("tree_track_nHitTOB",      &tree_track_nHitTOB);
    smalltree->Branch("tree_track_nHitTEC",      &tree_track_nHitTEC);
    smalltree->Branch("tree_track_nHitPXB",      &tree_track_nHitPXB);
    smalltree->Branch("tree_track_nHitPXF",      &tree_track_nHitPXF);
    smalltree->Branch("tree_track_isHitPixel",   &tree_track_isHitPixel);
    smalltree->Branch("tree_track_nLayers",      &tree_track_nLayers);
    smalltree->Branch("tree_track_nLayersPixel", &tree_track_nLayersPixel);

    
    smalltree->Branch("tree_track_x",            &tree_track_x );
    smalltree->Branch("tree_track_y",            &tree_track_y );
    smalltree->Branch("tree_track_z",            &tree_track_z );
    smalltree->Branch("tree_track_firstHit",     &tree_track_firstHit);
    smalltree->Branch("tree_track_firstHit_x",   &tree_track_firstHit_x);
    smalltree->Branch("tree_track_firstHit_y",   &tree_track_firstHit_y);
    smalltree->Branch("tree_track_firstHit_z",   &tree_track_firstHit_z);
    smalltree->Branch("tree_track_firstHit_phi", &tree_track_firstHit_phi);

    //
    smalltree->Branch("tree_track_ddxy", &tree_track_ddxy);
    smalltree->Branch("tree_track_ddxyError", &tree_track_ddxyError);
    smalltree->Branch("tree_track_dpt",&tree_track_dpt);
    smalltree->Branch("tree_track_deta",&tree_track_deta);
    smalltree->Branch("tree_track_dphi",&tree_track_dphi);
    smalltree->Branch("tree_track_dNChi2",&tree_track_dNChi2);
    smalltree->Branch("tree_track_dvx",&tree_track_dvx);
    smalltree->Branch("tree_track_dvy",&tree_track_dvy);
    smalltree->Branch("tree_track_dvz",&tree_track_dvz);
    smalltree->Branch("tree_track_ddz",&tree_track_ddz);
    smalltree->Branch("tree_track_ddzError",&tree_track_ddzError);
    
//$$$$
    smalltree->Branch("tree_track_extraTrue_dx",  &tree_track_extraTrue_dx);
    smalltree->Branch("tree_track_extraTrue_dy",  &tree_track_extraTrue_dy);
    smalltree->Branch("tree_track_extraTrue_dz",  &tree_track_extraTrue_dz);
    smalltree->Branch("tree_track_extraTrue_dphi",&tree_track_extraTrue_dphi);
    smalltree->Branch("tree_track_fromLayer_dx",  &tree_track_fromLayer_dx);
    smalltree->Branch("tree_track_fromLayer_dy",  &tree_track_fromLayer_dy);
    smalltree->Branch("tree_track_fromLayer_dz",  &tree_track_fromLayer_dz);
    smalltree->Branch("tree_track_fromLayer_dphi",&tree_track_fromLayer_dphi);
//$$$$
    smalltree->Branch("tree_track_iJet",         &tree_track_iJet);
    smalltree->Branch("tree_track_ntrk10",       &tree_track_ntrk10);
    smalltree->Branch("tree_track_ntrk20",       &tree_track_ntrk20);
    smalltree->Branch("tree_track_ntrk30",       &tree_track_ntrk30);
    smalltree->Branch("tree_track_MVAval",       &tree_track_MVAval);
    smalltree->Branch("tree_track_Hemi",         &tree_track_Hemi);
    smalltree->Branch("tree_track_Hemi_dR",      &tree_track_Hemi_dR);
    smalltree->Branch("tree_track_Hemi_mva_NChi2", &tree_track_Hemi_mva_NChi2);
    smalltree->Branch("tree_track_Hemi_LLP",     &tree_track_Hemi_LLP);

    smalltree->Branch("tree_track_ptError",           &tree_track_ptError );
    smalltree->Branch("tree_track_outerPt",           &tree_track_outerPt );
    smalltree->Branch("tree_track_isLoose",           &tree_track_isLoose );
    smalltree->Branch("tree_track_isTight",           &tree_track_isTight );
    smalltree->Branch("tree_track_numberOfLostHits",  &tree_track_numberOfLostHits );
    smalltree->Branch("tree_track_originalAlgo",      &tree_track_originalAlgo );
    smalltree->Branch("tree_track_algo",              &tree_track_algo );
    smalltree->Branch("tree_track_stopReason",        &tree_track_stopReason );
    smalltree->Branch("tree_track_isSimMatched",      &tree_track_isSimMatched );
    smalltree->Branch("tree_track_stripTECLayersWithMeasurement", &tree_track_stripTECLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTIBLayersWithMeasurement", &tree_track_stripTIBLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTIDLayersWithMeasurement", &tree_track_stripTIDLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTOBLayersWithMeasurement", &tree_track_stripTOBLayersWithMeasurement);
    smalltree->Branch("tree_track_recoVertex_idx",    &tree_track_recoVertex_idx);
    smalltree->Branch("tree_track_recoAK4PFJet_idx",  &tree_track_recoAK4PFJet_idx);
    smalltree->Branch("tree_track_reco08Jet_idx",     &tree_track_reco08Jet_idx);
    smalltree->Branch("tree_track_recoCaloJet_idx",   &tree_track_recoCaloJet_idx);
    smalltree->Branch("tree_track_drSig",              &tree_track_drSig );
    //!!!!
    smalltree->Branch("tree_track_GeoBarrel_firsthit_x",&tree_track_GeoBarrel_firsthit_x);
    smalltree->Branch("tree_track_GeoBarrel_firsthit_y",&tree_track_GeoBarrel_firsthit_y);
    smalltree->Branch("tree_track_GeoBarrel_firsthit_z",&tree_track_GeoBarrel_firsthit_z);
    smalltree->Branch("tree_track_RECOvsMINI_GeoBarrel_firsthit_x",&tree_track_RECOvsMINI_GeoBarrel_firsthit_x);
    smalltree->Branch("tree_track_RECOvsMINI_GeoBarrel_firsthit_y",&tree_track_RECOvsMINI_GeoBarrel_firsthit_y);
    smalltree->Branch("tree_track_RECOvsMINI_GeoBarrel_firsthit_z",&tree_track_RECOvsMINI_GeoBarrel_firsthit_z);

    smalltree->Branch("tree_track_PropBarrel_firsthit_x",&tree_track_PropBarrel_firsthit_x);
    smalltree->Branch("tree_track_PropBarrel_firsthit_y",&tree_track_PropBarrel_firsthit_y);
    smalltree->Branch("tree_track_PropBarrel_firsthit_z",&tree_track_PropBarrel_firsthit_z);
    smalltree->Branch("tree_track_RECOvsMINI_PropBarrel_firsthit_x",&tree_track_RECOvsMINI_PropBarrel_firsthit_x);
    smalltree->Branch("tree_track_RECOvsMINI_PropBarrel_firsthit_y",&tree_track_RECOvsMINI_PropBarrel_firsthit_y);
    smalltree->Branch("tree_track_RECOvsMINI_PropBarrel_firsthit_z",&tree_track_RECOvsMINI_PropBarrel_firsthit_z);

    smalltree->Branch("tree_track_PropBarrel_firsthit_opti_x",&tree_track_PropBarrel_firsthit_opti_x);
    smalltree->Branch("tree_track_PropBarrel_firsthit_opti_y",&tree_track_PropBarrel_firsthit_opti_y);
    smalltree->Branch("tree_track_PropBarrel_firsthit_opti_z",&tree_track_PropBarrel_firsthit_opti_z);
    smalltree->Branch("tree_track_RECOvsMINI_PropBarrel_firsthit_opti_x",&tree_track_RECOvsMINI_PropBarrel_firsthit_opti_x);
    smalltree->Branch("tree_track_RECOvsMINI_PropBarrel_firsthit_opti_y",&tree_track_RECOvsMINI_PropBarrel_firsthit_opti_y);
    smalltree->Branch("tree_track_RECOvsMINI_PropBarrel_firsthit_opti_z",&tree_track_RECOvsMINI_PropBarrel_firsthit_opti_z);

    smalltree->Branch("tree_track_GeoDisk_firsthit_x",&tree_track_GeoDisk_firsthit_x);
    smalltree->Branch("tree_track_GeoDisk_firsthit_y",&tree_track_GeoDisk_firsthit_y);
    smalltree->Branch("tree_track_GeoDisk_firsthit_z",&tree_track_GeoDisk_firsthit_z);
    smalltree->Branch("tree_track_RECOvsMINI_GeoDisk_firsthit_x",&tree_track_RECOvsMINI_GeoDisk_firsthit_x);
    smalltree->Branch("tree_track_RECOvsMINI_GeoDisk_firsthit_y",&tree_track_RECOvsMINI_GeoDisk_firsthit_y);
    smalltree->Branch("tree_track_RECOvsMINI_GeoDisk_firsthit_z",&tree_track_RECOvsMINI_GeoDisk_firsthit_z);
                            
    smalltree->Branch("tree_track_GeoDisk_firsthit_opti_x",&tree_track_GeoDisk_firsthit_opti_x);
    smalltree->Branch("tree_track_GeoDisk_firsthit_opti_y",&tree_track_GeoDisk_firsthit_opti_y);
    smalltree->Branch("tree_track_GeoDisk_firsthit_opti_z",&tree_track_GeoDisk_firsthit_opti_z);
    smalltree->Branch("tree_track_RECOvsMINI_GeoDisk_firsthit_opti_x",&tree_track_RECOvsMINI_GeoDisk_firsthit_opti_x);
    smalltree->Branch("tree_track_RECOvsMINI_GeoDisk_firsthit_opti_y",&tree_track_RECOvsMINI_GeoDisk_firsthit_opti_y);
    smalltree->Branch("tree_track_RECOvsMINI_GeoDisk_firsthit_opti_z",&tree_track_RECOvsMINI_GeoDisk_firsthit_opti_z);

    smalltree->Branch("tree_track_PropDisk_SLCC_firsthit_x",&tree_track_PropDisk_SLCC_firsthit_x);
    smalltree->Branch("tree_track_PropDisk_SLCC_firsthit_y",&tree_track_PropDisk_SLCC_firsthit_y);
    smalltree->Branch("tree_track_PropDisk_SLCC_firsthit_z",&tree_track_PropDisk_SLCC_firsthit_z);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_x",&tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_x);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_y",&tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_y);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_z",&tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_z);
                                            
    smalltree->Branch("tree_track_PropDisk_SLCC_firsthit_opti_x",&tree_track_PropDisk_SLCC_firsthit_opti_x);
    smalltree->Branch("tree_track_PropDisk_SLCC_firsthit_opti_y",&tree_track_PropDisk_SLCC_firsthit_opti_y);
    smalltree->Branch("tree_track_PropDisk_SLCC_firsthit_opti_z",&tree_track_PropDisk_SLCC_firsthit_opti_z);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_x",&tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_x);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_y",&tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_y);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_z",&tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_z);

                                
    smalltree->Branch("tree_track_PropDisk_SLPC_firsthit_x",&tree_track_PropDisk_SLPC_firsthit_x);
    smalltree->Branch("tree_track_PropDisk_SLPC_firsthit_y",&tree_track_PropDisk_SLPC_firsthit_y);
    smalltree->Branch("tree_track_PropDisk_SLPC_firsthit_z",&tree_track_PropDisk_SLPC_firsthit_z);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_x",&tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_x);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_y",&tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_y);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_z",&tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_z);

    smalltree->Branch("tree_track_PropDisk_SLPC_firsthit_opti_x",&tree_track_PropDisk_SLPC_firsthit_opti_x);
    smalltree->Branch("tree_track_PropDisk_SLPC_firsthit_opti_y",&tree_track_PropDisk_SLPC_firsthit_opti_y);
    smalltree->Branch("tree_track_PropDisk_SLPC_firsthit_opti_z",&tree_track_PropDisk_SLPC_firsthit_opti_z);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_x",&tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_x);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_y",&tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_y);
    smalltree->Branch("tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_z",&tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_z);

    smalltree->Branch("tree_track_theta",&tree_track_theta);
    smalltree->Branch("tree_track_surf",&tree_track_surf);
    smalltree->Branch("tree_track_hitpattern", &tree_track_hitpattern);
    //!!!!

    // info about the simulated track matched to the reco track
    smalltree->Branch("tree_track_sim_LLP",      &tree_track_sim_LLP );
    smalltree->Branch("tree_track_sim_isFromB",  &tree_track_sim_isFromB );
    smalltree->Branch("tree_track_sim_isFromC",  &tree_track_sim_isFromC );
    smalltree->Branch("tree_track_sim_pt",       &tree_track_sim_pt );
    smalltree->Branch("tree_track_sim_eta",      &tree_track_sim_eta  );
    smalltree->Branch("tree_track_sim_phi",      &tree_track_sim_phi  );
    smalltree->Branch("tree_track_sim_charge",   &tree_track_sim_charge );
    smalltree->Branch("tree_track_sim_pdgId",  	 &tree_track_sim_pdgId );
    smalltree->Branch("tree_track_sim_mass",     &tree_track_sim_mass   );
    smalltree->Branch("tree_track_sim_x",        &tree_track_sim_x );
    smalltree->Branch("tree_track_sim_y",        &tree_track_sim_y );
    smalltree->Branch("tree_track_sim_z",        &tree_track_sim_z );
    
    smalltree->Branch("tree_track_sim_longLived",             &tree_track_sim_longLived );
    smalltree->Branch("tree_track_sim_numberOfTrackerHits",   &tree_track_sim_numberOfTrackerHits   );
    smalltree->Branch("tree_track_sim_numberOfTrackerLayers", &tree_track_sim_numberOfTrackerLayers );
    smalltree->Branch("tree_track_sim_status",                &tree_track_sim_status );
    smalltree->Branch("tree_track_sim_isFromBC_mother_pdgId", &tree_track_sim_isFromBC_mother_pdgId );
    smalltree->Branch("tree_track_sim_dV",                    &tree_track_sim_dV );
    smalltree->Branch("tree_track_sim_dist",                  &tree_track_sim_dist );

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
    
    // gen info
    smalltree->Branch("tree_GenPVx" ,  &tree_GenPVx);
    smalltree->Branch("tree_GenPVy" ,  &tree_GenPVy);
    smalltree->Branch("tree_GenPVz" ,  &tree_GenPVz);
    
    smalltree->Branch("tree_genParticle_pt"  ,          &tree_genParticle_pt);
    smalltree->Branch("tree_genParticle_eta" ,          &tree_genParticle_eta);
    smalltree->Branch("tree_genParticle_phi" ,          &tree_genParticle_phi);
    smalltree->Branch("tree_genParticle_charge" ,       &tree_genParticle_charge);
    smalltree->Branch("tree_genParticle_pdgId" ,        &tree_genParticle_pdgId);
    smalltree->Branch("tree_genParticle_x"  ,           &tree_genParticle_x);
    smalltree->Branch("tree_genParticle_y" ,            &tree_genParticle_y);
    smalltree->Branch("tree_genParticle_z" ,            &tree_genParticle_z);
    smalltree->Branch("tree_genParticle_mass" ,         &tree_genParticle_mass);
    smalltree->Branch("tree_genParticle_statusCode",    &tree_genParticle_statusCode);
    smalltree->Branch("tree_genParticle_mother_pdgId" , &tree_genParticle_mother_pdgId);
    smalltree->Branch("tree_genParticle_LLP" ,          &tree_genParticle_LLP);

    smalltree->Branch("tree_nFromC",                 &tree_nFromC,  "tree_nFromC/I");
    smalltree->Branch("tree_genFromC_pt"  ,          &tree_genFromC_pt);
    smalltree->Branch("tree_genFromC_eta" ,          &tree_genFromC_eta);
    smalltree->Branch("tree_genFromC_phi" ,          &tree_genFromC_phi);
    smalltree->Branch("tree_genFromC_charge" ,       &tree_genFromC_charge);
    smalltree->Branch("tree_genFromC_pdgId" ,        &tree_genFromC_pdgId);
    smalltree->Branch("tree_genFromC_x"  ,           &tree_genFromC_x);
    smalltree->Branch("tree_genFromC_y" ,            &tree_genFromC_y);
    smalltree->Branch("tree_genFromC_z" ,            &tree_genFromC_z);
    smalltree->Branch("tree_genFromC_mother_pdgId" , &tree_genFromC_mother_pdgId);
    smalltree->Branch("tree_genFromC_generation" ,   &tree_genFromC_generation);
    smalltree->Branch("tree_genFromC_LLP" ,          &tree_genFromC_LLP);

    smalltree->Branch("tree_nFromB",                 &tree_nFromB,  "tree_nFromB/I");
    smalltree->Branch("tree_genFromB_pt"  ,	     &tree_genFromB_pt);
    smalltree->Branch("tree_genFromB_eta" ,	     &tree_genFromB_eta);
    smalltree->Branch("tree_genFromB_phi" ,	     &tree_genFromB_phi);
    smalltree->Branch("tree_genFromB_charge" ,	     &tree_genFromB_charge);
    smalltree->Branch("tree_genFromB_pdgId" ,	     &tree_genFromB_pdgId);
    smalltree->Branch("tree_genFromB_x"  ,	     &tree_genFromB_x);
    smalltree->Branch("tree_genFromB_y" ,	     &tree_genFromB_y);
    smalltree->Branch("tree_genFromB_z" ,	     &tree_genFromB_z);
    smalltree->Branch("tree_genFromB_mother_pdgId" , &tree_genFromB_mother_pdgId);
    smalltree->Branch("tree_genFromB_generation" ,   &tree_genFromB_generation);
    smalltree->Branch("tree_genFromB_LLP" ,          &tree_genFromB_LLP);
    
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
    
    smalltree->Branch("tree_nLLP",          &tree_nLLP);
    smalltree->Branch("tree_LLP",           &tree_LLP);
    smalltree->Branch("tree_LLP_pt" ,       &tree_LLP_pt);
    smalltree->Branch("tree_LLP_eta",       &tree_LLP_eta);
    smalltree->Branch("tree_LLP_phi",       &tree_LLP_phi);
    smalltree->Branch("tree_LLP_x",         &tree_LLP_x);
    smalltree->Branch("tree_LLP_y",         &tree_LLP_y);
    smalltree->Branch("tree_LLP_z",         &tree_LLP_z);
    smalltree->Branch("tree_LLP_nTrks",     &tree_LLP_nTrks);
    smalltree->Branch("tree_LLP_Vtx_nTrks", &tree_LLP_Vtx_nTrks);
    smalltree->Branch("tree_LLP_Vtx_NChi2", &tree_LLP_Vtx_NChi2);
    smalltree->Branch("tree_LLP_Vtx_dx",    &tree_LLP_Vtx_dx);
    smalltree->Branch("tree_LLP_Vtx_dy",    &tree_LLP_Vtx_dy);
    smalltree->Branch("tree_LLP_Vtx_dz",    &tree_LLP_Vtx_dz);

    smalltree->Branch("tree_Hemi",           &tree_Hemi);
    smalltree->Branch("tree_Hemi_njet",      &tree_Hemi_njet);
    smalltree->Branch("tree_Hemi_eta",       &tree_Hemi_eta);
    smalltree->Branch("tree_Hemi_phi",       &tree_Hemi_phi);
    smalltree->Branch("tree_Hemi_dR",        &tree_Hemi_dR);
    smalltree->Branch("tree_Hemi_nTrks",     &tree_Hemi_nTrks);
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
    // do anything here that needs to be done at destruction time
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
  bool showlog = false;
//$$
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  runNumber = iEvent.id().run();
  //std::cout << "runNumber = " << runNumber << std::endl;
  eventNumber = iEvent.id().event();
  //std::cout << "eventNumber = "<< eventNumber <<std::endl;
  lumiBlock = iEvent.luminosityBlock();
  
//$$$$
    bool Central[2000];
    float Layer[2000];
    for (unsigned int i = 0; i<2000; i++) {
      Central[i] = false;
      Layer[i] = 0;
    }  
    Central[1160] = true;
    Central[1168] = true;
    Central[1176] = true;
    Central[1184] = true;
    Central[1416] = true;
    Central[1420] = true;
    Central[1424] = true;
    Central[1428] = true;
    Central[1432] = true;
    Central[1440] = true;
    Central[1672] = true;
    Layer[1160] =   3.0;
    Layer[1168] =   6.8;
    Layer[1176] =  10.9;
    Layer[1184] =  16.0;
    Layer[1288] =  32.3;
    Layer[1296] =  39.3;
    Layer[1304] =  48.9;
    Layer[1416] =  24.0;
    Layer[1420] =  27.0;
    Layer[1424] =  32.4;
    Layer[1428] =  35.3;
    Layer[1432] =  41.7;
    Layer[1440] =  49.7;
    Layer[1544] =  77.9;
    Layer[1548] =  80.8;
    Layer[1552] =  90.4;
    Layer[1556] =  94.2;
    Layer[1560] = 102.7;
    Layer[1564] = 107.2;
    Layer[1672] =  60.5;
    Layer[1800] = 131.6;
    Layer[1804] = 129.5;
    Layer[1808] = 145.7;
    Layer[1812] = 143.0;
    Layer[1816] = 160.2;
    Layer[1820] = 157.8;
    Layer[1824] = 174.4;
    Layer[1828] = 173.1;
    Layer[1832] = 188.7;
    Layer[1836] = 186.0;
//$$$$


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
  
  iEvent.getByToken(ak4slimmedJetToken_,ak4slimmedJets);
  iEvent.getByToken(ak4PFJetToken_,ak4PFJets);
  iEvent.getByToken(CaloJetToken_,CaloJets);
  iEvent.getByToken(ak8jetToken_,ak8jets);
  
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
  const pat::PackedCandidateCollection* pc = pfcands.product();    
  
  edm::Handle<pat::ElectronCollection> electronsPAT;
  iEvent.getByToken(electronPATToken_,electronsPAT);
  
  edm::Handle<pat::MuonCollection> slimmedmuons;
  iEvent.getByToken(slimmedmuonToken_,slimmedmuons);
  
  edm::Handle<vector<reco::CompositeCandidate> > dimuons;

  edm::Handle<pat::PackedCandidateCollection> pcs;
  iEvent.getByToken(pfcandsToken_, pcs);
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
    
  for (auto const & vertex : *vertices) {
      tree_vtx_PosX.push_back(vertex.x());
      tree_vtx_PosY.push_back(vertex.y());
      tree_vtx_PosZ.push_back(vertex.z());
      tree_vtx_NChi2.push_back(vertex.normalizedChi2());
      tree_vtx_PosXError.push_back(vertex.xError());
      tree_vtx_PosYError.push_back(vertex.yError());
      tree_vtx_PosZError.push_back(vertex.zError());
      tree_vtx_ndf.push_back(vertex.ndof());
  }

  tree_nPV = vertices->size();
  if ( !vertices->empty() ) {
    tree_PV_x     = tree_vtx_PosX[0]; // l'index 0 donne le PV!
    tree_PV_y     = tree_vtx_PosY[0];
    tree_PV_z     = tree_vtx_PosZ[0];
    tree_PV_ez    = tree_vtx_PosZError[0];
    tree_PV_NChi2 = tree_vtx_NChi2[0];
    tree_PV_ndf   = tree_vtx_ndf[0];
  }
  const reco::Vertex &PV = vertices->front();

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
  tree_njet = 0;
  float HT_val = 0;
  float jet_pt_min = 20.;
  for (int ij=0;ij<int(ak4slimmedJets->size());ij++)
  {
    const Jet& jet = ak4slimmedJets->at(ij);
  if ( jet.pt() < jet_pt_min ) continue;
    tree_jet_E.push_back(jet.energy());
    tree_jet_pt.push_back(jet.pt());
    tree_jet_eta.push_back(jet.eta());
    tree_jet_phi.push_back(jet.phi());
    tree_njet++;
    if ( abs(jet.eta()) < 2.4 ) HT_val += jet.pt(); // used in HT filter !
//$$
//       TLorentzVector thejet, track;
//       thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
//       for (unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++)
//       {
//   	  track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
//   	  if(thejet.DeltaR(track) < 0.4 ) tree_jet_idxTrack.push_back(iTrack);
//       }
//$$
  }
  
  for (int ij=0;ij<int(ak4PFJets->size());ij++)
  {
      const Jet& jet = ak4PFJets->at(ij);
  if ( jet.pt() < jet_pt_min ) continue;
      tree_AK4PFjet_E.push_back(jet.energy());
      tree_AK4PFjet_pt.push_back(jet.pt());
      tree_AK4PFjet_eta.push_back(jet.eta());
      tree_AK4PFjet_phi.push_back(jet.phi());
//$$
//       TLorentzVector thejet, track;
//       thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
//       for (unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++)
//       {
//   	  track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
//   	  if(thejet.DeltaR(track) < 0.4 ) tree_AK4PFjet_idxTrack.push_back(iTrack);
//       }
//$$
  }
  
  for (int ij=0;ij< int(CaloJets->size());ij++)//plus simpliste, utilisete just eles cellules calorimetriques
  {
      const Jet& jet = CaloJets->at(ij);
  if ( jet.pt() < jet_pt_min ) continue;
      tree_CaloJet_E.push_back(jet.energy());
      tree_CaloJet_pt.push_back(jet.pt());
      tree_CaloJet_eta.push_back(jet.eta());
      tree_CaloJet_phi.push_back(jet.phi());
//$$
//       TLorentzVector thejet, track;
//       thejet.SetPtEtaPhiE(jet.pt(),jet.eta(), jet.phi(), jet.energy() );
//       for (unsigned int iTrack=0; iTrack< tree_track_pt.size() ; iTrack++) {
//   	  track.SetPtEtaPhiM(tree_track_pt[iTrack],tree_track_eta[iTrack], tree_track_phi[iTrack], 0 );
//   	  if(thejet.DeltaR(track) < 0.4 ) tree_CaloJet_idxTrack.push_back(iTrack);
//       }
//$$
  }
  
  for (int ij=0;ij< int(ak8jets->size() );ij++)
  {
      const Jet& jet = ak8jets->at(ij);
  if ( jet.pt() < jet_pt_min ) continue;
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
  if ( electron.pt() < 5. ) continue;
      tree_electron_pt.push_back(electron.pt());
      tree_electron_eta.push_back(electron.eta());
      tree_electron_phi.push_back(electron.phi());
      tree_electron_x.push_back(electron.vx());
      tree_electron_y.push_back(electron.vy());
      tree_electron_z.push_back(electron.vz());
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
      tree_muon_pt.push_back(       muon.pt());
      tree_muon_eta.push_back(      muon.eta());
      tree_muon_phi.push_back(      muon.phi());
      tree_muon_x.push_back(        muon.vx());
      tree_muon_y.push_back(        muon.vy());
      tree_muon_z.push_back(        muon.vz());
      tree_muon_energy.push_back(   muon.energy());
      tree_muon_dxy.push_back(	    muon.bestTrack()->dxy(PV.position()));
      tree_muon_dxyError.push_back( muon.bestTrack()->dxyError());
      tree_muon_dz.push_back(	    muon.bestTrack()->dz(PV.position()));
      tree_muon_dzError.push_back(  muon.bestTrack()->dzError());
      tree_muon_charge.push_back(   muon.charge());
      tree_muon_isLoose.push_back(  muon.passed(reco::Muon::CutBasedIdLoose));
      tree_muon_isTight.push_back(  muon.passed(reco::Muon::CutBasedIdTight));
      tree_muon_isGlobal.push_back( muon.passed(muon.isGlobalMuon()));
      
      tree_muon_PFisoVeryTight.push_back(         muon.passed(reco::Muon::PFIsoVeryTight));
      tree_muon_PFisoTight.push_back(             muon.passed(reco::Muon::PFIsoTight) );
      tree_muon_PFisoMedium.push_back(            muon.passed(reco::Muon::PFIsoMedium));
      tree_muon_PFisoLoose.push_back(             muon.passed(reco::Muon::PFIsoLoose ));
      tree_muon_MVAisoLoose.push_back(            muon.passed(reco::Muon::MvaLoose )  );
      tree_muon_MVAisoMedium.push_back(           muon.passed(reco::Muon::MvaMedium)  );
      tree_muon_MVAisoTight.push_back(            muon.passed(reco::Muon::MvaTight )  );
      tree_muon_isStandAloneMuon.push_back(       muon.passed(muon.isStandAloneMuon()));
      tree_muon_CutBasedIdMedium.push_back(       muon.passed(reco::Muon::CutBasedIdMedium));
      tree_muon_CutBasedIdMediumPrompt.push_back( muon.passed(reco::Muon::CutBasedIdMediumPrompt));
      
      if ( !muon.passed(muon.isGlobalMuon()) ) continue; //Need global muons
      mupt1  = muon.pt();
      if ( mupt1 < 10. ) continue; //Zmu filter
      dVr = TMath::Sqrt( (muon.vx()-tree_PV_x)*(muon.vx()-tree_PV_x) + (muon.vy()-tree_PV_y)*(muon.vy()-tree_PV_y) );
      dVz = muon.vz()-tree_PV_z;
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
  	  if ( !tree_muon_isGlobal[muon2_idx] ) continue;
  	  if ( tree_muon_charge[muon_idx] == tree_muon_charge[muon2_idx] ) continue;
  	  dVr = TMath::Sqrt( (tree_muon_x[muon2_idx]-tree_PV_x)*(tree_muon_x[muon2_idx]-tree_PV_x) + (tree_muon_y[muon2_idx]-tree_PV_y)*(tree_muon_y[muon2_idx]-tree_PV_y) );
  	  dVz = tree_muon_z[muon2_idx]-tree_PV_z;
  	  if ( dVr > 0.1 || abs(dVz) > 0.2 ) continue;
  	  mupt2  = tree_muon_pt[muon2_idx];
  	  if ( mupt2 < 10. ) continue;
  	  if ( mupt1 < 28. && mupt2 < 28. ) continue; //Zmu FIlter
  	  mueta2 = tree_muon_eta[muon2_idx];
  	  muphi2 = tree_muon_phi[muon2_idx];
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
  //  std::cout<< "number of muons first method : " << tree_muon_pt.size() << std::endl;

  if ( tree_muon_pt[imu2] > tree_muon_pt[imu1] ) {
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
  tree_nFromC = 0; 
  tree_nFromB = 0;
      
  // Gen Information  for event axis //
  float  Gen_neu1_eta=-10, Gen_neu1_phi=-10;
  float  Gen_neu2_eta=-10, Gen_neu2_phi=-10;
  float dRneuneu = 0.;
  int  neu[2],  nneu = 0;
  TLorentzVector vneu[2];

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
	    LLP2_dist = TMath::Sqrt( (LLP2_x - tree_GenPVx)*(LLP2_x - tree_GenPVx) 
				   + (LLP2_y - tree_GenPVy)*(LLP2_y - tree_GenPVy) 
				   + (LLP2_z - tree_GenPVz)*(LLP2_z - tree_GenPVz) ); 
	  }
	}
	if ( nllp == 0 ) {
	  nllp = 1;
	  LLP1_x = genParticle.vx();
	  LLP1_y = genParticle.vy();
	  LLP1_z = genParticle.vz();
	  LLP1_dist = TMath::Sqrt( (LLP1_x - tree_GenPVx)*(LLP1_x - tree_GenPVx) 
				 + (LLP1_y - tree_GenPVy)*(LLP1_y - tree_GenPVy) 
				 + (LLP1_z - tree_GenPVz)*(LLP1_z - tree_GenPVz) ); 
	}
      }
      
      float dV0 = (genParticle.vx() - tree_GenPVx)*(genParticle.vx() - tree_GenPVx)
        	+ (genParticle.vy() - tree_GenPVy)*(genParticle.vy() - tree_GenPVy)
        	+ (genParticle.vz() - tree_GenPVz)*(genParticle.vz() - tree_GenPVz);
      float dV1 = (genParticle.vx() - LLP1_x)     *(genParticle.vx() - LLP1_x)
        	+ (genParticle.vy() - LLP1_y)     *(genParticle.vy() - LLP1_y)
        	+ (genParticle.vz() - LLP1_z)     *(genParticle.vz() - LLP1_z);
      float dV2 = (genParticle.vx() - LLP2_x)     *(genParticle.vx() - LLP2_x)
        	+ (genParticle.vy() - LLP2_y)     *(genParticle.vy() - LLP2_y)
        	+ (genParticle.vz() - LLP2_z)     *(genParticle.vz() - LLP2_z);
      int fromLLP = -1;
      if      ( dV1 < dV2 && dV1 < 0.01 ) fromLLP = 1;
      else if ( dV2 < dV1 && dV2 < 0.01 ) fromLLP = 2;
      else if ( dV0 < 0.01 )		  fromLLP = 0;

    if (genParticle.pt() < 0.9 || fabs(genParticle.eta()) > 4.0) continue;
      
      tree_genParticle_pt.push_back(           genParticle.pt());
      tree_genParticle_eta.push_back(          genParticle.eta());
      tree_genParticle_phi.push_back(          genParticle.phi());
      tree_genParticle_charge.push_back(       genParticle.charge());
      tree_genParticle_pdgId.push_back(        genParticle.pdgId());
      tree_genParticle_x.push_back(            genParticle.vx());
      tree_genParticle_y.push_back(            genParticle.vy());
      tree_genParticle_z.push_back(            genParticle.vz());
      tree_genParticle_mass.push_back(         genParticle.mass());
      tree_genParticle_statusCode.push_back(   genParticle.status());
      tree_genParticle_mother_pdgId.push_back( mom ? mom->pdgId() :  -1 );
      tree_genParticle_LLP.push_back(          fromLLP);

    } // end loop on gen particles
    
    tree_nLLP = nllp;
    
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
              tree_nFromC++;
              tree_genFromC_pt.push_back(    gen1->pt());
              tree_genFromC_eta.push_back(   gen1->eta());
              tree_genFromC_phi.push_back(   gen1->phi());
              tree_genFromC_charge.push_back(gen1->charge());
              tree_genFromC_pdgId.push_back( gen1->pdgId());
              tree_genFromC_x.push_back(    gen1->vx());
              tree_genFromC_y.push_back(    gen1->vy());
              tree_genFromC_z.push_back(    gen1->vz());
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
                  tree_nFromC++;
                  tree_genFromC_pt.push_back(    gen2->pt());
                  tree_genFromC_eta.push_back(   gen2->eta());
                  tree_genFromC_phi.push_back(   gen2->phi());
                  tree_genFromC_charge.push_back(gen2->charge());
                  tree_genFromC_pdgId.push_back( gen2->pdgId());
                  tree_genFromC_x.push_back(    gen2->vx());
                  tree_genFromC_y.push_back(    gen2->vy());
                  tree_genFromC_z.push_back(    gen2->vz());
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
        	      tree_nFromC++;
                      tree_genFromC_pt.push_back(    gen3->pt());
                      tree_genFromC_eta.push_back(   gen3->eta());
                      tree_genFromC_phi.push_back(   gen3->phi());
                      tree_genFromC_charge.push_back(gen3->charge());
                      tree_genFromC_pdgId.push_back( gen3->pdgId());
                      tree_genFromC_x.push_back(    gen3->vx());
                      tree_genFromC_y.push_back(    gen3->vy());
                      tree_genFromC_z.push_back(    gen3->vz());
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
        	          tree_nFromC++;
                          tree_genFromC_pt.push_back(	 gen4->pt());
                          tree_genFromC_eta.push_back(	 gen4->eta());
                          tree_genFromC_phi.push_back(	 gen4->phi());
                          tree_genFromC_charge.push_back(gen4->charge());
                          tree_genFromC_pdgId.push_back( gen4->pdgId());
                          tree_genFromC_x.push_back(	 gen4->vx());
                          tree_genFromC_y.push_back(	 gen4->vy());
                          tree_genFromC_z.push_back(	 gen4->vz());
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
        	              tree_nFromC++;
                              tree_genFromC_pt.push_back(    gen5->pt());
                              tree_genFromC_eta.push_back(   gen5->eta());
                              tree_genFromC_phi.push_back(   gen5->phi());
                              tree_genFromC_charge.push_back(gen5->charge());
                              tree_genFromC_pdgId.push_back( gen5->pdgId());
                              tree_genFromC_x.push_back(    gen5->vx());
                              tree_genFromC_y.push_back(    gen5->vy());
                              tree_genFromC_z.push_back(    gen5->vz());
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
        	                  tree_nFromC++;
                                  tree_genFromC_pt.push_back(    gen6->pt());
                                  tree_genFromC_eta.push_back(   gen6->eta());
                                  tree_genFromC_phi.push_back(   gen6->phi());
                                  tree_genFromC_charge.push_back(gen6->charge());
                                  tree_genFromC_pdgId.push_back( gen6->pdgId());
                                  tree_genFromC_x.push_back(    gen6->vx());
                                  tree_genFromC_y.push_back(    gen6->vy());
                                  tree_genFromC_z.push_back(    gen6->vz());
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
        	                      tree_nFromC++;
                                      tree_genFromC_pt.push_back(    gen7->pt());
                                      tree_genFromC_eta.push_back(   gen7->eta());
                                      tree_genFromC_phi.push_back(   gen7->phi());
                                      tree_genFromC_charge.push_back(gen7->charge());
                                      tree_genFromC_pdgId.push_back( gen7->pdgId());
                                      tree_genFromC_x.push_back(    gen7->vx());
                                      tree_genFromC_y.push_back(    gen7->vy());
                                      tree_genFromC_z.push_back(    gen7->vz());
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
              tree_nFromB++;
              tree_genFromB_pt.push_back(    gen1->pt());
              tree_genFromB_eta.push_back(   gen1->eta());
              tree_genFromB_phi.push_back(   gen1->phi());
              tree_genFromB_charge.push_back(gen1->charge());
              tree_genFromB_pdgId.push_back( gen1->pdgId());
              tree_genFromB_x.push_back(    gen1->vx());
              tree_genFromB_y.push_back(    gen1->vy());
              tree_genFromB_z.push_back(    gen1->vz());
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
                  tree_nFromB++;
                  tree_genFromB_pt.push_back(    gen2->pt());
                  tree_genFromB_eta.push_back(   gen2->eta());
                  tree_genFromB_phi.push_back(   gen2->phi());
                  tree_genFromB_charge.push_back(gen2->charge());
                  tree_genFromB_pdgId.push_back( gen2->pdgId());
                  tree_genFromB_x.push_back(    gen2->vx());
                  tree_genFromB_y.push_back(    gen2->vy());
                  tree_genFromB_z.push_back(    gen2->vz());
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
        	      tree_nFromB++;
                      tree_genFromB_pt.push_back(    gen3->pt());
                      tree_genFromB_eta.push_back(   gen3->eta());
                      tree_genFromB_phi.push_back(   gen3->phi());
                      tree_genFromB_charge.push_back(gen3->charge());
                      tree_genFromB_pdgId.push_back( gen3->pdgId());
                      tree_genFromB_x.push_back(    gen3->vx());
                      tree_genFromB_y.push_back(    gen3->vy());
                      tree_genFromB_z.push_back(    gen3->vz());
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
        	          tree_nFromB++;
                          tree_genFromB_pt.push_back(    gen4->pt());
                          tree_genFromB_eta.push_back(   gen4->eta());
                          tree_genFromB_phi.push_back(   gen4->phi());
                          tree_genFromB_charge.push_back(gen4->charge());
                          tree_genFromB_pdgId.push_back( gen4->pdgId());
                          tree_genFromB_x.push_back(    gen4->vx());
                          tree_genFromB_y.push_back(    gen4->vy());
                          tree_genFromB_z.push_back(    gen4->vz());
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
        	              tree_nFromB++;
                              tree_genFromB_pt.push_back(    gen5->pt());
                              tree_genFromB_eta.push_back(   gen5->eta());
                              tree_genFromB_phi.push_back(   gen5->phi());
                              tree_genFromB_charge.push_back(gen5->charge());
                              tree_genFromB_pdgId.push_back( gen5->pdgId());
                              tree_genFromB_x.push_back(    gen5->vx());
                              tree_genFromB_y.push_back(    gen5->vy());
                              tree_genFromB_z.push_back(    gen5->vz());
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
        	                  tree_nFromB++;
                                  tree_genFromB_pt.push_back(    gen6->pt());
                                  tree_genFromB_eta.push_back(   gen6->eta());
                                  tree_genFromB_phi.push_back(   gen6->phi());
                                  tree_genFromB_charge.push_back(gen6->charge());
                                  tree_genFromB_pdgId.push_back( gen6->pdgId());
                                  tree_genFromB_x.push_back(    gen6->vx());
                                  tree_genFromB_y.push_back(    gen6->vy());
                                  tree_genFromB_z.push_back(    gen6->vz());
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
        	                      tree_nFromB++;
                                      tree_genFromB_pt.push_back(    gen7->pt());
                                      tree_genFromB_eta.push_back(   gen7->eta());
                                      tree_genFromB_phi.push_back(   gen7->phi());
                                      tree_genFromB_charge.push_back(gen7->charge());
                                      tree_genFromB_pdgId.push_back( gen7->pdgId());
                                      tree_genFromB_x.push_back(    gen7->vx());
                                      tree_genFromB_y.push_back(    gen7->vy());
                                      tree_genFromB_z.push_back(    gen7->vz());
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
        	                          tree_nFromB++;
                                          tree_genFromB_pt.push_back(    gen8->pt());
                                          tree_genFromB_eta.push_back(   gen8->eta());
                                          tree_genFromB_phi.push_back(   gen8->phi());
                                          tree_genFromB_charge.push_back(gen8->charge());
                                          tree_genFromB_pdgId.push_back( gen8->pdgId());
                                          tree_genFromB_x.push_back(    gen8->vx());
                                          tree_genFromB_y.push_back(    gen8->vy());
                                          tree_genFromB_z.push_back(    gen8->vz());
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
        	                              tree_nFromB++;
                                              tree_genFromB_pt.push_back(    gen9->pt());
                                              tree_genFromB_eta.push_back(   gen9->eta());
                                              tree_genFromB_phi.push_back(   gen9->phi());
                                              tree_genFromB_charge.push_back(gen9->charge());
                                              tree_genFromB_pdgId.push_back( gen9->pdgId());
                                              tree_genFromB_x.push_back(    gen9->vx());
                                              tree_genFromB_y.push_back(    gen9->vy());
                                              tree_genFromB_z.push_back(    gen9->vz());
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
    

  //////////////////////////////////
  //////////////////////////////////
  //////// HT FILTER CHECK /////////
  //////////////////////////////////
  //////////////////////////////////
  
  bool MuonSup10GeV = false;
  bool MuonSup28GeV = false;
  
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
  
  if ( tree_NbrOfZCand > 0 && HT_val > 180. && MuonSup10GeV && MuonSup28GeV && invMassGood ) {
    tree_passesHTFilter = true;
  }
  else tree_passesHTFilter = false;
//$$$$  tree_nTracks = -1; 
    
//$$$$  if ( tree_passesHTFilter ) {

    //////////////////////////////////
    //////////////////////////////////
    //////////   Tracks  /////////////
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
      if ( jet.pt() < jet_pt_min ) continue;
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
      if ( jet.pt() < jet_pt_min ) continue;
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
      if ( jet.pt() < jet_pt_min ) continue;
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
      if ( jet.pt() < jet_pt_min ) continue;  
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
      
      idxTrack++;
    }
    
    //reco::BeamSpot vertexBeamSpot= *recoBeamSpotHandle;
    edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);
    
    /////////////////
    //prepare association to tracks by hit
    reco::RecoToSimCollection recSimColl ;
    if ( !runOnData_ ) recSimColl= associatorByHits.associateRecoToSim(trackRefs, tpCollection);
    //fake rate determination : when a reco track has no matched simtrack
    
    //    int nUnmatchTrack_fromPU = 0;
    //    int nUnmatchTrack_fromPV = 0;
    //    int nUnmatchTrack= 0;
    //    int nPUTrack= 0;
    
    //---------------
    // loop on tracks
    //---------------

    //----!!!!-----//
    vector<reco::TransientTrack> BestTracks;
    int count =0;//can be useful if a TT cannot be created out of a pat:candidate
    std::vector<std::pair<uint16_t,float> > Players;

    //----------------DataBase for the layers of the tracker that are cocnerned by the analysis--//
    //-uint16_t is the hit pattern information and the float is the mean radius/z of the corresponding layer-//
        Players.push_back(make_pair(1160,2.959));//PIXBL1
        Players.push_back(make_pair(1168,6.778));//PIXBL2
        Players.push_back(make_pair(1176,10.89));//PIXBL3
        Players.push_back(make_pair(1184,16));//PIXBL4
        Players.push_back(make_pair(1416,23.83));//TIBL1
        Players.push_back(make_pair(1420,27.02));//TIBL1stereo
        Players.push_back(make_pair(1424,32.23));//TIBL2
        Players.push_back(make_pair(1428,35.41));//TIBL2stereo
        Players.push_back(make_pair(1432,41.75));//TIBL3
        Players.push_back(make_pair(1440,49.71));//TIBL4
        Players.push_back(make_pair(1672,60.43));//TOBL1
        Players.push_back(make_pair(1288,32.35));//PXFdisk1
        Players.push_back(make_pair(1296,39.41));//PXFdisk2
        Players.push_back(make_pair(1304,48.96));//PXFdisk3
        Players.push_back(make_pair(1544,77.9));//TIDWHeel1
        Players.push_back(make_pair(1548,80.71));//TIDWHeel1stereo
        Players.push_back(make_pair(1552,90.4));//TIDWHeel2
        Players.push_back(make_pair(1556,94.02));//TIDWHeel2stereo
        Players.push_back(make_pair(1560,102.7));//TIDWHeel3
        Players.push_back(make_pair(1564,107.0));//TIDWHeel3stereo
        Players.push_back(make_pair(1800,131.6));//TECWHeel1
        Players.push_back(make_pair(1804,129.4));//TECWHeel1stereo
        Players.push_back(make_pair(1808,145.5));//TECWHeel2
        Players.push_back(make_pair(1812,142.8));//TECWHeel2stereo
        Players.push_back(make_pair(1816,160.1));//TECWHeel3
        Players.push_back(make_pair(1820,157.1));//TECWHeel3stereo
        Players.push_back(make_pair(1824,174.2));//TECWHeel4
        Players.push_back(make_pair(1828,172.7));//TECWHeel4stereo
        Players.push_back(make_pair(1832,188.4));//TECWHeel5
        Players.push_back(make_pair(1836,186.3));//TECWHeel5stereo
        Players.push_back(make_pair(1840,203.2));//TECWHeel6
        Players.push_back(make_pair(1844,203.8));//TECWHeel6stereo
        Players.push_back(make_pair(1848,222.2));//TECWHeel7   
    //-----------------------------------------------------------------//
    //----!!!!-----//
    tree_nTracks = 0; 

    int count3_mini=0;
    for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack) {
      const auto& itTrack = trackRefs[iTrack];
      tree_nTracks++; 
      //------------------------
      //general track properties
      //------------------------
      // std::cout << " pt : " << itTrack->pt()  << std::endl;
      // std::cout << " eta : " << itTrack->eta()  << std::endl;
      // std::cout << " phi : " << itTrack->phi()  << std::endl;
      // std::cout << " dxy : " << itTrack->dxy(PV.position())  << std::endl;
      // std::cout << " dz : " << itTrack->dz()  << std::endl;
      // std::cout << " dxyError : " << itTrack->dxyError()  << std::endl;
      // std::cout << " dzError : " << itTrack->dzError()  << std::endl;
      tree_track_pt.push_back(           itTrack->pt());
      tree_track_outerPt.push_back(      itTrack->outerPt());
      tree_track_ptError.push_back(      itTrack->ptError());
      tree_track_pz.push_back(itTrack->pz());
      tree_track_eta.push_back(          itTrack->eta());
      tree_track_phi.push_back(          itTrack->phi());
      tree_track_charge.push_back(       itTrack->charge());
      tree_track_NChi2.push_back(        itTrack->normalizedChi2());
      tree_track_x.push_back(            itTrack->vx());
      tree_track_y.push_back(            itTrack->vy());
      tree_track_z.push_back(            itTrack->vz());
      // std::cout<<"recoTrack::dxy PV : "<<itTrack->dxy(PV.position()) <<"recoTrack::dxyError PV : "<<itTrack->dxyError(PV.position(),PV.covariance()) <<std::endl;
      // std::cout<<"recoTrack::dxy : "<<itTrack->dxy() <<"recoTrack::dxyError : "<<itTrack->dxyError() <<std::endl;
      // loop on PF Candidates 
    
      tree_track_firstHit_x.push_back(   itTrack->innerPosition().X());
      tree_track_firstHit_y.push_back(   itTrack->innerPosition().Y());
      tree_track_firstHit_z.push_back(   itTrack->innerPosition().Z());
      tree_track_firstHit_phi.push_back( itTrack->innerPosition().phi());
           
      if( itTrack->quality(reco::TrackBase::highPurity) ){tree_track_isHighPurity.push_back(true);}
      else {tree_track_isHighPurity.push_back(false);}
      if( itTrack->quality(reco::TrackBase::loose) )	 {tree_track_isLoose.push_back(true);}
      else {tree_track_isLoose.push_back(false);}
      if( itTrack->quality(reco::TrackBase::tight))	 {tree_track_isTight.push_back(true);}
      else {tree_track_isTight.push_back(false);}
      
      tree_track_dxy.push_back( 	 itTrack->dxy(PV.position()));
      tree_track_dxyError.push_back(	 itTrack->dxyError());
      tree_track_dz.push_back(           itTrack->dz(PV.position()));
      tree_track_dzError.push_back(	 itTrack->dzError());

      tree_track_numberOfLostHits.push_back( itTrack->numberOfLostHits());
      tree_track_originalAlgo.push_back(itTrack->originalAlgo());
      tree_track_algo.push_back(itTrack->algo());
      tree_track_stopReason.push_back(itTrack->stopReason());
      
      //--------------------------------
      //general hit properties of tracks
      //--------------------------------
      
      const reco::HitPattern& hp = itTrack->hitPattern();
      
      tree_track_nHit.push_back(         itTrack->numberOfValidHits());
      tree_track_nHitPixel.push_back(    hp.numberOfValidPixelHits());
//       tree_track_numberOfValidStripHits.push_back(hp.numberOfValidStripHits());
      tree_track_nHitTIB.push_back(      hp.numberOfValidStripTIBHits());
      tree_track_nHitTID.push_back(      hp.numberOfValidStripTIDHits());
      tree_track_nHitTOB.push_back(      hp.numberOfValidStripTOBHits());
      tree_track_nHitTEC.push_back(      hp.numberOfValidStripTECHits());
      tree_track_nHitPXB.push_back(      hp.numberOfValidPixelBarrelHits());
      tree_track_nHitPXF.push_back(      hp.numberOfValidPixelEndcapHits());
      tree_track_nLayers.push_back(      hp.trackerLayersWithMeasurement());
      tree_track_nLayersPixel.push_back( hp.pixelLayersWithMeasurement());
      tree_track_stripTECLayersWithMeasurement.push_back(hp.stripTECLayersWithMeasurement() );
      tree_track_stripTIBLayersWithMeasurement.push_back(hp.stripTIBLayersWithMeasurement());
      tree_track_stripTIDLayersWithMeasurement.push_back(hp.stripTIDLayersWithMeasurement());
      tree_track_stripTOBLayersWithMeasurement.push_back(hp.stripTOBLayersWithMeasurement());
      
      int hitPixelLayer = 0;
      if ( hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) )  hitPixelLayer += 1;
      if ( hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2) )  hitPixelLayer += 10;
      if ( hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3) )  hitPixelLayer += 100;
      if ( hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 4) )  hitPixelLayer += 1000;
      if ( hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1) )  hitPixelLayer += 2;
      if ( hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 2) )  hitPixelLayer += 20;
      if ( hp.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 3) )  hitPixelLayer += 200;
      
      tree_track_isHitPixel.push_back(hitPixelLayer);

      //----!!!!-----//
              //----------------MINIAOD_RECO_COMPARISON-----------//
                  //-----------------IMPORTANT----------------//
                  // BestTracks seems to be destructed somehow//
                  // and therefore cannot be used after-------//
                  // TSOS is said to be better for the -------//
                  // propagators (see Propagator.h)...--------//
                  //------------------------------------------//
      //-----hitpattern -> Database ---/
      uint16_t firsthit = hp.getHitPattern(HitPattern::HitCategory::TRACK_HITS,0);
      tree_track_hitpattern.push_back(firsthit);

      //---Creating State to propagate from  TT---//
      const reco::Track* RtBTracks = trackRefs[iTrack].get();
      BestTracks.push_back(theTransientTrackBuilder->build(RtBTracks));
      const MagneticField* B = BestTracks[count].field();//3.8T
      reco::TransientTrack TT (*RtBTracks,BestTracks[count].field());
      // const FreeTrajectoryState Freetraj = TT.initialFreeState();//Propagator in the barrel can also use FTS
      GlobalPoint vert (itTrack->vx(),itTrack->vy(),itTrack->vz());//Point where the propagation will start
      const TrajectoryStateOnSurface Surtraj = TT.stateOnSurface(vert);//TSOS of this point
      AnalyticalPropagator* Prop = new AnalyticalPropagator(B);//Propagator that will be used for barrel, crashes in the disks when using Plane

      float radius = sqrt(itTrack->innerPosition().X()*itTrack->innerPosition().X()+itTrack->innerPosition().Y()*itTrack->innerPosition().Y());
      //needed for the comparison
      
              //----------Initialitsation for the propagtion------------------//

      float rad = 0;
      float theta_prop=2*atan(exp(-itTrack->eta()));//=> needed for propagation
      tree_track_theta.push_back(theta_prop);
      float zlayers=0;
      TkRotation<float> rot(1,0,0,0,1,0,0,0,1);//Cylinder/Plane are already well-orientated => along/normal to the z-axis
      TrajectoryStateOnSurface PropTSOS;//TSOS for the Barrel
      TrajectoryStateOnSurface PropTSOStrueR;//TSOS for the Barrel for the comparison 

              //---------------------------------------------------------------//

              //--------------StraightLineCylinderCrossing (Propagator) for disks--------------------------//
        GlobalVector GVp (itTrack->vx()+itTrack->px(),itTrack->vy()+itTrack->py(),itTrack->vz()+itTrack->pz());
        GlobalPoint GPv (itTrack->vx(),itTrack->vy(),itTrack->vz());
        LocalVector LVp (itTrack->px(),itTrack->py(),itTrack->pz());
        LocalPoint LPv (0,0,0);
        StraightLineCylinderCrossing SLCC(LPv,LVp, PropagationDirection::alongMomentum);//(localPoint, localvector, )
        GloballyPositioned<float>::PositionType P3D(0.,0.,0.);//GLobal 0
        std::pair<bool,double> spair;
        std::pair<bool,double> spair2;
        std::pair<bool,Basic3DVector<float>> spairPlane;
        std::pair<bool,Basic3DVector<float>> spairPlanetrue;
            //--------------------------------------------------------------------------------------------//
     
            //-------------------------------------------Ease the use for TBrowser----------------------------------------------//
      if (firsthit==1288 || firsthit==1296 || firsthit==1304 || firsthit==1544 || firsthit==1548 || firsthit==1552 || firsthit==1556 || firsthit==1560 || firsthit==1564 || firsthit==1800 || firsthit==1804 || firsthit==1808 || firsthit==1812 || firsthit==1816 || firsthit==1820 || firsthit==1824 || firsthit==1828 || firsthit==1832 || firsthit==1836 || firsthit==1840 ||firsthit== 1844 || firsthit==1848)//supposed to be plane
      { tree_track_surf.push_back(1);}//plane
      else
      { tree_track_surf.push_back(0);}//cylinder
            //------------------------------------------------------------------------------------//    

            //-----------------------------------------Propagation----------------------------------------------//
     
 for (int i=0; i<33;i++)
            {
            if (Players[i].first==firsthit )//layers
                {
//------------------------------------------------Cylinder----------------------------------------------------------//
                    rad = Players[i].second;
                    Cylinder Cylind(rad);//miniaod case
                    PropTSOS = Prop->propagate(Surtraj,Cylind);//works well for barrel
                    Cylinder CylindtrueR(radius);//best case
                    PropTSOStrueR = Prop->propagate(Surtraj,CylindtrueR);//works well for barrel
                    if (i<11 && firsthit!=0 )//cylinder
                        {   
                          //--------------With geom vectors-is a little bit worse than with TSOS => worse case : ~std dev *2----------------------------------//
                            float z0 = (rad+itTrack->vz()*tan(theta_prop))/tan(theta_prop);
                            float x0 = rad*cos(itTrack->phi());
                            float y0 = rad*sin(itTrack->phi());  
                            tree_track_GeoBarrel_firsthit_x.push_back(x0);
                            tree_track_GeoBarrel_firsthit_y.push_back(y0);
                            tree_track_GeoBarrel_firsthit_z.push_back(z0);
                            tree_track_RECOvsMINI_GeoBarrel_firsthit_x.push_back(itTrack->innerPosition().X()-x0);
                            tree_track_RECOvsMINI_GeoBarrel_firsthit_y.push_back(itTrack->innerPosition().Y()-y0);
                            tree_track_RECOvsMINI_GeoBarrel_firsthit_z.push_back(itTrack->innerPosition().Z()-z0);
                            //-----------------------With Propagator TSOS-mean R-----------------------//
                            if (PropTSOS.isValid())
                                {
                            tree_track_PropBarrel_firsthit_x.push_back(PropTSOS.globalPosition().x());
                            tree_track_PropBarrel_firsthit_y.push_back(PropTSOS.globalPosition().y());
                            tree_track_PropBarrel_firsthit_z.push_back(PropTSOS.globalPosition().z());
                            tree_track_RECOvsMINI_PropBarrel_firsthit_x.push_back(itTrack->innerPosition().X()-PropTSOS.globalPosition().x());
                            tree_track_RECOvsMINI_PropBarrel_firsthit_y.push_back(itTrack->innerPosition().Y()-PropTSOS.globalPosition().y());
                            tree_track_RECOvsMINI_PropBarrel_firsthit_z.push_back(itTrack->innerPosition().Z()-PropTSOS.globalPosition().z());
                                }

                            if(PropTSOStrueR.isValid())
                            {
                            //-----------------------With Propagator TSOStrueR-----------------------------//
                            tree_track_PropBarrel_firsthit_opti_x.push_back(PropTSOStrueR.globalPosition().x());
                            tree_track_PropBarrel_firsthit_opti_y.push_back(PropTSOStrueR.globalPosition().y());
                            tree_track_PropBarrel_firsthit_opti_z.push_back(PropTSOStrueR.globalPosition().z());
                            tree_track_RECOvsMINI_PropBarrel_firsthit_opti_x.push_back(itTrack->innerPosition().X()-PropTSOStrueR.globalPosition().x());
                            tree_track_RECOvsMINI_PropBarrel_firsthit_opti_y.push_back(itTrack->innerPosition().Y()-PropTSOStrueR.globalPosition().y());
                            tree_track_RECOvsMINI_PropBarrel_firsthit_opti_z.push_back(itTrack->innerPosition().Z()-PropTSOStrueR.globalPosition().z());
                            }

                        }
 //-------------------------------------------------disks/wheel-------------------------------------------------------//
                            //------------With geom aspects (mean Z)----------------------------//
                            zlayers = Players[i].second;
                            if (itTrack->pz()<0){zlayers=-zlayers;}
                            float deltaZ =itTrack->innerPosition().Z()-zlayers;
                            
                            if (abs(deltaZ)>30 ){zlayers=-zlayers;}//30 being a arbitrary value, could be mini~5, max ~55

                            float R = (zlayers-itTrack->vz())*tan(theta_prop);
                            float Rreal =  sqrt(itTrack->innerPosition().X()*itTrack->innerPosition().X()+itTrack->innerPosition().Y()+itTrack->innerPosition().Y());
                            float x0 = R*cos(itTrack->phi()); 
                            float y0 = R*sin(itTrack->phi() );  
                            float x0real = Rreal*cos(itTrack->phi());
                            float y0real = Rreal*sin(itTrack->phi());    
                            float zreal = (Rreal+itTrack->vz()*tan(theta_prop))/tan(theta_prop);

                            tree_track_GeoDisk_firsthit_x.push_back(x0);
                            tree_track_GeoDisk_firsthit_y.push_back(y0);
                            tree_track_GeoDisk_firsthit_z.push_back(zlayers);
                            tree_track_RECOvsMINI_GeoDisk_firsthit_x.push_back(itTrack->innerPosition().X()-x0);
                            tree_track_RECOvsMINI_GeoDisk_firsthit_y.push_back(itTrack->innerPosition().Y()-y0);
                            tree_track_RECOvsMINI_GeoDisk_firsthit_z.push_back(itTrack->innerPosition().Z()-zlayers);
                            
                            tree_track_GeoDisk_firsthit_opti_x.push_back(x0real);
                            tree_track_GeoDisk_firsthit_opti_y.push_back(y0real);
                            tree_track_GeoDisk_firsthit_opti_z.push_back(zreal);
                            tree_track_RECOvsMINI_GeoDisk_firsthit_opti_x.push_back(itTrack->innerPosition().X()-x0real);
                            tree_track_RECOvsMINI_GeoDisk_firsthit_opti_y.push_back(itTrack->innerPosition().Y()-y0real);
                            tree_track_RECOvsMINI_GeoDisk_firsthit_opti_z.push_back(itTrack->innerPosition().Z()-zreal);

                            //------------------WIth Propagators mean Z SLCC // true R
                            //------------------SLCC Mean Z------------------------------//
                            SimpleCylinderBounds SCB(0.,R,-zlayers,zlayers);
                            Cylinder CylindBounds(P3D,rot,SCB);
                            spair = SLCC.pathLength(CylindBounds);
                            
                            //     Versus the Optimal resolution for this method   //
                            //------------------SLCC true R-----------------------------//
                            float Ropti = sqrt(itTrack->innerPosition().X()*itTrack->innerPosition().X()+itTrack->innerPosition().Y()*itTrack->innerPosition().Y());
                            SimpleCylinderBounds SCB2(0.,Ropti,-abs(itTrack->innerPosition().Z()),abs(itTrack->innerPosition().Z()));
                            Cylinder CylindBounds2(P3D,rot,SCB2);
                            spair2 = SLCC.pathLength(CylindBounds2);
                            //---------------------------Positions-------------------------------//

                            double s  = spair.second;
                            double s2 = spair2.second;
                            LocalPoint GP = SLCC.position(s);
                            LocalPoint GPopti = SLCC.position(s2);

                            //------------------SLPC with mean Z from hitpattern----------------------://
                            GloballyPositioned<float>::PositionType P3D_(0.,0.,zlayers);
                            Plane P(P3D_,rot);
                            Basic3DVector<float> P3D2(itTrack->vx(),itTrack->vy(),itTrack->vz());//global frame
                            Basic3DVector<float> B3DV (itTrack->px(),itTrack->py(),itTrack->pz());//global frame 
                            StraightLinePlaneCrossing SLPC(P3D2,B3DV);
                            spairPlane = SLPC.position(P);
                            Basic3DVector<float> sPlane(-1000.,-1000.,-1000.);
                            if (spairPlane.first){ sPlane = spairPlane.second;}

                            //------------------SLPC with true Z from reco----------------------://
                            GloballyPositioned<float>::PositionType P3D_true(0.,0.,itTrack->innerPosition().Z());
                            Plane Ptrue(P3D_true,rot);
                            Basic3DVector<float> P3D2true(itTrack->vx(),itTrack->vy(),itTrack->vz());//global frame
                            Basic3DVector<float> B3DVtrue (itTrack->px(),itTrack->py(),itTrack->pz());//global frame => both vectors have to be given in the same frame
                            StraightLinePlaneCrossing SLPCtrue(P3D2true,B3DVtrue);
                            spairPlanetrue = SLPCtrue.position(Ptrue);
                            Basic3DVector<float> sPlanetrue (-1000.,-1000.,-1000.);
                            if (spairPlanetrue.first){sPlanetrue = spairPlanetrue.second;}


                                            //-------SLCC with miniaod data-------------//
                                            //The trasnformation in the "push_back" is used to go from local to global frame as GP is given in local coordinates
                                tree_track_PropDisk_SLCC_firsthit_x.push_back(2*(GP.x()+itTrack->vx()/2));
                                tree_track_PropDisk_SLCC_firsthit_y.push_back(2*(GP.y()+itTrack->vy()/2));
                                tree_track_PropDisk_SLCC_firsthit_z.push_back(2*(GP.z()+itTrack->vz()/2));
                                tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_x.push_back(itTrack->innerPosition().X()-2*GP.x()-itTrack->vx());//-itTrack->vx()
                                tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_y.push_back(itTrack->innerPosition().Y()-2*GP.y()-itTrack->vy());//-itTrack->vy()
                                tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_z.push_back(itTrack->innerPosition().Z()-2*GP.z()-itTrack->vz());//-itTrack->vz()
                                            //----SLCC with RECO data-----//
                                tree_track_PropDisk_SLCC_firsthit_opti_x.push_back(2*(GPopti.x()+itTrack->vx()/2));//Go to GLobal Coordiantes coordiantes... The transformation is not well understood atm
                                tree_track_PropDisk_SLCC_firsthit_opti_y.push_back(2*(GPopti.y()+itTrack->vy()/2));
                                tree_track_PropDisk_SLCC_firsthit_opti_z.push_back(2*(GPopti.z()+itTrack->vz()/2));
                                tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_x.push_back(itTrack->innerPosition().X()-2*GPopti.x()-itTrack->vx());
                                tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_y.push_back(itTrack->innerPosition().Y()-2*GPopti.y()-itTrack->vy());
                                tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_z.push_back(itTrack->innerPosition().Z()-2*GPopti.z()-itTrack->vz());

                                //-----------------SLPC--------------------//
                                tree_track_PropDisk_SLPC_firsthit_x.push_back(sPlane.x());
                                tree_track_PropDisk_SLPC_firsthit_y.push_back(sPlane.y());
                                tree_track_PropDisk_SLPC_firsthit_z.push_back(sPlane.z());
                                tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_x.push_back(itTrack->innerPosition().X()-sPlane.x());
                                tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_y.push_back(itTrack->innerPosition().Y()-sPlane.y());
                                tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_z.push_back(itTrack->innerPosition().Z()-sPlane.z());

                                tree_track_PropDisk_SLPC_firsthit_opti_x.push_back(sPlanetrue.x());
                                tree_track_PropDisk_SLPC_firsthit_opti_y.push_back(sPlanetrue.y());
                                tree_track_PropDisk_SLPC_firsthit_opti_z.push_back(sPlanetrue.z());
                                tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_x.push_back(itTrack->innerPosition().X()-sPlanetrue.x());
                                tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_y.push_back(itTrack->innerPosition().Y()-sPlanetrue.y());
                                tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_z.push_back(itTrack->innerPosition().Z()-sPlanetrue.z());
                            // }        
                }
            }
        //TBrowser example of lines to compare REOC and MINIAOD information, pt, NChi2 and drSig cut can also be added but don't chage much when the LLP flag is given:
        // ttree->Draw("tree_track_RECOvsMINI_GeoDisk_firsthit_y","tree_track_surf==1 && abs(tree_track_RECOvsMINI_GeoDisk_firsthit_y)<20 && tree_track_sim_LLP>0")
        //ttree->Draw("tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_x","tree_track_surf==1 && abs(tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_x)<20 && tree_track_sim_LLP>0")
        //ttree->Draw("tree_track_RECOvsMINI_GeoBarrel_firsthit_opti_z","tree_track_surf==0 && abs(tree_track_RECOvsMINI_GeoBarrel_firsthit_opti_z)<20 && tree_track_sim_LLP>0")



        count+=1;


      //----!!!!-----//

      //----------------------------
      //matching to simulated tracks
      //----------------------------
      //matching par hits
      
      if ( !runOnData_ ) 
      {
        int nSimHits = 0;
        bool isSimMatched = false;
        std::vector<int> tpIdx;
        std::vector<float> sharedFraction;
        std::vector<float> tpChi2;
        
        //initialized values for trackingParticle
        float sim_vx = -1000;
        float sim_vy = -1000;
        float sim_vz = -1000;
        
        int   sim_charge      =-1000;
        float sim_pt	      =-1;
        float sim_eta	      =-1000;
        float sim_phi	      =-1000;
        bool  sim_longLived   = false;
        //int	sim_matchedHit	= 0;
        int   sim_pdgId	      = -1000;
        int   sim_numberOfTrackerHits   = 0;
        int   sim_numberOfTrackerLayers = 0;
        float sim_mass		= 0.;
        int   sim_status	= -1000;

        int   sim_fromLLP               = -1;
	bool  matchB                    = false;
	bool  matchC                    = false;
        int   sim_isFromBC_mother_pdgId = -1000;
        float sim_dV	                = 1000.;
        float sim_dist	                = -10.;
        
        auto foundTPs = recSimColl.find(itTrack);
        if ( foundTPs != recSimColl.end() ) 
        {
          //if (!foundTPs->val.empty()) {
          isSimMatched = true;
          TrackingParticleRef tpr = foundTPs->val[0].first;
          nSimHits = tpr->numberOfTrackerHits();
          
          sim_charge	         = tpr->charge();
          sim_pt		 = tpr->pt();
          sim_eta  	         = tpr->eta();
          sim_phi  	         = tpr->phi();
          sim_longLived	         = tpr->longLived();
          // sim_matchedHit	    = tpr->matchedHit();
          sim_pdgId	         = tpr->pdgId();
          sim_numberOfTrackerHits   = tpr->numberOfTrackerHits();
          sim_numberOfTrackerLayers = tpr->numberOfTrackerLayers();
          sim_mass 	         = tpr->mass();
          sim_status	         = tpr->status();
          
          //determine x,y,z position of the genVertex which produced the associated simtrack
          sim_vx	 = tpr->vx();
          sim_vy	 = tpr->vy();
          sim_vz	 = tpr->vz();
          
          // std::cout<< "nllp : "<< nllp<<std::endl;
          float dV0 = (sim_vx - tree_GenPVx)*(sim_vx - tree_GenPVx)
                    + (sim_vy - tree_GenPVy)*(sim_vy - tree_GenPVy)
                    + (sim_vz - tree_GenPVz)*(sim_vz - tree_GenPVz);
	  float dV1 = 1000.;
          if ( nllp >= 1 ) {
            dV1 = (sim_vx - LLP1_x)*(sim_vx - LLP1_x)
                + (sim_vy - LLP1_y)*(sim_vy - LLP1_y)
                + (sim_vz - LLP1_z)*(sim_vz - LLP1_z);
            // std::cout<< "genVertexPosX : "<<sim_vx<<" & LLP1_X : "<< LLP1_x<< " & dV1 : "<< dV1<<std::endl;
          }
	  float dV2 = 1000.;
          if ( nllp >= 2 ) {
            dV2 = (sim_vx - LLP2_x)*(sim_vx - LLP2_x)
            	+ (sim_vy - LLP2_y)*(sim_vy - LLP2_y)
            	+ (sim_vz - LLP2_z)*(sim_vz - LLP2_z);
          }
          if      ( dV1 < dV2 && dV1 < 0.01 ) sim_fromLLP = 1;
          else if ( dV2 < dV1 && dV2 < 0.01 ) sim_fromLLP = 2;
          else if ( dV0 < 0.01 )	      sim_fromLLP = 0;

          sim_dist = TMath::Sqrt( (sim_vx-tree_GenPVx)*(sim_vx-tree_GenPVx) 
	                        + (sim_vy-tree_GenPVy)*(sim_vy-tree_GenPVy)
	                        + (sim_vz-tree_GenPVz)*(sim_vz-tree_GenPVz) );
          if ( dV1 < dV2 ) {
	    sim_dV = TMath::Sqrt(dV1);
	    if ( sim_dist < LLP1_dist ) sim_dV = -sim_dV;
	  }
          else {
	    sim_dV = TMath::Sqrt(dV2);
	    if ( sim_dist < LLP2_dist ) sim_dV = -sim_dV;
	  }
	
          // match simtrack to genparticle from B hadron
	  matchB = false;
          if ( tree_nFromB > 0 ) {
            for (int k=0; k<tree_nFromB; k++) { // loop on genFromB
	    if ( tree_genFromB_pdgId[k] != sim_pdgId ) continue;
	      float dpt  = sim_pt / tree_genFromB_pt[k] - 1.;
	      float deta = sim_eta - tree_genFromB_eta[k];
  	      float dphi = Deltaphi( sim_phi , tree_genFromB_phi[k] );
              if ( abs(deta) < 0.01 && abs(dphi) < 0.01 && abs(dpt) < 0.01 ) {
	        matchB = true;
		sim_isFromBC_mother_pdgId = tree_genFromB_mother_pdgId[k];
              }
	    } // end loop on genFromB
	  }
	 
          // match simtrack to genparticle from C hadron
	  matchC = false;
          if ( tree_nFromC > 0 ) {
            for (int k=0; k<tree_nFromC; k++) { // loop on genFromC
	    if ( tree_genFromC_pdgId[k] != sim_pdgId ) continue;
	      float dpt  = sim_pt / tree_genFromC_pt[k] - 1.;
	      float deta = sim_eta - tree_genFromC_eta[k];
  	      float dphi = Deltaphi( sim_phi , tree_genFromC_phi[k] );
              if ( abs(deta) < 0.01 && abs(dphi) < 0.01 && abs(dpt) < 0.01 ) {
	        matchC = true;
		if ( !matchB ) {
		  sim_isFromBC_mother_pdgId = tree_genFromC_mother_pdgId[k];
		}
              }
	    } // end loop on genFromC
          }
        }
        
        tree_track_nSimHits	             .push_back(nSimHits);
        tree_track_isSimMatched              .push_back(isSimMatched);
        
        tree_track_sim_pt  	             .push_back(sim_pt);
        tree_track_sim_eta 	             .push_back(sim_eta);
        tree_track_sim_phi 	             .push_back(sim_phi);
        tree_track_sim_charge	             .push_back(sim_charge);
        tree_track_sim_longLived             .push_back(sim_longLived);
        tree_track_sim_pdgId	             .push_back(sim_pdgId);
        tree_track_sim_numberOfTrackerHits   .push_back(sim_numberOfTrackerHits);
        tree_track_sim_numberOfTrackerLayers .push_back(sim_numberOfTrackerLayers);
        tree_track_sim_mass	             .push_back(sim_mass);
        tree_track_sim_status	             .push_back(sim_status);
        tree_track_sim_x	             .push_back(sim_vx);
        tree_track_sim_y	             .push_back(sim_vy);
        tree_track_sim_z	             .push_back(sim_vz);
        tree_track_sim_LLP	             .push_back(sim_fromLLP);
        tree_track_sim_isFromB	             .push_back(matchB);
        tree_track_sim_isFromC	             .push_back(matchC);
        tree_track_sim_isFromBC_mother_pdgId .push_back(sim_isFromBC_mother_pdgId);
        tree_track_sim_dV	             .push_back(sim_dV);
        tree_track_sim_dist                  .push_back(sim_dist);
      }
      
      /*if( !isSimMatched &&  trackToVextexMap[iTrack] != 0 ) nUnmatchTrack_fromPU++;
       if( !isSimMatched &&  trackToVextexMap[iTrack] == 0 ) nUnmatchTrack_fromPV++;
       if( !isSimMatched		    ) nUnmatchTrack++;
       if( trackToVextexMap[iTrack] != 0	    ) nPUTrack++;*/
      
      tree_track_recoVertex_idx.push_back(trackToVextexMap[iTrack]);
      tree_track_iJet.push_back(trackToAK4SlimmedJetMap[iTrack]);
      tree_track_recoAK4PFJet_idx.push_back(trackToAK4PFJetMap[iTrack]);
      tree_track_reco08Jet_idx.push_back(trackTo08JetMap[iTrack]);
      tree_track_recoCaloJet_idx.push_back(trackToCaloJetMap[iTrack]);
      
    } //end loop on tracks
    std::cout<<" tracks size : pcref sze : "<<trackRefs.size()<<"  // "<<pc->size()<<std::endl;
     std::cout<<" mini : 1/2/3 : "<<count3_mini<<std::endl;
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
      if ( imu1 >= 0 ) deltaR1 = Deltar( jet_eta, jet_phi, tree_muon_eta[imu1], tree_muon_phi[imu1] );
      if ( imu2 >= 0 ) deltaR2 = Deltar( jet_eta, jet_phi, tree_muon_eta[imu2], tree_muon_phi[imu2] );
      if ( deltaR1 < 0.4 || deltaR2 < 0.4 )
      {
        if ( deltaR1 < 0.4 )
        { //if muon is inside, we remove the muons infomation from the jet
          v1.SetPtEtaPhiM( tree_muon_pt[imu1],
          		  tree_muon_eta[imu1],
          		  tree_muon_phi[imu1],
          		  0 );
          v -= v1; //v TLorentzFactor being just above, defined by jet data
        }
        if ( deltaR2 < 0.4 )
        {
          v2.SetPtEtaPhiM( tree_muon_pt[imu2],
          		  tree_muon_eta[imu2],
          		  tree_muon_phi[imu2],
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
    float dRcut_hemis  = 1.5; // subjective choice
    float dRcut_tracks = 10.; // no cut is better (could bias low track pT and high LLP ct) 
     
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
    
//$$
//     // force the axes to the true LLP
//     vaxis1 = vneu[0];
//     vaxis2 = vneu[1];
//$$

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

    ///// MVA for track selection coming from displaced tops
    
    //ajoute par Paul /*!*/
    float drSig, isinjet;
    // float dptSig; /*!*/
    int jet; /*!*/
    float ntrk10, ntrk20, ntrk30; /*!*/
    float firsthit_X, firsthit_Y, firsthit_Z, dxy, dxyError, pt, eta,phi, NChi2, nhits;
    // float algo;
    // float  track_dR;
    // float track_dRmax ; /*!*/
    double bdtval = -100.;

    LLP1_nTrks = 0;
    LLP2_nTrks = 0;

    int nTrks_axis1 = 0, nTrks_axis1_sig=0, nTrks_axis1_bad=0;
    int nTrks_axis2 = 0, nTrks_axis2_sig=0, nTrks_axis2_bad=0;
    int nTrks_axis1_sig_mva=0, nTrks_axis1_bad_mva=0;
    int nTrks_axis2_sig_mva=0, nTrks_axis2_bad_mva=0;
    
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
    // reader->AddVariable( "mva_track_algo", &algo);
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
    double bdtcut = -0.1456; // optimal value w/o track association to axis: -0.0401
    //optimal value : -0.0815 for TMVAClassification_BDTG50sansalgo.weights.xml
    // optimal value 0.0057 for TMVAClassification_BDTG50cm.weights.xml
    // optimal value : -0.0401: TMVAbgctau50withnhits.xml
    // -0.1456 for TMVAClassification_BDTG50cm_HighPurity.weights.xml
    // 0.0327 for TMVAClassification_BDTG50cm_NewSignal.weights.xml
    // -0.1083 for TMVAClassification_BDTG_FromBC.weights.xml
//$$    double bdtcut = -10.; // if NO BDT cut !
//$$

    int counter_track = -1;
    
    for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack)  // Loop on all the tracks
    {
      counter_track++;
      const auto& itTrack = trackRefs[iTrack];
      firsthit_X = tree_track_firstHit_x[counter_track];
      firsthit_Y = tree_track_firstHit_y[counter_track];
      firsthit_Z = tree_track_firstHit_z[counter_track];
      dxy	 = tree_track_dxy[counter_track];
      dxyError   = tree_track_dxyError[counter_track];
      pt	 = tree_track_pt[counter_track];
      // ptError    = tree_track_ptError[counter_track];
      eta	 = tree_track_eta[counter_track];
      phi	 = tree_track_phi[counter_track];
      NChi2	 = tree_track_NChi2[counter_track];
      nhits	 = tree_track_nHit[counter_track];
      // algo	 = tree_track_algo[counter_track];
      // float dptSig=-1;
      // if (pt>0) dptSig=ptError/pt;
      
      //Ajoute par Paul /*!*/
      drSig = -1.;
      if ( dxyError > 0 ) drSig = abs(dxy) / dxyError; /*!*/
      tree_track_drSig.push_back(drSig);
      ntrk10 = 0;
      ntrk20 = 0;
      ntrk30 = 0;
      bdtval = -10.;
      dR = -1.;
      int tracks_axis = 0; // flag to check which axis is the closest from the track

//$$
      if ( pt > pt_Cut && NChi2 < NChi2_Cut && drSig > drSig_Cut ) // preselection : pt > 1. && NChi2 < 5. && drSig > 5.
//$$
      { // On regarde les track_selec[i] qui sont True donc de potentielles tracks secondaires
        
        jet = tree_track_iJet[counter_track]; /*!*/
        isinjet = 0.;
        if ( jet >= 0 ) isinjet = 1.; /*!*/
        int isFromLLP    = tree_track_sim_LLP[counter_track];

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
          if ( counter_othertrack == counter_track ) continue;
          float pt2  = tree_track_pt[counter_othertrack];
          float dr2 = abs(tree_track_dxy[counter_othertrack]);
          float drSig2 = -1.;
          if ( tree_track_dxyError[counter_othertrack] > 0 ) drSig2 = dr2 / tree_track_dxyError[counter_othertrack];
          NChi2 = tree_track_NChi2[counter_othertrack];
//$$
        if ( !(pt2 > pt_Cut && NChi2 < NChi2_Cut && drSig2 > drSig_Cut ) ) continue; // On regarde les autres track_selec[i] qui sont True donc de potnetielles tracks secondaires
//$$
          float x2 = tree_track_firstHit_x[counter_othertrack];
          float y2 = tree_track_firstHit_y[counter_othertrack];
          float z2 = tree_track_firstHit_z[counter_othertrack];
          float dist = TMath::Sqrt( (firsthit_X-x2)*(firsthit_X-x2) + (firsthit_Y-y2)*(firsthit_Y-y2) + (firsthit_Z-z2)*(firsthit_Z-z2) );//pour chaque reconstruite, on regarde les autres tracks,
          if ( dist < 10. )	     {ntrk10++;} // les sctocker les 3 , on teste sur une seule couche quand on regarde vers l'avant
          if ( dist < 20. )	     {ntrk20++;}
          if ( dist < 30. )	     {ntrk30++;}
        }  // end Loop on other Tracks

        if ( dR < dRcut_tracks ) 
	{
	  if ( isFromLLP == 1 ) LLP1_nTrks++;
	  if ( isFromLLP == 2 ) LLP2_nTrks++;
	
          bdtval = reader->EvaluateMVA( "BDTG" ); //default value = -10 (no -10 observed and -999 comes from EvaluateMVA)
          //cout << "BDT VAL " << bdtval <<endl;

          if ( tracks_axis == 1 ) {
	    nTrks_axis1++;
            if ( isFromLLP == iLLPrec1 ) nTrks_axis1_sig++;
            else if ( isFromLLP >= 1 )   nTrks_axis1_bad++;
          }
        
          if ( tracks_axis == 2 ) {
	    nTrks_axis2++;
            if ( isFromLLP == iLLPrec2 ) nTrks_axis2_sig++;
            else if ( isFromLLP >= 1 )   nTrks_axis2_bad++;
          }
        
          if ( bdtval > bdtcut )      //optimal cut for 50 cm ,trained on 10k events //Paul, expected to change
          { 
            ////--------------Control tracks-----------------////
            if ( isFromLLP == 1 )
            {
              displacedTracks_llp1_mva.push_back(theTransientTrackBuilder->build(*itTrack));
            }
            if ( isFromLLP == 2 )
            {
              displacedTracks_llp2_mva.push_back(theTransientTrackBuilder->build(*itTrack));
            }
          
            if ( tracks_axis == 1 )
            {
              displacedTracks_Hemi1_mva.push_back(theTransientTrackBuilder->build(*itTrack));
              if ( isFromLLP == iLLPrec1 ) nTrks_axis1_sig_mva++;
              else if ( isFromLLP >= 1 )   nTrks_axis1_bad_mva++;
            }
 
            if ( tracks_axis == 2 )
            {
              displacedTracks_Hemi2_mva.push_back(theTransientTrackBuilder->build(*itTrack));
              if ( isFromLLP == iLLPrec2 ) nTrks_axis2_sig_mva++;
              else if ( isFromLLP >= 1 )   nTrks_axis2_bad_mva++;
            }
          }
        }
      }
      
      tree_track_ntrk10.push_back(ntrk10);
      tree_track_ntrk20.push_back(ntrk20);
      tree_track_ntrk30.push_back(ntrk30);
      tree_track_MVAval.push_back(bdtval);
      tree_track_Hemi.push_back(tracks_axis);
      tree_track_Hemi_dR.push_back(dR);
      if      ( tracks_axis == 1 ) tree_track_Hemi_LLP.push_back(iLLPrec1);
      else if ( tracks_axis == 2 ) tree_track_Hemi_LLP.push_back(iLLPrec2);
      else		           tree_track_Hemi_LLP.push_back(0);
      
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
//$$$$  } // endif ZMu and HT selections
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
    
    tree_vtx_PosX.clear();
    tree_vtx_PosY.clear();
    tree_vtx_PosZ.clear();
    tree_vtx_NChi2.clear();
    tree_vtx_PosXError.clear();
    tree_vtx_PosYError.clear();
    tree_vtx_PosZError.clear();
    tree_vtx_ndf.clear();
    
    tree_jet_E.clear();
    tree_jet_pt.clear();
    tree_jet_eta.clear();
    tree_jet_phi.clear();
//$$    tree_jet_idxTrack.clear();
    
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
    
    tree_electron_pt.clear();
    tree_electron_eta.clear();
    tree_electron_phi.clear();
    tree_electron_x.clear();
    tree_electron_y.clear();
    tree_electron_z.clear();
    tree_electron_energy.clear();
    tree_electron_charge.clear();
    
    tree_muon_pt.clear();
    tree_muon_eta.clear();
    tree_muon_phi.clear();
    tree_muon_x.clear();
    tree_muon_y.clear();
    tree_muon_z.clear();
    tree_muon_energy.clear();
    tree_muon_dxy.clear();
    tree_muon_dxyError.clear();
    tree_muon_dz.clear();
    tree_muon_dzError.clear();
    tree_muon_charge.clear();
    tree_muon_isLoose.clear();
    tree_muon_isTight.clear();
    tree_muon_isGlobal.clear();
    tree_muon_PFisoTight.clear();
    tree_muon_PFisoVeryTight.clear();
    tree_muon_PFisoTight.clear();
    tree_muon_PFisoMedium.clear();
    tree_muon_PFisoLoose.clear();
    tree_muon_MVAisoLoose.clear();
    tree_muon_MVAisoMedium.clear();
    tree_muon_MVAisoTight.clear();
    tree_muon_isStandAloneMuon.clear();
    tree_muon_CutBasedIdMedium.clear();
    tree_muon_CutBasedIdMediumPrompt.clear();
    
    //-----------------------
    //fill the tree per track
    //tree_track_nclusters.clear();
    tree_track_pz.clear();
    tree_track_pt.clear();
    tree_track_ptError.clear();
    tree_track_eta.clear();
    tree_track_phi.clear();
    tree_track_charge.clear();
    tree_track_NChi2.clear();
    tree_track_isHighPurity.clear();
    tree_track_dxy.clear();
    tree_track_dxyError.clear();
    tree_track_dz.clear();
    tree_track_dzError.clear();
    tree_track_nHit.clear();
    tree_track_nHitPixel.clear();
    tree_track_nHitTIB.clear();
    tree_track_nHitTID.clear();
    tree_track_nHitTOB.clear();
    tree_track_nHitTEC.clear();
    tree_track_nHitPXB.clear();
    tree_track_nHitPXF.clear();
    tree_track_isHitPixel.clear();
    tree_track_nLayers.clear();
    tree_track_nLayersPixel.clear();
    tree_track_x.clear();
    tree_track_y.clear();
    tree_track_z.clear();
    tree_track_firstHit.clear();
    tree_track_firstHit_x.clear();
    tree_track_firstHit_y.clear();
    tree_track_firstHit_z.clear();
    tree_track_firstHit_phi.clear();

    //
    tree_track_ddxy.clear();
    tree_track_ddxyError.clear();
                        tree_track_dpt.clear();
          tree_track_deta.clear();
                    tree_track_dphi.clear();
          tree_track_dNChi2.clear();
                    tree_track_dvx.clear();
          tree_track_dvy.clear();
                    tree_track_dvz.clear();
          tree_track_ddz.clear();
                    tree_track_ddzError.clear();
    //
//$$$$
    tree_track_extraTrue_dx.clear();
    tree_track_extraTrue_dy.clear();
    tree_track_extraTrue_dz.clear();
    tree_track_extraTrue_dphi.clear();
    tree_track_fromLayer_dx.clear();
    tree_track_fromLayer_dy.clear();
    tree_track_fromLayer_dz.clear();
    tree_track_fromLayer_dphi.clear();
//$$$$
    tree_track_iJet.clear();
    tree_track_ntrk10.clear();
    tree_track_ntrk20.clear();
    tree_track_ntrk30.clear();
    tree_track_MVAval.clear();
    tree_track_Hemi.clear();
    tree_track_Hemi_dR.clear();
    tree_track_Hemi_mva_NChi2.clear();
    tree_track_Hemi_LLP.clear();

    tree_track_outerPt.clear();
    tree_track_isLoose.clear();
    tree_track_isTight.clear();
    tree_track_numberOfLostHits.clear();
    tree_track_originalAlgo.clear();
    tree_track_algo.clear();
    tree_track_stopReason.clear();

    tree_track_stripTECLayersWithMeasurement .clear();
    tree_track_stripTIBLayersWithMeasurement.clear();
    tree_track_stripTIDLayersWithMeasurement.clear();
    tree_track_stripTOBLayersWithMeasurement.clear();
    tree_track_recoVertex_idx.clear();
    tree_track_recoAK4PFJet_idx.clear();
    tree_track_reco08Jet_idx.clear();
    tree_track_recoCaloJet_idx.clear();

    tree_track_theta.clear();
    tree_track_surf.clear();
    tree_track_hitpattern.clear();

    tree_track_GeoBarrel_firsthit_x.clear();
    tree_track_GeoBarrel_firsthit_y.clear();
    tree_track_GeoBarrel_firsthit_z.clear();
    tree_track_RECOvsMINI_GeoBarrel_firsthit_x.clear();
    tree_track_RECOvsMINI_GeoBarrel_firsthit_y.clear();
    tree_track_RECOvsMINI_GeoBarrel_firsthit_z.clear();

    tree_track_PropBarrel_firsthit_x.clear();
    tree_track_PropBarrel_firsthit_y.clear();
    tree_track_PropBarrel_firsthit_z.clear();
    tree_track_RECOvsMINI_PropBarrel_firsthit_x.clear();
    tree_track_RECOvsMINI_PropBarrel_firsthit_y.clear();
    tree_track_RECOvsMINI_PropBarrel_firsthit_z.clear();

    tree_track_PropBarrel_firsthit_opti_x.clear();
    tree_track_PropBarrel_firsthit_opti_y.clear();
    tree_track_PropBarrel_firsthit_opti_z.clear();
    tree_track_RECOvsMINI_PropBarrel_firsthit_opti_x.clear();
    tree_track_RECOvsMINI_PropBarrel_firsthit_opti_y.clear();
    tree_track_RECOvsMINI_PropBarrel_firsthit_opti_z.clear();

    tree_track_GeoDisk_firsthit_x.clear();
    tree_track_GeoDisk_firsthit_y.clear();
    tree_track_GeoDisk_firsthit_z.clear();
    tree_track_RECOvsMINI_GeoDisk_firsthit_x.clear();
    tree_track_RECOvsMINI_GeoDisk_firsthit_y.clear();
    tree_track_RECOvsMINI_GeoDisk_firsthit_z.clear();
      
    tree_track_GeoDisk_firsthit_opti_x.clear();
    tree_track_GeoDisk_firsthit_opti_y.clear();
    tree_track_GeoDisk_firsthit_opti_z.clear();
    tree_track_RECOvsMINI_GeoDisk_firsthit_opti_x.clear();
    tree_track_RECOvsMINI_GeoDisk_firsthit_opti_y.clear();
    tree_track_RECOvsMINI_GeoDisk_firsthit_opti_z.clear();

    tree_track_PropDisk_SLCC_firsthit_x.clear();
    tree_track_PropDisk_SLCC_firsthit_y.clear();
    tree_track_PropDisk_SLCC_firsthit_z.clear();
    tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_x.clear();
    tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_y.clear();
    tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_z.clear();

    tree_track_PropDisk_SLCC_firsthit_opti_x.clear();
    tree_track_PropDisk_SLCC_firsthit_opti_y.clear();
    tree_track_PropDisk_SLCC_firsthit_opti_z.clear();
    tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_x.clear();
    tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_y.clear();
    tree_track_RECOvsMINI_PropDisk_SLCC_firsthit_opti_z.clear();

          
    tree_track_PropDisk_SLPC_firsthit_x.clear();
    tree_track_PropDisk_SLPC_firsthit_y.clear();
    tree_track_PropDisk_SLPC_firsthit_z.clear();
    tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_x.clear();
    tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_y.clear();
    tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_z.clear();

    tree_track_PropDisk_SLPC_firsthit_opti_x.clear();
    tree_track_PropDisk_SLPC_firsthit_opti_y.clear();
    tree_track_PropDisk_SLPC_firsthit_opti_z.clear();
    tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_x.clear();
    tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_y.clear();
    tree_track_RECOvsMINI_PropDisk_SLPC_firsthit_opti_z.clear();

    tree_track_sim_LLP.clear();
    tree_track_sim_isFromB.clear();
    tree_track_sim_isFromC.clear();
    tree_track_sim_pt.clear();
    tree_track_sim_eta.clear();
    tree_track_sim_phi.clear();
    tree_track_sim_charge.clear();
    tree_track_sim_pdgId.clear();
    tree_track_sim_mass.clear();
    tree_track_sim_x.clear();
    tree_track_sim_y.clear();
    tree_track_sim_z.clear();

    tree_track_sim_status.clear();
    tree_track_sim_longLived.clear();
    tree_track_isSimMatched.clear();
//     tree_track_sim_matchedHit .clear();
    tree_track_nSimHits.clear();
    tree_track_sim_numberOfTrackerHits.clear();
    tree_track_sim_numberOfTrackerLayers.clear();
    tree_track_sim_isFromBC_mother_pdgId.clear();
    tree_track_sim_dV.clear();
    tree_track_sim_dist.clear();
    
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
    
    tree_genParticle_pt.clear();
    tree_genParticle_eta.clear();
    tree_genParticle_phi.clear();
    tree_genParticle_charge.clear();
    tree_genParticle_pdgId.clear();
    tree_genParticle_x.clear();
    tree_genParticle_y.clear();
    tree_genParticle_z.clear();
    tree_genParticle_mass.clear();
    tree_genParticle_statusCode.clear();
    tree_genParticle_mother_pdgId.clear();
    tree_genParticle_LLP.clear();

    tree_genFromC_pt.clear();
    tree_genFromC_eta.clear();
    tree_genFromC_phi.clear();
    tree_genFromC_charge.clear();
    tree_genFromC_pdgId.clear();
    tree_genFromC_x.clear();
    tree_genFromC_y.clear();
    tree_genFromC_z.clear();
    tree_genFromC_mother_pdgId.clear();
    tree_genFromC_generation.clear();
    tree_genFromC_LLP.clear();

    tree_genFromB_pt.clear();
    tree_genFromB_eta.clear();
    tree_genFromB_phi.clear();
    tree_genFromB_charge.clear();
    tree_genFromB_pdgId.clear();
    tree_genFromB_x.clear();
    tree_genFromB_y.clear();
    tree_genFromB_z.clear();
    tree_genFromB_mother_pdgId.clear();
    tree_genFromB_generation.clear();
    tree_genFromB_LLP.clear();

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

    tree_track_drSig.clear();
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



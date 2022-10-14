# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase1_2018_realistic -n 10 --era Run2_2018 --eventcontent RECOSIM,MINIAODSIM,DQM --runUnscheduled -s RAW2DIGI,L1Reco,RECO,RECOSIM,EI,PAT,VALIDATION:@standardValidation+@miniAODValidation,DQM:@standardDQM+@ExtraHLT+@miniAODDQM --datatier GEN-SIM-RECO,MINIAODSIM,DQMIO --geometry DB:Extended --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

#$$
# for new tracking
from Configuration.ProcessModifiers.displacedTracking_cff import displacedTracking 
#$$

process = cms.Process('RECO',Run2_2018)
# process = cms.Process('MINIAODSIM',Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
#$$
    # input = cms.untracked.int32(200)
   input = cms.untracked.int32(50)
#$$
)

# Input source
process.source = cms.Source("PoolSource",
#$$
    fileNames = cms.untracked.vstring(
     'file:step2HLT_1.root'
    #  'file:MINIAODSIM_v16_L1v1.root'
    # 'file:/opt/sbg/cms/ui2_data1/blochd/CMSSW_10_2_16_UL/src/step2HLT/step2HLT_10evts.root'
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_1.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_10.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_11.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_12.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_13.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_14.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_15.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_16.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_17.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_18.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_19.root",
    ##"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_2.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_20.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_21.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_22.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_23.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_24.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_25.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_26.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_27.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_28.root",
   #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_29.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_3.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_30.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_31.root",
   #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_32.root",
   # "/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_33.root",
   # "/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_34.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_35.root",
   # "/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_36.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_37.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_38.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_39.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_4.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_40.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_41.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_42.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_43.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_44.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_45.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_46.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_47.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_48.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_49.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_5.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_50.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_6.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_7.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_8.root",
    #"/store/user/blochd/CMSSW_10_6_20_LLP/MC/UDD_bgctau50_smu275_snu225/2018_step2HLT/220124_090507/0000/step2HLT_9.root"
    ), 
#$$
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(      
    # SkipEvent = cms.untracked.vstring('ProductNotFound'),

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:RECO_v16_L1v1.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('file:MINIAODSIM_v16_L1v1.root'),
    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    overrideBranchesSplitLevel = cms.untracked.VPSet(
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedCandidates_packedPFCandidates__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenParticles_prunedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patTriggerObjectStandAlones_slimmedPatTrigger__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patPackedGenParticles_packedGenParticles__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoVertexs_offlineSlimmedPrimaryVertices__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoCaloClusters_reducedEgamma_reducedESClusters_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEBRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedEERecHits_*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('recoGenJets_slimmedGenJets__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('patJets_slimmedJetsPuppi__*'),
            splitLevel = cms.untracked.int32(99)
        ), 
        cms.untracked.PSet(
            branch = cms.untracked.string('EcalRecHitsSorted_reducedEgamma_reducedESRecHits_*'),
            splitLevel = cms.untracked.int32(99)
        )
    ),
    overrideInputFileSplitLevels = cms.untracked.bool(True),
    splitLevel = cms.untracked.int32(0)
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:DQM_v16_L1v1.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")

from Configuration.AlCa.GlobalTag import GlobalTag
#$$
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')
#$$

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_hfNoisyHitsFilter = cms.Path(process.hfNoisyHitsFilter)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_ecalBadCalibFilter = cms.Path(process.ecalBadCalibFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.Flag_BadPFMuonDzFilter = cms.Path(process.BadPFMuonDzFilter)
process.prevalidation_step = cms.Path(process.prevalidation)
process.prevalidation_step1 = cms.Path(process.prevalidationMiniAOD)
process.validation_step = cms.EndPath(process.validation)
process.validation_step1 = cms.EndPath(process.validationMiniAOD)
process.dqmoffline_step = cms.EndPath(process.DQMOffline)
process.dqmoffline_1_step = cms.EndPath(process.DQMOfflineExtraHLT)
process.dqmoffline_2_step = cms.EndPath(process.DQMOfflineMiniAOD)
process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
process.dqmofflineOnPAT_1_step = cms.EndPath(process.PostDQMOffline)
process.dqmofflineOnPAT_2_step = cms.EndPath(process.PostDQMOfflineMiniAOD)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)
process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

from CommonTools.RecoAlgos.trackingParticleRefSelector_cfi import trackingParticleRefSelector as _trackingParticleRefSelector
process.trackingParticlesIntime = _trackingParticleRefSelector.clone(
    signalOnly = False,
    intimeOnly = False,
    chargedOnly = True,
    tip = 1e5,
    lip = 1e5,
    minRapidity = -2.4,
    maxRapidity =  2.4,
    ptMin = 1,
)

process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimGeneral.MixingModule.trackingTruthProducerSelection_cfi")

process.trackingParticles.simHitCollections = cms.PSet( )

############## CHANGED THIS BECAUSE ERROR MESSAGE ABOUT MIXING  
#process.mix.playback = cms.untracked.bool(False)
#process.mix.digitizers = cms.PSet(
#     mergedtruth = cms.PSet(process.trackingParticles)
#)
for a in process.aliases: delattr(process, a)

process.trackingPerf = cms.EDAnalyzer('TrackingPerf',     
#$$
weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50cm.weights.xml" ),
#weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG.weights.xml" ),
    # weightFileMVA = cms.untracked.string( "TMVAbgctau50withnhits.xml" ),
    # weightFileMVA = cms.untracked.string( "TMVAClassification_BDT_10cm_10k.weights.xml" ),
#$$

      tracks                   = cms.untracked.InputTag('packedPFCandidates'),#'generalTracks'ou 'packedPFCandidates'
      trackLabel               = cms.InputTag('packedPFCandidates'),#not used,'generalTracks'ou 'packedPFCandidates'
      trackingParticles        = cms.untracked.InputTag('trackingParticlesIntime'),
      trackingParticlesRef     = cms.untracked.bool(True),
    #   trackAssociator         = cms.untracked.InputTag('trackingParticleRecoTrackAsssociation'),
      trackAssociator          = cms.untracked.InputTag('quickTrackAssociatorByHits'),
      pileup                   = cms.InputTag('slimmedAddPileupInfo'),
      beamSpot                 = cms.untracked.InputTag('offlineBeamSpot'),
      vertices                 = cms.untracked.InputTag('offlineSlimmedPrimaryVertices'),
      ak4slimmedJetInput       = cms.InputTag('slimmedJets'),
    #   ak4PFJetInput            = cms.InputTag('ak4PFJets'), #collection not in miniaod
    #   caloJetInput             = cms.InputTag('ak4CaloJets'),#collection not in miniaod
      ak8jetInput              = cms.InputTag('slimmedJetsAK8'),
    #  ak8CaloJetInput          = cms.InputTag('ak8CaloJets'),#collection not in miniaod
#      pfmetInput               = cms.InputTag('pfMet'),
      genParticles             = cms.InputTag('genParticles'),
      genJetInput              = cms.InputTag("slimmedGenJets"),
      ak8GenJetInput           = cms.InputTag("ak8GenJetsNoNu"),
      genEventInfoInput        = cms.InputTag("generator"),#not used
      LHEEventProductInput     = cms.InputTag("externalLHEProducer"),#not used
      pfcands                  = cms.InputTag("packedPFCandidates"),
      parametersDefiner        = cms.untracked.string('LhcParametersDefinerForTP'),
      electronInput            = cms.InputTag("slimmedElectrons"),
      #not me electronInput            = cms.InputTag("patElectronsPFlow"),
      slimmedmuonInput         = cms.InputTag("slimmedMuons"),
# notme      recomuonInput            = cms.InputTag("muons"), 
      metInput                 = cms.InputTag("slimmedMETs"),
      TTRHBuilder              = cms.string('WithTrackAngle'),
      filterTriggerNames       = cms.untracked.vstring("HLT_Mu*", "HLT_IsoMu*"),
      Zmumucand                = cms.InputTag('dimuonsZMuSkim'),
      useCluster               = cms.untracked.bool(False),
      runOnData                = cms.untracked.bool(False),
      KVFParameters = cms.PSet(
        maxDistance = cms.double(1),
        maxNbrOfIterations = cms.int32(100)
      )
    #   ,
    #   pruned = cms.InputTag("prunedGenParticles")#byme
    #   packed = cms.InputTag("packedGenParticles")#byme 
      
    #   GSFParameters = cms.PSet( #AVF
    #     maxshift = cms.double(0.0001), #Max transverse distance bertween vertex i and i-1
    #     maxstep = cms.int32(30), # Criterion for the relineariation of the tracks
    #     maxlpshift = cms.double(0.1), # Max number of iterations to perform
    #     weightthreshold = cms.double(0.001) # Min track weight for a track to be considered significant
    #   )
    #   GSFParameters = cms.PSet( #TrimmedVertexFitter
    #     minPt = cms.double(1.), #pT cut to apply to the tracks
    #     trackCompatibilityCut = cms.double(0.05), # proba below which a track is considered incompatible with the 1st vertex candidate formed
    #     vtxFitProbCut = cms.double(0.01) # Proba below which a vertex is rejected
    #   )
     
)
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideTrimmedVertexFitter
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAdaptiveVertexFitter
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideKinematicVertexFit

#Pset(1,100)default
process.trackingPerf_step = cms.EndPath(process.trackingPerf)

# process.load("RecoJets.JetProducers.ak8CaloJets_cfi")

process.load('DPGAnalysis.Skims.ZMuSkim_cff')

process.TrackingPerfNtuple = cms.Path( 
#    process.LHCTransport*
#     process.g4SimHits*
    process.mix *
#    process.simHitTPAssocProducer*
    process.tpClusterProducer*
    process.quickTrackAssociatorByHits*
    process.trackingParticleRecoTrackAsssociation*
    process.trackingParticlesIntime*
#    process.ZMuHLTFilter *
    process.looseMuonsForZMuSkim *
    process.ConcretelooseMuonsForZMuSkim *
    process.tkIsoDepositTk *
    process.allPatTracks *
    process.looseIsoMuonsForZMuSkim * 
    process.tightMuonsCandidateForZMuSkim *
    process.tightMuonsForZMuSkim *
    process.dimuonsZMuSkim *
#    process.ak8CaloJets*
    process.trackingPerf
)

########## output of ntuple

#$$
process.TFileService = cms.Service("TFileService", fileName = cms.string("Ntuple_v16_L1v1.root") )
#$$

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,
process.eventinterpretaion_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,
process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,
process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,
process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_ecalBadCalibFilter,
process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,
process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,
process.Flag_BadPFMuonFilter,process.Flag_BadPFMuonDzFilter,process.Flag_hfNoisyHitsFilter,process.Flag_BadChargedCandidateSummer16Filter,
process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,
process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,
process.prevalidation_step,process.prevalidation_step1,process.validation_step,process.validation_step1,
process.dqmoffline_step,process.dqmoffline_1_step,process.dqmoffline_2_step,
process.dqmofflineOnPAT_step,process.dqmofflineOnPAT_1_step,process.dqmofflineOnPAT_2_step,
#$$ 
# process.RECOSIMoutput_step,
# process.MINIAODSIMoutput_step,
# process.DQMoutput_step,
process.TrackingPerfNtuple)
#$$ 

process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
process = setCrossingFrameOn(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

#$$ 
process.options.numberOfThreads=cms.untracked.uint32(4)#16 when using ui3, 4 when using ui2 and crab
#$$ 

# End of customisation functions

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


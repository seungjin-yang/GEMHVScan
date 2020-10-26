# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: runAnalyzerMC --mc --era Run3 --conditions auto:phase1_2021_cosmics --geometry DB:Extended --step NONE --filein file:/hdfs/store/user/seyang/store/data/Commissioning2020/Cosmics/RAW-RECO/CosmicSP-PromptReco-v1/000/337/234/00000/7D93F41F-50A7-474B-9ED7-A184BAEDC2D7.root --fileout file:output.root -n -1 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('NONE',Run3)

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/hdfs/store/user/seyang/store/data/Commissioning2020/Cosmics/RAW-RECO/CosmicSP-PromptReco-v1/000/337/234/00000/7D93F41F-50A7-474B-9ED7-A184BAEDC2D7.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
)

# Additional output definition
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('file:out.root')
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_cosmics', '')

# Path and EndPath definitions
process.GEMCosmicAnalyzer = cms.EDAnalyzer('GEMCosmicAnalyzer', 
    process.MuonServiceProxy, 
    muonTag = cms.InputTag('muons'),
    recHitTag = cms.InputTag('gemRecHits'),
    segmentTag = cms.InputTag('gemSegments'),
)

process.p = cms.Path(process.GEMCosmicAnalyzer)

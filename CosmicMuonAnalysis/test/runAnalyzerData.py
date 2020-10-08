# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: runAnalyzerData --data --era Run3 --conditions auto:run3_data_promptlike --geometry DB:Extended --step NONE --filein file:/hdfs/store/user/seyang/store/data/Commissioning2020/Cosmics/RAW-RECO/CosmicSP-PromptReco-v1/000/337/234/00000/7D93F41F-50A7-474B-9ED7-A184BAEDC2D7.root --fileout file:output.root -n -1 --no_exec
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
from FWCore.ParameterSet import VarParsing
import glob
import os.path

options = VarParsing.VarParsing ('analysis')
options.register(
    "inPath",
    "/hdfs/store/user/seyang/store/data/Commissioning2020/Cosmics/RAW-RECO/CosmicSP-PromptReco-v1/000/337/234/00000/7D93F41F-50A7-474B-9ED7-A184BAEDC2D7.root",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,
    "a path to a RECO file or a directory having such files")
options.parseArguments()

fileNames = cms.untracked.vstring('file:' + options.inPath)


process.source = cms.Source("PoolSource",
    fileNames = fileNames,
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('runAnalyzerData nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Additional output definition
outFileName = cms.string('file:' + '_'.join(options.inPath.split('/')[7:]))
process.TFileService = cms.Service("TFileService",
    fileName = outFileName)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data_promptlike', '')

# Path and EndPath definitions
process.GEMCosmicAnalyzer = cms.EDAnalyzer('GEMCosmicAnalyzer', 
    process.MuonServiceProxy, 
    muonTag = cms.InputTag('muons'),
    recHitTag = cms.InputTag('gemRecHits'),
    segmentTag = cms.InputTag('gemSegments'),
)


process.p = cms.Path(process.GEMCosmicAnalyzer)

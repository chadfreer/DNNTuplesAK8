import FWCore.ParameterSet.Config as cms

# ---------------------------------------------------------

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.outputFile = 'output.root'
options.inputFiles = '/store/mc/RunIISummer16MiniAODv2/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/48E7598C-A7E6-E611-8092-002590DE6E76.root'
options.maxEvents = 500

options.register('inputScript', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, "input Script")
options.register('skipEvents', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "skip N events")
options.register('job', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "job number")
options.register('nJobs', 1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "total jobs")
options.register('fjKeepFlavors', [], VarParsing.multiplicity.list, VarParsing.varType.int, "Types of fatjet to keep in this sample")
# options.register('gluonReduction', 0.0, VarParsing.multiplicity.singleton, VarParsing.varType.float, "gluon reduction")
options.register('inputDataset',
                 '',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input dataset")

options.setupTags(tag='%d', ifCond='nJobs > 1', tagArg='job')
options.parseArguments()

# ---------------------------------------------------------

process = cms.Process("DNNFiller")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
if not options.inputScript:  # this is probably for testing
	process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.options = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True),  
   wantSummary=cms.untracked.bool(False)
)

print ('Using output file ' + options.outputFile)

process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

process.source = cms.Source('PoolSource',
    fileNames=cms.untracked.vstring(options.inputFiles),
    skipEvents=cms.untracked.uint32(options.skipEvents)
)

if options.inputScript:
    process.load(options.inputScript)

numberOfFiles = len(process.source.fileNames)
numberOfJobs = options.nJobs
jobNumber = options.job

process.source.fileNames = process.source.fileNames[jobNumber:numberOfFiles:numberOfJobs]
if options.nJobs > 1:
    print ("running over these files:")
    print (process.source.fileNames)

# ---------------------------------------------------------

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# ---------------------------------------------------------

bTagInfos = [
	'pfImpactParameterTagInfos',
	'pfInclusiveSecondaryVertexFinderTagInfos',
	'pfDeepCSVTagInfos',
	'pfBoostedDoubleSVAK8TagInfos'
]

bTagDiscriminators = [
	'softPFMuonBJetTags',
	'softPFElectronBJetTags',
	'pfJetBProbabilityBJetTags',
	'pfJetProbabilityBJetTags',
	'pfCombinedInclusiveSecondaryVertexV2BJetTags',
	'pfDeepCSVJetTags:probudsg',
	'pfDeepCSVJetTags:probb',
	'pfDeepCSVJetTags:probc',
	'pfDeepCSVJetTags:probbb',
	'pfDeepCSVJetTags:probcc',
	'pfBoostedDoubleSecondaryVertexAK8BJetTags'
]

jetCorrectionsAK8 = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
        process,
        labelName = "DeepFlavour",
#        jetSource = cms.InputTag('slimmedJetsAK8'),
        jetSource = cms.InputTag('slimmedJets'),
        jetCorrections = jetCorrectionsAK4,
        pfCandidates = cms.InputTag('packedPFCandidates'),
        pvSource = cms.InputTag("offlineSlimmedPrimaryVertices"),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        muSource = cms.InputTag('slimmedMuons'),
        elSource = cms.InputTag('slimmedElectrons'),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        explicitJTA = False
)

if hasattr(process,'updatedPatJetsTransientCorrectedDeepFlavour'):
	process.updatedPatJetsTransientCorrectedDeepFlavour.addTagInfos = cms.bool(True) 
	process.updatedPatJetsTransientCorrectedDeepFlavour.addBTagInfo = cms.bool(True)
else:
	raise ValueError('I could not find updatedPatJetsTransientCorrectedDeepFlavour to embed the tagInfos, please check the cfg')
# ---------------------------------------------------------

# DeepNtuplizer
process.load("DeepNTuples.NtupleAK8.DeepNtuplizerAK8_cfi")
#process.deepntuplizer.jets = cms.InputTag('selectedUpdatedPatJetsDeepFlavour')
process.deepntuplizer.jets = cms.InputTag("slimmedJets")
process.deepntuplizer.bDiscriminators = bTagDiscriminators 
#process.deepntuplizer.bDiscriminators.append('pfCombinedMVAV2BJetTags')
process.deepntuplizer.LooseSVs = cms.InputTag("looseIVFinclusiveCandidateSecondaryVertices")
process.deepntuplizer.fjKeepFlavors = cms.untracked.vuint32(options.fjKeepFlavors)
process.deepntuplizer.isQCDSample = '/QCD_' in options.inputDataset

process.deepntuplizer.jetR = cms.double(0.4)

# process.p = cms.Path(process.QGTagger + process.genJetSequence * process.deepntuplizer)
process.p = cms.Path(process.deepntuplizer)

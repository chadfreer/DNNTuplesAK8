import FWCore.ParameterSet.Config as cms

deepntuplizer = cms.EDAnalyzer('DeepNtuplizerAK8',
                                vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                puInfo     = cms.InputTag("slimmedAddPileupInfo"),
                                rhoInfo    = cms.InputTag("fixedGridRhoFastjetAll"),    
                                jets       = cms.InputTag("slimmedJets"),
                                jetR       = cms.double(0.4),
                                SVs        = cms.InputTag("slimmedSecondaryVertices"),
                                LooseSVs   = cms.InputTag("inclusiveCandidateSecondaryVertices"),
#                                 genJetMatchWithNu = cms.InputTag("patGenJetMatchWithNu"),
#                                 genJetMatchRecluster = cms.InputTag("patGenJetMatchRecluster"),
                                genParticles = cms.InputTag("prunedGenParticles"),
                                muons        = cms.InputTag("slimmedMuons"),
                                electrons    = cms.InputTag("slimmedElectrons"),
                                jetPtMin     = cms.untracked.double(30),
                                jetPtMax     = cms.untracked.double(-1),
                                jetAbsEtaMax = cms.untracked.double(2.4),
                                tagInfoName = cms.string('pfDeepCSV'),
                                fjTagInfoName = cms.string('pfBoostedDoubleSVAK8'),
                                bDiscriminators = cms.vstring(),
                                fjKeepFlavors = cms.untracked.vuint32(),
                                isQCDSample   = cms.untracked.bool(False)
                                )

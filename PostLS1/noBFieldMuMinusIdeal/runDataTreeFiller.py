import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.TFileService = cms.Service("TFileService",
                                   	fileName=cms.string('L1MuonDataTreeNoBMuMinusIdeal.root'),
                                   )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000))


import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile('cmsswSourceFiles')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(*mylist)
)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
process.load('Configuration.StandardSequences.MagneticField_0T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')


from TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi import *
#process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi')

from TrackingTools.TrackAssociator.default_cfi import TrackAssociatorParameterBlock

#L1Extra
process.load('L1Trigger.Configuration.L1Extra_cff')

#horeco
process.load('Configuration.StandardSequences.Reconstruction_cff')

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = "DESIGN72_V2::All"
print process.GlobalTag.globaltag

parameters = TrackAssociatorParameterBlock.TrackAssociatorParameters
parameters.useEcal = False
parameters.useHcal = False
parameters.useMuon = False

process.options = cms.untracked.PSet(

)

process.hoDataTreeFiller = cms.EDAnalyzer(
    'DataTreeFiller',
    genSrc = cms.InputTag("genParticles"),
    l1MuonSrc=cms.InputTag("l1extraParticles"),
    horecoSrc = cms.InputTag("horeco"),
    hoEnergyThreshold = cms.double(0.2),
	maxDeltaR = cms.double(0.3),
	debug = cms.bool(True),
	TrackAssociatorParameters=parameters,
	hoDigiSrc = cms.InputTag('simHcalDigis'),
	hoAdcThreshold = cms.int32(60)
    )

process.genfilter = cms.EDFilter("MCSingleParticleFilter",
#	Status = cms.untracked.vint32(1,1),
	MinPt = cms.untracked.vdouble(5.0,5.0),
	MinEta = cms.untracked.vdouble(-0.8,-0.8),
	MaxEta = cms.untracked.vdouble(0.8,0.8),
	ParticleID = cms.untracked.vint32(13,-13),
	 )

#Try using different source for hoReco
process.horeco.digiLabel = cms.InputTag('simHcalDigis')

#Path definitions
process.genFilter_step = cms.Path(process.genfilter)
process.horeco_step = cms.Path(process.horeco)
process.analyzer_step = cms.Path(process.hoDataTreeFiller)

process.p = cms.Path(process.genfilter*
					process.horeco*
					process.hoDataTreeFiller)

#Schedule Definition
process.schedule = cms.Schedule(
	process.p
	)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions

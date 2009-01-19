import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD2")

process.load("MagneticField.Engine.volumeBasedMagneticField_1103l_cfi")

process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

#process.load("Geometry.CommonDetUnit.bareGlobalTrackingGeometry_cfi")
#process.load("Geometry.DTGeometry.dtGeometry_cfi")
#process.load("Geometry.CSCGeometry.cscGeometry_cfi")
#process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
process.load("RecoMuon.TrackingTools.MuonServiceProxy_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'IDEAL_V5::All'


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    debugFlag = cms.untracked.bool(False),
    debugVebosity = cms.untracked.uint32(10),
    fileNames = cms.untracked.vstring(
#'rfio:/castor/cern.ch/user/s/sckao/scale1000_210_2a.root',
'rfio:/castor/cern.ch/user/s/sckao/scale10_210_1a.root'
    )
)

process.MuonNumberingInitialization = cms.ESProducer("MuonNumberingInitialization")

process.muonseed = cms.EDFilter("MuonSeedValidator",
    process.MuonServiceProxy,
    debug = cms.untracked.bool(False),
    muonseedLabel = cms.untracked.string('MuonSeed'),
    glbTrackLabel = cms.InputTag("globalMuons"),
    staTrackLabel = cms.InputTag("standAloneMuons"),
    simHitLabel = cms.untracked.string('g4SimHits'),
    simTrackLabel = cms.untracked.string('g4SimHits'),
    recHitLabel = cms.untracked.string('csc2DRecHits'),
    dtrecHitLabel = cms.untracked.string('dt1DRecHits'),
    dtSegmentLabel = cms.untracked.string('dt4DSegments'),
    cscSegmentLabel = cms.untracked.string('cscSegments'),
    dfMax = cms.untracked.double(0.5),
    dtMax = cms.untracked.double(0.35),
    scope = cms.untracked.bool(False),
    pTCutMin = cms.untracked.double(0.0),
    pTCutMax = cms.untracked.double(30.0),
    eta_Low = cms.untracked.double(0.0),
    eta_High = cms.untracked.double(2.5),
    rootFileName = cms.untracked.string('h_mp10_210a.root')
)

process.p = cms.Path(process.muonseed)


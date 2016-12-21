// -*- C++ -*-
// 
/**

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
 */
//
// Original Author:  Christopher Anelli  
//         Created:  Fri, 16 May 2014 04:20:05 GMT
// $Id$
//
//

//hoMuonAnalyzer header file
#include "../interface/hoMuonAnalyzer.h"

#include <DataFormats/CaloRecHit/interface/CaloRecHit.h>
#include <DataFormats/Candidate/interface/LeafCandidate.h>
#include <DataFormats/Common/interface/HandleBase.h>
#include <DataFormats/Common/interface/Ref.h>
#include <DataFormats/Common/interface/RefToBase.h>
#include <DataFormats/Common/interface/SortedCollection.h>
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/DetId/interface/DetId.h>
#include <DataFormats/GeometryVector/interface/GlobalPoint.h>
#include <DataFormats/GeometryVector/interface/GlobalVector.h>
#include "DataFormats/GeometryVector/interface/GlobalTag.h"
#include <DataFormats/GeometryVector/interface/Phi.h>
#include <DataFormats/GeometryVector/interface/PV3DBase.h>
#include <DataFormats/HcalDetId/interface/HcalDetId.h>
#include <DataFormats/HcalDetId/interface/HcalSubdetector.h>
#include <DataFormats/HcalRecHit/interface/HORecHit.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include <DataFormats/L1Trigger/interface/L1MuonParticle.h>
#include <DataFormats/Math/interface/deltaR.h>
#include <DataFormats/Provenance/interface/EventID.h>
#include <DataFormats/TrajectorySeed/interface/PropagationDirection.h>
#include <FWCore/Common/interface/EventBase.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/EventSetupRecord.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <FWCore/Framework/src/Factory.h>
#include <FWCore/Framework/src/WorkerMaker.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/ParameterSet/interface/ParameterSetDescription.h>
#include <FWCore/ParameterSet/interface/ParameterSetDescriptionFiller.h>
#include <FWCore/ParameterSet/interface/ParameterSetDescriptionFillerPluginFactory.h>
#include <FWCore/PluginManager/interface/PluginFactory.h>
#include <Geometry/CaloGeometry/interface/CaloGeometry.h>
#include <Geometry/Records/interface/CaloGeometryRecord.h>
#include <MagneticField/Records/interface/IdealMagneticFieldRecord.h>
#include <Math/GenVector/DisplacementVector3D.h>
#include <RecoMuon/MuonIdentification/interface/MuonHOAcceptance.h>
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSimpleRecAlgo.h"

#include <SimDataFormats/CaloHit/interface/PCaloHit.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>
#include <TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h>
#include <TrackingTools/TrackAssociator/interface/TrackDetMatchInfo.h>
#include <TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h>
#include <TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h>
#include <TROOT.h>
#include <TTree.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <utility>
#include <math.h>

#include "../interface/FilterPlugin.h"

#include "CalibFormats/HcalObjects/interface/HcalDbService.h"

#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
using namespace::std;

const float hoMuonAnalyzer::MAX_ETA = 0.8;
const float hoMuonAnalyzer::L1PHI_OFFSET = 3.1415926535897931/144.;

hoMuonAnalyzer::hoMuonAnalyzer(const edm::ParameterSet& iConfig){
	coutPrefix = "[hoMuonAnalyzer] ";
	//now do what ever initialization is needed

	//Get Input Tags from the Configuration
	isData = iConfig.getParameter<bool>("isData");

	//That stuff is only available if we use MC events
	if(!isData){
		_genInput = iConfig.getParameter<edm::InputTag>("genSrc");
		_l1MuonGenMatchInput = iConfig.getParameter<edm::InputTag>("l1MuonGenMatchSrc");

		assoc.useDefaultPropagator();
		edm::ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
		edm::ConsumesCollector iC = consumesCollector();
		assocParams.loadParameters( parameters, iC );
	}

	_l1MuonInput = iConfig.getParameter<edm::InputTag>("l1MuonSrc");
	_horecoInput = iConfig.getParameter<edm::InputTag>("horecoSrc");
	_hltSumAODInput = iConfig.getParameter<edm::InputTag>("hltSumAODSrc");
	deltaR_Max = iConfig.getParameter<double>("maxDeltaR");
	threshold = iConfig.getParameter<double>("hoEnergyThreshold");
	debug = iConfig.getParameter<bool>("debug");
	deltaR_L1MuonMatching = iConfig.getParameter<double>("maxDeltaRL1MuonMatching");
	useArtificialPrimaryVertex = iConfig.getParameter<bool>("useArtificialPrimaryVertex");

	singleMuOpenTrigName = "L1_SingleMuOpen";
	singleMu3TrigName = "L1_SingleMu3";
	singleMu7TrigName = "L1_SingleMu7";
	singleMu10TrigName = "L1_SingleMu10";

	doubleMu0TrigName = "L1_DoubleMu_10_Open";
	doubleMu5TrigName = "L1_DoubleMu5 ";

	defineTriggersOfInterest();

	hoMatcher = new HoMatcher(iConfig);
	functionsHandler = new CommonFunctionsHandler(iConfig);

	firstRun = true;

}


hoMuonAnalyzer::~hoMuonAnalyzer()
{
	delete hoMatcher;
	delete functionsHandler;
	hoMatcher = 0;
	functionsHandler = 0;
}

// ------------ method called for each event  ------------
void
hoMuonAnalyzer::analyze(const edm::Event& iEvent, 
		const edm::EventSetup& iSetup){
	if(firstRun){
		//for now we do not need the channel qualities any more
		//Uncomment again, if needed
//		printChannelQualities(iSetup);
		firstRun = false;
	}

	if(!isData){
		iEvent.getByLabel(_genInput,truthParticles);
		iEvent.getByLabel(_l1MuonGenMatchInput,l1MuonGenMatches);
	}


	iEvent.getByLabel(_l1MuonInput, l1Muons);
	iEvent.getByLabel(_horecoInput, hoRecoHits);
	iEvent.getByLabel(_l1MuonInput,l1MuonView);
	iEvent.getByLabel(edm::InputTag("muons"),recoMuons);

	iEvent.getByLabel(edm::InputTag("selectedPatMuons"),patMuons);
	iEvent.getByLabel(edm::InputTag("gtDigis","DT"),dtRegionalCands);

	iSetup.get<IdealMagneticFieldRecord>().get(theMagField );
	iSetup.get<CaloGeometryRecord>().get(caloGeo);

	hoMatcher->getEvent(iEvent,iSetup);
	functionsHandler->getEvent(iEvent);

	//Try getting the event info for weights
	edm::Handle<GenEventInfoProduct> genEventInfo;
	iEvent.getByLabel(edm::InputTag("generator"), genEventInfo);

	if (!MuonHOAcceptance::Inited()) MuonHOAcceptance::initIds(iSetup);

	/*
	 * Set Up Level 1 Global Trigger Utility
	 */

	bool useL1EventSetup = true;
	bool useL1GtTriggerMenuLite = true;

	m_l1GtUtils.getL1GtRunCache(iEvent, iSetup, useL1EventSetup, useL1GtTriggerMenuLite);

	/*
	 * Get the GMT readout record
	 */
	iEvent.getByLabel(edm::InputTag("gtDigis"), pCollection);

	if ( !pCollection.isValid() ) {
		edm::LogError("DataNotFound") << "can't find L1MuGMTReadoutCollection";
	}

	/*
	 *
	 *  Start of Analysis
	 *
	 */

	histogramBuilder.fillCountHistogram("ProcessedEvents");

	//Use this variable to store whether the event has muons in acceptance
	bool hasMuonsInAcceptance = false;
	if(isData){
		auto recoMuon = recoMuons->begin();
		for(; recoMuon!= recoMuons->end(); ++recoMuon){
			if(fabs(recoMuon->eta()) <= MAX_ETA){
				hasMuonsInAcceptance = true;
			}
		}
	}
	//Assume, that we simulated muons only in our preferred acceptance
	else{
		hasMuonsInAcceptance = true;
	}

	if(!hasMuonsInAcceptance){
		if(debug){
			std::cout << coutPrefix << "No L1Muons to process" << std::endl;
		}
		return;
	}

	/**
	 * Debug stuff. Look into DTTF cands to see if there are at all
	 * any cands with eta fine bit set
	 */
	// get GMT readout collection
	L1MuGMTReadoutCollection const* gmtrc = pCollection.product();
	std::vector<L1MuGMTReadoutRecord> gmt_records = gmtrc->getRecords();
	std::vector<L1MuGMTReadoutRecord>::const_iterator RRItr;
	int fineCounter = 0;
	for ( RRItr = gmt_records.begin(); RRItr != gmt_records.end(); ++RRItr ) {
		std::vector<L1MuRegionalCand> dttfCands = RRItr->getDTBXCands();
		std::vector<L1MuRegionalCand>::iterator dttfCand;
		for( dttfCand = dttfCands.begin(); dttfCand != dttfCands.end();	++dttfCand ) {
			if(dttfCand->empty()) continue;
			if(fabs(dttfCand->etaValue()) <= MAX_ETA ){
				histogramBuilder.fillEtaPhiGraph(dttfCand->etaValue(),dttfCand->phiValue(),"DttfCands");
				if(dttfCand->isFineHalo()){
					fineCounter++;
					histogramBuilder.fillEtaPhiGraph(dttfCand->etaValue(),dttfCand->phiValue(),"DttfCands_Fine");
				} else {
					histogramBuilder.fillEtaPhiGraph(dttfCand->etaValue(),dttfCand->phiValue(),"DttfCands_NotFine");
				}
			}
		}
	}
	histogramBuilder.fillMultiplicityHistogram(fineCounter,"fineBitsPerEvent");

	analyzeL1Sources();

    iEvent.getByLabel("hcalDigis", hoTPDigis);
    if(hoTPDigis.isValid()){
    	analyzeHoTriggerPrimitives();
    }
	iEvent.getByLabel(edm::InputTag("offlinePrimaryVertices"), vertexColl);
	if(vertexColl.isValid()){
		histogramBuilder.fillMultiplicityHistogram(vertexColl->size(),"primaryVertexSize");
	}
	iEvent.getByLabel(edm::InputTag("offlineBeamSpot"),recoBeamSpotHandle);

	histogramBuilder.fillCountHistogram("Events");
	if(!isData)
		processGenInformation(iEvent,iSetup);
	processRecoInformation(iEvent,iSetup);
	analyzeEnergyDeposit(iEvent,iSetup);
    analyzeTimingSupport();

	//###############################
	// Loop over L1MuonObjects DONE
	//###############################

	//###############################
	// BEGIN Loop over HO Rec hits only
	//###############################
	string horeco_key = "horeco";
	string horecoT_key ="horecoAboveThreshold";
	/**
	 * Collect information of all HO Rec hits when there are
	 * l1 muons in the acceptance region
	 */
	int recHitAbThrCounter = 0;
	auto hoRecoIt = hoRecoHits->begin();
	for( ; hoRecoIt != hoRecoHits->end() ; hoRecoIt++){
		double ho_eta = caloGeo->getPosition(hoRecoIt->id()).eta();
		double ho_phi = caloGeo->getPosition(hoRecoIt->id()).phi();
		histogramBuilder.fillTimeHistogram(hoRecoIt->time(),"hoRecHits");
		if(hoRecoIt->energy() >= threshold){
			histogramBuilder.fillTimeHistogram(hoRecoIt->time(),"hoRecHitsAboveThr");
			histogramBuilder.fillEtaPhiHistograms(ho_eta, ho_phi,"hoRecHitsAboveThr");
			int hoIEta = hoRecoIt->id().ieta();
			int hoIPhi = hoRecoIt->id().iphi();
			histogramBuilder.fillIEtaIPhiHistogram(hoIEta,hoIPhi,"hoRecHitsAboveThr");
			recHitAbThrCounter++;
		}
		histogramBuilder.fillCountHistogram(horeco_key);
		histogramBuilder.fillEnergyVsIEta(hoRecoIt->energy(),hoRecoIt->id().ieta(),horeco_key);
		histogramBuilder.fillEtaPhiHistograms(ho_eta, ho_phi, horeco_key);
	}
	histogramBuilder.fillMultiplicityHistogram(recHitAbThrCounter,horecoT_key);
	histogramBuilder.fillMultiplicityHistogram(hoRecoHits->size(),horeco_key);
	makeHoRecHitThresholdScan();
	/*
	 * L1 Trigger Decisions
	 */
	singleMuOpenTrig = processTriggerDecision(singleMuOpenTrigName,iEvent);
//	singleMu3Trig = processTriggerDecision(singleMu3TrigName,iEvent);
//	singleMu7Trig = processTriggerDecision(singleMu7TrigName,iEvent);
//	singleMu10Trig = processTriggerDecision(singleMu10TrigName,iEvent);
	doubleMu0Trig = processTriggerDecision(doubleMu0TrigName,iEvent);
	//	processTriggerDecision(doubleMu5TrigName,iEvent);
	//###############################
	// Loop over HO Rec hits only DONE
	//###############################

	//###############################
	// Loop over L1 Muons and match to HO
	//###############################

	int countGenMatches = 0;
	string l1MuonWithHoMatch_key = "L1MuonWithHoMatch";
	histogramBuilder.fillMultiplicityHistogram(l1Muons->size(),"L1MuonPresent");
	if(l1Muons->size() > 0){
		histogramBuilder.fillCountHistogram("L1MuonPresent");
	}
	for( unsigned int i = 0 ; i < l1Muons->size(); i++ ){
		const l1extra::L1MuonParticle* bl1Muon = &(l1Muons->at(i));
		float l1Muon_eta = bl1Muon->eta();
		float l1Muon_phi = bl1Muon->phi() + L1PHI_OFFSET;
		if(fabs(l1Muon_eta) > MAX_ETA){
			continue;
		}
		histogramBuilder.fillCountHistogram("L1Muon");
		histogramBuilder.fillBxIdHistogram(bl1Muon->bx(),"L1MuonPresent");
		histogramBuilder.fillBxIdVsPt(bl1Muon->bx(),bl1Muon->pt(),"L1MuonPresent");
		histogramBuilder.fillL1MuonPtHistograms(bl1Muon->pt(),"L1MuonPresent");
		histogramBuilder.fillPdgIdHistogram(bl1Muon->pdgId(),"L1MuonPresent");
		histogramBuilder.fillVzHistogram(bl1Muon->vz(),"L1MuonPresent");
		histogramBuilder.fillEtaPhiGraph(l1Muon_eta, l1Muon_phi, "L1MuonPresent");

		// ####################################################################################
		// Fill ho energy in direction opposite to L1Muon when there is no reco muon found
		// Trying to measure the HO Noise
		// ####################################################################################
		const pat::Muon* patMuon = getBestPatMatch(l1Muon_eta, l1Muon_phi + M_PI);
		if(!patMuon){//Do not want muons to interfere
			for(auto hoRecHitIt = hoRecoHits->begin() ; hoRecHitIt != hoRecoHits->end(); hoRecHitIt++){
				if(hoMatcher->isRecHitInGrid(l1Muon_eta,l1Muon_phi + M_PI,&*hoRecHitIt,5)){
					histogramBuilder.fillEnergyVsIEta(hoRecHitIt->energy(),hoRecHitIt->id().ieta(),"L1MuonPlusPi");
				}
			}
		}

		/**
		 * Fill eta phi graphs with coordinates of L1Muons depending on their detector/quality.
		 * No bias by matching to HO introduced here
		 */
		const L1MuRegionalCand* l1RegCand = findBestCandMatch(bl1Muon);
		switch (bl1Muon->gmtMuonCand().quality()) {
			case 7:
				histogramBuilder.fillEtaPhiGraph(l1Muon_eta, l1Muon_phi, "L1MuonPresent_DT_RPC");
				if(l1RegCand){
					if(l1RegCand->isFineHalo()){
						histogramBuilder.fillEtaPhiGraph(bl1Muon->eta(), bl1Muon->phi() + L1PHI_OFFSET,"L1MuonPresent_DT_RPC_Fine");
					} else {
						histogramBuilder.fillEtaPhiGraph(bl1Muon->eta(), bl1Muon->phi() + L1PHI_OFFSET,"L1MuonPresent_DT_RPC_NotFine");
					}
				}
				break;
			case 6:
				histogramBuilder.fillEtaPhiGraph(l1Muon_eta, l1Muon_phi, "L1MuonPresent_DT");
				break;
			default:
				histogramBuilder.fillEtaPhiGraph(l1Muon_eta, l1Muon_phi, "L1MuonPresent_Other");
				break;
		}

		fillAverageEnergyAroundL1Direction(bl1Muon,"L1MuonPresent");
		//For variable binning
//		listL1MuonPt.push_back(bl1Muon->pt());
		/*
		 * Fill histogram for different pt thresholds
		 * CAUTION!! THIS IS NOT A REAL RATE YET!!
		 */
		for (int j = 0; j < 200; j+=2) {
			if(bl1Muon->pt() >= j){
				histogramBuilder.fillTrigRateHistograms(j,"L1MuonPresent");
			}
		}
		//Look for matches in grid around L1
		calculateGridMatchingEfficiency(&*bl1Muon,bl1Muon->pt(),"L1Muon");
		fillGridMatchingQualityCodes(&*bl1Muon,bl1Muon->pt(),"L1Muon");

		if(MuonHOAcceptance::inGeomAccept(l1Muon_eta,l1Muon_phi)&& !hoMatcher->isInChimney(l1Muon_eta,l1Muon_phi)){
			histogramBuilder.fillCountHistogram("L1MuonInGA_L1Dir");
			if(bl1Muon->bx() == 0)
				histogramBuilder.fillCountHistogram("L1MuoninGaBx0");
		}
		const HORecHit* matchedRecHit = hoMatcher->matchByEMaxInGrid(l1Muon_eta,l1Muon_phi,2);
		if(matchedRecHit){
			double hoEta,hoPhi;
			hoEta = caloGeo->getPosition(matchedRecHit->detid()).eta();
			hoPhi = caloGeo->getPosition(matchedRecHit->detid()).phi();
			histogramBuilder.fillCountHistogram("L1MuonPresentHoMatch");
			histogramBuilder.fillTimeHistogram(matchedRecHit->time(),"L1MuonPresentHoMatch");
			histogramBuilder.fillDeltaTimeHistogram(matchedRecHit->time(),bl1Muon->bx(),"L1MuonPresentHoMatch");
			histogramBuilder.fillBxIdHistogram(bl1Muon->bx(),"L1MuonPresentHoMatch");
			histogramBuilder.fillEnergyVsIEta(matchedRecHit->energy(),matchedRecHit->id().ieta(),l1MuonWithHoMatch_key);
			histogramBuilder.fillEtaPhiHistograms(hoEta, hoPhi,l1MuonWithHoMatch_key);
			histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1Muon_eta,hoEta,l1Muon_phi, hoPhi,l1MuonWithHoMatch_key);
			histogramBuilder.fillEnergyVsPosition(hoEta,hoPhi,matchedRecHit->energy(),l1MuonWithHoMatch_key);
			if (MuonHOAcceptance::inGeomAccept(l1Muon_eta,l1Muon_phi/*,deltaR_Max,deltaR_Max*/)&& !hoMatcher->isInChimney(l1Muon_eta,l1Muon_phi)){
				histogramBuilder.fillCountHistogram("L1MuonPresentHoMatchInAcc");
				histogramBuilder.fillBxIdHistogram(bl1Muon->bx(),"L1MuonPresentHoMatchInAcc");

				//This energy check is done to test, whether the results depend on the order of the cuts applied
				//So far, the answer is no
				if(matchedRecHit->energy() >= threshold ){
					histogramBuilder.fillCountHistogram("L1MuonPresentHoMatchInAccThr");

				}

				if (MuonHOAcceptance::inNotDeadGeom(l1Muon_eta,l1Muon_phi/*,deltaR_Max,deltaR_Max*/)){
					histogramBuilder.fillCountHistogram("L1MuonPresentHoMatchInAccNotDead");
					histogramBuilder.fillBxIdHistogram(bl1Muon->bx(),"L1MuonPresentHoMatchInAccNotDead");
					//This one is filled for the sake of completeness. The SiPM regions are hard-coded in the class!!
					if (MuonHOAcceptance::inSiPMGeom(l1Muon_eta,l1Muon_phi/*,deltaR_Max,deltaR_Max*/)){
						histogramBuilder.fillCountHistogram("L1MuonPresentHoMatchInAccNotDeadInSipm");
					}
				}
			}
			/**
			 * Dump events where delta i phi is 1
			 * Maybe cms show can help for the systematic shift
			 */
			if(hoMatcher->getDeltaIphi(l1Muon_phi,matchedRecHit) == 1){
				ofstream myfile;
				myfile.open ("deltaPhiOneEvents.txt",std::ios::app);
				myfile << iEvent.id().event() << std::endl;
			}
			//Pseudo trigger rate
			for (int i = 0; i < 200; i+=2) {
				if(bl1Muon->pt() >= i)
					histogramBuilder.fillTrigRateHistograms(i,"L1MuonWithHoNoThr");
			}
			if(!isData){
				const reco::GenParticle* bestGenMatch = getBestGenMatch(l1Muon_eta,l1Muon_phi);
				if(bestGenMatch){
					countGenMatches++;
					fillEfficiencyHistograms(bl1Muon->pt(),bestGenMatch->pt(),"L1MuonHoReco");
				}
			}
			//###########################################################
			//###########################################################
			//Now fill information for hits above threshold
			//###########################################################
			//###########################################################
			if(matchedRecHit->energy() >= threshold){
				//Fill some counting histograms. Can be used for cut flow in efficiency
				histogramBuilder.fillCountHistogram("L1MuonAboveThr");
				histogramBuilder.fillBxIdHistogram(bl1Muon->bx(),"L1MuonAboveThr");
				histogramBuilder.fillDeltaTimeHistogram(matchedRecHit->time(),bl1Muon->bx(),"L1MuonAboveThr");
				histogramBuilder.fillTimeHistogram(matchedRecHit->time(),"L1MuonAboveThr");
				histogramBuilder.fillBxIdVsPt(bl1Muon->bx(),bl1Muon->pt(),"L1MuonAboveThr");
				histogramBuilder.fillEnergyVsIEta(matchedRecHit->energy(),matchedRecHit->id().ieta(),"L1MuonWithHoMatchAboveThr");
				histogramBuilder.fillEtaPhiHistograms(hoEta,hoPhi,"L1MuonWithHoMatchAboveThr_HO");
				histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1Muon_eta,hoEta,l1Muon_phi, hoPhi,"L1MuonWithHoMatchAboveThr");
				histogramBuilder.fillL1MuonPtHistograms(bl1Muon->pt(),"L1MuonWithHoMatchAboveThr");
				histogramBuilder.fillEnergyVsPosition(hoEta,hoPhi,matchedRecHit->energy(),"L1MuonWithHoMatchAboveThr");


				TH2D* hist = new TH2D("hoEnergyVsTime","HO Energy vs. Time;Time / ns;E_{Rec} / GeV",201,-100.5,100.5,2100, -5.0, 100.0);
				histogramBuilder.fillCorrelationHistogram(matchedRecHit->time(),matchedRecHit->energy(),"hoEnergyVsTime",hist);
				delete hist;

				//Make time correlation plots depending on the different detector subsystems
				switch (bl1Muon->gmtMuonCand().detector()) {
					//RPC
					case 1:
						histogramBuilder.fillDeltaTimeHistogram(matchedRecHit->time(),bl1Muon->bx(),"RpcHoAboveThr");
						break;
					//DT
					case 2:
						histogramBuilder.fillDeltaTimeHistogram(matchedRecHit->time(),bl1Muon->bx(),"DtHoAboveThr");
						break;
					//RPC/DT
					case 3:
						histogramBuilder.fillDeltaTimeHistogram(matchedRecHit->time(),bl1Muon->bx(),"RpcAndDtHoAboveThr");
						break;
					default:
						break;
				}
				double hoEta,hoPhi;
				hoEta = caloGeo->getPosition(matchedRecHit->detid()).eta();
				hoPhi = caloGeo->getPosition(matchedRecHit->detid()).phi();
				//Fill the HO information
				//Fill the counters
				if (MuonHOAcceptance::inGeomAccept(l1Muon_eta,l1Muon_phi/*,deltaR_Max,deltaR_Max*/)&& !hoMatcher->isInChimney(l1Muon_eta,l1Muon_phi)){
					histogramBuilder.fillCountHistogram("L1MuonAboveThrInAcc");
					histogramBuilder.fillBxIdHistogram(bl1Muon->bx(),"L1MuonAboveThrInAcc");
					if (MuonHOAcceptance::inNotDeadGeom(l1Muon_eta,l1Muon_phi/*,deltaR_Max,deltaR_Max*/)){
						histogramBuilder.fillTimeHistogram(matchedRecHit->time(),"L1MuonAboveThrInAccNotDead");
						histogramBuilder.fillDeltaTimeHistogram(matchedRecHit->time(),bl1Muon->bx(),"L1MuonAboveThrInAccNotDead");
						histogramBuilder.fillCountHistogram("L1MuonAboveThrInAccNotDead");
						histogramBuilder.fillBxIdHistogram(bl1Muon->bx(),"L1MuonAboveThrInAccNotDead");
						histogramBuilder.fillTrigHistograms(caloGeo->present(matchedRecHit->id()),"caloGeoPresent_L1MuonHoMatchAboveThrFilt");
						histogramBuilder.fillEtaPhiHistograms(hoEta,hoPhi,"L1MuonWithHoMatchAboveThrFilt_HO");
						histogramBuilder.fillDeltaEtaDeltaPhiHistograms(l1Muon_eta,hoEta,l1Muon_phi, hoPhi,"L1MuonWithHoMatchAboveThrFilt");
						histogramBuilder.fillL1MuonPtHistograms(bl1Muon->pt(),"L1MuonWithHoMatchAboveThrFilt");
						histogramBuilder.fillEnergyVsPosition(hoEta,hoPhi,matchedRecHit->energy(),"L1MuonWithHoMatchAboveThrFilt");
						//Make time correlation plots depending on the different detector subsystems
						//This time for HO in GA only
						switch (bl1Muon->gmtMuonCand().detector()) {
							//RPC
							case 1:
								histogramBuilder.fillDeltaTimeHistogram(matchedRecHit->time(),bl1Muon->bx(),"RpcHoAboveThrFilt");
								break;
							//DT
							case 2:
								histogramBuilder.fillDeltaTimeHistogram(matchedRecHit->time(),bl1Muon->bx(),"DtHoAboveThrFilt");
								break;
							//RPC/DT
							case 3:
								histogramBuilder.fillDeltaTimeHistogram(matchedRecHit->time(),bl1Muon->bx(),"RpcAndDtHoAboveThrFilt");
								break;
							default:
								break;
						}
						if(!isData){
							const reco::GenParticle* bestGenMatch = getBestGenMatch(l1Muon_eta,l1Muon_phi);
							if(bestGenMatch){
								fillEfficiencyHistograms(bl1Muon->pt(),bestGenMatch->pt(),"L1MuonHoRecoAboveThrFilt");
							}
						}
					}
				}

				//Make the pseudo trig rate plot
				for (int i = 0; i < 200; i+=2) {
					if(bl1Muon->pt() >= i)
						histogramBuilder.fillTrigRateHistograms(i,l1MuonWithHoMatch_key);
				}
				histogramBuilder.fillL1MuonPtHistograms(bl1Muon->pt(), l1MuonWithHoMatch_key);
				histogramBuilder.fillEtaPhiGraph(l1Muon_eta, l1Muon_phi, "L1MuonWithHoMatchAboveThr_L1Mu");

				if(!isData){
					const reco::GenParticle* bestGenMatch = getBestGenMatch(l1Muon_eta,l1Muon_phi);
					if(bestGenMatch){
						fillEfficiencyHistograms(bl1Muon->pt(),bestGenMatch->pt(),"L1MuonHoRecoAboveThr");
						histogramBuilder.fillBxIdVsPt(bl1Muon->bx(),bl1Muon->pt(),"L1MuonAboveThrGenMatch");
					}
					//Try to find a corresponding Gen Muon
					edm::RefToBase<l1extra::L1MuonParticle> l1MuonCandiateRef(l1MuonView,i);
					reco::GenParticleRef ref = (*l1MuonGenMatches)[l1MuonCandiateRef];

					if(ref.isNonnull()){
						string l1MuonAndHoRecoAndGenref_key = "L1MuonandHORecowithMipMatchAndGenMatch";
						//Make the pseudo trig rate plot
						for (int i = 0; i < 200; i+=2) {
							if(bl1Muon->pt() >= i)
								histogramBuilder.fillTrigRateHistograms(i,l1MuonAndHoRecoAndGenref_key);
						}
						histogramBuilder.fillPdgIdHistogram(ref->pdgId(),l1MuonWithHoMatch_key);
					} else{
						histogramBuilder.fillPdgIdHistogram(0,l1MuonWithHoMatch_key);
					}
				}
			}// E > thr.
		}
	}//<-- For loop over all l1muons
	if(!isData)
		histogramBuilder.fillMultiplicityHistogram(countGenMatches,"nL1WithGenMatch");
	//################################
	//################################
	//	This is for the case where no L1Muon was found
	//################################
	//################################
	if(l1Muons->size() == 0){
		histogramBuilder.fillCountHistogram("NoL1Muon");
		int recHitAbThrNoL1Counter = 0;
		auto hoRecoIt = hoRecoHits->begin();
		for( ; hoRecoIt != hoRecoHits->end() ; hoRecoIt++){
			if(hoRecoIt->energy() >= threshold){
				recHitAbThrNoL1Counter++;
				double ho_eta = caloGeo->getPosition(hoRecoIt->id()).eta();
				double ho_phi = caloGeo->getPosition(hoRecoIt->id()).phi();
				histogramBuilder.fillCountHistogram("NoL1");
				histogramBuilder.fillEnergyHistograms(hoRecoIt->energy(),"NoL1");
				histogramBuilder.fillEtaPhiHistograms(ho_eta, ho_phi, "NoL1");
				histogramBuilder.fillEnergyVsPosition(ho_eta,ho_phi,hoRecoIt->energy(),"NoL1");
			}
		}
		histogramBuilder.fillMultiplicityHistogram(recHitAbThrNoL1Counter,"NoL1");

		if (!isData) {
			for(reco::GenParticleCollection::const_iterator genIt = truthParticles->begin();
					genIt != truthParticles->end(); genIt++){
				//Check for muons in Full barrel only
				//Try to find a corresponding Gen Muon
				float genEta = genIt->eta();
				float genPhi = genIt->phi();

				/**
				 * #############################################
				 * # Use TrackDetMatchInfo to see where the gen particle ends up in HO. When selecting only events, where gen is in
				 * # the geometric acceptance of HO, this gives a more realistic estimation for the number of recoverable
				 * # Triggers
				 * #############################################
				 */
				TrackDetMatchInfo * muMatch = getTrackDetMatchInfo(*genIt,iEvent,iSetup);
				double muMatchPhi = muMatch->trkGlobPosAtHO.phi();
				double muMatchEta = muMatch->trkGlobPosAtHO.eta();
				delete muMatch;
				histogramBuilder.fillEtaPhiGraph(genEta,genPhi,"NoL1GenAny");
				histogramBuilder.fillEtaPhiGraph(muMatchEta,muMatchPhi,"NoL1TdmiAny");
				if(MuonHOAcceptance::inGeomAccept(muMatchEta,muMatchPhi)
					&& MuonHOAcceptance::inNotDeadGeom(muMatchEta,muMatchPhi)
					&& !hoMatcher->isInChimney(muMatchEta,muMatchPhi)
				){
					histogramBuilder.fillEtaPhiGraph(genEta,genPhi,"NoL1GenInGA");
					histogramBuilder.fillEtaPhiGraph(muMatchEta,muMatchPhi,"NoL1TdmiInGA");
				}
				const HORecHit* matchedRecHit = 0;
				matchedRecHit = hoMatcher->matchByEMaxInGrid(muMatchEta,muMatchPhi,5);
				if(matchedRecHit){
					histogramBuilder.fillDeltaEtaDeltaPhiHistogramsWithWeights(muMatchEta
							,float(hoMatcher->getRecHitEta(matchedRecHit))	,muMatchPhi
							,float(hoMatcher->getRecHitPhi(matchedRecHit))	,matchedRecHit->energy()
							,"averageEMaxAroundPoint_NoL1_Tdmi");
					histogramBuilder.fillDeltaEtaDeltaPhiEnergyHistogram(muMatchEta
							,float(hoMatcher->getRecHitEta(matchedRecHit))	,muMatchPhi
							,float(hoMatcher->getRecHitPhi(matchedRecHit))	,matchedRecHit->energy()
							,"averageEMaxAroundPoint_NoL1_Tdmi");
				}
			}
		}

	}// <-- l1muons size == 0
	else{
		/**
		 * #################################
		 * # L1 Muon objects contain data
		 * # Loop over l1 and try to find gens. This way, the direction information of the strange l1 is already available
		 * # Loop over gens and inspect also the case, where the matching from gen to l1 failed
		 * # This might be a hint on ghost reduction
		 * #################################
		 */
		if (!isData) {
			for(reco::GenParticleCollection::const_iterator genIt = truthParticles->begin();
					genIt != truthParticles->end(); genIt++){
				float genEta = genIt->eta();
				float genPhi = genIt->phi();

				TrackDetMatchInfo * muMatch = getTrackDetMatchInfo(*genIt,iEvent,iSetup);
				double muMatchPhi = muMatch->trkGlobPosAtHO.phi();
				double muMatchEta = muMatch->trkGlobPosAtHO.eta();

				if(MuonHOAcceptance::inGeomAccept(muMatchEta,muMatchPhi)&& !hoMatcher->isInChimney(muMatchEta,muMatchPhi)){
					histogramBuilder.fillCountHistogram("TdmiInGA_TdmiDir");
					std::vector<const HORecHit*> crossedHoRecHits = muMatch->crossedHORecHits;
					calculateGridMatchingEfficiency(muMatchEta,muMatchPhi,genIt->pt(),"tdmi");
				}
				delete muMatch;

				const l1extra::L1MuonParticle* l1Part = functionsHandler->getBestL1MuonMatch(muMatchEta,muMatchPhi);
				if(l1Part){
					double deltaEta = muMatchEta - l1Part->eta();
					double deltaPhi = FilterPlugin::wrapCheck(muMatchPhi, l1Part->phi() + L1PHI_OFFSET);
					histogramBuilder.fillGraph(deltaEta,deltaPhi,"deltaEtaDeltaPhiTdmiL1");
				}
				const HORecHit* matchedRecHit = hoMatcher->matchByEMaxInGrid(genEta,genPhi,2);
				if(matchedRecHit){
					double hoEta = caloGeo->getPosition(matchedRecHit->id()).eta();
					double hoPhi = caloGeo->getPosition(matchedRecHit->id()).phi();
					histogramBuilder.fillDeltaEtaDeltaPhiEnergyHistogram(genEta,hoEta,genPhi,hoPhi,matchedRecHit->energy(),"WithSingleMu");
					if(matchedRecHit->energy() >= threshold){
						histogramBuilder.fillDeltaEtaDeltaPhiHistograms(genEta,hoEta,genPhi,hoPhi,"WithSingleMuAboveThr");
						if (MuonHOAcceptance::inGeomAccept(genEta,genPhi/*,deltaR_Max,deltaR_Max*/)&& !hoMatcher->isInChimney(genEta,genPhi)){
							histogramBuilder.fillDeltaEtaDeltaPhiHistograms(genEta,hoEta,genPhi,hoPhi,"WithSingleMuGeomAcc");
							if (MuonHOAcceptance::inNotDeadGeom(genEta,genPhi/*,deltaR_Max,deltaR_Max*/)){
								histogramBuilder.fillDeltaEtaDeltaPhiHistograms(genEta,hoEta,genPhi,hoPhi,"WithSingleMuNotDead");
							}
						}
					}
				}//Matched rec hit
			}//Gen particle loop
		}
	}
	//#############################################################
	//#############################################################
	// NO SINGLE MU TRIGGER
	//#############################################################
	//#############################################################


	if(!singleMuOpenTrig){
		histogramBuilder.fillMultiplicityHistogram(l1Muons->size(),"NoSingleMu_L1Muon");
		if(!isData){
			analyzeNoSingleMuEventsL1Loop(iEvent,iSetup);
			analyzeNoSingleMuEventsGenLoop(iEvent,iSetup);
		}
		//################################
		//################################
		// Match Ho to gen info and try to recover mu trigger
		//################################
		//################################
		if (!isData) {
			for(reco::GenParticleCollection::const_iterator genIt = truthParticles->begin();
					genIt != truthParticles->end(); genIt++){
				//Check for muons in Full barrel only
				//Try to find a corresponding Gen Muon
				float genEta = genIt->eta();
				float genPhi = genIt->phi();

				/**
				 * #############################################
				 * # Use TrackDetMatchInfo to see where the gen particle ends up in HO. When selecting only events, where gen is in
				 * # the geometric acceptance of HO, this gives a more realistic estimator for the number of recoverable
				 * # Triggers
				 * #############################################
				 */
				TrackDetMatchInfo * muMatch = getTrackDetMatchInfo(*genIt,iEvent,iSetup);
				double muMatchPhi = muMatch->trkGlobPosAtHO.phi();
				double muMatchEta = muMatch->trkGlobPosAtHO.eta();

				//Studies depending on whether there are L1 Objects or not
				if(l1Muons->size()>0){
					histogramBuilder.fillCountHistogram("NoTriggerButL1Muons");
				}
				else{
					histogramBuilder.fillCountHistogram("NoTrgNoL1Any");
					histogramBuilder.fillEtaPhiGraph(genEta,genPhi,"NoTrgNoL1GenAny");
					histogramBuilder.fillEtaPhiGraph(muMatchEta,muMatchPhi,"NoTrgNoL1TdmiAny");
					if(MuonHOAcceptance::inGeomAccept(muMatchEta,muMatchPhi) && !hoMatcher->isInChimney(muMatchEta,muMatchPhi)){
						histogramBuilder.fillEtaPhiGraph(genEta,genPhi,"NoTrgNoL1GenInGA");
						histogramBuilder.fillEtaPhiGraph(muMatchEta,muMatchPhi,"NoTrgNoL1TdmiInGA");
					}
				}

				histogramBuilder.fillEtaPhiGraph(genEta,genPhi,"NoTrgGenAny");
				histogramBuilder.fillEtaPhiGraph(muMatchEta,muMatchPhi,"NoTrgTdmiAny");
				//The muon needs to hit the HO geometric acceptance
				if(MuonHOAcceptance::inGeomAccept(muMatchEta,muMatchPhi) && !hoMatcher->isInChimney(muMatchEta,muMatchPhi)){
					histogramBuilder.fillEtaPhiGraph(muMatchEta,muMatchPhi,"NoTrgTdmiInGA");
					histogramBuilder.fillCountHistogram("NoTrgTdmiInGA");
					histogramBuilder.fillEnergyVsPosition(muMatchEta,muMatchPhi,muMatch->hoCrossedEnergy(),"NoTrgTdmiXedE");
					const HORecHit* matchedRecHit = hoMatcher->matchByEMaxInGrid(muMatchEta,muMatchPhi,2);
					//Where is the Rec hit in a delta R cone with the largest E?
					//Did we find any?
					if(matchedRecHit){
						histogramBuilder.fillCountHistogram("NoTrgTdmiInGAHoMatch");
						double hoEta = caloGeo->getPosition(matchedRecHit->id()).eta();
						double hoPhi = caloGeo->getPosition(matchedRecHit->id()).phi();
						histogramBuilder.fillDeltaEtaDeltaPhiHistograms(muMatchEta,hoEta,muMatchPhi,hoPhi,"NoTrgTdmi");
						//Apply energy cut on the matched RecHit
						if(matchedRecHit->energy() >= threshold ){
							histogramBuilder.fillCountHistogram("NoTrgTdmiInGAHoAboveThr");
							histogramBuilder.fillDeltaEtaDeltaPhiHistograms(muMatchEta,hoEta,muMatchPhi,hoPhi,"NoTrgTdmiAboveThr");
							histogramBuilder.fillEtaPhiGraph(muMatchEta,muMatchPhi,"NoTrgTdmiAboveThr");
							histogramBuilder.fillEtaPhiPtHistogram(muMatchEta,muMatchPhi,genIt->pt(),"NoTrgTdmiAboveThr");
							histogramBuilder.fillEtaPhiGraph(hoEta,hoPhi,"NoTrgTdmiAboveThrHoCoords");
							histogramBuilder.fillDeltaEtaDeltaPhiHistograms(genEta,hoEta,genPhi,hoPhi,"NoTrgGenAboveThr");
							histogramBuilder.fillEnergyVsPosition(muMatchEta,muMatchPhi,muMatch->hoCrossedEnergy(),"NoTrgTdmiAboveThrXedE");
							histogramBuilder.fillEnergyVsPosition(hoEta,hoPhi,matchedRecHit->energy(),"NoTrgTdmiAboveThrHoE");
							if( (muMatchEta > -0.35 && muMatchEta < -0.185) || (muMatchEta > 0.16 && muMatchEta < 0.3) ){
								if( (muMatchPhi > 0.7 && muMatchPhi < 1.36) || (muMatchPhi > 1.2 && muMatchPhi < 1.9) ){
									ofstream myfile;
									myfile.open ("eventNumbers.txt",std::ios::app);
									myfile << iEvent.id().event() << std::endl;
								}
							}
						//inspect the crossed energy, when the matched Rec hit in the cone was below threshold
						} else{
							histogramBuilder.fillEnergyVsPosition(muMatchEta,muMatchPhi,muMatch->hoCrossedEnergy(),"NoTrgTdmiBelowThrXedE");
							histogramBuilder.fillEtaPhiGraph(muMatchEta,muMatchPhi,"NoTrgTdmiBelowThr");
							histogramBuilder.fillEtaPhiPtHistogram(muMatchEta,muMatchPhi,genIt->pt(),"NoTrgTdmiBelowThr");
						}
					//Count the events, where we could not match a Rec hit in the delta dR cone
					} else{
						histogramBuilder.fillCountHistogram("NoTrgHoMatchFail");
						histogramBuilder.fillEtaPhiGraph(muMatchEta,muMatchPhi,"NoTrgHoMatchFail");
						histogramBuilder.fillEtaPhiPtHistogram(muMatchEta,muMatchPhi,genIt->pt(),"NoTrgHoMatchFail");
					}
				}//<-- in GA
				else{
					histogramBuilder.fillEtaPhiGraph(muMatchEta,muMatchPhi,"NoTrgTdmiNotInGA");
					histogramBuilder.fillEtaPhiPtHistogram(muMatchEta,muMatchPhi,genIt->pt(),"NoTrgTdmiNotInGA");
				}
				delete muMatch;

			}//Loop over gen particles
		}//Only in case of MC events
	}//<-- Not single mu trg
	//#############################################################
	//#############################################################
	// SINGLE MU TRIGGER FIRED
	//#############################################################
	//#############################################################
	else {
		/**
		 * #######################
		 * ######################
		 * The following part serves for checking cases where there is a L1 Muon trigger
		 * It's possible that some code is similar or the same as above. But for clarity reasons
		 * the code below is used
		 * #######################
		 * ######################
		 */
		if (!isData) {
			for(reco::GenParticleCollection::const_iterator genIt = truthParticles->begin();
					genIt != truthParticles->end(); genIt++){
				float genEta = genIt->eta();
				float genPhi = genIt->phi();

				TrackDetMatchInfo * muMatch = getTrackDetMatchInfo(*genIt,iEvent,iSetup);
				double muMatchPhi = muMatch->trkGlobPosAtHO.phi();
				double muMatchEta = muMatch->trkGlobPosAtHO.eta();
				delete muMatch;

				//Require the muon to hit the HO area
				if(MuonHOAcceptance::inGeomAccept(muMatchEta,muMatchPhi) && !hoMatcher->isInChimney(muMatchEta,muMatchPhi)){

					histogramBuilder.fillCountHistogram("SMuTrgTdmiInGA");
					const l1extra::L1MuonParticle* l1Ref = functionsHandler->getBestL1MuonMatch(genEta,genPhi);
					if(l1Ref){
						histogramBuilder.fillCountHistogram("SMuTrgFoundL1Match");
						float l1Muon_eta = l1Ref->eta();
						float l1Muon_phi = l1Ref->phi() + L1PHI_OFFSET;
						fillEfficiencyHistograms(l1Ref->pt(),genIt->pt(),"SMuTrgL1AndGenMatch");
						const HORecHit* matchedRecHit = hoMatcher->matchByEMaxInGrid(l1Muon_eta,l1Muon_phi,2);
						//Check whether an HO match could be found by the delta R method
						if(matchedRecHit){
							histogramBuilder.fillCountHistogram("SMuTrgL1AndFoundHoMatch");
							//inspect the data where HO was above the threshold
							if(matchedRecHit->energy() >= threshold){
								histogramBuilder.fillCountHistogram("SMuTrgL1AndHoAboveThr");
								fillEfficiencyHistograms(l1Ref->pt(),genIt->pt(),"SMuTrgL1AndGenMatchHoAboveThr");
							}
						}
					}//<-- Found L1 ref
					else{
						const HORecHit* matchedRecHit = hoMatcher->matchByEMaxInGrid(muMatchEta,muMatchPhi,2);
						//Check whether an HO match could be found by the delta R method
						if(matchedRecHit){
							histogramBuilder.fillCountHistogram("SMuTrgAndFoundHoMatch");
							//inspect the data where HO was above the threshold
							if(matchedRecHit->energy() >= threshold){
								histogramBuilder.fillCountHistogram("SMuTrgAndHoAboveThr");
							}
						}
					}
				}
			}//Gen loop
		}


		//Loop over L1 Objects
		//##################################################
		//##################################################
		// L1 Muons for the "Oliver Style" efficiency
		//##################################################
		//##################################################
		for( unsigned int i = 0 ; i < l1Muons->size(); i++  ) {
			const l1extra::L1MuonParticle* bl1Muon = &(l1Muons->at(i));
			double l1Muon_eta = bl1Muon->eta();
			double l1Muon_phi = bl1Muon->phi() + L1PHI_OFFSET;
			if(MuonHOAcceptance::inGeomAccept(l1Muon_eta,l1Muon_phi)&& !hoMatcher->isInChimney(l1Muon_eta,l1Muon_phi)){
				histogramBuilder.fillCountHistogram("L1MuonSMuTrgInGA_L1Dir");
				if(hoMatcher->hasHoHitInGrid(l1Muon_eta,l1Muon_phi,0)){
					histogramBuilder.fillCountHistogram("L1MuonSMuTrgCentral");
				}
				if(hoMatcher->hasHoHitInGrid(l1Muon_eta,l1Muon_phi,1)){
					histogramBuilder.fillCountHistogram("L1MuonSMuTrg3x3");
				}
			}
		}
	}//<-- End of Single mu trg
}
//#############################
// End of Analyze function
//#############################


/**
 * Small helper function to print the number of triggers for a certain algorithm name
 */
bool hoMuonAnalyzer::processTriggerDecision(std::string algorithmName,const edm::Event& iEvent){
	// Select on events that pass a specific L1Trigger Decision
	int iErrorCode = -1;
	bool trigDecision = m_l1GtUtils.decision(iEvent, algorithmName, iErrorCode);
	if(iErrorCode == 0){
		histogramBuilder.fillTrigHistograms(trigDecision,algorithmName);
		if(trigDecision){
			histogramBuilder.fillCountHistogram(algorithmName);
		}

	} else if (iErrorCode == 1) {
		cout<< coutPrefix << "trigger " << algorithmName << " does not exist in the L1 menu" << endl;
	} else {
		// error - see error code
		cout << coutPrefix << "Error Code " << iErrorCode << std::endl;
		cout << coutPrefix << "Algorithm name: " << algorithmName << std::endl;
	}
	return trigDecision;
}


// ------------ method called once each job just before starting event loop  ------------
void 
hoMuonAnalyzer::beginJob()
{
	if(debug){
		std::cout << coutPrefix << "Begin Job HoMuonAnalyzer." << std::endl;
	}
}

// ------------ method called once each job just after ending the event loop  ------------
void 
hoMuonAnalyzer::endJob() 
{
	if(debug){
		std::cout << coutPrefix << "End Job HoMuonAnalyzer." << std::endl;
	}
}

// ------------ method called when starting to processes a run  ------------

void 
hoMuonAnalyzer::beginRun(const edm::Run& iRun, 
		const edm::EventSetup& evSetup)
{

	bool useL1EventSetup = true;
	bool useL1GtTriggerMenuLite = true;
	std::cout << coutPrefix << "getL1GtRunCache" << std::endl;
	std::cout << coutPrefix << "UseL1EventSetup: " << useL1EventSetup << std::endl;
	std::cout << coutPrefix << "UseL1GtTriggerMenuLite: " << useL1GtTriggerMenuLite << std::endl;
	m_l1GtUtils.getL1GtRunCache(iRun, evSetup, useL1EventSetup, useL1GtTriggerMenuLite);

}


// ------------ method called when ending the processing of a run  ------------

void 
hoMuonAnalyzer::endRun(const edm::Run& iRun, const edm::EventSetup& evSetup)
{

	//Only interested in unique values
	listL1MuonPt.sort();
	listL1MuonPt.unique();  //NB it is called after sort
	if(debug){
		std::cout << coutPrefix << "The list contains " << listL1MuonPt.size() << "unique entries:";
		std::list<float>::iterator it;
		for (it=listL1MuonPt.begin(); it!=listL1MuonPt.end(); ++it){
			std::cout << ' ' << *it;
		}
		std::cout << endl;
	}
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
hoMuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

/*############################
 * Helper Functions
 *############################
 */
//Look for an L1 object in the direction of the HLT object
const l1extra::L1MuonParticle* hoMuonAnalyzer::getMatchedL1Object(trigger::TriggerObject hltObject
		,edm::Handle<l1extra::L1MuonParticleCollection> l1muons){
	double hltPhi,hltEta;
	double l1Phi,l1Eta;
	hltEta = hltObject.eta();
	hltPhi = hltObject.phi();
	for(unsigned int i = 0; i < l1muons->size(); i++){
		const l1extra::L1MuonParticle* l1muon = &(l1muons->at(i));
		l1Eta = l1muon->eta();
		l1Phi = l1muon->phi() + L1PHI_OFFSET;
		if(FilterPlugin::isInsideDeltaR(hltEta,l1Eta,hltPhi,l1Phi,deltaR_Max))
			return l1muon;
	}
	return NULL;
}

/**
 * Returns a pointer to the closest gen particle of all particles that are closer
 * than delta R < delta R max
 */
const reco::GenParticle* hoMuonAnalyzer::getBestGenMatch(float eta, float phi){
	const reco::GenParticle* bestGen = 0;
	float bestDR = 999.;
	reco::GenParticleCollection::const_iterator genIt = truthParticles->begin();
	reco::GenParticleCollection::const_iterator genEnd = truthParticles->end();
	for(; genIt!=genEnd; ++genIt) {
		if (abs(genIt->pdgId()) == 13 ) {
			float genPhi = genIt->phi();
			float genEta = genIt->eta();
			float dR = deltaR(eta,phi,genEta,genPhi);
			if (dR < deltaR_L1MuonMatching && dR < bestDR) { // CB get it from CFG
				bestDR = dR;
				bestGen = &(*genIt);
			}
		}
	}
	return bestGen;
}

/**
 * Returns a pointer to the closest gen particle of all particles that are closer
 * than delta R < delta R max
 *
 * Could do this and the gen function with a template
 *
 */
const pat::Muon* hoMuonAnalyzer::getBestPatMatch(float eta, float phi) {
	const pat::Muon* bestReco = 0;
	float bestDR = 999.;
	pat::MuonCollection::const_iterator recoIt = patMuons->begin();
	pat::MuonCollection::const_iterator recoEnd = patMuons->end();
	for (; recoIt != recoEnd; ++recoIt) {
		if (abs(recoIt->pdgId()) == 13) {
			float genPhi = recoIt->phi();
			float genEta = recoIt->eta();
			float dR = deltaR(eta, phi, genEta, genPhi);
			if (dR < deltaR_L1MuonMatching && dR < bestDR) { // CB get it from CFG
				bestDR = dR;
				bestReco = &(*recoIt);
			}
		}
	}
	return bestReco;
}


void hoMuonAnalyzer::defineTriggersOfInterest(){

	/*
	 * HLT Triggers
	 */

	string hltIsoMu24_key = "hltIsoMu24";
	hltNamesOfInterest.insert(pair<string, string>(hltIsoMu24_key,"HLT_IsoMu24_v18"));
	hltFiltersOfInterest.insert(pair<string, edm::InputTag>(hltIsoMu24_key,
			edm::InputTag("hltL3crIsoL1sMu16L1f0L2f16QL3"
					"f24QL3crIsoRhoFiltered0p15",
					"","HLT")));

	string hltMu5_key = "hltMu5";
	hltNamesOfInterest.insert(pair<string, string>(hltMu5_key, "HLT_Mu5_v21"));
	hltFiltersOfInterest.insert(pair<string, edm::InputTag>(hltMu5_key,
			edm::InputTag("hltL3fL1sMu3L3Filtered5",
					"","HLT")));

	string l1SingleMuOpen_key = "hlt_l1SingleMuOpen";
	hltNamesOfInterest.insert(std::pair<std::string,std::string>(l1SingleMuOpen_key,"HLT_L1SingleMuOpen_v7"));
	hltFiltersOfInterest.insert(std::pair<std::string,edm::InputTag>(l1SingleMuOpen_key,edm::InputTag("hltL1MuOpenL1Filtered0","","HLT")));

}

/**
 * Gets the TrackDetMatchInfo for a Gen Particle by using its vertex, momentum and charge information
 * Needs the magnetic field, the edm::event and the edm::setup of an event
 */
TrackDetMatchInfo* hoMuonAnalyzer::getTrackDetMatchInfo(reco::GenParticle genPart,const edm::Event& iEvent,
		const edm::EventSetup& iSetup){
	//Create the Track det match info
	GlobalPoint vertexPoint(genPart.vertex().X(),genPart.vertex().Y(),genPart.vertex().Z());
	GlobalVector mom (genPart.momentum().x(),genPart.momentum().y(),genPart.momentum().z());
	int charge = genPart.charge();
	const FreeTrajectoryState *freetrajectorystate_ = new FreeTrajectoryState(vertexPoint, mom ,charge , &(*theMagField));
	TrackDetMatchInfo* info = new TrackDetMatchInfo(assoc.associate(iEvent, iSetup, *freetrajectorystate_, assocParams));
	delete freetrajectorystate_;
	return info;
}

/**
 * Automatically call the histogram builder to fill the efficiency objects with different cut thresholds
 * Needs the pt for evaluation of the efficiency (ptMeasured) and the real pT
 */
void hoMuonAnalyzer::fillEfficiencyHistograms(double ptMeasured, double ptReal,std::string key){
	histogramBuilder.fillEfficiency(ptMeasured >= 5, ptReal, key + "Pt5");
	histogramBuilder.fillEfficiency(ptMeasured >= 10, ptReal, key + "Pt10");
	histogramBuilder.fillEfficiency(ptMeasured >= 15, ptReal, key + "Pt15");
	histogramBuilder.fillEfficiency(ptMeasured >= 20, ptReal, key + "Pt20");
	histogramBuilder.fillEfficiency(ptMeasured >= 25, ptReal, key + "Pt25");
}

/**
 * Prints all channel qualities for HO into histogram which itself is stored
 * in a root file
 */
void hoMuonAnalyzer::printChannelQualities(const edm::EventSetup& iSetup){
	std::cout << coutPrefix << "Printing Channel qualities" << std::endl;
	edm::ESHandle<HcalChannelQuality> p;
	iSetup.get<HcalChannelQualityRcd>().get(p);
	HcalChannelQuality *myqual = new HcalChannelQuality(*p.product());
	edm::ESHandle<HcalSeverityLevelComputer> mycomputer;
	iSetup.get<HcalSeverityLevelComputerRcd>().get(mycomputer);
	const HcalSeverityLevelComputer *mySeverity = mycomputer.product();
	int ieta, iphi;
	TFile* channelStatusfile = new TFile("channelStatus.root","RECREATE");
	TH2D* channelStatusHist = new TH2D("channelStatusHist","Channel status of HO",34,-16.5,16.5,73,-0.5,73.5);
	for (ieta=-15; ieta <= 15; ieta++) {
		if (ieta != 0) {
			for (iphi = 1; iphi <= 72; iphi++) {
				HcalDetId did(HcalOuter,ieta,iphi,4);
				const HcalChannelStatus *mystatus = myqual->getValues(did.rawId());
				channelStatusHist->SetBinContent(channelStatusHist->FindBin(ieta,iphi),mystatus->getValue());
				if (mySeverity->dropChannel(mystatus->getValue())) {
				}
			}
		}
	}
	channelStatusHist->Write();
	channelStatusfile->Write();
	channelStatusfile->Close();
}

/**
 * Analyzer function for events with no single muon trigger.
 * Tries to find l1 muons and match them to ho information
 */
void hoMuonAnalyzer::analyzeNoSingleMuEventsL1Loop(const edm::Event& iEvent,const edm::EventSetup& iSetup){
	for(unsigned int i = 0; i < l1Muons->size() ; i++){
		const l1extra::L1MuonParticle* l1Muon = &(l1Muons->at(i));
		edm::RefToBase<l1extra::L1MuonParticle> l1MuonCandiateRef(l1MuonView,i);
		reco::GenParticleRef ref = (*l1MuonGenMatches)[l1MuonCandiateRef];
		if(ref.isNonnull()){
			histogramBuilder.fillEfficiency(true,l1Muon->pt(),"L1GenRefNoSingleMu");
			//Once there is a gen ref, get the Track det match info
			TrackDetMatchInfo * muMatch = getTrackDetMatchInfo(*ref,iEvent,iSetup);
			double muMatchEta = muMatch->trkGlobPosAtHO.eta();
			double muMatchPhi = muMatch->trkGlobPosAtHO.phi();
			delete muMatch;
			histogramBuilder.fillCountHistogram("L1GenRefNoSingleMu");
			if(MuonHOAcceptance::inGeomAccept(muMatchEta,muMatchPhi)
			&& MuonHOAcceptance::inNotDeadGeom(muMatchEta,muMatchPhi)
			&& !hoMatcher->isInChimney(muMatchEta,muMatchPhi)){
				histogramBuilder.fillCountHistogram("L1GenRefNoSingleMuInGa");
				calculateGridMatchingEfficiency(&*l1Muon, ref->pt(),"L1GenRefNoSingleMuInGa");
			}

		} else {
			histogramBuilder.fillEfficiency(false,l1Muon->pt(),"L1GenRefNoSingleMu");
			histogramBuilder.fillEtaPhiGraph(l1Muon->eta(),l1Muon->phi() + L1PHI_OFFSET,"L1GenRefNoSingleMuFail");
		}
	}
}

/**
 * Analyzer function for events with no single muon trigger.
 * loops over gen muons and match them to ho information
 */
void hoMuonAnalyzer::analyzeNoSingleMuEventsGenLoop(const edm::Event& iEvent,const edm::EventSetup& iSetup){
	for(reco::GenParticleCollection::const_iterator genIt = truthParticles->begin();
			genIt != truthParticles->end(); genIt++){

		histogramBuilder.fillPtHistogram(genIt->pt(),"NoSingleMu");

		//Once there is a gen ref, get the Track det match info
		TrackDetMatchInfo * muMatch = getTrackDetMatchInfo(*genIt,iEvent,iSetup);
		double muMatchEta = muMatch->trkGlobPosAtHO.eta();
		double muMatchPhi = muMatch->trkGlobPosAtHO.phi();
		delete muMatch;
		histogramBuilder.fillCountHistogram("NoSingleMu");
		if(MuonHOAcceptance::inGeomAccept(muMatchEta,muMatchPhi)
		&& MuonHOAcceptance::inNotDeadGeom(muMatchEta,muMatchPhi)
		&& !hoMatcher->isInChimney(muMatchEta,muMatchPhi)){
			histogramBuilder.fillCountHistogram("NoSingleMuInGa");
			calculateGridMatchingEfficiency(genIt->eta(),genIt->phi(), genIt->pt(),"NoSingleMuInGa");
		}
	}
}

/**
 * Use this function to make the efficiency plots with root's TEfficiency.
 */
void hoMuonAnalyzer::analyzeWithGenLoop(const edm::Event& iEvent,const edm::EventSetup& iSetup){
	for(reco::GenParticleCollection::const_iterator genIt = truthParticles->begin();
			genIt != truthParticles->end(); genIt++){
		float genEta = genIt->eta();
		float genPhi = genIt->phi();
		const l1extra::L1MuonParticle* l1Part = 0;
		l1Part = functionsHandler->getBestL1MuonMatch(genEta,genPhi);
		histogramBuilder.fillCountHistogram("Gen");
		if(l1Part){

			double l1Eta, l1Phi;
			l1Eta = l1Part->eta();
			l1Phi = l1Part->phi() + L1PHI_OFFSET;

			fillEfficiencyHistograms(l1Part->pt(),genIt->pt(),"GenAndL1Muon");
			/**
			 * Fill Correlation between gen pt and l1 pt.
			 * Here, only the L1s best matching to the Gen are used. Above, Every L1 is
			 * used
			 */
			histogramBuilder.fillCorrelationHistogram(genIt->pt(),l1Part->pt(),"L1MuonPtGenLoop");
			histogramBuilder.fillCountHistogram("GenAndL1Muon");
			calculateGridMatchingEfficiency(&*l1Part,l1Part->pt(),"L1MuonTruth");
			fillGridMatchingQualityCodes(&*l1Part,genIt->pt(),"L1MuonTruth");
			fillAverageEnergyAroundL1Direction(l1Part,"L1MuonTruth");
			/**
			 * Find a rec hit that can be matched to the l1 particle. Use this information for the efficiency
			 * plots. This time it is ensured that only as many entries as there are gen particles is used
			 * This fixes double counting of ghost l1 muons
			 */
			const HORecHit* matchedRecHit = 0;
			matchedRecHit = hoMatcher->matchByEMaxInGrid(l1Eta,l1Phi,2);
			if(matchedRecHit){
				if(matchedRecHit->energy() > threshold){
					fillEfficiencyHistograms(l1Part->pt(),genIt->pt(),"GenAndL1MuonAndHoAboveThr");
					histogramBuilder.fillEnergyVsIEta(matchedRecHit->energy(),matchedRecHit->id().ieta(),"l1TruthAndHoMatch");
					histogramBuilder.fillCountHistogram("GenAndL1MuonAndHoAboveThr");
					double hoPhi = hoMatcher->getRecHitPhi(matchedRecHit);
					double hoIPhi = matchedRecHit->id().iphi();
					histogramBuilder.fillCorrelationGraph(hoPhi,l1Phi,"l1PhiVsHoPhi");
					histogramBuilder.fillCorrelationGraph(hoIPhi,l1Phi,"l1PhiVsHoIPhi");
					histogramBuilder.fillCorrelationGraph(hoIPhi,hoPhi,"hoPhiVsHoIPhi");

					TH2D* hist = new TH2D("hoTruthEnergyVsTime","HO Energy vs. Time;Time / ns;E_{Rec} / GeV",201,-100.5,100.5,2100, -5.0, 100.0);
					histogramBuilder.fillCorrelationHistogram(matchedRecHit->time(),matchedRecHit->energy(),"hoTruthEnergyVsTime",hist);
					delete hist;

					//Implement efficiency analysis for time window
					if(matchedRecHit->time() > -12.5 && matchedRecHit->time() < 12.5){
						//Might probably be removed. Fixed timing window won't work in data
					}
				}
			}
		}
	}
}

/**
 * Fill a histogram with the measured energy around a given L1.
 * For now the grid size is hard-coded
 */
void hoMuonAnalyzer::fillAverageEnergyAroundL1Direction(const l1extra::L1MuonParticle* l1Muon,std::string key){
	int gridSize = 5;

	double l1Eta, l1Phi;
	l1Eta = l1Muon->eta();
	l1Phi = l1Muon->phi() + L1PHI_OFFSET;

	for(auto recHitIt = hoRecoHits->begin(); recHitIt != hoRecoHits->end(); recHitIt++){
		if(hoMatcher->isRecHitInGrid(l1Eta, l1Phi,&*recHitIt,gridSize)){
			double hoEta = hoMatcher->getRecHitEta(&*recHitIt);
			double hoPhi = hoMatcher->getRecHitPhi(&*recHitIt);

			histogramBuilder.fillAverageEnergyHistograms(l1Eta,hoEta, l1Phi,hoPhi,recHitIt->energy(),"averageEnergyAroundPoint" + key);
			histogramBuilder.fillDeltaEtaDeltaPhiEnergyHistogram(l1Eta ,hoEta ,l1Phi ,hoPhi ,recHitIt->energy()
					,"averageEnergyAroundPoint" + key);//Use this function for the 1D distributions for each delta eta and delta phi

			double deltaPhi;
			deltaPhi = FilterPlugin::wrapCheck(l1Phi,hoMatcher->getRecHitPhi(&*recHitIt));
			double deltaEta = hoEta - l1Eta;

			double deltaPhiLocalTest;
			deltaPhiLocalTest = FilterPlugin::wrapCheck(l1Phi,hoPhi);
			histogramBuilder.fillGraph(deltaEta,deltaPhiLocalTest,"eAvAboveThrCounter" + key);
			if(hoMatcher->isRecHitInGrid(l1Eta, l1Phi,&*recHitIt,2))
				histogramBuilder.fillGraph(deltaEta,deltaPhiLocalTest,"eAvAboveThrCounter5x5" + key);
			if(recHitIt->energy() >= 0.2){
				if(hoMatcher->isRecHitInGrid(l1Eta, l1Phi,&*recHitIt,0)){
					histogramBuilder.fillGraph(deltaEta,deltaPhiLocalTest,"eAvAboveThrCentral" + key);
				}
				if (hoMatcher->isRecHitInGrid(l1Eta, l1Phi,&*recHitIt,1)) {
					histogramBuilder.fillGraph(deltaEta,deltaPhiLocalTest,"eAvAboveThr3x3" + key);
				}
				if (hoMatcher->isRecHitInGrid(l1Eta, l1Phi,&*recHitIt,2)) {
					histogramBuilder.fillGraph(deltaEta,deltaPhiLocalTest,"eAvAboveThr5x5" + key);
					if(l1Muon->pt() > 80){
						histogramBuilder.fillGraph(deltaEta,deltaPhiLocalTest,"eAvAboveThrPt805x5" + key);
					}
				}
			}

			TH1D* hist1D = new TH1D(("deltaPhi" + key).c_str(),"#Delta#phi;#Delta#phi;N Entries",81,-40*HoMatcher::HALF_HO_BIN/2. - HoMatcher::HALF_HO_BIN/4.
					,40*HoMatcher::HALF_HO_BIN/2. + HoMatcher::HALF_HO_BIN/4.);
			histogramBuilder.fillHistogram(deltaPhi,"deltaPhi" + key,hist1D);
			delete hist1D;

			TH1D* histL1Phi = new TH1D(("averageEnergyL1Phi" + key).c_str(),"L1 #phi;#phi;N Entries",145
					,-36*HoMatcher::HO_BIN - HoMatcher::HALF_HO_BIN/2.,36 * HoMatcher::HO_BIN + HoMatcher::HALF_HO_BIN/2.);
			histogramBuilder.fillHistogram(l1Phi,"averageEnergyL1Phi" + key,histL1Phi);
			delete histL1Phi;

			TH1D* histHoPhi = new TH1D(("averageEnergyHoPhi" + key).c_str(),"HO #phi;#phi;N Entries",145
					,-36*HoMatcher::HO_BIN - HoMatcher::HALF_HO_BIN/2.,36 * HoMatcher::HO_BIN + HoMatcher::HALF_HO_BIN/2.);
			histogramBuilder.fillHistogram(hoPhi,"averageEnergyHoPhi" + key,histHoPhi);
			delete histHoPhi;

		}
	}//for loop

	//Filling the average energy only for the highest energetic particle
	const HORecHit* matchedRecHit = 0;
	matchedRecHit = hoMatcher->matchByEMaxInGrid(l1Eta,l1Phi,5);
	if(matchedRecHit){
		histogramBuilder.fillDeltaEtaDeltaPhiHistogramsWithWeights(l1Eta
				,float(hoMatcher->getRecHitEta(matchedRecHit))	,l1Phi
				,float(hoMatcher->getRecHitPhi(matchedRecHit))	,matchedRecHit->energy()
				,("averageEMaxAroundPoint" + key).c_str());
		histogramBuilder.fillDeltaEtaDeltaPhiEnergyHistogram(l1Eta
				,float(hoMatcher->getRecHitEta(matchedRecHit))	,l1Phi
				,float(hoMatcher->getRecHitPhi(matchedRecHit))	,matchedRecHit->energy()
				,("averageEMaxAroundPoint" + key).c_str());
	}
}

/**
 * Fills an eta phi graph in case the gen particle is in the geometric
 * acceptance of HO.
 * Use a function for this to keep the code clean
 */
void hoMuonAnalyzer::fillHoGeomAcceptanceGraph(reco::GenParticle genPart){
	if(MuonHOAcceptance::inGeomAccept(genPart.eta(),genPart.phi())
				&& MuonHOAcceptance::inNotDeadGeom(genPart.eta(),genPart.phi())
				&& !hoMatcher->isInChimney(genPart.eta(),genPart.phi())){
		histogramBuilder.fillEtaPhiGraph(genPart.eta(),genPart.phi(),"HoGeomAcceptance");
	}
}

/**
 * Overloaded function to make the use for l1 objects easier
 */
void hoMuonAnalyzer::fillGridMatchingQualityCodes(const l1extra::L1MuonParticle* l1muon, float truePt, std::string key){

	double l1Eta, l1Phi;
	l1Eta = l1muon->eta();
	l1Phi = l1muon->phi() + L1PHI_OFFSET;

	//#####
	// Central tile
	//#####
	double variableBinArray[] = {0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,100,120,140,180};
	int l1MuonQuality = l1muon->gmtMuonCand().quality();
	histogramBuilder.fillMultiplicityHistogram(l1MuonQuality,key+"AllQualityCodes");
	if(hoMatcher->hasHoHitInGrid(l1Eta,l1Phi,0)){
		histogramBuilder.fillMultiplicityHistogram(l1MuonQuality,key + "QualityCodesCentral" );
		fillEfficiencyHistograms(l1muon->pt(),truePt,key + "GenPtCentral");
	} else{
		histogramBuilder.fillMultiplicityHistogram(l1MuonQuality,key + "QualityCodesCentralFail" );
		TH2D* hist = new TH2D((key + "pTvsQCCentralFail").c_str(),"p_{T} vs. QC (Central);QC;p_{T} / GeV",7,1.5,8.5,33,variableBinArray);
		histogramBuilder.fillCorrelationHistogram(l1MuonQuality,l1muon->pt(),key + "pTvsQCCentralFail",hist);
		delete hist;
	}
	//#####
	// 3 x 3
	//#####
	if(hoMatcher->hasHoHitInGrid(l1Eta,l1Phi,1)){
		histogramBuilder.fillMultiplicityHistogram(l1MuonQuality,key + "QualityCodes3x3");
		fillEfficiencyHistograms(l1muon->pt(),truePt,key + "GenPt3x3");
	} else {
		histogramBuilder.fillMultiplicityHistogram(l1MuonQuality,key + "QualityCodes3x3Fail");
		TH2D* hist = new TH2D((key + "pTvsQC3x3Fail").c_str(),"p_{T} vs. QC (3x3);QC;p_{T} / GeV",7,1.5,8.5,33,variableBinArray);
		histogramBuilder.fillCorrelationHistogram(l1MuonQuality,l1muon->pt(),key + "pTvsQC3x3Fail",hist);
		delete hist;
	}
	//#####
	// 5 x 5
	//#####
	if(hoMatcher->hasHoHitInGrid(l1Eta,l1Phi,2)){
		histogramBuilder.fillMultiplicityHistogram(l1MuonQuality,key + "QualityCodes5x5");
		fillEfficiencyHistograms(l1muon->pt(),truePt,key + "GenPt5x5");
	} else {
		histogramBuilder.fillMultiplicityHistogram(l1MuonQuality,key + "QualityCodes5x5Fail");
		TH2D* hist = new TH2D((key + "pTvsQC5x5Fail").c_str(),"p_{T} vs. QC (5x5);QC;p_{T} / GeV",7,1.5,8.5,33,variableBinArray);
		histogramBuilder.fillCorrelationHistogram(l1MuonQuality,l1muon->pt(),key + "pTvsQC5x5Fail",hist);
		delete hist;
	}
	histogramBuilder.fillCorrelationGraph(l1muon->pt(),l1MuonQuality,key);
}

/**
 * Overloaded function. Some additional stuff is done with the BX ID information from the L1 Object
 */
void hoMuonAnalyzer::calculateGridMatchingEfficiency(const l1extra::L1MuonParticle* l1muon, float pt, std::string key){

	/**
	 * Implement a quality code for the 3x3 matching
	 *
	 * 10 * Bx ID + L1 QC + HO time window information
	 *
	 */

	double l1Eta, l1Phi;
	l1Eta = l1muon->eta();
	l1Phi = l1muon->phi() + L1PHI_OFFSET;

	//Restrict L1 information to be within HO
	if(fabs(l1Eta) > MAX_ETA){
		return;
	}

	//Fill n multiple possbile HORecHits
	int nHits3x3 = hoMatcher->countHoHitsAboveThr(l1Eta,l1Phi,1);
	int nHits5x5 = hoMatcher->countHoHitsAboveThr(l1Eta,l1Phi,2);
	histogramBuilder.fillMultiplicityHistogram(nHits3x3,key + "_nHoHits3x3");
	histogramBuilder.fillMultiplicityHistogram(nHits5x5,key + "_nHoHits5x5");
	histogramBuilder.fillGraph(l1muon->pt(),nHits3x3,key + "_nHoHitsVsPt3x3");
	histogramBuilder.fillGraph(l1muon->pt(),nHits5x5,key + "_nHoHitsVsPt5x5");

	calculateGridMatchingEfficiency(l1Eta, l1Phi,pt, key);
	//Analyze the BX ID of L1 objects that do not have a match in the grid
	const HORecHit* recHit = hoMatcher->matchByEMaxInGrid(l1Eta,l1Phi,2);//getClosestRecHitInGrid(l1Eta,l1Phi,2);
	for(int i = 0; i < 3 ; i++){
		std::string gridString = CommonFunctionsHandler::getGridString(i);
		if(!recHit){
			histogramBuilder.fillBxIdVsPt(l1muon->bx(),l1muon->pt(),key + gridString + "Fail");
		} else {
			int qualityCode = 0;
			qualityCode = l1muon->gmtMuonCand().quality();
			qualityCode += l1muon->bx()*10;
			qualityCode += isInTimeWindow(recHit->time()) ? 0 : 100;
			qualityCode *= (recHit->time() > 0 ? 1 : -1);

			if(hoMatcher->isRecHitInGrid(l1Eta,l1Phi, recHit,i)){
				histogramBuilder.fillBxIdVsPt(l1muon->bx(),l1muon->pt(),key + gridString + "Match");
				histogramBuilder.fillQualityCodeVsPt(qualityCode,l1muon->pt(),key + gridString + "Match");
				int hoIEta = recHit->id().ieta();
				int hoIPhi = recHit->id().iphi();
				histogramBuilder.fillIEtaIPhiHistogram(hoIEta,hoIPhi,key + gridString);
				histogramBuilder.fillEtaPhiHistograms(hoMatcher->getRecHitEta(recHit),hoMatcher->getRecHitPhi(recHit),key + gridString);
			} else {
				histogramBuilder.fillBxIdVsPt(l1muon->bx(),l1muon->pt(),key + gridString + "Fail");
				histogramBuilder.fillQualityCodeVsPt(qualityCode,l1muon->pt(),key + gridString + "Fail");
			}
		}
	}
}

void hoMuonAnalyzer::fillTriggerRatesForQualityCodes(){
	//For L1 pt use the l1 binning
	//Normalize per bin width
	//Do this as stacked plot?
	//Use tight information only?
	//Quality codes as bit field
	// 0 - 7 	L1 Quality		-> Bit 0 - 2
	// 8, 16 	HO time window	-> Bit 3 - 4
	//32,64,xx 	L1 BX ID		-> Bit 5 - 11

	//MSB first
	//L1 Quality best 	xxx111
	//	Second best		xxx110
	//	third best		xxx101
	// ...
#define HO_time_good 11xxxb
#define HO_time_lo	01xxx
#define HO_time_hi	10xxx
#define L1BXID0		0001000xxxx....
#define L1BXID_P_1		0010000xxxx....
#define L1BXID_P_2		0100000xxxx....
#define L1BXID_P_3		1000000xxxx....
#define L1BXID_M_1		0000100xxxx....
#define L1BXID_M_2		0000010xxxx....
#define L1BXID_M_3		0000001xxxx....

	/**
	 * Start development by running over all L1 Objects.
	 * Also cases for no HO match can be depicted by having all 0 for the Ho time part.
	 * In future versions, more studies can be performed, using only the tight L1 information,
	 * or having the L1 algorithm result evaluated
	 */
	for(auto l1MuonIt = l1Muons->begin(); l1MuonIt != l1Muons->end(); l1MuonIt++){

	}
}

/**
 * Automatically fill efficiency and count histograms for the grid matching for grid sizes
 * central, 3x3 and 5x5. Also store the position information
 */
void hoMuonAnalyzer::calculateGridMatchingEfficiency(double eta, double phi, float pt, std::string key){
	const HORecHit* recHit = hoMatcher->matchByEMaxInGrid(eta,phi,2);//getClosestRecHitInGrid(eta,phi,2);
		for(int i = 0; i < 3 ; i++){
			if(!recHit){
				fillGridMatchingHistograms(false,i,pt,999,key,eta,phi);
			} else{
				fillGridMatchingHistograms(hoMatcher->isRecHitInGrid(eta,phi,recHit,i),i,pt,recHit->time(),key,eta,phi);
			}
		}
}

/**
 * This function automatically fills the corresponding histograms for the grid matching efficiency and the time window
 */
void hoMuonAnalyzer::fillGridMatchingHistograms(bool passed, int grid, double pt, double time, std::string key, double eta, double phi){
	std::string gridString = CommonFunctionsHandler::getGridString(grid);
	if(passed){
		histogramBuilder.fillCountHistogram(key + gridString);
		histogramBuilder.fillEfficiency(true,pt,key + gridString);
		histogramBuilder.fillEtaPhiGraph(eta,phi,key + gridString);
	} else{
		histogramBuilder.fillEfficiency(false,pt,key + gridString);
		histogramBuilder.fillEtaPhiGraph(eta,phi,key + gridString + "Fail");
	}
	if(passed && isInTimeWindow(time)){
		histogramBuilder.fillEfficiency(true,pt,key + "TimeWindow" + gridString);
	} else{
		histogramBuilder.fillEfficiency(false,pt,key + "TimeWindow" + gridString);
	}
}


void hoMuonAnalyzer::fillL1ResolutionPlots(const l1extra::L1MuonParticle* l1Part, const pat::Muon* patMuon,std::string label){
	float l1Eta = l1Part->eta();
	float l1Phi = l1Part->phi() + L1PHI_OFFSET;
	bool isTight = patMuon->isTightMuon(getPrimaryVertex());
	const HORecHit* matchedRecHit = 0;
	matchedRecHit = hoMatcher->matchByEMaxInGrid(l1Eta,l1Phi,2);
	histogramBuilder.fillL1ResolutionHistogram(l1Part->pt(), patMuon->pt(),label);
	if (matchedRecHit) {
		histogramBuilder.fillL1ResolutionHistogram(l1Part->pt(),patMuon->pt(),label +  "HoMatch");
	} else {
		histogramBuilder.fillL1ResolutionHistogram(l1Part->pt(),patMuon->pt(),label +  "NotHoMatch");
	}
	if(isTight){
		histogramBuilder.fillL1ResolutionHistogram(l1Part->pt(), patMuon->pt(), label + "Tight");
		if(matchedRecHit){
			histogramBuilder.fillL1ResolutionHistogram(l1Part->pt(), patMuon->pt(),label +  "TightHoMatch");
		} else {
			histogramBuilder.fillL1ResolutionHistogram(l1Part->pt(), patMuon->pt(),label +  "TightNotHoMatch");
		}
	}
}

/**
 * Create plots for pt resolution with x axis pt info coming from pat muons
 */
void hoMuonAnalyzer::analyzeL1Resolution(){
	for(auto patMuonIt = patMuons->begin(); patMuonIt != patMuons->end(); patMuonIt++){
		const l1extra::L1MuonParticle* l1Part = 0;
		l1Part = functionsHandler->getBestL1MuonMatch(patMuonIt->eta(),patMuonIt->phi());
		if(l1Part){
			if(fabs(l1Part->eta()) > MAX_ETA){
				continue;
			}
			fillL1ResolutionPlots(l1Part,&*patMuonIt,"patToL1Muon");
			const L1MuRegionalCand* dttfCand = findBestCandMatch(l1Part);
			if(dttfCand){
				if(dttfCand->isFineHalo()){
					histogramBuilder.fillGraph(patMuonIt->eta(),l1Part->eta(),"l1EtaVsPatEtaFine");
					if(patMuonIt->isTightMuon(getPrimaryVertex())){
						histogramBuilder.fillGraph(patMuonIt->eta(),l1Part->eta(),"l1EtaVsPatEtaFineTight");
					}
				} else{
					histogramBuilder.fillGraph(patMuonIt->eta(),l1Part->eta(),"l1EtaVsPatEtaNotFine");
					if(patMuonIt->isTightMuon(getPrimaryVertex())){
						histogramBuilder.fillGraph(patMuonIt->eta(),l1Part->eta(),"l1EtaVsPatEtaNotFineTight");
					}
				}
			}
		}
	}// End loop over pat

	for (auto l1Muon = l1Muons->begin(); l1Muon != l1Muons->end(); l1Muon++) {
		float l1Eta = l1Muon->eta();
		float l1Phi = l1Muon->phi() + L1PHI_OFFSET;
		if (fabs(l1Eta) > MAX_ETA) {
			continue;
		}
		const pat::Muon* patMuon = getBestPatMatch(l1Eta, l1Phi);
		if (patMuon) {
			fillL1ResolutionPlots(&*l1Muon,patMuon,"L1Muon");
		}
	}
}

void hoMuonAnalyzer::recoControlPlots(){
	for(auto recoIt = recoMuons->begin(); recoIt != recoMuons->end(); recoIt++){
		histogramBuilder.fillPtHistogram(recoIt->pt(),"recoMuons");
		histogramBuilder.fillEtaPhiGraph(recoIt->eta(), recoIt->phi(),"recoMuons");
		histogramBuilder.fillEtaPhiHistograms(recoIt->eta(),recoIt->phi(),"recoMuons");
	}
	int nTightMuons = 0;
	for(auto patMuonIt = patMuons->begin(); patMuonIt != patMuons->end(); ++patMuonIt){
		histogramBuilder.fillPtHistogram(patMuonIt->pt(),"patMuons");
		histogramBuilder.fillEtaPhiGraph(patMuonIt->eta(), patMuonIt->phi(),"patMuons");
		histogramBuilder.fillEtaPhiHistograms(patMuonIt->eta(),patMuonIt->phi(),"patMuons");
		if(patMuonIt->isTightMuon(getPrimaryVertex())){
			nTightMuons++;
			histogramBuilder.fillPtHistogram(patMuonIt->pt(),"patMuonsTight");
			histogramBuilder.fillEtaPhiGraph(patMuonIt->eta(), patMuonIt->phi(),"patMuonsTight");
			histogramBuilder.fillEtaPhiHistograms(patMuonIt->eta(),patMuonIt->phi(),"patMuonsTight");
			const HORecHit* recHit = hoMatcher->matchByEMaxInGrid(patMuonIt->eta(), patMuonIt->phi(),2);
			if(recHit){
				histogramBuilder.fillEtaPhiHistograms(hoMatcher->getRecHitEta(recHit),hoMatcher->getRecHitPhi(recHit),"patMuonsTight_HO");
				histogramBuilder.fillEnergyVsIEta(recHit->energy(),recHit->id().ieta(),"patMuonsTight");
			}
		}
	}
	//Count number of pat muons and tight pat muons in an event
	histogramBuilder.fillMultiplicityHistogram(patMuons->size(),"patMuonsSize");
	histogramBuilder.fillMultiplicityHistogram(nTightMuons,"tightPatMuonsSize");
}

/**
 * Try to find the best match in dt regional candidates for a given
 * L1 Muon particle
 */
const L1MuRegionalCand* hoMuonAnalyzer::findBestCandMatch(const l1extra::L1MuonParticle* l1Muon){
	const L1MuRegionalCand* dtCand = 0;
	float bestDeltaR = 999;
	float bestQuality = 0;

	// get GMT readout collection
	L1MuGMTReadoutCollection const* gmtrc = pCollection.product();
	std::vector<L1MuGMTReadoutRecord> gmt_records = gmtrc->getRecords();
	std::vector<L1MuGMTReadoutRecord>::const_iterator RRItr;

	for ( RRItr = gmt_records.begin(); RRItr != gmt_records.end(); ++RRItr ) {
		std::vector<L1MuRegionalCand> dttfCands = RRItr->getDTBXCands();
		std::vector<L1MuRegionalCand>::iterator dttfCand;
		for( dttfCand = dttfCands.begin(); dttfCand != dttfCands.end();	++dttfCand ) {
			if(dttfCand->empty()) continue;
			//phiValue is lower edge of bin -> + half L1Bin
			float dR = deltaR(l1Muon->eta(),l1Muon->phi()+ L1PHI_OFFSET,dttfCand->etaValue(),dttfCand->phiValue() + L1PHI_OFFSET);
			if(dR < 0.5){
				if(dR < bestDeltaR){
					if(dttfCand->quality() > bestQuality){
						dtCand = &*dttfCand;
						bestDeltaR = dR;
					}
				}
			}
		}// each dttf cand
	}
	return dtCand;
}

void hoMuonAnalyzer::fillTimingHistograms(const l1extra::L1MuonParticle* l1Muon, const HORecHit* hoRecHit, bool isTight, std::string extraId){
	std::string nameTrunk = "timingSupport_";
	nameTrunk += extraId;
	if(isTight){
		nameTrunk += "tight_";
	}
	histogramBuilder.fillCountHistogram(nameTrunk);
	double hoTime;
	if(hoRecHit){
		hoTime = hoRecHit->time();
	} else {
		hoTime = -999;
	}

	/**
	 * Fill graph containing eta and bxID information of L1Muons
	 */
	histogramBuilder.fillGraph(l1Muon->eta(),l1Muon->bx(),nameTrunk + "bxidVsEta");

	//No test for non-null needed. Starting from L1 information so there should be a regional
	//candidate somewhere
	const L1MuRegionalCand* l1RegCand = findBestCandMatch(l1Muon);

	switch(l1Muon->gmtMuonCand().quality()){
	case 7:
		//Matched DT-RPC
		if(hoTime == -999){
			histogramBuilder.fillCountHistogram(nameTrunk + "MatchedDtRpc");
			histogramBuilder.fillBxIdHistogram(l1Muon->bx(),nameTrunk + "MatchedDtRpc");
		}
		else{
			histogramBuilder.fillCountHistogram(nameTrunk + "MatchedDtRpcHo");
			histogramBuilder.fillBxIdHistogram(l1Muon->bx(),nameTrunk + "MatchedDtRpcHo");
			histogramBuilder.fillDeltaTimeHistogram(hoTime,l1Muon->bx(),nameTrunk + "MatchedDtRpcHo");
			histogramBuilder.fillTimeHistogram(hoTime,nameTrunk + "MatchedDtRpcHo");
			histogramBuilder.fillGraph(hoRecHit->id().ieta(),hoRecHit->time(),nameTrunk + "MatchedDtRpcHoTimeGraph");
			histogramBuilder.fillGraphPerIeta(hoRecHit->id().ieta(),hoRecHit->id().iphi(),hoRecHit->time(),nameTrunk + "MatchedDtRpcHoIphiTime");

		}
		break;
	case 6:
		//Do a crosscheck on the detector index and the quality code
		histogramBuilder.fillMultiplicityHistogram(l1Muon->gmtMuonCand().detector(),nameTrunk + "_detectorIndexInUnmatchedDt");
		histogramBuilder.fillGraph(l1Muon->eta(),l1Muon->bx(),nameTrunk + "dtOnly_bxidVsEta");
		if(hoTime == -999){
			histogramBuilder.fillCountHistogram(nameTrunk + "UnmatchedDt");
			histogramBuilder.fillBxIdHistogram(l1Muon->bx(),nameTrunk + "UnmatchedDt");
			histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDt");
			if(l1RegCand->isFineHalo()){
				histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDtFine");
			} else {
				histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDtNotFine");
			}
			if(l1Muon->bx() != 0){
				histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDtBxNot0");
				if(l1RegCand->isFineHalo()){
					histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDtBxNot0Fine");
				} else{
					histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDtBxNot0NotFine");
				}
			}
		}
		else{
			histogramBuilder.fillCountHistogram(nameTrunk + "UnmatchedDtHo");
			histogramBuilder.fillBxIdHistogram(l1Muon->bx(),nameTrunk + "UnmatchedDtHo");
			histogramBuilder.fillDeltaTimeHistogram(hoTime,l1Muon->bx(),nameTrunk + "UnmatchedDtHo");
			histogramBuilder.fillTimeHistogram(hoTime,nameTrunk + "UnmatchedDtHo");
			histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDtHo");
			if(l1RegCand->isFineHalo()){
				histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDtHoFine");
			} else {
				histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDtHoNotFine");
			}
			histogramBuilder.fillGraph(hoRecHit->id().ieta(),hoRecHit->time(),nameTrunk + "UnmatchedDtHoTimeGraph");
			histogramBuilder.fillGraphPerIeta(hoRecHit->id().ieta(),hoRecHit->id().iphi(),hoRecHit->time(),nameTrunk + "UnmatchedDtHoIphiTime");
			if(l1Muon->bx() != 0){
				histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDtHoBxNot0");
				histogramBuilder.fillGraph(hoRecHit->id().ieta(),hoRecHit->time(),nameTrunk + "UnmatchedDtHoBxNot0TimeGraph");
				histogramBuilder.fillGraphPerIeta(hoRecHit->id().ieta(),hoRecHit->id().iphi(),hoRecHit->time(),nameTrunk + "UnmatchedDtHoBxNot0IphiTime");

				if(l1RegCand->isFineHalo()){
					histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDtHoBxNot0Fine");
				} else {
					histogramBuilder.fillEtaPhiGraph(l1Muon->eta(), l1Muon->phi() + L1PHI_OFFSET,nameTrunk + "UnmatchedDtHoBxNot0NotFine");
				}
			}
		}
		//Unmatched DT
		break;
	default:
		if(hoTime == -999){
			histogramBuilder.fillCountHistogram(nameTrunk + "OtherCodes");
			histogramBuilder.fillBxIdHistogram(l1Muon->bx(),nameTrunk + "OtherCodes");
		}
		else{
			histogramBuilder.fillCountHistogram(nameTrunk + "OtherCodesHo");
			histogramBuilder.fillBxIdHistogram(l1Muon->bx(),nameTrunk + "OtherCodesHo");
			histogramBuilder.fillDeltaTimeHistogram(hoTime,l1Muon->bx(),nameTrunk + "OtherCodesHo");
			histogramBuilder.fillTimeHistogram(hoTime,nameTrunk + "OtherCodesHo");
		}
		//Other codes
		break;
	}
}


/**
 * Study the timing information when RPC info is not available
 */
void hoMuonAnalyzer::analyzeTimingSupport(){
	//switch for later use
	bool isDtOnlyEvent = false;
	for(auto l1Muon = l1Muons->begin(); l1Muon != l1Muons->end(); l1Muon++){
		float l1Eta = l1Muon->eta();
		float l1Phi = l1Muon->phi() + L1PHI_OFFSET;
		if( fabs(l1Eta) > MAX_ETA){
			continue;
		}
		//Count the number of L1 objects in a dt only l1 muon event
		//But only count once!
		if(!isDtOnlyEvent){
			if(l1Muon->gmtMuonCand().quality() == 6){
				histogramBuilder.fillMultiplicityHistogram(l1Muons->size(),"nL1InDtOnly");
				isDtOnlyEvent = true;
			}
		}
		const pat::Muon* patMuon = getBestPatMatch(l1Eta,l1Phi);
		if(patMuon){

//			const reco::Vertex primVertex = getPrimaryVertex();
//			std::cout << "DEBUG #############################################################" << std::endl;
//			std::cout << primVertex.x() << std::endl;
//			std::cout << primVertex.y() << std::endl;
//			std::cout << primVertex.z() << std::endl;
//						std::cout << "GlobalMuon: " << (patMuon->isGlobalMuon() ? "OK" : "\033[91mFAIL\033[0m") << std::endl;
//			std::cout << "PFMuon: " << (patMuon->isPFMuon() ? "OK" : "\033[91mFAIL\033[0m") << std::endl;
//			if (patMuon->globalTrack().isNonnull()){
//				std::cout << "Normalized Chi^2/ndof: " << (patMuon->globalTrack()->normalizedChi2() < 10? "OK" : "\033[91mFAIL\033[0m") << std::endl;
//				std::cout << "n muon hits: " << (patMuon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 ? "OK" : "\033[91mFAIL\033[0m") << std::endl;
//				std::cout << "Matched Stations: " << (patMuon->numberOfMatchedStations() > 1 ? "OK" : "\033[91mFAIL\033[0m") << std::endl;
//				std::cout << "dxy : " << ( fabs(patMuon->muonBestTrack()->dxy(primVertex.position())) < 0.2  ? "OK" : "\033[91mFAIL\033[0m") << std::endl;
//				std::cout << "dz : " << ( fabs(patMuon->muonBestTrack()->dz(primVertex.position())) < 0.5 ? "OK" : "\033[91mFAIL\033[0m") << std::endl;
//				std::cout << "pixel hits: " << (patMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 ? "OK" : "\033[91mFAIL\033[0m") << std::endl;
//				std::cout << "tracker layers : " << ( patMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 ? "OK" : "\033[91mFAIL\033[0m") << std::endl;
//				std::cout << "dz value: " << patMuon->muonBestTrack()->dz(primVertex.position()) << std::endl;
//			}
//			std::cout << "DEBUG #############################################################" << std::endl;

			bool isTight = patMuon->isTightMuon(getPrimaryVertex());
			const HORecHit* hoRecHit = hoMatcher->matchByEMaxInGrid(l1Eta,l1Phi,1);
			fillTimingHistograms(&*l1Muon,hoRecHit,false);
			if(isTight){
				fillTimingHistograms(&*l1Muon,hoRecHit,true);
			}
			/**
			 * Fill histograms with L1 pT cut
			 */
			for(int i = 10; i<=25; i+=5){
				if(l1Muon->pt() <= i){
					break;
				}
				bool isTight = patMuon->isTightMuon(getPrimaryVertex());
				const HORecHit* hoRecHit = hoMatcher->matchByEMaxInGrid(l1Eta,l1Phi,1);
				fillTimingHistograms(&*l1Muon,hoRecHit,false,Form("pt%d_",i));
				if(isTight){
					fillTimingHistograms(&*l1Muon,hoRecHit,true,Form("pt%d_",i));
				}
			}
		}
	}
}

void hoMuonAnalyzer::analyzeHoTriggerPrimitives(){
	histogramBuilder.fillMultiplicityHistogram(hoTPDigis->size(),"hoTpSize");
	for( auto tpIt = hoTPDigis->begin(); tpIt != hoTPDigis->end(); tpIt++){
		histogramBuilder.fillIEtaIPhiHistogram(tpIt->ieta(),tpIt->iphi(),"hoTpCollection");
	}
}

void hoMuonAnalyzer::analyzeGridMatching(){
	for(auto l1Muon = l1Muons->begin(); l1Muon != l1Muons->end(); l1Muon++){
		float l1Eta = l1Muon->eta();
		float l1Phi = l1Muon->phi() + L1PHI_OFFSET;
		if(fabs(l1Eta) > MAX_ETA){
			continue;
		}
		const pat::Muon* patMuon = getBestPatMatch(l1Eta,l1Phi);
		if(patMuon){
			calculateGridMatchingEfficiency(&*l1Muon,patMuon->pt(),"gridMatching_loose");
			calculateGridMatchingEfficiency(&*l1Muon,l1Muon->pt(),"gridMatching_L1pT_loose");
			if(patMuon->isTightMuon(getPrimaryVertex())){
				calculateGridMatchingEfficiency(&*l1Muon,patMuon->pt(),"gridMatching_tight");
				calculateGridMatchingEfficiency(&*l1Muon,l1Muon->pt(),"gridMatching_L1pT_tight");
			}
		}
	}

	/**
	 * Is probably not the correct order of matching. Should start from L1
	 * Can be removed
	 * 12. 05. 2016: Keep this order for clean signal L1!
	 */
	for(auto patMuonIt = patMuons->begin(); patMuonIt != patMuons->end(); ++patMuonIt){
		const l1extra::L1MuonParticle* l1Part = 0;
		l1Part = functionsHandler->getBestL1MuonMatch(patMuonIt->eta(),patMuonIt->phi());
		if(l1Part){
			float l1Eta = l1Part->eta();
			//FIXME: Need to perform a wrapcheck on statements like this one
			float l1Phi = l1Part->phi() + L1PHI_OFFSET;
			//Restrict the L1 information to Ho range
			if(fabs(l1Eta) > MAX_ETA){
				continue;
			}
			histogramBuilder.fillEtaPhiGraph(l1Eta, l1Phi,"patToL1Muons");
			histogramBuilder.fillCountHistogram("patToL1Muons");
			calculateGridMatchingEfficiency(&*l1Part,l1Part->pt(),"patToL1Muons_L1pT");
			calculateGridMatchingEfficiency(&*l1Part,patMuonIt->pt(),"patToL1Muons");
			fillAverageEnergyAroundL1Direction(&*l1Part,"patToL1Muons");
			if(patMuonIt->isTightMuon(getPrimaryVertex())){
				histogramBuilder.fillEtaPhiGraph(l1Eta, l1Phi,"patTightToL1Muons");
				histogramBuilder.fillCountHistogram("patTightToL1Muons");
				calculateGridMatchingEfficiency(&*l1Part,l1Part->pt(),"patTightToL1Muons_L1pT");
				calculateGridMatchingEfficiency(&*l1Part,patMuonIt->pt(),"patTightToL1Muons");
				fillAverageEnergyAroundL1Direction(&*l1Part,"patTightToL1Muons");
			}
		}
	}
}

/**
 * Call all function that process information coming from RECO
 */
void hoMuonAnalyzer::processRecoInformation(const edm::Event& iEvent, const edm::EventSetup& iSetup){
	analyzeL1Resolution();
	recoControlPlots();
	analyzeGridMatching();
}

void hoMuonAnalyzer::processGenInformation(const edm::Event& iEvent,const edm::EventSetup& iSetup){
	//Some Gen stuff
		std::string gen_key = "gen";
		int genMuonCounter = 0;
		for(reco::GenParticleCollection::const_iterator genIt = truthParticles->begin();
				genIt != truthParticles->end(); genIt++){
			//Check for muons in Full barrel only
			if( abs(genIt->pdgId()) == 13 ){
				genMuonCounter++;
				histogramBuilder.fillPtHistogram(genIt->pt(),gen_key);
				histogramBuilder.fillEtaPhiGraph(genIt->eta(),genIt->phi(),gen_key);
				//For Ho geometric acceptance
				fillHoGeomAcceptanceGraph(*genIt);
				for (int i = 0; i < 200; i+=2) {
					if(genIt->pt() >= i){
						histogramBuilder.fillTrigRateHistograms(i,gen_key);
					}
				}
			}
		}
		histogramBuilder.fillMultiplicityHistogram(genMuonCounter,gen_key);
		analyzeWithGenLoop(iEvent,iSetup);

		/*
		 * Level 1 Muons
		 */
		string l1muon_key = "L1Muon";

		/**
		 * first loop over all L1 Muon objects. The contents of this Loop
		 * may be moved to the larger loop over l1 objects later in the code
		 */
		int successfulMatches = 0;
		int failedMatches = 0;
		//Define iterators
		l1extra::L1MuonParticleCollection::const_iterator bl1Muon = l1Muons->begin();
		l1extra::L1MuonParticleCollection::const_iterator el1Muon = l1Muons->end();
		for( unsigned int i = 0 ; i < l1Muons->size(); i++  ) {
			const l1extra::L1MuonParticle* bl1Muon = &(l1Muons->at(i));
			float l1Eta = bl1Muon->eta();
			float l1Phi = bl1Muon->phi() + L1PHI_OFFSET;
			const reco::GenParticle* genMatch = getBestGenMatch(l1Eta, l1Phi);
			if(genMatch){
				successfulMatches++;
				histogramBuilder.fillDeltaVzHistogam( (genMatch->vz() - bl1Muon->vz()) ,l1muon_key);
				histogramBuilder.fillCorrelationHistogram(genMatch->pt(),bl1Muon->pt(),"L1MuonPt");
				histogramBuilder.fillEtaPhiGraph(genMatch->eta(),genMatch->phi(),"L1ToGen");
				histogramBuilder.fillEtaPhiPtHistogram(genMatch->eta(), genMatch->phi(),genMatch->pt(),"L1ToGen");
				fillEfficiencyHistograms(bl1Muon->pt(),genMatch->pt(),"L1Muon");
				if(bl1Muon->gmtMuonCand().detector() == 2 /*DT only*/){
					histogramBuilder.fillBxIdHistogram(bl1Muon->bx(),"BxDtOnly");
				}
				if(bl1Muon->bx() != 0){
					histogramBuilder.fillPtHistogram(genMatch->pt(),"BxWrongGen");
					histogramBuilder.fillEtaPhiGraph(genMatch->eta(),genMatch->phi(),"BxWrongGen");
					histogramBuilder.fillEtaPhiPtHistogram(genMatch->eta(), genMatch->phi(),genMatch->pt(),"BxWrongGen");
					/**
					 * Fill a multiplicity histogram with the detector index of the underlying GMT Cand
					 * From L1MuGmtExtendedCand:
					 * 1 RPC, 2 DT, 3 DT/RPC, 4 CSC, 5 CSC/RPC
					 *
					 * if (quality() == 7) // matched ?
						return isFwd() ? 5 : 3;
					   else
						return isRPC() ? 1 : ( isFwd()? 4 : 2);
					 */
					histogramBuilder.fillMultiplicityHistogram(bl1Muon->gmtMuonCand().detector(),"detectorIndexBxWrong");
				} else {
					histogramBuilder.fillPtHistogram(genMatch->pt(),"BxRightGen");
					histogramBuilder.fillEtaPhiGraph(genMatch->eta(),genMatch->phi(),"BxRightGen");
					histogramBuilder.fillEtaPhiPtHistogram(genMatch->eta(), genMatch->phi(),genMatch->pt(),"BxRightGen");
					histogramBuilder.fillMultiplicityHistogram(bl1Muon->gmtMuonCand().detector(),"detectorIndexBxRight");
				}
				/* Built this to fix the strange behavior of the efficiency plots.
				 * Did not yet help completely. The reason for the strange behavior is probably the fact,
				 * that there may be more than one l1 muons that can be matched to the Gen particle
				 */
				const HORecHit* matchedRecHit = hoMatcher->matchByEMaxInGrid(l1Eta, l1Phi,2);
				if(matchedRecHit){
					if(matchedRecHit->energy() > threshold)
						fillEfficiencyHistograms(bl1Muon->pt(),genMatch->pt(),"L1MuonAndHoAboveThr");
				}
			} else{
				failedMatches++;
			}
			edm::RefToBase<l1extra::L1MuonParticle> l1MuonCandiateRef(l1MuonView,i);
			reco::GenParticleRef ref = (*l1MuonGenMatches)[l1MuonCandiateRef];
			if(ref.isNonnull())
				histogramBuilder.fillPdgIdHistogram(ref->pdgId(),l1muon_key);
			else
				histogramBuilder.fillPdgIdHistogram(0,l1muon_key);
		}
		/**
		 * Count how often the matches to Gen were successful and how often they failed per Event.
		 * Also Count how often all matchings failed
		 */
		histogramBuilder.fillMultiplicityHistogram(failedMatches,"failedL1ToGenMatches");
		histogramBuilder.fillMultiplicityHistogram(successfulMatches,"successfulL1ToGenMatches");
		if(failedMatches && !successfulMatches){
			histogramBuilder.fillCountHistogram("allL1ToGenFailed");
		}
}

/**
 * Get an artificial primary vertex.
 * Trying to use this in the PU sample. Might fix the low tight muon output.
 */
const reco::Vertex hoMuonAnalyzer::getArtificialPrimaryVertex(){
	// =================================================================================
	// Look for the Primary Vertex (and use the BeamSpot instead, if you can't find it):
	reco::Vertex::Point posVtx;
	reco::Vertex::Error errVtx;

	reco::BeamSpot bs = *recoBeamSpotHandle;

	posVtx = bs.position();
	errVtx(0,0) = bs.BeamWidthX();
	errVtx(1,1) = bs.BeamWidthY();
	errVtx(2,2) = bs.sigmaZ();

	const reco::Vertex vtx(posVtx,errVtx);
	return vtx;
}

/**
 * Get the primary vertex for the current event
 * Used to identify tight muons
 */
const reco::Vertex hoMuonAnalyzer::getPrimaryVertex(){
	// =================================================================================
	// Look for the Primary Vertex (and use the BeamSpot instead, if you can't find it):
	reco::Vertex::Point posVtx;
	reco::Vertex::Error errVtx;
	reco::Vertex vtx;
	unsigned int theIndexOfThePrimaryVertex = 999.;

	if (vertexColl.isValid()){
		for (unsigned int ind=0; ind<vertexColl->size(); ++ind) {
			if ( (*vertexColl)[ind].isValid() && !((*vertexColl)[ind].isFake()) ) {
				theIndexOfThePrimaryVertex = ind;
				break;
			}
		}
	}

	if ((theIndexOfThePrimaryVertex<100) && !useArtificialPrimaryVertex) {
		posVtx = ((*vertexColl)[theIndexOfThePrimaryVertex]).position();
		errVtx = ((*vertexColl)[theIndexOfThePrimaryVertex]).error();
		//vtx = reco::Vertex(posVtx,errVtx);
		vtx = (reco::Vertex(posVtx,errVtx));
	}
	else {
		vtx = getArtificialPrimaryVertex();
	}
	return vtx;
}

/**
 * Make threshold scan for the number of passing rec hits for different energy thresholds
 */
void hoMuonAnalyzer::makeHoRecHitThresholdScan(){
	int iterationCounter = 0;
	for(double eThr = 0.025; eThr <= 5; eThr += 0.025){
		int thrCounter = 0;
		for(auto recHitIt = hoRecoHits->begin(); recHitIt != hoRecoHits->end(); ++recHitIt){
			if(recHitIt->energy() >= eThr){
				thrCounter++;
			}
		}
		histogramBuilder.fillMultiplicityHistogram(thrCounter,Form("recHitThrScan%d",iterationCounter));
		iterationCounter++;
	}
}

/**
 * Try to determine the origin of the additional L1 Objects in the PU sample
 */
void hoMuonAnalyzer::analyzeL1Sources(){
	for(auto patIt = patMuons->begin(); patIt != patMuons->end(); patIt++){
		const reco::Candidate* mother = patIt->mother(0);
		if(mother){
			std::cout << "Mother Id: " << mother->pdgId() << std::endl;
		}
	}
}

/**
 * Do the filling of the Energy histograms and at the same time analyze the efficiencies for matching L1 to reco,
 * to HO, and so on...
 */
void hoMuonAnalyzer::analyzeEnergyDeposit(const edm::Event& iEvent,const edm::EventSetup& iSetup){
	for (auto l1It = l1Muons->begin(); l1It != l1Muons->end();l1It++){
		float l1Eta = l1It->eta();
		float l1Phi = l1It->phi() + L1PHI_OFFSET;
		if( fabs(l1Eta) > MAX_ETA ){
			continue;
		}
		histogramBuilder.fillCountHistogram("energyDeposit_L1");
		const pat::Muon* patMuon = getBestPatMatch(l1Eta,l1Phi);
		if(patMuon){
			histogramBuilder.fillCountHistogram("energyDeposit_L1Reco");
			const HORecHit* hoRecHit = hoMatcher->matchByEMaxInGrid(l1Eta,l1Phi,1);
			if(hoRecHit){
				histogramBuilder.fillCountHistogram("energyDeposit_L1RecoHo");
				histogramBuilder.fillEnergyVsIEta(hoRecHit->energy(),hoRecHit->id().ieta(),"energyDeposit_L1RecoHo");
				if(patMuon->isTightMuon(getPrimaryVertex())){
					histogramBuilder.fillCountHistogram("energyDeposit_L1RecoHoTight");
					histogramBuilder.fillEnergyVsIEta(hoRecHit->energy(),hoRecHit->id().ieta(),"energyDeposit_L1RecoHoTight");
				}
			}
			hoRecHit = 0;
			/**
			 * Do inverted cut order as well
			 */
			if(patMuon->isTightMuon(getPrimaryVertex())){
				histogramBuilder.fillCountHistogram("energyDeposit_L1RecoTight");
				const HORecHit* hoRecHit = hoMatcher->matchByEMaxInGrid(l1Eta,l1Phi,1);
				if(hoRecHit){
					histogramBuilder.fillCountHistogram("energyDeposit_L1RecoTightHo");
					histogramBuilder.fillEnergyVsIEta(hoRecHit->energy(),hoRecHit->id().ieta(),"energyDeposit_L1RecoTightHo");
				}
			}
			/**
			 * Do the GA acceptance test with reco information and no threshold
			 */
			const HORecHit* hoRecHitNoThr = hoMatcher->matchByEMaxInGrid(l1Eta,l1Phi,1,true);
			if(hoRecHitNoThr){
				histogramBuilder.fillCountHistogram("energyDeposit_L1RecoHoNoThr");
				histogramBuilder.fillEnergyVsIEta(hoRecHitNoThr->energy(),hoRecHitNoThr->id().ieta(),"energyDeposit_L1RecoHoNoThr");
			}
			hoRecHitNoThr = 0;
			if(MuonHOAcceptance::inGeomAccept(patMuon->eta(),patMuon->phi())){
				histogramBuilder.fillCountHistogram("energyDeposit_L1RecoGa");
				const HORecHit* hoRecHitNoThr = hoMatcher->matchByEMaxInGrid(l1Eta,l1Phi,1,true);
				if(hoRecHitNoThr){
					histogramBuilder.fillCountHistogram("energyDeposit_L1RecoGaHoNoThr");
					histogramBuilder.fillEnergyVsIEta(hoRecHitNoThr->energy(),hoRecHitNoThr->id().ieta(),"energyDeposit_L1RecoGaHoNoThr");
				}
			}
		}
	}
}

//define this as a plug-in
DEFINE_FWK_MODULE(hoMuonAnalyzer);

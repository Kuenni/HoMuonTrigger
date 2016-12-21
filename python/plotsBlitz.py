#!/usr/bin/env python
import argparse
from types import FunctionType

parser = argparse.ArgumentParser()
parser.add_argument('--scripts','-s', dest='scripts', nargs='+',
                   help='The subscript(s) to be called for plotting')

parser.add_argument('--data','-d'
					,dest='data'
					,action="store_true",default=False
					,help='Skip anything that has to do with GEN')

parser.add_argument('--list','-l'
					,dest='list'
					,action="store_true",default=False
					,help='Show list of available scripts')

parser.add_argument('--fileScheme','-f'
					,dest='source'
					,help='File name scheme for the source of events')

parser.add_argument('--moduleName','-m'
					,dest='moduleName'
					,help='Set a different CMSSW module name. Necessary if the analyzer plugin was not '
						'hoMuonAnalyzer')

parser.add_argument('--DEBUG'
					,dest='DEBUG',action="store_true",default=False
					,help='Enable more verbose output in modules.')

parser.add_argument('--non-interactive'
					,dest='nonInteractive',action="store_true",default=False
					,help='Do not ask for user input to end the script. Plots will close automatically.')

parser.add_argument('--all','-a'
					,dest='plotAll',action="store_true",default=False
					,help='Produces all plots that are available. Implies non interactive mode.')

args = parser.parse_args()

if args.DEBUG:
	print "Running in DEBUG mode"

if args.plotAll:
	args.nonInteractive = True
	import sys
	#Need to do this before importing pyROOT
	#to enable batch mode
	sys.argv.append( '-b' )
	print "Running in batch mode"

def updateModuleName(lib):
	if args.moduleName:
		lib.setModuleName(args.moduleName)

def checkUserInput():
	if not args.nonInteractive:
		raw_input('Continue with <Enter>')
	return

from comparison.Comparison import Comparison
from dataQuality.ControlPlots import ControlPlots
from dataQuality.Energy	import Energy
from dataQuality.EvsEtaPhi import EvsEtaPhi
from phishift.DeltaPhi import DeltaPhi
from efficiency.Counters import Counters
from efficiency.GridMatching import GridMatching
from efficiency.HoThresholdScan import HoThresholdScan 
from efficiency.TimeWindow import TimeWindow
from dataQuality.PtResolution import PtResolution
from dataQuality.QualityCode import QualityCode
from dataQuality.Timing import Timing
from dataQuality.DtOnlyCoordinates import DtOnlyCoordinates
from dataQuality.HoTimeVsEta import HoTimeVsEta
from dataQuality.DtRpc import DtRpc
from dataQuality.DttfCands import DttfCands
from dataQuality.NoL1Muon import NoL1Muon
from dataQuality.Pileup import Pileup

def plotPileup():
	lib = Pileup(filename=args.source,data=args.data,debug = args.DEBUG)
	res = lib.plotPileup()
	checkUserInput()
	return

def plotNoL1Muons():
	lib = NoL1Muon(filename=args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	resEmax = lib.plotEMaxNoL1Muon()
	lib.printNoL1Info()
	checkUserInput()
	return

def plotDttfCands():
	lib = DttfCands(filename=args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	resDttf = lib.plotDttfCands()
	resDttfFine = lib.plotDttfCandsFine()
	resDttfNotFine = lib.plotDttfCandsNotFine()
	checkUserInput()
	return

def plotDtRpc():
	lib = DtRpc(filename=args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	resDtRpc = lib.plotDtRpc()
	resDtRpcFine = lib.plotDtRpcFine()
	resDtRpcNotFine = lib.plotDtRpcNotFine()
	resDt = lib.plotDt()
	checkUserInput()
	return

def plotHoTimeVsEta():
	lib = HoTimeVsEta(filename=args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	res = lib.plotTimeVsPhi()
	res2 = lib.plotTimeVsPhiTight()
	res3 = lib.plotTimeVsPhiDtRpc()
 	resCombined = lib.plotCombined()
 	resCombinedTight = lib.plotCombinedTight()
 	resL1TimeVsEta = lib.plotL1TimeVsEta()
 	resHoTimeEta = lib.plotHoTimeVsEta()
 	resHoTimeEtaBxWrng = lib.plotHoTimeVsEtaBxWrong()
 	resTightHoTimeEtaBxWrng = lib.plotTightHoTimeVsEtaBxWrong()
 	resTightDtRpc = lib.plotHoTimeVsEtaDtRpcTight()
	checkUserInput()
	return

def plotDtOnlyCoordinates():
	lib = DtOnlyCoordinates(filename=args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	resDtEtaFine = lib.doAllEtaFineBitPlots()
	checkUserInput()
	resDtOnlyCoordinates = lib.plotDtOnlyCoordinates()
	resDtOnlyTightCoordinates = lib.plotDtOnlyTightCoordinates()
	resDtOnlyBxWrong = lib.plotDtOnlyBxWrongCoordinates()
	resDtOnlyTightBxWrong = lib.plotDtOnlyTightBxWrongCoordinates()
	resDtOnlyAndHoBxWrong = lib.plotDtOnlyAndHoBxWrongCoordinates()
	resDtOnlyTightAndHoBxWrong = lib.plotDtOnlyTightAndHoBxWrongCoordinates()
	checkUserInput()
	lib.printFractionsForDtOnly()
	return

def plotTiming():
	lib = Timing(filename=args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
 	resEvsTime = lib.plotHoEnergyVsTime()
 	resDeltaTime = lib.plotDeltaTime()
 	resBxId = lib.plotL1BxId()
 	resTightBxId = lib.plotL1BxId(tight=True)
 	resTimeHo = lib.plotMatchedHoTime()
 	res = lib.plotHoTime()
 	res6 = lib.plotHoTimeLog()
 	res4 = lib.plotImprovementInDt()
 	resTightImpr = lib.plotImprovementInTightDt()
 	if not args.data:
 		lib.plotEtaOfWrongBxId()
 		lib.plotEtaPhiOfWrongBxId()
 		lib.plotFractionsOfBxId()
 		res2 = lib.plotDetectorContributionsToTiming()
 		res5 = lib.plotPtAndPhiOfWrongBxId()
 		res4 = lib.plotImprovementInDt()
 		res3 = lib.plotPtAndEtaOfWrongBxId()
	checkUserInput()
	return

def plotEVsEtaPhi():
	lib = EvsEtaPhi(filename = args.source, data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	lib.calculateCentralFractionInTight()
	reseMaxCountsTight = lib.plotEMaxCountsForTightMuons()
	reseMaxCounts = lib.plotEMaxCounts()
	res = lib.plotAverageEnergyAroundL1()
	res2 = lib.plotAverageEMaxAroundL1()
	res3 = lib.plot1DEnergyAroundL1()
	res4 = lib.plot1DEMaxAroundL1()
	res5 = lib.compareHistogramMethods()
	res6 = lib.plotEavForTightMuons()
	res7 = lib.plotEavPerWheelForTightMuons()
	resWheelwise = lib.plotEAveragePerWheel()
	resEtaPhiTight = lib.plotEtaPhiForTightL1()
	checkUserInput()
	return

def plotControlPlots():
	lib = ControlPlots(filename=args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	resEtaFine = lib.plotL1EtaVsPatEtaFine()
	resEtaNotFine = lib.plotL1EtaVsPatEtaNotFine()
	resEtaFineTight = lib.plotL1EtaVsPatEtaFineTight()
	resEtaNotFineTight = lib.plotL1EtaVsPatEtaNotFineTight()
	resNHits = lib.plotNHitsPerL1()
	resNHitsTight = lib.plotNHitsPerTightL1()
	resHitsVsPt3x3 = lib.plotNHitsVsPt()
	resTightHitsVsPt3x3 = lib.plotTightNHitsVsPt()
	res1 = lib.plotL1PerPt()
	res2 = lib.plotHoIEtaIPhi()
	res3 = lib.plotHoEtaPhiMatchedToL1()
	resTightAndHoIetaIphi = lib.plotHoIEtaIPhiMatchedToTightL1()
	resL1AndHoIetaIphi = lib.plotHoIEtaIPhiMatchedToL1()
	resPatTightHoEtaPhi = lib.plotHoEtaPhiMatchedToTightL1()
	resSameScale = lib.plotIEtaIPhiOnSameScales()
	resL1Multiplicity = lib.plotL1PresentMultiplicity()
	#res2 = lib.plotHoDigiMatchesPerDetId()
	if not args.data:
		res4 = lib.plotEfficiencyCountCheck()
		res5 = lib.plotGenEtaPhi()
	checkUserInput()
	return

def plotTimeWindow():
	lib = TimeWindow(filename=args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	resAllL1 = lib.plotAllL1Together()
	resBxidFail = lib.plotBxidVsPtFails()
	resBxidSucc = lib.plotBxidVsPtMatch()
	if not args.data:
		resTimeWindowAlone = lib.plotTimeWindowAlone()
		resTruthL1 = lib.plotTruthL1Together()
	checkUserInput()
	return

def plotPtResolution():
	lib = PtResolution(filename=args.source,truth=False,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	res4 = lib.plotL1()
	res5 = lib.plotL1Tight()
	res6 = lib.plotL1AndHo()
	res7 = lib.plotL1TightAndHo()
	res8 = lib.plotL1NotHo()
	res9 = lib.plotL1TightNotHo()
	res1 = lib.plotPtResolutionHistograms()
	res2 = lib.plotTightPtResolution()
	res3 = lib.plotLoosePtResolution()
	resTogether = lib.printThresholds()
	checkUserInput()
	return

def plotTruthPtResolution():
	libTruth = PtResolution(filename=args.source,truth=True,data=args.data,debug = args.DEBUG)
	updateModuleName(libTruth)
	resTruth4 = libTruth.plotL1()
	resTruth5 = libTruth.plotL1Tight()
	resTruth6 = libTruth.plotL1AndHo()
	resTruth7 = libTruth.plotL1TightAndHo()
	resTruth8 = libTruth.plotL1NotHo()
	resTruth9 = libTruth.plotL1TightNotHo()
	resTruth1 = libTruth.plotPtResolutionHistograms()
	resTruth2 = libTruth.plotTightPtResolution()
	resTruth3 = libTruth.plotLoosePtResolution()
		
	resTogether = libTruth.printThresholds()

	checkUserInput()

def plotQualityCodes():
	lib = QualityCode(filename=args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	res1 = lib.plot3x3MatchQualityCodesVsPt()
	res2 = lib.plot3x3FailQualityCodesVsPt()
	res3 = lib.plotQualityCodesStacked(1)
	res4 = lib.plotAllQualitiyCodes()
	checkUserInput()
	return

def plotCounters():
	lib = Counters(filename=args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	res = lib.plotL1AndTightL1Counters()
	res2 = lib.plotTightL1EtaPhiRatio()
	checkUserInput()
	return

def plotThresholdScan():
	lib = HoThresholdScan(filename=args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	res = lib.plotHoThresholdScan()
	checkUserInput()
	return

def plotEfficiency():
	lib = GridMatching(filename=args.source, data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	#res = lib.plotL1GridMatchingEfficiency()
	res2 = lib.plotL13x3AndL1Tight3x3()
	res3 = lib.plotL13x3AndL1Tight3x3L1Coordinates()
	res4 = lib.plotL13x3AndL1Tight3x3FromPatL1Pt()
	res5 = lib.plotL13x3AndL1Tight3x3FromPat()
	c = ROOT.TCanvas('123','asdf')
	c.cd()
	res2[3].SetMarkerStyle(23)
	temp = res3[3].Draw()
	temp2 = res4[3].Draw('same')
	c.Update()
	if not args.data:
		resL1Truth = lib.plotL1TruthGridMatchingPlot()
		res3x3Together = lib.plot3x3GridTogether()
		resN3x3 = lib.plotNtotalGridMatching3x3()
		res5x5Together = lib.plot5x5GridTogether()
	checkUserInput()
	return

def plotEnergy():
	lib = Energy(filename = args.source, data = args.data,debug = args.DEBUG)
	updateModuleName(lib)
	resEnergy = lib.plotEnergy()
	resEnergyNorm = lib.plotEnergyNormalized()
	resMipNorm = lib.plotEnergyNormalizedToMip()
	resEPerWheelHo = lib.plotHoEnergyPerWheel()
	resEPerWheelMatchedHo = lib.plotMatchedHoEnergyPerWheel()
	resEPerWheelTogether = lib.plotMatchedAndNotMatchedPerWheel()
	checkUserInput()
	return

def plotCompareEnergy():
	lib = Comparison(data = args.data,debug = args.DEBUG)
	updateModuleName(lib)
	resL1Count = lib.compareL1Count()
	resEPerWheel = lib.compareEnergyPerWheel()
	resEAbsolute = lib.compareEnergyAbsolute()
	resEAbsoluteMatched = lib.compareEnergyAbsoluteHoMatched()
	resEIntegralNorm = lib.compareEnergyNormalizedToIntegral()
	resEIntegralNormTight = lib.compareEnergyTightNormalizedToIntegral()
	checkUserInput()
	
	return

def plotPhiShift():
	lib = DeltaPhi(filename = args.source,data=args.data,debug = args.DEBUG)
	updateModuleName(lib)
	resDeltaPhi = lib.plotDeltaPhiHistogram()
	resEta = lib.plotDeltaPhiVsL1Eta()
	resdPhiVsPtL1 = lib.plotDeltaPhiVsL1Pt()
	resdPhiVsPtL1Tight = lib.plotDeltaPhiVsL1TightPt()
	resAllEtaPhi = lib.plotEtaPhiForAllL1()
	resEtaPhiMap = lib.plotEtaPhiForDeltaPhiOne()
	if not args.data:
		#res3 = lib.plotDeltaPhiVsGenPt()
		res4 = lib.plotL1PhiVsHoPhi()
		res5 = lib.plotL1PhiVsHoIPhi()
		#res = lib.plotDeltaPhiVsL1Phi()
		#res2 = lib.plotDeltaPhiVsL1Pt()
	checkUserInput()
	return

from ROOT import gROOT
import ROOT
print ROOT.__file__
gROOT.ProcessLine(".L $HOMUONTRIGGER_BASE/python/loader.C+");

scripts = ['controlPlots','eVsEtaPhi','timeWindow','ptResolution','ptResolutionTruth','qualityCodes',
			'counters','thresholdScan','efficiency','energy','compareEnergy','timing','phiShift',
			'dtOnly', 'hoTimeVsEta','dtRpc','dttfCands','noL1Muons','pilup']

if args.scripts:
	for script in args.scripts:
		if(script == 'controlPlots'):
			plotControlPlots()
		elif(script == 'eVsEtaPhi'):
			plotEVsEtaPhi()
		elif(script == 'phiShift'):
			plotPhiShift()
		elif(script == 'timeWindow'):
			plotTimeWindow()
		elif(script == 'ptResolution'):
			plotPtResolution()
		elif(script == 'ptResolutionTruth'):
			plotTruthPtResolution()
		elif(script == 'qualityCodes'):
			plotQualityCodes()
		elif (script=='timing'):
			plotTiming()
		elif (script == 'counters'):
			plotCounters()
		elif (script == 'thresholdScan'):
			plotThresholdScan()
		elif (script == 'efficiency'):
			plotEfficiency()
		elif (script == 'energy'):
			plotEnergy()
		elif (script == 'comparison'):
			plotCompareEnergy()
		elif (script == 'dtOnly'):
			plotDtOnlyCoordinates()
		elif (script == 'hoTimeVsEta'):
			plotHoTimeVsEta()
		elif (script == 'dtRpc'):
			plotDtRpc()
		elif (script == 'dttfCands'):
			plotDttfCands()
		elif (script == 'noL1Muons'):
			plotNoL1Muons()
		elif (script == 'pileup'):
			plotPileup()
		else:
			print 'Unknown script requested: %s' % (script)
			print "Available Scripts:"
			for script in scripts:
				print '\t',script
else:
	if args.plotAll:
		[t() for f,t in locals().items() if type(t) == FunctionType and (f.find('plot') == 0)]
		print "#########################"
		print "# And that's it, folks! #"
		print "#########################"

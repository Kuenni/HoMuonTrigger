from ROOT import TCanvas,TFile,ROOT,TPaveText,TGraph,TLegend,TMarker
import os,sys
sys.path.append(os.path.abspath("../../python"))

import PlotStyle

file = TFile.Open('L1MuonHistogram.root')

if( not os.path.exists('plots')):
	os.mkdir('plots')
if( not os.path.exists('plots/graphsEtaPhi')):
	os.mkdir('plots/graphsEtaPhi')

hTotalEvents = file.Get('hoMuonAnalyzer/count/Events_Count')

totalEvents = hTotalEvents.GetBinContent(2)
print totalEvents

hNoSingleMuEventsInGa = file.Get("hoMuonAnalyzer/count/NoTrgTdmiInGA_Count")
noSingleMuEventsInGa = hNoSingleMuEventsInGa.GetBinContent(2)
print noSingleMuEventsInGa

PlotStyle.setPlotStyle()

def plotAllEventsInAcceptance():
	c = TCanvas("c","c",1200,1200)
	gTdmiInGa = file.Get("hoMuonAnalyzer/graphs/tdmiInGA")
	gTdmiInGa = PlotStyle.convertToHcalCoords(gTdmiInGa)
	gTdmiInGa.GetXaxis().SetTitle("i#eta / a.u.")
	gTdmiInGa.GetYaxis().SetTitle("i#phi / a.u.")
	gTdmiInGa.SetMarkerStyle(6)
	gTdmiInGa.SetMarkerColor(PlotStyle.colorRwthDarkBlue)
	gTdmiInGa.Draw("AP")
	
	pText = TPaveText(0.7,0.85,0.9,0.9,'NDC')
	pText.AddText('Total Events: %d' % (totalEvents))
	pText.AddText('Events in Plot: %d' % (gTdmiInGa.GetN()))
	pText.SetBorderSize(1)
	pText.Draw()
	
	chimney1 = PlotStyle.chimney1
	chimney2 = PlotStyle.chimney2
	
	chimney1Converted = PlotStyle.convertToHcalCoords(chimney1)
	chimney2Converted = PlotStyle.convertToHcalCoords(chimney2)
	chimney1Converted.SetLineColor(PlotStyle.colorRwthMagenta)
	chimney2Converted.SetLineColor(PlotStyle.colorRwthMagenta)
	chimney1Converted.Draw("same,l")
	chimney2Converted.Draw("same,l")
	
	labelCmsPrivateSimulation = PlotStyle.getLabelCmsPrivateSimulation()
	labelCmsPrivateSimulation.Draw()
	
	legend = TLegend(0.1,0.87,0.3,0.9)
	legend.AddEntry(chimney2Converted,"chimney","l")
	legend.Draw()
	
	c.Update()
	
	c.SaveAs("plots/graphsEtaPhi/gTdmiInGa.pdf")
	c.SaveAs("plots/graphsEtaPhi/gTdmiInGa.png")

def plotHoMatchFail():
	c2 = TCanvas("c2","c2",1200,1200)
	
	gNoTrgNoL1TdmiInGa = file.Get("hoMuonAnalyzer/graphs/NoTrgNoL1TdmiInGA")
	gNoTrgNoL1TdmiInGa.GetXaxis().SetTitle("#eta")
	gNoTrgNoL1TdmiInGa.GetYaxis().SetTitle("#phi")
	gNoTrgNoL1TdmiInGa.SetMarkerStyle(6)
	gNoTrgNoL1TdmiInGa.Draw("AP")
	chimney1.Draw("same,l")
	chimney2.Draw("same,l")
	pText3 = TPaveText(0.7,0.85,0.9,0.9,'NDC')
	pText3.AddText('Total Events: %d' % (totalEvents))
	pText3.AddText('Events in Plot: %d' % (gNoTrgNoL1TdmiInGa.GetN()))
	pText3.SetBorderSize(1)
	pText3.Draw()
	
	c2.SaveAs("plots/graphsEtaPhi/gNoTrgNoL1TdmiInGa.pdf")

def plotHoMatchFailInHoCoords():
	c3 = TCanvas("c3","c3",1200,1200)
	  
	gHoMatchFail = file.Get("hoMuonAnalyzer/graphs/NoTrgHoMatchFail")
	gHoMatchFail.GetXaxis().SetTitle("#eta")
	gHoMatchFail.GetYaxis().SetTitle("#phi")
	gHoMatchFail.SetMarkerStyle(6)
	gHoMatchFail.Draw("AP")
	chimney1.Draw("same,l")
	chimney2.Draw("same,l")
	labelCmsPrivateSimulation.Draw()
	pText2 = TPaveText(0.7,0.85,0.9,0.9,'NDC')
	pText2.AddText('Total Events: %d' % (totalEvents))
	pText2.AddText('Events in Plot: %d' % (gHoMatchFail.GetN()))
	pText2.SetBorderSize(1)
	pText2.Draw()
	
	c3.SaveAs("plots/graphsEtaPhi/gHoMatchFail.pdf")

def plotNoL1AndTdmiInAcceptance():
	c4 = TCanvas("c4","c4",1200,1200)
	gConvertedToiEtaiPhi = PlotStyle.convertToHcalCoords(gNoTrgNoL1TdmiInGa)
	gConvertedToiEtaiPhi.SetTitle('No Single #mu Trg, in HO Acceptance')
	gConvertedToiEtaiPhi.SetMarkerColor(PlotStyle.colorRwthDarkBlue)
	gConvertedToiEtaiPhi.Draw("AP")
	chimney1Converted = PlotStyle.convertToHcalCoords(chimney1)
	chimney2Converted = PlotStyle.convertToHcalCoords(chimney2)
	chimney1Converted.SetLineColor(PlotStyle.colorRwthMagenta)
	chimney2Converted.SetLineColor(PlotStyle.colorRwthMagenta)
	chimney1Converted.Draw('Same,l')
	chimney2Converted.Draw('Same,l')
	legend = TLegend(0.7,0.85,0.9,0.9)
	legend.AddEntry(chimney2Converted,"chimneys","le")
	legend.Draw()
	
	label = PlotStyle.getLabelCmsPrivateSimulation()
	label.Draw()
	c4.Update()
	#boxes = PlotStyle.drawHcalBoxesHcalCoords(c4)
	c4.SaveAs("plots/graphsEtaPhi/gNoTrgNoL1TdmiInGaHocoords.pdf")
	c4.SaveAs("plots/graphsEtaPhi/gNoTrgNoL1TdmiInGaHocoords.png")
    
    
def plotEventsInAcceptance():
	c = TCanvas("c","c",1200,1200)
	
	gTdmiInGaNotConverted = file.Get("hoMuonAnalyzer/graphs/tdmiInGaNotDead")
	gTdmiInGa = PlotStyle.convertToHcalCoords(gTdmiInGaNotConverted)
	gTdmiInGa.GetXaxis().SetTitle("i#eta")
	gTdmiInGa.GetYaxis().SetTitle("i#phi")
	gTdmiInGa.SetMarkerStyle(6)
	gTdmiInGa.SetMarkerColor(PlotStyle.colorRwthDarkBlue)
	gTdmiInGa.SetTitle("#eta #phi plot of HO geom. Acceptance and not dead channels")
	gTdmiInGa.Draw("AP")
	
	pText = TPaveText(0.7,0.85,0.9,0.9,'NDC')
	pText.AddText('Total Events: %d' % (totalEvents))
	pText.AddText('Events in Plot: %d' % (gTdmiInGa.GetN()))
	pText.SetBorderSize(1)
	pText.Draw()
	
	chimney1 = PlotStyle.chimney1
	chimney2 = PlotStyle.chimney2
	labelCmsPrivateSimulation = PlotStyle.getLabelCmsPrivateSimulation()
	chimney1Converted.Draw("same,l")
	chimney2Converted.Draw("same,l")
	labelCmsPrivateSimulation.Draw()
	legend = TLegend(0.1,0.87,0.3,0.9)
	legend.AddEntry(chimney2Converted,"chimney","l")
	legend.Draw()
	c.Update()
	c.SaveAs("plots/graphsEtaPhi/gTdmiInGaNotDead.png")
   

#####
# HO > 0.2 GeV; No Single Muon Trigger
#####
def plotHoAboveThr():
	c5 = TCanvas("c5","c5",1200,1200)
	
	noTrgTdmiAboveThrNotConverted = file.Get("hoMuonAnalyzer/graphs/NoTrgTdmiAboveThr")
	noTrgTdmiAboveThr = PlotStyle.convertToHcalCoords(noTrgTdmiAboveThrNotConverted)
	noTrgTdmiAboveThr.GetXaxis().SetTitle("i#eta / a.u.")
	noTrgTdmiAboveThr.GetYaxis().SetTitle("i#phi / a.u.")
	noTrgTdmiAboveThr.SetMarkerStyle(6)
	noTrgTdmiAboveThr.SetMarkerColor(PlotStyle.colorRwthDarkBlue)
	noTrgTdmiAboveThr.SetTitle("#eta #phi plot of HO > 0.2 GeV in no Single #mu Trg. events")
	noTrgTdmiAboveThr.Draw("AP")
	
	pText = TPaveText(0.7,0.85,0.9,0.9,'NDC')
	pText.AddText('No Single #mu in GA: %d' % (noSingleMuEventsInGa))
	pText.AddText('Events in Plot: %d' % (noTrgTdmiAboveThr.GetN()))
	pText.SetBorderSize(1)
	pText.Draw()
	
	chimney1 = PlotStyle.chimney1
	chimney2 = PlotStyle.chimney2
	labelCmsPrivateSimulation = PlotStyle.getLabelCmsPrivateSimulation()
	chimney1Converted.Draw("same,l")
	chimney2Converted.Draw("same,l")
	labelCmsPrivateSimulation.Draw()
	legend = TLegend(0.1,0.87,0.3,0.9)
	legend.AddEntry(chimney2Converted,"chimney","l")
	legend.Draw()
	#c5.SetGridY(0)
	#c5.SetGridX(0)
	c5.Update()
	c5.SaveAs("plots/graphsEtaPhi/gNoTrgTdmiHoAboveThr.png")

def plotFailedHoMatchesNoTrg():
	c = TCanvas("cFailedHoMatchesNoTrg","cFailedHoMatchesNoTrg",1200,1200)
	c.cd().SetRightMargin(0.25)
	c.cd().SetLeftMargin(0.08)
	
	#Graph for in events not in geometric acceptance
	grNotInGaNC = file.Get("hoMuonAnalyzer/graphs/NoTrgTdmiNotInGA")
	grNotInGa = PlotStyle.convertToHcalCoords(grNotInGaNC)
	grNotInGa.GetYaxis().SetTitle("i#phi / a.u.")
	grNotInGa.GetXaxis().SetTitle("i#eta / a.u.")
	grNotInGa.GetYaxis().SetTitleFont(62)
	grNotInGa.GetYaxis().SetLabelFont(62)
	grNotInGa.SetMarkerStyle(6)
	grNotInGa.SetMarkerColor(PlotStyle.colorRwthDarkBlue)
	grNotInGa.SetTitle("#eta #phi plot failed HO matches in no Single #mu Trg. events")
	grNotInGa.Draw("AP")
	
	#Graph of events with HO match below threshold
	grHoBelowThrNC = file.Get("hoMuonAnalyzer/graphs/NoTrgTdmiBelowThr")
	grHoBelowThr = PlotStyle.convertToHcalCoords(grHoBelowThrNC)
	grHoBelowThr.SetMarkerStyle(20)
	grHoBelowThr.SetMarkerSize(1.2)
	grHoBelowThr.SetMarkerColor( PlotStyle.colorRwthGruen )
	grHoBelowThr.Draw("samep")
	
	#Graph for events where HO matching failed
	grHoMatchFailNC = file.Get("hoMuonAnalyzer/graphs/NoTrgHoMatchFail")
	grHoMatchFail = PlotStyle.convertToHcalCoords(grHoMatchFailNC)
	grHoMatchFail.SetMarkerStyle(21)
	grHoMatchFail.SetMarkerSize(1)
	grHoMatchFail.SetMarkerColor( PlotStyle.colorRwthRot )
	grHoMatchFail.Draw("samep")
	
	#Draw chimneys
	chimney1Converted = PlotStyle.convertToHcalCoords(PlotStyle.chimney1)
	chimney2Converted = PlotStyle.convertToHcalCoords(PlotStyle.chimney2)
	chimney1Converted.SetLineColor(PlotStyle.colorRwthMagenta)
	chimney2Converted.SetLineColor(PlotStyle.colorRwthMagenta)
	chimney1Converted.Draw('same')
	chimney2Converted.Draw('same')
	
	#cms private label
	label = TPaveText(PlotStyle.getLabelCmsPrivateSimulation(x1ndc=0.5,x2ndc=0.75))
	label.Draw()
	
	#create extra marker for the legend
	marker = TMarker(1,1,2)
	marker.SetMarkerSize(3)
	marker.SetMarkerColor(PlotStyle.colorRwthDarkBlue)
	
	#legend
	legend = TLegend(0.75,0.8,0.99,0.9)
	legend.AddEntry(chimney2Converted,"chimney","l")
	legend.AddEntry(marker,'Not in GA','p')
	legend.AddEntry(grHoMatchFail,'HO match fail','p')
	legend.AddEntry(grHoBelowThr,'HO match < 0.2 GeV','p')
	legend.Draw()
	
	nNotMatching = grHoMatchFail.GetN()
	nNotInGa = grNotInGa.GetN()
	nBelowThr = grHoBelowThr.GetN()
	nTotal = nNotMatching + nNotInGa + nBelowThr
	
	print 80*'#'
	print 'Not Matching:\t%5d/%d\t=> %5.2f%% +- %f%%' % (nNotMatching,nTotal,nNotMatching/float(nTotal)*100,PlotStyle.calcSigma(nNotMatching,float(nTotal)))
	print 'Not in GA:\t%5d/%d\t=> %5.2f%% +- %f%%' % (nNotInGa,nTotal,nNotInGa/float(nTotal)*100,PlotStyle.calcSigma(nNotInGa,float(nTotal)))
	print 'Below Thr:\t%5d/%d\t=> %5.2f%% +- %f%%' % (nBelowThr,nTotal,nBelowThr/float(nTotal)*100,PlotStyle.calcSigma(nBelowThr,float(nTotal)))
	print 80*'#'
	
	c.Update()
	c.SaveAs('plots/graphsEtaPhi/gNoTrgHoMatchingFailed.png')
	c.SaveAs('plots/graphsEtaPhi/gNoTrgHoMatchingFailed.pdf')
	return c,grNotInGa,label,chimney1Converted,chimney2Converted,legend,grHoMatchFail,grHoBelowThr


res = plotFailedHoMatchesNoTrg()

raw_input('--> Enter')

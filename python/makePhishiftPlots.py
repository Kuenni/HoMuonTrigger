#!/usr/bin/python

from phishift.DeltaPhi import *

resWheelwise = plotEAveragePerWheel()
raw_input('-->')
resEtaPhiMap = plotEtaPhi()
resDeltaPhi = plotDeltaPhiHistogram()
resEta = plotDeltaPhiVsL1Eta()
res4 = plotL1PhiVsHoPhi()
res5 = plotL1PhiVsHoIPhi()
res = plotDeltaPhiVsL1Phi()
res2 = plotDeltaPhiVsL1Pt()
res3 = plotDeltaPhiVsGenPt()

raw_input('-->')

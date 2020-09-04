import TASKit
import numpy as np

import sys

##Test load
test_load = TASKit.TrajHandle.TrajLoader(prmtop="test/Peptide.prmtop", DCDS="test/Peptide.dcd", stride=1 )
print ("Test Load: {}".format(test_load))

##Test EP timeseries generation


EPTSCorr = test_load.GetPotentialEnergies(ProgressInt=10000)
print ("EPTS Corr: {}".format(EPTSCorr))
print ("EPTS Corr length: {}".format(EPTSCorr.shape))


##Decorrelate EP TS
EPTSDecorr = TASKit.Decorrelation.DecorrelateTimeseries(EPTSCorr, stride=5)
print ("EPTS Decorr: {}".format(EPTSDecorr))

##Generate Plot

#EPTSDecorr.GenerateNEFFPlot()

##Print decorr frames
DecorrFrames = EPTSDecorr.GetDecorrSamples()
#print ("Decorr Frames: {}".format(DecorrFrames))

##Generate new DCD with Decorr Frames
TrimmedTraj = test_load.TrimTraj(DecorrFrames)
print ("TrimmedTraj: {}".format(TrimmedTraj))

RMSDMat = TASKit.RMSDAvA(TrimmedTraj, RMSDSele="backbone")
print (RMSDMat)

RMSDMat.RMSDAvAHeat
import TASKit
import numpy as np

import sys

##Test load
test_load1 = TASKit.TrajHandle.TrajLoader(prmtop="test/Peptide.prmtop", DCDS="test/Peptide_01.dcd", stride=1)
test_load2 = TASKit.TrajHandle.TrajLoader(prmtop="test/Peptide.prmtop", DCDS="test/Peptide_02.dcd", stride=1)

##Test EP timeseries generation
EPTSCorr1 = test_load1.GetPotentialEnergies(ProgressInt=10000)
EPTSCorr2 = test_load2.GetPotentialEnergies(ProgressInt=10000)


##Decorrelate EP TS
EPTSDecorr1 = TASKit.Decorrelation.DecorrelateTimeseries(EPTSCorr1, stride=5)
EPTSDecorr2 = TASKit.Decorrelation.DecorrelateTimeseries(EPTSCorr2, stride=5)

##Generate Plot

##Print decorr frames
DecorrFrames1 = EPTSDecorr1.GetDecorrSamples()
DecorrFrames2 = EPTSDecorr2.GetDecorrSamples()

##Generate new DCD with Decorr Frames
TrimmedTraj1 = test_load1.TrimTraj(DecorrFrames1)
TrimmedTraj2 = test_load2.TrimTraj(DecorrFrames2)

RMSDMat = TASKit.RMSDAvA.Calc(TrajIn=[TrimmedTraj1, TrimmedTraj2], RMSDSele="backbone")
print ("RMSD FrameIntervals are:\n{}".format(RMSDMat.frameIntervals))
print ("RMSD Matrix:\n{}".format(RMSDMat.RMSD_Matrix))

RMSDMat.RMSDAvAHeat()

Clusters = TASKit.Cluster.ClusterizeMatrix(RMSDMat.RMSD_Matrix, 3.25)
Clusters.GenerateHeatmap()

print ("If this message has been printed, that means the tests have been passed succesfully. \
You have all necesarry modules installed. Have a nice day!")

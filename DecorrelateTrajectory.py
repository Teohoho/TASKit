import simtk.openmm.app
import simtk.openmm
from simtk.unit import *

import mdtraj as md
import numpy as np

import pymbar.timeseries
import matplotlib.pyplot as plt

class DecorrelateTrajectory:

    """A class that takes a prmtop and any number of DCDs and generates 
        a potential energy timeseries, using openMM and MDTraj.
         
    Parameters
    ----------
    
    prmtop: str
             A list object containing the PRMTOPs of the systems for which RMSDS
             will be calculated
             
    DCDS: list or str
             A DCD or a list of DCDs to be loaded into the prmtop

    stride: int, default=1
             Use this stride when loading the trajectories to sparsify the trajectories
             
    """

    def __init__(self, prmtop, DCDS, stride=1):

        DCDsToLoad = []     ##This list serves to make the final list of systems more clear
            
        MDtrajTrajectoryObjectPRMTOP = md.load_prmtop(prmtop)
   
       ## If the user provides only one DCD:
        if (isinstance(DCDS,str) == True):
            MDtrajTrajectoryObject = md.load_dcd(DCDS, top= MDtrajTrajectoryObjectPRMTOP, stride=stride)
            print ("Loaded system {} ({} frames)".format(prmtop, MDtrajTrajectoryObject.n_frames))
            
        ## If the user provides multiple DCDs, in the form of a list or a tuple:
        elif (isinstance(DCDS, list) == True):
            for DCD in DCDS:
                MDtrajTrajectoryObject = md.load_dcd(DCD, top= MDtrajTrajectoryObjectPRMTOP, stride=stride)
                print ("Loaded DCD {} ({} frames)".format(DCD, MDtrajTrajectoryObject.n_frames))
                DCDsToLoad.append(MDtrajTrajectoryObject)
            
            MDtrajTrajectoryObject = md.join(DCDsToLoad)
            print ("Loaded system {} ({} frames)".format(prmtop, MDtrajTrajectoryObject.n_frames))
            
        else:
            print ("Please supply DCDS either as one string, if one DCD, or as a List of strings if multiple. You supplied {}".format(type(DCDS)))
            
        ##Create OpenMM system
        self.OpenMMPrmtop = simtk.openmm.app.AmberPrmtopFile(prmtop)
        
        self.MDtrajTrajectoryObject = MDtrajTrajectoryObject

                            
    def GetPotentialEnergies (self, verbose=True, ProgressInt=100, TimeseriesOut=None):
       
        """Generate the Potential Energy timeseries of the loaded trajectory.
       
        Parameters
        ----------
        
        verbose: bool
                 Print progress of energy calculation (useful for long trajectories
                 to be sure the proccess didn't hang).
                 
        ProgressInt: int
                How many frames' energy to calculate before printing a progress report.
        """
    
    
        self.T        = 0*kelvin
        self.timestep = 0.002*picoseconds       
        
        ##Create OpenMM system
        
        OpenMMSystem = self.OpenMMPrmtop.createSystem(implicitSolvent=None,nonbondedMethod=simtk.openmm.app.NoCutoff) ##Vacuum

        
        NOFrames =  self.MDtrajTrajectoryObject.n_frames
        
        ##Create array to store energies
        EnergiesArray = np.zeros(NOFrames)
        
        ##Construct OpenMMSimulation Object
        OpenMMIntegrator = simtk.openmm.LangevinIntegrator (self.T, 1/picosecond, self.timestep)
        OpenMMSimulation = simtk.openmm.app.Simulation(self.OpenMMPrmtop.topology, OpenMMSystem, OpenMMIntegrator, 
                                        platformProperties={'DisablePmeStream':'true'})
        
        ##Iterate over all frames in MDtrajTrajectoryObject, and copy positions to 
        ##OpenMM system
        
        ##A small bit of code that prints progress messages every 10%-th of progress
        progressCheck = 0
        
        for FrameNumber in range(NOFrames):
            
            if ((progressCheck % ProgressInt == 0) and (verbose == True)):
                print ("Energy calculation is {:4.1f}% done!".format(progressCheck*100/NOFrames))
            
            OpenMMSimulation.context.setPositions(self.MDtrajTrajectoryObject.openmm_positions(FrameNumber))
            
            ##Check for unitcell vectors, and apply them if they exist
            if (self.MDtrajTrajectoryObject.openmm_boxes(FrameNumber) is not None):
                 OpenMMSimulation.context.setPeriodicBoxVectors(self.MDtrajTrajectoryObject.openmm_boxes(FrameNumber)[0], 
                                                                self.MDtrajTrajectoryObject.openmm_boxes(FrameNumber)[1],
                                                                self.MDtrajTrajectoryObject.openmm_boxes(FrameNumber)[2])
                                                                
            EnergiesArray[FrameNumber] = OpenMMSimulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
            progressCheck += 1
        
        self.EnergiesArray = EnergiesArray
        
        if (TimeseriesOut is not None):
            print ("Saving PotEn Timeseries to file {}…".format(TimeseriesOut))
            np.save(TimeseriesOut, self.EnergiesArray)
        
        return EnergiesArray
              
    def TimeseriesDecorrelation(self, timeseries=None):
        """A method that takes a timeseries and detects the equilibration point, calculates 
        the statistical inefficiency and the effective sample size, using PyMBAR.
         
        Parameters
        ----------
        
        timeseries: numpy.ndarray
            A numpy array that contains the timeseries to be used. If None, then the one generated by the
                getPotentialEnergy will be used.
            
        """
        
        if (timeseries is None):
            self.timeseries_Corr = self.EnergiesArray
        else:
            self.timeseries_Corr = timeseries
            
        [self.Eq_point, self.statIneff, self.Neff_Max] = pymbar.timeseries.detectEquilibration(self.timeseries_Corr, nskip=1)
        self.timeseries_Decorr = self.timeseries_Corr[self.Eq_point:]
        
        print ("The equilibration is done at sample number {}, the statistical inefficiency for the post-Eq region is {} and the" +
                        "number of Decorrelated samples in the post-Eq region is {}".format(self.Eq_point, self.statIneff, self.Neff_Max))
         
    def GenerateNEFFPlot(self, fig=None):
    
        """A method that generates the Effective Sample Size vs. Equilibration Point plot. Uses matplotlib.pyplot
             
        Parameters
        ----------
        
        fig: str    
            Save the figure to this filename. If None, then it will be showed on-screen.
        
        """
       
        NOfSamples = self.timeseries_Corr.shape[0]
        g_t    = np.ones(NOfSamples - 1)
        Neff_t = np.ones(NOfSamples - 1)
        for t in range(NOfSamples - 1):
            g_t[t] = pymbar.timeseries.statisticalInefficiency(self.timeseries_Corr[t:NOfSamples], fast=True)

            Neff_t[t] = (NOfSamples - t +1)/ g_t[t]

        ##Plot the Neff
        plt.plot(Neff_t, color="k")

        ##Plot a long red line where NeffMax is
        plt.plot([self.Eq_point, self.Eq_point], [0, self.Neff_Max+(self.Neff_Max*0.1)], color="r", linewidth=2)
        
        plt.ylim(0, self.Neff_Max+(self.Neff_Max*0.1))
        
        plt.xlabel("Step",         fontsize=18, labelpad=22)
        plt.ylabel("Neff samples", fontsize=18, labelpad=22)

        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        
        if (fig is not None):
            plt.tight_layout(pad=0.5)
            plt.savefig(fig, dpi=800)
        else:
            plt.show()
                         
    def GetDecorrSamples(self, log=None):
    
        """A method that generates a list containing the Decorrelated Sample Indices, given the decorrelated timeseries and 
            the statistical inefficiency.
             
        Parameters
        ----------
        
        log: str    
            Save the log to this file.
        
        """
     
        decorrelatedSampleIndices = pymbar.timeseries.subsampleCorrelatedData(self.timeseries_Decorr, g=self.statIneff)
        decorrelatedSamples = self.timeseries_Decorr[decorrelatedSampleIndices]
        
        if (log is not None):
            with open(log, "w") as OutFN:
                for i in range(decorrelatedSamples.shape[0]):
                    OutFN.write(str(decorrelatedSampleIndices[i] + self.Eq_point))
                    OutFN.write(" ")
                OutFN.close()
        
        decorrelatedSampleIndicesAll = np.add(decorrelatedSampleIndices, self.Eq_point)
        
        self.decorrelatedSampleIndicesAll = decorrelatedSampleIndicesAll
        
        return decorrelatedSampleIndicesAll     
   
    def TrimTraj(self, dcdOut=None):

        """A method that generates a dcd that contains only the decorrelated frames, computed by  
            by the GetDecorrSamples function.
             
        Parameters
        ----------
        
        dcdOut: str
            Save the decorrelated trajectory to this file. (generates a DCD file).
           
        """
             
        ##Check if the GetDecorrSamples has been used
        if (self.decorrelatedSampleIndicesAll is None):
            print ("Run the GetDecorrSamples first, to compute the Decorrelated Sample Indices")
            return
             
        TrimmedTraj = MDtrajTrajectoryObject.slice(self.decorrelatedSampleIndicesAll)
        print ("The trimmed trajectory has {} frames".format(TrimmedTraj.n_frames))
        
        if (dcdOut is not None):
            f = md.formats.DCDTrajectoryFile(dcdOut, mode='w', force_overwrite=False)
            f.write(TrimmedTraj.xyz)
            f.close()
        
        
        return (TrimmedTraj)
 

def QuickTrim(prmtop, dcdIn, dcdOut=None, stride=1, plot=True):

    """A shortcut function that does what most users will want to do with this
        class. Uses most common options.

    Parameters
    ----------
    
    prmtop: str
             A list object containing the PRMTOPs of the systems for which RMSDS
             will be calculated
             
    dcdIn: list or str
             A DCD or a list of DCDs to be loaded into the prmtop

    stride: int, default=1
             Use this stride when loading the trajectories to sparsify the trajectories
             
    onlyProt: bool, default=True   
             If true, then load only the protein atoms. This saves memory, since most
             DCDs also contain more water atoms than protein atoms and often we want 
             the energy of the protein, without solvent contribution.
             
    plot: bool, default=True
            If True, then it will display the plot that would be generated by
            GenerateNEFFPlot.
            
    """

    InitTraj = DecorrelateTrajectory(prmtop=prmtop, DCDS=dcdIn, stride=stride)
    InitTraj.GetPotentialEnergies()
    InitTraj.TimeseriesDecorrelation()
    
    if (plot == True):
        InitTraj.GenerateNEFFPlot()
     
    InitTraj.GetDecorrSamples()
     
    TrimmedTraj = InitTraj.TrimTraj(dcdOut=dcdOut)
    
    return(TrimmedTraj)
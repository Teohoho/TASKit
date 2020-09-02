import simtk.openmm.app
import simtk.openmm
import simtk.unit

import mdtraj as md
import numpy as np

import matplotlib.pyplot as plt

class TrajHandle:

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
    
    
        T        = 0*simtk.unit.kelvin
        timestep = 0.002*simtk.unit.picoseconds       
        
        ##Create OpenMM system
        
        OpenMMSystem = self.OpenMMPrmtop.createSystem(implicitSolvent=None,nonbondedMethod=simtk.openmm.app.NoCutoff) ##Vacuum

        
        NOFrames =  self.MDtrajTrajectoryObject.n_frames
        
        ##Create array to store energies
        EnergiesArray = np.zeros(NOFrames)
        
        ##Construct OpenMMSimulation Object
        OpenMMIntegrator = simtk.openmm.LangevinIntegrator (T, 1/simtk.unit.picosecond, timestep)
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
                                                                
            EnergiesArray[FrameNumber] = OpenMMSimulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(simtk.unit.kilocalories_per_mole)
            progressCheck += 1
        
        self.EnergiesArray = EnergiesArray
        
        if (TimeseriesOut is not None):
            print ("Saving PotEn Timeseries to file {}…".format(TimeseriesOut))
            np.save(TimeseriesOut, self.EnergiesArray)
        
        return EnergiesArray
              
    def TrimTraj(self, FramesToKeep dcdOut=None):

        """A method that generates a DCD that contains only the frames described in
            FramesToKeep. Written to be used with the indices computed by the
            GetDecorrSamples function, but can be used with any numpy array of ints.
             
        Parameters
        ----------
        
        FramesToKeep: numpy.ndarray,
            Array of ints, coressponding to what frames to export into a new DCD file.
        
        dcdOut: str
            Save the decorrelated trajectory to this file. (generates a DCD file).
           
        """
        
        if (FramesToKeep.dtype != "int64"):
            print ("Warning! The type of the FramesToKeep array is not int64. This may cause errors!")
         
        
        TrimmedTraj = MDtrajTrajectoryObject.slice(self.decorrelatedSampleIndicesAll)
        print ("The trimmed trajectory has {} frames".format(TrimmedTraj.n_frames))
        
        if (dcdOut is not None):
            f = md.formats.DCDTrajectoryFile(dcdOut, mode='w', force_overwrite=False)
            f.write(TrimmedTraj.xyz)
            f.close()
        
        
        return (TrimmedTraj)

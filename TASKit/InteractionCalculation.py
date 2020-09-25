from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

##Define some handy functions beforehand

def NonBondedCalculator(distance,chargeProd,Epsilon,Sigma, NBArray):

    """
    Function that calculates the NB Energies array for a given frame,
    given the distance, chargeProd, Epsilon and Sigma matrices for 
    all atom pairs in a given system.

    """
    
	##Generate an array to store pair-wise NB energy values, with shape depending on n.of interactions
	ChargeProdOverDistance      = np.divide(chargeProd,distance)
	SigmaOverDistance6          = np.power(np.divide(Sigma, distance),6)

	SigmaOverDistance12         = np.power(SigmaOverDistance6, 2)
	SigmaOverDistanceDifference = np.subtract(SigmaOverDistance12, SigmaOverDistance6)

	
	coulomb_part = np.zeros((MDtrajTrajectoryObject.n_atoms, MDtrajTrajectoryObject.n_atoms))
	LJ_part      = np.zeros((MDtrajTrajectoryObject.n_atoms, MDtrajTrajectoryObject.n_atoms))
	
	coulomb_part  = np.multiply(one_4PI_EPS0, ChargeProdOverDistance)
	LJ_part       = np.multiply(np.multiply(4, Epsilon), SigmaOverDistanceDifference)
	AllNBEnergies = np.add(NBArray, (np.multiply(np.add(coulomb_part,LJ_part), 0.239006)))
	return (AllNBEnergies)

def MultiplyListByEachElement(ListIn, parameter):


    """
    Function that generates the prodCharge, Sigma and Epsilon matries
    given a list of each atom's charge, Sigma or Epsilon value respectively.
    Uses the Lorentz-Berthelot combining rule for Sigma and Epsilon.

    """
    
	##Generate an array to store all values
	ArrayFromList = np.zeros((ListIn.shape[0], ListIn.shape[0]))
	if (parameter == "prodCharge"):
		for i in range(ListIn.shape[0]):
                    ArrayFromList[i] = np.multiply(ListIn,ListIn[i])
	elif (parameter == "Sigma"):
		for i in range(ListIn.shape[0]):
		    ArrayFromList[i] = np.add(ListIn,ListIn[i])
		    ArrayFromList[i] = np.divide(ArrayFromList[i],2)
	elif (parameter == "Epsilon"):
		for i in range(ListIn.shape[0]):
		    ArrayFromList[i] = np.multiply(ListIn,ListIn[i])
		    ArrayFromList[i] = np.sqrt(ArrayFromList[i])
	else:
		print ("MultiplyListByEachElement function was given unknown parameter: {}".format(parameter))
	
	
	return (ArrayFromList)


class InteractionCalculator:

    """
    Class that calculates the Interaction Energies between all pairs of atoms.
    Currently only NB energies are computed, but in future Bonded, Angle and 
    Dihedral energies will also be included.
         
    Parameters
    ----------
    
    MDTrajTrajectory: mdtraj.Trajectory object
        MDTraj Trajectory object from which distances will be taken. Loadable
        with TrajHandle class from this package, or suppliable by user.
    
    AmberPRMTOPFile: str
        AMBER parameter file. Needed to initialize a system in openMM.
        
    Frames: numpy.ndarray
        Frames to use in calculation of interaction energies. The final energies
        reported are averaged from the energies of the individual frames
        
    OutFileRoot: str
        Each energy matrix will be saved to its own file, with the appropriate
        suffix. This will be their common root name. Saved as Binary Numpy file (.npy)
        
    OutHeatmap: str
        Save Heatmap to this file, same 
        
    OutHeatmapValues: list of float
        Two values that represent the minimum and maximum values for the heatmap. 
        Optional
             
    """
    
    def __init__ (MDTrajTrajectory, AmberPRMTOPFile, Frames=None, 
                  OutFileRoot=None, OutHeatmap=None, OutHeatmapValues=[None,None]):
                  
        print ("There are {} frames in the chosen DCD.".format(MDtrajTrajectoryObject.n_frames))
        print ("There are {} particles in the system, so a total of {} interactions will be computed".format((MDtrajTrajectoryObject.n_atoms, MDtrajTrajectoryObject.n_atoms**2)))
        
        
        ##Create OpenMM system, from which we will extract 

    
    
##Licence notices
#    TASKit: A Python-based kit containing Post-MD analysis tools.
#    Copyright (C) 2020  Șulea Teodor Asvadur
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>

#imports
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import sys

class Calc:
   
    def __init__ (self, TrajIn, RMSDSele, Verbose=True, RMSDMatrixOut=None, VerbosityCoarseness = 50):
    
        """
        The function that actually calculates the RMSD of the structures
        loaded.
        
        Parameters
        ----------
        
        TrajIn:   MDTrajTrajectory Object
                    Trajectories whose RMSD will be calculated.
        
        
        RMSDSele: str
                    Selection of atoms to use in calculating the RMSD between
                    structures. Note: It must contain the same number of atoms 
                    across all structures. Note2: Atoms and Residues use a 
                    0-based index.
        
        Verbose: bool
                  Print progress for RMSD calculation, useful for long DCDs, to 
                  make sure the process doesn't hang.
                  
        RMSDMatrixOut: str
                    Save the RMSD Matrix to this file. This is useful for doing a series
                    of clusterizations on the same RMSD Matrix, when the time necessary for the 
                    calculation of said matrix is the limiting factor.
                    Note: The Matrix will be saved as Binary Numpy file (.npy).
        VerbosityCoarseness: int
                    If Verbose=True, this sets how frequent the progress prints should be.
                    So, the script only prints when (current_Frame % VerbosityCoarseness) == 0
                    
        
        """

        ## If only 1 traj loaded, it's best to insert it into a list.

        if (isinstance(TrajIn,md.core.trajectory.Trajectory)):
            Trajs = []
            Trajs.append(TrajIn)

            TrajIn = Trajs
        ## Define frameIntervals
        frameIntervals = np.zeros((len(TrajIn), 2), dtype=int)
        total_frames = 0
        for TrajIx in range(len(TrajIn)):
            frameIntervals[TrajIx][0] = total_frames
            total_frames+=TrajIn[TrajIx].n_frames
            frameIntervals[TrajIx][1] = total_frames
            
        print ("{} systems have been loaded, with a total of {} frames".format(len(TrajIn), total_frames))    
        
        self.frameIntervals = frameIntervals
        self.total_frames   = total_frames
        
        
        ##Determine atom_indices for each system:
        atom_indicesRMSD = []
        for syst in TrajIn:
            atom_indicesRMSD.append(syst.topology.select(RMSDSele))
        
        atom_indicesRMSD = np.array(atom_indicesRMSD)
    
        ##Do a check if the RMSDSele containts the same number of atoms in all prmtops
        atom_numbers = np.empty(len(TrajIn), dtype=int)
        for i in range(len(atom_indicesRMSD)):
            atom_numbers[i] = len(atom_indicesRMSD[i])
          
        if (np.amin(atom_numbers) != np.amax(atom_numbers)):
            print ("The number of atoms selected isn't the same in all systems:")
            for i in range(len(TrajIn)):
                print (TrajIn[i] + ": ", atom_numbers[i])
            sys.exit(0)
        
        
        ##Now on to calculate the actual RMSD Matrix
        
        RMSD_ij  = np.zeros(self.total_frames)
        RMSD_Tot = np.full(self.total_frames, -1, dtype=float)
        curr_frame = 1
        
        for refIx in range(len(TrajIn)):
            for frameNum in range(TrajIn[refIx].n_frames):
                if (Verbose == True):
                    if (curr_frame % VerbosityCoarseness == 0):
                        print ("All-vs-All calculation is {:4.1f}% done!".format((curr_frame*100/self.total_frames)))
                    curr_frame +=1
                for targIx in range(len(TrajIn)):
        
                    RMSD_ij[self.frameIntervals[targIx][0]:self.frameIntervals[targIx][1]] = md.rmsd(TrajIn[targIx], TrajIn[refIx], frame=frameNum, atom_indices = atom_indicesRMSD[targIx], ref_atom_indices = atom_indicesRMSD[refIx])
                
                RMSD_Tot = np.vstack((RMSD_Tot, RMSD_ij))
          
        RMSD_Tot = np.delete(RMSD_Tot, 0, 0)
        
        ##Make all diagonals 0 (MDTraj has an issue with this, making the diags a very low number
        ##which is dfferent to 0).
        
        for i in range (RMSD_Tot.shape[0]):
            RMSD_Tot[i][i] = 0
        
        ##Convert from nm to Angstrom:
        
        RMSD_Tot = np.multiply (RMSD_Tot,10)

        print ("RMSD Matrix has been computed.")
        
        ##Save RMSD Array to Binary Numpy file (.npy), to be able to clusterize after generating initial RMSD matrix (Save time)
        
        if (RMSDMatrixOut is not None):
            print ("Saving RMSD matrix to file {}…".format(RMSDMatrixOut))
            np.save(RMSDMatrixOut, RMSD_Tot)

        self.RMSD_Matrix = RMSD_Tot

    def RMSDAvAHeat(self, vmax=None, SaveFigName=None, show=True, saveHTML=None):
        
        """
        This function generates (and optionally saves) a heatmap of
        the RMSD All-vs-All matrix, for easier visualization.
        
        Parameters
        ----------
        
        vmax: int
            Maximum value on the colorbar.        
        
        SaveFigName: str
            Save Generated heatmap to this file 
            
        show: bool
            Show the heatmap
         
        saveHTML: str
            If not None, saves RMSD Heatmap as an interactive HTML object, to this file.
                            
            
        Notes
        -----
        
        frameIntervals:
            An array of shape(NumberOfSystems, 2), which delimits 
            the frames to each system. Useful for visualzing differences 
            between structures. Generated when this class is instantiated.
           
        """

        heatmap = plt.imshow(self.RMSD_Matrix, cmap="jet", vmax=vmax)

        colormap = plt.colorbar(heatmap)
        colormap.ax.tick_params(labelsize=14)

    ##Use lines to separate the prmtops

        for i in range(self.frameIntervals.shape[0]-1):
            limit = self.frameIntervals[i][1]-0.5
        ##Vertical dividing lines
            plt.plot([limit, limit], [0, self.total_frames], color="k", linewidth=2)        
        ##Horizontal dividing lines
            plt.plot([0, self.total_frames], [limit, limit], color="k", linewidth=2)        

        plt.tight_layout(pad=0.5)

    ##Set labels, titles etc
        
        plt.title("All-vs-All RMSD Heatmap ($\AA\,$)", fontsize=20)


        plt.xlabel("Frame #")
        plt.ylabel("Frame #")

        plt.xlim(0,self.total_frames-0.5)
        plt.ylim(self.total_frames-0.5,0)
        plt.tight_layout()
        if (SaveFigName is not None):
            plt.savefig(SaveFigName, dpi=800)
        if (show == True):
            plt.show()
            
    ##Save as HTML:
        if (saveHTML is not None):
            fig = px.imshow(self.RMSD_Matrix, cmap="jet", vmax=vmax)
            fig.write_html(saveHTML)
        

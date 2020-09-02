class RMSDAvA:

    """A class that calculate an All-vs-All type RMSD matrix, save
       it to disk for future use, and generate a detailed heatmap of said matrix
       
    Parameters
    ----------
    
    prmtops: list
             A list object containing the PRMTOPs of the systems for which RMSDS
             will be calculated
    DCDS: list
             A nested list which contains the DCDS for the PRMTOPs. Its length needs
             to be the same as prmtops, since each PRMTOP needs its own DCD list

    stride: int, default=1
             Use this stride when loading the trajectories
    onlyProt: bool, default=True   
             If true, then load only the protein atoms. This saves memory, since most
             DCDs also contain more water atoms than protein atoms, and they are often
             not taken into account when calculating RMSD
             
    Notes
    -----
    
    An example of a valid pair of PRMTOPs and DCDS is:
    prmtops = [prmtop1, prmtop2]
    DCDS = [[DCD11, DCD12, DCD13], [DCD21, DCD22]]
     
    Also, this uses MDTraj to load parameters and trajectories.
    
    """
    
    
    def __init__(self, prmtops,DCDS,stride=1,onlyProt=True): 
        parmlist,trajlist = isinstance(prmtops, list) , isinstance(DCDS, list)
        
        ##Check if prmtops and DCDS are lists
        if (parmlist and trajlist):
            pass
        else:
            print ("Prmtops and DCDS need to be formated as lists.")
            sys.exit(0)
        
        ##Check if they're the same length
        if (len(prmtops) != len(DCDS)):
            print ("The number of PRMTOPs (%d) is different than the number of DCD lists (%d)." % (len(prmtops), len(DCDS)))
            sys.exit(0)   

         ##We assume that we will use many sets of prmtops and DCDs
        listOfTrajs    = []
        frameIntervals = np.zeros((len(prmtops), 2), dtype=int)
        total_frames   = 0   
        prmtopList     = []
        
        
        for PRMTOPIx in range (len(prmtops)):
            DCDsToLoad = []     ##This list serves to maek the final list of systems more clear
            framesLoaded = 0    ##This is used for delimitation of the systems
            
            MDtrajTrajectoryObjectPRMTOP = md.load_prmtop(prmtops[PRMTOPIx])
            if (onlyProt == True):
                AtomsToLoad = MDtrajTrajectoryObjectPRMTOP.select("protein")
            else:
                AtomsToLoad = None
            prmtopList.append(prmtops[PRMTOPIx])
        
            for DCDIx in range(len(DCDS[PRMTOPIx])):
                MDtrajTrajectoryObject = md.load_dcd(DCDS[PRMTOPIx][DCDIx], 
                    top= MDtrajTrajectoryObjectPRMTOP, stride=stride, atom_indices=AtomsToLoad)
                
                print ("Loaded DCD {} ({} frames)".format(DCDS[PRMTOPIx][DCDIx], MDtrajTrajectoryObject.n_frames))
                
                DCDsToLoad.append(MDtrajTrajectoryObject)
                framesLoaded = framesLoaded + MDtrajTrajectoryObject.n_frames
            
            listOfTrajs.append(md.join(DCDsToLoad))
            
            
            frameIntervals[PRMTOPIx][0] = total_frames
            total_frames = total_frames + framesLoaded
            frameIntervals[PRMTOPIx][1] = total_frames
                
            print ("Molecule {} loaded!".format(prmtops[PRMTOPIx]))

        print("{} frames have been loaded!".format(total_frames))
        
        self.listOfTrajs = listOfTrajs
        self.frameIntervals = frameIntervals
        self.total_frames = total_frames
        self.prmtopList = prmtopList

    def CalcRMSD(self, RMSDSele, Verbose=True, RMSDMatrixOut=None, VerbosityCoarseness = 50):
    
        """
        The function that actually calculates the RMSD of the structures
        loaded.
        
        Parameters
        ----------
        
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
        
        ##Determine atom_indices for each system:
        atom_indicesRMSD = []
        for syst in self.listOfTrajs:
            atom_indicesRMSD.append(syst.topology.select(RMSDSele))
        
        atom_indicesRMSD = np.array(atom_indicesRMSD)
        print (atom_indicesRMSD)
    
        ##Do a check if the RMSDSele containts the same number of atoms in all prmtops
        atom_numbers = np.empty(len(self.listOfTrajs), dtype=int)
        for i in range(len(atom_indicesRMSD)):
            atom_numbers[i] = len(atom_indicesRMSD[i])
          
        if (np.amin(atom_numbers) != np.amax(atom_numbers)):
            print ("The number of atoms selected isn't the same in all systems:", end="\n\n")
            for i in range(len(self.listOfTrajs)):
                print (self.prmtopList[i] + ": ", atom_numbers[i])
            sys.exit(0)
        
        
        ##Now on to calculate the actual RMSD Matrix
        
        RMSD_ij  = np.zeros(self.total_frames)
        RMSD_Tot = np.full(self.total_frames, -1, dtype=float)
        curr_frame = 1
        
        for refIx in range(len(self.listOfTrajs)):
            for frameNum in range(self.listOfTrajs[refIx].n_frames):
                if (Verbose == True):
                    if (curr_frame % VerbosityCoarseness == 0):
                        print ("All-vs-All calculation is {:4.1f}% done!".format((curr_frame*100/self.total_frames)))
                    curr_frame +=1
                for targIx in range(len(self.listOfTrajs)):
        
                    RMSD_ij[self.frameIntervals[targIx][0]:self.frameIntervals[targIx][1]] = md.rmsd(self.listOfTrajs[targIx], self.listOfTrajs[refIx], frame=frameNum, atom_indices = atom_indicesRMSD[targIx], ref_atom_indices = atom_indicesRMSD[refIx])
                
                RMSD_Tot = np.vstack((RMSD_Tot, RMSD_ij))
          
        RMSD_Tot = np.delete(RMSD_Tot, 0, 0)
        
        ##Make all diagonals 0 (MDTraj has an issue with this, making the diags a very low
        ##which is dfferent to 0.
        
        for i in range (RMSD_Tot.shape[0]):
            RMSD_Tot[i][i] = 0
        
        ##Convert from nm to Angstrom:
        
        RMSD_Tot = np.multiply (RMSD_Tot,10)

        print ("RMSD Matrix has been computed.")
        
        ##Save RMSD Array to Binary Numpy file (.npy), to be able to clusterize after generating initial RMSD matrix (Save time)
        
        if (RMSDMatrixOut is not None):
            print ("Saving RMSD matrix to file {}…".format(RMSDMatrixOut))
            np.save(RMSDMatrixOut, RMSD_Tot)

        return(RMSD_Tot)
        

    def RMSDAvAHeat(self, RMSDMatrixIn, vmax=None, SaveFigName=None, show=True):
        
        """
        This function generates (and optionally saves) a heatmap of
        the RMSD All-vs-All matrix, for easier visualization.
        
        Parameters
        ----------
        
        RMSDMatrixIn: numpy.ndarray
           The RMSD Matrix to represent.    
                    
        vmax: int
            Maximum value on the colorbar.        
        
        SaveFigName: str
            Save Generated heatmap to this file 
            
        show: bool
            Show the heatmap
            
            
        Notes
        -----
        
        frameIntervals:
            An array of shape(NumberOfSystems, 2), which delimits 
            the frames to each system. Useful for visualzing differences 
            between structures. Generated when this class is instantiated.
            
        """

        heatmap = plt.imshow(RMSDMatrixIn, cmap="jet", vmax=vmax)

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
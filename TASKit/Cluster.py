#imports
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform

import numpy as np
import matplotlib.pyplot as plt



class ClusterizeMatrix:

    """A class that will take a RMSD matrix and clusterize its values,
        keeping as attributes the cluster numbers and what frames are in
        what cluster and the most representative frame for each cluster.
       
    Parameters
    ----------
    
    RMSDMatrix: numpy.ndarray
        The RMSD All-vs-All array, formatted as a numpy.ndarray.

    FClusterCutoff: float
        Value to be used cutoff for fcluster function. (scipy.cluster.hierarchy.fcluster)
    
    verbosity: bool, default=True  
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

    def __init__(self, RMSDMatrix, FClusterCutoff, verbosity=True):
        reduced_distances = squareform(RMSDMatrix, checks=False) # Memory error
        Z = sch.linkage(reduced_distances, method='weighted') # 'single', 'complete', 'weighted', and 'average'
        flatClusts = sch.fcluster(Z, FClusterCutoff, criterion='distance')

        flatClusts = np.array(flatClusts)
        s = set()
        for a in range (flatClusts.shape[0]):
            s.add(flatClusts[a])
        
        print ("Number of Clusters: " + str(len(s)))            
            
        self.framesInClust = len(s) * [None]
        for j in range(len(s)):
            if (verbosity == True):
                print ("In cluster " + str(j+1) + " are the following frames: " , end='')
            self.framesInClust[j] = []
            for i in range (flatClusts.shape[0]):
                if (flatClusts[i] == list(s)[j]):
                    if (verbosity == True):
                        print (i, end = ", ")
                    self.framesInClust[j].append(i)
            if (verbosity == True):
                print()


        # Iterate through clusters
        self.mins = np.empty((len(self.framesInClust)), dtype=int)
        for clustIx in range(len(self.framesInClust)):
            # Iterate through frames matrix
            distM = np.empty((len(self.framesInClust[clustIx]), len(self.framesInClust[clustIx])))
            for Ix in range(len(self.framesInClust[clustIx])):
                for Jx in range(len(self.framesInClust[clustIx])):
                    i = self.framesInClust[clustIx][Ix]
                    j = self.framesInClust[clustIx][Jx]
                    distM[Ix][Jx] = RMSDMatrix[i][j]
            # Teodor score
            avgs = [np.mean(distM[i]) for i in range(distM.shape[0])]
            self.mins[clustIx] = self.framesInClust[clustIx][np.argmin(avgs)]

        print ("The min of each cluster, RMSD-wise:")
        print (self.mins)
        self.RMSDMatrix = RMSDMatrix
        
    def GenerateHeatmap(self, FigOut=None, show=True, vmax=None):
        clusterColors = ['g', 'c', 'y', 'k', 'w']
        clusterMarkers = ['o', 'v', 's', 'p', '*', 'h', 'D']

        ## First we have to plot the RMSD AvA
        heatmap = plt.imshow(self.RMSDMatrix, cmap="viridis", vmax=vmax)
        colormap = plt.colorbar(heatmap)
 
        ## Then we add colored points to mark the clusters

        for i in range (len(self.framesInClust)):
            for ClusterPoint in (self.framesInClust[i]):
                plt.plot(ClusterPoint, ClusterPoint, clusterMarkers[int(i/5)] + clusterColors[i%5], label="Cluster " + str(i+1))
    
        for ClusterCenterIx in range(len(self.mins)):
            plt.plot(self.mins[ClusterCenterIx], self.mins[ClusterCenterIx], "xm")

        plt.title("All-vs-All RMSD Heatmap ($\AA\,$) + Cluster Data", fontsize=20)
       
        ## Now we add a legend, so we know what point is in what cluster. Note:
## in order to have a non-redundant legend, we use a dict (since there can't
## be 2 of the same key in a dict object

        handles,labels = plt.gca().get_legend_handles_labels()
        test = dict(zip(labels,handles))
        plt.legend(test.values(),test.keys())
        
        
        ## Then just save and/or display
 
        if (FigOut is not None):
            plt.savefig(FigOut + ".png", dpi=800)   
    
        if (show == True):
            plt.show() 

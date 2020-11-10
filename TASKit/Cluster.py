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
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform

import numpy as np
import matplotlib.pyplot as plt

import plotly.express as px
import plotly.graph_objects as go

import pandas as pd

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
        If true, then print cluster-frame info, not just the mins of all frames 
    
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

    def GenerateHeatmap(self, FrameIntervals, show=True, saveHTML=None):


        """
        This function generates (and optionally saves) a heatmap of
        the RMSD All-vs-All matrix, for easier visualization.
        
        Parameters
        ----------

        FrameIntervals: numpy.ndarray
            Numpy array of shape (NumberOfSystems, 2), to visualise frame
            distribution between clusters. It's an attr of the Calc Object
            (RMSDAvA Module)

		  FigOut: str
			Save Generated heatmap to this file       
 
        vmax: int
            Maximum value on the colorbar.        
        
        show: bool
            Show the heatmap
         
        saveHTML: str
            If not None, saves RMSD Heatmap as an interactive HTML object, to this file. 

        PrmtopList: str
            If not None, in the Cluster TreeMap, instead of "structure #", the corresponding
            "structure name" will be used
             
        """
        


        ##Generate a PANDAS Dataframe, to be able to use Plotly easily

        Column_names_clusters = ["No. Of Frames"]

        clusters = []
        clusterSize = []
        for i in range (len(self.framesInClust)):			
            clusters.append("Cluster " + str(i+1)) 
            clusterSize.append(len(self.framesInClust[i]))
 
        Pandas_Data = {"Cluster": clusters, "Number Of Frames": clusterSize, "Frames": self.framesInClust, "Representative Frames": self.mins}

        Pandas_clusters = pd.DataFrame(data=Pandas_Data)
 
        print (Pandas_clusters)
 
        ##HTML Part
       
        fig = px.treemap(Pandas_clusters,
              names="Cluster",
              path=[clusters],
              values="Number Of Frames"
              )

        if (saveHTML is not None):
            fig.write_html(saveHTML)
        if (show == True):
            fig.show()

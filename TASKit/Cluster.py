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

    def __init__(self, RMSDMatrix, FClusterCutoff, ClusterGoal, MaxTrials=10, verbosity=True):
        reduced_distances = squareform(RMSDMatrix, checks=False) # Memory error
        Z = sch.linkage(reduced_distances, method='weighted') # 'single', 'complete', 'weighted', and 'average'
        if (ClusterGoal is not None):
            ClusterNumber = None
            Trials = 0
            Increment = 1
            while (Increment > 0.001):
                while (ClusterNumber != ClusterGoal and Trials < MaxTrials):
                    
                    flatClusts = sch.fcluster(Z, FClusterCutoff, criterion='distance')
    
                    flatClusts = np.array(flatClusts)
                    listOfClusters = set()
                    for a in range (flatClusts.shape[0]):
                        listOfClusters.add(flatClusts[a])
                    ClusterNumber = len(listOfClusters)
                
                    #print ("Number of Clusters: " + str(len(listOfClusters)))
                    if (len(listOfClusters) > ClusterGoal):
                       	NewClusterCutoff = FClusterCutoff + Increment 
                        if (NewClusterCutoff < 0):
                            Increment = Increment/10
                            NewClusterCutoff = FClusterCutoff - Increment
                          
                        print("Number of clusters ({}) larger than Cluster Goal ({}). Switching to Cluster Cutoff: {}".format(len(listOfClusters),ClusterGoal,NewClusterCutoff))	 
                        Trials += 1  
                        FClusterCutoff = NewClusterCutoff
                    if (len(listOfClusters) < ClusterGoal):
                       	NewClusterCutoff = FClusterCutoff - Increment 
                        print("Number of clusters ({}) smaller than Cluster Goal ({}). Switching to Cluster Cutoff: {}".format(len(listOfClusters),ClusterGoal,NewClusterCutoff))	 
                        Trials += 1              
                        FClusterCutoff = NewClusterCutoff
         
                    #print(Increment)
                    #print(Trials)

                    if (len(listOfClusters) == ClusterGoal):
                        Increment = 0
                        print ("Cluster Number reached! ({}). FClusterCutoff value is: {}".format(len(listOfClusters), FClusterCutoff))
                #print(Increment)
                Trials = 0
                Increment = Increment/10
            if (len(listOfClusters) != ClusterGoal):
                print ("Increment limit reached. Closest number of clusters reached: {}. \
Restart with this FClusterCutoff to reach desired number of clusters.".format((len(listOfClusters))))

        else:
            flatClusts = sch.fcluster(Z, FClusterCutoff, criterion='distance')
            
            flatClusts = np.array(flatClusts)
            listOfClusters = set()
            for a in range (flatClusts.shape[0]):
                listOfClusters.add(flatClusts[a])
            ClusterNumber = len(listOfClusters)
            
            print ("Number of Clusters: " + str(len(listOfClusters)))
            

        self.framesInClust = len(listOfClusters) * [None]
        for j in range(len(listOfClusters)):
            if (verbosity == True):
                print ("In cluster " + str(j+1) + " are the following frames: " , end='')
            self.framesInClust[j] = []
            for i in range (flatClusts.shape[0]):
                if (flatClusts[i] == list(listOfClusters)[j]):
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

#            print (distM)
            
            # MDTraj Score
#            beta = 1
#            index = (np.exp(-beta*distM[i] / distM[i].std()).sum(axis=1) for i in range (distM.shape[0]))
#            print (index)
#            highestSim = np.argmax(index)
#            self.mins[clustIx] = index
 

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

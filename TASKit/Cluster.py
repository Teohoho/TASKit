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

def GetLimitsOfConsequentialNumbers(ListIn,Stride=1):
    """ A function that gets limits of consequential number sequences from a large, non-continuous
        sequence of ordered numbers i.e. :

        [0,1,2,3,6,7,8,9] => [[0,3],[6,9]]
 
    Parameters
    ----------

    ListIn: list 
        List to be parsed
    
    Stride: int
        Minimum value that separates two consecutive elements of the list


    Returns 
    ----------

    ListOut:
        Resulting nested list

    """

#    for elementIx in range(len(ListIn)):
#        if ListIn[elementIx] == ListIn[elementIx+1]-Stride

    return "None"

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
        
    def GenerateHeatmap(self, FrameIntervals, FigOut=None, show=True, vmax=None, saveHTML=None):


        """
        This function generates (and optionally saves) a heatmap of
        the RMSD All-vs-All matrix, for easier visualization.
        
        Parameters
        ----------

        FrameIntervals: numpy.ndarray
            Numpy array of shape (NumberOfSystems, 2), to visualise frame
            distribution between clusters.

		  FigOut: str
			Save Generated heatmap to this file       
 
        vmax: int
            Maximum value on the colorbar.        
        
        show: bool
            Show the heatmap
         
        saveHTML: str
            If not None, saves RMSD Heatmap as an interactive HTML object, to this file. 

             
        """
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



        ##HTML Part
        if (saveHTML is not None):
            fig = px.imshow(self.RMSDMatrix, color_continuous_scale="viridis", zmax=vmax)

            for i in range (len(self.framesInClust)):
#               for ClusterPoint in (self.framesInClust[i]):
                print (self.framesInClust[i])
                fig.add_trace(go.Scatter(x=[self.framesInClust[i], self.framesInClust[i]], 
                                     y=[self.framesInClust[i], self.framesInClust[i]], 
                                     name="Cluster " + str(i+1),
                                     marker_symbol=100,
                                     mode="markers",
                                     fill="toself"))

            for ClusterCenterIx in range(len(self.mins)):
                fig.add_trace(go.Scatter(x=[self.mins[ClusterCenterIx]], 
                                         y=[self.mins[ClusterCenterIx]], 
                                         marker=dict(size=12, line=dict(width=4, color='DarkSlateGrey')),
                                         marker_symbol=4))

			

            ##Add frameIntervals lines
            if (FrameIntervals is not None):
	            for i in range(self.frameIntervals.shape[0]-1):
	                CutoffPoint = self.frameIntervals[i][1]-0.5
	            ##Vertical Lines
	                fig.add_trace(go.Scatter(x=[CutoffPoint, CutoffPoint],
	                                         y=[0,self.frameIntervals[-1][1]],
	                                         hoverinfo=None,
	                                         showlegend=False,
	                                         mode='lines',
	                                         line=dict(color="black", width=4)))
	            ##Horizontal Lines
	                fig.add_trace(go.Scatter(x=[0,self.frameIntervals[-1][1]],
	                                         y=[CutoffPoint, CutoffPoint],
	                                         hoverinfo=None,
	                                         showlegend=False,
	                                         mode='lines',
	                                         line=dict(color="black", width=4)))
	 
            fig.update_layout(legend=dict(
                yanchor="top",
                xanchor="right"))
            
            #fig.write_html(saveHTML)
            fig.show()

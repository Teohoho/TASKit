#__init__.py

##Import non-TASKit modules

import simtk.openmm.app
import simtk.openmm
import simtk.unit

import mdtraj as md
import numpy as np
import sys

import matplotlib.pyplot as plt

import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform

import pymbar.timeseries

##Print hello message

print ("TASKit: A Python Toolkit for Post-MD analysis. Written by Șulea Teodor Asvadur")

##Import all TASKit modules

from TASKit import Decorrelation
from TASKit import Cluster
from TASKit import RMSDAvA
from TASKit import TrajHandle

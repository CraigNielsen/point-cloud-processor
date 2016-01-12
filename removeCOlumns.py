'''
Created on 26 Mar 2015

@author: Craig
'''
import CraigIOonly as cio                                                                                                                                                                         
import combinedLeastSqr as cbl  
import numpy as np                                                                                                                                                                                            
import scipy.spatial as ss         
from numpy import arccos

if __name__ == '__main__':
    pointCloud=cio.readIn("c:/scripts/PCAfield.xyz",6,False) #can set origin within this function
    cio.WritePointstoTxt("C:/scripts/PCAfieldML.xyz", pointCloud)  
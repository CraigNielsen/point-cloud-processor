'''
Created on 26 Mar 2015

@author: Craig Ferguson
'''

#!/usr/bin/env python                                                                                                                                                                                         
#coding:utf-8                                                                                                                                                                                                 
                                                                                                                                                                                                                
from itertools import combinations
from itertools import permutations
import time
import sys
from PyQt4 import QtGui
from PyQt4 import *
from sympy import *

import CraigIOonly as cio                                                                                                                                                                         
import combinedLeastSqr as cbl  
import numpy as np                                                                                                                                                                                            
import scipy.spatial as ss         
from numpy import arccos


def getUnitVector(vector):
    u=0
    for i in vector:
        u+=i**2 
    u=sqrt(u)
    vector=vector/u

    c=0
    for i in vector:
        if (i<1e-15 and i>-1e-15):
            vector[c]=0
        c+=1
    return vector
def getKDtree(matrix):

    x=[]
    y=[]
    z=[]
    for i in Range(len(matrix)):
        x.append(matrix[i,0])
        y.append(matrix[i,1])
        z.append(matrix[i,2])
        
    points = zip(x,y,z)
    tree = ss.KDTree(points)
    return tree,points
def getPlanePointsGrid():
    x,y,z = np.mgrid[0:10, 0:10, 0:1]
    points = zip(x.ravel(), y.ravel(), z.ravel())
    pointsmat=np.matrix(points)
    
    return pointsmat


def main():

    timer=True
    
    areaForNormals=1
    areaIncrement=2
    testNighbourCount=True
    showNeighboursCount=False
    
    
    np.set_printoptions(precision=6, suppress=True,linewidth=500)
#     pointCloud=getPlanePointsGrid()

    pointCloud=cio.readIn("c:/scripts/urbancloud.csv",3,True) #can set origin within this function
    size=pointCloud.shape[0]
    print "Point Cloud Size: "+str(size)+" Points"
#     rotation=cbl.rotationMatrix(0, pi/4., 0)
#     pointCloud= pointCloud*rotation

    tree,points=getKDtree(pointCloud)
    # get index of each point which distance from (x=0, y=0) is under 1
    
#     a = tree.query_ball_point([1,1,1], 2)
    
    normals=np.zeros((len(points),3))
    theta=np.zeros((len(points),1))
    s=0
    
    average=0

    
    for point in pointCloud:
        start0=time.clock()
        a= point[0,0]
        b= point[0,1]
        c= point[0,2]
        vec=[a,b,c]

        #get nearest points

        localPoints = tree.query_ball_point(vec, areaForNormals)#gives index to local points close to vec
        sz=len (localPoints)
        i=0
        if testNighbourCount:
            while sz<20:
                localPoints = tree.query_ball_point(vec, areaForNormals+(i*areaIncrement))#gives index to local points close to vec
                sz=len (localPoints)
                i+=1
        if showNeighboursCount:
            print sz
            
            
        mat=[]
        for i in localPoints:
            mat.append(points[i])
     
        #calculate the local planes normal
#         s1=time.clock()
#         normal=cbl.solvePlane(mat, [0.1,0.1,1.1])
        normal=cbl.solvePlanePCA(mat)
#         f1=time.clock()
#         print s1-f1
        
#         s1=time.clock()
        normal=getUnitVector(normal)
#         f1=time.clock()
#         print s1-f1
        #associate the normal with that point

        normals[s,0]=normal[0]
        normals[s,1]=normal[1]
        normals[s,2]=normal[2]
        a=(np.dot(normal,[0,0,1]))
        a=float(a)
        theta[s]=(np.arccos(a))*180/pi
        
        s+=1
        if timer:
            if (s<100) :
                
                finish0= time.clock()
                average=(average+(finish0-start0))/2
                
            else:
                print "time left:"
                print int( (size-s)*(average))
                

#         finish0= time.clock()
#         print finish0-start0
    pointCloud=np.append(pointCloud, normals, 1)
    pointCloud=np.append(pointCloud, theta, 1)

    #print point cloud to file
    
    cio.WritePointstoTxt("C:/scripts/PCAfield.xyz", pointCloud)    
#     normal=cbl.solvePlane(pointCloud, [0.1,0.1,1.1])
    
    
#     print sqrt(normal[0]**2+normal[1]**2+normal[2]**2)
    
    
if __name__ == "__main__":
    main()
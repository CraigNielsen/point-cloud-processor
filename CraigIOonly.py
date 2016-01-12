import csv

from numpy import float_
from sympy import Range

import numpy as np


def readIn(path,upToColumn,changeOrigin):
    if changeOrigin:
        xOrigin=56502
        yOrigin=3754100
        zOrigin=-100
    else:
        xOrigin=0
        yOrigin=0
        zOrigin=0
    with open(path, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        data = []
        for row in reader:
            data.append(row)
        
        MAT=np.matrix(data,dtype=float_)
        safe=0
        while (safe<100) :
            try:
                MAT=np.delete(MAT, upToColumn, 1)
            except:
                print "outa while loop"
                break
            safe+=1



    for i in Range((MAT.shape[0])):
        MAT[i,0]=MAT[i,0]+xOrigin
        MAT[i,1]=MAT[i,1]+yOrigin
        MAT[i,2]=MAT[i,2]+zOrigin
#         pprint (MAT)
    print MAT
    
    return MAT

def WritePointstoTxt(path,Mat):
    '''
    given a Points file,
    this function will write to file path given as path 
    '''
    np.savetxt(path,Mat, delimiter=" ")
    return True
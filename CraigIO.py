'''
Created on 19 Feb 2015

@author: Craig
IsInt function checks for an int
'''
import csv
import numpy as np
from sympy import *
import sympy as sp
from numpy import double, dtype, float_
import ObsObject as ooj
from ObsObject import obsObject
import time

    
    
    
    
def readIn(path):

    with open(path, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        data = []
        for row in reader:
            data.append(row)
        MAT=np.matrix(data,dtype=float_)
#         pprint (MAT)
 
    return MAT
def createPMatrix(size,refv_OVER_obsVariance):
    '''
    given an observations file,
    this function will generate an I matrix 
    '''

    P=np.zeros((size,size))
    for row in Range(size):
        for column in Range(size):
            if row==column:
                P[row,column]=refv_OVER_obsVariance
            else:
                P[row,column]=0

    
#         pprint (MAT)
    P=np.asmatrix(P,dtype=float)
    return P

def writeVforObs(path,obsVector):
        with open(path, 'wb') as csvfile2:               
            writer = csv.writer(csvfile2,delimiter=',',quotechar='|')
            for row in obsVector:
                writer.writerow([row[0,0],])
def writeXforObs(path,obsVector):
        with open(path, 'wb') as csvfile2:               
            writer = csv.writer(csvfile2,delimiter=',',quotechar='|')
            for row in obsVector:
                writer.writerow([row[0,0],])

def readInOBS(path,symbolsArray):
    
    with open(path, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        first_row = next(reader)
        num_cols = len(first_row)
        csvfile.seek(0)
        row_count = sum(1 for row in reader)
#         print num_cols
#         print row_count
        csvfile.seek(0)
        
        data=np.zeros(shape=(row_count,num_cols),dtype=object)
#         pprint (data)
        i=0
        for row2 in reader:
            k=0
            for column in row2:
#                 print data[i,k]
#                 t=symbols(str(symbolsArray[k])+str(i))
#                 s=ooj.obsObject(symbolsArray[k],column)  #create obs object  for each x , y  and z
                
                data[i,k]=ooj.obsObject(symbolsArray[k],column)
                k=k+1
            

            i+=1
        
#         print data[0][0].value
#         print data[1][1].symbol
        return data
def copyIntoOBSObjectMatrix(Mat,symbolsArray):
        row_count=len(Mat)
        num_cols=len(Mat[0])
        data=np.zeros(shape=(row_count,num_cols),dtype=object)
#         pprint (data)
        i=0
        for row2 in Mat:
            k=0
            for column in row2:
#                 print data[i,k]
#                 t=symbols(str(symbolsArray[k])+str(i))
#                 s=ooj.obsObject(symbolsArray[k],column)  #create obs object  for each x , y  and z
                
                data[i,k]=ooj.obsObject(symbolsArray[k],column)
                k=k+1
            

            i+=1
        
        return data
 
    

def Linearize2WRT(function,x_,y_,initialX,initialY):
    x, y= symbols(x_+' '+y_+' ')
    
    return diff(function,x)+diff(function,y)+function.subs([(x,initialX),(y,initialY)])

def Linearize3WRT(function,x,y,z,initialX,initialY,initialZ):
   
    return diff(function,x)+diff(function,y)+diff(function,z)+function.subs([(x,initialX),(y,initialY),(z,initialZ)])

def Linearize4WRT(function,x,y,z,a,initialX,initialY,initialZ,initialA):
    dx,dy,dz,da=symbols('d'+x.name + ' d'+y.name + ' d'+z.name +' d'+a.name)
    return diff(function,x)*dx + diff(function,y)*dy+diff(function,z)*dz+diff(function,a)*da+function.subs([(x,initialX),(y,initialY),(a,initialA),(z,initialZ)])

def Linearize5WRT(function,x,y,z,a,b,initialX,initialY,initialZ,initialA,initialB):
    
    dx,dy,dz,da,db=symbols('d'+x.name + ' d'+y.name + ' d'+z.name +' d'+a.name+' d'+b.name)
    return diff(function,x)*dx + diff(function,y)*dy+diff(function,z)*dz+diff(function,a)*da+function.subs([(x,initialX),(y,initialY),(z,initialZ),(a,initialA),(b,initialB)])
def GetSymbolicMatrices(func,obsVariablesArray,unknownsArray,rowsForB):
    B=sp.zeros(rowsForB*len(func), len(obsVariablesArray)*rowsForB)
#     print B
    A=sp.zeros(rowsForB*len(func), len(unknownsArray))
    V=sp.zeros(len(obsVariablesArray)*rowsForB,1)
    X=sp.zeros(len(unknownsArray),1)
    #increment by size of function array (ADDING MULTIPLE FUNCTION FUNCTIONALITY)
    #loop through fun array with counter, fill in diff/wrt for each funtion while count is less than length of function array
    increment=len(func)
    vcount=0
    #range loop doesnt work so making a variable to use for indexing correct row.. i will be as if func had only 1 funtion
    row=0
    for i in Range(rowsForB):
        for j in Range(len(obsVariablesArray)):
            for level in Range(len(func)):
                B[row+level,vcount]=diff(func[level],obsVariablesArray[j])
            t=symbols('v'+obsVariablesArray[j].name+str(i+1))              #not sure if must be row as well
            V[vcount,0]=t
            vcount+=1
        row+=increment
    #             pprint (B)
    row=0   
    for i in Range(rowsForB):
        for j in Range(len(unknownsArray)):
            for level in Range(len(func)):
                A[row+level,j]=diff(func[level],unknownsArray[j])
            t=symbols('d'+unknownsArray[j].name)
            X[j,0]=t
        row+=increment
    return B,V,A,X

def isINT(x):
    try:
        int(x)
        
        return true
    except TypeError:
        return false
    
def fillMatrixB(obs,b,subsGroup,functionsArrayLength):    #groupForSubs=[[m_p,b_p],bSymbols,aSymbols]
    '''Given the Observations 2D array, B symbolic matrix, A symbols, x symbols, subsgroup
    This will return the numeric B Matrix'''
    
#     B=np.zeros(b.shape[0], b.shape[1])
    i=0
    j=0
    increment=functionsArrayLength
    rowCount=0
    #read fillMatrixA for explaination on things here:
    totalRows=(b.shape[0])
    c=b[:,:]                                #'''check this test'''
    print "total rows: "+ str(totalRows)
#     start12=time.clock()
    for row in Range(b.shape[0]/increment):
        
        
#         print "iteration: " + str(row)
#         print start12

        
#         pprint (b)                                                                         # C is a temperary row vector copy of symbolic vector B. (can change this if full b matrix is given into function
#         print B.shape[0]                                                                    #/copy of b to use for subs
        for column in Range(b.shape[1]):
            h=0
            for item in subsGroup[1]:#B SYmbols obs
#                 print c[row,j]
                if (isINT(c[row,j])):                                                       #if symbol in symbolic matrix is a number, skip over
                    continue
#                 print item
#                 print obs[i,h]
                k=obs[i,h]
                for level in Range(increment):
                    c[rowCount+level,j]=c[rowCount+level,j].subs([(item,k)])#                                        JUST DEBUG For an update of where you are>>> next function is likely iterate through A
#                 print c[row,j]
                h+=1
            h=0
            for item in subsGroup[2]:#A symbols unknowns/provisionals
#                 print c[row,j]
                if (isINT(c[row,j])):                                                       #if symbol in symbolic matrix is a number, skip over
                    continue
#                 print item
#                 print subsGroup[0][h]
                k=subsGroup[0][h]
                for level in Range(increment):
                    c[rowCount+level,j]=c[rowCount+level,j].subs([(item,k)])                                        
#                 print c[row,j]
                h+=1
            j+=1

        
        i+=1
        j=0
        rowCount+=increment
        fin12=time.clock()
        
#         print "ETA: " + str( (fin12 -start12)*totalRows-fin12 )
#     for j1 in Range(B.shape[1]):
#             for level in Range(increment): 
#                 B[rowCount+level,j1]=c[rowCount+level,j1]
#     t3=time.clock()
#     print "Fill B took: " + str(t3-start12)
    print type(c)
    c=np.matrix(c)
    c=c.astype(float_)
    print c.shape
    print type(c)
    print type(c[1,1])
    return c.astype(float_)
            

            
def fillMatrixA(obs,b,v,a,x,subsGroup,functionarrayLength):    
    #____________________________________
    
    A=np.zeros((a.shape[0], a.shape[1]))
#     pprint (A)
    i=0
    j=0
    #again, create a new row variable called increment that increments the row for larger b matrix
    #row now will just be a count for obs, each obs has to be subbed for each function hence (no of functions * obs) for num of rows
    increment=functionarrayLength
    rowCount=0
    for row in Range(A.shape[0]/increment):#/2 because we are using the row as if it were a single function (but it isnt)
        c=a[:,:]     
#         pprint (c)                                                                       # C is a temperary row vector copy of symbolic vector B. (can change this if full b matrix is given into function
#         print A.shape[0]                                                                    #/copy of b to use for subs
        for column in Range(A.shape[1]):
            h=0
            for item in subsGroup[1]:#B SYmbols
#                 print c[row,j]
                if (isINT(c[row,j])):                                                       #if symbol in symbolic matrix is a number, skip over
                    continue
#                 print item
#                 print obs[i,h]
                k=obs[i,h]
                #for each function, sub in values
                for level in Range(functionarrayLength):
                    c[rowCount+level,j]=c[rowCount+level,j].subs([(item,k)])#    sub B symbol for observation                                    JUST DEBUG For an update of where you are>>> next function is likely iterate through A
#                     print c[rowCount+level,j]
        
                h+=1
            h=0
            for item in subsGroup[2]:#A symbols
#                 print c[row,j]
                if (isINT(c[row,j])):                                                       #if symbol in symbolic matrix is a number, skip over
                    continue
#                 print item
#                 print subsGroup[0][h]
                k=subsGroup[0][h]
                for level in Range(functionarrayLength):
                    c[rowCount+level,j]=c[rowCount+level,j].subs([(item,k)])                                        
#                     print c[rowCount+level,j]
                h+=1
            j+=1
        
        for j1 in Range(A.shape[1]):
            for level in Range(functionarrayLength): 
                A[rowCount+level,j1]=c[rowCount+level,j1]
#         pprint(A)
        i+=1
        j=0
        rowCount+=increment
    A=A.astype(float_)
    print A.shape
    pprint (A)
    print dtype(A[1,1])
    return A      
    
def getWeightMat(covarianceMatrixB,Matrix_B):
        W=Matrix_B*covarianceMatrixB*Matrix_B.T
#         print W**(-1)
        return W**(-1)

      
def misclosureMatrix(obs,function,provObject_Vector): 
    W=[]
    i=0
    increment=len(function)
    for i in Range(len(obs)):
        for j in Range(increment):
            W.append(function[j])

#     pprint (W)
    #go through each observation
    #for each row, replace the function symbols with the provisional values,(iterate through prov values) 
    #and the observed value (iterate through obs values)
    #then add that symbolically replaced value into the W vector and repeat 
    #need a copy of function to use for replacements, and this is refreshed upon each iteration
    i=0
    j=0
    #create a fulll copy of w, and substitute like before where
    rowCount=0 
    for row in Range(len(W)/increment):
#         print row
        c=[]
        for level in Range(increment):
            c.append(W[rowCount+level])  
#         pprint (func2)                                                                       # C is a temperary row vector copy of symbolic vector B. (can change this if full b matrix is given into function
#         print len(W)  
        
        for object2 in provObject_Vector:  # prov values                              #/copy of b to use for subs
            for level in Range(increment):
                c[level]=c[level].subs(object2.symbol,object2.value)
#             print c
#             print object2.symbol
#             print object2.value
        h=0
#         print c
        
        for object in obs[row]:
            for level in Range(increment):
                c[level]=c[level].subs(object.symbol,object.value)
#             print c
#             print object.symbol
#             print object.value
 
#         print c
        
        
        for level in Range(increment):
            W[rowCount+level]=c[level]
#         i+=1    
#     pprint (W)
        rowCount+=increment
#         pprint(W)
    tempW=np.matrix(W)
    
    
    tempW=(tempW.T).astype(float_)
    return tempW
'''
-0.099
0.458
-0.484
0.131
'''
         
def chi_test(SampleVariance,PopulationVariance,dof,significance):
        from scipy.stats import chi2
        ChiSqr=dof * (SampleVariance/PopulationVariance)
        testChi=chi2.ppf(1-significance,dof)
        print "test stat: " + str(testChi)
        print "Chi2 :" + str(ChiSqr)
        if ChiSqr>testChi:
            print "reject Ho:\n    :    S>V \n"
        if ChiSqr<=testChi:
            print "Accept Ho:\n    :    S<=V \n"
    
    
    
    
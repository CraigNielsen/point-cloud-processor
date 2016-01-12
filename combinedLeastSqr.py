'''
Created on 10 Feb 2015

@author: Craig
'''
import CraigIO as Craig
from CraigIO import *
from sympy import *
from sympy import latex
import sympy as sp
import provisionalObject
from provisionalObject import provObject
import pprint
import time
from numpy import asmatrix
# from ctypes.test.test_python_api import c_py_ssize_t
    
#     Step 1: Compute the initial approximations.
#     Step 2: Develop the appropriate matrices.
#     Step 3: Solve the system.
#     Step 4: Apply the corrections to m0 and b0 to update the initial approximations for the second iteration
#             {-    Second Iteration. During the second iteration, only the unknown parameters are                -}
#             {-    updated, and thus only B, We, and K matrices differ from their first iteration counterparts.  -}
#     find equivalent residual vector
#     find the observational residuals
def getW_Plane(Matrix,provs):
    W=np.zeros((Matrix.shape[0],1))
    a=provs[0]
    b=provs[1]
    c=provs[2]
    
    for i in Range(Matrix.shape[0]):
        x=Matrix[i,0]
        y=Matrix[i,1]
        z=Matrix[i,2]
        W[i]=-1-(a*x + b*y + c*z)
    W=asmatrix(W,dtype=float_)
    return W 
def rotationMatrix(x,y,z):
    Rx=np.matrix([[1,0,0],[0,cos(x),-sin(x)],[0,sin(x),cos(x)]])
    Ry=np.matrix([[cos(y),0,sin(y)],[0,1,0],[-sin(y),0,cos(y)]])
    Rz=np.matrix([[cos(z),-sin(z),0],[sin(z),cos(z),0],[0,0,1]])
    
    R=Rx*Ry*Rz
    return R
def solvePlanePCA(mat):
    mat=np.matrix(mat)
    cov_mat=np.cov(mat.T)
    eigen_value,eigen_vector=np.linalg.eig(cov_mat)
#     print eigen_value,eigen_vector
#     for i in range(len(eigen_value)):
#         print('Eigenvector {}: \n{}'.format(i, eigen_vector[i]))
#         print('Eigenvalue {}: \n{}'.format(i, eigen_value[i]))
    eig_pairs = [(np.abs(eigen_value[i]), eigen_vector[:,i]) for i in range(len(eigen_value))]
    eig_pairs.sort()
    
    a=(np.dot(eig_pairs[0][1],[0,0,1]))
    a=float(a)
    theta=(np.arccos(a))*180/pi

    if theta>90:
        return -eig_pairs[0][1]
    return eig_pairs[0][1]

def solvePlane(MatofPoints,provsArray): # in linear plane, matofPoints is the A matrix in a quasi parametric case, see notes :7-13
    MatofPoints=np.matrix(MatofPoints)
    A=MatofPoints

    l=getW_Plane(MatofPoints,provsArray)
    try:
        X=(A.T*A)**-1*(A.T*l)
    except:
        
        print "encounted singular..setting x to zero"
        return np.asarray(provsArray, dtype=float_)
    for i in Range(len(provsArray)):
        provsArray[i]=provsArray[i]+X[i,0]



    provsArray=np.asarray(provsArray, dtype=float_)
    return provsArray
def process(MatOfPoints,provsArray):
    
    '''if not quasi, need to disable the check int function in the FillBmatrix function'''
    '''matOfPoints is a matrix with each obs variable in its own column,eg a point per row'''
    '''DONT FORGET TO CHANGE 'GROUPSFORSUBS' IN ITERATOR
    '''
    
    
    pp = pprint.PrettyPrinter(width=80,indent=4)    

    quasi=False
#     quasi=False
    obsVariance=0.09**2
    refVariance=0.09**2
    #obsvariance assumes equal weights for all obs
    

##______________________________________________________________________   Transform  ___________________________________________________________________________
    '''______________________________________________________________________________________________________________________________________________'''
    
#     a,b,c,d,x,y,X,Y=symbols('a b c d x y X Y')
#     func=[a*x-b*y+c-X,b*x+a*y+d-Y]   
#     OBS=Craig.readIn('F:/Scripts/transform.csv')
#     OBS2=Craig.readInOBS('F:\\scripts\\transform.csv', [x,y,X,Y])   #the updated read for an obs object, can access data without issues in the obs container..chilled
#     a_p=25.386458
#     b_p=-0.8158708
#     c_p=-137.216
#     d_p=-150.600
#     bSymbols=[x,y,X,Y]
#     aSymbols=[a,b,c,d]  #unknowns (with provisionals)
#     provs=[provObject(a, a_p),provObject(b, b_p),provObject(c, c_p),provObject(d, d_p)]
#     groupForSubs=[[a_p,b_p,c_p,d_p],bSymbols,aSymbols]                           #provisionals....B symbols(obs substituts)>>>>>A Symbols(provisonal substitutes)
#     P=Craig.createIdentityMatrix('F:\\scripts\\transform.csv',bSymbols,len(func))                    #only because data file had root values that needed to be sqaured)**-1
#     
    
##______________________________________________________________________   LINE  ___________________________________________________________________________
    '''______________________________________________________________________________________________________________________________________________'''
    
#     x,y,m,b=symbols('x y m b')
#     func=[(y) -(m*(x)) - b  ] 
#     OBS=Craig.readIn('C:/Scripts/test.csv')
#     OBS2=Craig.readInOBS('C:\\scripts\\test.csv', [x,y])   #the updated read for an obs object, can access data without issues in the obs container..chilled
#     m_p=0.246
#     b_p=3.663
#     bSymbols=[x,y]
#     aSymbols=[m,b]#unknowns (with provisionals)
#     provs=[provObject(m, m_p),provObject(b, b_p)]
#     groupForSubs=[[m_p,b_p],bSymbols,aSymbols]                           #provisionals....B symbols(obs substituts)>>>>>A Symbols(provisonal substitutes)
#     covarMatB=Craig.readIn('C:\\scripts\\covarianceB.csv')
#     P=(covarMatB*covarMatB.T )                      #only because data file had root values that needed to be sqaured)**-1
#     P=np.matrix(P)
#     P=P.astype(float_)
    
##______________________________________________________________________   PLANE  ___________________________________________________________________________
    '''______________________________________________________________________________________________________________________________________________'''
    
    x,y,z,a,b,c,d=symbols('x y z a b c d')
    func=[a*(x)+b*(y)+c*(z) - d ] 
    OBS=MatOfPoints#numpy matrix of columns of unknown observations. x in column 1, y in column 2 etc. 
    OBS2=Craig.copyIntoOBSObjectMatrix(OBS, [x,y,z])   #the updated read for an obs object, can access data without issues in the obs container..chilled
    
    a_p=provsArray[0]
    b_p=provsArray[1]
    c_p=provsArray[2]
    d_p=provsArray[3]
    bSymbols=[x,y,z]
    aSymbols=[a,b,c,d]#unknowns (with provisionals)
    provs=[provObject(a, a_p),provObject(b, b_p),provObject(c, c_p),provObject(d, d_p)]
    groupForSubs=[[a,b,c,d],bSymbols,aSymbols]                           #provisionals....B symbols(obs substituts)>>>>>A Symbols(provisonal substitutes)
#     covarMatB=Craig.readIn('C:\\scripts\\covarianceB.csv')
    size=len(OBS)
    P=Craig.createPMatrix(size,1)
    
##______________________________________________________________________ CIRCLE ___________________________________________________________________________
    '''__________________________________________________________________________________________________________________________________________________'''

#     x,y,r,xo,yo=symbols('x y r xo yo')
#     func=[(x-xo)**2 + (y-yo)**2 - (r)**2] 
# #     covarMatB=Craig.readIn('C:\\scripts\\covarianceB.csv')
#     OBS=Craig.readIn('C:\\scripts\\CircleValues.csv')
#     OBS2=Craig.readInOBS('C:\\scripts\\CircleValues.csv', [x,y])   #the updated read for an obs object, can access data without issues in the obs container..chilled
#     x_p=200
#     y_p=300
#     r_p=50        
#     bSymbols=[x,y]
#     aSymbols=[xo,yo,r]
#     provs=[provObject(xo, x_p),provObject(yo, y_p),provObject(r,r_p)]
#     groupForSubs=[[x_p,y_p,r_p],bSymbols,aSymbols]                           #provisionals....B symbols(obs substituts)>>>>>A Symbols(provisonal substitutes)
#     P=Craig.createIdentityMatrix('C:\\scripts\\CircleValues.csv',bSymbols,len(func))
#     P=np.matrix(P)
#     P=P.astype(float_)

##______________________________________________________________________ SPHERE ___________________________________________________________________________
    '''__________________________________________________________________________________________________________________________________________________'''

#     x,y,z,r,xo,yo,zo=symbols('x y z r xo yo zo')
#     func=[(x-xo)**2 + (y-yo)**2 +(z-zo)**2 - (r)**2] 
# #     covarMatB=Craig.readIn('F:\\scripts\\covarianceB.csv')
#     OBS=Craig.readIn('C:\\scripts\\sphereObsx.csv')
#     OBS2=Craig.readInOBS('C:\\scripts\\sphereObsx.csv', [x,y,z])   #the updated read for an obs object, can access data without issues in the obs container..chilled
#     x_p=-0.000218
#     y_p=0.00206703
#     z_p=0.03214974
#     r_p=3.08396034    
# 
# #     x_p=20000
# #     y_p=30000
# #     z_p=5000
# #     r_p=50       
#     bSymbols=[x,y,z]
#     aSymbols=[xo,yo,zo,r]
#     provs=[provObject(xo, x_p),provObject(yo, y_p),provObject(zo,z_p),provObject(r,r_p)]
#     groupForSubs=[[x_p,y_p,z_p,r_p],bSymbols,aSymbols]                           #provisionals....B symbols(obs substituts)>>>>>A Symbols(provisonal substitutes)
#     P=Craig.createPMatrix('C:\\scripts\\sphereObsx.csv',bSymbols,len(func),refVariance/obsVariance)
#     P=np.matrix(P)
#     P=P.astype(float_)
#     pp.pprint (P)

##__________________________________________________________________ Compute Symbolic Matrices _____________________________________________________________
    '''__________________________________________________________________________________________________________________________________________________'''

    B_,V_,A_,X_= GetSymbolicMatrices(func, bSymbols, aSymbols,OBS.shape[0]) #get Matrix. hand in a function, and tell what unknowns are needed etc
    pp.pprint(B_)
    pp.pprint(A_)

    start9=time.time()
    A = fillMatrixA(OBS, B_, V_, A_, X_,groupForSubs,len(func))
    fin9=time.time()
    print 'fillAmatrix took: '
    print fin9-start9
#     pp.pprint(A)
#   ______________________________________ . . .  ITERATE . . . ____________________________________________________________________________________
    '''__________________________________________________________________________________________________________________________________________________'''
    print "CALCULATING..."
    start0=time.clock()
    
    
    
    
    for k in Range(1):
        print "iteration: "+ str(k)
#         m_p=provs[0].value
#         b_p=provs[1].value
        try:
            groupForSubs=[[a,b,c,d],bSymbols,aSymbols]          #Plane
#             groupForSubs=[[m_p,b_p],bSymbols,aSymbols] #line
#             groupForSubs=[[x_p,y_p,r_p],bSymbols,aSymbols]    #circle
#             groupForSubs=[[x_p,y_p,z_p,r_p],bSymbols,aSymbols]     #sphere
    #         groupForSubs=[[a_p,b_p,c_p,d_p],bSymbols,aSymbols]     #transform
        except(NameError):
            print "ERROR: please change Groupforsubs in iterator..."
        print "start B matrix Filling"

        B = fillMatrixB(OBS, B_,groupForSubs,len(func)) #can optimaze further by not sub for every value in B, just update the provisionals
#         pp.pprint(B)

        W=misclosureMatrix(OBS2,func,provs)
#         pp.pprint(W)
        BT=B.T

        if (quasi==True):
            print B.shape
#             print "qausi started.."
#             st=time.time()
#             WL=(B*BT)
#             ft=time.time()
#             print "time taken for * is " + str(ft-st)
# #             st1=time.time()
# #             KL=np.dot(B,BT)
# #             ft1=time.time()
# #             print "time taken for dot is " + str(ft1-st1)
#             
#             
#             
#             for row in Range( WL.shape[0]):
#                 print "row " + str(row) + " / " + str(WL.shape[0])
#                 for column in Range(WL.shape[1]):
#                         if row==column:
#                             WL[row,column]=1/WL[row,column]
            print "optimized for sphere"
            WL = np.zeros((B.shape[0],B.shape[0]))
            j=0
            Pob=refVariance/obsVariance
            for k in Range(B.shape[0]):
                WL[k,k]=1./(B[k,j]**2/Pob+B[k,j+1]**2/Pob+B[k,j+2]**2/Pob)
                j+=3

        else:
            WL=(B*P**-1*B.T)**-1
        
        AT=A.T
        WL=np.asmatrix(WL)
#         pp.pprint(WL)
        X=-(AT*WL*A)**-1*(AT*WL*W)

        i=0
        for object in provs:
            object.value=object.value+X[i]
            i+=1
# 
#         Vl=A*X+W
# 
#         Vobs=BT*WL*Vl



        
#         pp.pprint (Vobs)
        i=0
        k=0


        
    finish0= time.clock()
#_________________________________________ . . .
    print "time elapsed:"
    print finish0-start0
    print "______\n"
        
  
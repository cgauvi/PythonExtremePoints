'''
Created on May 29, 2017

@author: charles

'''
import math as m
import numpy as np
from scipy.special import comb



def convertCodifiedToNumeric(codifiedList):
    '''
    Convert the points codified with  {1,2,3,4} to the actual varrho value
    '''
    dim         = len( codifiedList[0,:] )
    numPoints   = len( codifiedList[:,0] )
    epsilon =1 
    
    numericList = np.zeros((numPoints,dim))
    
    for p in range(numPoints):
        for i in range(dim):
            if codifiedList[p,i] == 0 : 
                numericList[p,i] = -    m.sqrt(dim*epsilon)
            elif codifiedList[p,i] == 1 : 
                numericList[p,i] =      m.sqrt(dim*epsilon)
            elif codifiedList[p,i] == 2 : 
                numericList[p,i] = -    m.sqrt(dim*epsilon) * (  m.sqrt(dim) - m.floor( m.sqrt(dim)) )
            elif codifiedList[p,i] == 3 : 
                numericList[p,i] =      m.sqrt(dim*epsilon) * (  m.sqrt(dim) - m.floor( m.sqrt(dim)) )
            else :          #codifiedList[p,i] == 4 
                numericList[p,i] = 0
                
                
    return numericList
                
                
                
                
def checkIfInBudgetedSet(codifiedList):
    '''
    Check if the generated point actually respects the norm 1 and norm inf constraints
    '''
    
    tolerance = m.pow(10, -6)
    
    numericList     =  convertCodifiedToNumeric(codifiedList)
    dim             = len(numericList[0,:])
    numPoints       = len( numericList[:,0] )
    epsilon         = 1
    
    numPointsInFacets = 0
     
    #First check the norm inf constraints
    #|x_i| <= m.sqrt( dim * epsilon) 
    for p in range(numPoints):
        oneOfNormInfBinding = 0
        for i in range(dim):
            if m.fabs( numericList[p,i] ) > m.sqrt( dim * epsilon) + tolerance :
                print("\nError, point{}, abs value of element: {} > {} (norm inf violated)\n".format(p,
                                                               m.fabs( numericList[p,i] ), 
                                                               m.sqrt( dim * epsilon))      \
                      )
                return 0
            
            if m.fabs( numericList[p,i] ) == m.sqrt( dim * epsilon):
                oneOfNormInfBinding =1 
                
        if oneOfNormInfBinding:
            numPointsInFacets += 1

    #Next, check the norm 1 constraint holds
    #sum_i |x_i| <= dim * m.sqrt(epsilon) 
    for p in range(numPoints):
        partialSum =0
        for i in range(dim):
            partialSum += m.fabs( numericList[p,i] ) 
        if partialSum > dim * m.sqrt(epsilon) + tolerance:
            print("\nError, point{},sum of abs values: {} > {} (norm 1 violated)\n".format(p,
                                                                                           partialSum, 
                                                                                           dim * m.sqrt(epsilon))      \
                  )
            return 0

    print("Ok all the points are within the budget (up to tolerance {})\n".format(tolerance));
    print("There are {} points that lie in a facet of the polyhedron out of the total {} \n".format(numPointsInFacets,numPoints) );
    
    
    return 1

def returnSubList(myList,mySize):
    '''
    Recursive method to get all n_C_k combinations of size k from a list of size n
    Return a ->sublist<- (not a scalar) of myList of size mySize
    '''
    
    for i in range(len(myList)):
        if mySize == 1:
            yield (myList[i],)
        else:
            for next in returnSubList(myList[i+1:len(myList)], mySize-1):
                yield (myList[i],) + next

    
   
   
   

def createVectorCodifiedExtPoints(listOfListsAtBounds, numCols):
    '''
    For each of the L_C_floor(sqrt(l)) combinations 
        take 1) floor(sqrt(l)) at ub or lb [+- sqrt(L eps)]
        &    2) 1 at +-(L sqrt(eps) - floor(sqrt(l)) sqrt(L eps) )
    for some give eps > 0
    Codification:   0 - LB         (< 0)
                    1 - UB         (> 0)
                    2 - Fractional (< 0)
                    3 - Fractional (> 0)
                    4 - 0
    Hence, there should be a total of: (this is bullshit, not smart enough to get exact #)
                                        n C k * 2^k * (n-k) * 2 
    Where n=dim and k = floor(sqrt(dim))
    '''
     
 
    
    n=numCols
    k=m.floor(m.sqrt(numCols))
    numRows= np.int64(  comb(n, k,exact=False) * m.pow(2, k) * (n-k) * 2     )
    
    codifiedList=np.ones((numRows,numCols) ) *4     #Defalt value of 4 means ith dim of ext point is fixed at 0
    codifiedListCount=0;
    
    posNegBounds=(0,1)
    posNegFrac  =(2,3)
    
    bigList=[i for i in range(numCols)]
    
    for l in listOfListsAtBounds :          #Go through each of the n C k combos
        
        elem1   = l[0]
        elem2   = l[1]
        
        for v1 in posNegBounds:
            for v2 in posNegBounds:
                listNot = [i for i in bigList if i not in l ]     
                for elemp in listNot:  #len(listNot) = n-k
                    for pp in posNegFrac:
                        codifiedList[codifiedListCount,elem1] = v1  #UB/LB
                        codifiedList[codifiedListCount,elem2] = v2  #Elem at frac value
                        codifiedList[codifiedListCount,elemp] = pp  #Elem at frac value
                        codifiedListCount += 1
 
    print("Number of ext points: {}".format(codifiedListCount) )
    
    return codifiedList
 

def printAMPLFormat(codExtPoints):
    '''
    Print codified matrix to AMPL format
    '''
    numRows= len( codExtPoints[:,0] )
    numCols= len( codExtPoints[0,:] )
    
    outFile = open("C:/Users/charles/git/ExtPoints/ExtPoints/ExtPoints/outputFil.txt", 'w')
    
    
    outFile.write('param p_codifiedMatrixExtPoints : 1  2 3 4 5 6 7   := \n' );
    for e in range(numRows):
        outFile.write( "{}\t".format( np.int64(e+1) ))
        for i in range(numCols):
             outFile.write( "{}\t".format( np.int64(codExtPoints[e,i]) ) )
        outFile.write( "\n")   
    outFile.write( ";\n")          

    
    #outFile.write("\n ".join(map(str,codExtPoints)))
    
    
    outFile.close()

if __name__ == '__main__':
    
    dim=7;
    k=m.floor(m.sqrt(dim))
    
    entireList=[i for i in range(dim)]
    
    #Save the values else magic generator object will erase them or i don't know what
    lOfL=[]
    for i in returnSubList(entireList,k):
        lOfL.append(i)
    
    
    print (", ".join(map(str,  lOfL ) ))
   
    
    codExtPoints=createVectorCodifiedExtPoints(lOfL, dim)
  
  
    #Print the vector
    printAMPLFormat(codExtPoints)
    
    
    #Check if the points generated are really within the budgeted uncertainty set
    checkIfInBudgetedSet(codExtPoints)
    
     
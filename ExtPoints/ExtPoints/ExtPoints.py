'''
Created on May 29, 2017

@author: charles

'''
import math as m
import numpy as np
from scipy.special import comb

 

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
    
    outFile = open("C:/Users/charles/GIT/ExtPointsV2GIT/ExtPointsV2/ExtPointsV2/outputFil.txt", 'w')
    
    
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
  
  
    #print the vector
    printAMPLFormat(codExtPoints)
    
     
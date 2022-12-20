#import itertools
import math
# import time, random, sys, os
# from   pathlib import Path
# from builtins import True False
# import math, time, random, sys, os
# import numpy as np, pandas as pd
# import pickle

# from printf import printf
import F2P, CEDAR

# printAllVals (cntrSize=8, hyperSize=2, hyperMaxSize=None, mode='F2P', verbose=[]):
# if __name__ == '__main__':

    # offsetOfExpVal = [0, 8, 24]
    # expVal  = list(offsetOfExpVal).index(0)
    # print ('offsetOfExpVal={}, expVal={}' .format (offsetOfExpVal, expVal))
    
# Below is an example of using cntrMaster
myCntrMaster = F2P.CntrMaster(cntrSize=6, hyperSize=1, hyperMaxSize=1, mode='F3P', numCntrs=2)
for i in range (5):
    cntrVal = print (myCntrMaster.incCntr(cntrIdx=0, delta=1))
for i in range (10):
    cntrVal = print (myCntrMaster.incCntr(cntrIdx=0, delta=10))

myCntrMaster = CEDAR.CntrMaster(cntrSize=6, delta=0.2, numCntrs=2)
for i in range (5):
    cntrVal = print (myCntrMaster.incCntr(cntrIdx=0, value=1))
for i in range (10):
    cntrVal = print (myCntrMaster.incCntr(cntrIdx=0, value=10))





# myF2PCntrMaster.num2cntr (33)
# myF2PCntrMaster.cntr2num ("11001111")
# printAllCntrMaxVals (mode='F3P')
# myF2PCntrMaster.num2cntr (num=2.3)
import math
import pickle
import random

import numpy as np

from printf import printf

# verbose constants
VERBOSE_COUT_CNTRLINE = 1
VERBOSE_DEBUG = 2
VERBOSE_RES = 3
VERBOSE_PCL = 4
VERBOSE_DETAILS = 5
VERBOSE_NOTE = 6

# Generates a strings that details the counter's settings (param vals).
genSettingsStr = lambda cntrSize, delta: 's{}_d{}'.format(cntrSize, delta)



class CntrMaster(object):
    """
    Generate, check and parse counters
    """
    # This is the CEDAR formula to calculate the diff given the delta and the sum of the previous diffs
    calc_diff = lambda self, sum_of_prev_diffs: (1 + 2 * self.delta ** 2 * sum_of_prev_diffs) / (1 - self.delta ** 2)

    # print the details of the counter in a convenient way
    printCntrLine = lambda self, cntrSize, delta, numCntrs, mantVal, cntrVal: print(
        'cntrSize={}, delta={}'
        .format(cntrSize, delta))

    cntr2num = lambda self, i: self.counterscntr[i]

    def __init__(self, cntrSize=8, delta=0.1, numCntrs=1, verbose=[]):
        """
        Initialize an array of cntrSize counters. The cntrs are initialized to 0.
        Inputs:
        cntrSize  - num of bits in each counter.
        Delta - the max relative error
        numCntrs - number of counters in the array.
        """
        self.cntrSize = cntrSize
        self.delta = delta
        self.numCntrs = numCntrs
        self.counters = [0 for i in range(self.numCntrs)]
        self.shared_estimators, self.different = self.cedar()

    def cedar(self) -> tuple:
        shared_estimators = []
        different = []
        shared_estimators.append(0)
        i = 0
        while i < 2 ** self.cntrSize:
            # using the cedar's formula
            different.append(self.calc_diff(sum_of_prev_diffs=shared_estimators[i]))
            i += 1
            shared_estimators.append(shared_estimators[i - 1] + different[i - 1])
        shared_estimators.pop()
        return shared_estimators, different

    def incCntr(self, cntrIdx: int, value: int = 1, verbose: list = []) -> dict:
        """
                Add value to the current value of the counter.
                Input:
                cntrIdx - the index in the counters array.
                value - can be either positive/negative (default: 1).
                verbose - determines which data will be written to the screen.
                Output:
                cntrDict: a dictionary representing the modified counter where:
                    - cntrDict['cntr'] is the counter's binary representation; cntrDict['val'] is its value.
                Operation:
                If cntr+delta > maximum cntr's value, return the cntr representing the max possible value.
                If cntr+delta < 0, return a cntr representing 0.
                If cntr+delta can be represented correctly by the counter, return the exact representation.
                Else, use probabilistic cntr's modification.
                """
        for i in range(value):
            # The probability to increment is calculate according to the diff
            if random.uniform(0, 1) < 1 / self.different[self.counters[cntrIdx]]:
                self.counters[cntrIdx] += 1

        return {'cntr': np.binary_repr(self.counters[cntrIdx], self.cntrSize), 'val': self.counters[cntrIdx]}

    def queryCntr(self, cntrIdx=0) -> dict:
        """
        Query a cntr.
        Input:
        cntrIdx - the counter's index.
        Output:
        cntrDic: a dictionary, where:
            - cntrDict['cntr'] is the counter's binary representation; cntrDict['val'] is its value.
        """
        self.checkCntrIdx(cntrIdx)
        return {'cntr': np.binary_repr(self.counters[cntrIdx], self.cntrSize), 'val': self.counters[cntrIdx]}


def printAllVals(cntrSize=8, delta=0.1, verbose=[]):
    """
    Loop over all the binary combinations of the given counter size.
    For each combination, print to file the respective counter, and its value.
    The prints are sorted in an increasing order of values.
    """
    listOfVals = []
    myCntrMaster = CntrMaster(cntrSize=cntrSize,delta=delta, numCntrs=1)
    for i in range(2 ** cntrSize):
        cntr = np.binary_repr(i, cntrSize)  # cntr = str("{0:{fill}8b}".format(i, fill='0'))
        val = myCntrMaster.incCntr(cntrIdx=0, value=1)
        listOfVals.append({'cntr': cntr, 'val': val})

    # listOfVals = sorted(listOfVals, key=lambda item: item['val'])
    if VERBOSE_RES in verbose:
        outputFile = open(
            '../res/{}.res'.format(genSettingsStr(myCntrMaster.cntrSize,myCntrMaster.delta), 'w'))
        for item in listOfVals:
            printf(outputFile, '{}={}\n'.format(item['cntr'], item['val']))


def printAllValsOfDiff(ced):
    print(len(ced.different))
    for x in ced.different:
        print(x)


def printAllValsOfSharedEstimators(ced):
    print(len(ced.shared_estimators))
    for x in ced.shared_estimators:
        print(x)


# if __name__ == '__main__':
#     ced = CntrMaster(4, 0.5)
#     printAllValsOfSharedEstimators(ced)
#     ced.incCntr(0, 5)

# import itertools
import math
# import time, random, sys, os
# from   pathlib import Path
# from builtins import True False
import math, time, random, sys, os
from printf import printf
import numpy as np, pandas as pd
import pickle

# verbose constants
VERBOSE_COUT_CNTRLINE = 1
VERBOSE_DEBUG = 2
VERBOSE_RES = 3
VERBOSE_PCL = 4
VERBOSE_DETAILS = 5
VERBOSE_NOTE = 6
# Generates a strings that details the counter's settings (param vals).    
genSettingsStr = lambda mode, cntrSize, hyperSize=None, hyperMaxSize=None: '{}_n{}_h{}'.format(mode, cntrSize,
                                                                                               hyperSize if mode == 'F2P' else hyperMaxSize)

# returns the value of a number given its offset, exp and mant
valOf = lambda offset, mantVal, expVal: int(offset + mantVal * 2 ** expVal)


def printBinVec(binVec, grp=4):
    """
    format-print a binary-vec. Bits are grouped into grp-sized groups.
    Example:  
    > printBinVec ("000011110", 4)
    Will print:
    0000 1111 0
    """
    L = []
    for i, b in enumerate(binVec):
        if b == "1":
            L.append("1")
        else:
            L.append("0")
        if (i + 1) % grp == 0:
            L.append(" ")
    print('{}'.format("".join(L)))


class CntrMaster(object):
    """
    Generate, check and parse counters
    """

    # Given an exponent E, calculate the exponent range to which this exponent belongs
    calc_rangeOfExpVal = lambda self, expVal: max(
        [j for j in range(len(self.expRange)) if self.expRange[j] <= expVal]) if expVal > 0 else 0

    # Calculate the maximum feasible hyper-exp size
    calcHyperMaxSize = lambda self: math.floor((self.cntrSize - 1) / 2)

    # print the details of the counter in a convenient way
    printCntrLine = lambda self, cntr, expVec, expVal, mantVal, cntrVal: print(
        'hyperVec={}, expVec={}, bias={}, exp={}, mantVec={}, mant={} \nmantMinSize={}, offset={}, val={}'
        .format(cntr[0:self.hyperSize], expVec, self.bias, self.mantVec, mantVal, self.mantMinSize,
                self.offsetOfExpVal[expVal], cntrVal))

    # Given the values of the mantissa and the exponent, returns the binary cntr representing them - when the mode is F3P.
    mantNexpVals2cntr = lambda self, mantVal, expVal: self.mantNexpVals2cntrF2P(mantVal, expVal) if (
                self.mode == 'F2P') else self.mantNexpVals2cntrF3P(mantVal, expVal)

    # Given the vector of the exponent, calculate the value it represents 
    expVec2expVal = lambda self, expVec: self.biasOfExpSize[self.expSize] - (
        int(expVec, base=2) if self.expSize > 0 else 0)

    # Given the value of the exponent, return the exponent vector representing this value 
    expVal2expVec = lambda self, expVal, expSize: np.binary_repr(num=self.biasOfExpSize[expSize] - expVal,
                                                                 width=expSize) if expSize > 0 else ""

    # Returns the maximum value that may be represented by this cntr. 
    calcCntrMaxVal = lambda self: self.calcCntrMaxValF2P() if (self.mode == 'F2P') else self.calcCntrMaxValF3P()

    # Returns the maximum value of the counter with its current params 
    cntrMaxVal = lambda self: self.cntrMaxVal

    # Returns the counter that reaches the max value  
    cntrMaxVec = lambda self: self.cntrMaxVec

    def calcExpRanges(self):
        """
        Calculate the ranges of the exponent (E_0, E_1, ...)
        """
        self.expRange = np.zeros(self.expMaxSize + 1, 'uint32')
        for j in range(0, self.expMaxSize):
            self.expRange[j + 1] = int(sum([2 ** (i) for i in range(self.expMaxSize - j, self.expMaxSize + 1)]))

    def mantNexpVals2cntrF2P(self, mantVal, expVal):
        """
        Given the values of the mantissa and the exponent, returns the binary cntr representing them - when the mode is F2P.
        """

        mantSize = self.mantSizeOfExpVal[expVal]
        expSize = self.cntrSize - self.hyperSize - mantSize
        return np.binary_repr(num=expSize, width=self.hyperSize) + self.expVal2expVec(expVal=expVal,
                                                                                      expSize=expSize) + np.binary_repr(
            num=mantVal, width=mantSize)

    def mantNexpVals2cntrF3P(self, mantVal, expVal):
        """
        Given the values of the mantissa and the exponent, returns the binary cntr representing them - when the mode is F3P.
        """

        hyperSize = self.hyperMaxSize - self.rangeOfExpVal[expVal]
        expSize = hyperSize

        return '1' * hyperSize + \
               ('0' if hyperSize < self.hyperMaxSize else '') + \
               self.expVal2expVec(expVal=expVal, expSize=expSize) + \
               np.binary_repr(num=mantVal, width=self.mantSizeOfExpVal[expVal])  # mantissa
        # print ('expVal={}, self.expRange={}, self.rangeOfExpVal={}, expSize={}, cntr={}' .format (expVal, self.expRange, self.rangeOfExpVal, expSize, cntr))

    def calcOffsets(self):
        """
        Pre-calculate all the offsets to be added to the counters, according to its exponent values.
        """
        self.calcExpRanges()
        self.calcMantSizeOfExpVal()
        self.offsetOfExpVal = np.zeros(self.bias + 1,
                                       'uint64')  # self.offsetOfExpVal[j] will hold the value to be added to the counter when the exponent is j
        for expVal in range(1, self.bias + 1):  # for each potential exponent value
            for i in range(expVal):
                self.offsetOfExpVal[expVal] += 2 ** (i + self.mantSizeOfExpVal[i])

    def calcMantSizeOfExpVal(self):
        """
        Calculate M(E), namely, the size of the mantissa implied by having each given value of exponent. 
        In particular, this function fills the array
        self.mantSizeOfExpVal, where self.mantSizeOfExpVal[e] is the size of mantissa implied by a given exponent value.
        """
        self.rangeOfExpVal = np.zeros(self.bias + 1, dtype='uint8')
        for expVal in range(self.bias + 1):
            self.rangeOfExpVal[expVal] = self.calc_rangeOfExpVal(expVal)
        if (self.mode == 'F2P'):
            self.mantSizeOfExpVal = np.array(
                [self.mantMinSize + self.rangeOfExpVal[expVal] for expVal in range(self.bias + 1)], 'uint32')
        elif (self.mode == 'F3P'):
            self.mantSizeOfExpVal = np.array(
                [self.mantMinSize + 2 * self.rangeOfExpVal[expVal] - 1 for expVal in range(self.bias + 1)], 'uint32')
            self.mantSizeOfExpVal[:self.expRange[1]] = self.mantMinSize
            # self.expSizeOfExpVal  = np.array ([self.mantMinSize + 2*self.rangeOfExpVal[exp]-1 for exp in range (self.bias+1)], 'uint32')

    def calcParams(self):
        """
        Calc the basics param, which are depended upon the counter size, and the hyper-exp' size.
        """
        self.mantMinSize = self.cntrSize - self.hyperMaxSize - self.expMaxSize
        if (self.mantMinSize < 1):
            print(
                'cntrSize={} and hyperSize={} implies min mantissa size={}. Mantissa size should be at least 1. Please use a smaller hyperSize'.format(
                    self.cntrSize, self.hyperSize, self.mantMinSize))
            exit()
        self.bias = sum([2 ** i for i in range(1, self.expMaxSize + 1)])
        self.biasOfExpSize = np.ones(self.expMaxSize + 1,
                                     'uint32')  # self.biasOfExpSize[j] will hold the bias to be added when the exp size is j
        for j in range(self.expMaxSize + 1):
            self.biasOfExpSize[j] = int(self.bias - sum([2 ** i for i in range(j)]))
        self.calcOffsets()

    def __init__(self, cntrSize=8, hyperSize=1, hyperMaxSize=None, mode='F2P', numCntrs=2, verbose=[]):

        """
        Initialize an array of cntrSize counters at the given mode. The cntrs are initialized to 0.
        Inputs:
        cntrSize  - num of bits in each counter.
        hyperSize - size of the hyper-exp field, in bits. Relevant only for F2P counters. 
        hyperMaxSize - maximal size of the hyper-exp field, in bits. Relevant only for F3P counters.
        mode - either 'F2P', or 'F3P'.
        numCntrs - number of counters in the array.
        verbose - can be either:
            VERBOSE_COUT_CNTRLINE - print to stdout details about the concrete counter and its fields.
            VERBOSE_DEBUG         - perform checks and debug operations during the run. 
            VERBOSE_RES           - print output to a .res file in the directory ../res
            VERBOSE_PCL           = print output to a .pcl file in the directory ../res/pcl_files
            VERBOSE_DETAILS       = print to stdout details about the counter
            VERBOSE_NOTE          = print to stdout notes, e.g. when the target cntr value is above its max or below its min.
        """

        if (cntrSize < 3):
            print('error: cntrSize requested is {}. However, cntrSize should be at least 3.'.format(cntrSize))
            exit()
        self.cntrSize = int(cntrSize)
        self.numCntrs = numCntrs
        self.mode = mode
        self.verbose = verbose
        if (self.mode == 'F2P'):
            self.setHyperSizeF2P(hyperSize)
        elif (self.mode == 'F3P'):
            self.setHyperMaxSize(hyperMaxSize)
        else:
            print('error: mode {} is not supported'.format(self.mode))
        self.calcParams()
        self.calcCntrMaxVal()

        self.cntrs = [self.cntrZero for i in range(self.numCntrs)]

    def cntr2num(self, cntr, hyperSize=None, hyperMaxSize=None, verbose=None):
        """
        Convert a counter, given as a binary vector (e.g., "11110"), to an integer num.
        """

        if (verbose != None):  # if a new verbose was given, it overrides the current verbose
            self.verbose = verbose
        if (len(cntr) < 4):
            print('error: cntr given is {}. However, F3P cntr size should be at least 4.'.format(cntr))
            exit()
        if (
                len(cntr) != self.cntrSize):  # if the cntr's size differs from the default, we have to update the basic params
            self.cntrSize = len(cntr)

        if (self.mode == 'F2P'):
            return self.cntr2numF2P(cntr, hyperSize)
        elif (self.mode == 'F3P'):
            return self.cntr2numF3P(cntr, hyperMaxSize)

    def calcNprintCntr(self, cntr, expVec):
        """
        Perform the final calculation (which are common for F2P, F3P modes); calculate the counter; and print the res (if requested by the user's verbose).
        Returns the value of the cntr (as int). 
        """
        expVal = self.expVec2expVal(expVec)
        if (VERBOSE_DEBUG in self.verbose):
            if (expVec != self.expVal2expVec(expVal, self.expSize)):
                print('error: expVec={}, expVal={}, expSize={}, Back to expVec={}'.format(expVec, expVal, self.expSize,
                                                                                          self.expVal2expVec(expVal,
                                                                                                             self.expSize)))
                exit()
        mantVal = int(self.mantVec, base=2)
        cntrVal = int(self.offsetOfExpVal[expVal] + mantVal * 2 ** expVal)
        if (VERBOSE_COUT_CNTRLINE in self.verbose):
            self.printCntrLine(cntr=cntr, expVec=expVec, expVal=expVal, mantVal=mantVal, cntrVal=cntrVal)
        return cntrVal

    def idxOfLeftmostZero(self, cntr, maxIdx):
        """
        if the index of the leftmost 0 in the cntr >= maxIdx, return maxIdx.
        else, return the index of the leftmost 0 in the cntr.
        """
        if (cntr == '1' * len(cntr)):
            return maxIdx
        return min(cntr.index('0'), maxIdx)

    def setHyperSizeF2P(self, hyperSize):
        """
        Sets the size of the hyper-exponent field in F2P counters as follows.
        - Check whether the hyper-exponent field size is feasible.
        - If yes - assign the relevant "self" fields (exponent's field max-size).
        - If not - print an error msg and exit.
        """
        if (hyperSize < 1 or hyperSize > self.cntrSize - 2):
            print('error: requested hyperSize {} is not feasible for counter size {}'.format(hyperSize, self.cntrSize))
            exit()
        self.hyperSize = hyperSize
        self.hyperMaxSize = hyperSize
        self.expMaxSize = 2 ** (
            self.hyperSize) - 1  # the maximum value that can be represented by self.hyperSize bits, using standard binary representation.
        if (self.expMaxSize + self.hyperSize > self.cntrSize - 1):
            print('error: requested hyperSize {} is not feasible for counter size {}'.format(hyperSize, self.cntrSize))
            exit()

    def setHyperMaxSize(self, hyperMaxSize):
        """
        Sets the maximal size of the hyper-exponent field in F3P counters as follows.
        - Check whether the hyper-exponent field size is feasible.
        - If yes - assign the relevant "self" fields (exponent's field max-size, which is identical to hyperMaxSize).
        - If not - print an error msg and exit.
        """
        hyperMaxSize = self.calcHyperMaxSize() if (hyperMaxSize == None) else hyperMaxSize
        if (2 * hyperMaxSize > self.cntrSize - 1):
            print(
                'error: requested hyperSize {} is not feasible for counter size {}'.format(hyperMaxSize, self.cntrSize))
            exit()
        self.hyperMaxSize = hyperMaxSize
        self.expMaxSize = self.hyperMaxSize

    def cntr2numF3P(self, cntr, hyperMaxSize=None):
        """
        Convert an F3P counter, given as a binary vector (e.g., "11110"), to an integer num.
        Inputs:
        cntr - the counter, given as a binary vector. E.g., "0011"
        hyperMaxSize - maximum size of the hyper-exp field.
        Output:
        the integer value of the given cntr.    
        """
        self.updateHyperMaxSize(hyperMaxSize)

        # Extract the hyper-exponent field, and value
        self.hyperSize = self.idxOfLeftmostZero(cntr=cntr, maxIdx=self.hyperMaxSize)
        self.expSize = self.hyperSize
        if (
                self.hyperSize < self.hyperMaxSize):  # if the # of trailing max < hyperMaxSize, the cntr must have a a delimiter '0'
            expVecBegin = self.hyperSize + 1
        else:
            expVecBegin = self.hyperMaxSize
        self.mantVec = cntr[expVecBegin + self.expSize:]

        return self.calcNprintCntr(cntr=cntr, expVec=cntr[expVecBegin: expVecBegin + self.expSize])

    def cntr2numF2P(self, cntr, hyperSize):
        """
        Convert an F3P counter, given as a binary vector (e.g., "11110"), to an integer num.
        Inputs:
        cntr - the counter, given as a binary vector. E.g., "0011"
        hyperSize - size of the hyper-exponent field.
        """
        if (len(cntr) < 3):
            print('error: cntr given is {}. However, F2P cntr size should be at least 4.'.format(cntr))
            exit()
        self.updateSelfHyperSize(
            hyperSize)  # if a new hyperSize was given, override the previous self.hyperSize and update the relevant params

        self.hyperVec = cntr[0:self.hyperSize:1]
        self.expSize = int(self.hyperVec, base=2)
        self.mantVec = cntr[self.hyperSize + self.expSize:]
        return self.calcNprintCntr(cntr=cntr, expVec=cntr[self.hyperSize:self.hyperSize + self.expSize])

    def checkCntrIdx(self, cntrIdx):
        """
        Check if the given cntr index is feasible.
        If not - print error msg and exit.
        """
        if (cntrIdx < 0 or cntrIdx > (self.numCntrs - 1)):
            print('error: wrong cntrIdx. Please select cntrIdx between 0 and {}'.format(self.numCntrs - 1))
            exit()

    def queryCntr(self, cntrIdx=0):
        """
        Query a cntr.
        Input: 
        cntrIdx - the counter's index. 
        Output:
        cntrDic: a dictionary, where: 
            - cntrDict['cntr'] is the counter's binary representation; cntrDict['val'] is its value.        
        """
        self.checkCntrIdx(cntrIdx)
        return {'cntr': self.cntrs[cntrIdx], 'val': self.cntr2num(self.cntrs[cntrIdx])}

    def incCntr(self, cntrIdx=0, delta=1, verbose=[]) -> dict:
        """
        Add delta to the current value of the counter.
        Input:
        cntr - given as binary.
        delta - can be either positive/negative (default: 1).
        verbose - determines which data will be written to the screen.
        Output:
        cntrDict: a dictionary representing the modified counter where: 
            - cntrDict['cntr'] is the counter's binary representation; cntrDict['val'] is its value.
        Operation:
        If cntr+delta > maximum cntr's value, return the a cntr representing the max possible value. 
        If cntr+delta < 0, return a cntr representing 0.
        If cntr+delta can be represented correctly by the counter, return the exact representation.
        Else, use probabilistic cntr's modification.
        
        If verbose==VERBOSE_DETAILS, the function will print to stdout:
        - the target value (the cntr's current value + delta)
        - optionalModifiedCntr - an array with entries, representing the counters closest to the target value from below and from above.
          If the target value can be accurately represented by the counter, then optionalModifiedCntr will include 2 identical entries. 
          Each entry in optionalModifiedCntr is a cntrDict that consists of: 
          - cntrDict['cntr'] - the binary counter.
          - cntrDict['val']  - the counter's value.
        """
        self.checkCntrIdx(cntrIdx)
        targetVal = self.cntr2num(self.cntrs[cntrIdx]) + delta
        optionalModifiedCntr = self.num2cntr(targetVal)
        if (VERBOSE_DETAILS in verbose):
            print('targetVal={}\n, optionalModifiedCntr={}'.format(targetVal, optionalModifiedCntr))
        if (
                len(optionalModifiedCntr) == 1):  # there's a single option to modify the cntr -- either because targetVal is accurately represented, or because it's > maxVal, or < 0.
            self.cntrs[cntrIdx] = optionalModifiedCntr[0]['cntr']
        else:
            probOfFurtherInc = float(targetVal - optionalModifiedCntr[0]['val']) / float(
                optionalModifiedCntr[1]['val'] - optionalModifiedCntr[0]['val'])
            # print ('probOfFurtherInc={}' .format (probOfFurtherInc))
            self.cntrs[cntrIdx] = optionalModifiedCntr[1]['cntr'] if (random.random() < probOfFurtherInc) else \
            optionalModifiedCntr[0]['cntr']
        return {'cntr': self.cntrs[cntrIdx], 'val': self.cntr2num(self.cntrs[cntrIdx])}

    def updateSelfHyperSize(self, hyperSize):
        """
        Sets self.hyperSize, and the relevant fields (self.expMaxSize) to the input hyperSize.
        The function works as follows: 
        If a new hyperSize was given - update self.hyperSize and all the relevant parameters accordingly.
        If no new hyperSize was given, but a hyperSize was stated already upon init - return.
        If no new hyperSize was given, and no hyperSize is already stored in self --> print error msg, saying that hyperSize is not know, and exit. 
        """

        if (hyperSize == None):
            if (self.hyperSize == None):
                print(
                    'error in F2P mode: hyperSize was nor specified upon init of CntrMaster, neither when calling a function')
                exit()
            # now we know that hyperSize was already set, and the relevant params were already calculated during __init__.
        else:  # hyperSize is not None --> need to re-calculate params accordingly
            self.setHyperSizeF2P(hyperSize)
            self.calcParams()

    def updateHyperMaxSize(self, hyperMaxSize):
        """
        Sets self.hyperMaSize, and the relevant fields (self.expMaxSize) to the input hyperSize.
        If a new hyperMaxSize was given - update self.hyperMaxSize.
        If no new hyperMaxSize was given, but a hyperMaxSize was stated already upon init - return.
        If no new hyperMaxSize was given, and no hyperMaxSize is already stored in self --> set hyperMaxSize to the maximum feasible value.
        If self.hyperMaxSize was now set, set the relevant parameters (e.g., expMaxSize) accordingly. 
        """

        if (hyperMaxSize == None):
            if (self.hyperMaxSize == None):
                print(
                    'note: hyperMaxSize was nor specified upon init of CntrMaster, neither when calling a function, so I am using the largest feasible maxHyperSize')
                self.setHyperMaxSize(hyperMaxSize=self.calcHyperMaxSize())
                return
                # No new hyperMaxSize was given, but self.hyperMaxSize and the relevant params were set already upon init --> do nothing
        else:
            self.setHyperMaxSize(hyperMaxSize)

    def num2cntr(self, targetVal, hyperSize=None, hyperMaxSize=None, verbose=None):
        """
        given a target value, find the closest counters to this targetVal from below and from above.
        Output:
        - A list of dictionaries, where, at each entry, 'cntr' is the binary counter, 'val' is its integer value.
        - If an exact match was found (the exact targetVal can be represented), the list contains a single dict entry: the cntr representing this targetVal. 
        - If targetVal <= 0, the list has a single dict entry: the cntr representing 0 
        - If targetVal > maxVal that this cntr can represent, the list has a single dict entry: the cntr repesenting maxVal
        - Else, 
            The first entry in the list is the dict of the max cntr value that is < targetVal.
            The second entry is the dict of min cntr val that is > targetVal.
        """
        if (verbose != None):  # if a new verbose was given, it overrides the current verbose
            self.verbose = verbose
        if (targetVal > self.cntrMaxVal):
            print('Note: the requested cntr value {} is above the max feasible cntr for this configuratiohn'.format(
                targetVal))
            return [{'cntr': self.cntrMaxVec, 'val': self.cntrMaxVal}]
        if (targetVal < 0):
            print('Note: the requested cntr value {} is negative'.format(targetVal))
            return [{'cntr': self.cntrZero, 'val': 0}]

        if (self.mode == 'F2P'):
            self.updateSelfHyperSize(hyperSize)
        elif (self.mode == 'F3P'):
            self.updateHyperMaxSize(hyperMaxSize)

        offset = max([offset for offset in self.offsetOfExpVal if offset <= targetVal])
        expVal = list(self.offsetOfExpVal).index(offset)
        # print ('list(self.offsetOfExpVal)={}, expVal={}' .format (list(self.offsetOfExpVal), expVal))
        mantVal = math.floor(float(targetVal - offset) / float(2 ** expVal))
        cntr = self.mantNexpVals2cntr(mantVal, expVal)
        cntrVal = valOf(offset, mantVal, expVal)
        numVal = self.cntr2num(cntr=cntr, hyperSize=self.hyperSize)
        if (VERBOSE_DEBUG in self.verbose):
            if (cntrVal != numVal):
                print('error in num2cntr: cntrVal={}, but the val of the generated cntr={}'.format(cntrVal,
                                                                                                   self.cntr2num(cntr)))
                exit()
        if (numVal == targetVal):  # found a cntr that accurately represents the target value
            return [{'cntr': cntr, 'val': cntrVal}]

        # now we know that the counter found is < the target value
        if (mantVal < 2 ** self.mantSizeOfExpVal[expVal] - 1):
            cntrpp = self.mantNexpVals2cntr(mantVal + 1, expVal)
            cntrppVal = valOf(offset, mantVal + 1, expVal)
        else:
            cntrpp = self.mantNexpVals2cntr(mantVal=0, expVal=expVal + 1)
            cntrppVal = valOf(offset=self.offsetOfExpVal[expVal + 1], mantVal=0, expVal=expVal + 1)
        return [{'cntr': cntr, 'val': cntrVal}, {'cntr': cntrpp, 'val': cntrppVal}]

    def calcCntrMaxValF3P(self):
        """
        sets self.cntrMaxVal to the maximum value that may be represented by this F3P cntr. 
        """

        self.cntrZero = np.binary_repr(2 ** self.cntrSize - 2 ** (self.cntrSize - 2 * self.hyperMaxSize), self.cntrSize)
        self.cntrMaxVec = np.binary_repr(2 ** (self.cntrSize - 1) - 1,
                                         self.cntrSize)  # the cntr that reaches the highest value
        self.cntrMaxVal = self.cntr2numF3P(self.cntrMaxVec)

        if (VERBOSE_DEBUG in self.verbose):
            cntrMaxValByFormula = 2 ** self.bias * (2 ** (self.cntrSize - 1) - 1)
            for i in range(self.bias):
                cntrMaxValByFormula += 2 ** (i + self.mantSizeOfExpVal[i])

            if (cntrMaxValByFormula != self.cntrMaxVal):
                print('error: cntrMaxValByFormula={}, cntrMaxValByCnt={}'.format(cntrMaxValByFormula, self.cntrMaxVal))
                exit()

    def calcCntrMaxValF2P(self):
        """
        sets self.cntrMaxVal to the maximum value that may be represented by this F2P cntr. 
        """

        self.cntrZero = np.binary_repr(2 ** self.cntrSize - 2 ** (self.cntrSize - self.hyperSize - self.expMaxSize),
                                       self.cntrSize)  # the cntr that reaches the lowest value (zero)
        self.cntrMaxVec = np.binary_repr(2 ** (self.cntrSize - self.hyperSize) - 1,
                                         self.cntrSize)  # the cntr that reaches the highest value
        cntrMaxValByCntr = self.cntr2numF2P(self.cntrMaxVec, hyperSize=self.hyperSize)
        self.cntrMaxVal = cntrMaxValByCntr

        if (VERBOSE_DEBUG in self.verbose):
            cntrMaxValByFormula = 2 ** self.bias * (2 ** (self.cntrSize - self.hyperSize) - 1)
            for i in range(self.bias):
                cntrMaxValByFormula += 2 ** (i + self.mantSizeOfExpVal[i])

            if (cntrMaxValByFormula != cntrMaxValByCntr):
                print('error: cntrMaxValByFormula={}, cntrMaxValByCnt={}'.format(cntrMaxValByFormula, cntrMaxValByCntr))
                exit()


def printAllVals(cntrSize=8, hyperSize=2, hyperMaxSize=None, mode='F2P', verbose=[]):
    """
    Loop over all the binary combinations of the given counter size. 
    For each combination, print to file the respective counter, and its value. 
    The prints are sorted in an increasing order of values.
    """
    print('running printAllVals. mode={}'.format(mode))
    myCntrMaster = CntrMaster(cntrSize=cntrSize, hyperSize=hyperSize, hyperMaxSize=hyperMaxSize, mode=mode)
    listOfVals = []
    if (mode == 'F2P'):
        for i in range(2 ** cntrSize):
            cntr = np.binary_repr(i, cntrSize)  # cntr = str("{0:{fill}8b}".format(i, fill='0'))
            val = myCntrMaster.cntr2num(cntr, hyperSize=hyperSize)
            listOfVals.append({'cntr': cntr, 'val': val})
    elif (mode == 'F3P'):
        for i in range(2 ** cntrSize):
            cntr = np.binary_repr(i, cntrSize)  # cntr = str("{0:{fill}8b}".format(i, fill='0'))
            val = myCntrMaster.cntr2num(cntr, hyperMaxSize=hyperMaxSize)
            listOfVals.append({'cntr': cntr, 'val': val})
    else:
        print('sorry, mode {} that you chose is not supported yet'.format(mode))
        exit()
    listOfVals = sorted(listOfVals, key=lambda item: item['val'])

    if (VERBOSE_RES in verbose):
        outputFile = open(
            '../res/{}.res'.format(genSettingsStr(mode, cntrSize, hyperSize=hyperSize, hyperMaxSize=hyperMaxSize)), 'w')
        for item in listOfVals:
            printf(outputFile, '{}={}\n'.format(item['cntr'], item['val']))

    if (VERBOSE_PCL in verbose):
        with open('../res/pcl_files/{}.pcl'.format(genSettingsStr(mode, cntrSize, hyperSize=hyperSize)),
                  'wb') as pclOutputFile:
            pickle.dump(listOfVals, pclOutputFile)


def printAllCntrMaxVals(mode='F3P', hyperSizeRange=[], hyperMaxSizeRange=[], cntrSizeRange=[], verbose=[VERBOSE_RES]):
    """
    print the maximum value a cntr reach for several "configurations" -- namely, all combinations of cntrSize and hyperSize. 
    """

    if (VERBOSE_RES in verbose):
        outputFile = open('../res/maxCntrVals.txt', 'a')
    if (mode == 'F2P'):
        for hyperSize in range(1, 3):
            for cntrSize in cntrSizeRange:
                myCntrMaster = CntrMaster(mode=mode, cntrSize=cntrSize, hyperSize=hyperSize)
                printf(outputFile,
                       '{} maxCntrVal={}\n'.format(genSettingsStr(mode, cntrSize, hyperSize), myCntrMaster.cntrMaxVal))
    elif (mode == 'F3P'):
        for hyperMaxSize in hyperMaxSizeRange:
            for cntrSize in cntrSizeRange:
                myCntrMaster = CntrMaster(mode='F3P', cntrSize=cntrSize, hyperMaxSize=hyperMaxSize)
                printf(outputFile, '{} maxCntrVal={}\n'.format(
                    genSettingsStr(mode='F3P', cntrSize=cntrSize, hyperMaxSize=hyperMaxSize),
                    myCntrMaster.CntrMaxVal()))
    else:
        print('Sorry, mode {} is not supported yet'.format(mode))


# hyperSizeRange = [2, 3]
# hyperMaxSizeRange = [2, 3]
# cntrSizeRange = range(8, 9)
# printAllCntrMaxVals(mode='F2P', hyperSizeRange=hyperSizeRange, hyperMaxSizeRange=hyperMaxSizeRange,
#                     cntrSizeRange=cntrSizeRange)
# if __name__ == '__main__':

# offsetOfExpVal = [0, 8, 24]
# expVal  = list(offsetOfExpVal).index(0)
# print ('offsetOfExpVal={}, expVal={}' .format (offsetOfExpVal, expVal))

# myCntrMaster = CntrMaster(cntrSize=5, hyperMaxSize=1, mode='F3P')
# for i in range (15):
#     myCntrMaster.incCntr(cntrIdx=0, delta=3)
# myCntrMaster.num2cntr (33)
# myCntrMaster.cntr2num ("11001111")
# printAllCntrMaxVals (mode='F3P')
# myCntrMaster.num2cntr (num=2.3)

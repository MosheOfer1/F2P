import math
import random

import numpy as np


class Counter:
    def __int__(self, cntr_num: int, method: str, cnt_size: int, exp_size: int):
        """Initialise an array of counters"""
        pass

    def print_all_vals(self):
        """Prints all the values that the counter capable to count"""
        pass

    def print_max_values(self, array_of_counters):
        """Prints all the max values that all the counters in the array capable to count"""
        pass

    def read(self, current_val: int) -> int:
        """Read the current value in the counter"""
        pass

    def update(self, added_val: int, current_val: int):
        """Add the value to the counter"""
        pass


class SeadCounter(Counter):
    def __init__(self, cntr_num, method=None, cnt_size=8, exp_size=3):
        self.cnt_array = []
        self.stage = None
        self.expansion_array = None
        self.cnt_size = cnt_size
        self.exp_size = exp_size
        self.capacity = 0
        # begin with zero
        self.value_resided = 0
        if method is None:
            self.static_pre_calculate_stages()
            self.static_capacity()
        elif method == 'dynamic':
            self.dynamic_pre_calculate_stages()
            self.dynamic_capacity()

    """
    Generate, check and parse counters
    """

    # Given an exponent E, calculate the exponent range to which this exponent belongs
    calc_rangeOfExpVal = lambda self, expVal: max(
        [j for j in range(len(self.expRange)) if self.expRange[j] <= expVal]) if expVal > 0 else 0

    # Calculate the maximum feasible hyper-exp size
    calcHyperMaxSize = lambda self: math.floor((self.cntrSize - 1) / 2)

    # print the details of the counter in a convenient way
    printCntrLine = lambda self, cntr, expVec, expVal, mantVal, cntrlLine: print(
        'hyperVec={}, expVec={}, bias={}, exp={}, mantVec={}, mant={} \nmantMinSize={}, offset={}, val={}'.format(
            cntr[0:self.hyperSize], expVec, self.bias, self.mantVec, mantVal, self.mantMinSize,
            self.offsetOfExpVal[expVal], cntrVal))

    # Given the values of the mantissa and the exponent, returns the binary cntr representing them - when the mode is
    # F3P.
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

    def static_capacity(self):
        for i in range(0, (2 ** self.exp_size)):
            self.capacity += self.expansion_array[i] * (2 ** (self.cnt_size - self.exp_size))

    def dynamic_capacity(self):
        for i in range(0, self.cnt_size - 1):
            self.capacity += self.expansion_array[i] * (2 ** (self.cnt_size - i - 1))

    def static_pre_calculate_stages(self):
        self.expansion_array = [(2 ** exp) for exp in range(0, 2 ** self.exp_size)]
        expansion_array_sum = [1]
        for i in range(1, len(self.expansion_array)):
            expansion_array_sum.append(expansion_array_sum[i - 1] + self.expansion_array[i])
        self.stage = [((2 ** (self.cnt_size - self.exp_size)) * expansion_array_sum[j]) for j in
                      range(0, len(expansion_array_sum) - 1)]
        self.stage.insert(0, 0)

    def dynamic_pre_calculate_stages(self):
        expansion_array = [(2 ** exp) for exp in range(0, self.cnt_size - 1)]
        stage = [0]
        for j in range(1, self.cnt_size - 1):
            stage.append(stage[j - 1] + (expansion_array[j] * 2 ** (self.cnt_size - 1 - j)))
        self.stage = stage
        self.expansion_array = expansion_array

    def read(self, l: int):
        value_resided_str = np.binary_repr(l, self.cnt_size)
        value_of_sign = int(value_resided_str[:self.exp_size], 2)
        value_of_count = int(value_resided_str[self.exp_size:], 2)
        return self.stage[value_of_sign] + (self.expansion_array[value_of_sign] * value_of_count)

    def update(self, v: int, l: int):
        if self.read(l) + v > self.capacity:  # is it the real capacity?
            self.value_resided = 2 ** self.cnt_size - 1
        else:
            cl = np.binary_repr(l, self.cnt_size)
            if l < 0:
                cl = np.binary_repr(~l + 1, self.cnt_size)
            s0 = int(cl[:self.exp_size], 2)
            c0 = int(cl[self.exp_size:], 2)
            q = v / self.expansion_array[s0]
            r = v % self.expansion_array[s0]
            if r != 0:
                p = random.uniform(0, 1)
                if p < r / self.expansion_array[s0]:
                    self.value_resided += 1
            if 0 < q < self.stage[s0] - c0:
                self.value_resided += q
            else:
                v1 = v - (self.stage[s0] - c0) * self.expansion_array[s0]
                self.value_resided = self.stage[s0]
                self.update(v1, self.value_resided)


if __name__ == '__main__':
    sCounter = SeadCounter(2)
    sCounter.update(5, sCounter.value_resided)
    print(sCounter.read(sCounter.value_resided))


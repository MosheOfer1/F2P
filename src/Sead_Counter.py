import random

import numpy as np


class SeadCounter:

    def __init__(self):
        self.stage = None
        self.expansion_array = None
        self.__int__('static', 8, 3)

    def __int__(self, method, cnt_size, exp_size):
        self.cnt_size = cnt_size
        self.exp_size = exp_size
        self.capacity = 0
        # begin with zero
        self.value_resided = 205
        if method == 'static':
            self.static_pre_calculate_stages()
            self.static_capacity()
        elif method == 'dynamic':
            self.dynamic_pre_calculate_stages()
            self.dynamic_capacity()

    def static_capacity(self):
        for i in range(0, (2 ** self.exp_size) - 1):
            self.capacity += self.expansion_array[i] * (2 ** (self.cnt_size - self.exp_size))

    def dynamic_capacity(self):
        for i in range(0, self.cnt_size - 2):
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


sCounter = SeadCounter()
sCounter.update(5, sCounter.value_resided)
print(sCounter.read(sCounter.value_resided))

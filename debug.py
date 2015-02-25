#!/usr/bin/env python

import graphcontraction
import numpy as np

np.set_printoptions(precision=4,suppress=True,linewidth=200)
x1 = abs(np.random.rand(5,5))*.5
x1[-1,:]*.1
x2 = abs(np.random.rand(5,5))*.5+.5
X = np.concatenate((x1,x2),axis=1)
print(X)

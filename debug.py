#!/usr/bin/env python

import graphcontraction
import numpy as np

np.set_printoptions(precision=4,suppress=True,linewidth=200)
h = 2
w = 2
c = 3
x1 = abs(np.random.randn(h,w,c))*.5
x1[-1,:]*.1
x2 = abs(np.random.randn(h,w,c))*.5+.5
X = np.concatenate((x1,x2),axis=1)

print("Input data:")
for ci in range(c):
    print(X[:,:,ci])
    print("")

print(X.flatten())
print("")
gc = graphcontraction.GC_D3(.1)
gc.init_grid_adjacency(h,2*w)
gc.fit(X.reshape(h*2*w,c))

print("")
print("Representer:")
res = gc.get_representer().reshape(h,2*w,c)
for ci in range(c):
    print(res[:,:,ci])

print("")
print("Labels:")
print(gc.get_labels().reshape(h,2*w))

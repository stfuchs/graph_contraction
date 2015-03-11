#!/usr/bin/env python

import graphcontraction
import numpy as np

n = np.load("normals.npy")
h,w,c = n.shape
gc = graphcontraction.QGC_F3(.004,7)
gc.fit(h,w,n.reshape(h*w,c))
C = gc.get_representer().reshape(h,w,3)

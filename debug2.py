#!/usr/bin/env python

import graphcontraction
import numpy as np
import cv2

np.set_printoptions(precision=4,suppress=True,linewidth=200)

files = ["img160x120.npy","img320x240.npy","img640x480.npy"]

img = map(np.load, files)

def image_info(pic):
    print("Res: %s, range: [ %s, %s ]" % (pic.shape, np.min(pic), np.max(pic)))
    print(pic.dtype)
    print(pic.flags)
    print()

map(image_info, img)

def normals(pic, l=5):
    x = cv2.GaussianBlur(pic,(l*2+1,l*2+1),10.)
    diffx = (x[:,l:] - x[:,:-l])[:-l,:]
    diffy = (x[l:,:] - x[:-l,:])[:,:-l]
    diffz = np.ones_like(diffx)/742.*float(l)
    info = (pic.shape, np.min(diffx), np.max(diffx), np.min(diffy), np.max(diffy))
    denom = 1./np.linalg.norm(np.dstack([diffx,diffy,diffz]),axis=2)
    return (np.dstack([diffx*denom,diffy*denom,diffz*denom])+1.)*.5

def gradients(z,l):
    #z = cv2.GaussianBlur(z,(l*2+1,l*2+1),10.)
    yy,xx=z.shape
    diffx = np.zeros([yy,xx-2*l,2*l], z.dtype)
    diffy = np.zeros([yy-2*l,xx,2*l], z.dtype)
    for i in range(1,l+1):
        dx = 1./float(i)*(z[:,i:] - z[:,:-i])[:,:xx-l]
        diffx[:,:,(i-1)*2  ] = dx[:,  l:  ]
        diffx[:,:,(i-1)*2+1] = dx[:,l-i:-i]

        dy = 1./float(i)*(z[i:,:] - z[:-i,:])[:yy-l,:]
        diffy[:,:,(i-1)*2  ] = dy[l  :  ,:]
        diffy[:,:,(i-1)*2+1] = dy[l-i:-i,:]
    return diffx[l:-l,:,:],diffy[:,l:-l,:]

def medNormals(z,l):
    diffx,diffy = gradients(z,l)
    dx = np.median(diffx,axis=2)
    dy = np.median(diffy,axis=2)
    dz = np.ones_like(dx)/742.
    info = (z.shape, np.min(dx), np.max(dx), np.min(dy), np.max(dy))
    print("Res: %s, x:[%s, %s], y:[%s, %s]"%info)
    denom = 1./np.linalg.norm(np.dstack([dx,dy,dz]),axis=2)
    #return (np.dstack([dx*denom,dy*denom,dz*denom])+1.)*.5
    return np.dstack( [(dx*denom+1.)*.5, (dy*denom+1.)*.5, dz*denom] )

def meanNormals(z,l):
    diffx,diffy = gradients(z,l)
    dx = np.mean(diffx,axis=2)
    dy = np.mean(diffy,axis=2)
    dz = np.ones_like(dx)/742.
    info = (z.shape, np.min(dx), np.max(dx), np.min(dy), np.max(dy))
    print("Res: %s, x:[%s, %s], y:[%s, %s]"%info)
    denom = 1./np.linalg.norm(np.dstack([dx,dy,dz]),axis=2)
    return (np.dstack([dx*denom,dy*denom,dz*denom])+1.)*.5

"""
X = (np.abs(np.random.randn(6,6))*100).astype(np.int)
print(X)
N = mednormals(X,2)
print(N[:,:,0])
print(N[:,:,1])
"""

"""
nimg = map(normals, img, [1,3,5])
for x,y in zip(nimg,img):
    scaled1 = cv2.resize(x, (800,600), interpolation=cv2.INTER_NEAREST)
    scaled2 = cv2.resize(y, (800,600), interpolation=cv2.INTER_NEAREST)
    cv2.imshow("org", scaled2)
    cv2.imshow("normals",scaled1)
    cv2.moveWindow("org", 10, 50)
    cv2.moveWindow("normals", 1000, 50)
    cv2.waitKey()
"""

def pair_comp(d,th):
    return lambda p: np.abs(d.flat[p[0]] - d.flat[p[1]])<th

def adjacency(h,w,mask=None,condition=None):
    aand = lambda i,j: i and j # elementwise array and
    idx = np.arange(h*w).reshape(h,w)
    xpair = np.array( zip(idx[:,:-1].flat, idx[:,1:].flat) )
    ypair = np.array( zip(idx[:-1,:].flat, idx[1:,:].flat) )
    if mask is None and condition is None:
        return np.vstack( [xpair,ypair] )
    
    if condition is not None:
        cxpair = np.array( map(condition, xpair) )
        cypair = np.array( map(condition, ypair) )
        if mask is None:
            return np.vstack( [xpair[cxpair], ypair[cypair]] )
    if mask is not None:
        mxpair = np.array( map(aand, mask[:,:-1].flat, mask[:,1:].flat) )
        mypair = np.array( map(aand, mask[:-1,:].flat, mask[1:,:].flat) )
        if condition is None:
            return np.vstack( [xpair[mxpair], ypair[mypair]] )

    mxpair = np.array(map(aand, cxpair, mxpair))
    mypair = np.array(map(aand, cypair, mypair))
    return np.vstack( [xpair[mxpair], ypair[mypair]] )

l = 5
#N1 = normals(img[2],5)
N2 = medNormals(img[2],l)
nan_mask = (img[2] != 1)[l:-l,l:-l]
#N3 = meanNormals(img[2],5)
Nsc = cv2.resize(N2, (4*2**7,3*2**7))#, interpolation=cv2.INTER_NEAREST)
h,w,c = Nsc.shape
pairs = adjacency(h,w)#,nan_mask,pair_comp((img[2])[l:-l,l:-l],8./742.))
gc = graphcontraction.GC_F3(.0015)
#gc.init_grid_adjacency(h,w)
gc.init_adjacency(pairs,h*w)
gc.fit(Nsc.reshape(h*w,3))
C = gc.get_representer().reshape(h,w,3)
#from sklearn.cluster import KMeans, MiniBatchKMeans
#mbkm = MiniBatchKMeans(8)
#mbkm.fit(Nsc.reshape(h*w,c))
#C = (mbkm.cluster_centers_[mbkm.labels_,:]).reshape(h,w,c)

sc1 = Nsc#cv2.resize(Nsc, (640,480), interpolation=cv2.INTER_NEAREST)
sc2 = C#cv2.resize(C, (640,480), interpolation=cv2.INTER_NEAREST)
cv2.imshow("N2", sc1)
cv2.imshow("C",sc2)
cv2.moveWindow("N2", 10, 50)
cv2.moveWindow("C", 850, 50)
cv2.waitKey()
"""
cv2.imshow("N1", N1)
cv2.imshow("N2", N2)
cv2.imshow("N3", N3)
cv2.moveWindow("N1", 10, 50)
cv2.moveWindow("N2", 650, 50)
cv2.moveWindow("N3", 650, 550)
cv2.waitKey()
"""

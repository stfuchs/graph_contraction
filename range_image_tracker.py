#!/usr/bin/env python

import rospy
from sensor_msgs.msg import Image

import numpy as np
import cv2
from cv_bridge import CvBridge

from graphcontraction import GC_F3

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
    #print("Res: %s, x:[%s, %s], y:[%s, %s]"%info)
    denom = 1./np.linalg.norm(np.dstack([dx,dy,dz]),axis=2)
    return (np.dstack([dx*denom,dy*denom,dz*denom])+1.)*.5


class RangeImageTracker:
    def __init__(self):
        self.sub_img = rospy.Subscriber("/camera/depth/image", Image,
                                        self.callback, queue_size=1, buff_size=2**24)
        self.pub_img1 = rospy.Publisher("/tracking/range_image/pseudo_normals",Image,queue_size=1)
        self.pub_img2 = rospy.Publisher("/tracking/range_image/image2",Image,queue_size=1)
        self.bridge = CvBridge()
        print("Default input /camera/depth/image [sensor_msgs::Image]")
        print("Default output /tracking/range_image/image [sensor_msgs::Image]")
        self.N = np.zeros([])

    def callback(self, img_msg):
        start = rospy.Time.now()
        img = self.bridge.imgmsg_to_cv2(img_msg).copy()
        img = 1. - 348./(742.*img)
        img[np.isnan(img)] = 1.
        N = medNormals(img,5)
        if True: #self.N.shape != N.shape:
            self.N = N
        else:
            self.N = (self.N+N)*.5
        msg = self.bridge.cv2_to_imgmsg( (self.N*255.).astype(np.uint8),'bgr8')
        msg.header = img_msg.header
        self.pub_img1.publish(msg)

        #msg = self.bridge.cv2_to_imgmsg( (img_cartoon*255.).astype(np.uint8),'bgr8')
        #self.pub_img2.publish(msg)
        #self.img_prev = img_mean
        print("Process took %s sec"% (rospy.Time.now() - start).to_sec())

if __name__ == '__main__':
    import sys
    rospy.init_node('range_image_tracker')
    node = RangeImageTracker()
    rospy.spin()

#!/usr/bin/env python

import rospy
from sensor_msgs.msg import Image

import numpy as np
import cv2
from cv_bridge import CvBridge

class RangeImageTracker:
    def __init__(self):
        self.sub_img = rospy.Subscriber("/camera/depth/image", Image,
                                        self.callback, queue_size=1)
        self.pub_img1 = rospy.Publisher("/tracking/range_image/image1",Image,queue_size=1)
        self.pub_img2 = rospy.Publisher("/tracking/range_image/image2",Image,queue_size=1)
        self.bridge = CvBridge()
        print("Default input /camera/depth/image [sensor_msgs::Image]")
        print("Default output /tracking/range_image/image [sensor_msgs::Image]")
        self.img_prev = np.zeros([])

    def callback(self, img_msg):
        start = rospy.Time.now()
        img = self.bridge.imgmsg_to_cv2(img_msg).copy()
        img = 1090. - 348./img
        nan_mask = np.isnan(img)
        img[nan_mask] = 0
        #diffx1 = img[:,3:] - img[:,:-3]
        #diffy1 = img[3:,:] - img[:-3,:]
        #diff11 = 1.-np.exp(-.1*np.abs(img[1:,:] - img[:-1,:]))[:,:-1]
        #diff21 = 1.-np.exp(-.1*np.abs(img[:,1:] - img[:,:-1]))[:-1,:]
        #diff13 = 1.-np.exp(-.1*np.abs(img[3:,:] - img[:-3,:]))[:,:-3]
        #diff23 = 1.-np.exp(-.1*np.abs(img[:,3:] - img[:,:-3]))[:-3,:]
        diff15 = np.exp(-.3*np.abs(img[5:,:] - img[:-5,:]))[:,:-5]
        diff25 = np.exp(-.3*np.abs(img[:,5:] - img[:,:-5]))[:-5,:]
        #diff15[nan_mask[:-5,:-5]] = 0
        #diff25[nan_mask[:-5,:-5]] = 0
        img_new = np.dstack([np.zeros_like(diff15), diff15, diff25 ])
        #diffx = np.abs(img_new[5:,:] - img_new[:-5,:])[:,:-5]
        #diffy = np.abs(img_new[:,5:] - img_new[:,:-5])[:-5,:]
        #img_new = .5*(diffx+diffy)
        #print(img_new.min(),img_new.max())
        #img_new = 1.-np.amax(np.dstack([diff11[:-4,:-4],diff21[:-4,:-4],diff15,diff25]),axis=2)
        #diffx = 1.-np.exp(-.1*np.abs(np.diff(img,axis=0)))[:,:-1]
        #diffy = 1.-np.exp(-.1*np.abs(np.diff(img,axis=1)))[:-1,:]
        img_blur = cv2.GaussianBlur(img_new,(7,7),1.)
        img_blur[nan_mask[:-5,:-5],:] = 0
        if self.img_prev.shape != img_blur.shape:
            img_mean = img_blur
        else:
            img_mean = .5*(img_blur + self.img_prev)

        msg = self.bridge.cv2_to_imgmsg( (img_mean*255.).astype(np.uint8),'bgr8')
        msg.header = img_msg.header
        self.pub_img1.publish(msg)
        img_cartoon = cv2.bilateralFilter(img_mean,5,300,300)
        msg = self.bridge.cv2_to_imgmsg( (img_cartoon*255.).astype(np.uint8),'bgr8')
        self.pub_img2.publish(msg)
        self.img_prev = img_mean
        #print(img.min(),img.max())
        #print(diffx.shape, diffx.min(), diffx.max())
        #print(diffy.shape, diffy.min(), diffy.max())
        print("Process took %s sec"% (rospy.Time.now() - start).to_sec())

if __name__ == '__main__':
    import sys
    rospy.init_node('range_image_tracker')
    node = RangeImageTracker()
    rospy.spin()

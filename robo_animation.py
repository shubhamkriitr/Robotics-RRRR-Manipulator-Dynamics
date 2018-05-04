#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:04:08 2018

@author: shubham
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation

def create_animation (X, Y, Z, colors=["red", "blue","green", "yellow"],figsize=(10,10),
                    lim = [[-20,20],[-20,20],[-20,20]],view = (30,180)):
    """
    Args:
        X (`class` numpy array): x-coordinate value for link end
                                    points(including origin)
    """
    N = X.shape[0]-1
    fig = plt.figure(figsize=figsize)
    ax = plt.axes(projection='3d')

    ax.set_xlim3d(lim[0])
    ax.set_xlabel('X')

    ax.set_ylim3d(lim[1])
    ax.set_ylabel('Y')

    ax.set_zlim3d(lim[2])
    ax.set_zlabel('Z')

    ax.set_title('3D Test')

    ax.view_init(view[0],view[1])

    ax_arr = []
    for i in range(N):
        clr = colors[i%len(colors)]
        ax_arr.append((ax.plot([0,0],[0,0],[0,0],'o-',color=clr))[0])

    def animate(i,ax_arr, xd, yd, zd):
        if i<xd.shape[1]:
            for j in range(N):
                ax_arr[j].set_data([xd[j,i],xd[j+1,i]],[yd[j,i],yd[j+1,i]])
                ax_arr[j].set_3d_properties([zd[j,i],zd[j+1,i]])
        return ax_arr

    anim = animation.FuncAnimation(fig, animate, fargs = (ax_arr, X, Y, Z),
                               frames=100, interval=100, blit=True)
    return anim



if __name__ == "__main__":
    X = np.linspace(0,20,800)
    Y = 20*np.sin(np.linspace(0,20,800))
    Z = 20*np.cos(np.linspace(0,20,800))
#    X = np.zeros(800)
#    Y = np.zeros(800)
#    Z = np.zeros(800)
#    for i in range(4):
#        X[i*200:(i+1)*200] = i*5
#        Y[i*200:(i+1)*200] = i*5
#        Z[i*200:(i+1)*200] = i*7

    X = X.reshape((8,100))
    Y = Y.reshape((8,100))
    Z = Z.reshape((8,100))
    #anim1 = create_animation(X,Y,Z)
    anim2 = create_animation(Y,Z,X,view=(60,-90))
    plt.show()
    anim2.save('/home/shubham/Dropbox/COURSES/Robotics Works/robo_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

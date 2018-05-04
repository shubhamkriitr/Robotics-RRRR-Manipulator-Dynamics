#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 22:26:33 2018

@author: shubham
"""

from __future__ import division
from sympy import *
from IPython.display import display
import robo_utils as ru
get_tf_mat = ru.get_tf_mat

#init_printing(use_latex='mathjax')


PI = 3.14159265
d2r = PI/180
L = symbols('dummy L1 L2 L3 L4')
Lc = symbols('dummy L_c1 L_c2 L_c3 L_c4')
th = symbols('dummy theta1 theta2 theta3 theta4')
#display(L)
#display(th)

DHT = [[L[1], -90*d2r, 0, th[1]],
       [L[2], +90*d2r, 0, th[2]],
       [L[3], -90*d2r, 0, th[3]],
       [L[4],    0   , 0, th[4]]]




T01 = get_tf_mat(DHT[0])
T12 = get_tf_mat(DHT[1])
T02 = T01*T12
T23 = get_tf_mat(DHT[2])
T03 = T02*T23
T34 = get_tf_mat(DHT[3])
T04 = T03*T34
#print("T01")
#display(T01)
#print("T02")
#display(T02)
#print("T03")
#display(T03)
#print("T04")
#display(T04)


Z0 = Matrix([0,0,1])
Z1 = T01[0:3,2]
Z2 = T02[0:3,2]
Z3 = T03[0:3,2]
Z4 = T04[0:3,2]


P0 = Matrix([0,0,0])
P1 = T01[0:3,3]
P2 = T02[0:3,3]
P3 = T03[0:3,3]
P4 = T04[0:3,3]

Zs = [Z0, Z1, Z2, Z3, Z4]
Ps = [P0, P1, P2, P3, P4]

#COMPUTE JACOBIANS
J1 = zeros(6,4)
J2 = zeros(6,4)
J3 = zeros(6,4)
J4 = zeros(6,4)

Js = [None, J1, J2, J3, J4]
Jvc = [None]
Jw


for v in range(1,N_links+1):
    for i in range(1,v+1):
        Js[v][0:3,i-1] = Zs[i-1].cross(Ps[v]-Ps[i-1])
        Js[v][3:6,i-1] = Zs[i-1]
    Jvc.append((Js[v][0:3,0:N_links]).subs({L[v]:Lc[v]}))

#ROTATIONAL INERTIA
w = symbols('w')
a = -w/2
b = +w/2
I1 = ru.get_inertia_tensor([0, L[1]],[a,b],[a,b])
I2 = ru.get_inertia_tensor([0, L[2]],[a,b],[a,b])
I3 = ru.get_inertia_tensor([0, L[3]],[a,b],[a,b])
I4 = ru.get_inertia_tensor([0, L[4]],[a,b],[a,b])

I1 = I1.subs({'w':1,L[1]:6,"rho":10})










# coding: utf-8

# In[41]:


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 13 22:26:33 2018

@author: shubham
"""

from __future__ import division
from sympy import *
from IPython.display import display

init_printing(use_latex='mathjax')

N_links = 4
PI = 3.14159265
d2r = PI/180
L = symbols('duh L1 L2 L3 L4')

Lc = symbols('dummy L_c1 L_c2 L_c3 L_c4')
th = symbols('duh  theta1 theta2 theta3 theta4')

display(L)
display(th)

DHT = [[L[1], -90*d2r, 0, th[1]],
       [L[2], +90*d2r, 0, th[2]],
       [L[3], -90*d2r, 0, th[3]],
       [L[4],    0   , 0, th[4]]]

def get_tf_mat (paramlist):
    a = paramlist[0]
    ap = paramlist[1]
    d = paramlist[2]
    th = paramlist[3]
    cosap = cos(ap)
    sinap = sin(ap)
    if (abs(cos(ap))<1e-5):
        cosap = 0
        if ap<0:
            sinap = -1
        else:
            sinap = 1
    T = [[cos(th), -sin(th)*cosap, sin(th)*sinap, a*cos(th)],
         [sin(th), cos(th)*cosap, -cos(th)*sinap, a*sin(th)],
         [  0,        sinap,           cosap,        d     ],
         [  0,           0,                 0,           1]]
    T = Matrix(T)
    return T

#display(DHT)
#display(get_tf_mat(DHT[0]))
T01 = get_tf_mat(DHT[0])
T12 = get_tf_mat(DHT[1])
T02 = T01*T12
T23 = get_tf_mat(DHT[2])
T03 = T02*T23
T34 = get_tf_mat(DHT[3])
T04 = T03*T34
print("T01")
display(T01)
print("T02")
display(T02)
print("T03")
display(T03)
print("T04")
display(T04)

Tmat = [None, T01, T02, T03, T04]


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

J1 = zeros(6,4)
J2 = zeros(6,4)
J3 = zeros(6,4)
J4 = zeros(6,4)

Js = [None, J1, J2, J3, J4]

"""for v in range(1,N_links+1):
    for i in range(1,v+1):
        Js[v][0:3,i-1] = Zs[i-1].cross(Ps[v]-Ps[i-1])
        Js[v][3:6,i-1] = Zs[i-1]"""


Jvc = [None]
Jwc = [None]
for v in range(1,N_links+1):
    for i in range(1,v+1):
        Js[v][0:3,i-1] = Zs[i-1].cross(Ps[v]-Ps[i-1])
        Js[v][3:6,i-1] = Zs[i-1]
    Jvc.append((Js[v][0:3,0:N_links]).subs({L[v]:Lc[v]}))
    Jwc.append(Js[v][3:6,0:N_links])
    
    
print("Jacobian")
for i in range(1,N_links+1):
    print("J"+str(i))
    display(Js[i])
    print("Jvc"+str(i))
    display(Jvc[i])
    print("Jwc"+str(i))
    display(Jwc[i])









# In[28]:


def get_inertia_tensor (limx, limy, limz):
    #INERTIA TENSORS
    x,y,z,rho =  symbols('x y z rho')

    ixx = rho*(y*y + z*z)
    iyy = rho*(x*x + z*z)
    izz = rho*(x*x + y*y)
    ixy = -rho*x*y
    iyz = -rho*y*z
    izx = -rho*z*x

    Ixx = integrate(ixx,(x, limx[0], limx[1]),(y,limy[0],
                         limy[1]),(z,limz[0], limz[1]))
    Iyy = integrate(iyy,(x, limx[0], limx[1]),(y,limy[0],
                         limy[1]),(z,limz[0], limz[1]))
    Izz = integrate(izz,(x, limx[0], limx[1]),(y,limy[0],
                         limy[1]),(z,limz[0], limz[1]))
    Ixy = integrate(ixy,(x, limx[0], limx[1]),(y,limy[0],
                         limy[1]),(z,limz[0], limz[1]))
    Iyz = integrate(iyz,(x, limx[0], limx[1]),(y,limy[0],
                         limy[1]),(z,limz[0], limz[1]))
    Izx = integrate(izx,(x, limx[0], limx[1]),(y,limy[0],
                         limy[1]),(z,limz[0], limz[1]))

    #Ixx,Iyy,Izz,Ixy,Iyz,Izx = symbols("I_xx I_yy I_zz I_xy I_yz I_zx")

    I = Matrix([[Ixx, Ixy, Izx],
                [Ixy, Iyy, Iyz],
                [Izx, Iyz, Izz]])
    return I


# In[54]:


w = symbols('w')
I = [None]
for i in range(1,N_links+1):
    I.append(get_inertia_tensor([-L[i]/2, +L[i]/2],[-w/2,+w/2],[-w/2,+w/2]))

display(I[1])
#I1 = I[1].subs({'w':1,L[1]:6,"rho":10})
#display(I1)


# In[66]:


q = Matrix([th[i] for i in range(1,5)])
display(q)
display(Tmat[4])
display(Tmat[4][2,0].diff(q[2]))


# In[56]:


M = symbols("dummy M1 M2 M3 M4")
for i in range(1,N_links+1):
    Ri = Tmat[i][0:3,0:3]
    if i==1:
        D = M[i]*Transpose(Jvc[i])*Jvc[i] + Transpose(Jwc[i])*Ri*I[i]*Transpose(Ri)*Jwc[i]
    else:
        D = D + M[i]*Transpose(Jvc[i])*Jvc[i] + Transpose(Jwc[i])*Ri*I[i]*Transpose(Ri)*Jwc[i]

print("D")
#display(D)


# In[57]:


D.shape


# In[60]:


display(I[1])
s = I[1][0,0]
sw  =  s.diff('w')
sL  = s.diff('L1')
display(sw)
display(sL)


# In[73]:


#CALCULATE C
C = []
for k in range(N_links):
    C.append(zeros(N_links,N_links))
    for i in range(N_links):
        for j in range(i, N_links):
            C[k][i,j] = (D[k,j].diff(q[i])+D[k,i].diff(q[k]) - D[i,j].diff(q[j]))
            C[k][j,i] = C[k][i,j]


# In[76]:


C[0].shape


# In[87]:


from sympy.physics.quantum import TensorProduct
mat1 = Matrix([[1,2],[3,1]])
u = Matrix([1,1])
display(Transpose(u))
display(mat1)
ut = Matrix(Transpose(u))
tp = TensorProduct(mat1,u)
display(tp)
display(ut*tp)


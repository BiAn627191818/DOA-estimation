# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 09:34:10 2021

@author: Administrator
"""

import math
import numpy as np
import matplotlib.pyplot as plt


derad=math.pi/180
radeg=180/math.pi
twpi=2*math.pi
kelm=8
dd=0.5
d=np.arange(0,kelm*dd,dd)
d=d.reshape(1,kelm)
iwave=3
theta=np.array([[10,30,60]])
snr=10
n=500
A=np.exp(-1j*twpi*np.dot(d.T,np.sin(theta*derad)))
S=np.random.randn(iwave,n)
X=np.dot(A,S)


def wgn(x,snr,outtype):
    snr=10**(snr/10.0)
    xsum=0
    for i ,p in enumerate(x):
        xsum=xsum+np.abs(p)**2
    l=len(x)
    xpower=xsum/len(x)
    npower=xpower/snr
    if(outtype=='complex'):
        z=np.random.randn(l)
        z=z.reshape(l,1)
        W=np.dot(z,np.sqrt(npower/2))+np.dot(1j,np.dot(z,np.sqrt(npower/2)))
        W=np.array(W)
        W=W.reshape(l,1)
    else:
        W=np.random.randn(l)*np.sqrt(npower)
    return x+W


X0=X.reshape(kelm*n,1)
X1=wgn(X0,snr,'complex')

X1=X1.reshape(kelm,n)
X1_mat=np.mat(X1)
Rxx=np.dot(X1,X1_mat.H)/n
InvS=np.linalg.inv(Rxx)
EV,D1=np.linalg.eigh(Rxx)
EVA=np.sort(EV)
D=np.fliplr(D1)


SP=np.array([[]]*1)
for iang in range(0,360):
    angle=(iang-180)/2
    phim=derad*angle
    a=np.exp(-1j*twpi*np.dot(d.T,np.sin(phim)))
    a_mat=np.mat(a)
    L=iwave
    En=D[:,L:kelm]
    En_mat=np.mat(En)
    SP_iang=np.dot(a_mat.H,a)/np.dot(a_mat.H,np.dot(En,np.dot(En_mat.H,a)))
    SP=np.hstack((SP,SP_iang))
    
SP=np.abs(SP)
SPmax=np.max(SP)
SP=10*np.log(10*(SP/SPmax))
SP=np.array(SP)
SP_list=SP.tolist()


angle=np.arange(-90,90,0.5)
angle=angle.reshape(1,360)
angle_list=angle.tolist()

fig,ax=plt.subplots(1,1)
ax.plot(angle_list[0],SP_list[0],label='trend')
plt.xlim(-90,90)
plt.ylim(-100,50)
plt.title("")#标题
plt.xlabel("angle(degree)")#坐标轴标签
plt.ylabel("magnitude(dB)")

plt.show()





























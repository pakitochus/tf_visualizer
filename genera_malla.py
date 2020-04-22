#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 13:59:45 2020

@author: pakitochus
"""


import numpy as np
import bpy
N = 75
X = np.linspace(-3,3,N)
Y = np.linspace(-3,3,N)
x,y = np.meshgrid(X,Y)

def hs(x,y, th=5):
    s = (x+1j*y)
    z = abs((s**2)/(s**2+2*0.5*1*s+1))
    z[z>th] = th
    return z

z = hs(x,y)

verts = [tuple(el) for el in np.c_[x.flatten(), y.flatten(),z.flatten()]]

def makeEdges(Nr,Nc):
    r = np.arange(Nr*Nc).reshape(Nr,Nc)
    cand1 = np.c_[r.flatten()[:-1], r.flatten()[1:]]
    cand2 = np.c_[r[:-1,1:].flatten(), r[1:,1:].flatten()]
    return np.vstack((cand1,cand2))


def MakeFacesVectorized1(Nr,Nc):
    out = np.empty((Nr-1,Nc-1,2,3),dtype=int)
    r = np.arange(Nr*Nc).reshape(Nr,Nc)
    out[:,:, 0,0] = r[:-1,:-1]
    out[:,:, 1,0] = r[:-1,1:]
    out[:,:, 0,1] = r[:-1,1:]
    out[:,:, 1,1] = r[1:,1:]
    out[:,:, :,2] = r[1:,:-1,None]
    out.shape =(-1,3)
    return out

faces = MakeFacesVectorized1(N,N)
faces = [tuple(el) for el in faces]

edges = makeEdges(N,N)
edges = [tuple(el) for el in edges]

mesh = bpy.data.meshes.new("hs")
object = bpy.data.objects.new("hs",mesh)
 
#set mesh location
object.location = bpy.context.scene.cursor.location
bpy.context.collection.objects.link(object)
 
#create mesh from python data
mesh.from_pydata(verts, [], faces)#]]],[]) #[], faces)
mesh.update(calc_edges=True)
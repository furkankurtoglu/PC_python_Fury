# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 18:00:01 2020

@author: fkurtog
"""

from pyMCDS import pyMCDS
import numpy as np
import pandas as pd
from fury import window, actor, ui
import itertools
import vtk
import glob
import time

#reading data
output = 'initial.xml'

#function to create sphere actors, this function works for only CRC organoid project
mcds=pyMCDS(output)
data=mcds.get_cell_df()
Types = np.array([mcds.data['discrete_cells']['cell_type']])
Fibro = np.where(Types == 2)
Organoid = np.where(Types == 1)
#Cells
C_xpos = np.array([pd.DataFrame.to_numpy(data['position_x'][:])])
C_ypos = np.array([pd.DataFrame.to_numpy(data['position_y'][:])])
C_zpos = np.array([pd.DataFrame.to_numpy(data['position_z'][:])])
C_xyz = np.concatenate((C_xpos,C_ypos,C_zpos),axis=0)
C_xyz=C_xyz.transpose()
# Whole Cell
# Cell Radius Calculation
C_volume = np.array([pd.DataFrame.to_numpy(data['total_volume'][:])])
C_volume = C_volume/4/np.pi*3
C_radii = np.power(C_volume,1/3).transpose()
# Coloring
C_R = np.array([np.ones(len(C_radii))]).transpose()
C_G = np.array([np.ones(len(C_radii))]).transpose()
C_B = np.array([np.ones(len(C_radii))]).transpose()
C_O = np.array([np.ones(len(C_radii))]).transpose()*0.6
# Type 1 (Organoid)
C_R[Fibro[1]] = 0
C_G[Fibro[1]] = 0.353
C_B[Fibro[1]] = 1
# Type 2 (Fibroblast)
C_R[Organoid[1]] = 1
C_G[Organoid[1]] = 1
C_B[Organoid[1]] = 0
C_colors = np.concatenate((C_R,C_G,C_B,C_O),axis=1)
# Nucleus
N_xyz=C_xyz
# Nucleus Radii
N_volume = np.array([pd.DataFrame.to_numpy(data['nuclear_volume'][:])])
N_volume = N_volume/4/np.pi*3
N_radii = np.power(N_volume,1/3).transpose()
N_R = np.array([np.ones(len(N_radii))]).transpose()*0.35
N_G = np.array([np.ones(len(N_radii))]).transpose()*0.2
N_B = np.array([np.ones(len(N_radii))]).transpose()*0.1
N_O = np.array([np.ones(len(N_radii))]).transpose()*0.9
N_colors = np.concatenate((N_R,N_G,N_B,N_O),axis=1)
# Concatenations
xyz = np.concatenate((C_xyz,N_xyz),axis=0)
colors = np.concatenate((C_colors,N_colors),axis=0)
radii = np.concatenate((C_radii,N_radii),axis=0)
# Creating Sphere Actor for one time-point
sphere_actor = actor.sphere(centers=xyz,colors=colors,radii=radii)

#%%

#Creating Scene and showmanager
scene = window.Scene()
scene.set_camera(position=(71.47, 19.45, 3469.03), focal_point=(6.70, -22.29, -18.06),
                 view_up=(0.00, 1.00, -0.01))

showm = window.ShowManager(scene,
                           size=(1080, 720), reset_camera=True,
                           order_transparent=True)

#Adding Sphere actor to scene
scene.add(sphere_actor)

#xdomain_size=


#Drawing Domain Boundaries
lines = [np.array([[-2880.,-500.,-2880.],[2880.,-500.,-2880.],[2880.,500.,-2880.],[-2880.,500.,-2880.],[-2880.,-500.,-2880.],[-2880.,-500.,2880.],[-2880.,500.,2880.],[-2880.,500.,-2880.],[-2880.,500.,2880.],[2880.,500.,2880.],[2880.,500.,-2880.],[2880.,500.,-2880.],[2880.,-500.,-2880.],[2880.,-500.,2880.],[2880.,500.,2880.],[2880.,-500.,2880.],[-2880.,-500.,2880.]])]
colors = np.random.rand(1,3)
c = actor.line(lines, colors)
scene.add(c)

#Adding Dimension Labels
x_label = actor.text_3d(text='x axis (micron)',position=(-2500,-500,-3000),font_size=50,justification='left')
scene.add(x_label)
y_label = actor.text_3d(text='y axis (micron)',position=(-2666,0,-2666),font_size=50,justification='left')
scene.add(y_label)
z_label = actor.text_3d(text='z axis (micron)',position=(-3300,-500,-2600),font_size=50,justification='left')
scene.add(z_label)

#Render
showm.render()

#Start
showm.start()
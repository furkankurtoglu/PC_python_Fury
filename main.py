# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 18:00:01 2020

@author: fkurtog
"""

from pyMCDS import pyMCDS
import numpy as np
import pandas as pd
from fury import window, actor, utils, primitive, io, ui
from fury.data import read_viz_textures, fetch_viz_textures
import itertools
import vtk
import glob
import time
import random
import os
#reading data

os.chdir(r"C:\Users\fkurt\Documents\GitHub\CRC-Organoids-Multiscale-Model\PhysiCell_Model\output")
output = "output00000024.xml"

scene = window.Scene()


showm = window.ShowManager(scene,
                           size=(1600, 900), reset_camera=True,
                           order_transparent=True)



#%%
#function to create sphere actors, this function works for only CRC organoid project
mcds=pyMCDS(output)
data=mcds.get_cell_df()




Types = np.array([mcds.data['discrete_cells']['cell_type']])
Fibro = np.where(Types == 1)
#Organoid = np.where(Types == 1)
KRAS_positive = np.where(Types == 3)
KRAS_negative = np.where(Types == 0)

#%%

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
C_O = np.array([np.ones(len(C_radii))]).transpose()*1
#%%
# Type 1 KRAS Positive
C_R[KRAS_positive[1]] = 1
C_G[KRAS_positive[1]] = 0
C_B[KRAS_positive[1]] = 0

# Type 2 KRAS Negative
C_R[KRAS_negative[1]]= 1
C_G[KRAS_negative[1]] = 0
C_B[KRAS_negative[1]] = 0
# Type 3 (Fibroblast)
C_R[Fibro[1]] = 0
C_G[Fibro[1]]= 0
C_B[Fibro[1]] = 1


C_colors = np.concatenate((C_R,C_G,C_B,C_O),axis=1)


# # Nucleus
# N_xyz=C_xyz
# # Nucleus Radii
# N_volume = np.array([pd.DataFrame.to_numpy(data['nuclear_volume'][:])])
# N_volume = N_volume/4/np.pi*3
# N_radii = np.power(N_volume,1/3).transpose()
# N_R = np.array([np.ones(len(N_radii))]).transpose()*0.35
# N_G = np.array([np.ones(len(N_radii))]).transpose()*0.2
# N_B = np.array([np.ones(len(N_radii))]).transpose()*0.1
# N_O = np.array([np.ones(len(N_radii))]).transpose()*0.9
# N_colors = np.concatenate((N_R,N_G,N_B,N_O),axis=1)
# Concatenations
xyz = C_xyz
radii = C_radii
colors = C_colors
# xyz = np.concatenate((C_xyz,N_xyz),axis=0)
# colors = np.concatenate((C_colors,N_colors),axis=0)
# radii = np.concatenate((C_radii,N_radii),axis=0)





###### Method 1 Start 
# Creating Sphere Actor for one time-point
sphere_actor = actor.sphere(centers=xyz,colors=colors,radii=radii)
scene.add(sphere_actor)
###### Method 1 End

#%%

####### Method 2 Start

# # fetch_viz_textures()
# filename = read_viz_textures("Yellow_Cell.png")
# image = io.load_image(filename)

# counter = 0;
# for pos in xyz:
#     rand = random.uniform(0, 90)
#     cell_actor = actor.texture_on_sphere(image)
#     cell_actor.SetScale(radii[counter],radii[counter],radii[counter])
#     cell_actor.SetPosition(pos[0],pos[1],pos[2])
#     utils.rotate(cell_actor,(rand, 1, 0, 0))
#     counter+=1
#     scene.add(actor.texture_on_sphere(image))
#     # print(counter)

#%%


x_ring_slider = ui.RingSlider2D(center=(30, 50), initial_value=0,
                              text_template="{angle:5.1f}°",font_size=12,slider_inner_radius=20,slider_outer_radius=24,handle_outer_radius=5)


y_ring_slider = ui.RingSlider2D(center=(30, 110), initial_value=0,
                              text_template="{angle:5.1f}°°",font_size=12,slider_inner_radius=20,slider_outer_radius=24,handle_outer_radius=5)

z_ring_slider = ui.RingSlider2D(center=(30, 170), initial_value=0,
                              text_template="{angle:5.1f}°°",font_size=12,slider_inner_radius=20,slider_outer_radius=24,handle_outer_radius=5)

scene.add(x_ring_slider)
scene.add(y_ring_slider)
scene.add(z_ring_slider)
x_ring_slider.set_visibility(True)
y_ring_slider.set_visibility(True)
z_ring_slider.set_visibility(True)



def change_slice_x(slider):
    angle = slider.value
    previous_angle = slider.previous_value
    rotation_angle = angle - previous_angle
    showm.scene.elevation(rotation_angle)
    
x_ring_slider.on_change = change_slice_x

def change_slice_y(slider):
    angle = slider.value
    previous_angle = slider.previous_value
    rotation_angle = angle - previous_angle
    showm.scene.azimuth(rotation_angle)
    
y_ring_slider.on_change = change_slice_y



def change_slice_z(slider):
    angle = slider.value
    previous_angle = slider.previous_value
    rotation_angle = angle - previous_angle
    showm.scene.roll(rotation_angle)
    
z_ring_slider.on_change = change_slice_z

xlabel = ui.TextBlock2D(position=(70,40),text="x-axis")
ylabel = ui.TextBlock2D(position=(70,100),text="y-axis")
zlabel = ui.TextBlock2D(position=(70,160),text="z-axis")
scene.add(xlabel)
scene.add(ylabel)
scene.add(zlabel)


ax = actor.axes(scale=(100, 100, 100))
# scene.add(ax)

# mesh_structure = mcds.get_mesh()
# X_domain = np.unique(mesh_structure[0][0])
# Y_domain = np.unique(mesh_structure[1][0])
# Z_domain = np.unique(mesh_structure[2][0])
# dx = X_domain[1]-X_domain[0]
# dy = Y_domain[1]-Y_domain[0]
# dz = Z_domain[1]-Z_domain[0]
# xlims = np.array([X_domain[0]-dx/2, X_domain[-1]+dx/2])
# ylims = np.array([Y_domain[0]-dy/2, Y_domain[-1]+dy/2])
# zlims = np.array([Z_domain[0]-dz/2, Z_domain[-1]+dz/2])

#Drawing Domain Boundaries
lines = [np.array([[-2880.,-500.,-2880.],[2880.,-500.,-2880.],[2880.,500.,-2880.],[-2880.,500.,-2880.],[-2880.,-500.,-2880.],[-2880.,-500.,2880.],[-2880.,500.,2880.],[-2880.,500.,-2880.],[-2880.,500.,2880.],[2880.,500.,2880.],[2880.,500.,-2880.],[2880.,500.,-2880.],[2880.,-500.,-2880.],[2880.,-500.,2880.],[2880.,500.,2880.],[2880.,-500.,2880.],[-2880.,-500.,2880.]])]
colors = np.array([0.5, 0.5, 0.5])
c = actor.line(lines, colors)
scene.add(c)

#Adding Dimension Labels
x_label = actor.text_3d(text='x axis (micron)',position=(-750.0,-700.0,3000.0),font_size=200,justification='left',color=(0,0,0))
scene.add(x_label)
y_label = actor.text_3d(text='y axis (micron)',position=(3000,0,3000.0),font_size=200,justification='left',color=(0,0,0))
scene.add(y_label)
z_label = actor.text_3d(text='z axis (micron)',position=(3000,-700.0,0.0),font_size=200,justification='left',color=(0,0,0))
scene.add(z_label)




#scene.add(actor.texture_on_sphere(image))

showm.initialize()
showm.scene.reset_camera()
scene.background((1, 1, 1))
scene.set_camera(position=(5026.62, 2766.0, 9293.52), focal_point=(221.95, -75.04, -77.73), view_up=(-0.06, 0.963, -0.26))

#showm.render()

#showm.start()
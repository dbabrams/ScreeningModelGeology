# -*- coding: utf-8 -*-
"""
This script intends to create a "toy model" that figures out how to implement 
the lake package in FloPy.

Created on Tue Feb 11 10:58:14 2020

@author: alljones
"""

#%% IMPORTS
import flopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import os

# import Shelby Specials
from modules.postProcessing import *

#%% INITIALIZE IMPORTANT DIRECTORIES AND THE MODEL OBJECT

# define the model name
model_name = "lake_toy_model"

# directory of MODFLOW executable
exe_dir = './modflowdir/mf2005'

# directory for saving model files
from datetime import datetime
folder_date   = datetime.now().strftime('%Y%m%d')
modelfile_dir = ('../TOY_MODEL-Model Files/v{}/'.format(folder_date)) 

# create the model object
m = flopy.modflow.Modflow( model_name, exe_name=exe_dir, model_ws=modelfile_dir )


#%% INITIALIZE THE GRID AND DISCRETIZE MODEL - SIMPLE

#assign discretization variables 
Lx = 100. #x dimension (m)
Ly = 100. #y dimension (m)
ztop = 0. #top of model (m)
zbot = -10. #bottom of model (m)
nlay = 1 #number of model layers
nrow = 25 #number model rows (cuts up model along Ly dimension)
ncol = 25 #number of model columns (cuts up model along Lx dimension)
dx = Lx/ncol #grid spacing along Lx dimension
dy = Ly/nrow #grid spacing along Ly dimension
dz = (ztop - zbot) / nlay #grid spacing between layers

#specify number of stress periods
nper = 1

#specify if stress period is transient (False) or steady-state (True)
steady = [True]

#create flopy discretization object, length and time are meters (2) and days (4)
dis = flopy.modflow.ModflowDis(model=m, 
                               nlay=nlay, 
                               nrow=nrow, 
                               ncol=ncol, 
                               delr=dx, 
                               delc=dy, 
                               top=ztop, 
                               botm=zbot, 
                               itmuni = 4, 
                               lenuni = 2, 
                               nper=nper, 
                               steady=steady)


#%% ESTABLISH RIVER CELLS - WESTERN BOUNDARY
#create list to store river stress period data
rivers = []
stage_riv1 = 5 #set west river stage to 5 m
cond = 1 #set sediment conductance as 1 m^2/d
rbot = 0 #set location of river bottom (set to top of model, 0 m)

#add rivers 
for i in range(1,nrow-1):
   rivers.append([0,i,0,stage_riv1,cond,rbot]) #set left side as column of river cells

riv_spd = {0: rivers} #appropriate dictionary setup for flopy stress period data
riv = flopy.modflow.ModflowRiv(m, stress_period_data=riv_spd, ipakcb=1) #assign river package


#%% ADD A WELL IN THE SW QUADRANT

#add well
#Create Single Well at center of domain with [lay, row, col, flux] list
pumping_rate = -50 #in m^3/d, negative for pumping/positive for injection
well_1 = [0, 3*nrow/4, ncol/4, pumping_rate]
# well_2 = [0,nrow/1.3,ncol/4,0]
wel_spd = {0: [well_1]}
#Create flopy wel object 
wel = flopy.modflow.ModflowWel(model=m, stress_period_data=wel_spd)


#%% ESTABLISH INACTIVE CELLS

#create ibound as array (1: active, 0: inactive, -1: constant head)
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
#assign top and bottom as inactive
ibound[:,0,:] = 0
ibound[:,-1,:] = 0

#!!! MAKE LAKE CELLS INACTIVE!!! - method?
lakerc = np.zeros( (nrow,ncol), dtype=bool )
lakerc[int(np.floor(nrow/3)-2):int(np.floor(nrow/3)+2), 
       int(np.floor(ncol/2)-2):int(np.floor(ncol/2)+2)] = True
ibound[0, lakerc] = 0

#designate starting heads as array of floats = 1.0 m
strt = np.ones((nlay, nrow, ncol), dtype=np.float32)

#create flopy bas object 
bas = flopy.modflow.ModflowBas(m, ibound=ibound, strt=strt)


#%% CHECK THE MODEL BOUNDARY CONDITIONS - PRE-LAKE

# # inifitalize figure
# fig = plt.figure( figsize=(8,6) )
# #use flopy to check location of ibound, and rivers
# modelmap  = flopy.plot.PlotMapView(model=m, layer=0)
# grid      = modelmap.plot_grid()
# ib        = modelmap.plot_ibound()
# riv_cell  = modelmap.plot_bc(ftype='RIV',color='cyan')
# well_cell = modelmap.plot_bc(ftype='WEL',color='red')
# #add labels and legend
# plt.xlabel('Lx (m)',fontsize = 14)
# plt.ylabel('Ly (m)',fontsize = 14)
# plt.title('Boundary Conditions (pre Lake)', fontsize = 15, fontweight = 'bold')
# plt.legend(handles=[mp.patches.Patch(color='white',label='Active Cell',ec='black'),
#                     mp.patches.Patch(color='black',label='Inactive Cell',ec='black'),
#                     mp.patches.Patch(color='cyan',label='Rivers',ec='black'),
#                     mp.patches.Patch(color='red', label='Well',ec='black')
#                     ], bbox_to_anchor=(1.01,0.5), loc='center left')
# plt.tight_layout()
# plt.show()


#%% ADDING THE LAKE BC

# array of lake values -> !!! opportunity for method/module
# array (nlay, nrow, ncol) of lake locations
# numbers in cell correspond to lake number
#
# NOTE:
# Lake cells must be inactive cells (IBOUND = 0) and should not be
# convertible to active cells (WETDRY = 0).
lakeid = np.zeros( (nlay, nrow, ncol) )
lakeid[:, lakerc] = 1

# bed leakances
bed_leaks = np.ones( (nlay, nrow, ncol) )*0.01

# sublake lists -> not including for now... (AEJ 2020-02-11)
# keyed - to stress period
# list of tuples - ( # sub lakes, [id numbers in lakeid for sublakes] )
sublake = { 0: [(0, [])] }

# stress period data - flux_data
# keyed to stress period
# list of lists - [[ PRCPLK EVAPLK RNF WTHDRW [SSMN]** [SSMX]** ], ...]
#
# ** [SSMN] [SSMX] variables are only necessary if you want to update/modify 
# stage_range variable for later steady-state stress periods 
lak_spd= { 0: [[1.0, 2.4, -1.0, -100.0],]}

# NOTES ON RUNOFF:
# RNF : float
#             Overland runoff from an adjacent watershed entering the lake.
#             If RNF > 0, it is specified directly as a volumetric rate, or
#             flux (L3 /T). If RNF < 0, its absolute value is used as a
#             dimensionless multiplier applied to the product of the lake
#             precipitation rate per unit area (PRCPLK) and the surface area
#             of the lake at its full stage (occupying all layer 1 lake
#             cells). When RNF is entered as a dimensionless multiplier
#             (RNF < 0), it is considered to be the product of two
#             proportionality factors. The first is the ratio of the area of
#             the basin contributing runoff to the surface area of the lake
#             when it is at full stage. The second is the fraction of the
#             current rainfall rate that becomes runoff to the lake. This
#             procedure provides a means for the automated computation of
#             runoff rate from a watershed to a lake as a function of
#             varying rainfall rate. For example, if the basin area is 10
#             times greater than the surface area of the lake, and 20 percent
#             of the precipitation on the basin becomes overland runoff
#             directly into the lake, then set RNF = -2.0.

# NOTES ON WITHDRAWAL:
# Negative values are supposed to reflect anthropogenic discharge 
# ("augmentation") of the lake with a specific volume flow rate.

# initiating lake package
lak = flopy.modflow.ModflowLak(model=m,
                               nlakes=1,  # number of lakes
                               ipackcb=1, # write cell-by-cell flows
                               lwrt=[0], # printout from lake package
                               theta=1.0,
                               nssitr=25., # newton solver iterations
                               sscncr=1.0, # convergence criteria
                               surfdepth=0.01, # height of undulations in lake-bottom
                               stages=5.0, # initial stage of lake
                               stage_range=[ (0.5, 20.0) ], # list of stage range for lakes -> steady-state solutions ONLY; remove in transient
                               lakarr=lakeid, # lake locations populated with lake number
                               bdlknc=bed_leaks, # bed leakances array
                               flux_data=lak_spd # dictionary w/lists of lake data -> stress_period_data
                               )
# OTHER INPUTS:
# sill_data=sublake, # dictionary w/lists of sublakes

# creating a stress period data file
# FOR PLOTTING PURPOSES ONLY!!
lak_loc = [(0, int(row), int(col)) for row, col in zip( np.where(lakerc)[0], np.where(lakerc)[1] )]
lak_loc = np.array( lak_loc, dtype=[('k', '<i4'), ('i', '<i4'), ('j', '<i4')] )
lak.stress_period_data = {0: lak_loc }


#%% DEFINE LAY PEROPERTIES - LPF PACKAGE
#define horizontal and vertical hydraulic conductivity
hk = np.ones((nlay,nrow,ncol), dtype=np.float32) #horizontal conductivity per cell = 1m/d
vk = np.ones((nlay,nrow,ncol), dtype=np.float32) #vertical conductivity per cell = 1m/d

#define specific storage
ss = np.ones((nlay,nrow,ncol), dtype=np.float)
ss[:,:,:] = 1e-5 #per cell

#define layer type as confined (flag: confined = 0, unconfined > 0)
laytyp = np.zeros((nlay,), dtype=np.int32)

#create flopy layer property flow object
lpf = flopy.modflow.ModflowLpf(model=m, hk=hk, vka=vk, ss=ss, laytyp=laytyp, ipakcb=1)


#%% INITIATE RECHARGE RATE
rch_val = .005 #rch value per cell in m^3/d
rch_array = np.zeros((nrow,ncol)) #rch array to hold values per cell

#apply recharge to all cells except inactive and river cell boundary conditions
rch_array[1:-1,1:-1] = rch_val
rch_array[ lakerc ]  = 0.0

#create recharge object
rch = flopy.modflow.ModflowRch(model=m,rech=rch_val)


#%% CREATE OUTPUT CONTROL (OC) PACKAGE
#create oc stress period data. 
oc_spd = {(0, 0): ['print head', 'print budget', 'save head', 'save budget']}

#create flopy output control object
oc = flopy.modflow.ModflowOc(model=m, stress_period_data=oc_spd, compact=True)


#%% ASIGN PCG AS SOLVER
#assign groundwater flow solver
pcg = flopy.modflow.ModflowPcg(model=m)


#%% WRITE THE MODEL FILES
#write MODFLOW input files
m.write_input()


#%% CHECK THE MODEL BOUNDARY CONDITIONS - POST-LAKE

# inifitalize figure
fig = plt.figure( figsize=(8,6) )
#use flopy to check location of ibound, and rivers
modelmap  = flopy.plot.PlotMapView(model=m, layer=0)
grid      = modelmap.plot_grid()
ib        = modelmap.plot_ibound()
lak_cell  = modelmap.plot_bc(ftype='LAK',color='navy')
riv_cell  = modelmap.plot_bc(ftype='RIV',color='cyan')
well_cell = modelmap.plot_bc(ftype='WEL',color='red')
#add labels and legend
plt.xlabel('Lx (m)',fontsize = 14)
plt.ylabel('Ly (m)',fontsize = 14)
plt.title('Boundary Conditions (POST-Lake)', fontsize = 15, fontweight = 'bold')
plt.legend(handles=[mp.patches.Patch(color='white',label='Active Cell',ec='black'),
                    mp.patches.Patch(color='black',label='Inactive Cell',ec='black'),
                    mp.patches.Patch(color='navy',label='Lakes',ec='black'),
                    mp.patches.Patch(color='cyan',label='Rivers',ec='black'),
                    mp.patches.Patch(color='red', label='Well',ec='black')
                    ], bbox_to_anchor=(1.01,0.5), loc='center left')
plt.tight_layout()
plt.show()




#%% RUN THE MODEL
# Run the model
success, mfoutput = m.run_model(pause=False, report=True)
if not success:
    raise Exception('MODFLOW did not terminate normally.')
else:
    print('\n'+'*****************************\n')


#%% GET DATA FOR VISUALIZATION!

# identify the first time period
totim=1
    
# obtain the head and flow budget objects
headobj = flopy.utils.binaryfile.HeadFile(m._model_ws+model_name+'.hds')
budgobj = flopy.utils.binaryfile.CellBudgetFile(m._model_ws+model_name+'.cbc')

# get the data for the specific stress period (?)
head = headobj.get_data(totim=totim)
frf = budgobj.get_data(text='flow right face', totim=totim)
fff = budgobj.get_data(text='flow front face', totim=totim)

#%% PLOT 3D SURFACE

# from Shelby's postProcessing
figure = headSurfacePlot(m, Lx, Ly, head, title='Head Surface')


#%% STILL NEED...
# - Layer properties - LPF package
# - Recharge - Rch package
# - Output control - OC package
# - assign PCG as solver
# - run the model
# - show outputs...? -> postProcessing



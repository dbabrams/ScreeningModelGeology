# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 10:13:29 2018

@author: shelbya2
"""
#import packages
from mpl_toolkits.mplot3d import Axes3D
import flopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import os

'''
see accompanying 'Notes on Post-Processing' PDF
in 'Documentation' folder
'''
def getMODFLOWResults(modelname,totim,get_ch=False, ws=''):
    try:
        headobj = flopy.utils.binaryfile.HeadFile(modelname+'.hds')
        budgobj = flopy.utils.binaryfile.CellBudgetFile(modelname+'.cbc')

        head = headobj.get_data(totim=totim)
        frf = budgobj.get_data(text='flow right face', totim=totim)
        fff = budgobj.get_data(text='flow front face', totim=totim)
        
        if get_ch is not False:
            ch = budgobj.get_data(text='constant head', totim=totim)
            return(head, frf, fff, ch)
        else:
            return(head,frf,fff)
    
    except FileNotFoundError as e:
        print(e, '\n \n Change working directory to model datafile location \n (Current directory:',
                os.getcwd(), ') \n'
              '\n If file does not exist:',
              '\n Check OC Package save-data designation & other MODFLOW Pkg save-budget flags (ipackb > 1) ')


def headContourPlot(model, head, frf, fff, 
                    title, subtitle, grid=True, cmap='viridis', 
                    plot_disch=True, save_plot=False):
    #data, color-plotting setup
    masked = np.ma.masked_where(model.bas6.ibound.array == 0 , head) #mask inactive cells cells
    ticker = mp.ticker.MaxNLocator(nbins=model.dis.ncol)
    levels = ticker.tick_values(np.min(masked), np.max(masked)) #set contour levels
    cmap = plt.get_cmap(cmap) #get colormap
    norm = mp.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True) #normalize colormap to dataset levels
    
    #create figure and plots
    fig, ax = plt.subplots(figsize=(8,7)) #create figure
    modelmap = flopy.plot.PlotMapView(model=model, layer=0, ax=ax) #use modelmap to attach plot to model
    pcolor = modelmap.plot_array(masked, norm = norm, animated=True) #create contours
    if grid is True:
        grid = modelmap.plot_grid() #plot grid if True
    if plot_disch is True:
        modelmap.plot_discharge(frf[0], fff[0]) #create discharge arrows if true
    
    #figure display parameters
    plt.xlabel('Lx (m)',fontsize = 13)
    plt.ylabel('Ly (m)',fontsize = 13)
    plt.suptitle(title, fontsize = 20, fontweight = 'bold',x=.44,y=.98) #bold title
    plt.title(subtitle, fontsize = 14,y=1.015) #subtitle
    cb = plt.colorbar(pcolor)
    cb.set_label('Head (m)', fontsize = 13, labelpad = 12) #add colorbar for head data
    plt.show()
    return(fig)

##NOTE! 3d plots utilize a meshgrid which puts the head on the corners so isn't exatly realistic
    #but gives good idea for head topography (especially large models)

def headSurfacePlot(model, Lx, Ly, head, 
                    title = '', cmap='viridis'):
    #create 3d figure
    fig_3d = plt.figure(figsize=(12,5))
    ax = fig_3d.gca(projection='3d')
    
    #set up data space
#    masked = np.where(model.bas6.ibound.array !=0, head, np.nan)
    masked = head[:, 1:-1, :]
    x = np.linspace(0,Lx,model.dis.ncol) #set x domain using model dis pkg
    y = np.linspace(0,Ly,model.dis.nrow) #set y domain using model dis pkg
#    x, y = np.meshgrid(x, y) #mesh domain
    x, y = np.meshgrid(x, y[1:-1]) #mesh domain
    
    #set color levels
    levels = np.linspace(np.nanmin(masked),
                         np.nanmax(masked),
                         model.dis.ncol-1)
    cmap = plt.get_cmap(cmap) #get colormap
    norm = mp.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True) #normalize colormap to dataset values
    surf = ax.plot_surface(x,y,np.flipud(masked[0]), 
                           cmap = cmap, linewidth=0, antialiased=False, 
                           label='head', norm=norm) #plot active head surface
    surf.set_facecolor((0,0,0,0))
    cb = fig_3d.colorbar(surf,shrink=0.5,aspect=5)
    cb.set_label('Head (m)',fontsize=15,labelpad = 10) #add colorbar 
    ax.set_xlabel('Lx (m)', fontsize=13, labelpad = 10)
    ax.set_ylabel('Ly (m)', fontsize=13, labelpad = 10)
    ax.set_title(title, fontsize=15, y=1.05)
    #suppress matplotlib warning (3d plot is angry at nan values but still works)
    import warnings
    warnings.filterwarnings("ignore") 
    plt.show()
    return(fig_3d)



def chError(model, head_flow, head_ch, trial, grid=True, SP=0):
    #calculate difference and mask inactive grid cells
    diff = head_flow-head_ch
    masked = np.ma.masked_where(model.bas6.ibound.array == 0 , diff)
    
    #plot difference
    fig, ax = plt.subplots(figsize=(8,7))
    modelmap = flopy.plot.PlotMapView(model=model, layer=0, ax=ax) #use modelmap to attach plot to model
    dif_plot = modelmap.plot_array(masked)
    if grid is True:
        modelmap.plot_grid()
    
    #display parameters
    plt.xlabel('Lx (m)',fontsize = 13)
    plt.ylabel('Ly (m)',fontsize = 13)
    plt.suptitle(' Head Differences (m) SP: %s'%(SP), 
                 fontsize = 20, fontweight = 'bold',x=.44,y=.98)
    plt.title('Flow Spec - Head Spec (Iter %s)'%(trial-1), 
              fontsize = 14,y=1.02)
    cb = plt.colorbar(dif_plot)
    cb.set_label('Head (m)', fontsize = 13, labelpad = 12) 
    plt.show()
    return(diff,fig)
    
    
def chPercentDiff(model, head_flow, head_ch, trial, grid=True, SP=0):
    #calculate percent difference and mask inactive grid cells
    perc_diff = np.empty((model.dis.nlay,model.dis.nrow,model.dis.ncol))
    for i in range(model.dis.nrow):
        for j in range(model.dis.ncol):
            perc_diff[0,i,j] = abs(head_flow[0,i,j]-head_ch[0,i,j])/head_flow[0,i,j]*100
    masked = np.ma.masked_where(model.bas6.ibound.array == 0 , perc_diff)
    
    #plot difference
    fig, ax = plt.subplots(figsize=(8,7))
    modelmap = flopy.plot.PlotMapView(model=model, layer=0, ax=ax) #use modelmap to attach plot to model
    dif_plot = modelmap.plot_array(masked)
    if grid is True:
        modelmap.plot_grid()
    
    #display parameters
    plt.xlabel('Lx (m)',fontsize = 13)
    plt.ylabel('Ly (m)',fontsize = 13)
    plt.suptitle('Head Percent Diff, SP: %s'%(SP), 
                 fontsize = 20, fontweight = 'bold',x=.44,y=.98)
    plt.title('Head vs Flow Spec. (Iter %s)'%(trial-1), 
              fontsize = 14,y=1.02)
    cb = plt.colorbar(dif_plot)
    cb.set_label('Head (m)', fontsize = 13, labelpad = 12) 
    plt.show()
    return(perc_diff,fig)
    
        #plot results (mason county only)
def plotMasonCounty(model, cell_size, head, 
                    title = 'Head Spec. Mason County', subtitle = ''):

    row_cut = 250
    col_cut = 175
    x = np.linspace(0,col_cut*cell_size,col_cut)
    y = np.linspace(0,row_cut*cell_size,row_cut)
    x, y = np.meshgrid(x, y)
    masked_head = np.ma.masked_where(model.bas6.ibound.array == 0 , head)
    fig, ax = plt.subplots(figsize = (10,9))
    flipped = np.flipud(masked_head[0])
    f = ax.contourf(x,y,flipped[0:row_cut,0:col_cut],30)
    #display parameters
    plt.xlabel('Lx (m)',fontsize = 14, labelpad = 10)
    plt.ylabel('Ly (m)',fontsize = 14, labelpad = 10)
    plt.suptitle(title, fontsize = 20, fontweight = 'bold',x=.44,y=.96)
    plt.title(subtitle, fontsize = 14,y=1.015)
    plt.colorbar(f).set_label('Head (ft)', fontsize = 14, labelpad = 15)
    plt.show()
    return(fig)
    
def masonAnimation(model, cell_size, head, 
                    title = 'Head Spec. Mason County', subtitle = ''):

    row_cut = 250
    col_cut = 175
    x = np.linspace(0,col_cut*cell_size,col_cut)
    y = np.linspace(0,row_cut*cell_size,row_cut)
    x, y = np.meshgrid(x, y)
    masked_head = np.ma.masked_where(model.bas6.ibound.array == 0 , head)
    fig, ax = plt.subplots(figsize = (10,9))
    flipped = np.flipud(masked_head[0])
    f = ax.contourf(x,y,flipped[0:row_cut,0:col_cut],30)
    #display parameters
    plt.xlabel('Lx (m)',fontsize = 14, labelpad = 10)
    plt.ylabel('Ly (m)',fontsize = 14, labelpad = 10)
    plt.suptitle(title, fontsize = 20, fontweight = 'bold',x=.44,y=.96)
    plt.title(subtitle, fontsize = 14,y=1.015)
    plt.colorbar(f).set_label('Head (ft)', fontsize = 14, labelpad = 15)
    plt.show()
    return(fig)
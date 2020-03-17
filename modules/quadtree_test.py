# -*- coding: utf-8 -*-
"""
Created on Tues Jan 21 2020

@author: krasows2
"""
# "Supreme executive power derives from a mandate from the masses, not from some farcical aquatic ceremony."

# copied after https://github.com/modflowpy/flopy/blob/develop/examples/Notebooks/flopy3_gridgen.ipynb

#%% IMPORTS
import flopy
from flopy.utils.gridgen import Gridgen
import geopandas as geopd
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import time

# Set path to location of gridgen.exe
os.environ["PATH"] = "./modflowdir/"

#%% READ IN FILES
mod_dom_fn = 'ESL_model_domain'
mod_dom_fp = 'A:/CGS/ESLgisData/Shapefiles/ESL_model_domain.shp'
mod_dom = geopd.read_file(mod_dom_fp)

act_dom_fn = 'ESL_Active_Domain'
act_dom_fp = 'A:/CGS/ESLgisData/Shapefiles/ESL_Active_Domain.shp'
act_dom = geopd.read_file(act_dom_fp)

ambot_fn = 'AmericanBottoms'
ambot_fp = 'A:/CGS/ESLgisData/Shapefiles/AmericanBottoms.shp'
ambot = geopd.read_file(ambot_fp)

refline_rb_fn = 'ESL_RefLines_riverbank'
refline_rb_fp = 'A:/CGS/ESLgisData/Shapefiles/ESL_RefLines_riverbank.shp'
refline_rb = geopd.read_file(refline_rb_fp)

refline_bl_fn = 'ESL_RefLines_bluff'
refline_bl_fp = 'A:/CGS/ESLgisData/Shapefiles/ESL_RefLines_bluff.shp'
refline_bl = geopd.read_file(refline_bl_fp)

miss_riv_fn = 'ESL_MississippiRv'
miss_riv_fp = 'A:/CGS/ESLgisData/Shapefiles/ESL_MississippiRv.shp'
miss_riv = geopd.read_file(miss_riv_fp)

idot_wells_fn = 'IDOT_dewat_wells'
idot_wells_fp = 'A:/CGS/ESLgisData/Shapefiles/IDOT_dewat_wells.shp'
idot_wells = geopd.read_file(idot_wells_fp)


#%% Setup Base MODFLOW Grid
Lx = mod_dom.bounds.loc[0, 'maxx'] - mod_dom.bounds.loc[0, 'minx']
Ly = mod_dom.bounds.loc[0, 'maxy'] - mod_dom.bounds.loc[0, 'miny']
nlay = 3
delr = 2640
delc = 2640
nrow = math.ceil(Ly/delc)
ncol = math.ceil(Lx/delr)
h0 = 10
h1 = 5
top = h0
botm = np.zeros((nlay, nrow, ncol), dtype=np.float32)
botm[1, :, :] = -10.

#defining upper left reference for model discretization
xul = mod_dom.bounds.loc[0, 'minx']
yul = mod_dom.bounds.loc[0, 'maxy']


#%% INITIALIZE MODFLOW AND DISCRETIZATION OBJECTS
ms = flopy.modflow.Modflow()  # rotation! ms = flopy.modflow.Modflow(rotation=-20.)
dis = flopy.modflow.ModflowDis(ms, nlay=nlay, nrow=nrow, ncol=ncol, delr=delr,
                               delc=delc, top=top, botm=botm, xul=xul, yul = yul)


#%% CREATE GRIDGEN OBJECT
model_ws ='A:/CGS/ESLgisData/Shapefiles' # os.path.join('.', 'data')
g = Gridgen(dis, model_ws=model_ws)


#%% ADD ACTIVE DOMAIN (OPTIONAL)
g.add_active_domain(act_dom_fn, range(nlay))


#%% SPECIFY REFINEMENTS

# for point/IDOT wells
g.add_refinement_features(idot_wells_fn, 'point', 5, range(nlay))

# for polyline/Mississippi River bank
g.add_refinement_features(refline_rb_fn, 'line', 4, range(nlay))

# for polyline/American Bottoms eastern bluffs
g.add_refinement_features(refline_bl_fn, 'line', 2, range(nlay))

# for polygon/Mississippi River
g.add_refinement_features(miss_riv_fn, 'polygon', 1, range(nlay))


#%% PLOT GRIDGEN INPUTS (PRE-REFINEMENT)
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
mm = flopy.plot.PlotMapView(model=ms)
mm.plot_grid()
flopy.plot.plot_shapefile( miss_riv_fp, ax=ax, facecolor='none', edgecolor='blue', alpha=1 )
flopy.plot.plot_shapefile( mod_dom_fp, ax=ax, facecolor='none', edgecolor='yellow', alpha=1 )
flopy.plot.plot_shapefile( idot_wells_fp, ax=ax, facecolor='red', edgecolor='red', alpha=1 )
flopy.plot.plot_shapefile( ambot_fp, ax=ax, facecolor='none', edgecolor='green', alpha=1 )


#%% Build the Grid
# note: more refinement means longer run time
start = time.time()
print('Started "Build the Grid" cell.')

g.build(verbose=False)

end = time.time()
print('"Build the Grid" cell finished after ' + 
      '{0:0.2f} minutes'.format((end-start)/60))


#%% PLOT GRIDGEN OUTPUTS (POST-REFINEMENT)
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
g.plot(ax, linewidth=0.5)
flopy.plot.plot_shapefile( miss_riv_fp, ax=ax, facecolor='none', edgecolor='blue', alpha=1 )
flopy.plot.plot_shapefile( mod_dom_fp, ax=ax, facecolor='none', edgecolor='yellow', alpha=1 )
flopy.plot.plot_shapefile( idot_wells_fp, ax=ax, facecolor='red', edgecolor='red', alpha=1 )
flopy.plot.plot_shapefile( ambot_fp, ax=ax, facecolor='none', edgecolor='green', alpha=1 )


#%% CREATE A FLOPY ModflowDisu OBJECT
mu = flopy.modflow.Modflow(model_ws=model_ws, modelname='mfusg')
disu = g.get_disu(mu)
disu.write_file()










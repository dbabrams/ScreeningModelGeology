# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 14:16:02 2019

This file will take the previous iteration of the East St. Louis python 
notebook and translate it into an improved script and module component.

@author: Illinois State Water Survey
"""

# "You're all individuals." "We're all individuals."

#%% IMPORTS

# import model packages
import flopy
import geopandas as geopd
import matplotlib.pyplot as plt
import numpy as np
import time

# interpolation
from pykrige.uk import UniversalKriging

# import American Bottoms function module
from modules.amerBot_module import *

# setting start time (just for kicks)
start_time = time.time()


#%% INITIALIZE IMPORTANT DIRECTORIES AND MODEL OBJECTS

# define the model name
model_name = "esl_steady"

# directory of MODFLOW executable
exe_dir = './modflowdir/MODFLOW-NWT.exe'

# directory for saving model files
from datetime import datetime
folder_date   = datetime.now().strftime('%Y%m%d')
modelfile_dir = ('../MODFLOW_FloPy-Model Files/{}_ESL_RGWFM/'
                 .format(folder_date)) 

### ------------------------------ model object ------------------------------

mf_obj = flopy.modflow.Modflow( model_name,
                                version="mfnwt",
                                exe_name=exe_dir,
                                model_ws=modelfile_dir )

### ----------------------------- specify solver -----------------------------

nwt = flopy.modflow.ModflowNwt( mf_obj )

### ----------------------------- output control -----------------------------

# specify number of stress periods and if it is steady state
nper = 1
perlen = [1]
nstp = [1]
steady = [True]
# Output control
stress_period_data = {}
for kper in range(nper):
    for kstp in range(nstp[kper]):
        stress_period_data[(kper, kstp)] = ['save head',
                                            'save drawdown',
                                            'save budget',
                                            'print head',
                                            'print budget']

oc = flopy.modflow.ModflowOc( mf_obj,
                              stress_period_data=stress_period_data,
                              compact=True)


#%% READ IN SHAPEFILES

# the model and active domains, and refinement zone are all incorporated into 
# the __init() function of the SimpleGridRefine class.
# current shapefiles are:
# - 'ESL_model_domain.shp'
# - 'ESL_Active_Domain.shp'
# - 'ESL_RefZone1.shp'
# in the following filepath:
# - 'A:/CGS/ESLgisData/Shapefiles/'

# largely used for plotting
ambot_fn = 'AmericanBottoms'
ambot_fp = 'A:/CGS/ESLgisData/Shapefiles/AmericanBottoms.shp'
ambot = geopd.read_file(ambot_fp)

# largely used for plotting at the moment, future, boundary condition
miss_riv_fn = 'ESL_MississippiRv'
miss_riv_fp = 'A:/CGS/ESLgisData/Shapefiles/ESL_MississippiRv.shp'
miss_riv = geopd.read_file(miss_riv_fp)

# largely used for plotting at the moment, future, boundary condition
idot_wells_fn = 'IDOT_dewat_wells'
idot_wells_fp = 'A:/CGS/ESLgisData/Shapefiles/IDOT_dewat_wells.shp'
idot_wells = geopd.read_file(idot_wells_fp)


#%% GRID REFINEMENT BASED ON SHAPEFILES

# assigns refined grid characteristics to object (see simpleFloPyGridRefine.py
#   for more information and to specify shapefiles to be used)
from modules.simpleFloPyGridRefine import SimpleGridRefine 
refined = SimpleGridRefine()


#%% READ IN ELEVATION RASTERS AND ASSIGN TO MODELGRID

# pass the file path to the folder with your elevation rasters (.tif) and the 
#   SimpleGridRefine() output to assign_prop_by_raster module
elevs_fp = 'A:/CGS/ESLgisData/Elevation/ModelSurfaces/'
elevs_array = assign_prop_by_raster(elevs_fp, refined)
nlay = elevs_array.shape[0] -1 


#%% PUT LAYERS IN ORDER (BY MEAN ELEVATION)
# (works well for 9 layer model, not sure about others. For others, surely
#   want to make sure cell by cell that layers don't intersect.)

mean_elev = np.zeros( nlay + 1 )
for iii in range( len( mean_elev ) ):
    mean_elev[iii] = elevs_array[iii][elevs_array[iii] > 0 ].mean()

# set no_data_value pixels to np.nan
elevs_array[ np.where( elevs_array < -10**10 ) ] = np.nan

# make sure layers are in order by comparing means
check_bool = False
while check_bool == False:
    check_bool = True
    for iii in range( len( mean_elev ) ):
        if iii < len(mean_elev)-1:
            if mean_elev[iii] >= mean_elev[iii+1]:
                print('Surface {0} is above surface {1}'.format( iii, iii+1 ))
            else:
                print('Surface {0} is BELOW surface {1}'.format( iii, iii+1 ))
                check_bool = False
                problem = iii+1
                
                place_holder = mean_elev[iii]*1
                mean_elev[iii] = mean_elev[problem]*1
                mean_elev[problem] = place_holder*1
                
                ph2 = elevs_array[iii]*1
                elevs_array[iii] = elevs_array[problem]*1
                elevs_array[problem] = ph2*1
                
                print('Surfaces reordered')
                
print('Surfaces ordered by mean elevation.')


#%% READING IN "PERCENT COARSE" FROM THE 9-LAYER MODEL

# pass the file path to the folder with your rasters (.tif) and the 
#   SimpleGridRefine() output
hk_fp = 'A:/CGS/ESLgisData/HydCond/'
pct_array = assign_prop_by_raster(hk_fp, refined)


#%% IF COVERAGE OF RASTERS ISN'T PERFECT

# This is an ESL specific solution. The coverage isn't perfect along the west
#   edge (MS River) and therefore data sets are extended by filling in with
#   the western-most value extended westward.

# Filling in elevations
for iii in range( elevs_array.shape[0] ):
    print('Surface elevation extended:', iii)

    #filling in the rows...
    for jjj in range( refined.nrow_ref ):
        elevs_array[iii][jjj, np.isnan(elevs_array[iii][jjj,:])] = \
            elevs_array[iii][jjj, ~np.isnan(elevs_array[iii][jjj,:])][0]


# Filling in pct values
pct_array[np.where(pct_array < -10**10 )] = np.nan
for iii in range ( pct_array.shape[0] ):
    print('Layer hk extended:', iii)

    #filling in the rows...
    for jjj in range( refined.nrow_ref ):
        pct_array[iii][jjj, np.isnan(pct_array[iii][jjj,:])] = \
            pct_array[iii][jjj, ~np.isnan(pct_array[iii][jjj,:])][0]
        

#%% ASSIGNING K VALUES FROM PERCENT COARSE - BINNED

# create array for hk values to be read into
hk_array = np.zeros( pct_array.shape )
hk_array = pct_array*1 # fight the mirror

# set the number of bins you'd like
num_of_bins = 8
# set the k values for the bins
k_by_bin = np.asarray( ( 10**(-3), 10**(-1), 1, 25, 50, 100, 400, 500 ) )
# set bin width
bin_width = (np.max(pct_array) - np.min(pct_array))/num_of_bins

# breaks percent coarse into bins and assigns hyd.cond. by those bins
for bin_num in range(num_of_bins):

    if bin_num == 0:
        hk_array[ np.where( pct_array <= np.nanmin( pct_array ) + bin_width ) ] \
            = k_by_bin[ bin_num ]

    elif bin_num < num_of_bins - 1:
        hk_array[ np.where( np.bitwise_and( 
        pct_array > ( np.min( pct_array ) + bin_width*bin_num ), 
        pct_array <= ( np.min( pct_array ) + bin_width*(bin_num+1) ) ) ) ] = \
                                                            k_by_bin[bin_num]

    elif bin_num == num_of_bins - 1:
        hk_array[ np.where( pct_array > ( np.min( pct_array ) +
                           bin_width*bin_num ) ) ]  = k_by_bin[bin_num]

    else:
        print('something tricky is afoot...')

# plotting percent coarse vs hk assigned relationship visualization
plt.scatter(np.linspace(np.min(pct_array), np.max(pct_array), num_of_bins), k_by_bin)


#%% STATS ON PCT DISTRIBUTIONS
# if you want to look at what different bin numbers does to bin groupings

pct_array[np.where(pct_array < -10**10 )] = np.nan
fig, axes = plt.subplots(9,1, tight_layout=True, figsize=(5,10))

dummy = [4, 5, 8, 10]
for iii in dummy:
    for jjj in range( pct_array.shape[0] ):
        axes[jjj].hist(pct_array[jjj].reshape(-1), bins=np.linspace(np.min(pct_array),np.max(pct_array), (iii+1)))
        axes[jjj].set_xlim([0,100])


#%% RIVER CORRIDOR ADJUSTMENTS (ELEVATIONS AND CONDUCTIVITIES)
# We have good topo-bathymerty data. The top layer elevation raster right now
#   is a mix of t-b  and 9-layer
        
# store the bottom/bedrock elevation for use below
elev_bedrock = elevs_array[9]*1

# With the top raster being edited from the 9-layer model, the elevations
#   intersect the lower layers. Where this occurs, reassign the lower layer to
#   to one foot below the layer above. (also solves points with thickness = 0.)
#   Also reassigns hk where layer elevation is belpw bottom of 9-layer.
for iii in range( elevs_array.shape[0] ):
    if iii > 0:
        ouofor_idx = np.where( elevs_array[iii-1 ] <= elevs_array[iii] + 1 )
        elevs_array[iii][ouofor_idx] = elevs_array[iii-1][ouofor_idx] - 1
        
        if np.any( elevs_array[iii-1] <= elev_bedrock):
            blw_bdrk_idx = np.where( elevs_array[iii-1] <= elev_bedrock )
            hk_array[iii-1][blw_bdrk_idx] = 10**(-4)


#%% ASSIGNING THE RATIO BETWEEN hk AND vk

# from docs: "...the ratio of horizontal to vertical hydraulic conductivity..."

# create array for vka values (ratio of hk to vk)
vk_array = np.zeros( pct_array.shape )
# assign ratios. note ratio > 1 means hk > vk (length is +1 over hk_by_bin
#   b/c of added bedrock value)
vk_by_bin = np.asarray((1, 1, 2, 10, 10, 10, 2, 2, 1))

# get list of values of hk (including added bedrock value)
hk_unique = np.unique(hk_array)

# assign vka values to array by location of hk values
for iii in hk_unique:
    vk_array[ np.where( hk_array == iii ) ] = vk_by_bin[ np.where(hk_unique == iii) ]
    
vk_actual = hk_array/vk_array

# plotting percent coarse vs vk assigned relationship visualization
plt.scatter(np.linspace(np.min(pct_array), np.max(pct_array), num_of_bins), np.unique(vk_actual)[1:])
    

#%% INITIALIZE DISCRETIZATION OBJECT

dis = flopy.modflow.ModflowDis(mf_obj, 
                               nlay=nlay, 
                               nrow=refined.nrow_ref, 
                               ncol=refined.ncol_ref, 
                               delr=refined.delr, 
                               delc=refined.delc, 
                               top=elevs_array[0], 
                               botm=elevs_array[1:10], 
                               xul=refined.cbds.minx[0], 
                               yul=refined.cbds.maxy[0])

#%% INITIALIZE LAYER PROPERTIES FILE OBJECT

# layvka = 1 means values assigned to vka are hk:vk, not literal vk values
lpf = flopy.modflow.ModflowLpf(mf_obj, hk=hk_array, layvka=1, vka=vk_array)


#%% SETTING IN/ACTIVE CELLS

# This requires using the flopy.modflow.ModflowBas method where in/active is 
#   set by ibound which is an array of int where 1 is active, -1 is fixed 
#   heads (lake example), and 0 is inactive.

# determine geometries of model cells and add grid offset back in
row_coor, col_coor,_ = dis.get_node_coordinates()
# meshgrid of coordinates
colmesh, rowmesh = np.meshgrid( (col_coor +  refined.cbds.minx[0]), 
                                (row_coor + (refined.cbds.maxy[0]-refined.Ly)) )
# initializing ibound_ad (active domain)
ibound_ad = np.zeros( colmesh.shape )

# determine the points inside polygons for active domain
refined.adgdf.loc[:,'Inside Indices'] =  -9999 # create place holder

for idx in range( len(refined.adgdf['geometry']) ):
    # for each polygon determine cells inside (bespoke function)
    in_idx = find_cells_within_polygon(refined.adgdf.loc[idx, 'geometry'],
                                        colmesh, rowmesh)
    # store the determined indices
    refined.adgdf.loc[idx, 'Inside Indices'] = [in_idx]
    
    # assign inside to ibound as active (1)
    ibound_ad[ refined.adgdf.loc[idx, 'Inside Indices'] ] = 1
            
# pass the modflow object and the ibound array
bas = flopy.modflow.ModflowBas(mf_obj, ibound=ibound_ad)


#%% SETTING BC - HORSESHOE LAKE TEST CASE

# other shapefiles to plot in Layer 1
pathofpaths = 'A:/CGS/ESLgisData/Shapefiles/'
# Horseshoe Lake 
horseshoe = pathofpaths + 'Horseshoe_lake.shp'
# Canals
canals = pathofpaths + 'ESL_canals.shp'
# south lakes
souths = pathofpaths + 'ESL_southlakes.shp'
# north lakes
norths = pathofpaths + 'ESL_northlakes.shp'

# !!! this value goes into the bc work and is currently a dummy value
# extremely approximate estimate of stage/bottom of rivers
h0 = 500

# create list of data for each lake input
# initialize an empty dictionary of lake data
lrcd = {0:[]}
lay  = 0 # top?
k_lakbott = 1 #lake bottom hydraulic conductivity in ft/d
sed_thick = 1 #thickness of riverbed sediment in ft
numlks = 1 # number of lakes in model

# calculate regional conductivity from meshgrid
cond = np.zeros( (nlay, refined.nrow_ref, refined.ncol_ref) )
spacex, spacey = np.meshgrid( refined.delr, refined.delc )
cond[:,...] = (spacex*spacey*k_lakbott)/(sed_thick)

# import shapefile data
horse  = geopd.read_file( horseshoe )
southl = geopd.read_file( souths )
northl = geopd.read_file( norths )
# # find GIS data for Horseshoe lake proper
# horse = horse.loc[horse.Name == 'Horseshoe Lake', :]

for df in [horse, southl, northl, miss_riv]:
    for iii in range( len(df) ):
        # for each polygon determine cells inside (bespoke function)
        in_idx = find_cells_within_polygon( df.loc[iii,'geometry'], 
                                            colmesh, rowmesh )
        # # create lake array
        # # 0 -> not a lake
        # # # -> lake # 
        # lake_id_arr = np.zeros( botm_ref.shape )
        # lake_id_arr[2, in_idx] = numlks
        
        # if indices were found for a given shapefile
        if len(in_idx) != 0:
            for rrr, ccc in zip( np.where(in_idx)[0], np.where(in_idx)[1]): 
                # assuming list format: 
                # [layer, row, col, water elevation, conductivity, lake bed elevation]
                lrcd[0].append( [lay, rrr, ccc, 5+h0, cond[lay,rrr,ccc], h0] )
                # horseshoe lake is 4 or 5 ft deep typically, according to state park author
# attach lake package
# lakes = flopy.modflow.ModflowLak(mf_obj, 
#                                  nlakes=numlks,
#                                  stages=[5.],
#                                  stage_range=[(3., 7.)], 
#                                  lakarr=lake_id_arr, 
#                                  bdlknc=cond, 
#                                  stress_period_data=lrcd,
#                                  ) 

# model lake as River package??
lakes = flopy.modflow.ModflowRiv(mf_obj, stress_period_data=lrcd)


#%% INITIALIZE UPSTREAM WEIGHTING

# for use with NWT only. Has to be called after some undetermined/unspecified
#   modification to the modflow model object ( mf_obj ) that happens above
if "NWT" in exe_dir:
    upw = flopy.modflow.ModflowUpw( mf_obj )
    

#%% WRITING THE INPUT FILES

mf_obj.write_input()


#%% PLOTTING, PLOTTING, PLOTTING - LAYER PROPERTIES

show_plots = True

if show_plots == True:
    # create layer titles list
    layer_titles = ['Shallow 1\n', 'Shallow 2\n', 'Shallow 3\n', 
                    'Middle 1\n','Middle 2\n','Middle 3\n',
                    'Deep 1\n', 'Deep 2\n', 'Deep 3\n']
    
    import matplotlib as mpl
    
    # find vk_actual equivalent according to the model (lpf object)
    vk_in_model = lpf.hk.array/lpf.vka.array
    
    # plotting all layers
    fig, axes = plt.subplots(1,nlay, figsize=(36,8))
    fig2, axes2 = plt.subplots(1,nlay, figsize=(36,8))
    fig3, axes3 = plt.subplots(1,nlay, figsize=(36,8))
    for lll in range( nlay ):
        
        ### ---------------------- Surface Elevations ----------------------
        
        plt_mm = flopy.plot.PlotMapView( model=mf_obj, layer=lll, ax=axes[lll] )
        plt_mm.plot_grid( linewidth=0.2 )
        
        if lll == 0:
            xkcd = np.where( ibound_ad > 0, mf_obj.modelgrid.top, np.nan )
        else:
            xkcd = np.where( ibound_ad > 0, mf_obj.modelgrid.botm[lll], np.nan )
            
        elev_plot = plt_mm.plot_array( xkcd, cmap='gist_earth')
        
        elev_plot.set_clim( np.nanpercentile(elevs_array, 1), np.nanpercentile(elevs_array,99))
        
        # if lll==0:
        #     lake_plot = plt_mm.plot_bc( package=lakes, color='cyan' )
        #     flopy.plot.plot_shapefile( horseshoe, ax=axes[lll], facecolor='none', 
        #                             edgecolor='darkcyan', alpha=1 ) 
        #     flopy.plot.plot_shapefile( canals, ax=axes[lll], facecolor='none', 
        #                             edgecolor='darkgoldenrod', alpha=1 ) 
        #     flopy.plot.plot_shapefile( souths, ax=axes[lll], facecolor='none', 
        #                             edgecolor='teal', alpha=1 ) 
        #     flopy.plot.plot_shapefile( norths, ax=axes[lll], facecolor='none', 
        #                             edgecolor='teal', alpha=1 )
        
        plt_mm.plot_inactive()
        
        flopy.plot.plot_shapefile( miss_riv_fp, ax=axes[lll], facecolor='none', 
                                  edgecolor='blue', alpha=1 )
        axes[lll].plot( idot_wells.loc[:, 'geometry'].x, idot_wells.loc[:, 'geometry'].y, 
                  color='red', marker='x', linestyle='')
        flopy.plot.plot_shapefile( ambot_fp, ax=axes[lll], facecolor='none', 
                                  edgecolor='green', alpha=1 )
        axes[lll].set_title( '{} Hydrogeologic Zone'.format( layer_titles[lll]) )
        axes[lll].set_xticklabels( [str(int(tick)) for tick in axes[lll].get_xticks()], 
                                   rotation=-45, ha='left')
        
        
        ### ------------------ Horizontal Hydraulic Conductivity ------------------
        
        plt_nn = flopy.plot.PlotMapView( model=mf_obj, layer=lll, ax=axes2[lll] )
        plt_nn.plot_grid( linewidth=0.2 )
        
        xkcd_2 = np.where( ibound_ad > 0, lpf.hk.array[lll], np.nan )
        
        # hk_plot = nn.plot_array( xkcd_2, cmap='jet')
        
        hk_plot = plt_nn.plot_array( xkcd_2, 
                                norm=mpl.colors.LogNorm(vmin=np.nanmin(lpf.hk.array), 
                                                        vmax=np.nanmax(lpf.hk.array)), 
                                cmap='jet')
        
        hk_plot.set_clim( np.nanmin(lpf.hk.array), np.nanmax(lpf.hk.array) )
        
        plt_nn.plot_inactive()
        
        ### ------------------ Vertical Hydraulic Conductivity -------------------
        
        plt_oo = flopy.plot.PlotMapView( model=mf_obj, layer=lll, ax=axes3[lll] )
        plt_oo.plot_grid( linewidth=0.2 )
        
        xkcd_2 = np.where( ibound_ad > 0, vk_in_model[lll], np.nan )
        
        # hk_plot = nn.plot_array( xkcd_2, cmap='jet')
        
        vk_plot = plt_oo.plot_array( xkcd_2, 
                                norm=mpl.colors.LogNorm(vmin=np.nanmin(vk_in_model), 
                                                        vmax=np.nanmax(vk_in_model)), 
                                cmap='jet')
        
        vk_plot.set_clim( np.nanmin(vk_in_model), np.nanmax(vk_in_model) )
        
        plt_oo.plot_inactive()
        
    # add a bit of white space between subplots    
    plt.subplots_adjust( wspace=0.2 )
    plt.tight_layout()
    
    
    #%% PLOTTING, PLOTTING, PLOTTING - CROSS SECTIONS
    
    row_xs = 90
    col_xs = 70
    
    ### ------------------------------ BY ROW X-SEC ------------------------------
    # creating the cross-section object 
    fig4 = plt.figure(figsize=(14, 8))
    
    # deprecated
    # xsec=flopy.plot.crosssection.ModelCrossSection(model= mf_obj, dis=dis, line={"row":row_xs})
    
    xsec = flopy.plot.PlotCrossSection(model=mf_obj, line={'Row': row_xs})
    
    grd_plots = xsec.plot_grid(linewidths=0.25)
        
    xsec.plot_array(pct_array)
    xsec.plot_ibound()
    
    # # defining the data to displayed in the cross-section
    # z = xsec.plot_surface(head_2[1]/3.28) # the /3.28 converts to meters
    
    # labeling the plot
    plt.xlabel('Row width (ft)',fontsize = 12)
    plt.ylabel('Elevation (ft)',fontsize = 12)
    plt.title("Cross section at Row " + str(row_xs))
    
    # adding East/West labels to plot
    axes4  = fig4.add_subplot(111)
    plt.text( axes4.get_xticks()[0]+1000, axes4.get_yticks()[-2], 'West', 
              va="center", ha='left', size=16 )
    plt.text( axes4.get_xticks()[-2], axes4.get_yticks()[-2], 'East', 
              va="center", ha='right', size=16 )
    
    ### ------------------------------ BY COL X-SEC ------------------------------
    # creating the cross-section object 
    fig5 = plt.figure(figsize=(18, 8))
    
    # deprecated:
    # xsec=flopy.plot.crosssection.ModelCrossSection(model= mf_obj, dis=dis, line={"column":col_xs})
    
    xsec = flopy.plot.PlotCrossSection(model=mf_obj, line={'Column': col_xs})
    
    grd_plots = xsec.plot_grid(linewidths=0.25)
    
    xsec.plot_array(hk_array)
    xsec.plot_ibound()
    
    # # defining the data to displayed in the cross-section
    # z = xsec.plot_surface(head_2[1]/3.28) # the /3.28 converts to meters
    
    # labeling the plot
    plt.xlabel('Column length (ft)',fontsize = 12)
    plt.ylabel('Elevation (ft)',fontsize = 12)
    plt.title("Cross section at Column " + str(col_xs))
    
    # adding East/West labels to plot
    axes5  = fig5.add_subplot(111)
    plt.text( axes5.get_xticks()[0]+1000, axes5.get_yticks()[-2], 'North', 
              va="center", ha='left', size=16 )
    plt.text( axes5.get_xticks()[-2], axes5.get_yticks()[-2], 'South', 
              va="center", ha='right', size=16 )


#%% HOW LONG IS THIS GONNA TAKE ME?? (INPUTS)

# how long did it take?
end_time_inputs = time.time()
print("\nWriting ALL the inputs took {0:0.2f} minutes\n"
      .format((end_time_inputs - start_time)/60))

#%% TIME TO RUN THIS

success, mfoutput = mf_obj.run_model(silent=False, pause=False)
if not success:
    raise Exception('MODFLOW did not terminate normally.')

import flopy.utils.binaryfile as bf

headobj = bf.HeadFile('../MODFLOW_FloPy-Model Files/{}_ESL_RGWFM/'
                 .format(folder_date)+model_name+'.hds')
head = headobj.get_data()

# plotting
fig6, axes6 = plt.subplots(1,nlay, figsize=(36,8))
for lll in range( nlay ):

    ### -------------------------------- HEADS --------------------------------
    
    modelmap = flopy.plot.ModelMap(model=mf_obj, layer=lll, ax=axes6[lll])
    
    lc = modelmap.plot_grid( linewidth=0.2 )
    # qm = modelmap.plot_bc('GHB', alpha=0.5)
    cs = modelmap.contour_array(head)
    plt.clabel(cs, inline=1, fontsize=10, fmt='%1.1f')
    qm = modelmap.plot_ibound()
    
# add a bit of white space between subplots    
plt.subplots_adjust( wspace=0.2 )
plt.tight_layout()




#%% HOW LONG IS THIS GONNA TAKE ME?? (RUNNING)

# how long did it take?
end_time_running = time.time()
print("\nRunning the model took {0:0.2f} minutes\n"
      .format((end_time_running - end_time_inputs)/60))


















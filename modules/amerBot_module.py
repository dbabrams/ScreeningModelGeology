 # -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 14:29:33 2019

This is the module of functions for the East St. Louis area regional 
groundwater flow model (RGWFM).

@author: alljones
"""

import numpy as np
import os
import shapely
import rasterio

#%% CREATE FUNCTION TO DETERMINE POINTS WITHIN PROVIDED POLYGON

def find_cells_within_polygon( polygon, gridx, gridy ):
    # create a boolean array of 'False' values (i.e., 0)
    pts_in = np.zeros( gridx.shape, dtype='bool')

    # get bounds of polygon <- [west, south, east, North]
    bbox = polygon.bounds
    # store the index of the array location where gridx, gridy are in bounds
    bbox_idx = np.where( (gridx >= bbox[0]) & (gridx <= bbox[2]) & \
                         (gridy >= bbox[1]) & (gridy <= bbox[3]) )

    for rowi, coli in zip( bbox_idx[0], bbox_idx[1] ):
        # check if polygon contains point -> boolean
        bool_val = polygon.contains(
                shapely.geometry.Point(gridx[rowi,coli],
                                       gridy[rowi,coli]) )
        if bool_val: # if true:
            # store True value
            pts_in[rowi, coli] = bool_val

    # ensuring points are found.
    if np.any( pts_in ):
        return pts_in
    else:
        return []
        print('Function did not find any cells within polygon.')


#%% FUNCTION TO ASSIGN ELEVATIONS TO THE MODEL

def assign_prop_by_raster( filepath_to_rasters, sgr_output):
    
    # build get_node_coordinates analog for locations at which to read elevations
    row_coor_analog = np.zeros((np.shape(sgr_output.delc)))
    for iii in range(len(sgr_output.delc)):
        if iii == 0:
            row_coor_analog[iii] =  sgr_output.cbds['miny'] + (sgr_output.delc[iii]/2) # sgr_output.cbds['miny'] +
        else:
            row_coor_analog[iii] = ( row_coor_analog[iii-1] + 
                                    (np.flipud(sgr_output.delc)[iii-1]/2) +
                                    (np.flipud(sgr_output.delc)[iii]/2))
    row_coor_analog = np.flipud(row_coor_analog)
    
    col_coor_analog = np.zeros((np.shape(sgr_output.delr)))
    for iii in range(len(sgr_output.delr)):
        if iii == 0:
            col_coor_analog[iii] = sgr_output.cbds['minx'] + (sgr_output.delr[iii]/2)  # sgr_output.cbds['minx'] + 
        else:
            col_coor_analog[iii] = ( col_coor_analog[iii-1] + 
                                    (sgr_output.delr[iii-1]/2) +
                                    (sgr_output.delr[iii]/2))
    
    colmesh_analog, rowmesh_analog = np.meshgrid( col_coor_analog, 
                                                  row_coor_analog  )
    
    # count the number of files that end in '.tif'
    num_surfaces = 0
    for item in os.listdir(filepath_to_rasters):
        fn, fext = os.path.splitext(filepath_to_rasters + item)
        if fext == '.tif':
            num_surfaces = num_surfaces + 1
    
    # define 3D array for botm elevations
    # numbers of rows and columns
    nrow_ref = len(sgr_output.delc)
    ncol_ref = len(sgr_output.delr)
    elevs_ref = np.zeros((num_surfaces, nrow_ref, ncol_ref), dtype=np.float32)
    
    # Bring in the elevation rasters and assign the values to an array
    counter = 0
    for item in os.listdir(filepath_to_rasters):
        fn, fext = os.path.splitext(filepath_to_rasters + item)
        if fext == '.tif':
            
            print(item)
            
            # obtain the array from the raster file
            with rasterio.open( filepath_to_rasters + item ) as elev_rast:
                elev_array   =  np.flipud( elev_rast.read(1) ) # elev_rast.read(1)   # np.flipud( elev_rast.read(1) )
                elev_affine  = elev_rast.profile["transform"]
                
            # find nearest raster value to grid cell center
            elev_array_bounds = rasterio.transform.array_bounds(
                                elev_array.shape[0], elev_array.shape[1],
                                elev_affine)
    
            # creating the coordinates for the elevation array
            xes  = np.linspace( elev_array_bounds[0], elev_array_bounds[2],
                        num=elev_array.shape[1] )
            yes  = np.linspace( elev_array_bounds[1], elev_array_bounds[3],
                        num=elev_array.shape[0] )
            
            xxx,yyy = np.meshgrid( xes, yes, )
            
            radius = np.zeros( np.shape( xxx ) )
            
            for rrr2, ccc2 in zip( rowmesh_analog.reshape(-1), colmesh_analog.reshape(-1) ):
                
                # find indicies of raster pixels near row col grid locations
                near_idx = np.where( (xxx >= (ccc2-(elev_affine[0]/(2**0.5)))) & \
                                     (xxx <= (ccc2+(elev_affine[0]/(2**0.5)))) & \
                                     (yyy >= (rrr2-(elev_affine[0]/(2**0.5)))) & \
                                     (yyy <= (rrr2+(elev_affine[0]/(2**0.5)))) ) 
                # distance (radius) from near raster pixels to grid cell
                radius = (((xxx[near_idx]-ccc2)**2)+((yyy[near_idx]-rrr2)**2))**0.5
                
                # now find the correct layer elevation and assign it to the right
                #   spot within the elevs_ref array
                # note, "if" statements here are for the case in which there 
                #   happen to be multiple identical distances
                _, radius_counts = np.unique(radius, return_counts=True)
                if any in radius_counts == 2:
    
                    elevs_ref[ counter, np.where(row_coor_analog == rrr2)[0][0], \
                                       np.where(col_coor_analog == ccc2)[0][0] ]=\
                                       ( np.min( elev_array[ near_idx[0][ np.where(\
                                         radius == np.min( radius ) )[0][0] ],   \
                                                     near_idx[1][ np.where(      \
                                         radius == np.min( radius ) )[0][0] ] ] ) )
                                           
                elif any in radius_counts > 2:
                    raise Exception( ('You never thought it would happen '
                                      'but here we are...') )
                    
                else:
                    elevs_ref[ counter, np.where(row_coor_analog == rrr2)[0][0], \
                                       np.where(col_coor_analog == ccc2)[0][0] ]=\
                                       ( elev_array[ near_idx[0][ np.where(      \
                                         radius == np.min( radius ) )[0][0] ],   \
                                                     near_idx[1][ np.where(      \
                                         radius == np.min( radius ) )[0][0] ] ] )
            
            counter = counter + 1
    
    return elevs_ref


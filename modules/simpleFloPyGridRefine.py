# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:58:10 2020

@author: alljones
"""

# imports
import geopandas as gpd
import numpy as np
import pandas as pd

# "Simple" Grid refinement class
class SimpleGridRefine:
    
    def __get_shapefiles(self):
        
        # read in the information
        self.mdgdf = gpd.read_file( self.__mdfl ) # total model domain - GeoDataFrame
        self.adgdf = gpd.read_file( self.__adfl ) # active model domain - GeoDataFrame
        self.rsgdf = gpd.read_file( self.__rsfl ) # refinement shape - GeoDataFrame


    def __IL_RGWFM_compatible(self):        
        # intializing variables
        xoff=self.__grid_dict['IGWFM_x_offset']
        yoff=self.__grid_dict['IGWFM_y_offset']
        cs  =self.__grid_dict['IGFWM_cellsize']
        
        # determine bounds of active domain
        bds=self.adgdf.bounds
        
        # extent of model area compat. w/statewide model grid
        #   may need to be adjusted (floor/ceil, add/subtract) depending on study area 
        #   location relative to IGWFM offset location
        compatible_xmin = ( xoff + ( np.floor( ( bds.minx.min() - xoff ) / cs ) ) * cs)
        compatible_xmax = ( xoff + ( np.ceil(  ( bds.maxx.max() - xoff ) / cs ) ) * cs)
        compatible_ymin = ( yoff + ( np.ceil(  ( bds.miny.min() - yoff ) / cs ) ) * cs)
        compatible_ymax = ( yoff + ( np.floor( ( bds.maxy.max() - yoff ) / cs ) ) * cs)
        
        # create compatible bounds object
        self.cbds = pd.DataFrame( {'minx': [compatible_xmin],
                                   'maxx': [compatible_xmax],
                                   'miny': [compatible_ymin],
                                   'maxy': [compatible_ymax]} )
    
    
    def __grid_refine(self):
        
        # model width x direction
        self.Lx = self.cbds.maxx[0] - self.cbds.minx[0] # ft
        # model "height" y direction
        self.Ly = self.cbds.maxy[0] - self.cbds.miny[0] # ft
        
        # spatial discretization of each cell by row (delr_) and column (delc_)
        self.delr_init = self.__grid_dict['delta_row_initial'] # ft
        self.delc_init = self.__grid_dict['delta_col_initial'] # ft
        
        # create calculation lambda for refine_shape points
        rfshp1 = lambda a, b, c: a + ( np.floor( (b - a) / c ) * c )
        rfshp2 = lambda a, b, c: a + ( np.ceil( (b - a) / c ) * c ) 
        
        ###-----ROW REFINEMENTS-----
        
        # beginning and ending values for refinement zone compat. w/statewide model
        rz_xmin_rnd_dwn = rfshp1( self.cbds.minx[0], self.rsgdf.bounds.minx[0], self.delr_init )
        rz_xmax_rnd_up  = rfshp2( self.cbds.minx[0], self.rsgdf.bounds.maxx[0], self.delr_init )
        
        # refinment zones need to add up to these numbers
        rowref_left   = rz_xmin_rnd_dwn - self.cbds.minx[0]
        rowref_center = rz_xmax_rnd_up - rz_xmin_rnd_dwn
        rowref_right  = self.cbds.maxx[0] - rz_xmax_rnd_up
        
        # delr_left will be a variable with cell-by-cell deltas across a row 
        #   (so column widths) to the left (West) of the refinement zone
        delr_left = self.__refine_rowcol(rowref_left, side='left')
        
        # delr_right will be a variable with cell-by-cell deltas across a row 
        #   (so column widths) to the right (East) of the refinement zone
        delr_right = self.__refine_rowcol(rowref_right, side='right')
        
        # delr_center will be a variable with cell-by-cell deltas across a row 
        #   (so column widths) in the refinement zone
        if delr_right.min() == delr_left.min():
            delr_center = np.ones( int(rowref_center/delr_right.min()) )*delr_right.min()
        else:
            raise Exception( ('Minimum delta refinement not the same along eastern'+'\n'
                              'and western stretches of model domain (i.e., along rows).') )
        
        # append delr_left, delr_center, and delr_right to create one array
        self.delr = np.append(delr_left, np.append(delr_center, delr_right))
        self.ncol_ref = len(self.delr)
        
        ###-----COLUMN REFINEMENTS-----
        
        # setting refinement zone span (min and max)
        rz_ymin_rnd_dwn = rfshp1( self.cbds.miny[0], self.rsgdf.bounds.miny[0], self.delc_init )
        rz_ymax_rnd_up = rfshp2( self.cbds.miny[0], self.rsgdf.bounds.maxy[0], self.delc_init )
        
        # refinment zones need to add up to these numbers
        colref_bottom = rz_ymin_rnd_dwn - self.cbds.miny[0]
        colref_center = rz_ymax_rnd_up - rz_ymin_rnd_dwn
        colref_top    = self.cbds.maxy[0] - rz_ymax_rnd_up
        
        # reset ref min to del _init to start refinements from coarsest cell
        delc_below = self.__refine_rowcol( colref_bottom, side='bot' )
        
        # reset ref min to del _init to start refinements from coarsest cell
        delc_above = self.__refine_rowcol( colref_top, side='top' )
        
        # delc_center will be a variable with cell-by-cell deltas across a column 
        #   (so row widths) in the refinement zone
        if delc_above.min() == delc_below.min():
            delc_center = np.ones( int(colref_center/delc_above.min()) )*delc_above.min()
        else:
            raise Exception( ('Minimum delta refinement not the same along northern'+'\n'
                              'and southern stretches of model domain (i.e., along columns).') )
        
        # append delr_below, delr_center, and delr_above to create one array
        delc_ref = np.append(delc_below, np.append(delc_center, delc_above))
        
        # flip becasue flopy.modflow.ModflowDis method takes the upper left corner as
        #   the origin and the variable has been defined with the lower left as origin
        self.delc = np.flipud(delc_ref)
        self.nrow_ref = len(self.delc)
        
        # report to user the lowest refinement size
        print("\nLowest refinement size = " + str(np.min(delr_left)) + " ft\n")
        

        
    # private - refine the rows and columns within the specified shapefile
    def __refine_rowcol(self, refine_length, deltas=None, side='left'):
        """
        The following function returns the delta (i.e., change in spacing) over a 
        provided model length, given a refinement factor (refine_factor=2) and 
        initial grid spacing (indiv_delta=1250).
        """
        #----- pull pertinent information from grid dictionary
        refine_factor    =self.__grid_dict['refine_factor']
        target_refinement=self.__grid_dict['target_refinement']
        
        # start with current IL RGWFM cells size to refine from
        cur_refine_min=self.__grid_dict['IGFWM_cellsize']
        
        # use initial grid spacing
        indiv_delta=self.delc_init
        
        # determine which side of the bounding box you are refining
        if 'l' in side.lower() or 'b' in side.lower(): # checking for 'left or 'bottom'/'below'/'bot'
            side = -1
        else:
            side = 0
            
        # create the initial array of delta values
        if deltas is None:
            deltas = np.ones( int(refine_length/indiv_delta) )*indiv_delta
        
        # loop through calculations and append/adjust arrays while 
        # splitting/refining row, or column deltas.
        while cur_refine_min >= target_refinement:
            # calculate the next level of refinement
            deltas[side]   = deltas[side]/refine_factor
            # perform analysis depending upon side
            if side ==-1:
                deltas = np.append(deltas, deltas[side])
            else:
                deltas = np.insert(deltas, 0, deltas[0])
            # update check variable
            cur_refine_min = np.min(deltas)
            
        if np.sum( deltas ) != refine_length:
            raise Exception('AEJ: Refinement deltas do not equal specified refinement length.')
            
        # after successful delta refinement and array creation,
        return deltas
        

   
    # function that initializes the object
    def __init__(self, 
                 model_domain_filename='ESL_model_domain.shp', 
                 model_domain_filepath='A:/CGS/ESLgisData/Shapefiles/',
                 active_domain_filename='ESL_Active_Domain3.shp',
                 active_domain_filepath='A:/CGS/ESLgisData/Shapefiles/',
                 refine_shape_filename='ESL_RefZone1.shp',
                 refine_shape_filepath='A:/CGS/ESLgisData/Shapefiles/',
                 grid_dictionary=None):
        # establishing file location
        self.__mdfl = model_domain_filepath + model_domain_filename
        self.__adfl = active_domain_filepath + active_domain_filename
        self.__rsfl = refine_shape_filepath + refine_shape_filename
        
        # initiate the dictionary of grid variables
        if grid_dictionary is None:
            self.__grid_dict = {'IGWFM_x_offset':2440250., 
                                'IGWFM_y_offset':2394750., 
                                'IGFWM_cellsize':2500.,
                                'delta_row_initial':1250., 
                                'delta_col_initial':1250.,
                                'refine_factor':2,
                                'target_refinement':500} # < specified ft
        else:
            self.__grid_dict = grid_dictionary
        
        # grab the data from the shapefiles
        self.__get_shapefiles()
        
        # set grid compatible with the larger IL RGWFM
        self.__IL_RGWFM_compatible()
        
        # begin grid refinement
        self.__grid_refine()
        

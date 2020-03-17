# ScreeningModelGeology

Create an approximation of the geologic framework in any area of Illinois for purposes of developing a MODFLOW screening model

## ESL approach to GitHub so far...

At the start of the ESL GitHub repository Allan did some great digging on different approaches to managing branches (great resource here: https://nvie.com/posts/a-successful-git-branching-model/). The gist of it is to maintain a "master" branch, "develop" branch, and "feature" branches. The develop branch is updated with additions made via feature branches and then the master is updated periodically from the develop branch. An example: Daniel created this GitHub repository "ScreeningModelGeology" and I populated the master branch with files from the ESL repository. That gave us a place to start the master. I then made a develop branch from the master. Going forward, if I want to make a change to the model I will create a new branch and name it something like "new-wells/mike" this is the nomenclature we've been using for feature branches ('thing-being-added/name'). It helps people to see what others are working on at a glance. Once I've completed the work intended for the feature branch I'll then merge "new-wells/mike" to develop. The changes made to develop can be merged to master at any time, but often there are good checkpoints (for instance,  I could have two feature branches called "pumping-wells/mike" and "observation-wells/mike" and merge those to develop individually, but only merge develop to master once they're both done.). If you have questions about GitHub and some of the terms used here, there are a numer of people in the office who have experience with it. If you have questions about the branching structure explained here, feel free to ask Allan or I.

## Changes that will need to be made at the start

Right off the bat there will be a number of changes that need to be made to get this up and running in terms of workflow and the script. I will update the kanban project board here (https://github.com/dbabrams/ScreeningModelGeology/projects/1 it's also the 5th tab in the header on this page). I'm going to leave some of the default project cards in case people find them helpful. They can be deleted later. Things that will need to be addressed (and are in the kanban board:
 * common filepaths and folder structures for input files that are accessible by everyone working on the project
 * gathering the input files and placing them in the right place in the right format (so far this is all .tif and .shp)
 * deleting the lakes portion of the script (it is still ESL-specific) and getting the model to run successfully
 * Coverage of the 9-layer model rasters may not perfectly coincide with the desired modeling area. This is generally a small area around the perimeter of the state. In ESL this was solved by extending the westernmost values of the 9-layer model rasters westward. For the entire state a similar method could be used to extend the rasters.

## Where to find things

There are a number of filepaths and inputs that will need to be reset for the statewide 9-layer model. Here are the important things and where to find them:
 * filepaths for shapefiles for modeling: simpleFloPyGridRefine.py, under init
 * filepaths for shapefiles for plotting: americanbottoms_RGWFM_beta.py, "read in shapefiles"
 * filepaths for elevation rasters: americanbottoms_RGWFM_beta.py, "read in elevation rasters..."
 * filepaths for percent coarse rasters: americanbottoms_RGWFM_beta.py, "reading in percent coarse..."
 * inputs for refinement/grid creation: simpleFloPyGridRefine.py, under init
 * 

## Further things of note

A few thoughts come to mind:
 * the "simpleFloPyGridRefine.py" script and function can still be used even if you're not including a refinement zone. In theory you should just be able to feed it the same shapefile as you did for the active domain.
 

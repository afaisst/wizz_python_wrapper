#### THIS SCRIPT RUNS THE WIZZ #####
# This script runs the WiZZ from beginning to the end.
# The script has to be executed in the DOCKER environment and it runs
# the WiZZ scripts in the command line.
#
#
#
# NOTE: need to install these python packages in the container:
#   - matplotlib
#
#
####################################

####### IMPORTS AND PYTHON DEFINITIONS ###########

import os, sys
import numpy as np
import time

import stomp

from astropy.io import fits, ascii
from astropy.table import Table, Column, MaskedColumn, hstack, vstack

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

## Plotting stuff
mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelpad'] = 10
mpl.rcParams['xtick.major.pad'] = 7
mpl.rcParams['ytick.major.pad'] = 7
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['xtick.minor.top'] = True
mpl.rcParams['xtick.minor.bottom'] = True
mpl.rcParams['ytick.minor.left'] = True
mpl.rcParams['ytick.minor.right'] = True
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rc('text', usetex=True)
mpl.rc('font', family='serif')
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['hatch.linewidth'] = 1

def_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']

###########################

#### FUNCTIONS ###########
## To be out-sourced at some point...


def plotStompMap(mapname):
    '''
    This function plots a STOMP map
    INPUT:
        - mapname: file name of map
    OUTPUT:
        - creates PDF figure of map. Saved at location of input map
    '''

    # Read map
    stomp_table = ascii.read( mapname.replace(".map" , ".tab") )
    stomp_table["RAwidth"] = stomp_table["RAmax"] - stomp_table["RAmin"]
    stomp_table["DECwidth"] = stomp_table["DECmax"] - stomp_table["DECmin"]

    # Plot
    fig = plt.figure(figsize=(6,6))
    ax1 = fig.add_subplot(1,1,1)

    ax1.plot(stomp_table["RAcenter"][0],stomp_table["DECcenter"][0], "o" , markersize=0.1 , color= "white" , zorder=-1)

    recs = [ Rectangle(xy = (tt["RAcenter"]-(tt["RAwidth"])/2 , tt["DECcenter"]-(tt["DECwidth"])/2 ),
                    width = tt["RAwidth"],
                    height = tt["DECwidth"]
                    ) for tt in stomp_table]
    pc = PatchCollection(recs, facecolor=def_cols[0], alpha=0.5,
                            edgecolor=def_cols[0])
    ax1.add_collection(pc)
        
    ax1.set_title(mapname.split("/")[-1] , fontsize=10)
    ax1.set_xlabel(r"R.A. [deg]")
    ax1.set_ylabel(r"Decl. [deg]")
    ax1.set_xlim(150.9 , 149.3)
    ax1.set_ylim(1.55 , 2.85)

    plt.savefig(mapname.replace(".map" , ".pdf") , bbox_inches="tight" )
    plt.close()



def CreateNiceSTOMPTable(mapname):
    '''
    Creates a nice table from a STOMP map.

    INPUT:
        - mapname: file name of STOMP map
    
    OUTPUT: Table of STOMP map saved at the same location of input STOMP map
    '''


    # Load stomp map
    stomp_map = stomp.Map( mapname )

    # create a pixel vector
    pix_vect = stomp.PixelVector()

    # load pixel vector of map
    stomp_map.Pixels(pix_vect)

    # save table
    fileout = mapname.replace(".map" , ".tab")
    f = open(fileout, 'w')
    f.write("# Map table for file {}\n".format(mapname.split("/")[-1]) )
    f.write("# File created using createNiceStompTable.py\n")
    f.write("pixid \t RAcenter \t DECcenter \t RAmin \t RAmax \t DECmin \t DECmax \t Weight \t UnitSphereX \t UnitSphereY \t UnitSphereZ \t GalLat \t GalLon\n")
    for ii in range(0 , pix_vect.size() ):
        #print(pix_vect[ii].RA() )
        f.write("%g \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n" % (ii+1, pix_vect[ii].RA() ,  pix_vect[ii].DEC() ,  pix_vect[ii].RAMin() ,  pix_vect[ii].RAMax() ,  pix_vect[ii].DECMin() ,  pix_vect[ii].DECMax() ,  pix_vect[ii].Weight() , pix_vect[ii].UnitSphereX() , pix_vect[ii].UnitSphereY() , pix_vect[ii].UnitSphereZ() , pix_vect[ii].GalLat() , pix_vect[ii].GalLon() ) )

    f.close()


def SmoothSTOMPMap(mapname , res_out , smooth_out):
    '''
    Smooths a STOMP map.

    INPUT:
        - mapname: name of STOMP map to be smoothed
        - res_out: output resolution
        - smooth_out: smoothing (coverage in percent)

    OUTPUT: New smoothed map saved at the same location of input map
        - path to smoothed map is output as first argument
    '''

    # output map 
    outputmap_name = "{}_SmoothTo{}-{}.map".format(mapname , res_out ,  smooth_out)

    # smooth map
    stomp_map = stomp.Map(mapname)
    pix_vect = stomp.PixelVector()
    stomp_map.Coverage(pix_vect, res_out, True)
    output_map = stomp.Map()
    for pix in pix_vect:
        if pix.Weight() > smooth_out/100.0:
            output_map.AddPixel(pix)
    output_map.Initialize()
    output_map.Write(outputmap_name)

    # count area of original map
    area1 = []
    pix_vect1 = stomp.PixelVector()
    stomp_map.Pixels(pix_vect1)

    for pix in pix_vect1:
        area1.append(pix.Area(pix))

    print("Area of original map %f" % (sum(area1)) )


    # count area of smoothed map
    area2 = []
    pix_vect2 = stomp.PixelVector()
    output_map.Pixels(pix_vect2)

    for pix in pix_vect2:
        area2.append(pix.Area(pix))

    print("Area of smoothed map %f" % (sum(area2)) )

    return(outputmap_name)


##########################




### DEFINITIONS ####

## PATHS AND OTHERS ------
# note that all the paths have to be given relative to the
# root in the DOCKER container.

# Main directory for scripts
wizz_main_dir = "/home/the-wizz"

# Input catalogs (these are reference and  unknown catalogs)
# Later, we can also modify them in this script.
input_catalog_main_dir = "/home/work/test_data"


# General output directory
# Here, all the files are saved in sub-directories
output_main_dir = "/home/work/wizz_output"

# Process ID
# This is a unique process ID. For now
# this is entered by hand. A sub-directory is created
process_id = "test"


## INPUT CATALOGS ---------
ref_catalog_name = "REFERENCE_size50000.fits"
unknown_catalog_name = "UNKNOW_z1T1.5_size1000.fits"

## CATALOG COLUMN NAMES -------
ref_cat_ra_name = "ALPHA_J2000"
ref_cat_dec_name = "DELTA_J2000"
ref_cat_z_name = "lp_zPDF"
ref_cat_id_name = "ID_UNIQUE_REF"


unknown_cat_ra_name = "ALPHA_J2000"
unknown_cat_dec_name = "DELTA_J2000"
unknown_cat_id_name = "ID_UNIQUE_UK"


## STOMP -------
stomp_map_resolution = 1024
stomp_smoothing_res_out = 512
stomp_smoothing_percentage = 25


## Pair Maker ------
n_regions = 8
min_scale = "3,30"
max_scale = "30,300"
n_randoms = 10
z_range = [0.1,3]

#########################

###### RUN #################


## PREPARATION ----------
# This step prepares file names, directories, etc

## Directories
work_dir = os.path.join(output_main_dir , process_id)
output_stomp_dir = os.path.join(work_dir , "stompmaps")
output_catalog_dir = os.path.join(work_dir , "catalogs")
output_logs_dir = os.path.join(work_dir , "logs")

## New file names
output_stomp_map_name = "{}_r{}.map".format(ref_catalog_name.replace(".fits",""),
                                        stomp_map_resolution)

output_masked_unknown_catalog_name = "{}_MASKED_nregions{}.fits".format(unknown_catalog_name.replace(".fits",""),n_regions)
output_pair_hdf5_name = "{}_PAIR_nregion{}.hdf5".format(unknown_catalog_name.replace(".fits","") , n_regions)

## Paths to files (only the frequently used ones)
ref_catalog_path = os.path.join(input_catalog_main_dir,ref_catalog_name)
unknown_catalog_path = os.path.join(input_catalog_main_dir,unknown_catalog_name)
output_stomp_map_path = os.path.join(output_stomp_dir , output_stomp_map_name)
output_masked_unknwon_catalog_path = os.path.join(output_catalog_dir , output_masked_unknown_catalog_name)
output_pair_hdf5_path = os.path.join(output_catalog_dir , output_pair_hdf5_name)

## Create directory
if ( not os.path.exists(work_dir) ):
    print("Creating directory: {}".format(work_dir) )
    os.mkdir(work_dir)
    print("Creating directory: {}".format(output_stomp_dir))
    os.mkdir(output_stomp_dir)
    print("Creating directory: {}".format(output_catalog_dir))
    os.mkdir(output_catalog_dir)
    print("Creating directory: {}".format(output_logs_dir))
    os.mkdir(output_logs_dir)





## RUN STOMP MAP -------------
# Creates a stomp map from the input reference catalog
# For now, we simply create a mask from the reference catalog. Later
# we might want to create one for the reference *and* unknown catalog
# and then merge them.


## Create command and then run it ##
stomp_script = os.path.join(wizz_main_dir , "utility_programs" , "stomp_adapt_map.py")
log_stomp_path = os.path.join(output_logs_dir , "stomp.log")

stomp_cmd = "python {} --input_fits_file={} --output_stomp_map={} --ra_name={} --dec_name={} --resolution={} > {}".format( stomp_script ,
        ref_catalog_path,
        output_stomp_map_path,
        ref_cat_ra_name,
        ref_cat_dec_name,
        stomp_map_resolution,
        log_stomp_path
)


## Run ##
print("Running: {}".format(stomp_cmd))
time1 = time.time()
os.system(stomp_cmd)
print("Done in {:2.2f} minutes".format( (time.time()-time1)*60 ) )


## Original Map: create table, create figure
print("Original Map: create table, create figure...")
_ = CreateNiceSTOMPTable(output_stomp_map_path)
_ = plotStompMap(output_stomp_map_path)

## Smoothed Map: smooth map, create table, create figure
print("Smoothed Map: smooth, create table, create figure")
output_stomp_map_smooth_path = SmoothSTOMPMap(mapname = output_stomp_map_path,
                                                res_out = stomp_smoothing_res_out, 
                                                smooth_out=stomp_smoothing_percentage)
_ = CreateNiceSTOMPTable(output_stomp_map_smooth_path)
_ = plotStompMap(output_stomp_map_smooth_path)



print("All done.")


## MASK CATALOG ----------------------
# In this step, we mask the catalog using the STOMP
# map that we created above.

## Create command and then run it ##
masking_script = os.path.join(wizz_main_dir , "utility_programs" , "stomp_mask_catalog.py")
log_masking_path = os.path.join(output_logs_dir , "masking.log")

masking_cmd = "python {} --stomp_map={} --n_regions={} --input_fits_file={} --output_fits_file={} --ra_name={} --dec_name={} > {}".format(masking_script,
output_stomp_map_smooth_path,
n_regions,
unknown_catalog_path,
output_masked_unknwon_catalog_path,
unknown_cat_ra_name,
unknown_cat_dec_name,
log_masking_path
)

## Run ##
print("Running: {}".format(masking_cmd))
time1 = time.time()
os.system(masking_cmd)
print("Done in {:2.2f} minutes".format( (time.time()-time1)*60 ) )


## RUN PAIR MAKER -----------------------------
# Here we run the pair maker on the reference and unknown sample 
# using the STOMP mask created in the previous step.


## Create command and then run it ##
pairmaker_script = os.path.join(wizz_main_dir , "pair_maker.py")
log_pair_path = os.path.join(output_logs_dir , "pair.log")


pair_cmd = "python {} --stomp_map={} --reference_sample_file={} --reference_redshift_name={} --reference_ra_name={} --reference_dec_name={} --reference_index_name={} --unknown_sample_file={} --unknown_ra_name={} --unknown_dec_name={} --unknown_index_name={} --min_scale={} --max_scale={} --n_randoms={} --z_min={} --z_max={} --n_regions={} --output_pair_hdf5_file={} > {}".format(pairmaker_script,
output_stomp_map_smooth_path,
ref_catalog_path,
ref_cat_z_name,
ref_cat_ra_name,
ref_cat_dec_name,
ref_cat_id_name,
output_masked_unknwon_catalog_path,
unknown_cat_ra_name,
unknown_cat_dec_name,
unknown_cat_id_name,
min_scale,
max_scale,
n_randoms,
z_range[0],
z_range[1],
n_regions,
output_pair_hdf5_path,
log_pair_path
)

## Run ##
print("Running: {}".format(pair_cmd))
time1 = time.time()
os.system(pair_cmd)
print("Done in {:2.2f} minutes".format( (time.time()-time1)*60 ) )


print("ALL DONE!!!")
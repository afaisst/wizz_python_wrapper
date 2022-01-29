### THIS IS AN EXAMPLE SCRIPT ON HOW TO RUN THE-WIZZ ####
# By Andreas Faisst (afaisst@caltech.edu)
# Versions:
# - Jan 28, 2022: First version


### IMPORT #####
import os, sys
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

from WIZZclasses import *

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
mpl.rc('font', family='serif')
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['hatch.linewidth'] = 1

def_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']


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
process_id = "test3"


## INPUT CATALOGS ---------
ref_catalog_name = "REFERENCE_size50000.fits"
unknown_catalog_name = "UNKNOW_z1T1.5_size1000.fits"


## STOMP -------
stomp_map_resolution = 1024
stomp_smoothing_res_out = 512
stomp_smoothing_percentage = 25


## Pair Maker ------
n_regions = 8
min_scale = "3,30"
max_scale = "30,300"
n_randoms = 10
z_range = [0.1,6]


## PDF Maker -------
z_binning_type = "logspace"  # linear | logspace
z_n_bins = 40
n_processes = 10
unknown_stomp_region_name = "None"
use_inverse_weighting = True
unknown_weight_name = "None"
n_bootstrap = 1000
bootstrap_samples = "None"
output_bootstrap_name = "None"
n_reference_load_size = 100000

plot_z_bin_width = 0.25 

#########################

###### RUN #################



## PREPARATION ----------
# This step prepares file names, directories, etc

## Directories
work_dir = os.path.join(output_main_dir , process_id)
ref_catalog_path = os.path.join(input_catalog_main_dir , ref_catalog_name)
unknown_catalog_path = os.path.join(input_catalog_main_dir , unknown_catalog_name)


## INITIALIZE WIZZ -------

wizz = WIZZ(output_main_dir=output_main_dir,
            process_id=process_id,
            wizz_main_dir=wizz_main_dir)


## RUN STOMP MAP -------------
# Creates a stomp map from the input reference catalog
# For now, we simply create a mask from the reference catalog. Later
# we might want to create one for the reference *and* unknown catalog
# and then merge them.

sm = CreateStompMap(work_dir=work_dir,
            stomp_map_resolution=stomp_map_resolution)
_ = sm.run(catalog_path=ref_catalog_path)

_ = sm.CreateNiceSTOMPTable(sm.output_stomp_map_path)
_ = sm.plotStompMap(sm.output_stomp_map_path)

_ = sm.SmoothSTOMPMap(res_out=stomp_smoothing_res_out , smooth_out=stomp_smoothing_percentage)
_ = sm.CreateNiceSTOMPTable(sm.output_smooth_map_path)
_ = sm.plotStompMap(sm.output_smooth_map_path)

print("All done.")


## MASK CATALOG ----------------------
# In this step, we mask the catalog using the STOMP
# map that we created above.

mk = MaskCatalog(work_dir=work_dir)
mk.run(stomp_map_path=sm.output_smooth_map_path ,
        catalog_path=unknown_catalog_path)





## RUN PAIR MAKER -----------------------------
# Here we run the pair maker on the reference and unknown sample 
# using the STOMP mask created in the previous step.


pm = RunPairMaker(work_dir=work_dir , z_range=z_range)
pm.run(ref_catalog_path = ref_catalog_path ,
        unknown_catalog_path=mk.output_masked_catalog_path ,
        stomp_map_path = sm.output_smooth_map_path)
print(pm.z_range)

print("Output HDF5 file: {}".format(pm.output_pair_hdf5_path))



## RUN PDF MAKER -----------------------------
# This will create the redshift PDF of for the unknown sources.
# Later, when we make the code more modular, we can select also
# a sub-sample of the catalog (for example different stellar
# masses or redshifts).
pdfm = RunPDFMaker(work_dir=work_dir, z_range=z_range)
_ = pdfm.run(unknown_catalog_path=mk.output_masked_catalog_path , pair_hdf5_path=pm.output_pair_hdf5_path)
pdfm.PlotPDF(plot_z_range=z_range , plot_z_bin_width=plot_z_bin_width)



print("ALL DONE!!!")
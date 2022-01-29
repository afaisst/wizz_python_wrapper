#### THESE ARE CLASSES AND FUNCTIONS FOR THE-WIZZ PYTHON WRAPPER ####
# By Andreas Faisst (afaisst@caltech.edu)
# Versions:
# - Jan 28, 2022: First version

import os, sys
import numpy as np
import time

import stomp

from astropy.io import fits, ascii
from astropy.table import Table, Column, MaskedColumn, hstack, vstack

import matplotlib
matplotlib.use('Agg') # add this so that plotting works
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
mpl.rc('font', family='serif')
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['hatch.linewidth'] = 1

def_cols = plt.rcParams['axes.prop_cycle'].by_key()['color']


def bin_data(x,y,yerr , bin_centers , bin_width):
    '''
    General function to bin data.

    Parameters
    ----------
    x : `float, array`
        X-values
    y : `float, array`
        Y-values
    yerr : `float, array`
        Uncertainties in y-values
    bin_centers : `float, array`
        Centers of the bins
    bin_width : `float`
        Width of the bins

    Output
    ------
    Returns a table with binned values for x, y, and yerr.

    '''
    
    tab = Table(names=["x","y","yerr"])
    for bin_center in bin_centers:
        
        sel = np.where( np.abs(x - bin_center) < bin_width/2 )[0]
        if len(sel) > 0:
            
            y_med = np.nanmedian(y[sel])
            y_err_med = np.sqrt( np.nansum( yerr[sel]**2 ) ) / np.sqrt(len(sel))**2
            
            tab.add_row([bin_center , y_med , y_err_med])
            
            
        else:
            pass
    
    return(tab)

class WIZZ(object):
    '''
    This class includes all general parameters for the WiZZ, such as directories or general variables.

    Parameters
    -----------
    output_main_dir : `string`
        Main directory for outputs
    process_id : `string`
        Process ID. This creates a sub-directory called `process_id` in the main output directory (`output_main_dir`)
    wizz_main_dir : `string`
        Main directory to where the-WiZZ code lives.
    '''
    def __init__(self,
                output_main_dir,
                process_id,
                wizz_main_dir,
                #ref_catalog_path,
                #unknown_catalog_path
                ):
        self.output_main_dir = output_main_dir
        self.process_id = process_id
        self.wizz_main_dir = wizz_main_dir
        #self.ref_catalog_path = ref_catalog_path
        #self.unknown_catalog_path = unknown_catalog_path
        
        ## Create directory names
        self.work_dir = os.path.join(self.output_main_dir , self.process_id)
        self.output_stomp_dir = os.path.join(self.output_main_dir , self.process_id , "stompmaps")
        self.output_catalog_dir = os.path.join(self.output_main_dir , self.process_id , "catalogs")
        self.output_logs_dir = os.path.join(self.output_main_dir , self.process_id , "logs")

        ## Create directory work directories
        if ( not os.path.exists(self.work_dir) ):
            print("Creating directory: {}".format(self.work_dir) )
            os.mkdir(self.work_dir)
            print("Creating directory: {}".format(self.output_stomp_dir))
            os.mkdir(self.output_stomp_dir)
            print("Creating directory: {}".format(self.output_catalog_dir))
            os.mkdir(self.output_catalog_dir)
            print("Creating directory: {}".format(self.output_logs_dir))
            os.mkdir(self.output_logs_dir)
        else:
            print("Warning: All directories exist - will overwrite")


class CreateStompMap(object):
    '''
    This class is to create a STOMP map from a catalog.

    Parameters
    -----------
    wizzobject : `object`
            The WiZZ object (contains general information)
    work_dir : `string`
        Directory where work is done. Ideally this is something like os.path.join(maindir , process_id)
    wizz_main_dir : `string`
        Main directory where the-WiZZ code lives in 
    stomp_map_resolution : `int`
        Resolution of stomp map to be created
    cat_ra_name : `string`
        Column name in catalog for R.A
    cat_dec_name : `string`
        Column name in catalog for Decl.


    List of Functions
    -----------------
    run(catalog_path) :
        runs the creation of a STOMP map from an input catalog
    CreateNiceSTOMPTable(stomp_map_path) :
        creates a nice table with RA/DEC centers and widths from a STOMP map. This table can be used for plotting the STOMP map.

    '''
    def __init__(self,
    work_dir,
    wizz_main_dir = "/home/the-wizz",
    stomp_map_resolution = 1024,
    cat_ra_name = "ALPHA_J2000",
    cat_dec_name = "DELTA_J2000"):
        #self.wizzobject = wizzobject
        self.work_dir = work_dir
        self.wizz_main_dir = wizz_main_dir
        self.stomp_map_resolution = stomp_map_resolution
        self.cat_ra_name = cat_ra_name
        self.cat_dec_name = cat_dec_name

        self.output_logs_dir = os.path.join(self.work_dir , "logs")
        self.output_stomp_dir = os.path.join(self.work_dir , "stompmaps")
        self.output_catalog_dir = os.path.join(self.work_dir , "catalogs")

        self.script_path = os.path.join(self.wizz_main_dir , "utility_programs" , "stomp_adapt_map.py")
        self.log_path = os.path.join(self.output_logs_dir , "stomp.log")
                
    def run(self , catalog_path):
        '''
        Main function to create a STOMP map from a catalog

        Parameters
        -----------
        catalog_path : `string`
            Full path to catalog from which a STOMP map should be created

        Output
        ------
        The STOMP map is saved in a directory `stompmaps` in the work directory (`work_dir`) defined in the class.
        '''
        
        self.output_stomp_map_name = "{}_r{}.map".format(catalog_path.split("/")[-1].replace(".fits",""),
                                        self.stomp_map_resolution)
        self.output_stomp_map_path = os.path.join(self.output_stomp_dir , self.output_stomp_map_name)

       
        self.cmd = "python {} --input_fits_file={} --output_stomp_map={} --ra_name={} --dec_name={} --resolution={} > {}".format( self.script_path ,
        catalog_path,
        self.output_stomp_map_path,
        self.cat_ra_name,
        self.cat_dec_name,
        self.stomp_map_resolution,
        self.log_path
        )


        ## Run ##
        print("Running: {}".format(self.cmd))
        time1 = time.time()
        os.system(self.cmd)
        print("Done in {:2.2f} minutes".format( (time.time()-time1)/60 ) )

    def CreateNiceSTOMPTable(self , stomp_map_path):
        '''
        Creates a nice human/program readable table from a STOMP map. This table can be used for plotting the STOMP map. The table includes RA/DEC centers and width for each STOMP-square.

        Parameters
        -----------
        stomp_map_path : `string`
            Full path to the STOMP map from which a table should be created.


        Output
        ------
        The table is saved in the same location as the original STOMP map.
        '''


        # Load stomp map
        stomp_map = stomp.Map( stomp_map_path )

        # create a pixel vector
        pix_vect = stomp.PixelVector()

        # load pixel vector of map
        stomp_map.Pixels(pix_vect)

        # save table
        self.table_path = stomp_map_path.replace(".map" , ".tab")
        f = open(self.table_path, 'w')
        f.write("# Map table for file {}\n".format(stomp_map_path.split("/")[-1]) )
        f.write("# File created using createNiceStompTable.py\n")
        f.write("pixid \t RAcenter \t DECcenter \t RAmin \t RAmax \t DECmin \t DECmax \t Weight \t UnitSphereX \t UnitSphereY \t UnitSphereZ \t GalLat \t GalLon\n")
        for ii in range(0 , pix_vect.size() ):
            #print(pix_vect[ii].RA() )
            f.write("%g \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n" % (ii+1, pix_vect[ii].RA() ,  pix_vect[ii].DEC() ,  pix_vect[ii].RAMin() ,  pix_vect[ii].RAMax() ,  pix_vect[ii].DECMin() ,  pix_vect[ii].DECMax() ,  pix_vect[ii].Weight() , pix_vect[ii].UnitSphereX() , pix_vect[ii].UnitSphereY() , pix_vect[ii].UnitSphereZ() , pix_vect[ii].GalLat() , pix_vect[ii].GalLon() ) )

        f.close()

    def plotStompMap(self , stomp_map_path , ra_limits=[150.9 , 149.3] , dec_limits=[1.55 , 2.85]):
        '''
        Creates a figure of a STOMP map (specifically of its table form created by `CreateNiceSTOMPTable()`)

        Parameters
        -----------
        stomp_map_path : `string`
            Full path to the STOMP map from which a table should be created.


        Output
        ------
        The PDF figure is saved in the same location as the original STOMP map.
        '''

        # Read map
        stomp_table = ascii.read( stomp_map_path.replace(".map" , ".tab") )
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
            
        ax1.set_title(stomp_map_path.split("/")[-1] , fontsize=10)
        ax1.set_xlabel(r"R.A. [deg]")
        ax1.set_ylabel(r"Decl. [deg]")
        ax1.set_xlim(ra_limits)
        ax1.set_ylim(dec_limits)

        plt.savefig(stomp_map_path.replace(".map" , ".pdf") , bbox_inches="tight" )
        plt.close()


    def SmoothSTOMPMap(self , res_out , smooth_out):
        '''
        Smooths the STOMP map creatd by `CreateStompMap.run()`.

        Parameters
        -----------
        res_out : `int`
            New resolution of the map (= output resolution)
        smooth_out : `float`
            Smoothing of map (coverage in percent)


        Output
        ------
        The smoothed STOMP map is saved in the same location as the original STOMP map.
        '''
        # output map 
        self.output_smooth_map_path = "{}_SmoothTo{}-{}.map".format(self.output_stomp_map_path.replace(".map","") , res_out ,  smooth_out)

        # smooth map
        stomp_map = stomp.Map(self.output_stomp_map_path)
        pix_vect = stomp.PixelVector()
        stomp_map.Coverage(pix_vect, res_out, True)
        output_map = stomp.Map()
        for pix in pix_vect:
            if pix.Weight() > smooth_out/100.0:
                output_map.AddPixel(pix)
        output_map.Initialize()
        output_map.Write(self.output_smooth_map_path)

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


    
class MaskCatalog(object):
    '''
    This class is used to mask an input catalog with a STOMP map.

    Parameters
    ----------
    work_dir : `string`
        Directory where work is done. Ideally this is something like os.path.join(maindir , process_id)
    wizz_main_dir : `string`
        Main directory where the-WiZZ code lives in 
    n_regions : `int`
        Number of regions to separate the area
    cat_ra_name : `string`
        Column name in catalog for R.A
    cat_dec_name : `string`
        Column name in catalog for Decl.
    
    List of Functions
    -----------------
    run(stomp_map_path , catalog_path) :
        Creates the masked catalog

    '''
    def __init__(self,
    work_dir,
    wizz_main_dir = "/home/the-wizz",
    n_regions = 8,
    cat_ra_name = "ALPHA_J2000",
    cat_dec_name = "DELTA_J2000"):
        self.work_dir = work_dir
        self.wizz_main_dir = wizz_main_dir
        self.n_regions = n_regions
        self.cat_ra_name = cat_ra_name
        self.cat_dec_name = cat_dec_name

        self.output_logs_dir = os.path.join(self.work_dir , "logs")
        self.output_stomp_dir = os.path.join(self.work_dir , "stompmaps")
        self.output_catalog_dir = os.path.join(self.work_dir , "catalogs")

        self.script_path = os.path.join(self.wizz_main_dir , "utility_programs" , "stomp_mask_catalog.py")
        self.log_path = os.path.join(self.output_logs_dir , "masking.log")




    def run(self, stomp_map_path , catalog_path):
        '''
        Creates the masked catalog.

        Parameters
        ----------
        stomp_map_path : `string`
            Full path to STOMP map
        catalog_path : `string`
            Full path to catalog from which a masked catalog should be created


        Output
        -------
        The masked catalog is saved in the directory called `catalogs` in the working directory (`work_dir` defined in the class).
        
        '''

        self.output_masked_catalog_name = "{}_MASKED_nregions{}.fits".format(catalog_path.split("/")[-1].replace(".fits",""),self.n_regions)
        self.output_masked_catalog_path = os.path.join(self.output_catalog_dir , self.output_masked_catalog_name)
        
        self.cmd = "python {} --stomp_map={} --n_regions={} --input_fits_file={} --output_fits_file={} --ra_name={} --dec_name={} > {}".format(self.script_path,
        stomp_map_path,
        self.n_regions,
        catalog_path,
        self.output_masked_catalog_path,
        self.cat_ra_name,
        self.cat_dec_name,
        self.log_path
        )

        ## Run ##
        print("Running: {}".format(self.cmd))
        time1 = time.time()
        os.system(self.cmd)
        print("Done in {:2.2f} minutes".format( (time.time()-time1)/60 ) )



class RunPairMaker(object):
    '''
    This class is used to create the pair HDF5 file (the cross-correlation between reference and unknown catalog).

    Parameters
    ----------
    work_dir : `string`
        Directory where work is done. Ideally this is something like os.path.join(maindir , process_id)
    wizz_main_dir : `string`
        Main directory where the-WiZZ code lives in 
    n_regions : `int`
        Number of regions to separate the area
    ref_cat_z_name : `string`
        Column name in reference catalog for redshift
    ref_cat_ra_name : `string`
        Column name in reference catalog for R.A
    ref_cat_dec_name : `string`
        Column name in reference catalog for Decl.
    ref_cat_id_name : `string`
        Column name in reference catalog for unique ID
    unknown_cat_ra_name : `string`
        Column name in unknown catalog for R.A
    unknown_cat_dec_name : `string`
        Column name in unknown catalog for Decl.
    unknown_cat_id_name : `string`
        Column name in unknown catalog for unique ID
    min_scale : `string'
        Comma-separated minimum correlation scale (e.g., "3,30")
    max_scale : `string`
        Comma-separated maximium correlation scale (e.g., "30,300")
    n_randoms: `int`
        Number of random samples
    z_range: `float, list[2]`
        List with redshift range like [minimum , maximum]


    List of Functions
    -----------------
    run(ref_catalog_path , unknown_catalog_path , stomp_map_path) : 
        Main function to run the cross-correlation

    '''
    def __init__(self,
    work_dir,
    wizz_main_dir = "/home/the-wizz",
    n_regions = 8,
    ref_cat_z_name = "lp_zPDF",
    ref_cat_id_name = "ID_UNIQUE_REF",
    ref_cat_ra_name = "ALPHA_J2000",
    ref_cat_dec_name = "DELTA_J2000",
    unknown_cat_ra_name = "ALPHA_J2000",
    unknown_cat_dec_name = "DELTA_J2000",
    unknown_cat_id_name = "ID_UNIQUE_UK",
    min_scale = "3,30",
    max_scale = "30,300",
    n_randoms = 10,
    z_range = [0.1,6]
    ):
        self.work_dir = work_dir
        self.wizz_main_dir = wizz_main_dir
        self.n_regions = n_regions
        self.ref_cat_z_name = ref_cat_z_name
        self.ref_cat_id_name = ref_cat_id_name
        self.ref_cat_ra_name = ref_cat_ra_name
        self.ref_cat_dec_name = ref_cat_dec_name
        self.unknown_cat_ra_name = unknown_cat_ra_name
        self.unknown_cat_dec_name = unknown_cat_dec_name
        self.unknown_cat_id_name = unknown_cat_id_name
        self.min_scale = min_scale
        self.max_scale = max_scale
        self.n_randoms = n_randoms
        self.z_range = z_range

        self.output_logs_dir = os.path.join(self.work_dir , "logs")
        self.output_stomp_dir = os.path.join(self.work_dir , "stompmaps")
        self.output_catalog_dir = os.path.join(self.work_dir , "catalogs")

        self.script_path = os.path.join(self.wizz_main_dir , "pair_maker.py")
        self.log_path = os.path.join(self.output_logs_dir , "pair.log")

    def run(self, ref_catalog_path , unknown_catalog_path , stomp_map_path):
        '''
        Creates the cross-correlation catalog.

        Parameters
        ----------
        ref_catalog_path : `string`
            Full path to reference catalog (e.g., including spectroscopic redshifts)
        unknown_catalog_path : `string`
            Full path to the unknown (masked!) catalog (e.g., containing no redshifts). Note that this can be a large catalog and later in PDFMaker it can be specified for which subset of galaxies the PDF should be calculated. Note: here, the masked catalog has to be used!
        stomp_map_path : `string`
            Full path to the STOMP map to be used.

        Output
        ------
        The cross-correlation is saved as HDF5 file in a directory called `catalogs` in the work directory (`work_dir` as defined in the class).

        '''

        self.output_pair_hdf5_name = "{}_PAIR_nregion{}.hdf5".format(unknown_catalog_path.split("/")[-1].replace(".fits","") , self.n_regions)
        self.output_pair_hdf5_path = os.path.join(self.output_catalog_dir , self.output_pair_hdf5_name)

        self.cmd = "python {} --stomp_map={} --reference_sample_file={} --reference_redshift_name={} --reference_ra_name={} --reference_dec_name={} --reference_index_name={} --unknown_sample_file={} --unknown_ra_name={} --unknown_dec_name={} --unknown_index_name={} --min_scale={} --max_scale={} --n_randoms={} --z_min={} --z_max={} --n_regions={} --output_pair_hdf5_file={} > {}".format(self.script_path,
        stomp_map_path,
        ref_catalog_path,
        self.ref_cat_z_name,
        self.ref_cat_ra_name,
        self.ref_cat_dec_name,
        self.ref_cat_id_name,
        unknown_catalog_path,
        self.unknown_cat_ra_name,
        self.unknown_cat_dec_name,
        self.unknown_cat_id_name,
        self.min_scale,
        self.max_scale,
        self.n_randoms,
        self.z_range[0],
        self.z_range[1],
        self.n_regions,
        self.output_pair_hdf5_path,
        self.log_path
        )

        ## Run ##
        print("Running: {}".format(self.cmd))
        time1 = time.time()
        os.system(self.cmd)
        print("Done in {:2.2f} minutes".format( (time.time()-time1)/60 ) )


class RunPDFMaker(object):
    '''
        This class is used to compute the redshift probability distribution functions for sources.

        Parameters
        ----------
        work_dir : `string`
            Directory where work is done. Ideally this is something like os.path.join(maindir , process_id)
        wizz_main_dir : `string`
            Main directory where the-WiZZ code lives in 
        n_regions : `int`
            Number of regions to separate the area
        min_scale : `string'
            Comma-separated minimum correlation scale (e.g., "3,30")
        max_scale : `string`
            Comma-separated maximium correlation scale (e.g., "30,300")
        unknown_cat_id_name : `string`
            Column name in unknown catalog for unique ID. Note, the redshift PDF is created for all objects in the input unknown catalog by matching this ID to the original masked input unknown catalog.
        unknown_weight_name : `string`
            ????
        unknown_stomp_region_name : `string`
            ????
        z_range : `float, list[2]`
            List with redshift range like [minimum , maximum]
        z_n_bins : `int`
            Number of redshift bins
        z_binning_type : `string`
            How to make the binning. Options are `linear` or `logspace`
        bootstrap_samples : `string`
            ????
        n_bootstrap : `int`
            Number of bootstrap samples to computed uncertainties
        output_bootstrap_name : `string`
            Output (full) path of bootstrap file (????)
        n_processes : `int`
            Number of processors to use
        use_inverse_weighting : `bool`
            Set to `True` if inverse weighting should be use (else set to `False`)
        n_reference_load_size: `int`
            How many reference source should be loaded.        


        List of Functions
        -----------------
        run(unknown_catalog_path, pair_hdf5_path) :
            Main function to create redshift PDF from HDF5 pair file and unknown catalog.
        PlotPDF(plot_z_range , plot_z_bin_width) :
            Function to plot the redshift PDF.
        _plotpdf(pdfname , z_range , z_bin_width) :
            Helper function which actually performs the plotting.
    
    '''
    def __init__(self,
    work_dir,
    wizz_main_dir = "/home/the-wizz",
    n_regions = 8,
    min_scale = "3,30",
    max_scale = "30,300",
    unknown_cat_id_name = "ID_UNIQUE_UK",
    unknown_weight_name = "None",
    unknown_stomp_region_name = "None",
    z_range = [0.1,6],
    z_n_bins = 40,
    z_binning_type = "logspace", # linear | logspace
    bootstrap_samples = "None",
    n_bootstrap = 1000,
    output_bootstrap_name = "None",
    n_processes = 10,
    use_inverse_weighting = True,
    n_reference_load_size = 100000
    ):
        self.work_dir = work_dir
        self.wizz_main_dir = wizz_main_dir
        self.n_regions = n_regions
        self.min_scale = min_scale
        self.max_scale = max_scale
        self.unknown_cat_id_name = unknown_cat_id_name
        self.unknown_weight_name = unknown_weight_name
        self.unknown_stomp_region_name = unknown_stomp_region_name
        self.z_range = z_range
        self.z_n_bins = z_n_bins
        self.z_binning_type = z_binning_type
        self.bootstrap_samples = bootstrap_samples
        self.n_bootstrap = n_bootstrap
        self.output_bootstrap_name = output_bootstrap_name
        self.n_processes = n_processes
        self.use_inverse_weighting = use_inverse_weighting
        self.n_reference_load_size = n_reference_load_size

        self.output_logs_dir = os.path.join(self.work_dir , "logs")
        self.output_stomp_dir = os.path.join(self.work_dir , "stompmaps")
        self.output_catalog_dir = os.path.join(self.work_dir , "catalogs")

        self.script_path = os.path.join(self.wizz_main_dir , "pdf_maker.py")
        self.log_path = os.path.join(self.output_logs_dir , "pdf.log")

        # set inverse weighting flag
        if (self.use_inverse_weighting):
            self.use_inverse_weighting_flag = "--use_inverse_weighting"
        else:
            self.use_inverse_weighting_flag = ""

        # set bootstrap file path
        if (self.output_bootstrap_name != "None"):
            self.output_bootstrap_path = os.path.join(self.output_catalog_dir , self.output_bootstrap_name)
        else:
            self.output_bootstrap_path = "None"

        # create array with pair scale names
        self.pair_scale_names = ["kpc{}t{}".format(a,b) for a,b in zip( self.min_scale.split(",") , self.max_scale.split(",") ) ]

    def run(self, unknown_catalog_path, pair_hdf5_path):
        '''
        Main function to create the redshift PDF.

        Parameters
        -----------
        unknown_catalog_path : `string`
            Full path to unknown (masked) catalog for which the redshift PDF should be created
        pair_hdf5_path : `string`
            Full path to HDF5 cross-correlation file created by the PairMaker.

        Output
        ------
        A set of redshift PDF text files for each cross-correlation scale. The files are saved at the same location of the masked unknown catalog. (This is usually the `catalogs` directory in the work directory defined for this class.)
        
        '''

        # create array with output PDF names and paths
        self.output_pdf_names = ["{}_PDF_pairscale{}_{}.dat".format(unknown_catalog_path.replace(".fits",""), pp , self.z_binning_type) for pp in self.pair_scale_names] 
        self.output_pdf_paths = [ os.path.join(self.output_catalog_dir , oo) for oo in self.output_pdf_names ]

        # Now run PDF maker for each pair scale
        for ii in range(len(self.pair_scale_names)):

            # create command
            self.cmd = "python {} --input_pair_hdf5_file={} --pair_scale_name={} --unknown_sample_file={} --unknown_index_name={} --unknown_weight_name={} --unknown_stomp_region_name={} --output_pdf_file_name={} --z_min={} --z_max={} --z_n_bins={} --z_binning_type={} --bootstrap_samples={} --n_bootstrap={} --output_bootstraps_file={} --n_processes={} {} --n_reference_load_size={} >> {}".format(self.script_path,
            pair_hdf5_path,
            self.pair_scale_names[ii], 
            unknown_catalog_path, # must be the masked one!
            self.unknown_cat_id_name,
            self.unknown_weight_name,
            self.unknown_stomp_region_name,
            self.output_pdf_paths[ii],
            self.z_range[0],
            self.z_range[1],
            self.z_n_bins,
            self.z_binning_type,
            self.bootstrap_samples,
            self.n_bootstrap,
            self.output_bootstrap_path,
            self.n_processes,
            self.use_inverse_weighting_flag,
            self.n_reference_load_size,
            self.log_path
            )


            # Run
            print("Running: {}".format(self.cmd))
            time1 = time.time()
            os.system(self.cmd)
            print("Done in {:2.2f} minutes".format( (time.time()-time1)/60 ) )


    def _plotpdf(self, pdfname , z_range , z_bin_width):
        '''
        Helper function to plot the redshift PDF.

        Parameters
        ----------
        pdfname : `string`
            Full path of the redshift PDF file produced by PDFMaker
        z_range : `float, list[2]`
            Redshift range to plot
        z_bin_width : `float`
            Width of redshift bins for binning

        Output
        ------
        Saves the figure at the same location as the redshift PDF file        
        '''

        ## load PDF catalog
        pdfcat = ascii.read(pdfname , 
                            names=["z","over_density","over_density_err","n_points","n_random","area","av_res"])


        ## Plot
        fig = plt.figure(figsize = (6,4))
        ax1 = fig.add_subplot(1,1,1)

        # original
        x = pdfcat["z"]
        y = pdfcat["over_density"]
        yerr = pdfcat["over_density_err"]

        # bin
        bin_centers = np.arange(z_range[0],z_range[1],z_bin_width)
        tab_bins = bin_data(x,y,yerr, bin_centers, z_bin_width)


        norm = np.trapz(y,x)
        #yerr = yerr / norm
        #y = y / norm


        ax1.axhline(0 , color="black" , linewidth=0.5 , linestyle=":")

        ax1.errorbar(x,y , fmt="o" , yerr=yerr , color="lightgray" , zorder=-1,  markersize=4)

        ax1.errorbar(tab_bins["x"],tab_bins["y"] , yerr=tab_bins["yerr"] , xerr=z_bin_width/2 ,fmt="o" , 
                    markersize=9, markeredgecolor="black" , markerfacecolor="white" , ecolor="black" , capsize=0)


        ax1.set_xlabel("Redshift")
        ax1.set_ylabel("normalized P(z)")
        ax1.set_title(pdfname.split("/")[-1] , fontsize=8)

        plt.savefig(pdfname.replace(".dat" , ".pdf") , bbox_inches="tight" )
        plt.close()


    def PlotPDF(self , plot_z_range , plot_z_bin_width):
        '''
        Function to create the redshift PDF plot (using the helper function `_plotpdf`).

        Parameters
        ----------
        plot_z_range : `float, list[2]`
            Redshift range to plot.
        plot_z_bin_width : `float`
            Redshift bin size for binning.
        '''
        for ii in range(len(self.pair_scale_names)):
            _ = self._plotpdf(pdfname=self.output_pdf_paths[ii],  z_range=plot_z_range , z_bin_width=plot_z_bin_width)

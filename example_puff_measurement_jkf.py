# -*- coding: utf-8 -*-
import windtunnel as wt
import time

# This is an example script for the use of a PuffConcentration object.
# The functionality of the PuffConcentration class is shown, as well
# as some of the functionality inherited from pandas.Dataframe.
# The PuffConcentration class returns a DataFrame, using the standard
# measurement txt-file output as input. The columns in the input file
# are expected to be time, wtref, slow FID, fast ID, release signal and
# open_rate.

start = time.time()
# Path to your data
path = '\\\\ewtl2\\work\\Johannes\Puff_Beispiele\\'
#edit 02/18/2020: new variable to specify name of csv file which contains ambient conditions data. If given dataset
#is not found in the given file, the program resosrts to the default values specified below. 
csv_file='Q2_Ambient_Conditions.csv'

# Name of your measurement
namelist = ['Q2_170_P09.txt.ts#0']
            
#Define user input variables
#Set theshold peak concentration (ppm, model scale). All puffs with a peak concentration
#less than this concentration will be ignored in the analysis.  
threshold_concentration=0#
#Set theshold dosage (ppmvs, model scale). All puffs with a total dose less
#than this concentration will be ignored in the analysis.  
threshold_dosage=0
#edit 03/18/2020: added variable 'n_exclude,' which specified how many outliers 
#to remove at the top of the measurements. Setting 'n_exclude' to None (without qutation marks) will 
#automatically select number of outliers to remove based on sample size. To turn
#off the removal of outliers, set 'n_exclude' to zero.  
n_exclude=None
#edit 02/04/2020: added variable 'time_threshold' to control dosage threshold to use for computing characteristic
#puff times. Note that a 'default' agreed-upon value for this variable 5%, however,this fails to properly
#capture the start and end times of several puffs.	
#edit 09/19/2019: added variable full_scale to determine whether to perform analysis at full scale or model scale. 
#Set full_scale to 'fs' to perform analysis at full scale, 'ms' for model scale, or 'both' for both full scale
#and model scale. 
time_threshold=0.05
if time_threshold != 0.05:
    print('Warning: threshold dosage used to compute characteristic start and end times set to '+str(100*time_threshold)+'%, which does not equal the default value of 5%. Consider using default value!')
#edit 02/04/2020: added non-dimensional mode
full_scale='nd'
#edit 09/19/2019: added a priori information necessary for full scale analysis. Potential for GUI usage
#at a futuretime
#edit 02/21/2020: added seperate input of source location (x_source, y_source, z_source) and measurement location (x_measure, y_measure, z_measure)
#edit 05/20/2020: the proposed GUI is well under development, but has been moved to a separate
#script, titled "PAPE_GUI_code_point.py."
#edit 07/23/2020: GUI has been restructured. Standalone script "PAPE_GUI_code_point.py" has 
#been abandoned in favor of the function "standard_point_analysis.py" in the windtunnel
#package. This avoids having to use the insecure method with the open and exec functions
#in the GUI script. 
x_source=0
y_source=0
z_source=0
x_measure=855.16
y_measure=176.29
z_measure=162
pressure=1009.38
temperature=23.5
#full_scale_wtref=0
wdir=0
#edit 10/21/2019: fix spelling error (calibration is not spelled with two ls)
calibration_curve=0.3 #0.3 oder 3
mass_flow_controller='X'
calibration_factor=1
scaling_factor=0.637
scale=250
ref_length=1/250
ref_height=None
gas_name='C12'
mol_weight=28.97
gas_factor=0.5
full_scale_wtref=6
full_scale_flow_rate=0.5

#edit 02/25/2020: added ability to run script in basic mode (if 'functions_mode' variable is set to 'basic'), which runs only the core features of the scipr,
#but much faster than full mode (if 'functions_mode' variable is set to 'full'), which runs all functions of the script. 
functions_mode='full'

#edit 05/31/2020: added abillity to determine y-axis range in puff plots using variable axis_range. Current options include 'auto' (whih determines y-axis limits
#automatically for each individual puff seperately), and 'same' (which sets y-axis limits from 0 to 1.05 times the maximum concentration in the time series). Potentially 
#add option to manually specify axis limits in the future. 
axis_range='auto'

#todo: add units (09/18/2019)

wt.standard_puff_analysis(path,csv_file,namelist,threshold_concentration,threshold_dosage,
n_exclude,time_threshold,full_scale,x_source,y_source,z_source,x_measure,y_measure,z_measure,
pressure,temperature,wdir,calibration_curve,mass_flow_controller,calibration_factor,
scaling_factor,scale,ref_length,ref_height,gas_name,mol_weight,gas_factor,full_scale_wtref,
full_scale_flow_rate,functions_mode,axis_range)

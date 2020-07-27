# -*- coding: utf-8 -*-

import windtunnel as wt
import numpy as np

# This is an example script for the use of a PointConcentration object.
# The functionality of the PointConcentration class is shown, as well
# as some of the functionality inherited from pandas.Dataframe.
# The PointConcentration class returns a DataFrame, using the standard
# measurement txt-file output as input. The columns in the input file
# are expected to be time, wtref, slow FID, fast ID, release signal and
# open_rate, where release signal will be ignored.

# Path to your data
path = '\\\\ewtl2\\work\\Johannes\Puff_Beispiele\\'
#edit 05/20/2020: new variable to specify name of csv file which contains ambient conditions data. If given dataset
#is not found in the given file, the program resosrts to the default values specified below. 
csv_file='Q2_Ambient_Conditions.csv'

# Name of your measurement
namelist = ['Q2_170_P09.txt.ts#0']
            
#edit 07/23/2020: added variable full_scale. Reserved for future implementation of full_scale, model scale, and
#non-dimensional mode, analogously to puff mode             
full_scale='nd'            
#edit 09/19/2019: moved a priori information to beginning of script. Potential for GUI usage
#at a futuretime
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
#edit 07/23/2020: added variable wdir for wind direction. To be implemented in future. 
wdir=0
#edit 07/23/2020: fix spelling error (calibration is not spelled with two ls)
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

#edit 07/23/2020: added variable axis_range. Reserved for future implementation of axis range specification, 
#analogously to puff mode
axis_range='auto'

wt.standard_point_analysis(path,csv_file,namelist,full_scale,x_source,y_source,z_source,x_measure,y_measure,z_measure,
pressure,temperature,wdir,calibration_curve,mass_flow_controller,calibration_factor,
scaling_factor,scale,ref_length,ref_height,gas_name,mol_weight,gas_factor,full_scale_wtref,
full_scale_flow_rate,axis_range)
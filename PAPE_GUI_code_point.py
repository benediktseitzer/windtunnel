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
#path = '//ewtl2/projects/CERN/concentration/'

# Name of your measurement
#namelist = ['Cern_C_ANE_LC_002_V2.txt.TS#0']
#edit 09/19/2019: moved a priori information to beginning of script. Potential for GUI usage
#at a futuretime
#x=855.16
#y=176.29
#z=162
#pressure=1009.38
#temperature=23.5
#calibration_curve=0.3 #0.3 oder 3
#mass_flow_controller='X'
#edit 05/13/2020: fix spelling error (calibration is not spelled with two ls´)
#calibration_factor=1
#scaling_factor=0.637
#scale=250
#ref_length=1/250
#ref_height=None
#gas_name='C12'
#mol_weight=28.97
#gas_factor=0.5
#full_scale_wtref=6
#full_scale_flow_rate=0.5

# Initialise dict object to store instances of PuffConcentration.
# If you only have one file to analyse you can remove the loops
# and replace them with:
# mydata = PuffConcentration.from(path + file)
# mydata.calc_net_concentration()
# etc.
#

conc_ts = {}
conc_ts.fromkeys(namelist)
data_dict = {}
data_dict.fromkeys(namelist)
for name in namelist:
    #edit 10/21/2019:added option to read ambient conditions from csv file    
    ambient_conditions=wt.PointConcentration.get_ambient_conditions(path=path,name=name,input_file=path+'Cern_Ambient_Conditions.csv')  
    if ambient_conditions is None:
      []
    else:
      x_source,y_source,z_source,x_measure,y_measure,z_measure,pressure,temperature,calibration_curve,mass_flow_controller,calibration_factor,scaling_factor,scale,ref_length,\
      ref_height,gas_name,mol_weight,gas_factor,full_scale_wtref,full_scale_flow_rate=wt.PointConcentration.read_ambient_conditions(ambient_conditions,name)    
    files = wt.get_files(path,name)
    conc_ts[name] = {}
    conc_ts[name].fromkeys(files)
    for file in files:
        conc_ts[name][file] = wt.PointConcentration.from_file(path + file) 
        #edit 09/19/2019: edited code to avvound for moving a priori information to beginning of script. 
        #conc_ts[name][file].ambient_conditions(x=855.16,y=176.29,z=162,pressure=1009.38,
                                               #temperature=23.5,
                                               #calibration_curve=0.3,#0.3 oder 3
                                               #mass_flow_controller='X',
                                               #calibration_factor=1)
        conc_ts[name][file].ambient_conditions(x_source=x_source,y_source=y_source,z_source=z_source,x_measure=x_measure,y_measure=y_measure,z_measure=z_measure,pressure=pressure,
                                               temperature=temperature,
                                               calibration_curve=calibration_curve,
                                               mass_flow_controller=mass_flow_controller,
                                               calibration_factor=calibration_factor)        
        #conc_ts[name][file].scaling_information(scaling_factor=0.637,scale=250,
                                              #ref_length=1/250,ref_height=None)                
        conc_ts[name][file].scaling_information(scaling_factor=scaling_factor,scale=scale,
                                              ref_length=ref_length,ref_height=ref_height)        
        #conc_ts[name][file].tracer_information(gas_name='C12',
                                               #mol_weight=28.97/1000,
                                               #gas_factor=0.5)
        conc_ts[name][file].tracer_information(gas_name=gas_name,
                                               mol_weight=mol_weight,
                                               gas_factor=gas_factor)        
        #conc_ts[name][file].full_scale_information(full_scale_wtref=6,
                                                   #full_scale_flow_rate=0.5)
        conc_ts[name][file].full_scale_information(full_scale_wtref=full_scale_wtref,
                                                   full_scale_flow_rate=full_scale_flow_rate)        
        conc_ts[name][file].convert_temperature()
        conc_ts[name][file].calc_wtref_mean()
        conc_ts[name][file].calc_model_mass_flow_rate()
        conc_ts[name][file].calc_net_concentration() 
        #edit 07/24/2019: clear all data points with a negative net_concentration. Function can be turned on or off
        conc_ts[name][file].clear_zeros()         
        conc_ts[name][file].calc_c_star()
        conc_ts[name][file].plot_hist_conc(path=path,name=name)       
        # Save full scale results in a variable.
        # to_full_scale() will only work if all
        # information necessary has already been
        # given and computed.
        data_dict[name] = conc_ts[name][file].to_full_scale()       
        # Save full scale results. Requires to_full_scale()
        #edit 10/21/2019. Save to path of data, not to installation path of windtunnel!
        conc_ts[name][file].save2file_fs(file,out_dir=path)
        # Save model scale results
        conc_ts[name][file].save2file_ms(file,out_dir=path)
        # Save average values. Requires to_full_scale()
        conc_ts[name][file].save2file_avg(file,out_dir=path)
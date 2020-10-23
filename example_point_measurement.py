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

conc_ts = {}
conc_ts.fromkeys(namelist)
conc_ts_fs = conc_ts
conc_ts_nd = conc_ts
dict_conc_ts = conc_ts
dict_conc_nd = conc_ts
dict_conc_fs = conc_ts
data_dict = {}
data_dict.fromkeys(namelist)
for name in namelist:
    # edit 10/21/2019:added option to read ambient conditions from csv file
    ambient_conditions = wt.PointConcentration.get_ambient_conditions(path=path, name=name, input_file=path + csv_file)
    if ambient_conditions is None:
        []
    else:
        x_source, y_source, z_source, x_measure, y_measure, z_measure, pressure, temperature, calibration_curve, mass_flow_controller, calibration_factor, scaling_factor, scale, ref_length, \
        ref_height, gas_name, mol_weight, gas_factor, full_scale_wtref, full_scale_flow_rate = wt.PointConcentration.read_ambient_conditions(
            ambient_conditions, name)
    files = wt.get_files(path, name)
    conc_ts[name] = {}
    conc_ts[name].fromkeys(files)
    for file in files:

        conc_ts[name][file] = wt.PointConcentration.from_file(path + file)
        # edit 09/19/2019: edited code to avvound for moving a priori information to beginning of script.
        # conc_ts[name][file].ambient_conditions(x=855.16,y=176.29,z=162,pressure=1009.38,
        # temperature=23.5,
        # calibration_curve=0.3,#0.3 oder 3
        # mass_flow_controller='X',
        # calibration_factor=1)
        conc_ts[name][file].ambient_conditions(x_source=x_source, y_source=y_source, z_source=z_source,
                                               x_measure=x_measure, y_measure=y_measure, z_measure=z_measure,
                                               pressure=pressure,
                                               temperature=temperature,
                                               calibration_curve=calibration_curve,
                                               mass_flow_controller=mass_flow_controller,
                                               calibration_factor=calibration_factor)
        # conc_ts[name][file].scaling_information(scaling_factor=0.637,scale=250,
        # ref_length=1/250,ref_height=None)
        conc_ts[name][file].scaling_information(scaling_factor=scaling_factor, scale=scale,
                                                ref_length=ref_length, ref_height=ref_height)
        # conc_ts[name][file].tracer_information(gas_name='C12',
        # mol_weight=28.97/1000,
        # gas_factor=0.5)
        conc_ts[name][file].tracer_information(gas_name=gas_name,
                                               mol_weight=mol_weight,
                                               gas_factor=gas_factor)
        # conc_ts[name][file].full_scale_information(full_scale_wtref=6,
        # full_scale_flow_rate=0.5)
        conc_ts[name][file].full_scale_information(full_scale_wtref=full_scale_wtref,
                                                   full_scale_flow_rate=full_scale_flow_rate)
        conc_ts[name][file].convert_temperature()
        conc_ts[name][file].calc_wtref_mean()
        conc_ts[name][file].calc_model_mass_flow_rate()
        conc_ts[name][file].calc_net_concentration()
        # edit 07/24/2019: clear all data points with a negative net_concentration. Function can be turned on or off
        conc_ts[name][file].clear_zeros()
        conc_ts[name][file].calc_c_star()

        # edit 07/27/2020: added options for outputting data in full-scale, model scale, and non-dimensionally.
        if full_scale == 'ms':
            dict_conc_ts = conc_ts
        elif full_scale == 'fs':
            dict_conc_ts = conc_ts_fs
            dict_conc_ts[name][file].to_full_scale()
        elif full_scale == 'nd':
            dict_conc_ts = conc_ts_nd
            dict_conc_ts[name][file].to_non_dimensional()
        else:
            print(
                "Error: invalid input for full_scale. Data can only be computed in model scale (full_scale='ms'), full scale (full_scale='fs'), or non-dimensionally (full_scale=nd).")
        dict_conc_ts[name][file].plot_hist_conc(path=path, name=name)
        # Save full scale results in a variable.
        # to_full_scale() will only work if all
        # information necessary has already been
        # given and computed.
        # data_dict[name] = conc_ts[name][file].to_full_scale()
        # Save full scale results. Requires to_full_scale()
        # edit 10/21/2019. Save to path of data, not to installation path of windtunnel!
        wt.check_directory(path + 'Point_Data\\' + name[:name.find('.')] + '\\')
        # dict_conc_ts[name][file].save2file_fs(file,out_dir=path+'Point_Data\\'+name[:name.find('.')]+'\\')
        if full_scale == 'ms':
            dict_conc_ts[name][file].save2file_ms(file, out_dir=path + 'Point_Data\\' + name[:name.find('.')] + '\\')
        elif full_scale == 'fs':
            dict_conc_ts[name][file].save2file_fs(file, out_dir=path + 'Point_Data\\' + name[:name.find('.')] + '\\')
        elif full_scale == 'nd':
            dict_conc_ts[name][file].save2file_nd(file, out_dir=path + 'Point_Data\\' + name[:name.find('.')] + '\\')
        else:
            print(
                "Error: invalid input for full_scale. Data can only be computed in model scale (full_scale='ms'), full scale (full_scale='fs'), or non-dimensionally (full_scale=nd).")
            # Save model scale results
        # conc_ts[name][file].save2file_ms(file,out_dir=path+'Puff_Data\\'+name[:name.find('.')]+'\\')
        # Save average values. Requires to_full_scale()
        # conc_ts[name][file].save2file_avg(file,out_dir=path+'Puff_Data\\'+name[:name.find('.')]+'\\')
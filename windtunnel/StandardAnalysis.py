# -*- coding: utf-8 -*-
#! /usr/bin/python3
import matplotlib.pyplot as plt
import pandas as pd
import windtunnel as wt
import numpy as np
import time
import logging

# This is an example script for the use of a PuffConcentration object.
# The functionality of the PuffConcentration class is shown, as well
# as some of the functionality inherited from pandas.Dataframe.
# The PuffConcentration class returns a DataFrame, using the standard
# measurement txt-file output as input. The columns in the input file
# are expected to be time, wtref, slow FID, fast ID, release signal and
# open_rate.

logger = logging.getLogger()
__all__ = [
'standard_puff_analysis',
'standard_point_analysis',
]

def standard_puff_analysis(path,csv_file,namelist,threshold_concentration,threshold_dosage,
n_exclude,time_threshold,full_scale,x_source,y_source,z_source,x_measure,y_measure,z_measure,
pressure,temperature,wdir,calibration_curve,mass_flow_controller,calibration_factor,
scaling_factor,scale,ref_length,ref_height,gas_name,mol_weight,gas_factor,full_scale_wtref,
full_scale_flow_rate,functions_mode,axis_range):

    """ Perform a routine data analysis of puff concentration data. Runs essentially all analysis routines 
    available in windtunnel package"""
    #edit 07/23/2020: new function, performs routine data analysis of puff data. Runs essentially all
    #analysis routines available in windtunnel package. Based largely on example_puff_analysis.py script.
    #Improved and more secure communication with GUI. Makes PAPE_GUI_code_puff.py redundant. Largely
    #replaces example_puff_analysis.py. 
    
    
    # Initialise dict object to store instances of PuffConcentration.
    # If you only have one file to analyse you can remove the loops
    # and replace them with:
    # mydata = PuffConcentration.from(path + file)
    # mydata.calc_net_concentration()
    # etc.
    #edit 09/19/2019: added dictionaries for full scale analysis
    #edit 01/14/2020: added dictionaries for non-dimensional data analysis
    conc_ts = {}
    conc_ts.fromkeys(namelist)
    conc_ts_fs=conc_ts
    conc_ts_nd=conc_ts
    dict_conc_ts=conc_ts
    dict_conc_nd=conc_ts
    #edit 08/05/2019: new dictionary called ensemble_ts, stores data for performing ensemble analysis
    ensemble_ts = {}
    ensemble_ts.fromkeys(namelist)
    ensemble_ts_fs=ensemble_ts
    ensemble_ts_nd=ensemble_ts
    #edit 08/08/2019: new dictionary to store statistics
    statistics = {}
    statistics.fromkeys(namelist)
    statistics_fs=statistics
    statistics_nd=statistics
    for name in namelist:
        #edit 10/21/2019:added option to read ambient conditions from csv file  
        ambient_conditions=wt.PuffConcentration.get_ambient_conditions(path=path,name=name,input_file=path+csv_file)  
        if ambient_conditions is None:
          []
        else:
          x_source,y_source,z_source,x_measure,y_measure,z_measure,pressure,temperature,wdir,calibration_curve,mass_flow_controller,calibration_factor,scaling_factor,scale,ref_length,\
          ref_height,gas_name,mol_weight,gas_factor,full_scale_wtref,full_scale_flow_rate=wt.PuffConcentration.read_ambient_conditions(ambient_conditions,name)
        files = wt.get_files(path, name)
        conc_ts[name] = {}
        conc_ts[name].fromkeys(files) 
        conc_ts_fs[name]=conc_ts[name]    
        #edit 08/05/2019: new dictionary called ensemble_ts, stores data for performing ensemble analysis.       
        ensemble_ts[name] = {}
        ensemble_ts[name].fromkeys(files) 
        ensemble_ts_fs[name]=ensemble_ts[name]     
        #edit 08/08/2019: new dictionary to store statistics
        statistics[name] = {}
        statistics[name].fromkeys(namelist)
        statistics_fs[name]=statistics[name]     
        for file in files:
            conc_ts[name][file] = wt.PuffConcentration.from_file(path + file)    
            #edit 09/19/2019: added calculations necessary for full scale analysis. 
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
            #edit 10/10/2019: moved offset correction to before clear zeros.
            #conc_ts[name][file].offset_correction()         
            #edit 07/24/2019: clear all data points with a negative net_concentration. Function can be turned on or off
            #conc_ts[name][file].clear_zeros()   
            conc_ts[name][file].calc_c_star()      
            
            conc_ts_fs[name][file]=conc_ts[name][file]
            #conc_ts_fs[name][file].to_full_scale()
            
            if full_scale == 'ms':           
               dict_conc_ts=conc_ts
               dict_ensemble_ts=ensemble_ts
               dict_statistics=statistics
            elif full_scale == 'fs':     
               dict_conc_ts=conc_ts_fs
               dict_conc_ts[name][file].to_full_scale()          
               dict_ensemble_ts=ensemble_ts_fs            
    #           dict_ensemble_ts[name][file].to_full_scale()             
               dict_statistics=statistics_fs
    #           dict_statistics[name][file].to_full_scale()
            elif full_scale == 'nd':     
               #edit 01/14/2020: added option to perform data anlysis in non-dimensional mode            
               dict_conc_ts=conc_ts_nd
               dict_conc_ts[name][file].to_non_dimensional()          
               dict_ensemble_ts=ensemble_ts_fs            
    #           dict_ensemble_ts[name][file].to_full_scale()             
               dict_statistics=statistics_fs
    #           dict_statistics[name][file].to_full_scale()           
            else:
               print("Error: invalid input for full_scale. Data can only be computed in model scale (full_scale='ms'), full scale (full_scale='fs'), or non-dimensionally (full_scale=nd).")
                
            dict_conc_ts[name][file].detect_begin_release_period()
            dict_conc_ts[name][file].begin_release_index_unmasked=dict_conc_ts[name][file].begin_release_index       
            dict_conc_ts[name][file].detect_end_release_period()
            dict_conc_ts[name][file].calc_release_length()
            dict_conc_ts[name][file].get_dosage()
            #dict_conc_ts[name][file].get_mean_puff()
            dict_conc_ts[name][file].detect_arrival_time(time_threshold=time_threshold) 
            dict_conc_ts[name][file].detect_leaving_time(time_threshold=time_threshold)
            dict_conc_ts[name][file].get_residence_time()
            dict_conc_ts[name][file].get_peak_concentration()          
            dict_conc_ts[name][file].get_peak_time()         
            dict_conc_ts[name][file].get_mask(threshold_concentration=threshold_concentration,threshold_dosage=threshold_dosage,n_exclude=n_exclude)   
            #edit 02/28/2020: moved get_mean_puff to after get_mask to allow referencing of mask variable in get_mean_puff variable
            dict_conc_ts[name][file].get_mean_puff()        
            dict_conc_ts[name][file].get_ascent_time()
            dict_conc_ts[name][file].get_descent_time()
            # Pass a threshold concentration to the results
            # in order to remove all puffs with a maximum
            # concentration beneath the threshold.                 
            dict_conc_ts[name][file].apply_threshold_concentration(threshold_concentration=threshold_concentration)        
           
             
                   
            # Pass a threshold dosage to the results
            # in order to remove all puffs with a total dosage
            #beneath the threshold.        
            dict_conc_ts[name][file].apply_threshold_dosage(threshold_dosage=threshold_dosage)        
            # Test each puff against the average puff of the
            # measurement. Save the results in a variable
            deviations = dict_conc_ts[name][file].check_against_avg_puff()       
            # Save output to a variable
            #edit 08/08/2019: renamed get_puff_statistics function to get_puffs. This is to avoid confusion with the 
            #new calc_puff_statistics function whicch calculates the actual statistics.  
            results = dict_conc_ts[name][file].get_puffs()
            results_keylist=results.keys().tolist()
            #edit 08/05/2019: new dictionary (same as above) and class to perform ensemble analysis
            dict_ensemble_ts[name][file] = {}
            dict_ensemble_ts[name][file].fromkeys(results_keylist)
            dict_statistics[name][file] = {}
            dict_statistics[name][file].fromkeys(results_keylist)        
            for key in results_keylist: 
                print(key)              
                dict_ensemble_ts[name][file][key]=wt.EnsembleAnalysis.from_results(results[key])  
                dict_ensemble_ts[name][file][key].ambient_conditions(x_source=x_source,y_source=y_source,z_source=z_source,x_measure=x_measure,y_measure=y_measure,z_measure=z_measure,pressure=pressure,
                                                   temperature=temperature,
                                                   calibration_curve=calibration_curve,
                                                   mass_flow_controller=mass_flow_controller,
                                                   calibration_factor=calibration_factor)        
    
                dict_ensemble_ts[name][file][key].scaling_information(scaling_factor=scaling_factor,scale=scale,
                                                  ref_length=ref_length,ref_height=ref_height)   
                dict_ensemble_ts[name][file][key].tracer_information(gas_name=gas_name,
                                                   mol_weight=mol_weight,
                                                   gas_factor=gas_factor)        
                #conc_ts[name][file].tracer_information(gas_name='C12',            
                dict_ensemble_ts[name][file][key].get_ensemble_min()         
                #edit 08/08/2019: added functions get_ensemble_max, get_ensemble_mean, get_ensemble_variance, and plot_convergence_ensemble. 
                dict_ensemble_ts[name][file][key].get_ensemble_max()                
                dict_ensemble_ts[name][file][key].get_ensemble_mean()  
                dict_ensemble_ts[name][file][key].get_ensemble_variance()              
                dict_ensemble_ts[name][file][key].plot_convergence_ensemble(key=key,name=name,path=path,full_scale=full_scale)
                
                #edit 08/12/2019: added calculation of classes. See Bachelor Thesis of Anne Philip (2010) for more details
                if functions_mode == 'basic':
                    []
                elif functions_mode == 'full':
                    dict_ensemble_ts[name][file][key].calc_class_width(n=5)
                    dict_ensemble_ts[name][file][key].calc_class_boundaries()
                    #edit 08/13/2019: added functions get_class_frequency and plot_class_statistics
                    dict_ensemble_ts[name][file][key].get_class_frequency() 
    
                    #edit 02/25/2020: added saving of full scale and non-dimensional data
                    #edit 07/23/2020: check to make sure that directory where data is to be saved exists
                    wt.check_directory(path+'Puff_Data\\'+name[:name.find('.')]+'\\')                
                    if full_scale == 'ms':           
                        dict_ensemble_ts[name][file][key].save2file_ms_ensemble(file,key,out_dir=path+'Puff_Data\\'+name[:name.find('.')]+'\\')
                    elif full_scale == 'fs':    
                        dict_ensemble_ts[name][file][key].save2file_fs_ensemble(file,key,out_dir=path+'Puff_Data\\'+name[:name.find('.')]+'\\')                
                    elif full_scale == 'nd':
                        dict_ensemble_ts[name][file][key].save2file_nd_ensemble(file,key,out_dir=path+'Puff_Data\\'+name[:name.find('.')]+'\\')  
                    else:
                        print("Error: invalid input for full_scale. Data can only be computed in model scale (full_scale='ms'), full scale (full_scale='fs'), or non-dimensionally (full_scale=nd).")
                    dict_ensemble_ts[name][file][key].plot_class_statistics(key=key,name=name,path=path,full_scale=full_scale)    
                else:
                    print("Error: invalid input for functions_mode. Program can only be run in basic mode (functions_mode='basic') or full mode (functions_mode='full').")                
                          
                #edit 08/08/2019: added calculation of statistical values
                dict_statistics[name][file][key]=wt.EnsembleAnalysis.from_results(results[key])  
                dict_statistics[name][file][key].calc_puff_statistics(x_source=x_source,y_source=y_source,z_source=z_source,x_measure=x_measure,y_measure=y_measure,z_measure=z_measure,pressure=pressure,temperature=temperature,wtref=full_scale_wtref,wdir=wdir) 
                
                             
                
            # Save DataFrame to txt file
            #dict_conc_ts[name][file].save2file(file)
            #edit 02/25/2020: added saving of full scale and non-dimensional data
            #edit 07/23/2020: check to make sure that directory where data is to be saved exists
            wt.check_directory(path+'Puff_Data\\'+name[:name.find('.')]+'\\')          
            if full_scale == 'ms':           
               dict_conc_ts[name][file].save2file_ms(file,out_dir=path+'Puff_Data\\'+name[:name.find('.')]+'\\')
            elif full_scale == 'fs':    
               dict_conc_ts[name][file].save2file_fs(file,out_dir=path+'Puff_Data\\'+name[:name.find('.')]+'\\')                
            elif full_scale == 'nd':
               dict_conc_ts[name][file].save2file_nd(file,out_dir=path+'Puff_Data\\'+name[:name.find('.')]+'\\')  
            else:
               print("Error: invalid input for full_scale. Data can only be computed in model scale (full_scale='ms'), full scale (full_scale='fs'), or non-dimensionally (full_scale=nd).")        
            # Save DataFrame to excel file
            writer = pd.ExcelWriter(path + 'test.xlsx')
            results.to_excel(writer, sheet_name='Puff Test')
            #edit 07/24/19: plot time series
            #edit 09/23/19: added proper plotting of full scale variables
            if functions_mode == 'basic':
                dict_conc_ts[name][file].plot_puff(path=path,name=name,full_scale=full_scale,n_puffs=5,axis_range=axis_range)
            elif functions_mode == 'full':
                dict_conc_ts[name][file].plot_puff(path=path,name=name,full_scale=full_scale,n_puffs='all',axis_range=axis_range)                 
            else:
                print("Error: invalid input for functions_mode. Program can only be run in basic mode (functions_mode='basic') or full mode (functions_mode='full').")          
    
            dict_conc_ts[name][file].plot_mean_puff(path=path,name=name,stats='on',dist='off',full_scale=full_scale)         
    
    # Preliminary hist plots of the results DataFrame.
    plt.figure(0)
    results['peak concentration'].plot.hist(title='Peak Concentration')
    plt.figure(1)
    results['peak time'].plot.hist(title='Peak Time')
    plt.figure(2)
    results['arrival time'].plot.hist(title='Arrival Time')
    plt.figure(3)
    results['leaving time'].plot.hist(title='Leaving Time')
    plt.figure(4)
    results['ascent time'].plot.hist(title='Ascent Time')
    plt.figure(5)
    results['descent time'].plot.hist(title='Descent Time')
    
def standard_point_analysis(path,csv_file,namelist,full_scale,x_source,y_source,z_source,x_measure,y_measure,z_measure,
pressure,temperature,wdir,calibration_curve,mass_flow_controller,calibration_factor,
scaling_factor,scale,ref_length,ref_height,gas_name,mol_weight,gas_factor,full_scale_wtref,
full_scale_flow_rate,axis_range):

    """ Perform a routine data analysis of point concentration data. Runs essentially all analysis routines 
    available in windtunnel package"""
    #edit 07/23/2020: new function, performs routine data analysis of puff data. Runs essentially all analysis
    #routines available in windtunnel package. Based largely on example_puff_analysis.py script. Improved and
    #more secure communication with GUI. Makes PAPE_GUI_code_puff.py redundant. Largely replaces
    #example_puff_analysis.py. 

    # TODO: add implementation of full_scale and axis_range variables, analogously to puff mode (07/23/2020)
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
        ambient_conditions=wt.PointConcentration.get_ambient_conditions(path=path,name=name,input_file=path+csv_file)  
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
# -*- coding: utf-8 -*-

################
'''
IMPORTS
'''
################

import numpy as np
import pandas as pd
from scipy import stats
import logging
import matplotlib.pyplot as plt
import warnings
import os

import windtunnel as wt

# supress SOURCE ID warnings by matplotlib backend
warnings.simplefilter("ignore")
# Create logger
logger = logging.getLogger()

################
'''
MAIN
'''
################

#%%#
# Specify path to data, path to wtref, output paths for plots and txt file, 
# file type for plots, name of files, scale and desired mode of analysis.
# Input paths for data and wtref with a list of names of the measurement files

# BL
# namelist = ['BA_BL_UW_001', 'BA_BL_UW_005', 'BA_BL_UW_008', 
#             'BA_BL_UW_007', 'BA_BL_UW_006', 'BA_BL_UW_009', 
#             'BA_BL_UW_004', 'BA_BL_UW_010']

# S5 long
# namelist = ['BA_S5_L_UW_007', 'BA_S5_L_UW_006', 'BA_S5_L_UW_005', 
#             'BA_S5_L_UW_004', 'BA_S5_L_UW_012', 'BA_S5_L_UW_003',
#             'BA_S5_L_UW_001', 'BA_S5_L_UW_011', 'BA_S5_L_UW_002']
# S5 short
# namelist = ['BA_S5_S_UW_001', 'BA_S5_S_UW_002', 'BA_S5_S_UW_003', 
#             'BA_S5_S_UW_004', 'BA_S5_S_UW_005', 'BA_S5_S_UW_006']
# S5 diag
# namelist = ['BA_S5_D_UW_001', 'BA_S5_D_UW_002', 'BA_S5_D_UW_003', 
#             'BA_S5_D_UW_004', 'BA_S5_D_UW_005', 'BA_S5_D_UW_006',
#             'BA_S5_D_UW_007']

# S6 long
# namelist = ['BA_S6_L_UW_001', 'BA_S6_L_UW_002', 'BA_S6_L_UW_003', 
#             'BA_S6_L_UW_004', 'BA_S6_L_UW_005', 'BA_S6_L_UW_006',
#             'BA_S6_L_UW_007']
# S6 short
# namelist = ['BA_S6_S_UW_003', 'BA_S6_S_UW_004', 'BA_S6_S_UW_005', 
#             'BA_S6_S_UW_006', 'BA_S6_S_UW_007', 'BA_S6_S_UW_008']
# S6 diag
# namelist = ['BA_S6_D_UW_001', 'BA_S6_D_UW_002', 'BA_S6_D_UW_003', 
#             'BA_S6_D_UW_004', 'BA_S6_D_UW_005', 'BA_S6_D_UW_006',
#             'BA_S6_D_UW_007']

# S7 long
# namelist = ['BA_S7_L_UW_001', 'BA_S7_L_UW_002', 'BA_S7_L_UW_003', 
#             'BA_S7_L_UW_004', 'BA_S7_L_UW_005', 'BA_S7_L_UW_006',
#             'BA_S7_L_UW_007']
# S7 short
# namelist = ['BA_S7_S_UW_001', 'BA_S7_S_UW_002', 'BA_S7_S_UW_003', 
#             'BA_S7_S_UW_004', 'BA_S7_S_UW_005', 'BA_S7_S_UW_006']
# S7 diag
# namelist = ['BA_S7_D_UW_001', 'BA_S7_D_UW_002', 'BA_S7_D_UW_003', 
#             'BA_S7_D_UW_004', 'BA_S7_D_UW_005', 'BA_S7_D_UW_006',
#             'BA_S7_D_UW_007']

# S8 UV-long
# namelist = ['BA_S8_L_UV_003', 'BA_S8_L_UV_004', 'BA_S8_L_UV_005', 
#             'BA_S8_L_UV_006', 'BA_S8_L_UV_007', 'BA_S8_L_UV_008',
#             'BA_S8_L_UV_009', 'BA_S8_L_UV_010', 'BA_S8_L_UV_011']
# S8 long
# namelist = ['BA_S8_L_UW_008', 'BA_S8_L_UW_002', 'BA_S8_L_UW_003', 
#             'BA_S8_L_UW_004', 'BA_S8_L_UW_005', 'BA_S8_L_UW_006',
#             'BA_S8_L_UW_007']
# S8 short
# namelist = ['BA_S8_S_UW_001', 'BA_S8_S_UW_002', 'BA_S8_S_UW_003', 
#             'BA_S8_S_UW_004', 'BA_S8_S_UW_005', 'BA_S8_S_UW_006']
# S8 diagonal
# namelist = ['BA_S8_D_UW_001', 'BA_S8_D_UW_002', 'BA_S8_D_UW_003', 
#             'BA_S8_D_UW_004', 'BA_S8_D_UW_005', 'BA_S8_D_UW_006',
#             'BA_S8_D_UW_007' ]

# Wiederholungsmessungen
# namelist = ['BA_compare_001', 'BA_compare_002', 'BA_compare_003']

# LONG
# namelist = ['BA_BL_UW_001', 'BA_BL_UW_long', 
#             'BA_S5_L_UW_long', 'BA_S5_S_UW_long', 'BA_S5_D_UW_long', 
#             'BA_S6_L_UW_long', 'BA_S6_S_UW_long', 'BA_S6_D_UW_long',
#             'BA_S7_L_UW_long', 'BA_S7_S_UW_long', 'BA_S7_D_UW_long', 
#             'BA_S8_L_UW_long', 'BA_S8_S_UW_long', 'BA_S8_D_UW_long']

# reference comparison
# namelist = ['BA_BL_UW_001', 'BA_BL_UW_010', 
#             'BA_S5_L_UW_002', 'BA_S5_S_UW_006', 'BA_S5_D_UW_007', 
#             'BA_S6_L_UW_007', 'BA_S6_S_UW_008', 'BA_S6_D_UW_007',
#             'BA_S7_L_UW_007', 'BA_S7_S_UW_006', 'BA_S7_D_UW_007',
#             'BA_S8_L_UW_007', 'BA_S8_S_UW_006', 'BA_S8_D_UW_007']
namelist = ['BA_BL_UW_001', 
            'BA_S5_L_UW_002',  
            'BA_S6_L_UW_007', 
            'BA_S7_L_UW_007', 
            'BA_S8_L_UW_007',
            'BA_BL_UW_010']

# namelist = ['HG_MR_M1_UV_021']
namelist = ['BA_BL_UW_001']

x_val_shift = 100.
print(x_val_shift)

# lateral
# namelist = ['BA_S6_S_UW_001', 'BA_S6_S_UW_002']

# set paths and files
experiment = 'balcony'
if len(namelist) == 1:
    config = namelist[0][3:5]
else:
    if namelist[1][-4:] == 'long':
        config = 'LONG'    
    else:
        config = namelist[0][3:5]
# config = 'CO_REF'
path = '../../experiments/balcony/{}/coincidence/timeseries/'.format(config) # path to timeseries folder
wtref_path = '../../experiments/balcony/{}/wtref/'.format(config)
mean_path = '../../experiments/balcony/{}/coincidence/mean/'.format(config)
txt_path = './postprocessed/'
#edit 06/20/19: set ref_path to none for unknown reference path
#ref_path = './plots/'
ref_path = '../wt_ref_path/'
file_type = 'png'
plot_path_0 = '../wt_results/{}/{}/'.format(experiment,config)

if os.path.exists('{}'.format(plot_path_0)):
    print('\n {} already exists \n'.format(plot_path_0))
else:
    os.mkdir('{}'.format(plot_path_0))

#set data_nd to 1 if using non-dimensional data
data_nd = 0

if data_nd == 1:
    u_err = 0.0169
    v_err = 0.0169
    u_v_err = u_err
    I_u_err = 0.0030
    I_v_err = 0.0016
    flux_err = 0.0001
    lux_err = 5.3421
elif data_nd == 0:
    u_err = 0.0169
    v_err = 0.0169
    u_v_err = u_err
    I_u_err = 0.0030
    I_v_err = 0.0016
    flux_err = 0.0001
    lux_err = 5.3421

# 1 = vertical profile
# 2 = lateral profile
# 5 = compare profiles
# 6 = local vector fields
# 7 = longitudinal profile
# 8 = statistical comparison (error evaluation) and data quality
# 9 = data quality assesment
mode = 1
if mode == 2:
    outdata_path = '../wt_outdata_lat/'# format in npz
else:
    outdata_path = '../wt_outdata/'# format in npz

# test weighting functions
mode_2 = False

# scale of the model. 
if mode == 9:
    scale = 1.
else:
    scale = 100.
#plot scatter
plot = True
scatter = True
save_data = False

# Check if all necessary output directories exist
wt.check_directory(plot_path_0)
wt.check_directory(txt_path)
time_series = {}
time_series.fromkeys(namelist)
# set plot style for all plots
plt.style.use('classic')

#edit 06/20/19: create seperate time series for equidistant data
time_series_eq = {}
time_series_eq.fromkeys(namelist)

# Gather all files into Timeseries objects, save raw timeseries
for name in namelist:
    # prepare environment
    if os.path.exists('{}{}'.format(plot_path_0,name)):
        print('\n {}{} already exists \n'.format(plot_path_0,name))
    else:
        os.mkdir('{}{}'.format(plot_path_0,name))

    files = wt.get_files(path,name)
    time_series[name] = {}
    time_series[name].fromkeys(files)
    time_series_eq[name] = {}
    time_series_eq[name].fromkeys(files)    
    for i,file in enumerate(files):
        ts = wt.Timeseries.from_file(path+file)            
        ts.get_wind_comps(path+file)
        ts.get_wtref(wtref_path,name,index=i)
        # edit 6/20/19: Assume that input data is dimensional, not non-dimensional
        if data_nd == 0:
           print('Warning: Assuming that data is dimensional. If using non-dimensional input data, set variable data_nd to 1')
           ts.nondimensionalise()
        else:
           if data_nd == 1:
              []
           else:
              print('Warning: data_nd can only be 1 (for non-dimensional input data) or 0 (for dimensional input data)')        
        #edit 06/20/19: added seperate functionto  calculate equidistant timesteps             
        ts.adapt_scale(scale)         
        ts.mask_outliers()
        ts_eq = ts
        ts_eq.calc_equidistant_timesteps()  
        ts.index = ts.t_arr         
        ts.weighted_component_mean
        ts_eq.weighted_component_mean
        ts.weighted_component_variance
        ts_eq.weighted_component_variance
        ts.mean_magnitude
        ts_eq.mean_magnitude
        ts.mean_direction
        ts_eq.mean_direction
        ts.save2file(file)     
        time_series[name][file] = ts
        time_series_eq[name][file] = ts_eq

    print('\n created Timeseries object for: {} \n'.format(name))

if files==[]:
   raise Exception('No Matching File Names Found. Please check namelist and/or path!') 
   
for name in namelist:
   # Check if positions in all files match for vertical profile
   files = wt.get_files(path, name)
   if not (mode == 2 or mode == 7):
       for i in range(np.size(files)-2):
           if time_series[name][files[i]].x != time_series[name][files[i+1]].x:
               raise Exception('Positions do not match! Check data file.')
           if time_series[name][files[i]].y != time_series[name][files[i+1]].y:
               raise Exception('Positions do not match! Check data file.')
   # Check if positions in all files match for lateral profile
   if mode == 2:
       for i in range(np.size(files)-2):
           if time_series[name][files[i]].x != time_series[name][files[i+1]].x:
               raise Exception('Positions do not match! Check data file.')
           if time_series[name][files[i]].z != time_series[name][files[i+1]].z:
               raise Exception('Positions do not match! Check data file.')
   # Check if positions in all files match for longitudinal profile
   if mode == 7:
       for i in range(np.size(files)-2):
           if time_series[name][files[i]].y != time_series[name][files[i+1]].y:
               raise Exception('Positions do not match! Check data file.')
           if time_series[name][files[i]].z != time_series[name][files[i+1]].z:
               raise Exception('Positions do not match! Check data file.')

if mode == 9:
    # # calculate particle arrival law
    for name in namelist:
        plot_path = '{}{}/'.format(plot_path_0,name)
        files = wt.get_files(path,name)
        t_arr = []
        for j,file in enumerate(files):
            print(file)
            # calculate particle arrival law
            t_arr = time_series[name][file].t_arr
            t_arr = t_arr[~np.isnan(t_arr)]
            data_rate = len(t_arr)/200.
            delta_t_arr, particle_arrival_law = wt.calc_theo_arrival_law(t_arr, 
                                                                        data_rate)
            binscenters, data_entries, popt = wt.calc_arrival_law(t_arr, data_rate)
            # plot particle arrival law
            plt.figure(len(t_arr))
            wt.plot_arrival_law(delta_t_arr, particle_arrival_law, 
                                binscenters, data_entries, popt)
            plt.savefig(plot_path + 'particleArrivalLaw_' + file[:-4] + '.' + file_type,
                                    bbox_inches='tight')
            plt.close(len(t_arr))
    print('calculated particle arrival law')

    # calculate histograms and skewness of transit-time
    for name in namelist:
        plot_path = '{}{}/'.format(plot_path_0,name)
        files = wt.get_files(path,name)
        for file in files:
            transit_time = time_series[name][file].t_transit
            plt.figure(len(transit_time))
            wt.plot_transit_time_distribution(transit_time, 
                                        wt.calc_transit_time_distribution(transit_time))
            plt.savefig(plot_path + 'hist_transit_time' + file[:-4] + '.' + file_type,
                                    dpi=150,bbox_inches='tight')
            plt.close(len(transit_time))
    print('calculated histogram of transit time')

    # calculate T_mes/t_int as additional measure of data quality 
    lux_ratio = {}
    lux_ratio.fromkeys(namelist)
    symbols = ['o', 'v', '+', '*', 'x','d','h','3']
    colors = ['r', 'g', 'b', 'c', 'm', 'k', 'y']
    for name in namelist:
        files = wt.get_files(path,name)
        T_mes = 200. * scale
        lux_ratio[name] = {}
        lux_ratio[name].fromkeys(files)
        for file in files:
            dt = time_series[name][file].t_eq[1] - \
                time_series[name][file].t_eq[0]        
            lux = wt.calc_lux_data_wght(time_series[name][file].t_transit, dt,
                            (time_series[name][file].u_eq.dropna().values*
                            time_series[name][file].wtref))
            lux_ratio[name][file] = T_mes * np.mean(time_series[name][file].u * time_series[name][file].wtref)/lux
            print('         calculation done for: {}'.format(file))    
        print('     calculation done for: {}'.format(name))
    plt.figure(3)
    for j,name in enumerate(namelist):
        plot_path = '{}{}/'.format(plot_path_0,name)    
        files = wt.get_files(path,name)
        for file in files:
            if file[-10:] == '000001.txt':
                plt.plot(time_series[name][file].z, lux_ratio[name][file],
                             marker = symbols[j], color = colors[j], label = 'x = {} m'.format(str(time_series[name][file].x*100 + x_val_shift)))
            else:
                plt.plot(time_series[name][file].z, lux_ratio[name][file],
                                        marker = symbols[j], color = colors[j])
    plt.grid()
    plt.legend()
    plt.ylabel(r'$T_{mes} \cdot \overline{u}/L_{u}^x$ (-)')
    plt.xlabel(r'$z$ (m)')
    plt.savefig(plot_path + 'ratio_lux_tmes' + file[:-4] + '.' + file_type,
                            dpi=150,bbox_inches='tight')
    plt.close(3)
    print('calculated Lux-ratios')
    
# Iniate first layer of dictionaries for results
wind_comps = {}
wind_comps.fromkeys(namelist)
wind_stats = {}
wind_stats.fromkeys(namelist)
wind_stats_nw = {}
wind_stats_nw.fromkeys(namelist)
turb_data = {}
turb_data.fromkeys(namelist)
lux_data = {}
lux_data.fromkeys(namelist)
spectra_data = {}
spectra_data.fromkeys(namelist)

if not (mode == 5 or mode == 8 or mode == 9):
    for name in namelist:

        # shift of coordinate system
        # L=>100, S=>75, D=>106, BL=>0
        if name[6] == 'D':
            x_val_shift = 106.
        elif name[6] == 'S':
            x_val_shift = 75.
        elif name[6] == 'L':
            x_val_shift = 100.
        else: 
            x_val_shift = 100.
        print(x_val_shift)

        plot_path = '{}{}/'.format(plot_path_0,name)

        if mode_2 == 'check_wght':
            heights = []
            u_mean_w = []
            u_mean_nw = []

        #edit 6/12/19: add variables n_outliers_u and n_outliers_v to keep track of number of outliers
        n_outliers_u=0
        n_outliers_v=0
        # Iniate second layer of dictionaries for results 
        files = wt.get_files(path,name)
        wind_comps[name] = {}
        wind_comps[name].fromkeys(files)
        if mode != 3:
            wind_stats[name] = {}
            wind_stats[name].fromkeys(files)
            wind_stats_nw[name] = {}
            wind_stats_nw[name].fromkeys(files)
            turb_data[name] = {}
            turb_data[name].fromkeys(files)
            lux_data[name] = {}
            lux_data[name].fromkeys(files)
            spectra_data[name] = {}
            spectra_data[name].fromkeys(files)
        for file in files:
            print('     Started data-processing for: {} \n'.format(file))        
            #edit 6/12/19: add variables n_outliers_u and n_outliers_v to keep track of number of outliers        
            n_outliers_u=n_outliers_u+time_series[name][file].n_outliers_u
            n_outliers_v=n_outliers_v+time_series[name][file].n_outliers_v
            #edit 6/12/19: add u and v outliers to determien total 
            wind_comps[name][file] = time_series[name][file].wind_comp1,\
                                    time_series[name][file].wind_comp2
        
            # Calculate mean wind quantities
            dt = time_series[name][file].t_eq[1] - time_series[name][file].t_eq[0]
            wind_stats[name][file] = wt.calc_wind_stats_wght(time_series[name][file].t_transit,time_series[name][file].u,
                                                        time_series[name][file].v)
            turb_data[name][file] = wt.calc_turb_data_wght(
                                        time_series[name][file].t_transit,
                                        time_series[name][file].u.dropna(),
                                        time_series[name][file].v.dropna())

            # compare weighted to non-weighted profiles
            if mode_2 == 'check_wght':
                heights.append((time_series[name][file].z))
                wind_stats_nw[name][file] = wt.calc_wind_stats(time_series[name][file].u, time_series[name][file].v)
                u_mean_nw.append(wind_stats_nw[name][file][1])
                u_mean_w.append(wind_stats[name][file][1])
                plt.figure(0)

                ax = plt.gca()
                for i in range(np.size(u_mean_w)):
                    if i == 1:             
                        ax.errorbar(u_mean_nw[i],heights[i], xerr=u_v_err, label = 'non-weighted', marker='o',
                                    color='navy')
                        ax.errorbar(u_mean_w[i],heights[i], xerr=u_v_err, label = 'weighted', marker='^',
                                color='darkorange')
                    else:
                        ax.errorbar(u_mean_nw[i],heights[i], xerr=u_v_err,  marker='o',
                                    color='navy')
                        ax.errorbar(u_mean_w[i],heights[i], xerr=u_v_err,  marker='^',
                                color='darkorange')
                ax.grid(True)
                ax.set_xlabel(r'$u/u_{ref}$ (-)')
                ax.set_ylabel(r'$z$ (m)')
                plt.yscale('log')
                plt.legend()
                plt.savefig(plot_path + 'wind_data_comp_wght' + name + '.' + file_type,
                            bbox_inches='tight')
                print('plotted weight-comparison')

            # edit 06/20/2019: changed script to adapt to dimensional and non-dimensional input data
            if data_nd==0:
                lux_data[name][file] = wt.calc_lux_data_wght(
                                            time_series[name][file].t_transit,
                                            dt,
                                            (time_series[name][file].u_eq.dropna().values*
                                            time_series[name][file].wtref))
            if data_nd==1:
                lux_data[name][file] = wt.calc_lux_data_wght(
                                            time_series[name][file].t_transit,
                                            dt,
                                            (time_series[name][file].u_eq.dropna().values))        
            
            if (mode == 1 or mode == 2) and scatter:
                
                # Plot scatter plot of raw data
                plt.figure(files.index(file)+100)
                wt.plots.plot_scatter(time_series[name][file].u,
                                    time_series[name][file].v)
                plt.savefig(plot_path + 'scatter_' + file[:-4] + '.' + file_type)
                
                # Plot histograms of each component
                plt.figure(files.index(file)+200)
                wt.plots.plot_hist(time_series[name][file].u)
                plt.savefig(plot_path + 'hist_u_' + file[:-4] + '.' + file_type)
                plt.figure(files.index(file)+300)
                wt.plots.plot_hist(time_series[name][file].v)
                plt.savefig(plot_path + 'hist_v_' + file[:-4] + '.' + file_type)

                
                #edit 4/7/19: changed script to use dimensional wind values for calculating spectra to make sure that frequency is dimensionless
                spectra_data[name][file] = wt.calc_spectra(
                                                    time_series_eq[name][file].u_eq.dropna()*time_series_eq[name][file].wtref,
                                                    time_series_eq[name][file].v_eq.dropna()*time_series_eq[name][file].wtref,
                                                    time_series_eq[name][file].t_eq[~np.isnan(time_series_eq[name][file].t_eq)],
                                                    time_series_eq[name][file].z)
                # Save spectra data
                np.savetxt(txt_path + 'spectra_' + file[:-4] + '.txt',
                        np.vstack((spectra_data[name][file][0],
                                    spectra_data[name][file][1],
                                    spectra_data[name][file][2],
                                    spectra_data[name][file][3])).transpose(),
                        fmt='%.8f',
                        header="dimensionless spectra - smoothed according to reduced frequency bins"+'\n'+\
                        "frequency=0 where no energy content"+'\n'+\
                        "format: standard numpy.genfromtxt()"+'\n'+\
                        "variables = \"f_sm\" \"S_uu_sm\" \"S_vv_sm\" \"S_uv_sm\" ")
                
                # Plot spectra
                plt.figure(files.index(file)+400)
                wt.plots.plot_spectra(spectra_data[name][file][0],
                                    spectra_data[name][file][1],
                                    spectra_data[name][file][2],
                                    spectra_data[name][file][3],
                                    spectra_data[name][file][4],
                                    spectra_data[name][file][5],
                                    spectra_data[name][file][6],
                                    wind_comps[name][file],
                                    time_series[name][file].z,ref_path=ref_path)
                
                plt.savefig(plot_path + 'spectra_' + file[:-4] + '.' + file_type)
                plt.close('all')
            
            print('     Finished data-processing for: {} \n'.format(file))
        # Initiate lists for all quantitites
        x = []
        y = []
        heights = []
        mean_mag = []
        u_mean = []
        u_mean_wght = []
        u_std = []
        u_std_wght = []
        v_mean = []
        v_mean_wght = []
        v_std = []
        v_std_wght = []
        I_u = []
        I_v = []
        fluxes = []
        lux = []
        wdir = []
        wtref = []
        
        # Gather all quantities for a complete profile        
        for file in files:
            mean_mag.append(time_series[name][file].mean_magnitude)
            u_mean.append(np.mean(time_series[name][file].u))
            u_std.append(np.std(time_series[name][file].u))
            wtref.append(time_series[name][file].wtref)
            if mode !=4:
                x.append((time_series[name][file].x))
                y.append((time_series[name][file].y))
                heights.append((time_series[name][file].z))
                u_mean_wght.append(time_series[name][file].weighted_component_mean[0])
                u_std_wght.append(np.sqrt(time_series[name][file].weighted_component_variance[0]))
                v_mean.append(np.mean(time_series[name][file].v))
                v_mean_wght.append(time_series[name][file].weighted_component_mean[1])
                v_std.append(np.std(time_series[name][file].v))
                v_std_wght.append(np.sqrt(time_series[name][file].weighted_component_variance[1]))
                wdir.append(time_series[name][file].mean_direction)
                I_u.append(turb_data[name][file][0])
                I_v.append(turb_data[name][file][1])
                fluxes.append(turb_data[name][file][2])
                lux.append(lux_data[name][file])
        if save_data:
            wt.check_directory(outdata_path)
            outfile = outdata_path + name + '.npz'
            np.savez(outfile, x=x,y=y,heights=heights,mean_mag=mean_mag,u_mean=u_mean,v_mean=v_mean
                    ,I_u=I_u,I_v=I_v,fluxes=fluxes,lux=lux,wtref=wtref)

        # calculate parameters z0 and alpha
        # estimate surface_height and BL_height from fluxes and mean windprofile
        if mode == 1:
            BL_height = 100.
            sfc_heights = [70., 60., 50., 40., 30., 20.]
            mode_predefine_sfc = 'boundingheights'

            for sfc_height in sfc_heights:
                # sfc_height = 100.
                sfc_l_bound = np.where(np.asarray(heights) <= 100.)
                xcen = np.mean(np.asarray(fluxes)[sfc_l_bound])
                xrange = np.abs(0.1*xcen)

                if mode_predefine_sfc == 'ten_percent':
                    sfc_layer_l = np.where(np.asarray(fluxes) > (xcen-xrange))
                    sfc_layer_u = np.where(np.asarray(fluxes) < (xcen+xrange))                    
                    sfc_layer = np.intersect1d(sfc_layer_l,sfc_layer_u)
                    BL_layer = []
                elif mode_predefine_sfc == 'between':
                    sfc_layer = np.intersect1d(sfc_layer_l,sfc_layer_u)
                    height_min = np.min(np.asarray(heights)[sfc_layer])
                    height_max = np.max(np.asarray(heights)[sfc_layer])
                    sfc_layer_l = np.where(np.asarray(heights) >= height_min)
                    sfc_layer_u = np.where(np.asarray(heights) <= height_max)   
                    sfc_layer = np.intersect1d(sfc_layer_l,sfc_layer_u)
                    BL_layer = []
                elif mode_predefine_sfc == 'boundingheights':
                    sfc_layer_l = np.where(np.asarray(heights) >= 7.)
                    sfc_layer_u = np.where(np.asarray(heights) <= sfc_height)
                    BL_layer_u = np.where(np.asarray(heights) <= 1000.)                    
                    sfc_layer = np.intersect1d(sfc_layer_l,sfc_layer_u)
                    BL_layer = np.intersect1d(sfc_layer_l,BL_layer_u)
                elif mode_predefine_sfc == 'lower':    
                    sfc_layer = []
                    BL_layer = []

                plotfluxes = False
                if plotfluxes:
                    # plot used points and check fit
                    plt.figure(8)

                    ax = plt.gca()
                    for i in range(np.size(heights)):
                        if i in sfc_layer:    
                            ax.errorbar(fluxes[i],heights[i], xerr=flux_err, 
                                        marker='o', color='blue')
                        else:
                            ax.errorbar(fluxes[i],heights[i], xerr=flux_err,  marker='o',
                                        color='grey')
                    # ax.plot(np.asarray(u_mean_wght), np.asarray(u_mean_wght)* z1+ z0)
                    ax.axvspan(xcen-xrange,xcen+xrange,facecolor='lightskyblue',
                            edgecolor='none', alpha=0.2,
                            label='10% range of low point mean')                    
                    ax.grid(True)
                    ax.set_xlabel(r'flux (-)')
                    ax.set_ylabel(r'$z$ (m)')
                    plt.yscale('log')
                    plt.legend()
                    plt.show()
                    plt.close(8)

                print('\n surface height: {}'.format(sfc_height))    
                z0, z0_err, z1 = wt.calc_z0(u_mean_wght, heights, 
                                    sfc_height=sfc_height, sfc_layer=sfc_layer)
                print('     z0 = ', z0, '+-', z0_err)
                
                alpha, alpha_err = wt.calc_alpha(u_mean_wght, heights, 
                                    BL_height = BL_height, BL=BL_layer)
                print('     a = ', alpha, ' +- ', alpha_err)
                print('     calculated roughness length and alpha. \n')
                # calculate and sort fitting functions
                z0_x, z0_y = zip(*sorted(zip(u_mean_wght, z1)))
                alpha_x2 = (np.asarray(heights)/170.)**alpha
                alpha_x2, alpha_y2 = zip(*sorted(zip(alpha_x2, heights)))

                # plot results and check fit
                plt.figure(0)
                ax = plt.gca()
                for i in range(np.size(heights)):
                    if i in sfc_layer:     
                        ax.errorbar(u_mean_wght[i],heights[i], xerr=u_v_err, 
                                    marker='o', color='blue')
                    else:
                        ax.errorbar(u_mean_wght[i],heights[i], xerr=u_v_err,  
                                    marker='o', color='grey')

                ax.plot(z0_x, z0_y, color='darkorange', linestyle='--', 
                        label=r'fit: $z_0$ = {}'.format(np.around(z0,decimals=4)))
                ax.plot(alpha_x2, alpha_y2, color='indianred', linestyle='-.',  
                        label=r'fit: $\alpha = {}$'.format(np.around(alpha,decimals=3)) 
                        + r' at $z_{ref}=$' + r'${}$'.format(170.))
                ax.grid(True)
                ax.set_xlabel(r'$u/u_{ref}$ (-)')
                ax.set_ylabel(r'$z$ (m)')
                ax.set_xlim(right=1.0)
                ax.set_ylim(1,1000.)                
                plt.yscale('log')
                plt.legend()
                # plt.show()
                plt.savefig(plot_path + 'wind_data_fit_' + str(sfc_height) + '_' + name + '.' + file_type,
                            bbox_inches='tight')
                plt.close(0)
            print(' calculated boundary layer parameters')

        # calculate empty boundary layer reference value from upstream position 
        # name and file are reference from boundary layer at 5m height (full scale)
        if mode  == 7:
            bl_wind = time_series['BA_BL_UW_001']['BA_BL_UW_001.000001.txt'].weighted_component_mean[0]
            bl_I_u, bl_I_v, bl_flux = wt.calc_turb_data_wght(
                                            time_series['BA_BL_UW_001']['BA_BL_UW_001.000001.txt'].t_transit,
                                            time_series['BA_BL_UW_001']['BA_BL_UW_001.000001.txt'].u.dropna(),
                                            time_series['BA_BL_UW_001']['BA_BL_UW_001.000001.txt'].v.dropna())
            dt = time_series['BA_BL_UW_001']['BA_BL_UW_001.000001.txt'].t_eq[1] - \
                time_series['BA_BL_UW_001']['BA_BL_UW_001.000001.txt'].t_eq[0]
            bl_lux =  wt.calc_lux_data_wght(
                            time_series['BA_BL_UW_001']['BA_BL_UW_001.000001.txt'].t_transit,
                            dt,
                            (time_series['BA_BL_UW_001']['BA_BL_UW_001.000001.txt'].u_eq.dropna().values
                            * time_series['BA_BL_UW_001']['BA_BL_UW_001.000001.txt'].wtref))

        # plot results
        if mode == 1 and plot:
            # Plot results of a vertical profile
            # Wind components
            plt.figure(0)
            ret, lgd = wt.plots.plot_winddata(mean_mag,u_mean_wght,v_mean_wght,heights,yerr=u_v_err)
            plt.savefig(plot_path + 'wind_data_' + name + '.' + file_type,
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
        
            # Wind components, logarithmic y-axis
            plt.figure(1)
            ret, lgd = wt.plots.plot_winddata_log(mean_mag,u_mean_wght,v_mean_wght,heights,xerr=u_v_err)
            plt.tight_layout()
            plt.savefig(plot_path + 'wind_data_log_' + name + '.' + file_type,
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
            
            # Turbulence intensity of the first component
            plt.figure(2)
            wt.plots.plot_turb_int(I_u, heights, xerr=I_u_err,component='I_'+ts.wind_comp1, ref_path=ref_path)
            plt.tight_layout()
            plt.savefig(plot_path + 'I_'+ts.wind_comp1+'_' + name + '.' + file_type)
        
            # Turbulence intensity of the second component
            plt.figure(3)
            wt.plots.plot_turb_int(I_v,heights,xerr=I_v_err,component='I_'+ts.wind_comp2,ref_path=ref_path)
            plt.tight_layout()
            plt.savefig(plot_path + 'I_'+ts.wind_comp2+'_' + name + '.' + file_type)
        
            # Profile of the fluxes
            plt.figure(4) 
            wt.plots.plot_fluxes(fluxes,heights,yerr=flux_err,component='w',sfc_height=40.)
            plt.tight_layout()
            plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)
        
            # Profiles of the fluxes, logarithmic y-axis
            plt.figure(5)
            wt.plots.plot_fluxes_log(fluxes,heights,yerr=flux_err,component='w',sfc_height=40.)
            plt.tight_layout()
            plt.savefig(plot_path + 'fluxes_log_' + name + '.' + file_type)
        
            # Double logarithmic profile of Lux data
            plt.figure(6)
            wt.plots.plot_lux(lux,heights, err=lux_err,component='w',ref_path=ref_path)
            plt.tight_layout()
            plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)
            plt.close('all')

        if mode == 2 and plot:
            # Results of a lateral profile
            # Wind components
            plt.figure(0)
            ret, lgd = wt.plots.plot_winddata(mean_mag,u_mean,v_mean,y,yerr=u_v_err,lat=True)
            plt.tight_layout()
            plt.savefig(plot_path + 'wind_data_' + name + '.' + file_type,
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
            
            # Turbulence intensity of the first component
            plt.figure(1)
            wt.plots.plot_turb_int(I_u, y, xerr=I_u_err,lat=True)
            plt.tight_layout()
            plt.savefig(plot_path + 'I_u_' + name + '.' + file_type)
        
            # Turbulence intensity of the second component
            plt.figure(2)
            wt.plots.plot_turb_int(I_v, y, xerr=I_v_err,component='I_w',lat=True)
            plt.tight_layout()
            plt.savefig(plot_path + 'I_v_' + name + '.' + file_type)
        
            # Profile of the fluxes
            plt.figure(3) 
            wt.plots.plot_fluxes(fluxes,y,yerr=flux_err,component='w',lat=True)
            plt.tight_layout()
            plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)
        
            # Lateral profile of Lux data
            plt.figure(4)
            wt.plots.plot_lux(lux,y,xerr=lux_err,component='w',lat=True)
            plt.tight_layout()
            plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)

        if mode == 7 and plot:
            # Results of a longitudinal profile
            x_val = np.array(x)
            x_val += x_val_shift
            # Wind components
            plt.figure(0)
            ret, lgd = wt.plots.plot_winddata(mean_mag, u_mean_wght, v_mean_wght, x_val, 
                                            yerr=u_v_err, lat=True, var_lat = 'x')
            plt.tight_layout()
            plt.hlines(bl_wind, linestyle = 'dashdot', xmin = min(x_val), xmax= max(x_val),
                             color='darkorange', label = 'Boundary layer reference')            
            plt.savefig(plot_path + 'wind_data_' + name + '.' + file_type,
                        bbox_extra_artists=(lgd,), bbox_inches='tight')
            plt.close(0)
            
            # Turbulence intensity of the first component
            plt.figure(1)
            wt.plots.plot_turb_int(I_u, x_val, xerr=I_u_err, lat=True, var_lat= 'x')
            plt.hlines(bl_I_u, linestyle = 'dashdot', xmin = min(x_val), xmax= max(x_val),
                    color='darkorange', label = 'Boundary layer reference')   
            plt.tight_layout()
            plt.savefig(plot_path + 'I_u_' + name + '.' + file_type)
            plt.close(1)

            # Turbulence intensity of the second component
            plt.figure(2)
            wt.plots.plot_turb_int(I_v, x_val, xerr=I_v_err, component='I_w', lat=True, var_lat = 'x')
            plt.hlines(bl_I_v, linestyle = 'dashdot', xmin = min(x_val), xmax= max(x_val),
                    color='darkorange', label = 'Boundary layer reference')
            plt.tight_layout()
            plt.savefig(plot_path + 'I_v_' + name + '.' + file_type)
            plt.close(2)

            # Profile of the fluxes
            plt.figure(3) 
            wt.plots.plot_fluxes(fluxes, x_val, yerr=flux_err, component='w', lat=True, var_lat = 'x')
            plt.hlines(bl_flux, linestyle = 'dashdot', xmin = min(x_val), xmax= max(x_val),
                    color='darkorange', label = 'Boundary layer reference')
            plt.tight_layout()
            plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)
            plt.close(3)

            # Lateral profile of Lux data
            plt.figure(4)
            wt.plots.plot_lux(lux, x_val, xerr=lux_err, component='w', lat=True, var_lat = 'x')
            plt.hlines(bl_lux, linestyle = 'dashdot', xmin = min(x_val), xmax= max(x_val),
                    color='darkorange', label = 'Boundary layer reference')            
            plt.tight_layout()
            plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)
            plt.close(4)

        print('\n Finished data-processing for: {} \n'.format(name))

# comparing mode for cummulative flux and wind-profile plots 
if mode == 5:
    plt.figure(7)
    # c_list = ['#00006f', '#000092', '#0000b4',
    #             '#003eda', '#007cff', '#3ebdff',
    #             '#7cffff', '#a8ebf4', '#d5d7e9']
    # c_list = ['#00006f', '#a8ebf4']    
    c_list = ['orangered', 'forestgreen', 'darkviolet', 'dodgerblue', 'darkorange', 'slategray']            
    component = 'w'
    # u_mean
    j=0
    for name in namelist:
        u_mean = []
        heights = []
        files = wt.get_files(path,name)
        turb_data[name] = {}
        turb_data[name].fromkeys(files)
        for file in files:
            u_mean.append(time_series[name][file].weighted_component_mean[0])
            heights.append((time_series[name][file].z))
            x_val = time_series[name][file].x
            x_val += x_val_shift
        ax = plt.gca()
        ax.grid(True)

        for i in range(np.size(u_mean)):
            xerror = 0.02*u_mean[i]
            if i == 1:
                # U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='o',
                #                 color=c_list[j], label=r'$u$ at $x={}$ m'.format(str(x_val)))
                if name[6] == 'L':
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='o',
                                    color=c_list[j], label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name[6] == 'S':                
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='d',
                                    color=c_list[j], label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name[6] == 'D':              
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='+',
                                    color=c_list[j], label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name == 'BA_BL_UW_010':
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='x',
                                    color=c_list[j], label=r'no balconies at x = {}'.format(x_val))
                elif name == 'BA_BL_UW_001':
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='^',
                                    color=c_list[j], label=r'boundary layer at x = {}'.format(-130))
            else:
                if name[6] == 'L':
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='o',
                                    color=c_list[j])
                elif name[6] == 'S':                
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='d',
                                    color=c_list[j])
                elif name[6] == 'D':
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='+',
                                    color=c_list[j])
                elif name == 'BA_BL_UW_010':
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='x',
                                    color=c_list[j])
                elif name == 'BA_BL_UW_001':
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='^',
                                    color=c_list[j])

            # lgd = ax.legend(U, 'u at ', bbox_to_anchor=(0.5,1.05), loc='lower center', 
            #                   borderaxespad=0., ncol=3, fontsize=16)
        j = j + 1
    ax.set_xlabel(r'$u/u_{ref}$ (-)')
    ax.set_ylabel(r'$z$ (m)')
    plt.legend()
    plt.yscale('log')    
    plt.tight_layout()
    plt.savefig(plot_path_0 + 'cummulative_wind_data_log_' + name[:-4] + '.' + file_type,
                bbox_inches='tight')
    plt.close(7)
    # flux
    plt.figure(8)
    j = 0
    for name in namelist:
        heights = []
        fluxes = []
        files = wt.get_files(path,name)
        turb_data[name] = {}
        turb_data[name].fromkeys(files)

        for file in files:
            heights.append((time_series[name][file].z))
            x_val = time_series[name][file].x
            x_val += x_val_shift
            turb_data[name][file] = wt.calc_turb_data_wght(
                                        time_series[name][file].t_transit,
                                        time_series[name][file].u.dropna(),
                                        time_series[name][file].v.dropna())            
            fluxes.append(turb_data[name][file][2])

        ax = plt.gca()
        ax.grid(True)

        for i in range(np.size(fluxes)):
            if i == 1:
                # l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='o',color=c_list[j],
                #                 label = r'fluxes at $x={}$ m'.format(str(x_val)))
                if name[6] == 'L':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='o',color=c_list[j],
                                label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name[6] == 'S':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='d',color=c_list[j],
                                label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name[6] == 'D':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='+',color=c_list[j],
                                label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='x',color=c_list[j],
                                label=r'no balconies at x = {}'.format(x_val))
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='^',color=c_list[j],
                                label=r'boundary layer at x = {}'.format(-130))
            else:
                if name[6] == 'L':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='o',color=c_list[j])
                if name[6] == 'S':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='d',color=c_list[j])
                if name[6] == 'D':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='+',color=c_list[j])
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='x',color=c_list[j])
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='^',color=c_list[j])
        j += 1
    ax.set_xlabel(r'u' + '\'' + component + '\' $\cdot$ $U_{0}^{-2}\ (-)$')
    ax.set_ylabel(r'$z$ (m)')
    plt.legend()
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(plot_path_0 + 'cummulative_flux_log_' + name[:-4] + '.' + file_type,
                bbox_inches='tight')
    plt.close(8)
    # I_u
    plt.figure(9)
    j = 0
    for name in namelist:
        heights = []
        I_u_list = []
        files = wt.get_files(path,name)
        turb_data[name] = {}
        turb_data[name].fromkeys(files)

        for file in files:
            heights.append((time_series[name][file].z))
            x_val = time_series[name][file].x
            x_val += x_val_shift
            turb_data[name][file] = wt.calc_turb_data_wght(
                                        time_series[name][file].t_transit,
                                        time_series[name][file].u.dropna(),
                                        time_series[name][file].v.dropna())            
            I_u_list.append(turb_data[name][file][0])

        ax = plt.gca()
        ax.grid(True)

        for i in range(np.size(I_u_list)):
            if i == 1:
                # l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='o',color=c_list[j],
                #                 label = r'fluxes at $x={}$ m'.format(str(x_val)))
                if name[6] == 'L':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='o',color=c_list[j],
                                    label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name[6] == 'S':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='d',color=c_list[j],
                                    label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name[6] == 'D':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='+',color=c_list[j],
                                    label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='x',color=c_list[j],
                                    label=r'no balconies at x = {}'.format(x_val))
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='^',color=c_list[j],
                                    label=r'boundary layer at x = {}'.format(-130))
            else:
                if name[6] == 'L':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='o',color=c_list[j])
                if name[6] == 'S':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='d',color=c_list[j])
                if name[6] == 'D':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='+',color=c_list[j])
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='x',color=c_list[j])
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='^',color=c_list[j])

                
        j += 1
    ax.set_xlabel(r'$I_u$ (-)')
    ax.set_ylabel(r'$z$ (m)')
    plt.legend()
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(plot_path_0 + 'cummulative_I_u_log_' + name[:-4] + '.' + file_type,
                bbox_inches='tight')
    plt.close(9)
    # I_v
    plt.figure(10)
    j = 0
    for name in namelist:
        heights = []
        I_v_list = []        
        files = wt.get_files(path,name)
        turb_data[name] = {}
        turb_data[name].fromkeys(files)

        for file in files:
            heights.append((time_series[name][file].z))
            x_val = time_series[name][file].x
            x_val += x_val_shift
            turb_data[name][file] = wt.calc_turb_data_wght(
                                        time_series[name][file].t_transit,
                                        time_series[name][file].u.dropna(),
                                        time_series[name][file].v.dropna())            
            I_v_list.append(turb_data[name][file][1])

        ax = plt.gca()
        ax.grid(True)

        for i in range(np.size(I_v_list)):
            if i == 1:
                # l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='o',color=c_list[j],
                #                 label = r'fluxes at $x={}$ m'.format(str(x_val)))
                if name[6] == 'L':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='o',color=c_list[j],
                                    label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name[6] == 'S':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='d',color=c_list[j],
                                    label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name[6] == 'D':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='+',color=c_list[j],
                                    label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='x',color=c_list[j],
                                    label=r'no balconies at x = {}'.format(x_val))
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='^',color=c_list[j],
                                    label=r'boundary layer at x = {}'.format(-130))
            else:
                if name[6] == 'L':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='o',color=c_list[j])
                if name[6] == 'S':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='d',color=c_list[j])
                if name[6] == 'D':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='+',color=c_list[j])
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='x',color=c_list[j])
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='^',color=c_list[j])
        j += 1
    ax.set_xlabel(r'$I_w$ (-)')
    ax.set_ylabel(r'$z$ (m)')
    plt.legend()
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(plot_path_0 + 'cummulative_I_w_log_' + name[:-4] + '.' + file_type,
                bbox_inches='tight')
    plt.close(10)
    # Lux
    plt.figure(11)
    j = 0
    for name in namelist:
        heights = []
        lux_list = []
        files = wt.get_files(path,name)
        turb_data[name] = {}
        turb_data[name].fromkeys(files)

        for file in files:
            heights.append((time_series[name][file].z))
            x_val = time_series[name][file].x
            x_val += x_val_shift
            dt = time_series[name][file].t_eq[1] - \
                time_series[name][file].t_eq[0]
            lux_list.append(wt.calc_lux_data_wght(
                                time_series[name][file].t_transit,
                                dt,
                                (time_series[name][file].u_eq.dropna().values
                                * time_series[name][file].wtref)))

        ax = plt.gca()
        ax.grid(True)

        for i in range(np.size(lux_list)):
            if i == 1:
                # l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='o',color=c_list[j],
                #                 label = r'fluxes at $x={}$ m'.format(str(x_val)))
                if name[6] == 'L':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='o',color=c_list[j],
                                    label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name[6] == 'S':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='d',color=c_list[j],
                                    label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name[6] == 'D':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='+',color=c_list[j],
                                    label=r'{} m spacing at x={}'.format(name[4], x_val))
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='x',color=c_list[j],
                                    label=r'no balconies at x = {}'.format(x_val))
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='^',color=c_list[j],
                                    label=r'boundary layer at x = {}'.format(-130))
            else:
                if name[6] == 'L':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='o',color=c_list[j])
                if name[6] == 'S':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='d',color=c_list[j])
                if name[6] == 'D':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='+',color=c_list[j])
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='x',color=c_list[j])
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='^',color=c_list[j])
        j += 1
        print('Lux done for {}'.format(name))
    ax.set_xlabel(r'$L_u^x$ (m)')
    ax.set_ylabel(r'$z$ (m)')
    plt.legend()
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(plot_path_0 + 'cummulative_Lux_log_' + name[:-4] + '.' + file_type,
                bbox_inches='tight')
    print(plot_path_0)
    plt.close(11)

    print('\n plotted cummulative profiles over length.\n')

# local eddy resoluting vector plots
if mode == 6:
    plt.figure(7)
    for name in namelist:
        u_mean = []
        v_mean = []
        heights = []
        files = wt.get_files(path,name)
        for file in files:
            u_mean.append(np.mean(time_series[name][file].u))
            v_mean.append(np.mean(time_series[name][file].v))
            heights.append((time_series[name][file].z))
            x_val = time_series[name][file].x + x_val_shift
            x_vals = np.full(len(heights),x_val)

        u_mean = np.asarray(u_mean)
        v_mean = np.asarray(v_mean)        
        heights = np.asarray(heights)
        stop_height = np.where(heights < 5.)
        ax = plt.gca()
        ax.grid(True)
        for i in range(np.size(u_mean)):
                U = ax.quiver(x_vals[stop_height], heights[stop_height], 
                                u_mean[stop_height], v_mean[stop_height], width = 0.005)
            # lgd = ax.legend(U, 'u at ', bbox_to_anchor=(0.5,1.05), 
            #                   loc='lower center', borderaxespad=0., ncol=3, fontsize=16)

    ax.set_xlabel(r'$x$ [m]')
    ax.set_ylabel(r'$z$ [m]')
    plt.legend()   
    plt.tight_layout()
    plt.savefig(plot_path_0 + 'local_mean_wind_' + name[:-4] + '.' + file_type,
                bbox_inches='tight')
    plt.close(7)
    print('\n plotted local-meanwind-profiles \n')

# calculate statistics from repetitive measurements
if mode == 8:
    print('\n statistical evaluation:')
    u_error = []
    u_err_val = []    
    x_val = []
    z_val = []
    lux_unc = []
    I_u_unc = []
    I_v_unc = []
    flux_unc = []
    for name in namelist:
        files = wt.get_files(path,name)
        # calculate errors
        for file in files:
            if file[15:] == '000001.txt':
                I_u1, I_v1, flux1 = wt.calc_turb_data_wght(
                                        time_series[name][file].t_transit,
                                        time_series[name][file].u.dropna(),
                                        time_series[name][file].v.dropna())
                I_u2, I_v2, flux2 = wt.calc_turb_data_wght(
                                        time_series[name][file[:14] + '.000002.txt'].t_transit,
                                        time_series[name][file[:14] + '.000002.txt'].u.dropna(),
                                        time_series[name][file[:14] + '.000002.txt'].v.dropna())
                diff_I_u = abs(I_u1-I_u2)
                diff_I_v = abs(I_v1-I_v2)
                diff_flux = abs(flux1-flux2)
                dt1 = time_series[name][file].t_eq[1] - \
                    time_series[name][file].t_eq[0]
                dt2 = time_series[name][file[:14] + '.000002.txt'].t_eq[1] - \
                    time_series[name][file[:14] + '.000002.txt'].t_eq[0]             
                u_error.append([time_series[name][file].x, time_series[name][file].z, 
                                np.mean(time_series[name][file].u)-np.mean(time_series[name][file[:14] + '.000002.txt'].u),
                                wt.calc_lux_data_wght(
                                    time_series[name][file].t_transit,
                                    dt1,
                                    (time_series[name][file].u_eq.dropna().values
                                    *time_series[name][file].wtref)) - 
                                wt.calc_lux_data_wght(
                                    time_series[name][file[:14] + '.000002.txt'].t_transit,
                                    dt2,
                                    (time_series[name][file[:14] + '.000002.txt'].u_eq.dropna().values
                                    *time_series[name][file[:14] + '.000002.txt'].wtref)),
                                diff_I_u, diff_I_v, diff_flux])
            if file[15:] == '000003.txt':
                I_u1, I_v1, flux1 = wt.calc_turb_data_wght(
                                        time_series[name][file].t_transit,
                                        time_series[name][file].u.dropna(),
                                        time_series[name][file].v.dropna())
                I_u2, I_v2, flux2 = wt.calc_turb_data_wght(
                                        time_series[name][file + '.000004.txt'].t_transit,
                                        time_series[name][file[:14] + '.000004.txt'].u.dropna(),
                                        time_series[name][file[:14] + '.000004.txt'].v.dropna())
                diff_I_u = abs(I_u1-I_u2)
                diff_I_v = abs(I_v1-I_v2)
                diff_flux = abs(flux1-flux2)                
                dt1 = time_series[name][file].t_eq[1] - \
                    time_series[name][file].t_eq[0]
                dt2 = time_series[name][file[:14] + '.000004.txt'].t_eq[1] - \
                    time_series[name][file[:14] + '.000004.txt'].t_eq[0]
                u_error.append([time_series[name][file].x, time_series[name][file].z, 
                                np.mean(time_series[name][file].u)-np.mean(time_series[name][file[:14] + '.000004.txt'].u),
                                wt.calc_lux_data_wght(
                                    time_series[name][file].t_transit,
                                    dt1,
                                    (time_series[name][file].u_eq.dropna().values
                                    *time_series[name][file].wtref)) - 
                                wt.calc_lux_data_wght(
                                    time_series[name][file + '.000004.txt'].t_transit,
                                    dt2,
                                    (time_series[name][file[:14] + '.000004.txt'].u_eq.dropna().values
                                    *time_series[name][file[:14] + '.000004.txt'].wtref)),
                                diff_I_u, diff_I_v, diff_flux])                            
            if file[15:] == '000005.txt':
                I_u1, I_v1, flux1 = wt.calc_turb_data_wght(
                                        time_series[name][file].t_transit,
                                        time_series[name][file].u.dropna(),
                                        time_series[name][file].v.dropna())
                I_u2, I_v2, flux2 = wt.calc_turb_data_wght(
                                        time_series[name][file[:14] + '.000006.txt'].t_transit,
                                        time_series[name][file[:14] + '.000006.txt'].u.dropna(),
                                        time_series[name][file[:14] + '.000006.txt'].v.dropna())
                diff_I_u = abs(I_u1-I_u2)
                diff_I_v = abs(I_v1-I_v2)
                diff_flux = abs(flux1-flux2)
                dt1 = time_series[name][file].t_eq[1] - \
                    time_series[name][file].t_eq[0]
                dt2 = time_series[name][file[:14] + '.000006.txt'].t_eq[1] - \
                    time_series[name][file[:14] + '.000006.txt'].t_eq[0]
                u_error.append([time_series[name][file].x, time_series[name][file].z, 
                                np.mean(time_series[name][file].u)-np.mean(time_series[name][file[:14] + '.000006.txt'].u),
                                wt.calc_lux_data_wght(
                                    time_series[name][file].t_transit,
                                    dt1,
                                    (time_series[name][file].u_eq.dropna().values
                                    *time_series[name][file].wtref)) - 
                                wt.calc_lux_data_wght(
                                    time_series[name][file[:14] + '.000006.txt'].t_transit,
                                    dt2,
                                    (time_series[name][file[:14] + '.000006.txt'].u_eq.dropna().values
                                    *time_series[name][file[:14] + '.000006.txt'].wtref)),
                                diff_I_u, diff_I_v, diff_flux])
    for lists in u_error:
        u_err_val.append(lists[2])
        x_val.append(lists[0])
        z_val.append(lists[1])
        lux_unc.append(lists[3])
        I_u_unc.append(lists[4])
        I_v_unc.append(lists[5])
        flux_unc.append(lists[6])
    u_err_std = np.std(u_err_val)
    print('     means = {}      lux-error = {}     x = {}       z = {}'.format(u_err_val, lux_unc, x_val, z_val))
    abs_u_err_val = [abs(number) for number in u_err_val]
    print('\n     mean deviation u = {}'.format(np.mean(abs_u_err_val)))
    print('     standard deviation u = {}'.format(u_err_std))
    print('     standard deviation abs_u = {}'.format(np.std(abs_u_err_val)))    
    abs_lux_unc = [abs(number) for number in lux_unc]
    print('\n     mean deviation lux = {}'.format(np.mean(abs_lux_unc)))
    print('     standard deviation lux = {}'.format(np.std(lux_unc)))
    print('     standard deviation abs_lux = {}'.format(np.std(abs_lux_unc)))

    print('\n     mean deviation I_u = {}'.format(np.mean(I_u_unc)))
    print('     standard deviation I_u = {}'.format(np.std(I_u_unc)))

    print('\n     mean deviation I_v = {}'.format(np.mean(I_v_unc)))
    print('     standard deviation I_v = {}'.format(np.std(I_v_unc)))

    print('\n     mean deviation flux_uw = {}'.format(np.mean(flux_unc)))
    print('     standard deviation flux_uw = {}'.format(np.std(flux_unc)))

    # number of bins calculated automatically using rounded value of n=sqrt(N)
    # plt.hist(u_err_val, density=True, bins='auto')
    # plt.ylabel('Probability')
    # plt.xlabel(r'$\Delta u$')
    # plt.show()

    print(' calculated error from repeated measurements \n')

# %%

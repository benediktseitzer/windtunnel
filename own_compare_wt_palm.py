# -*- coding: utf-8 -*-

################
'''
IMPORTS
'''
################

from matplotlib.pyplot import axes
import numpy as np
from numpy.lib.index_tricks import nd_grid
import pandas as pd
from scipy import stats
import logging
import matplotlib.pyplot as plt
import warnings
import os
import sys

sys.path.append('/home/bene/palm/palm_python/')
import windtunnel as wt
import palm_py as papy

# supress SOURCE ID warnings by matplotlib
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

# reference comparison
# namelist = ['BA_BL_UW_001', 'BA_BL_UW_010', 
#             'BA_S5_L_UW_002', 'BA_S5_S_UW_006', 'BA_S5_D_UW_007', 
#             'BA_S6_L_UW_007', 'BA_S6_S_UW_008', 'BA_S6_D_UW_007',
#             'BA_S7_L_UW_007', 'BA_S7_S_UW_006', 'BA_S7_D_UW_007',
#             'BA_S8_L_UW_007', 'BA_S8_S_UW_006', 'BA_S8_D_UW_007']

namelist = ['BA_BL_UW_001',
            'BA_S5_L_UW_002']

x_val_shift = 100.

# palm_python parameters
papy.globals.run_name = 'BA_BL_UW_001'
papy.globals.run_number = '.005'
papy.globals.run_numbers = ['.015', '.019']
# PHYSICS       
papy.globals.z0 = 0.021
papy.globals.alpha = 0.17
papy.globals.ka = 0.41
papy.globals.d0 = 0.

# set paths and files
experiment = 'balcony'
if len(namelist) == 1:
    config = namelist[0][3:5]
else:
    if namelist[1][-4:] == 'long':
        config = 'LONG'    
    else:
        config = namelist[0][3:5]
config = 'CO_REF'
path = '../../experiments/balcony/{}/coincidence/timeseries/'.format(config) # path to timeseries folder
wtref_path = '../../experiments/balcony/{}/wtref/'.format(config)
mean_path = '../../experiments/balcony/{}/coincidence/mean/'.format(config)
txt_path = './postprocessed/'
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
    flux_err = 0.001
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
# 3 = compare spectra
# 5 = compare profiles
# 6 = compare multiple palm to wind tunnel 
# 7 = longitudinal profile
mode = 3
calc_palm = True
outdata_path = '../wt_outdata/'# format in npz

# scale of the model.
scale = 100.

# to plot or not to plot:
plot = True

# Check if all necessary output directories exist
wt.check_directory(plot_path_0)
wt.check_directory(txt_path)

# allocate variables
time_series = {}
time_series.fromkeys(namelist)
time_series_eq = {}
time_series_eq.fromkeys(namelist)
files = []

# Gather all files into Timeseries objects, save raw timeseries
# prepare environment
for name in namelist:
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

        if data_nd == 0:
           print('Warning: Assuming that data is dimensional. If using non-dimensional input data, set variable data_nd to 1')
           ts.nondimensionalise()
        else:
           if data_nd == 1:
              []
           else:
              print('Warning: data_nd can only be 1 (for non-dimensional input data) or 0 (for dimensional input data)')        
          
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
   
# Check if positions in all files match for vertical profile
for name in namelist:
    files = wt.get_files(path, name)
    for i in range(np.size(files)-2):
        if time_series[name][files[i]].x != time_series[name][files[i+1]].x:
            raise Exception('Positions do not match! Check data file.')
        if time_series[name][files[i]].y != time_series[name][files[i+1]].y:
            raise Exception('Positions do not match! Check data file.')

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


# calculate palm measures
height_list = []
if calc_palm: 
    nc_file = '{}_pr{}.nc'.format(papy.globals.run_name, papy.globals.run_number)    
    nc_file_grid = '{}_pr{}.nc'.format(papy.globals.run_name, papy.globals.run_number)
    nc_file_path = '../../../../palm/palm/current_version/JOBS/{}/OUTPUT/'.format(papy.globals.run_name)
    mask_name_list = ['M01', 'M02', 'M03', 'M04', 'M05', 'M06', 'M07', 'M08', 'M09', 
                        'M10','M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20']
    height_list = [2., 4., 5., 7.5, 10., 15.,  20., 25., 30., 35., 40., 45., 50., 60.,
                        70., 80., 90., 100., 125., 150.]
    time_prof, time_prof_unit = papy.read_nc_time(nc_file_path, nc_file)
    
    # estimate palm_wtref
    wt_file_pr = '../../experiments/balcony/BL/coincidence/mean/BA_BL_UW_001.000001.txt'
    wt_file_ref = '../../experiments/balcony/BL/wtref/BA_BL_UW_001_wtref.txt'  
    wt_scale = 100.
    wt_u_pr, wt_u_ref, wt_z = papy.read_wt_ver_pr(wt_file_pr, wt_file_ref, wt_scale)
    print('\n wind tunnel profile loaded \n')
    z = np.linspace(0.,256.,65)
    u_pr, u_fric = papy.calc_input_profile(wt_u_pr, wt_z, z)
    palm_wtref = u_pr[-1]

    # calculate palm-velocity-profile
    grid_name = 'zu'
    z_u, z_unit = papy.read_nc_grid(nc_file_path, nc_file, grid_name)
    palm_u, palm_u_max, palm_u_unit = papy.read_nc_var_ver_pr(nc_file_path, nc_file, 'u')
    if data_nd == 0:
        palm_u = palm_u/palm_wtref
    palm_wtref = palm_u_max
    print('     palm_u_max = {}'.format(palm_u_max))
    print('     Calculated PLAM-u-profile')

    # calculate palm-fluxes
    grid_name = 'zw"u"'
    z_flux, z_unit = papy.read_nc_grid(nc_file_path, nc_file, grid_name)
    var1, var_max1, var_unit1 = papy.read_nc_var_ver_pr(nc_file_path, nc_file, 'w*u*')
    var2, var_max2, var_unit2 = papy.read_nc_var_ver_pr(nc_file_path, nc_file, 'w"u"')
    if data_nd == 0:
        palm_flux = var1/palm_wtref**2. + var2/palm_wtref**2.
        palm_flux1 = var1/palm_wtref**2.
        palm_flux2 = var2/palm_wtref**2.
    elif data_nd == 1:
        palm_flux = var1 + var2
        palm_flux1 = var1
        palm_flux2 = var2
    print('     Calculated PLAM-fluxes')

    # calculate palm-Lux
    palm_lux = np.zeros(len(height_list))
    var_name = 'u'
    for i,mask_name in enumerate(mask_name_list): 
        nc_file = '{}_masked_{}{}.nc'.format(papy.globals.run_name,mask_name,papy.globals.run_number)
        height = height_list[i]
        time, time_unit = papy.read_nc_var_ms(nc_file_path,nc_file,'time')        
        var, var_unit = papy.read_nc_var_ms(nc_file_path,nc_file,var_name)        
        palm_lux[i] = papy.calc_lux(np.abs(time[1]-time[0]),var)
        print('\n calculated integral length scale for {}'.format(str(height)))
    
    # calculate PALM turbulence intensities
    palm_Iu = np.zeros(len(height_list))
    palm_Iv = np.zeros(len(height_list))
    palm_Iw = np.zeros(len(height_list))
    u_variance_old = np.zeros(len(height_list))

    for i,mask_name in enumerate(mask_name_list): 
        nc_file = '{}_masked_{}{}.nc'.format(papy.globals.run_name, mask_name, papy.globals.run_number)
        height = height_list[i]
        
        var_u, var_unit_u = papy.read_nc_var_ms(nc_file_path, nc_file, 'u')
        var_v, var_unit_v = papy.read_nc_var_ms(nc_file_path, nc_file, 'v')
        var_w, var_unit_w = papy.read_nc_var_ms(nc_file_path, nc_file, 'w')
        u_variance_old[i] = np.std(var_u)
        turbint_dat = papy.calc_turbint(var_u, var_v, var_w)

        palm_Iu[i] = turbint_dat[0]
        palm_Iv[i] = turbint_dat[1]
        palm_Iw[i] = turbint_dat[2]

    # other try
    grid_name = 'zu'
    nc_file = '{}_pr{}.nc'.format(papy.globals.run_name, papy.globals.run_number)
    z, z_unit = papy.read_nc_grid(nc_file_path,nc_file_grid,grid_name)
    
    var_u, var_max_u, var_unit_u = papy.read_nc_var_ver_pr(nc_file_path, nc_file, 'u*2')
    
    fig, ax = plt.subplots()

    ax.semilogy(u_variance_old, height_list, label='np.std()')
    ax.semilogy(var_u[-1], z, label='pr: u*2')
    ax.grid()
    ax.legend()
    plt.show()
    plt.close()

    print('\n calculated turbulence intensities scale for {}'.format(str(height)))    

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

    plot_path = '{}{}/'.format(plot_path_0, name)

    n_outliers_u=0
    n_outliers_v=0

    # Iniate second layer of dictionaries for results 
    files = wt.get_files(path,name)
    wind_comps[name] = {}
    wind_comps[name].fromkeys(files)
    wind_stats[name] = {}
    wind_stats[name].fromkeys(files)
    wind_stats_nw[name] = {}
    wind_stats_nw[name].fromkeys(files)
    turb_data[name] = {}
    turb_data[name].fromkeys(files)
    lux_data[name] = {}
    lux_data[name].fromkeys(files)

    for file in files:
        print('     Started data-processing for: {} \n'.format(file))        
    
        n_outliers_u=n_outliers_u+time_series[name][file].n_outliers_u
        n_outliers_v=n_outliers_v+time_series[name][file].n_outliers_v
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
    print('\n Finished data-processing for: {} \n'.format(name))


# comparing spectra for single palm simulations and wind tunnel measurements
if mode == 3:
    c_list = ['orangered', 'forestgreen', 'gold', 'dodgerblue', 'darkorange', 'slategray']
    var_name_list = ['u', 'w']

    spectra_data_wt = {}
    spectra_data_wt.fromkeys(namelist)

    spectra_data_palm = {}
    spectra_data_palm.fromkeys(mask_name_list)

    for var_name in var_name_list:

        # calculation of spectra for wind tunnel data
        for j,name in enumerate(namelist):
            print('     Started spectra-processing for wind tunnel: {} \n'.format(name))            
            files = wt.get_files(path,name)         
            spectra_data_wt[name] = {}
            spectra_data_wt[name].fromkeys(files)

            for file in files:
                if var_name == 'u':
                    spectra_data_wt[name][file] = papy.calc_spectra(
                                                        time_series_eq[name][file].u_eq.dropna(), 
                                                        time_series_eq[name][file].t_eq[~np.isnan(time_series_eq[name][file].t_eq)],
                                                        time_series_eq[name][file].z,
                                                        palm_wtref)
                elif var_name == 'w':
                    spectra_data_wt[name][file] = papy.calc_spectra(
                                                        time_series_eq[name][file].v_eq.dropna(), 
                                                        time_series_eq[name][file].t_eq[~np.isnan(time_series_eq[name][file].t_eq)],
                                                        time_series_eq[name][file].z,
                                                        palm_wtref)

            print('     Finished spectra-processing for wind tunnel: {} \n'.format(name))

        # calculation of spectra for PALM-data
        print('     Compute PALM-spectra for masked output: \n')
        grid_name = 'zu'
        z, z_unit = papy.read_nc_grid(nc_file_path,nc_file_grid,grid_name)
        for i,mask_name in enumerate(mask_name_list):
            height = height_list[i]
            nc_file = '{}_masked_{}{}.nc'.format(papy.globals.run_name, mask_name, papy.globals.run_number)
            try:
                time, time_unit = papy.read_nc_var_ms(
                                            nc_file_path, 
                                            nc_file, 
                                            'time')   
            except: 
                print('\n Mask {} not in dataset. \n Check {} and the corresponding heights in the *_p3d-file'.format(mask_name, nc_file_path))
            dt, dt_unit = papy.read_nc_var_ts(
                            nc_file_path, 
                            '{}_ts{}.nc'.format(papy.globals.run_name, papy.globals.run_number), 
                            'dt')
            var, var_unit = papy.read_nc_var_ms(
                                            nc_file_path, 
                                            nc_file, 
                                            var_name)
            spectra_data_palm[mask_name] = papy.calc_spectra(
                                            var,
                                            time,
                                            height,
                                            palm_wtref)
        print('    calculated PALM-spectra for {}'.format(var_name))
        
        # plot all spectra cummulative
        for i,mask_name in enumerate(mask_name_list):
            plot_height = False
            plt.style.use('classic')
            fig, ax = plt.subplots()
            height_c = height_list[i]

            f_sm = spectra_data_palm[mask_name][0]
            S_uu_sm = spectra_data_palm[mask_name][1]
            comp1_aliasing = spectra_data_palm[mask_name][2]
            f_sm = [f_sm][np.argmin([np.nanmax(f_sm)])]
            f_sm = f_sm[:len(S_uu_sm)]
            gridsize = 2.
            filter_width_scaled = gridsize/(palm_wtref*np.mean(dt))
            h1 = ax.loglog(f_sm[:comp1_aliasing], S_uu_sm[:comp1_aliasing], 
                            marker='o', markersize=3, color='darkviolet',
                            label=r'PALM at $z={}$'.format(height_c))
            h2 = ax.loglog(f_sm[comp1_aliasing-1:], S_uu_sm[comp1_aliasing-1:], 
                            marker='o', markersize=3, color='violet',
                            fillstyle='none')
            hf = ax.axvline(x=filter_width_scaled, ls= '--', color='gray')

            l = 0
            for j,name in enumerate(namelist):
                files = wt.get_files(path,name)
                for file in files:
                    if height_c == time_series_eq[name][file].z:
                        l += 1
                        file_c = file
                        height = height_c
                        f_sm_wt = spectra_data_wt[name][file][0]
                        S_wt_sm = spectra_data_wt[name][file][1]
                        wt_aliasing = spectra_data_wt[name][file][2]
                        f_sm_wt = [f_sm_wt][np.argmin([np.nanmax(f_sm_wt)])]
                        f_sm_wt = f_sm_wt[:len(S_wt_sm)]
                        if name == 'BA_BL_UW_001':
                            h3 = ax.loglog(f_sm_wt[:wt_aliasing+1], S_wt_sm[:wt_aliasing+1],
                                            marker='x', markersize=3, color=c_list[j],
                                            label=r'boundary layer at $z={}$m'.format( 
                                            time_series_eq[name][file].z))
                        elif name == 'BA_BL_UW_010':
                            h3 = ax.loglog(f_sm_wt[:wt_aliasing+1], S_wt_sm[:wt_aliasing+1],
                                            marker='x', markersize=3, color=c_list[j],
                                            label=r'smooth wall at $z={}$m'.format( 
                                            time_series_eq[name][file].z))
                        else:
                            h3 = ax.loglog(f_sm_wt[:wt_aliasing+1], S_wt_sm[:wt_aliasing+1],
                                            marker='x', markersize=3, color=c_list[j],
                                            label=r'$s_b={}$m at $z={}$m'.format(name[4], 
                                            time_series_eq[name][file].z))
            try:
                f_refspecs = np.logspace(-4, 3, num=100, base = 10) 
                ref_specs = papy.get_reference_spectra(height_c,
                                '../../../../palm/palm_python/reference_data/')
                E_min, E_max = papy.calc_ref_spectra(f_refspecs, 
                                ref_specs, var_name)
                if var_name == 'u':
                    h5 = ax.fill_between(f_refspecs, E_min, E_max,
                                    facecolor=(1.,0.6,0.6),edgecolor='none',alpha=0.2,
                                    label=r'VDI-range $S _{uu}$')
                elif var_name == 'w':
                    h5 = ax.fill_between(f_refspecs, E_min, E_max,
                                    facecolor=(1.,0.6,0.6),edgecolor='none',alpha=0.2,
                                    label=r'VDI-range $S _{ww}$')
            except:
                print('\n There are no reference-spectra available for this flow \n')

            ax.set_xlim([10**-4,300])
            ax.set_ylim([10 ** -4, 1])
            ax.set_xlabel(r"$f\cdot z\cdot u_{ref}^{-1}$")
            ax.set_ylabel(r"$f\cdot S_{ij}\cdot (\sigma_i\sigma_j)^{-1}$")
            ax.legend(loc='lower left')
            ax.grid()
            if l == len(namelist):                       
                plot_height = True
            if plot_height:
                plt.savefig(plot_path + 'spectra_'+ var_name + '_' + mask_name + '.' + file_type, 
                            bbox_inches='tight')
            else:
                plt.close()


# comparing mode for single palm simulations and wind tunnel measurements
if mode == 5:
    plt.style.use('classic')
    plt.figure(7)
    # c_list = ['#00006f', '#000092', '#0000b4',
    #             '#003eda', '#007cff', '#3ebdff',
    #             '#7cffff', '#a8ebf4', '#d5d7e9']
    # c_list = ['#00006f', '#a8ebf4']    
    c_list = ['orangered', 'forestgreen', 'gold', 'dodgerblue', 'darkorange', 'slategray']
    component = 'w'

    y_max = 300.
    y_min = 1.

    #########
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
                                    color=c_list[j], label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name[6] == 'S':                
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='d',
                                    color=c_list[j], label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name[6] == 'D':              
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='+',
                                    color=c_list[j], label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name == 'BA_BL_UW_010':
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='x',
                                    color=c_list[j], label=r'smooth wall'.format(x_val))
                elif name == 'BA_BL_UW_001':
                    U = ax.errorbar(u_mean[i], heights[i], xerr=xerror, marker='^',
                                    color=c_list[j], label=r'boundary layer'.format(-130))
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
        j = j + 1

    for i in range(len(time_prof)-1,len(time_prof)):
        try:
            ax.plot(palm_u[i,:-1], z_u[:-1], label='PALM', 
                    color='darkviolet')
        except:
            print('Exception has occurred: StopIteration - plot_ver_profile')

    ax.set_xlabel(r'$u/u_{ref}$ (-)')
    ax.set_ylabel(r'$z$ (m)')
    ax.set_ylim(y_min, y_max) 
    if data_nd == 0:
        ax.set_xlim(0.4,1.)   
    elif data_nd == 1:
        axes.set_ylim(min(palm_u), max(palm_u))   
    plt.legend(loc= 'upper left', numpoints=1)
    plt.yscale('log')    
    plt.tight_layout()
    plt.grid(True,'both','both')
    plt.savefig(plot_path_0 + 'cummulative_wind_data_log_' + name[:-4] + '.' + file_type,
                bbox_inches='tight')
    plt.close(7)

    #######
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
                                label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name[6] == 'S':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='d',color=c_list[j],
                                label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name[6] == 'D':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='+',color=c_list[j],
                                label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='x',color=c_list[j],
                                label=r'smooth wall'.format(x_val))
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='^',color=c_list[j],
                                label=r'boundary layer')
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

    for i in range(len(time_prof)-1,len(time_prof)):
        try:
            ax.plot(palm_flux[i,:-1], z_flux[:-1], label='PALM', color = 'darkviolet')
            # ax.plot(palm_flux1[i,:-1], z_flux[:-1], label='$>\Delta$', color = 'plum')
            # ax.plot(palm_flux2[i,:-1], z_flux[:-1], label='$<\Delta$', color = 'magenta')
        except:
            print('Exception has occurred: StopIteration - plot_ver_profile')

    ax.set_xlabel(r'u' + '\'' + component + '\' $\cdot$ $u_{ref}^{-2}\ (-)$')
    ax.set_ylabel(r'$z$ (m)')
    ax.set_ylim(y_min, y_max)
    if data_nd == 0:
        ax.set_xlim(-0.004,0.)
    elif data_nd == 1:
        ax.set_xlim(-0.1,0.)
    plt.legend(loc= 'upper left', numpoints=1)
    plt.yscale('log')
    plt.tight_layout()
    plt.grid(True,'both','both')    
    plt.savefig(plot_path_0 + 'cummulative_flux_log_' + name[:-4] + '.' + file_type,
                bbox_inches='tight')
    plt.close(8)

    ######
    # I_u
    plt.figure(9)
    j = 0
    for name in namelist:
        heights = []
        I_u_list = []
        files = wt.get_files(path,name)
        turb_data[name] = {}
        turb_data[name].fromkeys(files)

        #wind tunnel data
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
                                    label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name[6] == 'S':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='d',color=c_list[j],
                                    label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name[6] == 'D':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='+',color=c_list[j],
                                    label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='x',color=c_list[j],
                                    label=r'smooth wall'.format(x_val))
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(I_u_list[i],heights[i],xerr=I_u_err,fmt='^',color=c_list[j],
                                    label=r'boundary layer')
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

    pl = ax.errorbar(palm_Iu, height_list, xerr=0.02*palm_Iu, markersize=3,
                fmt='o', color = 'darkviolet', label='PALM')
    ax.set_xlabel(r'$I_u$ (-)')
    ax.set_ylabel(r'$z$ (m)')
    ax.set_ylim(y_min, y_max)
    plt.legend(numpoints=1)
    plt.yscale('log')
    plt.tight_layout()
    plt.grid(True,'both','both')    
    plt.savefig(plot_path_0 + 'cummulative_I_u_log_' + name[:-4] + '.' + file_type,
                bbox_inches='tight')
    plt.close(9)

    ######
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
                                    label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name[6] == 'S':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='d',color=c_list[j],
                                    label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name[6] == 'D':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='+',color=c_list[j],
                                    label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='x',color=c_list[j],
                                    label=r'smooth wall'.format(x_val))
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(I_v_list[i],heights[i],xerr=I_v_err,fmt='^',color=c_list[j],
                                    label=r'boundary layer')
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

    pl = ax.errorbar(palm_Iw, height_list, xerr=0.02*palm_Iv, markersize=3,
                fmt='o', color = 'darkviolet', label='PALM')

    ax.set_xlabel(r'$I_w$ (-)')
    ax.set_ylabel(r'$z$ (m)')
    ax.set_ylim(y_min, y_max)
    plt.legend(numpoints=1)    
    plt.yscale('log')
    plt.tight_layout()
    plt.grid(True,'both','both')    
    plt.savefig(plot_path_0 + 'cummulative_I_w_log_' + name[:-4] + '.' + file_type,
                bbox_inches='tight')
    plt.close(10)

    ######
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
                                    label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name[6] == 'S':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='d',color=c_list[j],
                                    label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name[6] == 'D':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='+',color=c_list[j],
                                    label=r'$s_b={}$ m'.format(name[4], x_val))
                elif name == 'BA_BL_UW_010':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='x',color=c_list[j],
                                    label=r'smooth wall'.format(x_val))
                elif name == 'BA_BL_UW_001':
                    l = ax.errorbar(lux_list[i],heights[i],xerr=lux_err,fmt='^',color=c_list[j],
                                    label=r'boundary layer')
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
    h1 = ax.errorbar(palm_lux, height_list, xerr=lux_err, 
                fmt='o', markersize=3, color='darkviolet', label=r'PALM')
                
    Lux_10,Lux_1,Lux_01,Lux_001,Lux_obs_smooth,Lux_obs_rough = \
    wt.get_lux_referencedata(ref_path)
    ref1 = ax.plot(Lux_10[1,:],Lux_10[0,:],'k-',linewidth=1)
    ref2 = ax.plot(Lux_1[1,:],Lux_1[0,:],'k--',linewidth=1)
    ref3 = ax.plot(Lux_01[1,:],Lux_01[0,:],'k-.',linewidth=1)
    ref4 = ax.plot(Lux_001[1,:],Lux_001[0,:],'k:',linewidth=1)
    ref5 = ax.plot(Lux_obs_smooth[1,:],Lux_obs_smooth[0,:],'k+',
                    linewidth=1)
    ref6 = ax.plot(Lux_obs_rough[1,:],Lux_obs_rough[0,:],'kx',linewidth=1)

    ax.set_xlabel(r'$L_u^x$ (m)')
    ax.set_ylabel(r'$z$ (m)')
    ax.set_ylim(y_min, y_max)
    ax.set_xlim(10., 1000.)    
    plt.legend(loc= 'upper left', numpoints=1)
    plt.yscale('log')
    plt.xscale('log')
    plt.tight_layout()
    plt.grid(True,'both','both')    
    plt.savefig(plot_path_0 + 'cummulative_Lux_log_' + name[:-4] + '.' + file_type,
                bbox_inches='tight')
    print(plot_path_0)
    plt.close(11)

    print('\n plotted cummulative profiles over length.\n')


# comparing mode for multiple palm simulations and wind tunnel measurements
if mode == 6:
    plt.style.use('classic')
    c_list = ['orangered', 'forestgreen', 'gold', 'dodgerblue', 'darkorange', 'slategray']
    # c_list = ['orangered', 'darkorange', 'peru', 'chocolate']
    palm_data = {}
    palm_data.fromkeys(papy.globals.run_numbers)
    var_name_list = ['u', 'flux']
    #read palm-data and init 
    for run in papy.globals.run_numbers:
        print('     Start processing palm-run #{}'.format(run[-3:]))
        papy.globals.run_number = run
        palm_data[papy.globals.run_number] = {}
        palm_data[papy.globals.run_number].fromkeys(var_name_list)
        nc_file = '{}_pr{}.nc'.format(papy.globals.run_name,papy.globals.run_number)        
        
        # read variables for plot
        time, time_unit = papy.read_nc_time(nc_file_path,nc_file)
        # read wind tunnel profile
        wt_pr, wt_u_ref, wt_z = papy.read_wt_ver_pr(wt_file_pr, wt_file_ref , wt_scale)        
        
        for i,var_name in enumerate(var_name_list):
            if var_name == 'u':
                grid_name = 'z{}'.format(var_name)        
                var, var_max, var_unit = papy.read_nc_var_ver_pr(nc_file_path,nc_file,var_name)
                z, z_unit = papy.read_nc_grid(nc_file_path,nc_file,grid_name)
                palm_data[run][var_name] = var
                palm_wtref = var_max
                print(palm_wtref)
            elif var_name == 'flux':
                grid_name = 'zw*u*'
                var1, var_max1, var_unit1 = papy.read_nc_var_ver_pr(nc_file_path, nc_file, 'w*u*')
                var2, var_max2, var_unit2 = papy.read_nc_var_ver_pr(nc_file_path, nc_file, 'w"u"')
                if data_nd == 0:
                    var = var1/palm_wtref**2. + var2/palm_wtref**2.
                    var1 = var1/palm_wtref**2.
                    var2 = var2/palm_wtref**2.
                elif data_nd == 1:
                    var = var1 + var2
                    var1 = var1
                    var2 = var2
                palm_data[run][var_name] = var
            else:
                grid_name = 'z{}'.format(var_name)        
                var, var_max, var_unit = papy.read_nc_var_ver_pr(nc_file_path,nc_file,var_name)
                z, z_unit = papy.read_nc_grid(nc_file_path,nc_file,grid_name)
                palm_data[run][var_name] = var                
        print('     End processing palm-run #{} \n'.format(run[-3:]))

    print(' Start plotting of palm-runs of {}'.format(papy.globals.run_name))
    # plot flux data
    plot_flux = True
    if plot_flux:
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
            ax.grid(True, 'both', 'both')

            for i in range(np.size(fluxes)):
                if i == 1:
                    # l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='o',color=c_list[j],
                    #                 label = r'fluxes at $x={}$ m'.format(str(x_val)))
                    if name[6] == 'L':
                        l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='o',color=c_list[j],
                                    label=r'$s_b={} m$'.format(name[4]))
                    elif name[6] == 'S':
                        l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='d',color=c_list[j],
                                    label=r'$s_b={} m$'.format(name[4]))
                    elif name[6] == 'D':
                        l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='x',color=c_list[j],
                                    label=r'$s_b={} m$'.format(name[4]))
                    elif name == 'BA_BL_UW_010':
                        l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='+',color=c_list[j],
                                    label=r'smooth wall'.format(x_val))
                    elif name == 'BA_BL_UW_001':
                        l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='^',color=c_list[j],
                                    label=r'boundary layer')
                else:
                    if name[6] == 'L':
                        l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='o',color=c_list[j])
                    if name[6] == 'S':
                        l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='d',color=c_list[j])
                    if name[6] == 'D':
                        l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='x',color=c_list[j])
                    elif name == 'BA_BL_UW_010':
                        l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='+',color=c_list[j])
                    elif name == 'BA_BL_UW_001':
                        l = ax.errorbar(fluxes[i],heights[i],xerr=flux_err,fmt='^',color=c_list[j])
            j += 1

        color_list = ['orchid', 'hotpink']
        z0_list = [0.02, 0.06]

        for i in range(len(time)-1,len(time)):
            ax.fill_betweenx(z[:-1], palm_data[papy.globals.run_numbers[0]]['flux'][i,:-1], 
                    palm_data[papy.globals.run_numbers[1]]['flux'][i,:-1], color ='thistle',
                    label = 'PALM: $0.02 m<z_0<0.06 m$')
        for i in range(len(time_prof)-1,len(time_prof)):
            try:
                ax.plot(palm_flux[i,:-1], z_flux[:-1], label='PALM: $z_0=0.021$ $m$', color = 'darkviolet')
            except:
                print('Exception has occurred: StopIteration - plot_ver_profile')

        ax.set_xlabel(r'$u$' + '\'' + '$w$' + '\' $\cdot$ $u_{ref}^{-2}\ (-)$')
        ax.set_ylabel(r'$z$ $(m)$')
        ax.set_ylim(1., 265)
        if data_nd == 0:
            ax.set_xlim(-0.0045,0.)
        elif data_nd == 1:
            ax.set_xlim(-0.1,0.)
        plt.legend(loc= 'upper left', numpoints=1)
        # plt.yscale('log')
        ax.set_yscale('log', nonposy='clip')
        plt.tight_layout()
        plt.grid(True,'both','both')    
        plt.savefig(plot_path_0 + 'cummulative_flux_log_range' + '.' + file_type,
                    bbox_inches='tight')

    print(' End plotting of palm-runs of {} \n'.format(papy.globals.run_name))

# %%

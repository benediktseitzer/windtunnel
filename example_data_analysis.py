# -*- coding: utf-8 -*-

import numpy as np
import logging
import windtunnel as wt
import matplotlib.pyplot as plt


# Create logger
logger = logging.getLogger()




#%%#
# This is an example script. It can be modified to your personal needs.
# Specify path to data, path to wtref, output paths for plots and txt file, 
# file type for plots, name of files, scale and desired mode of analysis.
# Input paths for data and wtref with a list of names of the measurement files
path = '/path/to/your/data/' # path to timeseries folder
wtref_path = '/path/to/your/wtref/'

#namelist = ['HG_BL_UW_0' + str(x).zfill(2) for x in range(12,2,-1)]
namelist = ['name_of measurement_file']

txt_path = './postprocessed/'
#edit 06/20/19: set ref_path to none for unknown reference path
#ref_path = './plots/'
ref_path = None
file_type = 'png'
scale = 100
plot_path = '/path/to/your/plots'

#edit 08/08/2019: add errors for all quantities
u_err=0
v_err=0
#u_v_error is for plotting wind data, where u and v are plotted on one plot. 
#It is assumed in this case the error to be plotted is the u error. 
#TODO: modify  code to allow for seperate errors for u and v on the same plot. 
u_v_err=u_err
I_u_err=0
I_v_err=0
flux_err=0
lux_err=0



# 1 = vertical profile
# 2 = lateral profile
# 3 = convergence test
# 4 = Reynolds Number Independence
mode = 2
if mode == 2:
    outdata_path = '/path/to/your/outdata_lat/'# format in npz
else:
    outdata_path = '/path/to/your/outdata/'# format in npz

#plot scatter
plot = True
scatter = True
save_data = True
# Check if all necessary output directories exist
wt.check_directory(plot_path)
wt.check_directory(txt_path)

time_series = {}
time_series.fromkeys(namelist)

#edit 06/20/19: create seperate time series for equidistant data
time_series_eq = {}
time_series_eq.fromkeys(namelist)

#set data_nd to 1 if using non-dimensional data
data_nd=0
# Gather all files into Timeseries objects, save raw timeseries

for name in namelist:
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
        ts_eq=ts
        ts_eq.calc_equidistant_timesteps()  
        ts.index=ts.t_arr         
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




if files==[]:
   raise Exception('No Matching File Names Found. Please check namelist and/or path!') 
                   
for name in namelist:
   # Check if positions in all files match for vertical profile
   files = wt.get_files(path, name)
   if not mode == 2:
       for i in range(np.size(files)-2):
           if time_series[name][files[i]].x != time_series[name][files[i+1]].x:
               raise Exception('Positions do not match! Check data file.')
           if time_series[name][files[i]].y != time_series[name][files[i+1]].y:
               raise Exception('Positions do not match! Check data file.')
   # Check if positions in all files match for horizontal profile
   if mode == 2:
       for i in range(np.size(files)-2):
           if time_series[name][files[i]].x != time_series[name][files[i+1]].x:
               raise Exception('Positions do not match! Check data file.')
           if time_series[name][files[i]].z != time_series[name][files[i+1]].z:
               raise Exception('Positions do not match! Check data file.')


# Iniate first layer of dictionaries for results
wind_comps = {}
wind_comps.fromkeys(namelist)
wind_stats = {}
wind_stats.fromkeys(namelist)
turb_data = {}
turb_data.fromkeys(namelist)
lux_data = {}
lux_data.fromkeys(namelist)
spectra_data = {}
spectra_data.fromkeys(namelist)

for name in namelist:
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
        turb_data[name] = {}
        turb_data[name].fromkeys(files)
        lux_data[name] = {}
        lux_data[name].fromkeys(files)
        spectra_data[name] = {}
        spectra_data[name].fromkeys(files)
    for file in files:
        #edit 6/12/19: add variables n_outliers_u and n_outliers_v to keep track of number of outliers        
        n_outliers_u=n_outliers_u+time_series[name][file].n_outliers_u
        n_outliers_v=n_outliers_v+time_series[name][file].n_outliers_v
        #edit 6/12/19: add u and v outliers to determien total 
        wind_comps[name][file] = time_series[name][file].wind_comp1,\
                                 time_series[name][file].wind_comp2
       
        if mode == 3:
            # Perform convergence test and plot results, no
            # output saved as txt, as programme ends at "break"
            # Average u and v component for different (default) intervals
            # Set maximum interval size, interval and blocksize and initialise
            # list of qunatities.
            max_interval = int(np.size(time_series[name][file].u))
            interval =5000
            blocksize = 10000
            # time step
            dt = time_series[name][file].t_eq[1] - \
                 time_series[name][file].t_eq[0]
            # reference length for dimensionless time
            ref_length = 1
            # nondimensionalise time step using wtref and ref_length
            # dt = dt*time_series[name][file].wtref/ref_length
            # determine averaging intervals
            intervals = np.arange(interval, int(0.5 * max_interval), blocksize)
            # intervals = intervals*dt
            quantities = ['Magnitude', 'u_mean',
                          wind_comps[name][file][1] + '_mean', 'u_std',
                          wind_comps[name][file][1] + '_std', 'I_u',
                          'I_' + wind_comps[name][file][1], 'flux', 'Lux']

            # Initialise dictionary using list of quantities
            conv_data = {}
            conv_data.fromkeys(quantities)
            for quantity in quantities:
                conv_data[quantity] = {}
                conv_data[quantity].fromkeys(intervals)

            # Calculate convergence test of each quantity for intervals up to
            # the maximum interval size. Each iteration passes a larger array
            # interval to the wt.calc_* function (this is inefficient regarding
            # computing time but avoids very short arrays being passed to the
            # wt.calc_* functions). The data is saved in the conv_data
            # dictionary.
            while interval < max_interval / 2:
                M_list = []
                u_mean_list = []
                v_mean_list = []
                u_std_list = []
                v_std_list = []
                I_u_list = []
                I_v_list = []
                flux_list = []
                Lux_list = []
                for i in range(0, max_interval - interval, interval):
                    #                    dt = time_series[name][file].t_eq[i+interval] -\
                    #                         time_series[name][file].t_eq[i]
                    M, u_mean, v_mean, u_std, v_std, dd = wt.calc_wind_stats(
                        time_series[name][file].u.iloc[i:i + interval],
                        time_series[name][file].v.iloc[i:i + interval])
                    M_list.append(M)
                    u_mean_list.append(u_mean)
                    v_mean_list.append(v_mean)
                    u_std_list.append(u_std)
                    v_std_list.append(v_std)

                    I_u, I_v, flux = wt.calc_turb_data(
                        time_series[name][file].u.iloc[i:i + interval],
                        time_series[name][file].v.iloc[i:i + interval])
                    I_u_list.append(I_u)
                    I_v_list.append(I_v)
                    flux_list.append(flux)

                    Lux = wt.calc_lux_data(dt,
                                           time_series[name][file].u.dropna().values[i:i + interval])
                    Lux_list.append(Lux)

                conv_data['Magnitude'][interval] = np.asarray(M_list)
                conv_data['u_mean'][interval] = np.asarray(u_mean_list)
                conv_data[wind_comps[name][file][1] + '_mean'][interval] = \
                    np.asarray(v_mean_list)
                conv_data['u_std'][interval] = np.asarray(u_std_list)
                conv_data[wind_comps[name][file][1] + '_std'][interval] = \
                    np.asarray(v_std_list)
                conv_data['I_u'][interval] = np.asarray(I_u_list)
                conv_data['I_' + wind_comps[name][file][1]][interval] = \
                    np.asarray(I_v_list)
                conv_data['flux'][interval] = np.asarray(flux_list)
                conv_data['Lux'][interval] = np.asarray(Lux_list)

                interval = interval + blocksize

            # Plot convergence test results for each of the nine quanities
            # investigated. The plot is saved in plot_path, specified at the
            # beginning of this example programme.
            plt.figure(1001)
            wt.plots.plot_convergence(conv_data,
                                      wtref=time_series[name][file].wtref,
                                      ref_length=ref_length, scale=scale)
            plt.tight_layout()
            plt.savefig(plot_path + 'convergence_' + name + '.' + file_type,
                        dpi=1000,bbox_inches='tight')

            raise Exception('Convergence test complete!')
    
        # Calculate mean wind quantities
        dt = time_series[name][file].t_eq[1] - time_series[name][file].t_eq[0]
        wind_stats[name][file] = wt.calc_wind_stats_wght(time_series[name][file].t_transit,time_series[name][file].u,
                                                    time_series[name][file].v)
        turb_data[name][file] = wt.calc_turb_data(time_series[name][file].u.dropna(),
                                                  time_series[name][file].v.dropna())
        #edit 06/20/2019: changed script to adapt to dimensional and non-dimensional input data
        if data_nd==0:
           lux_data[name][file] = wt.calc_lux_data(dt,
                                                (time_series[name][file].u_eq.dropna().values*
                                                time_series[name][file].wtref))
        if data_nd==1:
           lux_data[name][file] = wt.calc_lux_data(dt,
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
                                  spectra_data[name][file][4],
                                  spectra_data[name][file][5],
                                  wind_comps[name][file],
                                  time_series[name][file].z,
                                  ref_path=ref_path)
            
            plt.savefig(plot_path + 'spectra_' + file[:-4] + '.' + file_type)
            plt.close('all')

            print('\n Start wavelet analysis for {}'.format(file))
            wavelet, scale = wt.calc_wavelet_transform(time_series[name][file].u_eq,
                                                    time_series[name][file].t_eq,
                                                    wavelet='morlet')
            # y_val = time_series[name][file].y-y_val_shift
            y_val = time_series[name][file].z            
            f_scale = y_val/(scale * np.mean(time_series[name][file].u_eq))

            # plot wavelet transform
            plt.figure(55)
            plt.style.use('classic')    
            wt.plots.plot_wavelet_transform(wavelet, 
                                    scale, 
                                    time_series[name][file].u_eq, 
                                    time_series[name][file].t_eq,
                                    y_val)

            plt.savefig(plot_path +  'wavelet_analysis_' + file +  '.' + file_type,
                        bbox_inches='tight', dpi=300)
            print(' Saved plot to:'+ plot_path + 'wavelet_analysis_' + file +  '.' + file_type)
            plt.close(55)
            print(' Finished wavelet analysis for {}'.format(file))



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
    
    for file in files:
        # Gather all quantities for a complete profile
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

    if mode == 4:
        # Perform and plot results of Reynolds Number Independence test, no
        # output saved as txt, as programme ends at "break"
        plt.figure(0)
        wt.plots.plot_Re_independence(mean_mag, wtref, yerr=u_err,ymin=0,ymax=1)
        #wt.plots.plot_Re_independence(u_std, wtref, yerr=u_err,ymin=0,ymax=0.3)
        plt.tight_layout()
        plt.savefig(plot_path + 'Re_u' + name + '.' + file_type)
        break

   # Save quantities for vertical and lateral profiles    
    #np.savetxt(txt_path + name + '_turb.txt',
               #np.vstack((x, y, heights, mean_mag, u_mean, u_mean_wght, v_mean,
                          #v_mean_wght, u_std, u_std_wght, v_std, v_std_wght,
                          #I_u, I_v, lux, fluxes, wdir, wtref)).transpose(),
               #fmt='%.8f', header=('flow and turbulence parameters\n'
                                   #'units: dimensionless!\n'
                                   #'format: standard numpy.genfromtxt()\n'
                                   #'variables = \"x\" \"y\" \"z\" \"M\" '
                                   #'\"{0}_mean\" \"{0}_mean_wght\" '
                                   #'\"{1}_mean\" \"{1}_mean_wght\" \"{0}_std\"'
                                   #' \"{0}_std_wght\" \"{1}_std\" '
                                   #'\"{1}_std_wght\" \"I_{0}\" \"I_{1}\" '
                                   #'\"L{0}x\" \"{0}\'{1}\'_flux\" \"wdir\" '
                                   #'\"wtref\"'.format(wind_comps[name][file][0],
                                                      #wind_comps[name][file][1])))  
    if mode == 1 and plot:
        # Plot results of a vertical profile
        # Wind components
        plt.figure(0)
        ret, lgd = wt.plots.plot_winddata(mean_mag,u_mean,v_mean,heights,yerr=u_v_err)
        plt.savefig(plot_path + 'wind_data_' + name + '.' + file_type,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
    
        # Wind components, logarithmic y-axis
        plt.figure(1)
        ret, lgd = wt.plots.plot_winddata_log(mean_mag,u_mean,v_mean,heights,xerr=u_v_err)
        plt.tight_layout()
        plt.savefig(plot_path + 'wind_data_log_' + name + '.' + file_type,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
        
        # Turbulence intensity of the first component
        plt.figure(2)
        wt.plots.plot_turb_int(I_u,heights,yerr=I_u_err,component='I_'+ts.wind_comp1,ref_path=ref_path)
        plt.tight_layout()
        plt.savefig(plot_path + 'I_'+ts.wind_comp1+'_' + name + '.' + file_type)
    
        # Turbulence intensity of the second component
        plt.figure(3)
        wt.plots.plot_turb_int(I_v,heights,yerr=I_v_err,component='I_'+ts.wind_comp2,ref_path=ref_path)
        plt.tight_layout()
        plt.savefig(plot_path + 'I_'+ts.wind_comp2+'_' + name + '.' + file_type)
    
        # Profile of the fluxes
        plt.figure(4) 
        wt.plots.plot_fluxes(fluxes,heights,yerr=flux_err,component='w')
        plt.tight_layout()
        plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)
    
        # Profiles of the fluxes, logarithmic y-axis
        plt.figure(5)
        wt.plots.plot_fluxes_log(fluxes,heights,yerr=flux_err,component='w')
        plt.tight_layout()
        plt.savefig(plot_path + 'fluxes_log_' + name + '.' + file_type)
    
        # Double logarithmic profile of Lux data
        plt.figure(6)
        wt.plots.plot_lux(lux,heights,yerr=lux_err,component='w',ref_path=ref_path)
        plt.tight_layout()
        plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)
        plt.close('all')
    if mode == 2:
        # Results of a lateral profile
        # Wind components
        plt.figure(0)
        ret, lgd = wt.plots.plot_winddata(mean_mag,u_mean,v_mean,y,yerr=u_v_err,lat=True)
        plt.tight_layout()
        plt.savefig(plot_path + 'wind_data_' + name + '.' + file_type,
                    bbox_extra_artists=(lgd,), bbox_inches='tight')
        
        # Turbulence intensity of the first component
        plt.figure(1)
        wt.plots.plot_turb_int(I_u,y,yerr=I_u_err,lat=True)
        plt.tight_layout()
        plt.savefig(plot_path + 'I_u_' + name + '.' + file_type)
    
        # Turbulence intensity of the second component
        plt.figure(2)
        wt.plots.plot_turb_int(I_v,y,yerr=I_v_err,component='I_w',lat=True)
        plt.tight_layout()
        plt.savefig(plot_path + 'I_v_' + name + '.' + file_type)
    
        # Profile of the fluxes
        plt.figure(3) 
        wt.plots.plot_fluxes(fluxes,y,yerr=flux_err,component='w',lat=True)
        plt.tight_layout()
        plt.savefig(plot_path + 'fluxes_' + name + '.' + file_type)
    
        # Lateral profile of Lux data
        plt.figure(4)
        wt.plots.plot_lux(lux,y,yerr=lux_err,component='w',lat=True)
        plt.tight_layout()
        plt.savefig(plot_path + 'Lux_' + name + '.' + file_type)
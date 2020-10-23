#! /usr/bin/python3
# -*- coding: utf-8 -*-
""" Statistical and calculation tools for boundary layer analysis.
"""

import numpy as np
import scipy.stats as sc
import math as m
from scipy.optimize import curve_fit

import windtunnel as wt

__all__ = [
    'calc_intervalmean',
    'calc_stats',
    'calc_exceedance_prob',
    'calc_wind_stats',
    'calc_wind_stats_wght',
    'calc_turb_data',
    'calc_turb_data_wght',
    'calc_lux_data',
    'calc_lux_data_wght',
    'calc_acorr',
    'calc_autocorr',
    'calc_spectra',
    'calc_ref_spectra',
    'convergence_test_1',
    'convergence_test_2',
    'convergence_test',
    'power_law',
    'calc_alpha',
    'calc_z0',
    'calc_alpha_profile',
    'calc_normalization_params',
    'calc_theo_arrival_law',
    'calc_arrival_law',
    'calc_transit_time_distribution'
]

def calc_intervalmean(indata,intervals,DD=False):    
    """ Calculates interval means of indata. If DD is set to True the means are 
    calculated for circular quantities. Returns a dictionary with 
    intervals as keys. If intervals has length 1 the function returns an array.
    @parameter: indata, type = any
    @parameter: intervals, type = list
    @parameter: DD, type = boolean"""
   
    outdata = {}
    outdata.fromkeys(intervals)
    
    if len(intervals) == 1:
        outdata = np.array([])
   
    for interval in intervals:
        avg_data = np.zeros(int(np.ceil(np.size(indata) / interval)))
        for n, i in enumerate(range(0, int(0.5*np.size(indata)), interval)):
            if DD:
                u_east = np.nanmean(np.sin(indata[i:i+interval] * np.pi/180))
                u_north = np.nanmean(np.cos(indata[i:i+interval] * np.pi/180))
                unit_WD = np.arctan2(u_east, u_north) * 180/np.pi
                avg_data[n] = (360 + unit_WD) % 360
            else:
                avg_data[n] = np.nanmean(indata[i:i+interval])            
        if len(intervals) > 1:    
            outdata[interval] = avg_data[:-1]
        else:
            outdata = avg_data[:-1]
            
    return outdata

def calc_stats(sets,DD=False):
    """Returns mean, standard deviation and variance of data in sets. If DD is 
    true then the circular equivalents are calculated. TO BE USED WITH CAUTION
    @parameter sets: iterable set of data
    @parameter DD: boolean"""
       
    means = np.array([])
    var = np.array([])
    stds = np.array([])
    for data in sets:
        mask = ~np.isnan(data)
        data = data[mask]
        if DD:
            u_east = np.nanmean(np.sin(data * np.pi/180))
            u_north = np.nanmean(np.cos(data * np.pi/180))
            tmp = np.arctan2(u_east, u_north) * 180/np.pi
            means = np.append(means,m.degrees(tmp)%360)
            var = np.append(var,sc.circvar(data))
            stds = np.append(stds,sc.circstd(data))
        else:
            means = np.append(means,np.mean(data))
            var = np.append(var,np.var(data))
            stds = np.append(stds,np.std(data))
        
    return means,var,stds

def calc_exceedance_prob(data,threshold):
    """ Calculates exceedance probability of threshold in data. Returns 
    threshold and exceedance probability in percent.
    @parameter data: 
    @parameter threshold: int """
    
    tmp = data[data>threshold]
    exceed_prob = (np.size(tmp)/np.size(data))*100
    
    return threshold, exceed_prob

def calc_wind_stats(u_comp,v_comp,wdir=0.):
    """ Calculate wind data from equidistant times series of u and 
    v components. wdir is a reference wind direction.
    @parameter: u_comp: np.array or list
    @parameter: v_comp: np.array or list
    @parameter: wdir: int"""
    
    mask = mask = np.logical_and(~np.isnan(u_comp),
                          ~np.isnan(v_comp))
    u = u_comp[mask]
    v = v_comp[mask]
    
    Magnitude = np.mean(np.sqrt(u**2+v**2))
    u_mean = np.mean(u)
    v_mean = np.mean(v)       
    Direction = wdir-np.arctan2(v_mean,u_mean)*180/np.pi
    if Direction>360: Direction-=360
    if Direction<0: Direction+=360    
    u_std = np.std(u)
    v_std = np.std(v)

    data = np.array([Magnitude,u_mean,v_mean,u_std,v_std,Direction])

    return data
       
def calc_wind_stats_wght(transit_time,u_comp,v_comp,wdir=0.):
    """ Calculate wind data from equidistant times series of u and 
    v components. wdir is a reference wind direction.
    @parameter: transit_time, type = np.array
    @parameter: u_comp, type = np.array
    @parameter: v_comp, type = np.array
    @parameter: wdir: int"""
    
    mask = mask = np.logical_and(~np.isnan(u_comp),
                          ~np.isnan(v_comp),~np.isnan(transit_time))
    u = u_comp[mask]
    v = v_comp[mask]
    tt = transit_time[mask]

    # TODO: test ways of applying TT weighting to Magnitude
    Magnitude = np.mean(np.sqrt(u**2+v**2))
    u_mean = wt.transit_time_weighted_mean(tt,u)
    v_mean = wt.transit_time_weighted_mean(tt,v)       
    Direction = wdir-np.arctan2(v_mean,u_mean)*180/np.pi
    if Direction>360: Direction-=360
    if Direction<0: Direction+=360    
    u_std = np.sqrt(wt.transit_time_weighted_var(tt,u))
    v_std = np.sqrt(wt.transit_time_weighted_var(tt,v))

    data = np.array([Magnitude,u_mean,v_mean,u_std,v_std,Direction])

    return data

def calc_turb_data(u_comp,v_comp):
    """ Calculate turbulence intensity and turbulent fluxes from equidistant
    times series of u and v components.
    @parameter: u_comp: np.array or list
    @parameter: v_comp: np.array or list""" 
    
    mask = mask = np.logical_and(~np.isnan(u_comp),
                          ~np.isnan(v_comp))    
    u = np.asarray(u_comp[mask])
    v = np.asarray(v_comp[mask])
    
    M = np.mean(np.sqrt(u_comp**2 +v_comp**2))
    u_mean = np.mean(u)
    v_mean = np.mean(v)
    u_dev = u - u_mean
    v_dev = v - v_mean
    u_std = np.std(u_comp)
    v_std = np.std(v_comp)
    ##  TURBULENCE INTENSITY
    I_u = u_std/np.mean(M)
    I_v = v_std/np.mean(M)
    ##  TURBULENT FLUXES
    flux = np.mean((u_dev*v_dev).round(7)).round(6)
    
    data = np.array([I_u,I_v,flux])
    
    return data

def calc_turb_data_wght(transit_time,u_comp,v_comp):
    """ Calculate turbulence intensity and turbulent fluxes from equidistant
    times series of u and v components using transit time weighted statistics.
    @parameter: transit_time. type = np.array
    @parameter: u_compy type = np.array
    @parameter: v_comp, type = np.array"""    
    mask = mask = np.logical_and(~np.isnan(u_comp),
                          ~np.isnan(v_comp))    
    u = u_comp[mask]
    v = v_comp[mask]
    
    # TODO: test ways of applying TT weighting to M
    M = np.mean(np.sqrt(u_comp**2 +v_comp**2))
    u_std = np.sqrt(wt.transit_time_weighted_var(transit_time,u))
    v_std = np.sqrt(wt.transit_time_weighted_var(transit_time,v))
    # TURBULENCE INTENSITY
    I_u = u_std/np.mean(M)
    I_v = v_std/np.mean(M)
    # Fluxes
    flux = wt.transit_time_weighted_flux(transit_time,u_comp,v_comp)
    
    data = np.array([I_u,I_v,flux])
    
    return data

def calc_lux_data(dt, u_comp):
    """ Calculates the integral length scale according to R. Fischer (2011) 
    from an equidistant time series of the u component using time step dt.
    @parameter: t_eq, type = int or float
    @parameter: u_comp, type = np.array or list """

    if np.size(u_comp) < 5:
        raise Exception('Too few value to estimate Lux!')

    mask = np.where(~np.isnan(u_comp))

    u = u_comp[mask]

    initial_guess = 1000 / dt  # number of values for the length of the first autocorrelation
    if initial_guess >= len(u):
        initial_guess = int(len(u) / 2)
    lag_eq = np.arange(1, np.size(u) + 1) * dt  # array of time lags

    for lag in range(int(initial_guess), len(lag_eq), int(initial_guess)):
        u_eq_acorr = calc_acorr(u, lag)  # autocorrelation (one sided) of time series u_eq
        if np.any(np.diff(u_eq_acorr) > 0):  # if a single value of the derivation autocorrelation is smaller than 0
            # the iteration of the autocorrelation stops
            break

    lag_eq = lag_eq[:len(u_eq_acorr)]

    Lux = 0.
    # see dissertation R.Fischer (2011) for calculation method
    for i in range(np.size(u_eq_acorr) - 2):

        autc1 = u_eq_acorr[i]

        autc2 = u_eq_acorr[i + 1]

        Lux = Lux + (autc1 + autc2) * 0.5

        if autc2 > autc1:
            acorr_fit = np.polyfit(lag_eq[:i], np.log(abs(u_eq_acorr[:i])), deg=1)
            acorr_fit = np.exp(acorr_fit[0] * lag_eq + acorr_fit[1])

            if np.min(acorr_fit) < 0.001:
                ix = np.where(acorr_fit < 0.001)[0][0]

            else:
                ix = acorr_fit.size

            Lux = Lux + (np.sum(acorr_fit[i + 1:ix]) +
                         np.sum(acorr_fit[i + 2:ix + 1])) * 0.5
            break

        elif autc1 <= 0:
            break

    Lux = abs(Lux * np.mean(u_comp) * dt)
    return Lux

def calc_lux_data_wght(transit_time, dt, u_comp):
    """ Calculates the integral length scale according to R. Fischer (2011)
    from an equidistant time series of the u component using time step dt.
    @parameter: t_eq, type = int or float
    @parameter: u_comp, type = np.array or list """

    if np.size(u_comp) < 5:
        raise Exception('Too few value to estimate Lux!')

    mask = np.where(~np.isnan(u_comp))

    u = u_comp[mask]

    initial_guess = 1000 / dt  # number of values for the length of the first autocorrelation
    if initial_guess >= len(u):
        initial_guess = int(len(u) / 2)
    lag_eq = np.arange(1, np.size(u) + 1) * dt  # array of time lags

    for lag in range(int(initial_guess), len(lag_eq), int(initial_guess)):
        u_eq_acorr = calc_acorr(u, lag)  # autocorrelation (one sided) of time series u_eq
        if np.any(np.diff(u_eq_acorr) > 0):  # if a single value of the derivation autocorrelation is smaller than 0
            # the iteration of the autocorrelation stops
            break

    lag_eq = lag_eq[:len(u_eq_acorr)]

    Lux = 0.
    # see dissertation R.Fischer (2011) for calculation method
    for i in range(np.size(u_eq_acorr) - 2):

        autc1 = u_eq_acorr[i]

        autc2 = u_eq_acorr[i + 1]

        Lux = Lux + (autc1 + autc2) * 0.5

        if autc2 > autc1:
            acorr_fit = np.polyfit(lag_eq[:i], np.log(abs(u_eq_acorr[:i])), deg=1)
            acorr_fit = np.exp(acorr_fit[0] * lag_eq + acorr_fit[1])

            if np.min(acorr_fit) < 0.001:
                ix = np.where(acorr_fit < 0.001)[0][0]

            else:
                ix = acorr_fit.size

            Lux = Lux + (np.sum(acorr_fit[i + 1:ix]) +
                         np.sum(acorr_fit[i + 2:ix + 1])) * 0.5
            break

        elif autc1 <= 0:
            break

    Lux = abs(Lux * wt.transit_time_weighted_mean(transit_time, u_comp) * dt)

    return Lux

def calc_acorr(timeseries, maxlags):
    """ Full autocorrelation of time series for lags up to maxlags.
    @parameter timeseries: np.array or list
    @parameter maxlags: int"""

    timeseries = timeseries[~np.isnan(timeseries)]
    acorr = np.asarray([1. if x == 0 else np.corrcoef(timeseries[x:], timeseries[:-x])[0][1] for x in range(maxlags)])
    return acorr

def calc_autocorr(timeseries, lag=1):
    """ Autocorrelation of time series with lag.
    @parameter tiemseries: np.array or list
    @parameter lag: int"""

    timeseries = timeseries[~np.isnan(timeseries)]
    autocorr = np.corrcoef(timeseries[0:np.size(timeseries) - lag],
                           timeseries[lag:])[1, 0]

    return autocorr

def calc_spectra(u_comp,v_comp,t_eq,height):
    """ Calculate dimensionless energy density spectra from an equidistant 
    time series.
    @parameter: u_comp, type = np.array or list
    @parameter: v_comp, type = np.array or list
    @parameter: t_eq, type = np.array or list """
    ## FREQUENCY
    freq = np.fft.fftfreq(np.size(u_comp),t_eq[1]-t_eq[0])

    ## FFT
    fft_u = np.fft.fft(u_comp)*1./np.size(u_comp)        #  normalized fft
    fft_v = np.fft.fft(v_comp)*1./np.size(v_comp) 
    
    ## DISCRETE SPECTRA
    nyquist_freq = np.size(t_eq)//2   #  nyquist frequency index

    #  S == E/dn ; dn == 1
    if np.size(u_comp)%2==0:
        E_uu = np.hstack((np.abs(fft_u[0:nyquist_freq])**2 *
                          2.,np.abs(fft_u[nyquist_freq])**2))
        
        # frequency transformation (n to f as referenced in Stull, 
        # p.314) P = N * dt; n=f/P
        S_uu = E_uu * len(t_eq)*(t_eq[1]-t_eq[0])
        E_vv = np.hstack((np.abs(fft_v[0:nyquist_freq])**2 *
                          2.,np.abs(fft_v[nyquist_freq])**2))
        
        # frequency transformation (n to f as referenced in Stull, p.314) 
        # P = N * dt; n=f/P
        S_vv = E_vv * len(t_eq)*(t_eq[1]-t_eq[0])
        E_uv = np.hstack((np.abs(fft_u[0:nyquist_freq])*
                          np.abs(fft_v[0:nyquist_freq]) * 2.,
                          np.abs(fft_u[nyquist_freq])*
                          np.abs(fft_v[nyquist_freq])))
        
        # frequency transformation (n to f as referenced in Stull, p.314) 
        # P = N * dt; n=f/P
        S_uv = E_uv * len(t_eq)*(t_eq[1]-t_eq[0])
    else:
        E_uu = np.abs(fft_u[0:nyquist_freq+1])**2 * 2.
        
        # frequency transformation (n to f as referenced in Stull, p.314) 
        # P = N * dt; n=f/P
        S_uu = E_uu * len(t_eq)*(t_eq[1]-t_eq[0])
        E_vv = np.abs(fft_v[0:nyquist_freq+1])**2 * 2.
        
        # frequency transformation (n to f as referenced in Stull, p.314) 
        # P = N * dt; n=f/P
        S_vv = E_vv * len(t_eq)*(t_eq[1]-t_eq[0])
        E_uv = np.abs(fft_u[0:nyquist_freq+1])*np.abs(fft_v[0:nyquist_freq+1])*2.
        
        # frequency transformation (n to f as referenced in Stull, p.314)
        # P = N * dt; n=f/P
        S_uv = E_uv * len(t_eq)*(t_eq[1]-t_eq[0])

    # dimensionless
    S_uu = np.abs(freq[0:nyquist_freq+1])*S_uu/(np.nanstd(u_comp)**2)
    S_vv = np.abs(freq[0:nyquist_freq+1])*S_vv/(np.nanstd(v_comp)**2)
    S_uv = np.abs(freq[0:nyquist_freq+1])*S_uv/np.abs(u_comp.std()*v_comp.std())

    ##  REDUCED FREQUENCY (PLOT and reference spectra)
    reduced_freq = freq*height/np.mean(u_comp)

    ##  SMOOTH SPECTRA
    # use exponents for equi-distant bins in log plot
    f_sm = np.hstack((np.array([0]),
                      np.log10(np.abs(reduced_freq[1:nyquist_freq]))))
    valcount, edges = np.histogram(f_sm,bins=np.arange(
                                   f_sm.min(),f_sm.max()+10**(-5),0.05))
    f_sm = np.zeros(valcount.shape)     
    S_uu_sm = np.zeros(valcount.shape)
    S_vv_sm = np.zeros(valcount.shape)
    S_uv_sm = np.zeros(valcount.shape)

    vc = 0  # counting values
    for i,n in enumerate(valcount):
        if n>0:
            f_sm[i] = np.mean(np.abs(reduced_freq)[vc:vc+n])
            S_uu_sm[i] = np.mean(S_uu[vc:vc+n])
            S_vv_sm[i] = np.mean(S_vv[vc:vc+n])
            S_uv_sm[i] = np.mean(S_uv[vc:vc+n])
            vc=vc+n   
    
    ##  SORTING
    f_arg = np.argsort(f_sm)
    f_sm = f_sm[f_arg]
    S_uu_sm = S_uu_sm[f_arg]
    S_vv_sm = S_vv_sm[f_arg]
    S_uv_sm = S_uv_sm[f_arg]
    
    ##  ALIASING
    u_aliasing = f_sm.size-9+np.hstack((np.where(
                                    	 np.diff(S_uu_sm[-10:])>=0.)[0],[9]))[0]
    
    v_aliasing = f_sm.size-9+np.hstack((np.where(
                                   np.diff(S_vv_sm[-10:])>=0.)[0],[9]))[0]
    
    uv_aliasing = f_sm.size-9+np.hstack((np.where(
                                         np.diff(S_uv_sm[-10:])>=0.)[0],[9]))[0]

    return f_sm,S_uu_sm,S_vv_sm,S_uv_sm,u_aliasing,v_aliasing,uv_aliasing

def calc_spectra_nc(u_comp, t_eq, height):
    """ Calculate dimensionless energy density spectra from an equidistant
    time series.
    @parameter: u_comp, type = np.array or list
    @parameter: t_eq, type = np.array or list """
    ## FREQUENCY
    freq = np.fft.fftfreq(np.size(u_comp), t_eq[1] - t_eq[0])

    ## FFT
    fft_u = np.fft.fft(u_comp) * 1. / np.size(u_comp)  # normalized fft

    ## DISCRETE SPECTRA
    nyquist_freq = np.size(t_eq) // 2  # nyquist frequency index

    #  S == E/dn ; dn == 1
    if np.size(u_comp) % 2 == 0:
        E_uu = np.hstack((np.abs(fft_u[0:nyquist_freq]) ** 2 *
                          2., np.abs(fft_u[nyquist_freq]) ** 2))

        # frequency transformation (n to f as referenced in Stull,
        # p.314) P = N * dt; n=f/P
        S_uu = E_uu * len(t_eq) * (t_eq[1] - t_eq[0])
    else:
        E_uu = np.abs(fft_u[0:nyquist_freq + 1]) ** 2 * 2.

        # frequency transformation (n to f as referenced in Stull, p.314)
        # P = N * dt; n=f/P
        S_uu = E_uu * len(t_eq) * (t_eq[1] - t_eq[0])

    # dimensionless
    S_uu = np.abs(freq[0:nyquist_freq + 1]) * S_uu / (np.nanstd(u_comp) ** 2)

    ##  REDUCED FREQUENCY (PLOT and reference spectra)
    reduced_freq = freq * height / np.mean(u_comp)

    ##  SMOOTH SPECTRA
    # use exponents for equi-distant bins in log plot
    f_sm = np.hstack((np.array([0]),
                      np.log10(np.abs(reduced_freq[1:nyquist_freq]))))
    valcount, edges = np.histogram(f_sm, bins=np.arange(
        f_sm.min(), f_sm.max() + 10 ** (-5), 0.05))
    f_sm = np.zeros(valcount.shape)
    S_uu_sm = np.zeros(valcount.shape)

    vc = 0  # counting values
    for i, n in enumerate(valcount):
        if n > 0:
            f_sm[i] = np.mean(np.abs(reduced_freq)[vc:vc + n])
            S_uu_sm[i] = np.mean(S_uu[vc:vc + n])
            vc = vc + n

    ##  SORTING
    f_arg = np.argsort(f_sm)
    f_sm = f_sm[f_arg]
    S_uu_sm = S_uu_sm[f_arg]

    ##  ALIASING
    u_aliasing = f_sm.size - 9 + np.hstack((np.where(
        np.diff(S_uu_sm[-10:]) >= 0.)[0], [9]))[0]

    return f_sm, S_uu_sm, u_aliasing


#    """ Calculate dimensionless energy density spectra from an equidistant 
#    time series.
#    @parameter: u_comp, type = np.array or list
#    @parameter: v_comp, type = np.array or list
#    @parameter: t_eq, type = np.array or list """
#    ## FREQUENCY
#    freq = np.fft.fftfreq(np.size(u_comp),t_eq[1]-t_eq[0])
#
#    ## FFT
#    fft_u = np.fft.fft(u_comp)*1./np.size(u_comp)        #  normalized fft
#    fft_v = np.fft.fft(v_comp)*1./np.size(v_comp) 
#    uv_param = np.sqrt(fft_u * fft_v) # This is used in place of a true Fourier 
#                                      # transform to calculate diagonal energy
#    
#    u_normalization_params = wt.calc_normalization_params(freq, fft_u, t_eq, 
#                                                       height, 
#                                                       np.nanmean(u_comp), 
#                                                       np.std(u_comp), 
#                                                       len(u_comp))
#    v_normalization_params = wt.calc_normalization_params(freq, fft_v, t_eq, 
#                                                       height, 
#                                                       np.nanmean(v_comp), 
#                                                       np.std(v_comp), 
#                                                       len(v_comp))
#    uv_normalization_params = wt.calc_normalization_params(freq, uv_param, 
#                                                       t_eq, 
#                                                       height, 
#                                                       np.nanmean(
#                                                       np.sqrt(v_comp*u_comp)), 
#                                                       np.std(
#                                                       np.sqrt(v_comp*u_comp)),
#                                                       len(v_comp))
#    
#    return u_normalization_params[0],u_normalization_params[1],\
#           v_normalization_params[1],uv_normalization_params[1],\
#           u_normalization_params[2],v_normalization_params[2],\
#           uv_normalization_params[2]


def calc_ref_spectra(reduced_freq,a,b,c,d,e):
   """ Calculate dimensionless reference spectra. ???
   @parameter: reduced_freq, type = ???
   @parameter: a, type = ???
   @parameter: b, type = ???
   @parameter: c, type = ???
   @parameter: d, type = ???
   @parameter: e, type = ???"""
   # TODO: figure what's going on here
   e=e+0j

   return a*reduced_freq/np.abs(e+b*reduced_freq**c)**d

def convergence_test_1(data,blocksize=100):
    """ Conducts a block-wise convergence test on non circular data using 
    blocksize for the size of each increment. Returns a dictionary block_data.
    Each entry is named after its respective interval. blocksize's default 
    value is 100.
    @parameter: data, type = np.array or list
    @parameter: blocksize, type = int or float"""
    
    if blocksize > 0.5*np.size(data):
        raise Exception('blocksize must be smaller than half of the length\
        of data in order to maintain independent values.')
        
    intervals = list(np.arange(1,int(0.5*np.size(data))-1,blocksize))
    block_data = wt.calc_intervalmean(data,intervals)
    
    return intervals, block_data

def convergence_test_2(data,interval=100,blocksize=100):
    """ Conducts a block-wise convergence test on non circular data using 
    blocksize for the size of each increment between intervals. Returns a 
    dictionary block_data. Each entry is named after its respective interval.
    blocksize's and interval's default values are 100.
    @parameter: data, type = np.array or list
    @parameter: interval, type = int
    @parameter: blocksize, type = int"""
    
    if blocksize > 0.5*np.size(data):
        raise Exception('blocksize must be smaller than half of the length\
        of data in order to maintain independent values.')

    max_interval = int(np.size(data))
    
    intervals = np.arange(interval,int(0.5*max_interval),blocksize)
    block_data = {}
    block_data.fromkeys(intervals)
    
    while interval < max_interval/2:
        tmp = []
        for i in range(0,max_interval-interval,interval):
            tmp.append(np.mean(data[i:i+interval]))
            block_data[interval] = np.asarray(tmp)

        interval = interval + blocksize
    
    return intervals, block_data

def convergence_test(values, min_interval = 1000, max_num_intervals = 100, calc_overlap = True,
                     max_overlap = 1, overlap_stepsize = 0.1):
    '''
    Conducts a convergence test on non circular data by sequentially splitting the data into smaller fractions. The
    convergence test begins by calculating the total mean (1/1), then calculating the next smaller fractions (1/2) and
    so on (1/n). This continues until the max_num_intervals (default 100) or a minimum interval length (default 1000) is
    reached. If an overlap between the blocks is desired the calc_overlap flag should be set to True (is default). With
    the overlap active the blocks are overlapping each other. A maximum overlap of 1 (default) doesnt result in a
    complete overlap but rather a maximum overlap of 1 - overlap_stepsize (default 0.1) = 0.9.

    Returns a dictionary mean_vals. Each entry is named after its respective blocksize.

    Parameters
    ----------
    values: ndarray
    min_interval: int
    max_num_intervals: int
    calc_overlap: bool
    max_overlap: float
    overlap_stepsize: float

    Returns
    -------
    mean_vals: dict

    '''
    vals_org = np.copy(values)
    overlaps = np.arange(0, max_overlap, overlap_stepsize)
    fracs = np.array([12, 24, 40, 60, 84])
    mean_vals = {len(vals_org): np.mean(vals_org)}
    n = 2
    length = len(vals_org)

    while (length > min_interval) and (n < max_num_intervals):
        # with an increase of n, the fractions in which vals is split up increases. For n = 2, vals is split up in
        # halfs. For n = 4 vals is split up in quarters and so on.
        #
        vals = np.copy(vals_org)
        chunks = np.array_split(vals, n)
        length = int(len(vals) / n)
        mean_vals[length] = ([np.mean(chunk) for chunk in chunks])
        # if an overlap is desired an overlap of the chunks is also used to calculate a mean. In case of halfs (0.5,
        # 50%) of values of the total number of values an initial overlap of 0.1 (10%) this results in a first chunk
        # that starts at 5% of the total number of values.
        # The last chunk (overlap 0.9, 90%) then begins at 45% of all values and ends on 95%.
        if calc_overlap:
            for overlap in overlaps:
                vals = vals_org[int((len(vals_org) / n) 
                    * overlap):-int((len(vals_org) / n) * (1 - overlap))]
                chunks = np.array_split(vals, n - 1)
                mean_vals[length].extend([np.mean(chunk) for chunk in chunks])
        # This section is a little experimental so far. It works properly but I'm not sure if it is needed...
        # Here the total number of values is not split up in chunks like quarters but into 5/12, 7/24 and so on. This
        # is implemented to fill the large gaps of the first chunk segments (halfs, thirds, quarters, etc.). These
        # uneven chunks also move over the total number of values. E.g. 5/12 starting on the first value to 5/12. Then
        # starting on 1/12 to 6/12 and so on.
        if np.any(n == fracs):
            numerator = int(n / (np.where(n == fracs)[0][0] + 2) - 1)
            bounds = []
            offset = 0
            while offset < int(numerator):
                lower_bounds = np.where((np.arange(0, n)[offset:] 
                            - offset) % numerator == 0)[0] + offset
                upper_bounds = np.where((np.arange(0, n)[offset:] 
                            - offset) % numerator == (numerator - 1))[0] + offset
                lower_bounds = lower_bounds[:len(upper_bounds)]
                bounds.extend(list(zip(lower_bounds, upper_bounds)))
                offset += 1
            vals = np.copy(vals_org)
            chunks = np.array_split(vals, n)
            length = int(len(vals) * numerator / n)
            mean_vals[length] = [np.mean([np.mean(chunk) for chunk in chunks[bound[0]:bound[1]]]) for bound in bounds]

        n += 1

    return mean_vals

def power_law(u_comp,height,u_ref,z_ref,alpha,d0=0):
    """ Estimate power law profile.
    @parameter: u_comp, type = int or float
    @parameter: height, type = int or float
    @parameter: u_ref, type = int or float
    @parameter: z_ref, type = int or float
    @parameter: alpha, type = int or float
    @parameter: d0, type = int or float """
   
    return np.abs(u_comp / u_ref - ((height-d0)/(z_ref-d0))**alpha)

def calc_alpha(u_mean,heights,d0=0.,sfc_height=120.,BL_height=600.):
    """Estimate the power law exponent alpha.
    @parameter: u_mean, type = list or np.array
    @parameter: heights, type = list or np.array
    @parameter: d0, type = float
    @parameter: sfc_height, type = float
    @parameter: BL_height, type = float"""
    #note 07/28/2020: extremely questionable assumptions for calculating
    #alpha. In particular, it is unknown why uref and zref are not constant
    #with height. Different approach used in 'calc_profile' function. To be
    #verified with Bernd Leitl and/or Frank Harms !!
    
    u_mean = np.asarray(u_mean)
    heights = np.asarray(heights)
    
    BL = np.where(heights<BL_height)
    if np.size(heights[BL]) < 6:
        print("   small sample - alpha estimation with high uncertainty")
    def tempf(x,ex):
        return x**ex
    explist=[]
    if np.size(heights[BL]) > 2:
        ref = 999.
        for ui,zi in enumerate(heights[BL]):
            zref=zi
            zref = 867.5
            uref = u_mean[BL][ui]
            uref = 1                
            B,covtemp = curve_fit(tempf, (heights[BL]-d0)/(zref-d0),
                                  u_mean[BL]/uref)        
            diff = power_law(u_mean[BL],heights[BL],uref,zref,B[0],d0)
            diff = np.sum(diff)
            explist.append(B[0])
            if diff < ref:
                if zi>800:
                    print("ATTENTION! ref height (full-scale) {0}m on the edge of the BL ".format(zi))
                alpha = B[0]
                ref = diff
        ref = (max(explist)-min(explist))/2.
    else:
         print('    too few points for alpha estimation')
         alpha = np.nan
         ref = np.nan
        
    return alpha, ref

def calc_z0(u_mean,heights,d0=0.,sfc_height=120.,BL_height=600.):
    """ Estimate the roughness length z0.
    @parameter: u_mean, type = list or np.array
    @parameter: heights, type = list or np.array
    @parameter: d0, type = float
    @parameter: sfc_height, type = float
    @parameter: BL_height, type = float """   
    u_mean = np.asarray(u_mean)
    heights = np.asarray(heights)
    
    ## z0
    sfc_layer = np.where(heights<sfc_height)
    if np.size(heights[sfc_layer]) > 2:
        z0 = np.polyfit(u_mean[sfc_layer],np.log(heights[sfc_layer]-d0),deg=1)
        err = np.mean(np.abs(np.exp(u_mean[sfc_layer]*z0[0]+z0[1]) -
                             heights[sfc_layer]+d0))
        #  if the error is too great, try deleting the worst point
        for i in heights[sfc_layer]:
            if err > 10:
                print("      z0 error - deleting one point ")
                print("      remaining points: ",np.size(heights[sfc_layer])-1)
                pt = np.argmax( np.abs(u_mean[sfc_layer]*z0[0]+z0[1] -
                                       np.log(heights[sfc_layer])) )

                sfc_layer = np.hstack((sfc_layer[:pt],sfc_layer[pt+1:]))
                if np.size(heights[sfc_layer]) > 0:
                    z0 = np.polyfit(u_mean[sfc_layer],
                                       np.log(heights[sfc_layer]-d0),deg=1)
                    err = np.mean(np.abs(np.exp(u_mean[sfc_layer]*z0[0]+z0[1])- 
                                         heights[sfc_layer]))        
            else:
                break       
        z0 = np.exp(z0[-1])
        z0_list = [z0]
        for i in range(1,np.size(heights[sfc_layer])//3):
            if i < np.size(heights[sfc_layer])//6+1:
                z00 = np.polyfit(u_mean[sfc_layer][i:],
                                 np.log(heights[sfc_layer][i:]),deg=1)
            else:
                z00 = np.polyfit(u_mean[sfc_layer][:-i+np.size(heights[sfc_layer])//6],
                                 np.log(heights[sfc_layer][:-i+np.size(heights[sfc_layer])//6]),
                                        deg=1)
            z0_list.append(np.exp(z00[-1]))
        err = (max(z0_list)-min(z0_list))/2.
    else:
        print('    too few points for z0 estimation')
        z0 = np.nan
        err = np.nan
    
    return z0, err

def calc_alpha_profile(mean_mag, heights, wtref, z_ref, d_0=0, mode='all',min_height=None,max_height=None,split_height=None):
    """Calculate profile exponent alpha and roughness length z_0, and save data to excel file. 
    Avaliable modes are 'all' (calculates profile between min_height and max_height), and 
    'seperate', which calulates two profiles, one above and one below split_height.
    Setting minimum and maximum height to 'None' calculates profile to bottom and top of
    data respectively. All heights assumed in m full-scale."""
    print("Warning: Assuming that wind data is non-dimensional, and that heights are in m full-scale. d_0 in mm.")
    #Note: heights should be in (m) full-scale. Change code if this is not the case!
    
    heights=np.asarray(heights)
    mean_mag=np.asarray(mean_mag)
    if min_height == None:
       min_height = np.min(heights)
    if max_height == None:
       max_height = np.max(heights)
    if mode == 'all':
       heights_mask=(heights>=min_height) & (heights<=max_height)
       heights=heights[heights_mask]
       mean_mag=mean_mag[heights_mask]
       z_norm = (np.array(heights)-d_0)/(z_ref-d_0)
       #Note: the wind data in this script is already non-dimensional.
       #If the code is edited to allow for dimensional plotting, change 
       #function to non-dimensionalise data if necessary!
       z_norm=np.vstack([np.log(z_norm),np.zeros(len(z_norm))]).transpose()
       print(z_norm)
       alpha=np.linalg.lstsq(z_norm,np.log(np.array(mean_mag)))
       print(alpha) 
       return alpha
    elif mode == 'seperate':
       heights_bottom=heights[(heights>=min_height) & (heights<=split_height)]
       mean_mag_bottom = mean_mag[(heights>=min_height) & (heights<=split_height)]
       z_norm_bottom = (np.array(heights_bottom)-d_0)/(z_ref-d_0)
       #Note: the wind data in this script is already non-dimensional.
       #If the code is edited to allow for dimensional plotting, change 
       #function to non-dimensionalise data if necessary!
       z_norm_bottom=np.vstack([np.log(z_norm_bottom),np.zeros(len(z_norm_bottom))]).transpose()
       print(z_norm_bottom)
       alpha_bottom=np.linalg.lstsq(z_norm_bottom,np.log(np.array(mean_mag_bottom)))
       print(alpha_bottom)
       heights_top=heights[(heights>=split_height) & (heights<=max_height)]
       mean_mag_top = mean_mag[(heights>=split_height) & (heights<=max_height)]
       z_norm_top = (np.array(heights_top)-d_0)/(z_ref-d_0)
       #Note: the wind data in this script is already non-dimensional.
       #If the code is edited to allow for dimensional plotting, change 
       #function to non-dimensionalise data if necessary!
       z_norm_top=np.vstack([np.log(z_norm_top),np.zeros(len(z_norm_top))]).transpose()
       print(z_norm_top)
       alpha_top=np.linalg.lstsq(z_norm_top,np.log(np.array(mean_mag_top)))
       print(alpha_top)  
       return alpha_top

def calc_normalization_params(freqs, transform, t, height, mean_x, sdev_x, 
                              num_data_points):
    pass
#    """ Calculate the normalized Fourier transform and frequency for the 
#    Fourier transform of x
#    Warning: A previous code version normalized segments, while this version 
#    normalizes the entire data set at once. The previous version also included 
#    a smoothing algorithm, which has been omitted for simplicity.
#    @parameter: freqs, type = list or np.array
#    @parameter: transform, type = list or np.array - this is the non-normalized
#                                                     Fourier transform
#    @parameter: t, type = float - this is time
#    @parameter: height, type = float - z in the Timeseries object
#    @parameter: mean_x, type = float - the mean of the parameter of the Fourier
#                                       transform F(x)
#    @parameter: sdev_x, type = float - the standard deviation of x
#    @parameter: num_data_points, type = int - the number of elements in x 
#                                              before the transform was found"""   
#
#    ## DISCRETE SPECTRA
#    transform = transform/num_data_points
#
#
#    E = transform ** 2
#    S = E * len(t)*(t[1]-t[0])
#   
#    ##  REDUCED FREQUENCY (PLOT and reference spectra)
#    reduced_freq = np.abs(freqs*height/mean_x)
#    reduced_transform = np.abs(np.meshgrid(S[0],
#                               freqs,sparse=True)[0]*S/sdev_x**2)
#    reduced_transform = reduced_transform[0]
#
#    ##  ALIASING
#    aliasing = reduced_freq.size - 9 + \
#               np.hstack((np.where(np.diff(
#                          reduced_transform[-10:])>=0.)[0],[9]))[0]
#    
#    return reduced_transform, reduced_freq, aliasing

def calc_theo_arrival_law(t_arr, data_rate):
    """ 
    calculate theoretical particle arrival law. 
    if exponential, there is temporally uniform seeding.
    
    @parameter: t_arr, type = list or np.array, arrival times
    @parameter: data_rate, type = float, data rate from mean file 
    """

    # allocate
    delta_t_arr = []
    # calculate inter arrival times for each burst
    delta_t_arr = [ t_arr[i+1] - t_arr[i] for i in range(len(t_arr)-1) ]
    delta_t_arr = np.asarray(delta_t_arr)
    # calculate particle arrival law: p(t_arr) = dN/dt * exp(- dN/dt * t_arr)
    particle_arrival_law = data_rate * np.exp(-data_rate * delta_t_arr)

    return delta_t_arr, particle_arrival_law

def calc_arrival_law(t_arr, data_rate):
    """ 
    calculate particle arrival law and fit the distribution. 
    if exponential, there is temporally uniform seeding.
    
    @parameter: t_arr, type = list or np.array, arrival times
    @parameter: data_rate, type = float, data rate N/T
    """

    # allocate
    delta_t_arr = []
    # calculate inter arrival times for each burst
    delta_t_arr = [ t_arr[i+1] - t_arr[i] for i in range(len(t_arr)-1) ]
    delta_t_arr = np.asarray(delta_t_arr)

    data_entries, bins = np.histogram(delta_t_arr, bins='auto',density=True)
    binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])

    def fit_function(x, A):
        return (A * np.exp(-x * A) )

    popt, pcov = curve_fit(fit_function, xdata=binscenters, ydata=data_entries)
    print('     fitted data rate = {}'.format(popt))
    print('     expected data rate = {}'.format(data_rate))
    # scale
    # data_entries = data_entries * popt

    return binscenters, data_entries, popt

def calc_transit_time_distribution(transit_time):
    """ 
    calculate particle arrival law. 
    if exponential, there is temporally uniform seeding.
    
    @parameter: transit_time, type = list or np.array, 
                                    time particle moves 
                                    through measurement volume
    """

    return sc.skew(transit_time, nan_policy='omit')


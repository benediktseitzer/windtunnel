PuffConcentration.py

07/24/2019:

in calc_net_concentration:
	changed 'net_concentration' from key to atribute,
 in alignment with the PointConcentration.py script
in apply_threshold_dosage:
	new function, nearly identical to 'apply_threshold_concentration' above, but 
here a threshold total dosage, in place of a threshold peak concentration, is used.
in clear_zeros:
	new function, nearly identical to 'clear_zeros' in PointConcentration.py, but
 the varibale 'signal' is also masked. Further, since no full scale concentration exists ts this point,
 the mask is based on the sign of net_concentration. Aditionally, masking operation is not performed on the variables, c*, full_scale concentration, and full_scale_time,
 which are (currently) not calculated in this script. 
in plot_puff: 
	new function, which plots time series of slected variables
	
07/26/2019:

in __init__:
	added attributes
	
in apply_threshold_concentration: 
	added logging of puffs below threshold concentration
in plot_puff:
	data confirms no puff before start of data aquisiton;
	display arrival and leaving time in plot;
	no puff before start of data aquisiton; 
	save plots	

07/29/2019:

in __init__:
	added attribute mask for proper plotting
in calc_net_concentration:
	changed clear_zeros to a standalone funciton, which is now called from 
outside calc_net_concentration
in get_dosage:
	added new algorithm to calcualte dosage to check the original algorithm. The questionable results of orignal algorithm noted on 07/06/2019 have been attributed to incorrect indexing 
	in plot_puff function. The algorithm for calcuaitng the dosage thus seems to be working correctly. 
in get_mask: 
	new fuction, returns the indeces of all data points which remain
 after applying threshold concentration & dosage.  	
in apply_threshold_concentration: 
	added logging of puffs below threshold concentration
in plot_puff:
	adjustments to account for masking of variables


08/01/2019:

in get_mean_puffs:
	new function, calculates mean puff, based on calculation in original C Program.
 See Bachelor Thesis of Anne Philipp (2010) for more details.
in apply_threshold_concentration: 
	added logging of puffs below threshold concentration
in plot_mean_puff:
	new function, which plots the time series of the mean puff, as calculated in get_mean_puff above, as a function of time. 

08/02/2019:

in __init__:
	added attributes signal_array, puffs_array, mean_array	, pct10_puff, pct10_signal, pct90_puff, pct90_signal
in get_mean_puff: 
	changed puffs_array to atribute, added computation of mean signal;
        compute 10th and 90th percentile of puff and signal in addition to mean
in plot_mean_puff:
	added plotting of mean signal. Also, added new input variables stats and dist. If stats is 
set to on, function also plots the mean dosage, arrival time, leaving time, and peak time.
        If dist is set to 
on, function also plots the 10th and 90th percentile of the puff and signal;
        added optional plotting of 10th and 90th percentile of puff and signal, in addition to mean puff and mean signal;
        added optional plotting of various statistics	

08/08/2019:

in get_puffs:
        renamed function to get_puffs (previously get_puff_statistics). This is mainly to avoid confusion with the new calc_puff_statistics function, and to clarify that the function get_puffs
        does not actually perform any statistical analysis (unlike 
calc_puff_statistics)

09/19/2019:

in __init__:
        added attibutes wtref_mean, c_star, calibration_curve, calibration_factor, full_scale_concentration, full_scale_flow_rate,
 full_scale_ref_length, full_scale_time, full_scale_wtref,
        gas_factor, gas_name, mol_weight, temperature, temperature_K, mass_flow_controller, mass_flow_rate, pressure, ref_height, ref_length, scaling_factor, scaling_factor, standard_pressure, 
        R (universal Gas constant), and __check_sum.
in full_scale:
       new function, based on to_full_scale in PointConcentration.py,
 which converts data to full scale. Note that this function, unlike the analogous function in PointConcentration.py,
       overwrites the variables in the dictionary.
in ambient_conditions:
       new function,based on ambient_conditions in PointConcentration.py,
 which collects ambient conditions during measurement. Pressure in [Pa]!
in scaling_information:
       new function, based on scaling_information in PointConcentration.py,
 which collects scaling data. Units (where applicable) is [m]
in tracer_information:
       new function, based on tracer_information in PointConcentration.py,
 which collects scaling data. Units (where applicable) is [m].
in full_scale_information:
       new function, based on full_scale_information in PointConcentration.py,
 which collects information on full scale information. full_scale_wtref is in [m/s], full_scale_flow_rate takes 
       input in [kg/s], and outputs the flow rate in m^3/s, adjusted
 to ambient conditions	
in convert_tempterature:
       new function, based on convert_temperature in PointConcentration.py,
 which converts temperature, Also edited code to account for removal of variable kelvin_temperature.  
in calc_model_mass_flow_rate:
       new function, based on calc_model_mass_flow_rate in PointConcentration.py,
 which calculates the model scale mass flow rate in [kg/s]	
in calc_full_scale_flow_rate:
       new function, based on calc_full_scale_mass_flow_rate in PointConcentration.py,
 which calculates the full scale mass flow rate in [m^3/s]
in calc_c_star:
       new function, based on calc_c_star in PointConcentration.py,
 which calculates the dimensionless concentration [unitless]	
in calc_full_scale_concentration:
       new function, based on calc_full_scale_concentration in
 PointConcentration.py, which calculates the full scale concentration [ppmV]
in calc_wtref_mean:
       new function, based on calc_wtref_mean in
 PointConcentration.py, which calculates the scaled wtref mean in [m/s]	
in calc_full_scale_time:
       new function, based on calc_full_scale_time in
 PointConcentration.py, which calculates the full scale time step in [s]
in calc_net_concentration:
       re-added __check_sum variable to allow for full scale analysis

09/23/2019:

in plot_puff:
       added proper plotting of full scale data   
in plot_mean_puff:
       added proper plotting of full scale data and revisions to plot formatting, including larger figure, larger text, larger markers, and specified tickmark locations 

09/26/2019:

in __init__: 	   
       added attribute begin_release_index_masked
in offset_correction:
       calculate and subtract offset individually for each puff,
 and perform analysis for all puffs, not just the first 200. Based on algorithm 
developed by Rasmus Fischer;
       see Dissertation of Rasmus Fischer (2011) for more
 details. 
in apply_threshold_concentration:
       added masked release index, to allow for proper masking of multiple variables. 
in apply_threshold_dosage:
       added masked release index, to allow for proper masking of multiple variables. 
in plot_puff:
       added units to dosage and time
in plot_mean_puff:
       added units to dosage and time 
	   
09/27/2019:

in __init__:
       added attribute dosage_unmasked, mask_full, dt
in get_dosage:
       compute dt, which is the sampling time interval, and multiply by dt to fix units of dosage. 
in detect_leaving_time:
       multiply by dt to fix units
in detect_arrival_time:
       multiply by dt to fix units
in get_mask:
       added seperate variables which return both the indeces of data points masked (self.mask_full)
	   and puffs masked (self.mask) sepeartely. Also added logging of msaked data points and puffs; 
       set data before first puff to nan. This ensures it is removed by the threshold dosage and threshold dosage.
in apply_threshold_concentration: 
       fixed logger, now outputs outliers based on size of applicable mask variables, which is equivalent to the
	   number of data points which are masked by applying the threshold concentration. Also outputs number of
	   ensembles which fall below the threshold dosage sepeartely from the total number of masked data points;
       set data before first puff to nan. This ensures it is removed by the threshold concentration. 
in apply_threshold_dosage:
       fixed logger, now outputs outliers based on size of applicable mask variables, which is equivalent to the
       number of data points which are masked by applying the threshold dosage. Also outputs number of ensembles
       which fall below the threshold dosage sepeartely from the total number of masked data points. Further, the
       logger now counts all data points/ensembles whcih fall below the threshold dosage (previously only data
       points which were below the threshold dosage but above the threshold concentration were counted);
       set data before first puff to nan. This ensures it is removed by the threshold dosage.  

10/01/2019:

in save2file_ms: 
       new function, similar to save2file_ms in PuffConcentration.py, saves model scale data to file, for
       (among other things),
 plotting the data in Tecplot. Generates a total of 2 different txt files, for
       puff data, and basic statistics.
 Note that data here is dimensional. 	  

10/04/2019:

in save2file_ms:
       added proper labeling of rows and columns in txt files to make them more readable.
	   
10/07/2019:

in save2file_ms:
       original _ms files depreciated in favor of puff files	
	   
10/14/2019:

in __init__:
       added attributes arrival_index and leaving_index
in detect_leaving_time:
       write index of leaving time to attribute leaving_index
in detect_arrival_time:
       write index of arrival time to attribute arrival_index
in get_peak_concentration:
       constrain peak concentration to between arrival time and leaving time	
in get_peak_time: 
       find peak time of peak concentration constrained to between arrival time and leaving time, consistent with
       updated algorithm for finding peak concentration (see function get_peak_concentration). Also added logging
       of puffs which have a concentration greater than the peak concentration before the arrival time or after
	   the leaving time. 
	   
10/18/2019:

in read_ambient_conditions:
       new function, which reads ambient conditions during from  seperated csv file. Assumes that the csv data is
       located in the same folder as the measurement data, and that each column represents an input variable, and
       each row represents a dataset. If no such file exists in the data directory, function does nothing and
       instead ambient conditions are read from input variables in example_puff_measurement.py.	 
in get_peak_time:
       set peak time to nan if peak concentration is nan. This occurs, among other scenarios, if the arrival time
       and the leaving time are identical. 	

10/21/2019:

in read_ambient_conditions.py:
        new function, which populates the individual variables representing the ambient conditions data based on
		data in ambient_conditions. Requires get_ambient_conditions to be called before calling function, further
		requires that get_ambient_conditions sucessfullly outputs the ambient_conditions array (i.e. requires the
		csv file containing the ambient conditions data to be in the proper format and location). function also
		assumes that csv file (and thus the ambient_conditions array) is in the correct format and contains the
		variables x,y,z,pressure,temperature,wdir,calibration_curve,mass_flow_controller, calibration_factor,
        scaling_factor,scale,ref_length,ref_height,gas_name,mol_weight,gas_factor,full_scale_wtref,and 
		full_scale_flow_rate.	 

10/25/2019:

in plot_puff:
        fix error in x-axis tickmark labeling. Previously only tickmark labels converted to np.int, which led to
        incorrect tickmark labeling.	
in plot_mean_puff:
        fix error in x-axis tickmark labeling. Previously only tickmark labels converted to np.int, which led to
        incorrect tickmark labeling.			   		

10/28/2019: 

in detect_end_release:
        ignore last puff if puff release extends beyond of the timeseries.	

11/12/2019:

in plot_puff:
        plot signal at height of individual peak concentration for each puff,
        not at height of mean peak concentration	
		 
01/10/2020:

in get_peak_time:
        delete puffs for which the arrival time is equal to the leaving time 

01/10/2020: 

in to_non_dimensional:
	new function, which converts data to non-dimensional values. Note that this function, like
	to_full_scale, overwrites the variables in the dictionary.
in calc_non_dimensional_flow_rate:
        new function, based on calc_full_scale_mass_flow_rate in PointConcentration.py,
        which calculates the non-dimensional mass flow rate in [-]	 
in calc_non_dimensional_time:
        new function, based on calc_full_scale_time, which
        calculates the non-dimensional time step [-]	 

01/28/2020:

in plot_puff:
        added plotting of non-dimensional data;
        added plotting of arrival and leaving time for first puff    
in plot_mean_puff:
        added plotting of non-dimensional data 

02/04/2020:

in detect_leaving_time:
        add variable 'time_threshold' to control dosage threshold for determinig leaving and arrival
        time to use for computing characteristic puff times. Note that a 'default' agreed-upon value
        for this variable 5%, however,this fails to properly capture the start and end times of several puffs.
in detect_arrival_time:
        add variable 'time_threshold' to control dosage threshold for determinig leaving and arrival
        time to use for computing characteristic puff times. Note that a 'default' agreed-upon value
        for this variable 5%, however,this fails to properly capture the start and end times of several puffs.
in plot_puff:
        set x_tick_step to 500 for non-dimensional data. 
in plot_mean_puff:
        set x_tick_step to 500 for non-dimensional data. 

02/18/2020:

in plot_puff:
        increase index by 1, even if puff is skipped. Previous version caused incorrect plotting of arrival and leaving times.      
        
Ensemble_Analysis.py

08/05/2019:

in calc_net_concentration:
	changed 'net_concentration' from key to atribute,
 in alignment with the PointConcentration.py script
in apply_threshold_dosage:
	new function, nearly identical to 'apply_threshold_concentration' above, but 
here a threshold total dosage, in place of a threshold peak concentration, is used.
in clear_zeros:
	new function, nearly identical to 'clear_zeros' in PointConcentration.py, but
 the varibale 'signal' is also masked. Further, since no full scale concentration exists as this
 
	point, the mask is based on the sign of net_concentration. Aditionally, masking operation is not performed on the variables, c*, full_scale concentration, and full_scale_time,
	which are 
(currently) not calculated in this script. 
in plot_puff: 
	new function, which plots time series of slected variables

08/08/2019: 

in __init__:
        added attributes ensemble_max, ensemble_mean, ensemble_var, and ensemble_std	
in get_ensemble_min:
        revised indexing algorithm, to ensure that correct minimum is written to the ensemble_min array
in get_ensemble_max: 
        new function, virtually the same as get_ensemble_min, but for maximum value of each ensemble
in get_ensemble_mean:
        new function, virtually the same as get_ensemble_min and get_ensemble_max, but for mean value of each ensemble
in get_ensemble_variance:
        new function, virtually the same as get_ensemble_min, get_ensemble_max, and get_ensemble_mean but for variance of each ensemble.
in calc_puff_statistics: 
        new function, similar to CalculatePuffStatistics in the original C program written by Anne Philip, which calculates basic statistical 
        parameters of puff variables from output of get_puffs. See Bachelor Thesis of Anne Philipp (2010) for more details. 
in plot_convergence_ensemble:
        new function, plots convergence analysis of puff data. Requires get_ensemble_mean to be called before calling function. Variable conv_step plots every conv_step ensemble numbers, 
        default configuration plots all ensemble numbers. 

08/12/2019:

in __init__:
        added attributes class_width, n_classes, class_min, class_max, class_center
in calc_n_classes:
        new function, which calculates the number of classes for each ensemble size. Similar to subsection of CalculationAndWritingFrequencyDistribution
 function from original C program
        written by Anne Philip. See Bachelor Thesis of Anne Philipp (2010) for more details. 
in calc_class_width:
        new function, which calculates the width of the classes. 
in calc_class_boundaries:
        new function, which computes the boundaries, and the center of each class, as a function
 of ensemble size and number. Requires calc_class_width to be called before calling function. 

08/13/2019:

in __init__:
        added attribute n_classes_raw, class_freq, class_freq_norm
in calc_class_width:
	odified code to also save variable n_classes raw to dictionary. Note that n_classes actually represents twice the number of classes,
 and n_classes raw represents the actual number of
        classes. This confusing convention was taken directly from Anne Philip's original C program, and may
 be changed in the future upon consultation with Bernd Leitl and/or Frank Harms;
        added multiplicative factor of 2 to class widhth, as in function mantioned above from Anne Philip's C program.  	
in calc_class_boundaries
	added error messages if n_classes is not an integer. Also changed code to iterate only to n_classes_raw instead of n_classes. This is because n_classes 
actually represents twice the
	number of classes, and n_classes raw represents the actual number of classes. This confusing convention was taken directly from Anne Philip's
 original C program, and may be simplified
	in the future upon consultation with Bernd Leitl and/or Frank Harms.    
in get_class_frequency:
	new function, which returns the number of data points inside each class. Similar to subsection of CalculationAndWritingFrequencyDistribution
 function from original C program written
	by Anne Philip. See Bachelor Thesis of Anne Philipp (2010) for more details. Requires calc_class_boundaries and calc_class_width 
to be called before calling function. 
in plot_class_statstics:
	new function, plots histogram using the class analysis for all data. Uses normalized class frequencies, and creates one
	plot for each ensemble size. Requires get_class_frequency to
	be called before calling function.

09/18/2019:

in plot_class_statistics:
        adjusted code to make sure data is properly plotted, formatted plots. Unable to specifiy title pad, presumably due to old matplotlib version.
		
09/23/2019:

in plot_convergence_ensemble:
        revisions to plot formatting, including larger figure, larger text, larger markers, and specified tickmark locations
in plot_class_statistics:
        revisions to plot formatting, specifically a larger figure size	

10/01/2019:

in __init__:
       added attributes x, y, z, scale, calibration_curve, calibration_factor,
 mass_flow_controller, ref_height, ref_length,
       scaling_factor, gas_factor, gas_name, mol_weight,
 temperature, standard_temp_K, and pressure
in ambient_conditions: new function, same as ambient_conditions in PuffConcentration.py,
 which collects ambient conditions 
       during measurement. Pressure in [Pa]!
in scaling_information: 
       new function, same as scaling_information in PuffConcentration.py,
 which collects scaling data. Units (where applicable) is [m]
in tracer_information: 
       new function, same as tracer_information in PuffConcentration.py,
 which collects tracer information. Units (where applicable) is [m].
in save2file_ms_ensemble: 
       new function, similar to save2file_ms in PuffConcentration.py, saves model scale data to file, for (among other things),

       plotting the data in Tecplot. Generates a total of 2 different txt files, for convergence analysis data and class data.
       Note that data here is dimensional. 

10/04/2019:

in save2file_ms_ensemble:
       added proper labeling of rows and columns in txt files to make them more readable.	
	   
10/07/2019:

in save2file_ms_ensemble:
       removed nan values, since these won't load properly in tecplot.	

01/28/2020:

in plot_convergence_ensemble:
       added support for non-dimensional data
in plot_class_statistics:
       added support for non-dimensional data

import numpy as np
import math
import logging
import os
import pandas as pd
import windtunnel as wt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import string

# Create logger
logger = logging.getLogger()
__all__ = ['PointConcentration']


# %%#
class PointConcentration(pd.DataFrame):
    """ PointConcentration is a class that holds data collected during
    a continuous release point concentration measurement. The class can hold
    the raw time series, the corresponding wtref and all other quantities
    necessary to analyse the time series. All the information in a
    PointConcentration object can be saved to a txt file.
    @parameter: time, type = np.array
    @parameter: wtref, type = np.array
    @parameter: fast_FID, type = np.array
    @parameter: slow_FID, type = np.array
    @parameter: open_rate, type = np.array"""

    def __init__(self, time, wtref, slow_FID, fast_FID, open_rate):
        #edit 09/19/2019: removed variable Kelvin_temp. This variable makes no sense, 
        #as it is simply a unit conversion factor equivalent to the 0°C in K used to
        #convert temperatures between °C and K, and has been replaced by adjustments
        #to the convert_temperature function.   Also removed duplicate variables
        #full_scale_time and scale. 
        """ Initialise PointConcentration object. """
        super().__init__()

        self['slow_FID'] = pd.Series(data=slow_FID)
        self['fast_FID'] = pd.Series(data=fast_FID)

        self.x = None
        self.y = None
        self.z = None
        self.scale = None
        self.wtref_mean = None
        self.open_rate = open_rate
        self.time = time
        self.wtref = wtref
        self.net_concentration = None
        self.c_star = None
        self.calibration_curve = None
        self.calibration_factor = None
        self.full_scale_concentration = None
        self.full_scale_flow_rate = None
        self.full_scale_ref_length = None
        self.full_scale_time = None
        self.full_scale_wtref = None
        self.gas_factor = None
        self.gas_name = None
        self.mol_weight = None
        self.temperature = None
        self.temperature_K = None
        self.mass_flow_controller = None
        self.mass_flow_rate = None
        self.pressure = None
        self.ref_height = None
        self.ref_length = None
        self.scaling_factor = None
        self.standard_temp_K = None
        self.standard_temp = 20  # [°C]
        self.standard_pressure = 101325  # [Pa]
        self.R = 8.3144621  # universal gas constant [kJ/kgK]
        self.__check_sum = 0

    def __repr__(self):
        """ Return the x, y and z coordinate of the PointConcentration
        object. """
        return 'PointConcentration (x={x}, y={y}, z={z})'.format(x=self.x,
                                                                 y=self.y,
                                                                 z=self.z)

    def __eq__(self, other):
        """ Two PointConcentration objects are considered equal, if their x, y
        and z coordinates are the same. """
        return self.x == other.x and self.y == other.y and self.z == other.z

    @classmethod
    def from_file(cls, filename):
        """ Create PointConcentration object from file. open_rate is converted
        to %."""
	    #edit 07/22/2019: data does not usually have six colums. 
        #TODO: open rate = flow rate?		
        time, wtref, slow_FID, fast_FID, open_rate = np.genfromtxt(filename,
                                                                   usecols=(0, 1, 2, 3, 4),
                                                                   unpack=True)

        return cls(time, wtref, slow_FID, fast_FID, open_rate * 10)

    def to_full_scale(self):
        """ Return all quantities to full scale. Requires XXXXXX to be
        specified."""
		#edit 07/26/2019: changed clear_zeros to a standalone funciton, which is now called from 
		#outside calc_net_concentration
        if self.__check_sum >= 8:

            quantities = ['x', 'y', 'z', 'time', 'concentration', 'flow rate']
            your_measurement = {}
            your_measurement.fromkeys(quantities)

            your_measurement['x'] = self.x = self.x * self.scale / 1000  # [m]
            your_measurement['y'] = self.y = self.y * self.scale / 1000  # [m]
            your_measurement['z'] = self.z = self.z * self.scale / 1000  # [m]
            your_measurement['flow rate'] = self.calc_full_scale_flow_rate()

            self.calc_full_scale_time()			
            self.calc_full_scale_concentration()

            your_measurement['time'] = self.full_scale_time

            your_measurement['concentration'] = \
                self.full_scale_concentration

            return your_measurement

        else:
            raise Exception('Please enter or calculate all full scale data '
                            'necessary!')
							
    def get_ambient_conditions(path=None,name=None,input_file=None):
        """Read ambient conditions from csv file. If no such file exists, function
		does nothing and instead ambient conditions are read from values in
		example_puff_measurement.py."""	
		#edit 10/21/2019: new function, similar to get_ambient_conditions in PuffConcentration.py
        #which reads ambient conditions during measurement from seperate csv file. Assumes that the
        #csv data is located in the same folder as the measurement data, and that each column represents
        #an input variable, and each row represents a dataset. If no such file exists in the data directory,
        #function does nothing and instead ambient conditions are read from input variables
		#in example_puff_measurement.py.
        if input_file==None:
           print('Warning: Input csv filename (for ambient conditions) not specified. Resorting to input data in example_puff_measurement.py')
           return		   
        elif name==None:
           print('Warning: Name of dataset not specified. Cannot attempt to locate csv file containing ambient conditions data. Resorting\
to input data in example_puff_measurement.py')
        elif path==None:
           print('Warning: Path of input csv file (for ambient conditions) not specified. Resorting to input data in example_puff_measurement.py')			   
           return
        elif not os.path.exists(input_file):
           print('Error: Cannont find csv file containing ambient conditions in specified directory. Check name and/or location of ambient \
conditions file. Resorting to input data in example_puff_measurement.py')	
           return		   
        else:	
           ambient_conditions=pd.read_csv(input_file,sep=',',index_col=0) 
        
        if name not in ambient_conditions.keys():
           print('Error: Dataset not found in csv file. Check to make sure that csv file to make sure that the csv file contains all necessary \
data and is properly formatted. Resorting to input data in example_puff_measurement.py')
           return	

        #list of all variables output by read_ambient_conditions fuction.  
        necessary_keys={'x','y','z','pressure','temperature','calibration_curve','mass_flow_controller','calibration_factor', \
        'scaling_factor','scale','ref_length','ref_height','gas_name','mol_weight','gas_factor','full_scale_wtref','full_scale_flow_rate' }
        if not all(name2 in ambient_conditions[name] for name2 in necessary_keys):
           print('Error: csv file does not contain all necessary ambient conditions data. Check to make sure that csv file to make sure that \
the csv file contains all necessary data and is properly formatted. Resorting to input data in example_puff_measurement.py')
           return			   
       		

        return ambient_conditions	

    def read_ambient_conditions(ambient_conditions,name):
        """Populate individual variables representing ambient conditions based on data
		in ambient_conditions array. """	
        #edit 10/21/2019: new function, similar to read_ambient_conditions in PuffConcentration.oy which
        #populates the individual variables representing the ambient conditions data based on data in
        #ambient_conditions. Requires get_ambient_conditions to be called before calling function, further
        #requires that get_ambient_conditions sucessfullly outputs the ambient_conditions array (i.e. requires
        #the csv file containing the ambient conditions data to be in the proper format and location). function
        #also assumes that csv file (and thus the ambient_conditions array) is in the correct format and contains
        #the variables x,y,z,pressure,temperature,calibration_curve,mass_flow_controller, calibration_factor,
        #scaling_factor,scale,ref_length,ref_height,gas_name,mol_weight,gas_factor,full_scale_wtref,and 
        #full_scale_flow_rate.
 	   
        x=None if ambient_conditions[name]['y'] =='None' else np.float(ambient_conditions[name]['x'])
        y=None if ambient_conditions[name]['y'] =='None' else np.float(ambient_conditions[name]['y'])
        z=None if ambient_conditions[name]['x'] =='None' else np.float(ambient_conditions[name]['z'])
        pressure=None if ambient_conditions[name]['pressure'] =='None' else np.float(ambient_conditions[name]['pressure'])		
        temperature=None if ambient_conditions[name]['temperature'] =='None' else np.float(ambient_conditions[name]['temperature'])
        calibration_curve=None if ambient_conditions[name]['calibration_curve'] =='None' else np.float(ambient_conditions[name]['calibration_curve'])
        mass_flow_controller=None if ambient_conditions[name]['mass_flow_controller'] =='None' else ambient_conditions[name]['mass_flow_controller']
        calibration_factor=None if ambient_conditions[name]['calibration_factor'] =='None' else np.float(ambient_conditions[name]['calibration_factor'])
        scaling_factor=None if ambient_conditions[name]['scaling_factor'] =='None' else np.float(ambient_conditions[name]['scaling_factor'])	
        scale=None if ambient_conditions[name]['scale'] =='None' else np.float(ambient_conditions[name]['scale'])
        ref_length=None if ambient_conditions[name]['ref_length'] =='None' else np.float(eval(ambient_conditions[name]['ref_length']))
        ref_height=None if ambient_conditions[name]['ref_height'] =='None' else np.float(ambient_conditions[name]['ref_height'])	
        gas_name=None if ambient_conditions[name]['gas_name'] =='None' else ambient_conditions[name]['gas_name']
        mol_weight=None if ambient_conditions[name]['mol_weight'] =='None' else np.float(ambient_conditions[name]['mol_weight'])
        gas_factor=None if ambient_conditions[name]['gas_factor'] =='None' else np.float(ambient_conditions[name]['gas_factor'])
        full_scale_wtref=None if ambient_conditions[name]['full_scale_wtref'] =='None' else np.float(ambient_conditions[name]['full_scale_wtref'])
        full_scale_flow_rate=None if ambient_conditions[name]['full_scale_flow_rate'] =='None' else np.float(ambient_conditions[name]['full_scale_flow_rate'])	
		
        return x,y,z,pressure,temperature,calibration_curve,mass_flow_controller,\
        calibration_factor, scaling_factor,scale,ref_length,ref_height,gas_name,mol_weight,\
        gas_factor,full_scale_wtref,full_scale_flow_rate							

    def ambient_conditions(self, x, y, z, pressure, temperature, calibration_curve,
                           mass_flow_controller, calibration_factor=0):
        """ Collect ambient conditions during measurement. pressure in [Pa],
        temperature in [°C]. """
        self.__check_sum = self.__check_sum + 1

        self.x = x
        self.y = y
        self.z = z
        self.pressure = pressure
        self.temperature = temperature
        self.calibration_curve = calibration_curve
        self.calibration_factor = calibration_factor
        self.mass_flow_controller = mass_flow_controller

    def scaling_information(self, scaling_factor, scale, ref_length, ref_height):
        """ Collect data necessary to scale the results. unit: [m], where
        applicable."""
        self.__check_sum = self.__check_sum + 1

        self.scaling_factor = scaling_factor
        self.scale = scale
        self.ref_length = ref_length
        self.ref_height = ref_height
        self.full_scale_ref_length = self.scale * self.ref_length

    def tracer_information(self, gas_name, mol_weight, gas_factor):
        """ Collect information on tracer gas used during measurement.
        Molecular weight in [kg/mol]. """
        #edit 10/24/2019: highly doubtful that weighth of tracer gas is in kg/mol.
        #Typical units for molecular weight are g/mol, 28.97 kg/mol seems outrageously 
        #large a value. 
        #TODO: verify units of molecular weight of gas		
        self.__check_sum = self.__check_sum + 1

        self.gas_name = gas_name
        self.mol_weight = mol_weight
        self.gas_factor = gas_factor

    def full_scale_information(self, full_scale_wtref, full_scale_flow_rate):
        """ Collect information on desired full scale information.
        full_scale_wtref in [m/s]. full_scale_flow_rate is automatically
        adjusted to standard atmosphere conditions.
        input in [kg/s], output in [m^3/s]. """
        self.__check_sum = self.__check_sum + 1

        self.full_scale_wtref = full_scale_wtref
        self.full_scale_flow_rate = full_scale_flow_rate

    def convert_temperature(self):
        """ Convert ambient temperature to °K. """
        #edit 09/19/2019: edited code to account for removal 
        #of variable kelvin_temperature. 		
        self.temperature_K = self.temperature + 273.15
        self.standard_temp_K = self.standard_temp + 273.15

    def calc_model_mass_flow_rate(self):
        """ Calculate the model scale flow rate in [kg/s]. """
        self.__check_sum = self.__check_sum + 1

        self.mass_flow_rate = self.gas_factor * (np.mean(self.open_rate) *
                                                 self.calibration_curve +
                                                 self.calibration_factor) * \
                              self.temperature_K * self.standard_pressure / \
                              (self.pressure * self.standard_temp_K)

        return self.mass_flow_rate

    def calc_full_scale_flow_rate(self):
        """ Convert flow rate to full scale flow rate in [m^3/s]. """
        self.full_scale_flow_rate = (self.full_scale_flow_rate * self.R *
                                     self.standard_temp_K) / \
                                    (self.standard_pressure * self.mol_weight)

        return self.full_scale_flow_rate

    def calc_net_concentration(self):
        """ Calculate net concentration in [ppmV]. """
        self.__check_sum = self.__check_sum + 1

        self.net_concentration = self.fast_FID - self.slow_FID

        return self.net_concentration

    def calc_c_star(self):
        """ Calculate dimensionless concentration. [-] """
        self.__check_sum = self.__check_sum + 1
        # TODO: calc_mass_flow_rate (for Point, Line and Area)
        self.c_star = self.net_concentration * self.wtref_mean * \
                      self.ref_length ** 2 / self.mass_flow_rate * 1000 * 3600

        return self.c_star

    def calc_full_scale_concentration(self):
        """ Calculate full scale concentration in [ppmV]. """
        self.full_scale_concentration = self.c_star * \
                                        self.full_scale_flow_rate / \
                                        (self.full_scale_ref_length ** 2 *
                                         self.full_scale_wtref)	
		
        return self.full_scale_concentration

    def calc_wtref_mean(self):
        """ Calculate scaled wtref mean in [m/s]. """
        self.__check_sum = self.__check_sum + 1

        self.wtref_mean = self.scaling_factor * np.mean(self.wtref)

        return self.wtref_mean

    def calc_full_scale_time(self):
        """ Calculate full scale timesteps in [s]. """
        #edit 10/24/2019: replaced full_scale_ref_length/ref_length with scale. This should be equivalent, and is likely to appear more logical	
        if self.wtref_mean is None:
            self.wtref_mean = PointConcentration.calc_wtref_mean()

        #self.full_scale_time = self.full_scale_ref_length / self.ref_length * \
                               #self.wtref_mean / self.full_scale_wtref * \
                               #self.time
        					   
        self.full_scale_time = self.scale * self.wtref_mean / self.full_scale_wtref * \
                               self.time

		

        return self.full_scale_time

    def clear_zeros(self):
	    
        """ Clear and count zeros in concentration measurements."""
		#edit 07/22/2019: convert mask to array to apply it to full_scale_time. 
		#Aditionally, mask c_star, and net concentration. Also changed script so 
		#that the negative values are removed from the time series, not the positive
		#ones as done previosuly, since negative concetration values don't make any sense. 
        concentration_size = np.size(self.net_concentration)

        # Mask zeros
        mask = self.net_concentration > 0
        
        self.time = self.time[np.asarray(mask)]			
        self.net_concentration = self.net_concentration[mask]			


        # Log outliers in console and to file
        logger.info('Values below 0: {} or {:.4f}%'.format(
            np.size(np.where(~mask)),
            np.size(np.where(~mask)) / concentration_size * 100
        ))
		
    def plot_hist_conc(self,n_classes=None,var='net_concentration',path=None,name=None):
        """Creates a historgram point concentration, i.e. continuous release, data."""
        #edit 10/17/2019: New function, largely based off of plot_hist in bl.py (for plotting wind data), 
        #as well as plot_class_statistics in EnsembleAnalysis.py (for plotting puff data), which plots
        #histogram of point concentration (continuous release) data.
        #TODO: add units to class mean label		

       

        data=getattr(self,var)
        data_size=np.size(data)		
        #if n_classes not specified, calculate the number of classes. Algorithm to determine number of
        #classes based on CalculationAndWritingFrequencyDistribution function from old C program written
        #by Anne Philip, and also used in calc_n_classes function in Enesemble_Analysis.py (for use with
        #puff data). 

        if n_classes == None:
           n_classes=np.int(1+math.log10(data_size)/math.log10(2))	

        #calculate class width, class min, class max, and class mean
        class_width=(np.max(data)-np.min(data))/n_classes
        class_min=[np.min(data)]
        class_max=[np.min(data)+class_width]		
        for i in range(n_classes-1):
            class_min=np.append(class_min,class_min[i]+class_width)
            class_max=np.append(class_max,class_max[i]+class_width)	
        class_mean=(class_min+class_max)/2
			
		
        #get class frequency, both absolute and normalized
        class_freq=np.zeros(np.shape(class_min),dtype=np.int)	
        class_freq_cum=np.zeros(np.shape(class_min),dtype=np.int)
        class_freq_norm=np.zeros(np.shape(class_min),dtype=np.int)	
        class_freq_cum_norm=np.zeros(np.shape(class_min),dtype=np.int)			
        for i in range(n_classes):		
            class_freq[i]=((class_min[i] <= data) & (data < class_max[i])).sum()	
            class_freq_cum[i]= (data < class_max[i]).sum()	
		
        class_freq_norm=class_freq/data_size	
        class_freq_cum_norm= class_freq_cum/data_size	
	
        #plot frequency distribution
        
        var_label=var.replace("_"," ")		
        ret=plt.figure(301)			
        plt.clf()	
        plt.bar(np.linspace(1,np.shape(class_freq_norm)[0],np.shape(class_freq_norm)[0]),class_freq_norm,width=1,align='center')		
        ax=plt.gca()		
        ax.set_title('Frequency Distribution (Model Scale) of '+string.capwords(var_label),fontsize=40)              			
                				
        ret.set_figwidth(26)	
        ret.set_figheight(16)								
        n=n_classes		
        plt.xticks(np.linspace(1,n,n),np.round(class_mean,1))	
        plt.tick_params(axis='both', labelsize=30)
        		
        plt.xlim(0.5,n+0.5)	
        plt.ylim(0,1)
        plt.xlabel('Class Mean',fontsize=40)
        plt.ylabel('Frequency',fontsize=40)		   		
        if path=='none':
           []
        else:	
           if name=='none':
              print('Name of dataset not specified. Plot of frequency distribution will not be saved to avoid confusion in the future.')	
           else:		
              plt.savefig(path + var + '.jpg')
			  
        #plot cumulative frequency distribution
		
        ret=plt.figure(302)			
        plt.clf()	
        plt.bar(np.linspace(1,np.shape(class_freq_norm)[0],np.shape(class_freq_norm)[0]),class_freq_cum_norm,width=1,align='center')		
        ax=plt.gca()		
        ax.set_title('Cumulative Frequency Distribution (Model Scale) of '+string.capwords(var_label),fontsize=40)              			
                				
        ret.set_figwidth(26)	
        ret.set_figheight(16)								
        n=n_classes		
        plt.xticks(np.linspace(1,n,n),np.round(class_mean,1))	
        plt.tick_params(axis='both', labelsize=30)
        		
        plt.xlim(0.5,n+0.5)	
        plt.ylim(0,1)
        plt.xlabel('Class Mean',fontsize=40)
        plt.ylabel('Cumulative Frequency',fontsize=40)		   		
        if path=='none':
           []
        else:	
           if name=='none':
              print('Name of dataset not specified. Plot of frequency distribution will not be saved to avoid confusion in the future.')	
           else:		
              plt.savefig(path + var + '_cumulative.jpg')		

        return ret		

    def save2file_ms(self, filename, out_dir=None):
        """ Save model scale data from PointConcentration object to txt file.
        filename must include '.txt' ending. If no out_dir directory is
        provided './' is set as standard.
        @parameter: filename, type = str
        @parameter: out_dir, type = str"""
        if out_dir is None:
            out_dir = './'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        output_file = out_dir + '_ms_' + filename		
        np.savetxt(output_file, np.vstack((self.time,
                                           self.c_star,
                                           self.net_concentration)
                                          ).transpose(),
                   fmt='%.4f',
                   header="General concentration measurement data:" + '\n' +
                          "" + '\n' +
                          "geometric scale: 1:{}".format(float(self.scale))
                          + "" + '\n' +
                          "Variables: x: {} [mm], y: {} [mm], z: {} [mm], "
                          "ambient temperature: {:.1f} [°C], ambient pressure: {:.2f} [Pa],"
                          " mass flow rate {:.4f} [kg/s], "
                          "reference length (model): {:.4f} [m], "
                          "Tracer gas: {}, mol. weight tracer: {:.4f} [mol/kg], "
                          "gas factor: {:.6f}, calibartion curve: {:.6f}, "
                          "wtref: {:.4f} [m/s]".format(self.x, self.y, self.z,
                                                       self.temperature,
                                                       self.pressure,
                                                       self.mass_flow_rate,
                                                       self.ref_length,
                                                       self.gas_name,
                                                       self.mol_weight,
                                                       self.gas_factor,
                                                       self.calibration_curve,
                                                       self.wtref_mean)
                          + "" + '\n' +
                          "\"time [ms]\" \"c_star [-]\" \"net_concentration [ppmV]\" ")

    def save2file_fs(self, filename, out_dir=None):
        """ Save full scale and model scale data from PointConcentration object
        to txt file. filename must include '.txt' ending. If no out_dir
        directory is provided './' is set as standard.
        @parameter: filename, type = str
        @parameter: out_dir, type = str"""
        if self.__check_sum < 8:
            raise Exception('Please enter or calculate all full scale data '
                            'necessary!')

        if out_dir is None:
            out_dir = './'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        output_file = out_dir + '_fs_' + filename		
        np.savetxt(output_file, np.vstack((self.full_scale_time,
                                           self.c_star,
                                           self.net_concentration,
                                           self.full_scale_concentration)
                                          ).transpose(),
                   fmt='%.4f',
                   header="General concentration measurement data:" + '\n' +
                          "" + '\n' +
                          "geometric scale: 1:{}".format(float(self.scale))
                          + "" + '\n' +
                          "Variables: x: {} [m], y: {} [m], z: {} [m], "
                          "ambient temperature: {:.1f} [°C], ambient pressure: {:.2f} [Pa],"
                          " mass flow rate {:.4f} [kg/s], "
                          "reference length (model): {:.4f} [m], "
                          "reference length (full-scale): {:.4f} [m], Tracer gas: {}, "
                          "mol. weight tracer: {:.4f} [mol/kg], "
                          "gas factor: {:.6f}, calibartion curve: {:.6f}, "
                          "wtref: {:.4f} [m/s], full scale flow rate: {:.4f} [m^3/s]".format(
                              self.x,
                              self.y,
                              self.z,
                              self.temperature,
                              self.pressure,
                              self.mass_flow_rate,
                              self.ref_length,
                              self.full_scale_ref_length,
                              self.gas_name,
                              self.mol_weight,
                              self.gas_factor,
                              self.calibration_curve,
                              self.wtref_mean,
                              self.full_scale_flow_rate)
                          + "" + '\n' +
                          "\"full scale time [s]\" \"c_star [-]\" "
                          "\"net_concentration [ppmV]\" \"full_scale_concentration [ppmV]\"")

    def save2file_avg(self, filename, out_dir=None):
        """ Save average full scale and model scale data from
        PointConcentration object to txt file. filename must include '.txt'
        ending. If no out_dir directory is provided './' is set as standard.
        @parameter: filename, type = str
        @parameter: out_dir, type = str"""
        if self.__check_sum < 8:
            raise Exception('Please enter or calculate all full scale data '
                            'necessary!')

        if out_dir is None:
            out_dir = './'
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

        output_file = out_dir + '_avg_' + filename

        np.savetxt(output_file, np.vstack((np.nanmean(self.c_star),
                                           np.nanmean(self.net_concentration),
                                           np.nanmean(
                                               self.full_scale_concentration))
                                          ).transpose(),
                   fmt='%.4f',
                   header="General concentration measurement data:" + '\n' +
                          "" + '\n' +
                          "geometric scale: 1:{}".format(float(self.scale))
                          + "" + '\n' +
                          "Variables: x: {} [m], y: {} [m], z: {} [m], ambient temperature: {:.1f} [°C], "
                          "ambient pressure: {:.2f} [Pa], mass flow rate {:.4f} [kg/s], "
                          "reference length (model): {:.4f} [m], "
                          "reference length (full-scale): {:.4f} [m], Tracer gas: {}, "
                          "mol. weight tracer: {:.4f} [mol/kg], "
                          "gas factor: {:.6f}, calibartion curve: {:.6f}, "
                          "wtref: {:.4f} [m/s], full scale flow rate: {:.4f} [m^3/s]".format(
                              self.x,
                              self.y,
                              self.z,
                              self.temperature,
                              self.pressure,
                              self.mass_flow_rate,
                              self.ref_length,
                              self.full_scale_ref_length,
                              self.gas_name,
                              self.mol_weight,
                              self.gas_factor,
                              self.calibration_curve,
                              self.wtref_mean,
                              self.full_scale_flow_rate)
                          + "" + '\n' +
                          "\"c_star [-]\" \"net_concentration [ppmV]\" "
                          "\"full_scale_concentration [ppmV]\"")

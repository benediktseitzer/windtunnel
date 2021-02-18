# Windtunnel 
The windtunnel python module is a processing library developed by the EWTL Group at the University of Hamburg. 


## 1. Installation
### 1.1 Requirements
 - python 3.0 
 - the following python libraries:
 - numpy v1.13.1
 - scipy v0.19.1
 - pandas v20.3
 - logging v0.5.1.2
 - matplotlib v2.0.2
 - skimage v0.13.0
 - openpyxl v2.4.8
 - tkinter v8.6
 - logging v0.5.1.2
 - csv v1.0
 - tsp_solver 
 - os, sys, time, fnmatch, math, mpl_toolkits (all built-in modules)
 
The listed versions represent minimum requirements. 

PAPE should run fine wih newer versions of the modules listed in Section 1.1. However, backwards compatability, especially with python v2.x, cannot be guaranteed. 
PAPE has been tested tested using Spyder 3.3.6 on a computer running Ubtuntu 20.04LTS, and using Spyder 3.2.3 on a computer running Windows 10, but should work on any distribution that supports the above requirements. Note that the purpose of the GUI installation is to run the program without the use of Spyder. 

### 1.2 Installation

To install PAPE, copy the program files to a directory of your choice. It is recommended to install the program in a directory to which you have read and write permission to avoid potential permission conflicts. If you want to run PAPE via the GUI, make sure that your system is set to allow the execution of python files as a program (varies depending on your OS, see below). 

#### 1.2.1 Windows

GUI mode should be enabled automatically. To verify this, double-click the file PAPE_GUI.pyw. If a warning along the lines of “Windows cannot open this kind of file” appears, select the command along the lines of “Select an app to open this file,” and select the file “pythonw.exe.” Note that this is not the same as “python.exe,” and if you don’t know the difference between the two, you should probably refrain from selecting “python.exe.” If you have trouble finding the “python.exe” file, type “python.exe” into the search bar next to the start menu. Select “open file location,” and the path of the appropriate file will open automatically. If windows cannot find the file “python.exe,” chances are that python is probably not or at least incorrectly installed. If you don't know how to fix this yourself, ask someone who knows! 

**Note**: on some windows systems, an error has been encountered in which python appears to be properly installed, but double clicking the PAPE_GUI files does nothing. If you encounter this issue, try opening python.exe  by itself (which should open the python interpreter). If a message along the lines of “Warning: this Python interpreter is in a conda environment, but the environment has not been activated” appears, chances are your python is not properly installed. Please contact your administrator for help on how to fix this! 

#### 1.2.2 Linux and MacOs
You can pip install all the packages by typing:
```sh
$ pip3 install numpy pandas scipy matplotlib tsp_solver
$ sudo apt-get install python3-tk
```

To allow PAPE to run in GUI mode, right-click the file called PAPE_GUI.pyw” in the PAPE root directory, and make sure that the option “Allow executing the file as program” is selected. Furthermore, under “Preferences → Behavior,” “Executable Text Files” should be set to either “Run them” or “Ask what to do.”

After following these steps with the required dependencies, PAPE should be installed and running on your computer. 

## 2. A priori and input values in windtunnel.py

This section briefly describes the key input and a priori variables of the functions included in the windtunnel.py package. The names of the variables are as given in the PAPE GUI, while  the variable names in the PAPE CLI are given in brackets. The following a priori variables are specified in both the puff and point concentration analysis: 
 - Source Location (x,y,z) [x_source, y_source, z_source]:  x-,y-, and z- position of the emission source, in mm (model scale). Note that this is not necessarily the same as the measurement location (see below). This quantity is a vector and hence consists of three different variables, one for each spatial direction. These variables are represented by three input lines in the GUI, while in CLI mode, the vector components are treated a three distinct variables. Note that while PAPE will allow negative z-coordinates, **any negative values for z_source should be approached with extreme skepticism, as any negative value for z_source implies that the emission source is below the wind tunnel surface.** 
 - Measurement Location (x,y,z) [x_measure, y_measure, z_measure]:  x-,y-, and z- position of the measurement source, in mm (model scale). Note that this is not necessarily the same as the source location (see above). This quantity is a vector and hence consists of three different variables, one for each spatial direction. These variables are represented by three input lines in the GUI, while in CLI mode, the vector components are treated a three distinct variables. Note that while PAPE will allow negative z-coordinates, **any negative values for z_measure should be approached with extreme skepticism, as any negative value for z_measure implies that measurements are taken below the wind tunnel surface.** 
 - Pressure [pressure]: reference air pressure outside wind tunnel during measurement, in Pa.
 - Temperature [temperature]: reference air temperature outside wind tunnel during measurement, in °C. 
 - Calibration Curve [calibration_curve]: tbd
 - Mass Flow Controller [mass_flow_controller]: tbd
 - Calibration Factor [calibration_factor]: tbd
 - Scaling Factor [scaling_factor]: tbd
 - Scale [scale]: scale of model. Note that the number here denotes the real-world units represented by one corresponding model unit; i.e. a scale of 1:250 would equal a value of 250 for the parameter ‘scale.’ Note further that this is not the same as the output scale-controlling variable “full_scale.”
 - Reference Length [ref_length]: reference length, in m (model scale), used to make data non-dimensional. 
 - Reference Height [ref_height]: reference height, in m (model scale). Currently not used.
 - Gas Name [gas_name]: Name of tracer gas. Currently only used for documentation purposes.
 - Molecular Weight [mol_weight]: molecular weight of gas, in kg/mol.  Units here are highly doubtful, **to be verified.** More likely units are g/mol. 
 - Gas Factor [gas_factor]: used to calculate model scale flow rate
 - Full Scale Reference Wind Speed [full_scale_wtref]: full scale reference wind velocity, in m/s, used to calculate full_scale concentration and time. 
 - Full scale mass flow rate [full_scale_flow_rate]: model scale flow rate out of emission source, in kg/s. Note that this is later converted to full scale flow rate, in m3/s. Convention of naming both the full scale and model scale flow rate full_scale_flow_rate taken from Benyamin Schliffke’s original script, the reason for this rather odd naming convention is at this point unknown. 
The following variables are located in the “advanced options” menu of the GUI: 
 - Number of Outliers to Exclude [n_exclude]: number of outliers to exclude. The default setting, None (without quotation marks of any kind), automatically selects the number of outliers to remove. It is explicitly noted that setting this variable to None does not mean that no outliers are removed from the datasets. Furthermore, **changing this variable can substantially compromise the integrity of the results of the program. Do not change the value of this variable unless you know what you are doing!** 
 - y-Axis Range [axis_range]:  determines how the y-axis on the puff plots are determined. Accepts ‘auto’ (Automatic in GUI mode) and ‘same’ as input values. Setting this variable to ‘auto’ will select the axis range individually for each puff plot, while setting this variable to ‘same’ will enforce the y-axis range to be the same for all plots.    
The puff analysis uses all of the a priori variables above, plus the following additional a priori variables:
 - Threshold Concentration [threshold_concentration]: the minimum peak concentration, in ppmv (model scale) of a puff to be included in the statistical analysis. All puffs with a peak concentration below this value will be filtered out of the analysis. 
 - Threshold Dosage [threshold_dosage]: the minimum dosage, in ppmv (model scale), of a puff to be included in the statistical analysis. All puffs with a peak concentration below this value will be filtered out of the analysis. 
 - Select Scale of Data Output [full_scale]: variable which controls which scale to perform the analysis in. Accepts ‘ms’ (for model scale data), ‘fs’ (for full scale data), and ‘nd’ (for non-dimensional data) as inputs. 
 - Full Mode, Basic Mode [functions_mode]: determines the scope of the analysis. Accepts ‘full’ for complete, usually time-consuming, analysis, and ‘basic’ for a trimmed-down, but usually faster, analysis. See Section 5.1.2 for more details.
 - Show Advanced Options Distribution Data with Mean Puff [dist]: Option to display 10th and 90th percentiles of puff time series in mean puff plot. In CLI mode, the variable should be set to ‘off’ (default) or ‘on.’   
 - Number of Puffs to Plot [n_puffs]: Set the number of puffs to plot. The default setting for this variable is 5 in basic mode and ‘all’ in full mode. In GUI mode, the defaults are pre-selected when selecting the program mode, while in CLI mode, the variable should be set to ‘default’ in order to use the default values. Note that in both GUI and CLI mode, this variable can be set to ‘all’ to plot all puffs, regardless of the size of the dataset. 
The following variables are located in the “advanced options” menu of the GUI, and are only used in puff mode analysis mode. 
 - Percentage to Use as Threshold for Computing Arrival and Leaving Times [time_threshold]: a number which determines at which percentages of the total dosage the arrival and leaving times should be computed. The number itself is the fraction of the dosage at which to define the arrival time, with the leaving time determined in a symmetric manner. For example, setting this variable to 0.01 will define the arrival time at 1% of the dosage, and the leaving time at 99% (or 100%-1%) of the dosage. The agreed-upon setting for this variable is 5%, but some puffs cannot be appropriately characterized by this convention. **Warning: changing this variable will substantially modify the results of the program, and can severely compromise the integrity of the resulting analysis. Do not change this variable unless you know what you are doing!!** 
 - signal_mode: determine the percentage of the mean release signal to be tolerated as noise. **Do not modify this variable unless you know what you are doing! Furthermore, there should be no reason to modify this variable unless you encounter problems with the puff detection.** This variable controls the threshold used for determining when the noisy release signal is read as 1 (puff is being released) or 0 (no puff release), as a percentage of the mean signal. In the default setting, this variable is set to 1, which means that whenever the puff release signal (usually in V) is greater than the mean puff release signal, the program interprets that puff gases are being released, while if the puff release signal is less than the mean release signal, the program interprets that no puff gases are being released. If this variable if set to 0.1, the program interprets that puff gases are being released if the release signal is more than 10% of the mean, and no gas is being released if the signal is less than 10% of the mean, etc. Effectively, this parameter sets the amount of permitted noise in the release signal. Since PAPE uses the binary form of the release signal to determine where the puffs start and end, **an incorrect value for this variable can severely compromise the integrity of the resulting analysis!** In particular, if this value if too low, the program may detect puffs which are not actually puffs but just a noisy release signal, while if this value if too high, a single puff release event may be split into multiple puffs.    

## 3. How to work with windtunnel.py

A key element of windtunnel.py is the ability to analyze both continuous release (point concentration) and puff release (puff concentration) timeseries separately. It is up to the user to determine which type of tracer gas release was used for a particular time series; it is also critical to determine this **before** running this program. The two analysis modes cannot be used interchangeably under any circumstances! **If the incorrect analysis mode is selected, all functions will appear to run and output normally, but the output will make as much sense as replacing the 16 anvils in Richard Wagner’s Das Rheingold with 25 Hautboys!** 

### 3.1 Continuous Release

The continuous release mode implies that the tracer gas was released over the entire measuring time at a steady, continuous rate into the wind tunnel. Analyzing this type of data is comparatively straightforward and the program will usually run very quickly (i.e. not more than a few seconds) unless a large amount of input data is specified (>100 mb).  

#### 3.2 Puff Release

The puff release mode is quite a bit more complicated, but in many cases also more realistic for many applications. Here, instead of being released continuously over time, the tracer gas is released only during short (usually not more than a few seconds) but repeated intervals. Each of these intervals is called a puff. For the purpose of diagnosing the characteristics of the individual puffs, the program plots, among other things, the release of the tracer gas over time into the puff plots. Note that the release signal is assumed to be binary; i.e at a point in time the gas can either be released or not released. As it has been found that the release signal is usually quite noisy, PAPE also contains a straightforward noise removal algorithm.  

#### 3.2.1 Issues with Puff Detection

If you encounter issues with the detection of puffs (i.e. multiple pluffs showing up in one puff plot), you can try changing the value of variable ‘signal_mode.’ Note that the default value (if signal_mode is set to ‘None’) is effectively 1, in short, lowering this value will increase the likelihood of noise being interpreted as a puff, while raising this value will increase the likelihood that puffs will be interpreted as noise and remain undetected. This variable is located in the advanced options menu in the GUI. **Do not modify this variable unless you know what you are doing!** In the default setting, the program will work properly, as long as the magnitude of the noise is less than the mean signal, and the assumption that the signal is binary holds true.  

### 3.2.2 Release Signal in Puff Mode

The mean release signal (plotted with the mean puff) can give crucial information about the release length in each puff; if the release lengths of the individual puffs are all the same, the mean release signal appears as a rectangle function (as do all the individual release signals for a particular puff). If the puff release lengths are all different, the mean release function will have a non-rectangular shape. 

## 4. Contributing to windtunnel.py

When contributing code, be sure to adhere to the following checklist:

1. **Pull** current version of windtunnel.py to your **local repository**.
2. Make sure everything works fine in your local environment.
    1. if not, check the **installation requirements** because it should work fine.
3. Create a **local branch** named ...
    1. ... unmistakeably after the **functionality** you want to implement.
    2. ... after an **issue** created on gitlab you want to solve.
4. **Make your change.** 
5. Ensure yourself your change work fine. 
6. Create a **Merge Request** with target branch **master** on the remote repository. 
    1. The merge requests are reviewed and accepted by the maintainers/owners. 
    2. Contributors cannot accept their own merge requests. This also holds for the maintainers/owners.

By following this checklist, we seek to minimize the number of merge conflicts and errors. Every contributor is asked to follow the [PEP8 Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/). 

Please also ensure that your **.gitignore**-file contains at least the following:

```glsl
# ignore pyc-files
*.pyc
# ignore generated files (log, txt) 
windtunnel.log
postprocessed*
*txt
# ignore own scripts
own*
example*
# ignore itself
.gitignore
```


## 5. PAPE

PAPE (**P**oint **A**nd **P**uff Concentration M**E**asurements Analysis Program) is the native software implementation of the windtunnel.py package. It comes in two modes (GUI and CLI), and should be sufficient for most data analysis applicatons of the windtunnel.py package. However, PAPE does not exhaustively use all of the functionalities provided by the windtunnel.py package. 

### 5.1. Program structure
The program can easily be operated and started by running one of several source code scripts, but also features a graphical user interface (GUI) for easier operation of the program. Note that most of the development occurs primarily on the source code; thus, not all current features may be available when operating the program through the GUI. Note that it is substantially easier to operate the program incorrectly through the source code than through the GUI. **It is recommended that users with little or no experience with python programming use the GUI version the program.** 

#### 5.1.1: Graphical User Interface (GUI)

To run PAPE using the GUI, follow the steps above to allow for the execution of python scripts as executables, then double click on the file "PAPE_GUI.py" in the PAPE root folder. If a dialog appears asking you what to do with the script, select “Run” (exact wording may vary slightly depending on the OS). If you want to see what the program is doing and what it outputs you can also select “Run in Terminal.” This option is most likely to be suitable only for experienced users, but has much less potential to do damage to the actual program than running and editing the program from a CLI. 

In PAPE, the program options and a priori valuesn are in the main window towards the left. **The "default values" which are loaded upon starting the program come from an arbitrary dataset and are by no means universally applicable!** They are just placeholders to make sure that there is “something” in these values, and **should be edited to match the dataset being analyzed.** A detailed explanation of the physical meaning of the options and a priori values is provided in Sections 3 and 4. Note that the two options “Select input data type” correspond to switching between the two functions standard_point_analysis and standard_puff_analysis (see Sections 2.2 and 3). 
	
To start the analysis, you must first select one or multiple files for the program to analyze. In GUI mode this is accomplished by clicking on the “path” button at the top left corner. After selecting the data location, a selectable list of files which the program recognizes for analysis pops up on the right hand side of the main window. **If none of the files are selected, PAPE does nothing but stays open.** Currently, only files with the extension .txt.ts#0 are supported. If you would like support for more data, please contact the developers. 

PAPE also supports specifying the a priori information through a csv file; this option is recommended for larger datasets with multiple files. You can specify the location of the file containing the a priori information through the “CSV file (optional)” button next in the top left of the window. The csv file has several formatting requirements, so a sample csv file containing the a priori information is located in the /windtunnel subfolder of your installation directory (currently in the root directory of PAPE). The sample file is called Q2_Ambient_Conditions.csv. 

The a priori values can be reset to their default values by clicking on the “Reset Parameters” button at the top right of the window. 

PAPE also has two options for regulating the scope of the data analysis performed.
- ”full mode,” which runs the entire set of analysis methods, including ensemble and class analysis. This mode also outputs a full set of plots, including a separate plot for the concentration time series during each individual puffs. The number of puffs to plot can be changed in the GUI by setting the parameter “Number of Puffs to Plot” (variable n_puffs in the CLI script) to the desired number of puffs. It should be noted that changing the number of puffs to plot defeats the purpose of "full mode." **Running PAPE in “full mode” can be very time-consuming, especially for large datasets.** 
- ”basic mode,” which leaves out the class analysis, and consequently the histograms, entirely. The default setting for this mode is to plot only the first five puffs. However, this can be changed in the GUI by setting the parameter “Number of Puffs to Plot” to the desired number of puffs to plot (variable n_puffs in the CLI script). Setting this variable to “all” will plot all puffs, regardless of puff size. **Running PAPE in “basic mode” is substantially faster than “full mode,” and is recommended for most purposes.**

These options can be set directly from the GUI. PAPEs default setting is to run in basic mode. 

The “Advanced Options” button opens a new subwindow, which pulls up a range of additional settings for the data analysis. **Do not modify these settings unless you know what you are doing.** Incorrectly modifying these parameters can several compromise the integrity of the output of this program. If you want to dive into these more advanced parameters, it is highly recommended to get started reading the description of the variables time_threshol, n_outliers, axis_range, and signal mode in section 4. **Warning: changing the variable signal_mode can substantially compromise the data analysis quality! Do not modify this variable unless you know what you are doing!** 

#### 5.1.2 Command-Line Interface (CLI)
To run PAPE in the command-line interface mode (CLI), start by opening the file example_point_measurement.py (or example_puff_measurement.py for puff measurement, if you wish), either in spyder or using a text editor. Enter the path of the dataset (variable ‘path’) and filename(s) (variable ‘namelist’) of the dataset you wish to analyze. You should also specify the a priori information located below the namelist variable. A detailed explanation of the physical meaning of the options and a priori values is provided in Sections 3 and 4. While the program will in theory run without any a priori values, correct a priori information is required for physically correct results. Of particular importance are the variables x_source, y_source, z_source, x_measure, y_measure, z_measure, pressure, temperature, scaling factor, scale, ref_length, and ref_height. Furthermore, you should specify the outputs scale of the data, by setting the variable full_scale to either “fs” (full-scale output), “ms” (model scale output), or “nd” (non-dimensional output). The variable ‘axis_range’ can also be set to control the range of the y-axis in all puff plots; setting this variable to ‘auto’ will automatically slecte the y-axis range individually for each puff, while setting this variable to ‘same’ will force the y-axis range to be equal for all puffs being plotted. Now you can run whichever file you just edited (in spyder just click ‘run’ in the top toolbar), and the program should output the results of your analysis. 

Note that if you want to analyze more than one dataset (i.e. if you have more than one entry in namelist), you may want to consider inputting the a priori information into a csv file. The csv file has several formatting requirements, so a sample csv file containing the a priori information is located in the /windtunnel subfolder of your installation directory (currently //ewtl2/work/Johannes/PAPE/windtunnel). The sample file is called Q2_Ambient_Conditions.csv; if you rename this file you must also enter the new name into the python script (variable csv_file, line 18 in the point measurement script and line 21 in the puff measurement script). In the csv file, the first line must contain the filename of each dataset, and this must be identical to the filenames of your input datasets (see sample file). If this is not the case, the program will resort to using the a priori values written in the respective python script. Note further that you must enter (proper) values for all variables located in the sample csv file, you cannot selectively specify some variables from the csv file and others directly from the python file. If the csv file does not contain values for all variables located in the sample csv file, the program will once again resort to using the values from the python script for all the a priori values (thus ignoring the values in the csv file). Everything else is taken care of by the program.

In puff mode, PAPE also has two options for regulating the scope of the data analysis performed. These are controlled by the variable “functions_mode”
 - “full mode,” (functions_mode=’full’) which runs the entire set of analysis methods, including. ensemble and class analayis. This mode also outputs a full set of plots, including a separate plot for the concentration time series during each individual puffs. The number of puffs to plot can be changed in the GUI by setting the parameter “Number of Puffs to Plot” (variable n_puffs in the CLI script) to the desired number of puffs. It should be noted that changing the number of puffs to plot defeats the purpose of "full mode." **Running PAPE in "full mode" can be very time-consuming, especially for large datasets.**
 - “basic mode,” (functions_mode=’basic’) which leaves out the class analysis, and consequently the histograms, entirely. The default setting of this mode also plots only the first five puffs. However, this can be changed in the GUI by setting the parameter “Number of Puffs to Plot” to the desired number of puffs to plot (variable n_puffs in the CLI script). Setting this variable to “all” will plot all puffs, regardless of puff size. **Running PAPE in “basic mode” is substantially faster than "full mode," and is recommended for most purposes.** 

The number of puffs to be plotted can also be specified using the variable “n_puffs.” This variable was added by request from Simon Michel. Setting “n_puffs” to ‘default’ will use the default values based on the mode selected in “functions_mode,” while specifying a number for “n_puffs” overrides the default settings. 

In puffs mode, there also exists a variable called ‘signal_mode.’ **Do not modify this variable unless you know what you are doing!** Changing the input for this variable can substantially compromise the validity of the data analysis. If you want to dive into modifying the variable signal_mode, it is highly recommended to read over the description of this variable in section 4 first. 

Note that the output is, per default, saved in the same location as your original dataset (i.e. under the path you entered for the variable ‘path’). While the point measurement script outputs mainly text files, the puff measurement script outputs several (i.e. many) plots in several subfolder. All necessary subfolders are automatically generated, provided that they don’t already exist.  


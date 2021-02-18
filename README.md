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

After following these steps with the required dependencies, should be running on your computer. 

## 2. Contributing to windtunnel.py

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


## 3. PAPE

PAPE (**P**oint **A**nd **P**uff Concentration M**E**asurements Analysis Program) is the native software implementation of the windtunnel.py package. It comes in two modes (GUI and CLI), and should be sufficient for most data analysis applicatons of the windtunnel.py package. However, PAPE does not exhaustively use all of the functionalities provided by the windtunnel.py package. 

### 3.1. Program structure
The program can easily be operated and started by running one of several source code scripts, but also features a graphical user interface (GUI) for easier operation of the program. Note that most of the development occurs primarily on the source code; thus, not all current features may be available when operating the program through the GUI. Note that it is substantially easier to operate the program incorrectly through the source code than through the GUI. **It is recommended that users with little or no experience with python programming use the GUI version the program.** 

#### 3.1.1: Graphical User Interface (GUI)

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

#### 3.1.2 Command-Line Interface (CLI)
To run PAPE in the command-line interface mode (CLI), start by opening the file example_point_measurement.py (or example_puff_measurement.py for puff measurement, if you wish), either in spyder or using a text editor. Enter the path of the dataset (variable ‘path’) and filename(s) (variable ‘namelist’) of the dataset you wish to analyze. You should also specify the a priori information located below the namelist variable. A detailed explanation of the physical meaning of the options and a priori values is provided in Sections 3 and 4. While the program will in theory run without any a priori values, correct a priori information is required for physically correct results. Of particular importance are the variables x_source, y_source, z_source, x_measure, y_measure, z_measure, pressure, temperature, scaling factor, scale, ref_length, and ref_height. Furthermore, you should specify the outputs scale of the data, by setting the variable full_scale to either “fs” (full-scale output), “ms” (model scale output), or “nd” (non-dimensional output). The variable ‘axis_range’ can also be set to control the range of the y-axis in all puff plots; setting this variable to ‘auto’ will automatically slecte the y-axis range individually for each puff, while setting this variable to ‘same’ will force the y-axis range to be equal for all puffs being plotted. Now you can run whichever file you just edited (in spyder just click ‘run’ in the top toolbar), and the program should output the results of your analysis. 

Note that if you want to analyze more than one dataset (i.e. if you have more than one entry in namelist), you may want to consider inputting the a priori information into a csv file. The csv file has several formatting requirements, so a sample csv file containing the a priori information is located in the /windtunnel subfolder of your installation directory (currently //ewtl2/work/Johannes/PAPE/windtunnel). The sample file is called Q2_Ambient_Conditions.csv; if you rename this file you must also enter the new name into the python script (variable csv_file, line 18 in the point measurement script and line 21 in the puff measurement script). In the csv file, the first line must contain the filename of each dataset, and this must be identical to the filenames of your input datasets (see sample file). If this is not the case, the program will resort to using the a priori values written in the respective python script. Note further that you must enter (proper) values for all variables located in the sample csv file, you cannot selectively specify some variables from the csv file and others directly from the python file. If the csv file does not contain values for all variables located in the sample csv file, the program will once again resort to using the values from the python script for all the a priori values (thus ignoring the values in the csv file). Everything else is taken care of by the program.

In puff mode, PAPE also has two options for regulating the scope of the data analysis performed. These are controlled by the variable “functions_mode”
 - “full mode,” (functions_mode=’full’) which runs the entire set of analysis methods, including. ensemble and class analayis. This mode also outputs a full set of plots, including a separate plot for the concentration time series during each individual puffs. The number of puffs to plot can be changed in the GUI by setting the parameter “Number of Puffs to Plot” (variable n_puffs in the CLI script) to the desired number of puffs. It should be noted that changing the number of puffs to plot defeats the purpose of "full mode." **Running PAPE in "full mode" can be very time-consuming, especially for large datasets.**
 - “basic mode,” (functions_mode=’basic’) which leaves out the class analysis, and consequently the histograms, entirely. The default setting of this mode also plots only the first five puffs. However, this can be changed in the GUI by setting the parameter “Number of Puffs to Plot” to the desired number of puffs to plot (variable n_puffs in the CLI script). Setting this variable to “all” will plot all puffs, regardless of puff size. **Running PAPE in “basic mode” is substantially faster than "full mode," and is recommended for most purposes.** 

The number of puffs to be plotted can also be specified using the variable “n_puffs.” This variable was added by request from Simon Michel. Setting “n_puffs” to ‘default’ will use the default values based on the mode selected in “functions_mode,” while specifying a number for “n_puffs” overrides the default settings. 

In puffs mode, there also exists a variable called ‘signal_mode.’ **Do not modify this variable unless you know what you are doing!** Changing the input for this variable can substantially compromise the validity of the data analysis. If you want to dive into modifying the variable signal_mode, it is highly recommended to read over the description of this variable in section 4 first. 

Note that the output is, per default, saved in the same location as your original dataset (i.e. under the path you entered for the variable ‘path’). While the point measurement script outputs mainly text files, the puff measurement script outputs several (i.e. many) plots in several subfolder. All necessary subfolders are automatically generated, provided that they don’t already exist.  


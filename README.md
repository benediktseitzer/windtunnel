# Windtunnel 
The windtunnel python module is a processing library developed by the EWTL Group at the University of Hamburg. 

## 1. Installation
To set up the windtunnel-package, simply place it in a directory of your choice. Executing scripts can, but must not be placed in the same directory as the windtunnel module. 

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

#### 1.2. Windows
Use the terminal or an Anaconda Prompt for the following steps:

    Create the environment from the environment.yml file:

    conda env create -f environment.yml

    The first line of the yml file sets the new environment's name. For details see Creating an environment file manually.

    Activate the new environment: conda activate myenv

    Verify that the new environment was installed correctly:

    conda env list

    You can also use conda info --envs.

For more information follow:
https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file    
#### 1.3. Linux and MacOs
You can pip install all the packages by simply typing:
```sh
$ pip3 install numpy pandas scipy matplotlib tsp_solver
$ sudo apt-get install python3-tk
```

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

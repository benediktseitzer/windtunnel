# Windtunnel 
The windtunnel python module is a processing library developed by the EWTL Laboratory at the University of Hamburg. 

## Installation
#### Requirements
 - python 3.0 
 - the following python libraries:
 -- numpy, pandas, os, sys, time, fnmatch, logging, math, scipy, matplotlib, csv, tsp_solver, mpl_toolkits, tkinter

#### Windows
somebody else can write this but it's most likely the same as below.

#### Linux and MacOs
You can pip install all the packages by typing:
```sh
$ pip3 install numpy pandas os sys time fnmatch logging math scipy matplotlib csv tsp_solver mpl_toolkits tkinter 
```

## Program structure
test

## How to work with windtunnel.py
test test

## Contributing to windtunnel.py
When contributing code, you’ll want to follow this checklist:
1. **Pull** current version of windtunnel.py to your **local repository**.
2. Make sure everything works fine in your local environment.
2.1 if not, check the **installation requirements** because it should work fine.
3. Create a **local branch** named ...
3.1 ... unmistakeably after the **functionality** you want to implement.
3.2 ... after an **issue** created on gitlab you want to solve.
4. **Make your change.** 
5. Ensure yourself your change work fine. 
6. Create a **Merge Request** with target branch **master** on the remote repository. 
6.1 The merge requests are reviewed and accepted by the maintainers/owners. 
6.2 Contributors cannot accept their own merge requests. This also holds for the maintainers/owners.

By following this checklist, we seek to minimize the number of merge conflicts and errors. Every contributor is asked to follow the [PEP8 Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/). Please also ensure that your **.gitignore**-file contains at least the following:
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

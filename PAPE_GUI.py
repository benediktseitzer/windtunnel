#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 17:04:48 2020

@author: johannes
"""
import tkinter as tk
import tkinter.filedialog
import tkinter.simpledialog
import tkinter.messagebox
import windtunnel as wt
import pandas as pd
import numpy as np
import os
import sys

class PAPE(tk.Frame):
    def __init__(self,master=None):
        super().__init__(master)
        self.master = master
        self.grid()
        self.create_widgets()
        self.path=None
        self.csv_file=None  
        self.namelist=[]
        self.var={}
        self.ambient_conditions={}
    def create_widgets(self):
        self.threshold_concentration=tk.DoubleVar()
        self.threshold_dosage=tk.DoubleVar()        
        self.n_exclude=tk.StringVar()  
        self.axis_range=tk.StringVar() 
        self.time_threshold=tk.DoubleVar()
        self.x_source=tk.DoubleVar() 
        self.y_source=tk.DoubleVar() 
        self.z_source=tk.DoubleVar()
        self.x_measure=tk.DoubleVar() 
        self.y_measure=tk.DoubleVar()
        self.z_measure=tk.DoubleVar() 
        self.pressure=tk.DoubleVar()  
        self.temperature=tk.DoubleVar()
        self.calibration_curve=tk.DoubleVar() 
        self.mass_flow_controller=tk.StringVar() 
        self.calibration_factor=tk.DoubleVar()
        self.scaling_factor=tk.DoubleVar()
        self.scale=tk.DoubleVar()
        self.ref_length=tk.DoubleVar()
        self.ref_height=tk.StringVar()
        self.gas_name=tk.StringVar()
        self.mol_weight=tk.DoubleVar()
        self.gas_factor=tk.DoubleVar()
        self.full_scale_wtref=tk.DoubleVar()
        self.full_scale_flow_rate=tk.DoubleVar()
        self.full_scale=tk.StringVar() 
        self.functions_mode=tk.StringVar() 
        self.dist=tk.IntVar()         
        self.puff=tk.StringVar()                             
        self.wdir=tk.DoubleVar()   
        self.f1=tk.Frame(self)
        self.f1.grid(row=0,column=0,columnspan=4,sticky='WE')                  
        self.path_button = tk.Button(self.f1)
        self.path_button["text"] = "Path"  
        self.path_button["command"] = self.open_path      
        self.path_button.grid(row=0,column=0)        
        self.csv_button = tk.Button(self.f1)
        self.csv_button["text"] = "CSV file (optional)"  
        self.csv_button["command"] = self.open_csv      
        self.csv_button.grid(row=0,column=1)
        self.n_exclude.set('Auto-Select')
        self.axis_range.set('auto') 
        self.time_threshold.set(5)          
        self.advanced_button = tk.Button(self.f1)
        self.advanced_button["text"] = "Advanced Options"  
        self.advanced_button["command"] = self.advanced_options    
        self.advanced_button.grid(row=0,column=2)  
        self.puff.set("puff")        
        self.f4=tk.Frame(self)
        self.f4.grid(row=1,column=1,columnspan=10,sticky='WE')        
        self.puff_button_label=tk.Label(self)
        self.puff_button_label["text"]="Select Input Data Type"
        self.puff_button_label.grid(row=1,column=0)
        self.puff_button_1=tk.Radiobutton(self.f4)
        self.puff_button_1["text"]="Puff Measurements"
        self.puff_button_1["variable"]=self.puff
        self.puff_button_1["value"]="puff" 
        self.puff_button_1["command"]=self.input_puff  
        self.puff_button_1.grid(row=0,column=0)  
        self.puff_button_2=tk.Radiobutton(self.f4)
        self.puff_button_2["text"]="Point Measurements"
        self.puff_button_2["variable"]=self.puff
        self.puff_button_2["value"]="point"   
        self.puff_button_2["command"]=self.input_point
        self.puff_button_2.grid(row=0,column=1)            
        self.threshold_concentration.set(0)
        self.e1=tk.Entry(self,textvariable=self.threshold_concentration)
        self.e1.grid(row=2,column=1,sticky='WE')
        self.e1_label=tk.Label(self)
        self.e1_label["text"]="Threshold Concentration"
        self.e1_label.grid(row=2,column=0)  
        self.threshold_dosage.set(0)        
        self.e2=tk.Entry(self,textvariable=self.threshold_dosage)
        self.e2.grid(row=3,column=1,sticky='WE')
        self.e2_label=tk.Label(self)
        self.e2_label["text"]="Threshold Dosage"        
        self.e2_label.grid(row=3,column=0)
        self.n_exclude.set('Auto-Select')    
        self.axis_range.set('auto')           
        self.time_threshold.set(5)
        self.x_source.set(0)
        self.y_source.set(0)
        self.z_source.set(0)           
        self.e5_label=tk.Label(self)
        self.e5_label["text"]="Source Location (x,y,z):"        
        self.e5_label.grid(row=4,column=0)         
        self.e5_1=tk.Entry(self,textvariable=self.x_source)
        self.e5_1.grid(row=4,column=1,sticky='WE')
        self.e5_2=tk.Entry(self,textvariable=self.y_source)
        self.e5_2.grid(row=5,column=1,sticky='WE') 
        self.e5_3=tk.Entry(self,textvariable=self.z_source)
        self.e5_3.grid(row=6,column=1,sticky='WE') 
        self.x_measure.set(855.16)
        self.y_measure.set(176.29)
        self.z_measure.set(162)
        self.e6_label=tk.Label(self)
        self.e6_label["text"]="Measurement Location (x,y,z):"        
        self.e6_label.grid(row=7,column=0)         
        self.e6_1=tk.Entry(self,textvariable=self.x_measure)
        self.e6_1.grid(row=7,column=1,sticky='WE')
        self.e6_2=tk.Entry(self,textvariable=self.y_measure)
        self.e6_2.grid(row=8,column=1,sticky='WE') 
        self.e6_3=tk.Entry(self,textvariable=self.z_measure)
        self.e6_3.grid(row=9,column=1,sticky='WE') 
        self.pressure.set(1009)
        self.e7=tk.Entry(self,textvariable=self.pressure)
        self.e7.grid(row=10,column=1,sticky='WE')
        self.e7_label=tk.Label(self)
        self.e7_label["text"]="Atmospheric Pressure"        
        self.e7_label.grid(row=10,column=0)  
        self.e7_label2=tk.Label(self)
        self.e7_label2["text"]="hPa"        
        self.e7_label2.grid(row=10,column=2,sticky='W')
        self.temperature.set(23.5)        
        self.e8=tk.Entry(self,textvariable=self.temperature)
        self.e8.grid(row=11,column=1,sticky='WE')
        self.e8_label=tk.Label(self)
        self.e8_label["text"]="Temperature"        
        self.e8_label.grid(row=11,column=0)  
        self.e8_label2=tk.Label(self)
        self.e8_label2["text"]="°C"        
        self.e8_label2.grid(row=11,column=2,sticky='W')
        self.wdir.set(0)        
        self.e9=tk.Entry(self,textvariable=self.wdir)
        self.e9.grid(row=12,column=1,sticky='WE')
        self.e9_label=tk.Label(self)
        self.e9_label["text"]="Wind Direction"        
        self.e9_label.grid(row=11,column=0,)  
        self.e9_label2=tk.Label(self)
        self.e9_label2["text"]="°"        
        self.e9_label2.grid(row=12,column=2,sticky='W') 
        self.calibration_curve.set(0.3)        
        self.e10=tk.Entry(self,textvariable=self.calibration_curve)
        self.e10.grid(row=13,column=1,sticky='WE')
        self.e10_label=tk.Label(self)
        self.e10_label["text"]="Calibration curve"        
        self.e10_label.grid(row=13,column=0)
        self.mass_flow_controller.set('X')
        self.e11=tk.Entry(self,textvariable=self.mass_flow_controller)
        self.e11.grid(row=14,column=1,sticky='WE')
        self.e11_label=tk.Label(self)
        self.e11_label["text"]="Mass flow controller"        
        self.e11_label.grid(row=14,column=0)
        self.calibration_factor.set(0.637)
        self.e12=tk.Entry(self,textvariable=self.calibration_factor)
        self.e12.grid(row=15,column=1,sticky='WE')
        self.e12_label=tk.Label(self)
        self.e12_label["text"]="Calibration factor"        
        self.e12_label.grid(row=15,column=0)
        self.scaling_factor.set(0.637)
        self.e13=tk.Entry(self,textvariable=self.scaling_factor)
        self.e13.grid(row=16,column=1,sticky='WE')
        self.e13_label=tk.Label(self)
        self.e13_label["text"]="Scaling factor"        
        self.e13_label.grid(row=16,column=0)
        self.scale.set(250)
        self.e14=tk.Entry(self,textvariable=self.scale)
        self.e14.grid(row=17,column=1,sticky='WE')
        self.e14_label=tk.Label(self)
        self.e14_label["text"]="Scale"        
        self.e14_label.grid(row=17,column=0)
        self.e14_label2=tk.Label(self)
        self.e14_label2["text"]="1:"        
        self.e14_label2.grid(row=17,column=0,stick='E')
        self.ref_length.set(1/250)        
        self.e15=tk.Entry(self,textvariable=self.ref_length)
        self.e15.grid(row=18,column=1,sticky='WE')
        self.e15_label=tk.Label(self)
        self.e15_label["text"]="Reference Length"        
        self.e15_label.grid(row=18,column=0) 
        self.ref_height.set(None)
        self.e16=tk.Entry(self,textvariable=self.ref_height)
        self.e16.grid(row=19,column=1,sticky='WE')
        self.e16_label=tk.Label(self)
        self.e16_label["text"]="Reference Height"        
        self.e16_label.grid(row=19,column=0)
        self.gas_name.set('C12')
        self.e17=tk.Entry(self,textvariable=self.gas_name)
        self.e17.grid(row=20,column=1,sticky='WE')
        self.e17_label=tk.Label(self)
        self.e17_label["text"]="Gas Name"        
        self.e17_label.grid(row=20,column=0)
        self.mol_weight.set(28.97)
        self.e18=tk.Entry(self,textvariable=self.mol_weight)
        self.e18.grid(row=21,column=1,sticky='WE')
        self.e18_label=tk.Label(self)
        self.e18_label["text"]="Molecular Mass"        
        self.e18_label.grid(row=21,column=0)
        self.e18_label2=tk.Label(self)
        self.e18_label2["text"]="g/mol"        
        self.e18_label2.grid(row=21,column=2,sticky='W')
        self.gas_factor.set(0.5)        
        self.e19=tk.Entry(self,textvariable=self.gas_factor)
        self.e19.grid(row=22,column=1,sticky='WE')
        self.e19_label=tk.Label(self)
        self.e19_label["text"]="Gas Factor"        
        self.e19_label.grid(row=22,column=0)
        self.full_scale_wtref.set(6)
        self.e20=tk.Entry(self,textvariable=self.full_scale_wtref)
        self.e20.grid(row=23,column=1,sticky='WE')
        self.e20_label=tk.Label(self)
        self.e20_label["text"]="Full Scale Referece Wind Speed"        
        self.e20_label.grid(row=23,column=0)
        self.e20_label2=tk.Label(self)
        self.e20_label2["text"]="m/s"        
        self.e20_label2.grid(row=23,column=2,sticky='W')
        self.full_scale_flow_rate.set(0.5)        
        self.e21=tk.Entry(self,textvariable=self.full_scale_flow_rate)
        self.e21.grid(row=24,column=1,sticky='WE')
        self.e21_label=tk.Label(self)
        self.e21_label["text"]="Full Scale Mass Flow Rate"        
        self.e21_label.grid(row=24,column=0)
        self.e21_label2=tk.Label(self)
        self.e21_label2["text"]="kg/s"        
        self.e21_label2.grid(row=24,column=2,sticky='W')
        self.full_scale.set("fs")
        self.f2=tk.Frame(self)
        self.f2.grid(row=25,column=1,columnspan=3,sticky='WE')
        self.scale_button_label=tk.Label(self)
        self.scale_button_label["text"]="Select Scale of Data Output"
        self.scale_button_label.grid(row=25,column=0)
        self.scale_button_1=tk.Radiobutton(self.f2)
        self.scale_button_1["text"]="Full Scale"
        self.scale_button_1["variable"]=self.full_scale
        self.scale_button_1["value"]="fs"       
        self.scale_button_1.grid(row=0,column=0)  
        self.scale_button_2=tk.Radiobutton(self.f2)
        self.scale_button_2["text"]="Model Scale"
        self.scale_button_2["variable"]=self.full_scale
        self.scale_button_2["value"]="ms"              
        self.scale_button_2.grid(row=0,column=1)        
        self.scale_button_3=tk.Radiobutton(self.f2)
        self.scale_button_3["text"]="Non-Dimensional"
        self.scale_button_3["variable"]=self.full_scale
        self.scale_button_3["value"]="nd"             
        self.scale_button_3.grid(row=0,column=2)
        self.functions_mode.set("basic")        
        self.f3=tk.Frame(self)
        self.f3.grid(row=26,column=1,columnspan=10,sticky='WE')
        self.mode_button_label=tk.Label(self)
        self.mode_button_label["text"]="Options:"
        self.mode_button_label.grid(row=26,column=0)
        self.mode_button_1=tk.Radiobutton(self.f3)
        self.mode_button_1["text"]="Full Mode: all functions and plots, \n but longer running time"
        self.mode_button_1["variable"]=self.functions_mode
        self.mode_button_1["value"]="full"       
        self.mode_button_1.grid(row=0,column=0)  
        self.mode_button_2=tk.Radiobutton(self.f3)
        self.mode_button_2["text"]="Basic Mode: no histograms or class analysis, and plot only \n first 5 puffs. Substantially faster than full mode."
        self.mode_button_2["variable"]=self.functions_mode
        self.mode_button_2["value"]="basic" 
        self.mode_button_2.grid(row=0,column=1) 
        self.reset = tk.Button(self.f1, text="Reset Parameters", command=self.reset_input)
        self.reset.grid(row=0,column=3)   
        self.dist.set(0)
        self.dist_button=tk.Checkbutton(self)
        self.dist_button["text"]="Show advanced distribution data with mean puff"
        self.dist_button["variable"]=self.dist
        self.dist_button.grid(row=27,column=1,sticky="W")                         
                  
        self.quit = tk.Button(self.f1, text="QUIT", fg="red", command=self.quit_pape)
        self.quit.grid(row=0,column=4)
        self.run_button = tk.Button(self, text="RUN", fg="red", command=self.run_script)
        self.run_button.grid(row=28,column=1)  
 
         
    def open_path(self):
        self.path=tkinter.filedialog.askdirectory()
        #self.namelist=os.scandir(path=self.path+"/*.txt.ts#0")
        for file in os.listdir(path=self.path):
            if file.endswith(".txt.ts#0"):
                self.namelist.append(file)
        self.namelist.sort()        
        for i in range(len(self.namelist)):
            self.var[i]=tk.IntVar()
            self.var[i].set(0)
            self.file_button=tk.Checkbutton(self)
            self.file_button["text"]=self.namelist[i]
            self.file_button["variable"]=self.var[i]           
            self.file_button.grid(row=i+3,column=3)
            self.path_label=tk.Label(self)
            self.path_label["text"]="Select Datasets to Analyze:"        
            self.path_label.grid(row=2,column=3)
    def advanced_options(self):
        dialog1=advanced_dialog(root)

                      
    def printSelection(self, i):
        print(self.var[i].get())            
    def open_csv(self):
        self.csv_file=tkinter.filedialog.askopenfile()	
        self.ambient_conditions=pd.read_csv(self.csv_file,sep=',',index_col=0) 
        necessary_keys={'x_source','y_source','z_source','x_measure','y_measure','z_measure','pressure','temperature','wdir','calibration_curve','mass_flow_controller','calibration_factor', \
        'scaling_factor','scale','ref_length','ref_height','gas_name','mol_weight','gas_factor','full_scale_wtref','full_scale_flow_rate' }
        name=self.ambient_conditions.keys()[0]        
        if not all(name2 in self.ambient_conditions[name] for name2 in necessary_keys):
           print('Error: csv file does not contain all necessary ambient conditions data. Check to make sure that csv file to make sure that \
the csv file contains all necessary data and is properly formatted. Resorting to input data in example_puff_measurement.py')
           return
        self.x_source.set(None) if self.ambient_conditions[name]['x_source'] =='None' else self.x_source.set(np.float(self.ambient_conditions[name]['x_source']))
        self.y_source.set(None) if self.ambient_conditions[name]['y_source'] =='None' else self.y_source.set(np.float(self.ambient_conditions[name]['y_source']))
        self.z_source.set(None) if self.ambient_conditions[name]['z_source'] =='None' else self.z_source.set(np.float(self.ambient_conditions[name]['z_source']))  
        self.x_measure.set(None) if self.ambient_conditions[name]['x_measure'] =='None' else self.x_measure.set(np.float(self.ambient_conditions[name]['x_measure']))
        self.y_measure.set(None) if self.ambient_conditions[name]['y_measure'] =='None' else self.y_measure.set(np.float(self.ambient_conditions[name]['y_measure']))
        self.z_measure.set(None) if self.ambient_conditions[name]['z_measure'] =='None' else self.z_measure.set(np.float(self.ambient_conditions[name]['z_measure']))          
        self.pressure.set(None) if self.ambient_conditions[name]['pressure'] =='None' else self.pressure.set(np.float(self.ambient_conditions[name]['pressure']))		
        self.temperature.set(None) if self.ambient_conditions[name]['temperature'] =='None' else self.temperature.set(np.float(self.ambient_conditions[name]['temperature']))
        self.wdir.set(None) if self.ambient_conditions[name]['wdir'] =='None' else self.wdir.set(np.float(self.ambient_conditions[name]['wdir']))	
        self.calibration_curve.set(None) if self.ambient_conditions[name]['calibration_curve'] =='None' else self.calibration_curve.set(np.float(self.ambient_conditions[name]['calibration_curve']))
        self.mass_flow_controller.set(None) if self.ambient_conditions[name]['mass_flow_controller'] =='None' else self.mass_flow_controller.set(self.ambient_conditions[name]['mass_flow_controller'])
        self.calibration_factor.set(None) if self.ambient_conditions[name]['calibration_factor'] =='None' else self.calibration_factor.set(np.float(self.ambient_conditions[name]['calibration_factor']))
        self.scaling_factor.set(None) if self.ambient_conditions[name]['scaling_factor'] =='None' else self.scaling_factor.set(np.float(self.ambient_conditions[name]['scaling_factor']))	
        self.scale.set(None) if self.ambient_conditions[name]['scale'] =='None' else self.scale.set(np.float(self.ambient_conditions[name]['scale']))
        self.ref_length.set(None) if self.ambient_conditions[name]['ref_length'] =='None' else self.ref_length.set(np.float(eval(self.ambient_conditions[name]['ref_length'])))
        self.ref_height.set(None) if self.ambient_conditions[name]['ref_height'] =='None' else self.ref_height.set(np.float(self.ambient_conditions[name]['ref_height']))	
        self.gas_name.set(None) if self.ambient_conditions[name]['gas_name'] =='None' else self.gas_name.set(self.ambient_conditions[name]['gas_name'])
        self.mol_weight.set(None) if self.ambient_conditions[name]['mol_weight'] =='None' else self.mol_weight.set(np.float(self.ambient_conditions[name]['mol_weight']))
        self.gas_factor.set(None) if self.ambient_conditions[name]['gas_factor'] =='None' else self.gas_factor.set(np.float(self.ambient_conditions[name]['gas_factor']))
        self.full_scale_wtref.set(None) if self.ambient_conditions[name]['full_scale_wtref'] =='None' else self.full_scale_wtref.set(np.float(self.ambient_conditions[name]['full_scale_wtref']))
        self.full_scale_flow_rate.set(None) if self.ambient_conditions[name]['full_scale_flow_rate'] =='None' else self.full_scale_flow_rate.set(np.float(self.ambient_conditions[name]['full_scale_flow_rate']))
        self.e1["state"]="disabled"
        self.e2["state"]="disabled"   
        #self.e3["state"]="disabled"
        #self.e4["state"]="disabled"  
        self.e5_1["state"]="disabled"
        self.e5_2["state"]="disabled"   
        self.e5_3["state"]="disabled"     
        self.e6_1["state"]="disabled"
        self.e6_2["state"]="disabled"   
        self.e6_3["state"]="disabled"
        self.e7["state"]="disabled"
        self.e8["state"]="disabled"   
        self.e9["state"]="disabled"
        self.e10["state"]="disabled"   
        self.e11["state"]="disabled"
        self.e12["state"]="disabled"   
        self.e13["state"]="disabled"
        self.e14["state"]="disabled"   
        self.e15["state"]="disabled"   
        self.e16["state"]="disabled"
        self.e17["state"]="disabled"   
        self.e18["state"]="disabled"
        self.e19["state"]="disabled"     
        self.e20["state"]="disabled"   
        self.e21["state"]="disabled"               
    def run_script(self):    
        
        namelist=np.array(self.namelist)
        functions_mode=self.functions_mode.get()
        path = self.path+'/'
        csv_file = None if self.csv_file==None else os.path.split(self.csv_file.name)[1]
        threshold_concentration=self.threshold_concentration.get()
        threshold_dosage=self.threshold_dosage.get() 
        #n_exclude=self.n_exclude.get()
        #print(self.n_exclude.get())
        n_exclude=None if self.n_exclude.get()=='Auto-Select' else np.float(self.n_exclude.get())         
        axis_range=self.axis_range.get()  
        time_threshold=self.time_threshold.get()/100        
        full_scale=self.full_scale.get() 
        x_source=self.x_source.get()
        y_source=self.y_source.get()
        z_source=self.z_source.get()
        x_measure=self.x_measure.get()
        y_measure=self.y_measure.get()
        z_measure=self.z_measure.get()
        pressure=self.pressure.get()
        temperature=self.temperature.get()
        wdir=self.wdir.get()
        calibration_curve=self.calibration_curve.get()
        mass_flow_controller=self.mass_flow_controller.get()
        calibration_factor=self.calibration_factor.get()
        scaling_factor=self.scaling_factor.get()
        scale=self.scale.get()
        ref_length=self.ref_length.get()
        ref_height=None if self.ref_height.get()=='None' else np.float(self.ref_height.get())
        gas_name=self.gas_name.get()
        mol_weight=self.mol_weight.get()
        gas_factor=self.gas_factor.get()
        full_scale_wtref=self.full_scale_wtref.get()
        full_scale_flow_rate=self.full_scale_flow_rate.get()   
        dist=self.dist.get()
        var=np.zeros(len(self.var),dtype=np.int)
        for i in range(len(self.var)):
            var[i]=self.var[i].get()
        namelist=namelist[var==1]
        if len(namelist)==0:
            tkinter.messagebox.showerror("Error","No file selected! Please select at least one file from the list to analyze.")
            return
        if self.puff.get()=="puff":
            wt.standard_puff_analysis(path,csv_file,namelist,threshold_concentration,threshold_dosage,
            n_exclude,time_threshold,full_scale,x_source,y_source,z_source,x_measure,y_measure,z_measure,
            pressure,temperature,wdir,calibration_curve,mass_flow_controller,calibration_factor,
            scaling_factor,scale,ref_length,ref_height,gas_name,mol_weight,gas_factor,full_scale_wtref,
            full_scale_flow_rate,functions_mode,axis_range) 
        elif self.puff.get()=="point":
            wt.standard_point_analysis(path,csv_file,namelist,full_scale,x_source,y_source,z_source,x_measure,y_measure,z_measure,
            pressure,temperature,wdir,calibration_curve,mass_flow_controller,calibration_factor,
            scaling_factor,scale,ref_length,ref_height,gas_name,mol_weight,gas_factor,full_scale_wtref, 
            full_scale_flow_rate,axis_range)
        
    def quit_pape(self):
        self.master.destroy()
        sys.exit
    def reset_input(self):
        self.threshold_concentration.set(0)
        self.threshold_dosage.set(0)             
        self.time_threshold.set(5)
        self.x_source.set(0)
        self.y_source.set(0)
        self.z_source.set(0)           
        self.x_measure.set(855.16)
        self.y_measure.set(176.29)
        self.z_measure.set(162)
        self.pressure.set(1009)
        self.temperature.set(23.5)        
        self.wdir.set(0)        
        self.calibration_curve.set(0.3)        
        self.mass_flow_controller.set('X')
        self.calibration_factor.set(0.637)
        self.scaling_factor.set(0.637)
        self.scale.set(250)
        self.ref_length.set(1/250)        
        self.ref_height.set(None)
        self.gas_name.set('C12')
        self.mol_weight.set(28.97)
        self.gas_factor.set(0.5)        
        self.full_scale_wtref.set(6)
        self.full_scale_flow_rate.set(0.5) 
        self.e1["state"]="normal"
        self.e2["state"]="normal"   
        #self.e3["state"]="normal"
        #self.e4["state"]="normal"  
        self.e5_1["state"]="normal"
        self.e5_2["state"]="normal"   
        self.e5_3["state"]="normal"     
        self.e6_1["state"]="normal"
        self.e6_2["state"]="normal"   
        self.e6_3["state"]="normal"
        self.e7["state"]="normal"
        self.e8["state"]="normal"   
        self.e9["state"]="normal"
        self.e10["state"]="normal"   
        self.e11["state"]="normal"
        self.e12["state"]="normal"   
        self.e13["state"]="normal"
        self.e14["state"]="normal"   
        self.e15["state"]="normal"   
        self.e16["state"]="normal"
        self.e17["state"]="normal"   
        self.e18["state"]="normal"
        self.e19["state"]="normal"     
        self.e20["state"]="normal"   
        self.e21["state"]="normal"          

    def input_puff(self):
        self.mode_button_1["state"]="normal"
        self.mode_button_2["state"]="normal"   
        self.dist_button["state"]="normal"
    def input_point(self):
        self.mode_button_1["state"]="disabled"
        self.mode_button_2["state"]="disabled" 
        self.dist_button["state"]="disabled"        
        
                                    
              
        
        

        #exec(open('/home/johannes/Desktop/Uni_Hamburg/Home_Office/Konzentrationen/PAPE_GUI_code.py').read(),globals(),locals())
           
class advanced_dialog(tkinter.simpledialog.Dialog):
    def __init__(self, master=None, title = None):

        tk.Toplevel.__init__(self, master)
        self.transient(master)

        if title:
            self.title(title)
            

     
        self.master = master
        
        self.create_widgets()

#        self.result = None
#
#        body = tk.Frame(self)
#        self.body(body)
#        body.pack(padx=5, pady=5)
#
        self.buttonbox()
    
#

#        #self.grab_set()
#
#       if not self.initial_focus:
#           self.initial_focus = self
#
#       self.protocol("WM_DELETE_WINDOW", self.cancel)
#
#       self.geometry("+%d+%d" % (self.master.winfo_rootx()+50,
#                                 self.master.winfo_rooty()+50))
#
#        self.initial_focus.focus_set()

        #self.wait_window(self)
    

    #
    # construction hooks    
    def create_widgets(self): 
        self.n_exclude=tk.StringVar()         
        self.n_exclude.set(app.n_exclude.get())
        self.axis_range=tk.StringVar()         
        self.axis_range.set(app.axis_range.get()) 
        self.time_threshold=tk.DoubleVar()        
        self.time_threshold.set(app.time_threshold.get())    
        self.e3=tk.Entry(self,textvariable=self.n_exclude)
        self.e3.grid(row=4,column=1,sticky='WE')
        self.e3_label=tk.Label(self)
        self.e3_label["text"]="Number of outliers to exclude"        
        self.e3_label.grid(row=4,column=0)    
        self.e4=tk.Entry(self,textvariable=self.time_threshold)
        self.e4.grid(row=5,column=1,sticky='WE')
        self.e4_label=tk.Label(self)
        self.e4_label["text"]="Percentage to use as threshold for \n computing arrival and leaving times"        
        self.e4_label.grid(row=5,column=0)  
        self.e4_label2=tk.Label(self)
        self.e4_label2["text"]="%"        
        self.e4_label2.grid(row=5,column=2,sticky='W')#
        self.f5=tk.Frame(self)
        self.f5.grid(row=6,column=1,columnspan=10,sticky='WE')        
        self.axis_options_button_1=tk.Radiobutton(self.f5)
        self.axis_options_button_1["text"]="Automatic"
        self.axis_options_button_1["variable"]=self.axis_range
        self.axis_options_button_1["value"]="auto" 
        self.axis_options_button_1.grid(row=6,column=1)         
        self.axis_options_button_2=tk.Radiobutton(self.f5)
        self.axis_options_button_2["text"]="Same"
        self.axis_options_button_2["variable"]=self.axis_range
        self.axis_options_button_2["value"]="same"     
        self.axis_options_button_2.grid(row=1,column=2)              
        self.axis_options_button_1.grid(row=1,column=1,sticky='WE')
        self.axis_options_button_1_label=tk.Label(self)
        self.axis_options_button_1_label["text"]="y-Axis Range:"        
        self.axis_options_button_1_label.grid(row=6,column=0)  

        
    def buttonbox(self):
        '''add standard button box.

        override if you do not want the standard buttons
        '''

        box = tk.Frame(self)

        w = tk.Button(box, text="OK", width=10, command=self.ok, default="active")
        w.pack(side="left", padx=5, pady=5)
        w = tk.Button(box, text="Cancel", width=10, command=self.destroy)
        w.pack(side="left", padx=5, pady=5)

        #self.bind("<Return>", self.ok)
        #self.bind("<Escape>", self.destroy)
        box.grid(column=5)
        
    def ok(self, event=None):

        if not self.validate():
            self.initial_focus.focus_set() # put focus back
            return

        self.withdraw()
        self.update_idletasks()

        try:
            self.apply()
        finally:
            self.destroy()        
        
    def apply(self):       
        self.n_exclude.set(self.e3.get())
        self.time_threshold.set(self.e4.get())   
        app.n_exclude.set(self.n_exclude.get())
        app.axis_range.set(self.axis_range.get())  
        app.time_threshold.set(self.time_threshold.get())

root = tk.Tk()
app = PAPE(master=root)
app.mainloop()
app.wait_window()



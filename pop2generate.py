#author Vishwas Verma
# date Oct 08,2016

import numpy as np
import scipy as sc
from scipy import integrate
import math
import random
import csv
import xlwt
import os
from functions2 import *

# an emperical relation is used to calculate deviation angle
class blade_section:
	
	def __init__(self,r_val):
		self.r_value = r_val
		self.alpha1 = 0
		self.beta1 = 0
		self.alpha2 = 0
		self.beta2 = 0
		self.Rx = 0
		self.C1 = 0
		self.C2 = 0
		self.V1 = 0
		self.V2 = 0
		self.de_haller_number = 0.72
		self.Cw1 = 0
		self.Cw2 = 0
		self.Ca1 = 0
		self.Ca2 = 0
		self.U = 0
		self.solidity_factor = 0   #c/s
		self.chord = 0
		self.spacing = 0
		self.stagger =0
		self.blade_angle1=0
		self.blade_angle2=0
		self.deviation_angle=0
		self.camber = 0
		self.diffusion_factor = 0
		self.T02 = 0
		self.T01 = 0
		self.delta_T_o = 0
		self.T1 = 0
		self.M_relative = 0
		self.cl_at_zero_incidence=0
		self.diffusion_factor_for_stator = 0
		

overall_pressure_ratio = 1.2
hub_tip_ratio = 0.5
tip_diameter = 0.75
hub_diameter = hub_tip_ratio * tip_diameter
mean_diameter = (hub_diameter + tip_diameter)/2

hub_radius = hub_diameter /2
tip_radius = tip_diameter/2
mean_radius = mean_diameter/2

design_speed = 10000     #rpm
operating_absolute_mach_number = 0.5
# General constants
gamma = 1.4
gas_constant = 287
stage_efficency = 0.9
Cp = 1005

R_tip = tip_diameter / mean_diameter

def pop_generator(Ca,Cw1,Cw2,T02,delta_T_o,n,aspect_ratio_blade,filename,saving_crit):
	
	a = (Cw1+Cw2)/(2*R_tip**n)
	b = ((Cw2-Cw1)*(R_tip))/2

	#print(a,b)
	blade_height = tip_radius - hub_radius

	#print(blade_height)
	total_blade_section = 20
	r_values = np.linspace(hub_radius,tip_radius,total_blade_section)

	blade_sections = np.empty(total_blade_section,dtype=object)
	for i,r in enumerate(r_values):
		blade_sections[i] = blade_section(r)
		R_val = blade_sections[i].r_value/mean_radius
		#print(blade_sections[i].r_value)
		blade_sections[i].Cw1 = a*R_val**n - b/R_val
		blade_sections[i].Cw2 = a*R_val**n + b/R_val	
		blade_sections[i].U = rotational_speed(design_speed,blade_sections[i].r_value*2)
		blade_sections[i].Rx =  Rx_calculator(blade_sections[i].Cw1,blade_sections[i].Cw2,blade_sections[i].U)
		if i<19 : 
			blade_sections[i].Ca1 = Ca_upstream( Ca,Cw1,blade_sections[i].Cw1,(tip_radius/mean_radius),R_val,a,b,n)
			blade_sections[i].Ca2 = Ca_downstream( Ca,Cw2,blade_sections[i].Cw2,(tip_radius/mean_radius),R_val,a,b,n)
		else : 
			blade_sections[i].Ca1 = Ca
			blade_sections[i].Ca2 = Ca
		blade_sections[i].alpha1 = radian_2_degree(math.atan(blade_sections[i].Cw1/blade_sections[i].Ca1) )
		blade_sections[i].alpha2 = radian_2_degree( math.atan(blade_sections[i].Cw2/blade_sections[i].Ca2) )
		blade_sections[i].C1 = np.sqrt(blade_sections[i].Cw1**2 + blade_sections[i].Ca1**2)
		blade_sections[i].C2 = np.sqrt(blade_sections[i].Cw2**2 + blade_sections[i].Ca2**2)
		blade_sections[i].beta1 = U_Ca_formula(blade_sections[i].U, blade_sections[i].Ca1 , blade_sections[i].alpha1)
		blade_sections[i].beta2 = U_Ca_formula(blade_sections[i].U, blade_sections[i].Ca2 , blade_sections[i].alpha2)
		blade_sections[i].V1 = blade_sections[i].Ca1 / np.cos(degree_2_radian(blade_sections[i].beta1))
		blade_sections[i].V2 = blade_sections[i].Ca2 / np.cos(degree_2_radian(blade_sections[i].beta2))
		blade_sections[i].de_haller_number = blade_sections[i].V2 / blade_sections[i].V1

	##### Creating camber and chord of airfoils ###############
	## All calculations at mean radius

	mean_count =0
	if total_blade_section % 2 ==0 :
		mean_count = int(total_blade_section/2 )
	else : 
		mean_count = int(total_blade_section/2)

	#rint(mean_count)
	#rint(mean_radius)
	#rint(blade_sections[mean_count].r_value)

	air_deflection_angle = blade_sections[mean_count].beta1 - blade_sections[mean_count].beta2
	air_outlet_angle = blade_sections[mean_count].beta2
	#print("air outlet angle is %f and air deflection angle is %f" %(air_outlet_angle,air_deflection_angle))

	solidity_factor = pitch_chord_ration (blade_sections[mean_count].V1,blade_sections[mean_count].V2)# from emperical data (s/c)
	#aspect_ratio_blade = input("Enter aspect ratio of blade")  #assumed value
	#aspect_ratio_blade = float(aspect_ratio_blade)
	chord_of_blade =  blade_height/aspect_ratio_blade
	blade_spacing = solidity_factor * chord_of_blade
	number_of_blades = (2*np.pi*mean_radius)/blade_spacing

	number_of_blades = int(number_of_blades) +1
	blade_spacing = (2*np.pi*mean_radius)/number_of_blades
	solidity_factor = blade_spacing/chord_of_blade

	#print(number_of_blades,solidity_factor)

	# after knowing the number of blades lets set solodity factor,spacing and chord_length for all blade section

	for i in range(total_blade_section):
		blade_sections[i].spacing = (2*np.pi*blade_sections[i].r_value)/number_of_blades
		#blade_sections[i].solidity_factor = pitch_chord_ration (blade_sections[i].V1,blade_sections[i].V2)
		blade_sections[i].chord = chord_of_blade
		blade_sections[i].solidity_factor = blade_sections[i].spacing / blade_sections[i].chord
		#print(blade_sections[i].solidity_factor,blade_sections[i].spacing/chord_of_blade)
		blade_sections[i].camber = ( blade_sections[i].beta1 - blade_sections[i].beta2 ) / \
										( 1- (0.25 + 0.1*(blade_sections[i].beta2/50))* np.sqrt(blade_sections[i].solidity_factor) )
		blade_sections[i].deviation_angle = (0.23 + 0.1*(blade_sections[i].beta2/50)) * np.sqrt(blade_sections[i].solidity_factor) * blade_sections[i].camber
		blade_sections[i].blade_angle1 = blade_sections[i].beta1
		blade_sections[i].blade_angle2 = blade_sections[i].beta2 - blade_sections[i].deviation_angle
		blade_sections[i].stagger = blade_sections[i].beta1- (blade_sections[i].camber/2)
		blade_sections[i].diffusion_factor = diffusion_factor(blade_sections[i].V1,blade_sections[i].V2,blade_sections[i].Cw1,blade_sections[i].Cw2,blade_sections[i].solidity_factor)
		blade_sections[i].T02 = T02
		blade_sections[i].delta_T_o = delta_T_o
		blade_sections[i].T01 = blade_sections[i].T02- delta_T_o
		blade_sections[i].T1= blade_sections[i].T01 - (0.5*blade_sections[i].C1**2)/(Cp)
		blade_sections[i].M_relative = blade_sections[i].V1/np.sqrt(gamma*gas_constant*blade_sections[i].T1)
		blade_sections[i].cl_at_zero_incidence = np.tan(degree_2_radian(blade_sections[i].camber) /4) / 0.1103
		blade_sections[i].diffusion_factor_for_stator = 1 - np.cos(degree_2_radian(blade_sections[i].alpha2)) + 0.5*9*np.sin(degree_2_radian(blade_sections[i].alpha2))* ( (np.cos(degree_2_radian(blade_sections[i].alpha2)))**2 - 0.433) 
		
	

	parameters = ['serial number','r_value', 'alpha1', 'beta1', 'alpha2', 'beta2', 'Rx', 'C1', 'C2', 'V1', 'V2', 'De haller number',\
					  'Cw1', 'Cw2', 'Ca1','Ca2','U','solidity factor(s/c)','chord','spacing','stagger','blade angle1','blade angle2','deviation angle','camber','diffusion_factor','T02','delta_T_0','T01','T1','M_relative','cl at zero incidence','stator diff factor']
	
	
	# constraint from stator diffusion factor and blade numbers
	#if saving_crit != 0 : print("Hello")
	if saving_crit ==0 :
		diffusion_counter = 1
		for i in range(20):
			if blade_sections[i].diffusion_factor_for_stator < 0.6 :
				diffusion_counter += 1
		if 0<number_of_blades <100 and diffusion_counter == 20:
			print("Number of blades %d" %number_of_blades)
			print(filename)
			if saving_crit !=0 : 
				save_data_excel(blade_sections,20,parameters,n,aspect_ratio_blade,number_of_blades,filename)
			return (number_of_blades)
	if saving_crit != 0 :
		save_data_excel(blade_sections,20,parameters,n,aspect_ratio_blade,number_of_blades,filename)
		print("YO")
#filename = "n%fAR%fTM%fdt%d.xls" %(0,4.5,0.8,27)
#pop_generator(178.8157,107.241,176.3396,481.4389,27,0,4.5,filename)

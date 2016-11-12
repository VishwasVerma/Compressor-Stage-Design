#author Vishwas Verma
# date Oct 08,2016

import numpy as np
import scipy as sc
from scipy import integrate
import math
import xlwt
import os
from functions import *

# Design Requirements
overall_pressure_ratio = 1.2
hub_tip_radius_ratio = 0.5
tip_diameter = 0.75
operating_absolute_mach_number = 0.5

# Calculate Geometry parameters 
hub_diameter = hub_tip_radius_ratio * tip_diameter
hub_radius = hub_diameter /2
tip_radius = tip_diameter/2
mean_radius = (hub_radius+tip_radius)/2
blade_height = tip_radius - hub_radius

# Assumed Values	convergence_function
stage_efficiency = 0.9 		# Assume Satge efficiency 
Max_relative_tip_mach_number = 1 #In order to aviod shock losses maximum relative tip mach number is set to 1

# soft constarints
design_speed = 10000 		#rpm  
de_haller_number = 0.72  # allowed to vary with penality over 0.76 and under 0.68

# other constants
gamma = 1.4		#Specific heat ratio
gas_constant = 287 	# Universal Gas Constant for air (R)
Cp = 1005   # Specific heat at constant pressure


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
		self.solidity_factor = 0   #s/c
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
		
	
# Value used in 3-D streching function   pitch_chord_relation
R_tip = tip_diameter / (mean_radius*2)

def pop_generator(Ca,Cw1,Cw2,T02,delta_T_o,n,aspect_ratio_blade,filename,saving_crit):
	# Calculation of a and b at tip of blade
	a = (Cw1+Cw2)/(2*R_tip**n)
	b = ((Cw2-Cw1)*(R_tip))/2
	total_blade_section = 20		# total number of sections in a blade
	r_values = np.linspace(hub_radius,tip_radius,total_blade_section)	# distance from central axis to blade section
	blade_sections = np.empty(total_blade_section,dtype=object)			# Blade sections at each section 
	for i,r in enumerate(r_values):
		blade_sections[i] = blade_section(r)	#initiliase each blade section with unique disatance from central axis
		R_val = blade_sections[i].r_value/mean_radius
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
		#print(blade_sections[i].Ca1,blade_sections[i].Ca2)
		blade_sections[i].alpha1 = radian_2_degree(math.atan(blade_sections[i].Cw1/blade_sections[i].Ca1) )
		blade_sections[i].alpha2 = radian_2_degree( math.atan(blade_sections[i].Cw2/blade_sections[i].Ca2) )
		blade_sections[i].C1 = np.sqrt(blade_sections[i].Cw1**2 + blade_sections[i].Ca1**2)
		blade_sections[i].C2 = np.sqrt(blade_sections[i].Cw2**2 + blade_sections[i].Ca2**2)
		blade_sections[i].beta1 = U_Ca_formula(blade_sections[i].U, blade_sections[i].Ca1 , blade_sections[i].alpha1)
		blade_sections[i].beta2 = U_Ca_formula(blade_sections[i].U, blade_sections[i].Ca2 , blade_sections[i].alpha2)
		blade_sections[i].V1 = blade_sections[i].Ca1 / np.cos(degree_2_radian(blade_sections[i].beta1))
		blade_sections[i].V2 = blade_sections[i].Ca2 / np.cos(degree_2_radian(blade_sections[i].beta2))
		blade_sections[i].de_haller_number = blade_sections[i].V2 / blade_sections[i].V1

	## All calculations at mean radius	cal_r_mean

	mean_count =0
	if total_blade_section % 2 ==0 :
		mean_count = int(total_blade_section/2 )
	else : 
		mean_count = int(total_blade_section/2)+1

	air_deflection_angle = blade_sections[mean_count].beta1 - blade_sections[mean_count].beta2
	air_outlet_angle = blade_sections[mean_count].beta2
	solidity_factor = pitch_chord_relation (blade_sections[mean_count].V1,blade_sections[mean_count].V2)	# from emperical data (s/c)
	chord_of_blade =  blade_height/aspect_ratio_blade
	blade_spacing = solidity_factor * chord_of_blade
	number_of_blades = (2*np.pi*mean_radius)/blade_spacing
	#print(type(number_of_blades))
	#print("Number of blades %f" %number_of_blades)
	number_of_blades = int(number_of_blades) +1
	blade_spacing = (2*np.pi*mean_radius)/number_of_blades
	solidity_factor = blade_spacing/chord_of_blade

	for i in range(total_blade_section): 
		blade_sections[i].spacing = (2*np.pi*blade_sections[i].r_value)/number_of_blades
		blade_sections[i].chord = chord_of_blade
		blade_sections[i].solidity_factor = blade_sections[i].spacing / blade_sections[i].chord
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
	if saving_crit ==0 :
		diffusion_counter = 0
		pitch = 0
		for i in range(20):
			if blade_sections[i].diffusion_factor_for_stator < 0.6 :
				diffusion_counter += 1
		if pitch_chord_relation(blade_sections[mean_count].C2,blade_sections[mean_count].Ca2) > 0.45 : pitch +=1
		if 39<number_of_blades <100 and diffusion_counter>18 and pitch==1 :
			print("Number of blades %d" %number_of_blades)
			print(filename)
			return (number_of_blades)
	if saving_crit != 0 :
		save_data_excel(blade_sections,20,parameters,n,aspect_ratio_blade,number_of_blades,filename)

## End of the programme
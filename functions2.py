import numpy as np
import scipy as sc
from scipy import integrate
import math
import random
import csv
import xlwt
import os

gamma,gas_constant,Cp,stage_efficency = 1.4,287,1005,0.85
overall_pressure_ratio = 1.2
operating_absolute_mach_number = 0.5
design_speed = 10000
hub_tip_ratio = 0.5
tip_diameter = 0.75
hub_diameter = hub_tip_ratio * tip_diameter
mean_diameter = (hub_diameter + tip_diameter)/2
hub_radius = hub_diameter /2
tip_radius = tip_diameter/2
mean_radius = mean_diameter/2

def deviation_angle(camber,solidity):
	x = 2 # experimental factor 
	deviation = (mc*camber)/math.sqrt(solidity) + x
	return deviation
	
def pitch_chord_ration (W1,W2):
	# Using McKenize Method we will calculate Solidity factor
	Cp = 1-(W2/W1)**2  # Cp = 1-dH^2
	solidity = 9*(0.567-Cp) # 1/solidity = S/C
	return solidity

def cal_r_mean(hub_tip_ratio,tip_dia): 
	hub_diameter = hun_tip_ratio * tip_dia
	mean_diameter = (hub_diameter + tip_dia)/2
	return ( mean_diameter)
	
def rotational_speed(N,m_d):
	return (np.pi*(N/60)*m_d)
'''
def convergence_function(delta_T_o):
	# parameters at tip of blade
	# Solving equation as stated above :
	# from calculated expression
	expression = ( 1+((gamma*gas_constant)/(8*Cp)) ) * (1/stage_efficency) * ( ( overall_pressure_ratio)**((gamma-1)/gamma) -1)
	#print(delta_T_o,expression)
	T1 = [dT/expression for dT in delta_T_o]
	#print(T1)
	count = 0
	k_val = np.linspace(0.65,0.95,31)
	for k in k_val:
		for i,T in enumerate(T1):
			#alpha1= (1/np.sqrt(gamma*gas_constant*T) ) * ( 125*np.pi - (T/2.51) ) 
			alpha1 = (1/np.sqrt(gamma*gas_constant*T) ) * ( 125*np.pi - T*( (k**2 -0.25)*( (gamma*gas_constant)/(125*np.pi)) ) ) 
			T02 =  delta_T_o[i] + T*(1+ ((gamma*gas_constant) /(8*Cp) ) )
			#print (alpha1,T02)
			if alpha1 <=0.7 and 460<T02<=500 :
				print("Got a value")
				#print(T02,alpha1,k,T,delta_T_o[i])
				return (T02,alpha1,k,T,delta_T_o[i])
'''

def convergence_function(delta_T_o):
	# parameters at tip of blade
	expression = ( 1+((gamma*gas_constant)/(8*Cp)) ) * (1/stage_efficency) * ( ( overall_pressure_ratio)**((gamma-1)/gamma) -1)
	T1 = [dT/expression for dT in delta_T_o]
	count = 0
	k_val = np.linspace(0.6,0.95,21)
	for k in k_val:
		for i,T in enumerate(T1):
			#alpha1= (1/np.sqrt(gamma*gas_constant*T) ) * ( 125*np.pi - (T/2.51) ) 
			alpha1 = (1/np.sqrt(gamma*gas_constant*T) ) * ( 125*np.pi - T*( (k**2 -0.25)*( (gamma*gas_constant)/(125*np.pi)) ) ) 
			T02 =  delta_T_o[i] + T*(1+ ((gamma*gas_constant) /(8*Cp) ) )
			C1 = operating_absolute_mach_number * np.sqrt(gamma*gas_constant*T)
			Ca = axial_velocity_given_total_vel(C1,radian_2_degree(alpha1))
			U = rotational_speed(design_speed,tip_diameter)
			V1 = np.sqrt(gamma*gas_constant*T) * k
			beta1 = U_Ca_formula(U,Ca,radian_2_degree(alpha1))
			Cw1 = tangential_velocity_given_total_velocity(C1,radian_2_degree(alpha1))
			Cw2 = Cw1 + (delta_T_o[i]*Cp)/U
			Rx = 1 - (Cw1+Cw2)/(2*U)
			beta2 = math.atan( (Rx*2*U)/Ca - math.tan(degree_2_radian(beta1)))	#radian
			alpha2 = U_Ca_formula(U,Ca,radian_2_degree(beta2))		#degree
			alpha2 = degree_2_radian(alpha2)
			diffusion_factor_constraint_from_stator = 1-np.cos(alpha2) + 0.5*np.sin(alpha2)*9*( (np.cos(alpha2))**2 - 0.433)
			solidity = 9*( (np.cos(alpha2))**2 - 0.433 )
			#print("HELLO HERE")
			#print(solidity,diffusion_factor_constraint_from_stator)
			
			if diffusion_factor_constraint_from_stator<0.6 and 420<T02<500 and 0.6<solidity :
				print("Got a value")
				#print(T02,alpha1,k,T,delta_T_o[i])
				return (T02,alpha1,k,T,delta_T_o[i],C1,Ca,U,V1,beta1,Cw1,Cw2,Rx,radian_2_degree(beta2),radian_2_degree(alpha2))

def total_temperature(T,velocity):
	return ( T+ velocity**2/(2*Cp) )

def axial_velocity_given_total_vel(velocity,angle):
	angle = angle*(np.pi/180)
	return (velocity*np.cos(angle))

def tangential_velocity_given_total_velocity(velocity,angle):
	angle = degree_2_radian(angle)
	return(velocity*np.sin(angle))

def degree_2_radian(degree):
	return( degree* (np.pi/180))

def radian_2_degree(radian):
	return (radian* (180/np.pi))

def U_Ca_formula (U,Ca,angle1):
	angle1 = degree_2_radian(angle1)
	angle2 = math.atan( (U/Ca) - np.tan(angle1) )
	angle2 = radian_2_degree(angle2)
	return (angle2)

def downstream_flows(beta1,V1,Ca,U,Cw1,delta_T_o):
	Cp=1005
	Cw2 = Cw1 + (delta_T_o*Cp)/U
	Rx = 1 - (Cw1+Cw2)/(2*U)
	beta2 = math.atan( (Rx*2*U)/Ca - math.tan(beta1))
	V2 = Ca/np.cos(beta2)
	de_haller_number = V2/V1
	return (Rx,radian_2_degree(beta2), de_haller_number,Cw2,V2)

def Rx_calculator(Cw1,Cw2,U):
	Rx = 1- (Cw1+Cw2)/(2*U)
	return (Rx)
			

# calculating Ca1 based on exponential law (upstream)
def Ca_upstream(Ca_t, Cw_t, Cw_r, R_t , R_r,a,b,n):
	fun = lambda R : (a*R**n - b/R)/R
	value = integrate.quad(fun, R_r,R_t)[0]
	Ca = np.sqrt( Ca_t**2 + (Cw_t**2 - Cw_r**2)  + value )
	return (Ca)
# Ca downstream
def Ca_downstream(Ca_t, Cw_t, Cw_r, R_t , R_r,a,b,n):
	fun = lambda R : (a*R**n + b/R)/R
	value = integrate.quad(fun, R_r,R_t)[0]
	Ca = np.sqrt( Ca_t**2 + (Cw_t**2 - Cw_r**2)  + value )
	return (Ca)

def diffusion_factor(V1,V2,Cw1,Cw2,solidity):
	d = 1 -(V2/V1) + (((Cw2-Cw1)/(2*V1))*solidity)
	return(d)			

## To save data into excel file
def save_data_excel(blade,total_number,parameters,n,aspect_ratio_blade,number_of_blades,filename):
	book = xlwt.Workbook(encoding="utf-8")
	sheet1 = book.add_sheet("Sheet 1")
	total_parameters = len(parameters)
	for i in range(total_parameters):
		sheet1.write(0,i,parameters[i])
	
	r_value,Rx,de_haller_number,chord,stagger = [],[],[],[],[]
	alpha1,alpha2,beta1,beta2 =  [],[],[],[]
	C1,C2,V1,V2,Cw1,Cw2,Ca1,Ca2,U = [],[],[],[],[],[],[],[],[]
	diffusion_factor = []
	solidity_factor,spacing,blade_angle1,blade_angle2,deviation_angle,camber = [],[],[],[],[],[]
	T02,T01,T1,delta_T_o,M_relative,cl_at_zero_incidence = [],[],[],[],[],[]
	diffusion_factor_for_stator = []
	for i in range(total_number):
		r_value.append(blade[i].r_value)
		alpha1.append(blade[i].alpha1)
		alpha2.append(blade[i].alpha2)
		beta1.append(blade[i].beta1)
		beta2.append(blade[i].beta2)
		Rx.append(blade[i].Rx)
		C1.append(blade[i].C1)
		C2.append(blade[i].C2)
		V1.append(blade[i].V1)
		V2.append(blade[i].V2)
		de_haller_number.append(blade[i].de_haller_number)
		Cw1.append(blade[i].Cw1)
		Cw2.append(blade[i].Cw2)
		Ca1.append(blade[i].Ca1)
		Ca2.append(blade[i].Ca2)
		U.append(blade[i].U)
		solidity_factor.append(blade[i].solidity_factor)
		spacing.append(blade[i].spacing)
		blade_angle1.append(blade[i].blade_angle1)
		blade_angle2.append(blade[i].blade_angle2)
		deviation_angle.append(blade[i].deviation_angle)
		camber.append(blade[i].camber)
		chord.append(blade[i].chord)
		stagger.append(blade[i].stagger)
		diffusion_factor.append(blade[i].diffusion_factor)
		T02.append(blade[i].T02)
		T01.append(blade[i].T01)
		T1.append(blade[i].T1)
		delta_T_o.append(blade[i].delta_T_o)
		M_relative.append(blade[i].M_relative)
		cl_at_zero_incidence.append(blade[i].cl_at_zero_incidence)
		diffusion_factor_for_stator.append(blade[i].diffusion_factor_for_stator)
	
	for i in range(1,(total_number+1)):
		sheet1.write(i,0,i)
		sheet1.write(i,1,r_value[i-1])
		sheet1.write(i,2,alpha1[i-1])
		sheet1.write(i,3,beta1[i-1])
		sheet1.write(i,4,alpha2[i-1])
		sheet1.write(i,5,beta2[i-1])
		sheet1.write(i,6,Rx[i-1])
		sheet1.write(i,7,C1[i-1])
		sheet1.write(i,8,C2[i-1])
		sheet1.write(i,9,V1[i-1])
		sheet1.write(i,10,V2[i-1])
		sheet1.write(i,11,de_haller_number[i-1])
		sheet1.write(i,12,Cw1[i-1])
		sheet1.write(i,13,Cw2[i-1])
		sheet1.write(i,14,Ca1[i-1])
		sheet1.write(i,15,Ca2[i-1])
		sheet1.write(i,16,U[i-1])
		sheet1.write(i,17,solidity_factor[i-1])
		sheet1.write(i,18,chord[i-1])
		sheet1.write(i,19,spacing[i-1])
		sheet1.write(i,20,stagger[i-1])
		sheet1.write(i,21,blade_angle1[i-1])
		sheet1.write(i,22,blade_angle2[i-1])
		sheet1.write(i,23,deviation_angle[i-1])
		sheet1.write(i,24,camber[i-1])
		sheet1.write(i,25,diffusion_factor[i-1])
		sheet1.write(i,26,T02[i-1])
		sheet1.write(i,27,delta_T_o[i-1])
		sheet1.write(i,28,T01[i-1])
		sheet1.write(i,29,T1[i-1])
		sheet1.write(i,30,M_relative[i-1])
		sheet1.write(i,31,cl_at_zero_incidence[i-1])
		sheet1.write(i,32,diffusion_factor_for_stator[i-1])
	sheet1.write(i+2,1,"aspect_ratio")
	sheet1.write(i+2,2,"n_value")
	sheet1.write(i+2,3,"number of blades")
	sheet1.write(i+3,1,aspect_ratio_blade)
	sheet1.write(i+3,2,n)
	sheet1.write(i+3,3,number_of_blades)
	try :
		os.remove(filename)
	except OSError:
		pass
	book.save(filename)
	print("Done Enjoy!!")
	
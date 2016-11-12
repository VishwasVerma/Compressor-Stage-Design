#Author Vishwas Verma
# Writen For course Project AE651

from functions import *      # functions file
from bladedesign import * 		# pop2generate filw
import numpy as np
from scipy import integrate
import math
import xlwt
import random
import os
import time

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

# Assumed Values
stage_efficiency = 0.9 		# Assume Satge efficiency 
Max_relative_tip_mach_number = 1 #In order to aviod shock losses maximum relative tip mach number is set to 1

# soft constarints
design_speed = 10000 		#rpm  
#de_haller_number = 0.72  # allowed to vary with penality over 0.76 and under 0.68

# other constants
gamma = 1.4		#Specific heat ratio
gas_constant = 287 	# Universal Gas Constant for air (R)
Cp = 1005   # Specific heat at constant pressure


# After analysing design variables, these are the independent variables which limits whole design
# Now set these values from Aerodynamic efficiency point of view
delta_T_o =  np.linspace(18,30,13)				#Temperature change across one stage
tip_mach_number = np.linspace(0.7,0.95,26)		#Relative Mach Number at tip
aspect_ratio = np.linspace(3,5,5) 				#Aspect Ratio of Blade, Since blade height is large so assume these values
n_values = np.linspace(0,1,11) 					#3-D streching parameter

# This function calculate design paramters at tip of blade, limiting tip mach number to be 1.
def upstream_flow(overall_pressure_ratio,operating_absolute_mach_number):
	expression = ( 1+((gamma*gas_constant)/(8*Cp)) ) * (1/stage_efficiency) * \
						( ( overall_pressure_ratio)**((gamma-1)/gamma) -1) 			# this term comes from pressure ratio and efficiency
	T1 = [dT/expression for dT in delta_T_o]
	Values = np.zeros((len(tip_mach_number),len(delta_T_o),19)) 				#To store data
	for j,k in enumerate(tip_mach_number):
		for i,T in enumerate(T1): 
			alpha1 = (1/np.sqrt(gamma*gas_constant*T) ) * \
						( 125*np.pi - T*( (k**2 -0.25)*( (gamma*gas_constant)/(125*np.pi)) ) )
			#print(alpha1)
			if alpha1>1:	pass
			else :
				alpha1 = radian_2_degree(math.asin(alpha1))
				T02 =  delta_T_o[i] + T*(1+ ((gamma*gas_constant) /(8*Cp) ) )
				C1 = operating_absolute_mach_number * np.sqrt(gamma*gas_constant*T)
				T01 = total_temperature(T,C1)
				Ca = axial_velocity_given_total_vel(C1,alpha1)
				U = rotational_speed(design_speed,tip_diameter)
				beta1 = U_Ca_formula(U,Ca,alpha1)
				beta1 = degree_2_radian(beta1)
				Cw1 = np.sin(degree_2_radian(alpha1)) * C1
				Cw2 = Cw1 + (delta_T_o[i]*Cp)/U
				Rx = 1 - (Cw2+Cw1)/(2*U)
				beta2 = math.atan( (Rx*2*U)/Ca - math.tan(beta1))
				alpha2 = math.atan(Cw2/Ca)
				V1 = Ca/np.cos(beta1)
				V2 = Ca/np.cos(beta2)
				de_haller_number = V2/V1
				solidity =  pitch_chord_relation (V1,V2)
				diffusion_factor = 1 - (V2/V1) + (((Cw2-Cw1)/(2*V1))*solidity)

				Values[j][i][2] = alpha1
				Values[j][i][3] = radian_2_degree(beta1)
				Values[j][i][4] = radian_2_degree(alpha2)
				Values[j][i][5] = radian_2_degree(beta2)
				Values[j][i][6] = T01
				Values[j][i][7] = T02
				Values[j][i][8] = de_haller_number
				Values[j][i][9] = Rx
				Values[j][i][1] = delta_T_o[i]
				Values[j][i][10] = C1
				Values[j][i][11] = Cw1
				Values[j][i][12] = Cw2
				Values[j][i][13] = Ca
				Values[j][i][14] = V1
				Values[j][i][15] = V2
				Values[j][i][16] = solidity
				Values[j][i][17] = diffusion_factor
				Values[j][i][18] = U
				Values[j][i][0] = k #tip_mach_number
			
	parameters = ['tip_mach_number','delta t0','alpha1', 'beta1', 'alpha2', 'beta2','T01','T02', 'De haller number', 												'Rx','C1', 'Cw1', 'Cw2', 'Ca', 'V1', 'V2','solidity','diffusion factor','U'] 		
	return(Values,parameters)
# the above function returns tipmachnumber,temperature chage,alpha1,alpha2,beta1,beta2,T01,T02,deHallerNumber,Rx.C1.Cw1.Cw2.Ca.V1.V2.solidity,diffusion factor and Rotational Speed at tip of blade.

Values,parameters= upstream_flow(overall_pressure_ratio,operating_absolute_mach_number)
'''
filename = 'first data'
try :
	os.remove("%s.xls" %filename)
except OSError:
	pass

#save_data(Values,tip_mach_number,parameters,delta_T_o,filename) 	# Saves the databin excel file.
'''
count = 0
blade_numbers = []
possible_tip_mach,possible_n_values,possible_AR,possible_delta_T_o = [],[],[],[]
population_size = 0
MAX_VELOCITY = [0.95,30,1,5]	# product of tip_mach_number,dT,n,AR
MIN_VELOCITY = [0.7,18,0,3]
#popsize=	# to have finite number of populations

for j in range(len(tip_mach_number)):
		for m in range(len(delta_T_o)):
			for k,n in enumerate(n_values):
				for l,AR in enumerate(aspect_ratio):
					filename = "data5/n%dAR%dc%d.xls" %(k,l,count)
					Ca,Cw1,Cw2,T02,dT = Values[j][m][13],Values[j][m][11],Values[j][m][12],Values[j][m][7],Values[j][m][1]
					blade_number = pop_generator(Ca,Cw1,Cw2,T02,dT,n,AR,filename,0) # dont want to save data thats why 0 to save data set1
					if blade_number != None : 	# total number of blades, its possible that we may get None because of the exit criterion set at stator, stator exit angle should be zero.
						blade_numbers.append(blade_number) 
						possible_tip_mach.append(tip_mach_number[j])
						possible_n_values.append(n)
						possible_AR.append(AR)
						possible_delta_T_o.append(delta_T_o[m])
						population_size +=1
					#if population_size > popsize : break
				#if population_size > popsize : break	
			count +=1
		#if population_size >popsize : break

popsize= (population_size)
possible_sol  = np.empty ( (4,population_size))

print("Total population size is %d" %population_size)
print([numbers for numbers in blade_numbers if 40<numbers<70])
print(blade_numbers)
for i in range(population_size):
	possible_sol[0][i] = possible_tip_mach[i]
	possible_sol[1][i] = possible_delta_T_o[i]
	possible_sol[2][i] = possible_n_values[i]
	possible_sol[3][i] = possible_AR[i]

def save_initials(parameters,blade_numbers,total_population,filename):
	book = xlwt.Workbook(encoding="utf-8")
	sheet1 = book.add_sheet("Sheet 1")
	header = [ 'n' , 'AR', 'tip_mach','dT','no. blades' ]
	for i in range(5):
		sheet1.write(0,i,header[i])
	
	for j in range(1,(total_population+1)):
		sheet1.write(j,2,parameters[0][j-1])	# tip_mach
		sheet1.write(j,3,parameters[1][j-1])	# dT
		sheet1.write(j,0,parameters[2][j-1])	# n
		sheet1.write(j,1,parameters[3][j-1]) 	#AR
		sheet1.write(j,4,blade_numbers[j-1])	#nblades
	
	try :
		os.remove(filename)
	except OSError:
		pass
	book.save(filename)
	print("Done Enjoy!!")
	
save_initials(possible_sol,blade_numbers,population_size,'data.xls')

# to calculate required parameters at tip of blade
def upstream_flow_at_tip(overall_pressure_ratio,operating_absolute_mach_number,dT,tip_mach_number):
	expression = ( 1+((gamma*gas_constant)/(8*Cp)) ) * (1/stage_efficiency) * \
						( ( overall_pressure_ratio)**((gamma-1)/gamma) -1)
	T1 = dT/expression	
	Values = np.zeros(7)	
	alpha1 = (1/np.sqrt(gamma*gas_constant*T1) ) * \
				( 125*np.pi - T1*( (tip_mach_number**2 -0.25)*( (gamma*gas_constant)/(125*np.pi)) ) )
	alpha1 = radian_2_degree(math.asin(alpha1))
	T02 =  dT + T1*(1+ ((gamma*gas_constant) /(8*Cp) ) )
	C1 = operating_absolute_mach_number * np.sqrt(gamma*gas_constant*T1)
	T01 = total_temperature(T1,C1)
	Ca = axial_velocity_given_total_vel(C1,alpha1)
	U = rotational_speed(design_speed,tip_diameter)
	beta1 = U_Ca_formula(U,Ca,alpha1)
	beta1 = degree_2_radian(beta1)
	Cw1 = np.sin(degree_2_radian(alpha1)) * C1
	Cw2 = Cw1 + (dT*Cp)/U
	Rx =  Rx_calculator(Cw1,Cw2,U)
	beta2 = math.atan( (Rx*2*U)/Ca - math.tan(beta1))
	alpha2 = math.atan(Cw2/Ca)
	#tip_diffusion = stator_diffusion_checker(radian_2_degree(alpha2))
	Values[0] = Ca
	Values[1] = Cw1
	Values[2] = Cw2
	Values[3] = T02
	Values[4] = dT
	Values[5] = radian_2_degree(beta1)
	Values[6] = radian_2_degree(beta2)
	#Values[7] = tip_diffusion
	
	return(Values)
#returns at tip Ca,Cw1,Cw2,T02,delta_T_o,beta1,beta2,tip_diffusion
# To calculate efficiencies facotors we need these vales at given blade section
# find at mean beta1,beta2, solidity_factor,mean_spacing
# find at tip camber
# find T02 and T01

#pop_generator(Ca,Cw1,Cw2,T02,delta_T_o,n,aspect_ratio_blade,filename)
# final blade generation
#returns at tip Ca,Cw1,Cw2,T02,delta_T_o,beta1,beta2
def final_blade_generation(filename):
	delta_T_o = input("Enter dT :")
	delta_T_o = float(delta_T_o)
	tip_mach_number = input("Enter tip mach number :")
	tip_mach_number = float(tip_mach_number)
	n = input("Enter n value :")
	n = float(n)
	AR = input("Enter AR value :")
	AR = float(AR)
	
	Ca,Cw1,Cw2,T02,delat_T_o,tip_beta1,tip_beta2 = upstream_flow_at_tip(overall_pressure_ratio,operating_absolute_mach_number,dT,tip_mach_number)
	
	number_of_blades = pop_generator(Ca,Cw1,Cw2,T02,delta_T_o,n,AR,filename,1)
	print("saved")
	
final_blade_generation("rotor_blade.xls")




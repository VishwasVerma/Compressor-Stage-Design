#author Vishwas Verma
# date Oct 08,2016
#Function file
import numpy as np
from scipy import integrate
import math
import random
import os
import xlwt
Cp,gamma,stage_efficency = 1005,1.4,0.9	

def pitch_chord_relation (W1,W2):
	# Using McKenize Method we will calculate Solidity factor
	Cp = 1-(W2/W1)**2  # Cp = 1-dH^2
	solidity = 9*(0.567-Cp) # solidity = S/C
	return solidity

def rotational_speed(N,m_d):
	return (np.pi*(N/60)*m_d)

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

def Rx_calculator(Cw1,Cw2,U):
	Rx = 1- (Cw1+Cw2)/(2*U)
	return (Rx)

def save_data(Values,tip_mach_number,parameters,delta_T_o,filename):
	book = xlwt.Workbook(encoding="utf-8")
	total_parameters = len(parameters)
	sheet1 = book.add_sheet("Sheet 1")
	count=0
	for i in range(total_parameters):
		sheet1.write(0,i,parameters[i])
	for i,t_mach in enumerate(tip_mach_number):
		#count +=1
		for j in range(len(delta_T_o)):
			if Values[i][j][7] == 0 : 
				continue
			else :
				count +=1
				for k in range(len(parameters)):
					sheet1.write( count,k,Values[i][j][k]) 	
	book.save("%s.xls" %filename)
	print("Done Enjoy!!")			

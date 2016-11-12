import xlwt
import xlrd
import numpy as np
import matplotlib.pyplot as plt
from functions import *

def read_data(filename,column):
	data = []
	book = xlrd.open_workbook(filename)
	sheet = book.sheet_by_index(0)
	for row in range(1,sheet.nrows):
		data.append(sheet.cell(row,column).value)
	#data = [float(d) for d  in data]
	return (data)

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

# other constants
gamma = 1.4		#Specific heat ratio
gas_constant = 287 	# Universal Gas Constant for air (R)
Cp = 1005   # Specific heat at constant pressure

# Total number of blade sections
total_blade_sections = 20

class stator_blades:
	
	def __init__(self,r_val):
		self.r_values = r_val
		self.C2 = 0
		self.C3 = 0
		self.Ca2 = 0
		self.Ca3 = 0
		self.Cw2 = 0
		self.Cw3 = 0
		self.alpha2 = 0
		self.alpha3 = 0
		self.T02 = 0
		self.T2 = 0
		self.absolute_mach = 0
		self.deflection = 0
		self.camber = 0
		self.stagger = 0
		self.chord =0
		self.spacing = 0
		self.deviation_angle = 0
		self.diffusion_factor = 0
		self.cl_at_zero_incidence = 0
		self.blade_angle1 = 0
		self.blade_angle2 = 0
		
	
	def set_component_vel(self):
		self.Cw2 = np.sin(degree_2_radian(self.alpha2))*self.C2
		self.Ca2 = np.cos(degree_2_radian(self.alpha2))*self.C2
		self.Ca3 = self.Ca2
		self.Cw3 = math.tan(degree_2_radian(self.alpha3))*self.Ca3
		
	def set_C3(self):
		self.C3 = np.sqrt(self.Cw3**2 + self.Ca3**2)
		
	def set_T2(self):
		self.T2 = self.T02 - (self.C2**2)/(2*Cp)

	def set_mach_number(self):
		self.absolute_mach = self.C2 / (np.sqrt(gamma*gas_constant*self.T2) )
		
	def set_deflection(self):
		self.deflection = degree_2_radian(self.alpha2) - degree_2_radian(self.alpha3)
	
	def set_blade_geometry(self,number_of_blades,AR):
		self.spacing = 2*np.pi*self.r_values/number_of_blades
		self.chord = blade_height/AR
		self.solidity_factor = self.spacing / self.chord
		self.camber = ( self.alpha2 - self.alpha3 ) / \
								( 1- (0.25 + 0.1*(self.alpha3/50))* np.sqrt(self.solidity_factor) )
		self.deviation_angle = (0.23 + 0.1*(self.alpha3/50)) * np.sqrt(self.solidity_factor) * self.camber
		self.blade_angle1 = self.alpha2
		self.blade_angle2 = self.alpha3 - self.deviation_angle
		self.stagger = self.alpha2 - self.camber/2
		self.diffusion_factor =  1 - (self.C3/self.C2) + ((self.Cw2-self.Cw3)/(2*self.C2))* 9*( (self.C3/self.C2)**2 - 0.433) 
		self.cl_at_zero_incidence = np.tan(degree_2_radian(self.camber) /4) / 0.1103

def stator_initilise(total_blade_sections):
	# data reading from rotor file
	filename = "rotor_blade.xls"#'blade_geometry.xls'
	
	C2 = read_data(filename,8)
	alpha2 = read_data(filename,4)
	r_values = read_data(filename,1)
	T02 = read_data(filename,26)
	# Design Stator
	stator_section = np.empty(total_blade_sections,dtype = object)
	for i in range(total_blade_sections):
		stator_section[i] = stator_blades(r_values[i])
		stator_section[i].C2 = C2[i]
		stator_section[i].alpha2 = alpha2[i]
		stator_section[i].T02 = T02[i]
		stator_section[i].alpha3 = 0
		stator_section[i].set_component_vel()
		stator_section[i].set_C3()
		stator_section[i].set_T2()
		stator_section[i].set_mach_number()
		stator_section[i].set_deflection()
	# Calculation at mean
	mean_count =0
	if total_blade_sections % 2 ==0 :
		mean_count = int(total_blade_sections/2 )
	else : 
		mean_count = int(total_blade_sections/2 + 1)
	
	solidity_factor = pitch_chord_relation (stator_section[mean_count].C2,stator_section[mean_count].C3)
	# Set Aspect Ratio of blade = 4.0 Same as for rotor design
	aspect_ratio = 4.0
	chord = blade_height/aspect_ratio
	spacing = solidity_factor * chord
	number_of_blades = 2*np.pi*mean_radius / spacing
	# make number of blades to integer value
	number_of_blades = int(number_of_blades) +1
	#blade_spacing = (2*np.pi*mean_radius)/number_of_blades
	#solidity_factor = blade_spacing/chord
	#print("Total number of blades %d , solidity is %f" %(number_of_blades,solidity_factor))
	
	for i in range(total_blade_sections):
		stator_section[i].set_blade_geometry(number_of_blades,aspect_ratio)
	
	return (stator_section,number_of_blades,aspect_ratio)
#diff_stator = 1 - np.cos(degree_2_radian(alpha2)) + 0.5*9*np.sin(degree_2_radian(alpha2))* ( (np.cos(degree_2_radian(alpha2)))**2 - 0.433)

	
def data_saving_stator(stator,total_blade_sections,parameters,filename,number_of_blades,AR):
	book = xlwt.Workbook(encoding="utf-8")
	sheet1 = book.add_sheet("Sheet 1")
	total_parameters = len(parameters)
	for i in range(total_parameters):
		sheet1.write(0,i,parameters[i])
	
	r_val,C2,C3,Cw2,Cw3,Ca2 = [],[],[],[],[],[]	#6
	alpha2,alpha3,T02,T2,absolute_mach,deflection,camber,stagger,chord,spacing,deviation_angle = [],[],[],[],[],[],[],[],[],[],[]	#11
	diffusion_factor,cl_zero_incidence, blade_angle1,blade_angle2 = [],[],[],[]	#4
	solidity_factor = []
	
	for i in range(total_blade_sections):
		r_val.append(stator[i].r_values)
		C2.append(stator[i].C2)
		C3.append(stator[i].C3)
		Ca2.append(stator[i].Ca2)
		#Ca3.append(stator[i].Ca3)
		Cw2.append(stator[i].Cw2)
		Cw3.append(stator[i].Cw3)
		alpha2.append(stator[i].alpha2)
		alpha3.append(stator[i].alpha3)
		T02.append(stator[i].T02)
		T2.append(stator[i].T2)
		absolute_mach.append(stator[i].absolute_mach)
		deflection.append(stator[i].deflection)
		camber.append(stator[i].camber)
		stagger.append(stator[i].stagger)
		chord.append(stator[i].chord)
		spacing.append(stator[i].spacing)
		deviation_angle.append(stator[i].deviation_angle)
		diffusion_factor.append(stator[i].diffusion_factor)
		cl_zero_incidence.append(stator[i].cl_at_zero_incidence)
		blade_angle1.append(stator[i].blade_angle1)
		blade_angle2.append(stator[i].blade_angle2)
		solidity_factor.append(stator[i].solidity_factor)
		
	for i in range(1,(total_blade_sections+1)):
		sheet1.write(i,0,i)
		sheet1.write(i,1,r_val[i-1])
		sheet1.write(i,2,alpha2[i-1])
		sheet1.write(i,3,alpha3[i-1])
		sheet1.write(i,4,C2[i-1])
		sheet1.write(i,5,Ca2[i-1])
		sheet1.write(i,6,Cw2[i-1])
		sheet1.write(i,7,C3[i-1])
		sheet1.write(i,8,Cw3[i-1])
		sheet1.write(i,9,solidity_factor[i-1])
		sheet1.write(i,10,chord[i-1])
		sheet1.write(i,11,spacing[i-1])
		sheet1.write(i,12,camber[i-1])
		sheet1.write(i,13,blade_angle1[i-1])
		sheet1.write(i,14,blade_angle2[i-1])
		sheet1.write(i,15,deviation_angle[i-1])
		sheet1.write(i,16,stagger[i-1])
		sheet1.write(i,17,T02[i-1])
		sheet1.write(i,18,T2[i-1])
		sheet1.write(i,19,absolute_mach[i-1])
		sheet1.write(i,20,deflection[i-1])
		sheet1.write(i,21,diffusion_factor[i-1])
		sheet1.write(i,22,cl_zero_incidence[i-1])
	
	sheet1.write(24,1,"number_of_blades")
	sheet1.write(25,1,number_of_blades)
	sheet1.write(24,2,"Aspect Ratio")
	sheet1.write(25,2,AR)
	
	try :
		os.remove(filename)
	except OSError:
		pass
	book.save(filename)
	print("Done Enjoy!!")

parameters = [ 'serial number', 'r_value','alpha2','alpha3','C2','Ca2','Cw2','C3','Cw3','solidity_factor','chord','spacing','camber','blade_angle1','blade_angle2','deviation_angle','stagger','T02','T2','absolute_mach','deflection','diffusion_factor','cl_at_zero_incidence']

stator,number_of_blades,AR = stator_initilise(total_blade_sections)
filename = "stator_blade.xls"
data_saving_stator(stator,total_blade_sections,parameters,filename,number_of_blades,AR)



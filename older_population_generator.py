from functions import *
from pop2generate import *
import numpy as np
from scipy import integrate
import math
import xlwt
import random
import os
import time
# hard constraints

overall_pressure_ratio = 1.2
hub_tip_radius_ratio = 0.5
tip_diameter = 0.75
hub_diameter = hub_tip_radius_ratio * tip_diameter

hub_radius = hub_diameter /2
tip_radius = tip_diameter/2
mean_radius = (hub_radius+tip_radius)/2
stage_efficiency = 0.9
operating_absolute_mach_number = 0.5
relative_tip_mach_number = 0.9999 # has to be alwasys less than one

# soft constarints
design_speed = 10000 #rpm  # appropriate penality is used if it goes below 9000 and higher than 11000
de_haller_number = 0.72  # allowed to vary with penality over 0.76 and under 0.68
#  penality is imposed if it exceeds 30 and less than 15
#alpha1 = 30 # penality is imposed if it exceeds above 30 degree
#low_Rx = 0.35
#high_Rx = 0.75 # appropriate penality is imposed violating these condition
#n_value = 0.5  # appropriate penality is imposed if it exceeds above 1 and goes below 0

# other constants
gamma = 1.4
gas_constant = 287
Cp = 1005
# independent parameters
delta_T_o =  np.linspace(20,30,11)
tip_mach_number = np.linspace(0.7,0.95,26)		#free parameter
aspect_ratio = np.linspace(3,5,5) # free parameter
n_values = np.linspace(0,1,10) #free parameter
#diffusion_factor_stator = np.linspace(0.5,0.6,10) # hard constraint

def upstream_flow(overall_pressure_ratio,operating_absolute_mach_number):
	
	delta_T_o =  np.linspace(20,30,11)		#free parameter
	expression = ( 1+((gamma*gas_constant)/(8*Cp)) ) * (1/stage_efficiency) * \
						( ( overall_pressure_ratio)**((gamma-1)/gamma) -1)
	T1 = [dT/expression for dT in delta_T_o]
	#print(T1)
	tip_mach_number = np.linspace(0.7,0.95,26)		#free parameter
	#aspect_ratio = np.linspace(3,5,5) # free parameter
	#n = np.linspace(0.05,1,10) #free parameter
	
	Values = np.zeros((26,16,19)) #alpha1*T02
	for j,k in enumerate(tip_mach_number):
		for i,T in enumerate(T1): 
			alpha1 = (1/np.sqrt(gamma*gas_constant*T) ) * \
						( 125*np.pi - T*( (k**2 -0.25)*( (gamma*gas_constant)/(125*np.pi)) ) )
			#print(alpha1)
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

Values,parameters= upstream_flow(overall_pressure_ratio,operating_absolute_mach_number)

filename = 'first data'
try :
	os.remove("%s.xls" %filename)
except OSError:
	pass
save_data(Values,tip_mach_number,parameters,delta_T_o,filename)

# independent parameters
delta_T_o =  np.linspace(20,30,11)
tip_mach_number = np.linspace(0.7,0.95,26)		#free parameter
aspect_ratio = np.linspace(3,5,5) # free parameter
n_values = np.linspace(0,1,10) #free parameter
#diffusion_factor_stator = np.linspace(0.5,0.6,10) # hard constraint
#pop_generator(Ca,Cw1,Cw2,T02,delta_T_o,n,aspect_ratio_blade,filename)

#print(len(delta_T_o))
count = 0
blade_numbers = []
blade_height = tip_radius - hub_radius
#print(delta_T_o[4])
possible_tip_mach,possible_n_values,possible_AR,possible_delta_T_o = [],[],[],[]
population_size = 0
MAX_VELOCITY = [0.95,30,1,5]	# product of tip_mach_number,dT,n,AR
MIN_VELOCITY = [0.7,20,0,3]
popsize=50

for j in range(len(tip_mach_number)):
		for m in range(len(delta_T_o)):
			for k,n in enumerate(n_values):
				for l,AR in enumerate(aspect_ratio):
					filename = "data5/n%dAR%dc%d.xls" %(k,l,count)
					Ca,Cw1,Cw2,T02,dT = Values[j][m][13],Values[j][m][11],Values[j][m][12],Values[j][m][7],Values[j][m][1]
					blade_number = pop_generator(Ca,Cw1,Cw2,T02,dT,n,AR,filename,0) # dont want to save data thats why 0
					if blade_number != None : 
						blade_numbers.append(blade_number) 
						possible_tip_mach.append(tip_mach_number[j])
						possible_n_values.append(n)
						possible_AR.append(AR)
						possible_delta_T_o.append(delta_T_o[m])
						population_size +=1
					if population_size > popsize : break
				if population_size > popsize : break	
			count +=1
		if population_size >popsize : break

possible_sol  = np.empty ( (4,population_size))
print("Total population size is %d" %population_size)
for i in range(population_size):
	possible_sol[0][i] = possible_tip_mach[i]
	possible_sol[1][i] = possible_delta_T_o[i]
	possible_sol[2][i] = possible_n_values[i]
	possible_sol[3][i] = possible_AR[i]

def stator_diffusion_checker(alpha2):
	alpha2= degree_2_radian(alpha2)
	D = 1 - np.cos(alpha2) + 0.5*9*np.sin(alpha2)* ( (np.cos(alpha2))**2 - 0.433)
	if D < 0.6 : return(1)
	else : return(0)

def upstream_flow_at_tip(overall_pressure_ratio,operating_absolute_mach_number,dT,tip_mach_number):
	expression = ( 1+((gamma*gas_constant)/(8*Cp)) ) * (1/stage_efficiency) * \
						( ( overall_pressure_ratio)**((gamma-1)/gamma) -1)
	T1 = dT/expression	
	Values = np.zeros(8)	
	alpha1 = (1/np.sqrt(gamma*gas_constant*T1) ) * \
				( 125*np.pi - T1*( (tip_mach_number**2 -0.25)*( (gamma*gas_constant)/(125*np.pi)) ) )
	#print(alpha1)
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
	tip_diffusion = stator_diffusion_checker(radian_2_degree(alpha2))
	Values[0] = Ca
	Values[1] = Cw1
	Values[2] = Cw2
	Values[3] = T02
	Values[4] = dT
	Values[5] = radian_2_degree(beta1)
	Values[6] = radian_2_degree(beta2)
	Values[7] = tip_diffusion
	
	
	return(Values)


	

#returns at tip Ca,Cw1,Cw2,T02,delta_T_o,beta1,beta2,tip_diffusion

# find at mean beta1,beta2, solidity_factor,mean_spacing
# find at tip camber
# find T02 and T01



def blade_parameters_for_optimiser(Ca,Cw1,Cw2,T02,delta_T_o,n,AR,tip_beta1,tip_beta2):
	#blade = blade_profile()
	R_tip = tip_diameter / mean_diameter
	R_mean = 1
	R_hub = hub_radius / mean_radius
	a = (Cw1+Cw2)/(2*R_tip**n)
	b = ((Cw2-Cw1)*(R_tip))/2
	#print(a,b)
	mean_Cw1 = a*R_mean**n - b/R_mean
	mean_Cw2 = a*R_mean**n + b/R_mean
	hub_Cw2 = a*R_hub**n + b/R_hub
	#print( "mean CW1 %f ,CW2 %f" %(mean_Cw1,mean_Cw2) )
	mean_U = rotational_speed(design_speed,mean_radius*2)
	hub_U = rotational_speed(design_speed,hub_radius*2)
	mean_Rx =  Rx_calculator(mean_Cw1,mean_Cw2,mean_U)
	mean_Ca1 = Ca_upstream( Ca,Cw1,mean_Cw1,(tip_radius/mean_radius),R_mean,a,b,n)
	mean_Ca2 = Ca_downstream( Ca,Cw2,mean_Cw2,(tip_radius/mean_radius),R_mean,a,b,n)
	hub_Ca2 = Ca_downstream( Ca,Cw2,hub_Cw2,(tip_radius/mean_radius),R_hub,a,b,n)
	mean_alpha1 = radian_2_degree(math.atan(mean_Cw1/mean_Ca1) )
	mean_alpha2 = radian_2_degree( math.atan(mean_Cw2/mean_Ca2) )
	hub_alpha2 = radian_2_degree( math.atan(hub_Cw2/hub_Ca2 ))
	mean_C1 = np.sqrt(mean_Cw1**2 + mean_Ca1**2)
	mean_C2 = np.sqrt(mean_Cw2**2 + mean_Ca2**2)
	mean_beta1 = U_Ca_formula(mean_U, mean_Ca1 , mean_alpha1) 	#mean beta1
	mean_beta2 = U_Ca_formula(mean_U, mean_Ca2 , mean_alpha2)	#mean beta2
	mean_V1 = mean_Ca1 / np.cos(degree_2_radian(mean_beta1))
	mean_V2 = mean_Ca2 / np.cos(degree_2_radian(mean_beta2))
	# AT MEAN
	air_deflection_angle = mean_beta1 - mean_beta2
	air_outlet_angle = mean_beta2
	solidity_factor = pitch_chord_ration (mean_V1,mean_V2)		# from emperical data (s/c)

	chord_of_blade =  blade_height/AR
	spacing = solidity_factor * chord_of_blade
	number_of_blades = (2*np.pi*mean_radius)/spacing

	number_of_blades = int(number_of_blades) +1
	mean_spacing = (2*np.pi*mean_radius)/number_of_blades		#mean spacing 
	mean_solidity_factor = mean_spacing/chord_of_blade				# mean  solidity factor
	# AT tip TIP CAMBER
	tip_spacing = (2*np.pi*tip_radius)/number_of_blades
	tip_solidity_factor = tip_spacing / chord_of_blade
	tip_camber = ( tip_beta1 - tip_beta2 ) / \
										( 1- (0.25 + 0.1*(tip_beta2/50))* np.sqrt(tip_solidity_factor) )
	hub_diffusion = stator_diffusion_checker(hub_alpha2)
	mean_diffusion = stator_diffusion_checker (mean_alpha2)
	total_diffusion = hub_diffusion + mean_diffusion
	return (mean_beta1,mean_beta2,mean_spacing,mean_solidity_factor,tip_camber,total_diffusion)
#returns mean_beta1,mean_beta2,mean_spacing,mean_solidity_factor,tip_camber


MAX_EPOCHS = 10000
POPULATION = population_size


# store at mean beta1,beta2, solidity_factor,mean_spacing
# store at tip camber
# store T02 and T01

class blade_profile:
	
	def __init__(self):
		self.mean_radius = ( hub_radius + tip_radius ) /2
		self.mean_beta1 = 0
		self.mean_beta2 = 0
		self.mean_solidity_factor = 0
		self.mean_spacing = 0
		self.tip_camber = 0
		self.T02 = 0
		self.T01 = 0
		self.total_efficiency = 0
		self.pressure_rise = 0
		self.entropy_change = 0
		self.pBest = 0
		#self.pVelocity = 0
		self.pVelocity = np.zeros(4)
		#self.pVelocity = np.zeros(4) # tip_mach,dT,n,AR
		#unique to each blade design
		self.tip_mach_number = 0
		self.dT = 0
		self.n = 0
		self.AR = 0

	def get_pBest(self): 
		return (self.pBest)

	def set_pBest(self, value) : 
		self.pBest = value

	def get_velocity(self):
		return self.pVelocity

	def set_pVelocity(self, velocityScore):
		for i in range(4):
			self.pVelocity[i] = velocityScore[i]

	def set_design_parameters(self,tip_mach,dT,n,AR):
		self.tip_mach_number,self.dT,self.n,self.AR = tip_mach,dT,n,AR
	
	def set_efficiency(self):
		#loss calculation
		#At zero incidence angle CDp = 0.018
		CDp = 0.018 # profile dr
		# requirement mean beta1,beta2,solidity_factor,pacing all at mean.
		alpha1,alpha2 = degree_2_radian(self.mean_beta1),degree_2_radian(self.mean_beta2)
		alpha_m = math.atan( 0.5* ( math.tan(alpha1) + math.tan(alpha2)) ) 
		CL = 2*self.mean_solidity_factor*(math.tan(alpha1) - math.tan(alpha2)) - ( CDp * math.tan(alpha_m) ) 
		CDs = 0.018*CL**2 #secondary losses
		CDa = 0.020 * (self.mean_spacing/blade_height) # annulus drag losses s:spacin ; h:blade height
		CD = CDp + CDs + CDa # total drag losses
		#print(self.mean_solidity_factor,alpha_m,alpha1,alpha2)
		#print( self.mean_solidity_factor * ( (np.cos(alpha_m) )**3 / (np.cos(alpha1))**2 ) )
		loss_coeff = CD / ( self.mean_solidity_factor * ( ((np.cos(alpha_m) )**3) / (np.cos(alpha1))**2 ) )
		actual_loss_coeff = 1 - ( (np.cos(alpha1))**2 / (np.cos(alpha2))**2 )
		#print(loss_coeff,actual_loss_coeff)
		blade_efficiency = 1 - (loss_coeff/actual_loss_coeff)
		#stator and rotor efficeincy are assumed to be same efficiency since at mean Rx is nearly 0.5
		self.total_efficiency = blade_efficiency
		self.pressure_rise = ( 1+ (self.total_efficiency*(self.T02-self.T01)) /self.T01 )  ** (gamma/(gamma-1))
	
	def set_entropy_change(self):
		Cv = gas_constant / (gamma-1)
		self.entropy_change = ( Cv * np.log(self.T02/self.T01 ) - gas_constant * np.log(1.2) ) # pressure ratio = 1.2

def population_initialise(POPULATION,possible_sol):
	blade = np.empty(POPULATION,dtype=object)
	for i in range(POPULATION):
		blade[i] = blade_profile()
		tip_Ca,tip_Cw1,tip_Cw2,blade[i].T02,dT,tip_beta1,tip_beta2,tip_diffusion = upstream_flow_at_tip(overall_pressure_ratio,operating_absolute_mach_number,possible_sol[1][i],possible_sol[0][i])
		#print(tip_beta1,tip_beta2,possible_sol[1][i],possible_sol[0][i],possible_sol[2][i],possible_sol[3][i])
		blade[i].T01 = blade[i].T02 - dT
		blade[i].mean_beta1,blade[i].mean_beta2,blade[i].mean_spacing,blade[i].mean_solidity_factor,blade[i].tip_camber,total_diffusion = blade_parameters_for_optimiser(tip_Ca,tip_Cw1,tip_Cw2,T02,dT,possible_sol[2][i],possible_sol[3][i],tip_beta1,tip_beta2)
		#print(blade[i].mean_beta1,blade[i].mean_beta2,blade[i].mean_spacing,blade[i].mean_solidity_factor)
		blade[i].set_efficiency()
		blade[i].set_entropy_change()
		if blade[i].total_efficiency - stage_efficiency > 0 : 
			value = ( blade[i].total_efficiency - stage_efficiency ) * 100 
		else : 
			value = ( -blade[i].total_efficiency  + stage_efficiency ) * 100
		blade[i].set_pBest(value)
		Values = np.zeros(4)
		for z in range(4):
			Values[z] = possible_sol[z][i]
		blade[i].set_pVelocity( Values )
		blade[i].set_design_parameters(possible_sol[0][i],possible_sol[1][i],possible_sol[2][i],possible_sol[3][i])
		
	return(blade)

def find_gbest(blade):
	gpbest,gnbest = [],[]
	gp,gn = -1,-1
	for i in range(POPULATION):
		if blade[i].total_efficiency <=0.9:
			gnbest.append(0.9 - blade[i].total_efficiency)
		else :
			gpbest.append( blade[i].total_efficiency - 0.9 )
	if len(gpbest)!=0 : gp = min(gpbest)
	if len(gnbest)!=0 : gn = min(gnbest)
	if gp>= 0 and gn >=0 :
		gp,gn = min(gpbest),min(gnbest)
		if gp>gn : global_best = 0.9 - gn
		else : global_best = 0.9 + gp 
	elif gp < 0 : global_best = 0.9 - gn
	elif gn <0 : global_best = 0.9 + gp 
	return (global_best)
	
def test_problem(blade,index):	#returns levenshtein distance of two arrays
	efficiency_difference = blade[index].total_efficiency - stage_efficiency
	pressure_difference = blade[index].pressure_rise - overall_pressure_ratio
	return ( efficency_difference , pressure_difference )

def update_population(blade):
	phi = 4.1
	c1,c2 = 2.05,2.05
	K = 0.729
	pgd = find_gbest(blade)
	for i in range(POPULATION):
		new_blade = blade_profile()
		new_set = np.zeros(4)
		new_velocity_set = np.zeros(4)
		random_index = random.randint(0,3) # gives the index which implies which variable has to be perturbed from tip_mach_number,dT,n,AR 
		#print(random_index)
		if random_index == 0 :
			new_set[0] = random.uniform(0.7,0.95) # tip_mach is changed
			xid = new_set[0]
			pid = possible_sol[0][i]
			new_set[1],new_set[2],new_set[3] = possible_sol[1][i],possible_sol[2][i],possible_sol[3][i]
			pVelocity = K *( blade[i].pVelocity[0] + c1*random.random()*(pid-xid) + c2*random.random()*(pgd-xid) )
			new_velocity_set[0] = pVelocity
			new_velocity_set[1],new_velocity_set[2],new_velocity_set[3] = possible_sol[1][i],possible_sol[2][i],possible_sol[3][i]
		elif random_index == 1 :
			new_set[1] = random.uniform(20,30)		#dT is changed
			xid = new_set[1]
			pid = possible_sol[1][i]
			new_set[0],new_set[2],new_set[3] = possible_sol[0][i],possible_sol[2][i],possible_sol[3][i]
			pVelocity = K *( blade[i].pVelocity[1] + c1*random.random()*(pid-xid) + c2*random.random()*(pgd-xid) )
			new_velocity_set[1] = pVelocity
			new_velocity_set[0],new_velocity_set[2],new_velocity_set[3] = possible_sol[0][i],possible_sol[2][i],possible_sol[3][i]
		elif random_index == 2 :
			new_set[2] = random.uniform(0,1) 	# n is changed
			xid = new_set[2]
			pid = possible_sol[2][i]
			new_set[0],new_set[1],new_set[3] = possible_sol[0][i],possible_sol[1][i],possible_sol[3][i]
			pVelocity = K *( blade[i].pVelocity[2] + c1*random.random()*(pid-xid) + c2*random.random()*(pgd-xid) )
			new_velocity_set[2] = pVelocity
			new_velocity_set[0],new_velocity_set[1],new_velocity_set[3] = possible_sol[0][i],possible_sol[1][i],possible_sol[3][i]
		else : 
			new_set[3] = random.uniform(3,5)	# AR is changed
			xid = new_set[3]
			pid = possible_sol[3][i]
			new_set[0],new_set[1],new_set[2] = possible_sol[0][i],possible_sol[1][i],possible_sol[2][i]
			pVelocity = K *( blade[i].pVelocity[3] + c1*random.random()*(pid-xid) + c2*random.random()*(pgd-xid) )
			new_velocity_set[3] = pVelocity
			new_velocity_set[0],new_velocity_set[1],new_velocity_set[2] = possible_sol[0][i],possible_sol[1][i],possible_sol[2][i]
		
		if np.any(new_velocity_set > MAX_VELOCITY) : pass
		elif np.any(new_velocity_set < MIN_VELOCITY) : pass
		else : 
			tip_Ca,tip_Cw1,tip_Cw2,new_blade.T02,dT,tip_beta1,tip_beta2,tip_diffusion = upstream_flow_at_tip(overall_pressure_ratio,operating_absolute_mach_number,new_set[1],new_set[0])
			#print(tip_beta1,tip_beta2,possible_sol[1][i],possible_sol[0][i],possible_sol[2][i],possible_sol[3][i])
			new_blade.T01 = new_blade.T02 - dT
			new_blade.mean_beta1,new_blade.mean_beta2,new_blade.mean_spacing,new_blade.mean_solidity_factor,new_blade.tip_camber,total_diffusion = blade_parameters_for_optimiser(tip_Ca,tip_Cw1,tip_Cw2,T02,dT,new_set[2],new_set[3],tip_beta1,tip_beta2)
			if tip_diffusion + total_diffusion == 3 :
				print("random index is %d" %random_index)
				blade[i] = new_blade
				#print(blade[i].mean_beta1,blade[i].mean_beta2,blade[i].mean_spacing,blade[i].mean_solidity_factor)
				blade[i].set_efficiency()
				blade[i].set_entropy_change()
				if blade[i].total_efficiency - stage_efficiency > 0 : 
					value = ( blade[i].total_efficiency - stage_efficiency ) * 100 
				else : 
					value = ( -blade[i].total_efficiency  + stage_efficiency ) * 100
				blade[i].set_pBest(value)
				blade[i].set_pVelocity(new_velocity_set)
				blade[i].set_design_parameters(new_set[0],new_set[1],new_set[2],new_set[3])
				#time.sleep(0.1)
				print("new tip_mach %f,dT %f,n is %f,AR is %f" %(new_set[0],new_set[1],new_set[2],new_set[3]))
			else : 
				pass
	
	return (blade)
	

breaking_point = 100

def PSO_algorthim():
	epoch = 0
	done = False
	Blades = population_initialise(POPULATION,possible_sol)
	breaker = 0
	parameters = np.zeros((7,breaking_point))
	while not done:
		# the loop will break in following cases
		# Case(1). Total number of iterations are as defined in initial
		# Case(2). Efficiency and pressure ratio are met
		# Case(3). Optional Entropy Change and Camber is approching to zero
		if epoch < MAX_EPOCHS:
			
			for j in range(POPULATION):
				if (0.9-0.005) <= Blades[j].total_efficiency <= (0.9+0.005):
					print("efficiency is %f,pressure rise %f,tip_camber %f \n" %(Blades[j].total_efficiency,Blades[j].pressure_rise,Blades[j].tip_camber))
					#print("tipmach %f ,dT %f ,n %f ,AR is %f " % (Blades[j].tip_mach_number,Blades[j].dT,Blades[j].n,Blades[j].AR) )
					#print(type(parameters[0][1]))
					parameters[0][breaker] =  Blades[j].total_efficiency
					parameters[1][breaker] =  Blades[j].pressure_rise
					parameters[2][breaker] =  Blades[j].tip_camber
					parameters[3][breaker] =  Blades[j].tip_mach_number
					parameters[4][breaker] =  Blades[j].dT
					parameters[5][breaker] =  Blades[j].n
					parameters[6][breaker] =  Blades[j].AR
					
					breaker += 1
					#old_tip_camber.append(Blades[j].tip_camber)
			Blades = update_population(Blades)
			#for j in range(POPULATION):
			#	print("tipmach %f ,dT %f ,n %f ,AR is %f " % (Blades[j].tip_mach_number,Blades[j].dT,Blades[j].n,Blades[j].AR) )
			print("Iteration_number %d" %epoch)	
			epoch  += 1
			if breaker ==breaking_point :
				#print(parameters)
				print("YO")
				return (parameters)
				break
		else:
			done = True
			return(parameters)

	return

parameters = PSO_algorthim()
header = ['total_efficiency','pressure_rise','tip_camber','tip_mach_number','dT','n','AR']
# saving data in excel sheet
def save_final_data(parameters,header,filename):
	book = xlwt.Workbook(encoding="utf-8")
	sheet1 = book.add_sheet("Sheet 1")
	total_parameters = len(header)
	for i in range(total_parameters):
		sheet1.write(0,i,header[i])
	
	for j in range(1,(breaking_point+1)):
		sheet1.write(j,0,parameters[0][j-1])
		sheet1.write(j,1,parameters[1][j-1])
		sheet1.write(j,2,parameters[2][j-1])
		sheet1.write(j,3,parameters[3][j-1])
		sheet1.write(j,4,parameters[4][j-1])
		sheet1.write(j,5,parameters[5][j-1])
		sheet1.write(j,6,parameters[6][j-1])
	
	try :
		os.remove(filename)
	except OSError:
		pass
	book.save(filename)
	print("Done Enjoy!!")

filename = "data.xls"
save_final_data(parameters,header,filename)

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
	
	Ca,Cw1,Cw2,T02,delat_T_o,tip_beta1,tip_beta2,tip_diffusion = upstream_flow_at_tip(overall_pressure_ratio,operating_absolute_mach_number,dT,tip_mach_number)
	
	number_of_blades = pop_generator(Ca,Cw1,Cw2,T02,delta_T_o,n,AR,filename,1)
	print("saved")
	
final_blade_generation("final_blade.xls")




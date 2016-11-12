#author Vishwas Verma
# Nov 05,2016
import csv
import numpy as np
import time
import os.path

def read_data(filename):
	data = np.empty(51,dtype=object )
	data1,data2 = [],[]
	with open(filename) as csvfile:
		readCSV = csv.reader(csvfile, delimiter=',')
		for i,row in enumerate(readCSV):
			#print(row[0])
			#print(i)
			if i < 26 : data1.append(row[0])#data[i][1]=row[0]#
			if i > 25 and i<52 : data2.append(row[0])#data[j][2] = row[0]	#
	print(len(data1),len(data2))
	count = len(data2)
	for i,d in enumerate(data1) :
		d = d +" " + '0'
		data[i] = d
		#print(i)
	for i,d in enumerate(data2) :
		if i >0:
			d = d + " " + '0'
			data[(count+25)-i] = d
	return data

fullpath = "/home/vishwas/Desktop/Aircraft Engine/AE310/modifying_velocity_function/to_submit/stator_blades/"
path = input("Enter path of file :")
file = "naca.gnu"
filename = fullpath + path+ '/' + file
print(filename)
data = read_data(filename)
savefile = fullpath + path + '/' + 'text.txt'
try :
		os.remove(savefile)
except OSError:
		pass
for d in data:
	#print(d)
	with open(savefile,'a') as out :
		out.write(d + '\n')
print(len(data))
		#with open('test.csv','w') as fp:
#	np.savetxt(fp,data,delimiter=' ') 

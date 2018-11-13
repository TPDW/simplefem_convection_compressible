#this plots the comparisons of the final temeperature view
#over the different approximations and parameters


import os
import numpy as np
import matplotlib.pyplot as plt

nnx=16

Ra_list = np.array([4,5])

Ra_dictionary = {4:0,5:1}

Di_list = [0.25,0.5,1.0]

approximation_dict = {1:"ALA",2:"TALA",3:"BA",4:"EBA"}

#iterate over the allowed combinations of values
#find a good way of storing the arrays
#maybe a multidimensional array?
#that could work, might be a bit weird though
#how else could I do it?
#should I care that much about how I do it?
#yes, I want this to be extendable to other sets of values
#and to figure out this problem for it's own sake
#and for the sake of learning

Nu_list = []
Temp_list = []
Vrms_list = []

for approximation in range(1,5):
	Nu_list.append([])
	Temp_list.append([])
	Vrms_list.append([])
	for Ra in Ra_list:
		if(approximation == 3):
			filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(nnx) + "/Nu.dat"
			Nu_list[2].append(np.loadtxt(filename))
			filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(nnx) + "/temp.dat"
			Temp_list[2].append(np.loadtxt(filename))
			filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(nnx) + "/vrms.dat"
			Vrms_list[2].append(np.loadtxt(filename))

		else:
			Nu_list[approximation-1].append([])
			Temp_list[approximation-1].append([])
			Vrms_list[approximation-1].append([])
			for Di in Di_list:
				filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(Di)+ "_" + str(nnx)  +"/Nu.dat"
				Nu_list[approximation-1][Ra_dictionary[Ra]].append(np.loadtxt(filename))
				filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(Di)+ "_" + str(nnx)  +"/temp.dat"
				Temp_list[approximation-1][Ra_dictionary[Ra]].append(np.loadtxt(filename))
				filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(Di)+ "_" + str(nnx)  +"/vrms.dat"
				Vrms_list[approximation-1][Ra_dictionary[Ra]].append(np.loadtxt(filename))


#so that's the data reading and storing done
#here plot the data
#see notes for plan of how these plots should look

nrows = 4 #ALA,TALA,BA,EBA

ncols = len(Ra_list)

f_Nu,   ax_Nu   = plt.subplots(nrows,ncols,figsize=(10,10))
f_vrms, ax_vrms = plt.subplots(nrows,ncols,figsize=(10,10))
f_temp, ax_temp = plt.subplots(nrows,ncols,figsize=(10,10))

for Ra in Ra_list:
	for approximation in range(1,5):
		if approximation==3:
			ax_temp[approximation-1][Ra_dictionary[Ra]].plot(Temp_list[approximation-1][Ra_dictionary[Ra]][:,0],Temp_list[approximation-1][Ra_dictionary[Ra]][:,1])
			ax_vrms[approximation-1][Ra_dictionary[Ra]].plot(Vrms_list[approximation-1][Ra_dictionary[Ra]][:,0],Vrms_list[approximation-1][Ra_dictionary[Ra]][:,1])
			ax_Nu[approximation-1][Ra_dictionary[Ra]].plot(Nu_list[approximation-1][Ra_dictionary[Ra]][:,0],Nu_list[approximation-1][Ra_dictionary[Ra]][:,1])
			
			ax_temp[approximation-1][Ra_dictionary[Ra]].set_xlabel("Time")
			ax_vrms[approximation-1][Ra_dictionary[Ra]].set_xlabel("Time")
			ax_Nu[approximation-1][Ra_dictionary[Ra]].set_xlabel("Time")

			ax_temp[approximation-1][Ra_dictionary[Ra]].set_ylabel("Temperature")
			ax_vrms[approximation-1][Ra_dictionary[Ra]].set_ylabel("$V_{rms}$")
			ax_Nu[approximation-1][Ra_dictionary[Ra]].set_ylabel("Nu")



		else:
			for Di_number in range(len(Di_list)):
				Di = Di_list[Di_number] #There has to be a better way to do this bit	
				ax_temp[approximation-1][Ra_dictionary[Ra]].plot(Temp_list[approximation-1][Ra_dictionary[Ra]][Di_number][:,0],Temp_list[approximation-1][Ra_dictionary[Ra]][Di_number][:,1],label="Di = {}".format(Di))
				ax_vrms[approximation-1][Ra_dictionary[Ra]].plot(Vrms_list[approximation-1][Ra_dictionary[Ra]][Di_number][:,0],Vrms_list[approximation-1][Ra_dictionary[Ra]][Di_number][:,1],label="Di = {}".format(Di))
				ax_Nu[approximation-1][Ra_dictionary[Ra]].plot(Nu_list[approximation-1][Ra_dictionary[Ra]][Di_number][:,0],Nu_list[approximation-1][Ra_dictionary[Ra]][Di_number][:,1],label="Di = {}".format(Di))
			ax_temp[approximation-1][Ra_dictionary[Ra]].legend()
			ax_vrms[approximation-1][Ra_dictionary[Ra]].legend()
			ax_Nu[approximation-1][Ra_dictionary[Ra]].legend()

			ax_temp[approximation-1][Ra_dictionary[Ra]].set_xlabel("Time")
			ax_vrms[approximation-1][Ra_dictionary[Ra]].set_xlabel("Time")
			ax_Nu[approximation-1][Ra_dictionary[Ra]].set_xlabel("Time")

			ax_temp[approximation-1][Ra_dictionary[Ra]].set_ylabel("Temperature")
			ax_vrms[approximation-1][Ra_dictionary[Ra]].set_ylabel("$V_{rms}$")
			ax_Nu[approximation-1][Ra_dictionary[Ra]].set_ylabel("Nu")





f_temp.tight_layout()
f_temp.subplots_adjust(left=0.175,top=0.9)

f_vrms.tight_layout()
f_vrms.subplots_adjust(left=0.175,top=0.9)

f_Nu.tight_layout()
f_Nu.subplots_adjust(left=0.175,top=0.9)

for approximation in range(1,5):
	f_temp.text(0.025,1.025-0.225*approximation,approximation_dict[approximation],fontsize=20)
	f_vrms.text(0.025,1.025-0.225*approximation,approximation_dict[approximation],fontsize=20)
	f_Nu.text(0.025,1.025-0.225*approximation,approximation_dict[approximation],fontsize=20)

for Ra in Ra_list:
	string = "Ra = $10^"+str(Ra)+"$"
	f_temp.text(0.15+0.45*(Ra_dictionary[Ra]+0.5),0.925,string,fontsize=20,horizontalalignment='center')
	f_vrms.text(0.15+0.45*(Ra_dictionary[Ra]+0.5),0.925,string,fontsize=20,horizontalalignment='center')
	f_Nu.text(0.15+0.45*(Ra_dictionary[Ra]+0.5),0.925,string,fontsize=20,horizontalalignment='center')


f_temp.savefig("temp.pdf")
f_vrms.savefig("vrms.pdf")
f_Nu.savefig("Nu.pdf")


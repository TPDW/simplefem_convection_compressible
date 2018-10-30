#this plots the comparisons of the final temeperature view
#over the different approximations and parameters


import os
import numpy as np
import matplotlib.pyplot as plt

nnx_list=np.array([16])


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

list_of_arrays = []

for approximation in range(1,5):
	list_of_arrays.append([])
	for Ra in Ra_list:
		if(approximation == 3):
			filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "/temp_final.dat"
			list_of_arrays[2].append(np.loadtxt(filename))
		else:
			list_of_arrays[approximation-1].append([])
			for Di in Di_list:
				filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(Di) + "/temp_final.dat"
				list_of_arrays[approximation-1][Ra_dictionary[Ra]].append(np.loadtxt(filename))


#so that's the data reading and storing done
#here plot the data
#see notes for plan of how these plots should look

nrows = 4 #ALA,TALA,BA,EBA

ncols = len(Ra_list)

for Di_number in range(len(Di_list)):
	f, ax = plt.subplots(nrows,ncols,figsize=(10,10))
	Di = Di_list[Di_number] #There has to be a better way to do this bit
	for approximation in range(1,5):
		for Ra in Ra_list:
			if approximation==3:
				#ax[approximation-1][Ra_dictionary[Ra]].axis('off')
				# ax[approximation-1][Ra_dictionary[Ra]].contour (list_of_arrays[approximation-1][Ra_dictionary[Ra]],cmap='coolwarm')
				ax[approximation-1][Ra_dictionary[Ra]].contourf(list_of_arrays[approximation-1][Ra_dictionary[Ra]],cmap='coolwarm')	
				#ax[approximation-1][Ra_dictionary[Ra]].spines['bottom'].set_visible(True)
				ax[approximation-1][Ra_dictionary[Ra]].tick_params(
			    axis='both',          # changes apply to the x-axis
			    which='both',      # both major and minor ticks are affected
			    bottom=False,      # ticks along the bottom edge are off
			    top=False,         # ticks along the top edge are off
			    left=False,
			    right=False,
			    labelbottom=False,
			    labelleft=False) # labels along the bottom edge are off
			else:	
				# ax[approximation-1][Ra_dictionary[Ra]].contour (list_of_arrays[approximation-1][Ra_dictionary[Ra]][Di_number],cmap='coolwarm')
				ax[approximation-1][Ra_dictionary[Ra]].contourf(list_of_arrays[approximation-1][Ra_dictionary[Ra]][Di_number],cmap='coolwarm')
				#ax[approximation-1][Ra_dictionary[Ra]].axis('off')
				ax[approximation-1][Ra_dictionary[Ra]].tick_params(
			    axis='both',          # changes apply to the x-axis
			    which='both',      # both major and minor ticks are affected
			    bottom=False,      # ticks along the bottom edge are off
			    top=False,         # ticks along the top edge are off
			    left=False,
			    right=False,
			    labelbottom=False,
			    labelleft=False) # labels along the bottom edge are off


	cs = ax[2][0].contourf(list_of_arrays[2][0],cmap='coolwarm')

	f.subplots_adjust(left=0.1,right=0.8)
	cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
	f.colorbar(cs, cax=cbar_ax)


	for approximation in range(1,5):
		f.text(0.025,1-0.2*approximation,approximation_dict[approximation],fontsize=20)

	for Ra in Ra_list:
		string = "Ra = $10^"+str(Ra)+"$"
		f.text(0.1+0.35*(Ra_dictionary[Ra]+0.5),0.9,string,fontsize=20,horizontalalignment='center')

	string = "Di = " + str(Di)
	f.text(0.025,0.9,string,fontsize=20)

	string = "Di_=_" + str(Di) + ".png"
	f.savefig(string)

import os
import numpy as np

nnx_list=np.array([32])


Ra_list = np.array([6])


Di_list = [0.25,0.5,1.0]

approximation_dict = {1:"ALA",2:"TALA",3:"BA",4:"EBA"}

approximation = 3
for nnx in nnx_list:
	for Ra in Ra_list:
		os.mkdir("OUT")
		os.mkdir("OUT/Paraview")
		string_to_run = "./simplefem " + str(nnx) + " " + str(approximation) + " " + str(Ra) + " 0.5 1 0"
		os.system(string_to_run)
		string_to_move = "mv OUT OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(nnx)
		os.system(string_to_move)




# for approximation in [1,2,4]:
# 	for nnx in nnx_list:
# 		for Ra in Ra_list:
# 			for Di in Di_list:
# 				os.mkdir("OUT")
# 				os.mkdir("OUT/Paraview")
# 				string_to_run = "./simplefem " + str(nnx) + " " + str(approximation) + " " + str(Ra) + " " + str(Di) + " 1 0"
# 				os.system(string_to_run)
# 				string_to_move = "mv OUT OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(Di) + "_" + str(nnx)
# 				os.system(string_to_move)



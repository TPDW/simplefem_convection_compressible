import os
import numpy as np

nnx_list=np.array([16])


Ra_list = np.array([4,5])
Di_list = np.array([0.25,0.5,1.0])
approximation = "1"

for nnx in nnx_list:
	for Ra in Ra_list:
		for Di in Di_list:
			print(approximation+str(nnx)+str(Ra)+str(Di))
			string_to_run = "./simplefem " + str(nnx) + " " + approximation + " " + str(Ra) +" " + str(Di) + " 0 0 >nul"
			print(string_to_run)
			os.system(string_to_run)
			string_to_move = "mv OUT/vrms.dat vrms_" + str(nnx) + "_ALA_" + str(Ra) + " " + str(Di) +".dat"
			os.system(string_to_move)
			string_to_move = "mv OUT/Nu.dat Nu_" + str(nnx) + "_ALA_" + str(Ra) + " " + str(Di) +".dat"
			os.system(string_to_move)
			string_to_move = "mv OUT/temp_final.dat temp_final" + str(nnx) + "_ALA_" + str(Ra) + " " + str(Di) +".dat"
			os.system(string_to_move)

approximation = "2"
for nnx in nnx_list:
	for Ra in Ra_list:
		for Di in Di_list:
			print(approximation+str(nnx)+str(Ra))
			string_to_run = "./simplefem " + str(nnx) + " " + approximation + " " + str(Ra) +" " + str(Di) + " 0 0 >nul"
			print(string_to_run)
			os.system(string_to_run)
			string_to_move = "mv OUT/vrms.dat vrms_" + str(nnx) + "_TALA_" + str(Ra) + " " + str(Di) +".dat"
			os.system(string_to_move)
			string_to_move = "mv OUT/Nu.dat Nu_" + str(nnx) + "_TALA_" + str(Ra) + " " + str(Di) +".dat"
			os.system(string_to_move)
			string_to_move = "mv OUT/temp_final.dat temp_final" + str(nnx) + "_TALA_" + str(Ra) + " " + str(Di) +".dat"
			os.system(string_to_move)

approximation = "4"
for nnx in nnx_list:
	for Ra in Ra_list:
		for Di in Di_list:
			print(approximation+str(nnx)+str(Ra))
			string_to_run = "./simplefem " + str(nnx) + " " + approximation + " " + str(Ra) +" " + str(Di) + " 0 0 >nul"
			print(string_to_run)
			os.system(string_to_run)
			string_to_move = "mv OUT/vrms.dat vrms_" + str(nnx) + "_EBA_" + str(Ra) + " " + str(Di) +".dat"
			os.system(string_to_move)
			string_to_move = "mv OUT/Nu.dat Nu_" + str(nnx) + "_EBA_" + str(Ra) + " " + str(Di) +".dat"
			os.system(string_to_move)
			string_to_move = "mv OUT/temp_final.dat temp_final" + str(nnx) + "_EBA_" + str(Ra) + " " + str(Di) +".dat"
			os.system(string_to_move)



approximation = "3"
Di=0.5
for nnx in nnx_list:
	for Ra in Ra_list:
		string_to_run = "./simplefem " + str(nnx) + " " + approximation + " " + str(Ra) +" " + str(Di) + " 0 0 >nul"
		print(string_to_run)
		os.system(string_to_run)
		string_to_move = "mv OUT/vrms.dat vrms_" + str(nnx) + "_BA_" + "4" + ".dat"
		os.system(string_to_move)
		string_to_move = "mv OUT/Nu.dat Nu_" + str(nnx) + "_BA_" + "4" + ".dat"
		os.system(string_to_move)
		string_to_move = "mv OUT/temp_final.dat temp_final" + str(nnx) + "_BA_" + "4" + ".dat"
		os.system(string_to_move)

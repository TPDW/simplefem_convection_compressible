import os
import numpy as np
import io

nnx_list=np.array([16])


Ra_list = np.array([4,5])

Ra_dictionary = {4:0,5:1}

Di_list = [0.25,0.5,1.0]

Di_dict = {0.25:0,0.5:1,1.0:2}

approximation_dict = {1:"ALA",2:"TALA",3:"BA",4:"EBA"}



#wait what am I doing?


file=io.open("stats_tables.txt","w")


BA_Nu = [4.8844,10.5340,21.9720]
BA_Vrms = [42.8650,193.2150,833.990]
BA_Temp = [0.500,0.500,0.500]


#Ok, I need to figure out a way to make these comprehensible
#How do I even do that?
#Maybe with dictionaries?
#or some kind of data structure
#now I know why compscis study data structures


#Just document in the comments
#Loop over Di, then Ra
EBA_Nu = [[4.09,3.38,2.19],[8.62,6.87,3.96],[16.47]]
EBA_Vrms = [[38.4,33.4,24.2],[173.5,152.5,107.7],[658.9]]
EBA_Temp = [[0.491,0.482,0.467],[0.504,0.502,0.482],0.519]

ALA_Nu = [[4.41,3.8,2.47],[9.2,7.5,3.88]]
ALA_Vrms = [[40.0,36.0,25],[179,155,85]]
ALA_Temp = [[0.515,0.522,0.512],[0.532,0.547,0.529]]

TALA_Nu = [[4.41,3.86,2.57],[9.22,7.61,3.93]]
TALA_Vrms = [[40.0,36.3,26.0],[178.0,155.9,84.9]]
TALA_Temp = [[0.513,0.519,0.509],[0.530,0.545,0.529]]


Nu_Array   = [ALA_Nu,TALA_Nu,BA_Nu,EBA_Nu]
Vrms_Array = [ALA_Vrms,TALA_Vrms,BA_Vrms,EBA_Vrms]
Temp_Array = [ALA_Temp,TALA_Temp,BA_Temp,EBA_Temp]




#TODO: Format the strings properly using the new {} formatting things
#So that we don't have stupid numbers of decimal places

for approximation in range(1,5):
	file.write("\\begin{table} \n")
	file.write("\\caption{"+approximation_dict[approximation]+"} \n")
	file.write("\\label{tab:"+approximation_dict[approximation]+"} \n")
	file.write("\\begin{center} \n")

	if(approximation == 3):
		file.write("\\begin{tabular}{|c|c|c|c|c|c|c|} \n")
		file.write("\\hline \n")
		file.write(" & \\multicolumn{6}{c|}{Code} \\\ \\hline \n")
		file.write("& \\multicolumn{3}{c|}{King et al.} & \\multicolumn{3}{c|}{SimpleFEM} \\\ \\hline \n")
		file.write("Ra & Nu & $V_{rms}$ & \\textless T \\textgreater & Nu & $V_{rms}$ & \\textless T \\textgreater \\\ \\hline \\hline \n")


		for Ra in Ra_list:
			filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "/temp_final.dat"
			mean_temp = np.mean(np.loadtxt(filename))
			King_Temp = Temp_Array[approximation-1][Ra_dictionary[Ra]]
			filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "/Nu.dat"
			Nu = np.loadtxt(filename)[-1,1]
			King_Nu = Nu_Array[approximation-1][Ra_dictionary[Ra]]
			filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "/vrms.dat"
			vrms = np.loadtxt(filename)[-1,1]
			King_Vrms = Vrms_Array[approximation-1][Ra_dictionary[Ra]]
			string_to_write = "$10^" + str(Ra) + "$ & " +str(King_Nu) +"& "+ str(King_Vrms)+"& "+ str(King_Temp)+"&" 
			#string_to_write +=  str(Nu) + "&" + str(vrms) + "&" + str(mean_temp) + "\n"
			string_to_write += "{:1.4f}&{:1.4f}&{:1.4f} \\\ \\hline \n".format(Nu,vrms,mean_temp)
			file.write(string_to_write)



	else:
		file.write("\\begin{tabular}{|c|c|c|c|c|c|c|c|} \n")
		file.write("\\hline \n")
		file.write("& & \\multicolumn{6}{c|}{Code} \\\ \\hline \n")
		file.write("& & \\multicolumn{3}{c|}{King et al.} & \\multicolumn{3}{c|}{SimpleFEM} \\\ \\hline \n")
		file.write("Ra & Di & Nu & $V_{rms}$ & \\textless T \\textgreater & Nu & $V_{rms}$ & \\textless T \\textgreater\\\ \\hline \\hline \n")
		for Ra in Ra_list:
			for Di in Di_list:
				filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(Di) + "/temp_final.dat"
				King_Temp = Temp_Array[approximation-1][Ra_dictionary[Ra]][Di_dict[Di]]
				mean_temp = np.mean(np.loadtxt(filename))
				filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(Di) + "/Nu.dat"
				Nu = np.loadtxt(filename)[-1,1]
				King_Nu = Nu_Array[approximation-1][Ra_dictionary[Ra]][Di_dict[Di]]
				filename = "OUT_" + approximation_dict[approximation] + "_" +  str(Ra) + "_" + str(Di) + "/vrms.dat"
				vrms = np.loadtxt(filename)[-1,1]
				King_Vrms = Vrms_Array[approximation-1][Ra_dictionary[Ra]][Di_dict[Di]]
				string_to_write = "$10^" + str(Ra) + "$&" + str(Di) +" & " +str(King_Nu) +"& "+ str(King_Vrms)+"& "+ str(King_Temp)+"&" 
				#string_to_write +=  str(Nu) + "&" + str(vrms) + "&" + str(mean_temp) + "\n"
				string_to_write += "{:1.4f}&{:1.4f}&{:1.4f} \\\ \\hline \n".format(Nu,vrms,mean_temp)
				file.write(string_to_write)


	file.write("\\end{tabular}\n")
	file.write("\\end{center}\n")
	file.write("\\end{table}\n")
	file.write("\n")


file.close()


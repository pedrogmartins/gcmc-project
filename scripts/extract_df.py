import os
import pandas as pd
import feather
import sys

#This script takes in a directory of LAMMPS simulations (argument #1)
#   and extracts the ouput for a GCCM calculation into a pandas data frame
#   saved as a compressed feather file in the /results directory with the
#   original simulation directory name.

#Input is folder name
folder = sys.argv[1]

os.chdir("./" + folder)
dir_path = os.getcwd()
print(dir_path)

#Loop over all simulations directories in folder
for subdir in os.listdir(dir_path):

    # get the full path of the subdirectory
    subdir_path = os.path.join(dir_path, subdir)
    # check if the path is a directory
    if os.path.isdir(subdir_path):
        # get the log file path
        log_path = os.path.join(subdir_path, "log.lammps")

        # check if the log file exists
        if os.path.exists(log_path):

            print(log_path)

            #Booleans to check progress
            found_top = False
            found_bottom = False
            steps_line = -1

            with open(log_path) as file:
                #Loop over lines to find initial and ending lines
                for i, line in enumerate(file):
                    if found_top == False:
                        if 'Step Temp' in line:
                            steps_line = i
                            found_top = True
                    if found_bottom == False:
                        if 'Loop time of' in line:
                            loop_line = i
                            found_bottom = True

            #Read extracted results as a pandas data frame
            if found_bottom == True:
                df = pd.read_csv(log_path, delim_whitespace=True, skiprows=steps_line, skipfooter=i-loop_line+1, comment='#', engine='python')
            else:
                df = pd.read_csv(log_path, delim_whitespace=True, skiprows=steps_line, skipfooter=1, comment='#', engine='python')

            feather.write_dataframe(df, "../results/" + subdir + ".feather")


#increase amount of data allowed on notebook if needed when going through
#    hundreds of directories. 
#jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000

#!/bin/bash

#This bash script creates atomatically LAMMPS GCMC input files from
#   a .txt file with inputs/parameters separated by " ", one simulations
#   per line. It populates a blank sample script lammps_inputs_p.sh with
#   all files, creates a directory for this simulation and submits the
#   slurm job.

while read p; do
   # Extracting individual fields
   field1=$(echo $p | cut -d' ' -f1)
   field2=$(echo $p | cut -d' ' -f2)
   field3=$(echo $p | cut -d' ' -f3)
   field4=$(echo $p | cut -d' ' -f4)
   field5=$(echo $p | cut -d' ' -f5)
   field6=$(echo $p | cut -d' ' -f6)
   field7=$(echo $p | cut -d' ' -f7)

   # Constructing directory name
   directory="all_183k_313k_${field1}_p${field2}_s${field3}_m${field4}_a${field5}_d${field6}_oc${field7}"

   # Creating directory
   mkdir "$directory"

   # Changing directory
   cd "$directory"

   # Running script with arguments
   bash ../lammps_inputs_p.sh $field1 $field2 $field3 $field4 $field5 $field6 $field7

   # Copying script
   cp ../run_lammps.sh ./

   # Submitting job
   sbatch run_lammps.sh

   # Returning to parent directory
   cd ..
done < inputs

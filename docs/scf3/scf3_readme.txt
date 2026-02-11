D(r) files are included for ScF3 with the format ScF3_xxxK_dofr.xye,
where xxx is the temperature the data was collected at.

Temperatures: 125, 140, 147, 152, 175, 200, 225, 250, 275, 325, 350, 375, 400, 425, 450 K

To run the "batch_run_scf3.py" script:
Save the script, along with the input file for scf3 and the scf3 isodistort cif in the Topas 6 directory.
This script uses the read_isodisplace_cif method, so copy that to your topas directory.
NOTE: The topas-6 directory should not have any spaces.
Save the data to a directory, again without any spaces.
Edit the input files to reflect the location of the data files.
Execute the script, e.g. by typing "python batch_run_scf3.py" in the console
Output files will be created in the topas directory for each irrep and each temperature.
After the script has finished executing, copy the output files to another directory for data analysis.



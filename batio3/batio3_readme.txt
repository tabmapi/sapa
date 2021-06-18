D(r) files are included for BaTiO3 with the format BaTiO3_xxxK.xy,
where xxx is the temperature the data was collected at.

Temperatures: 15, 150, 210, 250, 293, 350, 410, 500 K

To run the "batch_run_batio3.py" script:
Save the script, along with the input files for batio3 in the Topas 6 directory.
NOTE: The topas-6 directory should not have any spaces.
Save the data to a directory, again without any spaces.
Edit the input files to reflect the location of the data files.
Execute the script, e.g. by typing "python batch_run_batio3.py" in the console
Output files will be created in the topas directory for each irrep and each temperature.
After the script has finished executing, copy the output files to another directory for data analysis.


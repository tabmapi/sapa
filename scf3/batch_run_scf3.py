import os
import read_isodisplace_cif as r

##---------ScF3---------
#Run this file in the Topas directory
####----NOTE - EDIT INP FILES TO POINT TO CORRECT DATA LOCATION
temps = ["125","140","147","152","175","200","225","250","275","325","350","375","400","425","450"]
struc = r.iso_cif_file("scf3_alt_isodistort_v2.cif")
irreps = struc.irrep_list()
if ("GaM4-" in irreps):
  irreps.remove("GaM4-")
for i in irreps:
  for temp in temps:
    os.system("tc batch_modes_scf3.inp \" macro IRREP {%s} macro VAR {%s} #define %s \" "%(i,temp,i))
    
  

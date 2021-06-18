import os
import read_isodisplace_cif as r

####---separate inp files are included for 
####---different space groups
####---each outputs files with the same name format for ease of analysis

####----Run in Topas directory
####----NOTE - EDIT INP FILES TO POINT TO CORRECT DATA LOCATION


irreps = ["GM4-","GM5-","R3-","R4-","R5+","R5-","X1+","X2+","X3-","X5+","X5-","M1+","M2+","M2-","M3+","M4+","M5+","M5-"]
####----Rhombohedral-----
temps = ["15","150"]
for temp in temps:
  for i in irreps:
    os.system("tc batch_modes_batio3_rhomb.inp \" macro IRREP {%s} macro VAR {%s} #define %s \" "%(i,temp,i))
    
####----Orthorhombic
temps = ["210","250"]
for temp in temps:
  for i in irreps:
    os.system("tc batch_modes_batio3_ortho.inp \" macro IRREP {%s} macro VAR {%s} #define %s \" "%(i,temp,i))

####----Tetragonal
temps = ["293","350"]
for temp in temps:
  for i in irreps:
    os.system("tc batch_modes_batio3_tetragonal.inp \" macro IRREP {%s} macro VAR {%s} #define %s \" "%(i,temp,i))

####----Cubic
temps = ["410","500"]
for temp in temps:
  for i in irreps:
    os.system("tc batch_modes_batio3_cubic.inp \" macro IRREP {%s} macro VAR {%s} #define %s \" "%(i,temp,i))


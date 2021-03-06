# only takes single phase cif from isodisplace.
# blank space before cif file must be striped
import numpy as np
import pandas as pd
import re
# the class which can contain any cif file
# IMPORTANT - the cif file must have a line at the end starting with a # symbol
# This is included by default in cifs generated by ISODISTORT
# This is NOT included by default in cifs generated by TOPAS
class iso_cif_file:
    def __init__(self, cif_file):

        f = open(cif_file, "r")
        lines = f.readlines()
        count = 0
        count_cifloop = 0
        for n in range(0, len(lines)):
            line = lines[n]
            # reads single value entry in cif file and creats a class attribute with value stored as str
            if line[0] == "_" and len(line.split()) == 2:
                vars(self)['%s' % (line.split()[0].lstrip(
                    "_").strip())] = line.split()[1]

            # reads single line entries in cif file where data is on the next line
            # if line[0] == "_" and len(line.split())==1 and r'[\t]' not in line and lines[n+1][0] != '_':
            #   vars(self)['%s'%(line.lstrip("_").rstrip("\n").strip())] = lines[n+1].strip().sub(r'[\t]', '').rstrip("\n")

            # reads multiline entries in table format....
            if line[0:5] == "loop_":
                n = n + 1
                m = 0
                while "_" in lines[n+m][0]:
                    vars(self)['%s' %
                               (lines[n+m].lstrip("_").rstrip("\n"))] = []
                    m = m + 1
                    # work out out many columns the mutiline entreis have....
                q = 0
                while len(lines[n+m+q]) > 2 and "loop_" not in lines[n+m+q] and "#" not in lines[n+m+q]:
                    k = 0

                    for k in range(m):
                        # each multiline becomesa. a attribute of the class with the associated value set.
                        try:
                            vars(self)['%s' % (
                                lines[n+k].lstrip("_").rstrip("\n").strip())].append(lines[n+m+q].split()[k])
                        except:
                            print("cif file error in %s at line %d" %
                                  (cif_file, n))
                    q = q+1
                    count_cifloop = q + m

            # interate forwards so that I don't re-read the cif.
            n = n + count_cifloop

        f.close()

    # converts any matricies in cif loop format in matrix foramt
    def make_matrices(self):
        for class_attribute in dir(self):
            if "_row" in class_attribute:
                rows = np.array(vars(self)['%s' %
                                           class_attribute]).astype(np.int)
                cols = np.array(vars(self)['%s_col' %
                                           class_attribute[:-4]]).astype(np.int)
                values = vars(self)['%s_value' % class_attribute[:-4]]
                cif_loop_matrix = np.zeros([max(rows), max(cols)])
                # print rows
                # print cols
                for val_index in range(len(values)):
                    cif_loop_matrix[rows[val_index]-1,
                                    cols[val_index]-1] = values[val_index]

                vars(self)['%s_matrix' % class_attribute[:-4]
                           ] = np.matrix(cif_loop_matrix)

    # takes mode matrix
    def constraint_array(self, mode_number):
        mode_matrix = self.iso_displacivemodematrix_matrix
        if (len(mode_matrix[:, 0]) >= mode_number) and (mode_number > 0):
            vector = np.zeros([len(mode_matrix[:, 0]), 1])
            vector[mode_number-1, 0] = 1
            return mode_matrix*vector

        else:
            print("mode number is outside range")
            return 0

        # split constraints up into x, y and z array.



    def make_x_y_z_constrains(self, c_array):
        # make x, y, z array
        xyz_label = self.iso_deltacoordinate_label
        c_x = []  # np.zeros[len(xyz_label)/3,1]
        c_y = []  # np.zeros[len(xyz_label)/3,1]
        c_z = []  # np.zeros[len(xyz_label)/3,1]
        for i in range(len(c_array)/3):
            c_x.append(c_array[i*3, 0])
            c_y.append(c_array[(i*3)+1, 0])
            c_z.append(c_array[(i*3)+2, 0])

        return c_x, c_y, c_z

# method to return constraints all in one step.
    def constraints_out(self, mode_number):
        self.make_matrices()
        c_array = self.constraint_array(mode_number)
        c_x, c_y, c_z = self.make_x_y_z_constrains(c_array)
        # print c_x, c_y, c_z
        return c_x, c_y, c_z
    
    def irrep_list(self):
      irreps = []
      for i in self.iso_displacivemode_label:
        m = re.search("\](.+?)\(",i)
        if m:
          label = m.group(1)
          irreps.append(label)
      irreps = list(dict.fromkeys(irreps))
      return irreps       
        
        
        
    def write_inp(self,filename="batch_modes.inp"):
      f = open(filename,"w")
      irreps = self.irrep_list()
      f.write("r_wp 0.0 r_exp 0.0 r_p 0.0 r_wp_dash 0.0 r_exp_dash 0.0 weighted_Durbin_Watson 0.0 gof 0.0 \n")
      f.write("iters 10000000 \n")
      f.write("chi2_convergence_criteria 0.001 \n")
      f.write("continue_after_convergence")
      f.write(" \' Set fixed number of cycles, replace xxx with desired number \n")
      f.write("prm dummy  0.00000` val_on_continue = If(Cycle == xxx, Get(iters) = 0, 0); \n")
      f.write("\n")
      f.write(" \' File input macro example. The ##n## is replaced when called. \n")
      f.write("macro file_in_(n) \n")
      f.write("{ \n C:\Documents\Files\##n##K_pdf.xye \n } \n")
      f.write(" \' VAR is defined when running on command line \n")
      f.write("xdd file_in_(VAR) \n")
      f.write("   weighting = 1; \n")
      f.write(" \' Comment out the following line if using X-ray PDF \n")
      f.write("   neutron_data \n   pdf_data \n")
      f.write(" \' Following lines tell Topas to rebin the data into bins of size xxx, comment out if unnecessary \n")
      f.write("   rebin_start_x_at 0 \n")
      f.write("   rebin_with_dx_of xxx \n")
      f.write(" \' Define fitting range \n")
      f.write("   start_X xxx \n")
      f.write("   finish_X yyy \n")
      f.write(" \' The following are instrumental parameters and should be determined using a standard \n")
      f.write(" \' Only use one of dQ_lor_damping and dQ_damping \n")
      f.write(" \' convolute_alpha slows the refinement down considerably. Doesn't affect qualitative results but does improve the R_wp. \n")
      f.write(" \' dQ_lor_damping(!dQ,dq_val,!lor,lor_val) \n")
      f.write("   dQ_damping(!dQ,dq_val) \n")
      f.write("   convolute_Qmax_Sinc(!Qmax,qmax_val) \n")
      f.write(" \' convolute_alpha(!alpha,alpha_val) \n")
      f.write(" \n \n \n \n")
      f.write("\'=== PHASE INFORMATION === \n")
      f.write("\' If you want to constrain lattice parameters to be equal, edit the lprm_a to lprm for the values you wish to constrain \n")
      f.write("      str \n")
      f.write("         a lprm_a %s \n"%self.cell_length_a)
      f.write("         b lprm_b %s \n"%self.cell_length_b)
      f.write("         c lprm_c %s \n"%self.cell_length_c)
      f.write("         al %s \n"%self.cell_angle_alpha)
      f.write("         be %s \n"%self.cell_angle_beta)
      f.write("         ga %s \n"%self.cell_angle_gamma)
      f.write("         volume %s \n"%self.cell_volume)
      f.write("         space_group %s \n"%self.symmetry_Int_Tables_number)
      f.write("\' Topas and most PDF generating programs define bins differently. Replace xxx with half of the bin size. It needs to be negative \n")
      f.write("         pdf_zero -xxx \n")
      f.write("\n \n \n \n")
      f.write("\'====== MODE DEFINITIONS ===== \n")
      iso_label = self.iso_displacivemode_label
      for irrep in irreps:
        f.write("#ifndef %s \n"%irrep)
        for i in range(len(iso_label)):
          m = re.search("\](.+?)\(",iso_label[i])
          if m:
            label = m.group(1)
            if irrep == label:
              f.write("prm !a%s 0.00000 min -1 max 1 val_on_continue = Rand(-0.05,0.05); \' %s \n"%(i+1,iso_label[i]))
        f.write("#endif \n")        
        f.write("#ifdef %s \n"%irrep)
        for i in range(len(iso_label)):
          m = re.search("\](.+?)\(",iso_label[i])
          if m:
            label = m.group(1)
            if irrep == label:       
              f.write("prm a%s 0.00000 min -1 max 1 val_on_continue = Rand(-0.05,0.05); \' %s \n"%(i+1,iso_label[i]))
        f.write("#endif \n")      
      f.write("\n \n \n \n")
      f.write("\'===== MODE TO DELTAS TRANSFORMATION ===== \n")
      dm_rows = self.iso_displacivemodematrix_row
      dm_cols = self.iso_displacivemodematrix_col
      dm_vals = self.iso_displacivemodematrix_value
      dc_labels = self.iso_deltacoordinate_label
      
      for i in range(len(dc_labels)):
        dc_cols = []
        dc_vals = []
        f.write("prm %s ="%dc_labels[i])
        for j in range(len(dm_rows)):
          
          if int(dm_rows[j]) == i + 1:
            dc_cols.append(dm_cols[j])
            dc_vals.append(dm_vals[j])
        
        
        for j in range(len(dc_cols)):
          if j == len(dc_cols) - 1:
            if float(dc_vals[j]) > 0:
              f.write(" + %s*a%s ;: 0.0 \n"%(dc_vals[j],dc_cols[j]))
            else:
              f.write(" %s*a%s ;: 0.0\n"%(dc_vals[j],dc_cols[j]))
          else:
            if float(dc_vals[j]) > 0:
              f.write(" + %s*a%s "%(dc_vals[j],dc_cols[j]))
            else:
              f.write(" %s*a%s "%(dc_vals[j],dc_cols[j]))
      icl = self.iso_coordinate_label
      icf = self.iso_coordinate_formula
      f.write("\n")
      xs = self.atom_site_fract_x
      ys = self.atom_site_fract_y
      zs = self.atom_site_fract_z
      formula = []
      for i in range(len(xs)):
        formula.append(xs[i])
        formula.append(ys[i])
        formula.append(zs[i])
      for i in range(len(icl)):
        
        #formula = icf[i].strip("\"")
        f.write("prm %s = %s + %s;\n"%(icl[i],formula[i],dc_labels[i])) 
      atoms = self.atom_site_type_symbol
      f.write("\'Parameters for beq functions and scale. It is recommended to refine scale using average structure and fix within a range of that value \n")
      elements = set(atoms)
      elements = list(elements)
      elements = [x.strip("+-123456789") for x in elements]
      beq_vars = []
      for i in elements:
        bv = i + "_beq"
        beq_vars.append(bv)
      for i in beq_vars:
        f.write("prm %s 0.01 min 0 max 1 val_on_continue = Rand(0,0.1); \n"%i)
      f.write("prm rv 0.01 min 0 val_on_continue = Rand(0,0.1); \n")
      f.write("prm r2v 0.01 min 0 val_on_continue = Rand(0,0.1); \n")
      f.write("prm phase_scale 1 \n")
      sites = self.atom_site_label
      occs = self.atom_site_symmetry_multiplicity
      for i in range(len(sites)):
        f.write("site %s "%sites[i])
        f.write(" x = %s; "%icl[3*i])
        f.write(" y = %s; "%icl[3*i + 1])
        f.write(" z = %s; "%icl[3*i + 2])
        f.write(" occ %s %s "%(atoms[i].strip("+-123456789"),occs[i]))      
        for j in range(len(elements)):
          if atoms[i].strip("+-123456789") == elements[j]:
            bn = "beq_" + str(i) 
            f.write("beq_r_r2(%s,=%s;,d1,=rv;,d2,=r2v;) \n"%(bn,beq_vars[j]))
      f.write("scale = phase_scale; \n")
      f.write("\' Example macro to create output files. IRREP and VAR are replaced with a command line macro. \n")
      f.write("macro file_out_(m,n) \n { \n m##_xxxxxx_##n##_yyyyyy_out.txt \n } \n")
      f.write("out_prm_vals_on_convergence file_out_(IRREP,VAR) \n")
      
      

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import re
import matplotlib.markers
import matplotlib.patches as mpatches
import os
from math import ceil

import scipy.constants
from scipy.constants import hbar, k
import selenium
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time






class sapa_utils:
    """
    This class contains utilities for use with output from the sapa class for performing
    symmetry-adapted PDF analysis with Topas.

    Arguments
    ---------

    cif_file :  The sapa class takes an ISODISTORT-generated cif as it's only input. The filename is required to be a
        string and must have a line designating the end of the file (e.g. #End). Format - string

    temps : A list of temperatures or other variable to iterate over. The elements in this list replace the ##n##
            in the filenameformat argument. Format - list of strings

    sample : A sample name used for reference. It is used to read in the correct output file names. Format - string

    irrep : An irreducible representation label. NOTE: If irrep is required as a function argument, it should be entered
            in ISODISTORT format (i.e. endings with + or - preserved, not as described later for accessing class
            variables. Format - string

    sigma - a value that represents a noticeable different in fit quality for calculating bwma and bwmma. Format - string

    varname - The name of a variable refined during the SAPA refinements. Format - string

    Methods
    -------

    sapa_utils(cif_file) : This imports the ISODISTORT CIF into class variables, which is needed for other class functions

    import_files(sample, temps) : This processes the output files from SAPA

    plot_var(varname) :  This saves a plot of the value of the refined variable for each irrep and each temperature for
                            the lowest Rwp refinement.

    plot_all() : This executes plot_var(varname) for every varname in the output file (excluding Cycle and Iteration)

    calc_bwmas() : This calculates the bwma and bwmma for each irrep and temperatures and stores them as class variables.

    make_histograms() : This saves histograms of Rwp for each irrep and temperature

    generate_cif(irrep,temp) : For a given irrep (in isodistort format, e.g. X5+ not X5p) and temperature, this takes
                            the amplitudes for the lowest Rwp refinement and outputs a CIF of the structure.

    Variables
    ---------

    Accessing irrep variables - due to python, the irrep class variables cannot be stored if they end in a + or -
    When importing these variables, + -> p and - -> m, i.e. to access the information stored for the X5+ irrep (for a
    class object is named sapa) you use sapa.X5p - these variables are case sensitive.

    The irrep variables themselves are dictionaries. After import_files() has been run, the output file for each temp
    can be accessed as a Pandas dataframe by using sapa.irrep["temp"] and the values for a variable corresponding to
    the lowest Rwp for each temperature can be accessed with sapa.irrep["varname"]. Once calc_bwmas() has been ran,
    the bwmas and bwmmas can similarly be accessed using sapa.irrep["bwma"] and sapa.irrep["bwmma"], respectively.
    """
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
        #vars(self)["temps"] = temps
        f.close()
    def irrep_list(self):
      irreps = []
      for i in self.iso_displacivemode_label:
        m = re.search("\](.+?)\(",i)
        if m:
          label = m.group(1)
          irreps.append(label)
      irreps = list(dict.fromkeys(irreps))

      return irreps




    def import_files(self, sample, temps):
        self.irreps = self.irrep_list()
        self.sample = sample
        temps = [str(x) for x in temps]
        self.temps = temps

        irreps = self.irrep_list()
        irreps.append("nomodes")
        irrep_vars = []
        for irrep in irreps:
            irrep1 = irrep
            if irrep.endswith("+") == True:
                irrep1 = irrep.replace("+", "p")
            elif irrep.endswith("-") == True:
                irrep1 = irrep.replace("-", "m")
            irrep_vars.append(irrep1)
            vars(self)[irrep1] = {}
            for temp in self.temps:
                fn = irrep + "_" + sample + "_" + temp + "_out.txt"
                df = pd.read_csv(fn, sep="\s+", index_col=None)
                df = df.drop([0])
                df = df.sort_values("Rwp")
                df.index = range(len(df.index))
                nmfn = "nomodes_" + sample + "_" + temp + "_out.txt"
                df1 = pd.read_csv(nmfn, sep="\s+", index_col=None)
                df1 = df1.drop([0])
                df1 = df1.sort_values("Rwp")
                df1.index = range(len(df1.index))
                nm_rwp = df1.at[0,"Rwp"]
                df["delrwp"] = df["Rwp"] - nm_rwp
                vars(self)[irrep1]["%s"%temp] = df


        for irrep in irrep_vars:
            varnames = vars(self)[irrep][self.temps[0]].columns.values.tolist()
            for name in varnames:
                best = []
                for temp in self.temps:
                    best.append(vars(self)[irrep][temp].at[0,"%s"%name])
                vars(self)[irrep]["%s"%name] = best
        irrep_vars.remove("nomodes")
        self.irrep_vars = irrep_vars
        for irrep in irrep_vars:
            rwplist = vars(self)[irrep]["Rwp"]
            rwplist = [float(x) for x in rwplist]
            nomodesrwp = vars(self)["nomodes"]["Rwp"]
            nomodesrwp = [float(x) for x in nomodesrwp]
            deltas = []
            for i in range(len(rwplist)):
                drwp = rwplist[i] - nomodesrwp[i]
                deltas.append(drwp)
            vars(self)[irrep]["drwp"] = deltas
        best_deltas = []
        for i in range(len(self.temps)):
            deltalist = []
            for irrep in irrep_vars:
                delta = vars(self)[irrep]["drwp"][i]
                deltalist.append(delta)
            deltalist = [float(x) for x in deltalist]
            best_delta = min(deltalist)
            best_deltas.append(best_delta)
        vars(self)["nomodes"]["drwpmin"] = best_deltas
        #self.irreps = self.irrep_list().remove("nomodes")




    #def irrep_vars(self):
    #        #irreps = self.irrep_list()
    #        irrep_vars = []
    #        for irrep in self.irreps:
    #            irrep1 = irrep
    #            if irrep.endswith("+") == True:
    #                irrep1 = irrep.replace("+","p")
    #            elif irrep.endswith("-") == True:
    #                irrep1 = irrep.replace("-","m")
    #            irrep_vars.append(irrep1)
    #        return irrep_vars

    #irrep_vars = irrep_vars()

    def plot_var(self,varname):

        plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.gist_rainbow(np.linspace(0,1,len(self.irrep_vars))))
        markers = ["^","v","<",">","*","o","s"]
        markerrpt = ceil(len(self.irrep_vars)/len(markers))
        markerrpt = markerrpt + 1
        markerlist = np.tile(markers,int(markerrpt))

        fig = plt.figure()
        for i in range(len(self.irrep_vars)):
            irrep = self.irrep_vars[i]
            plt.plot(self.temps,vars(self)[irrep]["%s"%varname],linestyle="-",marker=markerlist[i],label=self.irreps[i],markersize=5,linewidth=2.0)
        plt.xlabel("Temperature (K)",fontsize=9)
        plt.ylabel("%s"%varname,fontsize=9)
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
        plt.tight_layout()
        plt.legend(loc="best",ncol=6)
        plt.savefig("%s_%s_plot.png"%(self.sample,varname),dpi=300)
        plt.close(fig)

    def plot_all(self,bwmas=False):

        ncols = len(vars(self)[self.irrep_vars[0]][self.temps[0]].columns)
        varnames = vars(self)[self.irrep_vars[0]][self.temps[0]].columns.values.tolist()
        varnames.remove("Cycle")
        varnames.remove("Iter")

        if bwmas == True:
            varnames.append("bwma")
            varnames.append("bwmma")
            varnames.remove("mode_amps")
            varnames.remove("weights")
            varnames.remove("weights*amps")
        print(varnames)
        regex = re.compile(r"^a[0-9]*")
        varnames = [x for x in varnames if not regex.search(x)]
        print(varnames)
        for i in varnames:
            self.plot_var(i)

    def get_indices(self,irrep):
      indices = []
      iso_label = self.iso_displacivemode_label
      for i in range(len(iso_label)):
        m = re.search("\](.+?)\(",iso_label[i])
        if m:
            label = m.group(1)
            if irrep == label:
                indices.append(i)
      return indices

    def get_rwglob(self):
      rwplist = []
      #irreps = self.irrep_vars()
      for temp in self.temps:
        for irrep in self.irrep_vars:
          rwplist.append(vars(self)[irrep][temp].at[0,"Rwp"])
      rwplist = [float(x) for x in rwplist]
      rwglob = min(rwplist)
      return rwglob

    def get_var_names(self,irrep):
      inds = self.get_indices(irrep)
      inds = [int(x) for x in inds]
      inds = [x + 1 for x in inds]
      varlist = ["a" + str(x) for x in inds]
      return varlist

    def get_norms(self,irrep):
      inds = self.get_indices(irrep)
      norms = self.iso_displacivemodenorm_value
      normlist = []
      for i in inds:
        normlist.append(norms[i])
      return normlist

    def get_dbwma(self,irrep,sigma):
        if irrep.endswith("+") == True:
            irrep1 = irrep.replace("+","p")
        elif irrep.endswith("-") == True:
            irrep1 = irrep.replace("-","m")
        else:
            irrep1 = irrep
        varlist = self.get_var_names(irrep)
        norms = self.get_norms(irrep)
        bwmas = []
        bwmmas = []
        for i in range(len(self.temps)):
            vars(self)[irrep1][self.temps[i]]["mode_amps"] = (vars(self)[irrep1][self.temps[i]][varlist[0]] / float(norms[0])) ** 2
            for j in range(1, len(norms)):
                vars(self)[irrep1][self.temps[i]]["mode_amps"] += (vars(self)[irrep1][self.temps[i]][varlist[j]] / float(norms[j])) ** 2
            vars(self)[irrep1][self.temps[i]]["weights"] = np.longdouble(np.exp((vars(self)["nomodes"]["drwpmin"][i] - vars(self)[irrep1][self.temps[i]]["delrwp"]) / sigma))
            sum_weights = np.sum(vars(self)[irrep1][self.temps[i]]["weights"])
            vars(self)[irrep1][self.temps[i]]["mode_amps"] = np.sqrt(vars(self)[irrep1][self.temps[i]]["mode_amps"])
            vars(self)[irrep1][self.temps[i]]["weights*amps"] = vars(self)[irrep1][self.temps[i]]["mode_amps"] * vars(self)[irrep1][self.temps[i]]["weights"]
            bwmma = np.sum(vars(self)[irrep1][self.temps[i]]["weights*amps"]) / sum_weights
            bwma = np.sum(vars(self)[irrep1][self.temps[i]]["weights*amps"])
            bwmas.append(bwma)
            bwmmas.append(bwmma)
        return bwmas, bwmmas

    def get_bwma(self,irrep,sigma):
        if irrep.endswith("+") == True:
            irrep1 = irrep.replace("+","p")
        elif irrep.endswith("-") == True:
            irrep1 = irrep.replace("-","m")
        else:
            irrep1 = irrep
        varlist = self.get_var_names(irrep)
        norms = self.get_norms(irrep)
        Rwglob = self.get_rwglob()
        bwmas = []
        bwmmas = []
        for temp in self.temps:
          vars(self)[irrep1][temp]["mode_amps"] = (vars(self)[irrep1][temp][varlist[0]]/float(norms[0]))**2
          for i in range(1,len(norms)):
            vars(self)[irrep1][temp]["mode_amps"] += (vars(self)[irrep1][temp][varlist[i]]/float(norms[i]))**2
          vars(self)[irrep1][temp]["weights"] = np.longdouble(np.exp((Rwglob - vars(self)[irrep1][temp]["Rwp"])/sigma))
          sum_weights = np.sum(vars(self)[irrep1][temp]["weights"])
          vars(self)[irrep1][temp]["mode_amps"] = np.sqrt(vars(self)[irrep1][temp]["mode_amps"])
          vars(self)[irrep1][temp]["weights*amps"] = vars(self)[irrep1][temp]["mode_amps"]*vars(self)[irrep1][temp]["weights"]
          bwmma = np.sum(vars(self)[irrep1][temp]["weights*amps"])/sum_weights
          bwma = np.sum(vars(self)[irrep1][temp]["weights*amps"])
          bwmas.append(bwma)
          bwmmas.append(bwmma)
        return bwmas, bwmmas

    def calc_bwmas(self,sigma):

        for irrep in self.irreps:
            bwma, bwmma = self.get_bwma(irrep,sigma)
            vars(self)[irrep]["bwma"] = bwma
            vars(self)[irrep]["bwmma"] = bwmma

    def calc_dbwmas(self,sigma):
        # amplitudes calculated here are RMS, i.e 1/sqrt(2)*amplitude of oscillation
        for irrep in self.irreps:
            if irrep.endswith("+") == True:
                irrep1 = irrep.replace("+", "p")
            elif irrep.endswith("-") == True:
                irrep1 = irrep.replace("-", "m")
            else:
                irrep1 = irrep
            bwma, bwmma = self.get_dbwma(irrep,sigma)
            vars(self)[irrep1]["bwma"] = bwma
            vars(self)[irrep1]["bwmma"] = bwmma

    def make_histograms(self):
        if not os.path.exists("./histograms"):
            os.makedirs("./histograms")
        for irrep in self.irrep_vars:
            for temp in self.temps:
                rwp = vars(self)[irrep][temp]["Rwp"].tolist()
                rwp = [float(x) for x in rwp]
                start = int(np.floor(min(rwp)))
                end = int(np.ceil(max(rwp)))
                binlist = list(range(start,end+1))
                f = plt.figure()
                plt.hist(rwp,bins=binlist)
                plt.xlabel("Rwp")
                plt.ylabel("Number")
                plt.grid(True)
                plt.savefig("./histograms/{0}_{1}_{2}_rwp_histogram.png".format(self.sample,irrep,temp),dpi=300)
                plt.close(f)
                
    def generate_cif(self,irrep,temp):
        if irrep.endswith("+") == True:
            irrep1 = irrep.replace("+","p")
        elif irrep.endswith("-") == True:
            irrep1 = irrep.replace("-","m")
        else:
            irrep1 = irrep
        df = vars(self)[irrep1][temp]
        varnames = self.get_var_names(irrep)
        amps = []
        for i in varnames:
            amps.append(df.at[0,"%s"%i])
        amps = [float(x) for x in amps]
        #amps = [0.05,0,0]
        print(amps)
        indices = self.get_indices(irrep)

        indices = [int(x) + 1 for x in indices]
        print(indices)
        file = open("%s_%s_%s.cif"%(self.sample,irrep,temp),"w")
        file.write("data_sapa-output \n \n")
        cla = df.at[0,"lprm_a"]
        file.write("_cell_length_a %s \n"%cla)
        if "lprm_b" and "lprm_c" in df.columns.values.tolist():
            clb = df.at[0,"lprm_b"]
            clc = df.at[0,"lprm_c"]
            file.write("_cell_length_b %s \n"%clb)
            file.write("_cell_length_c %s \n"%clc)
        elif "lprm_b" not in df.columns.values.tolist() and "lprm_c" in df.columns.values.tolist():
            clc = df.at[0,"lprm_c"]
            file.write("_cell_length_b %s \n"%cla)
            file.write("_cell_length_c %s \n"%clc)
        elif "lprm_b" and "lprm_c" not in df.columns.values.tolist():
            file.write("_cell_length_b %s \n"%cla)
            file.write("_cell_length_c %s \n"%cla)
        file.write("_cell_angle_alpha %s \n"%self.cell_angle_alpha)
        file.write("_cell_angle_beta %s \n"%self.cell_angle_beta)
        file.write("_cell_angle_gamma %s \n"%self.cell_angle_gamma)
        file.write("_cell_volume %s \n"%self.cell_volume)
        file.write("_symmetry_space_group_name_H-M \"P 1\" \n")
        file.write("_symmetry_Int_Tables_number 1 \n")
        file.write("_space_group.reference_setting \'001:P 1\' \n")
        file.write("_space_group.transform_Pp_abc a,b,c;0,0,0 \n \n")
        file.write("loop_ \n")
        file.write("_space_group_symop_id \n")
        file.write("_space_group_symop_operation_xyz \n")
        file.write("1 x,y,z \n \n")
        file.write("loop_ \n")
        file.write("_atom_site_label \n")
        file.write("_atom_site_type_symbol \n")
        file.write("_atom_site_symmetry_multiplicity \n")
        file.write("_atom_site_Wyckoff_label  \n")
        file.write("_atom_site_fract_x \n")
        file.write("_atom_site_fract_y  \n")
        file.write("_atom_site_fract_z \n")
        file.write("_atom_site_occupancy  \n")
        for i in range(len(self.atom_site_label)):
            print(self.atom_site_label[i])
            file.write("%s %s %s %s"%(self.atom_site_label[i],self.atom_site_type_symbol[i],self.atom_site_symmetry_multiplicity[i],self.atom_site_Wyckoff_label[i]))
            dx = float(0)
            dy = float(0)
            dz = float(0)
            for k in range(len(self.iso_displacivemodematrix_col)):
                val = self.iso_displacivemodematrix_row[k]
                val = int(val)
                index = int(self.iso_displacivemodematrix_col[k])
                delta = self.iso_displacivemodematrix_value[k]
                delta = float(delta)
                
                if (val == (3*i)+1) and (index in indices):
                    listindex = indices.index(index)
                    #print(index,val,i)
                    dx = dx + delta*amps[listindex]
                    #print(delta,amps[listindex],dx)
                elif (val == (3*i)+2) and (index in indices):
                    listindex = indices.index(index)
                    #print(index,val,i)
                    dy += delta*amps[listindex]
                elif (val == (3*i)+3) and (index in indices):
                    listindex = indices.index(index)
                    #print(index,val,i)
                    dz += delta*amps[listindex]
            fracx = np.longdouble(self.atom_site_fract_x[i]) + np.longdouble(dx)
            fracy = np.longdouble(self.atom_site_fract_y[i]) + np.longdouble(dy)
            fracz = np.longdouble(self.atom_site_fract_z[i]) + np.longdouble(dz)
            if fracx < 0:
                fracx = fracx + 1
            if fracy < 0:
                print(fracy)
                fracy = fracy + 1
                print(fracy)
            if fracz < 0:
                fracz = fracz + 1
            if fracx > 1:
                fracx = fracx - 1
            if fracy > 1:
                fracy = fracy - 1

            if fracz > 1:
                fracz = fracz - 1
            fracx = str(fracx)
            fracy = str(fracy)
            fracz = str(fracz)
            file.write(" %s %s %s "%(fracx,fracy,fracz))
            file.write("%s \n"%self.atom_site_occupancy[i])
        #file.write("\n \n")
        #file.write("#End")
        file.close()










            
            

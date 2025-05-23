# ==============================================================================
# Python script that reads config.ini with parameter settings
# for the sectional aerosol model in NorESM
# Output:
# namelist oslo_sectional_nl with settings for NorESM
# my_chem_mech.in: edited to contain aerosol tracers
# ==============================================================================

# TODO: Christina Brodowsky: Add documentation
# TODO: double underscores to make all members private
# TODO: Use only spaces, NO TABS!!

import numpy as np
import configparser
import sys, os
import argparse

#_CIMEROOT = os.environ.get("CIMEROOT") # TODO import in main fct instead
#if _CIMEROOT is None:
#    raise SystemExit("ERROR: must set CIMEROOT environment variable")

#_LIBDIR = os.path.join(_CIMEROOT, "scripts", "Tools")
#sys.path.append(_LIBDIR)

#from standard_script_setup          import * # TODO find out if we need these..
#from CIME.XML.standard_module_setup import *
from CIME.case                      import Case


def _parse_range(config, section, option):
    ''' Parse a range from the config file'''
    range_str = config.get(section, option).split(',')
    parsed_range = [float(i) for i in range_str]
    return parsed_range

# ==============================================================================
# Types for bin, range and species information
# ==============================================================================
class _AerosolSpecies:
    def __init__(self, config, species):
        self.active = config.getboolean(species, 'active', fallback=False)
        self.short_name = config.get(species, 'short_name')
        self.long_name = config.get(species, 'long_name')
        self.composition = config.get(species, 'composition')
        self.soluble = config.getboolean(species, 'soluble')
        self.range_bnds = _parse_range(config, species, 'range_bounds')
        self.range_idx = []
    def check_range_bnds(self, range_bnds):
        if any(i not in range_bnds for i in self.range_bnds):
            sys.exit('ERROR: Species range bound not equal to range bounds')
    def get_range_idx(self, range_bnds):
        # adjust species bounds to range bounds
        idx0 = (np.abs(range_bnds - self.range_bnds[0])).argmin()
        idx1 = (np.abs(range_bnds - self.range_bnds[1])).argmin()
        self.range_bnds = [range_bnds[idx0], range_bnds[idx1]]
        # get range indices for each species
        self.range_idx.extend([i for i in range(1, len(range_bnds))
                                if range_bnds[i] > self.range_bnds[0]
                                and range_bnds[i] <= self.range_bnds[1]])

# bin settings, maybe add 'method' to allow for other than
# volume ratio. e.g. 'custom' -> user defined bin bounds
class _BinSpecs:
    def __init__(self, config):
        self.N = config.getint('BIN SPECS', 'nbin')
        self.r_1 = config.getfloat('BIN SPECS', 'radius_1')
        self.r_N = config.getfloat('BIN SPECS', 'radius_N')
        self.r = None
        self.r_bnds = None
    def calc_bins(self):
        # Bin calculations - Volume ratio approach
        V_rat = (self.r_N / self.r_1) ** (3 / (self.N-1)) 		# calculate volume ratio
        v_0 = 3/4 * np.pi * (self.r_1) ** 3 # # calculate smallest volume
        v = v_0 * V_rat ** np.arange(self.N) # calculate all volumes for bins
        v_lo = (2*v) / (1+V_rat) # lower radius bounds
        v_hi = V_rat*v_lo[-1] # upper bnd for largest bin
        v_bnds = np.append(v_lo, v_hi) # append to get volume bounds

        self.r = (v*4/(3*np.pi))**(1/3) # calculate radii
        self.r_bnds = (v_bnds*4/(3*np.pi))**(1/3) # calculate radius bounds from volume

class _RangeSpecs:
    def __init__(self, config):
        self.ranges = config.getboolean('RANGE SPECS', 'ranges')
        self.range_bnds = np.asarray(_parse_range(config, 'RANGE SPECS', 'range_bounds'))
        self.range_bnd_bin_idx = np.array([1])
    def adjust_range_bnds(self, r_bnds):
        if self.ranges: # TODO: add explanation on what this computation is doing
            self.range_bnds[0] = r_bnds[0]
            self.range_bnds[-1] = r_bnds[-1]
            for i in range(1, len(self.range_bnds)-1):
                idx = (np.abs(r_bnds - self.range_bnds[i])).argmin()
                self.range_bnds[i] = r_bnds[idx]
                self.range_bnd_bin_idx = np.append(self.range_bnd_bin_idx, idx+1) # account for 1-indexing in Fortran
   			self.range_bnd_bin_idx = np.append(self.range_bnd_bin_idx, len(r_bnds)) # account for 1-indexing in Fortran

        else:
            # If no ranges (classic sectional scheme), set range bounds equal to bin bounds
            self.range_bnds = r_bnds

def bin_config(aerconf_file, chem_mech_file):
    CAM_CONFIG_OPTS     = case.get_value("CAM_CONFIG_OPTS")

    # ==============================================================================
    # Read input from config file
    # ==============================================================================
    # INPUT
    config = configparser.ConfigParser()
    config.read(aerconf_file)

    # ==============================================================================
    # Initialize
    # ==============================================================================´

    # initialize bins
    bin_specs = _BinSpecs(config)
    bin_specs.calc_bins()

    # initialize ranges
    range_specs = _RangeSpecs(config)

    # initialize aerosol species
    species_obj_list = []
    for section in config.sections():
    if section not in ['BIN SPECS', 'RANGE SPECS']:
        species_obj_list.append(_AerosolSpecies(config, section))

    [species_obj.check_range_bnds for species_obj in species_obj_list]  # check species range bnds

    # adjust range range_bnds
    range_specs.adjust_range_bnds(bin_specs.r_bnds)

    # get species range indices
    [species_obj.get_range_idx(range_specs.range_bnds) for species_obj in species_obj_list]

    nspecies = len(species_obj_list)
	# =====================================================================
	# Write to namelist
	# Change to write to e.g. atm_in namelist?
	# Write tracers instead to chem_mech file or my_chem_mech.in in case dir?
	# =====================================================================

	# TODO: Find out where to write out

    f = open(os.path.join(caseroot, "oslo_sectional_in"), "w")
    f.write("&oslo_sectional_properties_nl\n")
    f.write(" oslo_sectional_nspecies		=  ")
    f.write(f"{nspecies} \n")

    f.write(" oslo_sectional_bin_bounds		=  ")
    for i in range(bin_specs.N):
        f.write(f"'{bin_specs.r_bnds[i]:.3f}D0:{bin_specs.r_bnds[i+1]:.3f}D0'")
        if i != bin_specs.N-1:
            f.write(', ')
    f.write("\n")

    f.write(" oslo_sectional_bin_centers		=  ")
    for i in range(bin_specs.N):
        f.write(f"'{bin_specs.r[i]:.3f}D0'")
        if i != bin_specs.N-1:
            f.write(', ')
    f.write("\n")

    f.write(" oslo_sectional_range_bounds		=  ")
    for i in range(len(range_specs.range_bnds)-1):
        f.write(f"'{range_specs.range_bnd_bin_idx[i]}:{range_specs.range_bnd_bin_idx[i+1]}'")
        if i != len(range_specs.range_bnds)-2:
            f.write(', ')
    f.write("\n")

    for species in species_obj_list:
        if species.active:
            f.write("&oslo_sectional_properties_aerosol_nl\n")
            f.write(f" oslo_sectional_aerosol_name		=  '{species.short_name}' \n")
            f.write(" oslo_sectional_aerosol_range		=  ")
            f.write(f"'{species.range_idx[0]}:{species.range_idx[-1]}' \n")
            f.write(f" oslo_sectional_aerosol_soluble		=  .{species.soluble}. \n")

    f.write("/\n")
    f.close()

    # ==============================================================================
    # Prepare output for chem_mech.in file
    # The filenames and file dirs will need to be changed
    # e.g. the directory where the chem_mech file comes from. Probably can be retrieved
    # somehow from the compset.
    # Also the output, currently 'my_chem_mech.in' needs to go to the casedir
    # The code below reads the chem_mech.in file line by line and looks for keywords
    # "Solution" and "Implicit". Below these, the composition of the tracers and the
    # names of the tracers are added. The file is then written out to my_chem_mech.in
    # ==============================================================================

    # TODO: write output to a sensible place -> casefolder?

    config_opts = CAM_CONFIG_OPTS.split(' ')
    chem_index = config_opts.index('-chem')
    CHEM_OPT = config_opts[chem_index + 1]
    chemconf = os.path.join(srcroot, "src", "chemistry", "pp_" + CHEM_OPT, "chem_mech.in")
    chem_outfile = os.path.join(caseroot, 'my_chem_mech.in') # caseroot/...
    config_opts += ['-usr_mech_infile', chem_outfile]
    # TODO: Long-term find a better way to reference this chem_mech file than usr_mech_infile
    case.set_value("CAM_CONFIG_OPTS", config_opts)

    composition_list = []
    implicit_list = []

    for species in species_obj_list:
        if species.active:
            for i in range(len(species.range_idx)):
                composition_list.append(
                    f"{species.short_name}_R{species.range_idx[i]} -> {species.composition}"
                    ) # test
                implicit_list.append(
                    f"{species.short_name}_R{species.range_idx[i]}"
                    )

    for i in range(1,bin_specs.N+1):
        composition_list.append(
            f"num_{i} -> H"
            )
        implicit_list.append(
            f"num_{i}"
    )

    with open(chemconf, 'r') as chem_file:
    lines = chem_file.readlines()

    modified_chem = []

    for line in lines:
    # write lines from the old file for the new file
        modified_chem.append(line)
    # add species composition
        if 'Solution' in line and not 'End' in line and not 'Classes' in line:
            for i in range(len(composition_list)):
                modified_chem.append(f"{composition_list[i]}\n")
   	# add species for advection
        if 'Implicit' in line and not 'End' in line:
            for i in range(len(implicit_list)):
                modified_chem.append(f"{implicit_list[i]}\n")
    # write out to my_chem_mech.in
    with open(chem_outfile, 'w') as file:
        file.writelines(modified_chem)

def _main_func(): # TODO: import cimeroot and caseroot
    parser = argparse.ArgumentParser(description="Process aerosol configuration for the" \
    "sectional aerosol scheme in NorESM, write the namelist and add the tracers to chemistry.")
    parser.add_argument('--aerconf', required=True, help='Path to the aerosol configuration file')
    parser.add_argument('--chem_mech', required=True, help='Path to the chem_mech.in file')

    args = parser.parse_args()

    if not os.path.isfile(args.aerconf):
        sys.exit('Error: Specified aerosol configuration file does not exist')
    if not os.path.isfile(args.aerconf):
        sys.exit('Error: Specified chem_mech.in file does not exist')

    bin_config(args.aerconf, args.chem_mech)

if __name__ == "__main__":
    _main_func()
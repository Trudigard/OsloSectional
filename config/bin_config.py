# Python script that reads config.ini with parameter settings
# for the sectional aerosol model and writes a namelist file for
# NorESM

import numpy as np
import configparser
import sys

def check_species_bounds(species_bounds):
	''' Check if species bounds are equal to any range bounds'''
	if any(i not in range_bnds for i in species_bounds):
		sys.exit('ERROR: Species range bound not equal to range bounds')

def find_spec_bnds(range_bnds, spec_bnds):
	''' Set species bounds to the closest range bound'''
	range_bnds = np.asarray(range_bnds)
	idx0 = (np.abs(range_bnds - spec_bnds[0])).argmin()
	idx1 = (np.abs(range_bnds - spec_bnds[1])).argmin()
	spec_bnds[0] = range_bnds[idx0]
	spec_bnds[1] = range_bnds[idx1]

def get_spec_array(spec_range, spec_name):
	''' Add species name to array of species names.
	Each Species is named as speciesname_i where i is the index
	of the range it belongs to'''
	spec_list = [
		f"{spec_name}_R{i}"
		for i in range(1, len(range_bnds))
		if range_bnds[i] > spec_range[0] and range_bnds[i] <= spec_range[1]
	]
	return spec_list

def parse_range(config, section, option):
	''' Parse a range from the config file'''
	range_str = config.get(section, option).split(',')
	parsed_range = [float(i) for i in range_str]
	return parsed_range

pi = np.pi

# INPUT
config = configparser.ConfigParser()
config.read('config.ini')

# bin settings
N = config.getint('BIN SPECS', 'nbin')
r_1 = config.getfloat('BIN SPECS', 'radius_1')
r_N = config.getfloat('BIN SPECS', 'radius_N')

# Booleans for ranges, species
ranges = config.getboolean('RANGE SPECS', 'ranges')

# Dictionary for species information
# If species are added, add them here
# Species keys and their corresponding range keys
species_info = {
	'organic': ('organic_range', 'OA'),
	'dust': ('dust_range', 'DU'),
	'sulfate': ('sulfate_range', 'SO4'),
	'blackcarbon': ('blackcarbon_range', 'BC'),
	'blackcarbon_insoluble': ('blackcarbon_insoluble_range', 'BC_is'),
	'seasalt': ('seasalt_range', 'SS'),
	'nitrate': ('nitrate_range', 'NO3')
}

# Check what species are enabled
enabled_species = {species: config.getboolean('SPECIES', species) for species in species_info.keys()}
# get the range bounds for each species
range_bnds = parse_range(config, 'RANGE SPECS', 'range_bounds')

species_ranges = {species: parse_range(config, 'RANGE SPECS', range_key) for species, (range_key, _) in species_info.items()}

for species, bounds in species_ranges.items():
	if any(i not in range_bnds for i in bounds):
		sys.exit('ERROR: Species range bound not equal to range bounds')


# Bin calculations
# calculate volume ratio
V_rat = (r_N / r_1) ** (3 / (N-1))

# calculate smallest volume
v_0 = 3/4 * pi * (r_1) ** 3 # smallest

# calculate all volumes for bins
v = v_0 * V_rat ** np.arange(N)

# calculate radii
r = (v*4/(3*pi))**(1/3)

# volume widths / bin bounds
del_v = (2 * v * (V_rat - 1)) / (1+V_rat)

# lower radius bounds
v_lo = (2*v) / (1+V_rat)

# upper bnd for largest bin
v_hi = V_rat*v_lo[-1]

# append to get volume bounds
v_bnds = np.append(v_lo, v_hi)

# calculate radii from volume
r_bnds = (v_bnds*4/(3*pi))**(1/3)


# find bin bounds closest to range_bnds
if ranges:
	range_bnds[0] = r_bnds[0]
	range_bnds[-1] = r_bnds[-1]
	for i in range(1, len(range_bnds)-1):
		for j in range(len(r_bnds)):
   	     # check if range bnd is between adjacent range bnds
			if range_bnds[i] > r_bnds[j] and range_bnds[i] < r_bnds[j+1]:
			# chose the closest bnd
				if abs(r_bnds[j] - range_bnds[i]) < abs(r_bnds[j+1] - range_bnds[i]):
				# replace the value in range_bnds with the closest value
					range_bnds[i] = r_bnds[j]
				else:
					range_bnds[i] = r_bnds[j+1]
					break
else:
	# If no ranges (classic sectional scheme), set range bounds equal to bin bounds
	range_bnds = r_bnds

# set species bounds equal to range bounds
for species, spec_bnds in species_ranges.items():
	find_spec_bnds(range_bnds, spec_bnds)

spec_list = []
for species, enabled in enabled_species.items():
	if enabled:
		range_key, spec_name = species_info[species]
		spec_list.extend(get_spec_array(species_ranges[species], spec_name))

spec_array = np.array(spec_list)
nspec = len(spec_array)
ntrac = nspec + N

# =====================================================================
# Write to namelist
# Change to write to e.g. atm_in namelist?
# Write tracers instead to chem_mech file or my_chem_mech.in in case dir?
#
# =====================================================================

f = open("bin_nl", "w") # change to "a" later!!!
f.write("&oslo_sectional_nl\n")
f.write("nbin = " + str(N) + "\n")
f.write("nrange = " + str(len(range_bnds)-1) + "\n")
f.write("nspec = " + str(nspec) + "\n")
f.write("ntrac = " + str(ntrac) + "\n")
f.write("bin_list = ")
for i in range(len(r)):
	f.write(str(r[i]) + ",")
f.write("\n")
f.write("bin_bnd_list = ")
for i in range(len(r_bnds)):
	f.write(str(r_bnds[i]) + ",")
f.write("\n")
f.write("range_bnd_list = ")
for i in range(len(range_bnds)):
	f.write(str(range_bnds[i]) + ",")
f.write("\n")
#f.write("OA_active = ." + str(organic) + ". \n")
#f.write("DU_active = ." + str(dust) + ". \n")
#f.write("SO4_active = ." + str(sulfate) + ". \n")
#f.write("BC_active = ." + str(blackcarbon) + ". \n")
#f.write("BC_is_active = ." + str(blackcarbon_insoluble) + ". \n")
#f.write("SS_active = ." + str(seasalt) + ". \n")
#f.write("NO3_active = ." + str(nitrate) + ". \n")
f.write("bin_name_list = ")
for i in range(1, len(r)+1):
	f.write("'aer_" + str(i) + "',")
f.write("\n")
f.write("species_name_list = ")
for i in range(len(spec_array)):
	f.write("'" + str(spec_array[i]) + "',")
f.write("\n")
f.write("/\n")
f.close()
'''
f.write("species_names = 'OA', 'DU', 'SO4', 'BC', 'BC_is', 'SS', 'NO3'\n")

f.write("OA_range = " + str(OA_range[0]) + ", " + str(OA_range[1]) + "\n")
f.write("DU_range = " + str(DU_range[0]) + ", " + str(DU_range[1]) + "\n")
f.write("SO4_range = " + str(SO4_range[0]) + ", " + str(SO4_range[1]) + "\n")
f.write("BC_range = " + str(BC_range[0]) + ", " + str(BC_range[1]) + "\n")
f.write("BC_is_range = " + str(BC_is_range[0]) + ", " + str(BC_is_range[1]) + "\n")
f.write("SS_range = " + str(SS_range[0]) + ", " + str(SS_range[1]) + "\n")
f.write("NO3_range = " + str(NO3_range[0]) + ", " + str(NO3_range[1]) + "\n")'''

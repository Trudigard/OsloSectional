# Python script that reads config.ini with parameter settings
# for the sectional aerosol model and writes a namelist file for
# NorESM

import numpy as np
import configparser
import sys

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
# "range_bounds", "range_idx", and "enabled" are read from the config.ini file
# everything else is set here
species_info = {
	'organic': {
    	'short_name': 'OA',
		'long_name': 'Organic aerosol',
		'range_bounds': (0,0),
		'range_idx': [],
		'composition': 'C',
		'enabled': False,
  		'soluble': True
	},
	'dust': {
    	'short_name': 'DU',
		'long_name': 'Dust aerosol',
		'range_bounds': (0,0),
		'range_idx': [],
		'composition': 'AlSiO5',
		'enabled': False,
  		'soluble': True
	},
	'sulfate': {
		'short_name': 'SO4',
		'long_name': 'Sulfate aerosol',
		'range_bounds': (0,0),
		'range_idx': [],
		'composition': 'NH4SO4',
		'enabled': False,
  		'soluble': True
	},
	'blackcarbon': {
		'short_name': 'BC',
		'long_name': 'Black Carbon aerosol',
		'range_bounds': (0,0),
		'range_idx': [],
		'composition': 'C',
		'enabled': False,
  		'soluble': True
	},
	'blackcarbon_insoluble': {
		'short_name': 'BC_is',
		'long_name': 'Insoluble Black Carbon aerosol',
		'range_bounds': (0,0),
		'range_idx': [],
		'composition': 'C',
		'enabled': False,
  		'soluble': False
	},
	'seasalt': {
		'short_name': 'SS',
		'long_name': 'Seasalt aerosol',
		'range_bounds': (0,0),
      	'range_idx': [],
		'composition': 'NaCl',
		'enabled': False,
  		'soluble': True
	},
	'nitrate': {
		'short_name': 'NO3',
		'long_name': 'Nitrate aerosol',
		'range_bounds': (0,0),
       	'range_idx': [],
		'composition': 'NO3',
		'enabled': False,
  		'soluble': True
	}
}

# get the range bounds
range_bnds = np.asarray(parse_range(config, 'RANGE SPECS', 'range_bounds'))

for species in species_info.keys():
	''' Check if species are enabled.
	This will get the 'enabled' and 'range_bounds' keys for each species
	from the config file. If species are not found, they are set to False'''
	# read from file
	enabled = config.getboolean('SPECIES', species, fallback=False)
	species_range = parse_range(config, 'RANGE SPECS', str(species + '_range'))

	# set values in dictionary
	species_info[species]['enabled'] = enabled
	species_info[species]['range_bounds'] = species_range

	# check if species bounds are == range bounds
	if any(i not in range_bnds for i in species_range):
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
for species in species_info.keys():
	''' The species bounds are adjusted to the range bounds'''
	idx0 = (np.abs(range_bnds - species_info[species]['range_bounds'][0])).argmin()
	idx1 = (np.abs(range_bnds - species_info[species]['range_bounds'][1])).argmin()
	species_info[species]['range_bounds'] = np.asarray([range_bnds[idx0], range_bnds[idx1]])

spec_list = []
for species in species_info.keys():
	if species_info[species]['enabled']:
		spec_list.extend([f"{species_info[species]['short_name']}_R{i}"
					for i in range(1, len(range_bnds))
					if range_bnds[i] > species_info[species]['range_bounds'][0]
					and range_bnds[i] <= species_info[species]['range_bounds'][1]])
		species_info[species]['range_idx'].extend([i
					for i in range(1, len(range_bnds))
					if range_bnds[i] > species_info[species]['range_bounds'][0]
					and range_bnds[i] <= species_info[species]['range_bounds'][1]])

spec_array = np.array(spec_list)
nspec = len(spec_array)
ntrac = nspec + N

# =====================================================================
# Write to namelist
# Change to write to e.g. atm_in namelist?
# Write tracers instead to chem_mech file or my_chem_mech.in in case dir?
#
# =====================================================================

f = open("bin_nl", "w") # change to "a" later!!! -> atm_in?
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
f.write("bin_name_list = ")
for i in range(1, len(r)+1):
	f.write("'num_" + str(i) + "',")
f.write("\n")
f.write("species_name_list = ")
for i in range(len(spec_array)):
	f.write("'" + str(spec_array[i]) + "',")
f.write("\n")
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

# write composition
composition_list = []
implicit_list = []

for species in species_info.keys():
   for i in range(len(species_info[species]['range_idx'])):
      composition_list.append(
         species_info[species]['short_name'] +
         '_R' +
         str(species_info[species]['range_idx'][i]) +
         ' -> ' +
         species_info[species]['composition']
	  )
      implicit_list.append(
         species_info[species]['short_name'] +
         '_R' +
         str(species_info[species]['range_idx'][i])
	  )

for i in range(1,N+1):
   composition_list.append(
      'num_' + str(i) + ' -> H'
    )
   implicit_list.append(
      'num_' + str(i)
   )


# read chem_mech.in

with open('chem_mech.in', 'r') as chem_file:
   lines = chem_file.readlines()

modified_chem = []

for line in lines:
    # write lines from the old file for the new file
	modified_chem.append(line)
    # add species composition
	if 'Solution' in line and not 'End' in line and not 'Classes' in line:
		for i in range(len(composition_list)):
			modified_chem.append(composition_list[i] + '\n')
    # add species for advection
	if 'Implicit' in line and not 'End' in line:
		for i in range(len(implicit_list)):
			modified_chem.append(implicit_list[i] + '\n')
#      modified_lines.append()
# write out to my_chem_mech.in
with open('my_chem_mech.in', 'w') as file:
      file.writelines(modified_chem)
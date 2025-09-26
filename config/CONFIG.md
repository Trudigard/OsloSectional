# Oslo sectional configuration
User input and bin_config.py. These configuration files are used to specify the bin structure, as well as the aerosol species and their properties. bin_config.py is called by cime_config/buildnml during the build process when an oslo_sectional compset is being used.
The general workflow is:

### Input
- **configuration_file.ini** : Meant for user edits to adjust bin configuration or aerosol species and their properties
- **chem_mech.in** : Lists all gas species and their chemical reactions but *without* aerosol tracers

### Processing
- **bin_config.py** : Executed in buildnml, creates files read by the model and puts them in the CASEDIR

### Output
- **sectional_aerosol_properties_nl** : Namelist for the general aerosol properties, included in atm_in (CASEDIR)
- **sectional_aerosol_species_properties * nspecies** : Species specific namelists, one per aerosol species in atm_in (CASEDIR)
- **chem_mech.in** : Edited chem_mech.in including all aerosol tracers (CASEDIR)

## User input options
E.g. dust_oslo_sectional.ini.
Specifies the bin configuration, aerosol species and their properties.

### Bin specs
- **nbin**: Number of bins
- **radius_1**: *Center* radius of the smallest bin [nm]
- **radius_N**: *Center* radius of the largest bin [nm]

### Range specs
- **ranges**: Boolean indicating whether size resolution for chemical species should be lower than total aerosol number concentration
- **range_bounds**: The *boundary* radii for the ranges [nm]. These range bounds get adjusted to the closest bin bound, once they are calculated in bin_config.py

### Species
All species have their own section. bin_config.py will assume that every section (denoted by the [] parentheses) after [RANGE SPECS] is a new species.
Each species needs all the properties specified below:

- **[NAME]**                      : The name for the species, this is only relevant in bin_config.py
- **active** = True/False       : Is this species an active component or not. If set to true, the corresponding tracers are added to chem_mech.in and the oslo_sectional_aerosol_properties_nl namelist in atm_in
- **short_name** = XX           : A short string identifier for each species, e.g. *DU* or *SO4*
- **long_name** = some string   : A long name for the species, e.g. *Dust aerosol*
- **range_bounds** = 0.5, 10000 : The *boundary* radii for this species. These should be the same as two values in range_bounds, in order to avoid species ending up in unexpected bins.
- **composition** = X2Yy3       : Chemical formula for the compound, this is used for the tracer specification in chem_mech.in
- **density** = 0.0             : Density of the aerosol species [kg/m3]
- **molecular_weight** = 0.0    : Molecular weight of the aerosol species in [kg/kmol]
- **mixed** = True/False        : When True, the aerosol is internally mixed with all other aerosol of this kind. False means only externally mixed.


## bin_config.py

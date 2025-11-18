# OsloSectional
Oslo Sectional Aerosol Model

## Suggested names
Keywords: Oslo-Stockholm, Nordic, Flexible, Sectional, Aerosol, Bin, Adaptable
Name suggestions: FANCY (Full Aerosol And Not too much ChemistrY/Flexible Aerosol model of the Nordic CommunitY/...), BRO, CAM (Christinas Aerosol Model :D), TAM (THE Aerosol Model), GOAT (Greatest Of all Aerosol Tools), MIST, MOSS, OSCAR (Oslo-Stockholm Community Aerosol R..)

## config
Contains configuration files for the sectional aerosol model that are used/executed during the build process.
More information in the [config documentation](config/CONFIG.md)

## src_osloaero and src_cam
These folders contain files from OsloAero and CAM that have been slightly changed to make the oslo_sectional code run. However, this should only be a temporary fix to avoid changing CAM code at this stage.

## src
The oslo_sectional source code.
Files that have been significantly edited so far:
* sectional_aerosol_properties_mod.F90
* dust_model.F90: Handles dust emissions (INPUT)
* aero_model.F90: Main calls to objects, reads general namelists

## ../pp_dust_oslo_sectional
* chemistry.F90
* chem_mech.in

## Getting this thing running
Status when the code has not been merged properly :)
Note: This is temporary!! It will be more straight-forward in the (hopefully near) future
1. Clone Christina's CAM fork https://github.com/Trudigard/CAM.git
2. Switch branch (e.g. add-chem has the newest changes)
3. run ./bin/git-fleximod update
4. check .gitmodules, the last entry should be something about oslo_sectional pointing to the "trudigard" fork
5. There is some issue downloading the oslo_sectional repo, so you might have to go to src/chemistry/oslo_sectional change the url -> git remote set-url origin git@github.com:Trudigard/OsloSectional.git
6. Current development is on the add-bin-config branch
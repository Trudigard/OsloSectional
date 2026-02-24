# OsloSectional
Oslo Sectional Aerosol Model (Name to be decided!)

## config
Contains configuration files for the sectional aerosol model that are used/executed during the build process.
More information in the [config documentation](config/CONFIG.md)

## src_osloaero and src_cam
These folders contain files from OsloAero and CAM that have been changed to make the oslo_sectional code run. However, this should only be a temporary fix to avoid changing CAM code at this stage.

## src
The oslo_sectional source code.
Files that have been significantly edited so far:
* sectional_aerosol_properties_mod.F90
* sectional_aerosol_state_mod.F90
* dust_model.F90: Handles dust emissions (INPUT)
* aero_model.F90: Main calls to objects, reads general namelists

## ../pp_dust_oslo_sectional
* chemistry.F90
* chem_mech.in

## Getting this thing running
Note: This setup is temporary
1. Clone Christina's CAM fork https://github.com/Trudigard/CAM.git
2. Switch branch (e.g. sec-stable has the newest changes)
3. run ./bin/git-fleximod update
4. cd into src/chemistry/oslo_sectional
5. git switch sec-stable
6. Compset: SecDust

## branches
- sectional_develop : Reviewed code, currently very outdated
- sec-stable : NOT reviewed but short tests have passed
- stale-* : Old branches, to be deleted soon
- other branches are current working feature branches
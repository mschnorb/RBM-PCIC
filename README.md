# RBM-PCIC

## 1. Modifications to the original code

The original code was downloaded on the hydrology group’s website of the University of Washington (http://www.hydro.washington.edu/Lettenmaier/Models/RBM/download.shtml). The version downloaded is the VIC-RBM v2.2.

Several modifications to the code have been made, both technical modifications, but also modifications on the modeling approach.

### 1.1. Modeling modifications
- Reformulation of the energy fluxes modeling based on Wu et al. (2012) doi:10.1029/2012WR012082. The formulations seem to be more widely used in the community and they also make the conversion to the metric system easier. Indeed, in the original version of RBM, units were in the imperial system.
- The runoff and the baseflow contributions are now considered at each grid cell:
     - The runoff from the snow melt contributes at 0 degree Celsius.
     - The runoff contributes at the air temperature (dbt in the code).
     - The baseflow contributes at the bottom soil layer temperature (Tsoil in the code).
- Three options for the Headwater temperature are coded:
     - Mosheni statistical, original RBM approach.
     - Pseudo linear relation between average annual temperature (Tavg in the code) and Thead.
     - Thead = tsoil, used in the actual version.
- Width, depth and velocity : width and depth are modeled following the power functions designed in Kai et al. 2018. The power functions has been calibrated using the few observations available on the Fraser. Velocity is calculated using the streamflow.
		Width = 6,6588 * Q**0,4967
		Depth = 0,2307 * Q**0,5123
		Velocity = Streamflow / (Width * Depth)
- The depth used for the local contribution of runoff and baseflow is calculated with the same formula as for the main stream.

### 1.2. Technical modifications
- The Network file and the Forcing files are generated with two Python scripts instead of a Perl script. A configuration file is read to run these two scripts, an example is provided. The forcing files are still written in a text file. It is not the best solution as it can require a lot of memory for large watershed. The code could be modified to directly read VIC and RVIC files for example.
- Output data are written in a netcdf file. It is much faster for large watershed and long simulations.
- Added new modules to deal with dates, leap year and time related variables.


## 2. Run a simulation

We quickly describe here how to run a simulation with RBM assuming that you are already using VIC and RVIC. We include a brief description of the files required to set-up RBM.

### 2.1. Run VIC.
The output variables must include:
- OUT_RUNOFF
- OUT_BASEFLOW
- OUT_RUNOFF_SNOW (RUNOFF_SNOW is coded in myVIC version)
- OUT_AIR_TEMP
- OUT_SOIL_TEMP
- OUT_VP
- OUT_SHORTWAVE
- OUT_LONGWAVE
- OUT_DENSITY (not currently used)
- OUT_PRESSURE
- OUT_WIND

### 2.2. Run RVIC.
In this run, each grid cells has to be a pour point as simulated streamflow is required at each grid cell RBM.

### 2.3. Generate Network File.
Generate the network file (BASIN_Network) for RBM with the Python script generate_network.py. You need to run this script only once per watershed. You need a configuration file, an example is provided. You also have to provide a direction file (same as in RVIC) clipped to your watershed and in which the outlet is indicated by the number 9 and grid cells not in the watershed set to missing value.
 
### 2.4. Generate Forcing FIles.
the two forcing files (BASIN_flow and BASIN_heat) for RBM with the Python script generate_forcings.py. You need a configuration file, the same as before. It reads the VIC and RVIC output files. You can setup a subsample of the original simulation.

The first file is the flow forcing file, it includes the following variables:
- tstep : time step
- nnode : cell number
- Qout : output streamflow at tsetp from cell nnode (m**3.s**-1)
- Qrun : runoff that reaches the cell at tstep (m**3.s**-1)
- Qbas : baseflow that reaches the cell at tstep (m**3.s**-1)
- QrunSnow : snowmelt contribution to the runoff that reaches the cell at tstep (m**3.s**-1)
- Velocity : streamflow velocity (m.s**-1)

The second is the heat forcing file, it includes the following variables:
- nnode: cell number
- Air : air temperature (deg C)
- AirAvg : air temperature annual average (deg C). Used only for the estimation of the headwater temperature with the pseudo-linear relation.
- Soil : Botton layer soil temperature (deg C)
- VP : vapor pressure (hPa)
- Short : solar short-wave radiation (W.m**-2)
- Long : solar long-wave radiation (W.m**-2)
- Density : atmospheric density
- Pressure : atmospheric pressure (hPa)
- Wind : wind speed (m.s**-1)

### 2.5. Run RBM.
The Network file and the two forcing files have to be in the execution folder. The output files will be written in this folder.

## Example python configuration file.

RBM-PCIC_MOUTH.cfg :

```
[global]
dir:         /home/gdayon/project/group_writable/vicgl/rbm/fraser/
domain:      rvic.domain_fraser_v2.nc

[parameters]
dir:         /home/gdayon/project/group_writable/vicgl/rbm/fraser/
direction:   rvic.parameters_MOUTH_v2.nc
pour_points: FULL.rvic.prm.MOUTH.nc

[vic_output]
dir:         /scratch/gdayon/projections/output/fraser/MOUTH/ACCESS1-0_rcp45_r1i1p1/flux/
vic_output:  results.nc

[rvic_output]
dir:      /scratch/gdayon/projections/output/fraser/MOUTH/ACCESS1-0_rcp45_r1i1p1/flow_full/
flow:     hist/MOUTH.rvic.h0a.2100-01-01.nc
runoff:   flow_runoff/hist/MOUTH.rvic.h0a.2100-01-01.nc
baseflow: flow_baseflow/hist/MOUTH.rvic.h0a.2100-01-01.nc

[rbm_forcings]
dir:     /scratch/gdayon/projections/output/fraser/MOUTH/ACCESS1-0_rcp45_r1i1p1/rbm/
network: MOUTH_Network
flow:    MOUTH_flow
heat:    MOUTH_heat
dStart:  19450101
dEnd:    20991231
```

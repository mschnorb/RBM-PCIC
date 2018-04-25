# -*- coding: utf-8 -*-

import os
import sys
import math
import numpy as np
from netCDF4 import Dataset
import configparser
Config = configparser.ConfigParser()

# Set arrays for column/row moves using numbers 1-8
# from the direction file
# Work in Python, modified from the initial Perl script.
col_move = col_move=(0, 1, 1, 1, 0,-1,-1,-1)
row_move = row_move=(1, 1, 0,-1,-1,-1, 0, 1)

# QUAD is an array with the obverse directions, e.g.
# 5 comes from 1, 6 comes from 2 ...
quad=(5,6,7,8,1,2,3,4)
print("")

################################################
################################################
############ Open an read configuration file
cfg_file = str(sys.argv[1])

Config.read(cfg_file)
gpath       = Config.get('global', 'dir')
domain      = Config.get('global', 'domain')

ppath       = Config.get('parameters', 'dir')
direction   = Config.get('parameters', 'direction')
pour_points = Config.get('parameters', 'pour_points')

vpath    = Config.get('vic_output', 'dir')
energy   = Config.get('vic_output', 'vic_output')

rpath    = Config.get('rvic_output', 'dir')
flow     = Config.get('rvic_output', 'flow')
runoff   = Config.get('rvic_output', 'runoff')
baseflow = Config.get('rvic_output', 'baseflow')

opath    = Config.get('rbm_forcings', 'dir')
network  = Config.get('rbm_forcings', 'network')
fflow    = Config.get('rbm_forcings', 'flow')
fheat    = Config.get('rbm_forcings', 'heat')
############ End of Open an read configuration file
################################################
################################################

################################################
################################################
############ Open an read files
### Parameters
with Dataset(ppath+direction) as pfic:
   # Velocity
   Velocity = np.array(pfic.variables['velocity'])
   ncells   = np.count_nonzero(~np.isnan(Velocity)) # Number of active cells

   # Read the Velocity grid
   vlat     = np.array(pfic.variables['lat'])
   vlon     = np.array(pfic.variables['lon'])

### Drainage Area
with Dataset(ppath+pour_points) as afic:
   # Local Unit Hydrograph
   UH    = np.array(afic.variables['unit_hydrograph']) # Dims = (timesteps=100, sources, tracers=1)

   outlet_number = np.array(afic.variables['outlet_number'])
   source2outlet_ind = np.array(afic.variables['source2outlet_ind'])

   outlet_lat = np.array(afic.variables['outlet_lat'])
   outlet_lon = np.array(afic.variables['outlet_lon'])

   source_lat = np.array(afic.variables['source_lat'])
   source_lon = np.array(afic.variables['source_lon'])

   # Watershed area
   Area  = np.array(afic.variables['outlet_upstream_area']) # Same lat/lon as Streamflow
   Area  = Area / 1e6 # Convert from m2 to km2
   Depth = 0.93 * Area**0.490  # From Yearsley 2012 (in m)
   Width = 0.08 * Area**0.396  # From Yearsley 2012 (in m)

###### Forcings water
#### Streamflow from Runoff
#with Dataset(rpath+runoff) as wfic:
   #QRunoff = np.array(wfic.variables['streamflow']) # In m3.s-1
   #stime   = np.array(wfic.variables['time'])
   #slat    = np.array(wfic.variables['lat'])
   #slon    = np.array(wfic.variables['lon'])

#### Streamflow from Baseflow
#with Dataset(rpath+baseflow) as wfic:
   #QBaseflow = np.array(wfic.variables['streamflow'])
   
### Streamflow from Runoff and Baseflow
with Dataset(rpath+flow) as wfic:
   QFlow   = np.array(wfic.variables['streamflow']) # In m3.s-1
   stime   = np.array(wfic.variables['time'])
   slat    = np.array(wfic.variables['lat'])
   slon    = np.array(wfic.variables['lon'])


###### Forcings energy and local water fluxes
with Dataset(vpath+energy) as efic:
   # AirTemp, VapPress, Short, Long etc...
   AirTemp    = np.array(efic.variables['AIR_TEMP'])
   SoilTemp   = np.array(efic.variables['SOIL_TEMP'])
   #Ratio      = np.copy(AirTemp)# To delete later...
   VapPress   = np.array(efic.variables['VP'])
   Short      = np.array(efic.variables['SHORTWAVE'])
   Long       = np.array(efic.variables['LONGWAVE'])
   Density    = np.array(efic.variables['DENSITY'])
   Press      = np.array(efic.variables['PRESSURE'])
   Wind       = np.array(efic.variables['WIND'])
   elat       = np.array(efic.variables['lat'])
   elon       = np.array(efic.variables['lon'])
   etime      = efic.variables['time']
   ndays   = len(etime)

   VapPress    = VapPress * 10 # kPa to hPa

   SoilTemp[ SoilTemp == 1e+20 ] = np.nan
   SoilTemp = SoilTemp[:,0,:,:]

   # Runoff and Baseflow
   Runoff     = np.array(efic.variables['RUNOFF'])
   RunoffSnow = np.array(efic.variables['RUNOFF_SNOW'])
   Baseflow   = np.array(efic.variables['BASEFLOW'])
   
   Runoff[     Runoff > 1e+19]     = np.nan
   RunoffSnow[ RunoffSnow > 1e+19] = np.nan
   Baseflow[   Baseflow > 1e+19]   = np.nan

### Grid Area to runoff and baseflow convertion
with Dataset(gpath+domain) as dfic:
   GridArea   = np.array(dfic.variables['area']) # In m2
   glat       = np.array(dfic.variables['lat'])
   glon       = np.array(dfic.variables['lon'])

   ilat = np.where( (min(elat)<=glat) & (glat<=max(elat)) )[0]
   ilon = np.where( (min(elon)<=glon) & (glon<=max(elon)) )[0]
   myGrid     = GridArea[ ilat, :]
   myGrid     = myGrid[:, ilon]

   Runoff     = Runoff     * 1e-3 * myGrid / 86400 # In meters3.s-1
   RunoffSnow = RunoffSnow * 1e-3 * myGrid / 86400 # In meters3.s-1
   Baseflow   = Baseflow   * 1e-3 * myGrid / 86400 # In meters3.s-1

### Compute the Annual temperature
navg    = 365
AirYear = np.copy(AirTemp)

for i in range(ndays):
   if( i < round(navg/2)):
      AirYear[i,:,:] = np.nanmean(AirTemp[0:navg-1,:,:], axis = 0)
   elif( i > (ndays - round(navg/2)) ):
      AirYear[i,:,:] = np.nanmean(AirTemp[ndays-navg-1:ndays-1,:,:], axis = 0)
   else:
      AirYear[i,:,:] = np.nanmean(AirTemp[ i-round(navg/2):i+round(navg/2)+1,:,:], axis = 0)

### Network
nfic  = open(opath+network, 'r')

# The node to write with their lat/lon
allNode, allLat, allLon = [], [], []
Node, nlat, nlon = [], [], []
for lines in nfic:
   data = lines.strip().split()
   if data[0] == 'Node': # Check if we are on a node
      allNode.append(int(data[1]))
      allLat.append(float(data[7]))
      allLon.append(float(data[9]))

      if data[12] != '0.00': # Check if we are not at the end of the stream
         Node.append(int(data[1]))
         nlat.append(float(data[7]))
         nlon.append(float(data[9]))

# Find the Headwaters, Outlet and Contributing nodes
nfic  = open(opath+network, 'r')
lines = nfic.readlines()
HeadCell = [] # Will have nheadwater (= nstream) values at the end
OutCell  = []
OutLat   = []
OutLon   = []
TribCell = []

for l in range(len(lines)):
   line  = lines[l]
   split = line.strip().split()

   if len(split) > 1:
      if split[1] == 'Headwaters':
         # Find the Headwater
         headline = lines[l+2].strip().split()
         headnode = headline[1]
         HeadCell.append(int(headnode))
         
         # Find the Trib
         tribnode = split[4]
         TribCell.append(int(tribnode))

         # Find the cell just upstream Outlet
         outline = lines[l+int(split[0])].strip().split()
         outnode = outline[1]
         OutCell.append(int(outline[1]))
         OutLat.append(float(outline[7]))
         OutLon.append(float(outline[9]))

nfic.close()
nnodes = len(Node)
############ End of Open an read direction file
################################################
################################################

################################################
################################################
############ Write Flow forcing file
iinq  = []   # Indice for inflow and Area
ioutq = []   # Indice for outflows and Area
ivlat = []   # Indice for velocity
ivlon = []   # Indice for velocity
irlat = []   # Indice for Runoff and Baseflow
irlon = []   # Indice for Runoff and Baseflow
myUH  = []   # ...

print("Write Flow forcing file")
ofic = open(opath+fflow, "w")
for t in range(len(stime)):
   #print("Flow forcings file t:",t)
   for n,mynode,mylat,mylon in zip(range(nnodes),Node,nlat,nlon):
      iinq.append([])

      if t == 0: # We need to do that only once
         # Output streamflow
         tlat = np.isin(slat, mylat)
         tlon = np.isin(slon, mylon)
         new_sii  = np.where(tlat & tlon)[0]
         ioutq.append(new_sii[0])

         # Velocity
         new_ivlat = list(vlat).index(mylat)
         new_ivlon = list(vlon).index(mylon)
         ivlat.append(new_ivlat)
         ivlon.append(new_ivlon)
         
         # Runoff and Baseflow
         new_irlat = list(elat).index(mylat)
         new_irlon = list(elon).index(mylon)
         irlat.append(new_irlat)
         irlon.append(new_irlon)
         
         # Unit Hydrograph of the cell
         outlet_ilat = np.isin(outlet_lat, mylat)
         outlet_ilon = np.isin(outlet_lon, mylon)
         myoutlet = outlet_number[np.where(outlet_ilat & outlet_ilon)[0]] # The outlet we are working with

         sources  = np.isin(source2outlet_ind, myoutlet) # All source that contribute to the outlet
         isource  = np.where(sources)[0]

         sources_ilat = np.isin(source_lat[isource], mylat) # Source cell with same Lat/Lon as outlet
         sources_ilon = np.isin(source_lon[isource], mylon)
         iUH          = np.where(sources_ilat & sources_ilon)[0]

         myUH.append(UH[0,iUH,0][0]) # Append to the list of UH coef.


         if mynode in HeadCell: # We have a headwater : Qin = 0.
            ... # Do nothing, Qin = 0.
            
         elif mynode not in TribCell: # Only one stream, the previous one
            iinq[n].append(ioutq[n-1])

         elif mynode in TribCell: # Several streams arrive her : Qin = SUM(Qupstream)
            # First, the upstrean cell
            iinq[n].append(ioutq[n-1])

            # Then the Tributaries
            iTribs = (i for (i,n) in enumerate(TribCell) if n==mynode)
            for i in iTribs:
               ocell = OutCell[i]
               olat  = OutLat[i]
               olon  = OutLon[i]

               nolat = np.isin(slat, olat)
               nolon = np.isin(slon, olon)
               new_sii  = np.where(nolat & nolon)[0]
               iinq[n].append(new_sii[0])

      Qout = QFlow[t,ioutq[n]]
      #Qout = QRunoff[t,ioutq[n]] + QBaseflow[t,ioutq[n]]
         
      Qrun     = myUH[n] * Runoff[t,irlat[n],irlon[n]]
      QrunSnow = myUH[n] * RunoffSnow[t,irlat[n],irlon[n]]
      Qbas     = myUH[n] * Baseflow[t,irlat[n],irlon[n]]
      
      # Area, Depth, Width
      mydepth, mywidth = 0.2307 * Qout**0.4123, 6.6588 * Qout**0.4967
      mydepth, mywidth = Depth[ioutq[n]],Width[ioutq[n]]
      mydepth, mywidth = '{:6.1f}'.format(mydepth), '{:7.1f}'.format(mywidth)

      Test = Qout - Qrun - Qbas
      Test = '{:12.5f}'.format(Test)
      
      Qout,Qrun,QrunSnow,Qbas = '{:12.5f}'.format(Qout), '{:12.5f}'.format(Qrun), '{:12.5f}'.format(QrunSnow), '{:12.5f}'.format(Qbas)
      
      # Velocity
      Vel   = Velocity[ivlat[n],ivlon[n]]
      Vel   = format(Vel, "6.2f")

      tstep = '{:8d}'.format(t+1)
      nnode = '{:8d}'.format(mynode)

      ofic.write(str(tstep)+str(nnode)+str(Qout)+str(Qrun)+str(Qbas)+str(QrunSnow)+str(mydepth)+str(mywidth)+str(Vel)+" ")

ofic.close()
############# End of Write Flow forcing file
#################################################
#################################################

################################################
################################################
############ Write Heat forcing file
ihlat  = []   # Lat indice for energy
ihlon  = []   # Lon indice for energy

print("Write Heat forcing file")
ofic = open(opath+fheat, "w")
for t in range(len(stime)):
   #print("Heat forcings file t:",t)
   for n,mynode,mylat,mylon in zip(range(len(allNode)),allNode,allLat,allLon):
      
      if t == 0:
         # Lat and Lon indexes
         ielat = list(elat).index(mylat)
         ielon = list(elon).index(mylon)
         
         ihlat.append(ielat)
         ihlon.append(ielon)
      
      AirOut   = AirTemp[t,ihlat[n],ihlon[n]]
      AirAvg   = AirYear[t,ihlat[n],ihlon[n]]
      SoilOut  = SoilTemp[t,ihlat[n],ihlon[n]]
      VPOut    = VapPress[t,ihlat[n],ihlon[n]]
      ShortOut = Short[t,ihlat[n],ihlon[n]]
      LongOut  = Long[t,ihlat[n],ihlon[n]]
      DensOut  = Density[t,ihlat[n],ihlon[n]]
      PressOut = Press[t,ihlat[n],ihlon[n]]
      WindOut  = Wind[t,ihlat[n],ihlon[n]]
      
      
      AirOut   = '{:6.1f}'.format(AirOut)
      AirAvg   = '{:6.1f}'.format(AirAvg)
      SoilOut  = '{:6.1f}'.format(SoilOut)
      VPOut    = '{:6.1f}'.format(VPOut)
      ShortOut = '{:10.4f}'.format(ShortOut)
      LongOut  = '{:10.4f}'.format(LongOut)
      DensOut  = '{:6.3f}'.format(DensOut)
      PressOut = '{:7.1f}'.format(PressOut)
      WindOut  = '{:5.1f}'.format(WindOut)
      
      nnode = '{:5d}'.format(mynode)
      
      ofic.write(str(nnode)+str(AirOut)+str(AirAvg)+str(SoilOut)+str(VPOut)+str(ShortOut)+str(LongOut)+str(DensOut)+str(PressOut)+str(WindOut)+"   ")
ofic.close()
############ End of Write Heat forcing file
################################################
################################################

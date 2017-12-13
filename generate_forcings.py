# -*- coding: utf-8 -*-

import os
import numpy as np
from netCDF4 import Dataset  #, date2num, date2index, num2date
#import time

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
############ Open an read files
#ipath = "/home/programmer_analyst/Workspace/rbm/"  # On my laptop
iparam = "/storage/home/gdayon/Workspace/rbm/" # On lynx
iforcs = "/storage/home/gdayon/hydro/BAKER/" # On lynx
opath  = "/storage/home/gdayon/hydro/RBM-PCIC/" # On lynx

### Output files
oflow = opath+"Baker.DA_flow_TG"
oheat = opath+"Baker.DA_heat_TG"

### Parameters
fic   = iparam+"rvic.parameters_baker_v2.nc"        # Fraser file
pfic  = Dataset(fic)

# Velocity
Velocity = np.array(pfic.variables['velocity'])
Velocity = Velocity * 3.2808398950131 # meters.s-1 to feet.s-1
ncells   = np.count_nonzero(~np.isnan(Velocity)) # Number of active cells

# Read the Velocity grid
vlat     = pfic.variables['lat']
vlon     = pfic.variables['lon']

### Drainage Area
fic   = iparam+"FULL.rvic.prm.BAKER.20171106.nc"
afic  = Dataset(fic)
Area  = np.array(afic.variables['outlet_upstream_area']) # Same lat/lon as Streamflow
Area  = Area / 1e6 # Convert from m2 to km2
Depth = 0.93 * Area**0.490  # From Yearsley 2012 (in m)
Width = 0.08 * Area**0.396  # From Yearsley 2012 (in m)
Depth,Width = Depth * 3.2808398950131, Width * 3.2808398950131 # meters to feet

### Forcings water
fic   = iforcs+"flow_full/hist/BAKER.rvic.h0a.2006-01-01.nc"
wfic  = Dataset(fic)

# Streamflow
Streamflow = np.array(wfic.variables['streamflow'])
Streamflow = Streamflow * 35.3146665722226 # meters3.s-1 to feet3.s-1
stime      = wfic.variables['time']
slat       = wfic.variables['lat']
slon       = wfic.variables['lon']

### Forcings energy
fic   = iforcs+"flux/results.nc"
efic  = Dataset(fic)

# AirTemp, VapPress, Short, Long etc...
AirTemp    = np.array(efic.variables['AIR_TEMP'])
SoilTemp   = np.array(efic.variables['SOIL_TEMP'])
VapPress   = np.array(efic.variables['VP'])
Short      = np.array(efic.variables['SHORTWAVE'])
Long       = np.array(efic.variables['LONGWAVE'])
Density    = np.array(efic.variables['DENSITY'])
Press      = np.array(efic.variables['PRESSURE'])
Wind       = np.array(efic.variables['WIND'])
etime      = efic.variables['time']
elat       = efic.variables['lat']
elon       = efic.variables['lon']

VapPress    = VapPress * 10
Short, Long = Short * 0.238845897e-3, Long * 0.238845897e-3 # W.m-2 to Kcal.s-1.m-2

SoilTemp[ SoilTemp == 1e+20 ] = np.nan
SoilTemp    = np.nanmean(SoilTemp, axis = 1)

### Network
fic   = opath+"Baker_Network"
nfic  = open(fic, 'r')

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
nfic  = open(fic, 'r')
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

print("Write Flow forcing file")
ofic = open(oflow, "w")
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

         
         if mynode in HeadCell: # We have a headwater : Qin = Qout
            iinq[n].append(ioutq[n])
            
         elif mynode not in TribCell: # Only one stream, the previous one
            iinq[n].append(ioutq[n-1])

         elif mynode in TribCell: # Several streams arrive her : Qin = SUM(Qupstream)
            # First, the upstrean cell
            iinq[n].append(ioutq[n-1])

            # Then the Tributaries
            iTribs = (i for (i,n) in enumerate(TribCell) if n==mynode)
            for i,ocell,olat,olon in zip(iTribs,OutCell,OutLat,OutLon):
               nolat = np.isin(slat, olat)
               nolon = np.isin(slon, olon)
               new_sii  = np.where(nolat & nolon)[0]
               iinq[n].append(new_sii[0])

      Qin  = Streamflow[t,ioutq[n]]
      Qout = Streamflow[t,ioutq[n]]
      Qdif = 0.

      Qin,Qout,Qdif = '{:10.1f}'.format(Qin), '{:10.1f}'.format(Qout), '{:6.1f}'.format(Qdif)

      # Area, Depth, Width
      mydepth,mywidth  = Depth[ioutq],Width[ioutq]
      mydepth,mywidth  = '{:6.1f}'.format(mydepth[0]), '{:7.1f}'.format(mywidth[0])

      # Velocity
      Vel   = Velocity[ivlat[n],ivlon[n]]
      Vel   = format(Vel, "6.2f")

      tstep = '{:5d}'.format(t+1)
      nnode = '{:5d}'.format(mynode)

      ofic.write(str(tstep)+str(nnode)+str(Qin)+str(Qout)+str(Qdif)+str(mydepth)+str(mywidth)+str(Vel)+"     ")

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
ofic = open(oheat, "w")
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
      SoilOut  = SoilTemp[t,ihlat[n],ihlon[n]]
      VPOut    = VapPress[t,ihlat[n],ihlon[n]]
      ShortOut = Short[t,ihlat[n],ihlon[n]]
      LongOut  = Long[t,ihlat[n],ihlon[n]]
      DensOut  = Density[t,ihlat[n],ihlon[n]]
      PressOut = Press[t,ihlat[n],ihlon[n]]
      WindOut  = Wind[t,ihlat[n],ihlon[n]]
      
      #VPOut,PressOut   = VPOut*10, PressOut*10
      #ShortOut,LongOut = ShortOut * 0.238845897e-3, LongOut * 0.238845897e-3

      AirOut   = '{:6.1f}'.format(AirOut)
      SoilOut  = '{:6.1f}'.format(SoilOut)
      VPOut    = '{:6.1f}'.format(VPOut)
      ShortOut = '{:7.4f}'.format(ShortOut)
      LongOut  = '{:7.4f}'.format(LongOut)
      DensOut  = '{:6.3f}'.format(DensOut)
      PressOut = '{:7.1f}'.format(PressOut)
      WindOut  = '{:5.1f}'.format(WindOut)
      
      nnode = '{:5d}'.format(mynode)
      ofic.write(str(nnode)+str(AirOut)+str(SoilOut)+str(VPOut)+str(ShortOut)+str(LongOut)+str(DensOut)+str(PressOut)+str(WindOut)+"     ")
      #ofic.write(str(nnode)+str(AirOut)+str(VPOut)+str(ShortOut)+str(LongOut)+str(DensOut)+str(PressOut)+str(WindOut)+" ")

ofic.close()
############ End of Write Heat forcing file
################################################
################################################
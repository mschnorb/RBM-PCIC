# -*- coding: utf-8 -*-

import os
import numpy as np
from netCDF4 import Dataset  #, date2num, date2index, num2date

# Set arrays for column/row moves using numbers 1-8
# from the direction file
col_move=(0,1,1,1,0,-1,-1,-1)
row_move=(1,1,0,-1,-1,-1,0,1)

# QUAD is an array with the obverse directions, e.g.
# 5 comes from 1, 6 comes from 2 ...
quad=(5,6,7,8,1,2,3,4)
print("")

################################################
################################################
############ Open an read direction file
ipath = "/storage/home/gdayon/Workspace/rbm/"       # On lynx
ific  = ipath+"rvic.parameters_baker_v2.nc"       # Test file
fic   = Dataset(ific)

opath = "/storage/home/gdayon/hydro/RBM-PCIC/"       # On lynx
ofic  = opath+"Baker_Network"

# Flow direction
Flow_dir = fic.variables['Flow_Direction']
Flow_dir = np.matrix(Flow_dir)
Flow_dir[Flow_dir < 0] = np.nan
ncells   = np.count_nonzero(~np.isnan(Flow_dir)) # Number of active cells

# Flow distance
Flow_dis = fic.variables['Flow_Distance']
Flow_dis = np.matrix(Flow_dis)

# Read the grid
lat      = np.array( fic.variables['lat'] )
lon      = np.array( fic.variables['lon'] )
nlat     = len(lat)
nlon     = len(lon)
nall     = nlat * nlon # Total number of cells

if(lat[0] > lat[1]): # Check if lat is in ascending order
   lat      = lat[::-1]
   Flow_dir = np.flipud(Flow_dir)
   Flow_dis = np.flipud(Flow_dis)

############ End of Open an read direction file
################################################
################################################

################################################
################################################
############ Define some array
# Indexes of active cells
icells          = np.where(~np.isnan(Flow_dir))
cell_n01d       = np.arange(ncells)

cell_n0         = -1*np.ones(Flow_dir.shape, dtype=int)
cell_n0[icells] = np.arange(ncells)

# Rows and cols of active cells
grid     = np.indices(Flow_dir.shape, dtype=int)
cell_row = np.reshape(grid[0,icells[0],icells[1]], ncells)
cell_col = np.reshape(grid[1,icells[0],icells[1]], ncells)

# Lat and Lon of active cells
lat2d    = np.resize(lat, (Flow_dir.T).shape).T
lon2d    = np.resize(lon, Flow_dir.shape)
cell_lat = np.reshape(lat2d[icells[0],icells[1]], ncells)
cell_lon = np.reshape(lon2d[icells[0],icells[1]], ncells)

# Length of active cells
Length   = np.reshape(Flow_dis[icells[0],icells[1]], ncells)
Length   = np.array(Length.tolist()[0])
Length   = Length * 0.0006213712 # Meters to miles
############ End of Define some array
################################################
################################################

################################################
################################################
############ Find the headwaters and the outlet
# Look at all the cells in the domain
head_cell    = np.copy(Flow_dir) # Array to find headwater and outlet
ncell_incell = np.zeros(Flow_dir.shape, dtype=int) # Numbers of cells contributing to the cell
lcell_incell = [[] for x in range(ncells)] # List of cells contributing to the cell
for nr in range(1,nlat):
   for nc in range(1,nlon):
      if 8 >= Flow_dir[nr,nc] >= 1:
         # Direction and number
         direction = int(Flow_dir[nr,nc])
         n0        = cell_n0[nr,nc]

         # Dowstream (ds) cell
         nr_ds = grid[0,nr,nc] + row_move[direction-1]
         nc_ds = grid[1,nr,nc] + col_move[direction-1]
         n0_ds = cell_n0[nr_ds,nc_ds]
         
         # Contribution to the ds cell
         ncell_incell[nr_ds, nc_ds] += 1
         lcell_incell[n0_ds].append(n0)

         # Set to nan, this downstream cell is not a headwater
         head_cell[nr_ds,nc_ds] = np.nan
         lcell_incell

      elif Flow_dir[nr,nc] == 9:
         outlet_row = grid[0,nr,nc]
         outlet_col = grid[1,nr,nc]
         outlet_n0  = cell_n0[outlet_row,outlet_col]

# Identify rows and cols of headwaters
head_icells = np.where(~np.isnan(head_cell))
nhead       = len(head_icells[0])

head_row    = np.reshape(grid[0, head_icells[0],head_icells[1]], nhead)
head_col    = np.reshape(grid[1, head_icells[0],head_icells[1]], nhead)
############ End of Find the headwaters and the outlet
################################################
################################################

################################################
################################################
############ Read the network
############ Find the headwaters, the outlet and the main stream
stream_order  = np.ones(nhead)
stream_length = np.ones(nhead)
stream_cells  = [[] for x in range(nhead)] # List of cells along a stream

conf_row      = np.zeros(nhead, dtype=int) # Row of confluence cells
conf_col      = np.zeros(nhead, dtype=int) # Col of confluence cells
contrib_to    = np.zeros(nhead, dtype=int) # 0 for the main stream

nhead_incell  = np.zeros(Flow_dir.shape)    # Numbers of headwaters contributing to the cell
lhead_incell  = [[] for x in range(ncells)] # List of headwaters contributing to the cell

max_length = -9999.
for nh in range(nhead):
   nr       = head_row[nh]
   nc       = head_col[nh]
   flow_dir = int(Flow_dir[nr,nc])

   # Follow the network to the outlet
   length   = 0
   while 8 >= flow_dir >= 1:
      length += 1
      
      # Cell number
      n0 = cell_n0[nr,nc]

      # Dowstream (ds) cell
      nr_ds = grid[0,nr,nc] + row_move[flow_dir-1]
      nc_ds = grid[1,nr,nc] + col_move[flow_dir-1]
      
      # A new headwater contributes
      nhead_incell[nr,nc] += 1
      lhead_incell[n0].append(nh)
      
      # Append the list of cells in the stream
      stream_cells[nh].append(n0)
      
      # Let's move
      nr = nr_ds
      nc = nc_ds
      flow_dir = int(Flow_dir[nr,nc])
      
   stream_cells[nh].append(outlet_n0) # Add the outlet

   if length >= max_length:
      max_length  = length
      main_stream = nh

   stream_length[nh] = length

# The confluence of the main stream is the outlet
conf_row[main_stream] = outlet_row
conf_col[main_stream] = outlet_col
lorder_stream         = [[main_stream]] # The main stream order is 1. Push it in the list
print("Main stream :",main_stream,"Lenght:",max_length)
print("")
############ End of Read the network
################################################
################################################

################################################
################################################
############ Stream level, confluence nodes
remain_stream = nhead  # Init. (voir -1, a verifier)
stream_order  = 1      # Init.
while remain_stream > 0:
   remain_stream -= 1
   lorder_stream.append([])

   for nst in lorder_stream[stream_order-1]:

      start_row = conf_row[nst]
      start_col = conf_col[nst]
      start_n0  = cell_n0[start_row, start_col]
      in0_inseg = np.where(stream_cells[nst] == start_n0)

      mymainseg = stream_cells[nst][0:int(in0_inseg[0])]
      nseg      = len(mymainseg)

      for ns in range(nseg-1,0,-1):
         n0 = mymainseg[ns]
         nr = cell_row[n0]
         nc = cell_col[n0]
         
         nc_incell = ncell_incell[nr,nc] # Cells contributing to the actual cell

         if( nc_incell > 1): # If true, we found a confluence !
            myus_cell = mymainseg[ns-1]

            # Let's find the longest stream
            for new in (np.where(lcell_incell[n0] != myus_cell)[0] ): # Loop over stream different from the main

               us_cell = lcell_incell[n0][int(new)] # Cell upstream of the confluence
               us_row  = cell_row[us_cell]
               us_col  = cell_col[us_cell]

               nhd     = nhead_incell[us_row,us_col]      # Number of headwater contributing to the cell
               if(nhd > 1):
                  max_length = -9999
                  for nbr in lhead_incell[us_cell]:
                     mynew_stream = stream_cells[nbr]
                     istart       = int(np.where(mynew_stream == us_cell)[0])
                     length       = len(mynew_stream[0:istart+1])
                     
                     if(length >= max_length):
                        max_length = length
                        new_branch = nbr
                  
               else: # Only one headwater, it is our stream
                  # This is our new branch :)
                  new_branch = lhead_incell[us_cell][0]

               # Save the confluence cell and the 
               conf_row[new_branch] = nr
               conf_col[new_branch] = nc
               contrib_to[new_branch] = nst
               
               # Push the stream to the list
               lorder_stream[stream_order].append(new_branch)
                  
               remain_stream -= 1

   # Look for the upper level streams
   stream_order += 1
############ Stream level, confluence nodes
################################################
################################################

################################################
################################################
############ Go trhough the all network to reorder it
stream_n0_2write = np.zeros(nhead)
node_n0_2write   = np.zeros(ncells+nhead-1, dtype=int)
node_st_2write   = np.zeros(ncells+nhead-1, dtype=int)
ist = 0
ind = 0
for lv in range(len(lorder_stream)-1,-1,-1):

   for st in lorder_stream[lv]:
      stream_n0_2write[ist] = st
      ist += 1
      
      row = conf_row[st]
      col = conf_col[st]
      conf_cell = cell_n0[row,col]
      
      mystream  = stream_cells[st]
      istart    = int(np.where(mystream == conf_cell)[0])
      mystream  = mystream[0:istart+1]
      
      for nd in mystream[0:istart+1]:
         node_n0_2write[ind] = nd
         node_st_2write[ind] = st
         ind += 1
############ End of Go trhough the all network to reorder it
################################################
################################################

################################################
################################################
############ Write the Network file
print(ofic)
f = open(ofic, "w")

f.write("Networtk file for BAKER test case\n")
f.write("Baker.DA_flow_TG\n") # Forcing files
f.write("Baker.DA_heat_TG\n") # Forcing files
f.write("19890101".rjust(10)+"20051231".rjust(10)+"1".rjust(10)+"\n") # Start - End dates / Timesteps
f.write(str(nhead).rjust(10)+str(ncells-1).rjust(10)+str(ncells+nhead-1).rjust(10)+"FALSE".rjust(21)+"\n")

for lv in range(len(lorder_stream)-1,-1,-1):

   for st in lorder_stream[lv]:
      # My stream and its new number
      mystream  = stream_cells[st]      
      n0_stream = np.where(stream_n0_2write == st)[0]  + 1

      # Dowstream stream and its new number
      ctb       = contrib_to[st]
      n0_contri = np.where(stream_n0_2write == ctb)[0] + 1
      n0_contri = n0_contri[0]

      # Confluence cell and its new number
      row = conf_row[st]
      col = conf_col[st]
      conf_cell = cell_n0[row,col]
      trib_cell = list(set(np.where(node_n0_2write == conf_cell)[0]) & set(np.where(node_st_2write == ctb)[0]))

      # Cuidado for the main stream
      if lv == 0:
         n0_contri = 0
         ctb       = st
         trib_cell = list(set(np.where(node_n0_2write == conf_cell)[0]) & set(np.where(node_st_2write == ctb)[0]))

      trib_cell[0] = trib_cell[0] + 1

      # Length of my stream (in cell)
      istart    = int(np.where(mystream == conf_cell)[0])
      mystream  = mystream[0:istart+1]
      lmystream = str(len(mystream))

      # Length of my stream (in km)
      Total_length = sum(Length[mystream])

      f.write(lmystream.rjust(5)+" Headwaters"+str(n0_stream[0]).rjust(4)+" TribCell"+str(trib_cell[0]).rjust(6)+"  Headwaters"+str(n0_contri).rjust(8)+"R.M. =".rjust(15)+str(format(Total_length, "5.2f")).rjust(10)+"\n")
      f.write("   23.00    12.10   0.2500     0.50   0.1000\n")

      myLength = Total_length
      for nd in mystream[0:istart+1]:
         # My cell and its new number
         row = str(cell_row[nd] + 1)
         col = str(cell_col[nd] + 1)

         n0_cell = list(set(np.where(node_n0_2write == nd)[0]) & set(np.where(node_st_2write == st)[0]))
         n0_cell = str(n0_cell[0] + 1)

         # Lat - Lon of my cell
         lat = str(cell_lat[nd])
         lon = str(cell_lon[nd])

         # Length
         myLength = abs(myLength - Length[nd]) # In case it goes under zero because of truncation

         f.write("Node"+n0_cell.rjust(6)+" Row"+row.rjust(6)+" Column"+col.rjust(6)+"  Lat"+lat.rjust(9)+" Long"+lon.rjust(11)+" R.M. =")
         f.write(str(format(myLength, "5.2f")).rjust(10))
         f.write(str(2).rjust(5)+"\n")

f.close()
############ End of Write the Network file
################################################
################################################
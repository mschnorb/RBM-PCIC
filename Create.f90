SUBROUTINE CREATE  !!!(all_lat,all_lon)

use netcdf
use Block_Network
use Block_Netcdf

Implicit NONE

integer :: i,j,loc
real    :: mylat, mylon
real    :: res = 0.0625
character (len=10) :: yyyymmdd
character (len=21) :: time_units


call check( nf90_create( adjustl(out_file_ncdf), 0, ncid) )

ncells = heat_cells
nlat   = (maxval(llat) - minval(llat)) / res + 1
nlon   = (maxval(llon) - minval(llon)) / res + 1
nsegs  = maxval(ndelta)

allocate( cells(ncells) )
allocate( segs(nsegs) )

!!! Create lat and lon vectors (we assume a resolution here...) 
cells(:) = (/(i, i = 1,ncells, 1)/)
segs(:) = (/(i, i = 1,nsegs, 1)/)

!!! Define the dimensions
call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, t_dimid) )
call check( nf90_def_dim(ncid, 'cell', ncells, cell_dimid) )
call check( nf90_def_dim(ncid, 'seg', nsegs, seg_dimid) )

!!! Define the coordinate variables
call check( nf90_def_var(ncid, 'time', NF90_INT, t_dimid, tid) )
call check( nf90_def_var(ncid, 'cell', NF90_INT, cell_dimid, cellid) )
call check( nf90_def_var(ncid, 'seg', NF90_INT,  seg_dimid, segid) )

!!! Assign units attributes to coordinate variables
call check( nf90_put_att(ncid, cellid, 'units', '-') )
call check( nf90_put_att(ncid, cellid, 'long_name', 'Grid cell number') )

call check( nf90_put_att(ncid, segid, 'units', '-') )
call check( nf90_put_att(ncid, segid, 'long_name', 'Segments in the grid cell') )
call check( nf90_put_att(ncid, segid, 'infos', 'The smaller the integer, the more upstream you are') )

write(yyyymmdd, '(i4,A1,i2.2,A1,i2.2)') start_year, '-', start_month, '-', start_day
write(time_units, '(A11, A10)') 'days since ', yyyymmdd 
call check( nf90_put_att(ncid, tid, 'units', time_units) )
call check( nf90_put_att(ncid, tid, 'long_name', 'time') )

!!! The dimids array
dimids = (/ cell_dimid, seg_dimid, t_dimid /)
     
!!! Define the netCDF variables for the pressure and temperature data
call check( nf90_def_var(ncid, 'Ts', NF90_REAL, dimids, Tsid) )
call check( nf90_def_var(ncid, 'Qs', NF90_REAL, dimids, Qsid) )
call check( nf90_def_var(ncid, 'Tas', NF90_REAL, dimids, Tasid) )
call check( nf90_def_var(ncid, 'Tsoil', NF90_REAL, dimids, Tsoilid) )
call check( nf90_def_var(ncid, 'Width', NF90_REAL, dimids, Widthid) )
call check( nf90_def_var(ncid, 'Depth', NF90_REAL, dimids, Depthid) )
call check( nf90_def_var(ncid, 'Qrun', NF90_REAL, dimids, Qrunid) )
call check( nf90_def_var(ncid, 'Qbas', NF90_REAL, dimids, Qbasid) )
call check( nf90_def_var(ncid, 'Qrunsnow', NF90_REAL, dimids, Qrunsnowid) )


!!! Assign units attributes to the netCDF variables
call check( nf90_put_att(ncid, Tsid, 'units', 'degC') )
call check( nf90_put_att(ncid, Tsid, '_FillValue', 1.e+20) )
call check( nf90_put_att(ncid, Tsid, 'long_name', 'Stream temperature at grid cell') )

call check( nf90_put_att(ncid, Qsid, 'units', 'm3.s-1') )
call check( nf90_put_att(ncid, Qsid, '_FillValue', 1.e+20) )
call check( nf90_put_att(ncid, Qsid, 'long_name', 'Stream flow at grid cell') )

call check( nf90_put_att(ncid, Tasid, 'units', 'degC') )
call check( nf90_put_att(ncid, Tasid, '_FillValue', 1.e+20) )
call check( nf90_put_att(ncid, Tasid, 'long_name', 'Air temperature at grid cell') )

call check( nf90_put_att(ncid, Tsoilid, 'units', 'degC') )
call check( nf90_put_att(ncid, Tsoilid, '_FillValue', 1.e+20) )
call check( nf90_put_att(ncid, Tsoilid, 'long_name', 'Soil temperature at grid cell') )

call check( nf90_put_att(ncid, Widthid, 'units', 'm') )
call check( nf90_put_att(ncid, Widthid, '_FillValue', 1.e+20) )
call check( nf90_put_att(ncid, Widthid, 'long_name', 'River width at grid cell') )

call check( nf90_put_att(ncid, Depthid, 'units', 'm') )
call check( nf90_put_att(ncid, Depthid, '_FillValue', 1.e+20) )
call check( nf90_put_att(ncid, Depthid, 'long_name', 'River depth at grid cell') )

call check( nf90_put_att(ncid, Qrunid, 'units', 'm3.s-1') )
call check( nf90_put_att(ncid, Qrunid, '_FillValue', 1.e+20) )
call check( nf90_put_att(ncid, Qrunid, 'long_name', 'Runoff contribution at grid cell') )

call check( nf90_put_att(ncid, Qbasid, 'units', 'm3.s-1') )
call check( nf90_put_att(ncid, Qbasid, '_FillValue', 1.e+20) )
call check( nf90_put_att(ncid, Qbasid, 'long_name', 'Baseflow contribution at grid cell') )

call check( nf90_put_att(ncid, Qrunsnowid, 'units', 'm3.s-1') )
call check( nf90_put_att(ncid, Qrunsnowid, '_FillValue', 1.e+20) )
call check( nf90_put_att(ncid, Qrunsnowid, 'long_name', 'Runoff from snow and ice melt contribution at grid cell') )

!!! End define mode
call check( nf90_enddef(ncid) )

!!! Write the coordinate variable data
call check( nf90_put_var(ncid, cellid, cells) )
call check( nf90_put_var(ncid, segid, segs) )


!!! These settings tell netcdf to write one timestep of data
count = (/ ncells, nsegs, 1 /)
start = (/ 1, 1, 1 /)

END SUBROUTINE CREATE

Module Block_Netcdf

! Intermediate variables to write NETCDF file
   integer, dimension(:), allocatable  :: ivlon, ivlat   ! Size ncell
   real, dimension(:,:), allocatable   :: Qs_vec, Ts_vec, Tas_vec, Tsoil_vec, Width_vec, Depth_vec, Qrun_vec, Qbas_vec, Qrunsnow_vec ! Size ncell*nlvls
   real, dimension(:,:,:), allocatable :: Qs_2d, Ts_2d   ! Netcdf dimensions

! NETCDF Output file variables
   integer                            :: ncells, nlat, nlon, nsegs, nlvls
   integer, dimension(:), allocatable :: nnodes     ! Size ncell
   real, dimension(:), allocatable    :: llat, llon ! Size ncell
   integer, dimension(:), allocatable :: time, segs, lvls, cells ! Netcdf sizes
   real, dimension(:), allocatable    :: lats, lons        ! Netcdf sizes

! NETCDF Output file informations
   character (len=200)   :: out_file_ncdf
   integer, dimension(3) :: dimids, start, count
   integer :: ndims = 3
   integer :: ncid, Tsid, Qsid, Tasid, Tsoilid, Widthid, Depthid, Qrunid, Qbasid, Qrunsnowid
   integer :: cellid, latid, lonid, segid, lvlid, tid
   integer :: cell_dimid, lat_dimid, lon_dimid, seg_dimid, lvl_dimid, t_dimid
   integer :: csegs ! Compteur for levels writting

end Module Block_Netcdf
Subroutine BEGIN(param_file,spatial_file)

use netcdf
use Block_Energy
use Block_Hydro
use Block_Network
use Block_Netcdf

Implicit None
    
    character (len=8)   :: end_date,start_date     
    character (len=8)   :: clat
    character (len=10)  :: clon
    character (len=200) :: param_file,source_file,spatial_file

    integer :: Julian
    integer :: head_name,trib_cell
    integer :: jul_start,main_stem,nyear1,nyear2,nc,ncell,nseg
    integer :: ns_max_test,node,ncol,nrow,nr,cum_sgmnt

    logical :: first_cell,source
!
    real :: nndlta
    real :: rmile0,rmile1,xwpd
!
!   Mohseni parameters, if used
!
!
!
!     Card Group I
!
read(90,*) start_date,end_date
read(start_date,'(i4,2i2)') start_year,start_month,start_day
read(end_date,  '(i4,2i2)') end_year,end_month,end_day
nyear1 = start_year
nyear2 = end_year

write(*,'(2(2x,i4,2i2))')  &
 start_year,start_month,start_day,end_year,end_month,end_day
!
!     Establish the Julian day for which simulations begin
!
jul_start = Julian(start_year,start_month,start_day)
!
read(90,*) nreach,flow_cells,heat_cells,source
!
! Allocate dynamic arrays (messy)
!
 allocate(ndelta(heat_cells))
 allocate(mu(nreach))
 allocate(alphamu(nreach))
 allocate(beta(nreach))
 allocate(gmma(nreach))
 allocate (smooth_param(nreach))
 allocate(dx(heat_cells))
 allocate(no_celm(nreach))
 no_celm=0
 allocate(no_cells(nreach))
 no_cells=0
 allocate(no_tribs(heat_cells))
 no_tribs=0
 allocate(trib(heat_cells,10))
 trib=0
 allocate (conflnce(heat_cells,10))
 conflnce=0
 allocate(reach_cell(nreach,ns_max)) 
 allocate(head_cell(nreach))
 allocate(segment_cell(nreach,ns_max))
 allocate(x_dist(nreach,0:ns_max))
 
!!!!!! Allocate dynamic arrays (mine, much better :p)
! Network information
allocate(nnodes(heat_cells))
allocate(llat(heat_cells))
allocate(llon(heat_cells))
 
!
!     Start reading the reach date and initialize the reach index, NR
!     and the cell index, NCELL
!
ncell=0
!
ns_max_test=-1
!
!     Card Group IIb. Reach characteristics
!
do nr = 1,nreach
   ! Initialize NSEG, the total number of segments in this reach

   nseg = 0
   write(*,*) ' Starting to read reach ',nr

   ! Read the number of cells in this reach, the headwater #,
   ! the number of the cell where it enters the next higher order stream,
   ! the headwater number of the next higher order stream it enters, and
   ! the river mile of the headwaters.
   read(90,'(i5,11x,i4,10x,i5,15x,i5,15x,f15.2,i5)') no_cells(nr) &
      ,head_name,trib_cell,main_stem,rmile0

   ! Set boundary distance
   x_dist(nr,0)=rmile0

   ! If this is reach that is tributary to cell TRIB_CELL, give it the
   ! pointer TRIB(TRIB_CELL) the index of this reach for further use.
   ! Also keep track of the total number of tributaries for this cell
   ! Add main_stem condition (Feb. 2018, DAYON G.), watch outlet cell
   if (trib_cell.gt.0 .and. main_stem.ne.0) then
      no_tribs(trib_cell) = no_tribs(trib_cell) + 1
      trib(trib_cell,no_tribs(trib_cell)) = nr
   end if

   ! Reading Mohseni parameters for each headwaters (UW_JRY_2011/06/18)
   read(90,*) alphaMu(nr),beta(nr), &
              gmma(nr),mu(nr),smooth_param(nr)

   ! Reading Reach Element information
   first_cell=.true.
  
   do nc = 1,no_cells(nr)
      ncell = ncell+1

      ! Relate the cell number from the linear list to the local cell number
      reach_cell(nr,nc) = ncell

      !!! Read the data for point sources

      if (source) then
         ! Place holder for point source input
      end if 

      ! The headwaters index for each cell in this reach is given
      ! in the order the cells are read

      ! Card Type 3. Cell indexing #, Node # Row # Column Lat Long RM
      !
      ! Variable ndelta read in here.  At present, number of elements
      ! is entered manually into the network file (UW_JRY_2011/03/15)
      read(90,'(5x,i5,5x,i5,8x,i5,6x,a8,6x,a10,7x,f15.2,f5.0)')  &
              node,nrow,ncol,clat,clon,rmile1,ndelta(ncell)

      ! Keep some variables for netcdf
      nnodes(ncell) = node
      read(clat, *) llat(ncell)
      read(clon, *) llon(ncell)
      
      ! Set the number of segments of the default, if not specified
      if (ndelta(ncell).lt.1) ndelta(ncell)=n_default

      if(first_cell) then
         first_cell=.false.
         head_cell(nr)=ncell
         x_dist(nr,0)=rmile0
      end if

      ! Added variable ndelta (UW_JRY_2011/03/15)
      dx(ncell) = (rmile0 - rmile1) / ndelta(ncell)
      rmile0    = rmile1
      nndlta    = 0
      200 continue
      
      nndlta    = nndlta + 1
      nseg      = nseg + 1
      segment_cell(nr,nseg)   = ncell
      segment_cell(nr,nseg+1) = ncell
      
      write(*,*) 'nndlta --',nr,nndlta,nseg,ncell,segment_cell(nr,nseg),dx(ncell),rmile1
      
      x_dist(nr,nseg)=x_dist(nr,nseg-1)-dx(ncell)

      ! Write Segment List for mapping to temperature output (UW_JRY_2008/11/19)
      open(22,file=TRIM(spatial_file),status='unknown') ! (changed by WUR_WF_MvV_2011/01/05)
      write(22,'(4i6,1x,a8,1x,a10,f5.0)') nr,ncell,nrow,ncol,clat,clon,nndlta

      !  Added variable ndelta  (UW_JRY_2011/03/15)
      !
      if(nndlta.lt.ndelta(ncell)) go to 200  
      no_celm(nr)=nseg
      segment_cell(nr,nseg)=ncell
      x_dist(nr,nseg)=rmile1

   end do ! End of loop on segment (nc=1,no_cells(nr)

   !
   ! If this is a reach that is tributary to another, set the confluence cell to the previous 
   ! cell. This is necessary because the last cell in the reach has the same cell number
   ! as that of the cell it enters. This is to account for the half portion of the cell the
   ! the parcel traverses to the center of the grid.
   if (trib_cell .gt. 0) then
      conflnce(trib_cell,no_tribs(trib_cell)) = ncell
   end if

   if(ns_max_test.lt.nseg) ns_max_test = nseg

end do ! End of loop on reach (nr=1,nreach)

if(ns_max_test.gt.ns_max) then
  write(*,*) 'RBM is terminating because'
  write(*,*) 'NS_MAX exceeded. Change NS_MAX in Block_Network to: ',ns_max_test
  stop
end if

nwpd = 1
xwpd = nwpd
dt_comp = 86400./xwpd
!
!     ******************************************************
!                         Return to RMAIN
!     ******************************************************
!
900 continue
!
!
end subroutine BEGIN
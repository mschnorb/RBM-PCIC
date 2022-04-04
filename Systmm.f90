SUBROUTINE SYSTMM(temp_file,param_file)

use netcdf

use Block_Energy
use Block_Hydro
use Block_Network
use Block_Netcdf

use Date_Utility

Implicit None

character (len=200):: temp_file
character (len=200):: param_file

integer          :: i,j,k
integer          :: ncell,nncell,ncell0,nc_head,no_flow,no_heat
integer          :: nc,nd,ndd,nm,nr,ns
integer          :: nr_trib,ntribs
integer          :: nrec_flow,nrec_heat
integer          :: n1,n2,nnd,nobs,nyear,nd_year,ntmp
integer          :: npart,nseg,nx_s,nx_part,nx_head
integer          :: year, month, day, status

! Indices for lagrangian interpolation
integer              :: njb,npndx,ntrp
integer, dimension(2):: ndltp=(/-1,-2/)
integer, dimension(2):: nterp=(/2,3/)

real             :: dt_calc,dt_total,hpd,q_dot,q_surf,z,w,x_calc
real             :: Q_dstrb,Q_inflow,Q_outflow,Q_ratio,Q_trb,Q_trb_sum,Q_base,Q_runoff,Q_runoffsnow,Q_local
real             :: T_trb_load,T_local_load
!real             :: T_dstrb,T_dstrb_load,T_runoff_load,T_base_load
real             :: rminsmooth
real             :: T_0,T_dist,Thh,Tseuil,T_tribs,T_base,T_runoff,T_snow,T_mohs,T_local
real             :: x,xd_year,xwpd
real             :: tntrp
real             :: dt_ttotal
real,dimension(4):: ta,xa
real             :: local_w, local_z

real,dimension(:),allocatable     :: T_head,T_smth,T_trib,T_smth_local

logical:: DONE

! Allocate the arrays
allocate (temp(nreach,0:ns_max,2))
allocate (T_head(nreach))
allocate (T_smth(nreach))
allocate (T_trib(nreach))
allocate (T_smth_local(heat_cells))
allocate (depth(heat_cells))
allocate (Q_in(heat_cells))
allocate (Q_out(heat_cells))
allocate (Q_run(heat_cells))
allocate (Q_bas(heat_cells))
allocate (Q_runsnow(heat_cells))
allocate (Q_diff(heat_cells))
allocate (Q_trib(nreach))
allocate (width(heat_cells))
allocate (u(heat_cells))
allocate (dt(2*heat_cells))
allocate (dbt(heat_cells))
allocate (tavg(heat_cells))
allocate (tsoil(heat_cells))
allocate (ea(heat_cells))
allocate (Q_ns(heat_cells))
allocate (Q_na(heat_cells))
allocate (press(heat_cells))
allocate (wind(heat_cells))

! Initialize some arrays
dt_part   = 0.
x_part    = 0.
no_dt     = 0
nstrt_elm = 0
temp      = 0.5

! Initialize headwaters temperatures
T_head = 4.0

! Initialize smoothed air temperatures for estimating headwaters temperatures
T_smth = 4.0
T_smth_local = 4.0

! Create the netCDF file and define the dimensions
open(20,file = TRIM(temp_file), status = 'unknown')
write(out_file_ncdf, *) TRIM(temp_file), '.nc'
write(out_file_ncdf, *) TRIM( adjustl(out_file_ncdf) )

call CREATE

! My allocation to save in a netcdf file
allocate( Qs_vec(heat_cells, nsegs) )
allocate( Ts_vec(heat_cells, nsegs) )
allocate( Tas_vec(heat_cells, nsegs) )
allocate( Tsoil_vec(heat_cells, nsegs) )
allocate( Width_vec(heat_cells, nsegs) )
allocate( Depth_vec(heat_cells, nsegs) )
allocate( Qrun_vec(heat_cells, nsegs) )
allocate( Qbas_vec(heat_cells, nsegs) )
allocate( Qrunsnow_vec(heat_cells, nsegs) )

allocate( Qs_2d(nlon, nlat, nsegs) )
allocate( Ts_2d(nlon, nlat, nsegs) )

Qs_2d(:,:,:) = 1.e+20
Ts_2d(:,:,:) = 1.e+20

! Initialize dummy counters that facilitate updating simulated values
n1    = 1
n2    = 2
nobs  = 0
ndays = 0
csegs = 1
xwpd  = nwpd
hpd   = 1./xwpd

! Year loop starts
do nyear = start_year, end_year
   write(*,*) ' Simulation Year - ',nyear,start_year,end_year
   nd_year = 365
   if (Is_Leap_Year(nyear)) nd_year = 366
   
   ! Day loop starts
   do nd = 1, nd_year
      
      year    = nyear
      xd_year = nd_year
      
      ! Start the numbers of days-to-date counter
      ndays   = ndays+1
      
      ! Daily period loop starts
      do ndd = 1, nwpd
         
         status = Day_of_Year_to_Date(nd, year, Day_of_Month = day, Month = month)
         
         ! Read the hydrologic and meteorologic forcings
         call READ_FORCING
         
         !!! Begin reach computations
         ! Begin cycling through the reaches
         do nr = 1, nreach
!             if(time.eq.1989.0042) write(*,*) 'Begin of #reach in a stream loop'

            nc_head = segment_cell(nr,1)
            
            !!! Mohseni relation for T_head
            ! Determine smoothing parameters (UW_JRY_2011/06/21)
            rminsmooth = 1.0 - smooth_param(nr)
            T_smth(nr) = rminsmooth*T_smth(nr) + smooth_param(nr)*dbt(nc_head)
            ! Variable Mohseni parameters (UW_JRY_2011/06/16)
            T_mohs = mu(nr) + (alphaMu(nr) &
                   /(1.+exp(gmma(nr)*(beta(nr)-T_smth(nr))))) 

            !!! Simple pseudo-linear relation to calculate T_head
            Tseuil = 12.
            Thh = max(0., tavg(nc_head) + 1.2*( dbt(nc_head) - tavg(nc_head)) )
            if (Thh .gt. Tseuil) then
               Thh = min( Thh, Tseuil + 1.2*(dbt(nc_head) - Tseuil) )
            end if

            !!! Set the baseflow temperature (one of four choices)
!!! T_base based on Mohseni et al. regression
!            T_base = T_mohs
!!! T_base based on Pseudo-linear relation
!            T_base = Thh
!!! T_base equal to smoothed air temperature
!            T_base = T_smth(nr)
!!! T_base equal to soil temperature
            T_base = tsoil(nc_head)

            !!! Set headwater temperature, where local baseflow is at T_base and mix with
            ! runoff (at T_runoff) and runoff from snow and glacier melt (at T_snow)
            if (T_base .lt. 0.) T_base = 0.
            T_runoff = dbt(nc_head)
            if(T_runoff .le. 0.) T_runoff = 0.
            T_snow = 0.
            Q_local = Q_bas(nc_head) + Q_run(nc_head)
            T_head(nr) = Q_bas(nc_head)*T_base                          / Q_local &
                       + (Q_run(nc_head) - Q_runsnow(nc_head))*T_runoff / Q_local &
                       + Q_runsnow(nc_head)*T_snow                      / Q_local

            temp(nr,0,n1) = T_head(nr)
            temp(nr,1,n1) = T_head(nr)

            ! Begin cell computational loop
            do ns=1,no_celm(nr)
               
               DONE = .FALSE.
               
               ! Testing new code 8/8/2016
               !Establish particle tracks
               call Particle_Track(nr,ns,nx_s,nx_head)
               
               ncell = segment_cell(nr,ns)
               
               ! Now do the third-order interpolation to
               ! establish the starting temperature values
               ! for each parcel
               nseg = nstrt_elm(ns)
               
               !!! Perform polynomial interpolation
               ! Interpolation inside the domain
               npndx = 2
               
               ! Interpolation at the upstream boundary if the
               ! parcel has reached that boundary
               if(nx_head.eq.0) then
                  T_0 = T_head(nr)
               else 
                  ! Interpolation at the upstream or downstream boundary
                  if(nseg .eq. 1 .or. nseg .eq. no_celm(nr)) npndx=1
                  
                  do ntrp = nterp(npndx),1,-1
                     npart    = nseg + ntrp + ndltp(npndx)
                     xa(ntrp) = x_dist(nr,npart)
                     ta(ntrp) = temp(nr,npart,n1)
                  end do
                  
                  ! Start the cell counter for nx_s
                  x = x_part(nx_s)
                  
                  ! Call the interpolation function
                  T_0 = tntrp(xa,ta,x,nterp(npndx))
               end if ! End of if nx_head

               nncell = segment_cell(nr,nstrt_elm(ns))

               ! Initialize inflow, outflow, runoff and baseflow (usefull ?)
               !Q_inflow  = Q_in(nncell)
               !Q_outflow = Q_out(nncell)
!               Q_runoff = Q_run(nncell)
!               Q_base   = Q_bas(nncell)
!               Q_runoffsnow = Q_runsnow(nncell)

               ! Set NCELL0 for purposes of tributary input
               ncell0 = nncell
               dt_total = 0.0

               do nm=no_dt(ns),1,-1
                  dt_calc = dt_part(nm)
                  x_calc  = x_part(nm)
                  
                  z = depth(nncell)
                  w = width(nncell)
                  !if (z .eq. 0.) z = 0.01
                  !if (w .eq. 0.) w = 0.01
                  
                  ! Compute the net heat flux
                  call energy(T_0,q_surf,nncell)

                  T_0 = T_0 + q_surf*dt_calc / (rho*Cp*z)
                  if(T_0 .le. 0.) T_0 = 0.
                  
!!! Original version of RBM
                  !q_dot=(q_surf/(z*rfac))
                  !T_0=T_0+q_dot*dt_calc
                  !if(T_0 .lt. 0.0) T_0=0.0

                  !!!! Read inflow and outflow (assumed equal for the moment)
                  Q_inflow      = Q_out(nncell)
                  Q_outflow     = Q_out(nncell)
                  
                  !!!! Add distributed flows from baseflow and runoff
                  Q_base        = Q_bas(nncell)
                  Q_runoff      = Q_run(nncell)
                  Q_runoffsnow  = Q_runsnow(nncell)
                  Q_local = Q_base + Q_runoff
                  
                  !!! Calculate temperature of distributed inflow from mixing ratios
                  ! of local baseflow, runoff, and runoff from snow and glacier melt.
                  ! Note that these are taken directly form VIC output and have not
                  ! been subject to grid cell routing. Hence, their total, Q_local,
                  ! is not necessarily equal to Q_outflow - Q_inflow
 
                  !!! Mohseni relation for T_base
                  rminsmooth = 1.0 - smooth_param(nr)
                  T_smth_local(nncell) = rminsmooth*T_smth_local(nncell) &
                                       + smooth_param(nr)*dbt(nncell)
                  T_mohs = mu(nr) + (alphaMu(nr) &
                         / (1.+exp(gmma(nr)*(beta(nr)-T_smth_local(nncell))))) 
                  
                  !!! Simple pseudo-linear relation to calculate T_base
                  Tseuil = 12.
                  Thh = max(0., tavg(nncell) + 1.2*( dbt(nncell) - tavg(nncell)) )
                  if (Thh .gt. Tseuil) then
                     Thh = min( Thh, Tseuil + 1.2*(dbt(nncell) - Tseuil) )
                  end if
                  
                  !!! Set baseflow temperature (one of four choices)
!!! T_base from Mohseni et al. regression
!                  T_base = T_mohs
!!! T_base from Pseudo-linear relation
!                  T_base = Thh
!!! T_base equal to smoothed air temperature
!                  T_base = T_smth_local(nncell)
!!! T_head equal to soil temperature
                  T_base = tsoil(nncell) 

                  if(T_base .le. 0.) T_base = 0.
                  T_runoff      = dbt(nncell)
                  if(T_runoff .le. 0.) T_runoff = 0.
                  T_snow = 0.
                  T_local = Q_base * T_base                      / Q_local &
                          + (Q_runoff - Q_runoffsnow) * T_runoff / Q_local &
                          + Q_runoffsnow * T_snow                / Q_local
                  if(Q_local .gt. 0.) then
                     local_z   = 0.2307 * Q_local**0.4123
                     T_local    = T_local + q_surf*dt_calc / (rho*Cp*local_z)
                  end if
                  T_local_load = Q_local * T_local

                  !!!! Add distributed flows.
                  Q_dstrb       = Q_diff(nncell)
                  !T_dstrb_load  = Q_dstrb*T_dstrb

                  !!!! Look for a tributary.
                  ntribs = no_tribs(nncell)
                  Q_trb_sum   = 0.0
                  T_trb_load  = 0.0

                  if(ntribs.gt.0 .and. .not.DONE) then
                     do ntrb = 1,ntribs
                        nr_trib = trib(nncell,ntrb)
                        
                        if(Q_trib(nr_trib).gt.0.0) then
                           Q_trb        = Q_trib(nr_trib)
                           Q_trb_sum    = Q_trb_sum + Q_trb

                           ! Update water temperature with tributary input
                           if( .not.isnan(T_trib(nr_trib)) ) then
                              
                              T_tribs      = T_trib(nr_trib) + q_surf*dt_calc / (rho*Cp*z)
                              if(T_tribs .le. 0.) T_tribs = 0.

                              T_trb_load   = (Q_trb*T_tribs)  &
                                           +  T_trb_load
                           end if
                        end if ! End of if on existence of tributary

                     end do ! End of loop on number of tributaries
                     DONE = .TRUE.

                     ! Update Q_outflow
                     Q_inflow  = Q_outflow - Q_dstrb - Q_trb_sum
                  end if ! End of if on ntribs

                  !  Ratio (only when there are tributaries)
                  Q_ratio   = Q_inflow/Q_outflow
                  
                  !!! Original version of RBM
                  ! Do the mass/energy balance
!                  T_0  = T_0*Q_ratio                            &
!                       + (T_dstrb_load + T_trb_load)/Q_outflow  &
!                       + q_dot*dt_calc
                  
                  !!! Updated version
                  T_0  = T_0*(Q_inflow - Q_runoff - Q_base)  / Q_outflow &
                       + T_trb_load                          / Q_outflow &
                       + T_local_load                        / Q_outflow
                  
                  !T_0 minimum value is 0.5 (double check ?!?)
                  if (T_0 .lt. 0.) T_0 = 0.

                  ! Q_inflow becomes Q_outflow (for next time step ?)
                  Q_inflow = Q_outflow
                  
                  nseg   = nseg + 1
                  nncell = segment_cell(nr,nseg)

                  ! Reset tributary flag if this is a new cell
                  if(ncell0 .ne. nncell) then
                     ncell0   = nncell
                     Q_inflow = Q_in(nncell)
                     DONE     = .FALSE.
                  end if

                  ! ... ?
                  dt_total = dt_total + dt_calc
               end do ! End of loop on nm = no_dt(ns),1,-1 (?!?)

               temp(nr,ns,n2) = T_0
               T_trib(nr)     = T_0
               
               ! Write netcdf data in 1D variables
               Qs_vec(ncell, csegs)    = Q_outflow
               Ts_vec(ncell, csegs)    = T_0
               Tas_vec(ncell, csegs)   = dbt(ncell)
               Tsoil_vec(ncell, csegs) = tsoil(nncell)
               Width_vec(ncell, csegs) = width(ncell)
               Depth_vec(ncell, csegs) = depth(ncell)
               Qrun_vec(ncell, csegs)  = Q_runoff
               Qbas_vec(ncell, csegs)  = Q_base
               Qrunsnow_vec(ncell, csegs) = Q_runoffsnow
               
               csegs = csegs + 1
               if (csegs.gt.nsegs) csegs = 1
               
               !   Write all temperature output UW_JRY_11/08/2013
               !   The temperature is output at the beginning of the 
               !   reach.  It is, of course, possible to get output at
               !   other points by some additional code that keys on the
               !   value of ndelta (now a vector)(UW_JRY_11/08/2013)
               !call WRITE(year,month,day,nd,nr,ncell,nncell,T_0,T_base,dbt(ncell),Q_inflow,Q_outflow,Q_runoff,Q_base)

            end do ! End of computational element loop (ns=1,no_celm(nr))
         end do ! End of reach loop
         
         ntmp = n1
         n1 = n2
         n2 = ntmp
         
         4650 format(16x,12(6x,f6.0,6x))
         4700 format(f10.4,f6.0,15(f6.1,f8.3))
         4750 format(f10.4,10(i4,f8.0))
      end do ! End of weather period loop (NDD=1,NWPD)
      
      start(ndims) = ndays
      call check( nf90_put_var(ncid, Qsid,    Qs_vec,    start = start, count = count) )
      call check( nf90_put_var(ncid, Tsid,    Ts_vec,    start = start, count = count) )
      call check( nf90_put_var(ncid, Tasid,   Tas_vec,   start = start, count = count) )
      call check( nf90_put_var(ncid, Tsoilid, Tsoil_vec, start = start, count = count) )
      call check( nf90_put_var(ncid, Widthid, Width_vec, start = start, count = count) )
      call check( nf90_put_var(ncid, Depthid, Depth_vec, start = start, count = count) )
      call check( nf90_put_var(ncid, Qrunid,  Qrun_vec,  start = start, count = count) )
      call check( nf90_put_var(ncid, Qbasid,  Qbas_vec,  start = start, count = count) )
      call check( nf90_put_var(ncid, Qrunsnowid,  Qrunsnow_vec,  start = start, count = count) )
      
   end do ! End of day loop (ND=1,365/366)
end do ! End of year loop

!!! Create the time counter and write it:
allocate( time(ndays) )
time(:) = (/(i, i=0,ndays-1, 1)/)
call check( nf90_put_var(ncid, tid, time) )

!!! Close the netcdf file
write(*,*) 'ndays :', ndays
call check( nf90_close(ncid) )
!     ******************************************************
!                        return to rmain
!     ******************************************************
!
950 return
end SUBROUTINE SYSTMM

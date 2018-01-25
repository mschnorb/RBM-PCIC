SUBROUTINE SYSTMM(temp_file,param_file)

use Block_Energy
use Block_Hydro
use Block_Network

Implicit None

character (len=200):: temp_file
character (len=200):: param_file

integer          :: ncell,nncell,ncell0,nc_head,no_flow,no_heat
integer          :: nc,nd,ndd,nm,nr,ns
integer          :: nr_trib,ntribs
integer          :: nrec_flow,nrec_heat
integer          :: n1,n2,nnd,nobs,nyear,nd_year,ntmp
integer          :: npart,nseg,nx_s,nx_part,nx_head

! Indices for lagrangian interpolation
integer              :: njb,npndx,ntrp
integer, dimension(2):: ndltp=(/-1,-2/)
integer, dimension(2):: nterp=(/2,3/)

real             :: dt_calc,dt_total,hpd,q_dot,q_surf,z
real             :: Q_dstrb,Q_inflow,Q_outflow,Q_ratio,Q_trb,Q_trb_sum,Q_base,Q_runoff,myratio
real             :: T_dstrb,T_dstrb_load,T_trb_load,T_base_load, T_runoff_load
real             :: rminsmooth
real             :: T_0,T_dist,Thh,Tseuil,T_base,T_runoff
real(8)          :: time
real             :: x,xd,xdd,xd_year,xwpd,year
real             :: tntrp
real             :: dt_ttotal
real,dimension(4):: ta,xa

real,dimension(:),allocatable     :: T_head,T_smth,T_trib

logical:: DONE

! Allocate the arrays
allocate (temp(nreach,0:ns_max,2))
allocate (T_head(nreach))
allocate (T_smth(nreach))
allocate (T_trib(nreach))
allocate (depth(heat_cells))
allocate (Q_in(heat_cells))
allocate (Q_out(heat_cells))
allocate (Q_run(heat_cells))
allocate (Q_bas(heat_cells))
allocate (Q_diff(heat_cells))
allocate (Ratio(heat_cells))
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

! open the output file
open(20,file=TRIM(temp_file),status='unknown')

! Initialize dummy counters that facilitate updating simulated values
n1    = 1
n2    = 2
nobs  = 0
ndays = 0
xwpd  = nwpd
hpd   = 1./xwpd

! Year loop starts
do nyear=start_year,end_year
   write(*,*) ' Simulation Year - ',nyear,start_year,end_year
   nd_year = 365
   if (mod(nyear,4).eq.0) nd_year = 366 ! Not true for 1900. To correct.

   ! Day loop starts
   do nd=1,nd_year
      year=nyear
      xd=nd
      xd_year=nd_year
      
      ! Start the numbers of days-to-date counter
      ndays=ndays+1

      ! Daily period loop starts
      do ndd=1,nwpd
         xdd  = ndd
         time = year + (xd + (xdd-0.5)*hpd) / xd_year 


         ! Read the hydrologic and meteorologic forcings
         call READ_FORCING

         !!! Begin reach computations
         ! Begin cycling through the reaches
         do nr=1,nreach
            nc_head = segment_cell(nr,1)
            
            !!! Mohseni relation for T_head
            ! Determine smoothing parameters (UW_JRY_2011/06/21)
            rminsmooth = 1.0 - smooth_param(nr)
            T_smth(nr) = rminsmooth*T_smth(nr) + smooth_param(nr)*dbt(nc_head)

            ! Variable Mohseni parameters (UW_JRY_2011/06/16)
            T_head(nr) = mu(nr) + (alphaMu(nr) &
                        /(1.+exp(gmma(nr)*(beta(nr)-T_smth(nr))))) 
            
            !!! Other method to calculate T_head
            ! Simple pseudo-linear relation
            Tseuil = 12.
            Thh = max(0., tavg(nc_head) + 1.2*( dbt(nc_head) - tavg(nc_head)) )
            if (Thh .gt. Tseuil) then
               Thh = min( Thh, Tseuil + 1.2*(dbt(nc_head) - Tseuil) )
            end if
!             T_head(nr) = Thh
            
            !!! T_head = T_soil
            T_head(nr) = tsoil(nc_head)

            if (T_head(nr) .lt. 0.) T_head(nr) = 0.
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
                     npart    = nseg+ntrp+ndltp(npndx)
                     xa(ntrp) = x_dist(nr,npart)
                     ta(ntrp) = temp(nr,npart,n1)
                  end do

                  ! Start the cell counter for nx_s
                  x = x_part(nx_s)

                  ! Call the interpolation function
                  T_0=tntrp(xa,ta,x,nterp(npndx))

               end if ! End of if nx_head

               nncell = segment_cell(nr,nstrt_elm(ns))

               ! Initialize inflow
               Q_inflow = Q_in(nncell)
               Q_outflow = Q_out(nncell)
               myratio = Ratio(nncell)

               ! Set NCELL0 for purposes of tributary input
               ncell0 = nncell
               dt_total = 0.0

               do nm=no_dt(ns),1,-1
                  dt_calc = dt_part(nm)
                  z = depth(nncell)
                  
                  ! Compute some energy variables
                  call energy(T_0,q_surf,nncell)

                  q_dot = (q_surf/(z*rfac)) !* 0.5
                  T_0   = T_0 + q_dot*dt_calc
                  
                  ! T_0 minimum value is 0.0 (ice)
                   if(T_0.lt.0.0) T_0=0.0

                  !!!! Add distributed flows
                  T_dstrb_load  = 0.0
                  Q_dstrb       = Q_diff(nncell)
                  
                  !!!! Add baseflow
                  T_base        = tsoil(nncell)
                  if(T_base .le. 0.) T_base = 0.
                  Q_base        = Q_bas(nncell)
!                  T_base_load   = Q_base*(1-myratio)*T_base + Q_base*myratio*0.0 - Q_base*T_0
                  T_base_load   = Q_base*T_base - Q_base*T_0
                  
                  !!!! Add runflow
                  T_runoff      = dbt(nncell)
                  if(T_runoff .le. 0.) T_runoff = 0.
                  Q_runoff      = Q_run(nncell)
!                  T_runoff_load = Q_runoff*(1-myratio)*T_runoff + Q_base*myratio*0.0 - Q_runoff*T_0
                  T_runoff_load = Q_runoff*T_runoff - Q_runoff*T_0
                  
                  ! Temperature of distributed inflow assumed = 10.0 deg C
                  if(Q_dstrb.gt.0.001) then
                     T_dstrb  = 10.0
                  else
                     T_dstrb  = 10.0
                  end if
                  T_dstrb_load  = Q_dstrb*T_dstrb

                  !!! Add Baseflow and Runoff
                  T_dstrb_load  = T_dstrb_load + T_base_load  + T_runoff_load

                  !!! Look for a tributary.
                  ntribs = no_tribs(nncell)
                  Q_trb_sum   = 0.0
                  T_trb_load  = 0.0
                  
                  if(ntribs.gt.0 .and. .not.DONE) then
                     do ntrb=1,ntribs
                        nr_trib=trib(nncell,ntrb)
                        
                        if(Q_trib(nr_trib).gt.0.0) then
                           Q_trb        = Q_trib(nr_trib)
                           Q_trb_sum    = Q_trb_sum + Q_trb

                           ! Update water temperature with tributary input
                           T_trb_load   = (Q_trb*T_trib(nr_trib))       &
                                          +  T_trb_load
                        end if ! End of if on existence of tributary

                     end do ! End of loop on number of tributaries

                     DONE = .TRUE.

                  end if ! End of if on ntribs

                  !  Update inflow and outflow
                  Q_outflow = Q_inflow + Q_dstrb + Q_trb_sum
                  Q_ratio   = Q_inflow/Q_outflow

                  ! Do the mass/energy balance
                  T_0  = T_0*Q_ratio                            &
                       + (T_dstrb_load + T_trb_load)/Q_outflow  &
                       + q_dot*dt_calc

                  ! T_0 minimum value is 0.5 (double check ?!?)
                   if (T_0.lt.0.) T_0 = 0.

                  ! Q_inflow becomes Q_outflow (for next time step ?)
                  Q_inflow = Q_outflow

                  nseg   = nseg+1
                  nncell = segment_cell(nr,nseg)

                  ! Reset tributary flag if this is a new cell
                  if(ncell0.ne.nncell) then
                     ncell0 = nncell
                     Q_inflow = Q_in(nncell)
                     DONE = .FALSE.
                  end if

                  ! ... ?
                  dt_total = dt_total + dt_calc

               end do ! End of loop on nm=no_dt(ns),1,-1 (?!?)

               ! New check on minimum T_0 value... 0.5 again
                if (T_0.lt.0.) T_0 = 0.

               temp(nr,ns,n2) = T_0
               T_trib(nr) = T_0

               !   Write all temperature output UW_JRY_11/08/2013
               !   The temperature is output at the beginning of the 
               !   reach.  It is, of course, possible to get output at
               !   other points by some additional code that keys on the
               !   value of ndelta (now a vector)(UW_JRY_11/08/2013)
               !call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell),Q_inflow,Q_outflow)
               call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell),Q_runoff,Q_base)

            end do ! End of computational element loop (ns=1,no_celm(nr)) (segments ?)
         end do ! End of reach loop
      
         ntmp = n1
         n1 = n2
         n2 = ntmp

         4650 format(16x,12(6x,f6.0,6x))
         4700 format(f10.4,f6.0,15(f6.1,f8.3))
         4750 format(f10.4,10(i4,f8.0))
      
      end do ! End of weather period loop (NDD=1,NWPD)
   end do ! End of day loop (ND=1,365/366)
end do ! End of year loop
!
!
!     ******************************************************
!                        return to rmain
!     ******************************************************
!
950 return
end SUBROUTINE SYSTMM

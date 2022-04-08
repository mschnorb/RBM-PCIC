SUBROUTINE Read_Forcing
!
USE Block_Energy
USE Block_Hydro
USE Block_Network
! 
IMPLICIT NONE
!
integer :: nc,ncell,nnd,no_flow,no_heat,nr,nrec_flow,nrec_heat
real    :: Q_avg,Q_dmmy,rho_dmmy
logical, save :: NotReaded_RiverParam = .TRUE. !Shahab (2022)

no_flow = 0
no_heat = 0
if( NotReaded_RiverParam ) then !Shahab (2022)
    NotReaded_RiverParam = .FALSE. !Shahab (2022)
    read(300,"(4(f7.4))",end=100) Awidth,Bwidth,Adepth,Bdepth ! 
    write(*,*) 'Awidth  :', Awidth
    write(*,*) 'Bwidth  :', Bwidth
    write(*,*) 'Adepth  :', Adepth
    write(*,*) 'Bdepth  :', Bdepth
    100 close(300)
endif !Shahab (2022)

do nr=1,nreach ! Loop over reach in the network
  do nc=1,no_cells(nr)-1 ! Loop over cells in the reach
    no_flow = no_flow + 1
    no_heat = no_heat + 1

    nrec_flow = flow_cells*(ndays-1) + no_flow
    nrec_heat = heat_cells*(ndays-1) + no_heat

!     read(35,'(2i8,4f12.5,3f8.3)' &
!            ,rec=nrec_flow) nnd,ncell &
!            ,Q_out(no_heat),Q_run(no_heat),Q_bas(no_heat),Q_runsnow(no_heat) &  
!            ,depth(no_heat),width(no_heat),u(no_heat)
           
!    read(35,'(2i8,5f12.5)' &
!           ,rec=nrec_flow) nnd,ncell &
!           ,Q_out(no_heat),Q_run(no_heat),Q_bas(no_heat),Q_runsnow(no_heat) &
!           ,u(no_heat)
     read(35,'(2i8,5f12.5)' &
           ,rec=nrec_flow) nnd,ncell &
           ,Q_out(no_heat),Q_run(no_heat),Q_bas(no_heat),Q_runsnow(no_heat) &
           ,u(no_heat)  !S.Larabi (2019)
   
!    depth(no_heat) = 0.2307 * Q_out(no_heat)**0.4123
!    width(no_heat) = 6.6588 * Q_out(no_heat)**0.4967
    depth(no_heat) = Adepth * Q_out(no_heat)**Bdepth !S.Larabi (2019)
    width(no_heat) = Awidth * Q_out(no_heat)**Bwidth !S.Larabi (2019)
    !u(no_heat)     = 1.5
    u(no_heat)     = Q_out(no_heat) / (depth(no_heat) * width(no_heat))

    Q_diff(no_heat) = 0. ! We will deal with that later. If artificial lateral inflow.

 !   if(u(no_heat).lt.0.01) u(no_heat)=0.01
    if(u(no_heat).lt.0.45) u(no_heat)=0.45 !S.Larabi (2019) low min speed from YEARSLEY (2012)
    if(ncell.ne.no_heat) write(*,*) 'Flow file error',ncell,no_heat 

    read(36,'(i5,4f6.1,2f10.4,f6.3,f7.1,f5.1)' &
           ,rec=nrec_heat) ncell &
           ,dbt(no_heat),tavg(no_heat),tsoil(no_heat) &
           ,ea(no_heat),Q_ns(no_heat),Q_na(no_heat) &
           ,rho_dmmy,press(no_heat),wind(no_heat)

  if(ncell.ne.no_heat) write(*,*) 'Heat file error',ncell,no_heat
!
!  Added variable ndelta (UW_JRY_2011/03/15
!
    delta_n=ndelta(ncell)
! 
    Q_avg=0.5*(Q_in(no_heat)+Q_out(no_heat))
    dt(no_heat)=dx(no_heat)/u(no_heat)

!
!  Added check to see if travel time of parcel exceeds the
!  computational interval.  If so, it writes to file fort.45.
!  One should check to see if there are NaN's at that node.
!  (UW_JRY_2011/03/15)
!
    if(dt(no_heat).gt.dt_comp) write(45,*) &
           'Travel time=',dt(no_heat) &
            , '> dt_comp at node -',no_heat
  end do ! End of loop on cells (in a reach)
!
! Tributary flow is Q_out from the next to the last cell
! However, it will be updated in Water_Balance to account
! for one-half the runoff from the confluence cell.
!
!
!       Read the meteorology for the last cell, but not the flow
!
  no_heat = no_heat + 1

  Q_out(no_heat) = Q_out(no_heat-1)

  Q_run(no_heat)     = Q_run(no_heat-1)
  Q_bas(no_heat)     = Q_bas(no_heat-1)
  Q_runsnow(no_heat) = Q_runsnow(no_heat-1)

  Q_trib(nr) = Q_out(no_heat)
  nrec_heat  = heat_cells*(ndays-1) + no_heat

  read(36,'(i5,3f6.1,2f7.4,f6.3,f7.1,f5.1)' &
         ,rec=nrec_heat) ncell &
         ,dbt(no_heat),tsoil(no_heat),ea(no_heat) &   
         ,Q_ns(no_heat),Q_na(no_heat),rho_dmmy &
         ,press(no_heat),wind(no_heat)
!
!  The flow and hydraulics for the last cell has to be 
!  modified so they do not
!  take the values of the segment to which it is tributary
!

  Q_in(no_heat)   = Q_out(no_heat-1)
  u(no_heat)      = u(no_heat-1)
  depth(no_heat)  = depth(no_heat-1)
  width(no_heat)  = width(no_heat-1)
  dt(no_heat)     = dx(ncell) / u(no_heat)
end do ! End of loop on reach
! Call the water balance subroutine
!
call Water_Balance
!
END SUBROUTINE Read_Forcing

SUBROUTINE WRITE(year,month,day,nd,nr,ncell,ns,T_0,T_head,dbt,Q_inflow,Q_outflow,Q_runoff,Q_base)

use netcdf

Implicit NONE
!
integer :: year, month, day
integer :: nd,nr,ncell,ns 
real    :: Q_inflow,Q_outflow,Q_runoff,Q_base,Diff
real    :: T_0,T_dist
! real(8) :: time
real    :: T_head
real    :: dbt
!

Diff = Q_outflow - Q_runoff - Q_base

!!write (20,'(7i4,3f10.2,5f15.5,f8.5)') &
!!    year,month,day,nd,nr,ncell,ns,T_0,T_head,dbt,Q_inflow,Q_outflow,Q_runoff,Q_base,Diff
write (20,'(4i8,3f10.2,5f15.5)') &
    year,month,day,ncell,T_0,T_head,dbt,Q_outflow,Q_runoff,Q_base
end SUBROUTINE WRITE

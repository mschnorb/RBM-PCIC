SUBROUTINE WRITE(time,nd,nr,ncell,ns,T_0,T_head,dbt,Q_inflow,Q_outflow,Q_runoff,Q_base)
!
Implicit NONE
!
integer :: nd,nr,ncell,ns 
real    :: Q_inflow,Q_outflow,Q_runoff,Q_base,Diff
real    :: T_0,T_dist
real(8) :: time
real    :: T_head
real    :: dbt
!

Diff = Q_outflow - Q_runoff - Q_base

write (20,'(f12.4,4i4,3f10.2,5f15.5,f8.5)') &
    time,nd,nr,ncell,ns,T_0,T_head,dbt,Q_inflow,Q_outflow,Q_runoff,Q_base,Diff
end SUBROUTINE WRITE

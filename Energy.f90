SUBROUTINE Energy(T_surf,q_surf,ncell)
   use Block_Energy
   implicit none
   integer::i,ncell
   real::A,B,e0,q_surf,q_conv,q_evap,q_ws,T_surf
   real, dimension(2):: q_fit, T_fit
   real::SF,FW
   real::qs,ql,Le,Kl,esat,E,qe,qh
!-------------------------------------------------------------------------------!
! Subroutine to do the energy budget in a grid cell.
! The energy balance has been modified from Yearsley et al. (2009, 2012) to
! follow Wu et al. (2012): doi:10.1029/2012WR012082
!
! ** The formulation of the Net solar short-wave radiation considers a Shading
! Factor (SF) due to the vegetation.
! ** The turbulent exchange of water vapor also considers a Sheltering Factor (FW)
! due to the vegetation.
!     ---> Both are considered equal to 1. everywhere for the moment.

! ** The formulation of the Flux of sensible heat is slightly different and now
! use the atmospheric pressure.
!
! ** The conductive heat flux is assumed to be equal to 0. for the moment.
! 
! All variables are now in metric.
!-------------------------------------------------------------------------------!

   SF = 1. ! Shading Factor due to vegetation (Solar)
   FW = 1. ! Sheltering Factor due to vegetation (Wind)

   T_fit(1)=T_surf-1.0
   T_fit(2)=T_surf+1.0
   do i=1,2
      ! Net solar short-wave radiation (incoming - reflected)
      qs = 0.97 * SF * q_ns(ncell)

      ! Net solar long-wave radiation (incoming - reflected)
      ql = q_na(ncell) - 0.97 * Boltzmann * (T_fit(i)+273.15)**4

      ! Flux of latent heat (evaporation)
         Le   = 2499.64 - 2.51 * T_fit(i)               ! Latent heat of vaporization
         Kl   = 0.211 + 0.103 * wind(ncell) * FW        ! Turbulent exchange of water vapor
         Kl   = Kl / (86.4E6)                           ! mm.d**-1 to m.s**-1
         esat = 2.1718E8*EXP(-4157.0/(T_fit(i)+239.09)) ! Vapor pressure at saturation
         E    = Kl * (esat - ea(ncell))                 ! Evaporation rate
      qe = - rho * E * Le

      ! Flux of sensible heat
      qh = - gam * rho * Kl * Le * press(ncell)*10/1013. * (T_fit(i) - dbt(ncell))

      ! Conductive heat flux
      ! qc is not considered here. Assumed to be 5% of the net solar radiative
      ! flux in Wu et al. (2012).

      q_fit(i)= qs + ql + qe + qh
   end do

!
!     q=AT+B
!
!     Linear fit over the range of 2.0 deg C.
!     These results can be used to estimate the "equilibrium" 
!     temperature and linear rate constant.
!
   A=(q_fit(1)-q_fit(2))/(T_fit(1)-T_fit(2))
   q_surf=0.5*(q_fit(1)+q_fit(2))
   B=(q_surf/A)-(T_fit(1)+T_fit(2))/2.
!
!     ******************************************************
!               Return to Subroutine RIVMOD
!     ******************************************************
!
END Subroutine Energy

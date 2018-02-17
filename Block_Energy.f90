module Block_Energy
!
!   Energy budget variables
!
!   Incoming short wave radiation, kcal/m**2/sec
!
    real, dimension(:), allocatable::q_ns
!
!   Incoming atmospheric radiation, kcal/m**2/sec
!
    real, dimension(:), allocatable::q_na
!
!   Air temperature at surface, deg. C
!
    real, dimension(:), allocatable::dbt
!
!   Annual Air temperature at surface, deg. C
!
    real, dimension(:), allocatable::tavg
!
!   Soil temperature (bottom layeri of VIC), deg. C
!
    real, dimension(:), allocatable::tsoil
!  
!   Wind speed, m/sec
!
    real, dimension(:), allocatable::wind
!
!   Vapor pressure of air at surface, mb
!
    real, dimension(:), allocatable::ea
!
!   Air pressure at surface, mb
!
    real, dimension(:), allocatable::press 

!
    real, dimension (:), allocatable::mu,alphamu,beta,gmma,smooth_param

!   Some important constants
!
      real,parameter   :: Boltzmann = 5.67e-8  ! Stefan-Boltzmann constant (W.m**-2.K**-4)
      real,parameter   :: rho = 1000           ! Water density (kg.m**-3)
      real,parameter   :: Cp  = 4.187E+3       ! Specific heat of water J.kg**-1.degC**-1
      real,parameter   :: gam = 0.655          ! Psychrometric constant (hPa.degC**-1)
!
end module Block_Energy  

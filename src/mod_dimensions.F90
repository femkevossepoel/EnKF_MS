module mod_dimensions
#ifdef ZEROD
! nop: number of periods (1)
! nop=100
! dtro: time step length (1/60)
! dtro=1./60.
! eps: epsilon (0.0)
! eps=0.
! mu: coupling coefficient (2/3)
! mu=2./3.
! ac XXXXXXXXX (F)
! wsf XXXXXXXX (F)
! sstemp_0 (T_0: initial sea surface temperature (1.125/7.5)
! sstemp_0=1.125/7.5
! ocheight_0 (h_0): initial height (0)
! ocheight_0=0.
!
   integer, parameter :: nx=1                      ! zero-D, so number of gridpoints is 1
   integer, parameter :: nop=100
   real, parameter :: dtro=1./60.
   real, parameter ::  eps=0.
   real, parameter ::  mu=2./3.
   real, parameter ::  sstemp_0=1.125/7.5
   real, parameter ::  ocheight_0=0.
#else
   integer, parameter :: nx=1024                     ! Number of gridpoints in $x$-direction
#endif
   integer, parameter :: ndim=2*nx
end module mod_dimensions

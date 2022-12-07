module ro_constants
! CONSTANTS

! high end value of the coupling parameter
   real, parameter :: b_0=2.5

! feedback of the thermocline gradient on the SST gradient
   real, parameter :: gamma_a=0.75

! damping rate of SST anomalies
   real, parameter :: c_damp=1.0

! damping of upper heat content
   real, parameter :: r_damp=0.25

!relates enchanced easterly wind stress to the recharge of the ocean heat content
   real, parameter :: alpha=0.125

! frequency
   real, parameter :: omega_c=sqrt(3./32.)

! period (nondimensionalised)
   real, parameter :: pi=4*atan(1.0)
   real, parameter :: period = 2*pi/omega_c
end module ro_constants











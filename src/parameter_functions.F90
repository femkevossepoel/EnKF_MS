module parameter_functions
   use mod_dimensions
   use ro_constants

contains
!   function coupling(mu,time_step)
   function coupling(mu_coup)
   use mod_dimensions
   use ro_constants
   implicit none
!   real mu, time_step,coupling
   real mu_coup
   real coupling
! possibility to include annual cycle
!
!        mu = mu_0*(1 + mu_ann*cos((2*pi*time_step/tau) - (5*pi/6)))
!
!    relates stronger thermocline gradient to stronger easterly wind stress
    coupling = b_0 * mu_coup
   end function coupling

! here, I could include a subroutine to calculate xi in case I want to include an effect of wind-stress forcing
!begin subroutine xi(t, dt, wsf, f_ann = 0.02, f_ran = 0.2, eps = 0.1, mu_0 = 0.75, 
   function xi(time_step,wsf)
      use mod_dimensions
   real time_step,xi
   integer status
   logical wsf
!       mu_ann = 0.2
!       tau=(12/2)
!       tau_cor=1/(30*2))
! with wind-stress forcing
!    if wsf.eq.true.
!    endif
!
!    for now, without wind stress forcing
   if (wsf) then
!        W = -1 + (random()*2)
!        xi = f_ann*cos(2*pi*time_step/tau) + f_ran*W*(tau_cor/dt)
! need to define all constants
        print *,'Warning, this code needs to be completed'
        call exit(status)
     else
        xi = 0.0
     endif
!
        end function xi
end module

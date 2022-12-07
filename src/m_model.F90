module m_model
   use mod_dimensions
   use mod_state
   #ifdef ZEROD
   use ro_constants
   use parameter_functions
   #endif
   type(substate) u                      ! Model velocity
   type(substate) rh                     ! Decorrelation lengths of ocean and atmosphere
   real, parameter :: dx=1.0             ! Horizontal grid spacing
   real, parameter :: dt=1.0             ! Time step of atmospheric model

contains
#ifdef ZEROD

   subroutine recharge_osc(mem,old,leuler,time)
   use mod_dimensions
   use mod_state
   use ro_constants
   use parameter_functions
   implicit none
   type(state),     intent(inout):: mem
   type(state),     intent(inout):: old
   logical,         intent(in)   :: leuler
!   real T, k1a, k1o,k2a,k2o,k3a,k3o,k4a,k4o
   real k1a, k1o,k2a,k2o,k3a,k3o,k4a,k4o
   real dTdt,dhdt
!   integer N_t, k, i_count
   type(state) new
   real time

   interface dTdt
      module procedure dTdt
   end interface

   interface dhdt
      module procedure dhdt
   end interface

   interface coupling
      module procedure coupling
   end interface


   interface xi
      module procedure xi
   end interface



   if (leuler) then
       print '(a)', ' ONED does not work with euler'
    else


   k1a=0.
   k1o=0.
   k2a=0.
   k2o=0.
   k3a=0.
   k3o=0.
   k4a=0.
   k4o=0.
!    print '(a,f10.2,a,f10.2)', ' sstemp_0 = ',sstemp_0,' ocheight_0 ',ocheight_0

      ! at the moment, xi is constant, but kept time as argument and dt in imod_dimensions
      ! to allow for variable xi.
      !
        k1a=dTdt(mem%atmos,mem%ocean,time)
        k1o=dhdt(mem%atmos,mem%ocean, time)
!    print '(a,f10.2,a,f10.2,a,f10.2)', ' atmos = ',mem%atmos,' ocean ',mem%ocean,'time ',time
!!!!!!!!!!!!!!!!!!!!!!!!!
! FV line below is to store output, should now work via dumpsol
!
!    write(10,'(3f10.4)')time,mem%atmos,mem%ocean

        k2a= dTdt(mem%atmos+k1a* dtro/2,mem%ocean + k1a*dtro/2 , time + dtro/2)
        k2o= dhdt(mem%atmos+k1o* dtro/2,mem%ocean + k1o*dtro/2 , time + dtro/2)

        k3a= dTdt(mem%atmos+k2a* dtro/2,mem%ocean + k2a*dtro/2 , time + dtro/2)
        k3o= dhdt(mem%atmos+k2o* dtro/2,mem%ocean + k2o*dtro/2 , time + dtro/2)

        k4a= dTdt(mem%atmos+k3a* dtro/2,mem%ocean + k3a*dtro/2 , time + dtro)
        k4o= dhdt(mem%atmos+k3o* dtro/2,mem%ocean + k3o*dtro/2 , time + dtro)
        !
        new%atmos=mem%atmos+dtro*1/6*(k1a+2*(k2a*k3a*k4a))
        new%ocean=mem%ocean+dtro*1/6*(k1o+2*(k2o*k3o*k4o))

   old=mem
   mem=new
   endif

end subroutine


  function dTdt(seast, h, time_step)
!   make sure to calculate xi,b=coupling
  use mod_dimensions
  use parameter_functions
  use ro_constants
  implicit none
  real dTdt
  real R,seast,h,time_step
  logical wsf

!   R = gamma_a * coupling(mu,time_step)
   R = gamma_a * coupling(mu) - c_damp
   wsf=.false.
!   dTdt = R * seast + gamma_a * h - eps * (h + coupling(mu,time_step) * seast)**3 + gamma_a * xi
   dTdt = R * seast + gamma_a * h - eps * (h + coupling(mu) * seast)**3 + gamma_a * xi(time_step,wsf)

   end function dTdt

   function dhdt(seast, h, time_step)
!   function dhdt(seast, h, mu, dt)
   use mod_dimensions
   use parameter_functions
   use ro_constants
   implicit none
   real dhdt
   real seast, h, time_step
!   real seast, h
   logical wsf

   wsf=.false.
!    dhdt = - r_damp * h - alpha * coupling(mu,time_step) * seast - alpha * xi
    dhdt = - r_damp * h - alpha * coupling(mu) * seast - alpha * xi(time_step,wsf)

    end function dhdt

#else

   subroutine model(mem,old,leuler)
   use mod_dimensions
   use mod_state
   implicit none
   type(state),     intent(inout):: mem
   type(state),     intent(inout):: old
   logical,         intent(in)   :: leuler
   integer i,ia,ib
   type(state) new
   real, parameter :: o2a=0.01
   real, parameter :: a2o=0.00
   real, parameter :: b1=0.003
   real, parameter :: b2=0.001
   if (leuler) then
      do i=1,nx
         ia=mod(i-2+nx,nx)+1
         ib=mod(i,nx)+1
         !new%atmos(i) = mem%atmos(i) - dt*u%atmos*(mem%atmos(ib)-mem%atmos(ia))/(2.0*dx) &
         new%atmos(i) = mem%atmos(ia) &
                      + dt*o2a*mem%ocean(i) - dt*b1*mem%atmos(i)
         new%ocean(i) = mem%ocean(i) - dt*u%ocean*(mem%ocean(ib)-mem%ocean(ia))/(2.0*dx) &
                      + dt*a2o*mem%atmos(i) - dt*b2*mem%ocean(i)
      enddo
   else
      do i=1,nx
         ia=mod(i-2+nx,nx)+1
         ib=mod(i,nx)+1
         !new%atmos(i) =  old%atmos(i) - dt*u%atmos*(mem%atmos(ib)-mem%atmos(ia))/dx &
         new%atmos(i) =  mem%atmos(ia) &
                      +2.0*dt*o2a*mem%ocean(i) - 2.0*dt*b1*mem%atmos(i)
         new%ocean(i) = old%ocean(i) - dt*u%ocean*(mem%ocean(ib)-mem%ocean(ia))/dx &
                      +2.0*dt*a2o*mem%atmos(i) - 2.0*dt*b2*mem%ocean(i)
      enddo
   endif
   old=mem
   mem=new

end subroutine

#endif
end module

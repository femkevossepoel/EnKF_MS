module m_statistics
use mod_dimensions
use mod_state
real, allocatable :: fullave(:,:,:)
real, allocatable :: fullvar(:,:,:)
real, allocatable :: fullcor(:,:,:)
real, allocatable :: fullcoroa(:,:,:)
contains
subroutine statistics(full,nrt,nrens,nrf)
   integer, intent(in) :: nrt
   integer, intent(in) :: nrens
   integer, intent(in) :: nrf
   type(state), intent(in) :: full(0:nrt,nrens)
#ifdef ZEROD
   integer i,j,k,kc
#else
   integer i,j,k,ic,kc
#endif

   allocate(fullave(nx,0:nrt,2))
   allocate(fullvar(nx,0:nrt,2))
   allocate(fullcor(nx,0:nrt,2))
   allocate(fullcoroa(nx,0:nrt,2))


! mean
   do j=1,nrens
   do k=0,nrt
   do i=1,nx
      #ifdef ZEROD
      fullave(i,k,1)=fullave(i,k,1)+full(k,j)%ocean
      fullave(i,k,2)=fullave(i,k,2)+full(k,j)%atmos
      #else
      fullave(i,k,1)=fullave(i,k,1)+full(k,j)%ocean(i)
      fullave(i,k,2)=fullave(i,k,2)+full(k,j)%atmos(i)
      #endif
   enddo
   enddo
   enddo
   fullave=fullave/real(nrens)

! variance
   do j=1,nrens
   do k=0,nrt
   do i=1,nx
      #ifdef ZEROD
      fullvar(i,k,1)=fullvar(i,k,1)+(full(k,j)%ocean-fullave(i,k,1))**2
      fullvar(i,k,2)=fullvar(i,k,2)+(full(k,j)%atmos-fullave(i,k,2))**2
      #else
      fullvar(i,k,1)=fullvar(i,k,1)+(full(k,j)%ocean(i)-fullave(i,k,1))**2
      fullvar(i,k,2)=fullvar(i,k,2)+(full(k,j)%atmos(i)-fullave(i,k,2))**2
      #endif
   enddo
   enddo
   enddo
   fullvar=sqrt(fullvar/real(nrens-1))

! corariance oo and aa
#ifdef ZEROD
   kc=nrt/2
   do j=1,nrens
   do k=0,nrt
   do i=1,nx
      fullcor(i,k,1)=fullcor(i,k,1) + (full(k,j)%ocean-fullave(i,k,1)) * (full(kc,j)%ocean-fullave(nx,kc,1))
      fullcor(i,k,2)=fullcor(i,k,2) + (full(k,j)%atmos-fullave(i,k,2)) * (full(kc,j)%atmos-fullave(nx,kc,2))
   enddo
   enddo
   enddo
   fullcor=fullcor/real(nrens-1)
   do k=0,nrt
   do i=1,nx
      fullcor(i,k,1)=fullcor(i,k,1)/(fullvar(i,k,1)*fullvar(nx,kc,1))
      fullcor(i,k,2)=fullcor(i,k,2)/(fullvar(i,k,2)*fullvar(nx,kc,2))
   enddo
   enddo
#else
   ic=nx/2
   kc=nrt/2
   do j=1,nrens
   do k=0,nrt
   do i=1,nx
      fullcor(i,k,1)=fullcor(i,k,1) + (full(k,j)%ocean(i)-fullave(i,k,1)) * (full(kc,j)%ocean(ic)-fullave(ic,kc,1))
      fullcor(i,k,2)=fullcor(i,k,2) + (full(k,j)%atmos(i)-fullave(i,k,2)) * (full(kc,j)%atmos(ic)-fullave(ic,kc,2))
   enddo
   enddo
   enddo
   fullcor=fullcor/real(nrens-1)
   do k=0,nrt
   do i=1,nx
      fullcor(i,k,1)=fullcor(i,k,1)/(fullvar(i,k,1)*fullvar(ic,kc,1))
      fullcor(i,k,2)=fullcor(i,k,2)/(fullvar(i,k,2)*fullvar(ic,kc,2))
   enddo
   enddo
#endif

! corariance oa and ao
#ifdef ZEROD
   do j=1,nrens
   do k=0,nrt
   do i=1,nx
      fullcoroa(i,k,1)=fullcoroa(i,k,1) + (full(k,j)%atmos-fullave(i,k,2)) * (full(kc,j)%ocean-fullave(nx,kc,1))
      fullcoroa(i,k,2)=fullcoroa(i,k,2) + (full(k,j)%ocean-fullave(i,k,1)) * (full(kc,j)%atmos-fullave(nx,kc,2))
   enddo
   enddo
   enddo
   fullcoroa=fullcor/real(nrens-1)
   do k=0,nrt
   do i=1,nx
      fullcoroa(i,k,1)=fullcoroa(i,k,1)/(fullvar(i,k,2)*fullvar(nx,kc,1))
      fullcoroa(i,k,2)=fullcoroa(i,k,2)/(fullvar(i,k,1)*fullvar(nx,kc,2))
   enddo
   enddo
#else
   do j=1,nrens
   do k=0,nrt
   do i=1,nx
      fullcoroa(i,k,1)=fullcoroa(i,k,1) + (full(k,j)%atmos(i)-fullave(i,k,2)) * (full(kc,j)%ocean(ic)-fullave(ic,kc,1))
      fullcoroa(i,k,2)=fullcoroa(i,k,2) + (full(k,j)%ocean(i)-fullave(i,k,1)) * (full(kc,j)%atmos(ic)-fullave(ic,kc,2))
   enddo
   enddo
   enddo
   fullcoroa=fullcor/real(nrens-1)
   do k=0,nrt
   do i=1,nx
      fullcoroa(i,k,1)=fullcoroa(i,k,1)/(fullvar(i,k,2)*fullvar(ic,kc,1))
      fullcoroa(i,k,2)=fullcoroa(i,k,2)/(fullvar(i,k,1)*fullvar(ic,kc,2))
   enddo
   enddo
#endif



end subroutine
end module

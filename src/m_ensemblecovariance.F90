module m_ensemblecovariance
contains
subroutine ensemblecovariance(A,ave,cov,nrens)
   use mod_dimensions
   use mod_state
   implicit none
   integer, intent(in) :: nrens
   type(state), intent(in)    :: A(nrens)
   type(state), intent(in)    :: ave
   type(state), intent(inout) :: cov(2)
   integer j
   #ifndef ZEROD
   integer :: ix=nx/3
   #endif
   cov(1)=0.0
   cov(2)=0.0
   do j=1,nrens
      #ifdef ZEROD
      cov(1)=cov(1)+(A(j)%ocean-ave%ocean)*(A(j)-ave)
      cov(2)=cov(2)+(A(j)%atmos-ave%atmos)*(A(j)-ave)
      #else
      cov(1)=cov(1)+(A(j)%ocean(ix)-ave%ocean(ix))*(A(j)-ave)
      cov(2)=cov(2)+(A(j)%atmos(2*ix)-ave%atmos(2*ix))*(A(j)-ave)
      #endif
   enddo
   cov(1)=(1.0/real(nrens-1))*cov(1)
   cov(2)=(1.0/real(nrens-1))*cov(2)

end subroutine ensemblecovariance
end module m_ensemblecovariance

module m_gnuplot
contains
subroutine gnuplot(fname,variable,nrt,outdir)
   use mod_dimensions
   use mod_state
   character(len=*), intent(in) :: fname
   integer, intent(in) :: nrt
   type(state), intent(in) :: variable(0:nrt)
   character(len=25), intent(in) :: outdir
   integer k,i

   print *,trim(fname)//'o.dat'
   open(10,file=trim(outdir)//'/'//trim(fname)//'o.dat')
      write(10,'(1025i5)')nx,(i,i=1,nx)
      do k=0,nrt
         #ifdef ZEROD
         write(10,'(i5,1024g12.4)')k,variable(k)%ocean
         #else
         write(10,'(i5,1024g12.4)')k,variable(k)%ocean(:)
         #endif
      enddo
   close(10)

   print *,trim(fname)//'a.dat'
   open(10,file=trim(outdir)//'/'//trim(fname)//'a.dat')
      write(10,'(1025i5)')nx,(i,i=1,nx)
      do k=0,nrt
         #ifdef ZEROD
         write(10,'(i5,1024g12.4)')k,variable(k)%atmos
         #else
         write(10,'(i5,1024g12.4)')k,variable(k)%atmos(:)
         #endif
      enddo
   close(10)

end subroutine
end module

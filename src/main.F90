program main
   use mod_dimensions                      ! Defines state dimension
   use mod_state                           ! Defines model state
   use mod_observation                     ! Defines observation state
   use m_readinfile                        ! Reading infile.in
   use m_set_random_seed2                  ! Sets new random seed unless 'seed.dat' exists
   use m_model                             ! Model time stepping
!#ifdef ZEROD
!   use m_recharge_osc                      ! Zero-D recharge oscillator ! not needed, is in m_model
!#endif
   use m_pseudo1D                          ! Samples pseudo-random fields in 1-D
   use m_fixsample1D                       ! Ensures zero mean and unit variance of sampled ensembles
   use m_obsxloc                           ! Computes x position of measurements
   use m_obstloc                           ! Computes time of measurements
   use m_obscount                          ! Counts the number of measurements in a DA window
   use m_obspoints                         ! Saves the measurement locations in space and time for plotting
   use m_dumpsol                           ! Saves outputs for gnuplot animations (p.gnu)
   use m_ensemblemean                      ! Compute ensemble mean at current time (from mem)
   use m_ensemblevariance                  ! Compute ensemble variance at current time (from mem)
   use m_ensemblecovariance                ! Compute ensemble covariance at current time (from mem)
   use m_enkfprep                          ! Setting up predicted measurements and matrices for analysis
   use m_random                            ! Generate random normal numbers
   use m_tecfld                            ! Tecplot output (not used)
   use mod_shapiro                         ! Shapiro filter in case it is needed
   use m_windowstat                        ! Compute ensemble mean and variance in a DA window (from win)
   use m_covstat                           ! Computes space time covariances in case you have a lot of memory
   use m_gnuplot                           ! Generate space time plots for plotting in gnuplot (c.gnu)
#ifdef ZEROD
   use ro_constants                        ! Constants for recharge oscillator
   use parameter_functions                 ! Functions for recharge oscillator
#endif
   implicit none

   type(state), allocatable :: full(:,:)   ! The ensemble of realizations over the whole simulation
   type(state), allocatable :: win(:,:)    ! The ensemble of realizations over an assimilation window.
   type(state), allocatable :: winana(:)   ! The analytical solution over an assimilation window.
   type(state), allocatable :: mem(:)      ! The ensemble of realizations
   type(state), allocatable :: old(:)      ! The ensemble of realizations in previous time step used in Leapfrog
   type(state), allocatable :: sysnoise(:) ! The ensemble of sampled system noise updated every timestep
#ifdef ZEROD
   real, allocatable :: randat(:) ! randomized vector to add to intial conditions for initialisation atmosphere
   real, allocatable :: randoc(:) ! randomized vector to add to intial conditions for initialisation ocean
#endif
   type(state) ana                         ! analytical solution
   type(state) anaold                      ! analytical solution in previous time step used in Leapfrog
   type(state) ave                         ! ensemble average
   type(state) var                         ! ensemble variance
   type(state) cov(2)                      ! ensemble covariance for two points

#ifdef ZEROD
   real, allocatable :: samples(:)       ! work array used when sampling in random
#else
   real, allocatable :: samples(:,:)       ! work array used when sampling in pseudo1D
#endif

! Spacew time statistics diagnostic variables
   type(state), allocatable :: mean(:)     ! ensemble average as a function of space and time
   type(state), allocatable :: stdt(:)     ! ensemble std dev as a function of space and time
   type(state), allocatable :: covo(:)     ! ensemble std dev as a function of space and time
   type(state), allocatable :: cova(:)     ! ensemble std dev as a function of space and time

! Observation location variables
   integer, allocatable :: obsoloc(:)      ! location of ocean observations in space
   integer, allocatable :: obsaloc(:)      ! location of atmos observations in space
   integer, allocatable :: obsotimes(:)    ! location of ocean observations in time
   integer, allocatable :: obsatimes(:)    ! location of atmos observations in time

   type(observation), allocatable :: obs(:)! Stores all observation information in a DA window

! EnKF analysis variables
   real, allocatable :: S(:,:)             ! Predicted meaurements anomalies
   real, allocatable :: E(:,:)             ! Measurement perturbations
   real, allocatable :: D(:,:)             ! Innovations
   real, allocatable :: R(:,:)             ! Measurement error covariance matrix
   real, allocatable :: meanS(:)           ! Mean of predicted measurements
   real, allocatable :: innov(:)           ! Mean innovation 

! Shapiro filter variables
   real, allocatable :: x(:)               ! work array for shapiro filter (input)
   real, allocatable :: y(:)               ! work array for shapiro filter (output)


! other variables
   integer tini                            ! start time of an assimilation window
   integer tfin                            ! end   time of an assimilation window
   integer nrobs                           ! Total number of measurement per assimilation window
   logical leuler                          ! Euler timestepping if true

   integer j,k,l
   real time

   call readinfile()

   call set_random_seed2

   call system('rm -f eigenvalues.dat')
   call system('rm -f obsloc?.dat')
!   call system('touch obsloca.dat obsloco.dat')

! Now that we know all dimensions, allocate the main arrays
   if (lglobstat) allocate (full(0:nrt,nrens))
   allocate (win(0:nrw,nrens))
   allocate (winana(0:nrw))
   allocate (mean(0:nrt))
   allocate (stdt(0:nrt))
   allocate (covo(0:nrt))
   allocate (cova(0:nrt))
   allocate (mem(nrens))
   allocate (old(nrens))
   allocate (sysnoise(nrens))
#ifdef ZEROD
   allocate (samples(nrens))
   allocate (randat(nx))
   allocate (randoc(nx))
#else
   allocate (samples(nx,nrens))
#endif

   allocate (obsoloc(nro))
   allocate (obsaloc(nra))
   allocate (obsotimes(nrt))
   allocate (obsatimes(nrt))

! Shapiro filter
   if (nsh > 0) then
      allocate (x(nx))
      allocate (y(nx))
   endif

   time=0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Uniformly distributed measurements for ocean and atmosphere in space and time
! note, for ZEROD this nro and nra have to be 1
!
   print *,'Ocean observation locations'
   call obsxloc(nro,obsoloc)
   call obstloc(nrt,obst0o,obsdto,obsotimes)
   call obspoints(trim(outdir)//'/obsloco',obsotimes,nrt,obsoloc,nro)

   print *,'Atmos observation locations'
   call obsxloc(nra,obsaloc)
   call obstloc(nrt,obst0a,obsdta,obsatimes)
   call obspoints(trim(outdir)//'/obsloca',obsatimes,nrt,obsaloc,nra)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The true solution smooth pseudo random field drawn from  N(0,1,rh).
! FV Here, I initialize the truth (analytical solution)
!
#ifdef ZEROD
! CHECK?
   ana%ocean=ocheight_0
   ana%atmos=sstemp_0
   print *,'main: ana ok'
#else
   call pseudo1D(ana%ocean,nx,1,rh%ocean,dx,nx)
   call pseudo1D(ana%atmos,nx,1,rh%atmos,dx,nx)
   print *,'main: ana ok'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First guess solution stored in ave
#ifdef ZEROD
   call random(randat,nx)
   ! this only works if randat has nx=1
   ave%atmos=ave%atmos+randat(1)
   call random(randoc,nx)
   ave%ocean=ave%ocean+randoc(1)
   ! add random to truth
   ave=(ave + ana)*(1.0/sqrt(2.0))
!   print *,'need to change this', ave%ocean, nx, rh%ocean, dx
#else
   call pseudo1D(ave%ocean,nx,1,rh%ocean,dx,nx)
   call pseudo1D(ave%atmos,nx,1,rh%atmos,dx,nx)
   ! sqrt 2 to maintain correct variance
   ave=(ave + ana)*(1.0/sqrt(2.0))
   print *,'main: fg ok'
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open output file -obsolete
!
!#ifdef ZEROD
!!
!    print '(a,i4,a,f10.2,a,f10.2,a,f10.2,a,i10)','nop= ',nop,', dtro=',dtro,'epsilon=',eps,' mu=',mu,' nrw ',nrw
!    open(10,file='solution_ro.dat')
!    write(10,'(3a11)')'         t','  Height_Oc','  SST_Atmos'
!#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization of ensemble
!

#ifdef ZEROD
   call random(samples,nx)
#else
   call pseudo1D(samples,nx,nrens,rh%ocean,dx,nx)
   if (samp_fix) call fixsample1D(samples,nx,nrens)
#endif
   do j=1,nrens
#ifdef ZEROD
! ensemble mean should be some sample
      mem(j)%ocean=ocheight_0+samples(j)
#else
      mem(j)%ocean(1:nx)=samples(1:nx,j)
#endif
   enddo

#ifdef ZEROD
   call random(samples,nx)
#else
   call pseudo1D(samples,nx,nrens,rh%atmos,dx,nx)
   if (samp_fix) call fixsample1D(samples,nx,nrens)
#endif
   do j=1,nrens
#ifdef ZEROD
      mem(j)%atmos=sstemp_0+samples(j)
#else
      mem(j)%atmos(1:nx)=samples(1:nx,j)
#endif
   enddo

   do j=1,nrens
      mem(j)%ocean=ave%ocean + sqrt(inivar%ocean)*mem(j)%ocean
      mem(j)%atmos=ave%atmos + sqrt(inivar%atmos)*mem(j)%atmos
   enddo
   print *,'main: ensemble ok'

   nrobs=0
   call ensemblemean(mem,ave,nrens)
   print *,'main: ensemble mean ok'
   call ensemblevariance(mem,ave,var,nrens)
   print *,'main: ensemble variance ok'
   call ensemblecovariance(mem,ave,cov,nrens)
   print *,'main: ensemble covariance ok'
   call dumpsol(time,ana,ave,var,cov,nx,dx,obs,nrobs,mem,nrens,outdir,'I')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop over assimilation windows
   print *,'main: start time stepping'
   do l=1,nrwindows                                   ! Integrate model over the assimilation window
      winana(0)=ana                                   ! store reference solution at tini in winana
      win(0,:)=mem(:)                                 ! store ensemble at tini in win
      tini=(l-1)*nrw                                  ! time at beginning of assimilation window
      tfin=l*nrw                                      ! time at end of assimilation window
      print '(a,i3.3,a,i5.5,a,i5.5)','DA window number ',l,' running from ',tini,' to ',tfin
#ifdef ZEROD
      leuler=.false.                                   ! first timestep of each window is Euler step
#else
      leuler=.true.                                   ! first timestep of each window is Euler step
#endif
      do k=1,nrw                                      ! timestep loop over assimilation window
         if (k>1) leuler=.false.
         time=real((l-1)*nrw+k)
         print '(a,f10.2,a,i3,a,i4,a,l1)','time= ',time,', window=',l,', timestep=',k,' Euler step=',leuler

! Advection
         do j=1,nrens
#ifdef ZEROD
           call recharge_osc(mem(j), old(j), leuler, time)
#else
           call model(mem(j),old(j),leuler)
#endif
         enddo
#ifdef ZEROD
           call recharge_osc(ana, anaold, leuler, time)
#else
           call model(ana,anaold,leuler)
#endif

#ifdef ZEROD
! System noise
         if (sysvar%ocean > 0.0) then
            call pseudo1D(samples,nx,nrens,rh%ocean,dx,nx)
            if (samp_fix) call fixsample1D(samples,nx,nrens)
            do j=1,nrens
! FV do I need to use dtro for ZEROD?
               mem(j)%ocean=mem(j)%ocean+sqrt(2.0*sysvar%ocean*dt)*samples(j)
            enddo
         endif
         if (sysvar%atmos > 0.0) then
            call pseudo1D(samples,nx,nrens,rh%atmos,dx,nx)
            if (samp_fix) call fixsample1D(samples,nx,nrens)
            do j=1,nrens
! FV do I need to use dtro for ZEROD?
               mem(j)%atmos=mem(j)%atmos+sqrt(2.0*sysvar%atmos*dt)*samples(j)
            enddo
         endif

#else
! System noise
         if (sysvar%ocean > 0.0) then
            call pseudo1D(samples,nx,nrens,rh%ocean,dx,nx)
            if (samp_fix) call fixsample1D(samples,nx,nrens)
            do j=1,nrens
               mem(j)%ocean=mem(j)%ocean+sqrt(2.0*sysvar%ocean*dt)*samples(:,j)
! FV do I need to use dtro for ZEROD?
            enddo
         endif
         if (sysvar%atmos > 0.0) then
            call pseudo1D(samples,nx,nrens,rh%atmos,dx,nx)
            if (samp_fix) call fixsample1D(samples,nx,nrens)
            do j=1,nrens
               mem(j)%atmos=mem(j)%atmos+sqrt(2.0*sysvar%atmos*dt)*samples(:,j)
! FV do I need to use dtro for ZEROD?
            enddo
         endif

! Shapiro filter in case we need it
         if (nsh > 0) then
            x=mem(j)%ocean
            call  shfilt(nsh,sh,nx,x,1,y,1,nsh)
            mem(j)%ocean=y
            x=mem(j)%atmos
            call  shfilt(nsh,sh,nx,x,1,y,1,nsh)
            mem(j)%atmos=y
         endif
#endif

! Store ensemble and reference solution for window
         win(k,:)=mem(:)
         winana(k)=ana
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Close output file -obsolete
!
!#ifdef ZEROD
!      close(10)
!#endif


! Counting the number of measurements in the current DA window
      nrobs=obscount(nrt,tini,tfin,obsotimes,obsatimes,nro,nra)
      print *,'Total number of measurements in the DA window= ',nrobs

      if ((nrobs > 0) .and. (mode_analysis > 0)) then
! Assimilation step
         allocate(obs(nrobs))
         allocate(S(nrobs,nrens))
         allocate(E(nrobs,nrens))
         allocate(D(nrobs,nrens))
         allocate(meanS(nrobs))
         allocate(innov(nrobs))
         allocate(R(nrobs,nrobs))


! Assigning the active observations to obs and the predicted observations to S and setting up for analysis
         call enkfprep(mem,obs,S,E,D,meanS,R,innov,winana,win,nrobs,l,tini,tfin,obsoloc,obsaloc,obsotimes,obsatimes)

! Diagnostics at end of DA window
         call ensemblemean(mem,ave,nrens)
         call ensemblevariance(mem,ave,var,nrens)
         call ensemblecovariance(mem,ave,cov,nrens)
         call dumpsol(time,ana,ave,var,cov,nx,dx,obs,nrobs,mem,nrens,outdir,'F')

! ES update of DA window
         print '(a,i2)','Calling analysis with mode: ',mode_analysis
         call analysis(win, R, E, S, D, innov, ndim*(nrw+1), nrens, nrobs, .true., truncation, mode_analysis, &
                        lrandrot, lupdate_randrot, lsymsqrt, inflate, infmult, 1)

! Ensemble at end of DA window (used in continued integration)
         mem(:)=win(nrw,:)

! Diagnostics at end of DA window
         call ensemblemean(mem,ave,nrens)
         call ensemblevariance(mem,ave,var,nrens)
         call ensemblecovariance(mem,ave,cov,nrens)
         call dumpsol(time,ana,ave,var,cov,nx,dx,obs,nrobs,mem,nrens,outdir,'A')
         deallocate(obs,S,E,D,meanS,R,innov)
      endif

! Only needed to compute space-time covariances
      if (lglobstat) full(tini:tfin,:)=win(0:nrw,:)

! Compute ensemble mean and variance over current DA window
      call windowstat(win,nrw,nrens,mean(tini:tfin),stdt(tini:tfin))
   enddo

! Dumping mean and variance as function of space and time
   call gnuplot('gnu_ave',mean,nrt,outdir)
   call gnuplot('gnu_std',stdt,nrt,outdir)

   if (lglobstat) then
      print *,'calling full covariance statistics'
      call covstat(full,nrt,nrens,mean,stdt,covo,cova,outdir)
      call gnuplot('gnu_covo',covo,nrt,outdir)
      call gnuplot('gnu_cova',cova,nrt,outdir)

   endif

end program main


program H2O_vacf
   implicit none
   integer,parameter:: nstep=40000
   character :: b
   real(8) :: velocity(nstep,216*9),cor_dip(150000),Tmax,dt
   real(8) :: sumu,T=298.15d0,vacf
   integer :: i,j,k ,N,tau
   integer,parameter:: tau_nstep=4000
   real(8)::elec(3*216),vtmp(3),dotdip(tau_nstep,3),uacf(0:tau_nstep)
   open(14,file='mdvel')
   read(14,*) b
   do j=1,nstep
      read(14,*) velocity(j,:)
   enddo
   print *,'finish read velocity'

   do i=0,215
      elec(3*i+1)=-2.d0
      elec(3*i+2)=1.d0
      elec(3*i+3)=1.d0
   enddo
   ! calculate dot (dip) auto-correlation functions
   ! uacf(0:ncorr), vel(0:nsteps, ndof)
   call dotdipacf(uacf, velocity, elec, 216*3, 9, nstep, tau_nstep)

   open(13,file='vacf_test.txt')
   do j=0,tau_nstep
      write(13,*) uacf(j)
   enddo
   close(13)

end program

subroutine dotdipacf(uacf, vel, elec, natoms, ndof, nsteps, ncorr)
   integer :: ncorr
   integer :: natoms
   integer :: ndof
   real(8) :: vel(nsteps, ndof)
   real(8) :: elec(natoms)
   real(8) :: uacf(0:ncorr)
   real(8) :: acf(0:ncorr, 3)

   integer :: i, j, k
   real(8) :: dotdip(nsteps,3)
   real(8) :: vtmp(3)

   ! loop over all steps
   do i = 1, nsteps
      vtmp = 0d0
      k = 0
      ! loop over all atoms
      do j = 1, natoms
         vtmp(1) = vtmp(1) + elec(j) * vel(i, k+1)
         vtmp(2) = vtmp(2) + elec(j) * vel(i, k+2)
         vtmp(3) = vtmp(3) + elec(j) * vel(i, k+3)
         k = k + 3
      end do
      dotdip(i,1:3) = vtmp(1:3)
   end do

   do i = 1,3
      call autocorr(acf(0:ncorr,i), dotdip(1:nsteps,i), nsteps, ncorr)
   end do
   do i = 0, ncorr
      uacf(i) = acf(i,1) + acf(i,2) + acf(i,3)
   end do

   return
end subroutine dotdipacf

subroutine autocorr(acf, A, iTtotal, iTcorr)
   implicit none
   integer, intent(in)  ::  iTtotal       ! total number of A
   integer, intent(in)  ::  iTcorr        ! correlation length
   real(8), intent(in)  ::  A(1:iTtotal)  ! array of values of A
   real(8), intent(out) ::  acf(0:iTcorr) ! normalized auto-correlation function

   real(8) ::  nACF(0:iTcorr)             ! number of ACF, for normalization of ACF
   integer :: iT, iT0, iT0max, iTau       ! time index of t0+tau, t0, t0max and tau

   ! initialize all the auto-correlation functions to be calculated
   acf(0:iTcorr)  = 0d0
   nACF(0:iTcorr) = 0d0

   ! for each t0 from 1 to total, calculate auto-correlation function
   do iT0 = 1, iTtotal

      ! the max index of T0 must not exceed the last index
      iT0max = min(iTtotal, iT0+iTcorr)

      ! for all possible t0+tau, calculate a(t0)*a(t0+tau) and add them
      do iT = iT0, iT0max

         ! interval between t and t0
         iTau  = iT - iT0

         ! calculate auto-correlation function for current interval
         acf(iTau) = acf(iTau) + A(iT0) * A(iT)

         ! count the number of auto-correlation function
         nACF(iTau) = nACF(iTau) + 1d0
      end do

   end do

   ! normalize auto-correlation function
   do iT = 0, iTcorr
      acf(iT) = acf(iT) / nACF(iT)
   end do

   return
end subroutine autocorr

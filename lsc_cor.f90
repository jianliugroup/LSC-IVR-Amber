! ===========================================================
! Collect the LSC-IVR correlation function data
! 'NO pgf90 -fast' !!!
! optimization problem should be solved later
! ===========================================================
! -------------------------------------------
    Module lsccor_data
      implicit none
      save
      integer, parameter:: natom = 648
     integer, parameter :: nmolc = 216 
      integer, parameter:: nstep = 15000       ! edit according to settings
      real*8, parameter:: dtstep = 0.0002      ! edit according to settings
      real*8, parameter:: eH = 7.653366E+00    ! charge can be found 
      real*8, parameter:: eO = -1.5306732E+01  ! in prmtop file
      real*8, parameter:: mH = 1.0d0           ! mass can be found 
      real*8, parameter:: mO = 16.0d0          ! in prmtop file
      real*8, parameter:: beta_lsc = 1.68866887       ! AMBER internal unit 
      real*8, parameter:: kb =  8.31441d-3 / 4.184d0  ! AMBER internal unit
!       
    end Module lsccor_data
! -------------------------------------------
    Program lsc_cor
      use lsccor_data
      implicit none 
!      
      character*20:: force_field
      character*50:: line_space
      real*8:: pkubo(3*natom)
      integer:: ipkubo(3*natom) 
      real*8:: x(3*natom,0:nstep), v(3*natom,0:nstep)
      real*8:: xbox(3)
      real*8:: mass(3*natom)
      real*8:: jump(3*natom)
      real*8:: vcenter(3*natom/3,0:nstep)
      real*8:: v2_freq(3*natom), alpha_lsc(3*natom)
      real*8:: rho_a(3*natom)
      real*8:: dipole_0(3), dipole_t(3)
      real*8:: vkubo_sum_H0(3), v_sum_H_t(3,0:nstep), vstd_sum_H0(3)
      real*8:: vkubo_sum_water0(3), v_sum_water_t(3,0:nstep), vstd_sum_water0(3)
      real*8:: cor_H_std(0:nstep), cor_H_kubo(0:nstep)
      real*8:: cor_water_std(0:nstep), cor_water_kubo(0:nstep)
      real*8:: sum_cor_H_std(0:nstep), sum_cor_H_kubo(0:nstep)
      real*8:: sum_cor_water_std(0:nstep), sum_cor_water_kubo(0:nstep) 
      real*8:: cor_lsc(0:nstep)
      real*8:: cor_vv_O(0:nstep)
      real*8:: cor_vv_H(0:nstep)
      real*8:: cor_dipole(0:nstep)
      real*8:: cor_vv_center(0:nstep)
      real*8:: cor_vv_kubo(0:nstep)
      real*8:: sum_cor_lsc(0:nstep)
      real*8:: sum_cor_vv_O(0:nstep)
      real*8:: sum_cor_vv_H(0:nstep)
      real*8:: sum_cor_vv_center(0:nstep)
      real*8:: sum_cor_vv_kubo(0:nstep)
      real*8:: sum_cor_dipole(0:nstep)
      real*8:: time(0:nstep)
      integer:: itemp
      integer:: i, j, k, m
      integer:: ntraj
      real*8:: vvsum(0:nstep)
      real*8:: pc0_kubo(natom/3*3)
      real*8:: pct_kubo(natom/3*3)
      real*8:: vH_kubo0(3*natom/3*2), vH2_kubo0(3*natom/3)
      real*8:: vH_center(3*natom/3*2,0:nstep), vH2_center(3*natom/3,0:nstep)
      real*8:: cor_H_center_kubo(0:nstep), cor_H2_center_kubo(0:nstep)
      real*8:: sum_cor_H_center_kubo(0:nstep), sum_cor_H2_center_kubo(0:nstep)
      real*8:: dxOH, xH1, xH2
      real*8:: dyOH, yH1, yH2
      real*8:: dzOH, zH1, zH2
      real*8:: p2_kubo, sum_p2k
      real(8) :: pO_kubo(nmolc*3), vO(nmolc*3, 0:nstep) ! 216 O atoms, x,y,z      
      real(8) :: vO_acf_kubo(0:nstep), avg_vO_acf_kubo(0:nstep) 
      real(8) :: pH1_kubo(nmolc*3), vH1(nmolc*3, 0:nstep) ! 216 H1 atoms, x,y,z      
      real(8) :: pH2_kubo(nmolc*3), vH2(nmolc*3, 0:nstep) ! 216 H2 atoms, x,y,z      
      real(8) :: vH1_acf_kubo(0:nstep), vH2_acf_kubo(0:nstep) 
      real(8) :: vH_acf_kubo(0:nstep), avg_vH_acf_kubo(0:nstep) 
      
!
      open(11,file='mdcrd')
      open(12,file='mdvel')
      open(13,file='LSC_rhoa.dat')
      open(14,file='LSC_cor-216-vv.dat')
      open(15,file='LSC_screen.dat')

      ! read coordinates and velocities from 1 to nstep from mdcrd and mdvel
      read(11,*)force_field
      read(12,*)force_field
      do j = 1, nstep
         read(11,'(10f8.3)') x(1:3*natom,j)
         read(11,'(3f8.3)') xbox(1:3)
         read(12,'(10f8.3)') v(1:3*natom,j)

         write(15,*)j
      end do

      ! read random seed, initial coordinates and velocities from LSC_rhoa.dat
      read(13,*)i
      read(13,'(10f8.3)')x(:,0)
      read(13,'(10f8.3)')v(:,0)
      read(13,*)mass
      read(13,*)jump
      ! read pkubo from LSC_rhoa.dat
      do i = 1, 3*natom
         read(13,'(f10.3, I6)')pkubo(i), ipkubo(i)
      end do
      do i = 0, nstep
         time(i) = dble(i)*dtstep
      end do
      !
      read(14,*)ntraj
      if(ntraj == 0)then
         sum_cor_H_std = 0.0d0
         sum_cor_water_std = 0.0d0
         sum_cor_H_kubo = 0.0d0
         sum_cor_water_kubo = 0.0d0
         sum_cor_vv_kubo = 0.0d0
         sum_cor_H_center_kubo = 0.0d0
         sum_cor_H2_center_kubo = 0.0d0
         avg_vO_acf_kubo = 0d0
         avg_vH_acf_kubo = 0d0
         rewind(14)
      else
         do i = 0, nstep
            read(14,101)time(i), sum_cor_H_std(i),   &
               sum_cor_water_std(i), sum_cor_H_kubo(i), sum_cor_water_kubo(i), & 
               sum_cor_vv_kubo(i), sum_cor_H_center_kubo(i), sum_cor_H2_center_kubo(i), &
               avg_vO_acf_kubo(i), avg_vH_acf_kubo(i)
         end do
         rewind(14)
      endif
      !
      do i = 0, nstep
         time(i) = dble(i)*dtstep
      end do
      !
      pc0_kubo = 0.0d0
      vH_kubo0 = 0.0d0
      vH2_kubo0 = 0.0d0
      do i = 1, natom/3
         j = 3*i - 2
         pO_kubo(3*i-2) = pkubo(3*j-2)  ! xO
         pO_kubo(3*i-1) = pkubo(3*j-1)  ! yO
         pO_kubo(3*i-0) = pkubo(3*j-0)  ! zO

         pc0_kubo(3*i-2) = pc0_kubo(3*i-2) + pkubo(3*j-2) ! xO
         pc0_kubo(3*i-1) = pc0_kubo(3*i-1) + pkubo(3*j-1) ! yO
         pc0_kubo(3*i) = pc0_kubo(3*i) + pkubo(3*j)       ! zO
         !
         vH2_kubo0(3*i-2) = - 2.0d0*pkubo(3*j-2) / mO     ! qO*xO
         vH2_kubo0(3*i-1) = - 2.0d0*pkubo(3*j-1) / mO     ! qO*yO
         vH2_kubo0(3*i) = - 2.0d0*pkubo(3*j) / mO         ! qO*zO
         !
         k = 2*i -1
         vH_kubo0(3*k-2) = - pkubo(3*j-2) / mO 
         vH_kubo0(3*k-1) = - pkubo(3*j-1) / mO 
         vH_kubo0(3*k) = - pkubo(3*j) / mO
         k = 2*i
         vH_kubo0(3*k-2) = - pkubo(3*j-2) / mO 
         vH_kubo0(3*k-1) = - pkubo(3*j-1) / mO 
         vH_kubo0(3*k) = - pkubo(3*j) / mO
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         j = 3*i - 1
         pH1_kubo(3*i-2) = pkubo(3*j-2)  ! xH1
         pH1_kubo(3*i-1) = pkubo(3*j-1)  ! yH1
         pH1_kubo(3*i-0) = pkubo(3*j-0)  ! zH1

         pc0_kubo(3*i-2) = pc0_kubo(3*i-2) + pkubo(3*j-2)  ! xH1
         pc0_kubo(3*i-1) = pc0_kubo(3*i-1) + pkubo(3*j-1)  ! yH1
         pc0_kubo(3*i) = pc0_kubo(3*i) + pkubo(3*j)        ! zH1
         !
         vH2_kubo0(3*i-2) = vH2_kubo0(3*i-2) + pkubo(3*j-2) / mH  ! qH*xH1
         vH2_kubo0(3*i-1) = vH2_kubo0(3*i-1) + pkubo(3*j-1) / mH  ! qH*yH1
         vH2_kubo0(3*i) = vH2_kubo0(3*i) + pkubo(3*j) / mH        ! qH*zH1
         !
         k = 2*i -1
         vH_kubo0(3*k-2) =  vH_kubo0(3*k-2) + pkubo(3*j-2) / mH 
         vH_kubo0(3*k-1) = vH_kubo0(3*k-1) +  pkubo(3*j-1) / mH 
         vH_kubo0(3*k) = vH_kubo0(3*k) + pkubo(3*j) / mH
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         j = 3*i
         pH2_kubo(3*i-2) = pkubo(3*j-2)  ! xH2
         pH2_kubo(3*i-1) = pkubo(3*j-1)  ! yH2
         pH2_kubo(3*i-0) = pkubo(3*j-0)  ! zH2

         pc0_kubo(3*i-2) = pc0_kubo(3*i-2) + pkubo(3*j-2)  ! xH2
         pc0_kubo(3*i-1) = pc0_kubo(3*i-1) + pkubo(3*j-1)  ! yH2
         pc0_kubo(3*i) = pc0_kubo(3*i) + pkubo(3*j)        ! zH2
         !
         vH2_kubo0(3*i-2) = vH2_kubo0(3*i-2) + pkubo(3*j-2) / mH  ! qH*xH2
         vH2_kubo0(3*i-1) = vH2_kubo0(3*i-1) + pkubo(3*j-1) / mH  ! qH*yH2
         vH2_kubo0(3*i) = vH2_kubo0(3*i) + pkubo(3*j) / mH        ! qH*zH2
         !
         k = 2*i
         vH_kubo0(3*k-2) =  vH_kubo0(3*k-2) + pkubo(3*j-2) / mH 
         vH_kubo0(3*k-1) = vH_kubo0(3*k-1) +  pkubo(3*j-1) / mH 
         vH_kubo0(3*k) = vH_kubo0(3*k) + pkubo(3*j) / mH
      enddo

      vH2_center = 0.0d0    ! Single molecule
      vH_center = 0.0d0     ! Single atom
      vcenter = 0.0d0
      !
      do k = 0, nstep
         do i = 1, natom/3
            j = 3*i - 2
            vO(3*i-2,k) = v(3*j-2, k) 
            vO(3*i-1,k) = v(3*j-1, k)
            vO(3*i,k) = v(3*j, k)

            vH2_center(3*i-2,k) = vH2_center(3*i-2,k) - 2.0d0 * v(3*j-2,k)
            vH2_center(3*i-1,k) = vH2_center(3*i-1,k) - 2.0d0 * v(3*j-1,k)
            vH2_center(3*i,k) = vH2_center(3*i,k) - 2.0d0 * v(3*j,k)
            m = 2*i - 1
            vH_center(3*m-2,k) = vH_center(3*m-2,k) -  v(3*j-2,k)
            vH_center(3*m-1,k) = vH_center(3*m-1,k) -  v(3*j-1,k)
            vH_center(3*m,k) = vH_center(3*m,k) - v(3*j,k)
            m = 2*i
            vH_center(3*m-2,k) = vH_center(3*m-2,k) -  v(3*j-2,k)
            vH_center(3*m-1,k) = vH_center(3*m-1,k) -  v(3*j-1,k)
            vH_center(3*m,k) = vH_center(3*m,k) - v(3*j,k)
            !
            j = 3*i - 1
            vH1(3*i-2,k) = v(3*j-2, k) 
            vH1(3*i-1,k) = v(3*j-1, k)
            vH1(3*i,k) = v(3*j, k)

            vH2_center(3*i-2,k) = vH2_center(3*i-2,k) + v(3*j-2,k)
            vH2_center(3*i-1,k) = vH2_center(3*i-1,k) + v(3*j-1,k)
            vH2_center(3*i,k) = vH2_center(3*i,k) + v(3*j,k)
            m = 2*i - 1
            vH_center(3*m-2,k) = vH_center(3*m-2,k) +  v(3*j-2,k)
            vH_center(3*m-1,k) = vH_center(3*m-1,k) +  v(3*j-1,k)
            vH_center(3*m,k) = vH_center(3*m,k) + v(3*j,k)
            !
            j = 3*i
            vH2(3*i-2,k) = v(3*j-2, k) 
            vH2(3*i-1,k) = v(3*j-1, k)
            vH2(3*i,k) = v(3*j, k)

            vH2_center(3*i-2,k) = vH2_center(3*i-2,k) + v(3*j-2,k)
            vH2_center(3*i-1,k) = vH2_center(3*i-1,k) + v(3*j-1,k)
            vH2_center(3*i,k) = vH2_center(3*i,k) + v(3*j,k)
            m = 2*i
            vH_center(3*m-2,k) = vH_center(3*m-2,k) +  v(3*j-2,k)
            vH_center(3*m-1,k) = vH_center(3*m-1,k) +  v(3*j-1,k)
            vH_center(3*m,k) = vH_center(3*m,k) + v(3*j,k)
            !
            ! --------------------------------------------------------------
            j = 3*i - 2
            vcenter(3*i-2,k) = vcenter(3*i-2,k) + mO * v(3*j-2,k)
            vcenter(3*i-1,k) = vcenter(3*i-1,k) + mO * v(3*j-1,k)
            vcenter(3*i,k) = vcenter(3*i,k) + mO * v(3*j,k)
            !
            j = 3*i - 1
            vcenter(3*i-2,k) = vcenter(3*i-2,k) + mH * v(3*j-2,k)
            vcenter(3*i-1,k) = vcenter(3*i-1,k) + mH * v(3*j-1,k)
            vcenter(3*i,k) = vcenter(3*i,k) + mH * v(3*j,k)
            !
            j = 3*i
            vcenter(3*i-2,k) = vcenter(3*i-2,k) + mH * v(3*j-2,k)
            vcenter(3*i-1,k) = vcenter(3*i-1,k) + mH * v(3*j-1,k)
            vcenter(3*i,k) = vcenter(3*i,k) + mH * v(3*j,k)
         end do
      end do
      vcenter = vcenter/(mH + mH + mO)
      !
      cor_H_std = 0.0d0
      cor_water_std = 0.0d0
      cor_H_kubo = 0.0d0
      cor_water_kubo = 0.0d0
      cor_vv_kubo = 0.0d0
      cor_H_center_kubo = 0.0d0
      cor_H2_center_kubo = 0.0d0
      vO_acf_kubo = 0d0
      vH1_acf_kubo = 0d0
      vH2_acf_kubo = 0d0

      ! collective H-atom kubo-transformed velocity correlation function pkubo_sum_H0
      ! collective relative H-atom kubo-transformed velocity correlation function pkubo_sum_water0
      vkubo_sum_H0 = 0.0d0         ! only H-atom sum kubo
      vkubo_sum_water0 = 0.0d0     ! water molecules sum kubo  
      vstd_sum_H0 = 0.0d0          ! only H-atom sum std
      vstd_sum_water0 = 0.0d0      ! water molecules sum std
      v_sum_H_t = 0.0d0          
      v_sum_water_t = 0.0d0
      do i = 1, natom/3
         j = 3*i - 1
         vkubo_sum_H0(1) = vkubo_sum_H0(1) + pkubo(3*j-2) / mH
         vkubo_sum_H0(2) = vkubo_sum_H0(2) + pkubo(3*j-1) / mH
         vkubo_sum_H0(3) = vkubo_sum_H0(3) + pkubo(3*j) / mH
         vkubo_sum_water0(1) = vkubo_sum_water0(1) + pkubo(3*j-2) / mH
         vkubo_sum_water0(2) = vkubo_sum_water0(2) + pkubo(3*j-1) / mH
         vkubo_sum_water0(3) = vkubo_sum_water0(3) + pkubo(3*j) / mH
         !
         j = 3*i
         vkubo_sum_H0(1) = vkubo_sum_H0(1) + pkubo(3*j-2) / mH
         vkubo_sum_H0(2) = vkubo_sum_H0(2) + pkubo(3*j-1) / mH
         vkubo_sum_H0(3) = vkubo_sum_H0(3) + pkubo(3*j) / mH
         vkubo_sum_water0(1) = vkubo_sum_water0(1) + pkubo(3*j-2) / mH
         vkubo_sum_water0(2) = vkubo_sum_water0(2) + pkubo(3*j-1) / mH
         vkubo_sum_water0(3) = vkubo_sum_water0(3) + pkubo(3*j) / mH
         !
         j = 3*i - 2
         vkubo_sum_water0(1) = vkubo_sum_water0(1) - 2.0d0*pkubo(3*j-2) / mO
         vkubo_sum_water0(2) = vkubo_sum_water0(2) - 2.0d0*pkubo(3*j-1) / mO
         vkubo_sum_water0(3) = vkubo_sum_water0(3) - 2.0d0*pkubo(3*j) / mO        

      end do
      !
      do k = 0, nstep
         do i = 1, natom/3
            j = 3*i - 2
            v_sum_water_t(1,k) = v_sum_water_t(1,k) - 2.0d0 * v(3*j-2,k)
            v_sum_water_t(2,k) = v_sum_water_t(2,k) - 2.0d0 * v(3*j-1,k)
            v_sum_water_t(3,k) = v_sum_water_t(3,k) - 2.0d0 * v(3*j,k)
            !
            j = 3*i - 1
            v_sum_water_t(1,k) = v_sum_water_t(1,k) + v(3*j-2,k) !* mH
            v_sum_water_t(2,k) = v_sum_water_t(2,k) + v(3*j-1,k) !* mH
            v_sum_water_t(3,k) = v_sum_water_t(3,k) + v(3*j,k)   !* mH
            v_sum_H_t(1,k) = v_sum_H_t(1,k) + v(3*j-2,k)
            v_sum_H_t(2,k) = v_sum_H_t(2,k) + v(3*j-1,k)
            v_sum_H_t(3,k) = v_sum_H_t(3,k) + v(3*j,k)
            !
            j = 3*i
            v_sum_water_t(1,k) = v_sum_water_t(1,k) + v(3*j-2,k) !* mH
            v_sum_water_t(2,k) = v_sum_water_t(2,k) + v(3*j-1,k) !* mH
            v_sum_water_t(3,k) = v_sum_water_t(3,k) + v(3*j,k)   !* mH 
            v_sum_H_t(1,k) = v_sum_H_t(1,k) + v(3*j-2,k)
            v_sum_H_t(2,k) = v_sum_H_t(2,k) + v(3*j-1,k)
            v_sum_H_t(3,k) = v_sum_H_t(3,k) + v(3*j,k)
         end do
      end do
      vstd_sum_H0 = v_sum_H_t(:,0)
      !
      !
      vstd_sum_water0 = v_sum_water_t(:,0)  ! / (mH + mH + mO)

      do j = 0, nstep
         ! kubo-transform of single water molecules velocity correlation function
         cor_vv_kubo(j) = cor_vv_kubo(j) + dot_product(pc0_kubo, vcenter(:,j))   
         ! single vH-vO  kubo-transform        
         cor_H_center_kubo(j) = cor_H_center_kubo(j) + dot_product(vH_kubo0,vH_center(:,j))
         ! molecular vH-vO kubo-transform
         cor_H2_center_kubo(j) = cor_H2_center_kubo(j) + dot_product(vH2_kubo0,vH2_center(:,j))
         !
         ! sum vH  standard
         cor_H_std(j) = cor_H_std(j) + dot_product(vstd_sum_H0, v_sum_H_t(:,j))
         ! sum vH-vO  standard
         cor_water_std(j) = cor_water_std(j) + dot_product(vstd_sum_water0, v_sum_water_t(:,j))
         ! sum vH  Kubo-transform
         cor_H_kubo(j) = cor_H_kubo(j) + dot_product(vkubo_sum_H0, v_sum_H_t(:,j))
         ! sum vH-vO  Kubo-transform
         cor_water_kubo(j) = cor_water_kubo(j) + dot_product(vkubo_sum_water0, v_sum_water_t(:,j))

         vO_acf_kubo(j) = vO_acf_kubo(j) + dot_product(pO_kubo, vO(:,j))
         vH1_acf_kubo(j) = vH1_acf_kubo(j) + dot_product(pH1_kubo, vH1(:,j))
         vH2_acf_kubo(j) = vH2_acf_kubo(j) + dot_product(pH2_kubo, vH2(:,j))
      end do ! j  

      cor_H_std = cor_H_std / dble(2*natom/3) * 0.50d0 * mH / kb
      cor_water_std = cor_water_std / dble(2*natom/3) * 0.50d0 * mH / kb
      cor_H_kubo = cor_H_kubo / dble(2*natom/3) * mH 
      cor_water_kubo = cor_water_kubo / dble(2*natom/3) * mH

      cor_H_center_kubo = cor_H_center_kubo / dble(2*natom/3) * mH
      cor_H2_center_kubo = cor_H2_center_kubo / dble(2*natom/3) * mH
      cor_vv_kubo = cor_vv_kubo/dble(natom/3)

      vO_acf_kubo = vO_acf_kubo / dble(nmolc) * mH
      vH_acf_kubo = (vH1_acf_kubo + vH2_acf_kubo) / dble(2*nmolc) * mH
      !
      if(ntraj == 0)then
         sum_cor_H_std = cor_H_std
         sum_cor_water_std = cor_water_std
         sum_cor_H_kubo = cor_H_kubo
         sum_cor_water_kubo = cor_water_kubo
         sum_cor_vv_kubo = cor_vv_kubo
         sum_cor_H_center_kubo = cor_H_center_kubo
         sum_cor_H2_center_kubo = cor_H2_center_kubo
         avg_vO_acf_kubo = vO_acf_kubo
         avg_vH_acf_kubo = vH_acf_kubo
      else
         sum_cor_H_std = cor_H_std/dble(ntraj+1) +  &
            dble(ntraj)*sum_cor_H_std/dble(ntraj+1)
         sum_cor_water_std = cor_water_std/dble(ntraj+1) +  &
            dble(ntraj)*sum_cor_water_std/dble(ntraj+1)
         sum_cor_H_kubo = cor_H_kubo/dble(ntraj+1) +  &
            dble(ntraj)*sum_cor_H_kubo/dble(ntraj+1)
         sum_cor_water_kubo = cor_water_kubo/dble(ntraj+1) +  &
            dble(ntraj)*sum_cor_water_kubo/dble(ntraj+1)
         sum_cor_vv_kubo = cor_vv_kubo/dble(ntraj+1) +  &
            dble(ntraj)*sum_cor_vv_kubo/dble(ntraj+1) 
         sum_cor_H_center_kubo = cor_H_center_kubo/dble(ntraj+1) +  &
            dble(ntraj)*sum_cor_H_center_kubo/dble(ntraj+1)
         sum_cor_H2_center_kubo = cor_H2_center_kubo/dble(ntraj+1) +  &
            dble(ntraj)*sum_cor_H2_center_kubo/dble(ntraj+1)
         avg_vO_acf_kubo = vO_acf_kubo / dble(ntraj+1) &
            + dble(ntraj) * avg_vO_acf_kubo / dble(ntraj+1)
         avg_vH_acf_kubo = vH_acf_kubo / dble(ntraj+1) &
            + dble(ntraj) * avg_vH_acf_kubo / dble(ntraj+1)
      endif
      !
      ntraj = ntraj + 1
      write(14,*)ntraj
      !
      do i = 0, nstep
         write(14,101)time(i), &
            sum_cor_H_std(i), sum_cor_water_std(i), sum_cor_H_kubo(i), sum_cor_water_kubo(i), & 
            sum_cor_vv_kubo(i), sum_cor_H_center_kubo(i), sum_cor_H2_center_kubo(i), &
            avg_vO_acf_kubo(i), avg_vH_acf_kubo(i)
      end do  
      write(14,*)'END'

      101  format(2x,f16.8,9(x,e11.5))
      !
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
   end program lsc_cor

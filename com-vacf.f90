program com_vacf
    implicit none
    integer,parameter::nstep=40000
    integer,parameter::tau_step=4000
    integer,parameter::molecule=216
    integer,parameter::atom=molecule*3
    character ::b
    real(8),parameter::dt=0.0005d0      !time step= 0.5 fs
    real(8) ::velocity(nstep,atom*3),sum_v(3),v_center(molecule,3,nstep)
    real(8) ::vv_acf(0:tau_step,molecule)
    integer ::i,j
    real(8) :: mass(3),sumu

    open(14,file='mdvel')
    read(14,*) b
    do j=1,nstep
        read(14,*) velocity(j,:)
    enddo
    print *, 'finish read velocity'
    close(14)
    mass(1)=16.d0
    mass(2)=1.d0
    mass(3)=1.d0

    call center_velocity(velocity,mass,v_center,nstep,molecule)
    print *,'finish center_velocity'
    call acf(v_center,molecule,nstep,tau_step,vv_acf)
    print *,'finish calculate acf' 

    open(16,file='vv_correlation.txt')

    do i=0,tau_step
        sumu=0.d0
        do j=1,molecule  
            sumu=vv_acf(i,j)+sumu
        enddo
        write(16,*) sumu/molecule
    enddo
    close(16)
    print *,'finish write vv_acf'

end program com_vacf


!....................................................................
!To get the velocity of mass center
!we will get 216 H2O velocities of each directions in every timestep
!....................................................................
subroutine center_velocity(vel,mass,v_center,nstep,molecule)
    implicit none
    real(8),intent(in)::vel(nstep,molecule*9),mass(3)
    integer,intent(in) :: molecule,nstep
    real(8),intent(out) ::v_center(molecule,3,nstep)
    integer::i,j,k
    do i=1,nstep
        do j=0,molecule-1
            v_center(j+1,1,i)=vel(i,9*j+1)*mass(1)+vel(i,9*j+4)*mass(2)+vel(i,9*j+7)*mass(3)
            v_center(j+1,2,i)=vel(i,9*j+2)*mass(1)+vel(i,9*j+5)*mass(2)+vel(i,9*j+8)*mass(3)
            v_center(j+1,3,i)=vel(i,9*j+3)*mass(1)+vel(i,9*j+6)*mass(2)+vel(i,9*j+9)*mass(3)
        enddo
    enddo
    v_center(:,:,:)=v_center(:,:,:)/(mass(1)+mass(2)+mass(3))
    return
end subroutine center_velocity

!.............................................................
!To get the velocity-velocity correlation function of each H2O
!And Then we could average all 216 H2O.
!.............................................................
subroutine acf(v_center,molecule,nstep,tau_step,vv_acf)
    implicit none
    integer,intent(in) ::molecule,nstep,tau_step
    real(8),intent(in) ::v_center(molecule,3,nstep)
    real(8),intent(out) ::vv_acf(tau_step,molecule)
    integer ::it,it0,it0max,itau
    real(8)::uacf(0:tau_step,molecule),nacf(0:tau_step,molecule)
    integer ::i
    uacf(:,:)=0.d0
    nacf(:,:)=0.d0

    do i=1,molecule
        do it0=1,nstep
            it0max=min(nstep,it0+tau_step)
            do it=it0,it0max
                itau=it-it0
                uacf(itau,i)=uacf(itau,i)+dot_product(v_center(i,:,it0),v_center(i,:,it))
                nacf(itau,i)=nacf(itau,i)+1d0
            enddo
        enddo
    enddo

    do i=1,molecule
        do it=0,tau_step
            vv_acf(it,i)=uacf(it,i)/nacf(it,i)
        enddo
    enddo
    return
end subroutine acf

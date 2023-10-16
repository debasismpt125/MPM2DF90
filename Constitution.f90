
    !-------------------------------------------------------------------
    !- MPM3D - A Three Dimensional Explicit Material Point Method Code -
    !-                                                                 -
    !- Developed by                                                    -
    !-    Computational Dynamics Laboratory                            -
    !-    School of Aerospace, Tsinghua University                     -
    !-    Beijing 100084, China.                                       -
    !-                                                                 -
    !-    Professor Xiong Zhang                                        -
    !-    Email: xzhang@tsinghua.edu.cn                                -
    !-    Web  : http://www.comdyn.cn                                  -
    !-                                                                 -
    !-    Copyright (C) 2004 ~ 2012                                    -
    !-------------------------------------------------------------------

    module MaterialModel

    use MaterialData

    private
    public Constitution

    integer:: mid       ! material set id
    integer:: etype_    ! Type of EOS
    integer:: mtype_    ! Type of material model

    real(8):: young_, poisson_
    real(8):: yield0_, tangmod_

    real(8):: den0_     ! Initial density
    real(8):: den_      ! Current density

    real(8):: vold      ! volume of step n
    real(8):: vol_      ! Current volume
    real(8):: vol0_     ! initial volume
    real(8):: dvol      ! 0.5 * volume increment

    real(8):: dinc(6)   ! strain increment
    real(8):: sm        ! Mean stress
	real(8):: sm0       ! initial mean stress
    real(8):: sd(6)     ! deviatoric stress
    real(8):: sig(6)    ! stress components
    real(8):: sold(6)   ! deviatoric stress of step n
    real(8):: dsm
    real(8):: bqf       ! bulk viscosity force

    real(8):: seqv      ! Equivalent stress
    real(8):: epeff_    ! Effective plastic strain
    real(8):: sig_y_    ! Current yield stress
    real(8):: depeff    ! increment of equivalent plastic strain
    real(8):: ratio     ! for hardening caculation

    real(8):: G2, K3, PlaMod ! 2*G, 3*K, plastic hardening modulus

    real(8):: specheat_
    real(8):: tmprt     ! temperature

    real(8):: iener     ! internal energy
    real(8):: specen    ! internal energy per initial volume
    real(8):: ieinc     ! internal energy increment

    real(8):: mu        ! mu = den_/den0_ - 1
    real(8):: rv        ! relative volume rv = vol_/vol0_
    real(8):: bfac      ! burn fraction
    real(8):: cp        ! sound speed
	real(8):: pc        ! pre-consolidation pressure
	real(8):: evoid     ! void ratio
	real(8):: pore
	real(8):: KW
    real(8):: pwpa
	real(8):: pm
	real(8):: qd
	real(8):: ptot
    contains

    subroutine Constitution(de, vort, b, p)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Update stresses by using a constitution model             -
    !-  Inputs:                                                       -
    !-      de - strain increment                                     -
    !-             (D11, D22, D33, 2D23, 2D13, 2D12)*DT               -
    !-      vort  - vorticity * DT (W32, W13, W21)*DT                 -
    !-      b     - body index                                        -
    !-      p     - particle index                                    -
    !-    Outputs:                                                    -
    !-        stress component                                        -
    !-    Note:                                                       -
    !-        sd(i) and de(i) comply the Voigt rule                   -
    !------------------------------------------------------------------
    use ParticleData
    use FFI, only: iomsg

    implicit none

    integer, intent(in):: b, p
    real(8), intent(in):: de(6), vort(3)
    real(8):: mp_, ltim, elastic_strain(6), trial_elastic_strain(6), plastic_strain(6)
    logical:: failure
    type(Particle), POINTER :: pt

    ! Pick parameters
    dinc = de

    ! Volume at time step t
    pt => particle_list(p)
    vold = pt%VOL
    ! Current volume
    !pt%VOL = vold*(1d0+de(1)+de(2)+de(3))
    pt%VOL = pt%mass/mat_list(body_list(b)%mat)%Density*pt%Jacobian
    vol_ = pt%VOL
    dvol = 0.5d0 * (vol_ - vold)

    if (vol_ .le. 0d0) then
        write(*,*) '=== warning: negative volume at particle', p
        write(iomsg,*) '=== warning: negative volume at particle', p
        stop
    end if
    mid = body_list(b)%mat

    epeff_ = pt%epeff
    sig_y_ = pt%sig_y    ! SIGY
    seqv = pt%seqv       ! 2005-8-13

    ltim = pt%LT         ! lighting time
    tmprt = pt%celsius_t ! temperature

    iener = pt%ie           ! internal energy
    mp_ = pt%mass           ! particle mass

    failure = pt%failure

    mtype_ = mat_list(mid)%MatType
    etype_ = mat_list(mid)%EosType

    young_ = mat_list(mid)%Young       ! Young's Modulus
    poisson_ = mat_list(mid)%Poisson   ! Poisson Ratio
    yield0_ = mat_list(mid)%Yield0     ! Yield limit
    tangmod_ = mat_list(mid)%TangMod   ! Tangential modulus

    den0_ = mat_list(mid)%Density      ! Initial density
    vol0_ = mp_/den0_                  ! Initial volume
    den_ = mp_/vol_                    ! Current density
    mu = den_/den0_ - 1d0
    rv = vol_/vol0_                    ! Relative volume
    specen = iener / vol0_

    young_ = young_*(1.0d0-pt%DMG+epsilon(1.0d0))
    G2 = young_ / (1d0 + poisson_)
    K3 = young_ / (1d0 - 2d0*poisson_)
    PlaMod = young_ * tangmod_ / (young_ - tangmod_)

    sm = pt%SM           ! Mean stress
	sm0 = pt%SM0           ! initial mean stress  
    evoid=pt%evoid        ! Void ratio
	pc= pt%pc             !pre-consolidation pressure
	pore=pt%pore
	pwpa= pt%pwpa
	KW = pt%KW
	pm= pt%pm
	qd=pt%qd
	ptot=pt%ptot
    sd(1) = pt%SDxx      ! deviatoric stress
    sd(2) = pt%SDyy
    sd(3) = pt%SDzz
    sd(4) = pt%SDyz
    sd(5) = pt%SDxz
    sd(6) = pt%SDxy

    sig(1) = sd(1) + sm        ! cauchy stress
    sig(2) = sd(2) + sm
    sig(3) = sd(3) + sm
    sig(4) = sd(4)
    sig(5) = sd(5)
    sig(6) = sd(6)

    sold = sd    ! cauchy stress at time step t

    ! Select material model
    select case(mtype_)

    case(1)
        ! elas: elastic model
        call sigrot(vort, sig, sm, sd)    ! Rotate stress
        call M3DM1()
        call lieupd()

    case(2)
        ! pla1: elastic-perfectly plastic
        call sigrot(vort, sig, sm, sd)
        call M3DM2()
        call lieupd()

    case(3)
        ! pla2: isotropic hardening
        call sigrot(vort, sig, sm, sd)
        call M3DM3()
        call lieupd()

    case(4)
        ! john: Johnson-Cook plasticity model
        call sigrot(vort, sig, sm, sd)
        call M3DM4(mat_list(mid), DT, tmprt)
        call lieupd()
        ! call seleos()

        specheat_ = mat_list(mid)%SpecHeat
        pt%celsius_t = pt%celsius_t + &
            seqv*depeff/den_/specheat_

    case(5)
        ! sjc: Simplified Johnson-Cook plasticity model
        call sigrot(vort, sig, sm, sd)
        call M3DM5(mat_list(mid), DT)
        call seleos(failure)

    case(6)
        ! sjcf: Simplified Johnson-Cook plasticity model with
        !       failure
        call sigrot(vort, sig, sm, sd)
        if (.not.pt%failure) then
            call M3DM5(mat_list(mid), DT)
            call seleos(failure)
        else
            call M3DM7()
            call seleos(failure)
        end if

        if (epeff_ .gt. mat_list(mid)%epf) then
            epeff_ = mat_list(mid)%epf + 0.0000001
            pt%failure = .true.
        end if

    case(7)
        ! john: Johnson-Cook plasticity model with failure
        call sigrot(vort, sig, sm, sd)
        if (.not.pt%failure) then
            call M3DM8(mat_list(mid), DT, tmprt)
            call seleos(failure)
            specheat_ = mat_list(mid)%SpecHeat
            pt%celsius_t = pt%celsius_t &
                + seqv*depeff/den_/specheat_
        else
            call M3DM7()
            call seleos(failure)
        end if

        if (epeff_ .gt. mat_list(mid)%epf) then
            epeff_ = mat_list(mid)%epf + 0.0000001
            failure = .true.
        end if

        if (failure .and. sm < 0) then
            sm=0.0
            pt%VOL = vold
        end if
        pt%failure = failure

    case(8)
        ! hiex: High Explosive burn
        call M3DM6(ltim, CurrentTime, mat_list(mid)%D)
        call seleos(failure)

    case(9)
        ! null: used to model air
        call M3DM7()
        call seleos(failure)

    case(10)
        ! Drucker-Prager elastic-perfectly plastic
       ! call sigrot(vort, sig, sm, sd)
        call M3DM10(mat_list(mid), DT)
        call lieupd()

    case(11)
        ! Drucker-Prager elastic-perfectly plastic with piecewise cohesion softening
        !call sigrot(vort, sig, sm, sd)
        elastic_strain = pt%elastic_strain
        trial_elastic_strain = pt%trial_elastic_strain
        call M3DM11(mat_list(mid), DT, trial_elastic_strain, elastic_strain)
        call lieupd()
    
        ! update elastic strain
        pt%elastic_strain = elastic_strain
        
    case(12)
        ! Mohr-Coulomb elastic-perfectly plastic with tension cutting
        call sigrot(vort, sig, sm, sd)
        call M3DM12(mat_list(mid), DT)
        call lieupd()
        
    case(13)
        ! Mohr-Coulomb elastic-perfectly plastic with piecewise hardening
        elastic_strain = pt%elastic_strain
        trial_elastic_strain = pt%trial_elastic_strain

        ! call sigrot(vort, sig, sm, sd)
        call M3DM13(mat_list(mid), DT, trial_elastic_strain, elastic_strain)
        call lieupd()
        
        ! update elastic strain
        pt%elastic_strain = elastic_strain
    case(14)
        ! Modified cam-clay model
        plastic_strain = pt%plastic_strain
		

        call sigrot(vort, sig, sm, sd)
        call M3DM14(mat_list(mid), DT, plastic_strain)
        call lieupd()
        
        ! update elastic strain
        pt%plastic_strain = plastic_strain
      
    case default
        write(*, 10) mtype_
10      format(1x,'*** Stop *** material type ', i2, &
            ' has not been implemented !')
        stop

    end select

    ! Write stress result
    pt%SM0=sm0
    pt%SM = sm
    pt%Seqv = seqv

    pt%SDxx = sd(1)
    pt%SDyy = sd(2)
    pt%SDzz = sd(3)
    pt%SDyz = sd(4)
    pt%SDxz = sd(5)
    pt%SDxy = sd(6)

    pt%sig_y = sig_y_
    pt%epeff = epeff_
    pt%evoid = evoid
	pt%pc = pc
	pt%pore = pore
	pt%KW = KW
	pt%pm = pm
	pt%qd = qd
	pt%ptot = ptot
    pt%ie = iener
    pt%cp = cp

    end subroutine Constitution


    subroutine sigrot(vort, sig, sm, sd)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Rotate stresses, and then update the mean stress sm       -
    !-      and deviate stresses sd                                   -
    !-  Input                                                         -
    !-      vort - Vorticity increments (W32, W13, W21)*DT            -
    !-      sig  - Caucy stresses at time step t                      -
    !-             (S11, S22, S33, S23, S13, S12)                     -
    !-  Output                                                        -
    !-      sig  - Rotated caucy stresses                             -
    !-      sd   - Rotated deviate stress                             -
    !-      sm   - Rotated mean stress                                -
    !-  References                                                    -
    !-      Section 5.1                                               -
    !------------------------------------------------------------------
    implicit none
    real(8), intent(in) :: vort(3)
    real(8), intent(inout) :: sig(6)
    real(8), intent(out) :: sm, sd(6)

    real(8) :: rot(6), q(3)

    q(1) = 2d0*sig(6)*vort(3)
    q(2) = 2d0*sig(5)*vort(2)
    q(3) = 2d0*sig(4)*vort(1)

    rot(1) = - q(1) + q(2) ! (Eq: 5.4)
    rot(2) = + q(1) - q(3)
    rot(3) = - q(2) + q(3)
    rot(4) = vort(1)*(sig(2)-sig(3)) + vort(3)*sig(5) - &
        vort(2)*sig(6)
    rot(5) = vort(2)*(sig(3)-sig(1)) + vort(1)*sig(6) - &
        vort(3)*sig(4)
    rot(6) = vort(3)*(sig(1)-sig(2)) + vort(2)*sig(4) - &
        vort(1)*sig(5)

    sig = sig + rot ! First two terms in RHS of (Eq. 5.3)

    sm = (sig(1)+sig(2)+sig(3))/3d0

    sd(1) = sig(1) - sm
    sd(2) = sig(2) - sm
    sd(3) = sig(3) - sm
    sd(4) = sig(4)
    sd(5) = sig(5)
    sd(6) = sig(6)

    end subroutine sigrot


    subroutine elastic_devi()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Update devitoric stress by elastic relation               -
    !-  Inputs                                                        -
    !-      dinc - strain increment                                   -
    !-             (D11, D22, D33, 2D23, 2D13, 2D12)*DT               -
    !-  Outputs                                                       -
    !-    sd     - devitoric stress component                         -
    !-  References                                                    -
    !-      Section 5.2.1                                             -
    !------------------------------------------------------------------
    implicit none
    real(8):: dem, G

    dem = (dinc(1) + dinc(2) + dinc(3)) / 3.0d0

    G = 0.5d0*G2    ! Shear modulus

    sd(1) = sd(1) + G2*(dinc(1)-dem)    ! (Eq: 5.29)
    sd(2) = sd(2) + G2*(dinc(2)-dem)
    sd(3) = sd(3) + G2*(dinc(3)-dem)
    sd(4) = sd(4) + G*dinc(4)
    sd(5) = sd(5) + G*dinc(5)
    sd(6) = sd(6) + G*dinc(6)

    end subroutine elastic_devi


    subroutine elastic_p()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Update pressure by elastic relation                       -
    !-  Inputs                                                        -
    !       dinc - strain increment                                   -
    !   Outputs                                                       -
    !       sm   - mean stress (pressure)                             -
    !------------------------------------------------------------------
    implicit none
    real(8):: dem

    dem = (dinc(1) + dinc(2) + dinc(3)) / 3.0
    dsm = K3*dem

    sm = sm + dsm        ! (Eq: 5.30)

    end subroutine elastic_p
    
    subroutine compute_elastic_stress_state(elastic_strain)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Update pressure by elastic relation                       -
    !-  Inputs                                                        -
    !       dinc - strain increment                                   -
    !   Outputs                                                       -
    !       sm   - mean stress (pressure)                             -
    !------------------------------------------------------------------
    implicit none
    real(8), intent(in):: elastic_strain(6)
    real(8):: dem, G

    dem = (elastic_strain(1) + elastic_strain(2) + elastic_strain(3)) / 3d0
    sm = K3*dem

    G = 0.5d0*G2    ! Shear modulus

    sd(1) = G2*(elastic_strain(1)-dem)    ! (Eq: 5.29)
    sd(2) = G2*(elastic_strain(2)-dem)
    sd(3) = G2*(elastic_strain(3)-dem)
    sd(4) = G*elastic_strain(4)
    sd(5) = G*elastic_strain(5)
    sd(6) = G*elastic_strain(6)

    end subroutine compute_elastic_stress_state


    subroutine lieupd()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Update energy for material models without EOS             -
    !------------------------------------------------------------------
    use ParticleData, only: EngInternal
    implicit none

    real:: vavg

    vavg = vol_ + vold
    iener = iener + 0.25d0*ieinc*vavg ! (Eq. 5.10)

    end subroutine lieupd


    subroutine hieupd()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Update energy for material models calling a EOS           -
    !-  References                                                    -
    !-      Section 5.1                                               -
    !------------------------------------------------------------------
    implicit none
    real:: vavg

    vavg = vol_ + vold

    ieinc = dinc(1)*(sold(1)+sd(1)) + dinc(2)*(sold(2)+sd(2)) +  &
        dinc(3)*(sold(3)+sd(3)) + dinc(4)*(sold(4)+sd(4)) +     &
        dinc(5)*(sold(5)+sd(5)) + dinc(6)*(sold(6)+sd(6))

    iener = iener + 0.25d0*ieinc*vavg + dvol*sm    ! (Eq: 5.15)

    specen = iener / vol0_

    end subroutine hieupd


    real(8) function EquivalentStress()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Calculate the equivalent stress sqrt(3*J2)                -
    !-  Input                                                         -
    !-      sd    - the deviatoric stress components                  -
    !-  Return                                                        -
    !-      EquivalentStress - the equivalent stress                  -
    !------------------------------------------------------------------
    use ParticleData
    implicit none

    real(8) :: J2  ! the second stress invariant

    J2 = 0.5d0*(sd(1)**2 + sd(2)**2 + sd(3)**2) + sd(4)**2 + &
        sd(5)**2 + sd(6)**2
    J2 = J2*3.0d0
    EquivalentStress = sqrt(J2)

    return
    end function EquivalentStress


    subroutine M3DM1()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Elastic material model                                    -
    !-  References                                                    -
    !-      Section 5.2.1                                             -
    !------------------------------------------------------------------
    implicit none

    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +   &
        sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)

    call elastic_devi()
    call elastic_p()

    seqv = EquivalentStress()

    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) +   &
        (sd(3)+sm)*dinc(3) +sd(4)*dinc(4) + sd(5)*dinc(5) + &
        sd(6)*dinc(6)     ! (part of Eq: 5.10)

    end subroutine M3DM1


    subroutine M3DM2()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Elastic-perfectly plastic material model                  -
    !-  References                                                    -
    !-      Section 5.2.4                                             -
    !------------------------------------------------------------------
    implicit none

    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +   &
        sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)

    call elastic_devi()

    seqv = EquivalentStress()

    if (seqv .GT. yield0_) then
        depeff = (seqv - sig_y_) / (1.5d0*G2)    ! (Eq: 5.93)
        epeff_ = epeff_ + depeff                 ! (Eq: 5.98)

        ratio = yield0_/seqv    ! (Eq: 5.100)

        sd(1) = sd(1)*ratio     ! (Eq: 5.101)
        sd(2) = sd(2)*ratio
        sd(3) = sd(3)*ratio
        sd(4) = sd(4)*ratio
        sd(5) = sd(5)*ratio
        sd(6) = sd(6)*ratio

        seqv = seqv*ratio       ! (Eq: 5.102)
    end if

    call elastic_p()

    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) + &
        (sd(3)+sm)*dinc(3) + sd(4)*dinc(4) +  &
        sd(5)*dinc(5) + sd(6)*dinc(6)   ! (part of Eq: 5.10)

    end subroutine M3DM2


    subroutine M3DM3()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Isotropic hardening plastic material model                -
    !-      LS-DYNA theorectical manual : Material Model 10           -
    !-  References                                                    -
    !-      Section 5.2.4                                             -
    !------------------------------------------------------------------
    implicit none

    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +   &
        sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)

    call elastic_devi()

    seqv = EquivalentStress()

    if (seqv .GT. sig_y_) then
        depeff = (seqv - sig_y_) / (1.5e0*G2 + PlaMod)    ! (Eq: 5.93)
        epeff_ = epeff_ + depeff                          ! (Eq: 5.98)

        sig_y_ = sig_y_ + PlaMod*depeff    ! (Eq: 5.99)
        ratio = sig_y_/seqv                ! (Eq: 5.100)

        sd(1) = sd(1)*ratio                ! (Eq: 5.101)
        sd(2) = sd(2)*ratio
        sd(3) = sd(3)*ratio
        sd(4) = sd(4)*ratio
        sd(5) = sd(5)*ratio
        sd(6) = sd(6)*ratio

        seqv = seqv*ratio
    end if

    call elastic_p()

    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) +  &
        (sd(3)+sm)*dinc(3) + sd(4)*dinc(4) + &
        sd(5)*dinc(5) + sd(6)*dinc(6)

    end subroutine M3DM3


    subroutine M3DM4(mat, DT, tmprt)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Johnson-Cook plastic material model (john)                -
    !-              (G.R. Johnson 1988)                               -
    !-      LS-DYNA theorectical manual : Material Model 15           -
    !-  References                                                    -
    !-      Section 5.2.5                                             -
    !------------------------------------------------------------------
    implicit none

    real(8), intent(in) :: DT, tmprt
    type(material), intent(in) :: mat
    real(8):: Bjc, njc, Cjc, mjc, epso, tstar

    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +  &
        sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)

    Bjc = mat%B_jc
    njc = mat%n_jc
    Cjc = mat%C_jc
    mjc = mat%m_jc
    epso = mat%epso

    epeff_ = epeff_ + 0.0001
    ! Note:
    !    PlaMod=Bjc*njc*(epeff_**(njc-1))*
    !           (1 + Cjc*log(depeff/epso/DT))*(1-1 - tstar**mjc)
    ! simplied as follow:
    PlaMod = Bjc*njc*(epeff_**(njc-1))
    epeff_ = epeff_ - 0.0001

    call elastic_devi()
    seqv = EquivalentStress()

    if (seqv .GT. sig_y_) then
        depeff = (seqv - sig_y_) / (1.5e0*G2 + PlaMod)
        epeff_ = epeff_ + depeff

        tstar = (tmprt-mat%roomt)/(mat%melt-mat%roomt)

        ! (Eq: 5.107)
        sig_y_ = (yield0_ + Bjc*(epeff_**njc)) * &
            (1 + Cjc*log(depeff/epso/DT)) * (1 - tstar**mjc)

        ratio = sig_y_/seqv    ! (Eq: 5.100)

        sd(1) = sd(1)*ratio    ! (Eq: 5.101)
        sd(2) = sd(2)*ratio
        sd(3) = sd(3)*ratio
        sd(4) = sd(4)*ratio
        sd(5) = sd(5)*ratio
        sd(6) = sd(6)*ratio

        seqv = seqv*ratio      ! (Eq: 5.102)
    end if

    call elastic_p()

    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) + &
        (sd(3)+sm)*dinc(3) +  sd(4)*dinc(4) + &
        sd(5)*dinc(5) + sd(6)*dinc(6)  ! (part of Eq: 5.10)

    end subroutine M3DM4


    subroutine M3DM5(mat, DT)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Simplified Johnson-Cook plastic material model (sjc)      -
    !-                 Ignoring the temperature effects in            -
    !-                 Johnson-Cook model                             -
    !-      LS-DYNA theorectical manual : Material Model 98           -
    !-      should be used with EOS                                   -
    !-  References                                                    -
    !-      Section 5.2.4                                             -
    !------------------------------------------------------------------
    implicit none

    real(8), intent(in) :: DT
    type(material), intent(in) :: mat
    real(8):: Bjc, njc, Cjc, epso

    Bjc = mat%B_jc
    njc = mat%n_jc
    Cjc = mat%C_jc
    epso = mat%epso

    epeff_ = epeff_ + 0.0001
    PlaMod = Bjc*njc*(epeff_**(njc-1))
    epeff_ = epeff_ - 0.0001

    call elastic_devi()

    seqv = EquivalentStress()

    if (seqv .GT. sig_y_) then
        depeff = (seqv - sig_y_) / (1.5e0*G2 + PlaMod)
        epeff_ = epeff_ + depeff

        sig_y_ = (yield0_ + Bjc*(epeff_**njc)) * &
            (1 + Cjc*log(depeff/epso/DT))

        ratio = sig_y_/seqv

        sd(1) = sd(1)*ratio
        sd(2) = sd(2)*ratio
        sd(3) = sd(3)*ratio
        sd(4) = sd(4)*ratio
        sd(5) = sd(5)*ratio
        sd(6) = sd(6)*ratio

        seqv = seqv*ratio
    end if

    end subroutine M3DM5


    subroutine M3DM6(ltim,ctim,dvelo)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      High Explosive burn material                              -
    !-      an Equation of State must be defined                      -
    !-      LS-DYNA theorectical manual : Material Model 8            -
    !-  References                                                    -
    !-      Section 5.2.13                                            -
    !------------------------------------------------------------------
    implicit none
    real(8), intent(in) :: ltim, ctim, dvelo

    if (ctim.gt.ltim) then
        bfac = 1.0
    else
        bfac = 0.0
    end if

    end subroutine M3DM6


    subroutine M3DM7()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Null material model                                       -
    !------------------------------------------------------------------
    implicit none

    sd = 0.0

    end subroutine M3DM7


    subroutine M3DM8(mat, DT, tmprt)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Johnson-Cook plastic material model (john)                -
    !-      (G.R. Johnson 1988)                                       -
    !-      LS-DYNA theorectical manual : Material Model 15           -
    !-      should be used with EOS                                   -
    !-  References                                                    -
    !-      Section 5.2.5                                             -
    !------------------------------------------------------------------
    implicit none

    real(8), intent(in) :: DT, tmprt
    type(material), intent(in) :: mat
    real(8):: Bjc, njc, Cjc, mjc, epso, tstar

    Bjc = mat%B_jc
    njc = mat%n_jc
    Cjc = mat%C_jc
    mjc = mat%m_jc
    epso = mat%epso

    epeff_ = epeff_ + 0.0001
    PlaMod = Bjc*njc*(epeff_**(njc-1))
    epeff_ = epeff_ - 0.0001


    call elastic_devi()
    seqv = EquivalentStress()

    if (seqv .GT. sig_y_) then
        depeff = (seqv - sig_y_) / (1.5e0*G2 + PlaMod)
        epeff_ = epeff_ + depeff

        tstar = (tmprt-mat%roomt)/(mat%melt-mat%roomt)

        sig_y_ = (yield0_ + Bjc*(epeff_**njc)) * &
            (1 + Cjc*log(depeff/epso/DT)) * (1 - tstar**mjc)

        ratio = sig_y_/seqv

        sd(1) = sd(1)*ratio
        sd(2) = sd(2)*ratio
        sd(3) = sd(3)*ratio
        sd(4) = sd(4)*ratio
        sd(5) = sd(5)*ratio
        sd(6) = sd(6)*ratio

        seqv = seqv*ratio
    end if

    end subroutine M3DM8


    subroutine M3DM10(mat, DT)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      perfect elastic-plastic Drucker-Prager model for soil     -
    !-      LS-DYNA theorectical manual : Material Model 193          -
    !-      should be used without EOS                                -
    !-  References                                                    -
    !-      Section 5.2.6                                             -
    !------------------------------------------------------------------
    implicit none

    real(8), intent(in) :: DT
    type(material), intent(in) :: mat
    real(8):: qfai, kfai, qpsi, tenf, Gmod, Kmod,tenf_max
    real(8):: J2, Tau, dpFi, dpsig, dlamd, newTau, dem, G
    real(8):: dp_hfai, Taup,alphap, rprops(100),dsig(6),a,b,c
    integer :: iplas
    real(8), parameter :: sqrt3 = 1.73205080756888d0

    rprops = mat%mprops ! [density, young, poisson, q_fai, k_fai, q_psi, ten_f]
    qfai = rprops(4)
    kfai = rprops(5)
    qpsi = rprops(6)
    tenf = rprops(7)
    Gmod = young_ / (2d0* (1d0 + poisson_))
    Kmod =  young_ / (3d0*(1d0 - 2d0*poisson_))

    ! internal energy incremental
    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +  &
        sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)
  

	!KW=(3*(0.49-0.35)/((1-2.0d0*0.49)*(1+0.35)))*Kmod
	
    ! Give the tension stress value
    !* --- set fTension to cone apex if larger than apex --- *
    if (qfai == 0.0d0 )then
        tenf = 0.0d0
    else
        tenf_max =  kfai / qfai
        tenf = min(tenf,tenf_max)
    end if

    ! --- trial elastic stresses --
    call elastic_devi()
    call elastic_p()

    iplas = 0         ! elastic calculation
    J2 = 0.5d0*(sd(1)**2 + sd(2)**2 + sd(3)**2) + sd(4)**2 + &
        sd(5)**2 + sd(6)**2
    Tau = sqrt(J2)
    seqv = Tau*sqrt3
    dpFi = Tau + qfai*sm - kfai   ! D-P yield surface
    dpsig = sm - tenf             ! the spherical stress difference

    if (dpsig < 0.0d0) then
        if (dpFi >0.0d0) then
            iplas = 1 ! shear plastic flow
            ! plastic flow coefficient
            dlamd = dpFi/(Gmod + Kmod*qfai*qpsi)
            ! correct spherical stress
            sm = sm - Kmod*qpsi*dlamd
            ! correct shear stress
            newTau = kfai - qfai*sm
            ratio = newTau/Tau

            ! correct deviatoric stress
            sd(1) = sd(1)*ratio
            sd(2) = sd(2)*ratio
            sd(3) = sd(3)*ratio
            sd(4) = sd(4)*ratio
            sd(5) = sd(5)*ratio
            sd(6) = sd(6)*ratio

            seqv =  seqv*ratio     !correct the mises stress
            
            ! calculate the effective plastic strain
            depeff = dlamd*sqrt(1.0d0/3.0d0 + (2.0d0/9.0d0)*(qpsi**2))
            epeff_ = epeff_ + depeff
            

        end if
    else   !(dpsig >= 0.0)then
        alphap = sqrt(1d0 + qfai**2) - qfai
        Taup = kfai - qfai*tenf
        dp_hfai = Tau - Taup - alphap * dpsig

        if(dp_hfai > 0.0d0) then
            iplas = 1 ! shear plastic flow
            ! plastic flow coefficient
            dlamd = dpFi/(Gmod + Kmod*qfai*qpsi)
            ! correct spherical stress
            sm = sm - Kmod*qpsi*dlamd
            ! correct shear stress
            newTau = kfai - qfai*sm
            ratio = newTau/Tau

            ! correct deviatoric stress
            sd(1) = sd(1)*ratio
            sd(2) = sd(2)*ratio
            sd(3) = sd(3)*ratio
            sd(4) = sd(4)*ratio
            sd(5) = sd(5)*ratio
            sd(6) = sd(6)*ratio

            seqv =  seqv*ratio     !correct the mises stress
            ! calculate the effective plastic strain
            depeff = dlamd*sqrt(1.0d0/3.0d0 + (2.0d0/9.0d0)*(qpsi**2))
            epeff_ = epeff_ + depeff
             
        else
            iplas = 2 ! tension plastic flow
            dlamd = (sm - tenf)/Kmod
            sm = tenf
			
            depeff = dlamd*(1.0d0/3.0d0)*sqrt(2.0d0)
            epeff_ = epeff_ + depeff
             


        end if
    end if

    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) + &
        (sd(3)+sm)*dinc(3) + sd(4)*dinc(4) + &
        sd(5)*dinc(5) + sd(6)*dinc(6)

    end subroutine M3DM10


    subroutine M3DM11(mat, DT, trial_elastic_strain, elastic_strain)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      perfect elastic-plastic Drucker-Prager model for soil     -
    !-      LS-DYNA theorectical manual : Material Model 193          -
    !-      should be used without EOS                                -
    !-  References                                                    -
    !-      Section 5.2.6                                             -
    !------------------------------------------------------------------
    implicit none

    real(8), intent(in) :: DT, trial_elastic_strain(6)
    type(material), intent(in) :: mat
    real(8), intent(inout):: elastic_strain(6)
    real(8):: qfai, kfai, qpsi, tenf, Gmod, Kmod,tenf_max
    real(8):: J2, Tau, dpFi, dpsig, dlamd, newTau
    real(8):: dp_hfai, Taup, alphap, rprops(100), plas_strain_inc_vol
    real(8):: xi, cohesion, epbarn, plas_strain_inc(6), elastic_strain_inc(6)
    integer:: iplas, iphard = 9, nhard
    real(8), parameter :: sqrt3 = 1.73205080756888d0

    rprops = mat%mprops ! [density, young, poisson, q_fai, q_psi, ten_f, nhard, ep1, k_fai1, ...]
    qfai = rprops(4)
    xi = rprops(5)
    qpsi = rprops(6)
    tenf = rprops(7)
    nhard = int(rprops(8))
    cohesion = plfun(epeff_, nhard, rprops(iphard:))
    kfai = cohesion*xi
    Gmod = young_ / (2d0* (1d0 + poisson_))
    Kmod =  young_ / (3d0*(1d0 - 2d0*poisson_))

    ! internal energy incremental
    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +  &
        sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)

    ! Give the tension stress value
    !* --- set fTension to cone apex if larger than apex --- *
    if (qfai == 0.0d0 )then
        tenf = 0.0d0
    else
        tenf_max =  kfai / qfai
        tenf = min(tenf,tenf_max)
    end if

    ! trial elastic strain increment
    elastic_strain_inc = dinc
    plas_strain_inc = 0d0

    ! --- trial elastic stresses ---
    call compute_elastic_stress_state(trial_elastic_strain)

    iplas = 0         ! elastic calculation
    J2 = 0.5d0*(sd(1)**2 + sd(2)**2 + sd(3)**2) + sd(4)**2 + &
        sd(5)**2 + sd(6)**2
    Tau = sqrt(J2)
    seqv = Tau*sqrt3
    dpFi = Tau + qfai*sm - kfai   ! D-P yield surface
    dpsig = sm - tenf             ! the spherical stress difference

    if (dpsig < 0.0d0) then
        if (dpFi >0.0d0) then
            iplas = 1 ! shear plastic flow
            ! plastic flow coefficient
            dlamd = dpFi/(Gmod + Kmod*qfai*qpsi)
            ! correct spherical stress
            sm = sm - Kmod*qpsi*dlamd
            ! correct shear stress
            newTau = kfai - qfai*sm
            ratio = newTau/Tau

            ! correct deviatoric stress
            sd(1) = sd(1)*ratio
            sd(2) = sd(2)*ratio
            sd(3) = sd(3)*ratio
            sd(4) = sd(4)*ratio
            sd(5) = sd(5)*ratio
            sd(6) = sd(6)*ratio

            seqv =  seqv*ratio     !correct the mises stress

            ! calculate the effective plastic strain
            depeff = dlamd*sqrt(1.0d0/3.0d0 + (2.0d0/9.0d0)*(qpsi**2))
            epeff_ = epeff_ + depeff
            ! calculate elastic strain increment
            plas_strain_inc = dlamd*sd/newTau/2d0
            plas_strain_inc_vol = dlamd*qpsi/3d0
            plas_strain_inc(1:3) = plas_strain_inc(1:3) + plas_strain_inc_vol
            plas_strain_inc(4:6) = plas_strain_inc(4:6)*2 ! engineering strain
        end if
    else   !(dpsig >= 0.0)then
        alphap = sqrt(1d0 + qfai**2) - qfai
        Taup = kfai - qfai*tenf
        dp_hfai = Tau - Taup - alphap * dpsig

        if(dp_hfai > 0.0d0) then
            iplas = 1 ! shear plastic flow
            ! plastic flow coefficient
            dlamd = dpFi/(Gmod + Kmod*qfai*qpsi)
            ! correct spherical stress
            sm = sm - Kmod*qpsi*dlamd
            ! correct shear stress
            newTau = kfai - qfai*sm
            ratio = newTau/Tau

            ! correct deviatoric stress
            sd(1) = sd(1)*ratio
            sd(2) = sd(2)*ratio
            sd(3) = sd(3)*ratio
            sd(4) = sd(4)*ratio
            sd(5) = sd(5)*ratio
            sd(6) = sd(6)*ratio

            seqv =  seqv*ratio     !correct the mises stress
            ! calculate the effective plastic strain
            depeff = dlamd*sqrt(1.0d0/3.0d0 + (2.0d0/9.0d0)*(qpsi**2))
            epeff_ = epeff_ + depeff
            ! calculate elastic strain increment
            plas_strain_inc = dlamd*sd/newTau/2
            plas_strain_inc_vol = dlamd*qpsi/3d0
            plas_strain_inc(1:3) = plas_strain_inc(1:3) + plas_strain_inc_vol
            plas_strain_inc(4:6) = plas_strain_inc(4:6)*2 ! engineering strain
        else
            iplas = 2 ! tension plastic flow
            dlamd = (sm - tenf)/Kmod
            sm = tenf
            depeff = dlamd*(1.0d0/3.0d0)*sqrt(2.0d0)
            epeff_ = epeff_ + depeff
            ! calculate elastic strain increment
            plas_strain_inc_vol = dlamd/3d0
            plas_strain_inc = 0d0
            plas_strain_inc(1:3) = plas_strain_inc_vol
        end if
    end if
    
    elastic_strain = trial_elastic_strain - plas_strain_inc
    
    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) + &
        (sd(3)+sm)*dinc(3) + sd(4)*dinc(4) + &
        sd(5)*dinc(5) + sd(6)*dinc(6)

    end subroutine M3DM11
    
    
    subroutine M3DM12(mat, DT)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      perfect elastic-plastic Mohr-Coulomb model for soil       -
    !-      FLAC theorectical manual (Verison 5.00), page 715         -
    !-      should be used without EOS                                -
    !-  References                                                    -
    !-      FLAC theorectical manual (Verison 5.00), page 715         -
    !------------------------------------------------------------------
    implicit none

    real(8), intent(in) :: DT
    type(material), intent(inout) :: mat
	integer:: iplas, ii, jj, mm, i, nhard, iphard=8
    real(8):: rprops(100), cohesion, sinphi, tanphi, sinpsi, tenf, tenf_max
    real(8):: alpha1, alpha2, nphi, npusi, alphap, sigp, j3, theta
	real(8):: sig_mat(3,3)=0d0, sig_eig_trial(3,3), sig_eigv_trial(3,3)
	real(8):: pstrs(3), pstrs_new(3), pstrs_trial(3), pstrs_delta(3)
	real(8):: fa, fb, c1, c2, c3, c4, fnew, delta, lambda_a, lambda_b
	real(8):: dlamd, fs_mohr, ft_mohr, h_mohr, temp, epsp(3)
	real(8), parameter:: small = 1d-6

	! Material parameters
    rprops = mat%mprops ![density, young, poisson, nphi, npusi, tenf, nhard, ...]
    nphi = rprops(4) ![density, young, poisson, nphi, npusi, tenf, nhard, ...]
    npusi = rprops(5)
    tenf = rprops(6)
    nhard = int(rprops(7))
    cohesion = plfun(epeff_, nhard, rprops(iphard:))
    if (nphi .eq. 1d0) then
        tenf = 0d0
    else
        sinphi = (nphi - 1d0)/(nphi + 1d0)
        tanphi = sinphi/sqrt(1 - sinphi**2)
        tenf_max =  cohesion / tanphi;
        tenf = min(tenf,tenf_max);
    end if
    
	! variables for calculation
	alpha1 = K3/3d0 + G2*2d0/3d0
    alpha2 = K3/3d0 - G2*1d0/3d0
	alphap = sqrt(1d0 + nphi**2) + nphi
    sigp = tenf*nphi - 2d0*cohesion*sqrt(nphi)
    
    ! internal energy incremental
    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +  &
        sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)
    
    ! --- trial elastic stresses ---
    call elastic_devi()
    call elastic_p()
    
    ! transform stress vector into stress matrix
    do i = 1,3
	    sig_mat(i,i) = sm + sd(i)
	end do
	sig_mat(1,2) = sd(6)
	sig_mat(2,1) = sd(6)
	sig_mat(1,3) = sd(5)
	sig_mat(3,1) = sd(5)
	sig_mat(2,3) = sd(4)
	sig_mat(3,2) = sd(4)
    
    ! compute eigenvalues and eigenvectors of trial stresses
    call jacob(sig_mat,sig_eig_trial,sig_eigv_trial,3)
    
    ! identify minimum (pstrs1) and maximum (pstrs3) principal stresses
    ii = 1; jj = 1;
    pstrs_trial(1) = sig_eig_trial(ii,ii)
    pstrs_trial(3) = sig_eig_trial(jj,jj)
    do i = 2, 3
        if (sig_eig_trial(i,i)<pstrs_trial(1)) then
            ii = i
            pstrs_trial(1) = sig_eig_trial(ii,ii)
        end if
        if (sig_eig_trial(i,i)>=pstrs_trial(3)) then
            jj = i
            pstrs_trial(3) = sig_eig_trial(jj,jj)
        end if
    end do
    if (ii/=1 .and. jj/=1) mm = 1
    if (ii/=2 .and. jj/=2) mm = 2
    if (ii/=3 .and. jj/=3) mm = 3
    pstrs_trial(2) = sig_eig_trial(mm,mm)

    ! yield conditions
    fs_mohr = pstrs_trial(1) - pstrs_trial(3)*nphi + 2d0*cohesion*sqrt(nphi)
    ft_mohr = tenf - pstrs_trial(3)
    if (ft_mohr .ge. 0d0) then
		if (fs_mohr .gt. 0d0) then
		    ! elastic state, do nothing, ft>=0, fs>0
			iplas = 0
            pstrs_new(1) = pstrs_trial(1)
            pstrs_new(2) = pstrs_trial(2)
            pstrs_new(3) = pstrs_trial(3)
            fnew = fs_mohr
			depeff = 0d0
		else 
            ! plastic shear yielding, ft>=0, fs<=0
            iplas = 1
            dlamd = fs_mohr/((alpha1 - alpha2*npusi) - (alpha2 - alpha1*npusi)*nphi)
            pstrs_new(1) = pstrs_trial(1) - dlamd*(alpha1 - alpha2*npusi)
            pstrs_new(2) = pstrs_trial(2) - dlamd*alpha2*(1 - npusi)
            pstrs_new(3) = pstrs_trial(3) - dlamd*(alpha2 - alpha1*npusi)
            fnew = pstrs_new(1) - pstrs_new(3)*nphi + 2d0*cohesion*sqrt(nphi)

            ! calculate the effective plastic strain
            depeff = 2.0d0/3.0d0*abs(dlamd)*sqrt(1.0d0 + npusi + npusi**2)
		end if
    else
		!h_mohr = pstrs_trial(3) - tenf + alphap*(pstrs_trial(1) - sigp)
        h_mohr = pstrs_trial(3) + pstrs_trial(1)*nphi + 2*cohesion*nphi*sqrt(nphi) - (1+ nphi**2)*tenf
		if (h_mohr .lt. 0d0) then
			! plastic shear yielding, ft<0, h<0
		    iplas = 1
            dlamd = fs_mohr/((alpha1 - alpha2*npusi) - (alpha2 - alpha1*npusi)*nphi)
            pstrs_new(1) = pstrs_trial(1) - dlamd*(alpha1 - alpha2*npusi)
            pstrs_new(2) = pstrs_trial(2) - dlamd*alpha2*(1 - npusi)
            pstrs_new(3) = pstrs_trial(3) - dlamd*(alpha2 - alpha1*npusi)
            fnew = pstrs_new(1) - pstrs_new(3)*nphi + 2d0*cohesion*sqrt(nphi)

            ! calculate the effective plastic strain
            depeff = 2.0d0/3.0d0*abs(dlamd)*sqrt(1.0d0 + npusi + npusi**2)
		else
            ! tensile plastic yileding, ft<0, h>=0
            iplas = 2
            dlamd = ft_mohr/alpha1
            pstrs_new(1) = pstrs_trial(1) + dlamd*alpha2
            pstrs_new(2) = pstrs_trial(2) + dlamd*alpha2
            pstrs_new(3) = pstrs_trial(3) + dlamd*alpha1
            fnew = tenf - pstrs_new(3)
            
            ! calculate the effective plastic strain
            depeff = abs(dlamd)*2.0d0/3.0d0
        end if
    end if
    
    ! check whether multi-surface yielding corrector is needed
    if (iplas .eq. 1) then
        delta = maxval(abs(pstrs_new))*small
        if (pstrs_new(1) .le. pstrs_new(2) + delta .and. pstrs_new(2) .le. pstrs_new(3) + delta) then
            ! do nothing
		elseif (pstrs_new(1) .gt. pstrs_new(2) .and. pstrs_new(1) .le. pstrs_new(3) + delta) then
			iplas = 3;
		elseif (pstrs_new(1) .le. pstrs_new(3) + delta .and. pstrs_new(2) .gt. pstrs_new(3)) then
			iplas = 4;
		else
			iplas = 5;
        end if
    end if
    
	if (iplas .eq. 3) then
		! two yield surfaces considering sp1 <= sp2 <= sp3 and sp2 < sp1 <= sp3
		fa = fs_mohr
		fb = pstrs_trial(2) - pstrs_trial(3)*nphi + 2d0*cohesion*sqrt(nphi)
		c1 = (alpha1 - alpha2*npusi) - (alpha2 - alpha1*npusi)*nphi
		c2 = (alpha2 - alpha2*npusi) - (alpha2 - alpha1*npusi)*nphi
        c3 = c2; c4 = c1;
		lambda_a = (c4*fa - c2*fb)/(c1*c4 - c2*c3);
        lambda_b = (c1*fb - c3*fa)/(c1*c4 - c2*c3);
		epsp(1) = lambda_a
		epsp(2) = lambda_b
		epsp(3) = -(lambda_a + lambda_b)*npusi
		pstrs_new(1) = pstrs_trial(1) - alpha1*epsp(1) - alpha2*epsp(2) - alpha2*epsp(3)
		pstrs_new(2) = pstrs_trial(2) - alpha1*epsp(2) - alpha2*epsp(1) - alpha2*epsp(3)
		pstrs_new(3) = pstrs_trial(3) - alpha1*epsp(3) - alpha2*epsp(1) - alpha2*epsp(2)
		depeff = sqrt(2d0/9d0*((epsp(1) - epsp(2))**2 + (epsp(1) - epsp(3))**2 + (epsp(2) - epsp(3))**2))
        
	elseif (iplas .eq. 4) then
		! two yield surfaces considering sp1 <= sp2 <= sp3 and sp1 <= sp3 < sp2
		fa = fs_mohr;
		fb = pstrs_trial(1) - pstrs_trial(2)*nphi + 2*cohesion*sqrt(nphi);
		c1 = (alpha1 - alpha2*npusi) - (alpha2 - alpha1*npusi)*nphi;
		c2 = (alpha1 - alpha2*npusi) - (alpha2 - alpha2*npusi)*nphi;
        c3 = c2; c4 = c1;
		lambda_a = (c4*fa - c2*fb)/(c1*c4 - c2*c3);
        lambda_b = (c1*fb - c3*fa)/(c1*c4 - c2*c3);
		epsp(1) = lambda_a + lambda_b; 
		epsp(2) = -lambda_b*npusi; 
		epsp(3) = -lambda_a*npusi;
		pstrs_new(1) = pstrs_trial(1) - alpha1*epsp(1) - alpha2*epsp(2) - alpha2*epsp(3);
		pstrs_new(2) = pstrs_trial(2) - alpha1*epsp(2) - alpha2*epsp(1) - alpha2*epsp(3);
		pstrs_new(3) = pstrs_trial(3) - alpha1*epsp(3) - alpha2*epsp(1) - alpha2*epsp(2);
		depeff = sqrt(2d0/9d0*((epsp(1) - epsp(2))**2 + (epsp(1) - epsp(3))**2 + (epsp(2) - epsp(3))**2));

	elseif (iplas .eq. 5) then
		! return to apex
		if (nphi/=1d0) then
			pstrs_new(1) = 2d0*cohesion*sqrt(nphi)/(nphi - 1d0);
			pstrs_new(2) = pstrs_new(1);
			pstrs_new(3) = pstrs_new(1);
            pstrs_delta(1) = pstrs_trial(1) - pstrs_new(1);
            pstrs_delta(2) = pstrs_trial(2) - pstrs_new(2);
            pstrs_delta(3) = pstrs_trial(3) - pstrs_new(3);
            temp = (alpha1**2 + alpha1*alpha2 - 2d0*alpha2**2);
            epsp(1) = (pstrs_delta(1)*(alpha1 + alpha2))/temp - (alpha2*pstrs_delta(2))/temp - (alpha2*pstrs_delta(3))/temp;
            epsp(2) = (pstrs_delta(2)*(alpha1 + alpha2))/temp - (alpha2*pstrs_delta(1))/temp - (alpha2*pstrs_delta(3))/temp;
            epsp(3) = (pstrs_delta(3)*(alpha1 + alpha2))/temp - (alpha2*pstrs_delta(1))/temp - (alpha2*pstrs_delta(2))/temp;
            depeff = sqrt(2d0/9d0*((epsp(1) - epsp(2))**2 + (epsp(1) - epsp(3))**2 + (epsp(2) - epsp(3))**2));
		else
			pstrs_new(1) = 0d0;
			pstrs_new(2) = 0d0;
			pstrs_new(3) = 0d0;
			depeff = 0d0
        end if
	end if

	! update stress using corrected principal stresses
	sig_eig_trial(ii,ii) = pstrs_new(1)
	sig_eig_trial(jj,jj) = pstrs_new(3)
	sig_eig_trial(mm,mm) = pstrs_new(2)
	sig_mat = MATMUL(MATMUL(sig_eigv_trial,sig_eig_trial), transpose(sig_eigv_trial))
	do i = 1, 3
		sig(i) = sig_mat(i,i)
	end do
	sig(4) = sig_mat(2,3)
	sig(5) = sig_mat(1,3)
	sig(6) = sig_mat(1,2)
	
	! mean stress and deviate stresses
	sm = (sig(1)+sig(2)+sig(3))/3d0
	sd(1) = sig(1) - sm
	sd(2) = sig(2) - sm
	sd(3) = sig(3) - sm
	sd(4) = sig(4)
	sd(5) = sig(5)
	sd(6) = sig(6)
    epeff_ = epeff_ + depeff
	
    ! equivalent stress
    seqv = EquivalentStress()

    ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) + &
        (sd(3)+sm)*dinc(3) + sd(4)*dinc(4) + &
        sd(5)*dinc(5) + sd(6)*dinc(6)

    end subroutine M3DM12


    subroutine M3DM13(mat, DT, trial_elastic_strain, elastic_strain_state)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Mohr-Coulomb model for soil with piece-wise hardening     -
    !-      should be used without EOS                                -
    !-  References                                                    -
    !-      de Souza Neto et al., 2008. Computational Methods for     -
    !-      Plasticity: Theory and Applications                       -
    !------------------------------------------------------------------
    implicit none

    type(material), intent(in) :: mat
    real(8), intent(in) :: DT, trial_elastic_strain(6)
    real(8), intent(inout) :: elastic_strain_state(6)
    real(8):: J2, tau, dgam(2), sqrt3
    real(8):: rprops(100), epbarn
    logical:: lalgva(5)
    data sqrt3 /1.73205080756888d0/

    rprops = mat%mprops ![density, young, poisson, sinphi, cosphi, sinpsi, nhard]
    epbarn = epeff_

    ! internal energy incremental
    ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +  &
        sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)

    ! call Mohr-Coulomb return-mapping
    call sumc(dgam, lalgva, rprops, elastic_strain_state, trial_elastic_strain, epbarn, sig)

    ! accumulated plastic strain
    epeff_ = epbarn 
    
    ! deviatoric stress and mean stress
    sd = sig
    sm = (sig(1) + sig(2) + sig(3))/3.0d0
    sd(1) = sd(1) - sm
    sd(2) = sd(2) - sm
    sd(3) = sd(3) - sm

    ! von mises effective stress
    J2 = 0.5d0*(sd(1)**2 + sd(2)**2 + sd(3)**2) + sd(4)**2 + &
        sd(5)**2 + sd(6)**2
    tau = sqrt(J2)
    seqv = tau*sqrt3

    ! internal energy incremental
    ieinc = ieinc + sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +  &
        sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)

     end subroutine M3DM13

     subroutine M3DM14(mat, DT, plastic_strain)
	!-----------------------------------------------------------------------
	!Purpose                                                               -
    !     Modified Cam-clay model for granular materials should be used    -
	!     without EOS                                                      -  
    !References                                                            -  
    !    (1) D.M. Potts & A. Gens(1985) A critical assesment of methods    -
	!        for correcting drift from the yield surface in elasto-plastic - 
	!		 finite element analysis, IJNAMG                               -
    !    (2) R Lorenzo, R. P. Cunha & Manoel (2013) Material point method  -
	!        for geotechnical problems involving large deformation         - 
    !-----------------------------------------------------------------------
        implicit none

         real(8), intent(in) :: DT
         type(material), intent(inout) :: mat
	     real(8), intent(inout):: plastic_strain(6)
	     real(8):: rprops(100), phi,E, TOLF, TOLG,pcn, OCR
         real(8):: v,Pstr(6),po,K,T(3,3),DeltaE(3,3),ps 
         real(8):: lambda,kappa,M,nu,Lalambda,LaG,ELMOD(6,6)
         real(8):: p,q,cElas(3,3,3,3),xi,cEplas(3,3,3,3)
         real(8):: IdenT(3,3,3,3),IdenT1(3,3,3,3), Uni(3,3) 
         real(8):: NormTrDev, NormTrStr, TrF, Trp, Trq, theta
         real(8):: TrStr(3,3), TrDev(3,3), pc0, J2, J3
         real(8):: dphi, G, F, dGdpc, dFdpc, dTrStr(3,3), e0
         real(8):: ddphide(3,3), dirn(3,3), theta1, alph
         real(8):: dEplas3(3,3), dEplas6(6), TDev(3,3), NormTDev, NormdEplas 
         Integer:: i, l, n, j, maxiter
	     real(8), parameter:: PI=4.D0*DATAN(1.D0)
	     ! Material parameters
         rprops = mat%mprops ![density, lamda, kapa, poisson, phi, OCR, e0, pc]
	     lambda = rprops(2)
	     kappa =  rprops(3)
	     nu = rprops(4)
	     phi = rprops(5)
	     e0 = rprops(6)
	
	     M = 6.0d0 * sin(phi * PI / 180.0d0) / (3.0d0 - sin(phi * PI / 180.0d0));
	     v=1.0d0+evoid
	     ! internal energy incremental
         ieinc = sig(1)*dinc(1) + sig(2)*dinc(2) + sig(3)*dinc(3) +  &
                sig(4)*dinc(4) + sig(5)*dinc(5) + sig(6)*dinc(6)
	     !Plastic strains as state variable  
         Pstr(1)=plastic_strain(1)
         Pstr(2)=plastic_strain(2) 
	     Pstr(3)=plastic_strain(3) 
	     Pstr(4)=plastic_strain(6) 
	     Pstr(5)=plastic_strain(5)
	     Pstr(6)=plastic_strain(4)
         !Overconsolidated mean stress state variable 
         pc=abs(pc)
         CALL Initial( sm, sd, dinc, T, DeltaE)
        theta=(1.0d0+evoid)/(lambda-kappa)
        !Mean stress 
        p=(T(1,1)+T(2,2)+T(3,3))/3.0d0
        if (abs(p)<1000) then 
             p=1000
        endif
        if (abs(p)>abs(pc)) then
            pc=p
        endif
		 if (abs(pc)<1000) then 
             pc=1000
         endif
		
        !Material constants 
         !K=K3/3.0d0
		K= abs(v*p/kappa)
        if (abs(K)<1000.0d0) then 
            K=1000.0d0
        endif

        E=3*K*(1.0d0-2.0d0*nu) 
        !Lamé constants
        LaG=E/(2.0d0*(1.0d0+nu))
        LaLambda=(nu*E)/((1.0d0+nu)*(1.0d0-2.0d0*nu))

        !Elastic modulus and identity fourth order tensor 
         CALL TenIsoUni(IdenT)
         CALL TenUni1(IdenT1)
         CALL TenUni(Uni)
         CALL Emod(LaLambda,LaG,IdenT1,IdenT,cElas)

        ! Delta Trial Stress
         CALL dCon(cElas,DeltaE,dTrStr)
	     !	Trial Stress
         TrStr=T+dTrStr
         CALL normS(TrStr,NormTrStr) 
	     ! Trial p
         Trp=(TrStr(1,1)+TrStr(2,2)+TrStr(3,3))/3.0d0
         !Deviator Stress T*
         CALL dev(TrStr,TrDev,Uni)
	     !	Deviator Stress norm
         CALL normS(TrDev,NormTrDev)
	     !	Trial q
         Trq=sqrt(3.0d0/2.0d0)*NormTrDev
	     !	Trial F (Yielding surface)
         TrF=(Trq**2.0d0/M**2.0d0)+(abs(Trp))*(abs(Trp)-abs(pc))
         !Verificando elasticidad
         if (TrF<0.0d0) then
     
         !Elastic behavior
          do i = 1, 3
		  sig(i) = TrStr(i,i)
	      end do
	      sig(4) = TrStr(2,3)
	      sig(5) = TrStr(1,3)
	      sig(6) = TrStr(1,2)
	    ! mean stress and deviate stresses
	      sm = (sig(1)+sig(2)+sig(3))/3.0d0
	      sd(1) = sig(1) - sm
	      sd(2) = sig(2) - sm
	      sd(3) = sig(3) - sm
	      sd(4) = sig(4)
	      sd(5) = sig(5)
	      sd(6) = sig(6)
	
        ! Call Solution(NTENS, NDI, NSHR, TrStr, STRESS, cElas, DDSDDE)
        !Void ratio change
	     evoid=evoid+(1.0d0+evoid)*(DeltaE(1,1)+DeltaE(2,2)+DeltaE(3,3)) 

	    !Saving State variables 
		 seqv=EquivalentStress()
         depeff=0.0d0
        return
     else
      
       ! Plastic behavior 
	      pcn=abs(pc) 
	      dphi=0 
	      maxiter=100
    !	Tolerance is defined using a constant value and a variable
    !	that depends on functions F and G. In this case F depends on q and 
    !	q depends on the deviator stress tensor
         TolF=1.0E-5
    !	G depends on the overconsolidated stress 
         TolG=1.0E-5
         do i=0, maxiter, 1 
	     do j=0, maxiter, 1
         G=pcn*exp(theta*dphi*((2.0d0*abs(Trp)-abs(pc))/(1.0d0+2.0d0*K*dphi)))-abs(pc)
         if(abs(G)<TolG) then
	     q=Trq/(1.0d0+6.0d0*LaG*(dphi/M**2))
	     p=(abs(Trp)+dphi*K*abs(pc))/(1.0d0+2.0d0*dphi*K)
	     F=q**2.0d0/M**2.0d0+(abs(p))*(abs(p)-abs(pc))
	     exit
         else
         dGdpc=-(pcn*dphi*theta/(1.0d0+2.0d0*dphi*K))*exp(((2.0d0*abs(p)-abs(pc))*theta*dphi)/(1.0d0+2.0d0*dphi*K))-1.0d0
         pc=abs(pc)-G/dGdpc 
	     endif
         enddo
         if (abs(F)<TolF) then 
	     exit
         else
    !	 Derivative of F
         dFdpc=-K*((2.0d0*abs(p)-abs(pc))**2.0d0)/(1.0d0+dphi*(2.0d0*k+ &
           theta*abs(pc)))-(2.0d0*(q**2))/(M**2.0d0*(dphi+M**2.0d0/(6.0d0* &
           LaG)))-theta*abs(pc)*(((abs(p))*(2.0d0*abs(p)-abs(pc))/(1.0d0+ &
           (2.0d0*K+theta*abs(pc))*dphi)))
        dphi=dphi-F/dFdpc
   
        endif 
        enddo
     
   
    !   Flow rule. 
       dirn=TrDev/NormTrDev
    !	Reviewing equation for tensile space
       if (Trp<=0) then
         p=-p 
         pc=-pc 
       endif
	   ! if (Trp<=0) then
       !   dirn=-dirn
        !endif
    !	Plastic Strains
       if (NormTrDev==0) then 
       dEplas3=0.0d0
       CALL NormE(dEplas3, NormdEplas) 
       else
       dEplas3=dphi*(((1.0d0/3.0d0)*(2.0d0*p-pc))*Uni+(sqrt(3.0d0/2.0d0)*(2.0d0*q/(M**2.0d0)))*dirn)
       CALL NormE(dEplas3, NormdEplas) 
       endif
       evoid=evoid+(1.0d0+evoid)*(DeltaE(1,1)+DeltaE(2,2)+DeltaE(3,3)) 
       CALL R36(dEplas3,dEplas6)
       Pstr=Pstr+dEplas6
       CALL dCon(cElas,dEplas3,dTrStr)
       T=TrStr-dTrStr
       p=1.0d0/3.0d0*(T(1,1)+T(2,2)+T(3,3))
       CALL dev(T,TDev,Uni) 
       CALL normS(TDev,NormTDev)
       q=sqrt(3.0d0/2.0d0)*NormTDev 

        do i = 1, 3
           sig(i) = T(i,i)
	    end do
	       sig(4) = T(2,3)
	       sig(5) = T(1,3)
	       sig(6) = T(1,2)
	    ! mean stress and deviate stresses
	      sm = (sig(1)+sig(2)+sig(3))/3.0d0
	      sd(1) = sig(1) - sm
	      sd(2) = sig(2) - sm
	      sd(3) = sig(3) - sm
	      sd(4) = sig(4)
	      sd(5) = sig(5)
	      sd(6) = sig(6)

       ! Call Solution(NTENS, NDI, NSHR, TrStr, STRESS, cElas, DDSDDE)
   	
	    !Saving State variables 
	    plastic_strain(1)= Pstr(1)
        plastic_strain(2)= Pstr(2)
	    plastic_strain(3)= Pstr(3)
	    plastic_strain(5)= Pstr(4) 
	    plastic_strain(6)= Pstr(5)
	    plastic_strain(4)= Pstr(6)
        !print*, pc
		seqv=EquivalentStress()
		
        ps=(Pstr(1)+Pstr(2)+Pstr(3))/3
        epeff_= sqrt((2.0d0/3.0d0)*((pstr(1))**2.0d0+(pstr(2))**2.0d0+(pstr(3))**2.0d0 &
               +(pstr(4))**2.0d0+(pstr(5))**2.0d0+(pstr(6))**2.0d0))
	   
        return 
    endif

	   
        ieinc = ieinc + (sd(1)+sm)*dinc(1) + (sd(2)+sm)*dinc(2) + &
            (sd(3)+sm)*dinc(3) + sd(4)*dinc(4) + &
            sd(5)*dinc(5) + sd(6)*dinc(6)
    end subroutine M3DM14
	

    subroutine sumc(dgam, lalgva, rprops, elastic_strain, trial_elastic_strain, epbarn, stres)
    !***********************************************************************
    ! state update procedure for mohr-coulomb type elasto-plastic material
    ! with associative/non-associative flow rule and piece-wise linear
    ! isotropic hardening:
    ! implicit elastic predictor/return mapping algorithm (boxes 8.4-7).
    ! plane strain and axisymmetric implmentations.
    !
    ! reference: boxes 8.4-7
    !            section 8.2.2
    !***********************************************************************

    implicit double precision (a-h, o-z)
    parameter (iphard=8, mstre=6)
    ! arguments
    real(8), intent(in):: rprops(100), trial_elastic_strain(6)
    real(8), intent(inout):: elastic_strain(6), stres(6), epbarn
    logical, intent(out):: lalgva(5)
    real(8), intent(out):: dgam(2)
    ! local variables and arrays
    logical apex, edge, ifplas, right, sufail, twov
    dimension eigprj(mstre, 3), pstrs(3), strest(6)
    data r0, r1, r2, r3, r4, small, tol/0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 1.d-06, 1.d-10/
    data mxiter/50/

    ! initialize some algorithmic and internal variables
    dgama = r0
    dgamb = r0
    epbar = epbarn

    ! set some material properties
    young = rprops(2) ![density, young, poisson, sinphi, cosphi, sinpsi, nhard]
    poiss = rprops(3)
    sinphi = rprops(4)
    cosphi = rprops(5)
    sinpsi = rprops(6)
    nhard = int(rprops(7))

    ! set some constants
    gmodu = young/(r2*(r1+poiss))
    bulk = young/(r3*(r1-r2*poiss))
    r2g = r2*gmodu
    r4g = r4*gmodu
    r2bulk = r2*bulk
    r2cphi = r2*cosphi
    r1d3 = r1/r3

    ! compute elastic trial state
    ! ---------------------------
    ! elastic trial volumetric strain and pressure stress
    eetv = trial_elastic_strain(1) + trial_elastic_strain(2) + trial_elastic_strain(3)
    trial_pressure = bulk*eetv

    ! spectral decomposition of the elastic trial stress
    eevd3 = eetv*r1d3
    strest(1) = r2g*(trial_elastic_strain(1)-eevd3) + trial_pressure
    strest(2) = r2g*(trial_elastic_strain(2)-eevd3) + trial_pressure
    strest(3) = r2g*(trial_elastic_strain(3)-eevd3) + trial_pressure
    strest(4) = gmodu*trial_elastic_strain(4)
    strest(5) = gmodu*trial_elastic_strain(5)
    strest(6) = gmodu*trial_elastic_strain(6)
    call spdec3(eigprj, pstrs, strest)

    ! identify maximum (pstrs1) and minimum (pstrs3) principal stresses
    ii = 1
    jj = 1
    pstrs1 = pstrs(ii)
    pstrs3 = pstrs(jj)
    do i = 2, 3
        if (pstrs(i)>=pstrs1) then
            ii = i
            pstrs1 = pstrs(ii)
        end if
        if (pstrs(i)<pstrs3) then
            jj = i
            pstrs3 = pstrs(jj)
        end if
    end do
    if (ii/=1 .and. jj/=1) mm = 1
    if (ii/=2 .and. jj/=2) mm = 2
    if (ii/=3 .and. jj/=3) mm = 3
    pstrs2 = pstrs(mm)

    ! compute trial yield function and check for plastic consistency
    ! --------------------------------------------------------------
    ifplas = .false.
    sufail = .true.
    edge = .false.
    apex = .false.
    twov = .false.
    cohesion = plfun(epbarn, nhard, rprops(iphard:))
    smct = pstrs1 - pstrs3 + (pstrs1+pstrs3)*sinphi
    phia = smct - r2cphi*cohesion
    res = phia
    if (cohesion/=r0) res = res/abs(cohesion)
    if (res>tol) then
        ! plastic step: apply return mapping
        ! ==================================
        ifplas = .true.
    else
        sufail = .false.
    end if

    if (ifplas) then
        ! identify possible edge return: either right or left of main plane
        scaprd = pstrs1*(r1-sinpsi) + pstrs2*(-r2) + pstrs3*(r1+sinpsi)
        if (scaprd>=r0) then
            right = .true.
        else
            right = .false.
        end if

        ! apply one-vector return mapping first (return to main plane)
        ! ------------------------------------------------------------
        sphsps = sinphi*sinpsi
        consta = r4g*(r1+r1d3*sphsps) + r4*bulk*sphsps
        r4c2ph = r2cphi*r2cphi
        ! start newton-raphson iterations for dgama
        do nriter = 1, mxiter
            ! compute residual derivative
            denom = -consta - r4c2ph*dplfun(epbar, nhard, rprops(iphard:))

            ! compute newton-raphson increment and update variable dgama
            ddgama = -phia/denom
            dgama = dgama + ddgama

            ! compute new residual
            epbar = epbarn + r2cphi*dgama
            cohesion = plfun(epbar, nhard, rprops(iphard:))
            phia = smct - consta*dgama - r2cphi*cohesion

            ! check convergence
            resnor = abs(phia)
            if (smct/=r0) resnor = resnor/abs(smct)
            if (resnor<=tol) then
                ! check validity of 1-vector return (check sextant of converged stress)
                s1 = pstrs1 - (r2g*(r1+r1d3*sinpsi)+r2bulk*sinpsi)*dgama
                s2 = pstrs2 + (r4g*r1d3-r2bulk)*sinpsi*dgama
                s3 = pstrs3 + (r2g*(r1-r1d3*sinpsi)-r2bulk*sinpsi)*dgama
                delta = dmax1(abs(s1), abs(s2), abs(s3))*small
                if (s1+delta>=s2 .and. s2+delta>=s3) then
                    ! converged stress is in the same sextant as trial stress -> 1-vector
                    ! return is valid.
                    pressure = (s1+s2+s3)*r1d3
                    sufail = .false.
                    exit
                else
                    ! converged stress is not in the same sextant -> 1-vector result is
                    ! not valid. go to two-vector return map to edge
                    twov = .true.
                    exit
                end if
            end if
        end do
    end if

    if (twov) then
        ! apply two-vector return mapping to appropriate edge
        ! ---------------------------------------------------
        dgama = r0
        epbar = epbarn
        cohesion = plfun(epbarn, nhard, rprops(iphard:))
        smcta = pstrs1 - pstrs3 + (pstrs1+pstrs3)*sinphi
        if (right) then
            smctb = pstrs1 - pstrs2 + (pstrs1+pstrs2)*sinphi
        else
            smctb = pstrs2 - pstrs3 + (pstrs2+pstrs3)*sinphi
        end if
        phia = smcta - r2cphi*cohesion
        phib = smctb - r2cphi*cohesion
        if (right) then
            constb = r2g*(r1+sinphi+sinpsi-r1d3*sphsps) + r4*bulk*sphsps
        else
            constb = r2g*(r1-sinphi-sinpsi-r1d3*sphsps) + r4*bulk*sphsps
        end if

        ! start newton-raphson iterations for dgama and dgamb
        do nriter = 1, mxiter
            ! compute residual derivative matrix
            facta = r4c2ph*dplfun(epbar, nhard, rprops(iphard:))
            drvaa = -consta - facta
            drvab = -constb - facta
            drvba = -constb - facta
            drvbb = -consta - facta

            ! compute newton-raphson increment and update variables dgama and dgamb
            r1ddet = r1/(drvaa*drvbb-drvab*drvba)
            ddgama = (-drvbb*phia+drvab*phib)*r1ddet
            ddgamb = (drvba*phia-drvaa*phib)*r1ddet
            dgama = dgama + ddgama
            dgamb = dgamb + ddgamb

            ! compute new residual
            epbar = epbarn + r2cphi*(dgama+dgamb)
            cohesion = plfun(epbar, nhard, rprops(iphard:))
            phia = smcta - consta*dgama - constb*dgamb - r2cphi*cohesion
            phib = smctb - constb*dgama - consta*dgamb - r2cphi*cohesion

            ! check convergence
            resnor = (abs(phia)+abs(phib))
            factor = (abs(smcta)+abs(smctb))
            if (factor/=r0) resnor = resnor/factor
            if (resnor<=tol) then
                ! check validity of 2-vector return to edge
                aux1 = r2g*(r1+r1d3*sinpsi) + r2bulk*sinpsi
                aux2 = (r4g*r1d3-r2bulk)*sinpsi
                aux3 = r2g*(r1-r1d3*sinpsi) - r2bulk*sinpsi
                if (right) then
                    s1 = pstrs1 - aux1*(dgama+dgamb)
                    s2 = pstrs2 + aux2*dgama + aux3*dgamb
                    s3 = pstrs3 + aux3*dgama + aux2*dgamb
                else
                    s1 = pstrs1 - aux1*dgama + aux2*dgamb
                    s2 = pstrs2 + aux2*dgama - aux1*dgamb
                    s3 = pstrs3 + aux3*(dgama+dgamb)
                end if
                delta = dmax1(abs(s1), abs(s2), abs(s3))*small
                if (s1+delta>=s2 .and. s2+delta>=s3) then
                    ! converged stress is in the same sextant as trial stress -> 2-vector
                    ! return to edge is valid.
                    edge = .true.
                    pressure = (s1+s2+s3)*r1d3
                    sufail = .false.
                    exit
                else
                    ! converged stress is not in the same sextant -> 2-vector return to edge
                    ! is not valid. go to two-vector return map to apex
                    apex = .true.
                    exit
                end if
            end if
        end do
    end if

    if (apex) then
        ! apply multi-vector return mapping to apex
        ! ---------------------------------------
        ! check conditions for which return to apex does not make sense
        if (sinphi==r0) write(*,*)'ee0009'
        if (sinpsi==r0) write(*,*)'ee0010'

        ! set initial guess for volumetric plastic strain increment depv
        depv = r0
        epbar = epbarn
        cohesion = plfun(epbar, nhard, rprops(iphard:))
        cotphi = cosphi/sinphi
        res = cotphi*cohesion - trial_pressure

        ! newton-raphson iterations for depv
        do nriter = 1, mxiter
            denom = cosphi*cotphi/sinpsi*dplfun(epbar, nhard, rprops(iphard:)) + bulk
            ddepv = -res/denom
            depv = depv + ddepv
            epbar = epbarn + cosphi/sinpsi*depv
            cohesion = plfun(epbar, nhard, rprops(iphard:))
            pressure = trial_pressure - bulk*depv
            res = cotphi*cohesion - pressure

            ! check for convergence
            resnor = abs(res)
            if (trial_pressure/=r0) resnor = resnor/abs(trial_pressure)
            if (resnor<=tol) then
                dgama = depv
                dgamb = r0
                ! update principal stresses
                s1 = pressure
                s2 = pressure
                s3 = pressure
                sufail = .false.
                exit
            end if
        end do
    end if

    if (sufail) then
        write(*,*) 'WE0003'
    end if

    if (ifplas) then
        ! update internal variable epbar  and stress components
        ! -----------------------------------------------------
        epbarn = epbar
        pstrs(ii) = s1
        pstrs(jj) = s3
        pstrs(mm) = s2
        stres(1) = pstrs(1)*eigprj(1, 1) + pstrs(2)*eigprj(1, 2) + pstrs(3)*eigprj(1, 3)
        stres(2) = pstrs(1)*eigprj(2, 1) + pstrs(2)*eigprj(2, 2) + pstrs(3)*eigprj(2, 3)
        stres(3) = pstrs(1)*eigprj(3, 1) + pstrs(2)*eigprj(3, 2) + pstrs(3)*eigprj(3, 3)
        stres(4) = pstrs(1)*eigprj(4, 1) + pstrs(2)*eigprj(4, 2) + pstrs(3)*eigprj(4, 3)
        stres(5) = pstrs(1)*eigprj(5, 1) + pstrs(2)*eigprj(5, 2) + pstrs(3)*eigprj(5, 3)
        stres(6) = pstrs(1)*eigprj(6, 1) + pstrs(2)*eigprj(6, 2) + pstrs(3)*eigprj(6, 3)

        ! and elastic engineering strain
        eevd3 = pressure/bulk*r1d3
        elastic_strain(1) = (stres(1)-pressure)/r2g + eevd3
        elastic_strain(2) = (stres(2)-pressure)/r2g + eevd3
        elastic_strain(3) = (stres(3)-pressure)/r2g + eevd3
        elastic_strain(4) = stres(4)/gmodu
        elastic_strain(5) = stres(5)/gmodu
        elastic_strain(6) = stres(6)/gmodu
    else
        ! elastic step: update stress using linear elastic law
        ! ====================================================
        stres(1) = strest(1)
        stres(2) = strest(2)
        stres(3) = strest(3)
        stres(4) = strest(4)
        stres(5) = strest(5)
        stres(6) = strest(6)

        ! elastic engineering strain
        elastic_strain(1) = trial_elastic_strain(1)
        elastic_strain(2) = trial_elastic_strain(2)
        elastic_strain(3) = trial_elastic_strain(3)
        elastic_strain(4) = trial_elastic_strain(4)
        elastic_strain(5) = trial_elastic_strain(5)
        elastic_strain(6) = trial_elastic_strain(6)
    end if

    ! update algorithmic variables before exit
    ! ========================================
    dgam(1) = dgama
    dgam(2) = dgamb
    lalgva(1) = ifplas
    lalgva(2) = sufail
    lalgva(3) = edge
    lalgva(4) = right
    lalgva(5) = apex
    end subroutine sumc

    
    Subroutine jacob(a, d, v, n)
    Implicit Double Precision (A-H, O-Z)
    Parameter (mjiter=50, nmax=100)
    Dimension a(n, n), d(n, n), v(n, n), acopy(n, n)
    Dimension b(nmax), z(nmax)
    Data r0, rp2, rp5, r1, r100/0.0D0, 0.2D0, 0.5D0, 1.0D0, 100.0D0/
    Data toler/1.0D-12/
    !***********************************************************************
    ! JACOBI ITERATIVE PROCEDURE FOR SPECTRAL DECOMPOSITION OF A
    ! N-DIMENSIONAL SYMMETRIC MATRIX
    !
    ! REFERENCE: WH Press, SA Teukolsky, WT Vetting & BP Flannery. Numerical
    !            recipes in FORTRAN: The art of scientific computing. 2nd
    !            Edn., Cambridge University Press, 1992.
    !***********************************************************************

    Do ip = 1, n
        Do iq = 1, n
            v(ip, iq) = r0
            acopy(ip, iq) = a(ip, iq)
        End Do
        v(ip, ip) = r1
    End Do
    Do ip = 1, n
        b(ip) = acopy(ip, ip)
        d(ip, ip) = b(ip)
        z(ip) = r0
    End Do
    Do i = 1, mjiter
        sm = r0
        Do ip = 1, n - 1
            Do iq = ip + 1, n
                sm = sm + abs(acopy(ip,iq))
            End Do
        End Do
        If (sm<toler) return
        If (i<4) Then
            tresh = rp2*sm/dble(n**2)
        Else
            tresh = r0
        End If
        Do ip = 1, n - 1
            Do iq = ip + 1, n
                g = r100*abs(acopy(ip,iq))
                If ((i>4) .And. (abs(d(ip,ip))+g==abs(d(ip,ip))) .And. (abs(d(iq,iq))+g==abs(d(iq,iq)))) Then
                    acopy(ip, iq) = r0
                Else If (abs(acopy(ip,iq))>tresh) Then
                    h = d(iq,iq) - d(ip,ip)
                    If (abs(h)+g==abs(h)) Then
                        t = acopy(ip, iq)/h
                    Else
                        theta = rp5*h/acopy(ip, iq)
                        t = r1/(abs(theta)+sqrt(r1+theta**2))
                        If (theta<r0) t = -t
                    End If
                    c = r1/sqrt(r1+t**2)
                    s = t*c
                    tau = s/(r1+c)
                    h = t*acopy(ip, iq)
                    z(ip) = z(ip) - h
                    z(iq) = z(iq) + h
                    d(ip,ip) = d(ip,ip) - h
                    d(iq,iq) = d(iq,iq) + h
                    acopy(ip, iq) = r0
                    Do j = 1, ip - 1
                        g = acopy(j, ip)
                        h = acopy(j, iq)
                        acopy(j, ip) = g - s*(h+g*tau)
                        acopy(j, iq) = h + s*(g-h*tau)
                    End Do
                    Do j = ip + 1, iq - 1
                        g = acopy(ip, j)
                        h = acopy(j, iq)
                        acopy(ip, j) = g - s*(h+g*tau)
                        acopy(j, iq) = h + s*(g-h*tau)
                    End Do
                    Do j = iq + 1, n
                        g = acopy(ip, j)
                        h = acopy(iq, j)
                        acopy(ip, j) = g - s*(h+g*tau)
                        acopy(iq, j) = h + s*(g-h*tau)
                    End Do
                    Do j = 1, n
                        g = v(j, ip)
                        h = v(j, iq)
                        v(j, ip) = g - s*(h+g*tau)
                        v(j, iq) = h + s*(g-h*tau)
                    End Do
                End If
            End Do
        End Do
        Do ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip,ip) = b(ip)
            z(ip) = r0
        End Do
    End Do

    Return
    End Subroutine jacob
    
    
    subroutine spdec3(eigprj, eigx, x)
    implicit none
    integer:: mcomp, ndim,  idir, jdir
    parameter (mcomp = 6, ndim = 3)
    real(8):: eigprj(mcomp,ndim), eigx(ndim), x(6), eigmult, eigx_copy(ndim, ndim)
    real(8):: differ1, differ2, differ3, amxeig, small = 1.0d-8
    integer:: repeat, norepeat
    real(8):: auxmtx(ndim,ndim), eigvec(ndim,ndim), x_new(6), mat(6), xtemp(6)
    !***********************************************************************
    ! performs the closed form spectral decomposition of a
    ! symmetric 3-d tensor stored in vector form
    !
    ! reference: page 739
    !***********************************************************************

    auxmtx(1,1) = x(1)
    auxmtx(2,2) = x(2)
    auxmtx(3,3) = x(3)
    auxmtx(2,3) = x(4)
    auxmtx(1,3) = x(5)
    auxmtx(1,2) = x(6)
    auxmtx(2,1) = auxmtx(1,2)
    auxmtx(3,2) = auxmtx(2,3)
    auxmtx(3,1) = auxmtx(1,3)

    call jacob(auxmtx,eigx_copy,eigvec,ndim)

    eigx(1) = eigx_copy(1,1)
    eigx(2) = eigx_copy(2,2)
    eigx(3) = eigx_copy(3,3)
    differ1 = abs(eigx(1) - eigx(2))
    differ2 = abs(eigx(2) - eigx(3))
    differ3 = abs(eigx(1) - eigx(3))

    amxeig = dmax1(abs(eigx(1)),abs(eigx(2)),abs(eigx(3)))
    if (amxeig .ne. 0) then
        differ1 = differ1/amxeig
        differ2 = differ2/amxeig
        differ3 = differ3/amxeig
    end if
    repeat = 0
    norepeat = 0
    if (differ1 .lt. small) then
        repeat = 2
        norepeat = 3
    else if (differ2 .lt. small) then
        repeat = 2
        norepeat = 1
    else if (differ3 .lt. small) then
        repeat = 2
        norepeat = 2
    end if
    if (differ1 .lt. small .and. differ2 .lt. small) then
        repeat = 3
    end if

    eigprj = 0.0d0
    if (repeat == 0) then
        ! case-1: three distince eigvalues
        do idir = 1, ndim
            eigprj(1, idir) = eigvec(1,idir)*eigvec(1,idir)
            eigprj(2, idir) = eigvec(2,idir)*eigvec(2,idir)
            eigprj(3, idir) = eigvec(3,idir)*eigvec(3,idir)
            eigprj(4, idir) = eigvec(2,idir)*eigvec(3,idir)
            eigprj(5, idir) = eigvec(1,idir)*eigvec(3,idir)
            eigprj(6, idir) = eigvec(1,idir)*eigvec(2,idir)
        end do
    else if (repeat == 3) then
        ! case-2: three equal eigenvalues
        idir = 1
        eigprj(1, idir) = 1.0d0
        eigprj(2, idir) = 1.0d0
        eigprj(3, idir) = 1.0d0
        eigprj(4, idir) = 0d0
        eigprj(5, idir) = 0d0
        eigprj(6, idir) = 0d0
    else if (repeat == 2) then
        ! case-3: two equal eigenvalues
        idir = norepeat
        eigprj(1, idir) = eigvec(1,idir)*eigvec(1,idir)
        eigprj(2, idir) = eigvec(2,idir)*eigvec(2,idir)
        eigprj(3, idir) = eigvec(3,idir)*eigvec(3,idir)
        eigprj(4, idir) = eigvec(2,idir)*eigvec(3,idir)
        eigprj(5, idir) = eigvec(1,idir)*eigvec(3,idir)
        eigprj(6, idir) = eigvec(1,idir)*eigvec(2,idir)
        do idir = 1, ndim
            if (idir .ne. norepeat) then
                eigprj(1, idir) = 1d0 - eigprj(1, norepeat)
                eigprj(2, idir) = 1d0 - eigprj(2, norepeat)
                eigprj(3, idir) = 1d0 - eigprj(3, norepeat)
                eigprj(4, idir) = -eigprj(4, norepeat)
                eigprj(5, idir) = -eigprj(5, norepeat)
                eigprj(6, idir) = -eigprj(6, norepeat)
                exit
            end if
        end do
    end if

    !x_new(1) = eigx(1)*eigprj(1, 1) + eigx(2)*eigprj(1, 2) + eigx(3)*eigprj(1, 3)
    !x_new(2) = eigx(1)*eigprj(2, 1) + eigx(2)*eigprj(2, 2) + eigx(3)*eigprj(2, 3)
    !x_new(3) = eigx(1)*eigprj(3, 1) + eigx(2)*eigprj(3, 2) + eigx(3)*eigprj(3, 3)
    !x_new(4) = eigx(1)*eigprj(4, 1) + eigx(2)*eigprj(4, 2) + eigx(3)*eigprj(4, 3)
    !x_new(5) = eigx(1)*eigprj(5, 1) + eigx(2)*eigprj(5, 2) + eigx(3)*eigprj(5, 3)
    !x_new(6) = eigx(1)*eigprj(6, 1) + eigx(2)*eigprj(6, 2) + eigx(3)*eigprj(6, 3)

    return
    end subroutine spdec3


    real(8) function plfun(x, npoint, xfx)
    integer npoint, i
    real(8):: x,xfx(2,npoint)
    !***********************************************************************
    ! piecewise linear function defined by a set of npoint pairs
    ! {x,f(x)} stored in the matrix xfx (dim. 2*npoint).
    !***********************************************************************
    plfun=xfx(2,npoint)
    do i=1,npoint
        if (x.ge.xfx(1,i)) then
            cycle
        else
            if (i.eq.1) then
                !           -- x < x1 --> f(x)=f(x1) ---
                plfun=xfx(2,1)
                exit
            else
                !           -- x(i-1) <= x < x(i) ---
                plfun=xfx(2,i-1)+(x-xfx(1,i-1))*(xfx(2,i)-xfx(2,i-1))/(xfx(1,i)-xfx(1,i-1))
                exit
            endif
        endif
    end do
    !     ----  x >= x(npoint) --> f(x) = f(x(npoint))  ---
    return
    end function plfun

    real(8) function dplfun(x, npoint, xfx)
    integer npoint, i
    real(8) x, xfx(2,npoint), r0
    data r0 / 0.0d0 /
    !***********************************************************************
    ! derivative of the piecewise linear function 'plfun' defined by a set
    ! of npoint pairs {x,f(x)} stored in the matrix xfx (dim. 2*npoint).
    !***********************************************************************
    dplfun=r0
    do i=1,npoint
        if (x.ge.xfx(1,i)) then
            cycle
        else
            if (i.eq.1) then
                !           -- x < x1   --> f(x)=f(x1) --> df(x)/dx=0 ---
                dplfun=r0
                exit
            else
                !           -- x(i-1) <= x < x(i) ---
                dplfun=(xfx(2,i)-xfx(2,i-1))/(xfx(1,i)-xfx(1,i-1))
                exit
            endif
        endif
    end do
    !     ---- x >= x(npoint) --> f(x) = f(x(npoint)) --> df/dx=0 ---
    return
    end function dplfun

    subroutine seleos(failure)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Update pressure and energy using appropriate              -
    !-      equation of state                                         -
    !-  Input                                                         -
    !-      etype - Type of equation of state                         -
    !-      failure - failure state of particle                       -
    !------------------------------------------------------------------
    implicit none
    logical:: failure

    select case(etype_)

    case(1)
        call eos1(failure)
    case(2)
        call eos2(failure)
    case(3)
        call eos3(failure)
        case default
        stop "error eos id"
    end select

    end subroutine seleos


    subroutine eos1(failure)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Linear Equation Of State                                  -
    !-  Inputs                                                        -
    !-      mu      - ( = rho/rho0 - 1 )                              -
    !-      iener   - internal energy                                 -
    !-      den_    - density                                         -
    !-  Outputs                                                       -
    !-      sm      - mean stress (pressure)                          -
    !-  References                                                    -
    !-      Section 5.2.3                                             -
    !------------------------------------------------------------------
    use ParticleData, only: EngInternal
    implicit none
    real:: A, B, pnew, c0, c1, c2, c3, c4, c5, c6
    logical:: failure

    c0 = mat_list(mid)%cEos(1)
    c1 = mat_list(mid)%cEos(2)
    c2 = mat_list(mid)%cEos(3)
    c3 = mat_list(mid)%cEos(4)
    c4 = mat_list(mid)%cEos(5)
    c5 = mat_list(mid)%cEos(6)
    c6 = mat_list(mid)%cEos(7)

    A = c0 + mu * (c1 + mu * (c2 + mu * c3))
    B = c4 + mu * (c5 + mu * c6)

    call hieupd()

    pnew = (A + B * specen) / (1 + B * dvol / vol0_) !(Eq: 5.18)

    if(failure .and. pnew < 0) pnew=0.0

    iener = iener - dvol * pnew        ! (Eq: 5.14)

    sm = -pnew

    end subroutine eos1


    subroutine eos2(failure)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Mie-Gruneisen Equation Of State                           -
    !-  Inputs                                                        -
    !-      mu      - ( = rho/rho0 - 1 )                              -
    !-      iener   - internal energy                                 -
    !-      den_    - density                                         -
    !-  Outputs                                                       -
    !-      sm      - mean stress (pressure)                          -
    !-  References                                                    -
    !-      Section 5.3.5                                             -
    !------------------------------------------------------------------
    use ParticleData, only: EngInternal
    implicit none
    real:: A, B, pnew, c1, c2, c3, c4, c5
    logical:: failure

    c1 = mat_list(mid)%cEos(1)
    c2 = mat_list(mid)%cEos(2)
    c3 = mat_list(mid)%cEos(3)
    c4 = mat_list(mid)%cEos(4)
    c5 = mat_list(mid)%cEos(5)

    if (mu .gt. 0) then      ! compressed
        A = (c1*mu+c2*mu**2 + c3*mu**3) * (1-c4*mu*0.5d0/den_)
        B = c5
    else          ! expanded
        !mu = max(mu,-0.5)
        A = c1*mu
        B = 0
    end if

    call bulkq()
    call hieupd()

    pnew = (A + B * specen) / (1 + B * dvol / vol0_) !(Eq: 5.18)

    if(failure .and. pnew < 0) pnew=0.0

    iener = iener - dvol * pnew        ! (Eq: 5.14)

    sm = -pnew

    end subroutine eos2


    subroutine eos3(failure)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      JWL Equation Of State                                     -
    !-  Inputs                                                        -
    !-      rv    - Relative volume rv = vol_/vol0_                   -
    !-      Eng - specific internal energy                            -
    !-             (internal energy per unit mass)                    -
    !-  Outputs                                                       -
    !-      sm  - mean stress (pressure)                              -
    !-  References                                                    -
    !-      Section 5.2.4                                             -
    !------------------------------------------------------------------
    use ParticleData, only: EngInternal
    implicit none

    real(8):: r1v, r2v, wr1v, wr2v, wdr1v, wdr2v, er1v, er2v
    real(8):: A, B, pnew, c1, c2, r1, r2, w1
    logical:: failure

    ! initialize parameters
    c1 = mat_list(mid)%cEos(1)
    c2 = mat_list(mid)%cEos(2)
    r1 = mat_list(mid)%cEos(3)
    r2 = mat_list(mid)%cEos(4)
    w1 = mat_list(mid)%cEos(5)

    r1v = r1*rv
    r2v = r2*rv
    wr1v = c1*w1/r1v
    wr2v = c2*w1/r2v
    wdr1v = c1 - wr1v
    wdr2v = c2 - wr2v
    er1v = exp(-r1v)
    er2v = exp(-r2v)

    ! calculate bulk viscosity

    A = wdr1v*er1v + wdr2v*er2v + bqf
    B = w1/rv

    call hieupd()

    A = A * bfac
    B = B * bfac

    pnew = (A + B * specen) / (1 + B * dvol / vol0_) ! (Eq: 5.18)

    if(failure .and. pnew < 0) pnew=0.0
    iener = iener - dvol * pnew        ! (Eq: 5.14)

    sm = -pnew

    end subroutine eos3


    subroutine bulkq()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      bulk viscosity                                            -
    !-  Inputs                                                        -
    !-      bq1,bq2                                                   -
    !-  Outputs                                                       -
    !-      bqf                                                       -
    !------------------------------------------------------------------
    use GridData, only: DCell
    use ParticleData, only: DT
    implicit none
    real(8):: dd    ! bulk strain rate

    dd = (dinc(1) + dinc(2) + dinc(3)) / DT
    dd = min(dd,0.0)

    bqf = den_ * DCell * DCell * bq1 * dd * dd - &
        bq2 * den_ * DCell * cp * dd      !(Eq. 2.82)

    end subroutine bulkq

 ! subroutines for MCC model
 
SUBROUTINE normE(C,N)
real(8):: C(3,3),D(3,3),N
!returns the norm of a second order tensor A (strain)
 D=matmul(C,C)
 N=sqrt(D(1,1)+D(2,2)+D(3,3)+0.5d0*(D(1,2)+D(1,3)+D(2,3))) 
END SUBROUTINE normE

SUBROUTINE TenUni(Uni)
    real(8):: Uni(3,3) 
	integer:: i, j
	!Tensor unitario de 3*3 
	do i=1,3
    do j=1,3
    if (i==j) then 
	Uni(i,j)=1 
	else 
	Uni(i,j)=0 
	endif
    enddo
	enddo
END SUBROUTINE

SUBROUTINE TenIsoUni(IdenT)
     integer:: i, k, l, m
     real(8)::IdenT(3,3,3,3)
	 do i=1,3
     do k=1,3 
	 do l=1,3
	 do m=1,3
     if ((i==l).and.(k==m).and.(i==m).and.(k==l)) then 
	  IdenT(i,k,l,m)=1.0d0
    else if (((i/=l).or.(k/=m)).and.((i/=m).or.(k/=l)).or.((l/=i).or. &
           (m/=k)).and.((m/=i).or.(l/=k)))then 
		   IdenT(i,k,l,m)=0.0d0
    else 
	     IdenT(i,k,l,m)=0.5d0
	endif
    enddo
	enddo 
	enddo 
	enddo
END SUBROUTINE TenIsoUni

SUBROUTINE Initial( sm, sd, dinc, T, DeltaE)
integer :: i, j
real(8) :: T(3,3), DeltaE(3,3), sm, sd(6), dinc(6)
     
!Stress and delta strain in a 3x3-matrix 
    T=0.0d0
    DeltaE=0.0d0 

! transform stress vector into stress matrix
    do i = 1,3
	    T(i,i) = sm + sd(i)
	end do
	T(1,2) = sd(6)
	T(2,1) = sd(6)
	T(1,3) = sd(5)
	T(3,1) = sd(5)
	T(2,3) = sd(4)
	T(3,2) = sd(4)


    do i=1,3
    do j=1,3
    if (i==j) then 
	DeltaE(i,j)=dinc(i) 
    endif
	enddo 
	enddo
	
    DeltaE(1,2) = dinc(6)/2.0d0
	DeltaE(2,1) = dinc(6)/2.0d0
	DeltaE(1,3) = dinc(5)/2.0d0
	DeltaE(3,1) = dinc(5)/2.0d0
	DeltaE(2,3) = dinc(4)/2.0d0
	DeltaE(3,2) = dinc(4)/2.0d0
END SUBROUTINE Initial
	
SUBROUTINE TenUni1(IdenT1)
      integer :: i, l, k, m
      real(8) :: IdenT1(3,3,3,3) 
	  do i=1,3
      do l=1,3
	  do k=1,3
      do m=1,3
      if ((i==k).and.(l==m)) then 
	  IdenT1(i,k,l,m)=1.0d0
      else 
	  IdenT1(i,k,l,m)=0.0d0
	  endif
      enddo 
	  enddo 
	  enddo 
	  enddo
END SUBROUTINE TenUni1

SUBROUTINE Emod(LaLambda,LaG,IdenT1,IdenT,cElas)
     real(8):: LaLambda, LaG, IdenT1(3,3,3,3), IdenT(3,3,3,3), cElas(3,3,3,3)
     integer:: i, j, l, n
    ! Elastic modulus 
	 do i=1,3
     do j=1,3 
	 do l=1,3
	 do n=1,3
     cElas(i,j,l,n)=LaLambda*IdenT1(i,j,l,n)+2*LaG*IdenT(i,j,l,n)
	 enddo
     enddo 
	 enddo 
	 enddo
END SUBROUTINE Emod

SUBROUTINE dCon(A,B,Resp)
!	Returns 3*3 matrix
	real(8):: Resp(3,3), A(3,3,3,3), B(3,3) 
	integer :: i, j, k, l
	Resp=0.0d0 
	do i=1,3 
	do j=1,3
	do k=1,3
	do l=1,3
    Resp(i,j)=Resp(i,j)+A(i,j,k,l)*B(k,l)
	enddo
    enddo
	enddo
    enddo
END SUBROUTINE dCon

SUBROUTINE normS(C,N)
   real(8):: C(3,3),D(3,3),N
! returns the norm of a second order tensor A (stress)
   D=matmul(C,C)
   N=sqrt(D(1,1)+D(2,2)+D(3,3)+2.0d0*(D(1,2)+D(1,3)+D(2,3))) 
END SUBROUTINE normS

SUBROUTINE dev(C, devi, Uni)
real(8):: C(3,3), trC, devi(3,3), Uni(3,3)
!returns the deviator part of a second order tensor A 
trC=C(1,1)+C(2,2)+C(3,3)
devi=C-1.0d0/3.0d0*trC*Uni 
END SUBROUTINE dev

SUBROUTINE r36(B,A)
   ! R3-R6
	real(8):: A(6), B(3,3)
	integer i, j
	do i=1,3 
	do j=1,3
	if (i==j) then
        A(i)=B(i,j)
    !else 
      !  A(1+i+j)=2.0d0*B(i,j)
	endif
	A(6)=2.0d0*B(1,2)
	A(5)=2.0d0*B(1,3)
	A(4)=2.0d0*B(2,3)
	enddo 
	enddo
END SUBROUTINE r36
	 
end module MaterialModel
	

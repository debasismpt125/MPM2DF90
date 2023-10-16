
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
    ! ------------------------------------------------------------------
    ! -                                                                -
    ! -  Solution and particles information                            -
    ! -                                                                -
    ! ------------------------------------------------------------------  
    module ParticleData
#ifdef NDIM2
    type Particle
        real(8):: XX(2)     ! particle position at time step t+1
        real(8):: Xp(2)     ! particle position at time step t
        real(8):: VXp(2)    ! particle velocity
        real(8):: FXp(2)    ! load

        real(8):: VOL       ! current volume
        real(8):: sig_y     ! yield stress
        real(8):: SM, SM0, Seqv  ! mean stress and Mises stress
        ! deviatoric stress
        real(8):: SDxx, SDyy, SDzz, SDxy, SDyz, SDxz
        real(8):: epeff     ! effective plastic strain
        real(8):: celsius_t ! celsius temperature

        logical:: SkipThis  ! for plot less
        logical:: failure   ! failure
        integer:: icell     ! cell number
        real(8):: DMG       ! damage
        real(8):: LT        ! lighting time for explosive simulation
        real(8):: ie        ! internal energy
        real(8):: mass      ! particle mass
        real(8):: cp        ! sound speed
        
        integer:: nb_InflNode, InflNode(9)
        real(8):: SHP(9), DNDX(9), DNDY(9), DNDZ(9), rpg(9,2), SHP_QUAD(4)
        
        ! deformation gradient
        real(8):: DG(3,3), Jacobian, DFG(3,3)
        real(8):: strain_increment(6), vort(3)
        real(8):: elastic_strain(6), trial_elastic_strain(6)
		real(8):: plastic_strain(6)
		real(8):: evoid, pc, pore, pwpa, KW, pm, qd, ptot    !state variables for mcc
        ! for strain approximation
        real(8):: dilation, pressure
    end type Particle

    type Body
        integer:: mat       ! material set number (1 ~ nb_mat)
        integer:: comID     ! component set number (1 ~ nb_component)
        real(8):: Gravp(2)  ! gravity
        integer:: par_begin ! the beginning of particle of body
        integer:: par_end   ! the end of particle of body
    end type Body
    
#elif NDIM3
     type Particle
        real(8):: XX(3)     ! particle position at time step t+1
        real(8):: Xp(3)     ! particle position at time step t
        real(8):: VXp(3)    ! particle velocity
        real(8):: FXp(3)    ! load

        real(8):: VOL       ! current volume
        real(8):: sig_y     ! yield stress
        real(8):: SM, SM0, Seqv  ! mean stress and Mises stress
		real(8):: evoid, pc, pore, pwpa, KW, pm, qd, ptot    !state variables for mcc
        ! deviatoric stress
        real(8):: SDxx, SDyy, SDzz, SDxy, SDyz, SDxz
        real(8):: epeff     ! effective plastic strain
        real(8):: celsius_t ! celsius temperature

        logical:: SkipThis  ! for plot less
        logical:: failure   ! failure
        integer:: icell     ! cell number
        real(8):: DMG       ! damage
        real(8):: LT        ! lighting time for explosive simulation
        real(8):: ie        ! internal energy
        real(8):: mass      ! particle mass
        real(8):: cp        ! sound speed
        
        integer:: nb_InflNode, InflNode(27)
        real(8):: SHP(27), DNDX(27), DNDY(27), DNDZ(27), rpg(27,3), SHP_QUAD(9)
        
        ! deformation gradient
        real(8):: DG(3,3), Jacobian, DFG(3,3)
        real(8):: strain_increment(6), vort(3)
        real(8):: elastic_strain(6), trial_elastic_strain(6), plastic_strain(6)
        ! for strain approximation
        real(8):: dilation, pressure
    end type Particle

    type Body
        integer:: mat       ! material set number (1 ~ nb_mat)
        integer:: comID     ! component set number (1 ~ nb_component)
        real(8):: Gravp(3)  ! gravity
        integer:: par_begin ! the beginning of particle of body
        integer:: par_end   ! the end of particle of body
    end type Body
#endif
    
#ifdef NDIM2 
    integer, parameter:: ndim = 2, max_infl_nodes = 9
#elif NDIM3
    integer, parameter:: ndim = 3, max_infl_nodes = 27
#endif
    
    integer:: nb_particle = 0    ! number of particles
    type(Particle), target, allocatable:: particle_list(:) ! particles list

    ! number of components for contact simulating
    integer:: nb_component = 1
    integer:: nb_body = 0        ! number of bodies
    type(Body), target, allocatable:: body_list(:)  ! bodies list

    integer, parameter :: Version = 1, &   !  Version of MPM3D
    Release = 0      !  Release of MPM3D version

    character(256):: Title
    integer:: iTitle(64)
    equivalence (Title, iTitle)

    logical:: MUSL = .false.      ! use MUSL method?
    logical:: USL = .false.       ! use USL method
    logical:: USF = .false.       ! use USF method?
    logical:: GIMP = .false.      ! use cpGIMP method?
    logical:: contact = .false.   ! use contact method?
    logical:: Gravity = .false.   ! set gravity?
    logical:: IncLoad = .false.   ! use incremental loading
    logical:: DEBUG = .false.     ! debug option
    logical:: INITPART = .false.  ! initialize particles
    logical:: BBAR = .false.      ! bbar option
    logical:: CPDI = .false.      ! cpdi option
    logical:: StressApprox = .false.     ! stress approximation option
    logical:: StrainApprox = .false.     ! strain approximation option
    logical:: FiniteStrain = .false. ! finite strain option
    real(8):: ALPHA_PIC = 0.0d0   ! ratio of pic velocity in mixed velocity update

    integer:: istep = 0           ! Current time step
    real(8):: DT = 0              ! time step interval
    real(8):: CurrentTime = 0.0   ! Time for current step
    real(8):: EndTime = 0.0    ! Time for end of solution
    real(8):: DTScale = 0.9    ! time step size factor (<= 1.0)
    real(8):: LoadTime = 1.0d0    ! loading time for incremental load 
    
    real(8):: EngInternal = 0.0   ! Internal energy
    real(8):: EngKinetic = 0.0    ! Kinetic energy
    real(8):: Momentum(ndim)
    real(8):: Mombody1(ndim)
    real(8):: Mombody2(ndim)
    
    real(8), allocatable:: PreOutTime(:)

    contains

    subroutine InitParticle()
    !------------------------------------------------------------------
    !    purpose: initialize particle_list                            -
    !    used after particle_list space allocated                     -
    !------------------------------------------------------------------
    implicit none
    integer :: i, j
    
    particle_list%SM = 0.0d0
    particle_list%seqv = 0.0d0
    particle_list%SDxx = 0.0d0
    particle_list%SDyy = 0.0d0
    particle_list%SDzz = 0.0d0
    particle_list%SDxy = 0.0d0
    particle_list%SDyz = 0.0d0
    particle_list%SDxz = 0.0d0
	particle_list%SM0 =  0.0d0
    particle_list%epeff = 0.0d0
	particle_list%pc = 0.0d0
	particle_list%pore = 0.0d0
	particle_list%KW = 0.0d0
	particle_list%pwpa = 0.0d0
	particle_list%pm = 0.0d0
	particle_list%qd = 0.0d0
	particle_list%ptot = 0.0d0
    particle_list%celsius_t = 293.0

    do i = 1, ndim
        particle_list%VXp(i) = 0.0d0
        particle_list%FXp(i) = 0.0d0
    end do
    
    do i = 1, 3
        do j = 1, 3
            if (i == j) then
                particle_list%DG(i,j) = 1d0
            else
                particle_list%DG(i,j) = 0d0
            end if
        end do
    end do
    
    do i = 1, 6
        particle_list%elastic_strain(i) = 0d0
        particle_list%strain_increment(i) = 0d0
		particle_list%plastic_strain(i) = 0d0
	
    end do
    
    particle_list%DMG = 0.0
    particle_list%ie = 0.0

    particle_list%SkipThis = .false.
    particle_list%failure  = .false.
    particle_list%icell = 0
    
    particle_list%Jacobian = 1d0
    
     do i = 1, max_infl_nodes
        particle_list%SHP(i) = 0d0
        particle_list%DNDX(i) = 0d0
        particle_list%DNDY(i) = 0d0
        particle_list%DNDZ(i) = 0d0
        particle_list%InflNode(i) = 0
    end do
    
    end subroutine InitParticle

    subroutine InitBody
    !------------------------------------------------------------------
    !    purpose: initialize body_list                                -
    !    used after body_list space allocated                         -
    !------------------------------------------------------------------
    implicit none
    integer :: i

    do i = 1, ndim
        body_list%Gravp(i) = 0.0d0
    end do

    end subroutine InitBody

    subroutine calcEnergy()
    !-------------------------------------------------------------------
    !-   purpose: calculate kinematic energy  and internal energy      -
    !----------------------------------------------------------------- -
    implicit none
    real(8) :: mp   ! mass of a material point
    type(Particle), POINTER :: pt

    integer :: p    ! loop counter
    integer :: mat
    real(8)::delta_EngKinetic,iener

    EngKinetic  = 0.0      ! initial Kinetic energy
    EngInternal = 0.0      ! initial internal energy
    Momentum    = 0.0      ! set initial momentum equal to  zero
    do p = 1, nb_particle
        pt => particle_list(p)
        mp = pt%mass
        delta_EngKinetic = DOT_PRODUCT(pt%VXp, pt%VXp)*mp*0.5d0
        EngKinetic = EngKinetic + delta_EngKinetic
        EngInternal = EngInternal + pt%ie
        Momentum = Momentum + mp*pt%VXp
    end do

    end subroutine calcEnergy

    end module ParticleData

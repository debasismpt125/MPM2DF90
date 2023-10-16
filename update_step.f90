
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
    ! -  Calculation procedures                                        -
    ! -                                                                -
    ! ------------------------------------------------------------------
    
    subroutine GridMomentumInitial()
    !-------------------------------------------------------------------
    !-  Purpose                                                        -
    !-      1. The variables of particle is mapped to the grid node    -
    !-------------------------------------------------------------------

    use ParticleData
    use MaterialData
    use GridData
    
    implicit none

    integer:: b, p, n, c, parBegin, parEnd ! loop counter
    integer:: icell, inode, ix, iy, iz, mat_, comID = 1
    real(8):: sxx, syy, szz, sxy, syz, sxz
    real(8):: fx(ndim), f_int(ndim), f_ext(ndim), mp_, vol_, kw_
    real(8):: shm, SHPn, DNDXn, DNDYn, DNDZn

    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd

    ! Calculate the grid nodal masses, moemntum only
    ! Reset Grid data
    grid_list%Mg = 0.0d0         ! Grid nodal mass
    do n = 1, ndim
        grid_list%PXg(n) = 0.0d0    ! Nodal momentum
    end do
    grid_list%KI = 0.0d0         ! Grid nodal mass
	
    do b = 1, nb_body     ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID = body_list(b)%comID ! Get comID from body
        do p = parBegin, parEnd    ! Loop over all particles (1)
            pt => particle_list(p)
            pt%icell = InWhichCell(pt%Xp)
            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            vol_ = pt%VOL
            mp_ = pt%Mass
            kw_ = pt%KW
            ! Loop over the grid nodes of the hexhedron
            !   in which the particle is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                ! out of the computational grid
                if (inode .gt. nb_gridnode .or. &
                    inode .le. 0) cycle

                gd => grid_list(inode, comID)

                SHPn = pt%SHP(n)
                shm = SHPn*mp_

                gd%Mg = gd%Mg + shm            ! the nodal mass
                gd%PXg = gd%PXg + pt%VXp*shm   ! the nodal momentum
				gd%KI = gd%KI + SHPn*vol_/kw_  !
            end do !n

        end do !p
    end do    !b
    end subroutine GridMomentumInitial

    
    subroutine ParticleStressUpdate()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Update stresses by appropriate constitution law           -
    !-                                                                -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    use MaterialModel, only: Constitution
    use MaterialData
    implicit none

    integer:: b, p, parBegin, parEnd ! loop counter
    integer:: icell, comID = 1
    real(8):: de(6), vort(3)
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd

    ! Calculate the increment strain and vorticity
    ! de(i) comply the Voigt rule (d11, d22, d33, 2*d23, 2*d13, 2*d12)
    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        do p = parBegin, parEnd    ! Loop over all particles (4)
            pt => particle_list(p)
            icell = pt%icell    ! use old position
            ! Particle p is out of the computational region
            if (icell < 0) cycle
            
            de   = pt%strain_increment    ! Incremental strain
            vort = pt%vort    ! Incremental vorticity

            ! Update stress by constitution law
            call Constitution(de, vort, b, p)

        end do !p
    end do    !b

    end subroutine ParticleStressUpdate
    
    subroutine ParticlePositionUpdateOnly()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Update stresses by appropriate constitution law           -
    !-                                                                -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    use MaterialModel, only: Constitution
    use MaterialData
    implicit none

    integer:: b, p, parBegin, parEnd ! loop counter
    integer:: icell, comID = 1
    real(8):: de(6), vort(3)
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd

    ! Calculate the increment strain and vorticity
    ! de(i) comply the Voigt rule (d11, d22, d33, 2*d23, 2*d13, 2*d12)
    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        do p = parBegin, parEnd    ! Loop over all particles (4)
            pt => particle_list(p)
            icell = pt%icell    ! use old position
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            pt%Xp = pt%XX     ! the next particle position

        end do !p
    end do    !b

    end subroutine ParticlePositionUpdateOnly
    
    real(8) function MAT3DET(A)
    implicit none
    real(8):: A(3,3)
    !real(8):: detv

    MAT3DET = A(1,1)*A(2,2)*A(3,3)  &
        - A(1,1)*A(2,3)*A(3,2)  &
        - A(1,2)*A(2,1)*A(3,3)  &
        + A(1,2)*A(2,3)*A(3,1)  &
        + A(1,3)*A(2,1)*A(3,2)  &
        - A(1,3)*A(2,2)*A(3,1)
    return
    end function MAT3DET
    
#ifdef NDIM3
    subroutine ParticleStrainUpdate()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Calculate the strain rate and spin tensor                 -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    use MaterialModel, only: Constitution
    use MaterialData
    implicit none

    integer:: b, p, n, parBegin, parEnd ! loop counter
    integer:: icell, inode, ix, iy, iz, comID = 1
    real(8):: xx(3), vx(3), ax(3), vgx(3)
    real(8):: de(6), vort(3)
    real(8):: mp_, shm, SHPn, DNDXn, DNDYn, DNDZn
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd
    real(8):: dFG(3,3), DGnew(3, 3), Jacobian, MAT3DET

    ! Calculate the increment strain and vorticity
    ! de(i) comply the Voigt rule (d11, d22, d33, 2*d23, 2*d13, 2*d12)
    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID = body_list(b)%comID  ! Get comID from body
        do p = parBegin, parEnd    ! Loop over all particles (4)
            pt => particle_list(p)
            icell = pt%icell    ! use old position
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            de   = 0d0    ! Incremental strain
            vort = 0d0    ! Incremental vorticity
            dFG(:,:) = 0d0

            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if (inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle  ! out of the computational grid
                gd => grid_list(inode, comID)
                ! If the nodal mass is not too small
                if (gd%Mg > CutOff) then
                    vgx = gd%PXg / gd%Mg    ! Grid nodal velocity

                    DNDXn = pt%DNDX(n);  DNDYn = pt%DNDY(n);  DNDZn = pt%DNDZ(n)
                    de(1) = de(1) + DNDXn*vgx(1)            ! D11
                    de(2) = de(2) + DNDYn*vgx(2)            ! D22
                    de(3) = de(3) + DNDZn*vgx(3)            ! D33
                    ! 2*D23
                    de(4) = de(4) + (DNDYn*vgx(3) + DNDZn*vgx(2))
                    ! 2*D13
                    de(5) = de(5) + (DNDZn*vgx(1) + DNDXn*vgx(3))
                    ! 2*D12
                    de(6) = de(6) + (DNDXn*vgx(2) + DNDYn*vgx(1))

                    ! W32
                    vort(1) = vort(1) + (DNDYn*vgx(3) - DNDZn*vgx(2))
                    ! W13
                    vort(2) = vort(2) + (DNDZn*vgx(1) - DNDXn*vgx(3))
                    ! W21
                    vort(3) = vort(3) + (DNDXn*vgx(2) - DNDYn*vgx(1))
                    
                    ! Correct Vogit format
                    !! W23
                    !vort(1) = vort(1) + (DNDZn*vgx(2) - DNDYn*vgx(3))
                    !! W13
                    !vort(2) = vort(2) + (DNDZn*vgx(1) - DNDXn*vgx(3))
                    !! W12
                    !vort(3) = vort(3) + (DNDYn*vgx(1) - DNDXn*vgx(2))
                    
                    ! deformation gradient increment
                    dFG(1,1:ndim) = dFG(1,1:ndim) + DNDXn*vgx
                    dFG(2,1:ndim) = dFG(2,1:ndim) + DNDYn*vgx
                    dFG(3,1:ndim) = dFG(3,1:ndim) + DNDZn*vgx
                end if
            end do ! n

            de = de * DT
            vort = vort * DT / 2d0
            
            ! deformation gradient increment 
            dFG = dFG * DT
            dFG(1,1) = dFG(1,1) + 1d0
            dFG(2,2) = dFG(2,2) + 1d0
            dFG(3,3) = dFG(3,3) + 1d0
            
            ! new DG = deltaDG * DG
            DGnew = matmul(dFG, pt%DG)
            Jacobian = MAT3DET(DGnew)
            
            pt%DFG = dFG
            pt%DG = DGnew
            pt%Jacobian = Jacobian
            pt%strain_increment = de
            pt%vort = vort
            pt%trial_elastic_strain = pt%elastic_strain + de ! small strain
            
            ! compute trial elastic strain
            if (FiniteStrain) then
                ! finite strain
                !call compute_trial_elastic_log_strain(elastic_strain, pt%DFG, trial_elastic_strain)
            end if

        end do !p
    end do    !b

    end subroutine ParticleStrainUpdate
#endif
    
#ifdef NDIM3
    subroutine GridMomentumUpdate()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      1. calculate the background grid nodal force              -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    use MaterialData
    implicit none

    integer:: b, p, n, c, parBegin, parEnd ! loop counter
    integer:: icell, inode, ix, iy, iz, mat_, comID = 1
    real(8):: sxx, syy, szz, sxy, syz, sxz
    real(8):: fx(3), f_int(3), f_ext(3), f_dmp(3), mp_, vol_
    real(8):: shm, SHPn, DNDXn, DNDYn, DNDZn, load_factor

    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd
    type(ContactGridNodeProperty), POINTER :: CP

    ! Calculate the grid nodal forces only

    ! Reset nodal forces
    grid_list%FXg(1) = 0.0d0;    ! Nodal forces
    grid_list%FXg(2) = 0.0d0;
    grid_list%FXg(3) = 0.0d0;

    if(contact) then
        CP_list%ndir(1) = 0.0d0
        CP_list%ndir(2) = 0.0d0
        CP_list%ndir(3) = 0.0d0
    end if
    
    ! load factor
    load_factor = 1.0d0
    if (IncLoad) then
        load_factor = CurrentTime/LoadTime
        load_factor = min(load_factor, 1.0d0)
    end if

    do b = 1, nb_body             ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End

        if(contact) comID = body_list(b)%comID  ! Get comID from body
        do p = parBegin, parEnd    ! Loop over all particles
            pt => particle_list(p)

            icell = pt%icell        ! using old position

            ! Particle p is out of the computational region
            if (icell < 0) cycle

            vol_ = pt%VOL
            mp_ = pt%Mass
            sxx = pt%SM + pt%SDxx   ! Stresses
            syy = pt%SM + pt%SDyy
            szz = pt%SM + pt%SDzz
            sxy = pt%SDxy
            syz = pt%SDyz
            sxz = pt%SDxz
            
            ! External forces
            fx = pt%FXp*load_factor
            if (Gravity) then
                fx = fx + mp_ * (body_list(b)%Gravp) * load_factor
            end if

            ! Loop over the grid nodes of the hexhedron
            !  in which the particle is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if (inode .gt. nb_gridnode .or. inode .le. 0d0) &
                    cycle  ! out of the computational grid

                gd => grid_list(inode, comID)

                SHPn = pt%SHP(n)
                DNDXn = pt%DNDX(n)
                DNDYn = pt%DNDY(n)
                DNDZn = pt%DNDZ(n)

                f_int(1) = - (sxx*DNDXn + sxy*DNDYn + sxz*DNDZn)*vol_
                f_int(2) = - (sxy*DNDXn + syy*DNDYn + syz*DNDZn)*vol_
                f_int(3) = - (sxz*DNDXn + syz*DNDYn + szz*DNDZn)*vol_

                f_ext = fx*SHPn

                f_dmp = 0.0d0
               if (CurrentTime.lt.dampDuration) then
                   f_dmp = -dampCoef*abs(f_int + f_ext)*sign(1.0,pt%VXp)
               end if
                
                gd%FXg = gd%FXg + f_int + f_ext + f_dmp ! nodal force

                if(contact) then
                    CP => CP_list(comID, inode)
                    CP%ndir(1) = CP%ndir(1) + DNDXn*mp_
                    CP%ndir(2) = CP%ndir(2) + DNDYn*mp_
                    CP%ndir(3) = CP%ndir(3) + DNDZn*mp_
                end if

            end do !n

        end do !p
    end do    !b
    end subroutine GridMomentumUpdate
#endif

    subroutine IntegrateMomentum()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-    1. Integrating  the momentum equations on the               -
    !=        computational grid                                      -
    !-    2. Apply  boundary conditions                               -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    implicit none

    integer:: c, n, idim, fric_normal ! loop counter
    real(8):: ax_node(ndim), vx_node(ndim)
    type(GridNodeProperty), POINTER :: gd
    type(GridNode), POINTER :: node

    do c = 1, nb_component
        do n = 1, nb_gridnode
            gd => grid_list(n, c)
            node => node_list(n)
            
            ! Integrate momentum equation
            if(istep == 1) then
                gd%PXg = gd%PXg + gd%FXg * DT*0.5d0
            else
                gd%PXg = gd%PXg + gd%FXg * DT
            end if
            
            ! Find normal direction of frictional surface
            fric_normal = 0
            do idim = 1, ndim
                if (node%Fric(idim)) then
                    fric_normal = idim
                end if
            end do
            
            ! Frictional boundary
            if (fric_normal>0 .and. gd%Mg > CutOff) then
                
                call UpdateFrictionalForce(n, c, fric_normal)
            end if
            
            ! Apply boundary conditions on computational grid
            if (node%Fix_x) then  ! Grid node n is fixed in x direction
                gd%PXg(1) = 0.0d0
                gd%FXg(1) = 0.0d0
            end if

            if (node%Fix_y) then  ! Grid node n is fixed in y direction
                gd%PXg(2) = 0.0d0
                gd%FXg(2) = 0.0d0
            end if
            
#ifdef NDIM3
            if (node%Fix_z) then  ! Grid node n is fixed in z direction
                gd%PXg(3) = 0.0d0
                gd%FXg(3) = 0.0d0
            end if
#endif
            
        end do !n
    end do    !c
    end subroutine IntegrateMomentum
    
    subroutine GridVelocityAccelerationUpdate()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-    1. Integrating  the momentum equations on the               -
    !=        computational grid                                      -
    !-    2. Apply  boundary conditions                               -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    implicit none
    integer:: c, n, idim, fric_normal ! loop counter
    real(8):: ax_new(ndim), vx_old(ndim), vx_new(ndim), factor
    type(GridNodeProperty), POINTER :: gd
    type(GridNode), POINTER :: node

    if(istep == 1) then
        factor = 0.5d0
    else
        factor = 1d0
    end if
    do c = 1, nb_component
        do n = 1, nb_gridnode
            gd => grid_list(n, c)
            node => node_list(n)
            
            ! nodal acceleration and velocity
            vx_old = 0d0
            ax_new = 0d0
            if (gd%Mg > CutOff) then
                vx_old = gd%PXg/gd%Mg
                ax_new = gd%FXg/gd%Mg*factor
            end if
            

            ! apply friction boundary
            fric_normal = 0
            do idim = 1, ndim
                if (node%Fric(idim)) then
                    fric_normal = idim
                end if
            end do
            
            ! Frictional boundary
            if (fric_normal>0) then
                call UpdateAccelerationByFriction(ax_new, vx_old, fric_normal)
                ax_new(fric_normal) = 0d0
                vx_old(fric_normal) = 0d0
            end if
            
            ! Apply boundary conditions on computational grid
            if (node%Fix_x) then  ! Grid node n is fixed in x direction
                gd%PXg(1) = 0.0d0
                gd%FXg(1) = 0.0d0
                vx_old(1) = 0d0
                ax_new(1) = 0d0
            end if

            if (node%Fix_y) then  ! Grid node n is fixed in y direction
                gd%PXg(2) = 0.0d0
                gd%FXg(2) = 0.0d0
                vx_old(2) = 0d0
                ax_new(2) = 0d0
            end if
            
#ifdef NDIM3
            if (node%Fix_z) then  ! Grid node n is fixed in z direction
                gd%PXg(3) = 0.0d0
                gd%FXg(3) = 0.0d0
                vx_old(3) = 0d0
                ax_new(3) = 0d0
            end if
#endif
            
            vx_new = vx_old + ax_new*DT
            gd%AXg = ax_new
            gd%VXg = vx_new
            gd%VXg_old = vx_old
            gd%PXg = gd%Mg*vx_new
            gd%FXg = gd%Mg*ax_new
        end do !n
    end do    !c
    end subroutine GridVelocityAccelerationUpdate
    
    subroutine ParticleVelocityPositionUpdate()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      1. Update particle position and velocity                  -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    use MaterialModel, only: Constitution
    use MaterialData
    implicit none

    integer:: b, p, n, parBegin, parEnd ! loop counter
    integer:: icell, inode, ix, iy, iz, comID = 1
    real(8):: SHPn, xx(ndim), ax(ndim)
    real(8):: vx(ndim), vx_flip(ndim), vx_new(ndim)
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd

    ! Update particle position and velocity
    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        ! Get comID from body
        if(contact)  comID = body_list(b)%comID
        do p = parBegin, parEnd    ! Loop over all particles (2)
            pt => particle_list(p)

            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            xx = pt%Xp;  ! Particle position at time step k
            vx = 0d0
            ax = 0d0 ! PIC velocity 

            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if(inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle    ! out of the computational grid
                gd => grid_list(inode, comID)
                SHPn = pt%SHP(n)
                if (gd%Mg > CutOff) then
                    vx = vx + SHPn * gd%VXg
                    ax = ax + SHPn * gd%AXg
                end if
            end do ! n
            
            vx_flip = pt%VXp + DT*ax ! FLIP velocity
            vx_new = vx*ALPHA_PIC + vx_flip*(1d0 - ALPHA_PIC) ! mixed update 
            
            ! Time integration
            pt%XX = xx + vx * DT       ! Update particle position
            pt%VXp = vx_new            ! Update particle velocity
            if(USF)  pt%Xp = pt%XX     ! the next particle position

        end do ! p
    end do    ! b

    end subroutine ParticleVelocityPositionUpdate
    
    subroutine UpdateFrictionalForce(n, c, normal)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Calculate frictional acceleration                         -
    !-                                                                -
    !-                                                                -
    !------------------------------------------------------------------
    use ParticleData, only: ndim, DT, istep
    use GridData, only: grid_list, fric_coef_static, fric_coef_kinematic, GridNodeProperty
    implicit none
    integer:: tangent, normal, c, n
    real(8):: ax_node(ndim), vx_node(ndim), vx_old(ndim), vx_t, ax_n, ax_t, ax_delta, tol = 1d-5
    type(GridNodeProperty), POINTER :: gd
    
    ! grid nodal 
    gd => grid_list(n, c)
    ax_node = gd%FXg/gd%Mg
    vx_node = gd%PXg/gd%Mg
    if(istep == 1) then
        vx_old = (gd%PXg - gd%FXg * DT*0.5d0)/gd%Mg
    else
        vx_old = (gd%PXg - gd%FXg * DT)/gd%Mg
    end if
    
    ax_n = ax_node(normal) ! normal acceleration

    ! normal acceleration opposite to the face
    if (ax_n>0d0) then
        ! no friction
        return
    end if
    
    ! update acceleration in all tangent directions 
    do tangent = 1, ndim
        if (tangent==normal) then
            ! pass
            cycle
        else
            ! except for normal direction, 
            ! the other directions are tangent directions
            ax_t = ax_node(tangent)
            vx_t = vx_node(tangent)
        end if
        
        ! tangent acceleration under frictional effect
        if (dabs(vx_t)<tol) then
            ! static frictional condition
            if (dabs(ax_t)<fric_coef_static*dabs(ax_n)) then
                ! static
                ax_t = 0d0
                gd%FXg(tangent) = 0d0
                gd%PXg(tangent) = 0d0
            else
                ! gonna move immediately
                ax_delta = fric_coef_static*dabs(ax_n)*sign(1d0, ax_t)
                ax_t = ax_t - ax_delta
                gd%FXg(tangent) = gd%FXg(tangent) - ax_delta*gd%Mg
                gd%PXg(tangent) = gd%PXg(tangent) - ax_delta*gd%Mg*DT
            end if
    
        else
            ! kinematic frictional condition
            if (dabs(vx_t) <= fric_coef_kinematic*dabs(ax_n)*DT) then
                ! gonna stop
                ax_t = -vx_old(tangent)/DT
                gd%FXg(tangent) = -vx_old(tangent)*gd%Mg/DT
                gd%PXg(tangent) = 0d0
            else
                ! move under frictional effect
                ax_delta = fric_coef_kinematic*dabs(ax_n)*sign(1d0, vx_t)
                ax_t = ax_t - ax_delta
                gd%FXg(tangent) = gd%FXg(tangent) - ax_delta*gd%Mg
                gd%PXg(tangent) = gd%PXg(tangent) - ax_delta*gd%Mg*DT
                
            end if
        end if
        
        ! update tangent acceleration
        ax_node(tangent) = ax_t
    end do
    
    end subroutine UpdateFrictionalForce

    subroutine UpdateAccelerationByFriction(ax_try, vx_old, normal)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Calculate frictional acceleration                         -
    !-                                                                -
    !-                                                                -
    !------------------------------------------------------------------
    use ParticleData, only: ndim, DT, istep
    use GridData, only: grid_list, fric_coef_static, fric_coef_kinematic, GridNodeProperty
    implicit none
    integer, intent(in):: normal
    real(8), intent(inout):: ax_try(ndim)
    real(8), intent(in):: vx_old(ndim)
    real(8):: ax_n, ax_t, vx_t, ax_delta, tol = 1d-5, vx_new(ndim)
    integer:: tangent
    type(GridNodeProperty), POINTER :: gd
    
    ! do nothing if the acceleration point out of the surface
    if (ax_try(normal)>=0) then
        ! only applicable to the xmin, ymin surface 
        return
    end if
    
    ! compute new velocity and acceleration
    vx_new = vx_old + ax_try*DT;
    
    ! normal velocity and acceleration
    ax_n = ax_try(normal)
    do tangent = 1, ndim
        if (tangent .ne. normal) then
            ! except for normal direction,
            ! the other directions are tangent directions
            ax_t = ax_try(tangent)
            vx_t = vx_new(tangent)

            ! tangent acceleration under frictional effect
            if (dabs(vx_old(tangent))<tol) then
                ! static frictional condition
                if (dabs(ax_t)<fric_coef_static*dabs(ax_n)) then
                    ! static
                    ax_t = 0d0
                else
                    ! gonna move immediately
                    ax_delta = fric_coef_static*dabs(ax_n)*sign(1d0, ax_t)
                    ax_t = ax_t - ax_delta
                end if
            else
                ! kinematic frictional condition
                if (dabs(vx_t) <= fric_coef_kinematic*dabs(ax_n)*DT) then
                    ! gonna stop
                    ax_t = -vx_old(tangent)/DT
                else
                    ! move under frictional effect
                    ax_delta = fric_coef_kinematic*dabs(ax_n)*sign(1d0, vx_t)
                    ax_t = ax_t - ax_delta

                end if
            end if

            ! update tangent acceleration
            ax_try(tangent)= ax_t

        end if
    end do
    
    end subroutine 
    
    
#ifdef NDIM3
    subroutine Lagr_NodContact()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      1. Establishing the nodal contact criteria,          and  -
    !-      2. correct  the normal vectors                       and  -
    !-      3. Apply the contact force and adjust nodal velocities    -
    !------------------------------------------------------------------
    use  ParticleData
    use  GridData
    use  MaterialData

    implicit none

    integer:: p, n, c ! loop counter
    real(8):: nx, ny, nz, tt,tta,ttb,crit, crita, critb
    real(8):: nomforce, val_fslip, val_fstick, val_ffric, nodtolmg
    real(8):: fstick(3), fslip(3), cforce(3)
    integer:: abody,bbody

    type(GridNodeProperty), POINTER :: gd1
    type(GridNodeProperty), POINTER :: gd2
    type(ContactGridNodeProperty), POINTER :: CP1
    type(ContactGridNodeProperty), POINTER :: CP2
	type(GridNode), POINTER :: node

    tot_cont_for = 0.0 ! the total contact force between of bodies

    ! calculate contact force and adjust the nodal force and momentum
    do n = 1, nb_gridnode
        CP1 =>  CP_list(1,n)
        CP2 =>  CP_list(2,n)
        gd1 => grid_list(n, 1)
        gd2 => grid_list(n, 2)
		node => node_list(n)

        ! recalculate the nodal normal direction
        ! if normbody 0 then using average method;
        ! if 1,using abody; if 2,using bbody
        if(normbody == 0)then
            nx = CP1%ndir(1)  - CP2%ndir(1)
            ny = CP1%ndir(2)  - CP2%ndir(2)
            nz = CP1%ndir(3)  - CP2%ndir(3)
        end if

        if(normbody == 1)then
            nx = CP1%ndir(1)
            ny = CP1%ndir(2)
            nz = CP1%ndir(3)
        end if

        if(normbody == 2)then
            nx = - CP2%ndir(1)
            ny = - CP2%ndir(2)
            nz = - CP2%ndir(3)
        end if

        ! unitize normal vector
        tt = sqrt(nx*nx + ny*ny + nz*nz)
        if(tt > epsilon(tt)) then
            nx = nx / tt
            ny = ny / tt
            nz = nz / tt
        end if

        CP1%ndir(1) = nx; ! Nodal direction for contact
        CP1%ndir(2) = ny;
        CP1%ndir(3) = nz;
        CP2%ndir = -CP1%ndir

        crit = 0.0
        ! contact criteria using the unit normal vectors
        if ( gd1%Mg > CutOff .AND. gd2%Mg > CutOff) then
            ! Eq.(3.245)
            crit = (gd1%Pxg(1)*gd2%Mg - gd2%Pxg(1)*gd1%Mg)*nx +&
                (gd1%Pxg(2)*gd2%Mg - gd2%Pxg(2)*gd1%Mg)*ny +&
                (gd1%Pxg(3)*gd2%Mg - gd2%Pxg(3)*gd1%Mg)*nz
        end if

        if(crit > epsilon(crit)) then

            tt = (gd1%Mg + gd2%Mg)*Dt

            ! calculate the normal contact force
            nomforce =crit/tt   ! Eq.(3.252)

            ! for friction contact
            if(fricfa > epsilon(fricfa)) then

                ! calculate the contact force   Eq.(3.250)
                cforce = (gd1%Pxg*gd2%Mg - gd2%Pxg*gd1%Mg)/tt

                ! calculate the tangent contact force
                fstick = cforce - nomforce*CP1%ndir
                val_fstick = sqrt( fstick(1)*fstick(1) +  &
                    fstick(2)*fstick(2) + fstick(3)*fstick(3) )
                val_fslip = fricfa*abs(nomforce)
                if(val_fslip < val_fstick) then
                    cforce = nomforce*CP1%ndir + val_fslip*(fstick /val_fstick)
                end if

                ! for contact without friction
            else
                cforce = nomforce*CP1%ndir
            end if

            if (.not. node%Fix_x) then
				! add contact force to nodal force
				gd1%Fxg(1) = gd1%Fxg(1) - cforce(1)
				gd2%Fxg(1) = gd2%Fxg(1) + cforce(1)

				! adjust the nodal component by contact force
				gd1%Pxg(1) = gd1%Pxg(1) - cforce(1) * Dt
				gd2%Pxg(1) = gd2%Pxg(1) + cforce(1) * Dt
			end if

			if (.not. node%Fix_y) then
				gd1%Fxg(2) = gd1%Fxg(2) - cforce(2)
				gd2%Fxg(2) = gd2%Fxg(2) + cforce(2)

				gd1%Pxg(2) = gd1%Pxg(2) - cforce(2) * Dt
				gd2%Pxg(2) = gd2%Pxg(2) + cforce(2) * Dt
			end if

			if (.not. node%Fix_z) then
				gd1%Fxg(3) = gd1%Fxg(3) - cforce(3)
				gd2%Fxg(3) = gd2%Fxg(3) + cforce(3)

				gd1%Pxg(3) = gd1%Pxg(3) - cforce(3) * Dt
				gd2%Pxg(3) = gd2%Pxg(3) + cforce(3) * Dt
			end if

            tot_cont_for = tot_cont_for + cforce
        end if
    end do  !n

    end subroutine Lagr_NodContact
#endif
    subroutine ParticlePositionUpdate()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      1. Update particle position and velocity                  -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    use MaterialModel, only: Constitution
    use MaterialData
    implicit none

    integer:: b, p, n, parBegin, parEnd ! loop counter
    integer:: icell, inode, ix, iy, iz, comID = 1
    real(8):: xx(ndim), vx(ndim), ax(ndim), vgx(ndim), a_dmp(ndim)
    real(8):: de(6), vort(3)
    real(8):: mp_, shm, SHPn, DNDXn, DNDYn, DNDZn

    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd

    ! Update particle position and velocity
    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        ! Get comID from body
        if(contact)  comID = body_list(b)%comID

        do p = parBegin, parEnd    ! Loop over all particles (2)
            pt => particle_list(p)

            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            xx = pt%Xp;  ! Particle position at time step k
            vx = 0d0
            ax = 0d0

            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if(inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle    ! out of the computational grid

                gd => grid_list(inode, comID)
                if (gd%Mg > CutOff) then ! The nodal mass is not too small
                    SHPn = pt%SHP(n)
                    vx = vx + SHPn * (gd%PXg / gd%Mg)
                    ax = ax + SHPn * (gd%FXg / gd%Mg)
                end if

            end do ! n

            ! damping acceleration
            !a_dmp = 0d0
            !if (CurrentTime.lt.dampDuration) then
            !    a_dmp = -dampCoef*vx
            !end if
            !ax = ax + a_dmp

            ! Time integration
            pt%XX = xx + vx * DT       ! Update particle position
            if(istep == 1) then
                pt%VXp = pt%VXp + ax * DT * 0.5d0  ! Update particle velocity
            else
                pt%VXp = pt%VXp + ax * DT
            end if
            if(USF)  pt%Xp = pt%XX     ! the next particle position

        end do ! p
    end do    ! b

    end subroutine ParticlePositionUpdate

    subroutine GridMomentumMUSL()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      1. recalculate the grid node momentum by mapping          -
    !-         the updated particle information                       -
    !-      2. apply boundary condition                               -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    use MaterialModel, only: Constitution
    use MaterialData
    implicit none
    integer:: b, c, p, n, parBegin, parEnd ! loop counter
    integer:: icell, comID = 1, inode, fric_normal, idim
    real(8):: de(6), vort(3)
    real(8):: mp_, shm, SHPn
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd
    type(GridNode), POINTER :: node

    do n = 1, ndim
        grid_list%PXg(n) = 0.0d0
    end do
    
    ! Recalculate the grid node momentum
    do b = 1, nb_body        ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID = body_list(b)%comid   ! Get comID from body
        do p = parBegin, parEnd ! Loop over all particles (3)
            pt => particle_list(p)
            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            mp_ = pt%mass
            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if(inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle  ! out of the computational grid
                gd => grid_list(inode, comID)

                shm = pt%SHP(n)*mp_
                gd%PXg = gd%PXg + pt%VXp*shm
                
            end do ! n
        end do    ! p
    end do       ! b

    ! Applying essential boundary conditions
    do c = 1, nb_component
        do n = 1, nb_gridnode
            gd => grid_list(n, c)
            node => node_list(n)
            
            if (gd%Mg>CutOff) then
                gd%VXg = gd%PXg/gd%Mg
            else
                gd%VXg = 0d0
            end if
            
            ! apply friction boundary
            fric_normal = 0
            do idim = 1, ndim
                if (node%Fric(idim)) then
                    fric_normal = idim
                end if
            end do
            
            ! Frictional boundary
            if (fric_normal>0) then
                gd%PXg(fric_normal) = 0d0
                gd%VXg(fric_normal) = 0d0
            end if
            
            if (node%Fix_x) then
                gd%PXg(1) = 0.0d0
                gd%VXg(1) = 0.0d0
            end if

            if (node%Fix_y) then
                gd%PXg(2) = 0.0d0
                gd%VXg(2) = 0.0d0
            end if
            
#ifdef NDIM3
            if (node%Fix_z) then
                gd%PXg(3) = 0.0d0
                gd%VXg(3) = 0.0d0
            end if
#endif
        end do !n
    end do    ! c
    end subroutine GridMomentumMUSL

    subroutine ApplyBoundaryConditions()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      1. apply boundary condition                               -
    !------------------------------------------------------------------
    use GridData
    use ParticleData, only: nb_component
    implicit none

    integer:: n, c ! loop counter
    type(GridNodeProperty), POINTER :: gd
    type(GridNode), POINTER :: node

    do c = 1, nb_component
        do n = 1, nb_gridnode
            gd => grid_list(n, c)
            node => node_list(n)
            ! Apply boundary conditions on computational grid
            if (node%Fix_x) then  ! Grid node n is fixed in x direction
                gd%PXg(1) = 0.0d0
                gd%FXg(1) = 0.0d0
                gd%AXg(1) = 0.0d0
                gd%VXg(1) = 0.0d0
            end if

            if (node%Fix_y) then  ! Grid node n is fixed in y direction
                gd%PXg(2) = 0.0d0
                gd%FXg(2) = 0.0d0
                gd%AXg(2) = 0.0d0
                gd%VXg(2) = 0.0d0
            end if
            
#ifdef NDIM3
            if (node%Fix_z) then  ! Grid node n is fixed in z direction
                gd%PXg(3) = 0.0d0
                gd%FXg(3) = 0.0d0
                gd%AXg(3) = 0.0d0
                gd%VXg(3) = 0.0d0
            end if
#endif
        end do ! n
    end do    ! c

    end subroutine ApplyBoundaryConditions
   
    subroutine StrainApproximationByNode()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      strain approximation by node                              -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    implicit none
    
    integer:: b, c, p, n, parBegin, parEnd ! loop counter
    integer:: icell, comID = 1, inode
    real(8):: strain_new(6), r1d3
    real(8):: mp_, temp, dilation, dilation_new
    data r1d3 /0.333333333333333333d0/
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd
    
    ! reset nodal dilation summation
    grid_list(:,:)%dilation = 0d0
    
    ! compute nodal dilation summation
    do b = 1, nb_body        ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID = body_list(b)%comid   ! Get comID from body
        do p = parBegin, parEnd ! Loop over all particles (3)
            pt => particle_list(p)
            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            mp_ = pt%mass
            dilation = sum(pt%trial_elastic_strain(1:3))
            pt%dilation = dilation

            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if(inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle  ! out of the computational grid
                gd => grid_list(inode, comID)
                
                if (gd%Mg > CutOff) then
                    temp = dilation*pt%SHP(n)*mp_/gd%Mg
                    gd%dilation = gd%dilation + temp
                end if 
            end do ! n
        end do ! p
    end do ! b
 
    
    ! approximate point trial elastic strain
    do b = 1, nb_body        ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID = body_list(b)%comid   ! Get comID from body
        do p = parBegin, parEnd ! Loop over all particles (3)
            pt => particle_list(p)
            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            ! DO NOT approximate if the cell has only 1 or 2 particles
            !if (cell_list(icell)%NPInCell<2) then
            !    cycle
            !end if
            
            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            dilation_new = 0d0
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if(inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle  ! out of the computational grid
                gd => grid_list(inode, comID)
                if (gd%Mg > CutOff) then
                    dilation_new = dilation_new + gd%dilation*pt%SHP(n)
                end if
            end do ! n
            
            ! compute new trial elastic strain
            dilation = pt%dilation
            strain_new = pt%trial_elastic_strain
            strain_new(1) = strain_new(1) + (dilation_new - dilation)*r1d3
            strain_new(2) = strain_new(2) + (dilation_new - dilation)*r1d3
            strain_new(3) = strain_new(3) + (dilation_new - dilation)*r1d3
            
            ! update point dilation and trial elastic strain
            pt%dilation = dilation_new
            pt%trial_elastic_strain = strain_new
        end do ! p
    end do ! b
    
    end subroutine StrainApproximationByNode
    
    subroutine StressApproximationByNode()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      stress approximation by node                              -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    use MaterialData
    implicit none
    
    integer:: b, c, p, n, parBegin, parEnd ! loop counter
    integer:: icell, comID = 1, matid, inode
    real(8):: temp, pressure_new, density, pressure
    real(8):: r1d3
    data r1d3 /0.333333333333333333d0/
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd
    
    ! reset nodal pressure summation
    grid_list%pressure = 0d0
    
    ! compute nodal pressure summation
    do b = 1, nb_body        ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID = body_list(b)%comid   ! Get comID from body
        do p = parBegin, parEnd ! Loop over all particles (3)
            pt => particle_list(p)
            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            pressure = pt%SM
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if(inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle  ! out of the computational grid
                gd => grid_list(inode, comID)
                
                if (gd%Mg > CutOff) then
                    temp = pressure*pt%SHP(n)*pt%mass/gd%Mg
                    gd%pressure = gd%pressure + temp
                end if 
            end do ! n
        end do ! p
    end do ! b

    ! approximate point stress
    do b = 1, nb_body        ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        matid = body_list(b)%mat
        density = mat_list(matid)%Density
        if(contact) comID = body_list(b)%comid   ! Get comID from body
        do p = parBegin, parEnd ! Loop over all particles (3)
            pt => particle_list(p)
            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            pressure_new = 0d0
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if(inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle  ! out of the computational grid
                gd => grid_list(inode, comID)
                if (gd%Mg > CutOff) then
                    pressure_new = pressure_new + gd%pressure*pt%SHP(n)
                end if
            end do ! n
           

            ! update point pressure and deviatoric stress
            pt%SM = pressure_new
        end do ! p
    end do ! b

    end subroutine StressApproximationByNode

    subroutine StressApproximationByNodeQUAD()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      stress approximation by node                              -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    use MaterialData
    implicit none

    integer:: b, c, p, n, parBegin, parEnd ! loop counter
    integer:: icell, comID = 1, matid, inode
    real(8):: temp, pressure_new, density, pressure
    real(8):: r1d3
    data r1d3 /0.333333333333333333d0/
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd

    ! reset nodal pressure summation
    grid_list%pressure = 0d0

    ! compute nodal pressure summation
    do b = 1, nb_body        ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID = body_list(b)%comid   ! Get comID from body
        do p = parBegin, parEnd ! Loop over all particles (3)
            pt => particle_list(p)
            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            pressure = pt%SM*pt%VOL
            do n = 1, nCellNode
                inode = pt%InflNode(n)
                if(inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle  ! out of the computational grid
                gd => grid_list(inode, comID)

                if (gd%Mg > CutOff) then
                    temp = pressure*pt%SHP_QUAD(n)/gd%Mg
                    gd%pressure = gd%pressure + temp
                end if
            end do ! n
        end do ! p
    end do ! b

    ! approximate point stress
    do b = 1, nb_body        ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        matid = body_list(b)%mat
        density = mat_list(matid)%Density
        if(contact) comID = body_list(b)%comid   ! Get comID from body
        do p = parBegin, parEnd ! Loop over all particles (3)
            pt => particle_list(p)
            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle
            
            ! DO NOT approximate if the cell has only 1 or 2 particles
            !if (cell_list(icell)%NPInCell<2) then
            !    cycle
            !end if
            
            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            pressure_new = 0d0
            do n = 1, nCellNode
                inode = pt%InflNode(n)
                if(inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle  ! out of the computational grid
                gd => grid_list(inode, comID)
                if (gd%Mg > CutOff) then
                    pressure_new = pressure_new + gd%pressure*pt%SHP_QUAD(n)
                end if
            end do ! n
            pressure_new = pressure_new*pt%mass/pt%VOL
            
            ! update point pressure and deviatoric stress
            pt%SM = pressure_new
        end do ! p
    end do ! b
        
    end subroutine StressApproximationByNodeQUAD

    subroutine StrainIncrementApproximationByNode()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      strain approximation by node                              -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    implicit none

    integer:: b, c, p, n, parBegin, parEnd ! loop counter
    integer:: icell, comID = 1, inode
    real(8):: strain_new(6), r1d3
    real(8):: mp_, temp, dilation, dilation_new
    data r1d3 /0.333333333333333333d0/
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd

    ! reset nodal dilation summation
    grid_list(:,:)%dilation = 0d0

    ! compute nodal dilation summation
    do b = 1, nb_body        ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID = body_list(b)%comid   ! Get comID from body
        do p = parBegin, parEnd ! Loop over all particles (3)
            pt => particle_list(p)
            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            mp_ = pt%mass
            dilation = sum(pt%strain_increment(1:3))
            pt%dilation = dilation

            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if(inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle  ! out of the computational grid
                gd => grid_list(inode, comID)

                if (gd%Mg > CutOff) then
                    temp = dilation*pt%SHP(n)*mp_/gd%Mg
                    gd%dilation = gd%dilation + temp
                end if
            end do ! n
        end do ! p
    end do ! b

    ! approximate point trial elastic strain
    do b = 1, nb_body        ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID = body_list(b)%comid   ! Get comID from body
        do p = parBegin, parEnd ! Loop over all particles (3)
            pt => particle_list(p)
            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            ! DO NOT approximate if the cell has only 1 or 2 particles
            !if (cell_list(icell)%NPInCell<2) then
            !    cycle
            !end if

            ! Loop over all grid nodes of the hexhedron
            !  in which particle p is located
            dilation_new = 0d0
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if(inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle  ! out of the computational grid
                gd => grid_list(inode, comID)
                if (gd%Mg > CutOff) then
                    dilation_new = dilation_new + gd%dilation*pt%SHP(n)
                end if
            end do ! n

            ! compute new trial elastic strain
            dilation = pt%dilation
            strain_new = pt%strain_increment
            strain_new(1) = strain_new(1) + (dilation_new - dilation)*r1d3
            strain_new(2) = strain_new(2) + (dilation_new - dilation)*r1d3
            strain_new(3) = strain_new(3) + (dilation_new - dilation)*r1d3

            ! update point dilation and trial elastic strain
            pt%dilation = dilation_new
            pt%strain_increment = strain_new
        end do ! p
    end do ! b

    end subroutine StrainIncrementApproximationByNode

    subroutine ParticleSHPCompute()
    !-------------------------------------------------------------------
    !-  Purpose                                                        -
    !-      1. Find surface nodes % cells & particles                  -
    !-------------------------------------------------------------------

    use ParticleData
    use GridData
    use MaterialData
    implicit none
    integer:: b, p, n, parBegin, parEnd, inode, ix, iy ! loop counter
    integer:: icell, comID = 1, ThisCellNodes(nCellNode), node1
    logical:: surfaceCellNotFound
    real(8):: x(ndim), sx(nCellNode), sy(nCellNode)
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd
    type(GridCellProperty), POINTER :: gc

    ! Compute shape function for all particles
    cell_list%InOutProp = 0
    cell_list%NPInCell = 0
    grid_list(:,:)%InOutProp = 0
    do b = 1, nb_body     ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        do p = parBegin, parEnd    ! Loop over all particles (1)
            pt => particle_list(p)
            icell = InWhichCell(pt%Xp)
            pt%icell = icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            ! nodes of the cell
            if (GIMP) then
                call FindInflNode(p,icell)
                call NShape_GIMP(p)
                ! quad
                !call NShape_QUAD(CellsNode(icell,1),p,2)
            else
                pt%nb_InflNode = nCellNode
				pt%InflNode(1:nCellNode) = CellsNode(icell,:)
                call NShape(CellsNode(icell,1),p,2)
                pt%SHP_QUAD(1:nCellNode) = pt%SHP(1:nCellNode)
            end if

            ! 0 -- cell out of the body
            ! 1 -- cell in the body 
            ! 2 -- surface cell
            cell_list(icell)%InOutProp = 1 
            
            ! particle number in the cell
            cell_list(icell)%NPInCell = cell_list(icell)%NPInCell + 1 
            
        end do
    end do
    
    ! Find surface cells and nodes
    call FindSurfaceCN()
    
    end subroutine ParticleSHPCompute
    
    
    
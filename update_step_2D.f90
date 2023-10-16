
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
#ifdef NDIM2
    subroutine ParticleStrainUpdate()
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Calculate the strain rate and spin tensor                 -
    !-                                                                -
    !------------------------------------------------------------------
    use ParticleData
    use GridData
    use MaterialModel, only: Constitution
    use MaterialData
    implicit none

    integer:: b, p, n, parBegin, parEnd ! loop counter
    integer:: icell, inode, ix, iy, iz, comID = 1
    real(8):: xx(ndim), vx(ndim), ax(ndim), vgx(ndim)
    real(8):: de(6), vort(3)
    real(8):: mp_, shm, SHPn, DNDXn, DNDYn
    real(8):: B11, B12, B21, B22, B31, B32
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd
    real(8):: dFG(3,3), DGnew(3, 3), Jacobian

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
                    vgx = gd%VXg    ! Grid nodal velocity
                    
                    DNDXn = pt%DNDX(n)
                    DNDYn = pt%DNDY(n)
                    de(1) = de(1) + DNDXn*vgx(1)            ! D11
                    de(2) = de(2) + DNDYn*vgx(2)            ! D22
                    de(6) = de(6) + (DNDXn*vgx(2) + DNDYn*vgx(1)) ! D12
                    vort(3) = vort(3) + (DNDXn*vgx(2) - DNDYn*vgx(1)) ! omega12
                    
                    ! deformation gradient increment
                    dFG(1,1:ndim) = dFG(1,1:ndim) + DNDXn*vgx
                    dFG(2,1:ndim) = dFG(2,1:ndim) + DNDYn*vgx
                end if
            end do ! n
            !
            de = de * DT
            vort = vort * DT / 2d0
            
            ! deformation gradient increment 
            dFG = dFG * DT
            dFG(1,1) = dFG(1,1) + 1d0
            dFG(2,2) = dFG(2,2) + 1d0
            dFG(3,3) = 1d0
            
            ! new DG = deltaDG * DG
            DGnew = matmul(dFG, pt%DG)
            Jacobian = DGnew(1,1)*DGnew(2,2) - DGnew(1,2)*DGnew(2,1)
             
            pt%DFG = dFG
            pt%DG = DGnew
            pt%Jacobian = Jacobian
            pt%strain_increment = de
            pt%vort = vort
            ! small strain
            pt%trial_elastic_strain = pt%elastic_strain + de
            
            ! compute trial elastic strain
            if (FiniteStrain) then
                ! finite strain
                !call compute_trial_elastic_log_strain(pt%elastic_strain, pt%DFG, pt%trial_elastic_strain)
            end if

        end do !p
    end do    !b

    end subroutine ParticleStrainUpdate
	subroutine particlepwp()
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
    real(8):: vol_, SHPn, shv, pp, INVHI(3,3), xp_, yp_
	real(8):: pwp

    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd
    do n = 1, 3
        grid_list%AI(n) = 0.0d0    ! Nodal forces
    end do
	do n = 1, 3
        grid_list%AII(n) = 0.0d0    ! Nodal forces
    end do
    grid_list%HI(1,1) = 0.0d0 
	grid_list%HI(1,2) = 0.0d0 
	grid_list%HI(1,3) = 0.0d0 
	grid_list%HI(2,1) = 0.0d0 
	grid_list%HI(2,2) = 0.0d0 
	grid_list%HI(2,3) = 0.0d0 
	grid_list%HI(3,1) = 0.0d0 
	grid_list%HI(3,2) = 0.0d0 
	grid_list%HI(3,3) = 0.0d0 
	

	
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
            ! Loop over the grid nodes of the hexhedron
            !   in which the particle is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                ! out of the computational grid
                if (inode .gt. nb_gridnode .or. &
                    inode .le. 0) cycle

                gd => grid_list(inode, comID)

                SHPn = pt%SHP(n)
                gd%HI(1,1)= gd%HI(1,1)+1*SHPn*vol_
				gd%HI(1,2)= gd%HI(1,2)+pt%XX(1)*SHPn*vol_
				gd%HI(1,3)= gd%HI(1,3)+pt%XX(2)*SHPn*vol_
				
				gd%HI(2,1)= gd%HI(2,1)+pt%XX(1)*SHPn*vol_
				gd%HI(2,2)= gd%HI(2,2)+pt%XX(1)*pt%XX(1)*SHPn*vol_
				gd%HI(2,3)= gd%HI(2,3)+pt%XX(1)*pt%XX(2)*SHPn*vol_
				
				gd%HI(3,1)= gd%HI(3,1)+pt%XX(2)*SHPn*vol_
				gd%HI(3,2)= gd%HI(3,2)+pt%XX(1)*pt%XX(2)*SHPn*vol_
				gd%HI(3,3)= gd%HI(3,3)+pt%XX(2)*pt%XX(2)*SHPn*vol_

            end do !n			
        end do !p
    end do    !b
	
	do b = 1, nb_body     ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID = body_list(b)%comID ! Get comID from body
        do p = parBegin, parEnd    ! Loop over all particles (1)
            pt => particle_list(p)

            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            vol_ = pt%VOL
			pp =pt%KW*(pt%strain_increment(1)+pt%strain_increment(2)+pt%strain_increment(3))
            ! Loop over the grid nodes of the hexhedron
            !   in which the particle is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                ! out of the computational grid
                if (inode .gt. nb_gridnode .or. &
                    inode .le. 0) cycle

                gd => grid_list(inode, comID)
                
                SHPn = pt%SHP(n)
               
								
                gd%AII(1)=gd%AII(1)+SHPn*pp*vol_
				
				gd%AII(2)=gd%AII(2)+SHPn*pp*pt%XX(1)*vol_
                 
			    gd%AII(3)=gd%AII(3)+SHPn*pp*pt%XX(2)*vol_
            end do !n	
        end do !p
    end do    !b
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
            ! Loop over the grid nodes of the hexhedron
            !   in which the particle is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                ! out of the computational grid
                if (inode .gt. nb_gridnode .or. &
                    inode .le. 0) cycle

                gd => grid_list(inode, comID)
                
				call M33INV(gd%HI, INVHI)
                gd%AI(1)=gd%AII(1)*INVHI(1,1)+gd%AII(2)*INVHI(1,2)+gd%AII(3)*INVHI(1,3)
				
				gd%AI(2)= gd%AII(1)*INVHI(2,1)+gd%AII(2)*INVHI(2,2)+gd%AII(3)*INVHI(2,3)
                 
			    gd%AI(3)=gd%AII(1)*INVHI(3,1)+gd%AII(2)*INVHI(3,2)+gd%AII(3)*INVHI(3,3)
            end do !n	
        end do !p
    end do    !b
	
	do b = 1, nb_body     ! Loop over all bodies
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_End
        if(contact) comID = body_list(b)%comID ! Get comID from body
        do p = parBegin, parEnd    ! Loop over all particles (1)
            pt => particle_list(p)

            icell = pt%icell
            ! Particle p is out of the computational region
            if (icell < 0) cycle

            vol_ = pt%VOL
			pwp  = 0.0d0
            ! Loop over the grid nodes of the hexhedron
            !   in which the particle is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                ! out of the computational grid
                if (inode .gt. nb_gridnode .or. &
                    inode .le. 0) cycle

                gd => grid_list(inode, comID)
                
                SHPn = pt%SHP(n)
	            pwp= pwp+SHPn*(gd%AI(1)+gd%AI(2)*pt%XX(1)+gd%AI(3)*pt%XX(2))
            end do !n
			    pt%pwpa=pwp
        end do !p
    end do    !b
    end subroutine particlepwp

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
    real(8):: sxx, syy, szz, sxy
    real(8):: fx(ndim), f_int(ndim), f_ext(ndim), f_dmp(ndim), mp_, vol_
    real(8):: shm, SHPn, DNDXn, DNDYn, load_factor, B11, B12, B21, B22, B31, B32
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd
    type(ContactGridNodeProperty), POINTER :: CP

    ! Calculate the grid nodal forces only

    ! Reset nodal forces
    do n = 1, ndim
        grid_list%FXg(n) = 0.0d0    ! Nodal forces
    end do
    
    if(contact) then
        do n = 1, ndim
            CP_list%ndir(n) = 0.0d0
        end do
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
            sxy = pt%SDxy
            
            ! External forces
            fx = pt%FXp*load_factor
            if (Gravity) then
                fx = fx + mp_ * (body_list(b)%Gravp) * load_factor
            end if

            ! Loop over the grid nodes of the hexhedron
            !  in which the particle is located
            do n = 1, pt%nb_InflNode
                inode = pt%InflNode(n)
                if (inode .gt. nb_gridnode .or. inode .le. 0) &
                    cycle  ! out of the computational grid

                gd => grid_list(inode, comID)

                SHPn = pt%SHP(n)
                DNDXn = pt%DNDX(n)
                DNDYn = pt%DNDY(n)

                f_int(1) = - (sxx*DNDXn + sxy*DNDYn)*vol_
                f_int(2) = - (sxy*DNDXn + syy*DNDYn)*vol_

                f_ext = fx*SHPn

                f_dmp = 0.0d0
                if (CurrentTime.lt.dampDuration) then
                    f_dmp = -dampCoef*abs(f_int + f_ext)*sign(1.0d0,pt%VXp)
                end if
                
                gd%FXg = gd%FXg + f_int + f_ext + f_dmp ! nodal force

                if(contact) then
                    CP => CP_list(comID, inode)
                    CP%ndir(1) = CP%ndir(1) + DNDXn*mp_
                    CP%ndir(2) = CP%ndir(2) + DNDYn*mp_
                end if

            end do !n

        end do !p
    end do    !b
    end subroutine GridMomentumUpdate
    
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
    real(8):: fstick(ndim), fslip(ndim), cforce(ndim)
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
            !nz = CP1%ndir(3)  - CP2%ndir(3)
        end if

        if(normbody == 1)then
            nx = CP1%ndir(1)
            ny = CP1%ndir(2)
            !nz = CP1%ndir(3)
        end if

        if(normbody == 2)then
            nx = - CP2%ndir(1)
            ny = - CP2%ndir(2)
            !nz = - CP2%ndir(3)
        end if

        ! unitize normal vector
        !tt = sqrt(nx*nx + ny*ny + nz*nz)
        tt = sqrt(nx*nx + ny*ny)
        if(tt > epsilon(tt)) then
            nx = nx / tt
            ny = ny / tt
            !nz = nz / tt
        end if

        CP1%ndir(1) = nx; ! Nodal direction for contact
        CP1%ndir(2) = ny;
        !CP1%ndir(3) = nz;
        CP2%ndir = -CP1%ndir

        crit = 0.0
        ! contact criteria using the unit normal vectors
        if ( gd1%Mg > CutOff .AND. gd2%Mg > CutOff) then
            ! Eq.(3.245)
            !crit = (gd1%Pxg(1)*gd2%Mg - gd2%Pxg(1)*gd1%Mg)*nx +&
            !    (gd1%Pxg(2)*gd2%Mg - gd2%Pxg(2)*gd1%Mg)*ny +&
            !    (gd1%Pxg(3)*gd2%Mg - gd2%Pxg(3)*gd1%Mg)*nz
            !
            crit = (gd1%Pxg(1)*gd2%Mg - gd2%Pxg(1)*gd1%Mg)*nx +&
                (gd1%Pxg(2)*gd2%Mg - gd2%Pxg(2)*gd1%Mg)*ny
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
                !val_fstick = sqrt( fstick(1)*fstick(1) +  &
                !    fstick(2)*fstick(2) + fstick(3)*fstick(3) )
                
                val_fstick = sqrt( fstick(1)*fstick(1) +  &
                    fstick(2)*fstick(2) )
                
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

            tot_cont_for = tot_cont_for + cforce
        end if
    end do  !n

    end subroutine Lagr_NodContact
	!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!***********************************************************************************************************************************

    SUBROUTINE M33INV (A, AINV)

      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV


      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)
     IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         RETURN
      END IF


      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

   

      RETURN

      END SUBROUTINE M33INV
#endif
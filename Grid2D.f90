
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
    ! -  Background grid procedures                                    -
    ! -                                                                -
    ! -  NOTE: only 8 node cube element is available                   -
    ! -        only cuboid computational region is available           -
    ! -  computational region is:                                      -
    ! -   [SpanX(1),SpanX(2)]*[SpanY(1),SpanY(2)]*[SpanZ(1),SpanZ(2)]  -
    ! -  Grid,Cell(element) infomation will be computed by program     -
    ! -  GridNode Numbering is keypoint of searching algorithm         -
    ! ------------------------------------------------------------------
#ifdef NDIM2
    module GridData

    type GridNode
        real(8):: Xg(2)     ! grid node coordinate
        logical:: Fix_x, Fix_y    ! BC
        logical:: Fric(2)   ! frictional boundary
    end type GridNode

    type GridNodeProperty
        real(8):: Mg        ! mass on grid node
        real(8):: PXg(2)    ! momentum on grid node
        real(8):: FXg(2)    ! internal/external force on gride node
        ! for strain/stress approximation
        real(8):: AXg(2)    ! acceleration
        real(8):: VXg(2)    ! velocity
        real(8):: VXg_old(2)! old velocity
        real(8):: dilation, pressure
		! for pore-water pressure
		real(8):: HI(3,3)
		real(8)::AII(3)
		real(8):: AI(3)
        Integer:: InOutProp
        real(8):: KI        ! mass on grid node
    end type GridNodeProperty
    
    type GridCellProperty
        Integer:: NPInCell, InOutProp
    end type GridCellProperty

    type ContactGridNodeProperty
        ! the normal direction of contact grid node
        real(8):: ndir(2)
        ! the tangential unit vetors of contact grid node
        real(8):: sdir(2)
    end type ContactGridNodeProperty

    type(GridCellProperty), target, allocatable:: cell_list(:)
    type(GridNodeProperty), target, allocatable:: grid_list(:,:)
    type(GridNode), target, allocatable:: node_list(:)
    type(ContactGridNodeProperty), target, allocatable:: CP_list(:,:)
    real(8):: fricfa = 0.0    ! the frictional coefficient
    integer:: normbody = 0    ! the flag of computaional normal
    ! the flag of contact type: 0-no,1-langrange,2-penalty
    integer:: contact_type = 0
    ! the total contact force between of bodies
    real(8):: tot_cont_for(2)

    real(8):: SpanX(2) = 0.0d0    ! computational region
    real(8):: SpanY(2) = 0.0d0
    real(8):: SpanZ(2) = 0.0d0

    real(8):: DCell  = 0.0d0     ! Grid node interval

    real(8):: CutOff = 0.0d0     ! Grid mass cutoff value

    integer:: nb_gridnode = 0    ! number of gridnodes
    integer:: FixS(6) = 0

    ! Number of cells
    integer:: NumCell=0, NumCellx=0, NumCelly=0
    integer:: NGx, NGy
    ! NumCell * 8 - Define computational grid
    integer, allocatable:: CellsNode(:,:)
    integer, parameter:: nCellNode = 4
    !integer:: nb_InflNode = 4
    !integer:: InflNode(9)    ! influence node list for each particle
    !real(8):: rpg(9,2)       ! distance between particle and node

    real(8):: iJacobi, iJacobi4, iDCell ! 1/Jacobi
    ! shape function and its derivative
    !real(8), allocatable:: SHP(:), DNDX(:), DNDY(:)

    ! sign constant used by NShape
    integer,parameter:: SNX(4) = (/-1, 1, 1, -1/), &
        SNY(4) = (/-1, -1, 1, 1/)
    
    ! shape function and its spatial derivatives at element center
    real(8):: SHPC(4) = (/0.25d0, 0.25d0, 0.25d0, 0.25d0/)
    real(8):: DNDXC(4), DNDYC(4)
    
    real(8):: fric_coef_static = 0.2d0, fric_coef_kinematic = 0.2d0

    contains

    integer function InWhichCell(xx)
    ! -----------------------------------------------------------------
    ! - Purpose                                                       -
    ! -    Determine which cell the point (xx,yy,zz) is located in    -
    ! -                                                               -
    ! - Input                                                         -
    ! -    xx(3) - Coordinates of a point                             -
    ! -                                                               -
    ! - Return values                                                 -
    ! -    >0 : Cell number in which the point is located (MP)        -
    ! -    <0 : point out off the cell                                -
    ! -----------------------------------------------------------------
    implicit none
    real(8), intent(in):: xx(2)
    integer ix, iy

    if (xx(1)<SpanX(1) .or. xx(1)>SpanX(2) .or. xx(2)<SpanY(1) .or. &
        xx(2)>SpanY(2) ) then
    InWhichCell = -1
    return
    end if

    ix = int((xx(1)-SpanX(1))/DCell) + 1
    iy = int((xx(2)-SpanY(1))/DCell) + 1
    InWhichCell = (iy - 1)*NumCellx + ix

    if(InWhichCell.gt.NumCell .or. InWhichCell.le.0) then
        InWhichCell = -1
    end if

    end function InWhichCell


    subroutine SetGridData()
    ! -----------------------------------------------------------------
    ! - Purpose                                                       -
    ! -    Create computational grid                                  -
    ! -----------------------------------------------------------------
    use ParticleData
    use FFI

    implicit none

    integer:: i, j, ix, iy, icell, inode     ! loop counter
    integer:: NP2, Node1, Node5
    character(len=2) Velo
    real(8):: TempV, spx, spy
    real(8):: mat_, mp_

    spx = SpanX(2)-SpanX(1)
    spy = SpanY(2)-SpanY(1)

    if(DCell.eq.0) then
        stop '*** Error *** DCell must be defined !'
    end if

    if(spx.le.0 .or. spy.le.0 ) then
        stop '*** Error *** SPX/SPY must be defined !'
    end if

    NumCellx = int(spx/DCell + 0.5)
    NumCelly = int(spy/DCell + 0.5)

    SpanX(2) = SpanX(1) + NumCellx*DCell
    SpanY(2) = SpanY(1) + NumCelly*DCell
    NumCell = NumCellx * NumCelly

    NGx = NumCellx + 1
    NGy = NumCelly + 1

    nb_gridnode = NGx*NGy

    print *, 'Number of grid nodes = ', nb_gridnode
    write(iomsg,*)
    write(iomsg,"(a14,i10)") 'Number of grid nodes = ', nb_gridnode
    write(iomsg,"(a14,2i4)") 'cells (x,y) ', &
        numcellx,numcelly

    allocate(grid_list(nb_gridnode, nb_component))
    allocate(node_list(nb_gridnode))
    allocate(cell_list(NumCell))

    cell_list%InOutProp = 0
    cell_list%NPInCell = 0
    
    node_list%Fix_x = .false.
    node_list%Fix_y = .false.
    node_list%Fric(1) = .false.
    node_list%Fric(2) = .false.
    ! create grid node info

    do iy = 1, NGy
        do ix = 1, NGx
            i = (iy - 1)*NGx + ix
            node_list(i)%Xg(1) = (ix - 1)*DCell + SpanX(1)
            node_list(i)%Xg(2) = (iy - 1)*DCell + SpanY(1)

            if ((ix==1.and.FixS(1)==1).or.(ix==NGx.and.FixS(2)==1).or. &
                (iy==1.and.FixS(3)==1).or.(iy==NGy.and.FixS(4)==1)) then
                node_list(i)%Fix_x = .true.
                node_list(i)%Fix_y = .true.
            end if

            if ((ix==1.and.FixS(1)==2).or.(ix==NGx.and.FixS(2)==2)) then
                node_list(i)%Fix_x = .true.
            end if
            
            if ((iy==1.and.FixS(3)==2).or.(iy==NGy.and.FixS(4)==2)) then
                node_list(i)%Fix_y = .true.
            end if
            
            ! frictional surface ix = 1 
            if ((ix==1.and.FixS(1)==3)) then
                node_list(i)%Fric(1) = .true.
                node_list(i)%Fix_x = .true.
            end if

            ! frictional surface iy = 1 
            if ((iy==1.and.FixS(3)==3)) then
                node_list(i)%Fric(2) = .true.
                node_list(i)%Fix_y = .true.
            end if
                
        end do
    end do

    ! create CellsNode
    allocate(CellsNode(NumCell,4))

    ! loop over every cell to create CellsNode Matrix
    do iy = 1, NumCelly
        do ix = 1, NumCellx
            i = (iy - 1)*NumCellx + ix
            Node1 = (iy - 1)*NGx + ix
            CellsNode(i,1) = Node1
            CellsNode(i,2) = Node1 + 1
            CellsNode(i,3) = Node1 + 1 + NGx
            CellsNode(i,4) = Node1 + NGx
        end do
    end do

    ! set parameter used in NShape
    !if (GIMP) then
    !    allocate(SHP(9), DNDX(9), DNDY(9))
    !else
    !    allocate(SHP(4), DNDX(4), DNDY(4))
    !end if

    iJacobi  = 2d0/DCell    !Jacobi = DCell/2d0, iJacobi = 1/Jacobi
    iJacobi4 = iJacobi * 0.25d0    !iJacobi4 = iJacobi/4 for 2D
    iDCell = 1.0d0/DCell

    DNDXC = SNX*iJacobi4    ! derivative of shape function at cell center
    DNDYC = SNY*iJacobi4
    
    end subroutine SetGridData

    subroutine SetContact_GridNodeData()
    ! -----------------------------------------------------------------
    ! - Purpose                                                       -
    ! -    Create computational grid for simulating contact problem   -
    ! -----------------------------------------------------------------
    use ParticleData, only:nb_component
    use FFI

    implicit none

    allocate(CP_list(nb_component, nb_gridnode))

    end subroutine SetContact_GridNodeData

    subroutine NShape(node1, p, ider)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Evaluate the shape functions and/or their derivatives     -
    !-      at particle p associated with nodes of the cell in        -
    !-      which the particle p is located                           -
    !-  Inputs                                                        -
    !-      node1 - number of the first node of the cell in which the -
    !-              particle p is located                             -
    !-      p     - particle number                                   -
    !-      ider  - flag for shape function and derivative calculation-
    !-              0 - shape functions only                          -
    !-              1 - derivatives only                              -
    !-              2 - both shape functions and their derivatives    -
    !-  Outputs                                                       -
    !-      SHP(8)   - value of shape functions                       -
    !-      DNDX(8)  - value of derivative with respect to            -
    !-                 X of shape function                            -
    !-      DNDY(8)  - value of derivative with respect to            -
    !-                 Y of shape function                            -
    !-      DNDZ(8)  - value of derivative with respect to            -
    !-                 Z of shape function                            -
    !------------------------------------------------------------------
    use ParticleData
    implicit none

    integer, intent(in):: node1, p, ider
    real(8):: x(2) ! nature coordinate ( -1 < x,y,z < 1 )
    real(8):: sx(4), sy(4)
    type(Particle), POINTER :: pt
    type(GridNode), POINTER :: node

    node => node_list(node1)
    pt => particle_list(p)

    x = (pt%Xp - node%Xg)*iJacobi - 1d0

    sx = SNX*x(1) + 1d0    ! 1 + xi(i)  * xi
    sy = SNY*x(2) + 1d0    ! 1 + eta(i) * eta

    if (ider .NE. 1) pt%SHP(1:4) = sx * sy * 0.25d0

    if (ider .NE. 0) then
        pt%DNDX(1:4) = SNX * sy * iJacobi4
        pt%DNDY(1:4) = SNY * sx * iJacobi4
    end if

    end subroutine NShape

    subroutine NShape_QUAD(node1, p, ider)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Evaluate the shape functions and/or their derivatives     -
    !-      at particle p associated with nodes of the cell in        -
    !-      which the particle p is located                           -
    !-  Inputs                                                        -
    !-      node1 - number of the first node of the cell in which the -
    !-              particle p is located                             -
    !-      p     - particle number                                   -
    !-      ider  - flag for shape function and derivative calculation-
    !-              0 - shape functions only                          -
    !-              1 - derivatives only                              -
    !-              2 - both shape functions and their derivatives    -
    !-  Outputs                                                       -
    !-      SHP(8)   - value of shape functions                       -
    !-      DNDX(8)  - value of derivative with respect to            -
    !-                 X of shape function                            -
    !-      DNDY(8)  - value of derivative with respect to            -
    !-                 Y of shape function                            -
    !-      DNDZ(8)  - value of derivative with respect to            -
    !-                 Z of shape function                            -
    !------------------------------------------------------------------
    use ParticleData
    implicit none

    integer, intent(in):: node1, p, ider
    real(8):: x(2) ! nature coordinate ( -1 < x,y,z < 1 )
    real(8):: sx(4), sy(4)
    type(Particle), POINTER :: pt
    type(GridNode), POINTER :: node

    node => node_list(node1)
    pt => particle_list(p)

    x = (pt%Xp - node%Xg)*iJacobi - 1d0

    sx = SNX*x(1) + 1d0    ! 1 + xi(i)  * xi
    sy = SNY*x(2) + 1d0    ! 1 + eta(i) * eta

    if (ider .NE. 1) pt%SHP_QUAD(1:4) = sx * sy * 0.25d0

    !if (ider .NE. 0) then
    !    pt%DNDX_QUAD(1:4) = SNX * sy * iJacobi4
    !    pt%DNDY_QUAD(1:4) = SNY * sx * iJacobi4
    !end if

    end subroutine NShape_QUAD

    subroutine NShape_GIMP(p)
    !------------------------------------------------------------------
    !-  Purpose                                                       -
    !-      Evaluate the shape functions and/or their derivatives     -
    !-      at particle p associated with p's influence nodes (GIMP)  -
    !-  Inputs                                                        -
    !-      p     - particle number                                   -
    !-  Outputs                                                       -
    !-      SHP(27)   - value of shape functions                      -
    !-      DNDX(27)  - value of derivative with respect to           -
    !-                  X of shape function                           -
    !-      DNDY(27)  - value of derivative with respect to           -
    !-                  Y of shape function                           -
    !-      DNDZ(27)  - value of derivative with respect to           -
    !-                  X of shape function                           -
    !------------------------------------------------------------------
    use ParticleData, only:Particle, particle_list
    implicit none
    integer, intent(in):: p
    integer:: i, j, nb_InflNode, InflNode(9)
    real(8):: rpg(9,2), sh(2), dn(2), r(2)
    type(Particle), POINTER :: pt

    sh = 0.0d0! must initialize here instead of in variable statement block
    dn = 0.0d0
    pt => particle_list(p)
    rpg = pt%rpg
    nb_InflNode = pt%nb_InflNode
    InflNode = pt%InflNode
    do i = 1, nb_InflNode
        ! out of the computational grid
        if (InflNode(i) .gt. nb_gridnode .or. &
            InflNode(i) .le. 0) then
        cycle
        end if
        r = abs(rpg(i,:))
        do j = 1, 2
            if (r(j) <= 0.25d0) then
                sh(j) = (7d0 - 16d0 * r(j) * r(j)) * 0.125d0
                dn(j) = -4d0 * rpg(i,j)
            else if (r(j) <= 0.75d0) then
                sh(j) = 1d0 - r(j)
                dn(j) = -sign(real(1.0,8),rpg(i,j))
            else if (r(j) <= 1.25d0) then
                sh(j) = ((5d0 - 4d0 * r(j)) * 0.25d0) ** 2
                dn(j) = 2d0 * rpg(i,j) - sign(real(1.0,8),rpg(i,j)) * 2.5d0
            else
                sh(j) = 0.0d0
                dn(j) = 0.0d0
            end if
        end do
        pt%SHP(i) = sh(1)*sh(2)
        pt%DNDX(i) = dn(1)*sh(2)*iDCell
        pt%DNDY(i) = sh(1)*dn(2)*iDCell
    end do

    end subroutine NShape_GIMP

    subroutine FindInflNode(p, icell)
    !------------------------------------------------------------------
    !-  Purpose: Find Influence Nodes for particle p                  -
    !-  Inputs                                                        -
    !-      p     - particle number                                   -
    !-      icell - cell id containing particle p                     -
    !-  Outputs                                                       -
    !-      InflNode(27)                                              -
    !-      rpg(27,3)                                                 -
    !------------------------------------------------------------------
    use ParticleData, only:Particle, particle_list
    implicit none
    integer, intent(in):: p, icell
    integer:: nb_InflNode, InflNode(9), i
    real(8):: rpg(9,2), x(2) ! nature coordinate ( -1 < x,y < 1 )
    type(Particle), POINTER :: pt

    rpg(:,:) = 0d0
    pt => particle_list(p)
    do i = 1, 4
        InflNode(i) = CellsNode(icell,i)
        rpg(i,:) = (pt%Xp - node_list(InflNode(i))%Xg) * iDCell
    end do

    ! lower y
    if (rpg(1,2) < 0.25d0) then
        InflNode(5) = InflNode(1) - NGx
        InflNode(6) = InflNode(2) - NGx

        ! lower x (19)
        if (rpg(1,1) < 0.25d0) then
            InflNode(7) = InflNode(1) - 1
            InflNode(8) = InflNode(4) - 1
            InflNode(9) = InflNode(5) - 1
            nb_InflNode = 9

            ! upper x (20)
        else if (rpg(1,1) > 0.75d0) then
            InflNode(7) = InflNode(2) + 1
            InflNode(8) = InflNode(3) + 1
            InflNode(9) = InflNode(6) + 1
            nb_InflNode = 9

            ! middle x (21)
        else
            nb_InflNode = 6
        end if

        ! upper y
    else if (rpg(1,2) > 0.75d0) then
        InflNode(5) = InflNode(3) + NGx
        InflNode(6) = InflNode(4) + NGx

        ! lower x (22)
        if (rpg(1,1) < 0.25d0) then
            InflNode(7) = InflNode(1) - 1
            InflNode(8) = InflNode(4) - 1
            InflNode(9) = InflNode(6) - 1
            nb_InflNode = 9

            ! upper x (23)
        else if (rpg(1,1) > 0.75d0) then
            InflNode(7) = InflNode(2) + 1
            InflNode(8) = InflNode(3) + 1
            InflNode(9) = InflNode(5) + 1
            nb_InflNode = 9

            ! middle x (24)
        else
            nb_InflNode = 6
        end if

        ! middle y
    else

        ! lower x (25)
        if (rpg(1,1) < 0.25d0) then
            InflNode(5) = InflNode(1) - 1
            InflNode(6) = InflNode(4) - 1
            nb_InflNode = 6

            ! upper x (26)
        else if (rpg(1,1) > 0.75d0) then
            InflNode(5) = InflNode(2) + 1
            InflNode(6) = InflNode(3) + 1
            nb_InflNode = 6

            ! middle x (27)
        else
            nb_InflNode = 4
        end if
    end if

    do i = 5, nb_InflNode
        if (InflNode(i) .le. nb_gridnode .and. InflNode(i) .gt. 0) then
            rpg(i,:)= pt%Xp - node_list(InflNode(i))%Xg
        else    ! out of the computational grid
            rpg(i,:) = 100.0
        end if
        rpg(i,:) = rpg(i,:) * iDCell
    end do

    pt%nb_InflNode = nb_InflNode
    pt%InflNode = InflNode
    pt%rpg = rpg
    
    end subroutine FindInflNode

    subroutine FindSurfaceCN()
    !-------------------------------------------------------------------
    !-  Purpose                                                        -
    !-      1. Find surface nodes % cells                              -
    !-------------------------------------------------------------------

    implicit none
    integer:: n, inode, ix, iy ! loop counter
    integer:: icell, comID = 1
    logical:: surfaceCellNotFound
    type(GridNodeProperty), POINTER :: gd
    type(GridCellProperty), POINTER :: gc

    ! find surface cells from top to bottom in y direction
    do ix = 1, NumCellx
        surfaceCellNotFound = .true.
        do iy = NumCelly, 1, -1
            icell = (iy - 1)*NumCellx + ix
            gc => cell_list(icell)
            
            if (gc%InOutProp == 1) then
                do n = 1, nCellNode
                    inode = CellsNode(icell, n)
                    grid_list(inode, comID)%InOutProp = 1
                end do
            end if
            
            if (gc%InOutProp > 0 .and. surfaceCellNotFound) then
                surfaceCellNotFound = .false.
                gc%InOutProp = 2
                inode = CellsNode(icell, 3)
                grid_list(inode, comID)%InOutProp = 2
                inode = CellsNode(icell, 4)
                grid_list(inode, comID)%InOutProp = 2
            end if

        end do
    end do
    
    ! find surface cells from top to bottom in x direction
    do iy = 1, NumCelly
        surfaceCellNotFound = .true.
        do ix = NumCellx, 1, -1
            icell = (iy - 1)*NumCellx + ix
            gc => cell_list(icell)
            
            if (gc%InOutProp > 0 .and. surfaceCellNotFound) then
                surfaceCellNotFound = .false.
                gc%InOutProp = 2
                inode = CellsNode(icell, 2)
                grid_list(inode, comID)%InOutProp = 2
                inode = CellsNode(icell, 3)
                grid_list(inode, comID)%InOutProp = 2
            end if

        end do
    end do

    end subroutine FindSurfaceCN
    
    end module GridData
#endif


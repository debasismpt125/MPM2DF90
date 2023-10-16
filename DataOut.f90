
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

    !-------------------------------------------------------------------
    !-                                                                 -
    !-   Result output procedures                                      -
    !-                                                                 -
    !-------------------------------------------------------------------

    module DataOut
    use GridData
    use ParticleData
   

    integer:: nCvs = 0
    integer:: nAnm = 0

    integer, parameter:: nVariables = 24
    character(4), parameter:: OutputName(nVariables) = (/ &
        'seqv','epef','mat ','pres','volu',  &
        'engk','engi','velx','vely','velz',  &
        'cels','fail','sspd','damg','sxx ',  &
        'syy ','szz ','sxy ', 'void', 'pc  ',  &
		'pore', 'pm  ', 'qd  ', 'ptot' &
		
        /)

    logical:: WriteTecPlot  = .false.
    logical:: WriteParaView = .false.

    real(8):: OutTime = 0.0     ! time interval for plot data output
    real(8):: ReportTime = 0.0  ! time interval for status report

    ! Number of curves defined, default 2 curves
    integer:: nCurves  = 2
    ! Number of animation variables defined
    integer:: nAnimate = 1

    ! Limited 15 curves
    integer, parameter :: MaxCurves = 21
    ! Limited 10 animation variables
    integer, parameter :: MaxAnim = 16

    ! variable of each curve
    integer  CurveOption(MaxCurves) /6,7,19*0/
    ! animation variables
    integer  AnimOption(MaxAnim) /3,15*0/
    integer:: CurvePoint(MaxCurves) = 1

    real:: plot2d(6) = 0.0      ! plot less particles
    logical:: plot2dTrue = .false.

    contains

    subroutine OutCurve()
    !------------------------------------------------------------------
    !-  purpose: output result data for plotting time-history curve   -
    !------------------------------------------------------------------
    use FFI, only: iow2,iow03,iow04,iow05
    implicit none
    integer i

    character(20) istr

    type(Particle), pointer:: pt

    nCvs = nCvs + 1
    ! write variable name
    if (nCvs .eq. 1) then
        write(iow2,"(a5)",advance='no') 'Time '
        print *, 'ncurves=', nCurves
        do i = 1, nCurves
            write(istr,*) CurvePoint(i)
            istr = trim(adjustl(istr))
            write(iow2,"(a4,a8)",advance='no') &
                OutputName(CurveOption(i)), istr
        end do
        write(iow2,*)
    end if

    write(iow2,"(e12.4)", advance='no') CurrentTime

    do i = 1, nCurves

        pt => particle_list(CurvePoint(i))

        select case(CurveOption(i))

        case(1) !seqv von Mises equivalent stress
            write(iow2,"(e12.4)", advance='no') pt%seqv
        case(2) !epeff effective accumulative plastic strain
            write(iow2,"(e12.4)", advance='no') pt%epeff
        case(3) !no need to curve mat, so output epeff instead
            write(iow2,"(e12.4)", advance='no') pt%epeff
        case(4) !pressure
            write(iow2,"(e12.4)", advance='no') -pt%SM
        case(5) !volume
            write(iow2,"(e12.4)", advance='no') pt%vol
        case(6) !engk
            write(iow2,"(e12.4)", advance='no') EngKinetic
        case(7) !engi
            write(iow2,"(e12.4)", advance='no') EngInternal
        case(8) !vx
            write(iow2,"(e12.4)", advance='no') pt%VXp(1)
        case(9) !vy
            write(iow2,"(e12.4)", advance='no') pt%VXp(2)
#ifdef NDIM3
            write(iow2,"(e12.4)", advance='no') pt%VXp(3)
#endif
        case(11) !celsius_t
            write(iow2,"(e12.4)", advance='no') pt%celsius_t
        case(12) !fail
            write(iow2,"(i12)", advance='no') pt%failure
        case(13) !szz
            write(iow2,"(e12.4)", advance='no') pt%SDzz+pt%SM
        case(14) !damg
            write(iow2,"(e12.4)", advance='no') pt%DMG
		case(15) !sxx
            write(iow2,"(e12.4)", advance='no') pt%SDxx+pt%SM
        case(16) !syy
            write(iow2,"(e12.4)", advance='no') pt%SDyy+pt%SM
        case(17) !szz
            write(iow2,"(e12.4)", advance='no') pt%SDzz+pt%SM
        case(18) !sxy
        case(19) !void ratio
            write(iow2,"(e12.4)", advance='no') pt%evoid
        case(20) !pc
            write(iow2,"(e12.4)", advance='no') pt%pc
        case(21) !pore
            write(iow2,"(e12.4)", advance='no') pt%pore
		case(22) !pm
            write(iow2,"(e12.4)", advance='no') pt%pm
	    case(23) !qd
            write(iow2,"(e12.4)", advance='no') pt%qd
		case(24) !ptot
            write(iow2,"(e12.4)", advance='no') pt%ptot
        end select

    end do
    write(iow2,*)

    ! add code for plot energy and momentum

    write(iow03,"(4e12.4)") CurrentTime, EngKinetic+EngInternal, &
        EngKinetic, EngInternal
    
    
#ifdef NDIM3
    write(iow04,"(4e12.4)") CurrentTime, Mombody1(1), Mombody1(2), &
        Mombody1(3)
    if (contact) then
        ! output contact force
        write(iow05,"(4e12.4)") CurrentTime, tot_cont_for(1), &
            tot_cont_for(2), tot_cont_for(3)
    end if

#elif NDIM2
    write(iow04,"(4e12.4)") CurrentTime, Mombody1(1), Mombody1(2)
    if (contact) then
        ! output contact force
        write(iow05,"(4e12.4)") CurrentTime, tot_cont_for(1), &
            tot_cont_for(2)
    end if
#endif
    
    end subroutine OutCurve


    subroutine OutAnim()
    !------------------------------------------------------------------
    !-  purpose: output result data for animation                     -
    !------------------------------------------------------------------
    use FFI, only: iow1, iomsg
    implicit none

    type(Particle), POINTER :: pt
    integer p, i, b, parBegin, parEnd ! loop counter
    integer:: fCounter ! failure particle counter
    integer matID      ! material number

    fCounter = 0
    nAnm = nAnm + 1
    if (nAnm .eq. 1) then
        write(iow1,10) Title
        write(iow1,20) ('"', OutputName(AnimOption(i)), '"', &
            i=1,nAnimate)
    end if
10  format('TITLE = "', a60, '"')
#ifdef NDIM2
20 format('VARIABLES= "X"   "Y"   ', 50(A2, A4, A1))
#elif NDIM3
20  format('VARIABLES= "X"   "Y"   "Z"  ', 50(A2, A4, A1))
#endif
    write(iow1,"(a12,e13.5,a1)") 'ZONE T="time', CurrentTime,'"'

    do b = 1, nb_body
        parBegin = body_list(b)%par_begin
        parEnd = body_list(b)%par_end
        matID = body_list(b)%mat
        do p = parBegin, parEnd

            pt => particle_list(p)
            ! omit particles out of the computational region
           !if (pt%icell < 0) cycle
            if (pt%SkipThis) cycle

#ifdef NDIM2
            write(iow1, "(2e12.4)", advance='no') pt%Xp
#elif NDIM3
            write(iow1, "(3e12.4)", advance='no') pt%Xp            
#endif
            
            do i = 1, nAnimate
                select case(AnimOption(i))

                case(1) !seqv von Mises equivalent stress
                    write(iow1,"(e12.4)", advance='no') pt%seqv
                case(2) !epeff effective accumulative plastic strain
                    write(iow1,"(e12.4)", advance='no') pt%epeff
                case(3) !mat material id
                    write(iow1,"(i6)", advance='no') matID
                case(4) !pressure
                    write(iow1,"(e12.4)", advance='no') -pt%SM
                case(5) !vol
                    write(iow1,"(e12.4)", advance='no') pt%vol
				case(6) !engk
					write(iow1,"(e12.4)", advance='no') &
						0.5*pt%mass*(pt%VXp(1)**2+pt%VXp(2)**2)
                        !0.5*pt%mass*(pt%VXp(1)**2+pt%VXp(2)**2+pt%VXp(3)**2)
				case(7) !engi
					write(iow1,"(e12.4)", advance='no') pt%ie
                case(8) !vx
                    write(iow1,"(e12.4)", advance='no') pt%VXp(1)
                case(9) !vy
                    write(iow1,"(e12.4)", advance='no') pt%VXp(2)
                case(10) !vz
#ifdef NDIM3
                    write(iow1,"(e12.4)", advance='no') pt%VXp(3)
#endif
                case(11) !celsius_t
                    write(iow1,"(e12.4)", advance='no') pt%celsius_t
                case(12) !failure
                    if(pt%failure) then
                        write(iow1,"(i3)", advance='no') 1
                    else
                        write(iow1,"(i3)", advance='no') 0
                    end if
                case(13) !sound speed
                    write(iow1,"(e12.4)", advance='no') pt%cp
                case(14) !damage
                    write(iow1,"(e12.4)", advance='no') pt%DMG
                case(15) !sxx
                    write(iow1,"(e12.4)", advance='no') pt%SDxx+pt%SM
                case(16) !syy
                    write(iow1,"(e12.4)", advance='no') pt%SDyy+pt%SM
                case(17) !szz
                    write(iow1,"(e12.4)", advance='no') pt%SDzz+pt%SM
                case(18) !sxy
				case(19) !void ratio
				    write(iow1,"(e12.4)", advance='no') pt%evoid
			    case(20) !pc
				    write(iow1,"(e12.4)", advance='no') pt%pc
				case(21) !pore
				    write(iow1,"(e12.4)", advance='no') pt%pore
				case(22) !pm
				    write(iow1,"(e12.4)", advance='no') pt%pm
			    case(23) !qd
				    write(iow1,"(e12.4)", advance='no') pt%qd
				case(24) !ptot
				    write(iow1,"(e12.4)", advance='no') pt%ptot	
                end select
            end do ! i
            write(iow1,*)
            if(pt%failure) fCounter = fCounter + 1

        end do ! p
    end do    ! b
    if(fCounter.ne.0) then
        write(iomsg,"(a,i10)") 'number of failed particles = ', fCounter
        write(*,"(a,i10)") 'number of failed particles = ', fCounter
    end if

    end subroutine OutAnim


    subroutine OutAnimPV(iStepNo)
    !------------------------------------------------------------------
    !-  purpose: output result data for animation                     -
    !------------------------------------------------------------------
    use FFI, only: iow11, iow12, iomsg, FileAnim
    implicit none

    type(Particle), POINTER :: pt
    integer, parameter:: nb_elemnode = 1
    integer, parameter:: iElemType = 1
    integer:: p, i, b, parBegin, parEnd ! loop counter
    integer:: fCounter ! failure particle counter
    integer:: matID      ! material number
    integer:: rename, iStepNo, status, idx
    character*80:: FileAnimNext, stepc
    character*60, parameter:: gridtype = 'UnstructuredGrid'
    character*4:: indent = '    '

    fCounter = 0

    ! Identify particles
    parBegin = body_list(1)%par_begin
    parEnd = body_list(1)%par_end
    matID = body_list(1)%mat

    ! Prepare file names
    write(stepc,"(I0)") iStepNo

    idx = index(FileAnim,".")
    if (idx > 0) then
        FileAnim = FileAnim(1:index(FileAnim,".")-1)
    endif

    FileAnimNext = trim(FileAnim) // "_" // trim(stepc) // ".vtu"
    FileAnim = trim(FileAnim) // ".vtu"

    ! Open file
    open(iow11, file = trim(FileAnim), status = 'unknown')

    ! XML HEADER tags :: OPEN
    write(iow11,300)
    write(iow11,301) 'VTKFile', TRIM(gridtype)
    write(iow11,31)  '', TRIM(gridtype)
    write(iow11,"(A,I0,A,I0,A)") indent//'<Piece NumberOfPoints="', nb_particle, &
        '" NumberOfCells="', nb_particle,'">'

    ! 'POINTS' XML TAG :: OPEN
    write(iow11,31) indent//indent, 'Points'

    ! 'POINTS' Container
    write(iow11,*)  indent//indent//indent, &
        '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

    do p = 1, nb_particle
        pt => particle_list(p)
#ifdef NDIM3
        write(iow11,35) indent//indent//indent//indent, &
            dble(pt%Xp(1)), dble(pt%Xp(2)), dble(pt%Xp(3))
#elif NDIM2
        write(iow11,35) indent//indent//indent//indent, &
            dble(pt%Xp(1)), dble(pt%Xp(2))
#endif
    end do

    write(iow11,32) indent//indent//indent, 'DataArray'
    write(iow11,32) indent//indent, 'Points'

    ! 'POINTDATA' XML TAG :: OPEN
    write(iow11,31) indent//indent, 'PointData Scalars="scalars"'

    ! seqv
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'seqv'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, pt%seqv
    end do !p
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'seqv'

    ! epeff
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'epeff'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, pt%epeff
    end do !p
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'epeff'

    ! pressure
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'pressure'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, -pt%SM
    end do !p
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'pressure'

    ! vx
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'vx'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, -pt%VXp(1)
    end do !p
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'vx'

    ! vy
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'vy'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, -pt%VXp(2)
    end do !p
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'vy'

#ifdef NDIM3
    ! vz
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'vz'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, -pt%VXp(3)
    end do !p
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'vz'
#endif
    ! evoid
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'evoid'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, pt%evoid
    end do !p
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'evoid'
	! pc
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'pc'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, pt%pc
    end do !pc
	write(iow11,32) indent//indent//indent, 'DataArray' ! close 'pc'
	! pore
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'pore'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, pt%pore
    end do !pore
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'pore'
		! pm
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'pm'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, pt%pm
    end do !pm
	write(iow11,32) indent//indent//indent, 'DataArray' ! close 'pm'
    ! qd
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'qd'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, pt%qd
    end do !qd
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'qd'
	write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'ptot'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, pt%ptot
    end do !qd
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'ptot'
    ! vol
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'vol'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, pt%vol
    end do !p
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'vol'

    ! celsius_t
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'celsius_t'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, pt%celsius_t
    end do !p
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'celsius_t'

    ! damage
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'damage'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,e12.4)") indent//indent//indent//indent, pt%DMG
    end do !p
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'damage'

    ! 'POINTDATA' XML TAG :: CLOSE
    write(iow11,32) indent//indent, 'PointData'

    ! 'CELLS' XML TAG :: OPEN
    write(iow11,31) indent//indent, 'Cells'

    ! 'CONNECTIVITY'
    write(iow11,33) indent//indent//indent, &
        'DataArray', 'Int32', 'connectivity'
    do p = 1, nb_particle
        pt => particle_list(p)
        write(iow11, "(A,i5)") indent//indent//indent//indent, p-1
    end do
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'connectivity'

    ! 'OFFSETS'
    write(iow11,33) indent//indent//indent, 'DataArray', 'Int32', 'offsets'
    do p = 1, nb_particle
        write(iow11, "(A,i5)") indent//indent//indent//indent, p
    end do
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'offsets'

    ! 'TYPES'
    write(iow11,33) indent//indent//indent, 'DataArray', 'Float32', 'types'
    do p = 1, nb_particle
        write(iow11, "(A,i5)") indent//indent//indent//indent, iElemType
    end do
    write(iow11,32) indent//indent//indent, 'DataArray' ! close 'types'

    ! 'CELLS' XML TAG :: CLOSE
    write(iow11,32) indent//indent, 'Cells'

    ! XML HEADER tags :: CLOSE
    write(iow11,32) indent, 'Piece'
    write(iow11,32) '',TRIM(gridtype)
    write(iow11,32) '','VTKFile'
    close(iow11)

    ! Close file

    ! Rename for next animation step
    status = rename(trim(FileAnim), trim(FileAnimNext))
    if ( status .ne. 0 ) stop 'rename: error'

    ! Add file to animation collection
    write(iow12,"(A,e12.4,A)") indent//'<DataSet timestep="', CurrentTime, &
        '" group="" part="0" file="'//trim(FileAnimNext)//'"/>'

300 format('<?xml version="1.0"?>')
301 format('<', A, ' type="', A,'" version="0.1" byte_order="BigEndian">')

31  format(A,'<', A, '>')				! open
32  format(A,'</', A, '>')				! close
33	format(A,'<', A, ' type="', A,'" Name="', A,'" Format="ascii">') !DataArray


#ifdef NDIM3
35  format(A, e12.4, '  ', e12.4, '  ', e12.4) ! coordinate
#elif NDIM2
35  format(A, e12.4, '  ', e12.4) ! coordinate
#endif
39  format(a,i0,a,i0,a)

    end subroutine OutAnimPV


    subroutine OutGrid()
    ! -----------------------------------------------------------------
    ! - Purpose                                                       -
    ! -    Write computational grid to output file viewabla by PVW    -
    ! -----------------------------------------------------------------
    use FFI, only: iow10, fName
    implicit none

    integer:: iCell, iNode
#ifdef NDIM3
    integer:: iNodeOrder(8) = (/ 1, 2, 4, 3, 5, 6, 8, 7 /)
    integer, parameter:: nb_elemnode = 8
#elif NDIM2
    integer:: iNodeOrder(4) = (/ 1, 2, 3, 4 /)
    integer, parameter:: nb_elemnode = 4
#endif
    integer, parameter:: iElemType = 11
    character*60, parameter:: gridtype = 'UnstructuredGrid'
    character*4:: indent = '    '
    type(GridNode), POINTER:: node

    open(iow10, file = trim(fName) // "_Grid.vtu", status = 'unknown')

    ! XML HEADER tags :: OPEN
    write(iow10,300)
    write(iow10,301) 'VTKFile', TRIM(gridtype)
    write(iow10,31)  '',TRIM(gridtype)
    write(iow10,"(A,I0,A,I0,A)") indent// &
        '<Piece NumberOfPoints="', nb_gridnode,'" NumberOfCells="', NumCell,'">'

    ! 'POINTS' XML TAG :: OPEN
    write(iow10,31) indent//indent, 'Points'

    ! 'POINTS' Container
    write(iow10,*)  indent//indent//indent, &
        '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
    do iNode = 1, nb_gridnode
        node => node_list(iNode)
#ifdef NDIM3
        write(iow10,35) indent//indent//indent//indent, &
            dble(node%Xg(1)), dble(node%Xg(2)), dble(node%Xg(3))
#elif NDIM2
        write(iow10,35) indent//indent//indent//indent, &
            dble(node%Xg(1)), dble(node%Xg(2))
#endif
    end do
    write(iow10,32)  indent//indent//indent, 'DataArray'
    write(iow10,32) indent//indent, 'Points'

    ! 'CELLS' XML TAG :: OPEN
    write(iow10,31) indent//indent, 'Cells'

    ! 'CONNECTIVITY'
    write(iow10,33) indent//indent//indent, 'DataArray', 'Int32', 'connectivity'
    do iCell = 1, NumCell
        write(iow10, "(a)", advance='no') indent//indent//indent//indent
        do iNode = 1, nb_elemnode
            write(iow10,"(i0,A)", advance='no') CellsNode(iCell, iNodeOrder(iNode))-1, indent
        end do
        write(iow10, *)
    end do
    write(iow10,32) indent//indent//indent, 'DataArray' ! close 'connectivity'

    ! 'OFFSETS'
    write(iow10,33) indent//indent//indent, 'DataArray', 'Int32', 'offsets'
    do iCell = 1, NumCell
        write(iow10, "(a,i0)") indent//indent//indent//indent, nb_elemnode * iCell
    end do
    write(iow10,32) indent//indent//indent, 'DataArray' ! close 'offsets'
    ! 'TYPES'
    write(iow10,33) indent//indent//indent, 'DataArray', 'Int32', 'types'
    do iCell = 1, NumCell
        write(iow10, "(a,i5)") indent//indent//indent//indent, iElemType
    end do
    write(iow10,32) indent//indent//indent, 'DataArray' ! close 'types'
    ! 'CELLS' XML TAG :: CLOSE
    write(iow10,32) indent//indent, 'Cells'

    ! XML HEADER tags :: CLOSE
    write(iow10,32) indent,'Piece'
    write(iow10,32) '', TRIM(gridtype)
    write(iow10,32) '', 'VTKFile'

300 format('<?xml version="1.0"?>')
301 format('<', A, ' type="', A,'" version="0.1" byte_order="BigEndian">')

31  format(A,'<', A, '>')				! open
32  format(A,'</', A, '>')				! close
33  format(A,'<', A, ' type="', A,'" Name="', A,'" Format="ascii">') !DataArray
#ifdef NDIM3
35  format(A, e12.4, '  ', e12.4, '  ', e12.4) ! coordinate
#elif NDIM2
35  format(A, e12.4, '  ', e12.4) ! coordinate 
#endif
    
39  format(a,i0,a,i0,a)

    close(iow10)

    end subroutine OutGrid

#ifdef NDIM3
    subroutine OutPrint(timeFlag)
    use GridData
    use ParticleData
    implicit none
    integer:: i, comID = 1
    integer:: iow96 = 196, iow98 = 198, iow99 = 199
    logical:: timeFlag
    character(10) timestr
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd
    
    if (timeFlag) then
        write(timestr, '(F10.3)') CurrentTime
        open(iow98, file = "node_"// trim(adjustl(timestr)) // ".dat", status = 'replace')
        open(iow99, file = "particle_"// trim(adjustl(timestr)) // ".dat", status = 'replace')
    else
        open(iow96, file = "time.dat", status = 'replace')
        write(iow96,'(F10.3)') CurrentTime
        close(iow96)
        
        open(iow98, file = "node.dat", status = 'replace')
        open(iow99, file = "particle.dat", status = 'replace')
        
    end if
    
    write(iow98,"(A)") "   Xg(1)      Xg(2)       Xg(3)       Mg          PXg(1)      PXg(2)      PXg(3)"
    do i = 1, nb_gridnode
        write(iow98,"(1X,7(ES10.3,2X))") node_list(i)%Xg, grid_list(i, comID)%Mg,grid_list(i, comID)%PXg
    end do
    close(iow98)

    write(iow99,"(A)") "  Xp(1)       Xp(2)       Xp(3)       Mass        SM          SigX        SigY        SigZ        SigXZ"
    do i = 1, nb_particle
        pt => particle_list(i)
        write(iow99,"(1X,9(ES10.3,2X))") pt%Xp, pt%Mass, pt%SM, &
            pt%SDxx + pt%SM, pt%SDyy + pt%SM, &
            pt%SDzz + pt%SM, pt%SDxz
    end do
    close(iow99)

    end subroutine OutPrint
#endif
    
#ifdef NDIM2
    subroutine OutPrint(timeFlag)
    use GridData
    use ParticleData
    implicit none
    integer:: i, comID = 1
    integer:: iow96 = 196, iow98 = 198, iow99 = 199
    logical:: timeFlag
    character(10) timestr
    type(Particle), POINTER :: pt
    type(GridNodeProperty), POINTER :: gd
    
    if (timeFlag) then
        write(timestr, '(F10.3)') CurrentTime
        open(iow98, file = "node_"// trim(adjustl(timestr)) // ".dat", status = 'replace')
        open(iow99, file = "particle_"// trim(adjustl(timestr)) // ".dat", status = 'replace')
    else
        open(iow96, file = "time.dat", status = 'replace')
        write(iow96,'(F10.3)') CurrentTime
        close(iow96)
        
        open(iow98, file = "node.dat", status = 'replace')
        open(iow99, file = "particle.dat", status = 'replace')
        
    end if

    write(iow98,"(A)") "   Xg(1)      Xg(2)       Mg          PXg(1)      PXg(2)"
    do i = 1, nb_gridnode
        write(iow98,"(1X,5(ES10.3,2X))") node_list(i)%Xg, grid_list(i, comID)%Mg,grid_list(i, comID)%PXg
    end do
    close(iow98)

    write(iow99,"(A)") "  Xp(1)       Xp(2)       Mass        SM          SigX        &
	                      SigY        SigZ        SigXY        VX1        VX2"
    do i = 1, nb_particle
        pt => particle_list(i)
        write(iow99,"(1X,10(ES10.3,2X))") pt%Xp, pt%Mass, pt%SM, &
            pt%SDxx + pt%SM, pt%SDyy + pt%SM, &
            pt%SDzz + pt%SM, pt%SDxy, pt%VXp
    end do
    close(iow99)

    end subroutine OutPrint
#endif
    
    end module DataOut

    subroutine CheckMain
    use mesh
    use NRRMeshMod,only:NRRCheckPointer
    implicit none
    type(gridtyp),pointer::c
    integer::i

    !!Check the NRRpointer is pointing to a right cell.
    do i=1,total
        c=>cell(i)
        call NRRCheckPointer(c)
    enddo
    endsubroutine CheckMain
    !_________________________________________
    subroutine check_node(ccount)
    use mesh
    use outfile
    use typeMod
    implicit none
    real(pre)::nx,ny
    integer::i,j,k,repeat=0
    integer,intent(in)::ccount
    
    allocate(Node(ccount*2))

    call nodeinf
    do i=1,nnode
        nx=node(i)%x
        ny=node(i)%y
        do j=1+i,nnode	!!减少循环次数
            if(nx==node(j)%x.and.ny==node(j)%y.and.i/=j)then
            repeat=repeat+1
            print*,"repear node",i,j
            endif
        enddo
    enddo
    deallocate(Node)
    write(*,*) "Repeat Node count:", repeat
    endsubroutine
!_________________________________________
    subroutine stopall(from)
    implicit none
    character(*)::from
    
    WRITE(*,2) FROM
2   FORMAT(/' STOP_ALL called from routine ',A/)
    stop
    endsubroutine
!_________________________________________
    subroutine CheckLimiter
    use solution
    use inflow
    use typeMod
    implicit none

    if(nLUmax>0)then
        print*,"U_max   Limit over ",LimitUmax," for ",nLUmax," cells."
        nLUmax=0
    endif
    if(nLVmax>0)then
        print*,"V_max   Limit over ",LimitVmax," for ",nLVmax," cells."
        nLVmax=0
    endif
    if(nLRoumax>0)then
        print*,"Rou_max Limit over ",LimitRoumax," for ",nLRoumax," cells."
        nLRoumax=0
    endif
    if(nLRoumin>0)then
        print*,"Rou_min Limit over ",LimitRoumin," for ",nLRoumin," cells."
        nLRoumin=0
    endif
    if(nLTmax>0)then
        print*,"T_max   Limit over ",LimitTmax," for ",nLTmax," cells."
        nLTmax=0
    endif
    if(nLTmin>0)then
        print*,"T_min   Limit over ",LimitTmin," for ",nLTmin," cells."
        nLTmin=0
    endif
    if(nLPmax>0)then
        print*,"P_max   Limit over ",LimitPmax," for ",nLPmax," cells."
        nLPmax=0
    endif
    if(nLPmin>0)then
        print*,"P_min   Limit over ",LimitPmin," for ",nLPmin," cells."
        nLPmin=0
    endif
    endsubroutine CheckLimiter
!_________________________________________
!_________________________________________
!_________________________________________

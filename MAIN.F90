program main

  use precisionMod
  use typeMod
  use inflowGlobal
  use solution!,only:solution_initial,solute
  use outfile,only:output
  use flow_result,only:surface_data
  use NRRMeshMod,only:NRRMeshMain
  use neighbor
  use mesh
    implicit none

  integer i,j,cnum,alvx,alvy
  type(gridtyp),pointer::c
  real(pre) tx,ty,ttx,t1,t2,t3,t4  !
  integer a

  allocate(xd(dmax))
  allocate(yd(dmax))

  m=X0/h
  n=Y0/h
  total=m*n
  CellCount=total
  Rcenter(1)=x0/4
  Rcenter(2)=y0/2
  tstep=h/(2**lvlmax)*Courant
  allocate(cell(total))
 open(20,file='./Others/model.dat')
 open(21,file='./Others/modelnew.dat')
  do i=1,dmax
    read(20,*)xd(i),yd(i)

    xd(i)=xd(i)+Rcenter(1)
    yd(i)=yd(i)+Rcenter(2)
    write(21,*)xd(i),yd(i)
  end do
  close(20)
  close(21)

  call initialMesh

    if(NRR==1)then
        call NRRMeshMain
    else
        do i=1,total
          c=>cell(i)
          call gridrefm(c)
        end do
        print*,"initial mesh done"
        do j=1,lvltemp
          do i=1,total
            c=>cell(i)
            if(associated(c%son1)) call gridmodify(c)
          end do
        end do
        print*,"refine isotropic-mesh done"

        if(Anisotropic==1)then
            do i=1,total
              c=>cell(i)
            call gridgxyx(c)
            end do
            print*,"insitial anisotropic-mesh done"

            do j=1,lvlgxyx
                do i=1,total
                c=>cell(i)
                call init_s(c)
                end do
                do i=1,total
                c=>cell(i)
                call gxyxmodify_def(c)
                end do
                do i=1,total
                c=>cell(i)
                call gxyxmodify_ref(c)
                end do
            end do
            print*,"refine anisotropic-mesh done"
        endif
    endif
    
    if(check==1) call CheckMain
    
    do i=1,total
        c=>cell(i) 
        call solution_initial(c)
    end do
    
    call cellinf
    call output('Mesh',CellCount,.false.)
    if(MeshOnly==1)then
        goto 999
    endif
    
    call solute
    
    call output('ok',CellCount,.true.)
    call surface_data('ok')
999 stop
end program
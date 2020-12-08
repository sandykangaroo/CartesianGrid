module NRRMeshMod
    use neighbor
    use mesh
    use inflow
    use precisionMod
    use typeMod
    use inflowNRR
    implicit none
    type NRRray
        real(predd),dimension(4)::bbox,bboy  !!x1-4,y1-4
    endtype
    type(NRRray), dimension(NRRSeeds)::nrrrays
    real(predd),dimension(NRRSeeds,2)::SeedPoints
    integer::nrays  !!Corrent NRRSeeds count
    type(pointdd)::NRRnormal(NRRSeeds)
    
    private nrays
    private NRRBuildBox  !!private subroutine
    contains

!_________________________________________
    subroutine NRRMeshMain()
    type(gridtyp),pointer::c
    integer::i,j

    print*,"_______Normal Ray Refinement________"
    print*,"_______Xueliang Li  2020.11.9_______"
    print*,"Start NRR mesh generation"
    !!Initial cell%s1=0
    !!Locate NRR bounding box and generate mesh.
    if(NRRline==1)then
        call NRRBuildline
        do i=1,NRRSeeds
            nrays=i
            do j=1,total
                c=>cell(j)
                call NRRRefineLine(c)
            enddo
        enddo
    else
        do i=1,total
            c=>cell(i)
            call init_s(c)
        enddo
        call NRRBuildBox
        do i=1,NRRSeeds
            nrays=i
            do j=1,total
                c=>cell(j)
                call NRRRefineGlobal(c)
            enddo
        enddo
    endif
    !!!Special aera refinement.
    !do i=1,total
    !    c=>cell(i)
    !    call NRRspecialRefinement(c)
    !enddo
    !!Refine the cross-wall mesh.
    do i=1,total
        c=>cell(i)
        call NRRcrossMeshRefine(c)
    enddo
    !!Mesh optimization
    do j=1,NRRlvl
        do i=1,total
            c=>cell(i)
            if(associated(c%son1))  call gridmodify(c)
        enddo
    enddo
    !!NRR region mark for the cell in red, blue, and black.
    do i=1,total
        c=>cell(i)
        call NRRregionMark(c)
    enddo
    print*,"Done"
    endsubroutine NRRMeshMain
!_________________________________________
!_________________________________________
    recursive subroutine NRRRefineGlobal(c)
    type(gridtyp),pointer::c
    integer::i
    
    if(associated(c%son1))then
        if(.not.associated(c%son3))then
            call NRRRefineGlobal(c%son1)
            call NRRRefineGlobal(c%son2)
        else
            call NRRRefineGlobal(c%son1)
            call NRRRefineGlobal(c%son2)
            call NRRRefineGlobal(c%son3)
            call NRRRefineGlobal(c%son4)
        endif
        return
    else
        call NRRrefineSon(c)
    endif
    endsubroutine NRRRefineGlobal
!_________________________________________
    subroutine NRRBuildline
    integer::i
    real(predd)::ddeg
    !!       delta deg
    real(predd)::Length,deviant
    
    deviant=(h/2**(NRRlvl))*0.001
    ddeg=2*pi/(NRRSeeds)
    Length=NRRLength*(h/2**(NRRlvl))
    do i=1,NRRSeeds
        NRRnormal(i)%x=cos(ddeg*(i-1))
        NRRnormal(i)%y=sin(ddeg*(i-1))
        SeedPoints(i,1)=Rcenter(1)+NRRnormal(i)%x+deviant
        SeedPoints(i,2)=Rcenter(2)+NRRnormal(i)%y+deviant
        !!  Make a BBox
        nrrrays(i)%bbox(1)=SeedPoints(i,1)
        nrrrays(i)%bbox(2)=SeedPoints(i,1)+NRRnormal(i)%x*Length
        nrrrays(i)%bboy(1)=SeedPoints(i,2)
        nrrrays(i)%bboy(2)=SeedPoints(i,2)+NRRnormal(i)%y*Length
    enddo
    endsubroutine NRRBuildline
!_________________________________________
    recursive subroutine NRRRefineLine(c)
    type(gridtyp),pointer::c
    integer::i
    logical::cross

    
    if(associated(c%son1))then
        if(.not.associated(c%son3))then
            call NRRRefineLine(c%son1)
            call NRRRefineLine(c%son2)
        else
            call NRRRefineLine(c%son1)
            call NRRRefineLine(c%son2)
            call NRRRefineLine(c%son3)
            call NRRRefineLine(c%son4)
        endif
        return
    else
        cross=NRRcrossLine(c,nrrrays(nrays)%bbox,nrrrays(nrays)%bboy)
        if(c%lvl==NRRlvl)then
            c%NRRregion=1
            return
        elseif(cross==0)then
            return
        else
            c%typ=0
            call newson(c)
            call gridmodify(c)
            call NRRRefineLine(c%son1)
            call NRRRefineLine(c%son2)
            call NRRRefineLine(c%son3)
            call NRRRefineLine(c%son4)
        endif
    endif
    endsubroutine NRRRefineLine
!_________________________________________
    logical function NRRcrossLine(c,bbox,bboy)!!=1,cross; =0,detached
    real(predd),  DIMENSION(4), INTENT (IN) :: bbox,bboy
    type(gridtyp),pointer::c
    real(predd) k,B,xx1,xx2,yy1,yy2,maxx,minx,maxy,miny,x,y
    integer i,j,mm,j1,j2
    integer::counter
    
    !cell box: x1,x2,y1,y2
    xx1=c%center%x-h/2**(c%lvlx+1)
    xx2=c%center%x+h/2**(c%lvlx+1)
    yy1=c%center%y-h/(2**(c%lvly+1))
    yy2=c%center%y+h/(2**(c%lvly+1))
    !quick check
    if(bbox(1)>=bbox(2))then
        maxx=bbox(1)
        minx=bbox(2)
    else
        maxx=bbox(2)
        minx=bbox(1)
    endif
    if(bboy(1)>=bboy(2))then
        maxy=bboy(1)
        miny=bboy(2)
    else
        maxy=bboy(2)
        miny=bboy(1)
    endif
    if(xx2<minx.or.xx1>maxx.or.yy2<miny.or.yy1>maxy)then
        NRRcrossLine=0
        return
    endif
    NRRcrossLine=0
    !NRR line: y=kx+C
    !if(bbox(1)==bbox(2))then
    !    if(yy2<miny.or.yy1>maxy) NRRcrossLine=0
    !    return
    !endif
    !if(bboy(1)==bboy(2))then
    !    if(xx2<minx.or.xx1>maxx) NRRcrossLine=0
    !    return
    !endif
    k=(bboy(1)-bboy(2))/(bbox(1)-bbox(2))
    B=bboy(1)-k*bbox(1)
    !test in axis-x
    y=k*xx1+B
    if(yy1<y.and.y<yy2) NRRcrossLine=1
    y=k*xx2+B
    if(yy1<y.and.y<yy2) NRRcrossLine=1
    x=(yy1-B)/k
    if(xx1<x.and.x<xx2) NRRcrossLine=1
    x=(yy2-B)/k
    if(xx1<x.and.x<xx2) NRRcrossLine=1
    return
    endfunction NRRcrossLine
!_________________________________________
    subroutine NRRBuildBox() !!Locate NRR bounding box
    integer::i
    real(predd)::ddeg,sinddeg,cosddeg
    !!       delta deg
    real(predd)::Width,Length
    
    ddeg=2*pi/(NRRSeeds)
    Width=NRRWidth*(h/2**(NRRlvl))
    Length=NRRLength*(h/2**(NRRlvl))
    do i=1,NRRSeeds
        sinddeg=sin(ddeg*(i-1))
        cosddeg=cos(ddeg*(i-1))
        SeedPoints(i,1)=10+cosddeg
        SeedPoints(i,2)=10+sinddeg
    !!  Make a BBox
    nrrrays(i)%bbox(1)=SeedPoints(i,1)-sinddeg*Width/2
    nrrrays(i)%bbox(2)=SeedPoints(i,1)+sinddeg*Width/2
    nrrrays(i)%bbox(3)=SeedPoints(i,1)+sinddeg*Width/2+cosddeg*Length
    nrrrays(i)%bbox(4)=SeedPoints(i,1)-sinddeg*Width/2+cosddeg*Length
    nrrrays(i)%bboy(1)=SeedPoints(i,2)+cosddeg*Width/2
    nrrrays(i)%bboy(2)=SeedPoints(i,2)-cosddeg*Width/2
    nrrrays(i)%bboy(3)=SeedPoints(i,2)-cosddeg*Width/2+sinddeg*Length
    nrrrays(i)%bboy(4)=SeedPoints(i,2)+cosddeg*Width/2+sinddeg*Length
    enddo
    endsubroutine NRRBuildBox
!_________________________________________
    subroutine NRRNearestCell(x,y,c)
    !!Find the nearest cell(lvl=0) to p(x,y)
    real(predd)             :: x,y,bx1,bx2,by1,by2
    type(gridtyp),pointer  ::c
    integer::i
    real(predd)::dist,minDist
    
    minDist=9*h*h
    bx1=x-h; bx2=x+h; by1=y-h; by2=y+h
    do i=1,total
        if(cell(i)%center%x<bx1) cycle
        if(cell(i)%center%x>bx2) cycle
        if(cell(i)%center%y<by1) cycle
        if(cell(i)%center%y>by2) cycle
        dist=(cell(i)%center%x-x)**2+(cell(i)%center%y-y)**2
        if(dist<minDist)then
            minDist=dist
            c=>cell(i)
        endif
    enddo
    call NRRNearestCellSon(x,y,c)
    endsubroutine NRRNearestCell
!_________________________________________
    recursive subroutine NRRNearestCellSon(x,y,c)
    !!Find the nearest son cell to p(x,y)
    type(gridtyp),pointer  ::c,ctmp
    real(predd)::dist,minDist
    real(predd)             :: x,y
    
    if(associated(c%son1))then
        minDist=(c%son1%center%x-x)**2+(c%son1%center%y-y)**2
        ctmp=>c%son1
        dist=(c%son2%center%x-x)**2+(c%son2%center%y-y)**2
        if(dist<minDist)then
            minDist=dist
            ctmp=>c%son2
        endif
        dist=(c%son3%center%x-x)**2+(c%son3%center%y-y)**2
        if(dist<minDist)then
            minDist=dist
            ctmp=>c%son3
        endif
        dist=(c%son4%center%x-x)**2+(c%son4%center%y-y)**2
        if(dist<minDist)then
            minDist=dist
            ctmp=>c%son4
        endif
        c=>ctmp
        call NRRNearestCellSon(x,y,c)
    endif
    endsubroutine NRRNearestCellSon
!_________________________________________
    subroutine NRRRefine(x,y,c)
    type(gridtyp),pointer::c,ct1,ct2,ct3,ct4
    real(predd)             :: x,y
    integer::i
    
    do i=1,NRRlvl
        c%typ=0
        if(c%lvl>=NRRlvl) cycle
        call newson(c)
        call nullifyCell(c%son1)
        call nullifyCell(c%son2)
        call nullifyCell(c%son3)
        call nullifyCell(c%son4)
        call gridmodify(c)
        call NRRNearestCellSon(x,y,c)
    enddo
    !call NRRrefineBranch(c,0)
    !call NRRrefineUpBranch(c,0)
    !call NRRrefineDownBranch(c,0)
    !call NRRrefineLeftBranch(c,0)
    !call NRRrefineRightBranch(c,0)
    endsubroutine NRRRefine
!_________________________________________
    recursive subroutine NRRrefineSon(c)
    type(gridtyp),pointer::c
    integer::cross
    
    cross=NRRgridcross(c,nrrrays(nrays)%bbox,nrrrays(nrays)%bboy)
    if(cross==0)then
        return
    elseif(c%lvl>=NRRlvl)then
        c%NRRregion=nrays
        return
    else
        c%s1=1
        c%typ=0
        call newson(c)
        call gridmodify(c)
        call NRRrefineSon(c%son1)
        call NRRrefineSon(c%son2)
        call NRRrefineSon(c%son3)
        call NRRrefineSon(c%son4)
    endif
    endsubroutine NRRrefineSon
!_________________________________________
    recursive subroutine NRRrefineBranch(c,recurs)
    type(gridtyp),pointer::c,ct1,ct2,ct3,ct4
    integer::recurs
    
    if(associated(c%son1))then
        call NRRrefineBranch(c%son1,recurs)
        call NRRrefineBranch(c%son2,recurs)
        call NRRrefineBranch(c%son3,recurs)
        call NRRrefineBranch(c%son4,recurs)
    else
        if(c%s1==1)return
        c%s1=1
        if(c%lvl>=NRRlvl)return
        ct1=>up_neighbor(c)
        ct2=>right_neighbor(c)
        ct3=>down_neighbor(c)
        ct4=>left_neighbor(c)
        if(recurs/=0)then
            if(NRRgridcross(c,nrrrays(nrays)%bbox,nrrrays(nrays)%bboy)==0)return
        endif
        call NRRrefineSon(c)
        if(recurs/=3.and.associated(ct1)) call NRRrefineBranch(ct1,1)
        if(recurs/=4.and.associated(ct2)) call NRRrefineBranch(ct2,2)
        if(recurs/=1.and.associated(ct3)) call NRRrefineBranch(ct3,3)
        if(recurs/=2.and.associated(ct4)) call NRRrefineBranch(ct4,4)
    endif
    endsubroutine NRRrefineBranch
!_________________________________________
!_________________________________________
    integer function NRRgridcross(c,bbox,bboy)
    real(predd),  DIMENSION(4), INTENT (IN) :: bbox,bboy
    type(gridtyp),pointer::c
    integer t1,t2,t3,t4
    real(predd)::xx1,yy1,xx2,yy2,xx3,yy3,xx4,yy4

    xx1=c%center%x-h/2**(c%lvlx+1)
    yy1=c%center%y-h/(2**(c%lvly+1))
    xx2=c%center%x+h/2**(c%lvlx+1)
    yy2=c%center%y-h/(2**(c%lvly+1))
    xx3=c%center%x+h/2**(c%lvlx+1)
    yy3=c%center%y+h/(2**(c%lvly+1))
    xx4=c%center%x-h/2**(c%lvlx+1)
    yy4=c%center%y+h/(2**(c%lvly+1))
 
    t1=NRRgridsort(xx1,yy1,bbox,bboy)
    t2=NRRgridsort(xx2,yy2,bbox,bboy)
    t3=NRRgridsort(xx3,yy3,bbox,bboy)
    t4=NRRgridsort(xx4,yy4,bbox,bboy)
   
    if(t1==0.and.t2==0.and.t3==0.and.t4==0)then
    NRRgridcross=0     !!Outside bbox
    else
    NRRgridcross=1     !!Cross bbox
    end if
    end function NRRgridcross
!_________________________________________
    integer function NRRgridsort(xx,yy,bbox,bboy)
    real(predd),  DIMENSION(4), INTENT (IN) :: bbox,bboy
    integer,parameter::dnum=7
    real(predd) xx,yy,x,A,B,y
    real(predd)::d(dnum)=(/-3.0,-1.0,0.01,1.0,3.0,40.0,500.0/)
    integer i,j,mm,j1,j2
    integer::counter
    
    mm=0
    do i=1,dnum
        counter=0
        do j=1,4
            if(j<4)then
                j1=j
                j2=j+1
            else
                j1=4
                j2=1
            endif
            !!ray: y=dx+C
            !!line:y=Ax+B
            !!cross at point(x,y)
            !!if x1<x<x2, point is on the line
            if(bbox(j1)==bbox(j2))then
                y=d(i)*(bbox(j1)-xx)+yy
                if(bboy(j1)<=y.and.bboy(j2)>=y)then
                    counter=counter+1
                elseif(bbox(j1)>=y.and.bbox(j2)<=y)then
                    counter=counter+1
                endif
            else
                A=(bboy(j1)-bboy(j2))/(bbox(j1)-bbox(j2))
                B=bboy(j1)-A*bbox(j1)
                x=(yy-d(i)*xx-B)/(A-d(i))
                if(bboy(j1)<=x.and.bboy(j2)>=x)then
                    counter=counter+1
                elseif(bbox(j1)>=x.and.bbox(j2)<=x)then
                    counter=counter+1
                endif
            endif
        end do
        !if no cross point for any ray, is outside!
        if(counter==0)then
            NRRgridsort=0
            return
        endif
        if(mod(counter,2)==0)   mm=mm+1
    enddo
    if(mm>=(dnum/2))then
        NRRgridsort=0     !!outside bbox
    else
        NRRgridsort=1
    endif
    endfunction NRRgridsort
!_________________________________________
    !integer function NRRgridsort(xx,yy,bbox,bboy)
    !real(predd),  DIMENSION(4), INTENT (IN) :: bbox,bboy
    !integer,parameter::dnum=5
    !real(predd) xx,yy,a(dnum),t(dnum)
    !real(predd)::d(dnum)=(/1.123,-1.123,2.345,0.577,-3.49/)
    !integer i,j,mm,j1,j2,nn
    !integer::counter
    !
    !mm=0
    !nn=0
    !do i=1,dnum
    !    counter=0
    !    do j=1,4
    !        if(j<4)then
    !            j1=j
    !            j2=j+1
    !        else
    !            j1=4
    !            j2=1
    !        endif
    !        a(i)=(bboy(j2)-yy-(bbox(j2)-xx)*d(i))/((bbox(j1)-xd(j2))*d(i)-bboy(j1)+bboy(j2))
    !        if(a(i)>0.and.a(i)<=1)then
    !            t(i)=bbox(j2)-xx+a(i)*(bbox(j)-bbox(j1))  !分步判断，减少运算，提速
    !            if(t(i)>=0) counter=counter+1
    !        end if
    !    end do
    !    if(mod(counter,2)==0)   mm=mm+1
    !enddo
    !if(mm>=(dnum/2))then
    !    NRRgridsort=0     !!outside bbox
    !else
    !    NRRgridsort=1
    !endif
    !end function NRRgridsort
!_________________________________________
    recursive subroutine NRRregionMark(c)
    type(gridtyp),pointer::c
    real(predd)::Length
    real(predd)::theta
    
    if(associated(c%son1))then
        if(.not.associated(c%son3))then
        call NRRregionMark(c%son1)
        call NRRregionMark(c%son2)
        else
        call NRRregionMark(c%son1)
        call NRRregionMark(c%son2)
        call NRRregionMark(c%son3)
        call NRRregionMark(c%son4)
        endif
        return
    endif
    if(c%NRRregion>=-1) return
    !!inside cell,          c%NRRregion = -1
    if(c%sort==1) return
    !!red cell,             c%NRRregion =  1
    !if(c%lvl==(NRRlvl))then
    !    c%NRRregion=1
    !    return
    !endif
    !!blue cell,            c%NRRregion =  0
    Length=NRRLength*(h/2**(NRRlvl))
    if(WDist(c)<=Length)then
        theta=abs(Rtheta(c))
        if(theta<NRRrefineTheta1.or.theta>NRRrefineTheta2)then
            c%NRRregion=-1
        else
            c%NRRregion=0
            call NRRpointerFinder(c)    !!locate the pointer for the blue cell.
        endif
    else
    !!black cell
        c%NRRregion=-1
    endif

    endsubroutine NRRregionMark
!_________________________________________
    subroutine NRRpointerFinder(c)
    type(gridtyp),pointer::c
    real(predd)::dist
    integer::seed1,seed2
    type(pointdd)::position1,position2

    !c%NRRpointer1=>null()
    !c%NRRpointer2=>null()
    dist=WDist(c)
    call NRRnearestSeed(c,seed1,seed2)
    call NRRlocatePosition(seed1,dist,position1)!!point, dist, export point
    call NRRlocatePosition(seed2,dist,position2)
    call NRRNearestCell(position1%x,position1%y,c%NRRpointer1)
    call NRRNearestCell(position2%x,position2%y,c%NRRpointer2)
    endsubroutine NRRpointerFinder
!_________________________________________
    subroutine NRRnearestSeed(c,seed1,seed2)
    type(gridtyp),pointer::c
    real(predd)::dist1,dist2,mindist
    integer::seed1,seed2,a,b
    integer::i

    mindist=X0*X0+Y0*Y0
    do i=1,NRRSeeds
        dist1=(SeedPoints(i,1)-c%center%x)**2+(SeedPoints(i,2)-c%center%y)**2
        if(mindist>dist1)then
            mindist=dist1
            seed1=i
        endif
    enddo
    if(seed1==NRRSeeds)then
        a=1
        b=NRRSeeds-1
    elseif(seed1==1)then
        a=2
        b=NRRSeeds
    else
        a=seed1+1
        b=seed1-1
    endif
    dist1=(SeedPoints(a,1)-c%center%x)**2+(SeedPoints(a,2)-c%center%y)**2
    dist2=(SeedPoints(b,1)-c%center%x)**2+(SeedPoints(b,2)-c%center%y)**2
    if(dist1<=dist2)then
        seed2=a
    else
        seed2=b
    endif
    endsubroutine NRRnearestSeed
!_________________________________________
  !  subroutine NRRmakePointer(p,seeds,dist,mindist)
  !  type(gridtyp),pointer::p
  !  integer::seeds
  !  real(predd)::dist,mindist
  !  
  !  if(associated(p%son1))then
     !   if(associated(p%son3))then
     !       call NRRmakePointer(p%son1)
     !       call NRRmakePointer(p%son2)
     !       call NRRmakePointer(p%son3)
     !       call NRRmakePointer(p%son4)
        !else
        !    call NRRmakePointer(p%son1)
     !       call NRRmakePointer(p%son2)
  !      endif
  !      return
  !  endif
  !  if(p%NRRregion/=
  !      
  !  
  !  endsubroutine NRRmakePointer
!_________________________________________
    subroutine NRRlocatePosition(seed,dist,position)
    type(pointdd)::position
    integer::seed
    real(predd)::dist
    
    position%x=SeedPoints(seed,1)+dist*NRRnormal(seed)%x
    position%y=SeedPoints(seed,2)+dist*NRRnormal(seed)%y
    endsubroutine NRRlocatePosition
!_________________________________________
    recursive subroutine NRRCheckPointer(c)
    type(gridtyp),pointer::c
    integer::mlvl
    
    if(associated(c%son1))then
        if(associated(c%son3))then
            call NRRCheckPointer(c%son1)
            call NRRCheckPointer(c%son2)
            call NRRCheckPointer(c%son3)
            call NRRCheckPointer(c%son4)
        else
            call NRRCheckPointer(c%son1)
            call NRRCheckPointer(c%son2)
        endif
        return
    endif
    if(c%NRRregion/=0) return
    mlvl=NRRlvl
    if(associated(c%NRRpointer1))then
        if(c%NRRpointer1%lvl/=mlvl) print*, &
            '****error****    NRR Check: Wrong pointer1 level in cell',c%nelment
    else
        print*,'****error****    NRR Check: Not found pointer1.'
    endif
    if(associated(c%NRRpointer2))then
        if(c%NRRpointer2%lvl/=mlvl) print*, &
            '****error****    NRR Check: Wrong pointer2 level in cell',c%nelment
    else
        print*,'****error****    NRR Check: Not found pointer2.'
    endif
    endsubroutine NRRCheckPointer
!_________________________________________
    function Rtheta(c)
    real(predd)::Rtheta
    type(gridtyp),pointer::c
    real(predd)::x,y
    
    x=c%center%x-Rcenter(1)
    y=c%center%y-Rcenter(2)
    if(y==0)then
        if(x<0) Rtheta=0;
        if(x>0) Rtheta=180;
    elseif(x==0)then
        if(y>0) Rtheta=90
        if(y<0) Rtheta=-90
    else
        Rtheta=atan(abs(y/x))/pi*180  !!theta=acrtangent(x/y)
        if(x>0.and.y<0)then !4
            Rtheta=-180+Rtheta
        elseif(x<0.and.y<0)then !3
            Rtheta=-Rtheta
        elseif(x<0.and.y>0)then !2
            Rtheta=Rtheta
        elseif(x>0.and.y>0)then !1
            Rtheta=180-Rtheta
        endif
    endif
    !!turn to the upwind side as Rtheta=0
    !Rtheta=Rtheta+180
    !if(Rtheta>180) Rtheta=360-Rtheta
    endfunction Rtheta
!_________________________________________
    recursive subroutine NRRcrossMeshRefine(c)
    type(gridtyp),pointer::c
    real(predd)::theta

    if(c%cross/=1) return
    if(associated(c%son1))then
        if(associated(c%son3))then
            call NRRcrossMeshRefine(c%son1)
            call NRRcrossMeshRefine(c%son2)
            call NRRcrossMeshRefine(c%son3)
            call NRRcrossMeshRefine(c%son4)
        else
            call NRRcrossMeshRefine(c%son1)
            call NRRcrossMeshRefine(c%son2)
        endif
        return
    endif
    
!    if(c%NRRregion/=-1) return
    if(c%lvl<lvlmax)then
        theta=abs(Rtheta(c))
        if(theta<NRRrefineTheta1.or.theta>NRRrefineTheta2)then
            call newson(c)
            c%son1%NRRregion=-1
            c%son2%NRRregion=-1
            c%son3%NRRregion=-1
            c%son4%NRRregion=-1
            call NRRcrossMeshRefine(c%son1)
            call NRRcrossMeshRefine(c%son2)
            call NRRcrossMeshRefine(c%son3)
            call NRRcrossMeshRefine(c%son4)
        endif
    endif
    endsubroutine NRRcrossMeshRefine
!_________________________________________
    endmodule NRRMeshMod
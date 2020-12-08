module mesh  !变量信息
use inflow
use precisionMod
use typeMod
  implicit none
  integer       m,n,total   !水平竖直、总网格数
  type(gridtyp),pointer::cell(:)
  !外形数据点
  real(pre),allocatable::xd(:),yd(:) 
  !integer dmax
  !记录全区域旋度和散度标准差
  real(pre)::Sgra,Ssecgra,Srot,Sdiv
  integer       Ne     !全区域单元数
  integer       nconver !是否收敛
  integer::gxyxcount=0
  contains

!子函数，对点种类判断
integer function gridsort(xx,yy)
  implicit none
  integer,parameter::dnum=3
  real(pre) xx,yy,a(dnum),t(dnum)  !判断单元类型要用的a，t和射线方向d
  real(pre)::d(dnum)=(/1.1, -1.1, 3.5/)   !斜率 
  integer i,j,mm,nn   
  integer::counter     
      mm=0
      nn=0
   do i=1,dnum
      counter=0
   do j=1,dmax
    if(j<dmax)then
     a(i)=(yd(j+1)-yy-(xd(j+1)-xx)*d(i))/((xd(j)-xd(j+1))*d(i)-yd(j)+yd(j+1))
	 else
	 a(i)=(yd(1)-yy-(xd(1)-xx)*d(i))/((xd(j)-xd(1))*d(i)-yd(j)+yd(1))
	 endif
     if(a(i)>0.and.a(i)<=1)then
    if(j<dmax)then
       t(i)=xd(j+1)-xx+a(i)*(xd(j)-xd(j+1))  !分步判断，减少运算，提速
	   	 else
	 t(i)=xd(1)-xx+a(i)*(xd(j)-xd(1))
	 endif

       if(t(i)>=0)then
       counter=counter+1
       end if
     end if
   end do
 
if(mod(counter,2)==0)then
mm=mm+1
else
nn=nn+1
end if    
end do 

if(mm>nn)then
gridsort=0     !物面外部
else
gridsort=1
end if 

return
end function gridsort  


integer function gridcross(clx,cly,cx,cy) !单元中心 
!   !!程序不严谨，存在以下可能性：四个顶点都在物面外部，但物面侵入了网格某一条边
  implicit none
  type(gridtyp),pointer::c

  real(pre)::cx,cy   
  integer t1,t2,t3,t4,clx,cly
  real(pre)::xx1,yy1,xx2,yy2,xx3,yy3,xx4,yy4
  !四个点从左下到右上1234,
    xx1=cx-h/2**(clx+1)
    yy1=cy-h/(2**(cly+1))
    xx2=cx+h/2**(clx+1)
    yy2=cy-h/(2**(cly+1))
    xx3=cx+h/2**(clx+1)
    yy3=cy+h/(2**(cly+1))
    xx4=cx-h/2**(clx+1)
    yy4=cy+h/(2**(cly+1))
 
    t1=gridsort(xx1,yy1)
    t2=gridsort(xx2,yy2)
    t3=gridsort(xx3,yy3)
    t4=gridsort(xx4,yy4)
   
    if(t1==1.and.t2==1.and.t3==1.and.t4==1)then
    gridcross=2   !相交（1）判断，0是外部，1相交，2内部    
    else if(t1==0.and.t2==0.and.t3==0.and.t4==0)then
    gridcross=0 
    else
    gridcross=1
    end if
end function gridcross

!判断是否需要曲率加密的子函数，1加密，0不加密
!integer function gridcur(cl,cc,cx,cy)  !参数：单元层次，相交，中心坐标x,y
!  implicit none
!  integer i,cl,cc
!  real(pre)::mm,nn
!  real(pre)::cx,cy 
!if(cl>=lvlmax.and.cl<lvlcur.and.cc==1)then  !这里是外形加密，曲率加密
!gridcur=0      !开始先将cur置0，要不然循环完了还没有赋值给它
!    do i=1,dmax-1
!    !判断数据点是否在单元内部（或者右/上边上）
!
!        if(xd(i)>cx-h/(2**(cl+1)).and.xd(i)<=cx+h/(2**(cl+1))&
!        &.and.yd(i)>cy-h/(2**(cl+1)).and.yd(i)<=cy+h/(2**(cl+1)))then
!if(i==1)then
!mm=(xd(1)-xd(dmax-1))*(xd(2)-xd(1))+(yd(1)-yd(dmax-1))*(yd(2)-yd(1))
!nn=sqrt((xd(1)-xd(dmax-1))**2+(yd(1)-yd(dmax-1))**2)*sqrt((xd(2)-xd(1))**2+(yd(2)-yd(1))**2)    
!else
!mm=(xd(i)-xd(i-1))*(xd(i+1)-xd(i))+(yd(i)-yd(i-1))*(yd(i+1)-yd(i))
!nn=sqrt((xd(i)-xd(i-1))**2+(yd(i)-yd(i-1))**2)*sqrt((xd(i+1)-xd(i))**2+(yd(i+1)-yd(i))**2) 
!end if
!!若夹角弧度变化超过0.1（约6度）则曲率加密 
!if(acos(mm/nn)>0.5)then
!gridcur=1
!
!exit   !若判断需要加密则退出循环，而初次判定不需要加密循环仍然进行（针对单元内存在两个或以上数据点）
!else
!gridcur=0
!end if
!        end if  
!    end do
!else
!gridcur=0
!end if
!
!return
!end function gridcur

!基于外形加密，其中调用基于曲率加密
recursive subroutine gridrefm(cell)
  implicit none
  type(gridtyp),pointer::cell
!+++++++++++++++!
!加密方法：先对相交的单元加密三次（全单元加密），
!再对加密后相交的叶子单元加密两次
  if(cell%cross==1.and.cell%lvl<lvlmax)then
      call newson(cell)
      call gridrefm(cell%son1)
      call gridrefm(cell%son2)
      call gridrefm(cell%son3)
      call gridrefm(cell%son4)
  !else if(cell%lvl>0.and.cell%lvl<lvltemp.and.cell%cross/=2)then
  !    call newson(cell)
  !    call gridrefm(cell%son1)
  !    call gridrefm(cell%son2)
  !    call gridrefm(cell%son3)
  !    call gridrefm(cell%son4)
  !elseif(cell%lvl==lvltemp.and.cell%cross==1)then   !第二次统一加密，先找相交单元，再对其统一加密
  !    call newson(cell)                    
  !    call gridrefm(cell%son1)
  !    call gridrefm(cell%son2)
  !    call gridrefm(cell%son3)
  !    call gridrefm(cell%son4)
  !elseif((cell%center%x-X0/4)**2+(cell%center%y-Y0/2)**2<=(Rad+0.06)**2&
  !&.and.cell%lvl<lvlmax.and.cell%cross/=2)then 
  !    call newson(cell)
  !    call gridrefm(cell%son1)
  !    call gridrefm(cell%son2)
  !    call gridrefm(cell%son3)
  !    call gridrefm(cell%son4)
  ! elseif((cell%center%x-X0/4)**2+(cell%center%y-Y0/2)**2<=(Rad+1*h)**2&
  !&.and.cell%lvl<lvltemp.and.cell%sort==0)then 
  !    call newson(cell)
  !    call gridrefm(cell%son1)
  !    call gridrefm(cell%son2)
  !    call gridrefm(cell%son3)
  !    call gridrefm(cell%son4) 
  else
      nullify(cell%son1,cell%son2,cell%son3,cell%son4)
   end if
end subroutine gridrefm



subroutine newson(cell)
implicit none
type(gridtyp),pointer::cell
!每个指针需要给予地址allocate
allocate(cell%son1,cell%son2,cell%son3,cell%son4)
!!不同属性
      if(cell%typ==1)then
          CellCount=CellCount+1
	  cell%son1%lvlx=cell%lvlx
	  cell%son1%lvly=cell%lvly+1
	  cell%son1%spl=1
	  cell%son1%center%x=cell%center%x
	  cell%son1%center%y=cell%center%y-h/(2**(cell%son1%lvly+1))
	  cell%son2%lvlx=cell%lvlx
	  cell%son2%lvly=cell%lvly+1
	  cell%son2%spl=1
	  cell%son2%center%x=cell%center%x
	  cell%son2%center%y=cell%center%y+h/(2**(cell%son2%lvly+1))
      elseif(cell%typ==2)then
          CellCount=CellCount+1
	  cell%son1%lvlx=cell%lvlx+1
	  cell%son1%lvly=cell%lvly
	  cell%son1%spl=2
	  cell%son1%center%x=cell%center%x-h/(2**(cell%son1%lvlx+1))
	  cell%son1%center%y=cell%center%y
	  cell%son2%lvlx=cell%lvlx+1
	  cell%son2%lvly=cell%lvly
	  cell%son2%spl=2
	  cell%son2%center%x=cell%center%x+h/(2**(cell%son1%lvlx+1))
	  cell%son2%center%y=cell%center%y
      elseif(cell%typ==0)then
          CellCount=CellCount+3
	  cell%son1%lvlx=cell%lvlx+1
	  cell%son1%lvly=cell%lvly+1
	  cell%son1%spl=0
	  cell%son1%center%x=cell%center%x-h/(2**(cell%son1%lvlx+1))
	  cell%son1%center%y=cell%center%y-h/(2**(cell%son1%lvly+1))
	  cell%son2%lvlx=cell%lvlx+1
	  cell%son2%lvly=cell%lvly+1
	  cell%son2%spl=0
	  cell%son2%center%x=cell%center%x+h/(2**(cell%son1%lvlx+1))
	  cell%son2%center%y=cell%center%y-h/(2**(cell%son1%lvly+1))
	  cell%son3%lvlx=cell%lvlx+1
	  cell%son3%lvly=cell%lvly+1
	  cell%son3%spl=0
	  cell%son3%center%x=cell%center%x+h/(2**(cell%son1%lvlx+1))
	  cell%son3%center%y=cell%center%y+h/(2**(cell%son1%lvly+1))
	  cell%son4%lvlx=cell%lvlx+1
	  cell%son4%lvly=cell%lvly+1
	  cell%son4%spl=0
	  cell%son4%center%x=cell%center%x-h/(2**(cell%son1%lvlx+1))
	  cell%son4%center%y=cell%center%y+h/(2**(cell%son1%lvly+1))
	  
	  else
	  call stopall('newson')
	  endif
!!共同属性
	  cell%son1%num=cell%num
      cell%son1%lvl=cell%lvl+1
      cell%son1%father=>cell
      cell%son1%sp=1
	  cell%son1%typ=0
      cell%son1%s1=0
      cell%son1%s2=0
      cell%son1%s3=0
      cell%son1%s4=0
      !流场解传递
      cell%son1%u=cell%u
      cell%son1%v=cell%v
      cell%son1%rou=cell%rou
      cell%son1%T=cell%T
      cell%son1%p=cell%p
      cell%son1%grax=cell%grax
      cell%son1%secgrax=cell%secgrax
      cell%son1%gray=cell%gray
      cell%son1%secgray=cell%secgray
      cell%son1%graxp=cell%graxp
      cell%son1%grayp=cell%grayp
      cell%son1%rotx=cell%rotx
      cell%son1%roty=cell%roty
      cell%son1%divx=cell%divx
      cell%son1%divy=cell%divy
	  cell%son1%NRRregion=-999
      
      cell%son2%num=cell%num
      cell%son2%lvl=cell%lvl+1
      cell%son2%father=>cell
      cell%son2%sp=2
      cell%son2%typ=0
      cell%son2%s1=0
      cell%son2%s2=0
      cell%son2%s3=0
      cell%son2%s4=0
      !流场解传递
      cell%son2%u=cell%u
      cell%son2%v=cell%v
      cell%son2%rou=cell%rou
      cell%son2%T=cell%T
      cell%son2%p=cell%p
      cell%son2%grax=cell%grax
      cell%son2%secgrax=cell%secgrax
      cell%son2%gray=cell%gray
      cell%son2%secgray=cell%secgray
      cell%son2%graxp=cell%graxp
      cell%son2%grayp=cell%grayp
      cell%son2%rotx=cell%rotx
      cell%son2%roty=cell%roty
      cell%son2%divx=cell%divx
      cell%son2%divy=cell%divy
      cell%son2%NRRregion=-999
	  
	  cell%son1%sort=gridsort(cell%son1%center%x,cell%son1%center%y)
      cell%son1%cross=gridcross(cell%son1%lvlx,cell%son1%lvly,cell%son1%center%x,cell%son1%center%y)
	  cell%son2%sort=gridsort(cell%son2%center%x,cell%son2%center%y)
	  cell%son2%cross=gridcross(cell%son2%lvlx,cell%son2%lvly,cell%son2%center%x,cell%son2%center%y)
	  call nullifyCell(cell%son1)
	  call nullifyCell(cell%son2)
	  
	  
	  if(cell%typ==0)then
      cell%son3%num=cell%num
      cell%son3%lvl=cell%lvl+1
      cell%son3%father=>cell
      cell%son3%sp=3
	  cell%son3%typ=0
      cell%son3%s1=0
      cell%son3%s2=0
      cell%son3%s3=0
      cell%son3%s4=0
      !流场解传递
      cell%son3%u=cell%u
      cell%son3%v=cell%v
      cell%son3%rou=cell%rou
      cell%son3%T=cell%T
      cell%son3%p=cell%p
      cell%son3%grax=cell%grax
      cell%son3%secgrax=cell%secgrax
      cell%son3%gray=cell%gray
      cell%son3%secgray=cell%secgray
      cell%son3%graxp=cell%graxp!压力梯度
      cell%son3%grayp=cell%grayp
      cell%son3%rotx=cell%rotx!旋度散度
      cell%son3%roty=cell%roty
      cell%son3%divx=cell%divx
      cell%son3%divy=cell%divy
	  cell%son3%NRRregion=-999
      
      cell%son4%num=cell%num
      cell%son4%lvl=cell%lvl+1
      cell%son4%father=>cell
      cell%son4%sp=4
	  cell%son4%typ=0
      cell%son4%s1=0
      cell%son4%s2=0
      cell%son4%s3=0
      cell%son4%s4=0
      !流场解传递
      cell%son4%u=cell%u
      cell%son4%v=cell%v
      cell%son4%rou=cell%rou
      cell%son4%T=cell%T
      cell%son4%p=cell%p
      cell%son4%grax=cell%grax
      cell%son4%secgrax=cell%secgrax
      cell%son4%gray=cell%gray
      cell%son4%secgray=cell%secgray
      cell%son4%graxp=cell%graxp
      cell%son4%grayp=cell%grayp
      cell%son4%rotx=cell%rotx
      cell%son4%roty=cell%roty
      cell%son4%divx=cell%divx
      cell%son4%divy=cell%divy
	  cell%son4%NRRregion=-999
      
	  cell%son3%sort=gridsort(cell%son3%center%x,cell%son3%center%y)
	  cell%son3%cross=gridcross(cell%son3%lvlx,cell%son3%lvly,cell%son3%center%x,cell%son3%center%y)
	  cell%son4%sort=gridsort(cell%son4%center%x,cell%son4%center%y)
	  cell%son4%cross=gridcross(cell%son4%lvlx,cell%son4%lvly,cell%son4%center%x,cell%son4%center%y)
	  call nullifyCell(cell%son3)
	  call nullifyCell(cell%son4)
	  else
	  nullify(cell%son3,cell%son4)
	  endif
end subroutine newson



!---------------------------------------------------
recursive subroutine gridgxyx(cell)
implicit none
type(gridtyp),pointer::cell,c
if(associated(cell%son1))then
	    if(associated(cell%son3))then
	    call gridgxyx(cell%son1)
	    call gridgxyx(cell%son2)
	    call gridgxyx(cell%son3)
	    call gridgxyx(cell%son4)
		else
		call gridgxyx(cell%son1)
	    call gridgxyx(cell%son2)
		endif
!elseif(cell%lvl>=lvltemp.and.cell%lvl<lvlgxyx.and.cell%cross==1)then
!	call typ(cell)
!	call newson(cell)
!	    if(cell%typ==0)then
!	    call gridgxyx(cell%son1)
!	    call gridgxyx(cell%son2)
!	    call gridgxyx(cell%son3)
!	    call gridgxyx(cell%son4)
!		else
!		call gridgxyx(cell%son1)
!	    call gridgxyx(cell%son2)
!        endif
elseif(cell%lvl>=lvltemp.and.cell%lvl<lvlgxyx.and.cell%sort==0)then 
    call typ(cell)
	call newson(cell)
	    if(cell%typ==0)then
	    call gridgxyx(cell%son1)
	    call gridgxyx(cell%son2)
	    call gridgxyx(cell%son3)
	    call gridgxyx(cell%son4)
		else
		call gridgxyx(cell%son1)
	    call gridgxyx(cell%son2)
        endif  
endif
endsubroutine



subroutine typ(cell)  
  implicit none
  type(gridtyp),pointer::cell
  integer i,cl,cc
  real(pre)::mx,my,nn
  real(pre)::cx,cy 
!if(cell%lvl>=lvltemp.and.cell%lvl<lvlgxyx.and.cell%cross/=2)then  !这里是外形加密，曲率加密

    do i=1,dmax-1
    !判断数据点是否在单元内部（或者右/上边上）

        if(xd(i)>cell%center%x-h/(2**(cell%lvl+1)).and.xd(i)<=cell%center%x+h/(2**(cell%lvl+1))&
        &.and.yd(i)>cell%center%y-h/(2**(cell%lvl+1)).and.yd(i)<=cell%center%y+h/(2**(cell%lvl+1)))then
if(i==1)then
mx=abs((xd(1)-xd(dmax-1))*(xd(2)-xd(1)))
my=abs((yd(1)-yd(dmax-1))*(yd(2)-yd(1)))
nn=sqrt((xd(1)-xd(dmax-1))**2+(yd(1)-yd(dmax-1))**2)*sqrt((xd(2)-xd(1))**2+(yd(2)-yd(1))**2)     
else
mx=abs((xd(i)-xd(i-1))*(xd(i+1)-xd(i)))
my=abs((yd(i)-yd(i-1))*(yd(i+1)-yd(i)))
nn=sqrt((xd(i)-xd(i-1))**2+(yd(i)-yd(i-1))**2)*sqrt((xd(i+1)-xd(i))**2+(yd(i+1)-yd(i))**2) 
end if
        end if
    end do
!end if
     

!if(cell%lvl>lvltemp)then
!cell%typ=cell%father%typ
!return
!endif

!若夹角弧度变化超过0.1（约6度）则曲率加密 
if(abs(acos(mx/nn))>0.9.and.abs(acos(my/nn))<0.9)then
         cell%typ=2
        gxyxcount=gxyxcount+1  !x方向加密
      elseif(abs(acos(mx/nn))<0.9.and.abs(acos(my/nn))>0.9)then
        cell%typ=1
        gxyxcount=gxyxcount+1        !y方向加密
      elseif(abs(acos(mx/nn))>=0.9.and.abs(acos(my/nn))>=0.9)then
         cell%typ=0
      end if

return
end subroutine


!圆
!subroutine typ(c)
!implicit none
!type(gridtyp),pointer::c
!real(pre)::arc,xx,yy  !弧度，相对圆心的x/y坐标
!integer::quad   !象限
!!以30度夹角为分类
!xx=c%center%x-Rcenter(1)
!yy=c%center%y-Rcenter(2)
!if(xx>=0.and.yy>=0) quad=1
!if(xx<0.and.yy>=0)  quad=2
!if(xx<0.and.yy<0)   quad=3
!if(xx>=0.and.yy<0)  quad=4
!xx=abs(xx)
!yy=abs(yy)
!arc=asin(yy/sqrt(xx**2+yy**2))
!!print*,"arc:", arc
!if(c%lvl>lvltemp)then
!c%typ=c%father%typ
!return
!endif
!
!if(arc<=(pi/6))then
!c%typ=2
!gxyxcount=gxyxcount+1
!!write(*,*) "gxyx", gxyxcount, c%num, c%typ
!elseif(arc<=(pi/3))then
!c%typ=0
!elseif(arc<=(pi/2))then
!c%typ=1
!gxyxcount=gxyxcount+1
!!write(*,*) "gxyx", gxyxcount, c%num, c%typ
!endif
!endsubroutine

!!三角
!subroutine typ(c)
!implicit none
!type(gridtyp),pointer::c
!real(pre)::arc,xx,yy  !弧度，相对圆心的x/y坐标
!integer::quad   !象限
!!以30度夹角为分类
!if(c%lvl>lvltemp)then
!c%typ=c%father%typ
!return
!endif
!
!if(c%center%x>=5.8)then
!    if(c%center%y>9.7.and.c%center%y<10.3)then
!    c%typ=2
!    gxyxcount=gxyxcount+1
!	else
!	c%typ=0
!	endif
!!write(*,*) "gxyx", gxyxcount, c%num, c%typ
!elseif(c%center%x<4.2)then
!c%typ=0
!else
!c%typ=1
!gxyxcount=gxyxcount+1
!!write(*,*) "gxyx", gxyxcount, c%num, c%typ
!endif
!endsubroutine
!_________________________________________
    subroutine initialMesh
    integer::k1,k11,i
    real(pre)::x,y
    do i=1,total
        if(mod(i,m)/=0)then
        k1=i/m
        else
        k1=i/m-1
        end if
        k11=i-m*k1
        x=0.5*(2*k11-1)*h
        y=(k1+0.5)*h

        cell(i)%num=i
        cell(i)%lvl=0
        cell(i)%center%x=x
        cell(i)%center%y=y
        cell(i)%cur=0
        cell(i)%sp=0
        cell(i)%grax=0
        cell(i)%gray=0
        cell(i)%secgrax=0
        cell(i)%secgray=0
        cell(i)%graxp=0
        cell(i)%grayp=0
        cell(i)%rotx=0
        cell(i)%roty=0
        cell(i)%divx=0
        cell(i)%divy=0
        cell(i)%lvlx=0
        cell(i)%lvly=0
        cell(i)%typ=0
        cell(i)%NRRregion=-1
        cell(i)%spl=0
        cell(i)%father=>NULL()
        cell(i)%sort=gridsort(x,y)
        cell(i)%cross=gridcross(cell(i)%lvlx,cell(i)%lvly,x,y)
        nullify(cell(i)%son1,cell(i)%son2,cell(i)%son3,cell(i)%son4)
    enddo
    endsubroutine initialMesh
!_________________________________________
subroutine nullifyCell(cell)
type(gridtyp),pointer::cell

nullify(cell%son1,cell%son2,cell%son3,cell%son4)
endsubroutine
!_________________________________________
function WDist(c)
real(pre)::WDist
type(gridtyp),pointer::c
integer::i
real(pre)::mindist1,mindist2,dist1,dist2
real(pre)::crosspointx,crosspointy,k,a,b,normal
integer::geop1,geop2,p1,p2

mindist1=sqrt(X0*X0+Y0*Y0)
do i=1,dmax
    dist1=sqrt((xd(i)-c%center%x)**2+(yd(i)-c%center%y)**2)
    if(dist1<minDist1)then
        minDist1=dist1
        geop1=i
    endif
enddo
if(geop1==1)then
    p1=dmax
    p2=2
elseif(geop1==dmax)then
    p1=dmax-1
    p2=1
else
    p1=geop1-1
    p2=geop1+1
endif
dist1=sqrt((xd(p1)-c%center%x)**2+(yd(p1)-c%center%y)**2)
dist2=sqrt((xd(p2)-c%center%x)**2+(yd(p2)-c%center%y)**2)
if(dist1>dist2)then
    geop2=p2
    minDist2=dist2
else
    geop2=p1
    minDist2=dist1
endif
!dist from cell to line geop1-geop2
!normal direction
if(xd(geop1)==xd(geop2))then
    if(yd(geop1)<c%center%y.and.c%center%y<yd(geop2))then
        WDist=abs(xd(geop1)-c%center%x)
    elseif(yd(geop1)>c%center%y.and.c%center%y>yd(geop2))then
        WDist=abs(xd(geop1)-c%center%x)
    else
        WDist=min(mindist1,mindist2)
    endif
elseif(yd(geop1)==yd(geop2))then
    if(xd(geop1)<c%center%x.and.c%center%x<xd(geop2))then
        WDist=abs(yd(geop1)-c%center%y)
    elseif(xd(geop1)>c%center%x.and.c%center%x>xd(geop2))then
        WDist=abs(yd(geop1)-c%center%y)
    else
        WDist=min(mindist1,mindist2)
    endif
else
    !!yy=xx*normal+a
    !!y=kx+b
    k=(xd(geop1)-xd(geop2))/(yd(geop1)-yd(geop2))
    b=yd(geop1)-xd(geop1)*k
    normal=-1/k
    a=c%center%y-c%center%x*normal
    crosspointx=(a-b)/(n-k)
    if(xd(geop1)<crosspointx.and.crosspointx<xd(geop2))then
        crosspointy=crosspointx*k+b
        WDist=sqrt((crosspointx-c%center%x)**2+(crosspointy-c%center%y)**2)
    elseif(xd(geop1)>crosspointx.and.crosspointx>xd(geop2))then
        crosspointy=crosspointx*k+b
        WDist=sqrt((crosspointx-c%center%x)**2+(crosspointy-c%center%y)**2)
    else
        WDist=min(mindist1,mindist2)
    endif
endif
endfunction WDist
!---------------------------------------------------

end module mesh
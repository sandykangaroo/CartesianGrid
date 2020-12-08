module flow_result
use mesh
use precisionMod
use typeMod
implicit none

contains

subroutine surface_data(char)
integer i 
character(*)::char
character(40)::filename
real(pre)::phi,Cp,Vt

filename='./Others/'//trim(FileNameStr)//'-Surf-'//char//'.dat'
open(unit=88,file=filename) 
write(88,*)'"x","y","alfa","rou","Cp","Vt"'
do i=1,dmax
phi=360.0/(dmax-1)*i   
call get_data(xd(i),yd(i),phi,Cp,Vt)  !rou_t为相对来流
!数据点序号可以得到角度，针对圆
write(88,'(5F20.7)')xd(i),yd(i),phi,Cp,Vt   
end do
close(88)
print*,"Save to file:     ",filename
return
end subroutine surface_data
!确定数据点所在的叶子单元，由之得到相应的值

subroutine get_data(xx,yy,phi,Cp,Vt)
real(pre)::xx,yy,x1,x2,y1,y2,Cp,Vt,phi
type(gridtyp),pointer::ct,c
integer i

do i=1,total
  c=>cell(i)
  x1=c%center%x-h/2
  x2=c%center%x+h/2
  y1=c%center%y-h/2
  y2=c%center%y+h/2
  if(xx>=x1.and.xx<=x2.and.yy>=y1.and.yy<=y2)then
   call get_data_sub(c,xx,yy,phi,Cp,Vt)  !由大的子单元递归到具体的叶子单元，返回rou
  endif
end do
return
end subroutine get_data 

recursive subroutine get_data_sub(c,xx,yy,phi,Cp,Vt)
real(pre)::xx,yy,Cp,arc,phi,Vt
type(gridtyp),pointer::c

if(.not.associated(c%son1).and.xx>=c%center%x-h/2**(c%lvlx+1).&
&and.xx<=c%center%x+h/2**(c%lvlx+1).&
&and.yy>=c%center%y-h/2**(c%lvly+1).and.yy<=c%center%y+h/2**(c%lvly+1))then  
!这里没有讨论所在单元在物面内外情况，内部的话可以根据位置分别讨论其邻居的值代替
!rou=c%rou/rou0  !rou0远场的密度
!表面压力系数
arc=phi/360*2*pi
Cp=(c%p-R*rou0*T0)/(0.5*rou0*u0**2)
Vt=c%u*sin(arc)+c%v*cos(arc)
elseif(xx>=c%center%x-h/2**(c%lvlx+1).and.xx<=c%center%x+h/2**(c%lvlx+1).&
&and.yy>=c%center%y-h/2**(c%lvly+1).and.yy<=c%center%y+h/2**(c%lvly+1))then
    if(associated(c%son3))then  
        call get_data_sub(c%son1,xx,yy,phi,Cp,Vt)
        call get_data_sub(c%son2,xx,yy,phi,Cp,Vt)
        call get_data_sub(c%son3,xx,yy,phi,Cp,Vt)
        call get_data_sub(c%son4,xx,yy,phi,Cp,Vt)  
    else
        call get_data_sub(c%son1,xx,yy,phi,Cp,Vt)
        call get_data_sub(c%son2,xx,yy,phi,Cp,Vt)   
    endif
endif
return
end subroutine get_data_sub

function GetVt(c,phi)   !!切向速度，以沿切向为正
real(pre)::GetVt
type(gridtyp),pointer::c
real(pre)::phi

if(c%center%x<0.and.c%center%y>=0)then   !!2
    GetVt=c%u*sin(phi)+c%v*cos(phi)
elseif(c%center%x>=0.and.c%center%y>=0)then !!1
    
elseif(c%center%x>=0.and.c%center%y<0)then  !!4
    
elseif(c%center%x<0.and.c%center%y<0)then   !!3
    
endif

endfunction GetVt
end module
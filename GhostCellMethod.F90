module bc  !�߽紦��ģ��
use frame
use inflow
use precisionMod
use typeMod
implicit none
    contains
    
!���ⵥԪ����,�������ⵥԪ��������ⵥԪ��������
subroutine ghost_cell(cel)   
  implicit none              
  integer i,j,i0
  real(pre)::us,vs,ps,rous,Ts  !�ԳƵ�Ԫ��������
  type(gridtyp),pointer::celt,cel,cels,celsr,celsd,celsrd,celsu,celsru,cels1,cels2,cels3
  !���ⵥԪcell���ԳƵ�Ԫcells���ԳƵ�Ԫ�����ھ�cellsr�����ھӵ����ھ�cellsrd,cellt��ʱ
  real(pre)::xs,ys,xo,yo,xm,ym,hl,d,a,alfa  !�ԳƵ������,xo,yo
   !��ʱ��������t����������T��ͻ����a��

    xo=cel%center%x
    yo=cel%center%y  
    !hl=h/2**cel%lvl   !����λ�ÿɵ���
    hl=0
if(yo>=Y0/2)then    
  do j=1,dmax-1
   a=((xd(j)-xd(j+1))*(xd(j)-xo)+(yd(j)-yd(j+1))*(yd(j)-yo))/&
   &((xd(j)-xd(j+1))**2+(yd(j)-yd(j+1))**2)
   if(a>=0.and.a<1)then
     d=sqrt((a*xd(j+1)-a*xd(j)+xd(j)-xo)**2+(a*yd(j+1)-a*yd(j)+yd(j)-yo)**2)
     xs=(hl/d+2)*(a*xd(j+1)-a*xd(j)+xd(j)-xo)+xo    !!symmetry position
     ys=(hl/d+2)*(a*yd(j+1)-a*yd(j)+yd(j)-yo)+yo
     i0=cel%num-mod(cel%num,m)+1   !�ҶԳƵ����ڵ�Ԫ����ʼλ��
   !  nn=1
     exit
   end if
  end do
  !!find the nearest cell to the symmetry position.
  do i=i0,total                      
    if(abs(cell(i)%center%x-xs)<=h/2.and.abs(cell(i)%center%y-ys)<=h/2)then
      celt=>cell(i)
      exit
    end if
  end do
  call cellsym(xs,ys,celt,cels)   !�õ��ԳƵ�Ԫcels
  cels1=>left_neighbor(cels)
  cels2=>up_neighbor(cels)
  cels3=>right_neighbor(cels)
 else    !��������²������
 do j=dmax,2,-1
   a=((xd(j-1)-xd(j))*(xd(j-1)-xo)+(yd(j-1)-yd(j))*(yd(j-1)-yo))/&
   &((xd(j-1)-xd(j))**2+(yd(j-1)-yd(j))**2)
   if(a>=0.and.a<1)then
     d=sqrt((a*xd(j)-a*xd(j-1)+xd(j-1)-xo)**2+(a*yd(j)-a*yd(j-1)+yd(j-1)-yo)**2)
     xs=(hl/d+2)*(a*xd(j)-a*xd(j-1)+xd(j-1)-xo)+xo 
     ys=(hl/d+2)*(a*yd(j)-a*yd(j-1)+yd(j-1)-yo)+yo 
     i0=cel%num+m-mod(cel%num,m)
  !   nn=2   
     exit
   end if
 end do
  do i=i0,1,-1
    if(abs(cell(i)%center%x-xs)<=h/2.and.abs(cell(i)%center%y-ys)<=h/2)then
      celt=>cell(i)
      exit
    end if
  end do
  call cellsym(xs,ys,celt,cels)
  cels1=>left_neighbor(cels)
  cels2=>down_neighbor(cels)
  cels3=>right_neighbor(cels)
end if
!cels,cels1,cels2���������㵽�ԳƵ��ֵ�õ�
call interpose_sym(xs,ys,cels,cels1,cels2,cels3,us,vs,ps,rous,Ts)  !�����Ϊ������
!������������������ ST����
call ST_method(2*d+hl,xo,yo,xs,ys,us,vs,ps,rous,cel%ut,cel%vt,cel%rout,cel%pt,cel%Tt) 
                                            !������cel%u,cel%v,cel%rou�ȸ�ֵ
return  
end subroutine ghost_cell

!�����������ڳ�ʼ��Ԫ������Ҷ�ӽ��
recursive subroutine cellsym(xs,ys,celt,cels) !�ԳƵ㣬��ʱ��Ԫ�ͶԳƵ�Ԫ
  implicit none
  real(pre)::xs,ys
  type(gridtyp),pointer::cel,cels,celt 
  if((.not.associated(celt%son1)).and.abs(celt%center%x-xs)<=h/2**(celt%lvlx+1)&
     &.and.abs(celt%center%y-ys)<=h/2**(celt%lvly+1))then
     cels=>celt
     return
  else if(abs(celt%center%x-xs)<=h/2**(celt%lvlx+1)&
     &.and.abs(celt%center%y-ys)<=h/2**(celt%lvly+1))then
         if(associated(celt%son3))then
             call cellsym(xs,ys,celt%son1,cels)
             call cellsym(xs,ys,celt%son2,cels)
             call cellsym(xs,ys,celt%son3,cels)
             call cellsym(xs,ys,celt%son4,cels)
         else
             call cellsym(xs,ys,celt%son1,cels)
             call cellsym(xs,ys,celt%son2,cels)
         endif
         
  end if 
return
end subroutine cellsym  
!��ԳƵ������ֵ
subroutine interpose_sym(xs,ys,c1,c2,c3,c4,us,vs,ps,rous,Ts)
  implicit none
  real(pre)::xs,ys,us,vs,ps,rous,Ts
  type(gridtyp),pointer::c1,c2,c3,c4,c0
  real(pre)::e1,e2,e3,e4,d1,d2,d3,d4
  d1=sqrt((xs-c1%center%x)**2+(ys-c1%center%y)**2)
  d2=sqrt((xs-c2%center%x)**2+(ys-c2%center%y)**2)
  d3=sqrt((xs-c3%center%x)**2+(ys-c3%center%y)**2)
  d4=sqrt((xs-c4%center%x)**2+(ys-c4%center%y)**2)
  e1=exp(-d1)
  e2=exp(-d2)
  e3=exp(-d3)
  e4=exp(-d4)
  !�Բ�ֵ����ͨ���ӵ�Ԫ
  call interpose_slf(c2)
  call interpose_slf(c3)
  call interpose_slf(c4)
  us=(e1*c1%u+e2*c2%u+e3*c3%u+e4*c4%u)/(e1+e2+e3+e4)
  vs=(e1*c1%v+e2*c2%v+e3*c3%v+e4*c4%v)/(e1+e2+e3+e4) 
  rous=(e1*c1%rou+e2*c2%rou+e3*c3%rou+e4*c4%rou)/(e1+e2+e3+e4)
  Ts=(e1*c1%T+e2*c2%T+e3*c3%T+e4*c4%T)/(e1+e2+e3+e4)
  ps=R*rous*Ts
return
end subroutine interpose_sym

!�Ӻ������Բ�ֵ����ͨ���ӵ�Ԫ,�õ�cell�еĸ���������
recursive subroutine interpose_slf(c)
implicit none
  type(gridtyp),pointer::c
  real(pre)::e1,e2,e3,e4,d1,d2,d3,d4
 !û���ӽ������ 
  if(.not.associated(c%son1)) return 
  if(associated(c%son1)) then 
    if(associated(c%son3)) then 
        if(associated(c%son1%son1)) call interpose_slf(c%son1)
        if(associated(c%son2%son1)) call interpose_slf(c%son2)
        if(associated(c%son3%son1)) call interpose_slf(c%son3)
        if(associated(c%son4%son1)) call interpose_slf(c%son4)
    else
        if(associated(c%son1%son1)) call interpose_slf(c%son1)
        if(associated(c%son2%son1)) call interpose_slf(c%son2)
    endif
  endif
  !��ֵ��ʽ
  if(associated(c%son3)) then
  d1=sqrt((c%center%x-c%son1%center%x)**2+(c%center%y-c%son1%center%y)**2)
  d2=sqrt((c%center%x-c%son2%center%x)**2+(c%center%y-c%son2%center%y)**2)
  d3=sqrt((c%center%x-c%son3%center%x)**2+(c%center%y-c%son3%center%y)**2)
  d4=sqrt((c%center%x-c%son4%center%x)**2+(c%center%y-c%son4%center%y)**2)
  e1=exp(-d1)
  e2=exp(-d2)
  e3=exp(-d3)
  e4=exp(-d4)
  !��һ��ϵ������coefficient,�����ڲ�sort=1������ϵ��0
  c%u=(e1*c%son1%u+e2*c%son2%u+e3*c%son3%u+e4*c%son4%u)/(e1+e2+e3+e4)
  c%v=(e1*c%son1%v+e2*c%son2%v+e3*c%son3%v+e4*c%son4%v)/(e1+e2+e3+e4)
  c%rou=(e1*c%son1%rou+e2*c%son2%rou+e3*c%son3%rou+e4*c%son4%rou)/(e1+e2+e3+e4) 
  c%T=(e1*c%son1%T+e2*c%son2%T+e3*c%son3%T+e4*c%son4%T)/(e1+e2+e3+e4)
  c%p=R*c%rou*c%T
  else
      c%rou=(c%son1%rou+c%son2%rou)/2
      c%u=(c%son1%u+c%son2%u)/2
      c%v=(c%son1%v+c%son2%v)/2
      c%T=(c%son1%T+c%son2%T)/2
      c%p=R*c%rou*c%T   
   end if
  
return
end subroutine interpose_slf


!���ⵥԪ������ʽ
!subroutine ghost_cell_method(alfa,Rc,dcs,us,vs,ps,rous,Ts,u,v,rou,p,T,nn) 
!                            !�нǣ���������ʰ뾶������ԳƵ����
!  implicit none
!  integer nn
!  real(pre)::alfa,Rc,dcs,us,vs,ps,rous,Ts,us_t,vt,u,v,rou,p,T,alfa_u,alfa_v
!  !u��vΪˮƽ����ֱ������us��vsҲ�ǣ�
!  alfa_u=alfa-pi/2    !u������н�
!  alfa_v=pi-alfa
!  u=-us
!  v=-vs
!  if(nn==1)then  
!  us_t=vs*cos(alfa_v)+us*cos(alfa_u) !��pҪ�õ������ٶ�us_t
!  elseif(n==2)then
!  us_t=vs*cos(alfa_v)-us*cos(alfa_u)  
!  end if
!  p=ps-(2/3)*rous*(us_t**2)*dcs/Rc   
!  T=Ts               !������õ��±�ģ��
!  rou=p/(R*T)
!return  
!end subroutine ghost_cell_method  
 !�ԳƷ���߽�����
 subroutine ST_method(dcs,xo,yo,xs,ys,us,vs,ps,rous,u,v,rou,p,T)
 implicit none
 real(pre)::dcs,xo,yo,xs,ys,us,vs,ps,rous,u,v,rou,p,T,a
 
 
 p=ps
 rou=rous
 u=-us
 v=-vs
 T=p/(R*rou) 
 return
 end subroutine ST_method
 
end module bc
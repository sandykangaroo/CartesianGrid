!空间离散的NND和中心差分算法
module nndcenter
  use frame
use precisionMod
use typeMod
!  use sutherland 
  use mesh
  implicit none
contains

subroutine nnd_center(c)
implicit none
integer i
type(gridtyp),pointer::c
real(pre)::U10,U20,U30,U40,U11,U21,U31,U41
real(pre)::fu1,fu2,fu3,fu4!,tstep

!tstep=4D-6
!tstep=(h/2**c%lvl)*Courant
!记录下未时间推进前的分量，对应RK法里的U(n)，因为第一步会把c的值覆盖
U10=c%rou
U20=c%rou*c%u
U30=c%rou*c%v
U40=c%rou*(R*c%T/(gama-1)+(c%u*c%u+c%v*c%v)/2)
!i=1第一步更新,=2第二步,RK方法
!do i=1,2
call Res_nnd_center(c,fu1,fu2,fu3,fu4)   !残差的差分格式，返回四个分量各自差分结果
!if(i==1)then
U11=U10+tstep*fu1
U21=U20+tstep*fu2
U31=U30+tstep*fu3
U41=U40+tstep*fu4
!else
!U11=(U10+U11+tstep*fu1)/2
!U21=(U20+U21+tstep*fu2)/2
!U31=(U30+U31+tstep*fu3)/2
!U41=(U40+U41+tstep*fu4)/2
!end if
c%rout=U11
c%ut=U21/U11
c%vt=U31/U11
c%Tt=(U41/U11-((U21/U11)**2+(U31/U11)**2)/2)*(gama-1)/R
c%pt=R*c%rout*c%Tt
!end do 
return
end subroutine nnd_center

!应用nnd和中心格式求解
subroutine Res_nnd_center(c,fu1,fu2,fu3,fu4)
  implicit none
  type(gridtyp),pointer::c
  real(pre)::fu1,fu2,fu3,fu4
  real(pre)::dEv1,dEv2,dEv3,dEv4,dFv1,dFv2,dFv3,dFv4
  
  !call centerx(c,dEv1,dEv2,dEv3,dEv4)   !水平方向，输入c，返回四个差分后的量,
  !call centery(c,dFv1,dFv2,dFv3,dFv4)   !v代表miou，粘性
  
  fu1=(-1)*nndx(c,1)+(-1)*nndy(c,1) !+dEv1+dFv1
  fu2=(-1)*nndx(c,2)+(-1)*nndy(c,2) !+dEv2+dFv2
  fu3=(-1)*nndx(c,3)+(-1)*nndy(c,3) !+dEv3+dFv3
  fu4=(-1)*nndx(c,4)+(-1)*nndy(c,4) !+dEv4+dFv4
  return
end subroutine Res_nnd_center


!###############################
!粘性项中心差分，先看x方向
!subroutine centerx(c,dEv1,dEv2,dEv3,dEv4)
!  implicit none
!  type(gridtyp),pointer::c
!  real(pre)::dEv1,dEv2,dEv3,dEv4,vis,kk
!  character(6)::ua="u",va="v",xa="x",ya="y",xx="xx",yy="yy",xy="xy"  !加上a和module里的相同字母区分
!  call dynvis(c%T,vis)
!  call thermc(c%T,kk)
!  dEv1=0
!  dEv2=(1/3)*dvis(c,xa)*(4*duv(c,ua,xa)-2*duv(c,va,ya))+(1/3)*vis*&
!  &(4*duv(c,ua,xx)-2*duv(c,va,xy))
!  dEv3=dvis(c,xa)*(duv(c,va,xa)+duv(c,ua,ya))+vis*(duv(c,va,xx)+duv(c,ua,xy))
!  dEv4=(1/3)*vis*duv(c,ua,xa)*(4*duv(c,ua,xa)-2*duv(c,va,ya))+c%u*dEv2+duv(c,va,xa)*vis*&
!  &(duv(c,ua,ya)+duv(c,va,xa))+c%v*dEv3+dk(c,xa)*dT(c,xa)+kk*dT(c,xx)
!return
!end subroutine centerx
!!centery
!subroutine centery(c,dFv1,dFv2,dFv3,dFv4)
!  implicit none
!  type(gridtyp),pointer::c
!  real(pre)::dFv1,dFv2,dFv3,dFv4,vis,kk
!  character(6)::ua="u",va="v",xa="x",ya="y",xx="xx",yy="yy",xy="xy"
!  
!  call dynvis(c%T,vis)
!  call thermc(c%T,kk)
!  dFv1=0
!  dFv2=dvis(c,ya)*(duv(c,va,xa)+duv(c,ua,ya))+vis*(duv(c,va,xy)+duv(c,ua,yy))
!  dFv3=(1/3)*dvis(c,ya)*(4*duv(c,va,ya)-2*duv(c,ua,xa))+(1/3)*vis*&
!  &(4*duv(c,va,yy)-2*duv(c,ua,xy))
!  dFv4=vis*duv(c,ua,ya)*(duv(c,va,xa)+duv(c,ua,ya))+c%u*dFv2+(1/3)*duv(c,va,ya)*vis*&
!  &(4*duv(c,va,ya)-2*duv(c,ua,xa))+c%v*dFv3+dk(c,ya)*dT(c,ya)+kk*dT(c,yy)
!return
!end subroutine centery
!vis各向导数，2个
!function dvis(c,a)
!implicit none
!type(gridtyp),pointer::c
!character(6) a    
!real(pre)::dvis,dx,dy
!real(pre)::roui1,ui1,vi1,Ti1,pi1,roui_1,ui_1,vi_1,Ti_1,pi_1,&
!              &rouj1,uj1,vj1,Tj1,pj1,rouj_1,uj_1,vj_1,Tj_1,pj_1
!!dx=h/2**c%lvl
!!dy=dx
!dx=h/2**c%lvlx
!dy=h/2**c%lvly
!call Tip(c,roui1,ui1,vi1,Ti1,pi1)
!call Tim(c,roui_1,ui_1,vi_1,Ti_1,pi_1)
!call Tjp(c,rouj1,uj1,vj1,Tj1,pj1)
!call Tjm(c,rouj_1,uj_1,vj_1,Tj_1,pj_1)
!if(a=="x")then
!  dvis=dvisdt(c%T)*(Ti1-Ti_1)/(2*dx)
!elseif(a=="y")then
!  dvis=dvisdt(c%T)*(Tj1-Tj_1)/(2*dy)
!end if
!return
!end function dvis
!k的各向导数，类似粘性
!function dk(c,a)
!implicit none
!type(gridtyp),pointer::c
!character(6) a 
!real(pre)::dk,dx,dy
!real(pre)::roui1,ui1,vi1,Ti1,pi1,roui_1,ui_1,vi_1,Ti_1,pi_1,&
!              &rouj1,uj1,vj1,Tj1,pj1,rouj_1,uj_1,vj_1,Tj_1,pj_1
!!dx=h/2**c%lvl
!!dy=dx
!dx=h/2**c%lvlx
!dy=h/2**c%lvly
!call Tip(c,roui1,ui1,vi1,Ti1,pi1)
!call Tim(c,roui_1,ui_1,vi_1,Ti_1,pi_1)
!call Tjp(c,rouj1,uj1,vj1,Tj1,pj1)
!call Tjm(c,rouj_1,uj_1,vj_1,Tj_1,pj_1)
!if(a=="x")then
!  dk=dkdt(c%T)*(Ti1-Ti_1)/(2*dx)
!elseif(a=="y")then
!  dk=dkdt(c%T)*(Tj1-Tj_1)/(2*dy)
!end if
!return
!end function dk
!!T的各向导数，4个
!function dT(c,a)  !和时间步dt冲突，时间步改为tstep
!implicit none
!type(gridtyp),pointer::c
!character(6) a
!real(pre)::dT,dx,dy
!real(pre)::roui1,ui1,vi1,Ti1,pi1,roui_1,ui_1,vi_1,Ti_1,pi_1,&
!              &rouj1,uj1,vj1,Tj1,pj1,rouj_1,uj_1,vj_1,Tj_1,pj_1
!!dx=h/2**c%lvl
!!dy=dx
!dx=h/2**c%lvlx
!dy=h/2**c%lvly
!
!call Tip(c,roui1,ui1,vi1,Ti1,pi1)
!call Tim(c,roui_1,ui_1,vi_1,Ti_1,pi_1)
!call Tjp(c,rouj1,uj1,vj1,Tj1,pj1)
!call Tjm(c,rouj_1,uj_1,vj_1,Tj_1,pj_1)
!    if(a=="x")then     !dT/dx
!    dT=(Ti1-Ti_1)/(2*dx) 
!    elseif(a=="y")then  !dT/dy
!    dT=(Tj1-Tj_1)/(2*dy)
!    elseif(a=="xx")then
!    dT=(Ti1+Ti_1-2*c%T)/dx**2
!    elseif(a=="yy")then
!    dT=(Tj1+Tj_1-2*c%T)/dy**2
!    end if
!return
!end function dT
!u,v各向导数
function duv(c,a,b)  !a,b定类型，返回函数名
implicit none
type(gridtyp),pointer::c
character(6) a,b,u,v
real(pre)::duv,dx,dy
real(pre)::roui1,ui1,vi1,Ti1,pi1,roui_1,ui_1,vi_1,Ti_1,pi_1,&
              &rouj1,uj1,vj1,Tj1,pj1,rouj_1,uj_1,vj_1,Tj_1,pj_1,&
              &ma1i1,ma1i_1,ma1j1,ma1j_1
!dx=h/2**c%lvl
!dy=dx
dx=h/2**c%lvlx
dy=h/2**c%lvly

call Tip(c,roui1,ui1,vi1,Ti1,pi1)
call Tim(c,roui_1,ui_1,vi_1,Ti_1,pi_1)
call Tjp(c,rouj1,uj1,vj1,Tj1,pj1)
call Tjm(c,rouj_1,uj_1,vj_1,Tj_1,pj_1)



  if(a=="u")then
    if(b=="x")then     !du/dx
    duv=(ui1-ui_1)/(2*dx) 
    elseif(b=="y")then  !du/dy
    duv=(uj1-uj_1)/(2*dy)
    elseif(b=="xx")then
    duv=(ui1+ui_1-2*c%u)/dx**2
    elseif(b=="xy")then
    duv=(Ucn(c,ui1,uj1,c%u)+Ucn(c,ui_1,uj_1,c%u)-&
    &Ucn(c,ui_1,uj1,c%u)-Ucn(c,ui1,uj_1,c%u))/(4*dx*dy)
    elseif(b=="yy")then
    duv=(uj1+uj_1-2*c%u)/dy**2
    end if
  else if(a=="v")then
    if(b=="x")then     !dv/dx
    duv=(vi1-vi_1)/(2*dx) 
    elseif(b=="y")then  !dv/dy
    duv=(vj1-vj_1)/(2*dy)
    elseif(b=="xx")then
    duv=(vi1+vi_1-2*c%v)/dx**2
    elseif(b=="xy")then
    duv=(Ucn(c,vi1,vj1,c%v)+Ucn(c,vi_1,vj_1,c%v)-&
    &Ucn(c,vi_1,vj1,c%v)-Ucn(c,vi1,vj_1,c%v))/(4*dx*dy)
    elseif(b=="yy")then
    duv=(vj1+vj_1-2*c%v)/dy**2
    end if
   
   else if(a=="p")then
    if(b=="x")then     !dp/dx
    duv=(pi1-pi_1)/(2*dx) 
    elseif(b=="y")then  !dp/dy
    duv=(pj1-pj_1)/(2*dy)
    elseif(b=="xx")then
    duv=(pi1+pi_1-2*c%p)/dx**2
    elseif(b=="xy")then
    duv=(Ucn(c,pi1,pj1,c%p)+Ucn(c,pi_1,pj_1,c%p)-&
    &Ucn(c,pi_1,pj1,c%p)-Ucn(c,pi1,pj_1,c%p))/(4*dx*dy)
    elseif(b=="yy")then
    duv=(pj1+pj_1-2*c%p)/dy**2
    end if 
    
   else if(a=="ma1")then
       
    ma1i1=sqrt(ui1**2+vi1**2)/sqrt(gama*R*Ti1)!当地马赫数
    ma1i_1=sqrt(ui_1**2+vi_1**2)/sqrt(gama*R*Ti_1)
    ma1j1=sqrt(uj1**2+vj1**2)/sqrt(gama*R*Tj1)
    ma1j_1=sqrt(uj_1**2+vj_1**2)/sqrt(gama*R*Tj_1)
    
    if(b=="x")then     !dp/dx
    duv=(ma1i1-ma1i_1)/(2*dx) 
    elseif(b=="y")then  !dp/dy
    duv=(ma1j1-ma1j_1)/(2*dy)
    elseif(b=="xx")then
    duv=(ma1i1+ma1i_1-2*sqrt(c%u**2+c%v**2)/sqrt(gama*R*c%T))/dx**2
    elseif(b=="yy")then
    duv=(ma1j1+ma1j_1-2*sqrt(c%u**2+c%v**2)/sqrt(gama*R*c%T))/dy**2
    end if
    
    else if(a=="rou")then
    if(b=="x")then     !dp/dx
    duv=(roui1-roui_1)/(2*dx) 
    elseif(b=="y")then  !dp/dy
    duv=(rouj1-rouj_1)/(2*dy)
    elseif(b=="xx")then
    duv=(roui1+roui_1-2*c%rou)/dx**2
    elseif(b=="xy")then
    duv=(Ucn(c,roui1,rouj1,c%rou)+Ucn(c,roui_1,rouj_1,c%rou)-&
    &Ucn(c,roui_1,rouj1,c%rou)-Ucn(c,roui1,rouj_1,c%rou))/(4*dx*dy)
    elseif(b=="yy")then
    duv=(rouj1+rouj_1-2*c%rou)/dy**2
    end if
    
  end if   
return
end function duv
!针对用到四个角的情况ucorner
function Ucn(c,u1,u2,u0)
implicit none
type(gridtyp),pointer::c
real(pre)::u1,u2,Ucn,u0
 !Ucn=(d*u1+d*u2+1.414*d*u0)/(2*d+1.414*d)
  Ucn=(u1+u2+1.414*u0)/3.414
return
end function Ucn




!#######################
!nndx,水平的NND格式
function nndx(c,n)
implicit none
type(gridtyp),pointer::c
integer n
real(pre)::nndx,dx,H1,H2
real(pre)::rouip,uip,vip,TTip,pip,rouipp,uipp,vipp,TTipp,pipp,&
              &rouim,uim,vim,TTim,pim,rouimm,uimm,vimm,TTimm,pimm
real(pre)::Epi,Epi1,Epi2,Epi_1,Epi_2,Emi,Emi1,Emi2,Emi_1,Emi_2
              
dx=h/(2**c%lvlx)
!引用函数，对i-1，i+2等使用悬挂网格求解
call Tip(c,rouip,uip,vip,TTip,pip)     !注意Tip变量和函数名冲突
call Tipp(c,rouipp,uipp,vipp,TTipp,pipp)
call Tim(c,rouim,uim,vim,TTim,pim)
call Timm(c,rouimm,uimm,vimm,TTimm,pimm)

call Eplus(c%rou,c%u,c%v,c%T,c%p,Epi,n) 
call Eplus(rouip,uip,vip,TTip,pip,Epi1,n)
call Eplus(rouipp,uipp,vipp,TTipp,pipp,Epi2,n)
call Eplus(rouim,uim,vim,TTim,pim,Epi_1,n)
call Eplus(rouimm,uimm,vimm,TTimm,pimm,Epi_2,n)

call Emius(c%rou,c%u,c%v,c%T,c%p,Emi,n)
call Emius(rouip,uip,vip,TTip,pip,Emi1,n)
call Emius(rouipp,uipp,vipp,TTipp,pipp,Emi2,n)
call Emius(rouim,uim,vim,TTim,pim,Emi_1,n)
call Emius(rouimm,uimm,vimm,TTimm,pimm,Emi_2,n)

H1=Epi+Emi1+(minmod(Epi-Epi_1,Epi1-Epi)-minmod(Emi2-Emi1,Emi1-Emi))/2
H2=Epi_1+Emi+(minmod(Epi_1-Epi_2,Epi-Epi_1)-minmod(Emi1-Emi,Emi-Emi_1))/2
nndx=(H1-H2)/dx
return
end function nndx

!E+
subroutine Eplus(rou,u,v,T,p,Ep,n)
implicit none
real(pre)::rou,u,v,T,p,Ep,Et,a,Ma1 !,Ma2
integer n,nn

a=sqrt(abs(gama*p/rou))
Ma1=u/a
Et=rou*(R*T/(gama-1)+(u*u+v*v)/2)
!对于Ma>=1,E+ = E, E- = 0，Ma<=-1,反之
select case(n)
case(1)
  nn=1  
  if(abs(Ma1)<1)then
  Ep=rou*a*(1+Ma1)**2/4
  elseif(Ma1>=1)then
  Ep=rou*u
  else
  Ep=0
  endif
case(2)
  if(abs(Ma1)<1)then
  Ep=rou*a*(1+Ma1)**2/4*(2*a/gama)*((gama-1)*Ma1/2+1)
  elseif(Ma1>=1)then
  Ep=rou*u*u+p
  else
  Ep=0
  endif
case(3)
  if(abs(Ma1)<1)then
  Ep=rou*a*(1+Ma1)**2/4*v
  elseif(Ma1>=1)then
  Ep=rou*u*v
  else
  Ep=0
  endif
case(4)
  if(abs(Ma1)<1)then
  Ep=rou*a*(1+Ma1)**2/4*(2*a**2*((gama-1)*Ma1/2+1)**2/(gama**2-1)+v**2/2)
  elseif(Ma1>=1)then
  Ep=(Et+p)*u
  else
  Ep=0
  endif
end select

return
end subroutine Eplus

!E-
subroutine Emius(rou,u,v,T,p,Em,n)
implicit none
real(pre)::rou,u,v,T,p,Em,Et,a,Ma1 !,Ma2
integer n,nn

a=sqrt(abs(gama*p/rou))
Ma1=u/a
Et=rou*(R*T/(gama-1)+(u*u+v*v)/2)
select case(n)
case(1)
  if(abs(Ma1)<1)then
  Em=(-1)*rou*a*(1-Ma1)**2/4
  elseif(Ma1>=1)then
   Em=0
  else
  Em=rou*u
  end if
case(2)
  if(abs(Ma1)<1)then
  Em=(-1)*rou*a*(1-Ma1)**2/4*(2*a/gama)*((gama-1)*Ma1/2-1)
  elseif(Ma1>=1)then
  Em=0
  else
  Em=rou*u*u+p
  endif
case(3) 
  if(abs(Ma1)<1)then
  Em=(-1)*rou*a*(1-Ma1)**2/4*v
  elseif(Ma1>=1)then
  Em=0
  else
  Em=rou*u*v
  end if
case(4)
  if(abs(Ma1)<1)then
  Em=(-1)*rou*a*(1-Ma1)**2/4*(2*a**2*((gama-1)*Ma1/2-1)**2/(gama**2-1)+v**2/2)
  elseif(Ma1>=1)then
  Em=0
  else
  Em=(Et+p)*u
  endif
end select

return
end subroutine Emius

!nndy,竖直的NND格式
function nndy(c,n)
implicit none
type(gridtyp),pointer::c
integer n
real(pre)::nndy,dy,H1,H2
real(pre)::roujp,ujp,vjp,TTjp,pjp,roujpp,ujpp,vjpp,TTjpp,pjpp,&
              &roujm,ujm,vjm,TTjm,pjm,roujmm,ujmm,vjmm,TTjmm,pjmm
real(pre)::Fpj,Fpj1,Fpj2,Fpj_1,Fpj_2,Fmj,Fmj1,Fmj2,Fmj_1,Fmj_2
              
dy=h/(2**c%lvly)
!引用函数，对j-1，j+2等使用悬挂网格求解
call Tjp(c,roujp,ujp,vjp,TTjp,pjp)
call Tjpp(c,roujpp,ujpp,vjpp,TTjpp,pjpp)
call Tjm(c,roujm,ujm,vjm,TTjm,pjm)
call Tjmm(c,roujmm,ujmm,vjmm,TTjmm,pjmm)

call Fplus(c%rou,c%u,c%v,c%T,c%p,Fpj,n) 
call Fplus(roujp,ujp,vjp,TTjp,pjp,Fpj1,n)
call Fplus(roujpp,ujpp,vjpp,TTjpp,pjpp,Fpj2,n)
call Fplus(roujm,ujm,vjm,TTjm,pjm,Fpj_1,n)
call Fplus(roujmm,ujmm,vjmm,TTjmm,pjmm,Fpj_2,n)

call Fmius(c%rou,c%u,c%v,c%T,c%p,Fmj,n)
call Fmius(roujp,ujp,vjp,TTjp,pjp,Fmj1,n)
call Fmius(roujpp,ujpp,vjpp,TTjpp,pjpp,Fmj2,n)
call Fmius(roujm,ujm,vjm,TTjm,pjm,Fmj_1,n)
call Fmius(roujmm,ujmm,vjmm,TTjmm,pjmm,Fmj_2,n)


H1=Fpj+Fmj1+(minmod(Fpj-Fpj_1,Fpj1-Fpj)-minmod(Fmj2-Fmj1,Fmj1-Fmj))/2
H2=Fpj_1+Fmj+(minmod(Fpj_1-Fpj_2,Fpj-Fpj_1)-minmod(Fmj1-Fmj,Fmj-Fmj_1))/2
nndy=(H1-H2)/dy

return
end function nndy 

!F+
subroutine Fplus(rou,u,v,T,p,Fp,n)
implicit none
real(pre)::rou,u,v,T,p,Fp,Ft,a,Ma2 !,Ma1
integer n

a=sqrt(abs(gama*p/rou))
Ma2=v/a
Ft=rou*(R*T/(gama-1)+(u*u+v*v)/2)
select case(n)
case(1)
  if(abs(Ma2)<1)then
  Fp=rou*a*(1+Ma2)**2/4
  elseif(Ma2>=1)then
  Fp=rou*v
  else
  Fp=0
  endif
case(2)
  if(abs(Ma2)<1)then
  Fp=rou*a*(1+Ma2)**2/4*u
  elseif(Ma2>=1)then
  Fp=rou*u*v
  else
  Fp=0
  end if
case(3)
  if(abs(Ma2)<1)then
  Fp=rou*a*(1+Ma2)**2/4*(2*a/gama)*((gama-1)*Ma2/2+1)
  elseif(Ma2>=1)then
  Fp=rou*v*v+p
  else
  Fp=0
  end if
case(4)
  if(abs(Ma2)<1)then
  Fp=rou*a*(1+Ma2)**2/4*(2*a**2*((gama-1)*Ma2/2+1)**2/(gama**2-1)+u**2/2)
  elseif(Ma2>=1)then
  Fp=(Ft+p)*v
  else
  Fp=0
  endif
end select

return
end subroutine Fplus
!F-
subroutine Fmius(rou,u,v,T,p,Fm,n)
implicit none
real(pre)::rou,u,v,T,p,Fm,Ft,a,Ma2 !,Ma1
integer n
a=sqrt(abs(gama*p/rou))
Ma2=v/a
Ft=rou*(R*T/(gama-1)+(u*u+v*v)/2)
select case(n)
case(1)
  if(abs(Ma2)<1)then
  Fm=(-1)*rou*a*(1-Ma2)**2/4
  elseif(Ma2>=1)then
  Fm=0
  else
  Fm=rou*v
  endif
case(2)
  if(abs(Ma2)<1)then
  Fm=(-1)*rou*a*(1-Ma2)**2/4*u
  elseif(Ma2>=1)then
  Fm=0
  else
  Fm=rou*u*v
  endif
case(3)
  if(abs(Ma2)<1)then
  Fm=(-1)*rou*a*(1-Ma2)**2/4*(2*a/gama)*((gama-1)*Ma2/2-1)
  elseif(Ma2>=1)then
  Fm=0
  else
  Fm=rou*v*v+p
  endif
case(4)
  if(abs(Ma2)<1)then
  Fm=(-1)*rou*a*(1-Ma2)**2/4*(2*a**2*((gama-1)*Ma2/2-1)**2/(gama**2-1)+u**2/2)
  elseif(Ma2>=1)then
  Fm=0
  else
  Fm=(Ft+p)*v
  end if 
end select

return
end subroutine Fmius
!minmod
function minmod(x,y)
implicit none
real(pre)::x,y,minmod
real a
a=sign(1.0,x)+sign(1.0,y)  !注意函数参数和返回值类型
minmod=a*min(abs(x),abs(y))/2
return
end function minmod


end module
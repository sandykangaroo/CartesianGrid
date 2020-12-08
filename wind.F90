!一阶迎风格式，wind
module windd
use nndcenter
use mesh
use inflow
use precisionMod
use typeMod
implicit none

contains
subroutine wind(c)
implicit none
integer i,nn
type(gridtyp),pointer::c
real(pre)::U11,U21,U31,U41,U10,U20,U30,U40
real(pre)::fu1,fu2,fu3,fu4!,tstep

!tstep=h/(2**c%lvl)*Courant
!tstep=4D-6
U10=c%rou
U20=c%rou*c%u
U30=c%rou*c%v
U40=c%rou*(R*c%T/(gama-1)+((c%u)**2+(c%v)**2)/2)
!i=1第一步更新,=2第二步,即为RK方法
!do i=1,2

 call Res_wind_center(c,fu1,fu2,fu3,fu4)    !残差的差分格式，返回四个分量各自差分结果
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
end subroutine wind
!残差wind+center
subroutine Res_wind_center(c,fu1,fu2,fu3,fu4)
  implicit none
  type(gridtyp),pointer::c
  real(pre)::fu1,fu2,fu3,fu4,U11,a,b,wx1,wx2,wx3,wx4,wy1,wy2,wy3,wy4
  real(pre)::dEv1,dEv2,dEv3,dEv4,dFv1,dFv2,dFv3,dFv4,nn
  
  !call centerx(c,dEv1,dEv2,dEv3,dEv4)   !水平方向，输入c，返回四个差分后的量,
  !call centery(c,dFv1,dFv2,dFv3,dFv4)   !v代表miou，粘性
  !subroutine形式
  call windx(c,wx1,1)
  call windx(c,wx2,2)
  call windx(c,wx3,3)
  call windx(c,wx4,4)
  
  call windy(c,wy1,1)
  call windy(c,wy2,2)
  call windy(c,wy3,3)
  call windy(c,wy4,4)
  
  fu1=(-1)*wx1+(-1)*wy1 !+dEv1+dFv1
  fu2=(-1)*wx2+(-1)*wy2 !+dEv2+dFv2
  fu3=(-1)*wx3+(-1)*wy3 !+dEv3+dFv3
  fu4=(-1)*wx4+(-1)*wy4 !+dEv4+dFv4
  !function形式
  !fu1=(-1)*windx(c,1)+(-1)*windy(c,1)!+dEv1+dFv1
  !fu2=(-1)*windx(c,2)+(-1)*windy(c,2)!+dEv2+dFv2
  !fu3=(-1)*windx(c,3)+(-1)*windy(c,3)!+dEv3+dFv3
  !fu4=(-1)*windx(c,4)+(-1)*windy(c,4)!+dEv4+dFv4
  !只看粘性项
  !fu1=dEv1+dFv1 !+(-1)*windx(c,1)+(-1)*windy(c,1)
  !fu2=dEv2+dFv2 !+(-1)*windx(c,2)+(-1)*windy(c,2)
  !fu3=dEv3+dFv3 !+(-1)*windx(c,3)+(-1)*windy(c,3)
  !fu4=dEv4+dFv4 !+(-1)*windx(c,4)+(-1)*windy(c,4)
  
  !监视
  !U11=c%rou+tstep*fu1
  return
end subroutine Res_wind_center
!中心差分的直接引用，
!下面写对流向的wind格式
!水平windx
subroutine windx(c,wx,n)  !w是返回值
implicit none
type(gridtyp),pointer::c
integer n,nn
real(pre)::wx,dx,H1,H2
real(pre)::rouip,uip,vip,TTip,pip,&
              &rouim,uim,vim,TTim,pim
real(pre)::Epi,Epi1,Epi_1,Emi,Emi1,Emi_1            
dx=h/(2**c%lvlx)

call Tip(c,rouip,uip,vip,TTip,pip)
call Tim(c,rouim,uim,vim,TTim,pim)
!write(*,*)"c%rou,u,v,t,p",c%rou,c%u,c%v,c%t,c%p

call Eplus(c%rou,c%u,c%v,c%T,c%p,Epi,n) 
call Eplus(rouim,uim,vim,TTim,pim,Epi_1,n)

call Emius(c%rou,c%u,c%v,c%T,c%p,Emi,n)
call Emius(rouip,uip,vip,TTip,pip,Emi1,n)

!E+后差，E-前差
H1=Epi-Epi_1
H2=Emi1-Emi
!write(*,*)"Epi,Epi_1",Epi,Epi_1
!write(*,*)"Emi1-Emi",Emi1,Emi
!write(*,*)"H1,H2",H1,H2

!windx=(H1+H2)/dx
wx=(H1+H2)/dx
return
end subroutine windx
!windy,竖直的WIND格式
subroutine windy(c,wy,n)
implicit none
type(gridtyp),pointer::c
integer n
real(pre)::wy,dy,H1,H2
real(pre)::roujp,ujp,vjp,TTjp,pjp,&
              &roujm,ujm,vjm,TTjm,pjm
real(pre)::Fpj,Fpj1,Fpj_1,Fmj,Fmj1,Fmj_1
              
dy=h/(2**c%lvly)
!引用函数，对j-1，j+2等使用悬挂网格求解
call Tjp(c,roujp,ujp,vjp,TTjp,pjp)
call Tjm(c,roujm,ujm,vjm,TTjm,pjm)

call Fplus(c%rou,c%u,c%v,c%T,c%p,Fpj,n) 
call Fplus(roujm,ujm,vjm,TTjm,pjm,Fpj_1,n)

call Fmius(c%rou,c%u,c%v,c%T,c%p,Fmj,n)
call Fmius(roujp,ujp,vjp,TTjp,pjp,Fmj1,n)

H1=Fpj-Fpj_1
H2=Fmj1-Fmj 
!windy=(H1+H2)/dy
wy=(H1+H2)/dy
return
end subroutine windy

end module
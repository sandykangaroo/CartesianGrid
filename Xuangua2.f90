!悬挂网格处理方式
module frame
  use neighbor 
  use inflow
use precisionMod
use typeMod
  implicit none

 contains
!插值函数1，由四个子结点的信息得到父节点的信息，针对的是变量u 
recursive subroutine interposeA(c,rou,u,v,T,p)  
  implicit none
  type(gridtyp),pointer::c
  real(pre)::rou,u,v,T,p
  real(pre)::e1,e2,e3,e4,d1,d2,d3,d4,eTotal
  
  if(.not.associated(c%son1))then
       rou=c%rou
       u=c%u  
       v=c%v 
       T=c%T
       p=R*rou*T  
    return
  end if
  if(associated(c%son1%son1))then
    call interposeA(c%son1,c%son1%rou,c%son1%u,c%son1%v,c%son1%T,c%son1%p)
  end if
  if(associated(c%son2%son1))then
    call interposeA(c%son2,c%son2%rou,c%son2%u,c%son2%v,c%son2%T,c%son2%p)
  end if
  if(associated(c%son3).and.associated(c%son3%son1))then
    call interposeA(c%son3,c%son3%rou,c%son3%u,c%son3%v,c%son3%T,c%son3%p)
  end if
  if(associated(c%son3).and.associated(c%son4%son1))then
    call interposeA(c%son4,c%son4%rou,c%son4%u,c%son4%v,c%son4%T,c%son4%p)
  end if
  !插值公式
  if(associated(c%son1).and.(associated(c%son3)))then
d1=sqrt((c%center%x-c%son1%center%x)**2+(c%center%y-c%son1%center%y)**2)
d2=sqrt((c%center%x-c%son2%center%x)**2+(c%center%y-c%son2%center%y)**2)
d3=sqrt((c%center%x-c%son3%center%x)**2+(c%center%y-c%son3%center%y)**2)
d4=sqrt((c%center%x-c%son4%center%x)**2+(c%center%y-c%son4%center%y)**2)
e1=exp(-d1)
e2=exp(-d2)
e3=exp(-d3)
e4=exp(-d4)
eTotal=e1+e2+e3+e4
rou=(e1*c%son1%rou+e2*c%son2%rou+e3*c%son3%rou+e4*c%son4%rou)/eTotal
u=(e1*c%son1%u+e2*c%son2%u+e3*c%son3%u+e4*c%son4%u)/eTotal
v=(e1*c%son1%v+e2*c%son2%v+e3*c%son3%v+e4*c%son4%v)/eTotal
T=(e1*c%son1%T+e2*c%son2%T+e3*c%son3%T+e4*c%son4%T)/eTotal
p=R*rou*T    
elseif(associated(c%son1).and.(.not.associated(c%son3)))then
rou=(c%son1%rou+c%son2%rou)/2
u=(c%son1%u+c%son2%u)/2
v=(c%son1%v+c%son2%v)/2
T=(c%son1%T+c%son2%T)/2
p=R*rou*T   
end if
return
end subroutine interposeA

!插值函数2，由三个临近单元得到虚构单元的单元值,c2决定系数
recursive subroutine interposeB(c1,c2,c3,c4,n,m,rou,u,v,T,p)  !n控制序号-1 -2 1 2，m控制水平1竖直2
  implicit none
  integer n,m
  type(gridtyp),pointer::c1,c2,c3,c4,c0
  real(pre)::rou,u,v,T,p
  real(pre)::e1,e2,e3,e4,d1,d2,d3,d4,xx,yy,eTotal
  !通过单元c对虚构的单元中心点位置赋值
  call interposeA(c1,c1%rou,c1%u,c1%v,c1%T,c1%p)
  call interposeA(c2,c2%rou,c2%u,c2%v,c2%T,c2%p)
  call interposeA(c3,c3%rou,c3%u,c3%v,c3%T,c3%p)
  call interposeA(c4,c4%rou,c4%u,c4%v,c4%T,c4%p)
  
  if(m==1)then
      if(c1%spl==1)then
          xx=c1%center%x+n*h/2**(c1%lvlx)
          yy=c1%center%y
      else
          xx=c1%center%x+n*h/2**(c1%lvlx)
          yy=c1%center%y
      end if    
  else 
      if(c1%spl==2)then
          xx=c1%center%x
          yy=c1%center%y+n*h/2**(c1%lvly)
      else
          xx=c1%center%x
          yy=c1%center%y+n*h/2**(c1%lvly)
      end if       
  end if  
  d1=sqrt((xx-c1%center%x)**2+(yy-c1%center%y)**2)
  d2=sqrt((xx-c2%center%x)**2+(yy-c2%center%y)**2)
  d3=sqrt((xx-c3%center%x)**2+(yy-c3%center%y)**2)
  d4=sqrt((xx-c4%center%x)**2+(yy-c4%center%y)**2)
  
  e1=exp(-d1)
  e2=exp(-d2)
  e3=exp(-d3) 
  e4=exp(-d4)
  eTotal=e1+e2+e3+e4
  rou=(e1*c1%rou+e2*c2%rou+e3*c3%rou+e4*c4%rou)/eTotal
  u=(e1*c1%u+e2*c2%u+e3*c3%u+e4*c4%u)/eTotal
  v=(e1*c1%v+e2*c2%v+e3*c3%v+e4*c4%v)/eTotal
  T=(e1*c1%T+e2*c2%T+e3*c3%T+e4*c4%T)/eTotal
  p=R*rou*T
return  
end subroutine interposeB
 !T(i+1)为Tip，T(i+2)为Tipp，T(i-1)为Tim
recursive subroutine Tip(c,rou,u,v,T,p)
   implicit none
   integer::nn
   real(pre)::rou,u,v,T,p,roua,ua,va,Ta,pa,roub,ub,vb,Tb,pb
   type(gridtyp),pointer::c,cn_r,ct1,ct2,ct3,ccc,c1,pc,pc_r_nb,c_temp  !cn_r为c的右邻居,ct为右邻居的上邻居或者下邻居
   cn_r=>right_neighbor(c) 
   pc=>c%father
   if(associated(pc))then
        pc_r_nb=>right_neighbor(pc)
   else
        call interposeA(cn_r,rou,u,v,T,p)
        return
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !if(c%lvly-c%lvlx>1)then
 !    call Tip(c%father,rou,u,v,T,p)
 !    
 !    return
 !else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
if(c%spl==0)then   !分3+1种情况
    if(c%sp==1)then
        call interposeA(pc%son2,rou,u,v,T,p)
    elseif(c%sp==4)then
         call interposeA(pc%son3,rou,u,v,T,p)       
    elseif(c%sp==2)then
        !if(pc_r_nb%spl==0.and.(associated(pc_r_nb%son1)).and.(associated(pc_r_nb%son3)))then
        if(associated(pc_r_nb%son3))then
            call interposeA(pc_r_nb%son1,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb,ct1,ct2,1,1,rou,u,v,T,p)  
        endif
    elseif(c%sp==3)then
        !if(pc_r_nb%spl==0.and.(associated(pc_r_nb%son1)).and.(associated(pc_r_nb%son3)))then
        if(associated(pc_r_nb%son3))then
            call interposeA(pc_r_nb%son4,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb,ct1,ct2,1,1,rou,u,v,T,p) 
        endif
    else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb,ct1,ct2,1,1,rou,u,v,T,p)     
    endif
elseif(c%spl==1)then
    if(c%sp==1)then
        !if(pc_r_nb%spl==0.and.(associated(pc_r_nb%son1)).and.(associated(pc_r_nb%son3)))then
        if(associated(pc_r_nb%son3))then
            !call interposeA(pc_r_nb%son1,rou,u,v,T,p)
            !call interposeA(pc_r_nb%son2,rou,u,v,T,p)
            !rou=(pc_r_nb%son1%rou+pc_r_nb%son2%rou)/2
            !u=(pc_r_nb%son1%u+pc_r_nb%son2%u)/2
            !v=(pc_r_nb%son1%v+pc_r_nb%son2%v)/2
            !T=(pc_r_nb%son1%T+pc_r_nb%son2%T)/2
            !p=R*rou*T
            call interposeA(pc_r_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_r_nb%son2,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_r_nb%typ==1.and.associated(pc_r_nb%son1))then
            call interposeA(pc_r_nb%son1,rou,u,v,T,p)     
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb,ct1,ct2,1,1,rou,u,v,T,p) 
        endif
        
    else if(c%sp==2)then
        !if(pc_r_nb%spl==0.and.(associated(pc_r_nb%son1)).and.(associated(pc_r_nb%son3)))then
         if(associated(pc_r_nb%son3))then   
            call interposeA(pc_r_nb%son3,roua,ua,va,Ta,pa)
            call interposeA(pc_r_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_r_nb%typ==1.and.associated(pc_r_nb%son1))then
            call interposeA(pc_r_nb%son2,rou,u,v,T,p) 
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb,ct1,ct2,1,1,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==2)then
    if(c%sp==2)then
        !if(pc_r_nb%spl==0.and.(associated(pc_r_nb%son1)).and.(associated(pc_r_nb%son3)))then
        if(associated(pc_r_nb%son3))then
            call interposeA(pc_r_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_r_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_r_nb%typ==2.and.associated(pc_r_nb%son1))then
            call interposeA(pc_r_nb%son1,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb,ct1,ct2,1,1,rou,u,v,T,p)
        endif
        
    elseif(c%sp==1)then
        call interposeA(pc%son2,rou,u,v,T,p)
    endif
endif
return
!endif    
end subroutine Tip
!子函数 求T(i+2)的流场值
recursive subroutine Tipp(c,rou,u,v,T,p)      
  implicit none
  type(gridtyp),pointer::c,cn_r,cnn_r,ct1,ct2,ct3,ct4,pc,pc_r_nb,pc_r_nb_nb 
  real(pre)::rou,u,v,T,p,roua,ua,va,Ta,pa,roub,ub,vb,Tb,pb
  cn_r=>right_neighbor(c) 
  cnn_r=>right_neighbor(cn_r)
  pc=>c%father
  
   if(associated(pc))then
   pc_r_nb=>right_neighbor(pc)
   else
    ct1=>right_neighbor(c)
    ct2=>right_neighbor(ct1)
    call interposeA(ct2,rou,u,v,T,p)
    return
   endif
   
   if(associated(pc_r_nb))then
   pc_r_nb_nb=>right_neighbor(pc_r_nb) 
   endif
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !if(c%lvly-c%lvlx>1)then
 !    call Tipp(c%father,rou,u,v,T,p)
 !    
 !    return
 !else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  if(c%spl==0)then   !分3+1种情况
    if(c%sp==1)then
        if(associated(pc_r_nb%son3))then
            call interposeA(pc_r_nb%son1,rou,u,v,T,p)
        else
            ct1=>up_neighbor(pc%son2)!需不需要求son2的流场
            ct2=>down_neighbor(pc%son2)
            call interposeB(pc%son2,pc_r_nb,ct1,ct2,1,1,rou,u,v,T,p)
        endif    
    elseif(c%sp==4)then
        if(associated(pc_r_nb%son3))then
            call interposeA(pc_r_nb%son4,rou,u,v,T,p)
        else
            ct1=>up_neighbor(pc%son3)
            ct2=>down_neighbor(pc%son3)
            call interposeB(pc%son3,pc_r_nb,ct1,ct2,1,1,rou,u,v,T,p)
        endif
    elseif(c%sp==2)then
        if(associated(pc_r_nb%son3))then
            call interposeA(pc_r_nb%son2,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb,ct1,ct2,2,1,rou,u,v,T,p)
        endif
    elseif(c%sp==3)then
        if(associated(pc_r_nb%son3))then
            call interposeA(pc_r_nb%son3,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb,ct1,ct2,2,1,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==1)then
    if(c%sp==1)then
        if(pc_r_nb%spl==0)then
        if(associated(pc_r_nb_nb%son3))then
            call interposeA(pc_r_nb_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_r_nb_nb%son2,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_r_nb_nb%typ==1.and.associated(pc_r_nb_nb%son1))then
            call interposeA(pc_r_nb_nb%son1,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb_nb,ct1,ct2,2,1,rou,u,v,T,p)
        endif
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb,ct1,ct2,2,1,rou,u,v,T,p)
        endif
        
        
    elseif(c%sp==2)then
        if(pc_r_nb%spl==0)then
        if(associated(pc_r_nb_nb%son3))then
            call interposeA(pc_r_nb_nb%son3,roua,ua,va,Ta,pa)
            call interposeA(pc_r_nb_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_r_nb_nb%typ==1.and.associated(pc_r_nb_nb%son1))then
            call interposeA(pc_r_nb_nb%son2,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb_nb,ct1,ct2,2,1,rou,u,v,T,p)
        endif
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb,ct1,ct2,2,1,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==2)then
    if(c%sp==2)then
        if(associated(pc_r_nb%son3))then
            call interposeA(pc_r_nb%son2,roua,ua,va,Ta,pa)
            call interposeA(pc_r_nb%son3,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_r_nb%typ==1.and.associated(pc_r_nb%son1))then
            call interposeA(pc_r_nb%son2,rou,u,v,T,p) 
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_r_nb,ct1,ct2,2,1,rou,u,v,T,p) 
        endif
        
    elseif(c%sp==1)then
        if(associated(pc_r_nb%son3))then
            call interposeA(pc_r_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_r_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_r_nb%typ==1.and.associated(pc_r_nb%son1))then
            call interposeA(pc_r_nb%son1,rou,u,v,T,p)    
        else
            ct1=>up_neighbor(pc%son2)
            ct2=>down_neighbor(pc%son2)
            call interposeB(pc%son2,pc_r_nb,ct1,ct2,1,1,rou,u,v,T,p)
        endif
    endif
endif
  
return
!endif
end subroutine Tipp

!T(i-1)  Tim
recursive subroutine Tim(c,rou,u,v,T,p)
   implicit none
   real(pre)::rou,u,v,T,p,roua,ua,va,Ta,pa,roub,ub,vb,Tb,pb
   type(gridtyp),pointer::c,cn_r,ct1,ct2,ct3,ccc,c1,pc,pc_l_nb  
   cn_r=>left_neighbor(c) 
   pc=>c%father
   if(associated(pc))then
   pc_l_nb=>left_neighbor(pc)
   else
    call interposeA(cn_r,rou,u,v,T,p)
    return
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !if(c%lvly-c%lvlx>1)then
 !    call Tim(c%father,rou,u,v,T,p)
 !    
 !    return
 !else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(c%spl==0)then   !分3+1种情况
    if(c%sp==2)then
        call interposeA(pc%son1,rou,u,v,T,p)
    elseif(c%sp==3)then
        call interposeA(pc%son4,rou,u,v,T,p)
    elseif(c%sp==1)then
        if(associated(pc_l_nb%son3))then
            call interposeA(pc_l_nb%son2,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-1,1,rou,u,v,T,p)
        endif
    elseif(c%sp==4)then
        if(associated(pc_l_nb%son3))then
            call interposeA(pc_l_nb%son3,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-1,1,rou,u,v,T,p)
        endif
    else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-1,1,rou,u,v,T,p)
    endif
elseif(c%spl==1)then
    if(c%sp==1)then
        if(associated(pc_l_nb%son3))then
            call interposeA(pc_l_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_l_nb%son2,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_l_nb%typ==1.and.associated(pc_l_nb%son1))then
            call interposeA(pc_l_nb%son1,rou,u,v,T,p) 
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-1,1,rou,u,v,T,p)
        endif
        
    elseif(c%sp==2)then
        !if(pc_l_nb%spl==0.and.(associated(pc_l_nb%son1)).and.(associated(pc_l_nb%son3)))then
        if(associated(pc_l_nb%son3))then
            call interposeA(pc_l_nb%son3,roua,ua,va,Ta,pa)
            call interposeA(pc_l_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_l_nb%typ==1.and.associated(pc_l_nb%son1))then
            call interposeA(pc_l_nb%son2,rou,u,v,T,p) 
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-1,1,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==2)then
    if(c%sp==1)then
        !if(pc_l_nb%spl==0.and.(associated(pc_l_nb%son1)).and.(associated(pc_l_nb%son3)))then
        if(associated(pc_l_nb%son3))then
            call interposeA(pc_l_nb%son2,roua,ua,va,Ta,pa)
            call interposeA(pc_l_nb%son3,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_l_nb%typ==2.and.associated(pc_l_nb%son1))then
            call interposeA(pc_l_nb%son2,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-1,1,rou,u,v,T,p)
        endif    
    elseif(c%sp==2)then
        call interposeA(pc%son1,rou,u,v,T,p)
    endif
endif
   return 
!end if
end subroutine Tim
!子函数 求T(i-2)的流场值
recursive subroutine Timm(c,rou,u,v,T,p)      
  implicit none
  type(gridtyp),pointer::c,cn_r,cnn_r,ct1,ct2,ct3,ct4,pc,pc_l_nb,pc_l_nb_nb 
  real(pre)::rou,u,v,T,p,roua,ua,va,Ta,pa,roub,ub,vb,Tb,pb
  cn_r=>left_neighbor(c) 
  cnn_r=>left_neighbor(cn_r)
  pc=>c%father
  
  if(associated(pc))then
   pc_l_nb=>left_neighbor(pc)
   else
    ct1=>left_neighbor(c)
    ct2=>left_neighbor(ct1)
    call interposeA(ct2,rou,u,v,T,p)
    return
   endif
   if(associated(pc_l_nb))then
   pc_l_nb_nb=>left_neighbor(pc_l_nb)
   !else
   !    return
   endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !if(c%lvly-c%lvlx>1)then
 !    call Timm(c%father,rou,u,v,T,p)
 !    
 !    return
 !else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(c%spl==0)then   !分3+1种情况
    if(c%sp==2)then
        if(associated(pc_l_nb%son3))then
            call interposeA(pc_l_nb%son2,rou,u,v,T,p)
        else
            ct1=>up_neighbor(pc%son1)
            ct2=>down_neighbor(pc%son1)
            call interposeB(pc%son1,pc_l_nb,ct1,ct2,-1,1,rou,u,v,T,p)
        endif    
        
    elseif(c%sp==3)then
        if(associated(pc_l_nb%son3))then
            call interposeA(pc_l_nb%son3,rou,u,v,T,p)
        else
            ct1=>up_neighbor(pc%son4)
            ct2=>down_neighbor(pc%son4)
            call interposeB(pc%son4,pc_l_nb,ct1,ct2,-1,1,rou,u,v,T,p)
        endif
    elseif(c%sp==1)then
        if(associated(pc_l_nb%son3))then
            call interposeA(pc_l_nb%son1,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-2,1,rou,u,v,T,p)
        endif
    elseif(c%sp==4)then
        if(associated(pc_l_nb%son3))then
            call interposeA(pc_l_nb%son4,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-2,1,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==1)then
    if(c%sp==1)then
        if(pc_l_nb%spl==0)then
        if(associated(pc_l_nb_nb%son3))then
            call interposeA(pc_l_nb_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_l_nb_nb%son2,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_l_nb_nb%typ==1.and.associated(pc_l_nb_nb%son1))then
            call interposeA(pc_l_nb_nb%son1,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb_nb,ct1,ct2,-2,1,rou,u,v,T,p)
        endif
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-2,1,rou,u,v,T,p)
        endif
          
    elseif(c%sp==2)then
        if(pc_l_nb%spl==0)then
        if(associated(pc_l_nb_nb%son3))then
            call interposeA(pc_l_nb_nb%son3,roua,ua,va,Ta,pa)
            call interposeA(pc_l_nb_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_l_nb_nb%typ==1.and.associated(pc_l_nb_nb%son1))then
            call interposeA(pc_l_nb_nb%son2,rou,u,v,T,p)
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb_nb,ct1,ct2,-2,1,rou,u,v,T,p)
        endif
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-2,1,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==2)then
    if(c%sp==2)then
        if(associated(pc_l_nb%son3))then
            call interposeA(pc_l_nb%son2,roua,ua,va,Ta,pa)
            call interposeA(pc_l_nb%son3,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_l_nb%typ==2.and.associated(pc_l_nb%son1))then
            call interposeA(pc_l_nb%son2,rou,u,v,T,p) 
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-2,1,rou,u,v,T,p)
        endif
        
    elseif(c%sp==1)then
        if(associated(pc_l_nb%son3))then
            call interposeA(pc_l_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_l_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_l_nb%typ==2.and.associated(pc_l_nb%son1))then
            call interposeA(pc_l_nb%son1,rou,u,v,T,p)    
        else
            ct1=>up_neighbor(c)
            ct2=>down_neighbor(c)
            call interposeB(c,pc_l_nb,ct1,ct2,-2,1,rou,u,v,T,p)
        endif
    endif
endif
  
return
!endif
end subroutine Timm

!T(j+1) Tjp
recursive subroutine Tjp(c,rou,u,v,T,p)
   implicit none
   real(pre)::rou,u,v,T,p,roua,ua,va,Ta,pa,roub,ub,vb,Tb,pb
   type(gridtyp),pointer::c,cn_u,ct1,ct2,ct3,ccc,c1,pc,pc_u_nb  !cn_u为c的右邻居,ct为右邻居的上邻居或者下邻居
   cn_u=>up_neighbor(c) 
   pc=>c%father
   
   if(associated(pc))then
   pc_u_nb=>up_neighbor(pc)
   else
       call interposeA(cn_u,rou,u,v,T,p)
       return
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !if(c%lvlx-c%lvly>1)then
 !    call Tjp(c%father,rou,u,v,T,p)
 !    
 !    return
 !else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
if(c%spl==0)then   !分3+1种情况
    if(c%sp==1)then
        call interposeA(pc%son4,rou,u,v,T,p)
    elseif(c%sp==2)then
        call interposeA(pc%son3,rou,u,v,T,p)
    elseif(c%sp==3)then
        if(associated(pc_u_nb%son3))then    
            call interposeA(pc_u_nb%son2,rou,u,v,T,p)
        else
           ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb,ct1,ct2,1,2,rou,u,v,T,p) 
        endif
    elseif(c%sp==4)then 
        if(associated(pc_u_nb%son3))then
            call interposeA(pc_u_nb%son1,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb,ct1,ct2,1,2,rou,u,v,T,p)
        endif
    else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb,ct1,ct2,1,2,rou,u,v,T,p)
    endif
elseif(c%spl==2)then
    if(c%sp==1)then
        !if(pc_u_nb%spl==0.and.(associated(pc_u_nb%son1)).and.(associated(pc_u_nb%son3)))then
        if(associated(pc_u_nb%son3))then
            !call interposeA(pc_u_nb%son1,rou,u,v,T,p)
            !call interposeA(pc_u_nb%son4,rou,u,v,T,p)
            !rou=(pc_u_nb%son1%rou+pc_u_nb%son4%rou)/2
            !u=(pc_u_nb%son1%u+pc_u_nb%son4%u)/2
            !v=(pc_u_nb%son1%v+pc_u_nb%son4%v)/2
            !T=(pc_u_nb%son1%T+pc_u_nb%son4%T)/2
            !p=R*rou*T
            call interposeA(pc_u_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_u_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_u_nb%typ==2.and.associated(pc_u_nb%son1))then
            call interposeA(pc_u_nb%son1,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb,ct1,ct2,1,2,rou,u,v,T,p)
        endif
        
    elseif(c%sp==2)then
        !if(pc_u_nb%spl==0.and.(associated(pc_u_nb%son1)).and.(associated(pc_u_nb%son3)))then
        if(associated(pc_u_nb%son3))then
            !call interposeA(pc_u_nb%son3,rou,u,v,T,p)
            !call interposeA(pc_u_nb%son2,rou,u,v,T,p)
            !rou=(pc_u_nb%son3%rou+pc_u_nb%son2%rou)/2
            !u=(pc_u_nb%son3%u+pc_u_nb%son2%u)/2
            !v=(pc_u_nb%son3%v+pc_u_nb%son2%v)/2
            !T=(pc_u_nb%son3%T+pc_u_nb%son2%T)/2
            !p=R*rou*T
            call interposeA(pc_u_nb%son3,roua,ua,va,Ta,pa)
            call interposeA(pc_u_nb%son2,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_u_nb%typ==2.and.associated(pc_u_nb%son1))then
            call interposeA(pc_u_nb%son2,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb,ct1,ct2,1,2,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==1)then
    if(c%sp==2)then
        !if(pc_u_nb%spl==0.and.(associated(pc_u_nb%son1)).and.(associated(pc_u_nb%son3)))then
        if(associated(pc_u_nb%son3))then
            call interposeA(pc_u_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_u_nb%son2,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_u_nb%typ==1.and.associated(pc_u_nb%son1))then
            call interposeA(pc_u_nb%son1,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb,ct1,ct2,1,2,rou,u,v,T,p)
        endif
        
    elseif(c%sp==1)then
        call interposeA(pc%son2,rou,u,v,T,p)
    endif
endif

   return 
!endif   
end subroutine Tjp
!子函数 求T(j+2)的流场值
recursive subroutine Tjpp(c,rou,u,v,T,p)      
  implicit none
  type(gridtyp),pointer::c,cn_u,cnn_u,ct1,ct2,ct3,ct4,pc,pc_u_nb,pc_u_nb_nb 
  real(pre)::rou,u,v,T,p,roua,ua,va,Ta,pa,roub,ub,vb,Tb,pb
  cn_u=>up_neighbor(c) 
  cnn_u=>up_neighbor(cn_u)
  pc=>c%father
  
  
  if(associated(pc))then
   pc_u_nb=>up_neighbor(pc)
   else
    ct1=>up_neighbor(c)
    ct2=>up_neighbor(ct1)
    call interposeA(ct2,rou,u,v,T,p)
    return
   endif
   if(associated(pc_u_nb))then
   pc_u_nb_nb=>up_neighbor(pc_u_nb)
   !else
   !    return
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !if(c%lvlx-c%lvly>1)then
 !    call Tjpp(c%father,rou,u,v,T,p)
 !    
 !    return
 !else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  if(c%spl==0)then   !分3+1种情况
    if(c%sp==1)then
        if(associated(pc_u_nb%son3))then
            call interposeA(pc_u_nb%son1,rou,u,v,T,p)
        else
            ct1=>right_neighbor(pc%son4)
            ct2=>left_neighbor(pc%son4)
            call interposeB(pc%son4,pc_u_nb,ct1,ct2,1,2,rou,u,v,T,p)
        endif    
        
    elseif(c%sp==2)then
        if(associated(pc_u_nb%son3))then
            call interposeA(pc_u_nb%son2,rou,u,v,T,p)
        else
            ct1=>right_neighbor(pc%son3)
            ct2=>left_neighbor(pc%son3)
            call interposeB(pc%son3,pc_u_nb,ct1,ct2,1,2,rou,u,v,T,p)
        endif
    elseif(c%sp==3)then
        if(associated(pc_u_nb%son3))then
            call interposeA(pc_u_nb%son3,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb,ct1,ct2,2,2,rou,u,v,T,p)
        endif
    elseif(c%sp==4)then
        if(associated(pc_u_nb%son3))then
            call interposeA(pc_u_nb%son4,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb,ct1,ct2,2,2,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==2)then
    if(c%sp==1)then
        if(pc_u_nb%spl==0)then
        if(associated(pc_u_nb_nb%son3))then
            call interposeA(pc_u_nb_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_u_nb_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_u_nb_nb%typ==2.and.associated(pc_u_nb_nb%son1))then
            call interposeA(pc_u_nb_nb%son1,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb_nb,ct1,ct2,2,2,rou,u,v,T,p)
        endif
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb,ct1,ct2,2,2,rou,u,v,T,p)
        endif
        
        
    elseif(c%sp==2)then
        if(pc_u_nb%spl==0)then
        if(associated(pc_u_nb_nb%son3))then
            call interposeA(pc_u_nb_nb%son2,roua,ua,va,Ta,pa)
            call interposeA(pc_u_nb_nb%son3,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_u_nb_nb%typ==2.and.associated(pc_u_nb_nb%son1))then
            call interposeA(pc_u_nb_nb%son2,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb_nb,ct1,ct2,2,2,rou,u,v,T,p)
        endif
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb,ct1,ct2,2,2,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==1)then
    if(c%sp==2)then
        if(associated(pc_u_nb_nb%son3))then
            call interposeA(pc_u_nb_nb%son4,roua,ua,va,Ta,pa)
            call interposeA(pc_u_nb_nb%son3,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_u_nb_nb%typ==1.and.associated(pc_u_nb_nb%son1))then
            call interposeA(pc_u_nb_nb%son2,rou,u,v,T,p) 
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_u_nb,ct1,ct2,2,2,rou,u,v,T,p)
        endif
        
    elseif(c%sp==1)then
        if(associated(pc_u_nb_nb%son3))then
            call interposeA(pc_u_nb_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_u_nb_nb%son2,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_u_nb_nb%typ==1.and.associated(pc_u_nb_nb%son1))then
            call interposeA(pc_u_nb_nb%son1,rou,u,v,T,p)    
        else
            ct1=>right_neighbor(pc%son2)
            ct2=>left_neighbor(pc%son2)
            call interposeB(pc%son2,pc_u_nb,ct1,ct2,1,2,rou,u,v,T,p)
        endif
    endif
endif
  
return
!endif
end subroutine Tjpp
!T(j-1)
recursive subroutine Tjm(c,rou,u,v,T,p)
   implicit none
   real(pre)::rou,u,v,T,p,roua,ua,va,Ta,pa,roub,ub,vb,Tb,pb
   type(gridtyp),pointer::c,cn_u,ct1,ct2,ct3,ccc,c1,pc,pc_d_nb  !cn_u为c的右邻居,ct为右邻居的上邻居或者下邻居
   cn_u=>down_neighbor(c) 
   pc=>c%father
   
   if(associated(pc))then
   pc_d_nb=>down_neighbor(pc)
   else
       call interposeA(cn_u,rou,u,v,T,p)
       return
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !if(c%lvlx-c%lvly>1)then
 !    call Tjm(c%father,rou,u,v,T,p)
 !    
 !    return
 !else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(c%spl==0)then   !分3+1种情况
    if(c%sp==4)then
        call interposeA(pc%son1,rou,u,v,T,p)
    elseif(c%sp==3)then
        call interposeA(pc%son2,rou,u,v,T,p)
    elseif(c%sp==2)then
        if(associated(pc_d_nb%son3))then
            call interposeA(pc_d_nb%son3,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb,ct1,ct2,-1,2,rou,u,v,T,p)
        endif
    elseif(c%sp==1)then
        if(associated(pc_d_nb%son3))then
            call interposeA(pc_d_nb%son4,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb,ct1,ct2,-1,2,rou,u,v,T,p)
        endif
    else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb,ct1,ct2,-1,2,rou,u,v,T,p)
    endif
elseif(c%spl==2)then
    if(c%sp==1)then
        !if(pc_d_nb%spl==0.and.(associated(pc_d_nb%son1)).and.(associated(pc_d_nb%son3)))then
        if(associated(pc_d_nb%son3))then
            call interposeA(pc_d_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_d_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_d_nb%typ==2.and.associated(pc_d_nb%son1))then
            call interposeA(pc_d_nb%son1,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb,ct1,ct2,-1,2,rou,u,v,T,p)
        endif
        
    elseif(c%sp==2)then
        !if(pc_d_nb%spl==0.and.(associated(pc_d_nb%son1)).and.(associated(pc_d_nb%son3)))then
        if(associated(pc_d_nb%son3))then
            call interposeA(pc_d_nb%son3,roua,ua,va,Ta,pa)
            call interposeA(pc_d_nb%son2,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_d_nb%typ==2.and.associated(pc_d_nb%son1))then
            call interposeA(pc_d_nb%son2,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb,ct1,ct2,-1,2,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==1)then
    if(c%sp==1)then
        !if(pc_d_nb%spl==0.and.(associated(pc_d_nb%son1)).and.(associated(pc_d_nb%son3)))then
        if(associated(pc_d_nb%son3))then
            call interposeA(pc_d_nb%son3,roua,ua,va,Ta,pa)
            call interposeA(pc_d_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_d_nb%typ==1.and.associated(pc_d_nb%son1))then
            call interposeA(pc_d_nb%son2,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb,ct1,ct2,-1,2,rou,u,v,T,p)  
        endif
        
    elseif(c%sp==2)then
        call interposeA(pc%son1,rou,u,v,T,p)
    endif
endif

   return 
!endif
end subroutine Tjm
!子函数 求T(j-2)的流场值
recursive subroutine Tjmm(c,rou,u,v,T,p)      
  implicit none
  type(gridtyp),pointer::c,cn_u,cnn_u,ct1,ct2,ct3,ct4,pc,pc_d_nb,pc_d_nb_nb 
  real(pre)::rou,u,v,T,p,roua,ua,va,Ta,pa,roub,ub,vb,Tb,pb
  cn_u=>down_neighbor(c) 
  cnn_u=>down_neighbor(cn_u)
  pc=>c%father
  
  
  if(associated(pc))then
   pc_d_nb=>down_neighbor(pc)
   else
    ct1=>down_neighbor(c)
    ct2=>down_neighbor(ct1)
    call interposeA(ct2,rou,u,v,T,p)
    return
   endif
   if(associated( pc_d_nb))then
    pc_d_nb_nb=>down_neighbor(pc_d_nb)
   !else
   !    return
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !if(c%lvlx-c%lvly>1)then
 !    call Tjmm(c%father,rou,u,v,T,p)
 !    
 !    return
 !else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  if(c%spl==0)then   !分3+1种情况
    if(c%sp==4)then
        if(associated(pc_d_nb%son3))then
            call interposeA(pc_d_nb%son4,rou,u,v,T,p)
        else
            ct1=>right_neighbor(pc%son1)
            ct2=>left_neighbor(pc%son1)
            call interposeB(pc%son1,pc_d_nb,ct1,ct2,-1,2,rou,u,v,T,p)
        endif    
        
    elseif(c%sp==3)then
        if(associated(pc_d_nb%son3))then
            call interposeA(pc_d_nb%son3,rou,u,v,T,p)
        else
            ct1=>right_neighbor(pc%son2)
            ct2=>left_neighbor(pc%son2)
            call interposeB(pc%son2,pc_d_nb,ct1,ct2,-1,2,rou,u,v,T,p)
        endif
    elseif(c%sp==2)then
        if(associated(pc_d_nb%son3))then
            call interposeA(pc_d_nb%son2,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb,ct1,ct2,-2,2,rou,u,v,T,p)
        endif
    elseif(c%sp==1)then
        if(associated(pc_d_nb%son3))then
            call interposeA(pc_d_nb%son1,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb,ct1,ct2,-2,2,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==2)then
    if(c%sp==1)then
        if(pc_d_nb%spl==0)then
        if(associated(pc_d_nb_nb%son3))then
            call interposeA(pc_d_nb_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_d_nb_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_d_nb_nb%typ==2.and.associated(pc_d_nb_nb%son1))then
            call interposeA(pc_d_nb_nb%son1,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb_nb,ct1,ct2,-2,2,rou,u,v,T,p)
        endif
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb,ct1,ct2,-2,2,rou,u,v,T,p)
        endif
        
        
    elseif(c%sp==2)then
        if(pc_d_nb%spl==0)then
        if(associated(pc_d_nb_nb%son3))then
            call interposeA(pc_d_nb_nb%son2,roua,ua,va,Ta,pa)
            call interposeA(pc_d_nb_nb%son3,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_d_nb_nb%typ==2.and.associated(pc_d_nb_nb%son1))then
            call interposeA(pc_d_nb_nb%son2,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb_nb,ct1,ct2,-2,2,rou,u,v,T,p)
        endif
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb,ct1,ct2,-2,2,rou,u,v,T,p)
        endif
    endif
elseif(c%spl==1)then
    if(c%sp==1)then
        if(associated(pc_d_nb%son3))then
            call interposeA(pc_d_nb%son1,roua,ua,va,Ta,pa)
            call interposeA(pc_d_nb%son2,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_d_nb%typ==1.and.associated(pc_d_nb%son1))then
            call interposeA(pc_d_nb%son1,rou,u,v,T,p)
        else
            ct1=>right_neighbor(c)
            ct2=>left_neighbor(c)
            call interposeB(c,pc_d_nb,ct1,ct2,-2,2,rou,u,v,T,p) 
        endif
        
    elseif(c%sp==2)then
        if(associated(pc_d_nb%son3))then
            call interposeA(pc_d_nb%son3,roua,ua,va,Ta,pa)
            call interposeA(pc_d_nb%son4,roub,ub,vb,Tb,pb)
            rou=(roua+roub)/2
            u=(ua+ub)/2
            v=(va+vb)/2
            T=(Ta+Tb)/2
            p=R*rou*T
        elseif(pc_d_nb%typ==1.and.associated(pc_d_nb%son1))then
            call interposeA(pc_d_nb%son2,rou,u,v,T,p)    
        else
            ct1=>right_neighbor(pc%son1)
            ct2=>left_neighbor(pc%son1)
            call interposeB(pc%son1,pc_d_nb,ct1,ct2,-1,2,rou,u,v,T,p)  
        endif
    endif
endif
  
return
!endif
end subroutine Tjmm
!！对于求对称单元，判断插值单元或者
function coe(c0,c)
implicit none
integer coe
type(gridtyp),pointer::c,c0

if(c%sort==c0%sort)then
 coe=1
 else
 coe=0
end if 

return
end function coe

end module frame
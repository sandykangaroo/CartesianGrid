Module solution
!use mesh
use precisionMod
use inflowGlobal
use bc
use nndcenter
use windd
use outfile
use flow_result
use solution_aniso
use NRRSolveMod
use NRRMeshMod
use typeMod
implicit none
integer::nLUmax=0,nLVmax=0,nLRoumax=0,nLRoumin=0!!number for limiter has used
integer::nLTmax=0,nLTmin=0,nLPmax=0,nLPmin=0
integer::Residcount
real(pre)::Residu,Residv,Residrou,Residt,Residp
integer::maxRCu,maxRCv,maxRCrou,maxRCT,maxRCp
real(pre)::maxRu,maxRv,maxRrou,maxRT,maxRp
    contains

subroutine solute()
    implicit none
    integer i,j,it,nn
    type(gridtyp),pointer::c
    character(5)::char
    !首先初始化流场，先对于外部单元赋值，再对物面虚拟单元求解
    !外部单元初始化流场

    !!虚拟单元（物面边界赋值）
    do i=1,total
        c=>cell(i)
        call body_boundary(c)
    end do
    !!暂存值替代流场量，更新
    call update
    !
    write(*,*),"________Run calculation________"
    !!初始化完成之后，对每个单元递归，空间推进和时间推进更新信息
    do it=1,itmax
        !!Output calculation statue
        if(mod(it,20)==1) print'(//7(8XA5))',   &
            "Iter ","Limit","u    ","v    ","rou  ","T    ","p    "
        !!Initial
        Residcount=0
        Residu=0.0;Residv=0.0;Residrou=0.0;Residt=0.0;Residp=0.0;
        maxRu=0;   maxRv=0;   maxRrou=0;   maxRT=0;   maxRp=0
        !!Solute cell
        do i=m+1,total
            c=>cell(i)
            call solution_cell(c)
        end do
        do i=1,m
            c=>cell(i)
            call solution_cell(c)
        end do
        call update
        !!GhostCell
        do i=1,total
            c=>cell(i)
            call body_boundary(c)
        enddo
        !!Renew data
        call update
        !!Output calculation statue
        if(Limiter==1) call CheckLimiter
        print'(7XI6,7XI6,5(2XE11.4))', &
            it,OutofLimit,Residu,Residv,Residrou,Residt,Residp
        if(check==1)then
            print'(5XA17,4X,5(3XI10))',&
            'MaxResidual Cell:',maxRCu,maxRCv,maxRCrou,maxRCT,maxRCp
            print'(5XA18,3X,5(3XF10.4))',&
            'MaxResidual Value:',maxRu,maxRv,maxRrou,maxRT,maxRp
        endif
        !!AutoRefine
        if(AutoRefine==1)then
            if(it>=AutoRefineStartStep)then
                if(mod(it,AutoRefineStep)==1)then
                    if(Anisotropic==0)then
                            call gra_secgra
                            call solute_ref
                            call solute_gridmodify
                    else
                            call gra_secgra_aniso
                            call solute_ref_aniso
                            call solute_gridmodify
                    endif
                    if(NRR==1)then
                            do i=1,total
                                c=>cell(i)
                                call NRRregionMark(c)
                            enddo
                    endif
                    call cellinf
                    print*,'Mesh refinement to count: ', CellCount
                endif
            endif
        endif
        !!AutoSave
        if(SaveFile==1)then
            if(mod(it,SaveIterate)==0)then
                if(SaveStatue==0)then
                    call output('AutoSave',CellCount,.true.)
                else
                    write(char,"(I5.5)"),it
                    call output(char,CellCount,.true.)
                    if(SaveSurface==1) call surface_data(char)
                endif
            if(check==1) call check_node(CellCount)
            endif
        endif
    enddo
endsubroutine solute

!递归对于某个单元的所有叶子结点赋初值
recursive subroutine solution_initial(c)
implicit none
type(gridtyp),pointer::c

if(associated(c%son1))then
	!if(c%lvl<4)then
    if(associated(c%son3))then    
  call solution_initial(c%son1)
  call solution_initial(c%son2)
  call solution_initial(c%son3)
  call solution_initial(c%son4)
!elseif(c%lvl==4)then
  else
  call solution_initial(c%son1)
  call solution_initial(c%son2)
  endif
  else if(c%sort==0)then  !sort限定外部点
  c%u=u0
  c%v=v0
  c%T=T0
  c%rou=rou0
  c%p=R*c%rou*c%T
  
  c%ut=u0
  c%vt=v0
  c%Tt=T0
  c%rout=rou0
  c%pt=R*rou0*T0
  
  c%grax=0
  c%secgrax=0
  c%gray=0
  c%secgray=0
  c%graxp=0
  c%grayp=0
  c%rotx=0
  c%roty=0
  c%divx=0
  c%divy=0
  else                   
  c%u=0.00                
  c%v=0.00
  c%T=0.00
  c%rou=0.00
  c%p=0.00
  
  c%ut=0.00
  c%vt=0.00
  c%Tt=0.00
  c%rout=0.00
  c%pt=0.00
  
  c%grax=0
  c%secgrax=0
  c%gray=0
  c%secgray=0
  c%graxp=0
  c%grayp=0
  c%rotx=0
  c%roty=0
  c%divx=0
  c%divy=0
end if

return
end subroutine solution_initial

recursive subroutine body_boundary(c)
implicit none
type(gridtyp),pointer::c,cr,cu,cd

if(associated(c%son1))then
  if(associated(c%son3))then  
  call body_boundary(c%son1)
  call body_boundary(c%son2)
  call body_boundary(c%son3)
  call body_boundary(c%son4)
  else
  call body_boundary(c%son1)
  call body_boundary(c%son2)
  endif
else if(c%sort==1)then

   call ghost_cell(c)
 
end if

return
end subroutine body_boundary 
!对于某个初始单元迭代到叶子结点赋值或者离散求解
recursive subroutine solution_cell(c)
implicit none
integer::nn
type(gridtyp),pointer::c,cr,cu,cd,cn,cnn

if(associated(c%son1))then
    if(associated(c%son3))then
        call solution_cell(c%son1)
        call solution_cell(c%son2)
        call solution_cell(c%son3)
        call solution_cell(c%son4)
    else
        call solution_cell(c%son1)
        call solution_cell(c%son2)
    endif
    return
endif

!!NRRmod
if(NRR==1)then
    if(c%NRRregion==0)then
        call NRRinterposeBlue(c)
        return
    endif
endif

  !先对上下左赋值，上下左右四个边界也会加密,再调整
if(.not.associated(down_neighbor(c)))then 
  cn=>up_neighbor(c)
  cnn=>up_neighbor(cn)
  call interposeA(cn,cn%rou,cn%u,cn%v,cn%T,cn%p)
  call interposeA(cnn,cnn%rou,cnn%u,cnn%v,cnn%T,cnn%p)
  c%ut=2*cn%ut-cnn%ut
  c%vt=2*cn%vt-cnn%vt
  c%Tt=2*cn%Tt-cnn%Tt
  c%rout=2*cn%rout-cnn%rout
  c%pt=R*c%rout*c%Tt

else if(.not.associated(up_neighbor(c)))then 
  cn=>down_neighbor(c)
  cnn=>down_neighbor(cn)
  call interposeA(cn,cn%rou,cn%u,cn%v,cn%T,cn%p)
  call interposeA(cnn,cnn%rou,cnn%u,cnn%v,cnn%T,cnn%p)
  c%ut=2*cn%ut-cnn%ut
  c%vt=2*cn%vt-cnn%vt
  c%Tt=2*cn%Tt-cnn%Tt
  c%rout=2*cn%rout-cnn%rout
  c%pt=R*c%rout*c%Tt

!else if(mod(c%num,m)==1)then   !入口边界，左边界
else if(.not.associated(left_neighbor(c)))then  
  c%ut=u0
  c%vt=v0
  c%Tt=T0   
  c%rout=rou0 
  c%pt=R*rou0*T0

!else if(mod(c%num,m)==0)then   !右边界，前两点插值(单元值被更新的顺序)
else if(.not.associated(right_neighbor(c)))then 
  cn=>left_neighbor(c)
  cnn=>left_neighbor(cn)
  call interposeA(cn,cn%rou,cn%u,cn%v,cn%T,cn%p)
  call interposeA(cnn,cnn%rou,cnn%u,cnn%v,cnn%T,cnn%p)
  c%ut=2*cn%ut-cnn%ut
  c%vt=2*cn%vt-cnn%vt
  c%Tt=2*cn%Tt-cnn%Tt
  c%rout=2*cn%rout-cnn%rout
  c%pt=R*c%rout*c%Tt
 
!靠近边界的那一层运用一阶迎风+时间前差 
else if(c%num<2*m.or.c%num>(n-2)*m.or.mod(c%num,m)==2.or.mod(c%num+1,m)==0)then 
!输入要更新信息的单元，得到要更新的信息

   call wind(c)
!对剩余外部点运用NND+中心差分+RK方法 ,靠近物面边界也是一阶迎风
else if(c%sort==0)then
    !cl=>left_neighbor(c)
    !cr=>right_neighbor(c)
    !cu=>up_neighbor(c)
    !cd=>down_neighbor(c)
   !if(cl%sort==1.or.cr%sort==1.or.cd%sort==1.or.cu%sort==1)then
   !if(surfcell(c)==1)then   !贴着物面的网格
     !call wind(c)
     !else 
     call nnd_center(c)
     !call wind(c)
   !end if
 !else if(c%sort==0)then    
 !    call wind(c) 
endif
return
end subroutine solution_cell
!加个判定是否贴近物面需要用迎风格式的函数
integer function surfcell(c)
implicit none
type(gridtyp),pointer::c
real(pre)::x1,x2,x3,x4,y1,y2,y3,y4


if((c%center%x-X0/3)**2+(c%center%y-Y0/2)**2<=(2*Rad)**2)then
surfcell=0
else
surfcell=1 
endif
return
end function surfcell

!暂存值更新
subroutine update
implicit none
type(gridtyp),pointer::c
integer i
do i=1,total
  c=>cell(i)
  call update_sub(c)
end do
return
end subroutine update
!调用
recursive subroutine update_sub(c)
implicit none
type(gridtyp),pointer::c
real(pre)::Residut,Residvt,Residrout,Residtt,Residpt

if(associated(c%son1))then
  if(associated(c%son3))then  
  call update_sub(c%son1)
  call update_sub(c%son2)
  call update_sub(c%son3)
  call update_sub(c%son4) 
  else
  call update_sub(c%son1)
  call update_sub(c%son2) 
  endif
else
    if(residual==1)then
        Residut=abs((c%ut-c%u)/c%u);Residu=Residut+Residu
        if(c%v==0)then
            Residvt=0
        else
            Residvt=abs((c%vt-c%v)/c%v);Residv=Residvt+Residv
        endif
        Residrout=abs((c%rout-c%rou)/c%rou);Residrou=Residrout+Residrou
        Residtt=abs((c%tt-c%t)/c%t);Residt=Residtt+Residt
        Residpt=abs((c%pt-c%p)/c%p);Residp=Residpt+Residp
        Residcount=Residcount+1
    endif
    if(check==1)then
        if(maxRu<Residut)then;maxRu=Residut;  maxRCu=c%nelment; endif
        if(maxRv<Residvt)then;maxRv=Residvt;  maxRCv=c%nelment; endif
        if(maxRrou<Residrout)then;maxRrou=Residrout;  maxRCrou=c%nelment; endif
        if(maxRt<Residtt)then;maxRt=Residtt;  maxRCt=c%nelment; endif
        if(maxRp<Residpt)then;maxRp=Residpt;  maxRCp=c%nelment; endif
    endif
    c%u=c%ut
    c%v=c%vt
    c%rou=c%rout
    c%T=c%Tt
    c%p=c%pt
    if(Limiter==1)then  !!sort=0,outside object cells
        if(c%sort==0) return
        if(c%u>LimitUmax)then
            c%u=LimitUmax;  nLUmax=nLUmax+1;  OutofLimit=OutofLimit+1
        endif
        if(c%v>LimitVmax)then
            c%v=LimitVmax;  nLVmax=nLVmax+1;  OutofLimit=OutofLimit+1
        endif
        if(c%rou>LimitRoumax)then
            c%rou=LimitRoumax;nLRoumax=nLRoumax+1;  OutofLimit=OutofLimit+1
        elseif(c%rou<LimitRoumin)then
            c%rou=LimitRoumin;nLRoumin=nLRoumin+1;  OutofLimit=OutofLimit+1
        endif
        if(c%T>LimitTmax)then
            c%T=LimitTmax;  nLTmax=nLTmax+1;  OutofLimit=OutofLimit+1
        elseif(c%T<LimitTmin)then
            c%T=LimitTmin;  nLTmin=nLTmin+1;  OutofLimit=OutofLimit+1
        endif
        if(c%p>LimitPmax)then
            c%p=LimitPmax;  nLPmax=nLPmax+1;  OutofLimit=OutofLimit+1
        elseif(c%p<LimitPmin)then
            c%p=LimitPmin;  nLPmin=nLPmin+1;  OutofLimit=OutofLimit+1
        endif
        if(isnan(c%u))then
            c%u=u0;  nLUmax=nLUmax+1;  OutofLimit=OutofLimit+1
        endif
        if(isnan(c%v))then
            c%v=v0;  nLVmax=nLVmax+1;  OutofLimit=OutofLimit+1
        endif
        if(isnan(c%rou))then
            c%rou=rou0;nLRoumax=nLRoumax+1;  OutofLimit=OutofLimit+1
        endif
        if(isnan(c%t))then
            c%t=t0;  nLTmax=nLTmax+1;  OutofLimit=OutofLimit+1
        endif
        if(isnan(c%p))then
            c%p=R*rou0*T0;nLPmax=nLPmax+1;  OutofLimit=OutofLimit+1
        endif
    endif
endif
return
end subroutine update_sub
!子函数，求外部流场各自的速度梯度和二阶梯度，
subroutine gra_secgra
implicit none
type(gridtyp),pointer::c
integer i
real(pre):: Sgra,Ssecgra
  Sgra=0
  Ssecgra=0
  Srot=0
  Sdiv=0
  Ne=0  
do i=1,total
  c=>cell(i)  
  if(c%sort==0.and.c%num>m.and.(mod(c%num,m)/=1)&
  &.and.(mod(c%num,m)/=0).and.c%num<(n-1)*m)then  !去掉四个边界
 ! if(c%sort==0)then
  !速度梯度和二阶梯度
    call solute_secgra_gra(c)
  endif  
end do

return
end subroutine gra_secgra
!子函数，递归求子单元旋度，散度和单元计数叠加
recursive subroutine solute_secgra_gra(c)
implicit none
type(gridtyp),pointer::c
  if(.not.associated(c%son1).and.c%sort==0)then 
    if(NRRregionAutoRefine==0)then
      if(c%NRRregion>=0) return
    endif
  call Ugrax(c) 
  !call Ugray(c) 
  !call Usecgrax(c)
  !call Usecgray(c)
  
  !Sgra=Sgra+c%grax**2+c%gray**2
  !Ssecgra=Ssecgra+c%secgrax**2+c%secgray**2
  !Srot=Srot+c%rotx**2+c%roty**2
  !Sdiv=Sdiv+c%divx**2+c%divy**2
  Srot=Srot+c%rotx**2
  Sdiv=Sdiv+c%divx**2
  
  Ne=Ne+1
  elseif(associated(c%son1))then
      if(associated(c%son3))then                                 
  call solute_secgra_gra(c%son1) 
  call solute_secgra_gra(c%son2)
  call solute_secgra_gra(c%son3)
  call solute_secgra_gra(c%son4)
      else
  call solute_secgra_gra(c%son1) 
  call solute_secgra_gra(c%son2)
      endif
      
 endif 
return
end subroutine solute_secgra_gra
!子函数，解自适应
subroutine solute_ref
implicit none
type(gridtyp),pointer::c
integer i
real(pre)::asecgra,agra,arot,adiv
 
 !agra=sqrt(Sgra/(2*Ne))   !方差
 !asecgra=sqrt(Ssecgra/(2*Ne))
 
 !arot=sqrt(Srot/(2*Ne))  !方差
 !adiv=sqrt(Sdiv/(2*Ne))
  arot=sqrt(Srot/Ne)  !各向同性方差
 adiv=sqrt(Sdiv/Ne)
 
 
do i=1,total
   c=>cell(i)       

   if(c%cross/=2.and.c%num>m.and.mod(c%num,m)/=1&    
  &.and.mod(c%num,m)/=0.and.c%num<(n-1)*m)then 
 
    call solute_ref_sub_iso(c,agra,asecgra,arot,adiv)     !自适应加密的子函数 
  endif
  
   if(c%cross/=2.and.associated(c%son1))then     
   call solute_coarse_sub_iso(c,agra,asecgra,arot,adiv)   !粗化的子函数 
   endif
 
end do  
return
end subroutine solute_ref
!粗化单元子函数


recursive subroutine solute_coarse_sub_iso(c,agra,asecgra,arot,adiv)!各项同性笛卡尔网格粗化
implicit none
type(gridtyp),pointer::c   
real(pre)::asecgra,agra,arot,adiv
!外部非相交单元，且有子单元，且子单元为叶子单元
if(associated(c%son1))then!各项同性
if(c%lvl<solute_lvl_max.and.&
&((c%center%x-X0/4)**2+(c%center%y-Y0/2)**2>=(Rad+0.06)**2)&
&.and..not.associated(c%son1%son1).and..not.associated(c%son2%son1).and.&
&.not.associated(c%son3%son1).and..not.associated(c%son4%son1))then  

  !变系数的加粗判据
  
if((c%son1%rotx+c%son2%rotx+c%son3%rotx+c%son4%rotx<RotParameter2*arot).and.&
    &(c%son1%divx+c%son2%divx+c%son3%divx+c%son4%divx<DivParameter2*adiv))then
    if(NRRregionAutoRefine==0)then
        if(c%son1%NRRregion>=0)return
        if(c%son2%NRRregion>=0)return
        if(c%son3%NRRregion>=0)return
        if(c%son4%NRRregion>=0)return
    endif
    call gridmodify_coarse_iso(c)  !单元粗化具体操作的函数
  end if    
 elseif(c%sort==0.and.c%cross/=1)then
  if(associated(c%son1%son1))then
    call solute_coarse_sub_iso(c%son1,agra,asecgra,arot,adiv) 
  endif
  if(associated(c%son2%son1))then
    call solute_coarse_sub_iso(c%son2,agra,asecgra,arot,adiv)
  endif
  if(associated(c%son3%son1))then  
    call solute_coarse_sub_iso(c%son3,agra,asecgra,arot,adiv)
  endif
  if(associated(c%son4%son1))then
    call solute_coarse_sub_iso(c%son4,agra,asecgra,arot,adiv)
  end if
 endif      
endif

return
end subroutine solute_coarse_sub_iso

!粗化的具体操作
subroutine gridmodify_coarse_iso(c)
implicit none
type(gridtyp),pointer::c
integer i
!先获得粗化单元的流场值
call interposeA(c,c%rou,c%u,c%v,c%T,c%p) 
!再删除子单元
if(associated(c%son1).and.associated(c%son3))then
deallocate(c%son1,c%son2,c%son3,c%son4)
CellCount=CellCount-3
elseif(associated(c%son1).and..not.associated(c%son3))then
deallocate(c%son1,c%son2)   
CellCount=CellCount-1
endif
!再置空子单元
if(associated(c%son1).and.associated(c%son3))then
nullify(c%son1,c%son2,c%son3,c%son4)
elseif(associated(c%son1).and..not.associated(c%son3))then
nullify(c%son1,c%son2)
endif

return
end subroutine gridmodify_coarse_iso

recursive subroutine solute_ref_sub_iso(c,agra,asecgra,arot,adiv)
implicit none
type(gridtyp),pointer::c
real(pre)::asecgra,agra,arot,adiv

if(.not.associated(c%son1).and.c%sort==0.and.c%lvl<solute_lvl_max)then  
    if(NRRregionAutoRefine==0)then
        if(c%NRRregion>=0) return
    endif
if(c%rotx>RotParameter1*arot.or.c%divx>DivParameter1*adiv)then
    c%typ=0
    call newson(c)
    if(NRRregionAutoRefine==0)then
        if(NRR==1)then
            call NRRregionMark(c%son1)
            call NRRregionMark(c%son2)
            call NRRregionMark(c%son3)
            call NRRregionMark(c%son4)
        endif
    endif
    endif

elseif(associated(c%son1).and.c%sort==0.and.c%lvl<solute_lvl_max)then
  
  call solute_ref_sub_iso(c%son1,agra,asecgra,arot,adiv) 
  call solute_ref_sub_iso(c%son2,agra,asecgra,arot,adiv)
  call solute_ref_sub_iso(c%son3,agra,asecgra,arot,adiv)
  call solute_ref_sub_iso(c%son4,agra,asecgra,arot,adiv)
   
endif
return
end subroutine solute_ref_sub_iso


!解自适应后再消除悬差网格
subroutine solute_gridmodify
implicit none
type(gridtyp),pointer::c
integer i,j
 
    do j=1,solute_lvl_max
    !do j=1,2
            
		do i=1,total
		c=>cell(i) 
		call init_s(c) 
		end do
		do i=1,total
		c=>cell(i) 
        if(c%num>m.and.mod(c%num,m)/=1.and.mod(c%num,m)/=0.and.c%num<(n-1)*m)then
		call gxyxmodify_def(c) 
        endif
		end do
		do i=1,total
		c=>cell(i) 
		call gxyxmodify_ref2(c) 
        end do
  !      do i=1,total
		!c=>cell(i) 
  !      call gxyx_level_modify(c)
  !      end do 
	end do

return
end subroutine solute_gridmodify

!subroutine solute_gridmodify2
!implicit none
!type(gridtyp),pointer::c
!integer i,j
! 
!    do j=1,solute_lvl_max
!    !do j=1,2
!            
!		do i=1,total
!		c=>cell(i) 
!		call init_s(c) 
!		end do
!		do i=1,total
!		c=>cell(i) 
!		call gxyxmodify_def(c) 
!		end do
!		do i=1,total
!		c=>cell(i) 
!		call gxyxmodify_ref2(c) 
!        end do
!  !      do i=1,total
!		!c=>cell(i) 
!  !      call gxyx_level_modify(c)
!  !      end do 
!	end do
!
!return
!end subroutine solute_gridmodify2
!速度散度和旋度
subroutine Ugrax(c)
implicit none
type(gridtyp),pointer::c
real(pre)::hc
character(6)::ua="u",va="v",xa="x",xxa="xx",ya="y",yya="yy",pa="p",maa="ma1",roua="rou"

hc=h/2**c%lvl
c%rotx=abs(duv(c,va,xa)-duv(c,ua,ya))*hc**1.5
c%divx=abs(duv(c,ua,xa)+duv(c,va,ya))*hc**1.5

!c%grax=abs(duv(c,ua,xxa))*hc**1.5 !速度梯度
!c%graxp=abs(duv(c,pa,xa))*hc**1.5 !压力梯度
!c%grax=abs(duv(c,maa,xa))*hc**1.5 !马赫数梯度
!c%graxp=abs(duv(c,roua,xa))*hc**1.5 !密度梯度
return
end subroutine Ugrax
!散度
!subroutine Usecgrax(c)
!implicit none  
!type(gridtyp),pointer::c
!real(pre)::hc
!character(6)::ua="u",va="v",xxa="xx",ya="y"
!hc=h/2**c%lvlx
!c%secgrax=abs(duv(c,ua,xxa))*hc**3
!return
!end subroutine Usecgrax

subroutine Ugray(c)
implicit none
type(gridtyp),pointer::c
real(pre)::hc
character(6)::ua="u",va="v",xa="x",xxa="xx",ya="y",yya="yy",pa="p",maa="ma1",roua="rou"

!hc=h/2**c%lvl
!c%roty=abs(duv(c,va,xa)-duv(c,ua,ya))*hc**1.5
!c%divy=abs(duv(c,ua,xa)+duv(c,va,ya))*hc**1.5

!c%gray=abs(duv(c,va,yya))*hc**1.5 !速度梯度
!c%grayp=abs(duv(c,pa,ya))*hc**1.5 !压力梯度
!c%gray=abs(duv(c,maa,ya))*hc**1.5 !马赫数梯度
!c%grayp=abs(duv(c,roua,ya))*hc**1.5 !密度梯度
return
end subroutine Ugray
!散度
!subroutine Usecgray(c)
!implicit none  
!type(gridtyp),pointer::c
!real(pre)::hc
!character(6)::ua="u",va="v",xa="x",yya="yy"
!hc=h/2**c%lvly
!c%secgray=abs(duv(c,va,yya))*hc**3
!return
!end subroutine Usecgray



end module solution
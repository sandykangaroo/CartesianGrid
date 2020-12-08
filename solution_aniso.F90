!����������
Module solution_aniso
!use mesh
use inflowGlobal
use bc
use nndcenter
use windd
use outfile
use flow_result
use NRRMeshMod
use precisionMod
use typeMod

implicit none
contains

!�Ӻ��������ⲿ�������Ե��ٶ��ݶȺͶ����ݶȣ�
subroutine gra_secgra_aniso
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
  &.and.(mod(c%num,m)/=0).and.c%num<(n-1)*m)then  !ȥ���ĸ��߽�
 ! if(c%sort==0)then
  !�ٶ��ݶȺͶ����ݶ�
    call solute_secgra_gra_aniso(c)
  endif  
end do

return
end subroutine gra_secgra_aniso
!�Ӻ������ݹ����ӵ�Ԫ���ȣ�ɢ�Ⱥ͵�Ԫ��������
recursive subroutine solute_secgra_gra_aniso(c)
implicit none
type(gridtyp),pointer::c
  if(.not.associated(c%son1).and.c%sort==0)then 
    if(NRRregionAutoRefine==0)then
      if(c%NRRregion>=0) return
    endif
  call Ugrax_aniso(c) 
  call Ugray_aniso(c) 
  !call Usecgrax(c)
  !call Usecgray(c)
  
  Sgra=Sgra+c%grax**2+c%gray**2
  !Ssecgra=Ssecgra+c%secgrax**2+c%secgray**2
  Srot=Srot+c%rotx**2+c%roty**2
  Sdiv=Sdiv+c%divx**2+c%divy**2
  !Srot=Srot+c%rotx**2
  !Sdiv=Sdiv+c%divx**2
  
  Ne=Ne+1
  elseif(associated(c%son1))then
      if(associated(c%son3))then                                 
  call solute_secgra_gra_aniso(c%son1) 
  call solute_secgra_gra_aniso(c%son2)
  call solute_secgra_gra_aniso(c%son3)
  call solute_secgra_gra_aniso(c%son4)
      else
  call solute_secgra_gra_aniso(c%son1) 
  call solute_secgra_gra_aniso(c%son2)
      endif
      
 endif 
return
end subroutine solute_secgra_gra_aniso
!�Ӻ�����������Ӧ
subroutine solute_ref_aniso
implicit none
type(gridtyp),pointer::c
integer i
real(pre)::asecgra,agra,arot,adiv
 
 agra=sqrt(Sgra/(2*Ne))   !����
 !asecgra=sqrt(Ssecgra/(2*Ne))
 
 !arot=sqrt(Srot/(2*Ne))  !����
 !adiv=sqrt(Sdiv/(2*Ne))
  arot=sqrt(Srot/(2*Ne))  !����ͬ�Է���
  adiv=sqrt(Sdiv/(2*Ne))
 
 
do i=1,total
   c=>cell(i)       

   if(c%cross/=2.and.c%num>m.and.mod(c%num,m)/=1&    
  &.and.mod(c%num,m)/=0.and.c%num<(n-1)*m)then 
 
    call solute_ref_sub_aniso(c,agra,asecgra,arot,adiv)     !����Ӧ���ܵ��Ӻ��� 
  endif
  
   if(c%cross/=2.and.associated(c%son1))then     
   call solute_coarse_sub_aniso(c,agra,asecgra,arot,adiv)   !�ֻ����Ӻ��� 
   endif
 
end do  
return
end subroutine solute_ref_aniso
!�ֻ���Ԫ�Ӻ���
recursive subroutine solute_coarse_sub_aniso(c,agra,asecgra,arot,adiv)
implicit none
type(gridtyp),pointer::c   
real(pre)::asecgra,agra,arot,adiv
!�ⲿ���ཻ��Ԫ�������ӵ�Ԫ�����ӵ�ԪΪҶ�ӵ�Ԫ
if(associated(c%son1).and.associated(c%son3))then!����ͬ��
if(((c%center%x-X0/4)**2+(c%center%y-Y0/2)**2>=(Rad+0.1)**2)&
&.and..not.associated(c%son1%son1).and..not.associated(c%son2%son1).and.&
&.not.associated(c%son3%son1).and..not.associated(c%son4%son1))then  
    if(NRRregionAutoRefine==0)then
        if(c%son1%NRRregion>=0)return
        if(c%son2%NRRregion>=0)return
        if(c%son3%NRRregion>=0)return
        if(c%son4%NRRregion>=0)return
    endif

  !��ϵ���ļӴ��о�
  if((c%son1%grax+c%son2%grax+c%son3%grax+c%son4%grax<0.2*agra).or.&
  &(c%son1%gray+c%son2%gray+c%son3%gray+c%son4%gray<0.2*agra))then
       
  !  if((c%son1%rotx+c%son2%rotx+c%son3%rotx+c%son4%rotx+&
  !      &c%son1%roty+c%son2%roty+c%son3%roty+c%son4%roty<5.6*arot).or.&
  !  &(c%son1%divx+c%son2%divx+c%son3%divx+c%son4%divx&
  !&+c%son1%divy+c%son2%divy+c%son3%divy+c%son4%divy<5.6*adiv))then
        
     call gridmodify_coarse_aniso(c)  !��Ԫ�ֻ���������ĺ���
  end if    
 elseif(c%sort==0.and.c%cross/=1)then
  if(associated(c%son1%son1))then
    call solute_coarse_sub_aniso(c%son1,agra,asecgra,arot,adiv) 
  endif
  if(associated(c%son2%son1))then
    call solute_coarse_sub_aniso(c%son2,agra,asecgra,arot,adiv)
  endif
  if(associated(c%son3%son1))then  
    call solute_coarse_sub_aniso(c%son3,agra,asecgra,arot,adiv)
  endif
  if(associated(c%son4%son1))then
    call solute_coarse_sub_aniso(c%son4,agra,asecgra,arot,adiv)
  end if
 endif 
 
elseif(associated(c%son1).and..not.associated(c%son3))then!�����������
   if((c%center%x-X0/4)**2+(c%center%y-Y0/2)**2>=(Rad+0.1)**2&
&.and..not.associated(c%son1%son1).and..not.associated(c%son2%son1))then  
        
    if(((c%son1%grax+c%son2%grax<0.1*agra).and.c%son1%spl==2)&
        &.or.((c%son1%gray+c%son2%gray<0.1*agra).and.c%son1%spl==1))then
  !if((c%son1%rotx+c%son2%rotx+c%son1%roty+c%son2%roty<2.8*arot).or.&
  !  &(c%son1%divx+c%son2%divx+c%son1%divy+c%son2%divy<2.8*adiv))then     
       
     call gridmodify_coarse_aniso(c)  !��Ԫ�ֻ���������ĺ���
  end if    
 elseif(associated(c%son1).and.c%sort==0.and.c%cross/=1)then
  if(associated(c%son1%son1))then
    call solute_coarse_sub_aniso(c%son1,agra,asecgra,arot,adiv) 
  endif
  if(associated(c%son2%son1))then
    call solute_coarse_sub_aniso(c%son2,agra,asecgra,arot,adiv)
  endif
  
 endif   
       
endif

return
end subroutine solute_coarse_sub_aniso

subroutine gridmodify_coarse_aniso(c)
implicit none
type(gridtyp),pointer::c
integer i
!�Ȼ�ôֻ���Ԫ������ֵ
call interposeA(c,c%rou,c%u,c%v,c%T,c%p) 
!��ɾ���ӵ�Ԫ
if(associated(c%son1).and.associated(c%son3))then
deallocate(c%son1,c%son2,c%son3,c%son4)
CellCount=CellCount-3
elseif(associated(c%son1).and..not.associated(c%son3))then
deallocate(c%son1,c%son2)   
CellCount=CellCount-1
endif
!���ÿ��ӵ�Ԫ
if(associated(c%son1).and.associated(c%son3))then
nullify(c%son1,c%son2,c%son3,c%son4)
elseif(associated(c%son1).and..not.associated(c%son3))then
nullify(c%son1,c%son2)
endif

return
end subroutine gridmodify_coarse_aniso


!�Ӻ������ݹ����ӵ�Ԫ���ȣ�ɢ�Ⱥ͵�Ԫ��������
recursive subroutine solute_ref_sub_aniso(c,agra,asecgra,arot,adiv)
!use inflowGlobal,only:NRR
implicit none
type(gridtyp),pointer::c
real(pre)::asecgra,agra,arot,adiv
logical::newcell

newcell=0
if((.not.associated(c%son1).and.c%sort==0.and.c%lvlx<solute_lvl_max).and.&
    &(.not.associated(c%son1).and.c%sort==0.and.c%lvly<solute_lvl_max))then  
    if(NRRregionAutoRefine==0)then
        if(c%NRRregion>=0) return
    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ����һ�ֱ������ٶȡ�ѹ�����ܶȣ�����   
  if(c%grax>1.0*agra.and.c%gray<0.5*agra)then
     if(c%lvly>=c%lvlx)then
       c%typ=2; newcell=1
    else
        c%typ=0; newcell=1
    endif
  elseif(c%gray>1.0*agra.and.c%grax<0.5*agra)then    
     if(c%lvlx>=c%lvly)then
       c%typ=1; newcell=1
    else
        c%typ=0; newcell=1
    endif 
  elseif(c%gray>1.0*agra.and.c%grax>1.0*agra)then 
     c%typ=0; newcell=1
  end if 
  
  if(newcell==1)then
    call newson(c) 
    if(NRRregionAutoRefine==0)then
        if(NRR==1)then
            call NRRregionMark(c%son1)
            call NRRregionMark(c%son2)
            if(c%typ/=0) return
            call NRRregionMark(c%son3)
            call NRRregionMark(c%son4)
        endif
    endif
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  �����ٶȵ�ɢ�������ж��Ƿ����

!if((c%rotx>1.0*arot.and.c%roty>1.0*arot).or.(c%divx>1.0*adiv.and.c%divy>1.0*adiv))then
!       c%typ=0
!       call newson(c) 
!elseif((c%rotx>1.0*arot.and.c%roty<1.0*arot).or.(c%divx>1.0*adiv.and.c%divy<1.0*adiv))then
!    if(c%lvly>=c%lvlx)then
!       c%typ=2
!       call newson(c)
!    else
!        c%typ=0
!       call newson(c)
!    endif
!elseif((c%rotx<1.0*arot.and.c%roty>1.0*arot).or.(c%divx<1.0*adiv.and.c%divy>1.0*adiv))then
!    if(c%lvlx>=c%lvly)then
!       c%typ=1
!       call newson(c)
!    else
!        c%typ=0
!       call newson(c)
!    endif
!endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!if(c%rotx+c%roty>2.6*arot.or.c%divx+c%divy>2.6*arot)then
!    if(c%rotx>=c%roty.and.c%divx>=c%divy.and.c%grax>=c%gray)then!!.and.c%lvlx<=c%lvly
!       c%typ=2
!       call newson(c) 
!    elseif(c%rotx<=c%roty.and.c%divx<=c%divy.and.c%grax<=c%gray)then!!.and.c%lvlx>=c%lvly
!       c%typ=1
!       call newson(c)
!    else
!       c%typ=0
!       call newson(c)
!    endif
!endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    elseif((associated(c%son1).and.c%sort==0.and.c%lvlx<solute_lvl_max).and.&
    &(associated(c%son1).and.c%sort==0.and.c%lvly<solute_lvl_max))then
    if(associated(c%son3))then
  call solute_ref_sub_aniso(c%son1,agra,asecgra,arot,adiv) 
  call solute_ref_sub_aniso(c%son2,agra,asecgra,arot,adiv)
  call solute_ref_sub_aniso(c%son3,agra,asecgra,arot,adiv)
  call solute_ref_sub_aniso(c%son4,agra,asecgra,arot,adiv)
    else
  call solute_ref_sub_aniso(c%son1,agra,asecgra,arot,adiv) 
  call solute_ref_sub_aniso(c%son2,agra,asecgra,arot,adiv)
    endif
    
endif 
return
end subroutine solute_ref_sub_aniso


!�ٶ�ɢ�Ⱥ�����
subroutine Ugrax_aniso(c)
implicit none
type(gridtyp),pointer::c
real(pre)::hc
character(6)::ua="u",va="v",xa="x",xxa="xx",ya="y",yya="yy",pa="p",maa="ma1",roua="rou"

hc=h/2**c%lvlx
c%rotx=abs(duv(c,va,xa)-duv(c,ua,ya))*hc**1.5
c%divx=abs(duv(c,ua,xa)+duv(c,va,ya))*hc**1.5

!c%grax=abs(duv(c,ua,xxa))*hc**1.5 !�ٶ��ݶ�
c%grax=abs(duv(c,pa,xa))*hc**1.5 !ѹ���ݶ�
!c%grax=abs(duv(c,maa,xa))*hc**1.5 !�������ݶ�
!c%grax=abs(duv(c,roua,xa))*hc**1.5 !�ܶ��ݶ�
return
end subroutine Ugrax_aniso
!ɢ��
!subroutine Usecgrax(c)
!implicit none  
!type(gridtyp),pointer::c
!real(pre)::hc
!character(6)::ua="u",va="v",xxa="xx",ya="y"
!hc=h/2**c%lvlx
!c%secgrax=abs(duv(c,ua,xxa))*hc**3
!return
!end subroutine Usecgrax

subroutine Ugray_aniso(c)
implicit none
type(gridtyp),pointer::c
real(pre)::hc
character(6)::ua="u",va="v",xa="x",xxa="xx",ya="y",yya="yy",pa="p",maa="ma1",roua="rou"

hc=h/2**c%lvly
c%roty=abs(duv(c,va,xa)-duv(c,ua,ya))*hc**1.5
c%divy=abs(duv(c,ua,xa)+duv(c,va,ya))*hc**1.5

!c%gray=abs(duv(c,va,yya))*hc**1.5 !�ٶ��ݶ�
c%gray=abs(duv(c,pa,ya))*hc**1.5 !ѹ���ݶ�
!c%gray=abs(duv(c,maa,ya))*hc**1.5 !�������ݶ�
!c%gray=abs(duv(c,roua,ya))*hc**1.5 !�ܶ��ݶ�
return
end subroutine Ugray_aniso
!ɢ��
!subroutine Usecgray(c)
!implicit none  
!type(gridtyp),pointer::c
!real(pre)::hc
!character(6)::ua="u",va="v",xa="x",yya="yy"
!hc=h/2**c%lvly
!c%secgray=abs(duv(c,va,yya))*hc**3
!return
!end subroutine Usecgray



end module solution_aniso
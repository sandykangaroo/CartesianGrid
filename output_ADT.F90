module outfile
!use solution
use precisionMod
use mesh
use neighbor
use inflow
use typeMod
implicit none
integer elmentnum,lvldelta,temp
  !存储节点坐标的数组，大小可变，每个单元是piont类型
  type(point),allocatable::Node(:)
  !记录节点数目的值，以及单元数目
  integer,save::nNode,nElment

!integer::error
contains
subroutine output(char,ccount,stat)
implicit none
integer i
logical::stat   !!=0, write mesh; =1, write data.
character(*)::char
type(gridtyp),pointer::c
integer,intent(in)::ccount
character(40)::filename


filename='./Others/'//trim(FileNameStr)//'-'//char//'.dat'
!filename='./Others/'//FileNameStr//char//'.dat'
print*,"_______Start data output_______"

allocate(Node(ccount*2))

call nodeinf
write(*,"(A7,I10,5X,A5,I10)") "Elment:", nElment, "Node:", nNode
open(50,file=filename,STATUS='REPLACE',FORM='FORMATTED')
write(50,*)'TITLE="2D Results"'

if(stat==0)then
    write(50,*)'VARIABLES="X","Y","LVL","NRR"'
    write(50,*)'ZONE N=',nNode,'E=',nElment,'ZONETYPE=FEQUADRILATERAL'
    write(50,*)'DATAPACKING=BLOCK'
    write(50,*)'VARLOCATION=([1-2]=NODAL,[3-4]=CELLCENTERED)'
else
    write(50,*)'VARIABLES="X","Y","LVL","NRR","U","V","Rou","T","P","Ma"'
    write(50,*)'ZONE N=',nNode,'E=',nElment,'ZONETYPE=FEQUADRILATERAL'
    write(50,*)'DATAPACKING=BLOCK'
    write(50,*)'VARLOCATION=([1-2]=NODAL,[3-10]=CELLCENTERED)'
endif


    do i=1,nNode
        if(pre==my_r4)then
            write(50,"(F9.5)")Node(i)%x
        else
            write(50,"(F12.8)")Node(i)%x
        endif
    end do
    do i=1,nNode
        if(pre==my_r4)then
            write(50,"(F9.5)")Node(i)%y
        else
            write(50,"(F12.8)")Node(i)%y
        endif
    end do
    deallocate(Node)

if(stat==0)then
    do i=1,total
      c=>cell(i)
      call WriteLvL(c)
    end do
    do i=1,total
      c=>cell(i)
      call WriteNRR(c)
    end do
else
    do i=1,total
      c=>cell(i)
      call WriteLvL(c)
    end do
    do i=1,total
      c=>cell(i)
      call WriteNRR(c)
    end do
    do i=1,total
      c=>cell(i)
      call writeout_u(c)
    end do
    do i=1,total
      c=>cell(i)
      call writeout_v(c)
    end do
    do i=1,total
      c=>cell(i)
      call writeout_rou(c)
    end do
    do i=1,total
      c=>cell(i)
      call writeout_T(c)
    end do
    do i=1,total
      c=>cell(i)
      call writeout_P(c)
    end do
    do i=1,total
      c=>cell(i)
      call writeout_Ma(c)
    end do
endif
    !链接关系
    do i=1,total
      c=>cell(i)
      call link(c)
    end do
close (50)
print*,"Save to file:     ",filename
return
end subroutine output
!结点关系，递归到子单元
recursive subroutine link(c)
  implicit none
  type(gridtyp),pointer::c
       if(associated(c%son1)) then   
	       if(associated(c%son3))then
           call link(c%son1)
           call link(c%son2)
           call link(c%son3)
           call link(c%son4)
		   else
           call link(c%son1)
           call link(c%son2)
		   endif
         else
            write(50,*)c%cNode(1),c%cNode(2),c%cNode(3),c%cNode(4)
       end if        
return
end subroutine link

!!输出各个物理量的子函数
recursive subroutine WriteLvL(c)
  implicit none
  type(gridtyp),pointer::c
       if(associated(c%son1)) then 
           if(associated(c%son3))then
           call WriteLvL(c%son1)
           call WriteLvL(c%son2)
           call WriteLvL(c%son3)
           call WriteLvL(c%son4)
           else
           call WriteLvL(c%son1)
           call WriteLvL(c%son2)
           endif
       else
            write(50,"(I3)")c%lvl
       endif
return
end subroutine WriteLvL

recursive subroutine WriteNRR(c)
  implicit none
  type(gridtyp),pointer::c
       if(associated(c%son1)) then 
           if(associated(c%son3))then
           call WriteNRR(c%son1)
           call WriteNRR(c%son2)
           call WriteNRR(c%son3)
           call WriteNRR(c%son4)
           else
           call WriteNRR(c%son1)
           call WriteNRR(c%son2)
           endif
       else
           if(c%NRRregion<0)then
               c%NRRregion=-1
           endif
           write(50,*)c%NRRregion
       endif        
return
end subroutine WriteNRR


recursive subroutine writeout_u(c)
  implicit none
  type(gridtyp),pointer::c
       if(associated(c%son1)) then 
           if(associated(c%son3))then
           call writeout_u(c%son1)
           call writeout_u(c%son2)
           call writeout_u(c%son3)
           call writeout_u(c%son4)
           else
           call writeout_u(c%son1)
           call writeout_u(c%son2)
           endif
       else
           !write(50,"(2I4,7F20.7)")c%num,c%sort,c%center%x,c%center%y,c%u,c%v,c%rou,c%T,c%p
            if(c%sort==0)then
                if(isnan(c%u))then
                    write(50,*),'-999.9'
                else
                    write(50,*)c%u
                endif
            !write(103,*)c%center%x,c%center%y,c%u
            !if(c%lvlx==c%lvly.and.c%spl>0) then
            !    write(104,*)c%center%x,c%center%y
            !endif
            else
            write(50,"(F5.1)")0.0
            !write(*,*)0
            end if
       end if        
return
end subroutine writeout_u
!
recursive subroutine writeout_v(c)
  implicit none
  type(gridtyp),pointer::c
       if(associated(c%son1)) then
           if(associated(c%son3))then
           call writeout_v(c%son1)
           call writeout_v(c%son2)
           call writeout_v(c%son3)
           call writeout_v(c%son4)
           else
               call writeout_v(c%son1)
               call writeout_v(c%son2)
           endif
           
         else
           
            if(c%sort==0)then
                if(isnan(c%v))then
                    write(50,*),'-999.9'
                else
                    write(50,*)c%v
                endif
            
            else
            write(50,"(F5.1)")0.0
            end if
       end if        
return
end subroutine writeout_v
!
recursive subroutine writeout_rou(c)
  implicit none
  type(gridtyp),pointer::c
       if(associated(c%son1)) then
           if(associated(c%son3))then
           call writeout_rou(c%son1)
           call writeout_rou(c%son2)
           call writeout_rou(c%son3)
           call writeout_rou(c%son4)
           else
              call writeout_rou(c%son1)
              call writeout_rou(c%son2)
           endif
           
         else
            if(c%sort==0)then
                if(isnan(c%rou))then
                    write(50,*),'-999.9'
                else
                    write(50,*)c%rou
                endif
            !write(103,*)c%center%x,c%center%y,c%v
            else
            write(50,"(F5.1)")0.0
            end if
       end if        
return
end subroutine writeout_rou

recursive subroutine writeout_T(c)
  implicit none
  type(gridtyp),pointer::c
       if(associated(c%son1)) then 
           if(associated(c%son3))then
           call writeout_T(c%son1)
           call writeout_T(c%son2)
           call writeout_T(c%son3)
           call writeout_T(c%son4)
           else
              call writeout_T(c%son1)
              call writeout_T(c%son2) 
           endif
           
         else
            if(c%sort==0)then
                if(isnan(c%T))then
                    write(50,*),'-999.9'
                else
                    write(50,*)c%T
                endif
            else
            write(50,"(F5.1)")0.0
            end if
       end if        
return
end subroutine writeout_T

recursive subroutine writeout_P(c)
  implicit none
  type(gridtyp),pointer::c
       if(associated(c%son1)) then
           if(associated(c%son3))then
           call writeout_P(c%son1)
           call writeout_P(c%son2)
           call writeout_P(c%son3)
           call writeout_P(c%son4)
           else
               call writeout_P(c%son1)
               call writeout_P(c%son2)
           endif

         else
            if(c%sort==0)then
                if(isnan(c%p))then
                    write(50,*),'-999.9'
                else
                    write(50,*)c%p
                endif
            else
            write(50,"(F5.1)")0.0
            end if
       end if
return
end subroutine writeout_P
        
recursive subroutine writeout_Ma(c)
  implicit none
  type(gridtyp),pointer::c
  real(pre)::Ma_c,a_c,u_c  !当地马赫数，速度和a
       if(associated(c%son1)) then
           if(associated(c%son3))then
           call writeout_Ma(c%son1)
           call writeout_Ma(c%son2)
           call writeout_Ma(c%son3)
           call writeout_Ma(c%son4)
           else
                call writeout_Ma(c%son1)
                call writeout_Ma(c%son2)
           endif
           
         else
            if(c%sort==0)then
            !求Ma
                a_c=sqrt(abs(gama*c%p/c%rou))
                u_c=sqrt(c%u**2+c%v**2)
                Ma_c=u_c/a_c
                if(isnan(Ma_c))then
                    write(50,*),'-999.9'
                else
                    write(50,*)Ma_c
                endif
            else
            write(50,"(F5.1)")0.0
            end if
       end if                 
return
end subroutine writeout_Ma
!!++++++++++++++++++++++++++++++++++++++++++++!!
!ADT判断节点序号和对应坐标，按照初始单元顺序进行
subroutine nodeinf   !节点信息
implicit none
integer i,j
type(gridtyp),pointer::c
!按顺序将每个单元的节点加入Node(:)中
!先对第一个单元赋初值


nNode=4
nElment=1

cell(1)%cNode(1)=1
cell(1)%cNode(2)=2
cell(1)%cNode(3)=3
cell(1)%cNode(4)=4
Node(1)%x=0.0
Node(1)%y=0.0
Node(2)%x=h
Node(2)%y=0.0
Node(3)%x=h
Node(3)%y=h
Node(4)%x=0.0
Node(4)%y=h
!先将单元顶点序号都赋值为-1，方便判断未处理的点
cell(1)%nelment=nelment
do i=2,total
   c=>cell(i)
   call initial_node(c)
end do
!error=0
do i=2,total
   c=>cell(i)
   call addnode(c)
end do
return
end subroutine nodeinf
!初始化每个单元顶点
recursive subroutine initial_node(c)
implicit none
type(gridtyp),pointer::c

 if(.not.associated(c%son1))then
 c%cnode(1)=-1
 c%cnode(2)=-1
 c%cnode(3)=-1
 c%cnode(4)=-1
 else 
    if(c%son1%spl==1.or.c%son1%spl==2)then
    call initial_node(c%son1)
    call initial_node(c%son2)
    else
    call initial_node(c%son1)
    call initial_node(c%son2)
    call initial_node(c%son3)
    call initial_node(c%son4)
    endif
 endif
return
end subroutine initial_node

!每个单元，递归到叶子节点，再加到Node(:)里，注意1234的顺序
recursive subroutine addnode(c)
implicit none
type(point)::nt(4)
type(gridtyp),pointer::c,ct,cx,cxx,cxxx,cn,cnn,cnnn,ctemp1,ctemp2
real(pre)::ntx,nty,xmin,xmax,ymin,ymax
integer i,j,nn,nend,gxyxsd
if(associated(c%son1))then
    if(associated(c%son3))then
    call addnode(c%son1)
    call addnode(c%son2)
    call addnode(c%son3)
    call addnode(c%son4)
    else
    call addnode(c%son1)
    call addnode(c%son2)
    endif

else
   nElment=nElment+1       !递归到子单元，则将该单元计数
   c%nelment=nelment
   
!!错误位置指定
!if(nElment==PositionElement)then
!continue
!endif
!if(nNode==PositionNode)then
!continue
!endif
   
!邻居法确定重复点以及扩展点
if(c%sp==0)then     !无子单元的初始单元分类（主要是边界分类讨论），按遍历顺序
                  if(mod(c%num,m)==1)then     !左边界
                    cn=>down_neighbor(c)
                     if(associated(cn%son1))then
                     call corner(cn)
                     end if
    
                    c%cnode(1)=cn%cnode(4)
                    c%cnode(2)=cn%cnode(3)
    
                    nNode=nNode+1
                    c%cnode(3)=nNode
                    Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
                    Node(nNode)%y=c%center%y+h/2**(c%lvly+1)
   
                    nNode=nNode+1
                    c%cnode(4)=nNode
                    Node(nNode)%x=c%center%x-h/2**(c%lvlx+1)
                    Node(nNode)%y=c%center%y+h/2**(c%lvly+1)
                    elseif(c%num<=m)then    !下边界
                    cn=>left_neighbor(c)
					if(associated(cn%son1))then
					call corner(cn)
					end if
    
                    c%cnode(1)=cn%cnode(2)
                    c%cnode(4)=cn%cnode(3)
    
                    nNode=nNode+1
                    c%cnode(2)=nNode
                    Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
                    Node(nNode)%y=c%center%y-h/2**(c%lvly+1)
   
                    nNode=nNode+1
                    c%cnode(3)=nNode
                    Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
                    Node(nNode)%y=c%center%y+h/2**(c%lvly+1)
                   !elseif(mod(c%num,m)==0)then    !右和上边界不用讨论,直接其他点通用
                    else
                    cn=>down_neighbor(c)
                    cnn=>left_neighbor(c)
                    if(associated(cn%son1))then
                     call corner(cn)
                    end if
                    if(associated(cnn%son1))then
                     call corner(cnn)
                    end if
    
                    c%cnode(1)=cn%cnode(4)
                    c%cnode(2)=cn%cnode(3)
                    c%cnode(4)=cnn%cnode(3)
    
                    nNode=nNode+1
                    c%cnode(3)=nNode
                    Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
                    Node(nNode)%y=c%center%y+h/2**(c%lvly+1)
                  end if
!子单元位置1情况
elseif(c%sp==1)then
   cn=>down_neighbor(c)
   cnn=>left_neighbor(c) 
   cnnn=>right_neighbor(c)
   !无左邻居的单元
	if(.not.associated(cnn))then
	
            c%cnode(1)=cn%cnode(4)
			
            if(c%lvl==cn%lvl)then
            c%cnode(2)=cn%cnode(3)
			
            nNode=nNode+1
            c%cnode(3)=nNode
            Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
            Node(nNode)%y=c%center%y+h/2**(c%lvly+1)	   
			
            nNode=nNode+1
            c%cnode(4)=nNode
            Node(nNode)%x=c%center%x-h/2**(c%lvlx+1)
            Node(nNode)%y=c%center%y+h/2**(c%lvly+1)	  
            else
            nNode=nNode+1
            c%cnode(2)=nNode
            Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
            Node(nNode)%y=c%center%y-h/2**(c%lvly+1)
			
            nNode=nNode+1
            c%cnode(3)=nNode
            Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
            Node(nNode)%y=c%center%y+h/2**(c%lvly+1)	   
			
            nNode=nNode+1
            c%cnode(4)=nNode
            Node(nNode)%x=c%center%x-h/2**(c%lvlx+1)
            Node(nNode)%y=c%center%y+h/2**(c%lvly+1)	
            endif
	else   
	   
        if(associated(cn))then
        if(associated(cn%son1))then    !!若cn存在子节点，则通过子节点定义cn的顶点
        call corner(cn)
        end if
        endif  
        if(associated(cnn))then
        if(associated(cnn%son1))then
        call corner(cnn)
        end if
        endif
        if(associated(cnnn))then
        if(associated(cnnn%son1))then
        call corner(cnnn)
        end if
        endif
   
!顶点1_____________________________________________
   
    nNode=nNode+1
    c%cnode(1)=nNode 
    Node(nNode)%x=c%center%x-h/2**(c%lvlx+1)
    Node(nNode)%y=c%center%y-h/2**(c%lvly+1)

	if(cnn%cnode(2)>0)then
	   if(node(nNode)%x==node(cnn%cnode(2))%x.and.node(nNode)%y==node(cnn%cnode(2))%y)then
        Node(nNode)%x=-1
        Node(nNode)%y=-1
	    nNode=nNode-1
	    c%cnode(1)=cnn%cnode(2)
		goto 22
		endif
	endif
	if(associated(cn))then
	    if(cn%cnode(4)>0)then
	        if(node(nNode)%x==node(cn%cnode(4))%x.and.node(nNode)%y==node(cn%cnode(4))%y)then
            Node(nNode)%x=-1
            Node(nNode)%y=-1
	        nNode=nNode-1
	        c%cnode(1)=cn%cnode(4)
			endif
		endif
	endif
22  continue

!!Node2
    nNode=nNode+1
    c%cnode(2)=nNode
    Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
    Node(nNode)%y=c%center%y-h/2**(c%lvly+1)
	if(associated(cn))then
	if(cn%cnode(3)>0)then
	if(node(nNode)%x==node(cn%cnode(3))%x.and.node(nNode)%y==node(cn%cnode(3))%y)then
    Node(nNode)%x=-1
    Node(nNode)%y=-1
	nNode=nNode-1
	c%cnode(2)=cn%cnode(3)
	endif
	endif
	!elseif(.not.associated(cn))then
	!cn=>down_neighbor(c)
	endif
	
    nNode=nNode+1          !顶点3
    c%cnode(3)=nNode
    Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
    Node(nNode)%y=c%center%y+h/2**(c%lvly+1) 
	if(cnnn%cnode(4)>0)then
	if(node(nNode)%x==node(cnnn%cnode(4))%x.and.node(nNode)%y==node(cnnn%cnode(4))%y)then
    Node(nNode)%x=-1
    Node(nNode)%y=-1
	nNode=nNode-1
	c%cnode(3)=cnnn%cnode(4)
	endif
	endif

		   
    nNode=nNode+1
    c%cnode(4)=nNode 
    Node(nNode)%x=c%center%x-h/2**(c%lvlx+1)
    Node(nNode)%y=c%center%y+h/2**(c%lvly+1)
	if(cnn%cnode(3)>0)then
	if(node(nNode)%x==node(cnn%cnode(3))%x.and.node(nNode)%y==node(cnn%cnode(3))%y)then
    Node(nNode)%x=-1
    Node(nNode)%y=-1
	nNode=nNode-1
	c%cnode(4)=cnn%cnode(3)
	endif
	endif

endif

!子单元2情况
elseif(c%sp==2)then
   cn=>down_neighbor(c)
   cnn=>left_neighbor(c)
   cnnn=>right_neighbor(c)

  if(associated(cn))then 
    if(associated(cn%son1))then
     call corner(cn)
    end if
  endif  
  if(associated(cnn))then 
   if(associated(cnn%son1))then
     call corner(cnn)
   end if
   endif
   if(associated(cnnn))then
     if(associated(cnnn%son1))then
      call corner(cnnn)
     end if
   endif  
	
	!!Node1
	if(c%spl==0)then
	c%cnode(1)=cnn%cnode(2)
	elseif(c%spl==1)then
	c%cnode(1)=c%father%son1%cnode(4)
	elseif(c%spl==2)then
	c%cnode(1)=c%father%son1%cnode(2)
	endif
	
	!!Node2
    nNode=nNode+1          
    c%cnode(2)=nNode
    Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
    Node(nNode)%y=c%center%y-h/2**(c%lvly+1) 
    if(associated(cn).and.cn%cnode(3)>0)then
	if(node(nNode)%x==node(cn%cnode(3))%x.and.node(nNode)%y==node(cn%cnode(3))%y)then
    Node(nNode)%x=-1
    Node(nNode)%y=-1
	nNode=nNode-1
	c%cnode(2)=cn%cnode(3)
	elseif(temp==1)then      !temp啥作用
	nNode=nNode-1
	c%cnode(2)=cn%cnode(3)
	temp=0
	endif
	endif
	
	!!Node3
    nNode=nNode+1          
    c%cnode(3)=nNode
    Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
    Node(nNode)%y=c%center%y+h/2**(c%lvly+1) 
	
	if(associated(cnnn))then
    if(cnnn%cnode(4)>0)then
	if(node(nNode)%x==node(cnnn%cnode(4))%x.and.node(nNode)%y==node(cnnn%cnode(4))%y)then
    Node(nNode)%x=-1
    Node(nNode)%y=-1
	nNode=nNode-1
	c%cnode(3)=cnnn%cnode(4)
	endif
	endif
	endif

	
	!!Node4
    nNode=nNode+1
    c%cnode(4)=nNode 
	Node(nNode)%x=c%center%x-h/2**(c%lvlx+1)
	Node(nNode)%y=c%center%y+h/2**(c%lvly+1)
	if(associated(cnn))then
    if(cnn%cnode(3)>0)then
	if(node(nNode)%x==node(cnn%cnode(3))%x.and.node(nNode)%y==node(cnn%cnode(3))%y)then
    Node(nNode)%x=-1
    Node(nNode)%y=-1
	nNode=nNode-1
	c%cnode(4)=cnn%cnode(3)
    endif
	endif
	endif
	
!子单元3   
elseif(c%sp==3)then 
  cn=>down_neighbor(c)
  cnnn=>right_neighbor(c)

   if(associated(cn%son1))then
     call corner(cn)
   end if
   if(associated(cnnn))then
     if(associated(cnnn%son1))then
      call corner(cnnn)
     end if
   endif  
      
  c%cnode(1)=cn%cnode(4)
  c%cnode(2)=cn%cnode(3)

  nNode=nNode+1
  c%cnode(3)=nNode
  Node(nNode)%x=c%center%x+h/2**(c%lvlx+1)
  Node(nNode)%y=c%center%y+h/2**(c%lvly+1)
  
	if(associated(cnnn).and.cnnn%cnode(4)>0)then
	if(node(nNode)%x==node(cnnn%cnode(4))%x.and.node(nNode)%y==node(cnnn%cnode(4))%y)then
    Node(nNode)%x=-1
    Node(nNode)%y=-1
	nNode=nNode-1
	c%cnode(3)=cnnn%cnode(4)
	endif
	endif

	
  nNode=nNode+1
  c%cnode(4)=nNode 
  Node(nNode)%x=c%center%x-h/2**(c%lvlx+1)
  Node(nNode)%y=c%center%y+h/2**(c%lvly+1) 
  
	!if(cnnn%cnode(4)>0)then
	!if(node(nNode)%x==node(cnnn%cnode(4))%x.and.node(nNode)%y==node(cnnn%cnode(4))%y)then
 !   Node(nNode)%x=-1
 !   Node(nNode)%y=-1
	!nNode=nNode-1
	!c%cnode(3)=cnnn%cnode(4)
	!endif
	!endif

  
!-------------------------------------------
  !子单元4   
elseif(c%sp==4)then 
    cn=>down_neighbor(c)
    cnn=>left_neighbor(c)
    cnnn=>right_neighbor(c)

	   if(.not.associated(cnn))then
            if(associated(cn%son1))then
            call corner(cn)
            end if
            if(associated(cnnn%son1))then
            call corner(cnnn)
            end if
	           c%cnode(1)=cn%cnode(4)
	           c%cnode(2)=cn%cnode(3)
	           c%cnode(3)=cnnn%cnode(4)
	   
			nNode=nNode+1
            c%cnode(4)=nNode
            Node(nNode)%x=c%center%x-h/2**(c%lvlx+1)
            Node(nNode)%y=c%center%y+h/2**(c%lvly+1)	  
	   else
                if(associated(cn%son1))then
                call corner(cn)
                end if
                !if(associated(cnn))then
                if(associated(cnn%son1))then
                call corner(cnn)
                end if
                !endif 
                if(associated(cnnn%son1))then
                call corner(cnnn)
                end if
      
  c%cnode(1)=cn%cnode(4)
  c%cnode(2)=cn%cnode(3)
  c%cnode(3)=cnnn%cnode(4)


  nNode=nNode+1
  c%cnode(4)=nNode 
  Node(nNode)%x=c%center%x-h/2**(c%lvlx+1)
  Node(nNode)%y=c%center%y+h/2**(c%lvly+1)
  
	if(cnn%cnode(3)>0)then
	if(node(nNode)%x==node(cnn%cnode(3))%x.and.node(nNode)%y==node(cnn%cnode(3))%y)then
    Node(nNode)%x=-1
    Node(nNode)%y=-1
	nNode=nNode-1
	c%cnode(4)=cnn%cnode(3)
    endif
	endif

  
  endif
end if  
end if  
!以上是邻居方法,以下为遍历方法   
!if(statuePrint==0)then
!return
!else
!write(*,'(I5,2X,A4,I3,2X,A4,I1,2X,A4,I1,2X,A3,I1,2X,A5,I4,1X,I4,1X,I4,1X,I4)') nElment,"num:", c%num, "lvl:", c%lvl, "spl:", c%spl, "sp:", c%sp, "node:", c%cnode(1), c%cnode(2), c%cnode(3), c%cnode(4)
!!if(c%lvl>=lvltemp)then
!call location(c)
!!endif
!endif

!if(c%cnode(1)<=0.or.c%cnode(2)<=0.or.c%cnode(3)<=0.or.c%cnode(4)<=0)then
!    if(c%cnode(1)<=0.and.c%cnode(2)<=0.and.c%cnode(3)<=0.and.c%cnode(4)<=0)then
!	return
!	endif
!	error=error+1
!	print*,"error___________________________________________________________________________"
!endif
return
end subroutine addnode
!对于非叶子单元，通过子单元得到其单元的顶点坐标，递归
recursive subroutine corner(c)
implicit none
type(gridtyp),pointer::c

   if(associated(c%son1%son1))then
     call corner(c%son1)
   end if
   if(associated(c%son2%son1))then
     call corner(c%son2)
   end if
   if(associated(c%son3))then 
   if(associated(c%son3%son1))then  
     call corner(c%son3)
   end if
   endif
   if(associated(c%son4))then 
   if(associated(c%son4%son1))then 
     call corner(c%son4)
   end if 
   endif
if(c%son1%spl==0)then
  c%cnode(1)=c%son1%cnode(1)
  c%cnode(2)=c%son2%cnode(2) 
  c%cnode(3)=c%son3%cnode(3)
  c%cnode(4)=c%son4%cnode(4)
elseif(c%son1%spl==1)then
  c%cnode(1)=c%son1%cnode(1)
  c%cnode(2)=c%son1%cnode(2) 
  c%cnode(3)=c%son2%cnode(3)
  c%cnode(4)=c%son2%cnode(4)
elseif(c%son1%spl==2)then
  c%cnode(1)=c%son1%cnode(1)
  c%cnode(2)=c%son2%cnode(2) 
  c%cnode(3)=c%son2%cnode(3)
  c%cnode(4)=c%son1%cnode(4)
else
    call stopall("corner")
endif
end subroutine corner
!_________________________________________
subroutine cellinf
integer i,j
type(gridtyp),pointer::c

nElment=1
cell(1)%nelment=nelment
do i=2,total
   c=>cell(i)
   call cellinfsub(c)
end do
endsubroutine cellinf

recursive subroutine cellinfsub(c)
type(gridtyp),pointer::c

if(associated(c%son1))then
    if(associated(c%son3))then
    call cellinfsub(c%son1)
    call cellinfsub(c%son2)
    call cellinfsub(c%son3)
    call cellinfsub(c%son4)
    else
    call cellinfsub(c%son1)
    call cellinfsub(c%son2)
    endif
    return
endif
    nElment=nElment+1
    c%nelment=nelment
endsubroutine cellinfsub
!_________________________________________
end module outfile
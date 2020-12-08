!网格光顺,针对点和洞优化，是在解自适应完成后
module GS
use neighbor 
use frame
use precisionMod
use typeMod
implicit none

contains

subroutine guangshun
implicit none
integer i
type(gridtyp),pointer::c
!首先消除洞区域，即加密
do i=1,total
  c=>cell(i)
  if(c%cross==0)then
  call gs_hole(c)  
  end if 
end do
!再处理点区域，即粗化
do i=1,total
  c=>cell(i)
  if(c%cross==0.and.associated(c%son1))then
  call gs_point(c)  
  end if 
end do

return
end subroutine guangshun
!
recursive subroutine gs_hole(c)
implicit none
type(gridtyp),pointer::c,cleft,cr,cu,cd

if(associated(c%son1))then
     if (associated(c%son3)) then
  call gs_hole(c%son1)
  call gs_hole(c%son2)
  call gs_hole(c%son3)
  call gs_hole(c%son4)
     else
  call gs_hole(c%son1)
  call gs_hole(c%son2)  
     endif  
else 
!四个邻居指针赋值，考虑不存在的情况
  if(associated(left_neighbor(c)))then
  cleft=>left_neighbor(c)
  else
  cleft=>null()
  end if
  if(associated(right_neighbor(c)))then
  cr=>right_neighbor(c)
  else
  cr=>null()
  end if
  if(associated(up_neighbor(c)))then
  cu=>up_neighbor(c)
  else
  cu=>null()
  end if
  if(associated(down_neighbor(c)))then
  cd=>down_neighbor(c)
  else
  cd=>null()
  end if
!水平相邻的比较为例，左右邻居都存在，且左右邻居都存在子单元，则加密c，竖直类似
  if((associated(cleft).and.associated(cr).and.associated(cleft%son1).and.associated(cr%son1)).or.&
  &(associated(cu).and.associated(cd).and.associated(cu%son1).and.associated(cd%son1)))then
   call gridmodify_ref2(c)
  end if 
  
end if

return
end subroutine gs_hole 
!!!
recursive subroutine gs_point(c)
implicit none
type(gridtyp),pointer::c,cleft,cr,cu,cd

if(associated(c%son1).and.(.not.associated(c%son1%son1).and.&
&.not.associated(c%son2%son1).and..not.associated(c%son3%son1).and..not.associated(c%son4%son1)))then
!四个邻居指针赋值，考虑不存在的情况
  if(associated(left_neighbor(c)))then
  cleft=>left_neighbor(c)
  else
  cleft=>null()
  end if
  if(associated(right_neighbor(c)))then
  cr=>right_neighbor(c)
  else
  cr=>null()
  end if
  if(associated(up_neighbor(c)))then
  cu=>up_neighbor(c)
  else
  cu=>null()
  end if
  if(associated(down_neighbor(c)))then
  cd=>down_neighbor(c)
  else
  cd=>null()
  end if
!水平相邻的比较为例，左右邻居都存在，且左右邻居都为叶子单元，则粗化c，竖直类似
  if((associated(cleft).and.associated(cr).and.associated(cleft%son1).and.associated(cr%son1)).or.&
  &(associated(cu).and.associated(cd).and.associated(cu%son1).and.associated(cd%son1)))then
   
   call interposeA(c,c%rou,c%u,c%v,c%T,c%p) 
   deallocate(c%son1,c%son2,c%son3,c%son4)
!再置空子单元
   nullify(c%son1,c%son2,c%son3,c%son4)
  end if 

else if(associated(c%son1))then
  if(associated(c%son1%son1))then
    call gs_point(c%son1)
  endif
  if(associated(c%son2%son1))then
    call gs_point(c%son2)
  endif
  if(associated(c%son3%son1))then
    call gs_point(c%son3)
  endif
  if(associated(c%son4%son1))then
    call gs_point(c%son4)  
  endif 
end if

return
end subroutine gs_point 

end module GS
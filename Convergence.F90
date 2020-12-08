module conver
use mesh
use precisionMod
use typeMod
implicit none

contains

subroutine conver_criterion
implicit none
type(gridtyp),pointer::c
integer count,i

nconver=-1
count=0
do i=1,total
c=>cell(i)
if(c%cross/=2)then
 call criterion(c,count)
 if(count>0)then
   nconver=0
   exit
 end if 
endif
end do
if(count==0)then
  nconver=1  !收敛  !迭代程序里通过nconver判断是否收敛
endif

return
end subroutine conver_criterion
!
recursive subroutine criterion(c,count)
implicit none
type(gridtyp),pointer::c
integer count,i
real(pre)::eps
eps=0.001

if(.not.associated(c%son1).and.abs(c%ut-c%u)>eps)then
  count=count+1
elseif(associated(c%son1))then
    if(associated(c%son3))then
 call criterion(c%son1,count) 
 call criterion(c%son2,count)
 call criterion(c%son3,count)
 call criterion(c%son4,count)
    else
  call criterion(c%son1,count) 
  call criterion(c%son2,count) 
  endif
endif
return
end subroutine criterion

end module
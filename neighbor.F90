module neighbor
use mesh
use precisionMod
use typeMod
implicit none

    contains
!初始单元邻居
!up

function initial_up_neighbor(c)
  implicit none
  type(gridtyp),pointer::initial_up_neighbor,c
  integer i
  
   if(c%num>(n-1)*m)then
     initial_up_neighbor=>null()
    else
     i=c%num+m
     initial_up_neighbor=>cell(i)  
   end if  
return  
end function initial_up_neighbor
!初始邻居down
function initial_down_neighbor(c)
  implicit none
  type(gridtyp),pointer::initial_down_neighbor,c
  integer i
  
   if(c%num<=m)then
     initial_down_neighbor=>null()
    else
     i=c%num-m
     initial_down_neighbor=>cell(i)  
   end if  
return  
end function initial_down_neighbor
!初始邻居left
function initial_left_neighbor(c)
  implicit none
  type(gridtyp),pointer::initial_left_neighbor,c
  integer i
  
   if(mod(c%num-1,m)==0)then
     initial_left_neighbor=>null()
    else
     i=c%num-1
     initial_left_neighbor=>cell(i)
   end if  
return  
end function initial_left_neighbor
!初始邻居right
function initial_right_neighbor(c)
  implicit none
  type(gridtyp),pointer::initial_right_neighbor,c
  integer q
  
  if(mod(c%num,m)==0)then
     initial_right_neighbor=>null()
  else
     q=c%num+1 
     initial_right_neighbor=>cell(q)  
    end if  
return  
end function initial_right_neighbor

!邻居查询上up （下down左left右right）
recursive function up_neighbor(cell)    !返回的邻居是type类型的,celln是初始网格单元的
  implicit none    
  type(gridtyp),pointer::cell,cellf,cellfn
  type(gridtyp),pointer::up_neighbor
  
	up_neighbor=>null()
    cellf=>cell%father
    
	if(.not.associated(cell%father))then
	  up_neighbor=>initial_up_neighbor(cell)
	return
    endif
    
	cellfn=>up_neighbor(cellf)
	
	if(cell%spl==0)then
		if(cell%sp==1)then
			up_neighbor=>cellf%son4
		elseif(cell%sp==2)then
			up_neighbor=>cellf%son3
        elseif(cell%sp==3)then
            if(.not.associated(cellfn))then
			up_neighbor=>null()
			goto 21
			endif
			if(.not.associated(cellfn%son1))then
			up_neighbor=>cellfn
			goto 21
			endif
			if(cellfn%typ==0)then
			up_neighbor=>cellfn%son2
			elseif(cellfn%typ==1)then
			up_neighbor=>cellfn%son1
			elseif(cellfn%typ==2)then
			up_neighbor=>cellfn%son2
			endif
		elseif(cell%sp==4)then
            if(.not.associated(cellfn))then
			up_neighbor=>null()
			goto 21
			endif
			if(.not.associated(cellfn%son1))then
			up_neighbor=>cellfn
			goto 21
			endif
			if(cellfn%typ==0)then
			up_neighbor=>cellfn%son1
			elseif(cellfn%typ==1)then
			up_neighbor=>cellfn%son1
			elseif(cellfn%typ==2)then
			up_neighbor=>cellfn%son1
			endif
		endif
	elseif(cell%spl==1)then
		if(cell%sp==1)then
			up_neighbor=>cellf%son2
		elseif(cell%sp==2)then
            if(.not.associated(cellfn))then
			up_neighbor=>null()
			goto 21
			endif
			if(.not.associated(cellfn%son1))then
			up_neighbor=>cellfn
			goto 21
			endif
			if(cellfn%typ==0)then
			up_neighbor=>cellfn
			elseif(cellfn%typ==1)then
			up_neighbor=>cellfn%son1
			elseif(cellfn%typ==2)then
			up_neighbor=>cellfn
			endif
		endif
	elseif(cell%spl==2)then
		if(cell%sp==1)then
            if(.not.associated(cellfn))then
			up_neighbor=>null()
			goto 21
			endif
			if(.not.associated(cellfn%son1))then
			up_neighbor=>cellfn
			goto 21
			endif
			if(cellfn%typ==0)then
			up_neighbor=>cellfn%son1
			elseif(cellfn%typ==1)then
			up_neighbor=>cellfn%son1
			elseif(cellfn%typ==2)then
			up_neighbor=>cellfn%son1
			endif
		elseif(cell%sp==2)then
            if(.not.associated(cellfn))then
			up_neighbor=>null()
			goto 21
			endif
			if(.not.associated(cellfn%son1))then
			up_neighbor=>cellfn
			goto 21
			endif
			if(cellfn%typ==0)then
			up_neighbor=>cellfn%son2
			elseif(cellfn%typ==1)then
			up_neighbor=>cellfn%son1
			elseif(cellfn%typ==2)then
			up_neighbor=>cellfn%son2
			endif
		endif
	endif
21  return
end function up_neighbor

!下邻居down
recursive function down_neighbor(cell)    
  implicit none    
  type(gridtyp),pointer::cell,cellf,cellfn
  type(gridtyp),pointer::down_neighbor
  
	down_neighbor=>null()
    cellf=>cell%father
    
	if(.not.associated(cell%father))then
	   down_neighbor=>initial_down_neighbor(cell)
	return
    endif
    
	cellfn=>down_neighbor(cellf)
	
	if(cell%spl==0)then
		if(cell%sp==1)then
            if(.not.associated(cellfn))then
			down_neighbor=>null()
			goto 23
			endif
			if(.not.associated(cellfn%son1))then
			down_neighbor=>cellfn
			goto 23
			endif
			if(cellfn%typ==0)then
			down_neighbor=>cellfn%son4
			elseif(cellfn%typ==1)then
			down_neighbor=>cellfn%son2
			elseif(cellfn%typ==2)then
			down_neighbor=>cellfn%son1
			endif
		elseif(cell%sp==2)then
            if(.not.associated(cellfn))then
			down_neighbor=>null()
			goto 23
			endif
			if(.not.associated(cellfn%son1))then
			down_neighbor=>cellfn
			goto 23
			endif
			if(cellfn%typ==0)then
			down_neighbor=>cellfn%son3
			elseif(cellfn%typ==1)then
			down_neighbor=>cellfn%son2
			elseif(cellfn%typ==2)then
			down_neighbor=>cellfn%son2
			endif
		elseif(cell%sp==3)then
			down_neighbor=>cellf%son2
		elseif(cell%sp==4)then
			down_neighbor=>cellf%son1
		endif
	elseif(cell%spl==1)then
		if(cell%sp==1)then
            if(.not.associated(cellfn))then
			down_neighbor=>null()
			goto 23
			endif
			if(.not.associated(cellfn%son1))then
			down_neighbor=>cellfn
			goto 23
			endif
			if(cellfn%typ==0)then
			down_neighbor=>cellfn
			elseif(cellfn%typ==1)then
			down_neighbor=>cellfn%son2
			elseif(cellfn%typ==2)then
			down_neighbor=>cellfn
			endif
		elseif(cell%sp==2)then
			down_neighbor=>cellf%son1
		endif
	elseif(cell%spl==2)then
		if(cell%sp==1)then
            if(.not.associated(cellfn))then
			down_neighbor=>null()
			goto 23
			endif
			if(.not.associated(cellfn%son1))then
			down_neighbor=>cellfn
			goto 23
			endif
			if(cellfn%typ==0)then
			down_neighbor=>cellfn%son4
			elseif(cellfn%typ==1)then
			down_neighbor=>cellfn%son2
			elseif(cellfn%typ==2)then
			down_neighbor=>cellfn%son1
			endif
		elseif(cell%sp==2)then
            if(.not.associated(cellfn))then
			down_neighbor=>null()
			goto 23
			endif
			if(.not.associated(cellfn%son1))then
			down_neighbor=>cellfn
			goto 23
			endif
			if(cellfn%typ==0)then
			down_neighbor=>cellfn%son3
			elseif(cellfn%typ==1)then
			down_neighbor=>cellfn%son2
			elseif(cellfn%typ==2)then
			down_neighbor=>cellfn%son2
			endif
		endif
	endif
23  return
end function down_neighbor

!左邻居left
recursive function left_neighbor(cell)
  implicit none    
  type(gridtyp),pointer::cell,cellf,cellfn
  type(gridtyp),pointer::left_neighbor
  
	left_neighbor=>null()
    cellf=>cell%father
    
	if(.not.associated(cellf))then
	 left_neighbor=>initial_left_neighbor(cell)
	return
    endif
    
	cellfn=>left_neighbor(cellf)
	
	if(cell%spl==0)then
		if(cell%sp==1)then
            if(.not.associated(cellfn))then
			left_neighbor=>null()
			goto 25
			endif
			if(.not.associated(cellfn%son1))then
			left_neighbor=>cellfn
			goto 25
			endif
			if(cellfn%typ==0)then
			left_neighbor=>cellfn%son2
			elseif(cellfn%typ==1)then
			left_neighbor=>cellfn%son1
			elseif(cellfn%typ==2)then
			left_neighbor=>cellfn%son2
			endif
		elseif(cell%sp==2)then
			left_neighbor=>cellf%son1
		elseif(cell%sp==3)then
			left_neighbor=>cellf%son4
		elseif(cell%sp==4)then
            if(.not.associated(cellfn))then
			left_neighbor=>null()
			goto 25
			endif
			if(.not.associated(cellfn%son1))then
			left_neighbor=>cellfn
			goto 25
			endif
			if(cellfn%typ==0)then
			left_neighbor=>cellfn%son3
			elseif(cellfn%typ==1)then
			left_neighbor=>cellfn%son2
			elseif(cellfn%typ==2)then
			left_neighbor=>cellfn%son2
			endif
		endif
	elseif(cell%spl==1)then
		if(cell%sp==1)then
            if(.not.associated(cellfn))then
			left_neighbor=>null()
			goto 25
			endif
			if(.not.associated(cellfn%son1))then
			left_neighbor=>cellfn
			goto 25
			endif
			if(cellfn%typ==0)then
			left_neighbor=>cellfn%son2
			elseif(cellfn%typ==1)then
			left_neighbor=>cellfn%son1
			elseif(cellfn%typ==2)then
			left_neighbor=>cellfn%son2
			endif
		elseif(cell%sp==2)then
            if(.not.associated(cellfn))then
			left_neighbor=>null()
			goto 25
			endif
			if(.not.associated(cellfn%son1))then
			left_neighbor=>cellfn
			goto 25
			endif
			if(cellfn%typ==0)then
			left_neighbor=>cellfn%son3
			elseif(cellfn%typ==1)then
			left_neighbor=>cellfn%son2
			elseif(cellfn%typ==2)then
			left_neighbor=>cellfn%son2
			endif
		endif
	elseif(cell%spl==2)then
		if(cell%sp==1)then
            if(.not.associated(cellfn))then
			left_neighbor=>null()
			goto 25
			endif
			if(.not.associated(cellfn%son1))then
			left_neighbor=>cellfn
			goto 25
			endif
			if(cellfn%typ==0)then
			left_neighbor=>cellfn
			elseif(cellfn%typ==1)then
			left_neighbor=>cellfn
			elseif(cellfn%typ==2)then
			left_neighbor=>cellfn%son2
			endif
		elseif(cell%sp==2)then
			left_neighbor=>cellf%son1
		endif
	endif
25  return
end function left_neighbor

!右邻居right
recursive function right_neighbor(cell)    
  implicit none    
  type(gridtyp),pointer::cell,cellf,cellfn
  type(gridtyp),pointer::right_neighbor
  
	right_neighbor=>null()
    cellf=>cell%father
	
    if(.not.associated(cell%father))then
	right_neighbor=>initial_right_neighbor(cell)
	return
    endif
	
    cellfn=>right_neighbor(cellf)
	
	if(cell%spl==0)then
		if(cell%sp==1)then
			right_neighbor=>cellf%son2
		elseif(cell%sp==2)then
            if(.not.associated(cellfn))then
			right_neighbor=>null()
			goto 27
			endif
			if(.not.associated(cellfn%son1))then
			right_neighbor=>cellfn
			goto 27
			endif
			if(cellfn%typ==0)then
			right_neighbor=>cellfn%son1
			elseif(cellfn%typ==1)then
			right_neighbor=>cellfn%son1
			elseif(cellfn%typ==2)then
			right_neighbor=>cellfn%son1
			endif
		elseif(cell%sp==3)then
            if(.not.associated(cellfn))then
			right_neighbor=>null()
			goto 27
			endif
			if(.not.associated(cellfn%son1))then
			right_neighbor=>cellfn
			goto 27
			endif
			if(cellfn%typ==0)then
			right_neighbor=>cellfn%son4
			elseif(cellfn%typ==1)then
			right_neighbor=>cellfn%son2
			elseif(cellfn%typ==2)then
			right_neighbor=>cellfn%son1
			endif
		elseif(cell%sp==4)then
			right_neighbor=>cellf%son3
		endif
	elseif(cell%spl==1)then
		if(cell%sp==1)then
            if(.not.associated(cellfn))then
			right_neighbor=>null()
			goto 27
			endif
			if(.not.associated(cellfn%son1))then
			right_neighbor=>cellfn
			goto 27
			endif
			if(cellfn%typ==0)then
			right_neighbor=>cellfn%son1
			elseif(cellfn%typ==1)then
			right_neighbor=>cellfn%son1
			elseif(cellfn%typ==2)then
			right_neighbor=>cellfn%son1
			endif
		elseif(cell%sp==2)then
            if(.not.associated(cellfn))then
			right_neighbor=>null()
			goto 27
			endif
			if(.not.associated(cellfn%son1))then
			right_neighbor=>cellfn
			goto 27
			endif
			if(cellfn%typ==0)then
			right_neighbor=>cellfn%son4
			elseif(cellfn%typ==1)then
			right_neighbor=>cellfn%son2
			elseif(cellfn%typ==2)then
			right_neighbor=>cellfn%son1
			endif
		endif
	elseif(cell%spl==2)then
		if(cell%sp==1)then
			right_neighbor=>cellf%son2
		elseif(cell%sp==2)then
            if(.not.associated(cellfn))then
			right_neighbor=>null()
			goto 27
			endif
			if(.not.associated(cellfn%son1))then
			right_neighbor=>cellfn
			goto 27
			endif
			if(cellfn%typ==0)then
			right_neighbor=>cellfn
			elseif(cellfn%typ==1)then
			right_neighbor=>cellfn
			elseif(cellfn%typ==2)then
			right_neighbor=>cellfn%son1
			endif
		endif
	endif
27  return
endfunction
!
!针对网格层次差大于1，修正网格
recursive subroutine gridmodify(c)
implicit none
type(gridtyp),pointer::c,ct1,ct2,ct3,ct4

if(associated(c%son1))then
  call gridmodify(c%son1)
  call gridmodify(c%son2)
  call gridmodify(c%son3)
  call gridmodify(c%son4)
else
  ct1=>up_neighbor(c)
  ct2=>down_neighbor(c)
  ct3=>left_neighbor(c)
  ct4=>right_neighbor(c)
  if(associated(ct1))then
    if((c%lvl-ct1%lvl)>1)then
	call newson(ct1)
	call nullifyCell(ct1%son1)
	call nullifyCell(ct1%son2)
	call nullifyCell(ct1%son3)
	call nullifyCell(ct1%son4)
    end if
  end if
  
  if(associated(ct2))then
    if((c%lvl-ct2%lvl)>1)then
	call newson(ct2)
	call nullifyCell(ct2%son1)
	call nullifyCell(ct2%son2)
	call nullifyCell(ct2%son3)
	call nullifyCell(ct2%son4)
    end if
  end if
  
  if(associated(ct3))then  
    if((c%lvl-ct3%lvl)>1)then
	call newson(ct3)
	call nullifyCell(ct3%son1)
	call nullifyCell(ct3%son2)
	call nullifyCell(ct3%son3)
	call nullifyCell(ct3%son4)
    end if
  end if
  
  if(associated(ct4))then  
    if((c%lvl-ct4%lvl)>1)then
    call newson(ct4)
	call nullifyCell(ct4%son1)
	call nullifyCell(ct4%son2)
	call nullifyCell(ct4%son3)
	call nullifyCell(ct4%son4)
    end if
  end if
    
end if

return
end subroutine gridmodify
!同样针对解自适应
!针对网格层次差大于1，修正网格
recursive subroutine gridmodify2(c)
!!implicit none
type(gridtyp),pointer::c,ct1,ct2,ct3,ct4
!!
!!if(associated(c%son1))then
!!  call gridmodify2(c%son1)
!!  call gridmodify2(c%son2)
!!  call gridmodify2(c%son3)
!!  call gridmodify2(c%son4)
!!else
!!  ct1=>up_neighbor(c)
!!  ct2=>down_neighbor(c)
!!  ct3=>left_neighbor(c)
!!  ct4=>right_neighbor(c)
!!  if(associated(ct1))then
!!if((c%lvl-ct1%lvl)>1)then
!!
!!call gridmodify_ref2(ct1)
!!end if
!!  end if
!!  
!!  if(associated(ct2))then
!!    if((c%lvl-ct2%lvl)>1)then
!!
!!     call gridmodify_ref2(ct2)
!!    end if
!!  end if
!!  
!!  if(associated(ct3))then  
!!    if((c%lvl-ct3%lvl)>1)then
!!  
!!     call gridmodify_ref2(ct3)
!!    end if
!!  end if
!!  
!!  if(associated(ct4))then  
!!    if((c%lvl-ct4%lvl)>1)then
!!
!!     call gridmodify_ref2(ct4)
!!    end if
!!  end if   
!!end if
!!return
end subroutine gridmodify2
!加密需要加密的邻居单元
!!!subroutine gridmodify_ref(cell)
!!!implicit none
!!!type(gridtyp),pointer::cell
!!!  allocate(cell%son1,cell%son2,cell%son3,cell%son4)   !每个指针需要给予地址allocate
!!!      cell%son1%num=cell%num
!!!      cell%son1%lvl=cell%lvl+1
!!!      cell%son1%center%x=cell%center%x-h/(2**(cell%lvl+2))
!!!      cell%son1%center%y=cell%center%y-h/(2**(cell%lvl+2))
!!!      cell%son1%sort=gridsort(cell%son1%center%x,cell%son1%center%y)
!!!      cell%son1%cur=gridcur(cell%son1%lvl,cell%son1%cross,cell%son1%center%x,cell%son1%center%y)
!!!      cell%son1%father=>cell
!!!      cell%son1%sp=1
!!!	  cell%son1%spl=0
!!!      cell%son1%cross=gridcross(cell%son1)
!!!      nullify(cell%son1%son1,cell%son1%son2,cell%son1%son3,cell%son1%son4)   
!!!      !对每个子单元的四个子单元置空，要不然四个指针无指向，内存错误。
!!!      cell%son2%num=cell%num
!!!      cell%son2%lvl=cell%lvl+1
!!!      cell%son2%center%x=cell%center%x+h/(2**(cell%lvl+2))
!!!      cell%son2%center%y=cell%center%y-h/(2**(cell%lvl+2))
!!!      cell%son2%sort=gridsort(cell%son2%center%x,cell%son2%center%y)
!!!      cell%son2%cur=gridcur(cell%son2%lvl,cell%son2%cross,cell%son2%center%x,cell%son2%center%y)
!!!      cell%son2%father=>cell
!!!      cell%son2%sp=2
!!!	  cell%son2%spl=0
!!!      cell%son2%cross=gridcross(cell%son2)
!!!      nullify(cell%son2%son1,cell%son2%son2,cell%son2%son3,cell%son2%son4)
!!!      
!!!      cell%son3%num=cell%num
!!!      cell%son3%lvl=cell%lvl+1
!!!      cell%son3%center%x=cell%center%x+h/(2**(cell%lvl+2))
!!!      cell%son3%center%y=cell%center%y+h/(2**(cell%lvl+2))
!!!      cell%son3%sort=gridsort(cell%son3%center%x,cell%son3%center%y)
!!!      cell%son3%cur=gridcur(cell%son3%lvl,cell%son3%cross,cell%son3%center%x,cell%son3%center%y)
!!!      cell%son3%father=>cell
!!!      cell%son3%sp=3
!!!	  cell%son3%spl=0
!!!      cell%son3%cross=gridcross(cell%son3)
!!!      nullify(cell%son3%son1,cell%son3%son2,cell%son3%son3,cell%son3%son4)
!!!      
!!!      cell%son4%num=cell%num
!!!      cell%son4%lvl=cell%lvl+1
!!!      cell%son4%center%x=cell%center%x-h/(2**(cell%lvl+2))
!!!      cell%son4%center%y=cell%center%y+h/(2**(cell%lvl+2))
!!!      cell%son4%sort=gridsort(cell%son4%center%x,cell%son4%center%y)
!!!      cell%son4%cur=gridcur(cell%son4%lvl,cell%son4%cross,cell%son4%center%x,cell%son4%center%y)
!!!      cell%son4%father=>cell
!!!      cell%son4%sp=4
!!!	  cell%son4%spl=0
!!!      cell%son4%cross=gridcross(cell%son4)
!!!      nullify(cell%son4%son1,cell%son4%son2,cell%son4%son3,cell%son4%son4)
!!!	  
!!!	  cell%son1%gxyxsd=0
!!!	  cell%son2%gxyxsd=0
!!!	  cell%son3%gxyxsd=0
!!!	  cell%son4%gxyxsd=0
!!!	  
!!!return
!!!end subroutine gridmodify_ref
!针对解自适应加密需要加密的邻居单元
!主要是流场解的传递，区别于外形自适应加密
subroutine gridmodify_ref2(cell)
!!implicit none
type(gridtyp),pointer::cell
!!  allocate(cell%son1,cell%son2,cell%son3,cell%son4)   !每个指针需要给予地址allocate
!!      cell%son1%num=cell%num
!!      cell%son1%lvl=cell%lvl+1
!!      cell%son1%center%x=cell%center%x-h/(2**(cell%lvl+2))
!!      cell%son1%center%y=cell%center%y-h/(2**(cell%lvl+2))
!!      cell%son1%sort=gridsort(cell%son1%center%x,cell%son1%center%y)
!!      cell%son1%cross=gridcross(cell%son1%lvl,cell%son1%center%x,cell%son1%center%y)
!!      cell%son1%cur=gridcur(cell%son1%lvl,cell%son1%cross,cell%son1%center%x,cell%son1%center%y)
!!      cell%son1%father=>cell
!!      cell%son1%sp=1
!!      !流场解传递
!!      cell%son1%u=cell%u
!!      cell%son1%v=cell%v
!!      cell%son1%rou=cell%rou
!!      cell%son1%T=cell%T
!!      cell%son1%p=cell%p
!!      cell%son1%gra=cell%gra
!!      cell%son1%secgra=cell%secgra
!!      nullify(cell%son1%son1,cell%son1%son2,cell%son1%son3,cell%son1%son4)   
!!      !对每个子单元的四个子单元置空，要不然四个指针无指向，内存错误。
!!      cell%son2%num=cell%num
!!      cell%son2%lvl=cell%lvl+1
!!      cell%son2%center%x=cell%center%x+h/(2**(cell%lvl+2))
!!      cell%son2%center%y=cell%center%y-h/(2**(cell%lvl+2))
!!      cell%son2%sort=gridsort(cell%son2%center%x,cell%son2%center%y)
!!      cell%son2%cross=gridcross(cell%son2%lvl,cell%son2%center%x,cell%son2%center%y)
!!      cell%son2%cur=gridcur(cell%son2%lvl,cell%son2%cross,cell%son2%center%x,cell%son2%center%y)
!!      cell%son2%father=>cell
!!      cell%son2%sp=2
!!      cell%son2%u=cell%u
!!      cell%son2%v=cell%v
!!      cell%son2%rou=cell%rou
!!      cell%son2%T=cell%T
!!      cell%son2%p=cell%p
!!      cell%son2%gra=cell%gra
!!      cell%son2%secgra=cell%secgra
!!      nullify(cell%son2%son1,cell%son2%son2,cell%son2%son3,cell%son2%son4)
!!      
!!      cell%son3%num=cell%num
!!      cell%son3%lvl=cell%lvl+1
!!      cell%son3%center%x=cell%center%x+h/(2**(cell%lvl+2))
!!      cell%son3%center%y=cell%center%y+h/(2**(cell%lvl+2))
!!      cell%son3%sort=gridsort(cell%son3%center%x,cell%son3%center%y)
!!      cell%son3%cross=gridcross(cell%son3%lvl,cell%son3%center%x,cell%son3%center%y)
!!      cell%son3%cur=gridcur(cell%son3%lvl,cell%son3%cross,cell%son3%center%x,cell%son3%center%y)
!!      cell%son3%father=>cell
!!      cell%son3%sp=3
!!      cell%son3%u=cell%u
!!      cell%son3%v=cell%v
!!      cell%son3%rou=cell%rou
!!      cell%son3%T=cell%T
!!      cell%son3%p=cell%p
!!      cell%son3%gra=cell%gra
!!      cell%son3%secgra=cell%secgra
!!      nullify(cell%son3%son1,cell%son3%son2,cell%son3%son3,cell%son3%son4)
!!      
!!      cell%son4%num=cell%num
!!      cell%son4%lvl=cell%lvl+1
!!      cell%son4%center%x=cell%center%x-h/(2**(cell%lvl+2))
!!      cell%son4%center%y=cell%center%y+h/(2**(cell%lvl+2))
!!      cell%son4%sort=gridsort(cell%son4%center%x,cell%son4%center%y)
!!      cell%son4%cross=gridcross(cell%son4%lvl,cell%son4%center%x,cell%son4%center%y)
!!      cell%son4%cur=gridcur(cell%son4%lvl,cell%son4%cross,cell%son4%center%x,cell%son4%center%y)
!!      cell%son4%father=>cell
!!      cell%son4%sp=4
!!      cell%son4%u=cell%u
!!      cell%son4%v=cell%v
!!      cell%son4%rou=cell%rou
!!      cell%son4%T=cell%T
!!      cell%son4%p=cell%p
!!      cell%son4%gra=cell%gra
!!      cell%son4%secgra=cell%secgra
!!      nullify(cell%son4%son1,cell%son4%son2,cell%son4%son3,cell%son4%son4)
!!return
end subroutine gridmodify_ref2


!_____________________________________________
subroutine location(c)!输出网格单元位置信息
implicit none
type(gridtyp),pointer::c,tmp
integer,allocatable::p(:)
integer::lvl,i,j
tmp=>c
lvl=tmp%lvl
allocate(p(lvl))
do i=1,lvl
j=lvl-i+1
p(j)=tmp%sp
tmp=>tmp%father
enddo
write(*,"(5x,9I2)") p
deallocate(p)   !!再释放数组空间
endsubroutine
!_____________________________________________
recursive subroutine gxyxmodify_def(c)
implicit none
type(gridtyp),pointer::c,ct1,ct2,ct3,ct4

if(associated(c%son1))then
    if(associated(c%son3))then
    call gxyxmodify_def(c%son1)
    call gxyxmodify_def(c%son2)
    call gxyxmodify_def(c%son3)
    call gxyxmodify_def(c%son4)
	else
    call gxyxmodify_def(c%son1)
    call gxyxmodify_def(c%son2) 
	endif
	return
endif

    ct1=>up_neighbor(c)
    ct2=>down_neighbor(c)
    ct3=>left_neighbor(c)
    ct4=>right_neighbor(c)

    if(associated(ct1))then
    if((c%lvl-ct1%lvl)>1.and.c%spl/=1)then
    ct1%s1=.true.
	elseif((c%lvl-ct1%lvl)==1.and.ct1%spl==1.and.c%spl/=1)then
	ct1%s1=.true.
	end if
    end if

    if(associated(ct2))then
    if((c%lvl-ct2%lvl)>1.and.c%spl/=1)then
    ct2%s2=.true.
	elseif((c%lvl-ct2%lvl)==1.and.ct2%spl==1.and.c%spl/=1)then
	ct2%s2=.true.
	end if
    end if
  
    if(associated(ct3))then
    if((c%lvl-ct3%lvl)>1.and.c%spl/=2)then
    ct3%s3=.true.
	elseif((c%lvl-ct3%lvl)==1.and.ct3%spl==2.and.c%spl/=2)then
	ct3%s3=.true.
	end if
    end if
  
    if(associated(ct4))then
    if((c%lvl-ct4%lvl)>1.and.c%spl/=2)then
    ct4%s4=.true.
	elseif((c%lvl-ct4%lvl)==1.and.ct4%spl==2.and.c%spl/=2)then
	ct4%s4=.true.
	end if
    end if
endsubroutine
 

recursive subroutine init_s(c)
implicit none
type(gridtyp),pointer::c
if(associated(c%son1))then
    if(.not.associated(c%son3))then
    call init_s(c%son1)
    call init_s(c%son2)
    else
        
    call init_s(c%son1)
    call init_s(c%son2)
    call init_s(c%son3)
    call init_s(c%son4)
	endif
	return
else
	c%s1=.false.
    c%s2=.false.
    c%s3=.false.
    c%s4=.false.
endif
endsubroutine

recursive subroutine gxyxmodify_ref(c)
implicit none
type(gridtyp),pointer::c
if(associated(c%son1))then  !!确保不会对父单元加密
    if(.not.associated(c%son3))then
    call gxyxmodify_ref(c%son1)
    call gxyxmodify_ref(c%son2)
	else
    call gxyxmodify_ref(c%son1)
    call gxyxmodify_ref(c%son2)
    call gxyxmodify_ref(c%son3)
    call gxyxmodify_ref(c%son4)
	endif
else
	if((c%s1.or.c%s2).and.(c%s3.or.c%s4))then
	c%typ=0
	call newson(c)
	!call nullifyCell(c%son1)
	!call nullifyCell(c%son2)
	!call nullifyCell(c%son3)
	!call nullifyCell(c%son4)
    elseif((c%s1.or.c%s2).and.(.not.c%s3).and.(.not.c%s4))then
	!c%typ=2
        c%typ=0
	call newson(c)
	!call nullifyCell(c%son1)
	!call nullifyCell(c%son2)
    elseif((.not.c%s1).and.(.not.c%s2).and.(c%s3.or.c%s4))then
	!c%typ=1
        c%typ=0
	call newson(c)
	!call nullifyCell(c%son1)
	!call nullifyCell(c%son2)
	!elseif((c%s1.or.c%s2).and.(.not.c%s3).and.(.not.c%s4).and.c%spl==1)then
 !   call gxyxmodify_reverse(c)
	!elseif((.not.c%s1).and.(.not.c%s2).and.(c%s3.or.c%s4).and.c%spl==2)then
 !   call gxyxmodify_reverse(c)
	else
	return
	endif
endif

endsubroutine
  

recursive subroutine gxyxmodify_ref2(c)
implicit none
type(gridtyp),pointer::c
if(associated(c%son1))then  !!确保不会对父单元加密
    if(.not.associated(c%son3))then
    call gxyxmodify_ref2(c%son1)
    call gxyxmodify_ref2(c%son2)
	else
    call gxyxmodify_ref2(c%son1)
    call gxyxmodify_ref2(c%son2)
    call gxyxmodify_ref2(c%son3)
    call gxyxmodify_ref2(c%son4)
	endif
else
	if((c%s1.or.c%s2).and.(c%s3.or.c%s4))then
	c%typ=0
	call newson(c)
	!call nullifyCell(c%son1)
	!call nullifyCell(c%son2)
	!call nullifyCell(c%son3)
	!call nullifyCell(c%son4)
    elseif((c%s1.or.c%s2).and.(.not.c%s3).and.(.not.c%s4))then
	!c%typ=2
        c%typ=0
	call newson(c)
	!call nullifyCell(c%son1)
	!call nullifyCell(c%son2)
    elseif((.not.c%s1).and.(.not.c%s2).and.(c%s3.or.c%s4))then
	!c%typ=1
        c%typ=0
	call newson(c)
	!call nullifyCell(c%son1)
	!call nullifyCell(c%son2)
	!elseif((c%s1.or.c%s2).and.(.not.c%s3).and.(.not.c%s4).and.c%spl==1)then
 !   call gxyxmodify_reverse(c)
	!elseif((.not.c%s1).and.(.not.c%s2).and.(c%s3.or.c%s4).and.c%spl==2)then
 !   call gxyxmodify_reverse(c)
	else
	return
	endif
endif

endsubroutine

recursive subroutine gxyx_level_modify(c)
implicit none
type(gridtyp),pointer::c
if(associated(c%son1))then  !!确保不会对父单元加密
    if(.not.associated(c%son3))then
    call gxyx_level_modify(c%son1)
    call gxyx_level_modify(c%son2)
	else
    call gxyx_level_modify(c%son1)
    call gxyx_level_modify(c%son2)
    call gxyx_level_modify(c%son3)
    call gxyx_level_modify(c%son4)
	endif
else
	call gxyx_level_modify_sub(c)
endif

endsubroutine


recursive subroutine gxyx_level_modify_sub(c)
implicit none
type(gridtyp),pointer::c,c_right_neighbor,c_left_neighbor,c_up_neighbor,c_down_neighbor

    c_up_neighbor   =>up_neighbor(c)
    c_down_neighbor =>down_neighbor(c)
    c_left_neighbor =>left_neighbor(c)
    c_right_neighbor=>right_neighbor(c)
if(associated(c_up_neighbor).and.associated(c_down_neighbor)&
    &.and.associated(c_left_neighbor).and.associated(c_right_neighbor))then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!相邻网格层级相差不能大于1
if(c%lvly-c_up_neighbor%lvly>1)then
    c_up_neighbor%typ=1
    call newson(c_up_neighbor)
elseif(c_up_neighbor%lvly-c%lvly>1)then
    c%typ=1
    call newson(c)
endif
if(c%lvlx-c_up_neighbor%lvlx>1)then
    c_up_neighbor%typ=2
    call newson(c_up_neighbor)
elseif(c_up_neighbor%lvlx-c%lvlx>1)then
    c%typ=2
    call newson(c)
endif

if(c%lvly-c_down_neighbor%lvly>1)then
    c_down_neighbor%typ=1
    call newson(c_down_neighbor)
elseif(c_down_neighbor%lvly-c%lvly>1)then
    c%typ=1
    call newson(c)
endif
if(c%lvlx-c_down_neighbor%lvlx>1)then
    c_down_neighbor%typ=2
    call newson(c_down_neighbor)
elseif(c_down_neighbor%lvlx-c%lvlx>1)then
    c%typ=2
    call newson(c)
endif

if(c%lvlx-c_left_neighbor%lvlx>1)then
    c_left_neighbor%typ=2
    call newson(c_left_neighbor)
elseif(c_left_neighbor%lvlx-c%lvlx>1)then
    c%typ=2
    call newson(c)
endif
if(c%lvly-c_left_neighbor%lvly>1)then
    c_left_neighbor%typ=1
    call newson(c_left_neighbor)
elseif(c_left_neighbor%lvly-c%lvly>1)then
    c%typ=1
    call newson(c)
endif

if(c%lvlx-c_right_neighbor%lvlx>1)then
    c_right_neighbor%typ=2
    call newson(c_right_neighbor)
elseif(c_right_neighbor%lvlx-c%lvlx>1)then
    c%typ=2
    call newson(c)
endif
if(c%lvly-c_right_neighbor%lvly>1)then
    c_right_neighbor%typ=1
    call newson(c_right_neighbor)
elseif(c_right_neighbor%lvly-c%lvly>1)then
    c%typ=1
    call newson(c)
endif
endif
endsubroutine

!!!recursive subroutine gxyxmodify_reverse(c)
!!!implicit none
!!!type(gridtyp),pointer::c,cf
!!!
!!!	cf=>c%father
!!!	!!程序目前只能应对与该cf的子节点只有一层的情况
!!!	!!所以lvl设置过高会出现网格错误
!!!	 !   if(associated(cf%son2))then
!!!	 !   if(associated(cf%son2%son1))return
!!!	 !   endif
!!!	 !   if(associated(cf%son3))then
!!!	 !   if(associated(cf%son3%son1))return
!!!		!endif
!!!	 !   if(associated(cf%son4))then
!!!	 !   if(associated(cf%son4%son1))return
!!!		!endif
!!!		do while(cf%spl/=0)
!!!		cf=>cf%father
!!!		enddo
!!!	call nullifyCell(cf)
!!!	cf%spl=0
!!!	cf%typ=0
!!!	call newson(cf)
!!!	call nullifyCell(cf%son1)
!!!	call nullifyCell(cf%son2)
!!!	call nullifyCell(cf%son3)
!!!	call nullifyCell(cf%son4)
!!!	call init_s(cf) 
!!!
!!!endsubroutine



!_____________________________________________
end module
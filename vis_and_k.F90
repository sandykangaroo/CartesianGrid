!module sutherland      !用sutherland求粘性vis和热传导率k
!use inflow
!
!contains

!subroutine dynvis(T,vis)  !求粘性，输入T，输出vis粘性
!implicit none
!real(pre)::T,vis
!vis=vis0*T**1.5/(visT0+T)
!return
!end subroutine dynvis
!
!subroutine thermc(T,kk)  !求热传导率kk
!implicit none
!real(pre)::T,kk
!kk=kk0*T**1.5/(kkT0+T)
!return
!end subroutine thermc
!!k对T求导
!function dkdt(T)
!implicit none
!real(pre)::T,dkdt
!
!dkdt=kk0*(1.5*T**0.5*(kkT0+T)-T**1.5)/(kkT0+T)**2
!return
!end function dkdt
!!vis对T
!function dvisdt(T)
!implicit none
!real(pre)::T,dvisdt
!
!dvisdt=vis0*(1.5*T**0.5*(visT0+T)-T**1.5)/(visT0+T)**2
!return
!end function dvisdt
!
!end module 
    module NRRSolveMod
    use mesh
    use precisionMod
    use typeMod
    implicit none
    
    private NRRinterposeTwo
    contains
    
!_________________________________________
    subroutine NRRinterposeBlue(c)
    type(gridtyp),pointer::c
    real(pre)::dist1,dist2,d1t,d2t
    
    dist1=sqrt((c%center%x-c%NRRpointer1%center%x)**2+  &
        (c%center%y-c%NRRpointer1%center%y)**2)
    dist2=sqrt((c%center%x-c%NRRpointer2%center%x)**2+  &
        (c%center%y-c%NRRpointer2%center%y)**2)
    d1t=1/dist1; d2t=1/dist2
    c%ut=NRRinterposeTwo(c%NRRpointer1%u,c%NRRpointer2%u,d1t,d2t)
    c%vt=NRRinterposeTwo(c%NRRpointer1%v,c%NRRpointer2%v,d1t,d2t)
    c%rout=NRRinterposeTwo(c%NRRpointer1%Rou,c%NRRpointer2%Rou,d1t,d2t)
    c%Tt=NRRinterposeTwo(c%NRRpointer1%T,c%NRRpointer2%T,d1t,d2t)
    c%pt=NRRinterposeTwo(c%NRRpointer1%p,c%NRRpointer2%p,d1t,d2t)
    endsubroutine NRRinterposeBlue
!_________________________________________
    function NRRinterposeTwo(a,b,d1,d2)
    real(pre)::NRRinterposeTwo
    real(pre)::a,b,d1,d2
    NRRinterposeTwo=(a*d1+b*d2)/(d1+d2)
    endfunction NRRinterposeTwo
!_________________________________________
!_________________________________________
!_________________________________________
    endmodule NRRSolveMod
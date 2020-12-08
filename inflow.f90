!----------------------Global-------------------------------------------
    module inflowGlobal
        use precisionMod
        logical,parameter::NRR=0
        logical,parameter::MeshOnly=0
        logical,parameter::Anisotropic=0
        logical,parameter::AutoRefine=0
        logical,parameter::Limiter=1
        logical,parameter::Residual=1
        logical,parameter::Check=1
    endmodule inflowGlobal
    
module inflow
use precisionMod
!----------------------Mesh----------------------------------------------
integer,parameter::     dmax=100
integer,parameter::     lvlmax=0
integer,parameter::     lvltemp=4
integer,parameter::     solute_lvl_max=5
integer,parameter::     lvlgxyx=3
real(pre),parameter::   Rad=1.0
real(pre),parameter::   X0=50
real(pre),parameter::   Y0=50
real(pre),parameter::   h=0.1
!integer,parameter::     lvlcur=10
!----------------------AutoSave----------------------------------------------
character(20)::        FileNameStr='tmp'
!character(20)::        FileNameStr='R1NRR0L4AR4-'
logical,parameter    ::SaveFile=1   !!start auto save file.
logical,parameter    ::SaveStatue=1 !!=0,save into old file; =1,save into new file
integer,parameter    ::SaveIterate=1000
logical,parameter    ::SaveSurface=0
!----------------------Solute----------------------------------------------
integer,parameter     ::itmax=50000
integer,parameter::     AutoRefineStep=1000
integer,parameter::     AutoRefineStartStep=1000
real(pre),parameter::   Courant=1D-4  !!dt=dx/c, dx=Local mesh size.
real(pre),parameter::   u0=460.49
real(pre),parameter::   v0=0.0
real(pre),parameter::   T0=182.61
real(pre),parameter::   rou0=1.392
!real(pre),parameter::vis0=1.51e-6
!real(pre),parameter::visT0=124
!real(pre),parameter::kk0=2.5011e-3
!real(pre),parameter::kkT0=194.4
!integer::PositionNode=0,PositionElement=8485
!logical::statuePrint=0
!----------------------AutoRefine(Isotropic)----------------------------------------------
real(pre),parameter::DivParameter1=1.5
real(pre),parameter::RotParameter1=1.5
real(pre),parameter::DivParameter2=0.5
real(pre),parameter::RotParameter2=0.5
!----------------------Matrix----------------------------------------------
!integer,parameter    ::MGrids=100000     !Maximum number of grids.
!----------------------Limiter----------------------------------------------
real(pre),parameter    ::LimitUmax=   1E5
real(pre),parameter    ::LimitVmax=   1E5
real(pre),parameter    ::LimitRoumax= 1E3
real(pre),parameter    ::LimitRoumin= 1E-3
real(pre),parameter    ::LimitTmax=   1E5
real(pre),parameter    ::LimitTmin=   1
real(pre),parameter    ::LimitPmax=   5E10
real(pre),parameter    ::LimitPmin=   1
!----------------------BasicParameter----------------------------------------------
real(pre),parameter::pi=3.1415927
real(pre),parameter::gama=1.4
real(pre),parameter::R=287.0
!----------------------Others----------------------------------------------
!real(pre),parameter::tstep=4D-6
real(pre)::             tstep
real(pre)::             Rcenter(2)
integer::               CellCount   !Save the count of total cells.
integer::               OutofLimit=0  !Limiter used times.
    endmodule inflow
    
!----------------------NRR----------------------------------------------
    module inflowNRR
    use precisionMod
        logical,parameter::NRRline=1    !!Use NRRline or NRRBox
        integer,parameter::NRRlvl=4
        integer,parameter::NRRSeeds=16
        integer,parameter::NRRWidth=10  !!width count by maxlvl cell
        integer,parameter::NRRLength=12
        logical,parameter::NRRregionAutoRefine=0
        real(pre),parameter::NRRrefineTheta1=30.0   !!front degree
        real(pre),parameter::NRRrefineTheta2=110.0   !!rear degree
    endmodule inflowNRR
!----------------------Precision----------------------------------------------
    module precisionMod

    integer, parameter :: my_r4 = selected_real_kind(6 , 37)    !single kinddefs_module (IEEE 754)
    integer, parameter :: my_r8 = selected_real_kind(15, 307)   !double kinddefs_module (IEEE 754)
    integer, parameter :: pre   = my_r4
    integer, parameter :: predd = my_r8

    end module precisionMod
    
    
!----------------------TypeDefined----------------------------------------------
    module typeMod
    use precisionMod
    
    type point
        real(pre) x,y
    end type
    
    type pointdd
        real(predd) x,y
    end type
    
    type gridtyp
        integer      num,lvl     !初始单元序号，层次
        integer      sort        !网格单元类型 
        type(point)  center      !单元中心 
        type(gridtyp),pointer::father,son1,son2,son3,son4
        integer      cross,cur     !cur为1曲率加密，0则不加密
        integer      sp            !若为子单元，表明位置1,2,3,4，无子结点为0
        real(pre)::u,v,rou,T,p   !单元存储的信息，包括速度，密度，温度，压强
        !暂存的流场量,存储更新后的，再赋值给上面
        real(pre)::ut,vt,rout,Tt,pt
        integer      cnode(4)       !四个角顶点序号
        real(pre)::grax,gray,graxp,grayp,secgrax,secgray,rotx,roty,divx,divy       !速度梯度和二阶梯度
        integer::spl    !!spl=0 typ=0; spl=1 typ=1234; spl=2 typ=1324
        integer::lvlx,lvly !!gxyx的深度
        integer::nelment    !!网格单元编号
        integer      typ    !!typ=0 各向同性，typ=34 son3-son4; typ=12 son1-son2; typ=13 son1-son3; typ=24 son2-son4;
        logical::s1,s2,s3,s4
        integer::NRRregion    !!=-1:black cell. =0:blue cell. >0:red cell (NRR count).
        type(gridtyp),pointer::NRRpointer1,NRRpointer2
    end type
    endmodule typeMod

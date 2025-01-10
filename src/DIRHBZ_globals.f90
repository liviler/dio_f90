!==============================================================================!
! MODULE Globals                                                               !
!                                                                              !
! This module defines global variables that are shared by subroutines.         !                                                  !
!==============================================================================!
MODULE Globals

use Constants, only: r32,r64,i8,i16,i32,igfv,igfvbc,ngh,ngl,nghl,&
                     nt_max,ntb_max,nb_max,nz_max,nr_max,ml_max,nz_max_boson,nr_max_boson,&
                     nfgx, nkx, nhx, OUTPUT_PATH,u_config, nt3x, nh2x, NHHX, NBX,NB2X, &
                     MVX, NNNX, NHBX
implicit none

public

! --- input dio.par-------
type Input_Parameter
    character(len=10) :: force_name     ! Parameterset name of the Lagrangian
    character(len=3) :: nucleus_name    ! nucleus name
    integer(i16) :: nucleus_mass_number ! mass number
    integer(i8) :: basis_n0f            ! number of oscillator shells (fermions)
    integer(i8) :: basis_n0b            ! number of oscillator shells (bosons)
    real(r64)   :: basis_b0             ! b_0 oscillator parameter(fm**-1) of basis
    real(r64)   :: basis_q        
    real(r64)   :: basis_beta0          ! deformation parameter of basis
    real(r64)   :: woodssaxon_qs        ! WoodsSaxon qs = exp(beta2 * 3*sqrt(5/(16*pi)))
    real(r64)   :: woodssaxon_beta2     ! deformation beta2 of WoodsSaxon potential, beta2 = ln(qs)/(3*sqrt(5/(16*pi)))
    real(r64)   :: woodssaxon_beta3     ! deformation beta3 of WoodsSaxon potential
    integer(i8), dimension(2) :: pairing_ide    ! Pairing control: 1. no  2. Frozen  3. G   4. Delta
    real(r32), dimension(2)   :: pairing_dec    ! Delta,frozen gaps parameter(neutrons and protons)
    real(r32), dimension(2)   :: pairing_ga     ! Pairing-Constants(neutrons and protons)
    real(r64), dimension(2)   :: pairing_del    ! Average of single-particle pairing gaps(neutrons and protons)
    real(r64), dimension(2)   :: pairing_vpair  ! pairing strength for delta force(neutrons and protons)
    integer(i16) :: constraint_icstr ! quadratic constraint (no 0; beta2 1; b2+b3 2)
    real(r64)   :: constraint_cspr  ! spring constant
    real(r64)   :: constraint_cmax  ! cutoff for dE/db
    real(r32)   :: potential_mix    ! mixing parameter
    integer(i16) :: iteration_max    ! max number of iterations
    integer(i8) :: inin             ! Initialization of wavefunctions. 1(calc. from beginning);0(read saved pot.)
    integer(i8) :: option_iRHB  ! 0: BCS, 1: RHB
    integer(i8) :: option_iBlock ! 0: Non-blocking; 1: Block the given energy level; 2: Block according to K^\pi
    integer(i16),dimension(2) :: block_level ! block level for neutron and proton
    integer(i16),dimension(2) :: K ! block K for neutron and proton
    integer(i8),dimension(2) :: Pi ! block parity for neutron and proton
    integer(i16) :: option_blockMethod          ! 1: Block at the beginning of the iteration; 2:Block after non-blocking iteration converges; 3: Convergence blocking after convergence of non-blocking iterations
    integer(i16) :: option_Erot ! Rotation correction energy formula: 0: no; 1:Belyaev formula; 2: Nilsson formula; 3: Odd A formula
end type
type(Input_Parameter) :: input_par
!---- end input dio.par------


type Option_
    integer :: eqType ! 0: solve dirac equation, 1: solve RHB equation
    integer :: block_type !  0: Non-blocking; 1: Block the given energy level; 2: Block according to K^\pi
    integer :: block_method ! 1: blocking -> convergence;   2:convergence -> block; 3: convergence -> block -> convergence
    integer :: Erot_type  ! 0: no; 1:Belyaev formula; 2: Nilsson formula; 3: Odd A formula
end type
type(Option_) :: option


type Iteration_Parameters 
    integer(i16) :: iteration_max !maxi ! max number of iterations
    real(r32)   :: xmix  ! mixing parameter
    integer   :: ii ! iteration count
    real(r64) :: si ! maximum point of `V_out - V_in` in Broyden Mixing method
    real(r64) :: epsi = 0.00010
end type
type(Iteration_Parameters) :: iteration


type Constraint_Parameter
    integer(i32) :: length ! length of data store in betac,bet3c,clam2
    integer(i32) :: index ! index of betac,bet3c,clam2
    real(r64), dimension(:), allocatable :: betac ! Constraint beta2-values 
    real(r64), dimension(:), allocatable :: bet3c ! Constraint beta3-values
    real(r64), dimension(:), allocatable :: clam2
    ! quadratic constraint
    integer(i8) :: icstr ! quadratic constraint (no 0; beta2 1; b2+b3 2)
    real(r64) :: cspr ! spring constant
    real(r64) :: cmax ! cutoff for dE/db

    real(r64) :: c1x
    real(r64) :: c2x
    real(r64) :: c3x
    real(r64) :: c1xold
    real(r64) :: c2xold
    real(r64) :: c3xold

    real(r64) :: calq1
    real(r64) :: calq2
    real(r64) :: calq3

    real(r64), dimension(nghl,3) :: vc

end type
type(Constraint_Parameter) :: constraint


type WoodsSaxon_Parameters
    real(r64) :: qs !qs = exp(beta2 * 3*sqrt(5/(16*pi)))
    real(r64) :: beta2!betas ! deformation beta2 of WoodsSaxon potential, beta2 = ln(qs)/(3*sqrt(5/(16*pi)))
    real(r64) :: beta3!bet3s ! deformation beta3 of WoodsSaxon potential
    
    real(r64) :: v0
    real(r64) :: akv
    real(r64), dimension(2) :: r0v
    real(r64), dimension(2) :: av
    real(r64), dimension(2) :: vso
    real(r64), dimension(2) :: rso
    real(r64), dimension(2) :: aso
end type
type(WoodsSaxon_Parameters) :: woodssaxon


type Nucleus_Parameters
    character(len=2) :: name = ''!nucnam ! nucleus name
    integer(i16) :: mass_number_int!nama ! mass number
    real(r64) :: mass_number !amas ! float type of mass_number
    real(r64) :: proton_number = 0!npro ! proton number
    real(r64) :: neutron_number!nneu ! neutron number
end type
type(Nucleus_Parameters) :: nucleus_attributes


type Gauss_Integral
    ! GaussHermite
    integer :: nh ! points length
    real(r64),dimension(ngh) :: xh ! store guass points, actually is zeta
    real(r64),dimension(ngh) :: w_hermite ! store the weights of the related points
    real(r64),dimension(ngh) :: zb ! zb = xh * b_z, actually is z
    real(r64),dimension(ngh) :: wh ! wh = w_hermite * exp(xh**2)
    !GaussLaguerre
    integer :: nl
    real(r64),dimension(ngl) :: xl ! actually is eta
    real(r64),dimension(ngl) :: w_laguerre
    real(r64),dimension(ngl) :: rb ! rb = sqrt(xl) * b_\perp, actually is r_perp
    real(r64),dimension(ngl) :: wl ! wl = w_laguerre * exp(xl)

    real(r64),dimension(ngh,ngl) :: wdcor ! b_z (b_\perp)^2 \pi*wh*wl
end type
type(Gauss_Integral) :: gauss

type Cylindrical_Harmonic_Oscillator_Basis
    ! parameter
    integer(i8) :: n0f,n0b  !n0f n0b  ! number of oscillator shells
    real(r64)   :: b0       ! b_0 oscillator parameter(fm**-1) of basis
    real(r64)   :: beta0    !deformation parameter of basis
    real(r64)   :: q        !momentum?
    real(r64)   :: hb0      ! hbc/(2* amu)
    real(r64)   :: hom      ! hbar*omega_0
    real(r64)   :: bp       ! b_p/b_0
    real(r64)   :: bz       ! b_z/b_0
    ! -----------------------------------------------------------
    ! quantum numbers
    integer,dimension(nt_max) :: nz,nr,ml,ms
    character(len=8),dimension(nt_max) :: tb
    integer :: nzm,nrm,mlm ! max of nz,nr,ml
    integer,dimension(nb_max,2) :: ia  ! ia(b,1):begin of the large components of block b is ia(b,1)+1 ;ia(b,2):begin of the small components of block b is ia(b,2)+1
    integer,dimension(nb_max,2) :: id  ! id(b,1):dimension large components of block b;id(b,2):dimension small components of block b
    integer,dimension(nb_max) :: iag0  ! begin of the spurious components of block b is iag0(b)+1 
    integer,dimension(nb_max) :: idg0  ! number of spurious components of block b
    integer,dimension(nb_max) :: ikb   ! K-quantum number of each block (K+1/2) 
    integer,dimension(nb_max) :: ipb   ! parity of each block ( 1 for +, 2 for -)
    integer,dimension(nb_max) :: mb    ! 
    character(len=25),dimension(nb_max) :: txb
    integer :: nb ! number of K-parity-blocks 
    integer :: nt ! number of total levels

    integer,dimension(ntb_max) :: nzb,nrb
    integer :: nzbm,nrbm ! max of nzb,nrb
    integer :: ntb !no ! number of total levels for bosons
    ! -----------------------------------------------------------
    ! Basis in gauss points
    real(r64),dimension(:,:),allocatable :: qh  ! z-function at meshpoint xh(ih) with sqrt(wh)
    real(r64),dimension(:,:),allocatable :: qh1 ! the first-order derivative of z-funtion with sqrt(wh)
    real(r64),dimension(:,:,:),allocatable :: ql ! r-function at meshpoint xl(il) with sqrt(wl/2) 
    real(r64),dimension(:,:,:),allocatable :: ql1 ! the first-order derivative of r-funtion with r*sqrt(wl/2);NOTE: r d/dr psi_(nr,L)(rr) = QL1(nr,L,i) * sqrt( 2 / WL(i) )
    real(r64),dimension(:,:),allocatable :: qhb ! the z-function for the expansion of the mesonfields
    real(r64),dimension(:,:),allocatable :: qlb  ! the r-function for the expansion of the mesonfields
    !------------------------------------------------------------
    ! single particle states expansion coefficients
    ! real(r64),dimension(nhx,nkx,2) :: fg ! see dirac%fg ! energy state, fg(coefficients,the nth eigenstate,proton or neutron)
end type
type Spherical_Harmonic_Oscillator_Basis
    ! parameter
    integer(i8) :: n0f
    ! -----------------------------------------------------------
    ! quantum numbers
    integer,dimension(nt3x,4) :: nljm ! nr, l, j+1/2, m_j+1/2 
    integer,dimension(nt3x) :: nps    ! large( 1 ) or samll( 2 ) components
    integer,dimension(nb_max,2) :: iasp ! iasp(b,1): begin position of large components; iasp(b,2): begin position of small components
    integer,dimension(nb_max,2) :: idsp ! idsp(b,1): dimension of large components; idsp(b,2): dimension of small components
    character(len=14),dimension(nt3x) :: tts
    !-------------------------------------------------------------
    ! single particle states expansion coefficients
    real(r64), dimension(nh2x,nkx,2) :: fg
end type
type Basis
    type(Cylindrical_Harmonic_Oscillator_Basis) :: HO_cyl
    type(Spherical_Harmonic_Oscillator_Basis)   :: HO_sph
end type
type(Basis) BS


type Fields_
    integer(i8) :: inin  ! Initialization of fields. 1(calc. from beginning);0(read saved pot.)
    real(r64),dimension(ngh,ngl,2) :: vps !V+S
    real(r64),dimension(ngh,ngl,2) :: vms !V-S
    real(r64),dimension(ngh,ngl,2) :: vpstot
    real(r64),dimension(ngh,ngl,2) :: vmstot
 
    real(r64),dimension(ngh,ngl) :: coulomb
    real(r64),dimension(ngh,ngl) :: sigma
    real(r64),dimension(ngh,ngl) :: omega
    real(r64),dimension(ngh,ngl) :: delta
    real(r64),dimension(ngh,ngl) :: rho
end type
type(Fields_) :: fields

! type broyden_
!    real(r64),dimension(nn,mm) :: df ! \Delta F, but df(i,inex) sotre latest F 
!    real(r64),dimension(nn,mm) :: dv ! \Delta V, but dv(i,inex) sotre latest V_{in}
!    real(r64),dimension(nn,mm) :: bw0 ! \omega_0
!    real(r64),dimension(mm,mm) :: bbeta ! \beta_{kn}
!    real(r64),dimension(mm)    :: bwork ! c_k^m
!    real(r64),dimension(nn)    :: curv ! alpha * F^m
!    integer  ,dimension(mm)    :: ibwork
! end type
! type(broyden_) :: broyden


type Matrix_
    real(r64),dimension(nfgx,nb_max) :: sp
    real(r64),dimension(ntb_max*ntb_max,4) :: meson_propagators 
end type
type(Matrix_) matrices


type dirac_
    real(r64),dimension(nkx,2) :: ee ! energies
    real(r64),dimension(nhx,nkx,2) :: fg ! energy state, fg(coefficients,the nth eigenstate,proton or neutron)
    integer, dimension(2) :: nk ! dimension of ee, the number of energies
    integer, dimension(nb_max,2) :: ka ! ka(k,itx): begining of energies of k block in ee
    integer, dimension(nb_max,2) :: kd ! kd(k,itx): dimension of energies of k block in ee
    real(r64),dimension(nkx,2) :: nblo

    real(r64),dimension(NHHX,NB2X) :: hh ! store Hamiltonian matrix element for RHB equation
end type
type(dirac_) dirac


type Pairing_Parameters
    integer(i8), dimension(2) :: ide
    real(r32), dimension(2) :: dec !frozen gaps parameter(neutrons and protons)
    real(r32), dimension(2) :: ga  ! Pairing-Constants

    real(r64), dimension(2) :: del !average of single-particle pairing gaps : \Delta^{uv}\equiv \frac{\sum_{k} f_{k} u_kv_{k} \Delta_{k}}{\sum_{k} f_{k} u_kv_{k}}
    real(r64), dimension(2) :: gg ! gg = ga/A
    
    real(r64) :: pwi
    real(r64), dimension(2) :: vpair !pairing strength for delta force
    real(r64), dimension(nghl,2) :: delq ! pairing potential: \Delta_{tau}(\boldsymbol(r))
    
    real(r64), dimension(nkx,2) :: de  ! pairing gaps: \Delta_{k}
    real(r64), dimension(nkx,2) :: vv  ! 2*v_k^2 , v_k^2 is the occupied probability of a sigle particle state and 0< vv < 2
    real(r64), dimension(nkx,2) :: skk ! cutoff weight: f_k


    real(r64),dimension(2) :: ala ! fermi energy
    real(r64),dimension(2) :: ecut ! cut-off energy in pairing window
    real(r64),dimension(2) :: disp ! fluctuations of particle number: <(\Delta N) ^2> =  <N^2>-<N>^2 = 4*sum_(k>0) u^2_k * v^2_k
    ! real(r64),dimension(2) :: epair! pairing energy: E_{pair} = \sum_{k>0} f_{k} u_{k} v_{k} \Delta_{k}
    real(r64),dimension(2) :: spk ! trace of kappa: \sum_{k} f_k u_k v_k
    real(r64),dimension(2) :: dev2 ! average of single-particle pairing gaps : <\Delta>\equiv \frac{\sum_{k} f_{k} v_{k}^{2} \Delta_{k}}{\sum_{k} f_{k} v_{k}^{2}}

    ! Odd neutron/proton nucleus
    integer,dimension(2) :: block_level ! block level for BCS w.f.
    integer,dimension(2) :: block_K ! block K ! 1: 1/2;  2: 3/2; ......
    integer,dimension(2) :: block_Pi ! block parity, 1: + ; -1: -
    integer,dimension(50) :: ibk ! block level of K-Block
    logical :: allow_block ! control block action in BCS loop

end type
type(Pairing_Parameters) :: pairing


type densities_
    real(r64), dimension(nghl,2) :: rs ! scalar density; rs(i,1): \rho_s of neurons; rs(i,2): \rho_s of proton
    real(r64), dimension(nghl,2) :: rv ! vector density; rv(i,1): \rho_v of neurons; rv(i,2): \rho_v of proton
    real(r64), dimension(nghl,2) :: drs
    real(r64), dimension(nghl,2) :: drv
    real(r64), dimension(nghl,4) :: ro ! ro(i,1)=\rho_s; ro(i,2)=\rho_v ; ; ro(i,4)=\rho_3
    real(r64), dimension(nghl,4) :: dro
    real(r64), dimension(nghl,2) :: rkapp ! bcs \kappa
    real(r64), dimension(nghl) :: drvp
    
    real(r64),dimension(NHHX,NB2X) :: rosh ! rho_{n,n'}
    
end type
type(densities_) densities


type expectation_
    real(r64) :: rc ! charge radius
    real(r64) :: rms ! root mean square
    

    real(r64) :: beta2 
    real(r64) :: beta3
    
    real(r64) :: betg
    real(r64) :: beto

    real(r64) :: dd0 ! dipole moment: D0=e*(N/A*zp-Z/A*zn)
    real(r64) :: ddz 

    real(r64) :: qq2p
    real(r64) :: qq3p

    real(r64) :: etot ! total energy E
    real(r64) :: ecm
    real(r64) :: ea  ! E/A

    real(r64) :: Erot !rotational energy

end type
type(expectation_) expectations


type Output_FileName
    ! config file
    character(len=40) :: config = OUTPUT_PATH//'dio_config.out'
    integer :: u_config   = u_config
    ! standard output
    character(len=40) :: outputf
    integer :: u_outputf  = u_config + 1
    ! 
    character(len=40) :: outputw
    integer :: u_outputw  = u_config + 2
    !
    character(len=40) :: outdel
    integer :: u_outdel   = u_config + 3
    !
    character(len=40) :: outputwf
    integer :: u_outputwf = u_config + 4
    ! density
    character(len=40) :: outputd
    integer :: u_outputd = u_config + 5
    ! rotational correction energy
    character(len=40) :: rotationalE = OUTPUT_PATH//'E_rot.out'
    integer :: u_rotationalE = u_config + 11
    !
    character(len=40) :: outExpectation = OUTPUT_PATH//'Expectation.out'
    integer :: u_outExpectation = u_config + 12

end type
type(Output_FileName) :: outputfile

! -----define force parameters
type mass_parameters
    real(r64) :: amu
    real(r64) :: amsig
    real(r64) :: amome
    real(r64) :: amdel
    real(r64) :: amrho
    real(r64) :: ampi
end type
type couplg_parameters
    ! quadratic terms
    real(r64) :: ggsig
    real(r64) :: ggome
    real(r64) :: ggdel
    real(r64) :: ggrho
    integer(i8), dimension(4) :: lmes
end type
type coupld_parameters
    ! derivative terms
    real(r64) :: ddsig
    real(r64) :: ddome
    real(r64) :: dddel
    real(r64) :: ddrho
end type
type coupnl_parameters
    ! non-linear terms
    real(r64) :: ggbet
    real(r64) :: gggams
    real(r64) :: gggamv
end type
type couplm_parameters
    real(r64) :: gsig
    real(r64) :: gome
    real(r64) :: gdel
    real(r64) :: grho
    real(r64) :: gpi
end type
type nonlin_parameters
    real(r64) :: g2
    real(r64) :: g3
    real(r64) :: c3
end type
type dtypel_parameters
    real(r64) :: a_s
    real(r64) :: b_s
    real(r64) :: c_s
    real(r64) :: d_s
    real(r64) :: a_v
    real(r64) :: b_v
    real(r64) :: c_v
    real(r64) :: d_v
    real(r64) :: a_tv
    real(r64) :: b_tv
    real(r64) :: c_tv
    real(r64) :: d_tv
    real(r64) :: dsat
end type
type dpolyn_parameters
    real(r64), dimension(4) :: alf
    real(r64), dimension(4) :: bet
    real(r64), dimension(4) :: gam
end type
type option_parameters
    integer(i8) :: inl !non-linear meson-coupling
    integer(i8) :: ipc !density-dependent meson-coupling
    integer(i8) :: idd !point-coupling
end type

type Force_Parameter_Sets
    character(len=10) :: parname ! Parameterset name of the Lagrangian
    type(mass_parameters) :: masses
    type(couplg_parameters) :: couplg
    type(coupld_parameters) :: coupld
    type(coupnl_parameters) :: coupnl
    type(couplm_parameters) :: couplm
    type(nonlin_parameters) :: nonlin
    type(dtypel_parameters) :: dtypel
    type(dpolyn_parameters) :: dpolyn
    type(option_parameters) :: option
    real(r64) :: rosat
    real(r64), dimension(nghl,4,2) :: ff ! ff(i,m,1), ff(i,m,2) are the density-dependent coupling constant and its derivation respectively for meson m. 
end type
type(Force_Parameter_Sets) :: force
! -----end define force parameters


type math_gfv
    integer,dimension(0:igfv) :: iv
    real(r64),dimension(0:igfv) :: sq
    real(r64),dimension(0:igfv) :: sqi
    real(r64),dimension(0:igfv) :: sqh
    real(r64),dimension(0:igfv) :: shi
    real(r64),dimension(0:igfvbc,0:igfvbc) :: ibc
    real(r64),dimension(0:igfv) :: fak
    real(r64),dimension(0:igfv) :: fad
    real(r64),dimension(0:igfv) :: fi
    real(r64),dimension(0:igfv) :: fdi
    real(r64),dimension(0:igfv) :: wf
    real(r64),dimension(0:igfv) :: wfi
    real(r64),dimension(0:igfv) :: wfd
    real(r64),dimension(0:igfv) :: gm2
    real(r64),dimension(0:igfv) :: gmi
    real(r64),dimension(0:igfv) :: wg
    real(r64),dimension(0:igfv) :: wgi
end type
type(math_gfv) :: gfv

type RHB_Pairing_
    real(r64) :: ga ! a^2  
    real(r64), dimension(2) :: gl !G
    real(r64), dimension(NHHX,NB2X) :: delta
    real(r64), dimension(MVX,NNNX)  :: wnn
    real(r64), dimension(MVX,2) :: kappa !  pairing tensor
end type
type(RHB_Pairing_) :: RHB_pairing

type RHB_
    integer, dimension(NBX,4) :: ka ! ka(k,itx): begining of states of k block; itx=1:(neutron, positive state), itx=3:(neutron, negative), itx=2:(proton,positive), itx=4:(proton,negative) 
    integer, dimension(NBX,4) :: kd ! kd(k,itx): dimension of energies of k block; itx=1:(neutron, positive state), itx=3:(neutron, negative), itx=2:(proton,positive), itx=4:(proton,negative) 
    real(r64), dimension(nkx,4) :: equ ! energies
    real(r64), dimension(NHBX,nkx,4) :: fguv ! U,V
    ! real(r64),dimension(2) :: ala ! fermi energy
end type
type(RHB_) :: RHB

type OddA_
    integer :: qusiparticle_state_neutron ! create one quasiparticle state of neutron as odd-mass nuclei
end type
type(OddA_) OddA
END MODULE Globals
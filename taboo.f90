!
!
!   +------------------------------------------------------+
!   |                                                      |
!   |                                                      |
!   |                +---++---++---++------+               |
!   |                | T || A || B ||  OO  |               |
!   |                +---++---++---++------+               |
!   |                                                      |
!   |        ( a posT glAcial reBOund calculatOr )         |
!   |                                                      |
!   |                by Giorgio Spada et al. 	           |
!   |                                                      |
!   |                                                      |
!   |              Release 1.1 -  February 2017            |
!   |                                                      |
!   |                                                      |
!   |         http://github.com/danielemelini/TABOO        |
!   |                                                      |
!   |                                                      |
!   +------------------------------------------------------+
!                         
!
! v1.1 - DM, February 6, 2017
!        Ported to GNU gfortran and added two NV=3 models
!      - DM, February 7, 2017
!        Now TABOO computes the 'true' Earth and density for
!        M3-L70-V01. Adjusted the precision of PI values and
!        of the gravitational constant. Adjusted ice density to
!        931 kg/m**3 [r11]
!      - DM, June 16, 2020
!        Fixed some issues with External_Model
!      - DM, May 10, 2021
!        Added the VM5a model (NV=5, CODE=0) and the
!        VM7 model (NV=8, CODE=0)
!      - DM, Nov 29, 2022
!        Implementation of the degree 1 loading LN
!
 module STRATA 
! defines the kinds an PI 
 implicit NONE 
 save
 integer, parameter  :: i4b   = selected_int_kind(9)
 integer, parameter  :: sp    = kind(1.0)
 integer, parameter  :: dp    = kind(1.0d0)
 integer, parameter  :: qp    = kind(1.0q0)
 real(qp), parameter :: pi    = 3.141592653589793238462643383279502884197_qp
 end module STRATA 
!
!
!
!
!
module trigfcn
! defines trigonometric functions with input in degrees
  implicit none
!
  private
!
  integer, parameter :: sp = kind(1.e0), dp = kind(1.d0), qp = kind(1.q0)
  real(sp), parameter :: d2r_sp = 4 * atan(1.0_sp) / 180.0_sp
  real(dp), parameter :: d2r_dp = 4 * atan(1.0_dp) / 180.0_dp
  real(qp), parameter :: d2r_qp = 4 * atan(1.0_qp) / 180.0_qp
!
  public sind, cosd
!
  interface sind
    module procedure ssind ! Single precision 
    module procedure dsind ! Double precision
    module procedure qsind ! Quad precision
  end interface
!
  interface cosd
    module procedure scosd ! Single precision 
    module procedure dcosd ! Double precision
    module procedure qcosd ! Quad precision
  end interface
!  
  contains
!
    function ssind(x)
      implicit none
      real(sp) :: ssind
      real(sp) :: x
      ssind = sin(d2r_sp * x)
    end function ssind
!
    function dsind(x)
      implicit none
      real(dp) :: dsind
      real(dp) :: x
      dsind = sin(d2r_dp * x)
    end function dsind
!
    function qsind(x)
      implicit none
      real(qp) :: qsind
      real(qp) :: x
      qsind = sin(d2r_qp * x)
    end function qsind
!
    function scosd(x)
      implicit none
      real(sp) :: scosd
      real(sp) :: x
      scosd = cos(d2r_sp * x)
    end function scosd
!
    function dcosd(x)
      implicit none
      real(dp) :: dcosd
      real(dp) :: x
      dcosd = cos(d2r_dp * x)
    end function dcosd
!
    function qcosd(x)
      implicit none
      real(qp) :: qcosd
      real(qp) :: x
      qcosd = cos(d2r_qp * x)
    end function qcosd
!
end module trigfcn
!
!
!
!
!
module COMMON
!-----------------------------------------------------------
! Declares variables common to all of the three TABOO Tasks 
!-----------------------------------------------------------
use STRATA 
implicit none 
save 
!
 real(dp), parameter :: raggio    = 6.371D6    ! Earth radius (m)
 real(dp), parameter :: ref_emass = 5.97D24    ! Reference Earth mass (kg)
 real(dp), parameter :: ref_rhoea = 5.51157D3  ! Reference average density (kg/m3)
 real(dp), parameter :: ext_rhoea = 5.51467D3  ! Average density for the external model (kg/m3)
 real(dp), parameter :: rhoice    = 931.D0     ! Ice density (kg/m3)
!
 real(dp) :: emass                             ! Earth mass (kg)
 real(dp) :: rhoea                             ! Earth density (kg/m3)
!
integer(i4b), parameter :: nv_max     =  24  ! maximum allowed number of v.e. layers 
integer(i4b), parameter :: nroots_max =  96  ! maximum allowed number of modes  
integer(i4b), parameter :: llmin      =   0  ! minimum allowed degree 
integer(i4b), parameter :: llmax      = 128  ! maximum allowed degree 
!
INTEGER(i4b) :: IV          ! Verbose (1) or Silent (0) mode
!
integer(i4b) :: lmin, lmax  ! Current values of lmin e lmax, taken from TASK_*.DAT      
!
integer(i4b) :: NV          ! NV is the number of viscoelastic layers  
integer(i4b) :: NROOTS      ! NROOTS is the number of viscoelastic modes 
!
real (qp) :: r  (0:nv_max+1)   ! Radii 
real (qp) :: rho(0:nv_max+1)   ! Density 
real (qp) :: rmu(0:nv_max+1)   ! Shear Moduli 
real (qp) :: vis(0:nv_max+1)   ! Viscosity 
real (qp) :: pon(0:nv_max+1)   ! Gravity  
!
real (qp) :: s (llmin:llmax,0:nroots_max)      ! roots of the secular equation, in kyr**(-1)  
!
integer (sp) :: vec (llmin:llmax,nroots_max)   ! Marker for fake modes 
!
real (qp) :: h_e(llmin:llmax), l_e(llmin:llmax), k_e(llmin:llmax) ! Elastic modes
!
real (qp) :: h_f(llmin:llmax)  ! Fluid modes 
real (qp) :: l_f(llmin:llmax)  !   "     "
real (qp) :: k_f(llmin:llmax)  !   "     "
!
real (qp) :: h_v (llmin:llmax,1:nroots_max) ! Viscoelastic residues h
real (qp) :: l_v (llmin:llmax,1:nroots_max) !       "          "    l 
real (qp) :: k_v (llmin:llmax,1:nroots_max) !       "          "    k 
!
real (dp) :: SIGMA (llmin:llmax)    ! Harmonic coeff. for an axis-symmetric load
!
! _Sine and cosine coeff. for a rectangular load (real harmonics)_  
real (dp) :: FF(llmin:llmax,llmin:llmax), GG(llmin:llmax,llmin:llmax) 
!                             
! _Sine and Cosine Ocean function coefficients (real harmonics)_ 
real (dp) :: OC(llmin:llmax,llmin:llmax), OS(llmin:llmax,llmin:llmax) 
!     
! _Arrays containing the convolution between the time history & the Green functions_                      
real (sp) ::  conv_h(llmin:llmax),     conv_l(llmin:llmax),     conv_k(llmin:llmax)
real (sp) :: dconv_h(llmin:llmax),    dconv_l(llmin:llmax),    dconv_k(llmin:llmax)
!
REAL (sp) :: r_h (llmin:llmax,1:nroots_max)   ! Normalized residue for h
REAL (sp) :: r_l (llmin:llmax,1:nroots_max)   ! Normalized residue for l
REAL (sp) :: r_k (llmin:llmax,1:nroots_max)   ! Normalized residue for k
REAL (sp) :: tek (llmin:llmax,1:nroots_max)   ! Relaxation time in k-yrs
!
! _Statements for time-history l_h=0 (loading)_
   real(sp) :: tsl
!
! _Statements for time-history l_h=1 (unloading)_
   real(sp) :: tsu
!
! _Statements for time-history l_h=2 (loading and unloading)_
   real(sp) :: tau_lu
!  real(sp) :: tsu
!
! _Statements for time-history l_h=3 (simple deglaciation)_
   real(sp) :: tau_sd
   real(sp) :: tsbd
!
! _Statements for time-history l_h=4 (saw tooth)_
   integer (i4b) :: nr
   real(sp) :: tauc, dinc
   real(sp) :: tsbld
!
! _Statements for time-history l_h=5 (periodic)_
   real(sp) :: period
   REAL(sp) :: tsz
!
! _Statements for time-history l_h=6 (piecewise linear)_
  integer (i4b) , parameter :: nss_max = 100
  integer (i4b) :: nss
  REAL(sp) :: alf (0:nss_max), bet(0:nss_max)
  real(sp) :: tnn                             ! durata
  real(sp) :: a_hh    (0:nss_max)             ! altezze ai tempi t(k)
  real(sp) :: t_hh    (0:nss_max)             ! tempi t(k)
  REAL(sp) :: ri      (1:nss_max)
  REAL(sp) :: tsbplp
!
! _Statements for time-history l_h=7 (piecewise constant)_
  integer (i4b), parameter :: ns_max = 100
  real(sp), parameter :: dilta_min = 0.001  ! 1 yr is the *minimun dilta* allowed
  REAL(sp), parameter :: dilta_max = 1000.  ! 1000 kyr is the *maximum* allowed
  integer (i4b) :: ns                                         
  real(sp) :: tn, dilta                     ! durata e incremento
  real(sp) :: a_h(0:ns_max+1)               ! h del carico glaciale nel tempo
  REAL(sp) :: tsbpcp
!
! _Statements for time-history l_h=8 (piecewise constant with loading phase)_
  REAL(sp) :: tauf                         ! lenght of the finite loading phase
  real(sp), parameter :: tauf_min = 0.001  ! 1 yr is the *minimun lenght* allowed
  REAL(sp), parameter :: tauf_max = 1000.  ! 1000 kyr is the *maximum* allowed
! REAL(sp) :: tsbpcp
! real(sp), parameter :: dilta_min = 0.001  ! 1 yr is the *minimun dilta* allowed
! REAL(sp), parameter :: dilta_max = 1000.  ! 1000 kyr is the *maximum* allowed
! integer (i4b) :: ns
! real(sp) :: tn, dilta                     ! durata e incremento
! real(sp) :: a_h(0:ns_max+1)               ! h del carico glaciale nel tempo
! REAL(sp) :: tsbpcp
!
!
!_Switches for the definition of the variables to be computed_
  Integer (i4b) :: IR, dIR, it, dit, IG, dIG
!
!_Switch for the computation of the 'Pressure' (kg/m**2)_ 
  Integer (i4b) :: IPRESS 
!
!_Switch for the computation of Elastic Deformations_ 
  Integer (i4b) :: Only_Elastic  
!
end module COMMON 
!
!
!
!
module COMMON_FOR_SPECTRUM
!
! -------------------------------------------------------------------------+
! Variables in use by Sbr. 'spectrum' and in the routines depending on it  | 
! -------------------------------------------------------------------------+
!
use COMMON
use STRATA 
implicit none 
!
real (qp) :: pivo                 ! Largest coeff. of the secular polynomial 
real (qp) :: rubens               ! For the elastic part of the solution 
!
real (qp) :: ggg, xmass       ! gravity constant and Earth mass 
!
real (qp) :: aco (0:nv_max+1) ! <<A>> coefficients (bottom to top)  
real (qp) :: g (0:nv_max+1)   ! gravity at the interfaces (bottom to top) 
!
real (qp) :: a (6, 6), b (6, 6), c (6, 6), d (6, 6)  ! 6*6 propagators 
real (qp) :: ac(6, 6), ad(6, 6), bc(6, 6), bd(6, 6)  ! Outputs of MATPROD 
real (qp) :: k0 (6, 6), k1 (6, 6), k2 (6, 6)         ! 6*6 propagators 
real (qp) :: matela (6, 6)                           ! Elastic product 
real (qp) :: sinist_2 (3, 6), sinist_1 (3, 6)        ! Left products 
!
real (qp) :: rr (0:2 * nv_max, 3, 3)                     ! Matrix to be inverted 
real (qp) :: co (0:3 * 2 * nv_max), aa(nroots_max + 1)   ! Polynomial coefficints 
!
real (qp) :: rad (llmin:llmax,nroots_max) ! Roots in years 
!
real (qp) ::      cc (0:2, 1:nv_max, 6, 6) 	   ! Propagator 
real (qp) :: coefmat (0:2 * nv_max, 6, 6)          ! Propagator 
!
real (qp) :: ctmp (0:3 * 2 * nv_max - 1)         ! Derivative of the secular poly. 
!
real (qp) :: rrrr (0:2 * nv_max, 3, 6), qqqq (0:2 * nv_max, 3, 6)  ! Left products 
!
real (qp) :: cmb(6,3), bcs (3)    ! CMB and surface Boundary conditions 
!
real (qp) :: qq (0:2 * nv_max, 3, 3)  ! Q matrix 
!
real (qp) :: r_r (3, 3, 0:nroots_max), q_q (3, 3, 0:nroots_max)  ! Ausilium di R e Q  
real (qp) :: qr (3, 3, 0:nroots_max)                           ! Product Q*R 
real (qp) :: raggiu (3, 3, 0:nroots_max)                     ! Adjoint 
real (qp) :: derpo (1:nroots_max)                          ! Derivative of the secular poly. 
!
real (qp) :: rt1 (1:nroots_max), rt2 (1:nroots_max)  ! Real and Imag. parts of the roots, in kyears  
!
real (qp) ::  x_el (3), xx (3), xr (3, nroots_max)   ! Solution in vector form 
!
! _Projectors P1 and P2_ 
integer (i4b) :: io 
real (qp) :: pj1 (3, 6) = reshape ( (/ 1._qp, (0._qp,io = 2, 4), 1._qp,   &
         (0._qp, io = 6, 14), 1._qp, (0._qp, io = 16, 18) /), (/ 3, 6 /))
real (qp) :: pj2 (3, 6) = reshape ( (/ (0._qp,io = 1, 6), 1._qp,          & 
  (0._qp, io = 8, 10), 1._qp, (0._qp, io = 12, 17), 1._qp /), (/ 3, 6 /))
!
real (qp), external :: ded                           ! Useful to compute the determinant
!
!
END MODULE COMMON_FOR_SPECTRUM 
!
!
!
!
!
!
!
!+++++++++++++++++++++
      module for_task2
!+++++++++++++++++++++
!
! ----------------------------------------------------------------+
! Variables in use by Sbr. TASK_2 and in many of its dependencies |
! ----------------------------------------------------------------+
!
 USE COMMON
 implicit none
 save
!
 Real(qp) :: LT      ! Lithospheric Thickness
!
 REAL(dP) :: VECT(4), D_VECT(4)     ! Displacements & derivatives 
 REAL(dP) :: CORR(4), D_CORR(4)     ! Their oceanic corrections   
!
 REAL(dP) :: INER(3,3), D_INER(3,3) ! Inertia change & time-derivative 
 REAL(dP) :: W(2), D_W(2)           ! Cosine & sine Stokes coefficients 
!
 REAL(dP) :: JUNK                ! A junk
 Real(dp) :: HEIGHT              ! Maximum height of the load
 Real(dp) :: MASS		 ! Maximum mass of the load
 REAL(dP) :: PRESS               ! Maximum pressure at the pole for harmonic loads (kg/m**2)
 Real(dp) :: AMPLITUDE           ! Angular semi-amplitude of the load
 Real(dp) :: LONG_CEN, TETA_CEN  ! Long. and Colat. of the load center (centroid)
 REAL(dP) :: DL, DT              ! Delta long. and Delta Col. for RECT_LOAD0
 REAL(dP) :: CEN(2)              ! Center (or centroid) of the load (lon. & colat.)  
 Real(dp) :: OBS(2)              ! Longitude and colatitude of the observer
!
 REAL(sp) :: TIME             ! _Local_ time for the analysis
 REAL(SP) :: tmin, tmax       ! Time bounds
 REAL(SP) :: deltat           ! Time increment
!
 REAL(DP) :: given_time       ! A given time
 REAL(SP) :: load_history     ! Function which computes the TH
    EXTERNAL load_history
 REAL(SP) :: d_load_history   ! Function for the time-derivatives
    EXTERNAL d_load_history
 REAL(DP) :: FACT1, FACT2     ! Constant factors
 REAL(DP) :: MICE             ! Mass of the load as a function of time 
!
!   _These are useful to compute the elapsed time_
 Real(sp) :: txt(2), ela1, ela2   ! It was: Real(sp) :: etime, txt(2), ela1, ela2; EXTERNAL ETIME 
!
!        _Coordinates of the Observers_
    integer (i4b), parameter :: nobs_max=16200 ! Maximum accepted observers number
    integer (i4b) :: nobs                      ! Number of observers, as from input
    real(dp) :: long_obs(nobs_max), cola_obs(nobs_max)  ! Long. & Colat. of the obserbers 
    real(dp) :: L1, L2, C1, C2, EPL, EPC, CC, LL  ! Map bounds & increments in long. & colat.  
!
  integer (i4b), parameter :: max_num_righe=20001  ! Max. number of rows in task_2.dat
  INTEGER (i4b) :: ideja    ! checks for errors in the model KWs
  INTEGER (i4b) :: inone    ! checks for errors in the model KWs
  INTEGER (i4b) :: iarmo    ! detects the HArmonic_Degrees KW
  INTEGER (i4b) :: iload    ! detects the KW Load_Geometry
  INTEGER (I4b) :: ihist    ! detects the KW Load_History  
  INTEGER (i4b) :: n_load   ! counts the number of Load_Geometry KWs
  INTEGER (i4b) :: n_hist   ! counts the number of Load_History KWs
  integer (i4b) :: n_inte   ! counts the number of Local_Study KWs found 
  integer (i4b) :: n_este   ! counts the number of Global_Study KWs found  
  integer (i4b) :: NPUL, NPUC        ! points in long. & colat. for the maps    
  integer (i4b) :: CL                ! Load type 
  integer (i4b) :: LH    	     ! Time-History type 
  integer (i4b) :: type_inten  	     ! Type of local study 
  integer (i4b) :: type_esten        ! Type of global study 
  integer (i4b) :: IOC               ! Real Ocean Switch 
  integer (i4b) :: ISUB              ! Switch for AXIS_DISP0 & ESSA_IEZI0  
  integer (i4b) :: ITASK, JJ, L, M, I, J, K, IJUNK, JP   ! Various control variables 
  integer (i4b) :: LDE, MDE, LDE1, LDE2    ! Degrees & orders 
  integer (i4b) :: INDX; EXTERNAL indx     ! INDX = l*(l+1)/2+m+1
  Integer(I4B) :: CODE                 ! Codes to identify the model, given NV 
  Integer(I4B) :: ILM                  ! (0/1) An option for NV=7 and 9. 
  Integer(I4B) :: Inorm                ! (0/1) Switch for the Stokes coeff. normalization. 
!
  INTEGER(i4b) :: i_loading            ! loading/tidal switch
  integer(i4b) :: drop_modes           ! "drop modes" switch
!                                        
  CHARACTER*30 FILE_SPARSI     ! A filename 
  CHARACTER*30 KEYWORD          ! A Task#2 Keyword
  CHARACTER*20 date, timc         ! date and time
!
!+++++++++++++++++++++
  end module for_task2
!+++++++++++++++++++++
!
!



!
!
!
!+++++++++++++++
module FOR_TASK3  
!+++++++++++++++
!
! -----------------------------------------------------------------------------+
! Variables in use by Sbr. TASK_3 and in many of the routines depending on it  | 
! -----------------------------------------------------------------------------+
!
 USE COMMON
 implicit none
 save
!
 Integer(i4b)  Numero_di_elementi     ! Total number of elements of the aggregate                                                    
 Integer(i4b), parameter :: mel = 809 ! Maximum number of ice elements allowed
!
 REAL(DP):: TI_HI    ! Value of the TH at the current time 
 REAL(DP):: PRES     ! Pressure (kg/m**2) at point due to a load aggregate
 REAL(DP):: PRES_EL  ! Pressure (kg/m**2) at point due to a single element 
!
 Real(qp)  :: LT      ! Lithospheric Thickness for the call to sbr. SPEC
 Real(DP)  :: MICE    ! Mass of the Aggregate at a given time (kg)
 Real(DP)  :: Area_Oceans      ! Area of the Oceans 
!
 Real(sp)  :: Long_c (1:mel) ! Longitude of the center (centroid) of ice elements (deg)
 Real(sp)  :: Cola_c (1:mel) ! Colatitude   "        "      "      "             ( " )
 Real(sp)  :: Amplit (1:mel) ! Angular amplitude   (deg)
 Real(dp)  :: Grande (1:mel) ! MAXIMUM Mass in time of an ice element (kg)
 Real(sp)  :: Altezz (1:mel) ! MAXIMUM Thickness in time of an ice element (m)
 Real(sp)  :: D_long (1:mel) ! Amplitude in Long.
 Real(sp)  :: D_cola (1:mel) ! Amplitude in Colat.
!
! *** Block of statements useful for the time histories ***
!
 REAL(sp):: t_sl   (1:mel)  ! time Since Loading for lh=0
 REAL(sp):: t_su   (1:mel)  ! time Since Unloading for lh=1 & 2
 REAL(sp):: t_au   (1:mel)  ! time between loading and unloading for lh=2 & ...
!                           ! ... lenght of deglaciation for lh=3 ... &
!                           ! ... length of the loading phase for LH=8.
 REAL(sp):: t_sbd  (1:mel) ! time Since the Beginning of Deglaciation in lh=3
 REAL(sp):: t_sbld (1:mel) ! time Since the Beginning of Last Deglaciation in lh=4
 REAL(sp):: t_auc  (1:mel) ! length of a loading phase in lh=4
 REAL(sp):: d_inc  (1:mel) ! length of an un-loading phase in lh=4
 integer :: n_r    (1:mel) ! number of (loading+unloading) phases in addition to the
!                          ! las for lh=4
 REAL(sp):: t_stz  (1:mel) ! time Since t=zero in lh=5
 REAL(sp):: t_per  (1:mel) ! period of the sinusoid in lh=5
 REAL(sp):: t_sbplp(1:mel) ! time Since the Beginning of the Piecewise Linear
!                            Phase for Lh=6
 Integer(sp),parameter :: N_SS_MAX = 24  ! Max. number of linear phases
!                            within the piecewise linear phase for lh=6
 Real(sp) :: TMP (1:mel, 0:N_SS_MAX )   ! Times which mark the boundaries
!                            between adjacent linear phases in lh=6
 Real(sp) :: ALH (0:mel, 0:N_SS_MAX )   ! Thicknesses (or masses) correspondig to
!                             those times (lh=6)
 integer  :: n_ss (1:mel)   ! Number of linear phases which form the piecewise linear
!                             phase for lh=6
 Real(sp) :: T_NN (1:mel)   ! Lenght of the piecewise linear phase of LH=6.
!
 REAL(sp) :: t_sbpcp(1:mel) ! time Since the Beginning of the Piecewise Linear
!                             Phase for Lh=7 & 8
 Real(sp) :: D_ILT  (1:mel) ! Lenght of each time step  which form the piecewise
!                             constant phase for Lh=7 & 8
 integer  :: n_s (1:mel)    ! Number of steps which constitute the piecewise linear
!                             phase in Lh=7 & 8
 Real(sp) :: T_N (1:mel)    ! Lenght of the piecewise constant phase of LH=7 & 8.
!
 Integer(sp), parameter :: N_S_MAX  = 24 ! Max. number of steps within the
!                           ! pieceswise constant phase for  LH=7 & 8.
!
 Real(sp) :: alt    (1:mel, 0:N_S_MAX+1) ! Thicknesses (or masses) for each of the
!                           ! steps for lh=7 or 8.
 REAL(sp):: tau_8           ! length of the loading phase for the th7
!                             converted to th8.
!*** End of the block of statements useful for the time histories ***
!
 REAL(SP) :: load_history     ! Function which computes the TH
    EXTERNAL load_history
 REAL(SP) :: d_load_history   ! Function for the time-derivatives
    EXTERNAL d_load_history
!
  Real(sp) :: ajunk       ! a single precision junk   
  Real(sp) :: XR(0:18)    ! Helps reading file ICE3.DAT 
!
  REAL(DP) :: SSIGMA (1:mel,llmin:llmax) ! Spher. Arm. corfficients fon AX loads
!
  REAL(SP) :: T1, T2           ! Bounds on the time interval (kyr BP)
  REAL(SP) :: DTI              ! Time increment
  REAL(SP) :: GIVEN_TIME       ! A given time BP
!  
  REAL(DP) :: VECT(4), D_VECT(4)  ! 4-vector solution for a single load  
  REAL(DP) :: WECT(4), D_WECT(4)  ! 4-vector cumulative solution for many loads 
! 
  REAL(DP) :: HEIGHT      ! MAXIMUM height in time for a given load 
  REAL(DP) :: AMPLITUDE   ! Half-amplitude of the load (deg) 
  REAL(DP) :: MASS        ! MAXIMUM mass in time for a given load 
  REAL(DP) :: CEN(2)      ! Long. & Colat. of the centre (centroid) of the load  
  REAL(DP) :: OBS(2)      ! Long. & Colat for the observer 
  REAL(DP) :: DL, DT      ! Width in Long. & Colat. for a load of type 50
!
  REAL(SP) :: JUNK, junk3, junk4, junk5    ! A JUNK (sp)
  REAL(DP) :: JUNKD   ! A JUNK (DP) 
  real(dp) :: vjunkd(4) 
!
  REAL(SP) :: TFEA              ! A certain time feature for a given TH
  REAL(SP) :: TIME_BP           ! Time BP (kyr) 
  REAL(SP) :: LOCAL_TIME        ! LOCAL TIME for a given TH
!
  Real(sp) :: SCALE             ! Scaling factor for the Ice elements 
!
  Real(dp) :: Fact1, Fact2 ! Scaling factors for the Stokes Coefficients
  Real(dp) ::   UUU(2)     ! ==        (c_lm, s_lm) due to an aggregate
  Real(dp) :: D_UUU(2)     ! == (d/dt) (c_lm, s_lm) due to an aggregate
  Real(dp) ::   UU (2)     ! ==        (c_lm, s_lm) due to a single ice element
  Real(dp) :: D_UU (2)     ! == (d/dt) (c_lm, s_lm) due to a single ice element
!
  Real(dp) ::   Inerzia(3,3) ! (change of the) Inertia tensor 
  Real(dp) :: D_Inerzia(3,3) ! Time derivative of Inerzia(:,:)    
  Real(dp) ::       Ine(3,3) ! (change of the) Inertia tensor due to a single element 
  Real(dp) ::     D_Ine(3,3) ! Time-derivative of Ine(:,:) 
!
! _cpu time tools_
  REAL(sp) :: txt(2)       ! It was: REAL(sp) :: txt(2), etime; EXTERNAL ETIME 
  REAL(sp) :: ela4, ela1, ela0
!
! _observers coordinates_
 integer (i4b) :: nobs                         ! number of observers
 integer (i4b), parameter :: nobs_max=16200    ! maximum NOBS acceptes
 Real(dp) :: long_obs(nobs_max), cola_obs(nobs_max) ! observers long. & colat.
 Real(dp) :: L1, L2, C1, C2, EPL, EPC, CC, LL         ! maps bounds & increments
!
! _baselines declarations_
  integer (i4b), parameter :: max_sites= 300  ! Max. number of Couples of VLBI Sites
  Real(dp) :: Long_Site_1 (max_sites) ! Longitude of Site#1 
  Real(dp) :: Long_Site_2 (max_sites) ! Longitude of Site#2
  Real(dp) :: Cola_Site_1 (max_sites) ! Co-latitude of Site#1 
  Real(dp) :: Cola_Site_2 (max_sites) ! Co-latitude of Site#2
  Character*8 :: Name_1(max_sites), Name_2(max_sites)  ! A couple sites names
  CHARACTER*8 :: junk8                   ! a junk...
  REAL(DP):: D_WECT_1(4)   ! Solution vector (rad, cola, long, geoid) for Site#1
  REAL(DP):: D_WECT_2(4)   ! Solution vector (rad, cola, long, geoid) for Site#2
  REAL(DP):: OBS_1(2)      ! == (Long.,Colat.) for Observer #1
  REAL(DP):: OBS_2(2)      ! == (Long.,Colat.) for Observer #2
  REAL(DP):: D_VECT_1(4)   ! Solution at Site #1 due to a single element of the aggregate
  REAL(DP):: D_VECT_2(4)   ! Solution at Site #2 due to a single element of the aggregate
  REAL(DP):: L_RATE, T_RATE, V_RATE  ! L, T, & V rates
!
  Integer(i4b), parameter :: max_num_righe   = 20001 ! Max num. of rows in task_3.dat 
  Integer(i4b), parameter :: max_sub_agg     = 12    ! Max num. of subaggregates 
!  
  INTEGER(i4b), PARAMETER:: N_r_max=10  ! Max. number of ramps for TH=3 
  Integer(i4b) :: JLOC(max_sub_agg)     ! Number of elements of each sub-aggregate
  Integer(i4b) :: IICE                  ! Real Oceans Switch 
  Integer(i4b) :: C_L (1:mel)           ! Load type 
  Integer(i4b) :: L_H (1:mel)           ! Time-history type 
  Integer(i4b) :: IOC       (1:mel)     ! Real Oceans Switch for an ice element 
  Integer(i4b) :: type_inten            ! Type of local Study 
  Integer(i4b) :: type_esten            ! Type of global Study 
  Integer(i4b) :: FRTD        ! (0/1) Converts CL=50 to CL=10 (for ICE1 ed ICE2)
  Integer(i4b) :: F7T8        ! (0/1) Converte LH=7 to LH=8 (for ICE1, ICE2 & ICE3)
  Integer(i4b) :: KP                   ! Number fo subaggregates of ICE_*
  Integer(i4b) :: CR(0:18)             ! Ausilium for reading ICE1 & ICE2  
  Integer(i4b) :: Npul                 ! Punti in LONG. per mappe. 
  Integer(i4b) :: Npuc                 ! Punti in COLAT. per mappe. 
  Integer(i4b) :: Lde, Mde             ! Grado ed Ordine per Studi Estensivi 
  Integer(i4b) :: Lde1, Lde2           ! Gradi per studi estensivi
  Integer(I4B) :: CODE                 ! Codes to identify the model, given NV 
  Integer(I4B) :: ILM                  ! (0/1) An option for NV=7 and 9. 
  Integer(i4b) :: ijunk, file_type     ! Junks & indeces 
  Integer(i4b) :: I,J, KK, JJ, ITASK   !   "   "    "  
  INTEGER(I4B) :: JN, N, K             !   "   "    "
  Integer(i4b) :: l, m                 ! Junks & indeces
  Integer(i4b) :: Inorm                ! (0/1) Normalizzazione dei coeff. di stokes.  
  Integer(i4b) :: IGN0, IGN1, IGN2, IGN3  ! (0/1) Controlli per la formula del RSL
  Integer(i4b) :: n_study              ! Conteggio delle analisi
  INTEGER(i4b) :: iarmo, imake, iexte, iadhoc, iicest ! indices
  INTEGER(i4b) :: n_inte, n_este ! indices
  INTEGER(i4b) :: n_armo, n_make, n_exte, n_adhoc, n_iicest  ! indeces
  INTEGER(i4b) :: drop_modes           ! a switch
  INTEGER(i4b) :: iname(max_sub_agg)   ! checks for repetitions
  INTEGER(i4b) :: mel_eff              ! effectively counted ice elements
  INTEGER(i4b) :: ILB                  ! a baseline index
  INTEGER(i4b) :: I_THI                ! thickness switch for maps
!
  Character*30  :: filename           ! Filename
  Character*8   :: titre               ! Header 
  Character*30  :: keyword              ! A Keyword of Task#3 
  Character*8   :: Sub_Agg(max_sub_agg)  ! Name of a sub-aggregate of ICE*
  CHARACTER*30  :: FILE_SPARSI            ! A filename 
  CHARACTER*20  :: date, timc              ! date and time 
  CHARACTER*5   :: name_agg                 ! name of an aggregate
!  
!+++++++++++++++++++
end module FOR_TASK3
!+++++++++++++++++++
!
!
!



!
!
!
!
!  .......................... 
!  !                        !
 	program TABOO       !
!  !                        !
!  !........................!	
!
!
        use COMMON; implicit NONE 
!
        integer (i4b) :: oo(3)  
!	
        CHARACTER*6  ST
	CHARACTER*20 date, TIMC
!
        open(99,file='taboo.log',status='unknown')
!
!
!
        call DATE_AND_TIME (date,timc)      
        Write(99,*) '# ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
!
!
        Write(99,*)'Detecting the Active INPUT file !'
!
!
!
        OO(:)=0  
!
!
!
        OPEN(1,file='task_1.dat',status='unknown') 
	READ(1,'(A6)') ST
	IF( ST == 'Active' ) then
                          OO(1) = 1
                          Write(99,*)'task_1.dat is active '
                        else
                          Write(99,*)'task_1.dat is NOT active '
                        endif
        CLOSE(1)
!
!
        OPEN(1,file='task_2.dat',status='unknown') 
	READ(1,'(A6)') ST
	IF( ST == 'Active' ) then
                          OO(2) = 1
                          Write(99,*)'task_2.dat is active '
                        else
                          Write(99,*)'task_2.dat is NOT active '
                        endif
        CLOSE(1)	
!
        OPEN(1,file='task_3.dat',status='unknown') 
	READ(1,'(A6)') ST
	IF( ST == 'Active' ) then
                          OO(3) = 1
                          Write(99,*)'task_3.dat is active '
                        else
                          Write(99,*)'task_3.dat is NOT active '
                        endif
        CLOSE(1)
!
!
        IF(OO(1)+OO(2)+OO(3) == 0) Then 
	Write(99,*) 'ERROR in TABOO: No imput file is active '
	Write(99,*) '**** JOB ABORTED ********************** ';Stop
        Endif 
!
        IF(OO(1)+OO(2)+OO(3) >  1) Then 
	Write(99,*) 'ERROR in TABOO: The are more active input files  '
	Write(99,*) '**** JOB ABORTED ****************************** ';Stop
        Endif 	
!
!
!
        If (oo(1)==1) CALL TASK_1 
        If (oo(2)==1) CALL TASK_2 
        If (oo(3)==1) CALL TASK_3 
!
!
        call DATE_AND_TIME (date,timc) 
        Write(99,*) '# Closing this file (taboo.log) on', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
!
        close(99) 
!	
        STOP
!
!..........................
        end program TABOO
!..........................
!
!
!
!
!
!
!++++++++++++++++++++++++++++++
   SUBROUTINE TASK_3          !
!++++++++++++++++++++++++++++++
!
!
! o----------------------------------------o
! |                                        |
! |       Local & Global analyses for      |
! |                                        |
! |            Complex Ice Loads           |
! |                                        |
! o----------------------------------------o
!
!
!
!
 USE for_task3
 USE trigfcn
 implicit NONE
!
!
!
!    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
!    !  Beginning with the executable instructions !
!    ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
!
!
! # These integers are used to monitor the 
!     keyword sequence .... 
!
                     n_armo=0
                   n_make=0
                 n_exte=0
               n_adhoc=0
                 n_iicest=0
                   n_inte=0
                      n_este=0
                        iarmo=0
                        imake=0
                    iexte=0
                 iadhoc=0
               iicest=0
!
!
!
!
!
      Write(99,*) 'Opening file task_3.dat ...  '
!
      open(1,file='task_3.dat',status='unknown')
!
      Write(99,*) 'Reading file task_3.dat ...  '
! 
!
 DO 10001 ITASK = 1, max_num_righe      
!
!
!
!
!
   READ(1,'(a30)',END=10002) KEYWORD
!
!
!
!...............................................
   IF(KEYWORD(1:16) == 'Harmonic_Degrees') THEN
!...............................................
!
   iarmo=1
   n_armo = n_armo + 1
!
   READ(1,*)  LMIN, LMAX
!
   READ(1,*)  IV
!
     call m_3(30010) ! Checks the iv variable... 
!
     call m_3(30012) ! Harmonic_Degrees must appear once
!
        call DATE_AND_TIME (date,timc)
!
          IF    (IV==0) THEN
        Write(99,*) '# ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
          ELSEIF(IV==1) then
        Write(*,*) '# ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
        ENDIF
!
    IF(iv==1) then
     WRITE(*,*) '= = = = = = = = = = = = = = ='
     Write(*,*) '   Reading file task_3.dat   '
     WRITE(*,*) '= = = = = = = = = = = = = = ='
    ENDIF
!
                WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!

              Write(99,*) 'Lmin and Lmax ', Lmin, Lmax
    IF(IV==1) Write (*,*) 'Lmin and Lmax ', Lmin, Lmax
!
!
     call m_3(30020)  ! Checks l_min & l_max
!
    READ(1,*)  ONLY_ELASTIC
!
     call m_3(30030)  ! States the effect of only_elastic
!
    READ(1,*)  DROP_MODES 
!
     call m_3(30040)  ! States the effect of drop_modes
!

!..............................................
  ENDIF ! on KW Harmonic_Degrees
!..............................................
!
!
!
!
!
!
!* * * * * * * * * * * * * * * * * * * * * * 
   IF(KEYWORD(1:10) == 'Make_Model') THEN     
!* * * * * * * * * * * * * * * * * * * * * * 
!
   imake=1
   n_make = n_make + 1
!
   call m_3(30050) ! Harmonic_Degrees must be active
!
   call m_3(30052) ! Make_Model must appear once
!
!
                WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
                Write(99,*)    'Building the model '
    IF(IV==1)   Write(*, *)    'Building the model '
!
!
   READ(1,*) NV 
!  
   READ(1,*) CODE
!
   READ(1,*) LT 
!
   READ(1,*) ILM 
!      
!
   Write(99,*) 'Number of VE layers      ', NV  
   Write(99,*) 'Model CODE is            ', CODE 
!
IF(IV==1) then
   Write(*,*) 'Number of VE layers      ', NV
   Write(*,*) 'Model CODE is            ', CODE
ENDIF
!
!
    Write(99,'(a42,1x,f12.4,1x,a3)') &
    'From input, the lithospheric thickness is', LT,  ' km'
    Write(99,'(a47,1x,i3)')       &
    'The ILM parameter (for NV=7 and 9) is set to  ', ILM
!
   IF(IV==1) THEN
    Write(*, '(a42,1x,f12.4,1x,a3)') &
    'From input, the lithospheric thickness is', LT,  ' km'
    Write(*, '(a47,1x,i3)')       &
    'The ILM parameter (for NV=7 and 9) is set to  ', ILM
   ENDIF
!
!
               Write(99,*) 'Mantle viscosity from BOTTOM to TOP (/1E21) '
     IF(IV==1) WRITE(*,*)  'Mantle viscosity from BOTTOM to TOP (/1E21) '
   Do k=1, nv
      Read(1,*) vis(k)
	          Write(99, '(a19,i4,1x,a1F10.4)') &
                       'Viscosity of layer', k, '=', vis (k)
   IF(IV==1) Write(*,  '(a19,i4,1x,a1,F10.4)') &
                       'Viscosity of layer', k, '=', vis (k)
               vis(k) = vis(k) * 1.q21    	  
   Enddo 
! 
!
   CALL SPEC (Code, Ilm, LT)  ! Builds the Model ... 
!   
!
! +------------------------
       call SPECTRUM (1)      ! Computes the Spectral quantities        
! +------------------------
!
!
!+++++++++++++++++++++++++++
  If(only_elastic==0) then 
!+++++++++++++++++++++++++++  
!
!
            Write(99,*)'Writing the spectrum on file spectrum.dat'
 IF(IV==1)  Write(*,* )'Writing the spectrum on file spectrum.dat'
!
  open(3,file='spectrum.dat',status='unknown')
!
  call DATE_AND_TIME (date,timc)
!
  Write(3,*) '# File spectrum.dat, created by Task#3 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  
  Write(3,*) '# 1st column: l = Harmonic degree                '
  Write(3,*) '# 2nd   "   : LOG10 (l)                          '
  Write(3,*) '# 3th   "   : s, (kyrs**(-1))                    '
  Write(3,*) '# 4rd   "   : LOG10(-s), with s in kyrs**(-1)    '
  Write(3,*) '# 5th   "   : Relaxation time = -1000./s, (yrs)  '
  Write(3,*) '# 6th   "   : LOG10 (Relaxation time (yrs))      '
!
DO l=lmin, lmax 
! 
    write (3, '(/)')
      do k = 1, nroots
       write (3, '(i4,1x,5(e15.7,1x))') l,                  & 
                                        log10 (dfloat (l)), &
                                        s(l,k),            &
					log10(-s(l,k)),    &
				        (-1000.q0/s(l,k)), &
				        log10(-1000.q0/s(l,k))
      enddo
!
ENDDO; close(3)
!
!+++++++++++
   Endif 
!+++++++++++
!
          open(13,file='h.dat',status='unknown')       
          open(14,file='l.dat',status='unknown')
          open(15,file='k.dat',status='unknown')
!
!
  call DATE_AND_TIME (date,timc)
  Write(13,*) '# File h.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  call DATE_AND_TIME (date,timc)
  Write(14,*) '# File l.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  	      	  
  call DATE_AND_TIME (date,timc)
  Write(15,*) '# File k.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  	  
!
!
!
          if(only_elastic==0) then 
          WRITE(99,*) 'Writing the elastic, fluid, and v-elastic h ldc on h.dat '
          WRITE(99,*) 'Writing the elastic, fluid, and v-elastic l ldc on l.dat '
          WRITE(99,*) 'Writing the elastic, fluid, and v-elastic k ldc on k.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic, fluid, and v-elastic h ldc on h.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic, fluid, and v-elastic l ldc on l.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic, fluid, and v-elastic k ldc on k.dat '
          endif
!
          if(only_elastic==1) then 
          WRITE(99,*) 'Writing the elastic h ldc on h.dat '
          WRITE(99,*) 'Writing the elastic l ldc on l.dat '
          WRITE(99,*) 'Writing the elastic k ldc on k.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic h ldc on h.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic l ldc on l.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic k ldc on k.dat '
          endif
!	  
!
          DO l=lmin, lmax 
!
          if(only_elastic==0) then 
                    write (13, '((i3,1x,96(1x,e20.8)))') &
  	     l, h_e(l), h_f(l), (h_v (l,k), k = 1, nroots)
                    write (14, '((i3,1x,96(1x,e20.8)))') &
	     l, l_e(l), l_f(l), (l_v (l,k), k = 1, nroots)
                    write (15, '((i3,1x,96(1x,e20.8)))') &
	     l, k_e(l), k_f(l), (k_v (l,k), k = 1, nroots)
	  endif
          if(only_elastic==1) then 
                    write (13, '((i3,1x,96(1x,e20.8)))') &
  	     l, h_e(l) 
                    write (14, '((i3,1x,96(1x,e20.8)))') &
	     l, l_e(l) 
                    write (15, '((i3,1x,96(1x,e20.8)))') &
	     l, k_e(l) 
	  endif  
!
!++++++++++++++++++++++++++++++
    if(only_elastic==0) then 
!++++++++++++++++++++++++++++++ 
!
          if(drop_modes==0) then 
	                   do k=1, nroots
	                 r_h (l,k) = h_v(l,k)/s(l,k)      ! Normalized residue for h
                         r_l (l,k) = l_v(l,k)/s(l,k)      ! Normalized residue for l
                         r_k (l,k) = k_v(l,k)/s(l,k)      ! Normalized residue for k
                         tek (l,k) = -1.0/s(l,k)     
	                   enddo 
	  endif 
!
          if(drop_modes==1)  then
!	   
             do k=1, nroots
              r_h (l,k) = 0.0
              r_l (l,k) = 0.0
              r_k (l,k) = 0.0 
!
              IF(vec(l,k)==1._sp) then
!
                r_h (l,k) = h_v(l,k)/s(l,k)      ! Normalized residue for h
                r_l (l,k) = l_v(l,k)/s(l,k)      ! Normalized residue for l
                r_k (l,k) = k_v(l,k)/s(l,k)      ! Normalized residue for k
                tek (l,k) = -1.0/s(l,k)          ! Relaxation time in k-yrs
!
              endif
!
             enddo
!
             endif 
!
!+++++++++++++++++++ 
       ENDIF
!+++++++++++++++++++
!
          ENDDO
!
          close(13)
          CLOSE(14)
          CLOSE(15)
!
          WRITE(99,*) 'Done with the computation of the Spectral quantities '
IF(IV==1) WRITE(* ,*) 'Done with the computation of the Spectral quantities '
!
          WRITE(99,*) 'Looking now for KWs Ad_Hoc (or) Ice_*'
IF(IV==1) WRITE(* ,*) 'Looking now for KWs Ad_Hoc (or) Ice_*'
!
!
!
!* * * * * * * * * * * * * * * * * * 
   ENDIF  !  Endif on KW Make_Model
!* * * * * * * * * * * * * * * * * * 
!
!
!
!
!
! o o o o o o o o o o o o o o o o o o o o o o o 
   IF(KEYWORD(1:14) == 'External_Model') THEN
! o o o o o o o o o o o o o o o o o o o o o o o 
!
   iexte=1
!
   n_exte=n_exte+1
!
   rhoea = ext_rhoea
!
   call m_3(30060)  ! Checks if Make_Model is active
!
   call m_3(30062)  ! External_Model must appear once
!
               WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
   IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
   call m_3(30070)  ! States the model & a warning
!
    NV=1    
    NROOTS=3 
!    
 Open(67, file='external_spectrum.dat',status='unknown')   
   do J=0, llmax 
   do kk=1, nroots 
        read(67,'(i4,1x,2(e15.7,1x))',end=1901)    l,   s(l,kk) ,  s(l,kk)  
   enddo
   enddo
 1901 Close(67)
! stop 33411
!  
  Open(67, file='external_h.dat',status='unknown') 
          DO J=0, llmax 
             read (67, '((i3,1x,96(1x,e20.8)))', end=1902) &
	     l, h_e(l), h_f(l), (h_v (l,k), k = 1, nroots) 
	  ENDDO            
 1902 Close(67)   
  Open(68, file='external_l.dat',status='unknown')  
          DO J=0, llmax 
             read (68, '((i3,1x,96(1x,e20.8)))', end=1903) &
	     l, l_e(l), l_f(l), (l_v (l,k), k = 1, nroots)    
          ENDDO                             
 1903 close(68)                                     
      vec=1. 
!
!
!+++++++++++++++++++++++++++
  If(only_elastic==0) then 
!+++++++++++++++++++++++++++  
!
!
            Write(99,*)'Writing the spectrum on file spectrum.dat'
 IF(IV==1)  Write(*,* )'Writing the spectrum on file spectrum.dat'
!
 open(3,file='spectrum.dat',status='unknown')
!
  call DATE_AND_TIME (date,timc)
!
  Write(3,*) '# File spectrum.dat, created by Task#3 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  
  Write(3,*) '# 1st column: l = Harmonic degree                '
  Write(3,*) '# 2nd   "   : LOG10 (l)                          '
  Write(3,*) '# 3th   "   : s, (kyrs**(-1))                    '
  Write(3,*) '# 4rd   "   : LOG10(-s), with s in kyrs**(-1)    '
  Write(3,*) '# 5th   "   : Relaxation time = -1000./s, (yrs)  '
  Write(3,*) '# 6th   "   : LOG10 (Relaxation time (yrs))      '
!
!
DO l=lmin, lmax
!
! 
    write (3, '(/)')
      do k = 1, nroots
       write (3, '(i4,1x,5(e15.7,1x))') l,                  & 
                                        log10 (dfloat (l)), &
                                        s(l,k),            &
					log10(-s(l,k)),    &
				        (-1000.q0/s(l,k)), &
				        log10(-1000.q0/s(l,k))
      enddo
!
          if(drop_modes==0) then 
	                   do k=1, nroots
	                 r_h (l,k) = h_v(l,k)/s(l,k)      ! Normalized residue for h
                         r_l (l,k) = l_v(l,k)/s(l,k)      ! Normalized residue for l
                         r_k (l,k) = k_v(l,k)/s(l,k)      ! Normalized residue for k
                         tek (l,k) = -1.0/s(l,k)     
!
	                   enddo 
	  endif 
!
          if(drop_modes==1)  then
!	   
             do k=1, nroots
              r_h (l,k) = 0.0
              r_l (l,k) = 0.0
              r_k (l,k) = 0.0 
!
              IF(vec(l,k)==1._sp) then
!
                r_h (l,k) = h_v(l,k)/s(l,k)      ! Normalized residue for h
                r_l (l,k) = l_v(l,k)/s(l,k)      ! Normalized residue for l
                r_k (l,k) = k_v(l,k)/s(l,k)      ! Normalized residue for k
                tek (l,k) = -1.0/s(l,k)          ! Relaxation time in k-yrs
!
              endif
!
             enddo
!
             endif 
!
ENDDO; close(3)
!

!
!++++++++++
   Endif 
!++++++++++
!
!
!
!
!
          open(13,file='h.dat',status='unknown')
          open(14,file='l.dat',status='unknown')
!
  call DATE_AND_TIME (date,timc)
  Write(13,*) '# File h.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  call DATE_AND_TIME (date,timc)
  Write(14,*) '# File l.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  	      	  
!
          if(only_elastic==0) then 
          WRITE(99,*) 'Writing the elastic, fluid, and v-elastic h ldc on h.dat '
          WRITE(99,*) 'Writing the elastic, fluid, and v-elastic l ldc on l.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic, fluid, and v-elastic h ldc on h.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic, fluid, and v-elastic l ldc on l.dat '
          endif
!
          if(only_elastic==1) then 
          WRITE(99,*) 'Writing the elastic h ldc on h.dat '
          WRITE(99,*) 'Writing the elastic l ldc on l.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic h ldc on h.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic l ldc on l.dat '
          endif
!	  
          DO l=lmin, lmax 
!
          if(only_elastic==0) then 
                    write (13, '((i3,1x,96(1x,e20.8)))') &
  	     l, h_e(l), h_f(l), (h_v (l,k), k = 1, nroots)
                    write (14, '((i3,1x,96(1x,e20.8)))') &
	     l, l_e(l), l_f(l), (l_v (l,k), k = 1, nroots)
	  endif
          if(only_elastic==1) then 
                    write (13, '((i3,1x,96(1x,e20.8)))') &
  	     l, h_e(l) 
                    write (14, '((i3,1x,96(1x,e20.8)))') &
	     l, l_e(l) 
	  endif  
!
          ENDDO
!
          close(13)
          CLOSE(14)
!
!
          WRITE(99,*) 'Done with the computation of the Spectral quantities '
IF(IV==1) WRITE(* ,*) 'Done with the computation of the Spectral quantities '
!
          WRITE(99,*) 'Looking now for KWs Ad_Hoc (or) Ice_*'
IF(IV==1) WRITE(* ,*) 'Looking now for KWs Ad_Hoc (or) Ice_*'
!
!
! o o o o o o o o o o o o o o o o o o o o 
  ENDIF  ! Endif on KW External_Model
! o o o o o o o o o o o o o o o o o o o o         
!
!
!
!
!
! + + + + + + + + + + + + + + + + + +  
  IF (Keyword(1:6) == 'Ad_Hoc') THEN     
! + + + + + + + + + + + + + + + + + +  
!
  iadhoc = 1
!
  n_adhoc=n_adhoc+1
!
   call m_3(30075) ! One between Make_Model & External_Model must be active
   call m_3(30077) ! Ad_Hoc must appear once
!

                WRITE(99,*) '> found KEYWORD ', keyword
    IF(IV==1)   WRITE(* ,*) '> found KEYWORD ', keyword
!
!
  WRITE(99,*) 'Reading and verifying the parameters for the ''Ad_Hoc'' aggregate'
  IF(iv==1) &
  WRITE(* ,*) 'Reading and verifying the parameters for the ''Ad_Hoc'' aggregate'
!
!
  Read(1, '(a30 )') filename
!
!
  OPEN(8, file = filename, status='unknown')
!
                WRITE(99,*) 'Reading data from file: ', filename
    IF(IV==1)   WRITE(* ,*) 'Reading data from file: ', filename
!
!
!
  READ(8,*) Numero_di_elementi
!
  mel_eff = Numero_di_elementi
!
  call m_3(30080)  ! Check on Numero_di_elementi
!

!
!
!
!
!
!+---------------------------------+
  DO 666 J = 1, Numero_di_elementi !   # collecting the ice elements data
!+---------------------------------+
!
!
  READ(8, *) C_L(j), IOC(j), L_H(j)
!
!
  call m_3(30090)  ! Checks the range of c_l
  call m_3(30100)  ! Checks the range of ioc
  call m_3(30110)  ! Checks ioc for c_l=11 or 21 ( 1 element)
  call m_3(30120)  ! Checks ioc for c_l=11 or 21 (>1 element)
  call m_3(30130)  ! Checks the range of l_h
!
!
     WRITE(99,*) '+++++++++++++++'
     write(99,*) 'Element#', j
     WRITE(99,*) '+++++++++++++++'
     Write(99,*) 'Load  Type =',  C_L (j)
     Write(99,*) 'Real Ocean =',  IOC (j)
     Write(99,*) 'Type of TH =',  L_H (j)
!
    IF(iv==1) THEN
     WRITE(*,*)  '+++++++++++++++'
     write(*,*)  'Element#', j
     WRITE(*,*)  '+++++++++++++++'
     Write(*,*)  'Load  Type =',  C_L (j)
     Write(*,*)  'Real Ocean =',  IOC (j)
     Write(*,*)  'Type of TH =',  L_H (j)
    ENDIF
!
!
!
!
!
!
       If(c_l(j) <= 21) THEN
!
   		 Read(8,*) Long_c(j), Cola_c(j), Amplit(j)
		 Read(8,*) Altezz(j)
!
                 call m_3(30140) ! Checks the thickness Altezz(j)
                 call m_3(30145) ! Reports the parameters on taboo.log
                 call m_3(30150) ! A warning on Altezz(j) for l_h=6, 7, & 8
!
                        ENDIF
!
!
       If(c_l(j) == 30) THEN
!
                 Read(8,*) Long_c(j), Cola_c(j)
 		 Read(8,*) Grande(j)
!
                 call m_3(30160) ! Checks the thickness Altezz(j)
                 call m_3(30165) ! Reports the parameters on taboo.log
                 call m_3(30170) ! A warning on Grande(j) for l_h=6, 7, & 8
!
                        ENDIF
!
!
       If(c_l(j) == 50) THEN
!
                 Read(8,*) Long_c(j), Cola_c(j), D_long(j), D_cola(j)
		 Read(8,*) Altezz(j)
!
                 call m_3(30180) ! Checks the thickness Altezz(j)
                 call m_3(30185) ! Reports the parameters on taboo.log
                 call m_3(30190) ! A warning on Grande(j) for lh=6, 7, & 8
!
                        ENDIF
!
!
!
!
!
!
!  # L_h == 0: Instantaneous Unloading 
!    ---------------------------------
!
   If(L_H(j) == 0)      THEN
!
!
           WRITE(99,*) 'CHECKING THE USER INPUT for LH==0'
IF(iv==1)  WRITE(* ,*) 'CHECKING THE USER INPUT for LH==0'
!
!
   Read (8,*) T_sl(j)
!
                        ENDIF  !  On L_h == 0
!
!
!  # L_h == 1: Instantaneous Loading 
!    -------------------------------
!
   If(L_H(j) == 1)      THEN
!
           WRITE(99,*) 'CHECKING THE USER INPUT for LH==1'
IF(iv==1)  WRITE(* ,*) 'CHECKING THE USER INPUT for LH==1'
!
!
   Read (8,*) T_su(j)
!
                        ENDIF  !  On L_h == 1
!
!
!
!  # L_h == 2: Loading & Unloading 
!    -----------------------------
!
   If(L_H(j) == 2) THEN
!
           WRITE(99,*) 'CHECKING THE USER INPUT for LH==2'
IF(iv==1)  WRITE(* ,*) 'CHECKING THE USER INPUT for LH==2'
!
                       Read (8,*) T_su(j),   T_au(j)
!
                       call m_3(30200)  ! T_au must be > 0

                       ENDIF    ! On L_h == 2
!
!
!
!  # L_h == 3: Simple deglaciation 
!    -----------------------------
!
   IF(L_H(j) == 3) then
!
!
           WRITE(99,*) 'CHECKING THE USER INPUT for LH==3'
IF(iv==1)  WRITE(* ,*) 'CHECKING THE USER INPUT for LH==3'
!
                       READ (8,*) T_sbd(j),  T_au(j)
!
                       call m_3(30210)  ! T_au must be > 0
!
             WRITE(99,*) 'The deglaciation is ', T_au(j), ' kyrs long'
  IF(iv==1)  WRITE(* ,*) 'The deglaciation is ', T_au(j), ' kyrs long'
!
                       ENDIF   ! On L_h == 3
!
!
!
!  # L_h == 4: Saw tooth 
!    -------------------
!
   If(L_H(j) == 4) then
!
!
           WRITE(99,*) 'CHECKING THE USER INPUT for LH==4'
IF(iv==1)  WRITE(* ,*) 'CHECKING THE USER INPUT for LH==4'
!
!
                       Read (8,*) T_sbld(j), T_auc(j), D_inc(j), N_R(j)
!
                       call m_3(30220) ! It must be: 0<=n_r<=5
!
           WRITE(99,*)  &
           'There are', N_R(j)+1, 'glacial phases including the last'
 IF(iv==1) WRITE(* ,*)  &
           'There are', N_R(j)+1, 'glacial phases including the last'
!
                       call m_3(30230) ! T_auc must be > 0
!
                  WRITE(99,*)  'The loading phase is ', T_AUC(j), 'kyrs long'
        IF(iv==1) WRITE(* ,*) 'The loading phase is ', T_AUC(j), 'kyrs long'
!
                       call m_3(30240) ! D_inc must be > 0
!
                  WRITE(99,*)  'The un-loading phase is ',D_INC(j), 'kyrs long'
        IF(iv==1) WRITE(* ,*) 'The un-loading phase is ',D_INC(j), 'kyrs long'
!
!
                  ENDIF  ! on L_h == 4
!
!
!
!
!  # L_h == 5: Sinusoidal loading 
!    ----------------------------
!
   IF(L_H(j) == 5) then
!
           WRITE(99,*) 'CHECKING THE USER INPUT for LH==5'
IF(iv==1)  WRITE(* ,*) 'CHECKING THE USER INPUT for LH==5'
!
                      READ (8,*) T_stz(j),  T_per(j)
!
                      call m_3(30250) ! T_per must be > 0 kyrs
!
                  WRITE(99,*)  'The period is ', T_per(j) , 'kyrs'
        IF(iv==1) WRITE(* ,*) 'The period is ', T_per(j) , 'kyrs'
!
!
                  ENDIF   ! on L_h == 5
!
!
!
!
!
!  # L_h == 6  (piecewise linear)
!    ----------------------------
!
   IF(L_H(j) == 6) then
!
                   WRITE(99,*) 'CHECKING THE USER INPUT for LH==6'
        IF(iv==1)  WRITE(* ,*) 'CHECKING THE USER INPUT for LH==6'
!
!
                       Read (8,*) T_sbplp(j), N_ss(j)
!
                       Read (8,*) (Tmp(j,i), i=0, N_ss(j))
!
                        call m_3(30260)  ! n_ss(j) must be >=1
                        call m_3(30270)  ! n_ss must be <= max
                        call m_3(30280)  ! Tmp (j,0) must be = 0.0
                        call m_3(30290)  ! Tmp (j,i) must be > 0
                        call m_3(30300)  ! Tmp (j,i-1) must be < Tmp(j,i)
!
!
                                  T_nn(J) = Tmp(J,N_ss(j))

!
          Write(99,*) 'The piecewise linear phase is', T_nn(j), ' kyrs long'
IF(iv==1) Write(99,*) 'The piecewise linear phase is', T_nn(j), ' kyrs long'
!
!
!
!
                       Read (8,*) (Alh(j,i), i=0, N_ss(j))
!
!
                       call m_3(30320) ! Alh(j,i) must be > 0
!
!
           WRITE(99,*) 'Computing the maximum thickness (or mass) for LH==6'
IF(iv==1)  WRITE(* ,*) 'Computing the maximum thickness (or mass) for LH==6'
!
JUNK  = -1.e30
!
Do I = 0, N_ss(j)
!
IF(c_l(j)<=21.OR.c_l(j)==50)Write(99,*)'time (kyrs) & thickness (m)',Tmp(j,i),alh(j,i)
IF(c_l(j)==30              )Write(99,*)'time (kyrs) & mass (kg)'    ,Tmp(j,i),alh(j,i)
IF(iv==1) then
IF(c_l(j)<=21.OR.c_l(j)==50)Write(*,*)'time (kyrs) & thickness (m)',Tmp(j,i), alh(j,i)
IF(c_l(j)==30              )Write(*,*)'time (kyrs) & mass (kg)'    ,Tmp(j,i), alh(j,i)
endif
!
     If(Alh(j,i) >= JUNK) JUNK = Alh(j,i)
!
ENDDO
!
!
!
 IF(c_l(j)<=21.or.c_l(j)==50)Write(99,*)'Maximum thickness for l_h==6 (km) =',JUNK
IF((c_l(j)<=21.or.c_l(j)==50).and.iv==1)&
                             Write(*,*) 'Maximum thickness for l_h==6 (km) =',JUNK
!
  IF( c_l(j)==30)  Write(99,*) 'Maximum mass for l_h==6 (kg) =', JUNK
  IF( c_l(j)==30.and.iv==1)&
                   Write (*,*) 'Maximum mass for l_h==6 (kg) =', JUNK
!
!
!
                   call m_3(30330) ! the maximum thickness (mass) must be > 0
!
!
!
              WRITE(99,*)'Normalizing the thickness (mass) for l_h==6'
   IF(iv==1)  WRITE(*, *)'Normalizing the thickness (mass) for l_h==6'
!
IF( c_l(j)<=21.or.c_l(j)==50) Write(99,*) 'time(kyrs) & normalized thickness'
IF((c_l(j)<=21.or.c_l(j)==50).and.iv==1)&
              Write(*,*)  'time(kyrs) & normalized thickness'
!
  IF(c_l(j)==30)           &
              Write(99,*) 'time(kyrs) & normalized mass'
  IF(c_l(j)==30.and.iv==1)&
              Write(*,*)  'time(kyrs) & normalized mass'

!
!
   DO I =0, N_ss(j)
			   IF(JUNK  ==  0.0) ALH (J,I) = JUNK
			   If(JUNK  /=  0.0) ALH (J,I) = ALH(J,I)/JUNK
!
!
                  Write(99,*)  TMP(J,I), ALH(J,I)
       IF(iv==1)  Write(* ,*)  TMP(J,I), ALH(J,I)
!
   ENDDO
!
!
!
        If(C_L(j) <= 21 .or. C_L(j) == 50)    ALTEZZ(J) = JUNK
 	If(C_L(j) == 30)                      GRANDE(J) = JUNK
!
!
!
  ENDIF  ! on L_h == 6
!
!
!
!
!
!  # L_h == 7  (piecewise constant)
!    ------------------------------
!
   If(L_H(j) == 7) then
!
!
           WRITE(99,*) 'CHECKING THE USER INPUT for LH==7'
IF(iv==1)  WRITE(* ,*) 'CHECKING THE USER INPUT for LH==7'
!
!
!
             Read (8,*) T_sbpcp(j), D_ilt(j), N_s(j)
!
!
!
                        call m_3(30340) ! D_ilt(j) must be > 0
                        call m_3(30350) ! n_s(j) must be >= 1
                        call m_3(30360) ! n_s(j) msut be <= n_s_max
!
!
!
                        T_n(j) =  D_ilt (j) * Float(N_s(j))
!
!
!
           WRITE(99,*)'The piecewise constant phase is',T_n(j),'kyrs long for l_h=7'
  IF(iv==1)WRITE(* ,*)'The piecewise constant phase is',T_n(j),'kyrs long for l_h=7'
!
!
            WRITE(99,*) 'The duration of each time-step is', d_ilt(j), 'kyrs for l_h=7'
  IF(iv==1) WRITE(* ,*) 'The duration of each time-step is', d_ilt(j), 'kyrs for l_h=7'
!
!
!
             Read (8,*) (Alt(j,i), i=0, N_s(j))
!
!
!
                        call m_3(30370) ! Alt(j,0) must be > 0
                        call m_3(30380) ! Alt(j,i) must be > 0
!
!
!
           Alt(j,n_s(j)+1) = Alt(j,n_s(j))
!
!
           WRITE(99,*) 'Computing the maximum thickness (or mass) for LH==7'
IF(iv==1)  WRITE(* ,*) 'Computing the maximum thickness (or mass) for LH==7'
!
!
     JUNK  = -1.e30
!
!+++++++++++++++++++
  Do I =0, N_s(j)+1
!+++++++++++++++++++
!
     IF(c_l(j)<=21 .OR. c_l(j)==50)  Write(99,*) 'Step# & thickness (m)',i, alt(j,i)
     IF(c_l(j)==30                )  Write(99,*) 'Step# & mass (kg)'    ,i, alt(j,i)
!
     IF((c_l(j)<=21 .OR. c_l(j)==50) .and. iv==1)  &
                        Write(*,*) 'Step# & thickness (m)' ,i, alt(j,i)
     IF(c_l(j)==30                   .and. iv==1)  &
                        Write(*,*) 'Step# & thickness (m)' ,i, alt(j,i)
!
!
  If( ALT(J,I) >= JUNK) JUNK = ALT(J,I)
!
!
!+++++++
  ENDDO
!+++++++
!
!
!
!
  IF( c_l(j)<=21.or.c_l(j)==50) &
                 Write(99,*) 'Maximum thickness for l_h==7 (m) =',JUNK
  IF((c_l(j)<=21.or.c_l(j)==50).and.iv==1)&
                 Write(* ,*) 'Maximum thickness for l_h==7 (m) =',JUNK
!
  IF( c_l(j)==30)Write(99,*) 'Maximum mass for l_h==7 (kg) =', JUNK
  IF( c_l(j)==30.and.iv==1)&
                 Write (*,*) 'Maximum mass for l_h==7 (kg) =', JUNK
!
!
!
                 call m_3(30390) ! The maximum thickness (mass) must be > 0
!
!
!
              WRITE(99,*) 'Normalizing the thicknesses (or masses) for l_h==7'
   IF(iv==1)  WRITE(* ,*) 'Normalizing the thicknesses (or masses) for l_h==7'
!
  IF( c_l(j)<=21.or.c_l(j)==50) &
              Write(99,*) 'Step# & normalized thickness'
  IF((c_l(j)<=21.or.c_l(j)==50).and.iv==1)&
              Write(*,*)  'Step# & normalized  thickness'
!
  IF(c_l(j)==30)           &
              Write(99,*) 'Step# &  normalized mass'
  IF(c_l(j)==30.and.iv==1)&
              Write(*,*)  'Step# &  normalized mass'
!
!
!
!
!
  DO I = 0, n_s(j)+1
!
         Alt(j,I) = Alt(j,I)/JUNK
!
                    Write(99,*)  I, Alt(j,I)
         IF(iv==1)  Write(* ,*)  I, Alt(j,I)
!
  ENDDO
!
!
        If(C_L(j) <= 21 .or. C_L(j) == 50)    ALTEZZ(J) = JUNK
 	If(C_L(j) == 30)                      GRANDE(J) = JUNK
!
!
!
  ENDIF  ! on L_h == 7
!
!
!
!
!  # L_h == 8: Piecewise constant with loading phase   
!    ...............................................
!
   If(L_H(j) == 8) then
!
!
                   WRITE(99,*) 'CHECKING THE USER INPUT for LH==8'
        IF(iv==1)  WRITE(* ,*) 'CHECKING THE USER INPUT for LH==8'
!
!
                        Read (8,*) T_sbpcp(j), D_ilt(j), N_s(j), T_au(j)
!
!
                        call m_3(30400) ! T_au must be > 0
                        call m_3(30410) ! n_s must be >= 1
                        call m_3(30420) ! n_s must be <= n_s_max
!
!
                        T_n(j) =  D_ilt (j) * Float(N_s(j))
!
!
!
 WRITE(99,*) 'The piecewise constant phase is',T_n(j),'kyrs long for l_h=8'
  IF(iv==1)&
 WRITE(* ,*) 'The piecewise constant phase is',T_n(j),'kyrs long for l_h=8'
!
!
!
!
            WRITE(99,*) 'The duration of each time-step is', d_ilt(j), 'kyrs for l_h=8'
  IF(iv==1) WRITE(* ,*) 'The duration of each time-step is', d_ilt(j), 'kyrs for l_h=8'
!
!
             Read (8,*) (Alt(j,i), i=0, N_s(j))
!
!
             call m_3(30430) ! Alt(j,0) must be > 0
             call m_3(30440) ! Alt(j,i) must be > 0

!
           Alt(j,n_s(j)+1) = Alt(j,n_s(j))
!
!
           WRITE(99,*) 'Computing the maximum thickness (or mass) for LH==8'
IF(iv==1)  WRITE(* ,*) 'Computing the maximum thickness (or mass) for LH==8'
!
!
     JUNK  = -1.e30
!
!+++++++++++++++++++
  Do I =0, N_s(j)+1
!+++++++++++++++++++
!
     IF(c_l(j)<=21 .OR. c_l(j)==50)  Write(99,*) 'Step# & thickness (m)',i, alt(j,i)
     IF(c_l(j)==30                )  Write(99,*) 'Step# & mass (kg)'    ,i, alt(j,i)
!
     IF((c_l(j)<=21 .OR. c_l(j)==50) .and. iv==1)  &
                        Write(*,*) 'Step# & thickness (m)' ,i, alt(j,i)
     IF(c_l(j)==30                   .and. iv==1)  &
                        Write(*,*) 'Step# & thickness (m)' ,i, alt(j,i)
!
!
!
  If( ALT(J,I) >= JUNK) JUNK = ALT(J,I)
!
!
!+++++++
  ENDDO
!+++++++
!
!
!
!
  IF( c_l(j)<=21.or.c_l(j)==50) &
                 Write(99,*) 'Maximum thickness for l_h==8 (m) =',JUNK
  IF((c_l(j)<=21.or.c_l(j)==50).and.iv==1)&
                 Write(* ,*) 'Maximum thickness for l_h==8 (m) =',JUNK
!
  IF( c_l(j)==30) Write(99,*)'Maximum mass for l_h==8 (kg) =', JUNK
  IF( c_l(j)==30.and.iv==1)&
                  Write (*,*)'Maximum mass for l_h==8 (kg) =', JUNK
!
!
!
                  call m_3(30450) ! The maximum mass (thickness) must be > 0
!
!
!
              WRITE(99,*) 'Normalizing the thicknesses (or masses) for l_h==8'
  IF(iv==1)   WRITE(* ,*) 'Normalizing the thicknesses (or masses) for l_h==8'
!
  IF( c_l(j)<=21.or.c_l(j)==50) &
              Write(99,*) 'Step# & normalized thickness'
  IF((c_l(j)<=21.or.c_l(j)==50).and.iv==1)&
              Write(*,*)  'Step# & normalized  thickness'
!
  IF(c_l(j)==30)           &
              Write(99,*) 'Step# &  normalized mass'
  IF(c_l(j)==30.and.iv==1)&
              Write(*,*)  'Step# &  normalized mass'

!
!
!
!
  DO I = 0, n_s(j)+1
!
         Alt(j,I) = Alt(j,I)/JUNK
!
                    Write(99,*)  I, Alt(j,I)
         IF(iv==1)  Write(* ,*)  I, Alt(j,I)
!
  ENDDO
!
!
        If(C_L(j) <= 21 .or. C_L(j) == 50)    ALTEZZ(J) = JUNK
 	If(C_L(j) == 30)                      GRANDE(J) = JUNK
!
!

  ENDIF  ! on L_h == 8
!
!
!
!
666 CONTINUE
!
!

!
 ENDIF   !  on the kw 'Ad_Hoc'
!
!
!


!
!
! 
  IF (Keyword(1:5) == 'Ice_1' .OR. &         ! One of the kws ICE_* is active  
      Keyword(1:5) == 'Ice_2' .OR. &
      Keyword(1:5) == 'Ice_3')         THEN
!
      name_agg=keyword(1:5)
!
      iicest=1
!
      n_iicest=n_iicest+1
!
                call m_3(30456)  ! Make_Model OR External_Model must be active

                call m_3(30457)  ! Ad_Hoc & Ice_* are mutually exclusive
!
                call m_3(30458)  ! Ice_* must appear once
!
                WRITE(99,*) '> found KEYWORD   ', keyword
    IF(IV==1)   WRITE(*, *) '> found KEYWORD   ', keyword
!
!
!
!
    IF( name_agg == 'Ice_1') then
                WRITE(99,*) 'Reading and verifying the ICE1 parameters'
    IF(IV==1)   WRITE(*, *) 'Reading and verifying the ICE1 parameters'
    endif
    IF( name_agg == 'Ice_2') then
                WRITE(99,*) 'Reading and verifying the ICE2 parameters'
    IF(IV==1)   WRITE(*, *) 'Reading and verifying the ICE2 parameters'
    endif
    IF( name_agg == 'Ice_3') then
                WRITE(99,*) 'Reading and verifying the ICE3 parameters'
    IF(IV==1)   WRITE(*, *) 'Reading and verifying the ICE3 parameters'
    endif

!
!+++++++++++++++++
    READ(1,*) IICE            ! (0/1) Ocean switch
!+++++++++++++++++
              call m_3(30470) ! IICE = ocean switch must be =0 or =1
!
    IF(iice==1) then
         WRITE(99,*)'The ocean correction is set: IOC=1 for ALL of the ice elements'
IF(iv==1)WRITE(*, *)'The ocean correction is set: IOC=1 for ALL of the ice elements'
    ELSE
         WRITE(99,*)'The ocean correction is NOT set: IOC=0 for ALL of the ice elements'
IF(iv==1)WRITE(*, *)'The ocean correction is NOT set: IOC=0 for ALL of the ice elements'
    ENDIF
!
!
!+++++++++++++++++
    READ(1,*) FRTD            ! (0/1) Converts from rect (c_l=50) to disk (10)
!+++++++++++++++++
              call m_3(30480) ! FRTD must be =0 or =1
!
    IF(FRTD==1) then
              WRITE(99,*) 'Conversion from cl=50 to 10 (only for ICE1 & ICE2) = YES'
    IF(iv==1) WRITE(*, *) 'Conversion from cl=50 to 10 (only for ICE1 & ICE2) = YES'

!
    IF(iv==1) then
  WRITE(*, *)'Converting the Rectangles (CL=50) to discs (CL=10)'
  WRITE(*, *)'Thickness, mass, and centroid of the element remains the same...'
  WRITE(*, *)'... only the shape is changed !!!'
  WRITE(*, *)'Details of the conversion are on taboo.log'
    ENDIF
  WRITE(99,*)'Converting the Rectangles (CL=50) to discs (CL=10)'
  WRITE(99,*)'Thickness, mass, and centroid of the element remains the same...'
  WRITE(99,*)'... only the shape is changed !!!'
  WRITE(99,*)'Details of the conversion are on taboo.log'
!
    ENDIF
!
    IF (FRTD==1 .AND. name_agg=='Ice_3') call m_3(30570) ! warning for Ice_3
!
!
!++++++++++++++++++++++++
    READ(1,*) F7T8, TAU_8     ! (0/1) Converts from lh=7 to lh=8, & loading phase length
!++++++++++++++++++++++++
              call m_3(30500) ! F7T8 must be =0 or =1
              call m_3(30510) ! TAU_8 must be > 0
!
    IF(F7T8==1) then
                         WRITE(99,*) 'Conversion from lh=7 to 8 = YES'
               IF(iv==1) WRITE(*, *) 'Conversion from lh=7 to 8 = YES'
                         WRITE(99,*) 'The Loading phase is',tau_8,'kyrs long'
               IF(iv==1) WRITE(*, *) 'The Loading phase is',tau_8,'kyrs long'
    ENDIF

!++++++++++++++++++
    READ(1,*) Scale           ! Scales the heights of the elements
!++++++++++++++++++
              call m_3(30520) ! SCALE must be > 0
!
                         WRITE(99,*) 'The scaling factor is', SCALE
               IF(iv==1) WRITE(*, *) 'The scaling factor is', SCALE
!
!+++++++++++++++
    READ(1,*) KP              ! Number of sub-aggregates of ICE* to be loaded
!+++++++++++++++
              call m_3(30530) ! Check bounds on kp
!
!
!
!
!
!
!
! # Filters on the list of sub-aggregates for ICE*:
!
                  WRITE(99,*) 'checking the list of sub-aggregates '
        IF(iv==1) WRITE(*, *) 'checking the list of sub-aggregates '
!
! +++++++++++++++++++++++++++
!
  DO  KK = 1, KP                ! Recognition of the sub-aggregates names
    READ (1,'(a8)') Sub_Agg(kk)
    IF(Sub_Agg(kk)=='ice1.dat'.AND.kp>=2) call m_3(30550)  ! ice1.dat must stay alone
    IF(Sub_Agg(kk)=='ice2.dat'.AND.kp>=2) call m_3(30550)  ! ice2.dat must stay alone
    IF(Sub_Agg(kk)=='ice3.dat'.AND.kp>=2) call m_3(30550)  ! ice3.dat must stay alone
  ENDDO
!
  IF( name_agg=='Ice_1') CALL M_3(30540) ! The sub-agg. names must be spelled correctly
  IF( name_agg=='Ice_2') CALL M_3(30541) ! The sub-agg. names must be spelled correctly
  IF( name_agg=='Ice_3') CALL M_3(30542) ! The sub-agg. names must be spelled correctly
!
  DO KK=1, kp
    INAME(kk) = 0
      DO JJ=1, kp
        if(SUB_AGG(KK) == SUB_AGG(JJ)) INAME(kk) = INAME(kk) + 1
      ENDDO
    IF(INAME(kk)>=2) CALL M_3(30560) ! A repetition is found
  ENDDO
!
! +++++++++++++++++++++++++++
!
        IF(kp==1) THEN
                  WRITE(99,*) KP,'sub-aggregate will be loaded. It is:'
        IF(iv==1) WRITE(*, *) KP,'sub-aggregate will be loaded. It is:'
                  else
                  WRITE(99,*) KP,'sub-aggregate(s) will be loaded. They are: '
        IF(iv==1) WRITE(*, *) KP,'sub-aggregate(s) will be loaded. They are: '
                  ENDIF
  do kk=1, kp
                  WRITE(99,*) kk,'o--> ',sub_agg(kk)
     IF(iv==1)    WRITE(* ,*) kk,'o--> ',sub_agg(kk)
  end do
!
!
!
!
!
!
!
!
  mel_eff = 0              ! Counts the TOTAL number of elements
!
!
! ------------------+
  DO 7070 KK = 1, KP       ! Do-loop on the various sub- aggregates
! ------------------+
!
!
     OPEN(10,file= Sub_Agg(kk), status='OLD')
!
               Write(99,*) 'Reading file ',Sub_Agg(kk)
     IF(iv==1) Write(* ,*) 'Reading file ',Sub_Agg(kk)
!
!
      Read (10,'(a8)') Titre  ! File header
      Read (10,'(a8)') Titre  ! File header
!
!
       Jloc(kk) = 0             ! Counts the elements of the sub-aggregate
!
! +------------------------+
       Do 8080 JJ = 1, mel      ! do-loop the elements of the sub-aggregate
! +------------------------+
!
!
   IF( name_agg=='Ice_1') THEN
!
!
!
       Read (10,*,end=9090) junk3, junk4, &
                            cr(18), cr(14), cr(12), cr(10), cr(8), cr(6)
!
       mel_eff = mel_eff + 1    ! Counts the TOTAL number of elements
!
       j=mel_eff
!
       Cola_c(j)=junk3; Long_c(j)=junk4 
!       
       C_L (J) = 50             ! In ICE1, the load is 'rectangular' 
       Ioc(J)  = IICE            ! Realistic ocean correction (0/1)
       L_H (J)   = 7              ! In ICE1, the DEFAULT TH is by steps (l_h=7)
       T_SBPCP(J)  = 18.           ! Today is 18 kyrs since the beginning of the PCP
       D_long(J)    = 5.            ! The rectangle is 5 degrees wide in longitude
       D_cola(J)    = 5.             ! The rectangle is 5 degrees wide in co-latitude
       N_S(J)       = 9               ! There are 9 time steps for l_h=7 during the PCP
       D_ilt(J)     = 2.              ! Each is 2 kyrs long
       T_N(j)       = n_s(j)*d_ilt(j) ! Duration of the PCP
!
       Jloc(kk) = Jloc(kk) + 1        ! Counts the elements of the sub-aggregate
!
       cr(16) =  cr(14)    ! Ice thickness at 16 kyr BP = that of 14 kyrs BP
       cr(4)  =  cr(6)      ! From 6 kyrs BP to present the ...
       cr(2)  =  cr(4)       !  ... ice thickness
       cr(0)  =  cr(2)        !   ... remains constant in ICE1
!
!
       alt(j,0) = scale*cr(18)         ! Ice thickness before the PCP of ICE1 (step#0)
!
       do i=1, n_s(j)                  ! According to my conventions ...
         alt(J,I) = scale*cr(18-2*I)   ! Alt(j,1) = cr(16), ...,ALT(J,9)= cr(0)
       enddo                           ! Alt(j,i) = thickness at step#i
!
       alt(J,n_s(j)+1) = alt(J,n_s(j)) ! Ice thickness at step#i+1
!
   ENDIF ! on 'Ice_1'
!
!
   IF( name_agg=='Ice_2') THEN
!
!
       Read (10,*,end=9090) junk3, junk4, (cr(i),i = 18,4,-1)
!
       mel_eff = mel_eff + 1                ! Counts the TOTAL number of elements
!
       j=mel_eff
!
       Cola_c(j)=junk3; Long_c(j)=junk4 
!
       C_L (J) = 50             ! In ICE2, the load is 'rectangular'
       Ioc(J)  = IICE            ! Realistic ocean correction (0/1)
       L_H (J)   = 7              ! In ICE2, the DEFAULT TH is by steps (#7)
       T_SBPCP(J)  = 18.           ! Today is 18 kyrs since the beginning of the PCP
       D_long(J)    = 5.            ! The rectangle is 5 degrees wide in longitude
       D_cola(J)    = 5.             ! The rectangle is 5 degrees wide in co-latitude
       N_S(J)       = 18              ! There are 9 time steps for l_h=7 during the PCP
       D_ilt(J)     = 1.              ! Each is 1 kyrs long
       T_N(j)       = n_s(j)*d_ilt(j) ! Duration of the PCP
!
       Jloc(kk) = Jloc(kk) + 1  ! Counts the elements of the sub-aggregate
!
	cr(3)  =  cr(4)     ! From 4 kyrs BP to present the ...
	cr(2)  =  cr(3)     ! ... ice thickness
	cr(1)  =  cr(2)     !   ... remains constant in ICE2
	cr(0)  =  cr(1)
!
       alt(j,0) = scale*cr(18)         ! Ice thickness before the PCP of ICE2 (step#0)
!
       do i=1, n_s(j)                  ! According to my conventions ...
         alt(J,I) = scale*cr(18-I)     ! Alt(j,1) = cr(17), ...,ALT(J,18)= cr(0)
       enddo                           ! Alt(j,i) = thickness at step#i
!
       alt(J,n_s(j)+1) = alt(J,n_s(j)) ! Ice thickness at step#i+1
!
!
   ENDIF ! on 'Ice_2'
!
!
   IF( name_agg=='Ice_3') THEN
!

       Read (10,1111,END=9090) ijunk, junk3, junk4, junk5, ajunk, &
                               (xr(i),i = 18,11,-1)         ! 1st row
       Read (10,1112,END=9090) (xr(i),i = 10,0 ,-1)         ! 2nd row
1111   format(i4,2F9.3,2F6.3,1X,8F5.0)
1112   format(15F5.0) 
!
!
       mel_eff = mel_eff + 1   ! Counts the TOTAL number of elements
!
       j=mel_eff
!
!
       Cola_c(j)=junk3; Long_c(j)=junk4; amplit(j)=junk5  
!
       C_L (J) = 10             ! In ICE3, the load is a disc
       Ioc(J)  = IICE            ! Realistic ocean correction (0/1)
       L_H (J)   = 7              ! In ICE3, the DEFAULT TH is by steps (#7)
       T_SBPCP(J)  = 18.           ! Today is 18 kyrs since the beginning of the PCP
       N_S(J)       = 18              ! There are 9 time steps for l_h=7 during the PCP
       D_ilt(J)     = 1.              ! Each is 1 kyrs long
       T_N(j)       = n_s(j)*d_ilt(j) ! Duration of the PCP
!
       Jloc(kk) = Jloc(kk) + 1  ! Counts the elements of the sub-aggregate
!
       alt(j,0) = scale*xr(18)         ! Ice thickness before the PCP of ICE3 (step#0)
!
       do i=1, n_s(j)                  ! According to my conventions ...
         alt(J,I) = scale*xr(18-I)     ! Alt(j,1) = cr(17), ...,ALT(J,18)= cr(0)
       enddo                           ! Alt(j,i) = thickness at step#i
!
       alt(J,n_s(j)+1) = alt(J,n_s(j)) ! Ice thickness at step#i+1
!
!
  ENDIF  ! on 'Ice_3'
!
!++++++++++++
8080 CONTINUE     ! do-loop on the elements of the sub-aggregate
!++++++++++++
     CLOSE(10)    ! closes the sub-aggregate file
!
9090 CONTINUE
!
           WRITE(99,*) 'Read ', Jloc(kk), 'ice elements from file ', Sub_agg(kk)
IF(iv==1)  WRITE(* ,*) 'Read ', Jloc(kk), 'ice elements from file ', Sub_agg(kk)
!
!
!
!++++++++++++
7070 CONTINUE     ! do-loop on the sub-aggregates
!++++++++++++
!
!
!
!mel_eff = mel_eff - kp !1
!
!
!
WRITE(99,*)'Read ',mel_eff,'ice elements from the sub-aggregate(s) of ', name_agg
IF(iv==1) &
WRITE(*, *)'Read ',mel_eff,'ice elements from the sub-aggregate(s) of ', name_agg
!
          WRITE(99,*)'I have succesfully collected the ', name_agg, ' data'
IF(iv==1) WRITE(*, *)'I have succesfully collected the ', name_agg, ' data'
!
write (99,*)'The detail are given in the following ...'
IF(iv==1) &
WRITE (*, *)'I am going to report the details on file taboo.log'
!
!
!
!
!
!++++++++++++++
DO J=1, MEL_EFF    ! do-loop on the elements ...
!++++++++++++++
!
        WRITE(99,*) '++++++++++++++++++++++++++++'
        Write(99,*) 'Aggregate  ',name_agg, '  Element # ',J
        WRITE(99,*) '++++++++++++++++++++++++++++'
!
!
!
!
 IF(ioc(j)==1) WRITE(99,*) 'The element is compensated on a realistic ocean'
 IF(ioc(j)==0) WRITE(99,*) 'The element is NOT compensated on a realistic ocean'
!
 IF(FRTD==0)then
 WRITE(99,*) 'The element is of type CL = ',c_l(j)
 ELSE
 WRITE(99,*) 'Before conversion, the element is of type CL = ',c_l(j)
 ENDIF
!
 IF(F7T8==0)THEN
 WRITE(99,*) 'The element time-history is LH =', l_h(j)
 else
 WRITE(99,*) 'Before conversion, the element time-history is LH =', l_h(j)
 ENDIF
!
!
!
 IF(c_l(j)==10) THEN
   Write(99,*)'Longitude of the centre (deg)=', long_c(j)
   Write(99,*)'Colatitude of the centre (deg)=', cola_c(j)
   Write(99,*)'Half-amplitude of the disc (deg)=', Amplit(j)
 ENDIF
 IF(c_l(j)==50) THEN
   Write(99,*)'Longitude of the centroid (deg)=', long_c(j)
   Write(99,*)'Colatitude of the centroid (deg)=', cola_c(j)
 ENDIF
 IF(c_l(j)==50 .AND. FRTD==1) THEN
   Write(99,*)'Before conversion, the element width in longitude (deg)=', d_long(j)
   Write(99,*)'Before conversion, the lement width in colatitude (deg)=', d_cola(j)
 ENDIF
 IF(c_l(j)==50 .AND. FRTD==0) THEN
   Write(99,*)'Element width in longitude (deg)=', d_long(j)
   Write(99,*)'Element width in colatitude (deg)=', d_cola(j)
 ENDIF
!
!
IF (FRTD==1.AND.(name_agg=='Ice_1'.or.name_agg=='Ice_2')) THEN   
!
       C_L (j) = 10
!
       Amplit (j)    = acos( 1.0 - sind(cola_c(j))*sind(d_cola(j)/2.)* &
                              d_long(j)*pi/180.0/pi)   * 180.0 / pi   ! It was 'acosd'
!
   Write (99,*) &
            'After conversion, the element is of type CL = ',C_L(j)
   Write (99,*) &
            'After conversion, the half-amplitude is (deg)',amplit(j)
!
ENDIF ! --------------------------------------    Endif on FRTD==1
!
!
!
IF (F7T8==1) THEN

      L_H (J)   = 8
      T_AU(j)   = tau_8
!
   Write (99,*) &
            'After conversion, the time history is of type = ',l_h(j)
   write (99,*) 'The Loading phase is',t_au(j),'kyrs long'
!
ENDIF ! --------------------------------------    Endif on F7T8==1
!
!
   IF(F7T8==0) THEN
   Write(99,*) &
            'Element thickness before the PCP (step#0) (m)=',alt(j,0)
               ELSE
   Write(99,*) &
            'Element thickness at the end of loading (step#0) (m)=',alt(j,0)
   ENDIF
!
!
   WRITE(99,*) 'Step# & element thickness (m)'

   IF(name_agg=='Ice_3'.or.name_agg=='Ice_2') THEN
      WRITE(99,'(6(1x,i2,1x,f9.4))') (i, alt(j,i),i=1 , 6)
      WRITE(99,'(6(1x,i2,1x,f9.4))') (i, alt(j,i),i=7 ,12)
      WRITE(99,'(6(1x,i2,1x,f9.4))') (i, alt(j,i),i=13,18)
                       ELSE
      WRITE(99,'(3(1x,i2,1x,f9.4))') (i, alt(j,i),i=1,3)
      WRITE(99,'(3(1x,i2,1x,f9.4))') (i, alt(j,i),i=4,6)
      WRITE(99,'(3(1x,i2,1x,f9.4))') (i, alt(j,i),i=7,9)
                      ENDIF
!
   Write(99,*) &
            'Element thickness after the PCP (step#ns+1) (m)=',alt(j,n_s(j)+1)
!
!
!
!
        junk=alt(j,0)
!
        do i=1, n_s(j)+1
           if(alt(j,i) >= junk) junk=alt(j,i)
	enddo
!
        do i=0, n_s(j)+1
           alt(j,i) = alt(j,i)/junk
	enddo
!
        WRITE(99,*) 'Maximum thickness of the element ', junk
!
        Altezz(j) = junk 
!
   IF(F7T8==0) THEN
   Write(99,*) &
   'NORMALIZED element thickness before the PCP (step#0) = ',alt(j,0)
               ELSE
   Write(99,*) &
   'NORMALIZED element thickness at the end of loading (step#0) = ',alt(j,0)
   ENDIF
!
!
   WRITE(99,*) 'Step# & NORMALIZED element thickness'

   IF(name_agg=='Ice_3'.or.name_agg=='Ice_2') THEN
      WRITE(99,'(6(1x,i2,1x,f9.4))') (i, alt(j,i),i=1 , 6)
      WRITE(99,'(6(1x,i2,1x,f9.4))') (i, alt(j,i),i=7 ,12)
      WRITE(99,'(6(1x,i2,1x,f9.4))') (i, alt(j,i),i=13,18)
                       ELSE
      WRITE(99,'(3(1x,i2,1x,f9.4))') (i, alt(j,i),i=1,3)
      WRITE(99,'(3(1x,i2,1x,f9.4))') (i, alt(j,i),i=4,6)
      WRITE(99,'(3(1x,i2,1x,f9.4))') (i, alt(j,i),i=7,9)
                                              ENDIF
!
   Write(99,*) &
                'NORMALIZED Element thickness after the PCP (step#ns+1) (m)=',&
                alt(j,n_s(j)+1)
!
!
!
! <-------------- end of data acquisition & normalization --------------->
!
!
!+++++++
  ENDDO   ! enddo on the ice elements
!+++++++
!
!
!
!
ENDIF  ! on kw I_*
!
!
!
!
!
!
!-------------------------------------------
   IF(KEYWORD(1:11)=='Local_Study') THEN
!-------------------------------------------
!
   n_inte = n_inte+1
!
          call m_3(30580) ! Ad_Hoc or Ice_* must be active
!
          call m_3(30582) ! Local_Study must appear once
!
                WRITE(99,*) '> found KEYWORD  ', keyword
    IF(IV==1)   WRITE(* ,*) '> found KEYWORD  ', keyword
!
!
!
!
    READ(1,*) type_inten   
!
!
                call m_3(30590) ! Cheking bounds for type_inten
!
!
!
!++++++++++++++++++++++++
   IF(type_inten==1) THEN
!++++++++++++++++++++++++
!
!
            Write(99,*) 'Local Study of Type#1' 
 If (iv==0) Write(99,*) 'Local Study of Type#1' 
 If (iv==1) Write(* ,*) 'Local Study of Type#1' 
!
!
!
    READ(1,*) IR, IT, IG, dIR, DIT, dIG
!
!
                call m_3(30600) ! One of the switches must be == 1
!
!
              write(99,*) 'The following variables will be computed and printed:'
     if(iv==1)write(* ,*) 'The following variables will be computed and printed:'
! 
if( IR ==1) write(99,*) '-> Radial displacement'
if( it ==1) write(99,*) '-> Tangential displacement'
if( IG ==1) write(99,*) '-> Geoid Heigth'
if( dIR==1) write(99,*) '-> Radial velocity'
if( dit==1) write(99,*) '-> Tangential velocity'
if( dIG==1) write(99,*) '-> Rate of change of the geoid heigth'
!
if(iv==1 .and. IR ==1) write(*,*)  '-> Radial displacement'
if(iv==1 .and. it ==1) write(*,*)  '-> Tangential displacement'
if(iv==1 .and. IG ==1) write(*,*)  '-> Geoid Heigth'
if(iv==1 .and. dIR==1) write(*,*)  '-> Radial velocity'
if(iv==1 .and. dit==1) write(*,*)  '-> Tangential velocity'
if(iv==1 .and. dIG==1) write(*,*)  '-> Rate of change of the geoid heigth'
!
!
!
     READ(1,*) LONG_OBS(1), COLA_OBS(1)
!
     kk=1; call m_3(30625) ! Cheks the bounds on long. & colat.
!
!
  Write(99,*) 'Longitude of the observer (deg) =', LONG_OBS(1)
  Write(99,*) 'Colatitude of the observer (deg) =', COLA_OBS(1)     
 IF (iv==1) THEN
  Write(* ,*) 'Longitude of the observer (deg) =', LONG_OBS(1)
  Write(* ,*) 'Colatitude of the observer (deg) =', COLA_OBS(1)  
 ENDIF
!     
!
     READ(1,*) T1, T2, DTI
!
!
                   call m_3 (30610) ! T1 must be <= T2
                   call m_3 (30620) ! DTI must be > 0.0
!
     NOBS=1
!
!
!+++++++++++++++++++++++++++++++++++++++++
     ENDIF  ! endif on Local_Study == 1
!+++++++++++++++++++++++++++++++++++++++++
!
!
!
!
!
!
!
!
!
!רררררררררררררררררררררררר
   IF(type_inten==2) THEN
!רררררררררררררררררררררררר
!
!
            Write(99,*) 'Local Study of Type#2'
 If (iv==1) Write(* ,*) 'Local Study of Type#2'
!
!
!
    READ(1,*) IR, IT, IG, DIR, DIT, DIG
!
!
                call m_3(30600) ! One of the switches must be == 1
!
!
              write(99,*) 'The following variables will be computed and printed:'
     if(iv==1)write(* ,*) 'The following variables will be computed and printed:'
! 
if( IR ==1) write(99,*) '-> Radial displacement'
if( it ==1) write(99,*) '-> Tangential displacement'
if( IG ==1) write(99,*) '-> Geoid Heigth'
if( dIR==1) write(99,*) '-> Radial velocity'
if( dit==1) write(99,*) '-> Tangential velocity'
if( dIG==1) write(99,*) '-> Rate of change of the geoid heigth'
!
if(iv==1 .and. IR ==1) write(*,*)  '-> Radial displacement'
if(iv==1 .and. it ==1) write(*,*)  '-> Tangential displacement'
if(iv==1 .and. IG ==1) write(*,*)  '-> Geoid Heigth'
if(iv==1 .and. dIR==1) write(*,*)  '-> Radial velocity'
if(iv==1 .and. dit==1) write(*,*)  '-> Tangential velocity'
if(iv==1 .and. dIG==1) write(*,*)  '-> Rate of change of the geoid heigth'
!
!
!      			                   
     READ(1,'(A30)') file_sparsi
!
     j=0
     OPEN(2,file=file_sparsi,ERR=121,status='old')
     j=1
!
 121 IF(j==0) call m_3(30680) ! Error opening file_sparsi
!
!
            Write(99,*) 'Reading the coordinates of the observers from file ',file_sparsi  
 If (iv==1) Write(* ,*) 'Reading the coordinates of the observers from file ',file_sparsi 
!
     K=0 
     DO J=1, NOBS_MAX  
       READ(2,*,END=404) LL, CC 
       K=K+1
       LONG_OBS(K) = LL
       COLA_OBS(K) = CC
!
       kk=k; call m_3(30625) ! Cheks the bounds on long. & colat.
!
     ENDDO
!
 404 NOBS=K
     CLOSE(2)
!
       call m_3(30630) ! Too many observers
       call m_3(30640) ! Error on the user-supplied file

            Write(99,*) 'Read ', NOBS, ' observers from file ', file_sparsi  
 If (iv==1) Write(* ,*) 'Read ', NOBS, ' observers from file ', file_sparsi  
!
!
!
DO J=1, NOBS
 Write(99,'(a35,i4,1x,a7,1x,2(f10.4,1x))') 'Longitude & colatitude of observer', j, &
             '(deg) =',LONG_OBS(j),COLA_OBS(j)
 IF (iv==1) THEN
 Write(* ,'(a35,i4,1x,a7,1x,2(f10.4,1x))') 'Longitude & colatitude of observer', j, &
             '(deg) =',LONG_OBS(j),COLA_OBS(j)
 ENDIF
ENDDO
!     
!
     READ(1,*) T1, T2, DTI
!
!
                   call m_3 (30610) ! T1 must be <= T2
                   call m_3 (30620) ! DTI must be > 0.0
!
!
!
!
!+++++++++++++++++++++++++++++++++++++++++
     ENDIF  ! endif on Local_Study == 2
!+++++++++++++++++++++++++++++++++++++++++

!
!
!
!
!
!
!
!
!
!+-+-+-+-+-+-+-+-+-+-+-+-
   IF(type_inten==3) THEN
!-+-+-+-+-+-+-+-+-+-+-+-+
!
!
            Write(99,*) 'Local Study of Type#3_ On a map, at a given time '
 If (iv==1) Write(* ,*) 'Local Study of Type#3_ On a map, at a given time '
!
!
!
    READ(1,*) IR, IT, IG, DIR, DIT, DIG
!
!
                call m_3(30600) ! One of the above switches must be == 1
!
!
          READ(1,*)  L1, L2, EPL
!
             call m_3(30650) ! Bounds on L1, L2, EPL
!

	  READ(1,*)  C1, C2, EPC
!
             call m_3(30660) ! Bounds on C1, C2, EPC
!
!
          READ(1,*)  GIVEN_TIME       
!
          READ(1,*)  I_THI   ! Reads the primary load thickness switch
!
             call m_3(30661) ! Checks if the primary load thickness can
                             ! be computed
!
              write(99,*) 'The following variables will be computed and printed:'
     if(iv==1)write(* ,*) 'The following variables will be computed and printed:'
! 
if( IR ==1) write(99,*) '-> Radial displacement'
if( it ==1) write(99,*) '-> Tangential displacement'
if( IG ==1) write(99,*) '-> Geoid Heigth'
if( dIR==1) write(99,*) '-> Radial velocity'
if( dit==1) write(99,*) '-> Tangential velocity'
if( dIG==1) write(99,*) '-> Rate of change of the geoid heigth'
IF(i_thi==1)WRITE(99,*) '-> Thickness of the primary load'
!
if(iv==1 .and. IR ==1) write(*,*)  '-> Radial displacement'
if(iv==1 .and. it ==1) write(*,*)  '-> Tangential displacement'
if(iv==1 .and. IG ==1) write(*,*)  '-> Geoid Heigth'
if(iv==1 .and. dIR==1) write(*,*)  '-> Radial velocity'
if(iv==1 .and. dit==1) write(*,*)  '-> Tangential velocity'
if(iv==1 .and. dIG==1) write(*,*)  '-> Rate of change of the geoid heigth'
IF(iv==1 .and.i_thi==1)WRITE(*,*) '-> Thickness of the primary load'
!
!
!
!
          K=0
!     
          npul = int ((l2-l1)/epl) + 1
          npuc = int ((c2-c1)/epc) + 1
!
          DO I = 1, npuc 
           DO J	 = 1, npul 
!
             CC = C1 + (I-1)*EPC
             LL = L1 + (J-1)*EPL 
!
             if( CC>C2 .or. LL>L2 ) then 
	      goto 405
	     else 
	      k=k+1   	     		             
              LONG_OBS(K)=LL
              COLA_OBS(K)= CC
!
              kk=k; call m_3(30625) ! Cheks the bounds on long. & colat.
!
             ENDIF  
!
           ENDDO
          ENDDO 
!
 405 NOBS = K
!
          call m_3(30670) ! Too many observers
!
!
            Write(99,*) 'There are', NOBS, ' points on the selected map ' 
 If (iv==1) Write(* ,*) 'There are', NOBS, ' points on the selected map ' 
!
!
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
     ENDIF  ! endif on Local_Study == 3
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!
!
!
!
!
!
!+.+.+.+.+.+.+.+.+.+.+.+.+.
   IF(type_inten==4) THEN
!+.+.+.+.+.+.+.+.+.+.+.+.+.
!
!
!
            Write(99,*) 'Local Study of Type#4'
 If (iv==1) Write(* ,*) 'Local Study of Type#4'
!
!
!
  READ(1,'(A30)') file_sparsi
!
!
  READ(1,*)  GIVEN_TIME
!
!
!
  IR = 0; IT = 0; IG = 0
!  
  DIR= 1; DIT= 1; DIG= 0
!
              write(99,*) 'The following variables will be computed and printed:'
     if(iv==1)write(* ,*) 'The following variables will be computed and printed:'
! 
if( IR ==1) write(99,*) '-> Radial displacement'
if( it ==1) write(99,*) '-> Tangential displacement'
if( IG ==1) write(99,*) '-> Geoid Heigth'
if( dIR==1) write(99,*) '-> Radial velocity'
if( dit==1) write(99,*) '-> Tangential velocity'
if( dIG==1) write(99,*) '-> Rate of change of the geoid heigth'
!
if(iv==1 .and. IR ==1) write(*,*)  '-> Radial displacement'
if(iv==1 .and. it ==1) write(*,*)  '-> Tangential displacement'
if(iv==1 .and. IG ==1) write(*,*)  '-> Geoid Heigth'
if(iv==1 .and. dIR==1) write(*,*)  '-> Radial velocity'
if(iv==1 .and. dit==1) write(*,*)  '-> Tangential velocity'
if(iv==1 .and. dIG==1) write(*,*)  '-> Rate of change of the geoid heigth'
!
!
!
!
    j=0
    OPEN(2,file=file_sparsi,ERR=120,status='unknown')
    j=1
!
120 IF(j==0) call m_3(30680)  ! Error on opening file_sparsi
!
!
  READ(2,*) ILB
!
           call m_3(30690)    ! The first record of file_sparsi must be 1 or 2
!
!
!

!
!'''''''''''''''''''
  If (ilb==2 )  Then  ! Form #2 <- The file contains
!                                  long_1, cola_1, long_2, cola_2
!'''''''''''''''''''
!
          write(99,*) 'The user-supplied file is of the Form #2'
 if(iv==1)write(* ,*) 'The user-supplied file is of the Form #2'
!
!
 J=0
!------------------
 DO KK=1, max_sites
!------------------
!
  READ(2,*,END=474) Long_site_1(KK),Cola_site_1(KK),Long_site_2(KK),Cola_site_2(KK)
!
 J=J+1
!
                    call m_3(30700) ! The two sites cannot coincide
                    call m_3(30710) ! Checking their longitude & colatitude
!
!
Write(99,*) 'Couple of Sites # ', kk
Write(99,'(a14,F10.4,a7,1x,F10.4,a7)')&
'Site#1: Long=',Long_site_1(kk),'Cola =',Cola_site_1(kk),' (deg)'
Write(99,'(a14,F10.4,a7,1x,F10.4,a7)')&
'Site#2: Long=',Long_site_2(kk),'Cola =',Cola_site_2(kk),' (deg)'
        IF(iv==1) THEN
Write(* ,*) 'Couple of Sites # ', kk
Write(* ,'(a14,F10.4,a7,1x,F10.4,a7)')&
'Site#1: Long=',Long_site_1(kk),'Cola =',Cola_site_1(kk),' (deg)'
Write(* ,'(a14,F10.4,a7,1x,F10.4,a7)')&
'Site#2: Long=',Long_site_2(kk),'Cola =',Cola_site_2(kk),' (deg)'
        ENDIF

!---------
     ENDDO
!---------
!
 CLOSE(2)
!
 474 NOBS = J
!
 call m_3(30730)   ! NOBS must be >=1
!
           WRITE(99,*) 'Read ', NOBS, 'couples of sites from file ', file_sparsi
 IF(iv==1) WRITE(* ,*) 'Read ', NOBS, 'couples of sites from file ', file_sparsi
!
!
!'''''''''''''''''''''''
 ENDIF ! endif on Form#2
!'''''''''''''''''''''''
!
!
!
!
!
!++++++++++++++++++
  If (ilb==1)  Then  ! Form #1 <- The file contains Sitename#1, Sitename#2
!++++++++++++++++++
!

!
     J=0
     DO KK=1, max_sites
        READ(2,'(a8,2x,a8)',END=475) Name_1(kk), Name_2(kk)
!
        junk8 = Name_1(kk)
!
        IF(junk8(1:8) == '        ') GOTO 475
!
        J=J+1
!
        CALL m_3(30720)  ! Name_1 must be /= Name_2
!
!
 	CALL Find_Station (Name_1(kk), Long_site_1(kk), Cola_site_1(kk))
 	CALL Find_Station (Name_2(kk), Long_site_2(kk), Cola_site_2(kk))
!
!
Write(99,*) 'Couple of Sites # ', kk
Write(99,'(a8,a7,F10.4,a7,F10.4,a7)')&
Name_1(kk), 'Long =',Long_site_1(kk),'Cola =',Cola_site_1(kk),' (deg)'
Write(99,'(a8,a7,F10.4,a7,F10.4,a7)')&
Name_2(kk), 'Long =',Long_site_2(kk),'Cola =',Cola_site_2(kk),' (deg)'
        IF(iv==1) THEN
Write(* ,*) 'Couple of Sites # ', kk
Write(* ,'(1x,a8,a7,F10.4,a7,F10.4,a7)')&
Name_1(kk), 'Long =',Long_site_1(kk),'Cola =',Cola_site_1(kk),' (deg)'
Write(* ,'(1x,a8,a7,F10.4,a7,F10.4,a7)')&
Name_2(kk), 'Long =',Long_site_2(kk),'Cola =',Cola_site_2(kk),' (deg)'
        ENDIF
!
!---------
     ENDDO
!---------
!
 475 NOBS=J
!
     call m_3(30730)   ! NOBS must be >=1
!
 CLOSE(2)
!
           WRITE(99,*) 'Read ', NOBS, 'couples of sites from file ', file_sparsi
 IF(iv==1) WRITE(* ,*) 'Read ', NOBS, 'couples of sites from file ', file_sparsi
!
!
!'''''''''''''''''''''''
 ENDIF ! endif on Form#1
!'''''''''''''''''''''''
!
!
!
!+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+
     ENDIF  ! endif on Local_Study == 4
!.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.
!
!
!
!

!
!
!===========================================
!==========================================*
   ENDIF    ! Endif on kw Local_Study  *
!==========================================*
!===========================================
!
!
!
!
!
!--------------------------------------------
   IF (KEYWORD(1:12)=='Global_Study') THEN
!--------------------------------------------
!
   n_este = n_este+1
!
          call m_3(30581) ! Ad_Hoc or Ice_* must be active
!
          call m_3(30583) ! Global_Study must appear once

          call m_3(30584) ! If Global_Study is active, Local_Study must be not
!
                WRITE(99,*) '> found KEYWORD  ', keyword
    IF(IV==1)   WRITE(* ,*) '> found KEYWORD  ', keyword
!
!
   IR = 0;  IT=0;  IG=1
  DIR = 0; DIT=0; DIG=1
!
!
    READ(1,*) type_esten
!
!
               call m_3(30591) ! Cheking bounds for type_esten
!
!
!  ++++++++++++++++++++++
   IF(type_esten==1) THEN
!  ++++++++++++++++++++++
!
            Write(99,*) 'Global Study of Type#1: Stokes coefficients vs time '
 If (iv==0) Write(99,*) 'Global Study of Type#1: Stokes coefficients vs time '
 If (iv==1) Write(* ,*) 'Global Study of Type#1: Stokes coefficients vs time '
!
!
!
   READ(1,*) LDE, MDE
!
!
               call m_3(30770) ! A check LDE & MDE
!       
!
   Read(1,*) Inorm      
!
!
               call m_3(30780) ! A check for Inorm
!
!
   READ(1,*) T1, T2, DTI
!
!
               call m_3(30790) ! A check for T1, T2 & DTI
!
!
!  +++++++++++++++++++++++++++++++++
   ENDIF ! endif on type_estensive=1
!  +++++++++++++++++++++++++++++++++
!
!
!
!
!  ......................
   IF(type_esten==2) THEN
!  ......................
!
!
            Write(99,*) 'Global Study of Type#2: Inertia tensor vs time '
 If (iv==1) Write(* ,*) 'Global Study of Type#2: Inertia tensor vs time '
!
!
   READ(1,*) T1, T2, DTI
!
!
            call m_3(30790) ! A check for T1, T2 & DTI
!
            call m_3(30800) ! A check for l_min & l_max in kw Harmonic_Degrees
!
!  .................................
   ENDIF ! endif on type_estensive=2
!  .................................
!
!
!===========================================
!                                          *
   ENDIF    ! Endif on kw Global_Study  *
!                                          *
!===========================================
!
!
10001 CONTINUE    ! Loop on file task_3.dat
!
10002 CONTINUE 
!
!
! 
          Write(99,*) 'task_3.dat is being closed... done '
if(iv==1) Write(* ,*) 'task_3.dat is being closed... done '
!
CLOSE (1)
!
!
!
!
!
!
!
!
!
!
!
! ooooooooooooooooooooooooooooooooooooooooooooooooooo
! o                                                 o
! o     The input data  have been assimilated       o
! o                                                 o
! o     _Now I procede with the computations_       o
! o                                                 o
! ooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!
!  
!
!
!
  junkd = 0.d0 
!
         OPEN(31,file='coeff.tmp',status='UNKNOWN',form='unformatted')
!
          Write(99,*) 'Expanding the loads in spherical harmonics ...'
if(iv==1) Write(* ,*) 'Expanding the loads in spherical harmonics ...'
!
!++++++++++++++++++++++++++++
	DO 55011 j=1, Mel_Eff
!++++++++++++++++++++++++++++
!
             If(C_l(j) <= 21) Then      ! Axis-simmetric load with c_l<=21
                    HEIGHT  =  Altezz(j) ! *******************************
                 AMPLITUDE  =  Amplit(J)		  
                 CALL AXIS_LOAD (0, C_l(j), AMPLITUDE, HEIGHT, MASS)  		  
                 Grande(j)  =  MASS         
               SSIGMA(J,:)  =  SIGMA(:)
!
                         Endif  
!
             If(C_l(j) == 30) Then      ! Axis-simmetric load with c_l==30
                       MASS = Grande(j)  ! *******************************
                 CALL AXIS_LOAD (0, C_l(j), junkd, junkd, MASS)    
               SSIGMA(J,:)  = SIGMA(:)
!
                         Endif  
!
             If(C_l(j) == 50) Then      ! Non-Axis-simmetric load with c_l==50
	              HEIGHT = Altezz(J) ! ***********************************
		          DL = D_long(j)
                          DT = D_cola(j)
		      CEN(1) = Long_c(j) 
		      CEN(2) = Cola_c(j) 
                 CALL RECT_LOAD0 (0, CEN(1), CEN(2), DL, DT, HEIGHT, MASS)       
                   Grande(j)   = MASS   
		   Write(31) FF, GG
!
                         Endif 
!
!+++++++++++++++
 55011 CONTINUE 
!+++++++++++++++
!
!
  close(31)
!
!
!
!
!
!                                             On a point, or a group of points
!+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+
         IF ( TYPE_INTEN ==1 .OR. TYPE_INTEN==2 )  THEN
!.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.
!
!
  IF(type_inten==1)          WRITE(99,*) 'Local analysis of type #1'
  IF(type_inten==1.and.iv==1)WRITE(* ,*) 'Local analysis of type #1... WAIT'
  IF(type_inten==2)          WRITE(99,*) 'Local analysis of type #1'
  IF(type_inten==2.and.iv==1)WRITE(* ,*) 'Local analysis of type #1... WAIT'
!
   IF(IR ==1 .or.  IT==1 .or.  IG==1) Open(10,File='disp.his',Status='Unknown')
   IF(DIR==1 .or. DIT==1 .or. DIG==1) Open(11,File='rate.his',Status='Unknown')
!
   call read_ocean
!
!
!           ++++++++++++++++
	    DO 222 K=1, NOBS
!           ++++++++++++++++
!
            call etime(txt,ela0)       ! It was: ela0=etime(txt)
!
               OBS(1) = LONG_OBS(K)
	       OBS(2) = COLA_OBS(K)
!
               call printtre (101)   ! puts an header on disp.his & rate.his
!
!              ----------------------------
               DO 111 TIME_BP = T1, T2, DTI
!              ----------------------------
!
               OPEN(31,file='coeff.tmp',status='UNKNOWN',form='unformatted')
!
!
                  call etime(txt,ela1)     ! It was: ela1=ETIME(txt)
!
                  mice = 0.D0
!
                  WECT(:)=0.D0
 	        D_WECT(:)=0.D0
!
!               ...................
                DO 888 J=1, Mel_Eff
!               ...................
!
                  call printtre (12)
!
		  call change_time
!
                       LOCAL_TIME = TFEA - TIME_BP
!
                   call Convol_Load_History  (L_h(j), Local_Time)
!
                       mice =  mice + grande(j)*   Load_History (L_h(j), Local_Time)
!
                       If(C_L(j) <= 30) Then
 		                CEN(1)  = Long_c(j)
		                CEN(2)  = Cola_c(j)
	                        SIGMA(:) = SSIGMA(J,:)
                       Call AXIS_DISP0  ( IOC(j), OBS, CEN, VECT, D_VECT)
                                        Endif

                       If(C_L(j) == 50) Then
                                READ(31) FF, GG
	  	       CAll RECT_DISP0  ( IOC(J), OBS, VECT, D_VECT)
                                              Endif
!
                   WECT(:) =   WECT(:) +   VECT(:)
 	         D_WECT(:) = D_WECT(:) + D_VECT(:)
!
!
!
!..............
  888 CONTINUE  ! Enddo on the ice elements
!..............
!
!
   IF(IR ==1 .or. IT==1  .or.  IG ==1) &
   Write(10,'(f9.4,1x,4(1x,f9.4),1x,E14.5)') &
                 TIME_BP,   wect(1),   wect(2),  &
                                   wect(3),   wect(4),    Mice
!
   IF(DIR==1 .or. DIT==1 .or. DIG ==1) &
   Write(11,'(f9.4,1x,4(1x,f9.4),1x,E14.5)') &
                 TIME_BP, d_wect(1), d_wect(2),  &
                                 d_wect(3), d_wect(4),    Mice
!
   Close(31)
!
!-------------
 111 CONTINUE   ! Enddo on time
!-------------

  call etime(txt,ela4)      ! It was: ela4=etime(txt)
!
            WRITE(99,*) 'Time elapsed for an observer (s) =', ela4-ela0
!
!
!+++++++++++++
 222 CONTINUE   ! Enddo on the observers
!+++++++++++++
!
   IF(type_inten==1)           WRITE(99,*) 'Done'
   IF(type_inten==1.and.iv==1) WRITE(* ,*) 'Done'
   IF(type_inten==2)           WRITE(99,*) 'Done'
   IF(type_inten==2.and.iv==1) WRITE(* ,*) 'Done'
!
            Close(10)
            Close(11)
!
 IF(iv==1 .and. (IR==1 .or. it==1  .or.  IG==1)) &
        WRITE(*,*) '-> Output written on file disp.his'

 IF(iv==1 .and. (dIR==1.or.dit==1.or.dIG==1)) &
        WRITE(*,*) '-> Output written on file rate.his'
!
 IF(IR==1 .or. it==1  .or.  IG==1) &
        WRITE(99,*) '-> Output written on file disp.his'
!
 IF(dIR==1.or.dit==1.or.dIG==1) &
        WRITE(99,*) '-> Output written on file rate.his'
!
!.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+
        ENDIF  ! <--  On Type_inten =1 or =2
!+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.+.
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
            IF ( TYPE_INTEN ==3 )  THEN   ! On a map ...
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
!
!
!
            WRITE(99,*) 'Local analysis of type #3 (on a map)'
 IF(iv==1)  WRITE(* ,*) 'Local analysis of type #3 (on a map) please WAIT ...'
!
!
If(IR ==1)   Open(48,file= 'A_rad.dat', status='unknown')
If(IT ==1)   Open(49,file= 'A_lon.dat', status='unknown')
If(IT ==1)   Open(50,file= 'A_lat.dat', status='unknown')
If(IG ==1)   Open(51,file= 'A_geo.dat', status='unknown')
!
If(DIR==1)   Open(52,file='Ad_rad.dat', status='unknown')
If(DIT==1)   Open(53,file='Ad_lon.dat', status='unknown')
If(DIT==1)   Open(54,file='Ad_lat.dat', status='unknown')
If(DIG==1)   Open(55,file='Ad_geo.dat', status='unknown')
!
IF(i_thi==1) OPEN(56,file='load_thick.dat',status='unknown')
!
!
    Call READ_OCEAN
!
!
!========================
     DO 2227 K=1, NOBS
!========================
!
        call etime(txt,ela0)    ! It was: ela0=etime(txt)
!
         OBS(1) = LONG_OBS(K)
         OBS(2) = COLA_OBS(K)
!
         OPEN(31,file='coeff.tmp',status='UNKNOWN',form='unformatted')
!
         MICE  = 0.D0
!
           WECT(:)=0.D0
         D_WECT(:)=0.D0
!
         PRES     =0.D0    ! Mass/surface at the observer position
	 PRES_EL  =0.d0    ! MAss/surface due to a single ice element
!
!
!    '''''''''''''''''''''''''
        DO 8887 J=1, Mel_Eff
!    '''''''''''''''''''''''''
!
                call change_time
!
                  call printtre (3)
!
                       LOCAL_TIME = TFEA - GIVEN_TIME
!
                   call Convol_Load_History  (L_h(j), Local_Time)
!
                  IF(k==1) &
                  mice = mice   + grande(j)*  Load_History (L_h(j), Local_Time)
!
                  ti_hi =  Load_History (L_h(j), Local_Time)
!
                       If(C_L(j) <= 30) Then
 		                CEN(1)  = Long_c(j)
		                CEN(2)  = Cola_c(j)
	                        SIGMA(:) = SSIGMA(J,:)
                       Call AXIS_DISP0  (IOC(j), OBS, CEN, VECT, D_VECT)

                       If(I_THI==1) &
	               CALL AXIS_PRESSURE (OBS, CEN, TI_HI, Pres_El)

                                        Endif
!
                       If(C_L(j) == 50) Then
                                READ(31) FF, GG
	  	       CAll RECT_DISP0  (IOC(j), OBS, VECT, D_VECT)

                       If(I_THI==1) &
                       Call RECT_PRESSURE (OBS, TI_HI, PRES_EL)
                                              Endif
!
                   PRES    =   PRES    +   PRES_EL
                   WECT(:) =   WECT(:) +   VECT(:)
 	         D_WECT(:) = D_WECT(:) + D_VECT(:)
!
!
!
!''''''''''''''''''
  8887   CONTINUE  ! enddo on the ice elements
!''''''''''''''''''
!
!
         call printtre (105)  ! writes on the A*.dat files
!
!            
 Close(31) 
!
!
 call etime(txt,ela4)     ! It was: Ela4=Etime(txt)
!
!
 Write(99 ,'(a26,i5,a2,i6,a1,f14.5,a1)') &
      'Time elapsed for observer', k,'of', nobs,'=',ela4-ela0,'s'
!
!
!==================
  2227    Continue       ! Enddo on the observers
!==================
!
!
If(IR==1)   CLOSE(48);  If(DIR==1) CLOSE(52)
If(IT==1)   CLOSE(49);  If(DIT==1) CLOSE(53)
If(IT==1)   CLOSE(50);  If(DIT==1) CLOSE(54)
If(IG==1)   CLOSE(51);  If(DIG==1) CLOSE(55)
IF(i_thi==1)CLOSE(56)
!
            call printtre(31)  ! Print list of the output files...
!
!+-+-+-+-+
  ENDIF  ! <--  On Type_inten ==3
!-+-+-+-+-
!
!
!
!
!
!
!
!
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
            IF ( TYPE_INTEN ==4 )  THEN   ! For baselines ...
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!
!
           WRITE(99,*) 'Local analysis of type #4 (baselines)'
IF(iv==1)  WRITE(* ,*) 'Local analysis of type #4 (baselines) please WAIT ...'
!
!
  Open (55,file='ltv_rates.dat',status='unknown')  ! Output file for L, T, & V rates
!
  call printtre (129)   ! Writes an header on ltv_rates.dat
!
  Call READ_OCEAN
!
!
!------------------------------
	DO 2127 K=1, NOBS       ! do-loop on the couples
!------------------------------
!
         OPEN(31,file='coeff.tmp',status='UNKNOWN',form='unformatted')
!
         call etime(txt,ela0)    ! It was: ela0=etime(txt)
!
         OBS_1 (1) = LONG_Site_1 (K); OBS_1 (2) = COLA_Site_1 (K)
         OBS_2 (1) = LONG_Site_2 (K); OBS_2 (2) = COLA_Site_2 (K)
!
!
         D_WECT_1(:)=0D0    ! Solution vector for Site#1 
	 D_WECT_2(:)=0D0    !    "        "    "  Site#2          	        
!
!+===============================+
        DO 8187 J=1, Mel_Eff
!+===============================+
!
                    call change_time
!
                    call printtre (4)
!
                       LOCAL_TIME = TFEA - GIVEN_TIME
!
                    call Convol_Load_History  (L_h(j), Local_Time)
!
           If(c_l(j) <= 30) Then
 		   CEN(1)      = Long_c(j) 
		   CEN(2)      = Cola_c(j)  
  	           SIGMA(:)    = SSIGMA(J,:)                            
                   Call AXIS_DISP0  (ioc(J), OBS_1, CEN, vJUNKd, D_VECT_1)
                   Call AXIS_DISP0  (ioc(J), OBS_2, CEN, vJUNKd, D_VECT_2)
           Endif 
!	  		
           If(c_l(j) == 50) Then
!
	        READ(31) FF, GG
!
	  	CAll RECT_DISP0    ( IOC(J), OBS_1,      vJUNKd, D_VECT_1)
	  	CAll RECT_DISP0    ( IOC(J), OBS_2,      vJUNKd, D_VECT_2)		
           Endif	  
!
	   D_WECT_1(:) = D_WECT_1(:) + D_VECT_1(:)  ! Updating site#1
	   D_WECT_2(:) = D_WECT_2(:) + D_VECT_2(:)  ! Updating site#2
!  
!
!  
!===================
  8187   CONTINUE     ! enddo on the ice elements
!===================
!
!
 Call BASE_RATES ( Long_Site_1(k), Cola_Site_1(k), &
                   D_WECT_1(1), D_WECT_1(2), D_WECT_1(3), &
                   Long_Site_2(k), Cola_Site_2(k), &
                   D_WECT_2(1), D_WECT_2(2), D_WECT_2(3), & 
                   L_RATE, T_RATE, V_RATE) 
!
!
       IF(ilb==1) THEN
WRITE(55,'(2(a9,1x,2(f6.2)),3(f6.2))') &
   name_1(k), Long_Site_1(k), Cola_Site_1(k), &
   name_2(k), Long_Site_2(k), Cola_Site_2(k), &
   L_RATE, T_RATE, V_RATE
   ELSEIF(ilb==2) THEN
WRITE(55,'(2(a9,1x,2(f6.2)),3(f6.2))') &
   'site #1', Long_Site_1(k), Cola_Site_1(k), &
   'site #2', Long_Site_2(k), Cola_Site_2(k), &
   L_RATE, T_RATE, V_RATE
ENDIF
!
!
 CLOSE(31)
!
!
 call etime(txt,ela4)      ! It was: Ela4=Etime(txt)
 WRITE (99,'(a30,i5,a3,i6,a1,f14.5,a1)') &
            'Elapsed time for the couple #', k,' of', nobs,'=',ela4-ela0,'s'
!
!+------------------ 
  2127    Continue    ! enddo on the couples of sites
!+------------------
!
 CLOSE(55)
!    
!
 IF(iv==1) WRITE(* ,*) '-> Output written on file ltv_rates.dat'
           WRITE(99,*) '-> Output written on file ltv_rates.dat'
!
!
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
  ENDIF  ! <--  On Type_inten ==4
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!
!
!
!
!
!
!
!
!
!---------------------------------------------------------------
            IF ( TYPE_ESTEN ==1 )  THEN   ! Stokes vs. time ...
!---------------------------------------------------------------
!
!
!
           WRITE(99,*) 'Global analysis of type #1 (Stokes)'
IF(iv==1)  WRITE(* ,*) 'Global analysis of type #1 (Stokes) ...'
!
!
  open(10,file='stokes.his',status='unknown')
  open(11,file='stoked.his',status='unknown') 
!
!
  call DATE_AND_TIME (date,timc)
  do jj=10, 11 
  Write(JJ, *) '# done by Task#3 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  enddo 	     
!
!
  Do IJUNK=10, 11
     Write(IJUNK,*) '# Degree   =', LDE 
     Write(IJUNK,*) '# Order    =', MDE 
  End Do 
!
!
  If(Inorm == 1) Then 
  Do IJUNK=10, 11
     Write(IJUNK,*) '# The coefficients are Fully Normalized'  
  End Do   
      Else 
  Do IJUNK=10, 11
     Write(IJUNK,*) '# The coefficients are NOT Fully Normalized'  
  End Do   
      Endif 
! 
!
  Write(11,*) '# x x x x x x x x x x x x x x x x x x x x x x x x'
  Write(11,*) '#      Divide dc_lm and ds_lm by 10**11          '
  Write(11,*) '#    to obtain values in units of yr**(-1)       '
  Write(11,*) '#  (time is in kyrs, the load mass is in kg)     '
  WRITE(11,*) '#                                                '
!
  Write(10,*) '# x x x x x x x x x x x x x x x x x x x x x x x x x x x '
  Write(10,*) '# Divide c_lm and s_lm by 10**6 to obtain actual values '  
  Write(10,*) '#      (time is in kyrs, the load mass is in kg)        '
  WRITE(10,*) '# '
!  
  Write(10,*) '# time BP      c_lm          s_lm         mass   '
  Write(11,*) '# time BP     dc_lm         ds_lm         mass   '
!
!
  Call READ_OCEAN
!
!                           
!
!       +==========================+
	Do 131 TIME_BP = T1, T2, DTI
!       +==========================+
!
        call etime(txt,ela0)     ! It was: Ela0 = Etime(txt)
!
        MICE  = 0.D0
!
        OPEN(31,file='coeff.tmp',status='unknown',form='unformatted')
!
          UUU (:) = 0.D0
        D_UUU (:) = 0.D0
!	     	        
!       +.................+
        DO 838 J=1, Mel_Eff
!       +.................+
!
               call change_time
!
               call printtre (5)
!
               LOCAL_TIME = TFEA - TIME_BP
!
                   call Convol_Load_History  (L_h(j), Local_Time)
!
                  mice =   mice + grande(j)*   Load_History (L_h(j), Local_Time)
!
          If(c_l(j) <= 30) Then
 		CEN(1)      = Long_c(j) 
		CEN(2)      = Cola_c(j)   
	        SIGMA(:)    = SSIGMA(J,:)
                Call AXIS_STOK (Lde, Mde, IOC(j), Inorm, CEN, UU, D_UU)
          Endif 
!
          If(c_l(j) == 50) Then
                READ(31) FF, GG
	  	CAll RECT_STOK (Lde, Mde, IOC(J), Inorm,      UU, D_UU)
          Endif
!
!
            UUU(:) =   UUU(:) +   UU(:)   ! Updates the solution
	  D_UUU(:) = D_UUU(:) + D_UU(:)   !    "     "     "
!
!+...............+
  838    CONTINUE   
!+...............+
!
!
   call etime(txt,ela4)      ! It was: Ela4 = Etime(txt)
!
   Fact1= 1D6
   Fact2= 1D11
!
   Write(10,'(f9.4,3(e14.5))')    time_BP,   FACT1*  uuu(1), &
                                            FACT1*  uuu(2), &
                                           MICE
   Write(11,'(f9.4,3(e14.5))')    time_BP,   FACT2*d_uuu(1)/1000D0, &
                                            FACT2*d_uuu(2)/1000D0, &
                                           MICE
!
   Close(31)
!
             WRITE(99,*) 'time elapsed for this time step (s) = ', ela4-ela0
   IF(iv==1) WRITE(99,*) 'time elapsed for this time step (s) = ', ela4-ela0


!   ===========================
131 CONTINUE      !time DO-loop
!   ===========================
!
 IF(IV==1) then
   WRITE(* ,*) '-> Output written on file stokes.his'
   WRITE(* ,*) '-> Output written on file stoked.his'
 ENDIF
!
 IF(IV==0) then
   WRITE(99,*) '-> Output written on file stokes.his'
   WRITE(99,*) '-> Output written on file stoked.his'
 ENDIF
!
!
!
!----------------------------------
 ENDIF   ! endif on type_esten ==1
!----------------------------------
!
!
!
!
!
!
!
!
!
!---------------------------------------------------------------
            IF ( TYPE_ESTEN ==2 )  THEN   ! Inertia vs. time ...
!---------------------------------------------------------------
!
            Write(99,*) 'Computing for Global Study of Type#2 (inertia)'
 If (iv==1) Write(* ,*) 'Computing for Global Study of Type#2 (inertia)'
!
             Open (36, File ='load_mass.dat ', status='unknown')
!
  open(10,file='iner.his',status='unknown')
  open(11,file='ined.his',status='unknown') 
!
  call DATE_AND_TIME (date,timc)
  do jj=10, 11 
  Write(JJ, *) '# done by Task#3 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  enddo 	     
!
!
 Write(36,*) '# done by Task#3 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 Write(36,*) '##################################################'
 Write(36,*) '# Mass of the primary load as a function of time #'
 Write(36,*) '##################################################'
 Write(36,*) '#time_BP (kyrs)  &  primary load mass (kg)        '   
!
!
 Write(10,*) '# = = = = = = = = = = = = = = = = = = = = = = = = ='                                                     
 Write(10,*) '# The inertia tensor is in a non-dimensional form. ' 
 Write(10,*) '# To obtain the actual values in units of kg*m**2  ' 
 Write(10,*) '# multiply the following components by MR**2, with '      
 Write(10,*) '# M = Mass of the Earth, and with R = Earth radius ' 
 Write(10,*) '#       [n.b.: M R**2 = 2.42 x 10**32 kg*m**2]     ' 
 Write(10,*) '#                                                  '
 Write(10,*) '#   time_BP (column #1) is given in units of kyrs  '
 Write(10,*) '#= = = = = = = = = = = = = = = = = = = = = = = = = '
! 
!
 Write(11,*) '#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=' 
 Write(11,*) '# To obtain the actual derivative of the inertia tensor  '
 Write(11,*) '# in kg*m**2/yr, you must multiply the components  below '       
 Write(11,*) '# by MR**2, with M= Mass of the Earth, and R= its Radius '
 Write(11,*) '#        [n.b.: M R**2 = 2.42 x 10**32 kg*m**2]          ' 
 Write(11,*) '#                                                        '
 Write(11,*) '#     time_BP (column #1) is given in units of kyrs      '
 Write(11,*) '#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=' 
!  
 Write(10,*) '# time_BP  (xz)         (yz)         (zz)         (xy)         (yy)         (xx)'
 Write(10,*) '# ' 						      	                     
 Write(11,*) '# time_BP  (xz)         (yz)         (zz)         (xy)         (yy)         (xx)'
 Write(11,*) '# ' 				                                            
! 
!
 call READ_OCEAN
!
!
!
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
	Do 141 TIME_BP = T1, T2, DTI
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
!
        call etime(txt,ela0)     ! It was: Ela0 = Etime(txt)
!
        MICE = 0.D0 
!
        OPEN(31,file='coeff.tmp',status='unknown',form='unformatted')
!
           INERZIA (:,:) = 0.D0
	 D_INERZIA (:,:) = 0.D0
!
!+==========================+
        DO 848 J=1, Mel_Eff
!+==========================+
!
               call printtre (5)
!
               call change_time
!
               LOCAL_TIME = TFEA - TIME_BP
!
                  call Convol_Load_History  (L_h(j), Local_Time)
!
                  mice =   mice + grande(j)*   Load_History (L_h(j), Local_Time)
!
          If(c_l(j) <= 30) Then
 		CEN(1)      =  Long_c(j)
		CEN(2)      =  Cola_c(j)
	        SIGMA(:)    =  SSIGMA(J,:)
                Call AXIS_INER (IOC(j), CEN, INE, D_INE)
          Endif 

          If(c_l(j) == 50) Then
	        Read(31) FF, GG
	  	CAll RECT_INER (IOC(J),      INE, D_INE)
          Endif
!
            INERZIA(:,:) =   INERZIA(:,:) +   INE(:,:)   ! Updating the inertia ...
	  D_INERZIA(:,:) = D_INERZIA(:,:) + D_INE(:,:)   ! Updating the inertia ...
!	  
!=================
  848    CONTINUE   
!=================
!
  fact1=1000.D0
!   
!
        Write(10,'(f9.4,2x,7(e12.5,1x))') time_bp, &
	    inerzia(1,3), inerzia(2,3), inerzia(3,3), &
	    inerzia(1,2), inerzia(2,2), inerzia(1,1)
!
        Write(11,'(f9.4,2x,7(e12.5,1x))') time_bp,  &
            d_inerzia(1,3)/fact1, d_inerzia(2,3)/fact1, d_inerzia(3,3)/fact1, &
	    d_inerzia(1,2)/fact1, d_inerzia(2,2)/fact1, d_inerzia(1,1)/fact1
!
        Write(36,'(2x,f9.4,1x,E20.8)') time_bp, mice
!
	Close(31)
!
!
        WRITE(99,*) 'Elapsed time for this time step (s)', Ela4-Ela0
!
!
!
!.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
  141    CONTINUE    ! time do-loop
!-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
!
         Close(10)
         Close(11)
         CLOSE(36)
!
 IF(IV==1) then
   WRITE(* ,*) '-> Output written on file iner.his'
   WRITE(* ,*) '-> Output written on file ined.his'
   WRITE(*, *) '-> Output written on file load_mass.dat'
 ENDIF
!
 IF(IV==0) then
   WRITE(99,*) '-> Output written on file iner.his'
   WRITE(99,*) '-> Output written on file ined.his'
   WRITE(99,*) '-> Output written on file load_mass.dat'
 ENDIF
!
!
!
!----------------------------------
 ENDIF   ! endif on type_esten ==2
!----------------------------------
!

!
         call m_3(30810)  ! A warning if Local_ and Global_Study
                          ! are both inactive
!
!
!
!
  end subroutine TASK_3 
!
!
!
!
!
!
!

!
 subroutine change_time
! ------------------------------------------------------------------------
! Finds the time of occurrence af a specific feature of each time-history
! ------------------------------------------------------------------------
  USE for_task3
  implicit NONE
!
!
         IF (L_h(j)==0) THEN  
	           TSL = T_sl(j) 
		   TFEA = TSL
                   Endif 
!
         IF (L_h(j)==1) THEN
                   TSU = T_su(j)
                   TFEA = TSU
                   Endif
!
         IF (L_h(j)==2) THEN
	           TSU = T_su(j)
                   TAU_LU = T_au(j)
	           TFEA = TSU
                   Endif 
!
         IF (L_h(j)==3) THEN
                   TSBD = T_sbd(j)
                   TAU_SD = T_au(j)
                   TFEA = TSBD
                   Endif
!
         IF (L_h(j)==4) THEN
                   TSBLD   = T_sbld(j)
                   NR      = N_R(j)
                   TAUC  = T_auc(j)
                   DINC  = D_inc(j)
                   TFEA  = TSBLD
                   Endif
!
         IF (L_h(j)==5) THEN
                   Period  = T_per(j)
                   tsz  = T_stz(j)
                   TFEA  = TSZ
                   Endif
!
         IF (L_h(j)==6) THEN
                   NSS = N_ss(j)
                   TNN = T_nn(j)
                   do I=0,NSS
                       a_hh(I) = Alh(j,I)
                       t_hh(I) = Tmp(j,I)
                   end do
                   do i=1, nss
                       ri(i) = (a_hh(i) - a_hh(i-1))/ &
                       (t_hh(i) - t_hh(i-1))
                   end do
!
                   alf(0)   = (a_hh(1)-a_hh(0)) - ri(1)*t_hh(1)
                   bet(0)   = ri(1)
                   alf(nss) = ri(nss)*t_hh(nss)
                   bet(nss) = -ri(nss)
!
                   do I=1, nss-1
                   alf(i) = (a_hh(i+1)-a_hh(i)) + &
                          ri(i)*t_hh(i) - ri(i+1)*t_hh(i+1)
                   bet(i) = ri(i+1)-ri(i)
                   end do
!
                   TSBPLP = T_sbplp(j)
!
                   TFEA =  TSBPLP
                   END IF
!
          IF (L_h(j)==7) THEN
                    NS = N_s(j)
                    TN = T_N(j)
                    DILTA = D_ilt(j)
                    DO I=0, NS+1
                        A_H(I) = ALT(J,I)
                    ENDDO
                    TSBPCP = t_sbpcp(j)
!
                    TFEA =  TSBPCP
                    Endif
!
         IF (L_h(j)==8) THEN
                   NS = N_s(j)
                   TN = T_N(j)
!
                   DILTA = D_ilt(j)
!
                   DO I=0, NS+1
                       A_H(I) = ALT(J,I)
                   ENDDO
!
                   TAUF = TAU_8
                   TSBPCP = t_sbpcp(j)
                   TFEA =  TSBPCP
         Endif
!
 end subroutine change_time 
!
!
!
!
!
! +=====================+
   function gprod (l,m) ! 
! +=====================+   
! 
   use common; implicit none 
!
! # This routines computes (l-m)!/(l+m)! in double precision. 
!
   real(dp) gprod
!
   integer(i4b) :: l, m, k 
!
   If(m<0.or.l<0) Then
     Write(99,*) 'ERROR in gprod: negative arguments'
     Write(99,*) '**** JOB ABORTED *****************';Stop
   Endif
!
   If(m>l) Then 
     Write(99,*) 'ERROR in gprod: m must be <= l'
     Write(99,*) '**** JOB ABORTED *************';Stop
   endif
!
!
   If((m==0) .OR. (m==0.and.l==0)) Then
    gprod  =  1.d0
    Return 
   Endif 
!
   If(m<=l) then        
     gprod = 1._dp/dfloat(l+m)    
      do k = 1, 2*m -1 
        gprod =  gprod/dfloat(l+m-k) 
      enddo  
   Endif 
!
! +-------------------+
   end function gprod !
! +-------------------+
!
!
!
!
!
!
!
!
 SUBROUTINE Find_Station (Name, Long, Cola) 
!  
! +---------------------------------------------------------------------+
!  Finds the station with name 'Name' in the GSFC  site locations file, |
!  and returns the Longitude (Long) and Colatitude (Cola) of the site.  |
! [Done in Urbino on Aprile, 16, 2002 ]. 		                |
! +---------------------------------------------------------------------+
!    
  Use COMMON; Implicit NONE
!
!
!     Last change:  GS   16 Apr 2002    8:38 am
!
!
INTEGER, PARAMETER :: numstat=149   ! Number of sites on the GSFC file 
!
REAL*8 LONG, COLA        ! Long. and Colat. of the Station 'Name' (output) 
!
REAL*4 lats, lons        ! Lat. and long. in seconds
REAL*4 lat_deg, lat_min  ! Lat. in degrees & minutes
REAL*4 lon_deg, lon_min  ! Lon. in degrees & minutes 
!
INTEGER i                ! Do-loop index
!
INTEGER found
!
! # Auxiliary character variables for reading the GSFC file 
!
CHARACTER*1  j1        
CHARACTER*2  j2
CHARACTER*3  j3
CHARACTER*14 j14
!
CHARACTER*3 latd; CHARACTER*2 latm   
CHARACTER*3 lond; CHARACTER*2 lonm
!
!
CHARACTER*8 NAME           ! Input Station name 
CHARACTER*8 NAME_STAT      ! Station name as it appears in the GSFC file 
!
!
!
!
!
OPEN  (66,FILE='site_locations.2001cn.1',STATUS='unknown')
!
do i=1, 9
 READ (66, '(a1)') j1
enddo
!
!
!
found=0
!+-------------------+ 
   DO I=1, numstat   !     
!+-------------------+ 
!
!
READ   (66,'(1x,a8,a14,a2,a3,a2,a2,a2,f6.3, a2,a4,a2,a2,a2,f7.3)')  &
                         NAME_STAT, j14,                &
                         j2, latd, j2, latm, j2, lats,  &   
                         j3, lond, j2, lonm, j2, lons
!
!
If (NAME == NAME_STAT)       Then    ! <--- 
found=1
!
 OPEN(8,FILE='tmp.dat',STATUS='unknown')
!
write (8,'(a3)') latd
 IF(latm(1:1) == '0') WRITE(8,'(a1)')  latm(2:2)
 IF(latm(1:1) /= '0') WRITE(8,'(a2)')  latm
write (8,'(a3)') lond
 IF(lonm(1:1) == '0') WRITE(8,'(a1)')  lonm(2:2)
 IF(lonm(1:1) /= '0') WRITE(8,'(a2)')  lonm
!
 CLOSE(8)
!
!
 OPEN(8,FILE='tmp.dat',STATUS='unknown')
!
 READ (8,*)      lat_deg      ! lat. in Degrees
 READ (8,*)      lat_min      ! lat. in Minutes
 READ (8,*)      lon_deg      ! lon. in Degrees
 READ (8,*)      lon_min      ! lon. in Minutes
!
 CLOSE(8)
!
   LONG = lon_deg + lon_min/60. + lons/3600.
!
   If ( lat_deg <  0.0) COLA  = 90. -(lat_deg - lat_min/60. - lats/3600.)
   If ( lat_deg >= 0.0) COLA  = 90. -(lat_deg + lat_min/60. + lats/3600.) 
!
   ENDIF
!
!+----------+
  ENDDO     ! <--- enddo on the NASA file
!+----------+
!
  CLOSE(66)
!
IF(found==0) THEN
!
  IF(iv==0) THEN
    WRITE(99,*) 'ERROR from sbr TASK_3: Site name not found in the'
    WRITE(99,*) 'GSFC file site_locations.2001cn.1. JOB ABORTED ==';STOP
  ENDIF
  IF(iv==1) THEN
    WRITE(* ,*) 'ERROR from sbr TASK_3: Site name not found in the'
    WRITE(* ,*) 'GSFC file site_locations.2001cn.1. JOB ABORTED =='
    WRITE(99,*) 'ERROR from sbr TASK_3: Site name not found in the'
    WRITE(99,*) 'GSFC file site_locations.2001cn.1. JOB ABORTED ==';STOP
  ENDIF
!
ENDIF
!
 END SUBROUTINE FIND_Station 
!
!
!
!
!
!
!
!
!
  SUBROUTINE BASE_RATES ( Lon_1, Col_1, V_radi_1, V_Cola_1, V_long_1, & 
                          Lon_2, Col_2, V_radi_2, V_Cola_2, V_long_2, & 			                                      
			  L_RATE, T_RATE, V_RATE) 
!
!
! # Given the longitudes and colatitudes of two sites, and the
!   the velocity at each site in spherical coordinates, this 
!   routine computes the L, T, and V components of the velocity 
!   of site #2 w.r.t site #1. [done by GS in urbino on 18.02.2002] 
!
!
  Use COMMON; Use TRIGFCN
  Implicit NONE   
!
  Real(DP), PARAMETER :: EPS=1d-4    ! A small angle (degrees) 
!
  Real(DP) :: Lon_1     ! Longitude  of Site #1 (deg) 
  Real(DP) :: Col_1     ! Colatitude of Site #1 ( " )
!
  Real(DP) :: V_Radi_1  ! Velocity along the radius for Site #1  
  Real(DP) :: V_Cola_1  ! Velocity along Colatitude for Site #1 
  Real(DP) :: V_Long_1  ! Velocity along Longitude  for Site #1 
!  
  Real(DP) :: Lon_2     ! Longitude  of Site #2 (deg) 
  Real(DP) :: Col_2     ! Colatitude of Site #2 ( " )    
!
  Real(DP) :: V_Radi_2  ! Velocity along the radius for Site #2  
  Real(DP) :: V_Cola_2  ! Velocity along Colatitude for Site #2 
  Real(DP) :: V_Long_2  ! Velocity along Longitude  for Site #2    			             
!
  Real(DP) :: X(2), Y(2), Z(2) ! Cartesian Coordinates of the two Sites 
!
  Real(dp) :: lx, ly, lz  ! Cartesian components of versor 'L'
  Real(dp) :: tx, ty, tz  ! Cartesian components of versor 'T'  
  Real(dp) :: vx, vy, vz  ! Cartesian components of versor 'V'
!
  Real(dp) :: velo_1_x, velo_1_y, velo_1_z ! Cartesian components of the  
!					     absolute velocity at Site# 1
!
  Real(dp) :: velo_2_x, velo_2_y, velo_2_z ! Cartesian components of the  
!					     absolute velocity at Site# 2
!
!
  Real(dp) :: L_RATE  ! Velocity of #2 w.r.t. #1, L component      
  Real(dp) :: T_RATE  ! Velocity of #2 w.r.t. #1, T component     
  Real(dp) :: V_RATE  ! Velocity of #2 w.r.t. #1, V component        
!
  Real(dp) :: CC   ! Auxliary variable 
!
!  
 If(Lon_1<0.D0 .or. Lon_1>360d0) Then 
 Write(99,*) ' ERROR from Sbr. BASE_RATES: Longitude of  Site#1 out of Range '
 Write(99,*) ' #### JOB ABORTED ############################################ '; Stop 
 Endif 
!
 If(Col_1<0.D0 .or. Col_1>180d0) Then 
 Write(99,*) ' ERROR from Sbr. BASE_RATES: Colatitude of Site#1 out of Range '
 Write(99,*) ' #### JOB ABORTED ############################################ '; Stop 
 Endif   
! 
 If(Lon_2<0.D0 .or. Lon_2>360d0) Then 
 Write(99,*) ' ERROR from Sbr. BASE_RATES: Longitude of  Site#2 out of Range '
 Write(99,*) ' #### JOB ABORTED ############################################ '; Stop 
 Endif 
!
 If(Col_2<0.D0 .or. Col_2>180d0) Then 
 Write(99,*) ' ERROR from Sbr. BASE_RATES: Colatitude of Site#2 out of Range '
 Write(99,*) ' #### JOB ABORTED ############################################ '; Stop 
 Endif  
!
!
! # Cartesian coordinates of the two sites 
!  
        x(1) = Sind(Col_1) * Cosd(Lon_1)  
        y(1) = Sind(Col_1) * Sind(Lon_1)
	z(1) = Cosd(Col_1) 
!
        x(2) = Sind(Col_2) * Cosd(Lon_2)  
        y(2) = Sind(Col_2) * Sind(Lon_2)
	z(2) = Cosd(Col_2) 
!
!
   if ((raggio*sqrt((x(2)-x(1))**2 + (y(2)-y(1))**2 + (z(2)-z(1))**2) <= 1000.d0) .and. IV==0) THEN  
   Write(99,*) ' ERROR from Sbr. BASE_RATES:  The baseline must exceed 1000.0 m in lenght ' 
   Write(99,*) ' **** JOB ABORTED ******************************************************* '; Stop 
   Endif     	
   if ((raggio*sqrt((x(2)-x(1))**2 + (y(2)-y(1))**2 + (z(2)-z(1))**2) <= 1000.d0) .and. IV==1) THEN  
   Write(99,*) ' ERROR from Sbr. BASE_RATES:  The baseline must exceed 1000.0 m in lenght ' 
   Write(99,*) ' **** JOB ABORTED ******************************************************* '
   Write(* ,*) ' ERROR from Sbr. BASE_RATES:  The baseline must exceed 1000.0 m in lenght ' 
   Write(* ,*) ' **** JOB ABORTED ******************************************************* '; Stop  
   Endif  
!
!
! # Cartesian components of the unit vectors L, T, & V 
!
! # L 
!	
	lx = (x(2)-x(1)) 
	ly = (y(2)-y(1)) 
	lz = (z(2)-z(1)) 
!
        cc = Sqrt( lx*lx + ly*ly + lz*lz )
!
        lx = lx/cc 
	ly = ly/cc
	lz = lz/cc 
!
! # T         				
!
        tx = z(1)*y(2)-y(1)*z(2)
	ty = x(1)*z(2)-z(1)*x(2)
	tz = x(2)*y(1)-y(2)*x(1)
!	 
        cc = Sqrt( tx*tx + ty*ty + tz*tz )               		
!
        tx = tx/cc
	ty = ty/cc
	tz = tz/cc
!
! # V 
!
        vx = ly * tz - lz * ty 
	vy = lz * tx - lx * tz 
	vz = lx * ty - ly * tx  
!
!
! # Cartesian components of the velocity of site #1  
!
!
        velo_1_x  =  v_cola_1*cosd(Lon_1)*cosd(Col_1) - & 
	             v_long_1*sind(Lon_1)             + & 
		     v_radi_1*cosd(Lon_1)*sind(Col_1) 	
!		    		
        velo_1_y  =  v_cola_1*sind(Lon_1)*cosd(Col_1) + & 
	             v_long_1*cosd(Lon_1)             + & 
	   	     v_radi_1*sind(Lon_1)*sind(Col_1) 	 		
!
        velo_1_z  = -v_cola_1*sind(Col_1) + & 
	             v_radi_1*cosd(Col_1) 			
!
!
! # Cartesian components of the velocity of site #2  
!
!
        velo_2_x  =  v_cola_2*cosd(Lon_2)*cosd(Col_2) - & 
	             v_long_2*sind(Lon_2)             + & 
		     v_radi_2*cosd(Lon_2)*sind(Col_2) 	
!		    		
        velo_2_y  =  v_cola_2*sind(Lon_2)*cosd(Col_2) + & 
	             v_long_2*cosd(Lon_2)             + & 
		     v_radi_2*sind(Lon_2)*sind(Col_2) 	 		
!
        velo_2_z  = -v_cola_2*sind(Col_2) + & 
	             v_radi_2*cosd(Col_2) 			
!
! # Velocity of site#2 w.r.t. site#1, L, T, & V components 
!
        L_RATE = ( velo_2_x*lx + velo_2_y*ly + velo_2_z*lz ) - & 
	         ( velo_1_x*lx + velo_1_y*ly + velo_1_z*lz )
!
        T_RATE = ( velo_2_x*tx + velo_2_y*ty + velo_2_z*tz ) - & 
	         ( velo_1_x*tx + velo_1_y*ty + velo_1_z*tz )
!
        V_RATE = ( velo_2_x*vx + velo_2_y*vy + velo_2_z*vz ) - & 
	         ( velo_1_x*vx + velo_1_y*vy + velo_1_z*vz )	
!			
!  
  End Subroutine Base_Rates 
!
!
!
!
!
! ----------------------------+
  SUBROUTINE TASK_2           ! 
! ----------------------------+
!
! o--------------------------------------------------o
! | Computes the displacements and other quantities  |
! |    in the presence of a SINGLE glacial load      |
! o--------------------------------------------------o
!
!
 USE for_task2
 implicit NONE
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !  Beginning with the executable instructions  !
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!
! # The following integers are used to monitor the 
!     keyword sequence .... 
!
!
                 iarmo=0
              ideja=0
           inone=0
            iload=0
               ihist=0
             n_load=0
          n_hist=0
         n_inte=0
      n_este=0
!
!
!
!
!
Write(99,*) 'Opening file task_2.dat ...  '
!
!
      open(1,file='task_2.dat',status='unknown')
!
!
Write(99,*) 'Reading file task_2.dat ...  '
! 
!
!
 DO 10001 ITASK = 1, max_num_righe      
!
!
!
!
!
   READ(1,'(a30)',END=10002) KEYWORD
!
!
!
!
!
!...............................................
   IF(KEYWORD(1:16) == 'Harmonic_Degrees') THEN
!...............................................
!
   iarmo=1
!
   READ(1,*)  LMIN, LMAX
!
   READ(1,*)  IV
!
       call m_2(20000)  ! Check in IV=1 or 0
!
        call DATE_AND_TIME (date,timc)
!
          IF    (IV==0) THEN
        Write(99,*) '# ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
          ELSEIF(IV==1) then
        Write(*,*) '# ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
        ENDIF
!
    IF(iv==1) then
     WRITE(*,*) '= = = = = = = = = = = = = = ='
     Write(*,*) '   Reading file task_2.dat   '
     WRITE(*,*) '= = = = = = = = = = = = = = ='
    ENDIF
!
                WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
              Write(99,*) 'Lmin and Lmax ', Lmin, Lmax
    IF(IV==1) Write (*,*) 'Lmin and Lmax ', Lmin, Lmax
!
!
!
    i_loading=1
!
              call m_2(20010)   ! Check on l_min, l_max
!
!
    READ(1,*)  only_elastic
!
!
    IF(only_elastic==0.and.iv==1) &
           WRITE(*,*)  'Elastic AND viscous fields'
    IF(only_elastic==1.and.iv==1) &
           WRITE(*,*)  'ONLY Elastic fields'
    IF(only_elastic==0) WRITE(99,*) 'Elastic AND viscous fields'
    IF(only_elastic==1) WRITE(99,*) 'ONLY Elastic fields'
!
!
    READ(1,*)  DROP_MODES 
!
!
    IF(DROP_MODES==1.and.iv==1) &
           WRITE(*,*)  'The ''non physical'' modes will be killed'
    IF(DROP_MODES==0.and.iv==1) &
           WRITE(*,*)  'The ''non physical'' modes will NOT be killed' 
!
 IF(DROP_MODES==1) WRITE(99,*)'The ''non physical'' modes will be killed'
 IF(DROP_MODES==0) WRITE(99,*)'The ''non physical'' modes will NOT be killed'
!
!
!......................................
  ENDIF      ! on KW Harmonic_Degrees
!......................................
!
!
!
!
!
!
!
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
   IF(KEYWORD(1:10) == 'Make_Model') THEN     
!wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
!
   ideja=1
   inone=1
!
!
   call m_2(20020) ! status of Harmonic_Degrees
!
!
                WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
!
             Write(99,*)    'Building the model '
 IF(IV==1)   Write(*, *)    'Building the model '
!
   READ(1,*) NV 
!  
   READ(1,*) CODE
!
   READ(1,*) LT 
!
   READ(1,*) ILM 
!      
!
   Write(99,*) 'Number of VE layers      ', NV  
   Write(99,*) 'Model CODE is            ', CODE 
!
IF(IV==1) then
   Write(*,*) 'Number of VE layers      ', NV
   Write(*,*) 'Model CODE is            ', CODE
ENDIF
!
!
    Write(99,'(a42,1x,f12.4,1x,a3)') &
    'From input, the lithospheric thickness is', LT,  ' km'
    Write(99,'(a47,1x,i3)')       &
    'The ILM parameter (for NV=7 and 9) is set to  ', ILM
!
   IF(IV==1) THEN
    Write(*, '(a42,1x,f12.4,1x,a3)') &
    'From input, the lithospheric thickness is', LT,  ' km'
    Write(*, '(a47,1x,i3)')       &
    'The ILM parameter (for NV=7 and 9) is set to  ', ILM
   ENDIF
!
!
               Write(99,*) 'Mantle viscosity from BOTTOM to TOP (/1E21) '
     IF(IV==1) WRITE(*,*)  'Mantle viscosity from BOTTOM to TOP (/1E21) '
   Do k=1, nv
      Read(1,*) vis(k)
	          Write(99, '(a19,i4,1x,a1F10.4)') &
                       'Viscosity of layer', k, '=', vis (k)
   IF(IV==1) Write(*,  '(a19,i4,1x,a1,F10.4)') &
                       'Viscosity of layer', k, '=', vis (k)
               vis(k) = vis(k) * 1.q21    	  
   Enddo 
! 
!
   CALL SPEC (Code, Ilm, LT)

!
! +-----------------------              			      
       call SPECTRUM (1)
! +-----------------------
!
!+++++++++++++++++++++++++++
  If(only_elastic==0) then 
!+++++++++++++++++++++++++++  
!
             Write(99,*)'Writing the spectrum on file spectrum.dat'
  IF(IV==1)  Write(*,* )'Writing the spectrum on file spectrum.dat'
!
!
  open(3,file='spectrum.dat',status='unknown')
!
  Write(3,*) '# File spectrum.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  
  Write(3,*) '  # 1st column: l = Harmonic degree                '
  Write(3,*) '  # 2nd   "   : LOG10 (l)                          '
  Write(3,*) '  # 3th   "   : s, (kyrs**(-1))                    '
  Write(3,*) '  # 4rd   "   : LOG10(-s), with s in kyrs**(-1)    '
  Write(3,*) '  # 5th   "   : Relaxation time = -1000./s, (yrs)  '
  Write(3,*) '  # 6th   "   : LOG10 (Relaxation time (yrs))      '
!
DO l=lmin, lmax 
! 
    write (3, '(/)')
      do k = 1, nroots
       write (3, '(i4,1x,5(e15.7,1x))') l,                  & 
                                        log10 (dfloat (l)), &
                                        s(l,k),            &
					log10(-s(l,k)),    &
				        (-1000.q0/s(l,k)), &
				        log10(-1000.q0/s(l,k))
      enddo
!
ENDDO; close(3)
!
!++++++++++++
    endif
!++++++++++++
!
          open(13,file='h.dat',status='unknown')       
          open(14,file='l.dat',status='unknown')
          open(15,file='k.dat',status='unknown')
!
!
  call DATE_AND_TIME (date,timc)
  Write(13,*) '# File h.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  call DATE_AND_TIME (date,timc)
  Write(14,*) '# File l.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  	      	  
  call DATE_AND_TIME (date,timc)
  Write(15,*) '# File k.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  	  
!
!
!
          if(only_elastic==0) then 
          WRITE(99,*) 'Writing the elastic, fluid, and v-elastic h ldc on h.dat '
          WRITE(99,*) 'Writing the elastic, fluid, and v-elastic l ldc on l.dat '
          WRITE(99,*) 'Writing the elastic, fluid, and v-elastic k ldc on k.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic, fluid, and v-elastic h ldc on h.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic, fluid, and v-elastic l ldc on l.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic, fluid, and v-elastic k ldc on k.dat '
          endif
!
          if(only_elastic==1) then 
          WRITE(99,*) 'Writing the elastic h ldc on h.dat '
          WRITE(99,*) 'Writing the elastic l ldc on l.dat '
          WRITE(99,*) 'Writing the elastic k ldc on k.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic h ldc on h.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic l ldc on l.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic k ldc on k.dat '
          endif
!	  
!
          DO l=lmin, lmax 
!
          if(only_elastic==0) then 
                    write (13, '((i3,1x,96(1x,e20.8)))') &
  	     l, h_e(l), h_f(l), (h_v (l,k), k = 1, nroots)
                    write (14, '((i3,1x,96(1x,e20.8)))') &
	     l, l_e(l), l_f(l), (l_v (l,k), k = 1, nroots)
                    write (15, '((i3,1x,96(1x,e20.8)))') &
	     l, k_e(l), k_f(l), (k_v (l,k), k = 1, nroots)
	  endif
          if(only_elastic==1) then 
                    write (13, '((i3,1x,96(1x,e20.8)))') &
  	     l, h_e(l) 
                    write (14, '((i3,1x,96(1x,e20.8)))') &
	     l, l_e(l) 
                    write (15, '((i3,1x,96(1x,e20.8)))') &
	     l, k_e(l) 
	  endif  
!
!++++++++++++++++++++++++++++++
    if(only_elastic==0) then 
!++++++++++++++++++++++++++++++ 
!
          if(drop_modes==0) then 
	                   do k=1, nroots
	                 r_h (l,k) = h_v(l,k)/s(l,k)      ! Normalized residue for h
                         r_l (l,k) = l_v(l,k)/s(l,k)      ! Normalized residue for l
                         r_k (l,k) = k_v(l,k)/s(l,k)      ! Normalized residue for k
                         tek (l,k) = -1.0/s(l,k)     
	                   enddo 
	  endif 
!
          if(drop_modes==1)  then
!	   
             do k=1, nroots
              r_h (l,k) = 0.0
              r_l (l,k) = 0.0
              r_k (l,k) = 0.0 
!
              IF(vec(l,k)==1._sp) then
!
                r_h (l,k) = h_v(l,k)/s(l,k)      ! Normalized residue for h
                r_l (l,k) = l_v(l,k)/s(l,k)      ! Normalized residue for l
                r_k (l,k) = k_v(l,k)/s(l,k)      ! Normalized residue for k
                tek (l,k) = -1.0/s(l,k)          ! Relaxation time in k-yrs
!
              endif
!
             enddo
!
             endif 
!
!+++++++++++++++++++ 
       ENDIF
!+++++++++++++++++++
!
          ENDDO
!
          close(13)
          CLOSE(14)
          CLOSE(15)
!
!
          WRITE(99,*) 'Done with the computation of the Spectral quantities '
IF(IV==1) WRITE(* ,*) 'Done with the computation of the Spectral quantities '
!
          WRITE(99,*) 'Looking now for KWs Load_Geometry and Load_History'
IF(IV==1) WRITE(* ,*) 'Looking now for KWs Load_Geometry and Load_History'
!
!
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
  ENDIF  !  Endif on KW Make_Model
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!
!
!
!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
   IF(KEYWORD(1:14) == 'External_Model') THEN
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
   inone=1
!
           call m_2(20030) ! Checks the status of Make_Model
!
!
               WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
   IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
!
!
           call m_2(20040) ! Model description & a warning
!
!
    rhoea = ext_rhoea
!	
	NV=1    
    NROOTS=3 
!    
 Open(67, file='external_spectrum.dat',status='unknown')   
   do J=0, llmax 
   do k=1, nroots 
        read(67,'(i4,1x,2(e15.7,1x))',end=1901)    l,   s(l,k) ,  s(l,k)  
   enddo
   enddo
 1901 Close(67)
! stop 33411
!  
  Open(67, file='external_h.dat',status='unknown') 
          DO J=0, llmax 
             read (67, '((i3,1x,96(1x,e20.8)))', end=1902) &
	     l, h_e(l), h_f(l), (h_v (l,k), k = 1, nroots) 
	  ENDDO            
 1902 Close(67)   
  Open(68, file='external_l.dat',status='unknown')  
          DO J=0, llmax 
             read (68, '((i3,1x,96(1x,e20.8)))', end=1903) &
	     l, l_e(l), l_f(l), (l_v (l,k), k = 1, nroots)    
          ENDDO                             
 1903 close(68)                                     
      vec=1.
!
!      
!+++++++++++++++++++++++++++
  If(only_elastic==0) then 
!+++++++++++++++++++++++++++  
!
            Write(99,*)'Writing the spectrum on file spectrum.dat'
 IF(IV==1)  Write(*,* )'Writing the spectrum on file spectrum.dat'
!
!
  open(3,file='spectrum.dat',status='unknown')
!
  call DATE_AND_TIME (date,timc)
!
  Write(3,*) '# File spectrum.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  
!
  Write(3,*) '  # 1st column: l = Harmonic degree                '
  Write(3,*) '  # 2nd   "   : LOG10 (l)                          '
  Write(3,*) '  # 3th   "   : s, (kyrs**(-1))                    '
  Write(3,*) '  # 4rd   "   : LOG10(-s), with s in kyrs**(-1)    '
  Write(3,*) '  # 5th   "   : Relaxation time = -1000./s, (yrs)  '
  Write(3,*) '  # 6th   "   : LOG10 (Relaxation time (yrs))      '
!
!
DO l=lmin, lmax
! 
    write (3, '(/)')
      do k = 1, nroots
       write (3, '(i4,1x,5(e15.7,1x))') l,                  & 
                                        log10 (dfloat (l)), &
                                        s(l,k),            &
					log10(-s(l,k)),    &
				        (-1000.q0/s(l,k)), &
				        log10(-1000.q0/s(l,k))
      enddo
!
ENDDO
 close(3)
! 
 do l=lmin, lmax
!
      if(drop_modes==0) then 
	                   do k=1, nroots
	                     r_h (l,k) = h_v(l,k)/s(l,k)      ! Normalized residue for h
                         r_l (l,k) = l_v(l,k)/s(l,k)      ! Normalized residue for l
                         r_k (l,k) = k_v(l,k)/s(l,k)      ! Normalized residue for k
                         tek (l,k) = -1.0/s(l,k)     
	                   enddo 
	  endif 
!
          if(drop_modes==1)  then
!	   
             do k=1, nroots
              r_h (l,k) = 0.0
              r_l (l,k) = 0.0
              r_k (l,k) = 0.0 
!
              IF(vec(l,k)==1._sp) then
!
                r_h (l,k) = h_v(l,k)/s(l,k)      ! Normalized residue for h
                r_l (l,k) = l_v(l,k)/s(l,k)      ! Normalized residue for l
                r_k (l,k) = k_v(l,k)/s(l,k)      ! Normalized residue for k
                tek (l,k) = -1.0/s(l,k)          ! Relaxation time in k-yrs
!
              endif
!
             enddo
!
             endif 
!
 end do
!
!++++++++++
   Endif  
!++++++++++
!
          open(13,file='h.dat',status='unknown')
          open(14,file='l.dat',status='unknown')
!
  call DATE_AND_TIME (date,timc)
  Write(13,*) '# File h.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  call DATE_AND_TIME (date,timc)
  Write(14,*) '# File l.dat, created by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  	      	  
!
          if(only_elastic==0) then 
          WRITE(99,*) 'Writing the elastic, fluid, and v-elastic h ldc on h.dat '
          WRITE(99,*) 'Writing the elastic, fluid, and v-elastic l ldc on l.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic, fluid, and v-elastic h ldc on h.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic, fluid, and v-elastic l ldc on l.dat '
          endif
!
          if(only_elastic==1) then 
          WRITE(99,*) 'Writing the elastic h ldc on h.dat '
          WRITE(99,*) 'Writing the elastic l ldc on l.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic h ldc on h.dat '
IF(IV==1) &
          WRITE(* ,*) 'Writing the elastic l ldc on l.dat '
          endif
!	  
          DO l=lmin, lmax 
!
          if(only_elastic==0) then 
                    write (13, '((i3,1x,96(1x,e20.8)))') &
  	     l, h_e(l), h_f(l), (h_v (l,k), k = 1, nroots)
                    write (14, '((i3,1x,96(1x,e20.8)))') &
	     l, l_e(l), l_f(l), (l_v (l,k), k = 1, nroots)
	  endif
          if(only_elastic==1) then 
                    write (13, '((i3,1x,96(1x,e20.8)))') &
  	     l, h_e(l) 
                    write (14, '((i3,1x,96(1x,e20.8)))') &
	     l, l_e(l) 
	  endif  
!
          ENDDO
!
          close(13)
          CLOSE(14)
!
!
          WRITE(99,*) 'Done with the computation of the Spectral quantities '
IF(IV==1) WRITE(* ,*) 'Done with the computation of the Spectral quantities '
!
          WRITE(99,*) 'Looking now for KWs Load_Geometry and Load_History'
IF(IV==1) WRITE(* ,*) 'Looking now for KWs Load_Geometry and Load_History'
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ENDIF  ! Endif on KW External_Model
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
!
!
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IF(KEYWORD(1:13) =='Load_Geometry') THEN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
    iload=1
    n_load=n_load+1
!
!
              call m_2(20050)  ! Checks the kw chain
              call m_2(20060)  ! Load_Geometry must appear just one time
!
!
                WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
!
!
                Write(99,*) 'Reading the load geometry'
    IF(IV==1)   Write(*, *) 'Reading the load geometry'
!
!
!
   READ(1,*) CL 
!
!
              call m_2(20070)  ! check the range of CL
!
!
!
 If(CL <= 21) THEN
!
!
 IF( CL == 10 )      WRITE(99,*) 'Disk load (CL=10)'
 IF( CL == 10 .AND. IV==1 ) WRITE(* ,*) 'Disk load (CL=10)'
!
 IF( CL == 11 )      WRITE(99,*) 'Disk load + secondary disk load (CL=11)'
 IF( CL == 11 .AND. IV==1 ) WRITE(* ,*) 'Disk load + secondary disk load (CL=11)'
!
 IF( CL == 20 )      WRITE(99,*) 'Parabolic load (CL=20)'
 IF( CL == 20 .AND. IV==1 ) WRITE(* ,*) 'Parabolic load (CL=20)'
!
 IF( CL == 21 )      WRITE(99,*) 'Parabolic load + secondary disk load (CL=21)'
 IF( CL == 21 .AND. IV==1 ) WRITE(* ,*) 'Parabolic load + secondary disk load (CL=21)'
!
!
!	
 READ(1,*) IOC, LONG_CEN, TETA_CEN, AMPLITUDE     
!
!
!
 WRITE(99,*) 'Longitude of the centre (deg) =', long_cen
 WRITE(99,*) 'Colatitude of the centre (deg) =', teta_cen
 WRITE(99,*) 'Amplitude (deg) =', AMPLITUDE
!
 IF(IV==1) THEN
 WRITE(*,*) 'Longitude of the centre (deg) =', long_cen
 WRITE(*,*) 'Colatitude of the centre (deg) =', teta_cen
 WRITE(*,*) 'Amplitude (deg) =', AMPLITUDE
 endif
!
!
                     call m_2(20080)  ! IOC must be ==0 or ==1
                     call m_2(20090)  ! IOC=1 is not allowed for CL=11 or 21
!
!
 IF(iv==0 .AND. ioc==0) WRITE(99,*) 'The load is NOT compensated on a realistic ocean'
 IF(iv==1 .AND. ioc==0) THEN
           WRITE(99,*) 'The load is NOT compensated on a realistic ocean'
           WRITE(* ,*) 'The load is NOT compensated on a realistic ocean'
 ENDIF
 IF(iv==0 .AND. ioc==1) WRITE(99,*) 'The load IS compensated on a realistic ocean'
 IF(iv==1 .AND. ioc==1) THEN
           WRITE(99,*) 'The load IS compensated on a realistic ocean'
           WRITE(* ,*) 'The load IS compensated on a realistic ocean'
 ENDIF
!
 ENDIF   ! Endif on CL <= 21.
!   
!

!  
 IF (CL == 30) THEN
!
!
               WRITE(99,*) 'Point load (CL=30)'
   IF( IV==1 ) WRITE(* ,*) 'Point load (CL=30)'
!
!
!
   READ(1,*) IOC, LONG_CEN, TETA_CEN
!
!
!
 WRITE(99,*) 'Longitude of the point load (deg) =', long_cen
 WRITE(99,*) 'Colatitude of the point load (deg) =', teta_cen
!
 IF(IV==1) THEN
 WRITE(*,*) 'Longitude of the centre (deg) =', long_cen
 WRITE(*,*) 'Colatitude of the centre (deg) =', teta_cen
 endif
!
                   call m_2(20080)  ! IOC must be ==0 or ==1
!
!
!
 IF(iv==0 .AND. ioc==0) WRITE(99,*) 'The load is NOT compensated on a realistic ocean'
 IF(iv==1 .AND. ioc==0) THEN
           WRITE(99,*) 'The load is NOT compensated on a realistic ocean'
           WRITE(* ,*) 'The load is NOT compensated on a realistic ocean'
 ENDIF
 IF(iv==0 .AND. ioc==1) WRITE(99,*) 'The load IS compensated on a realistic ocean'
 IF(iv==1 .AND. ioc==1) THEN
           WRITE(99,*) 'The load IS compensated on a realistic ocean'
           WRITE(* ,*) 'The load IS compensated on a realistic ocean'
 ENDIF
!
 ENDIF   ! Endif on CL <= 30.
!
!
!
!
   If(CL == 40) THEN
!
                    WRITE(99,*) 'Harmonic load (CL=40)'
   IF( IV==1 ) WRITE(* ,*) 'Harmonic load (CL=40)'
!
!
!
   READ(1,*) IOC, LONG_CEN, TETA_CEN 
!
!
!
 WRITE(99,*) 'Longitude of the load pole (deg) =', long_cen
 WRITE(99,*) 'Colatitude of the load pole (deg) =', teta_cen
!
 IF(IV==1) THEN
 WRITE(*,*) 'Longitude of the load pole (deg) =', long_cen
 WRITE(*,*) 'Colatitude of the load pole (deg) =', teta_cen
 endif
!
!
                 call m_2(20100) ! l_min must be = l_max if cl=40
                 call m_2(20110) ! IOC /= 0 is not allowed for CL=40
!
!
 IF(iv==0 .AND. ioc==0) WRITE(99,*) 'The load is NOT compensated on a realistic ocean'
 IF(iv==1 .AND. ioc==0) THEN
           WRITE(99,*) 'The load is NOT compensated on a realistic ocean'
           WRITE(* ,*) 'The load is NOT compensated on a realistic ocean'
 ENDIF
 IF(iv==0 .AND. ioc==1) WRITE(99,*) 'The load IS compensated on a realistic ocean'
 IF(iv==1 .AND. ioc==1) THEN
           WRITE(99,*) 'The load IS compensated on a realistic ocean'
           WRITE(* ,*) 'The load IS compensated on a realistic ocean'
 ENDIF
!
 ENDIF   ! Endif on CL = 40.
!
!
!
!
   If(CL == 50)       THEN
!   
               WRITE(99,*) 'Rectangular load (CL=50)'
   IF( IV==1 ) WRITE(* ,*) 'Rectangular load (CL=50)'
!
!
!
     READ(1,*) IOC, LONG_CEN, DL, TETA_CEN, DT
!
!
  WRITE(99,*) 'Longitude  of the centroid (deg) =', long_cen
  WRITE(99,*) 'Colatitude of the centroid (deg) =', teta_cen
  WRITE(99,*) 'Amplitude in longitude  (deg) =', DL
  WRITE(99,*) 'Amplitude in colatitude (deg) =', DT
!
  IF(IV==1) THEN
  WRITE(*,*) 'Longitude  of the centroid (deg) =', long_cen
  WRITE(*,*) 'Colatitude of the centroid (deg) =', teta_cen
  WRITE(*,*) 'Amplitude in longitude  (deg) =', DL
  WRITE(*,*) 'Amplitude in colatitude (deg) =', DT
  ENDIF
!
!
  call m_2(20080)  ! IOC must be ==0 or ==1
!
!
 IF(iv==0 .AND. ioc==0) WRITE(99,*) 'The load is NOT compensated on a realistic ocean'
 IF(iv==1 .AND. ioc==0) THEN
           WRITE(99,*) 'The load is NOT compensated on a realistic ocean'
           WRITE(* ,*) 'The load is NOT compensated on a realistic ocean'
 ENDIF
 IF(iv==0 .AND. ioc==1) WRITE(99,*) 'The load IS compensated on a realistic ocean'
 IF(iv==1 .AND. ioc==1) THEN
           WRITE(99,*) 'The load IS compensated on a realistic ocean'
           WRITE(* ,*) 'The load IS compensated on a realistic ocean'
 ENDIF
!
 ENDIF   ! Endif on CL = 50.
!
!
!
!
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ENDIF   ! Endif on 'Load_Geometry'     
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!
!
!
!
!
!
!
!
!
!
!    
!******************************************
   IF(KEYWORD(1:13)=='Load_History') THEN
!******************************************
!
   n_hist=n_hist+1
!
   ihist = 1 
!
           call m_2(20120) ! Load_Geometry must be active
           call m_2(20130) ! Load_history must appear one time
!
                     WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
                     Write(99,*) 'Reading the load time-history '
    IF(IV==1)   Write(*, *) 'Reading the load time-history '
!
!
    READ  (1,*) lh
!     
!
          call m_2(20140) ! LH must be in the range [0:8]
!
!
     IF(lh == 0 .OR. lh == 1) THEN      
! 
  IF(lh == 0) WRITE(99,*) 'Time-history #0: Instantaneous loading at t=0 kyrs'
  IF(lh == 0 .and. IV==1) &
                   WRITE (*,*) 'Time-history #0: Instantaneous loading at t=0 kyrs'
  IF(lh == 1) WRITE(99,*) 'Time-history #1: Instantaneous Unloading at t=0 kyrs'
  IF(lh == 1 .and. IV==1) &
                   WRITE (*,*) 'Time-history #1: Instantaneous Unloading at t=0 kyrs'
!
!
!
        READ(1,*) JUNK    
!         
        If(CL <= 21 .or. CL == 50) HEIGHT   = JUNK 
	If(CL == 30)                        MASS   = JUNK 
	If(CL == 40)                       PRESS   = JUNK
!
!
!
        call m_2(20150) ! The thickness (mass, load pressure) must be > 0.0
!
!
  ENDIF  ! Endif on lh =0 or 1  <<<<<<<<<<<
!
!
!
!
     IF(lh==2) THEN 
! 
!
     WRITE (99,*)      'Time-history #2: Loading & Unloading'
     IF(IV==1) &
     WRITE  (*,*)      'Time-history #2: Loading & Unloading'
!
!         
     READ(1,*) TAU_LU
!
!
      call m_2(20160) ! The phase with load /= 0.0 must be > 0 kyrs long
!
!
        IF(iv==0) WRITE(99,*)  'The loading occurs at time', TAU_LU, 'kyrs'
        IF(iv==1) THEN
                   WRITE(99,*) 'The loading occurs at time', TAU_LU, 'kyrs'
                   WRITE(* ,*) 'The loading occurs at time', TAU_LU, 'kyrs'
                  ENDIF
        IF(iv==0) WRITE(99,*)  'The un-loading occurs at time', float(iv), 'kyrs'
        IF(iv==1) THEN
                   WRITE(99,*) 'The un-loading occurs at time', float(iv)-1., 'kyrs'
                   WRITE(* ,*) 'The un-loading occurs at time', float(iv)-1., 'kyrs'
                  ENDIF
!
     READ(1,*) JUNK
!
!
        If(CL <= 21 .or. CL == 50) HEIGHT   = JUNK
	If(CL == 30)                        MASS   = JUNK 
	If(CL == 40)                       PRESS   = JUNK
!
!
        call m_2(20150) ! The thickness (mass, load pressure) must be > 0.0
!
!
!
  ENDIF  ! Endif on NEW lh =2   <<<<<<<<<
!
!
!
!
!
!
     IF(lh==3) THEN
! 
!
     WRITE (99,*)      'Time-history #3: Simple melting'
     IF(IV==1) &
     WRITE  (*,*)      'Time-history #3: Simple melting'
!
!         
        READ(1,*) TAU_SD
!
!
        call m_2(20170)  ! Tau_sd must be > 0.0
!
!
        IF(iv==0) WRITE(99,*) 'The deglaciation is ', tau_sd, ' kyrs long'
        IF(iv==1) then
               WRITE(99,*) 'The deglaciation is ', tau_sd, ' kyrs long'
               WRITE(* ,*) 'The deglaciation is ', tau_sd, ' kyrs long'
        endif
!
        READ(1,*) JUNK
!
!
        call m_2(20150) ! The thickness (mass, load pressure) must be > 0.0
!
!
!
        If(CL <= 21 .or. CL == 50) HEIGHT   = JUNK
	If(CL == 30)                        MASS   = JUNK 
	If(CL == 40)                       PRESS   = JUNK
!
!
!
  ENDIF  ! Endif on lh =3    <<<<<<<<<
!
!
!
!	
!
     IF(lh==4) THEN
!
     WRITE (99,*)      'Time-history #4: Saw-tooth time history'
     IF(IV==1) &
     WRITE  (*,*)      'Time-history #4: Saw-tooth time history'
!
!          
        READ(1,*) NR
!
!
        call m_2(20180)  ! NR must be in the range [0:5]
!
!
        IF(iv==0) WRITE(99,*)  'There are', NR+1, ' glacial phases'
        IF(iv==1) THEN
                   WRITE(99,*) 'There are', NR+1, ' glacial phases'
                   WRITE(* ,*) 'There are', NR+1, ' glacial phases'
                  ENDIF
!
        READ(1,*) TAUC
!
!
        call m_2(20190)  ! TAUC must be > 0.0
!
!
        IF(iv==0) WRITE(99,*)  'The loading phase is ', TAUC, 'kyrs long'
        IF(iv==1) THEN
                   WRITE(99,*) 'The loading phase is ', TAUC, 'kyrs long'
                   WRITE(* ,*) 'The loading phase is ', TAUC, 'kyrs long'
                  ENDIF
!
!
        READ(1,*) DINC
!
!
        call m_2(20200)  ! DINC must be > 0.0
!
!
        IF(iv==0) WRITE(99,*)  'The un-loading phase is ', dinc, 'kyrs long'
        IF(iv==1) THEN
                   WRITE(99,*) 'The un-loading phase is ', dinc, 'kyrs long'
                   WRITE(* ,*) 'The un-loading phase is ', dinc, 'kyrs long'
                  ENDIF
!
!
        READ(1,*) JUNK    
!         
        If(CL <= 21 .or. CL == 50)   HEIGHT = JUNK
	If(CL == 30)                        MASS   = JUNK 
	If(CL == 40)                        PRESS  = JUNK    
!
!
!
        call m_2(20150) ! The thickness (mass, load pressure) must be > 0.0
!
!
  ENDIF  ! Endif on lh =4   <<<<<<<<<
!
!
!
!
     IF(lh==5) THEN
!
     WRITE (99,*)      'Time-history #5: Periodic Loading'
     IF(IV==1) &
     WRITE  (*,*)      'Time-history #5: Periodic Loading'
!
!          
        READ(1,*) Period
!
!
!
        call m_2(20210)  ! the period must be > 0.0
!
!
!
        IF(iv==0) WRITE(99,*)  'The period is ', Period, 'kyrs'
        IF(iv==1) THEN
                   WRITE(99,*) 'The period is ', Period, 'kyrs'
                   WRITE(* ,*) 'The period is ', Period, 'kyrs'
                  ENDIF
!
        READ(1,*) JUNK    
!         
        If(CL <= 21 .or. CL == 50)   HEIGHT = JUNK
	If(CL == 30)                        MASS   = JUNK 
	If(CL == 40)                        PRESS  = JUNK    
!
!
!
        call m_2(20150) ! The thickness (mass, load pressure) must be > 0.0
!
!
!
!
  ENDIF  ! Endif on lh =5  <<<<<<<<<
!
!    
!
!
!
!

     IF(lh==6) THEN
!    
     WRITE (99,*)      'Time-history #6: Piecewise linear'
     IF(IV==1) &
     WRITE  (*,*)      'Time-history #6: Piecewise linear'
!
!
     Open (2,FILE='timeh_6.dat',status='unknown')
!			 
         JUNK=0.0
!
     READ(2,*) ijunk
!
!
     call m_2(20220)  ! checking the first line of timeh_6.dat
!
!
   DO I=0, nss_max
	 READ  (2,*,END=8002)  T_HH(I), A_HH(I)
	 NSS=I
!
!
     call m_2(20230) ! t0 must be equal to zero
     call m_2(20240) ! t(i) must be > 0.0
     call m_2(20250) ! t(i) must form an increasing sequence
     call m_2(20260) ! the thickness (mass, load pressure) must be > 0.0
!

  IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) &
            write(99,*) 'Time (ka) and load thickness (m)', &
                                 t_hh(i),a_hh(i)
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) THEN
            write(* ,*) 'Time (ka) and load thickness (m)', &
                                 t_hh(i),a_hh(i)
            write(99,*) 'Time (ka) and load thickness (m)', &
                                 t_hh(i),a_hh(i)
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 30) &
            write(99,*) 'Time (ka) and load mass (kg)', &
                                  t_hh(i),a_hh(i)
  IF(iv==1 .AND. CL == 30)                       THEN
            write(*, *) 'Time (ka) and load mass (kg)', &
                                  t_hh(i),a_hh(i)
            write(99,*) 'Time (ka) and load mass (kg)', &
                                  t_hh(i),a_hh(i)
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 40) &
write(99,*) 'Time (ka) and load mass/surface (kg/m**2)', &
                              t_hh(i),a_hh(i)
  IF(iv==1 .AND. CL == 40)                       THEN
write(* ,*) 'Time (ka) and load mass/surface (kg/m**2)', &
                              t_hh(i),a_hh(i)
write(99,*) 'Time (ka) and load mass/surface (kg/m**2)', &
                              t_hh(i),a_hh(i)
                                                        ENDIF
!
!
!
       IF(A_HH(I) >= JUNK) JUNK=A_HH(I)
!
!
       ENDDO
!
 8002  CONTINUE
!
       CLOSE(2)
!
!
       call m_2(20270) ! At least two points must be given in file timeh_6.dat
       call m_2(20280) ! nss must be <= nss_max
       call m_2(20290) ! The maximum thickness (mass, load pressure) must be > 0.0
!
!
  IF(iv==0) WRITE(99,*) 'The th#6 is piecewise linear over ',nss, ' time interval(s)'
  IF(iv==1) THEN
            WRITE(* ,*) 'The th#6 is piecewise linear over ',nss, ' time interval(s)'
            WRITE(99,*) 'The th#6 is piecewise linear over ',nss, ' time interval(s)'
  ENDIF
!
!
!
       tnn = t_hh(nss)
!
!
!
  IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) THEN
    WRITE(99,*) 'For -oo<=time<0 kyr the load thickness is ', a_hh(0), ' m'
    WRITE(99,*) 'For time>=', tnn, 'kyr the load thickness is ', a_hh(nss), ' m'
    ENDIF
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) then
    WRITE(* ,*) 'For -oo<=time<0 kyr the load thickness is ', a_hh(0), ' m'
    WRITE(* ,*) 'For time>=', tnn, 'kyr the load thickness is ', a_hh(nss), ' m'
    WRITE(99,*) 'For -oo<=time<0 kyr the load thickness is ', a_hh(0), ' m'
    WRITE(99,*) 'For time>=', tnn, 'kyr the load thickness is ', a_hh(nss), ' m'
  endif
!
  IF(iv==0 .AND. CL == 30) THEN
    WRITE(99,*) 'For -oo<=time<0 kyr the load mass is ', a_hh(0), ' kg'
    WRITE(99,*) 'For time>=', tnn, 'kyr the load mass is ', a_hh(nss), ' kg'
    ENDIF
  IF(iv==1 .AND. CL == 30) then
    WRITE(* ,*) 'For -oo<=time<0 kyr the load mass is ', a_hh(0), ' kg'
    WRITE(* ,*) 'For time>=', tnn, 'kyr the load mass is ', a_hh(nss), ' kg'
    WRITE(99,*) 'For -oo<=time<0 kyr the load mass is ', a_hh(0), ' kg'
    WRITE(99,*) 'For time>=', tnn, 'kyr the load mass is ', a_hh(nss), ' kg'
  endif
!
  IF(iv==0 .AND. CL == 40) THEN
    WRITE(99,*) 'For -oo<=time<0 kyr the load mass is ', a_hh(0), ' kg'
    WRITE(99,*) 'For time>=', tnn, 'kyr the load mass is ', a_hh(nss), ' kg'
    ENDIF
  IF(iv==1 .AND. CL == 40) then
WRITE(* ,*) 'For -oo<=time<0 kyr the mass/surface is',a_hh(0),'kg/m**2'
WRITE(99,*) 'For -oo<=time<0 kyr the mass/surface is',a_hh(0),'kg/m**2'
WRITE(* ,*) 'For time>=', tnn, 'kyr the mass/surface is',a_hh(nss),'kg/m**2'
WRITE(99,*) 'For time>=', tnn, 'kyr the mass/surface is',a_hh(nss),'kg/m**2'
  endif
!
!
!
!
! 
  IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) &
            write(99,*) 'Maximum load thickness (m)', junk
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) THEN
            write(* ,*) 'Maximum load thickness (m)', junk
            write(99,*) 'Maximum load thickness (m)', junk
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 30) &
            write(99,*) 'Maximum mass (kg)', junk
  IF(iv==1 .AND. CL == 30)                       THEN
            write(* ,*) 'Maximum mass (kg)', junk
            write(99,*) 'Maximum mass (kg)', junk
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 40) &
write(99,*) 'Maximum mass/surface (kg/m**2)', junk
  IF(iv==1 .AND. CL == 40)                       THEN
write(* ,*) 'Maximum mass/surface (kg/m**2)', junk
write(99,*) 'Maximum mass/surface (kg/m**2)', junk
                                                        ENDIF
!
                          DO I=0, NSS 
			     IF(JUNK/=0.0) THEN
			     A_HH(I) = A_HH(I)/JUNK
			                      ELSE 
		             A_HH(I) = 0.0
			                      ENDIF
!
  IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) &
            write(99,*) &
            'Time (ka) and normalized load thickness (m)', t_hh(i),a_hh(i)
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) THEN
            write(* ,*) &
            'Time (ka) and normalized load thickness (m)',t_hh(i),a_hh(i)
            write(99,*) &
            'Time (ka) and normalized load thickness (m)',t_hh(i),a_hh(i)
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 30) &
            write(99,*) &
            'Time (ka) and normalized load mass (kg)', t_hh(i),a_hh(i)
  IF(iv==1 .AND. CL == 30)                       THEN
            write(* ,*) &
            'Time (ka) and normalized load mass (kg)', t_hh(i),a_hh(i)
            write(99,*) &
            'Time (ka) and normalized load mass (kg)', t_hh(i),a_hh(i)
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 40) write(99,*) &
       'Time (ka) and normalized mass/surface (kg/m**2)', t_hh(i),a_hh(i)
  IF(iv==1 .AND. CL == 40)                       THEN
write(* ,*) &
       'Time (ka) and normalized mass/surface (kg/m**2)', t_hh(i),a_hh(i)
write(99,*) &
       'Time (ka) and normalized mass/surface (kg/m**2)', t_hh(i),a_hh(i)
                                                        ENDIF
!
!
			  ENDDO 
!
!
        If(CL <= 21 .or. CL == 50) HEIGHT   = JUNK
	If(CL == 30)                        MASS   = JUNK 
	If(CL == 40)                        PRESS  = JUNK         
!     
!
!
! -------------------------------------------
! Auxiliary arrays useful for the convolution
! -------------------------------------------
!
  do i=1, nss
               ri(i) = (a_hh(i) - a_hh(i-1))/ (t_hh(i) - t_hh(i-1))
  end do
!
  alf(0)   = (a_hh(1)-a_hh(0)) - ri(1)*t_hh(1)
  bet(0)   = ri(1)
!
  alf(nss) = ri(nss)*t_hh(nss)
  bet(nss) = -ri(nss)
!
   do k=1, nss-1
      alf(k) = (a_hh(k+1)-a_hh(k)) + ri(k)*t_hh(k) - ri(k+1)*t_hh(k+1)
      bet(k) = ri(k+1)-ri(k)
   end do
!
!
  IF(iv==0) WRITE(99,*) 'Computing the auxiliary arrays... done'
  IF(iv==1) THEN
         WRITE(99,*) 'Computed the auxiliary arrays... done'
         WRITE(* ,*) 'Computed the auxiliary arrays... done'
  ENDIF
!
!
        ENDIF       ! <--- Endif on  lh = 6  <<<<<<<<<
	
!     
!     
!     
!          
!
!
     IF(lh==7) THEN
!
!
     WRITE (99,*)      'Time-history #7: Piecewise constant'
     IF(IV==1) &
     WRITE  (*,*)      'Time-history #7: Piecewise constant'
!
!
     Open (2,FILE='timeh_7.dat',status='unknown')
!			 
     JUNK=0.0
!
     READ(2,*) ijunk
!
!
         call m_2(20300) ! File timeh_7.dat must begin with '7'
!
     READ(2,*) dilta
!
!
         call m_2(20310) ! There are bounds on dilta
!
!
     JUNK =  -9999.
!
!
	DO L=0, 10001
!
           READ  (2,*,END=8001)  I, A_H(I)
!
           call m_2(20320) ! File timeh_7.dat should be written correctly
           call m_2(20330) ! The thickness (mass, load pressure) must be >= 0.0
!
           IF(A_H(I) >= JUNK) JUNK=A_H(I)
!
           ns = i
!
           ENDDO
!
 8001  continue
!
       CLOSE(2)
!
          call m_2(20340) ! At least two points must be given in file timeh_7.dat
          call m_2(20350) ! ns must be <= ns_max
          call m_2(20360) ! The maximum thickness (mass, load pressure) must be > 0.0
!
!
               a_h(ns+1)=a_h(ns)
!
               tn = float(ns)*dilta
!
!
!
!
  IF(iv==0) then
  WRITE(99,*) 'There are',ns,'steps between t = 0 kyrs and time t=',tn,'kyrs'
  ENDIF
  IF(iv==1) THEN
  WRITE(99,*) 'There are',ns,'steps between t = 0 kyrs and time t=',tn,'kyrs'
  WRITE(* ,*) 'There are',ns,'steps between t = 0 kyrs and time t=',tn,'kyrs'
  ENDIF
!
  IF(iv==0) THEN
write(99,*) 'Lenght of one these steps =', dilta, 'kyrs'
write(99,*) 'Lenght of the deglaciation phase=', tn, 'kyrs'
            ENDIF
  IF(iv==1) THEN
write(99,*) 'Lenght of one of these steps =', dilta, 'kyrs'
write(99,*) 'Lenght of the deglaciation phase=', tn, 'kyrs'
write(* ,*) 'Lenght of one of these steps =', dilta, 'kyrs'
write(* ,*) 'Lenght of the deglaciation phase=', tn, 'kyrs'
            ENDIF
!
IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) THEN
WRITE(99,*) 'For -oo<=t<0 kyr the load thickness is',  A_H(0), 'm (step#0)'
WRITE(99,*) 'For t>=',tn,'kyr the load thickness is', a_h(ns+1), 'm (step#',ns+1,')'
   ENDIF
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) then
WRITE(99,*) 'For -oo<=t<0 kyr the load thickness is',  A_H(0), 'm (step#0)'
WRITE(99,*) 'For t>=',tn,'kyr the load thickness is', a_h(ns+1), 'm (step#',ns+1,')'
WRITE(*,*)  'For -oo<=t<0 kyr the load thickness is',  A_H(0), 'm (step#0)'
WRITE(*,*)  'For t>=',tn,'kyr the load thickness is', a_h(ns+1), 'm (step#',ns+1,')'
  ENDIF
!
  IF(iv==0 .AND. CL == 30) THEN
WRITE(99,*) 'For -oo<=t<0 kyr the load mass is',  A_H(0), 'kg (step#0)'
WRITE(99,*) 'For t>=',tn,'kyr the load mass is', a_h(ns+1), 'm (step#',ns+1,')'
   ENDIF
  IF(iv==1 .AND. CL == 30) then
WRITE(*,*)  'For -oo<=t<0 kyr the load mass is',  A_H(0), 'kg (step#0)'
WRITE(*,*) 'For t>=',tn,'kyr the load mass is', a_h(ns+1), 'kg (step#',ns+1,')'
WRITE(99,*) 'For -oo<=t<0 kyr the load mass is',  A_H(0), 'kg (step#0)'
WRITE(99,*) 'For t>=',tn,'kyr the load mass is', a_h(ns+1), 'kg (step#',ns+1,')'
  ENDIF
!
  IF(iv==0 .AND. CL == 40) THEN
WRITE(99,*) 'For -oo<=t<0 kyr the mass/surface is', A_H(0),'kg/m**2 (step#0)'
WRITE(99,*) 'For t>=',tn,'kyr the mass/surface is', a_h(ns+1),'kg/m**2 (step#',ns+1,')'
   ENDIF
  IF(iv==1 .AND. CL == 40) THEN
WRITE(99,*) 'For -oo<=t<0 kyr the mass/surface is', A_H(0),'kg/m**2 (step#0)'
WRITE(99,*) 'For t>=',tn,'kyr the mass/surface is', a_h(ns+1),'kg/m**2 (step#',ns+1,')'
WRITE(* ,*) 'For -oo<=t<0 kyr the mass/surface is', A_H(0),'kg/m**2 (step#0)'
WRITE(* ,*) 'For t>=',tn,'kyr the mass/surface is', a_h(ns+1),'kg/m**2 (step#',ns+1,')'
  ENDIF
!
              WRITE(99,*) 'In detail:'
    IF(iv==1) WRITE(99,*) 'In detail:'
!
  do i=0,ns+1
!
  IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) &
         write(99,*) 'Step# and load thickness (m)', i,a_h(i)
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) THEN
         write(* ,*) 'Step# and load thickness (m)', i,a_h(i)
         write(99,*) 'Step# and load thickness (m)', i,a_h(i)
                                                        ENDIF
  IF(iv==0 .AND. CL == 30) &
             write(99,*) 'Step# and load mass (kg)', i,a_h(i)
  IF(iv==1 .AND. CL == 30)                       THEN
             write(* ,*) 'Step# and load mass (kg)', i,a_h(i)
             write(99,*) 'Step# and load mass (kg)', i,a_h(i)
                                                        ENDIF
  IF(iv==0 .AND. CL == 40) &
write(99,*) 'Step# and load mass/surface (kg/m**2)', i,a_h(i)
  IF(iv==1 .AND. CL == 40)                       THEN
write(* ,*) 'Step# and load mass/surface (kg/m**2)', i,a_h(i)
write(99,*) 'Step# and load mass/surface (kg/m**2)', i,a_h(i)
                                                        ENDIF
!
 end do
!
!
  IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) &
            write(99,'(a27,2(f10.4,1x))') 'Maximum load thickness (m)', junk
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) THEN
            write(* ,'(a27,2(f10.4,1x))') 'Maximum load thickness (m)', junk
            write(99,'(a27,2(f10.4,1x))') 'Maximum load thickness (m)', junk
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 30) &
            write(99,'(a18,2(f10.4,1x))') 'Maximum mass (kg)', junk
  IF(iv==1 .AND. CL == 30)                       THEN
            write(* ,'(a18,2(f10.4,1x))') 'Maximum mass (kg)', junk
            write(99,'(a18,2(f10.4,1x))') 'Maximum mass (kg)', junk
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 40) &
write(99,'(a31,2(f10.4,1x))') 'Maximum mass/surface (kg/m**2)', junk
  IF(iv==1 .AND. CL == 40)                       THEN
write(* ,'(a31,2(f10.4,1x))') 'Maximum mass/surface (kg/m**2)', junk
write(99,'(a31,2(f10.4,1x))') 'Maximum mass/surface (kg/m**2)', junk
                                                        ENDIF
!
!
!
!			  
                          DO I=0, NS+1
			     IF(JUNK/=0.0) THEN
			     A_H(I) = A_H(I)/JUNK
			                      ELSE 
		             A_H(I) = 0.0
			                      ENDIF
!
  IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) &
  write(99,*) 'Step# and normalized load thickness (m)',I,a_h(i)
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) THEN
  write(99,*) 'Step# and normalized load thickness (m)',I,a_h(i)
  write(* ,*) 'Step# and normalized load thickness (m)',I,a_h(i)
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 30) &
  write(99,*) 'Step# and normalized load mass (kg)', i,a_h(i)
  IF(iv==1 .AND. CL == 30)                       THEN
  write(99,*) 'Step# and normalized load mass (kg)', i,a_h(i)
  write(* ,*) 'Step# and normalized load mass (kg)', i,a_h(i)
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 40) &
  write(99,*)'Step# and normalized mass/surface (kg/m**2)',i,a_h(i)
  IF(iv==1 .AND. CL == 40)                       THEN
  write(99,*)'Step# and normalized mass/surface (kg/m**2)',i,a_h(i)
  write(* ,*)'Step# and normalized mass/surface (kg/m**2)',i,a_h(i)
                                                        ENDIF
!
!
			  ENDDO 
!
!
        If(CL <= 21 .or. CL == 50) HEIGHT   = JUNK
	If(CL == 30)                        MASS   = JUNK 
	If(CL == 40)                        PRESS  = JUNK         
!     
!
!     
        ENDIF       ! Endif On lh = 7  <<<<<<<<<
!
!	
!     
!
!
     IF(lh==8) THEN
!
!
     WRITE (99,*)  'Time-history #8: Piecewise constant with finite loading phase'
     IF(IV==1) &
     WRITE  (*,*)  'Time-history #8: Piecewise constant with finite loading phase'
!
!
     Open (2,FILE='timeh_8.dat',status='unknown')
!			 
     JUNK=0.0
!
     READ(2,*) ijunk
!
                call m_2(20370)   ! File timeh_8.dat must begin with '8'
!
!
!
     READ(2,*) dilta
!
                call m_2(20380)   ! There are bounds on dilta
!
!
     READ(2,*) TAUF
!
                call m_2(20390)   ! There are bounds on TAUF
!
!
     JUNK =  -9999.     
!
!
!			 
	DO L=0, 10001
		READ  (2,*,END=8003)  I, A_H(I)
!
                call m_2(20400)   ! File timeh_8.dat should be written correctly
                call m_2(20410)   ! The thickness (mass, load pressure) must be >= 0.0
!
                IF(A_H(I) >= JUNK) JUNK=A_H(I)
!
                ns = i
!
                ENDDO
!
 8003   continue
!
                      CLOSE(2)
!
!
         call m_2(20420)  ! At least two points must be given in file timeh_8.dat
         call m_2(20430)  ! ns must be <= ns_max
         call m_2(20440)  ! The maximum thickness (mass, load pressure) must be > 0.0
!
!
                      a_h(ns+1)=a_h(ns)
!
                      tn = float(ns)*dilta
!
!
!
!
  IF(iv==0) then
  WRITE(99,*) 'There are',ns,'deglaciation steps between t = 0 and t=',tn,'kyrs'
  ENDIF
  IF(iv==1) THEN
  WRITE(99,*) 'There are',ns,'deglaciation steps between t = 0 and t=',tn,'kyrs'
  WRITE(* ,*) 'There are',ns,'deglaciation steps between t = 0 and t=',tn,'kyrs'
  ENDIF
!
  IF(iv==0)                       THEN
write(99,*) 'Lenght of one step =', dilta, 'kyrs'
write(99,*) 'Lenght of the deglaciation phase=', tn, 'kyrs'
write(99,*) 'Lenght of the finite loading phase=', tauf, 'kyrs'
                                  ENDIF
  IF(iv==1)                       THEN
write(99,*) 'Lenght of one step =', dilta, 'kyrs'
write(99,*) 'Lenght of the deglaciation phase=', tn, 'kyrs'
write(99,*) 'Lenght of the finite loading phase=', tauf, 'kyrs'
write(* ,*) 'Lenght of one step =', dilta, 'kyrs'
write(* ,*) 'Lenght of the deglaciation phase=', tn, 'kyrs'
write(* ,*) 'Lenght of the finite loading phase=', tauf, 'kyrs'
                                  ENDIF
!
  IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) THEN
    write(99,*) 'At the end of the loading phase ... '
    WRITE(99,*) '... the load thickness is', a_h(0), '(m) (step#0)'
                 ENDIF
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) THEN
    write(99,*) 'At the end of the loading phase ... '
    WRITE(99,*) '... the load thickness is', a_h(0), '(m) (step#0)'
    write(* ,*) 'At the end of the loading phase ... '
    WRITE(* ,*) '... the load thickness is', a_h(0), '(m) (step#0)'
                 ENDIF
!
  IF(iv==0 .AND. CL == 30) THEN
    write(99,*) 'At the end of the loading phase ...'
    WRITE(99,*) '... the load mass is',a_h(0), 'kg (step#0)'
                 ENDIF
  IF(iv==1 .AND. CL == 30) THEN
    write(99,*) 'At the end of the loading phase ...'
    WRITE(99,*) '... the load mass is',a_h(0), 'kg (step#0)'
    write(* ,*) 'At the end of the loading phase ...'
    WRITE(* ,*) '... the load mass is',a_h(0), 'kg (step#0)'
                 ENDIF
!
  IF(iv==0 .AND. CL == 50) THEN
    write(99,*) 'At the end of the finite loading phase the mass/surface ...'
    WRITE(99,*) '... at the pole of the harmonic load is ',a_h(0), ' kg/m**2 (step#0)'
                 ENDIF
  IF(iv==1 .AND. CL == 50) THEN
    write(99,*) 'At the end of the finite loading phase the mass/surface...'
    WRITE(99,*) '...at the pole of the harmonic load is ',a_h(0), ' kg/m**2 (step#0)'
    write(* ,*) 'At the end of the finite loading phase the mass/surface...'
    WRITE(* ,*) '...at the pole of the harmonic load is ',a_h(0), ' kg/m**2 (step#0)'
                 ENDIF
!
  do i=0,ns+1
!
  IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) &
         write(99,*) 'Step# and load thickness (m)', i,a_h(i)
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) THEN
         write(* ,*) 'Step# and load thickness (m)', i,a_h(i)
         write(99,*) 'Step# and load thickness (m)', i,a_h(i)
                                                        ENDIF
  IF(iv==0 .AND. CL == 30) &
             write(99,*) 'Step# and load mass (kg)', i,a_h(i)
  IF(iv==1 .AND. CL == 30)                       THEN
             write(* ,*) 'Step# and load mass (kg)', i,a_h(i)
             write(99,*) 'Step# and load mass (kg)', i,a_h(i)
                                                        ENDIF
  IF(iv==0 .AND. CL == 40) &
write(99,*) 'Step# and load mass/surface (kg/m**2)', i,a_h(i)
  IF(iv==1 .AND. CL == 40)                       THEN
write(* ,*) 'Step# and load mass/surface (kg/m**2)', i,a_h(i)
write(99,*) 'Step# and load mass/surface (kg/m**2)', i,a_h(i)
                                                        ENDIF
!
 end do
!
!
IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) THEN
WRITE(99,*) 'For t>=',tn,'kyr the load thickness is', a_h(ns+1), 'm (step#',ns+1,')'
   ENDIF
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) then
WRITE(99,*) 'For t>=',tn,'kyr the load thickness is', a_h(ns+1), 'm (step#',ns+1,')'
WRITE(*,*)  'For t>=',tn,'kyr the load thickness is', a_h(ns+1), 'm (step#',ns+1,')'
  ENDIF
!
  IF(iv==0 .AND. CL == 30) THEN
WRITE(99,*) 'For t>=',tn,'kyr the load mass is', a_h(ns+1), 'm (step#',ns+1,')'
   ENDIF
  IF(iv==1 .AND. CL == 30) then
WRITE(*,*)  'For t>=',tn,'kyr the load mass is', a_h(ns+1), 'kg (step#',ns+1,')'
WRITE(99,*) 'For t>=',tn,'kyr the load mass is', a_h(ns+1), 'kg (step#',ns+1,')'
  ENDIF
!
  IF(iv==0 .AND. CL == 40) THEN
WRITE(99,*) 'For t>=',tn,'kyr the mass/surface is', a_h(ns+1),'kg/m**2 (step#',ns+1,')'
   ENDIF
  IF(iv==1 .AND. CL == 40) THEN
WRITE(99,*) 'For t>=',tn,'kyr the mass/surface is', a_h(ns+1),'kg/m**2 (step#',ns+1,')'
WRITE(* ,*) 'For t>=',tn,'kyr the mass/surface is', a_h(ns+1),'kg/m**2 (step#',ns+1,')'
  ENDIF
!
!
  IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) &
            write(99,'(a27,2(f10.4,1x))') 'Maximum load thickness (m)', junk
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) THEN
            write(* ,'(a27,2(f10.4,1x))') 'Maximum load thickness (m)', junk
            write(99,'(a27,2(f10.4,1x))') 'Maximum load thickness (m)', junk
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 30) &
            write(99,'(a18,2(f10.4,1x))') 'Maximum mass (kg)', junk
  IF(iv==1 .AND. CL == 30)                       THEN
            write(* ,'(a18,2(f10.4,1x))') 'Maximum mass (kg)', junk
            write(99,'(a18,2(f10.4,1x))') 'Maximum mass (kg)', junk
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 40) &
write(99,'(a31,2(f10.4,1x))') 'Maximum mass/surface (kg/m**2)', junk
  IF(iv==1 .AND. CL == 40)                       THEN
write(* ,'(a31,2(f10.4,1x))') 'Maximum mass/surface (kg/m**2)', junk
write(99,'(a31,2(f10.4,1x))') 'Maximum mass/surface (kg/m**2)', junk
                                                        ENDIF
!
!
!
                          DO I=0, NS+1
			     IF(JUNK/=0.0) THEN
			     A_H(I) = A_H(I)/JUNK
			                      ELSE 
		             A_H(I) = 0.0
			                      ENDIF
!
  IF(iv==0 .AND. (CL <=21 .OR. CL == 50)) &
  write(99,*) 'Step# and normalized load thickness (m)',I,a_h(i)
  IF(iv==1 .AND. (CL <=21 .OR. CL == 50)) THEN
  write(99,*) 'Step# and normalized load thickness (m)',I,a_h(i)
  write(* ,*) 'Step# and normalized load thickness (m)',I,a_h(i)
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 30) &
  write(99,*) 'Step# and normalized load mass (kg)', i,a_h(i)
  IF(iv==1 .AND. CL == 30)                       THEN
  write(99,*) 'Step# and normalized load mass (kg)', i,a_h(i)
  write(* ,*) 'Step# and normalized load mass (kg)', i,a_h(i)
                                                        ENDIF
!
  IF(iv==0 .AND. CL == 40) &
  write(99,*)'Step# and normalized mass/surface (kg/m**2)',i,a_h(i)
  IF(iv==1 .AND. CL == 40)                       THEN
  write(99,*)'Step# and normalized mass/surface (kg/m**2)',i,a_h(i)
  write(* ,*)'Step# and normalized mass/surface (kg/m**2)',i,a_h(i)
                                                        ENDIF
!
!
			  ENDDO 
!
!
        If(CL <= 21 .or. CL == 50) HEIGHT   = JUNK
	If(CL == 30)                        MASS   = JUNK 
	If(CL == 40)                        PRESS  = JUNK         
!     
!
!     
        ENDIF       ! Endif On lh = 8  <<<<<<<<<
!
!	
!     
!       
!
!
!
!
!
!
!
          WRITE(99,*) 'Done with Load_Geometry and Load_History'
IF(IV==1) WRITE(* ,*) 'Done with Load_Geometry and Load_History'
!
          WRITE(99,*) 'Looking now for KWs Local_Study or Global_Study'
IF(IV==1) WRITE(* ,*) 'Looking now for KWs Local_Study or Global_Study'

!     
!     
!     
!     
!<><><><><><>><><><><><><><><><><><><><><><><><>
        ENDIF     ! Endif su 'load_history' ***
!<><><><><><>><><><><><><><><><><><><><><><><><>
!
!
!
!
!

!
!
!
!
!
!'''''''''''''''''''''''''''''''''''''''''''''
   IF(KEYWORD(1:11)== 'Local_Study') THEN  
!'''''''''''''''''''''''''''''''''''''''''''''
!
   n_inte = n_inte+1
!
                     call m_2(20450) ! Load_history must be active
                     call m_2(20460) ! Local_Study must not appear more than once

!
                WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword

! 
!
!
    READ(1,*) type_inten   
!
!
                     call m_2(20470) ! Cheking bounds for type_inten
!
!
 READ(1,*) IR, it, IG, dIR, dit, dIG  
!
!
                     call m_2(20480) ! One of the switches must be == 1
!
!
              write(99,*) 'The following variables will be computed and printed:'
     if(iv==1)write(* ,*) 'The following variables will be computed and printed:'
! 
if( IR ==1) write(99,*) '-> Radial displacement'
if( it ==1) write(99,*) '-> Tangential displacement'
if( IG ==1) write(99,*) '-> Geoid Heigth'
if( dIR==1) write(99,*) '-> Radial velocity'
if( dit==1) write(99,*) '-> Tangential velocity'
if( dIG==1) write(99,*) '-> Rate of change of the geoid heigth'
!
if(iv==1 .and. IR ==1) write(*,*)  '-> Radial displacement'
if(iv==1 .and. it ==1) write(*,*)  '-> Tangential displacement'
if(iv==1 .and. IG ==1) write(*,*)  '-> Geoid Heigth'
if(iv==1 .and. dIR==1) write(*,*)  '-> Radial velocity'
if(iv==1 .and. dit==1) write(*,*)  '-> Tangential velocity'
if(iv==1 .and. dIG==1) write(*,*)  '-> Rate of change of the geoid heigth'
!
!
!
!
!
  IF(type_inten==1) THEN
!                     
!
!
            Write(99,*) 'Local Study of Type#1' 
 If (iv==1) Write(* ,*) 'Local Study of Type#1' 
!
!
!
     READ(1,*) LONG_OBS(1), COLA_OBS(1)  
!
!
  Write(99,*) 'Longitude  of the observer (deg) =', LONG_OBS(1)
  Write(99,*) 'Colatitude of the observer (deg) =', COLA_OBS(1)     
 If (iv==1) then 
  Write(* ,*) 'Longitude  of the observer (deg) =', LONG_OBS(1)
  Write(* ,*) 'Colatitude of the observer (deg) =', COLA_OBS(1)  
 Endif     
!     
!                   
     READ(1,*) TMIN, TMAX, DELTAT  
!
!
                   call m_2(20490) ! tmin must be <= tmax
                   call m_2(20500) ! deltat must be > 0.0

!
     READ(1,*) ISUB
!
!
!
 If (ISUB==0)             write(99,*) 'The Sbr. AXIS_DISP0 is employed (see the User Guide)' 
 If (ISUB==1)             write(99,*) 'The Sbr. ESSA_IEZI0 is employed (see the User Guide)'  
 If (iv==1 .and. ISUB==0) write(* ,*) 'The Sbr. AXIS_DISP0 is employed (see the User Guide)'         
 If (iv==1 .and. ISUB==1) write(* ,*) 'The Sbr. ESSA_IEZI0 is employed (see the User Guide)' 
!
!
     NOBS=1
!
!
!
     ENDIF    ! Endif on Local Study of Type#1	     					               
!
!
!
!
!
!
!
   IF(type_inten==2) THEN            
!
!
!
            Write(99,*) 'Local Study of Type#2' 
 If (iv==1) Write(* ,*) 'Local Study of Type#2' 
!
!
!     
!      			                   
     READ(1,'(A20)') file_sparsi
     OPEN(2,file=file_sparsi,status='old') 
!
            Write(99,*) 'Reading the coordinates of the observers from file ',file_sparsi  
 If (iv==1) Write(* ,*) 'Reading the coordinates of the observers from file ',file_sparsi 
!
     K=0 
     DO J=1, NOBS_MAX  
       READ(2,*,END=404) LL, CC 
       K=K+1
       LONG_OBS(K) = LL
       COLA_OBS(K) = CC 
     ENDDO
!
 404 NOBS=K
     CLOSE(2)
!
       call m_2(20510) ! Too many observers
       call m_2(20520) ! Error on the user-supplied file

            Write(99,*) 'Read ', NOBS, ' observers from file ', file_sparsi  
 If (iv==1) Write(* ,*) 'Read ', NOBS, ' observers from file ', file_sparsi  
!
!
     READ(1,*) TMIN, TMAX, DELTAT  
!
               call m_2(20530) ! tmin must be <= tmax
               call m_2(20540) ! deltat must be > 0
!
!
    READ(1,*) ISUB
!
!
!
 If (ISUB==0)             write(99,*) 'The Sbr. AXIS_DISP0 is employed (see the User Guide)' 
 If (ISUB==1)             write(99,*) 'The Sbr. ESSA_IEZI0 is employed (see the User Guide)'  
 If (iv==1 .and. ISUB==0) write(* ,*) 'The Sbr. AXIS_DISP0 is employed (see the User Guide)'         
 If (iv==1 .and. ISUB==1) write(* ,*) 'The Sbr. ESSA_IEZI0 is employed (see the User Guide)' 
!

     ENDIF    ! Endif on Local Study of Type#2	 
!
!
!
!
!
!
!
  IF(type_inten==3) THEN   
!
!
            Write(99,*) 'Local Study of Type#3' 
 If (iv==0) Write(99,*) 'Local Study of Type#3' 
 If (iv==1) Write(* ,*) 'Local Study of Type#3' 
!
!                    
!
       READ(1,*)  LONG_OBS(1), C1, C2, EPC
!       
!
       call m_2(20550) ! Bounds on Long_Obs
       call m_2(20560) ! Bounds on C1, C2, & EPC
!
!
       READ(1,*)  GIVEN_TIME 
!
!
       long_obs(:) = long_obs(1) 
!
       K=1
        
       CC = C1 
       COLA_OBS(K) = CC
 
 388   CC = CC + EPC         
             
       K=K+1
!       
       if(CC > C2) then 
        COLA_OBS(K) = C2	  	  	  
        goto 401	   
       else  		   
        COLA_OBS(K) = CC  	  
	goto 388  
       endif 
!       	  
 401   nobs=k       	  
!
!
       call m_2(20570) ! Too many observers
!
!
            Write(99,*) 'There are', NOBS, ' observers along the selected meridian ' 
 If (iv==1) Write(* ,*) 'There are', NOBS, ' observers along the selected meridian ' 
!
!
    READ(1,*) ISUB
!
!
!
 If (ISUB==0)             write(99,*) 'The Sbr. AXIS_DISP0 is employed (see the User Guide)' 
 If (ISUB==1)             write(99,*) 'The Sbr. ESSA_IEZI0 is employed (see the User Guide)'  
 If (iv==1 .and. ISUB==0) write(* ,*) 'The Sbr. AXIS_DISP0 is employed (see the User Guide)'         
 If (iv==1 .and. ISUB==1) write(* ,*) 'The Sbr. ESSA_IEZI0 is employed (see the User Guide)' 
!
!
!
     ENDIF    ! Endif on Local Study of Type#3	    	
!
!
!
!
!
!   
  IF(type_inten==4) THEN   
!
!
!
            Write(99,*) 'Local Study of Type#4' 
 If (iv==1) Write(* ,*) 'Local Study of Type#4' 
!	             
!
!
       READ(1,*)  COLA_OBS(1), L1, L2, EPL
!        
!
       call m_2(20580) ! Bounds on Cola_Obs
       call m_2(20590) ! Bounds on L1, L2, & EPL
!
!
       READ(1,*)  GIVEN_TIME        
!
!
!
       cola_obs(:) = cola_obs(1) 
!
       K=1
        
       LL = L1 
       LONG_OBS(K) = LL
 
 381   LL = LL + EPL       
             
       K=K+1
!       
       if(LL > L2) then 
        LONG_OBS(K) = L2	  	  	  
        goto 409	   
       else  		   
        LONG_OBS(K) = LL  	  
	goto 381  
       endif 
!       	  
 409   nobs=k      
!
       call m_2(20570) ! Too many observers
!
!
            Write(99,*) 'There are', NOBS, ' observers along the selected parallel ' 
 If (iv==1) Write(* ,*) 'There are', NOBS, ' observers along the selected parallel ' 
!
!
!
    READ(1,*) ISUB
!
!
!
 If (ISUB==0)             write(99,*) 'The Sbr. AXIS_DISP0 is employed (see the User Guide)' 
 If (ISUB==1)             write(99,*) 'The Sbr. ESSA_IEZI0 is employed (see the User Guide)'  
 If (iv==1 .and. ISUB==0) write(* ,*) 'The Sbr. AXIS_DISP0 is employed (see the User Guide)'         
 If (iv==1 .and. ISUB==1) write(* ,*) 'The Sbr. ESSA_IEZI0 is employed (see the User Guide)' 
!
!
!
!
     ENDIF    ! Endif on Local Study of Type#4	    	
!
!  
!

! 
!
!
  IF(type_inten==5) THEN        
!
!
!
            Write(99,*) 'Local Study of Type#5' 
 If (iv==1) Write(* ,*) 'Local Study of Type#5' 
!
!  			            
!
          READ(1,*)  L1, L2, EPL
!
             call m_2(20600) ! Bounds on L1, L2, EPL
!

	  READ(1,*)  C1, C2, EPC
!
             call m_2(20610) ! Bounds on C1, C2, EPC
!
!	  
          READ(1,*)  GIVEN_TIME        
!
!
!
          K=0 
!     
          npul = int ((l2-l1)/epl) + 1; npuc = int ((c2-c1)/epc) + 1
!
          DO I = 1, npuc 
           DO J	 = 1, npul 
!
             CC = C1 + (I-1)*EPC; LL = L1 + (J-1)*EPL 
!
             if( CC>C2 .or. LL>L2 ) then 
	      goto 405
	     else 
	      k=k+1   	     		             
              LONG_OBS(K)=LL;  COLA_OBS(K)= CC
             ENDIF  
!
           ENDDO
          ENDDO 
!
 405 NOBS = K
!
          call m_2(20570) ! Too many observers
!
!
            Write(99,*) 'There are', NOBS, ' points on the selected map ' 
 If (iv==1) Write(* ,*) 'There are', NOBS, ' points on the selected map ' 
!
!
!
    READ(1,*) ISUB
!
!
!
 If (ISUB==0)             write(99,*) 'The Sbr. AXIS_DISP0 is employed (see the User Guide)' 
 If (ISUB==1)             write(99,*) 'The Sbr. ESSA_IEZI0 is employed (see the User Guide)'  
 If (iv==1 .and. ISUB==0) write(* ,*) 'The Sbr. AXIS_DISP0 is employed (see the User Guide)'         
 If (iv==1 .and. ISUB==1) write(* ,*) 'The Sbr. ESSA_IEZI0 is employed (see the User Guide)' 
!
!
!
!
!
     ENDIF    ! Endif on Local Study of Type#5	
!
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     ENDIF    ! Endif on KW 'Local_Study' 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!////////////////////////////////////////////
   IF(KEYWORD(1:12)=='Global_Study') THEN
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!
!
   n_este = n_este + 1
!
      call m_2(20620) ! Load_History must be active
      call m_2(20630) ! Global_Study must appear once
      call m_2(20640) ! Global_Study & Local_Study are mutually exclusive
!
                WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword

! 
!
!
     IG  =1   
    dIG  =1    
!
!
!
    READ(1,*) type_esten  
!
!
              call m_2(20650) ! The range for type_esten is [0:3]
!
!
  IF(type_esten==1) THEN        
!
!             
!
            Write(99,*) 'Global Study of Type#1' 
 If (iv==1) Write(* ,*) 'Global Study of Type#1' 
!
!  
!
   READ(1,*) LDE, MDE
!
!
              call m_2(20660) ! check LDE & MDE
!       

!
!      
   Read(1,*) Inorm      
!
!
              call m_2(20670) ! check for Inorm
!

!               
   READ(1,*) TMIN, TMAX, DELTAT
!
!
              call m_2(20680) ! A check for TMIN, TMAX & DELTAT
!
!
  ENDIF      ! Endif on Global Study of Type#1
!
!
!
!
!
!
  IF(type_esten==2) THEN        
!             
!
!
!
!
!
            Write(99,*) 'Global Study of Type#2' 
 If (iv==1) Write(* ,*) 'Global Study of Type#2' 
!
!  
!
!
   READ(1,*) LDE1, LDE2 
!
!
             call m_2(20690) ! check for LDE1 & LDE2
!
!      
   Read(1,*) Inorm      
!
!
             call m_2(20670) ! check for Inorm
!
!
   READ(1,*) GIVEN_TIME
!
!    
   ENDIF   ! Endif on Global Study of Type#2   
!
!
!
!
  IF(type_esten==3)  Then 
!
!  
!
            Write(99,*) 'Global Study of Type#3' 
 If (iv==1) Write(* ,*) 'Global Study of Type#3' 
!
!  
!  
  READ(1,*) TMIN, TMAX, DELTAT
!
!
            call m_2(20680) ! A check for TMIN, TMAX & DELTAT
!
            call m_2(20700) ! check for l_min & l_max in Harmonic_Degrees
!
!
  ENDIF    ! Endif on Global Study of Type#3                      
!
!
!
!
!
!**************************************
  ENDIF  ! Endif on  'Global_Study' 
!**************************************
!
!
!
!
!
10001 CONTINUE    ! Loop sulla lettura di task_2.dat
!
10002 CONTINUE 
!
!
! 
          Write(99,*) 'task_2.dat is being closed... done '
if(iv==1) Write(* ,*) 'task_2.dat is being closed... done '
!
!
!
 CLOSE (1)                      
!
!
!
!
!
!
!
!
 write(99,*) 'The KWs Local_Study (Global_Study) are been configured correctly' 
if(iv==1) & 
 write(* ,*) 'The KWs Local_Study (Global_Study) are been configured correctly' 
!
!
 write(99,*) 'Writing the harmonic components of the load on file load_coeff.dat' 
if(iv==1) & 
 write(* ,*) 'Writing the harmonic components of the load on file load_coeff.dat' 
!
!
!
!
!
!
!
!
If (CL <= 21) Then
 CALL AXIS_LOAD (1, CL, amplitude, height, mass)
!
          Write(99,*) 'The Maximum Load MASS is ', MASS, 'kg'
If(iv==1) Write(* ,*) 'The Maximum Load MASS is ', MASS, 'kg'
!
 Endif  
!
!
!
!    
If (CL == 30) Then
 CALL AXIS_LOAD (1, CL, junk, junk, mass)
          Write(99,*) 'The Maximum Load MASS is ', MASS, 'kg'
If(iv==1) Write(* ,*) 'The Maximum Load MASS is ', MASS, 'kg'
Endif  
!  
!
!
!
If (CL == 40) Then
 CALL AXIS_LOAD (1, CL, junk, junk, press)
          Write(99,*) 'The Maximum Load MASS/surface is ', PRESS, 'kg/m/m'
If(iv==1) Write(* ,*) 'The Maximum Load MASS/surface is ', PRESS, 'kg/m/m'
Endif  
! 
!
!
!   
If (CL == 50) Then
 CALL RECT_LOAD0 (1, Long_cen, Teta_cen, dl, dt, height, mass)
          Write(99,*) 'The Maximum Load MASS is ', MASS, 'kg'
If(iv==1) Write(* ,*) 'The Maximum Load MASS is ', MASS, 'kg'
Endif  
!         
!
!
!
!
!
!
!
! ooooooooooooooooooooooooooooooooooooooooooooooooooo
! o                                                 o
! o     The input data  have been assimilated       o
! o                                                 o
! o     _Now I procede with the computations_       o
! o                                                 o
! ooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!
!
!
!
!
!
!
!----------------------------------------------
  IF( type_inten==1 .OR. type_inten==2)   THEN
!----------------------------------------------
!
!
!
 if(type_inten==1) then 
            Write(99,*) 'Computing for Local Study of Type#1' 
 If (iv==1) Write(* ,*) 'Computing for Local Study of Type#1' 
 endif
 if(type_inten==2) then 
            Write(99,*) 'Computing for Local Study of Type#2' 
 If (iv==1) Write(* ,*) 'Computing for Local Study of Type#2' 
 endif
! 
!
!
!
 IF(IR==1 .or. it==1  .or.  IG==1) Open(10,File='disp.his',Status='Unknown')
 IF(dIR==1.or.dit==1  .or. dIG==1) Open(11,File='rate.his',Status='Unknown')
!
!
!
      IF (IOC==1) then 
                  Call READ_OCEAN      
	                    Write(99,*) 'Read the ocean function from oceano.128'
		  if(iv==1) Write(* ,*) 'Read the ocean function from oceano.128'
      endif
!      
!
      CEN(1) = LONG_CEN
      CEN(2) = TETA_CEN 
!
      Do 222 K = 1, NOBS    
!
        OBS(1) = LONG_OBS(K)
        OBS(2) = COLA_OBS(K)
!
        call printdue (101)   ! Puts an header on disp.his & rate.his    
!

	Do 111 TIME = TMIN, TMAX, DELTAT 
!
	  Call CONVOL_load_history (lh, time)
!
          MICE = 0.d0
          If(CL<=30 .or. CL==50) MICE =  MASS  * load_history (lh, time)
          If(            CL==40) MICE = PRESS  * load_history (lh, time)
!
!
          If(CL <= 40 .and. Isub == 0) & 
	                  call axis_disp0 ( IOC, OBS, CEN, VECT, D_VECT ) 
!			  
          If(CL <= 40 .and. Isub == 1) & 
	                  call essa_iezi0 ( IOC, OBS, CEN, VECT, D_VECT ) 
!			  			  
          If(CL==50) CAll RECT_DISP0  ( IOC, OBS, VECT, D_VECT) 
!
!
!
   IF(IR ==1 .or. IT==1  .or.  IG ==1) &
   Write(10,'(5(1x,f9.4),1x,E14.5)') &
                 TIME,   vect(1),   vect(2),  &
                                   vect(3),   vect(4),    Mice
!
   IF(DIR==1 .or. DIT==1 .or. DIG ==1) &
   Write(11,'(5(1x,f9.4),1x,E14.5)') &
                 TIME, d_vect(1), d_vect(2),  &
                                 d_vect(3), d_vect(4),    Mice
!
!
!
111     Continue   ! Time do-loop
!
!
222     Continue   ! Observers do-loop
!
!
!
        IF(IV==1) THEN
!
 IF(IR==1 .or. it==1  .or.  IG==1) &
        WRITE(*,*) '-> Output written on file disp.his'

 IF(dIR==1.or.dit==1.or.dIG==1) &
        WRITE(*,*) '-> Output written on file rate.his'
!
        ENDIF
!
!
!
 IF(IR==1 .or. it==1  .or.  IG==1) &
        WRITE(99,*) '-> Output written on file disp.his'

 IF(dIR==1.or.dit==1.or.dIG==1) &
        WRITE(99,*) '-> Output written on file rate.his'
!
!
        Close(10)
        Close(11) 
!
!----------------------------------------------------
        ENDIF    ! Endif on Local_Study = 1 or 2
!----------------------------------------------------
!
!
!
!
!
!
!
!
!
!..............................................
  IF( type_inten==3 .OR. type_inten==4)   THEN 
!..............................................
!
 if(type_inten==3) then 
            Write(99,*) 'Computing for Local Study of Type#3' 
 If (iv==1) Write(* ,*) 'Computing for Local Study of Type#3' 
 endif
 if(type_inten==4) then 
            Write(99,*) 'Computing for Local Study of Type#4' 
 If (iv==1) Write(* ,*) 'Computing for Local Study of Type#4' 
 endif
! 
!
 IF(IR==1 .or. it==1  .or.  IG==1) &
        Open(10,File='disp.pro',Status='Unknown')

 IF(dIR==1.or.dit==1.or.dIG==1) &
        Open(11,File='rate.pro',Status='Unknown')
!
!
      IF (IOC==1) then 
                  Call READ_OCEAN      
	                    Write(99,*) 'Read the ocean function from oceano.128'
		  if(iv==1) Write(* ,*) 'Read the ocean function from oceano.128'
      endif
!      
     CEN(1) = LONG_CEN
     CEN(2) = TETA_CEN 
!
     TIME = GIVEN_TIME
!
	Call CONVOL_load_history (lh, TIME)
!
        MICE = 0.0
        If(CL<=30 .or. CL==50) MICE =  MASS  * load_history (lh, time)
        If(            CL==40) MICE = PRESS  * load_history (lh, time)
!
     Do 223 K = 1, NOBS    
!
        OBS(1) = LONG_OBS(K)
        OBS(2) = COLA_OBS(K)
!
!
          If(CL <= 40 .and. Isub == 0) & 
	                  call axis_disp0 ( IOC, OBS, CEN, VECT, D_VECT ) 
!			  
          If(CL <= 40 .and. Isub == 1) & 
	                  call essa_iezi0 ( IOC, OBS, CEN, VECT, D_VECT ) 
!			  			  
          If(CL==50) CAll RECT_DISP0  ( IOC, OBS, VECT, D_VECT) 
!
!
          call printdue (105)  ! Writes on the ".pro" files
!
!
 If(type_inten==3) then
  IF(IR==1 .or. it==1  .or.  IG==1) &
    Write(10,'(6(1x,f9.4))') obs(2),   vect(1),   vect(2),   vect(3),   vect(4)
  IF(dIR==1.or.dit==1.or.dIG==1) &
    Write(11,'(6(1x,f9.4))') obs(2), d_vect(1), d_vect(2), d_vect(3), d_vect(4)
 Endif     
! 
 If(type_inten==4) then
   IF(IR==1 .or. it==1  .or.  IG==1) &
     Write(10,'(6(1x,f9.4))') obs(1),   vect(1),   vect(2),   vect(3),   vect(4)
   IF(dIR==1.or.dit==1.or.dIG==1) &
     Write(11,'(6(1x,f9.4))') obs(1), d_vect(1), d_vect(2), d_vect(3), d_vect(4)
 Endif     
!
! 
!
223     Continue   ! Observers do-loop
!
!
!
 IF(IV==1) THEN
!
 IF(IR==1 .or. it==1  .or.  IG==1) &
        WRITE(*,*) '-> Output written on file disp.pro'

 IF(dIR==1.or.dit==1.or.dIG==1) &
        WRITE(*,*) '-> Output written on file rate.pro'
!
 ENDIF
!
 IF(IV==0) THEN
!
 IF(IR==1 .or. it==1  .or.  IG==1) &
        WRITE(99,*) '-> Output written on file disp.pro'

 IF(dIR==1.or.dit==1.or.dIG==1) &
        WRITE(99,*) '-> Output written on file rate.pro'
!
 ENDIF
!
        Close(10)
	Close(11) 
!
!(((((((((((((((((((((((((((((((((((((((((((((((((((
        ENDIF    ! Endif on Local_Study = 3 or 4 
!)))))))))))))))))))))))))))))))))))))))))))))))))))       
! 
!    
!    
!    
!
!
!
!:::::::::::::::::::::::::::
  IF(type_inten==5)  THEN 
!::::::::::::::::::::::::::: 
!
!
            Write(99,*) 'Computing for Local Study of Type#5' 
 If (iv==1) Write(* ,*) 'Computing for Local Study of Type#5' 
!
!
!
If(IR==1)   Open(48,file= 'A_rad.dat', status='unknown')  ! Vertical
If(it==1)   Open(49,file= 'A_lon.dat', status='unknown')  ! Along Longitude
If(it==1)   Open(50,file= 'A_lat.dat', status='unknown')  !   "   LATITUDE
If(IG==1)   Open(51,file= 'A_geo.dat', status='unknown')  ! Geoid
!
!
If(dIR==1)  Open(52,file='Ad_rad.dat', status='unknown')  ! Vertical
If(dit==1)  Open(53,file='Ad_lon.dat', status='unknown')  ! Along Longitude
If(dit==1)  Open(54,file='Ad_lat.dat', status='unknown')  !   "   LATITUDE
If(dIG==1)  Open(55,file='Ad_geo.dat', status='unknown')  ! Geoid
!
!
      IF (IOC==1) then 
                  Call READ_OCEAN      
	                    Write(99,*) 'Read the ocean function from oceano.128'
		  if(iv==1) Write(* ,*) 'Read the ocean function from oceano.128'
      endif
!      
     CEN(1) = LONG_CEN
     CEN(2) = TETA_CEN 
!
     TIME = GIVEN_TIME
!
!
	Call CONVOL_load_history (lh, TIME)
!
        MICE = 0.0
        If(CL<=30 .or. CL==50) MICE =  MASS  * load_history (lh, time)
        If(            CL==40) MICE = PRESS  * load_history (lh, time)
!
!	       
     Do 224 K = 1, NOBS    
!
        OBS(1) = LONG_OBS(K)
        OBS(2) = COLA_OBS(K)     
!    
!
!
          If(CL <= 40 .and. Isub == 0) & 
	                  call axis_disp0 ( IOC, OBS, CEN, VECT, D_VECT ) 
!			  
          If(CL <= 40 .and. Isub == 1) & 
	                  call essa_iezi0 ( IOC, OBS, CEN, VECT, D_VECT ) 
!			  			  
          If(CL==50)                   & 
	                  CAll RECT_DISP0  ( IOC, OBS, VECT, D_VECT) 
!
!
          call printdue(106)  ! Writes on the A*.dat files
!
!
224     Continue   ! Observers do-loop
!
!
!
If(IR==1) CLOSE(48);  If(dIR==1) CLOSE(52)
If(it==1) CLOSE(49);  If(dit==1) CLOSE(53)
If(it==1) CLOSE(50);  If(dit==1) CLOSE(54)
If(IG==1) CLOSE(51);  If(dIG==1) CLOSE(55)
!
!
 IF(IV==1)  THEN
            Write (*,*) 'Output files:'  
 If (IR==1) WRITE (*,*) ' > A_rad.dat'
 If (it==1) WRITE (*,*) ' > A_lon.dat, A_lat.dat'
 If (IG==1) WRITE (*,*) ' > A_geo.dat'
 If(dIR==1) WRITE (*,*) ' > A_rad.dat'
 If(dit==1) WRITE (*,*) ' > A_lon.dat, A_lat.dat'
 If(dIG==1) WRITE (*,*) ' > A_geo.dat'
 ENDIF
!
!
 IF (IV==0) THEN
            Write(99,*) 'Output files:'   
 If (IR==1) WRITE(99,*) ' > A_rad.dat'
 If (it==1) WRITE(99,*) ' > A_lon.dat, A_lat.dat'
 If (IG==1) WRITE(99,*) ' > A_geo.dat'
 If(dIR==1) WRITE(99,*) ' > A_rad.dat'
 If(dit==1) WRITE(99,*) ' > A_lon.dat, A_lat.dat'
 If(dIG==1) WRITE(99,*) ' > A_geo.dat'
 ENDIF
!
!
!::::::::::::::::::::::::::::::::::::::::::::::::::
        ENDIF     ! Endif on Local_Study = 5
!::::::::::::::::::::::::::::::::::::::::::::::::::	        
!
!
!
!
!
! 
!++++++++++++++++++++++++++++
  IF(type_esten == 1)  THEN 
!++++++++++++++++++++++++++++
!
!
            Write(99,*) 'Computing for Global Study of Type#1' 
 If (iv==1) Write(* ,*) 'Computing for Global Study of Type#1' 
!
!
  open(10,file='stokes.his',status='unknown') 
  open(11,file='stoked.his',status='unknown') 
!  
!

!
  Do IJUNK=10, 11
     call DATE_AND_TIME (date,timc)
     Write(IJUNK,*) '# done by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
     Write(IJUNK,*) '# Degree   =', LDE 
     Write(IJUNK,*) '# Order    =', MDE 
  End Do 
!
!
  If(Inorm == 1) Then 
  Do IJUNK=10, 11
     Write(IJUNK,*) '# The coefficients are Fully Normalized'  
  End Do   
      Else 
  Do IJUNK=10, 11
     Write(IJUNK,*) '# The coefficients are NOT Fully Normalized'  
  End Do   
      Endif 
!
!
  Write(11,*) '#      Divide dc_lm and ds_lm by 10**11             '
  Write(11,*) '#    to obtain values in units of yr**(-1)          '
  Write(11,*) '#  (time is in kyrs, the load mass is in kg)        '
  WRITE(11,*) '# '
!
  Write(10,*) '# Divide c_lm and s_lm by 10**6 to obtain actual values '  
  Write(10,*) '#      (time is in kyrs, the load mass is in kg)        '
  WRITE(10,*) '# '
!  
  Write(10,*) '#    time      c_lm          s_lm         mass   '
  Write(11,*) '#    time     dc_lm         ds_lm         mass   '
!
!
      IF (IOC==1) then 
                  Call READ_OCEAN      
	                    Write(99,*) 'Read the ocean function from oceano.128'
		  if(iv==1) Write(* ,*) 'Read the ocean function from oceano.128'
      endif
!      
     CEN(1) = LONG_CEN
     CEN(2) = TETA_CEN 
!
!
   Do 1111 TIME = TMIN, TMAX, DELTAT 
!
        MICE=0.d0 
!
	Call CONVOL_load_history (lh, time)
!
        If(CL<=30 .or. CL==50) MICE =  MASS  * load_history (lh, time)
        If(            CL==40) MICE = PRESS  * load_history (lh, time)
!
          If(CL<=40) Call AXIS_STOK (LDE, MDE, ioc, inorm, CEN, W, D_W)  
          If(CL==50) CAll RECT_STOK (LDE, MDE, ioc, inorm,      W, D_W) 
!
   Fact1= 1D6
   Fact2= 1D11
!
   Write(10,'(f9.4,3(e14.5))')    time,   FACT1*  w(1), &
                                            FACT1*  w(2), &
                                           MICE
   Write(11,'(f9.4,3(e14.5))')    time,   FACT2*d_w(1)/1000D0, &
                                            FACT2*d_w(2)/1000D0, &
                                           MICE
!
!
 1111   CONTINUE   ! Time do-loop
!
  Close(10)
  Close(11)
!
 IF(IV==1) then
   WRITE(* ,*) '-> Output written on file stokes.his'
   WRITE(* ,*) '-> Output written on file stoked.his'
 ENDIF
!
 IF(IV==0) then
   WRITE(99,*) '-> Output written on file stokes.his'
   WRITE(99,*) '-> Output written on file stoked.his'
 endif 
!
!
!+++++++++++++++++++++++++++++++++++++++++
  ENDIF    ! Endif on Global_Study = 1
!+++++++++++++++++++++++++++++++++++++++++    
!
!
!
!
!
!
!=============================
   IF(type_esten == 2)  THEN
!=============================
!
!
            Write(99,*) 'Computing for Global Study of Type#2' 
 If (iv==1) Write(* ,*) 'Computing for Global Study of Type#2' 
!
!
  open(10,file='stokes.pro',status='unknown') 
  open(11,file='stoked.pro',status='unknown') 
!
!
      IF (IOC==1) then 
                  Call READ_OCEAN      
	                    Write(99,*) 'Read the ocean function from oceano.128'
		  if(iv==1) Write(* ,*) 'Read the ocean function from oceano.128'
      endif
!
!
           CEN(1) = LONG_CEN
           CEN(2) = TETA_CEN 
!
        TIME = GIVEN_TIME
!
!
!
     DO IJUNK=10, 11 
     call DATE_AND_TIME (date,timc)
     Write(IJUNK,*) '# done by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
     Write(ijunk, *) '# time (kyr) = ', TIME
     ENDDO
!
!
        MICE=0.D0
        If(CL<=30 .or. CL==50) MICE =  MASS  * load_history (lh, time)
        If(            CL==40) MICE = PRESS  * load_history (lh, time)
!
!
  If(Inorm == 1) Then 
  Do IJUNK=10, 11
     Write(IJUNK,*) '# The coefficients are Fully Normalized'  
  End Do   
      Else 
  Do IJUNK=10, 11
     Write(IJUNK,*) '# The coefficients are NOT Fully Normalized'  
  End Do
      Endif 
!
 WRITE(11,*) '# *************************************************** '
 Write(11,*) '#         Divide dc_lm and ds_lm by 10**11            '
 Write(11,*) '#   to obtain actual values in units of yr**(-1)      '
!
 WRITE(10,*) '# ***************************************************** '
 Write(10,*) '# Divide c_lm and s_lm by 10**6 to obtain actual values '       
! 
!  
 Write(10,*) '#      l   m       c_lm         s_lm            '   
 Write(11,*) '#      l   m       dc_lm        ds_lm           '  
! 
         Call CONVOL_load_history (lh, time )
!
           Do 1113  L = LDE1, LDE2 
!
            Do 1114  M = 0, L
!
            If(CL<=40) Call AXIS_STOK (L, M, ioc, inorm, CEN, W, D_W)  
            If(CL==50) CAll RECT_STOK (L, M, ioc, inorm,      W, D_W) 	   
!
!
   Fact1= 1D6
   Fact2= 1D11
!   
 write (10, '(i4,1x,2(1x,i3),2(e14.5))') &
        indx(l,m), l,  m,  FACT1*w(1),   FACT1*w(2)
 write (11, '(i4,1x,2(1x,i3),2(e14.5))') &
        indx(l,m), l,  m,  FACT2*d_w(1), FACT2*d_w(2)
!
!    
  1114      CONTINUE
 1113      CONTINUE
!
!
 IF(IV==1) then
   WRITE(* ,*) '-> Output written on file stokes.pro'
   WRITE(* ,*) '-> Output written on file stoked.pro'
 ENDIF
!
 IF(IV==0) then
   WRITE(99,*) '-> Output written on file stokes.pro'
   WRITE(99,*) '-> Output written on file stoked.pro'
 ENDIF
        CLOSE(10)
	CLOSE(11) 	
!
!===========================================
  ENDIF     ! Endif on Global_Study = 2
!===========================================
!
!
!
!
!+++++++++++++++++++++++++++
  IF(type_esten == 3)  THEN
!+++++++++++++++++++++++++++
!
            Write(99,*) 'Computing for Global Study of Type#3' 
 If (iv==1) Write(* ,*) 'Computing for Global Study of Type#3' 
!
!
 if(cl <= 30 .or. cl ==50) Open( 12, File ='load_mass.dat ', status='unknown') 
 if(              cl ==40) Open( 12, File ='load_pole.dat ', status='unknown') 
!
  open(10,file='iner.his',status='unknown') 
  open(11,file='ined.his',status='unknown') 
!
! 
     DO IJUNK=10, 12, 1  	
     call DATE_AND_TIME (date,timc)
     Write(IJUNK,*) '# done by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
     ENDDO
!
!
 if(cl <= 30 .or. cl ==50) then
             Write(12,*) '##################################################'
             Write(12,*) '# Mass of the primary load as a function of time #'
             Write(12,*) '##################################################'
             Write(12,*) '#time (kyrs)   primary load mass (kg)             '
 Endif 
!
 if(              cl ==40) then 
       Write(12,*) '#######################################################'
       Write(12,*) '# Mass/surface at the load pole as a function of time #'
       Write(12,*) '#######################################################'
       Write(12,*) '#time (kyrs)    mass/surface (kg/m/m)                  '  
 Endif	     
!    
!
 Write(10,*) '# - - - - - - - - - - - - - - - - - - - - - - - -    '                                                     
 Write(10,*) '# The inertia tensor is in a non-dimensional form.   ' 
 Write(10,*) '# To obtain the actual values in units of kg*m**2    ' 
 Write(10,*) '# multiply the following components by MR**2, with   '      
 Write(10,*) '# M = Mass of the Earth, and with R = Earth radius   ' 
 Write(10,*) '#       [n.b.: M R**2 = 2.42 x 10**32 kg*m**2]       ' 
 Write(10,*) '#                                                    '
 Write(10,*) '#     time (column #1) is given in units of kyrs     '
 Write(10,*) '# - - - - - - - - - - - - - - - - - - - - - - - - -  '
! 
!
 Write(11,*) '# = = = = = = = = = = = = = = = = = = = = = = = = = = = =' 
 Write(11,*) '# To obtain the actual derivative of the  inertia tensor '
 Write(11,*) '# in kg*m**2/yr, you must multiply the components  below '       
 Write(11,*) '# by MR**2, with M= Mass of the Earth, and R= its Radius '
 Write(11,*) '#        [n.b.: M R**2 = 2.42 x 10**32 kg*m**2]          ' 
 Write(11,*) '#                                                        '
 Write(11,*) '#      time (column #1) is given in units of kyrs        '
 Write(11,*) '# = = = = = = = = = = = = = = = = = = = = = = = = = = = =' 
!  
 Write(10,*) '# time     (xz)         (yz)         (zz)         (xy)         (yy)         (xx)'
 Write(10,*) '# ' 						      	                     
 Write(11,*) '# time     (xz)         (yz)         (zz)         (xy)         (yy)         (xx)'
 Write(11,*) '# ' 				                                            
! 
!
      IF (IOC==1) then 
                  Call READ_OCEAN      
	                    Write(99,*) 'Read the ocean function from oceano.128'
		  if(iv==1) Write(* ,*) 'Read the ocean function from oceano.128'
      endif
!      
     CEN(1) = LONG_CEN
     CEN(2) = TETA_CEN   
!
	Do 1112 TIME = TMIN, TMAX, DELTAT  
!
        Mice=0.d0 
!
	  Call CONVOL_load_history (lh, time)
!
        If(CL<=30 .or. CL==50) MICE =  MASS  * load_history (lh, time)
        If(            CL==40) MICE = PRESS  * load_history (lh, time)
!
          If(CL<=40) Call AXIS_INER (ioc, CEN, INER, D_INER)  
          If(CL==50) CAll RECT_INER (ioc,      INER, D_INER) 
!
        fact1=1000.D0 
!   
!
        Write(10,'(f9.4,2x,7(e12.5,1x))') time, &
	    iner(1,3), iner(2,3), iner(3,3), & 
	    iner(1,2), iner(2,2), iner(1,1) 
!
        Write(11,'(f9.4,2x,7(e12.5,1x))') time,  &
            d_iner(1,3)/fact1, d_iner(2,3)/fact1, d_iner(3,3)/fact1, & 
	    d_iner(1,2)/fact1, d_iner(2,2)/fact1, d_iner(1,1)/fact1 
!
        Write(12,'(2x,f9.4,1x,E20.8)') time,   MICE
!
!
!
 1112   CONTINUE     ! Time do-loop
!
!
!
 IF(IV==1) then
   WRITE(* ,*) '-> Output written on file iner.his'
   WRITE(* ,*) '-> Output written on file ined.his'
   if(cl <= 30 .or. cl ==50) WRITE(*, *) '-> Output written on file load_mass.dat'
   if(              cl ==40) WRITE(*, *) '-> Output written on file load_pole.dat'
 ENDIF
!
 IF(IV==0) then
   WRITE(99,*) '-> Output written on file iner.his'
   WRITE(99,*) '-> Output written on file ined.his'
   if(cl <= 30 .or. cl ==50) WRITE(99,*) '-> Output written on file load_mass.dat'
   if(              cl ==40) WRITE(99,*) '-> Output written on file load_pole.dat'
 ENDIF
!
   CLOSE(10)
   CLOSE(11) 
   CLOSE(12)	
!	
!=========================================
  ENDIF    ! Endif on Global_Study = 3
!=========================================
!
!
!
!
!
!
!
!
!
!
!
!+==================================+
!                                   !
        END SUBROUTINE TASK_2       !
!                                   !
!+==================================+
!
!
!
!
!
!
!
!
  Function INDX(l,m) 
!
! +------------------------------------------------+
! | For l=degree and m=order, INDX = l*(l+1)/2+m+1 |
! +------------------------------------------------+
!
  Use STRATA; Implicit NONE    
!
  Integer(i4b) :: l, m     ! degree and order 
  Integer(i4b) :: INDX 
!
  INDX = L*(L+1)/2+M+1  
!
  END Function INDX 
!
!
!
!
!
!
!
!
!
!
!
! +--------------------------------------------------------------+					                       
   subroutine ESSA_IEZI0   (ICO, OBSE, CENT, VECT, D_VECT)       !                 
! +--------------------------------------------------------------+	
! 
!	This routine is an alternative to AXIS_DISP0.
! 
!	It computes the surface deformations and the geoid, as
! well as their time derivatives, for a given time-history
! and a glacial axissimmetric load placed anywhere. However,
! differently from AXIS_DISP0, ESSA_IEZI0 uses the formulas
! in which the axial symmetry of the load is NOT invoked to
! simplify the computations (and save CPU time).
!
!       Before calling ESSA_IEZI0, one must call sbr AXIS_LOAD,
! which prepares the harmonic coefficients of the load sigma(l),
! and sbr convol_load_history, which computes the convolutions
! between the load time-history and the load-deformation
! coefficients (ldcs) h, l, and k.
!
!       The routine ESSA_IEZI0 is much more slower than AXIS_DISP0,
! so that we recommnend the use of the latter. However, the formulas
! on which it is based are useful for computing the effect of
! the oceanic load when real continent bounaries are employed.
!
!
! # Fixed 'special case 5' (observer and load at the antipodes)
!   (DM Feb 20, 2017)
!
 use COMMON; use TRIGFCN; implicit NONE   
!
 real(dp) :: cos_delta
 real(dp), parameter :: eps = 1.D-200    ! Tolerance for antipode check
!
 real(dp) :: VECT(4), D_VECT(4) ! Solution vectors
 real(DP) :: CORR(4), D_CORR(4) ! Their time-derivatives
!!
 real(dp) :: obse(2), teta_obs, long_obs   ! Observer Long. & Colat.
 real(dp) :: cent(2), teta_cen, long_cen   ! Long. & Colat. of the load centroid.
!  
 real(dp) :: u_rad,   u_teta,  u_long,   geoid  ! Radial & Horizontal displacement,
!                                                 & Geoid.
 real(dp) :: du_rad, du_teta, du_long, d_geoid  ! Their time-derivatives.
!    
 real(dp) :: leg     ! Associated Legendre function
 REAL(dp) :: legdev1 ! Its time derivative w.r.t. colatitude
 real(dp) :: GPROD   ! FUNCTION gprod gives (l-m)!/(l+m)!
!    
 real(dp) :: CFA, COEFF, AUX, BRA1, BRA2  ! Auxiliary variables
! 
 real(SP) :: OBS_C(llmin:llmax)    ! cos(m*teta_obs)
 real(SP) :: OBS_S(llmin:llmax)    ! sin(m*teta_obs) 
 real(SP) :: CEN_C(llmin:llmax)    ! cos(m*teta_cen) 
 real(SP) :: CEN_S(llmin:llmax)    ! sin(m*teta_cen)         
!
 real(dp) :: delta(llmin:llmax)    ! Kronecker delta
!
 real(dp) :: PLM, PLMD             ! Current Legendre Functions and derivative
!
 integer (i4b) :: l, m             ! Degree & order
!
 integer (i4b) :: ICO              ! Ocean correction
!
 logical :: log_1, log_2, &        ! Auxiliary logical variables
            log_3, log_4
!
!
!
!
!
!
!
! +========================+ 
! | Executable statements  |
! +========================+
!
 teta_obs  = obse(2);      long_obs  = obse(1)  
 teta_cen  = cent(2);      long_cen  = cent(1)  
!
 VECT(:)=0.d0; D_VECT(:)=0.d0 
 CORR(:)=0.d0; D_CORR(:)=0.d0 
!
!                                   
!
  u_rad = 0._dp;   u_teta = 0._dp;   u_long = 0._dp;   geoid = 0._dp
 du_rad = 0._dp;  du_teta = 0._dp;  du_long = 0._dp; d_geoid = 0._dp 
!
!
!
! # these coefficients depend only on the order
! 
    Do m=0, lmax 
      obs_c(m) =  cosd(m*long_obs) 
      obs_s(m) =  sind(m*long_obs) 
      cen_c(m) =  cosd(m*long_cen) 
      cen_s(m) =  sind(m*long_cen) 
      If(m  ==  0) delta(m) = 1.d0 
      If(m /= 0) delta(m) = 0.d0 
    Enddo   
!
!
!
! +---------------------------------------+
! # Loop on the real spherical Harmonics  |
! +---------------------------------------+
!
!
    DO 100 L = LMIN, LMAX 
!
!
    CFA = (sigma(l)/(2d0*l + 1d0))*(3d0/rhoea) 
!
!
    DO 100 M = 0, L 
!
    PLM  =     Leg(l,m,cosd(teta_obs)) 
!
    COEFF = GPROD(l,m) * CFA * LEG (l,m,cosd(teta_cen)) 
!
!
    If(IR==1.or.dIR==1.or.IG==1.or.dIG==1) THEN 
!
!
    BRA1 = (2d0-delta(m))*COEFF*( OBS_C(m)*CEN_C(m) + OBS_S(m)*CEN_S(m)) * PLM     
!
    If(IR   ==1 )   u_rad =   u_rad +  conv_h(l) * BRA1 
    If(dIR ==1 )  du_rad =  du_rad + dconv_h(l) * BRA1
    If(IG   ==1 )   geoid =   geoid +  conv_k(l) * BRA1
    If(dIG ==1 ) d_geoid = d_geoid + dconv_k(l) * BRA1
!
                                           ENDIF
!
!
!						     
    If(it==1.or.dit==1)  THEN
!
    PLMD = Legdev1(l,m,cosd(teta_obs))     
!
    BRA1 =  (2d0-delta(m))*COEFF*(OBS_C(m)*CEN_C(m) + OBS_S(m)*CEN_S(m)) * PLMD
!
    BRA2 =   2d0* Float(m)*COEFF*(OBS_C(m)*CEN_S(m) - OBS_S(m)*CEN_C(m)) * PLM        
!    
    If(it==1   .and. sind(teta_obs) /= 0.d0) Then 
        u_teta =  u_teta + conv_l(l)* BRA1 
        u_long =  u_long + conv_l(l)* BRA2/Sind(teta_obs)
    Endif
!       						     
    If(dit==1 .and. sind(teta_obs) /= 0.d0) Then 
       du_teta = du_teta + dconv_l(l)* BRA1
       du_long = du_long + dconv_l(l)* BRA2/Sind(teta_obs)
    Endif
!
                         ENDIF
!
!
!
 100 CONTINUE    ! do-loop on the harmonic degree
!
!
! +----------------+
! | Special Cases  |
! +----------------+
!
! _Load at the pole & observer not at the pole_
!
log_1 = (teta_cen == 0._dp .or. teta_cen == 180._dp)
log_2 =  teta_obs /= 0._dp 
! 
if (log_1 .and. log_2) then 
u_long = 0._dp; du_long = 0._dp
                       Endif  
!                                                 
!
! _Load and observer at one of the poles_
!
log_1 = teta_cen == 0._dp   .AND. teta_obs ==   0._dp
log_2 = teta_cen == 0._dp   .AND. teta_obs == 180._dp
log_3 = teta_cen == 180._dp .AND. teta_obs ==   0._dp
log_4 = teta_cen == 180._dp .AND. teta_obs == 180._dp
!
if(log_1.or.log_2.or.log_3.or.log_4) then                                    
 u_teta = 0._dp ;  u_long = 0._dp 
du_teta = 0._dp ; du_long = 0._dp   
                                     endif
!
!
!
! _Load not at the pole, observer in one of the poles_
!
! The radial displacement & geoid are computed normally,
! whereas u_teta & u_long are set to zero conventionally
! (the direction of the colatitude & longitude unit vectors
! is not defined at the poles)
! 
log_1 = teta_cen /=    0._dp .and. teta_obs ==    0._dp 
log_2 = teta_cen /=  180._dp .and. teta_obs ==    0._dp 
log_3 = teta_cen /=    0._dp .and. teta_obs ==  180._dp 
log_4 = teta_cen /=  180._dp .and. teta_obs ==  180._dp 
!
if(log_1.or.log_2.or.log_3.or.log_4) then                                    
 u_teta = 0._dp ;  u_long = 0._dp  
du_teta = 0._dp ; du_long = 0._dp   
                                    endif 
!
!
!
! _Load not at the pole, and observer at its center_
!
log_1 =  teta_cen /=   0._dp .and. teta_obs == teta_cen   
log_2 =  teta_cen /= 180._dp .and. teta_obs == teta_cen 
log_3 =  long_cen == long_obs
!
if((log_1  .and. log_3) .or. (log_2 .and. log_3)) then 
 u_teta = 0._dp ;  u_long = 0._dp 
du_teta = 0._dp ; du_long = 0._dp 
                                                 endif
!
!
!
! _Load not at the pole, observer at its antipodes_
!
!log_1 = abs(long_cen - long_obs) == 180._dp 
!log_2 =     teta_cen + teta_obs  == 180._dp
!
!if( log_1 .or. log_2 ) then  
!
cos_delta = cosd(teta_obs)*cosd(teta_cen) + 			  & 
            sind(teta_obs)*sind(teta_cen)*cosd(long_obs-long_cen) 
log_1 = (cos_delta + 1.d0) < eps
!
if( log_1 ) then
  u_teta = 0._dp ;  u_long = 0._dp 
 du_teta = 0._dp ; du_long = 0._dp  
                      endif 
!
!
!
!+-------------------------------+
!  Assignements of vect & d_vect |
!+-------------------------------+
!
!
   vect(1) =   u_rad;   vect(2) =  u_teta
   vect(3) =  u_long;   vect(4) =   geoid
 d_vect(1) =  du_rad; d_vect(2) = du_teta
 d_vect(3) = du_long; d_vect(4) = d_geoid
!
!
!
!
!+----------------------------------+
! Correction for a realistic ocean  !
!+----------------------------------+ 
!
   If (ICO == 1) Then 
       Call OCEAN_CORRECTION (OBSE, CORR, D_CORR) 
         vect(:) =   vect(:) +   corr(:) 
       d_vect(:) = d_vect(:) + d_corr(:) 
   Endif 
!
!
 END SUBROUTINE ESSA_IEZI0 
!
!
!
!
!
!
!
!
!
!
!
!+----------------------------------------------------------------+ 
   SUBROUTINE AXIS_STOK  (Ld, Md, IOC, inorm, CEN, VECT, D_VECT)  !           
!+----------------------------------------------------------------+
    use COMMON; use TRIGFCN; implicit NONE
!
!	This routine computes the cosine (CI) and the sine STOKES
!   coefficients (SI) at degree LD and order MD for a load of
!   axissimmetric type with center at longitude = CEN(1) and
!   colatitude = CEN(2). It is to be called after the call to
!   AXIS_LOAD, which computes the coefficients sigma(l), and after
!   convol_load_history, which convolves the 'k' ldc with the load
!   time-history. For INORM ==0, the Stokes coefficients are computed
!   with the usual normalization, for INORM==1, they take the 'FULLY
!   NORMALIZED' form. The output is  VECT = (Ci, Si); vector D_VECT =
!   (D_Ci, D_Si) is the time-derivative of VECT. See the Theory section
!   for further details.
!
! +---------------+
! | Declarations  |
! +---------------+
!!
    real(dp) ::   VECT(2)   !   ==   (Ci,   Si) 
    real(dp) :: D_VECT(2)   !   ==  (D_Ci, D_Si)   
!
    real(dp)  :: CEN(2), long_cen, teta_cen ! Center of the load
    real(dp)  :: leg                        ! Legendre function
    real(dp)  :: GPROD        	            ! == (l-m)!/(l+m)!     
    real(dp)  :: AUX, dom                   ! Auxiliary variables
!    
    real(dp)  :: Ci, Si          ! Cosine & Sine Stokes coefficients
    real(dp)  :: D_Ci, D_Si      ! their time-derivatives
!  
    integer (i4b) :: IOC             ! Realistic ocean variable
    integer (i4b) :: Ld, Md, L, M    ! Degrees & Orders
    integer (i4b) :: Inorm           ! Normalization switch
!
! +-----------------------+
! | Executable statements |
! +-----------------------+
!
  long_cen  = CEN(1) 
  teta_cen  = CEN(2)   
!
!
! +---------------------------------------------------------------+
! | (Change in the) STOKES coefficients & their time-derivatives  |
! +---------------------------------------------------------------+
!
	L=Ld; M=Md   
!
        if(m==0) dom=1._dp
        if(m/=0) dom=0._dp
!
        IF (IOC == 0 .or. IOC ==1) THEN
!
         AUX = (4._DP*PI*RAGGIO*RAGGIO/EMASS)*(2._DP-DOM) *     &
               (SIGMA(L)/(2._DP*L+1._DP))*GPROD(L,M)      *     &
               LEG(L,M,COSD(TETA_CEN))
!
         Ci    =   conv_k(L) * AUX * COSD (M*LONG_CEN) 
         Si    =   conv_k(L) * AUX * SIND (M*LONG_CEN)
!
         D_Ci  =  Dconv_k(L) * AUX * COSD (M*LONG_CEN) 
         D_Si  =  Dconv_k(L) * AUX * SIND (M*LONG_CEN)
!
        ENDIF
!
! +---------------------------------------------------------+
! | For IOC==1, a correction is made for load compensation  |
! +---------------------------------------------------------+
!
!
        IF(IOC==1) THEN 
!
         AUX = (4._DP*PI*RAGGIO*RAGGIO/EMASS)*(SIGMA(0)/OC(0,0))
!
         Ci    =  Ci  -  conv_k (L) * AUX * OC(L,M) / dfloat(2*l+1)
         Si    =  Si  -  conv_k (L) * AUX * OS(L,M) / dfloat(2*l+1)
!
         D_Ci  =  D_Ci  - Dconv_k (L) * AUX * OC(L,M)  / dfloat(2*l+1)
         D_Si  =  D_Si  - Dconv_k (L) * AUX * OS(L,M)  / dfloat(2*l+1)
!
!
        ENDIF 
!
!
        IF( IOC < 0 .OR. IOC > 1) THEN
          WRITE(99,*) ' Error IN SBR. AXIS_STOK:   The ocean switch IOC'
          WRITE(99,*) ' can only assume values 0, 1, or 2.  JOB ABORTED'; Stop
        ENDIF
!
!
! +---------------------------------------------------+
! | For Inorm==1, we employ the fully normalized form |
! +---------------------------------------------------+
!
         If (Inorm == 1) Then 
!	 
         Ci = Ci /Sqrt(2d0-dom)/Sqrt(2d0*l+1d0)/Sqrt(gprod(l,m))          
         Si = Si /Sqrt(2d0-dom)/Sqrt(2d0*l+1d0)/Sqrt(gprod(l,m))     
!	 
         D_Ci = D_Ci /Sqrt(2d0-dom)/Sqrt(2d0*l+1d0)/Sqrt(gprod(l,m))          
         D_Si = D_Si /Sqrt(2d0-dom)/Sqrt(2d0*l+1d0)/Sqrt(gprod(l,m))   	 
!	 
	 Endif   
!
!
! +-----------------------------+
! | Definition of Vect & D_Vect |
! +-----------------------------+
!
        Vect(1) = Ci; Vect(2)=Si; D_Vect(1)=D_Ci; D_Vect(2)=D_Si
!
        END SUBROUTINE AXIS_STOK 
!
!
!
!
!
!
!
!
!
!
!
!
!+-------------------------------------------------+ 
    SUBROUTINE AXIS_INER  (IOC, CEN, INE, D_INE)   !           
!+-------------------------------------------------+
    use COMMON; use TRIGFCN; implicit NONE
!
!	This routine computes the change in the inertia tensor of
!   the Earth INE and its time--derivatine D_INE for an axis-simmetric
!   load with center at longitude = CEN(1) and colatitude = CEN(2).
!   It is to be called after the call to AXIS_LOAD, which computes the
!   coefficients sigma(l), and after  convol_load_history, which convolves
!   the 'k' ldc with the load time-history. See the Theory section
!   for further details.
!
! +----------------+
! |# Declarations  |
! +----------------+
!!
    real(dp)  ::   INE(3,3)   ! Change of the inertia tensor
    real(dp)  :: D_INE(3,3)   ! Its time-derivative
! 
    real(dp)  :: CEN(2), long_cen, teta_cen       ! Load center
    real(dp)  :: leg                              ! Legendre function
    real(dp)  :: AUX, dom                         ! Auxiliary variables
!
    integer (i4b) :: IOC      ! Realistic ocean switch
!
!
! +-----------------------+
! | Executable statements |
! +-----------------------+
!
  long_cen  = CEN(1); teta_cen  = CEN(2)   
!
!
!+---------------------------------------------------------+
!  Change in the inertia tensor and their time-derivatives |
!+---------------------------------------------------------+
!
!
       IF ( IOC ==0 .or. IOC ==1 ) THEN
!
       INE (1,1) = (4._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  & 
                       SIGMA(2)*conv_k(2)*LEG(2,0,COSD(TETA_CEN))   -  &
                       (2._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  &
                       SIGMA(2)*conv_k(2)*LEG(2,2,COSD(TETA_CEN))   *  &
                       COSD(2._DP*LONG_CEN)
!
       INE (2,2) = (4._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  & 
                       SIGMA(2)*conv_k(2)*LEG(2,0,COSD(TETA_CEN))   +  &
                       (2._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  &
                       SIGMA(2)*conv_k(2)*LEG(2,2,COSD(TETA_CEN))   *  &
                       COSD(2._DP*LONG_CEN)
!
       INE (3,3) = - (8._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)         *  & 
                       SIGMA(2)*conv_k(2)*LEG(2,0,COSD(TETA_CEN))      
!
       INE (1,3) = (4._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  & 
                       SIGMA(2)*conv_k(2)*LEG(2,1,COSD(TETA_CEN))   *  & 
                       COSD(1._DP*LONG_CEN)
!
       INE (2,3) = (4._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  & 
                       SIGMA(2)*conv_k(2)*LEG(2,1,COSD(TETA_CEN))   *  & 
                       SIND(1._DP*LONG_CEN)
!
       INE (1,2) = - (2._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)         *  & 
                       SIGMA(2)*conv_k(2)*LEG(2,2,COSD(TETA_CEN))   *  &     
                       SIND(2._DP*LONG_CEN)
!
       INE (2,1) = INE (1,2) 
       INE (3,1) = INE (1,3) 
       INE (3,2) = INE (2,3)
!
!
!  +--------------+
!  | Derivatives  |
!  +--------------+
!
!   
       D_INE (1,1) = (4._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  & 
                       SIGMA(2)*Dconv_k(2)*LEG(2,0,COSD(TETA_CEN))   -  &
                       (2._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  &
                       SIGMA(2)*Dconv_k(2)*LEG(2,2,COSD(TETA_CEN))   *  &
                       COSD(2._DP*LONG_CEN)
!
       D_INE (2,2) = (4._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  & 
                       SIGMA(2)*Dconv_k(2)*LEG(2,0,COSD(TETA_CEN))   +  &
                       (2._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  &
                       SIGMA(2)*Dconv_k(2)*LEG(2,2,COSD(TETA_CEN))   *  &
                       COSD(2._DP*LONG_CEN)
!
       D_INE (3,3) = - (8._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)         *  & 
                       SIGMA(2)*Dconv_k(2)*LEG(2,0,COSD(TETA_CEN))      
!
       D_INE (1,3) = (4._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  & 
                       SIGMA(2)*Dconv_k(2)*LEG(2,1,COSD(TETA_CEN))   *  & 
                       COSD(1._DP*LONG_CEN)
!
       D_INE (2,3) = (4._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)           *  & 
                       SIGMA(2)*Dconv_k(2)*LEG(2,1,COSD(TETA_CEN))   *  & 
                       SIND(1._DP*LONG_CEN)
!
       D_INE (1,2) = - (2._DP*PI*RAGGIO*RAGGIO/15._DP/EMASS)         *  & 
                       SIGMA(2)*Dconv_k(2)*LEG(2,2,COSD(TETA_CEN))   *  &     
                       SIND(2._DP*LONG_CEN)
!
       D_INE (2,1) = D_INE (1,2) 
       D_INE (3,1) = D_INE (1,3) 
       D_INE (3,2) = D_INE (2,3)
!
       ENDIF
!
!
! +---------------------------------------------------------+
! | For IOC==1, a correction is made for load compensation  |
! +---------------------------------------------------------+
!
!
        IF (IOC==1)                                            THEN 
!
!
        INE(1,1) = INE (1,1) -                                & 
            (4._DP*PI*RAGGIO*RAGGIO/3._DP/EMASS)*SIGMA(0)*conv_k(2)*  & 
            ( OC(2,0)/OC(0,0) - 6._DP*OC(2,2)/OC(0,0))/5._DP
!
        INE(2,2) = INE(2,2) -                                & 
            (4._DP*PI*RAGGIO*RAGGIO/3._DP/EMASS)*SIGMA(0)*conv_k(2)*  & 
            ( OC(2,0)/OC(0,0) + 6._DP*OC(2,2)/OC(0,0))/5._DP
!
        INE(3,3) = INE (3,3) + 		                      &
            (8._DP*PI*RAGGIO*RAGGIO/3._DP/EMASS)*SIGMA(0)*conv_k(2)*  & 
              OC(2,0)/OC(0,0)/5._DP
!
        INE(1,3) = INE (1,3) -                                &
            (4._DP*PI*RAGGIO*RAGGIO/1._DP/EMASS)*SIGMA(0)*conv_k(2)*  & 
              OC(2,1)/OC(0,0)/5._DP
!
        INE(2,3) = INE (2,3) - 			              & 
            (4._DP*PI*RAGGIO*RAGGIO/1._DP/EMASS)*SIGMA(0)*conv_k(2)*  & 
              OS(2,1)/OC(0,0)/5._DP
!
        INE(1,2) = INE (1,2) -                                & 
            (8._DP*PI*RAGGIO*RAGGIO/1._DP/EMASS)*SIGMA(0)*conv_k(2)*  & 
              OS(2,2)/OC(0,0)/5._DP
!
        INE (2,1) = INE (1,2) 
        INE (3,1) = INE (1,3) 
        INE (3,2) = INE (2,3) 
!
!
!  +--------------+
!  | Derivatives  |
!  +--------------+
!  
!
        D_INE(1,1) = D_INE (1,1) -                                & 
            (4._DP*PI*RAGGIO*RAGGIO/3._DP/EMASS)*SIGMA(0)*Dconv_k(2)*  & 
            ( OC(2,0)/OC(0,0) - 6._DP*OC(2,2)/OC(0,0))/5._DP
!
        D_INE(2,2) = D_INE(2,2) -                                & 
            (4._DP*PI*RAGGIO*RAGGIO/3._DP/EMASS)*SIGMA(0)*Dconv_k(2)*  & 
            ( OC(2,0)/OC(0,0) + 6._DP*OC(2,2)/OC(0,0))/5._DP
!
        D_INE(3,3) = D_INE (3,3) + 		                      &
            (8._DP*PI*RAGGIO*RAGGIO/3._DP/EMASS)*SIGMA(0)*Dconv_k(2)*  & 
              OC(2,0)/OC(0,0)/5._DP
!
        D_INE(1,3) = D_INE (1,3) -                                &
            (4._DP*PI*RAGGIO*RAGGIO/1._DP/EMASS)*SIGMA(0)*Dconv_k(2)*  & 
              OC(2,1)/OC(0,0)/5._DP
!
        D_INE(2,3) = D_INE (2,3) - 			              & 
            (4._DP*PI*RAGGIO*RAGGIO/1._DP/EMASS)*SIGMA(0)*Dconv_k(2)*  & 
              OS(2,1)/OC(0,0)/5._DP
!
        D_INE(1,2) = D_INE (1,2) -                                & 
            (8._DP*PI*RAGGIO*RAGGIO/1._DP/EMASS)*SIGMA(0)*Dconv_k(2)*  & 
              OS(2,2)/OC(0,0)/5._DP
!
        D_INE (2,1) = D_INE (1,2) 
        D_INE (3,1) = D_INE (1,3) 
        D_INE (3,2) = D_INE (2,3) 
!
        ENDIF
!
!
END SUBROUTINE AXIS_INER    
!
!
!
!
!
!
!
!
!
!   
!
! +---------------------------------------------------------+
    SUBROUTINE RECT_STOK (Ld, Md, IOC, Inorm, VECT, D_VECT) !
! +---------------------------------------------------------+
!
   use COMMON; implicit NONE
!
!	This routine computes the cosine (CI) and the sine STOKES
!   coefficients (SI) at degree LD and order MD for a load of
!   type 50 (i.e. a 'rectangular' load). It is to be called after
!   the call to RECT_LOAD0, which computes the harmonic coefficients
!   of the load, and after convol_load_history, which convolves the 'k'
!   ldc with the load time-history. For INORM ==0, the Stokes coefficients
!   are computed with the usual normalization, for INORM==1, they take
!   the 'FULLY NORMALIZED' form. The output is  VECT = (Ci, Si);
!   vector D_VECT = (D_Ci, D_Si) is the time-derivative of VECT. See
!   the Theory section for further details.
!
! +------------+
! | Statements |
! +------------+
!!
    real(dp) ::   VECT(2)   !   ==   (Ci,   Si) 
    real(dp) :: D_VECT(2)   !   ==  (D_Ci, D_Si)    
!
    real(dp)  :: leg                        ! Legendre function
    real(dp)  :: AUX, dom                   ! Auxiliary variables
!    
    real(dp)  ::   Ci,   Si      ! Cosine & Sine Stokes coefficients
    real(dp)  :: D_Ci, D_Si      ! Their time-derivatives
!  
    real(dp)  :: GPROD           ! Gprod == (l-m)!/(l+m)!
!
    integer (i4b) :: IOC             ! Realistic ocean switch
    integer (i4b) :: Ld, Md, L, M    ! Degrees & orders
    integer (i4b) :: Inorm           ! The normalization switch
!
!
! +-----------------------+
! | Executable statements |
! +-----------------------+
!
        L=Ld; M=Md 
!
        If(m == 0) Dom=1D0 
	If(m /= 0) Dom=0D0 
!

        IF( IOC == 0 .OR. IOC ==1 ) THEN
!
         AUX = (4._DP*PI*RAGGIO*RAGGIO/EMASS)/(2._DP*L+1._DP)
!
         Ci  =     conv_k(L) * AUX * FF (L,M)
         Si  =     conv_k(L) * AUX * GG (L,M)
!
         D_Ci  =  Dconv_k(L) * AUX * FF (L,M)
         D_Si  =  Dconv_k(L) * AUX * GG (L,M)
!
        ENDIF
!
! +---------------------------------------------------------+
! | For IOC==1, a correction is made for load compensation  |
! +---------------------------------------------------------+
!
        IF(IOC==1) THEN 
!
         AUX = (4._DP*PI*RAGGIO*RAGGIO/EMASS)*(FF(0,0)/OC(0,0))
!
         Ci = Ci  - conv_k (L) * AUX * OC(L,M)  / dfloat(2*l+1)
         Si = Si  - conv_k (L) * AUX * OS(L,M)  / dfloat(2*l+1)
!
         D_Ci  =  D_Ci  - Dconv_k (L) * AUX * OC(L,M) / dfloat(2*l+1)
         D_Si  =  D_Si  - Dconv_k (L) * AUX * OS(L,M) / dfloat(2*l+1)
!
        ENDIF
!
!
! +---------------------------------------------------+
! | For Inorm==1, we employ the fully normalized form |
! +---------------------------------------------------+
!
   
         If (Inorm == 1) Then 
!	 
         Ci = Ci /Sqrt(2d0-dom)/Sqrt(2d0*l+1d0)/Sqrt(gprod(l,m))          
         Si = Si /Sqrt(2d0-dom)/Sqrt(2d0*l+1d0)/Sqrt(gprod(l,m))     
!	 
         D_Ci = D_Ci /Sqrt(2d0-dom)/Sqrt(2d0*l+1d0)/Sqrt(gprod(l,m))          
         D_Si = D_Si /Sqrt(2d0-dom)/Sqrt(2d0*l+1d0)/Sqrt(gprod(l,m))   	 
!	 
	 Endif   
!
!
! +--------------------------------+
! |  Definition of VECT & D_VECT   !
! +--------------------------------+
!
!
        Vect(1) = Ci; Vect(2)=Si; D_Vect(1)=D_Ci; D_Vect(2)=D_Si 
!
!
	END SUBROUTINE RECT_STOK 
!
!
!
!
!
!
!
!
!
!
!+--------------------------------------------+ 
    SUBROUTINE RECT_INER  (IOC, INE, D_INE)   !           
!+--------------------------------------------+
    use COMMON; implicit NONE 
!
!
!	This routine computes the change in the inertia tensor of
!   the Earth INE and its time--derivatine D_INE for a 'rectangular'
!   load (type 50). It is to be called after the call to RECT_LOAD0,
!   which computes the coefficients sigma(l), and after convol_load_history,
!   which convolves the 'k' ldc with the load time-history. See the
!   Theory section for further details.
!
!
! # Statements
!!
    real(dp)  ::    INE(3,3) ! Change in the inertia tensor
    real(dp)  ::  D_INE(3,3) ! Its time-derivative
!
    real(dp)  ::  Ci20,  Ci21,  Ci22,  Si21,  Si22 ! Degree 2 Stokes coefficients
    real(dp)  :: DCi20, DCi21, DCi22, DSi21, DSi22 ! Theit time-derivatives
!    
    real(dp)  :: VECT(4), v1(2), v2(2)  ! Auxiliary arrays
!
    integer (i4b) :: IOC     ! Realistic ocean switch
!
!
!
! +----------------------------------------------------------------------+
! | Computing the degree 2 Stokes coefficients, & their time-derivatives |
! +----------------------------------------------------------------------+
!
       Call Rect_Stok (2,0, ioc, 0, v1, v2)    
        Ci20 = v1(1)    
       Dci20 = v2(1)    
!       
       Call Rect_Stok (2,1, ioc, 0, v1, v2)    
        Ci21 = v1(1)  
	Si21 = v1(2)  
       Dci21 = v2(1) 
       Dsi21 = v2(2)   
!
       Call Rect_Stok (2,2, ioc, 0, v1, v2)   
        Ci22 = v1(1) 
	Si22 = v1(2)  
       Dci22 = v2(1)  
       Dsi22 = v2(2) 
!
! +----------------+
! | Inertia change |
! +----------------+
!
       INE (1,1) = Ci20/3._DP - 2._DP*Ci22 
       INE (2,2) = Ci20/3._DP + 2._DP*Ci22 
       INE (3,3) = - 2._DP*Ci20/3._DP 
!
       INE (1,3) = Ci21 
       INE (2,3) = Si21 
       INE (1,2) = -2._DP*Si22 
!
       INE (2,1) = INE (1,2) 
       INE (3,1) = INE (1,3) 
       INE (3,2) = INE (2,3) 
!
! +----------------------------------------+
! | Time derivatives of the inertia change |
! +----------------------------------------+
!
!
       D_INE (1,1) = DCi20/3._DP - 2._DP*DCi22 
       D_INE (2,2) = DCi20/3._DP + 2._DP*DCi22 
       D_INE (3,3) = - 2._DP*DCi20/3._DP 
!
       D_INE (1,3) = DCi21 
       D_INE (2,3) = DSi21 
       D_INE (1,2) = -2._DP*DSi22 
!
       D_INE (2,1) = D_INE (1,2) 
       D_INE (3,1) = D_INE (1,3) 
       D_INE (3,2) = D_INE (2,3) 
!
!
END SUBROUTINE RECT_INER 
!
!
!
!
!
!
!
!
!
!
!
!
!
!
! +-------------------------------+
        subroutine TASK_1         !
! +-------------------------------+
!
 use COMMON; implicit NONE
!
!
 INTEGER (i4b) :: i_loading    ! loading/tidal switch
 INTEGER (i4b) :: ideja        ! checks for errors in the model KWs
 INTEGER (i4b) :: inone        ! checks for errors in the model KWs
 INTEGER (i4b) :: iarmo        ! detects the Harmonic_Degrees KW
 INTEGER (i4b) :: iexte        ! detects the kw External_Model  
 integer (i4b) :: CODE         ! Code to identify the models
 integer (i4b) :: ILM          ! (0/1) Controls the lower mantle stratification
 integer (i4b) :: ih, il, ik   ! Inputs for the (h/s) analysis
!
 integer (i4b) :: J, l, kk, i, k, m, ITASK      ! do-loops control variables
 integer (i4b) :: NUMD                          ! number of (harmonic) degree
 integer (i4b) :: NPUN                          ! number of time points
!
 integer (i4b), parameter :: max_num_righe = 20000  ! Max. n. of rows of file task_1.dat
 integer (i4b), parameter :: npun_max = 2001      ! Max number of time points
 integer (i4b), parameter :: numd_max = 6       ! Max. number of Harmonics
 integer (i4b) :: D(numd_max)                 ! Vector of harmonic degrees
!
 REAL(DP) :: T_MIN, T_MAX, P_MIN, P_MAX, ALFA, BETA, TIME   ! For the Heaviside TH
 REAL(DP) :: RESPH, RESPL, RESPK                            ! For the HEaviside TH
!
 real(sp) :: ajunk 
!
 Real(QP) :: LT 	 ! Lithospheric Thickness
!
 CHARACTER*30 KEYWORD    ! A keyword
 CHARACTER*20 date, timc ! date and time
!
!
!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !  Beginning with the executable instructions  !
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! # These integers are used to monitor the 
!     keywords sequence .... 
!
         iarmo=0
          ideja=0
        inone=0
      iexte=0
!
!
Write(99,*) 'Opening file  task_1.dat ...  '
!
! 
      open(1,file='task_1.dat',status='unknown')
!
!
Write(99,*) 'Reading file task_1.dat ...  '
! 
!
!
!
!
      DO 101 ITASK = 1, max_num_righe
!
!
      Read(1,'(a30)',End=102) KEYWORD
!
!
!
!
!
      If(KEYWORD(1:16) == 'Harmonic_Degrees') Then
!
!
       iarmo=1
!
!
       Read(1,*)  LMIN, LMAX
!
       READ(1,*)  IV
!
       IF(IV/=0.and.IV/=1) THEN
       WRITE(99,*) 'ERROR in sbr task_1.dat:     The VERBOSE switch '
       WRITE(99,*) 'can only assume values 0 and 1. ** JOB ABORTED* '; STOP
       endif
!
        call DATE_AND_TIME (date,timc)
          IF    (IV==0) THEN
    Write(99,*) '# Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
          ELSEIF(IV==1) then
    Write(99,*) '# Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
        ENDIF
!
    IF(iv==1) then
     WRITE(*,*) '= = = = = = = = = = = = = = ='
     Write(*,*) '   Reading file task_1.dat   '
     WRITE(*,*) '= = = = = = = = = = = = = = ='
    ENDIF
!
       READ(1,*)  i_loading
!
       IF((i_loading/=0.and.i_loading/=1).and.IV==0) THEN
       WRITE(99,*)'ERROR in sbr task_1.dat:     The i_loading switch '
       WRITE(99,*)'can only assume values 0 and 1. **** JOB ABORTED* '; STOP
       ENDIF
       IF((i_loading/=0.and.i_loading/=1).and.IV==1) THEN
       WRITE(99,*)'ERROR in sbr task_1.dat:     The i_loading switch '
       WRITE(99,*)'can only assume values 0 and 1. **** JOB ABORTED* '
       WRITE(*, *)'ERROR in sbr task_1.dat:     The i_loading switch '
       WRITE(*, *)'can only assume values 0 and 1. **** JOB ABORTED* '; STOP
       ENDIF
!
!
                     WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
!
    IF(IV==1 .AND. i_loading==1) &
                     WRITE(* ,*) 'Loading analysis'
    IF(IV==0 .AND. i_loading==1) &
                     WRITE(99,*) 'Loading analysis'
    IF(IV==1 .AND. i_loading==0) &
                     WRITE(* ,*) 'Tidal analysis'
    IF(IV==0 .AND. i_loading==0) &
                     WRITE(99,*) 'Tidal analysis'
!
!
                   Write(99,*) 'Lmin and Lmax ', Lmin, Lmax
    IF(IV==1) Write (*,*) 'Lmin and Lmax ', Lmin, Lmax
!
!
! # Checking Lmin and Lmax
!
IF(IV==0) then
If(lmin<llmin .or. lmax>llmax .or. lmax<lmin .or. lmin<0 .or. lmax<0) Then
Write(99,*)'ERROR IN SBR. TASK_1:     One of the following forbidden conditions'
Write(99,*)'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
Write(99,*)'Check the input file task_2.dat ****** JOB ABORTED ****************';Stop
Endif
Endif
IF(IV==1) then
If(lmin<llmin .or. lmax>llmax .or. lmax<lmin .or. lmin<0 .or. lmax<0) Then
Write(99,*)'ERROR IN SBR. TASK_1:     One of the following forbidden conditions'
Write(99,*)'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
Write(99,*)'Check the input file task_2.dat ****** JOB ABORTED ****************'
Write(*,*) 'ERROR IN SBR. TASK_1      One of the following forbidden conditions'
Write(*,*) 'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
Write(*,*) 'Check the input file task_2.dat ****** JOB ABORTED ****************';Stop
Endif
Endif
!
!
!
!
  Endif  ! on KW Harmonic_Degrees
!
!
!
!
!:::::::::::::::::::::::::::::::::::::::::
   IF(KEYWORD(1:10) == 'Make_Model') THEN   
!:::::::::::::::::::::::::::::::::::::::::
!
   ideja=1
   inone=1
!
    IF(iarmo==0 .and. IV==0) then
                 WRITE(99,*) ' ERROR IN SBR. TASK_1: the KW Harmonic_Degrees '
                 WRITE(99,*) ' is not active. **** JOB ABORTED ************* ';stop
    endif
    IF(iarmo==0 .and. IV==1) then
                 WRITE(99,*) ' ERROR IN SBR. TASK_1: the KW Harmonic_Degrees '
                 WRITE(99,*) ' is not active. **** JOB ABORTED ************* '
                 WRITE(*,*)  ' ERROR IN SBR. TASK_1: the KW Harmonic_Degrees '
                 WRITE(*,*)  ' is not active. **** JOB ABORTED ************* ';STOP
    endif
!
!
!
                     WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
!
                  Write(99,*)    'Building the model '
 IF(IV==1)   Write (*,*)    'Building the model '
!
   READ(1,*) NV 
!  
   READ(1,*) CODE
!
   READ(1,*) LT 
!
   READ(1,*) ILM 
!      
!
   Write(99,*) 'Number of VE layers      ', NV  
   Write(99,*) 'Model CODE is            ', CODE 
!
IF(IV==1) then
   Write(*,*) 'Number of VE layers      ', NV
   Write(*,*) 'Model CODE is            ', CODE
ENDIF
!
!
    Write(99,'(a42,1x,f12.4,1x,a3)') &
    'From input, the lithospheric thickness is', LT,  ' km'
    Write(99,'(a47,1x,i3)')       &
    'The ILM parameter (for NV=7 and 9) is set to  ', ILM
!
   IF(IV==1) THEN
    Write(*, '(a42,1x,f12.4,1x,a3)') &
    'From input, the lithospheric thickness is', LT,  ' km'
    Write(*, '(a47,1x,i3)')       &
    'The ILM parameter (for NV=7 and 9) is set to  ', ILM
   ENDIF
!
!
               Write(99,*) 'Mantle viscosity from BOTTOM to TOP (/1E21) '
IF(IV==1) WRITE(*,*)  'Mantle viscosity from BOTTOM to TOP (/1E21) '
   Do k=1, nv
      Read(1,*) vis(k)
	          Write(99, '(a19,i4,1x,a1F10.4)') &
                       'Viscosity of layer', k, '=', vis (k)
   IF(IV==1) Write(*,  '(a19,i4,1x,a1,F10.4)') &
                       'Viscosity of layer', k, '=', vis (k)
               vis(k) = vis(k) * 1.q21    	  
   Enddo 
! 
!
   CALL SPEC (Code, Ilm, LT)       

!
! +---------------------------------
       call SPECTRUM (i_loading)
! +---------------------------------
!
!
!
            Write(99,*)'Writing the spectrum on file spectrum.dat'
 IF(IV==1)  Write(*,* )'Writing the spectrum on file spectrum.dat'
!
!
!
  open(3,file='spectrum.dat',status='unknown')
!
  call DATE_AND_TIME (date,timc)
!
  Write(3,*) '# File spectrum.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  
  Write(3,*) '  # 1st column: l = Harmonic degree                '
  Write(3,*) '  # 2nd   "   : LOG10 (l)                          '
  Write(3,*) '  # 3th   "   : s, (kyrs**(-1))                    '
  Write(3,*) '  # 4rd   "   : LOG10(-s), with s in kyrs**(-1)    '
  Write(3,*) '  # 5th   "   : Relaxation time = -1000./s, (yrs)  '
  Write(3,*) '  # 6th   "   : LOG10 (Relaxation time (yrs))      '
!
DO l=lmin, lmax 
! 
    write (3, '(/)')
      do kk = 1, nroots  
       write (3, '(i4,1x,5(e15.7,1x))') l,                  & 
                                        log10 (dfloat (l)), &
                                        s(l,kk),            &
					log10(-s(l,kk)),    &
				        (-1000.q0/s(l,kk)), &
				        log10(-1000.q0/s(l,kk)) 
      enddo
!
ENDDO; close(3)
!
!
!::::::::::::::::::::::::::::::::::
   ENDIF  !  Endif on KW Make_Model
!::::::::::::::::::::::::::::::::::
!
!
!
!
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
   IF(KEYWORD(1:14) == 'External_Model') THEN
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
   inone=1
   iexte=1  
!   
   rhoea = ext_rhoea
!
!
!
    IF(ideja==1.and.IV==0) then
           WRITE(99,*)'ERROR IN SBR. TASK_1: The Kws Make_Model and External_Model '
           WRITE(99,*)'are mutually exclusive. If Make_Model is active,  External_ '
           WRITE(99,*)'model must be not. **** JOB ABORTED *********************** '
           STOP
    endif
    IF(ideja==1.and.IV==1) then
           WRITE(99,*)'ERROR IN SBR. TASK_1: The Kws Make_Model and External_Model '
           WRITE(99,*)'are mutually exclusive. If Make_Model is active,  External_ '
           WRITE(99,*)'model must be not. **** JOB ABORTED *********************** '
           WRITE(*,*) 'ERROR IN SBR. TASK_1: The Kws Make_Model and External_Model '
           WRITE(*,*) 'are mutually exclusive. If Make_Model is active,  External_ '
           WRITE(*,*) 'model must be not. **** JOB ABORTED *********************** '
           STOP
    endif
!
!
!
                     WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
!
 IF(i_loading==0.and.IV==0) THEN
 WRITE(99,*)'ERROR in sbr task_1.dat:     If the External_Model KW is  '
 WRITE(99,*)'set, the TIDAL analysis is not possible. **JOB ABORTED**  '; STOP
 ENDIF
 IF(i_loading==0.and.IV==1) THEN
 WRITE(99,*)'ERROR in sbr task_1.dat:     If the External_Model KW is  '
 WRITE(99,*)'set, the TIDAL analysis is not possible. **JOB ABORTED**  '
 WRITE(* ,*)'ERROR in sbr task_1.dat:     If the External_Model KW is  '
 WRITE(* ,*)'set, the TIDAL analysis is not possible. **JOB ABORTED**  '; STOP
 ENDIF
!
!
!
   Write(99,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '
   Write(99,*) ' The model is the one by   GIUNCHI & SPADA,  [2000] '
   Write(99,*) ' NV=1 and LT=100 km, see  GRL table (1) for details '
   Write(99,*) ' The model is Gravitating but not self -Gravitating '
   Write(99,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '    
!
   IF(IV==1) THEN
   Write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '
   Write(*,*) ' The model is the one by   GIUNCHI & SPADA,  [2000] '
   Write(*,*) ' NV=1 and LT=100 km, see  GRL table (1) for details '
   Write(*,*) ' The model is Gravitating but not self -Gravitating '
   Write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '
   ENDIF
!
    NV=1    
    NROOTS=3 
!    
 Open(67, file='external_spectrum.dat',status='unknown')   
   do J=0, llmax 
   do kk=1, nroots 
        read(67,'(i4,1x,2(e15.7,1x))',end=1901)    l,   s(l,kk) ,  s(l,kk)  
   enddo
   enddo
 1901 Close(67)
! stop 33411
!  
  Open(67, file='external_h.dat',status='unknown') 
          DO J=0, llmax 
             read (67, '((i3,1x,96(1x,e20.8)))', end=1902) &
	     l, h_e(l), h_f(l), (h_v (l,k), k = 1, nroots) 
	  ENDDO            
 1902 Close(67)   
  Open(68, file='external_l.dat',status='unknown')  
          DO J=0, llmax 
             read (68, '((i3,1x,96(1x,e20.8)))', end=1903) &
	     l, l_e(l), l_f(l), (l_v (l,k), k = 1, nroots)    
          ENDDO                             
 1903 close(68)                                     
      vec=1.
!    
!   
!
!
                 Write(99,*)'Writing the spectrum on file spectrum.dat'
 IF(IV==1)  Write(*,* )'Writing the spectrum on file spectrum.dat'
!
!
  open(3,file='spectrum.dat',status='unknown')
!
  call DATE_AND_TIME (date,timc)
!
  Write(3,*) '# File spectrum.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  
!
  Write(3,*) '  # 1st column: l = Harmonic degree                '
  Write(3,*) '  # 2nd   "   : LOG10 (l)                          '
  Write(3,*) '  # 3th   "   : s, (kyrs**(-1))                    '
  Write(3,*) '  # 4rd   "   : LOG10(-s), with s in kyrs**(-1)    '
  Write(3,*) '  # 5th   "   : Relaxation time = -1000./s, (yrs)  '
  Write(3,*) '  # 6th   "   : LOG10 (Relaxation time (yrs))      '
!
!
DO l=lmin, lmax
! 
   write (3, '(/)')
      do kk = 1, nroots  
       write (3, '(i4,1x,5(e15.7,1x))') l,                  & 
                                        log10 (dfloat (l)), &
                                        s(l,kk),            &
					log10(-s(l,kk)),    &
				        (-1000.q0/s(l,kk)), &
				        log10(-1000.q0/s(l,kk)) 
      enddo
!
ENDDO; close(3)
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ENDIF   ! Endif on KW External_Model
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++
      If(KEYWORD(1:19)=='Normalized_Residues') Then
!++++++++++++++++++++++++++++++++++++++++++++++++++++
!
    IF(inone==0.and.IV==0) then
           WRITE(99,*)'ERROR IN SBR. TASK_1:  One of the Kws Make_Model     and '
           WRITE(99,*)'External_Model MUST be active before Normalized_Residues '
           WRITE(99,*)'can be executed. ******* JOB ABORTED ******************* ';STOP
    endif
    IF(inone==0.and.IV==1) then
           WRITE(99,*)'ERROR IN SBR. TASK_1:  One of the Kws Make_Model     and '
           WRITE(99,*)'External_Model MUST be active before Normalized_Residues '
           WRITE(99,*)'can be executed. ******* JOB ABORTED ******************* '
           WRITE(*,*) 'ERROR IN SBR. TASK_1:  One of the Kws Make_Model     and '
           WRITE(*,*) 'External_Model MUST be active before Normalized_Residues '
           WRITE(*,*) 'can be executed. ******* JOB ABORTED ******************* ';STOP
    endif
!
                WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
                Write(99,*) 'Computing the normalized residues '
    IF(IV==1)   WRITE(*,*)  'Computing the normalized residues '
!
!
      Read(1,*) ih
      Read(1,*) il
      Read(1,*) ik 
!
      if(ik==1.and.iv==0.and.iexte==1) then 
        write(99,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
        write(99,*) 'WARNING in sbr TASK1: When the kw External_Model is active, '
	write(99,*) 'the k ldcs are zero by definition. I set IK==0 and continue '
        write(99,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
	IK=0
      endif 	 
      if(ik==1.and.iv==1.and.iexte==1) then 
        write(* ,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '      
        write(* ,*) 'WARNING in sbr TASK1: When the kw External_Model is active, '
	write(* ,*) 'the k ldcs are zero by definition. I set IK==0 and continue '
        write(* ,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
	IK=0
      endif 
!
!
        If(ih==0.and.il==0.and.ik==0.and.IV==0) Then
	Write(99,*) 'WARNING IN SBR. TASK_1:       No normalized residues '
	Write(99,*) 'will be printed: ih=il=ik =0 in input file task_1.dat'
	Endif
!
        If(ih==0.and.il==0.and.ik==0.and.IV==1) Then
	Write(* ,*) 'WARNING IN SBR. TASK_1:       No normalized residues '
	Write(* ,*) 'will be printed: ih=il=ik =0 in input file task_1.dat'
	Endif 
!
!
!
        IF(IH==1) THEN 
!
          open(3,file='ih.dat',status='unknown')
!
            WRITE(99,*) 'Writing on ih.dat the h normalized residues'
  IF(IV==1) WRITE(* ,*) 'Writing on ih.dat the h normalized residues'
!
  call DATE_AND_TIME (date,timc)
  Write(3,*) '# File ih.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  If(i_loading==1) write(3,*) '# These h normalized residues are of loading type'
  If(i_loading==0) write(3,*) '# These h normalized residues are of _tidal_ type'  	     
!
          DO l=lmin, lmax 
          write (3, '(/)')
           do kk = 1, nroots  
            write (3, '(i4,1x,e15.7)') l, - h_v(l,kk)/s(l,kk)
           enddo
          ENDDO; close(3) 
!
        ENDIF 
!
        IF(IL==1) THEN
!
          open(3,file='il.dat',status='unknown')
!
            WRITE(99,*) 'Writing on il.dat the l normalized residues'
  IF(IV==1) WRITE(* ,*) 'Writing on il.dat the l normalized residues'
!
!
  call DATE_AND_TIME (date,timc)
  Write(3,*) '# File il.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  If(i_loading==1) write(3,*) '# These l normalized residues are of loading type'
  If(i_loading==0) write(3,*) '# These l normalized residues are of _tidal_ type'  	     
!	  
!	  
          DO l=lmin, lmax 
          write (3, '(/)')
           do kk = 1, nroots  
            write (3, '(i4,1x,e15.7)') l, - l_v(l,kk)/s(l,kk)
           enddo
          ENDDO; close(3) 
!
        ENDIF 
!
        IF(IK==1) THEN 
!
          open(3,file='ik.dat',status='unknown') 
!
            WRITE(99,*) 'Writing on ik.dat the k normalized residues'
  IF(IV==1) WRITE(* ,*) 'Writing on ik.dat the k normalized residues'
!
!
  call DATE_AND_TIME (date,timc)
  Write(3,*) '# File ik.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  If(i_loading==1) write(3,*) '# These k normalized residues are of loading type'
  If(i_loading==0) write(3,*) '# These k normalized residues are of _tidal_ type'  	     
!	  
!
          DO l=lmin, lmax 
          write (3, '(/)')
           do kk = 1, nroots  
            write (3, '(i4,1x,e15.7)') l, - k_v(l,kk)/s(l,kk)
           enddo
          ENDDO; close(3) 
!
        ENDIF 
!
!+++++++++++++++++++++++++++++++++++++++++++++++
        Endif ! Endif on KW Normalized_Residues
!+++++++++++++++++++++++++++++++++++++++++++++++
!
! 
! 
!
!================================================
       If(KEYWORD(1:15)=='El_Fluid_Viscel') Then
!================================================
!
                WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
!       
      read(1,*) ih
      read(1,*) il
      read(1,*) ik
!
!
        If(ih==0.and.il==0.and.ik==0.and.IV==0) Then
	Write(99,*) 'WARNING IN SBR. TASK_1: No El., Fluid, Viscel amplitudes'
	Write(99,*) 'will be printed:    ih=il=ik =0 in input file task_1.dat'
	Endif
!
        If(ih==0.and.il==0.and.ik==0.and.IV==1) Then
	Write(* ,*) 'WARNING IN SBR. TASK_1: No El., Fluid, Viscel amplitudes'
	Write(* ,*) 'will be printed:    ih=il=ik =0 in input file task_1.dat'
	Endif 
!
      if(ik==1.and.iv==0.and.iexte==1) then 
        write(99,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
        write(99,*) 'WARNING in sbr TASK1: When the kw External_Model is active, '
	write(99,*) 'the k ldcs are zero by definition. I set IK==0 and continue '
        write(99,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
	IK=0
      endif 	 
      if(ik==1.and.iv==1.and.iexte==1) then 
        write(* ,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '      
        write(* ,*) 'WARNING in sbr TASK1: When the kw External_Model is active, '
	write(* ,*) 'the k ldcs are zero by definition. I set IK==0 and continue '
        write(* ,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
	IK=0
      endif 
!
        IF(IH==1) THEN 
          open(3,file='h.dat',status='unknown')
!
            WRITE(99,*) 'Writing on h.dat the elastic, fluid, and v-elastic h'
  IF(IV==1) WRITE(* ,*) 'Writing on h.dat the elastic, fluid, and v-elastic h'
!
!
  call DATE_AND_TIME (date,timc)
  Write(3,*) '# File h.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  If(i_loading==1) write(3,*) '# These h numbers are of loading type'
  If(i_loading==0) write(3,*) '# These h numbers are of _tidal_ type'  
!	  
!
          DO l=lmin, lmax 
          write (3, '((i3,1x,96(1x,e20.8)))') &
	  l, h_e(l), h_f(l), (h_v (l,k), k = 1, nroots)    
!
          ENDDO 
        ENDIF; close(3)  
!
        IF(IL==1) THEN
!
            WRITE(99,*) 'Writing on l.dat the elastic, fluid, and v-elastic l'
  IF(IV==1) WRITE(* ,*) 'Writing on l.dat the elastic, fluid, and v-elastic l'
!
          open(3,file='l.dat',status='unknown')
!
!
  call DATE_AND_TIME (date,timc)
  Write(3,*) '# File l.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  If(i_loading==1) write(3,*) '# These l numbers are of loading type'
  If(i_loading==0) write(3,*) '# These l numbers are of _tidal_ type'  	     
!	  
          DO l=lmin, lmax 
          write (3, '((i3,1x,96(1x,e20.8)))') &
	  l, l_e(l), l_f(l), (l_v (l,k), k = 1, nroots)      
!
          ENDDO 
        ENDIF; close(3)  
!
        IF(IK==1) THEN
!
            WRITE(99,*) 'Writing on k.dat the elastic, fluid, and v-elastic k'
  IF(IV==1) WRITE(* ,*) 'Writing on k.dat the elastic, fluid, and v-elastic k'
!
          open(3,file='k.dat',status='unknown') 
!
!
  call DATE_AND_TIME (date,timc)
  Write(3,*) '# File k.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  If(i_loading==1) write(3,*) '# These k numbers are of loading type'
  If(i_loading==0) write(3,*) '# These k numbers are of _tidal_ type'  
!	  
          DO l=lmin, lmax 
          write (3, '((i3,1x,96(1x,e20.8)))') &
	  l, k_e(l), k_f(l), (k_v (l,k), k = 1, nroots)      
!
          ENDDO 
        ENDIF; close(3)  
!
!============================================
        Endif  ! Endif on KW El_Fluid_Viscel
!============================================
!
!
!
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       If(KEYWORD(1:12)=='Heaviside_th')                                   Then
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
                WRITE(99,'(a16,2x,a30 )') '> found KEYWORD', keyword
    IF(IV==1)   WRITE(*, '(a16,2x,a30 )') '> found KEYWORD', keyword
!
Write(99,*) 'Love numbers as a function of time for an Heaviside time-history'
IF(IV==1) &
Write(* ,*) 'Love numbers as a function of time for an Heaviside time-history'
!
!
        READ(1,*) NUMD
!
!
   If(numd>numd_max.and.IV==0) then
   Write(99,*) 'ERROR IN SBR. TASK_1: The maximum number of degrees for '
   Write(99,*) 'the Heaviside Love numbers can be computed is NUMD = 6  '
   Write(99,*) 'Check the input file  task_1.dat  and submit a new task '
   Write(99,*) '******* JOB ABORTED *********************************** ';STOP
   Endif
   If(numd>numd_max.and.IV==1) then
   Write(99,*) 'ERROR IN SBR. TASK_1: The maximum number of degrees for '
   Write(99,*) 'the Heaviside Love numbers can be computed is NUMD = 6  '
   Write(99,*) 'Check the input file  task_1.dat  and submit a new task '
   Write(99,*) '******* JOB ABORTED *********************************** '
   Write(*,*)  'ERROR IN SBR. TASK_1: The maximum number of degrees for '
   Write(*,*)  'the Heaviside Love numbers can be computed is NUMD = 6  '
   Write(*,*)  'Check the input file  task_1.dat  and submit a new task '
   Write(*,*)  '******* JOB ABORTED *********************************** ';STOP
   Endif
!
!
!
        READ(1,*) (D(K),K=1,NUMD) 
!
        READ(1,*) T_MIN
!	
	READ(1,*) T_MAX
!	
	READ(1,*) NPUN 
!	
!
!
   IF(IV==0) THEN
   If((t_min<=0.q0.or.t_max<=0) .OR.  &
      (t_min==t_max)              .OR.  &
      (t_max<=t_min))  THEN
   Write(99,*)'ERROR IN SBR. TASK_1: One of the forbidden conditions  '
   Write(99,*)'t_min <=0; t_max <=0; t_min=t_max; t_max <= t_min  are '
   Write(99,*)'met. Check the input file task_1.dat and try again     '
   Write(99,*)'******* JOB ABORTED ********************************** ';STOP
   Endif
   ENDIF
!
   IF(IV==1) THEN
   If((t_min<=0.q0.or.t_max<=0) .OR.  &
      (t_min==t_max)              .OR.  &
      (t_max<=t_min))  THEN
   Write(99,*)'ERROR IN SBR. TASK_1: One of the forbidden conditions  '
   Write(99,*)'t_min <=0; t_max <=0; t_min=t_max; t_max <= t_min  are '
   Write(99,*)'met. Check the input file task_1.dat and try again     '
   Write(99,*)'******* JOB ABORTED ********************************** '
   Write(*,*) 'ERROR IN SBR. TASK_1: One of the forbidden conditions  '
   Write(*,*) 't_min <=0; t_max <=0; t_min=t_max; t_max <= t_min  are '
   Write(*,*) 'met. Check the input file task_1.dat and try again     '
   Write(*,*) '******* JOB ABORTED ********************************** ';STOP
   Endif    
   endif
!
!
!
!
!
   If(npun>npun_max.and.IV==0) then
   Write(99,*) 'ERROR IN SBR. TASK_1: The number of points between  '
   Write(99,*) 't_min and t_max must not excede ', npun_max
   WRITE(99,*) '***** JOB ABORTED ********************************* '; stop
   Endif
   If(npun>npun_max.and.IV==1) then
   Write(99,*) 'ERROR IN SBR. TASK_1: The number of points between  '
   Write(99,*) 't_min and t_max must not excede ', npun_max
   WRITE(99,*) '***** JOB ABORTED ********************************* '
   Write(*,*)  'ERROR IN SBR. TASK_1: The number of points between  '
   Write(*,*)  't_min and t_max must not excede ', npun_max
   WRITE(*,*)  '***** JOB ABORTED ********************************* '; stop
   ENDIF
!
!
IF(IV==0) THEN
do k=1, numd
if(d(k)<lmin.or.d(k)>lmax) then
Write(99,*) 'ERROR IN SBR. TASK_1: The degree chosen for the Heaviside computation '
Write(99,*) 'is out of range. Each degree must satisfy Lmin <= Degree <= Lmax      '
Write(99,*) 'Check the Input file task_1.dat and submit a new job                  '
Write(99,*) 'Here I have: '
Write(99,*) 'Degree   =   ', d(k)
Write(99,*) 'Lmin     =   ', lmin
Write(99,*) 'Lmax     =   ', lmax
Write(99,*) '********* JOB ABORTED ************* ';STOP
Endif
Enddo
endif
IF(IV==1) then
do k=1, numd
if(d(k)<lmin.or.d(k)>lmax) then
Write(99,*) 'ERROR IN SBR. TASK_1: The degree chosen for the Heaviside computation '
Write(99,*) 'is out of range. Each degree must satisfy Lmin <= Degree <= Lmax      '
Write(99,*) 'Check the Input file task_1.dat and submit a new job                  '
Write(99,*) 'Here I have: '
Write(99,*) 'Degree   =   ', d(k)
Write(99,*) 'Lmin     =   ', lmin
Write(99,*) 'Lmax     =   ', lmax
Write(99,*) '********* JOB ABORTED ************* '
Write(*,*) 'ERROR IN SBR. TASK_1: The degree chosen for the Heaviside computation '
Write(*,*) 'is out of range. Each degree must satisfy Lmin <= Degree <= Lmax      '
Write(*,*) 'Check the Input file task_1.dat and submit a new job                  '
Write(*,*) 'Here I have: '
Write(*,*) 'Degree   =   ', d(k)
Write(*,*) 'Lmin     =   ', lmin
Write(*,*) 'Lmax     =   ', lmax
Write(*,*) '********* JOB ABORTED ************* '; STOP
Endif
Enddo
!
!
endif
!
!
        read(1,*) ih
	read(1,*) il
        read(1,*) ik         
!
        If(ih==0.and.il==0.and.ik==0.and.IV==0) Then
	Write(99,*) 'WARNING IN SBR. TASK_1: No heaviside load-deformation'
	Write(99,*) 'will be printed: ih=il=ik =0 in input file task_1.dat' 
	Endif
        If(ih==0.and.il==0.and.ik==0.and.IV==1) Then
	Write(*,*) 'WARNING IN SBR. TASK_1: No heaviside load-deformation'
	Write(*,*) 'will be printed: ih=il=ik =0 in input file task_1.dat'
	Endif
!
      if(ik==1.and.iv==0.and.iexte==1) then 
        write(99,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
        write(99,*) 'WARNING in sbr TASK1: When the kw External_Model is active, '
	write(99,*) 'the k ldcs are zero by definition. I set IK==0 and continue '
        write(99,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
	IK=0
      endif 	 
      if(ik==1.and.iv==1.and.iexte==1) then 
        write(* ,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '      
        write(* ,*) 'WARNING in sbr TASK1: When the kw External_Model is active, '
	write(* ,*) 'the k ldcs are zero by definition. I set IK==0 and continue '
        write(* ,*) '+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + '
	IK=0
      endif 	 
!
        P_MIN=LOG10(T_MIN); P_MAX=LOG10(T_MAX)
!
        ALFA = (P_MAX-P_MIN)/(DFLOAT(NPUN)-1._DP); BETA = P_MIN-ALFA  
!
!
        if(ih==1)              open(13,file='h_heav.dat',status='unknown')
        if(il==1)              open(14,file='l_heav.dat',status='unknown')
        if(ik==1)              open(15,file='k_heav.dat',status='unknown')
!
!
  call DATE_AND_TIME (date,timc)

  if(ih==1) Write(13,*) '# File h_heav.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)
  If(i_loading==1) write(13,*) '# These h numbers are of loading type'
  If(i_loading==0) write(13,*) '# These h numbers are of _tidal_ type'  
!  	 
  if(il==1) Write(14,*) '# File l_heav.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  If(i_loading==1) write(14,*) '# These l numbers are of loading type'
  If(i_loading==0) write(14,*) '# These l numbers are of _tidal_ type'  
!	     	     
  if(ik==1) Write(15,*) '# File k_heav.dat, created by Task#1 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  If(i_loading==1) write(15,*) '# These k numbers are of loading type'
  If(i_loading==0) write(15,*) '# These k numbers are of _tidal_ type'  
!
!
        IF(ih==1.and.IV==1) WRITE(*,*)  &
                            'Computing the h heaviside ldc (file h_heav.dat)'
        IF(il==1.and.IV==1) WRITE(*,*)  &
                            'Computing the l heaviside ldc (file l_heav.dat)'
        IF(ik==1.and.IV==1) WRITE(*,*)  &
                            'Computing the k heaviside ldc (file k_heav.dat)'
!
        IF(ih==1.and.IV==0) WRITE(99,*)  &
                            'Computing the h heaviside ldc (file h_heav.dat)'
        IF(il==1.and.IV==0) WRITE(99,*)  &
                            'Computing the l heaviside ldc (file l_heav.dat)'
        IF(ik==1.and.IV==0) WRITE(99,*)  &
                            'Computing the k heaviside ldc (file k_heav.dat)'
!
          DO K =1, NUMD
!
          if(ih==1) WRITE(13,'(/)')
          if(il==1) WRITE(14,'(/)')
          if(ik==1) WRITE(15,'(/)')
!
!
              DO I =1, NPUN
! 
              TIME=10._DP**(ALFA*I+BETA)
!
              if (ih==1)               RESPH = h_f(D(K))
              if (il==1)               RESPL = l_f(D(K))
              if (ik==1)               RESPK = k_f(D(K))
!
                  DO M =1, NROOTS
!
if (ih==1)              RESPH = RESPH + (H_V(D(K),M)/S(D(K),M))*EXP( S(D(K),M)*TIME )  
if (il==1)              RESPL = RESPL + (L_V(D(K),M)/S(D(K),M))*EXP( S(D(K),M)*TIME )  
if (ik==1)              RESPK = RESPK + (K_V(D(K),M)/S(D(K),M))*EXP( S(D(K),M)*TIME )  
!
                  ENDDO
!
if (ih==1)              write(13, '(i4,1x,2(e15.7,1x))' ) d(k), time, resph/(2._dp*D(k)+1._dp)
if (il==1)              write(14, '(i4,1x,2(e15.7,1x))' ) d(k), time, respl
if (ik==1)              write(15, '(i4,1x,2(e15.7,1x))' ) d(k), time, respk
!
             ENDDO
!
             ENDDO
!
!
if (ih==1) CLOSE(13)
if (il==1) CLOSE(14)
if (ik==1) CLOSE(15)
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
           Endif ! Endif on KW Heaviside_th
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!
!
!
101 CONTINUE     ! End of the do-loop on the input file task_1.dat
!
102 CONTINUE
!
!

!
         CLOSE(1) 
!
!
        call DATE_AND_TIME (date,timc)
            IF(IV==0) THEN
     Write(99,*) '# Task#1 of TABOO closed on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
        ELSEIF(IV==1) then
     Write(99,*) '# Task#1 of TABOO closed on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
        ENDIF
!
         IF(IV==1) Write(*,*) 'For more details, see file taboo.log '	
!
!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
         END SUBROUTINE TASK_1
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!
!
!
!
!
!
!
!
!
!
!
subroutine SPECTRUM (i_loading)
!
use COMMON_FOR_SPECTRUM; implicit NONE
!
!
!       The sbr SPECTRUM computes the relaxation spectrum, as well as
! the ldcs (or the tLns) of the selected Earth model.
!
!
!
! # Local Variables and their meaning in short
!
!
integer (i4b) :: nl        ! nl = number of radii is always equal to nv+1
!
INTEGER (i4b) :: i_loading ! loading/tidal switch (0/1)
!
integer (i4b) :: indx, inu, i, ii, j, k, kk, l, ll, m ! do-loop control variables
!
!
!
!
nl = nv+1 
!
!
!
!
! # Normalization
!
 call DEFPA  
!
!
!
! # Loop on the harmonic degrees
!
do 1000 l = lmin, lmax
!

!
! # Products of fundamental matrices & inverses
!
do 111 k = 1, nl - 1  
!
call MATPROD  &  
                          (l, r (nl - k), r (nl - k - 1), & 
   rho(nl - k), aco (nl - k), g (nl - k), g (nl - k - 1)) 
!
k0 = ac + bd + bc * rmu (nl - k) + ad / rmu (nl - k)
k1 = ad / vis (nl - k)  
k2 = bc * ( - rmu (nl - k) **2 / vis (nl - k) )  
!
cc (0, k, :, :) = k1 * rmu (nl - k) / vis (nl - k)  
cc (1, k, :, :) = k0 * rmu (nl - k) / vis (nl - k) + k1 + k2
cc (2, k, :, :) = k0  
!
111 CONTINUE   
!

!
! # Elastic product
!
call diretta (l, r (nl),     rho (nl), aco (nl), g (nl),     a, b)  
call inversa (l, r (nl - 1), rho (nl), aco (nl), g (nl - 1), c, d)
!
matela = matmul (a + rmu (nl) * b, c + d / rmu (nl) )  
!
!
!Pj1 projector 
	pj1(:,:)=0.q0 
	pj2(:,:)=0.q0 
!
if(l==1) then 
	pj1(1,1)=1.q0
	pj1(2,2)=1.q0
	pj1(3,5)=1.q0
!
	pj2(1,3)=1.q0
	pj2(2,4)=1.q0
	pj2(3,5)=1.q0
else 
	pj1(1,1)=1.q0
	pj1(2,2)=1.q0
	pj1(3,5)=1.q0
!
	pj2(1,3)=1.q0
	pj2(2,4)=1.q0
	pj2(3,6)=1.q0
endif
!
sinist_2 = matmul (pj2, matela); sinist_1 = matmul (pj1, matela)  
!
! 
! # returns coefmat
! 
    call promat  
!
! # left products
!
    do inu = 0, 2 * nv  
     rrrr (inu, :, :) = matmul (sinist_2, coefmat (inu, :, :) )  
     qqqq (inu, :, :) = matmul (sinist_1, coefmat (inu, :, :) )  
    enddo  
!
! # Right products
!
    call corebo (l, r (0), rho (0), aco (0), cmb)  
!
    do inu = 0, 2 * nv  
     rr (inu, :, :) = matmul (rrrr (inu, :, :), cmb)  
     qq (inu, :, :) = matmul (qqqq (inu, :, :), cmb)  
    enddo  
!
!
! # Computes the determinant
!
    call det_tu
!
!
    do i = 1, nroots + 1  
    indx = 6 * nv + 1  
!
!
! # Coefficients of the secular polynomial
!
    aa (i) = co (indx - i)  
    enddo  
    pivo = aa (1)  
!
! # Roots of the secular polynomial: rt1 = REAL part, rt2 = IMAGINARY part
!
    call zrhqr (aa, nroots, rt1, rt2)  
!
!
!
! # VEC==0 for non physical modes, VEC==1 for 'good' modes
! 
    do i = 1, nroots 
    if(abs(rt2(i)) >= 1.q-30 .or. rt1(i) >= 0._qp) then
      vec(l,i)= 0._sp  
      else 
      vec(l,i)= 1._sp  
      endif 
    enddo  
!
!
!
do kk=1,nroots
!
!
if(abs(rt2(kk)) >= 1.q-30) then 
!
write(99,*) 'WARNING from sbr spectrum: A root with abs(imaginary part) <='
WRITE(99,*) '1E-30 kyrs**(-1)  has been found  in the relaxation spectrum '  
WRITE(99,*) 'degree = ', l, 'root # = ',kk
WRITE(99,*) 'R.P.= ', rt1(kk)
WRITE(99,*) 'I.P.= ', rt2(kk) 
!
if(iv==1) Write(*,*) 'There are''non physical modes''Check file taboo.log'
!
                   endif
!
if(rt1(kk) >= 0._qp) then
!
write(99,*) 'WARNING from sbr spectrum: A root with real part >= 0'
WRITE(99,*) 'has been found in the relaxation spectrum ********** '
WRITE(99,*) 'degree = ', l, 'root # = ',kk
WRITE(99,*) 'R.P.= ', rt1(kk)
WRITE(99,*) 'I.P.= ', rt2(kk) 
!
if(iv==1) Write(*,*) 'There are''non physical modes''Check file taboo.log'
!
                   endif  
!
!
!
!
enddo
!
!                                
s (l,0) = - 1.q30                     ! pole in -oo
do k = 1, nroots
s (l,k) = 1._qp / rt1 (k)             ! s is in kyrs^(-1) ( s<0 )
rad (l, k) = - 1000._qp / s (l,k)     ! rad is in years
enddo  
!
!
! 
! Arrays R and Q
!
do m = 0, nroots  
r_r (:, :, m) = rr (2 * nv, :, :)  
q_q (:, :, m) = qq (2 * nv, :, :)  
!	
do kk = 2 * nv - 1, 0, - 1  
r_r (:, :, m) = r_r (:, :, m) * s (l,m) + rr (kk, :, :)  
q_q (:, :, m) = q_q (:, :, m) * s (l,m) + qq (kk, :, :)  
enddo  
!	  	
enddo  
!
!
! Adjoint of R
!
do m = 0, nroots  
raggiu (1, 1, m) =   ded (r_r (2, 2, m), r_r (2, 3, m), r_r (3, 2, m), r_r (3, 3, m) )
raggiu (2, 1, m) = - ded (r_r (2, 1, m), r_r (2, 3, m), r_r (3, 1, m), r_r (3, 3, m) )
raggiu (3, 1, m) =   ded (r_r (2, 1, m), r_r (2, 2, m), r_r (3, 1, m), r_r (3, 2, m) )
raggiu (1, 2, m) = - ded (r_r (1, 2, m), r_r (1, 3, m), r_r (3, 2, m), r_r (3, 3, m) )
raggiu (2, 2, m) =   ded (r_r (1, 1, m), r_r (1, 3, m), r_r (3, 1, m), r_r (3, 3, m) )
raggiu (3, 2, m) = - ded (r_r (1, 1, m), r_r (1, 2, m), r_r (3, 1, m), r_r (3, 2, m) )
raggiu (1, 3, m) =   ded (r_r (1, 2, m), r_r (1, 3, m), r_r (2, 2, m), r_r (2, 3, m) )
raggiu (2, 3, m) = - ded (r_r (1, 1, m), r_r (1, 3, m), r_r (2, 1, m), r_r (2, 3, m) )
raggiu (3, 3, m) =   ded (r_r (1, 1, m), r_r (1, 2, m), r_r (2, 1, m), r_r (2, 2, m) )
enddo  
!
!
! Product QR
!
do m = 0, nroots  
qr (:, :, m) = matmul (q_q (:, :, m), raggiu (:, :, m) )  
enddo  
!
!
!
! # Surface boundary conditions
!
!   i_loading =1 --> Loading
!   i_loading =0 --> Tidal
!
call bf (l, i_loading)
!
!
!
!
!
! # Secular polynomial at s-> -oo
!
rubens = co (3 * 2 * nv)  
do k = 3 * 2 * nv - 1, 0, - 1  
rubens = rubens * s (l,0) + co (k)  
enddo  
!
! # Elastic part of the solution
!
x_el = matmul (qr (:, :, 0), bcs)  
x_el = x_el / rubens  
!
!
! # Derivative of the secular polynomial
!
do ii = 0, 3 * 2 * nv - 1  
ctmp (ii) = co (ii + 1) * float (ii + 1)  
enddo  
!
do m = 1, nroots  
derpo (m) = ctmp (3 * 2 * nv - 1)  
do k = 3 * 2 * nv - 2, 0, - 1  
derpo (m) = derpo (m) * s (l,m) + ctmp (k)  
enddo  
enddo  
!
!
! # Viscoelastic residues
!
do m = 1, nroots  
xx = matmul (qr (:, :, m), bcs)  
xr (1, m) = xx (1) / derpo (m)  
xr (2, m) = xx (2) / derpo (m)  
xr (3, m) = xx (3) / derpo (m)  
enddo  
!
!
!
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\
!   Loading Love number 'l'
! ///////////////////////////
!
l_e (l)   =  x_el (2) * xmass / r (nl)  
do k=1, nroots 
l_v (l,k) =  xr (2, k) * xmass / r (nl)  
enddo 
!
l_f (l) = l_e (l)  
do m = 1, nroots
l_f (l) = l_f (l) - l_v (l,m) / s (l,m)  
enddo
!
!
!
! ///////////////////////////
!   Loading Love number 'h'
! \\\\\\\\\\\\\\\\\\\\\\\\\\\
!
h_e (l)   = x_el (1) * xmass / r (nl)  
do k=1, nroots 
h_v (l,k) = xr (1, k) * xmass / r (nl)  
enddo 
!
h_f (l) = h_e (l) 
do m = 1, nroots
h_f (l) = h_f (l) - h_v (l,m) / s (l,m)  
enddo
!
!
!
! \\\\\\\\\\\\\\\\\\\\\\\\\\\
!   Loading Love number 'k'
! ///////////////////////////
!
k_e (l)    = - 1._qp - (xmass / r (nl) / g (nl) ) * x_el (3)  
do k=1, nroots 
k_v (l,k)  =         - (xmass / r (nl) / g (nl) ) * xr (3, k)  
enddo 
!
k_f (l) = k_e (l) 
do m = 1, nroots
k_f (l) = k_f (l) - k_v (l,m) / s (l,m)  
enddo
!
!
! 
! 
1000 CONTINUE ! enddo on the harmonic degrees
!
!
! 
RETURN
!
end subroutine SPECTRUM  
!
!
!
!
!
!
!
!
!
subroutine defpa  
!
! # Normalizes the variables and constants
!
use  COMMON_FOR_SPECTRUM; implicit none
integer (i4b) ::  j, k 
!
real (qp) :: massa_in (nv+1)
real (qp) :: ra0, rhor, rig, t0 
real (qp), parameter :: haskell = 1.q21   ! HASKELL unit 
!
!
!
    ggg =0.66732q-10    ! Newton constant
!
!
!
!
! * Earth mass
!
xmass = 4._qp * pi * rho (0) * r (0) **3 / 3._qp  
do k = 1, nv+1  
xmass = xmass + (4._qp / 3._qp) * pi * (r (k) **3 - r (k - 1) **3) &
 * rho (k)
enddo  
!
!
!
! * Surface gravity
!
g (nv+1) = ggg * xmass / r (nv+1) / r (nv+1)  
!
!
! * Gravity at the interfaces
!
g (0) = ggg * (4._qp / 3._qp) * pi * rho (0) * r (0)  
do k = 1, nv+1
massa_in (k) = g (0) * r (0) * r (0) / ggg  
do j = 1, k  
massa_in (k) = massa_in (k) + (4._qp / 3._qp) * pi * (r (j) **3 - &
 r (j - 1) **3) * rho (j)
enddo  
g (k) = massa_in (k) * ggg / r (k) / r (k)  
enddo  
!	
!
! * Real gravity at the interfaces
!
    pon(:) = g(:)  
!
!
 Write(99,*) '-----------------------------------'
 Write(99,*) ' Gravity at the interfaces (m/s/s) '
 Write(99,*) '       (from bottom to top)        '
 Write(99,*) '-----------------------------------'
 Do k=0, nv+1 
       Write(99,'(i4,2x,E14.8)') k, pon(k) 
 Enddo 
!
!	 
!
!
! * Normalization scheme
!
!
ra0 = r (nv+1)    ! reference radius
!
rhor = rho (1)    ! reference density
!
rig = rmu (nv+1)  ! reference shear modulus
!
t0 = 1000._qp * 365.25_qp * 24._qp * 3600._qp  ! seconds in 1 kyr
!
!
! * Normalized parameters
!
do k = 0, nv+1  
    r (k) =   r (k) / ra0
  rho (k) = rho (k) / rhor
  rmu (k) = rmu (k) / rig
  vis (k) = vis (k) / rig / t0
enddo  
!
! * Normalized ravity constant
!
ggg = ggg * rhor**2 * ra0**2 / rig  
!
!
! * Normalized gravity at the interfaces
!
g (0) = ggg * (4._qp / 3._qp) * pi * rho (0) * r (0)  
do k = 1, nv+1  
massa_in (k) = g (0) * r (0) * r (0) / ggg  
do j = 1, k  
massa_in (k) = massa_in (k) + (4._qp / 3._qp) * pi * (r (j) **3 - &
 r (j - 1) **3) * rho (j)
enddo  
g (k) = massa_in (k) * ggg / r (k) / r (k)  
enddo  
!
! * 'A' constants
!
do k = 0, nv+1  
aco (k) = (4._qp * pi / 3._qp) * rho (k) * ggg  
enddo  
!
! * Normalized mass of the Earth
!
xmass = 4._qp * pi * rho (0) * r (0) **3 / 3._qp  
do k = 1, nv+1  
xmass = xmass + (4._qp / 3._qp) * pi * (r (k) **3 - r (k - 1) **3) &
 * rho (k)
enddo  
!
end subroutine DEFPA 
!
!
!
!
!
!
subroutine COREBO (n, r, ro, a, g)
!
! # Boundary conditions at the CMB
!
use strata 
implicit none 
integer (i4b) :: n
real (qp) :: g (6, 3)  
real (qp) :: r, ro, a 
!
g = 0._qp  
g (1, 1) = - (r** (n - 1) ) / a  
g (1, 3) = 1._qp  
g (2, 2) = 1._qp  
g (3, 3) = ro * a * r  
g (5, 1) = r**n  
g (6, 1) = 2._qp * (float (n) - 1._qp) * r** (n - 1)  
g (6, 3) = 3._qp * a  
!
end subroutine COREBO 
!
!
!
!
!
!
subroutine det_tu
!
! # Secular determinant
!
use  COMMON_FOR_SPECTRUM
implicit none  
integer (i4b) :: inex, is, k, l, m, mu, nu
!
!
do inex = 0, 3 * 2 * nv  
co (inex) = 0._qp  
    do k = 0, 2 * nv  
        do l = 0, 2 * nv  
            do m = 0, 2 * nv  
                if (k + l + m == inex) then  
                   co (inex) = co (inex) + &
		   rr (k, 1, 1) * rr (l, 2, 2) * rr (m, 3, 3) - &
		   rr (k, 1, 1) * rr (l, 2, 3) * rr (m, 3, 2) - &
		   rr (k, 1, 2) * rr (l, 2, 1) * rr (m, 3, 3) + &
		   rr (k, 1, 2) * rr (l, 2, 3) * rr (m, 3, 1) + &
		   rr (k, 1, 3) * rr (l, 2, 1) * rr (m, 3, 2) - &
		   rr (k, 1, 3) * rr (l, 2, 2) * rr (m, 3, 1)
                 endif  
            enddo  
        enddo  
    enddo  
enddo
end subroutine det_tu
!
!
!
!
!
!
function ded (a, b, c, d)
!
! # Determinant od a 2x2 table
!
use strata
implicit none  
real (qp) :: a, b, c, d, ded  
ded = a * d-b * c  
end function ded
!
!
!
!
!
!
subroutine bf (l,char)  
!
! # Surface boundary conditions
!
use COMMON_FOR_SPECTRUM 
implicit none  
integer (i4b) :: char, l  
!
!
real (qp) :: gamma  
gamma = (2._qp * float (l) + 1._qp) / 4._qp / pi / r (nv+1) / r (nv+1)  
!
if (char == 1) then  
    bcs (1) = - g (nv+1) * gamma  
    bcs (2) = 0._qp  
	if(l.ne.1) bcs (3) = - 4._qp * pi * ggg * gamma  
	if(l.eq.1) bcs (3) = 0._qp  
endif  
!
if (char == 0) then
    bcs (1) = 0._qp  
    bcs (2) = 0._qp  
    bcs (3) = -4._qp * pi * ggg * gamma  
	if(l.eq.1) STOP '*** ERROR: Tidal conditions are not defined at degree 1'
endif  
end subroutine bf
!
!
!
!
!
!
subroutine inversa (n, r, rho, za, gra, a, b)  
!
! # Inverse of the fundamental matrix
!
use strata  
implicit none
integer (i4b) :: i, j, n  
real (qp) :: gra, r, rho, za
real (qp) :: a (6, 6), b (6, 6)  
real (qp) :: v(24)
real (qp) :: h1, h2, h3, h4, h5
real (qp) :: rn
!
!
a(:,:) = 0._qp 
b(:,:) = 0._qp 
!
!
rn = float (n)  
!
h1 = rn + 1._qp  
h2 = 1._qp + 2._qp * rn  
h3 = 2._qp * rn - 1._qp  
h4 = 3._qp + 2._qp * rn  
h5 = rn - 1._qp  
!
v (1) = h1 / h2  
v (2) = - 2._qp * h1 * (rn + 2._qp) / h2  
v (3) = rn * h1 / 2._qp / h3 / h2  
v (4) = - (rn**2 + 3._qp * rn - 1._qp) * rn / h2 / h3  
v (5) = rn / h2  
v (6) = 2._qp * rn * h5 / h2  
v (7) = rn * h1 / 2._qp / h2 / h4  
v (8) = h1 * (rn**2 - rn - 3._qp) / h2 / h4  
v (9) = - 2._qp * rn * h1 * (2._qp + rn) / h2  
v (10) = - h5 * rn * h1 * h1 / h3 / h2  
v (11) = - 2._qp * h5 * rn * h1 / h2  
v (12) = - rn**2 * h1 * (2._qp + rn) / h2 / h4  
v (13) = - v (1)  
v (14) = - v (3)  
v (15) = - v (5)  
v (16) = - v (7)  
v (17) = - rn * h1 / h2  
v (18) = (2._qp - rn) * h1 * rn / 2._qp / h3 / h2  
v (19) = rn * h1 / h2  
v (20) = rn * h1 * (3._qp + rn) / 2._qp / h2 / h4  
v (21) = - v (1)  
v (22) = - v (3)  
v (23) = - v (5)  
v (24) = - v (7)  
!
a (1, 1) = r**( - n - 1) * v (2)  
a (2, 1) = - r**( - n + 1) * v (4)  
a (3, 1) = 3._qp * za * r**( - n + 1) / (2._qp * float (n) + 1._qp)
a (4, 1) = r**(n) * v (6)  
a (5, 1) = - r**(n + 2) * v (8)  
a (6, 1) = - 3._qp * za * r**(n + 2) / (2._qp * float (n) + 1._qp)  
!
a (1, 2) = - r**( - n - 1) * v (9)  
a (2, 2) = r**( - n + 1) * v (10)  
a (4, 2) = - r**(n) * v (11)  
a (5, 2) = r**(n + 2) * v (12)  
!
a (6, 5) = - r**(n + 1)  
a (3, 6) = - r**( - n + 1) / (2._qp * float (n) + 1._qp)  
a (6, 6) = r**(n + 2) / (2._qp * float (n) + 1._qp)  
!
b (1, 1) = r**( - n) * rho * gra * v (1)  
b (2, 1) = - r**( - n + 2) * rho * gra * v (3)  
b (4, 1) = r**(n + 1) * rho * gra * v (5)  
b (5, 1) = - r**(n + 3) * rho * gra * v (7)  
!
b (1, 3) = r**( - n) * v (13)  
b (2, 3) = - r**( - n + 2) * v (14)  
b (4, 3) = r**(n + 1) * v (15)  
b (5, 3) = - r**(n + 3) * v (16)  
!
b (1, 4) = - r**( - n) * v (17)  
b (2, 4) = r**( - n + 2) * v (18)  
b (4, 4) = - r**(n + 1) * v (19)  
b (5, 4) = r**(n + 3) * v (20)  
!
b (1, 5) = - rho * r**( - n) * v (21)  
b (2, 5) = rho * r**( - n + 2) * v (22)  
b (4, 5) = - rho * r**(n + 1) * v (23)  
b (5, 5) = rho * r**(n + 3) * v (24)  
!
end subroutine inversa
!
!
!
!
!
!
subroutine diretta (n, r, rho, za, gra, a, b)  
!
! # Direct matrix
!
use strata 
implicit none
integer (i4b) :: i, j, n  
real (qp) :: a (6, 6), b (6, 6)
real (qp) :: gra, r, rho, za
real (qp) :: rn
real (qp) :: a1, a2, a3, a4, b1, b2, b3, b4 
real (qp) :: c1, c2, c3, c4, d1, d2, d3, d4
!
!
a(:,:) = 0._qp 
b(:,:) = 0._qp 
!
!
rn = float (n)  
!
a1 = rn / (2._qp * (2._qp * rn + 3._qp) )  
a2 = 1._qp  
a3 = (rn + 1._qp) / (2._qp * (2._qp * rn - 1._qp) )  
a4 = 1._qp  
!
b1 = (rn + 3._qp) / (2._qp * (2._qp * rn + 3._qp) * (1._qp * rn + 1._qp) )
b2 = 1._qp / rn  
b3 = ( - rn + 2._qp) / (2._qp * rn * (2._qp * rn - 1._qp) )  
b4 = - 1._qp / (rn + 1._qp)  
!
c1 = (rn * rn - rn - 3._qp) / (2._qp * rn + 3._qp)  
c2 = 2._qp * (rn - 1._qp)  
c3 = ( - rn * rn - 3._qp * rn + 1._qp) / (2._qp * rn - 1._qp)  
c4 = - 2._qp * (rn + 2._qp)  
!
d1 = rn * (rn + 2._qp) / ( (2._qp * rn + 3._qp) * (rn + 1._qp) )  
d2 = 2._qp * (rn - 1._qp) / rn  
d3 = (rn * rn - 1._qp) / (rn * (2._qp * rn - 1._qp) )  
d4 = 2._qp * (1._qp * rn + 2._qp) / (rn + 1._qp)  
!
a (1, 1) = a1 * r** (n + 1)  
a (1, 2) = a2 * r** (n - 1)  
a (1, 4) = a3 * r** ( - n)  
a (1, 5) = a4 * r** ( - n - 2)  
a (2, 1) = b1 * r** (n + 1)  
a (2, 2) = b2 * r** (n - 1)  
a (2, 4) = b3 * r** ( - n)  
a (2, 5) = b4 * r** ( - n - 2)  
a (3, 1) = a (1, 1) * rho * gra  
a (3, 2) = a (1, 2) * rho * gra  
a (3, 3) = - rho * r** (n)  
a (3, 4) = a (1, 4) * rho * gra  
a (3, 5) = a (1, 5) * rho * gra  
a (3, 6) = - rho * r** ( - n - 1)  
a (5, 3) = a (3, 3) / rho  
a (5, 6) = a (3, 6) / rho  
a (6, 1) = 3._qp * za * a (1, 1)  
a (6, 2) = 3._qp * za * a (1, 2)  
a (6, 3) = - (2._qp * float (n) + 1._qp) * r** (n - 1)  
a (6, 4) = 3._qp * za * a (1, 4)  
a (6, 5) = 3._qp * za * a (1, 5)  
!
b (3, 1) = c1 * r** (n)  
b (3, 2) = c2 * r** (n - 2)  
b (3, 4) = c3 * r** ( - n - 1)  
b (3, 5) = c4 * r** ( - n - 3)  
!
b (4, 1) = d1 * r** (n)  
b (4, 2) = d2 * r** (n - 2)  
b (4, 4) = d3 * r** ( - n - 1)  
b (4, 5) = d4 * r** ( - n - 3)  
!
end subroutine diretta
!
!
!
!
!
!
subroutine promat  
use COMMON_FOR_SPECTRUM 
implicit none
integer (i4b) :: i, j, k, l, n
!
! # Product od propagators
!
real (qp) :: ai (0:2 * nv, 6, 6), bi (0:2, 6, 6), ci (0:2 * nv, 6, 6)  
ai (0:2, :, :) = cc (:, 1, :, :)  
if (nv/=1) then  
    do i = 2, nv  
    n = 2 * i  
        do k = 0, n  
        ci (k, :, :) = 0._qp  
        bi (:, :, :) = cc (:, i, :, :)  
            do l = 0, n - 2  
                do j = 0, 2  
                    if (l + j == k) then  
                    ci (k, :, :) = ci (k, :, :) + &
                          matmul (ai (l, :, :), bi (j, :, :) )
                    endif  
                enddo  
            enddo  
        enddo  
    ai (0:n, :, :) = ci (0:n, :, :)  
    enddo  
    coefmat = ci  
else  
    coefmat (0:2, :, :) = ai (0:2, :, :)  
endif  
end subroutine promat
!
!
!
!
!
!
subroutine balanc (a, n, np)
use strata  
implicit none  
integer :: n, np  
real (qp) :: a (np, np), radix, sqrdx  
parameter (radix = 2._qp, sqrdx = radix**2)  
integer :: i, j, last  
real (qp) :: c, f, g, r, s  
    1 continue  
last = 1  
do i = 1, n  
    c = 0._qp  
    r = 0._qp  
    do j = 1, n  
        if (j/=i) then  
            c = c + abs (a (j, i) )  
            r = r + abs (a (i, j) )  
        endif  
    end do  
    if (c/=0._qp.and.r/=0._qp) then  
        g = r / radix  
        f = 1._qp  
        s = c + r  
    2         if (c < g) then  
            f = f * radix  
            c = c * sqrdx  
            goto 2  
        endif  
        g = r * radix  
    3         if (c > g) then  
            f = f / radix  
            c = c / sqrdx  
            goto 3  
        endif  
        if ( (c + r) / f < 0.95_qp * s) then  
            last = 0  
            g = 1._qp / f  
            do j = 1, n  
                a (i, j) = a (i, j) * g  
            end do  
            do j = 1, n  
                a (j, i) = a (j, i) * f  
            end do  
        endif  
    endif  
end do  
if (last == 0) goto 1  
return  
end subroutine balanc
!  (c) copr. 1986-92 numerical recipes software eslp3-1,1')!.
!
!
!
!
!
!
subroutine hqr (a, n, np, wr, wi)  
use strata
implicit none  
integer :: n, np  
real (qp) :: a (np, np), wi (np), wr (np)  
integer :: i, its, j, k, l, m, nn  
real (qp) :: anorm, p, q, r, s, t, u, v, w, x, y, z  
anorm = abs (a (1, 1) )  
do i = 2, n  
    do j = i - 1, n  
        anorm = anorm + abs (a (i, j) )  
    end do  
end do  
nn = n  
t = 0._qp  
    1 if (nn >= 1) then  
    its = 0  
2     do l = nn, 2, - 1  
        s = abs (a (l - 1, l - 1) ) + abs (a (l, l) )  
        if (s == 0._qp) s = anorm  
        if (abs (a (l, l - 1) ) + s == s) goto 3  
      end do  
    l = 1  
3     x = a (nn, nn)  
    if (l == nn) then  
        wr (nn) = x + t  
        wi (nn) = 0._qp  
        nn = nn - 1  
    else  
        y = a (nn - 1, nn - 1)  
        w = a (nn, nn - 1) * a (nn - 1, nn)  
        if (l == nn - 1) then  
            p = 0.5_qp * (y - x)  
            q = p**2 + w  
            z = sqrt (abs (q) )  
            x = x + t  
            if (q >= 0._qp) then  
                z = p + sign (z, p)  
                wr (nn) = x + z  
                wr (nn - 1) = wr (nn)  
                if (z/=0._qp) wr (nn) = x - w / z  
                wi (nn) = 0._qp  
                wi (nn - 1) = 0._qp  
            else  
                wr (nn) = x + p  
                wr (nn - 1) = wr (nn)  
                wi (nn) = z  
                wi (nn - 1) = - z  
            endif  
            nn = nn - 2  
        else  
            if (its == 30) pause 'too many iterations in hqr'  
            if (its == 10.or.its == 20) then  
                t = t + x  
                do i = 1, nn  
                    a (i, i) = a (i, i) - x  
                end do  
                s = abs (a (nn, nn - 1) ) + abs (a (nn - 1, nn - &
                 2) )
                x = 0.75_qp * s  
                y = x  
                w = - 0.4375_qp * s**2  
            endif  
            its = its + 1  
            do m = nn - 2, l, - 1  
                z = a (m, m)  
                r = x - z  
                s = y - z  
                p = (r * s - w) / a (m + 1, m) + a (m, m + 1)  
                q = a (m + 1, m + 1) - z - r - s  
                r = a (m + 2, m + 1)  
                s = abs (p) + abs (q) + abs (r)  
                p = p / s  
                q = q / s  
                r = r / s  
                if (m == l) goto 4  
                u = abs (a (m, m - 1) ) * (abs (q) + abs (r) )  
                v = abs (p) * (abs (a (m - 1, m - 1) ) + abs (z) &
                 + abs (a (m + 1, m + 1) ) )
                if (u + v == v) goto 4  
             end do  
    4             do i = m + 2, nn  
                     a (i, i - 2) = 0._qp  
                     if (i/=m + 2) a (i, i - 3) = 0._qp  
                  end do  
            do k = m, nn - 1  
                if (k/=m) then  
                    p = a (k, k - 1)  
                    q = a (k + 1, k - 1)  
                    r = 0._qp  
                    if (k/=nn - 1) r = a (k + 2, k - 1)  
                    x = abs (p) + abs (q) + abs (r)  
                    if (x/=0._qp) then  
                        p = p / x  
                        q = q / x  
                        r = r / x  
                    endif  
                endif  
                s = sign (sqrt (p**2 + q**2 + r**2), p)  
                if (s/=0._qp) then  
                    if (k == m) then  
                        if (l/=m) a (k, k - 1) = - a (k, k - 1)  
                    else  
                        a (k, k - 1) = - s * x  
                    endif  
                    p = p + s  
                    x = p / s  
                    y = q / s  
                    z = r / s  
                    q = q / p  
                    r = r / p  
                    do j = k, nn  
                        p = a (k, j) + q * a (k + 1, j)  
                        if (k/=nn - 1) then  
                            p = p + r * a (k + 2, j)  
                            a (k + 2, j) = a (k + 2, j) - p * z  
                        endif  
                        a (k + 1, j) = a (k + 1, j) - p * y  
                        a (k, j) = a (k, j) - p * x  
                    end do  
                    do i = l, min (nn, k + 3)  
                        p = x * a (i, k) + y * a (i, k + 1)  
                        if (k/=nn - 1) then  
                            p = p + z * a (i, k + 2)  
                            a (i, k + 2) = a (i, k + 2) - p * r  
                        endif  
                        a (i, k + 1) = a (i, k + 1) - p * q  
                        a (i, k) = a (i, k) - p  
                    end do  
                endif  
            end do  
            goto 2  
        endif  
    endif  
    goto 1  
endif  
return  
end subroutine hqr
!  (c) copr. 1986-92 numerical recipes software eslp3-1,1')!.
!
!
!
!
!
subroutine zrhqr (a, m, rtr, rti)
use strata; implicit none  
integer :: m, maxm  
real (qp) :: a (m + 1), rtr (m), rti (m)  
parameter (maxm = 500)  
!uses balanc,hqr
integer :: j, k  
real (qp) :: hess (maxm, maxm), xr, xi  
if (m > maxm.or.a (m + 1)  == 0._qp) pause 'bad args in zrhqr'  
do k = 1, m  
    hess (1, k) = - a (m + 1 - k) / a (m + 1)  
    do j = 2, m  
        hess (j, k) = 0._qp  
    end do  
    if (k/=m) hess (k + 1, k) = 1._qp  
end do  
call balanc (hess, m, maxm)  
call hqr (hess, m, maxm, rtr, rti)  
do j = 2, m  
    xr = rtr (j)  
    xi = rti (j)  
    do k = j - 1, 1, - 1  
        if (rtr (k)  <= xr) goto 1  
        rtr (k + 1) = rtr (k)  
        rti (k + 1) = rti (k)  
    end do  
    k = 0  
1   rtr (k + 1) = xr  
    rti (k + 1) = xi  
end do  
return  
end subroutine zrhqr
!
!
!
!
!
!
subroutine matprod (n, rb, rd, rhu, za, grx, gry) 
!
! # An inverse times a direct fundamental matrix
!
use common_for_spectrum 
implicit none
integer (i4b) :: i, j, n  
real (qp) :: grx, gry, rb, rd, rhu, za
real (qp) :: v(24)
real (qp) :: h1, h2, h3, h4, h5
real (qp) :: rn
real (qp) :: a1, a2, a3, a4, b1, b2, b3, b4 
real (qp) :: c1, c2, c3, c4, d1, d2, d3, d4
real (qp) :: x,y
ac = 0._qp  
ad = 0._qp  
bc = 0._qp  
bd = 0._qp  
!
rn = float (n)  
x=rb
y=rd
!
a1 = rn / (2._qp * (2._qp * rn + 3._qp) )  
a2 = 1._qp  
a3 = (rn + 1._qp) / (2._qp * (2._qp * rn - 1._qp) )  
a4 = 1._qp  
!
b1 = (rn + 3._qp) / (2._qp * (2._qp * rn + 3._qp) * (1._qp * rn + 1._qp) )
b2 = 1._qp / rn  
b3 = ( - rn + 2._qp) / (2._qp * rn * (2._qp * rn - 1._qp) )  
b4 = - 1._qp / (rn + 1._qp)  
!
c1 = (rn * rn - rn - 3._qp) / (2._qp * rn + 3._qp)  
c2 = 2._qp * (rn - 1._qp)  
c3 = ( - rn * rn - 3._qp * rn + 1._qp) / (2._qp * rn - 1._qp)  
c4 = - 2._qp * (rn + 2._qp)  
!
d1 = rn * (rn + 2._qp) / ( (2._qp * rn + 3._qp) * (rn + 1._qp) )  
d2 = 2._qp * (rn - 1._qp) / rn  
d3 = (rn * rn - 1._qp) / (rn * (2._qp * rn - 1._qp) )  
d4 = 2._qp * (1._qp * rn + 2._qp) / (rn + 1._qp)  
!
h1 = rn + 1._qp  
h2 = 1._qp + 2._qp * rn  
h3 = 2._qp * rn - 1._qp  
h4 = 3._qp + 2._qp * rn  
h5 = rn - 1._qp  
!
v (1) = h1 / h2  
v (2) = - 2._qp * h1 * (rn + 2._qp) / h2  
v (3) = rn * h1 / 2._qp / h3 / h2  
v (4) = - (rn**2 + 3._qp * rn - 1._qp) * rn / h2 / h3  
v (5) = rn / h2  
v (6) = 2._qp * rn * h5 / h2  
v (7) = rn * h1 / 2._qp / h2 / h4  
v (8) = h1 * (rn**2 - rn - 3._qp) / h2 / h4  
v (9) = - 2._qp * rn * h1 * (2._qp + rn) / h2  
v (10) = - h5 * rn * h1 * h1 / h3 / h2  
v (11) = - 2._qp * h5 * rn * h1 / h2  
v (12) = - rn**2 * h1 * (2._qp + rn) / h2 / h4  
v (13) = - v (1)  
v (14) = - v (3)  
v (15) = - v (5)  
v (16) = - v (7)  
v (17) = - rn * h1 / h2  
v (18) = (2._qp - rn) * h1 * rn / 2._qp / h3 / h2  
v (19) = rn * h1 / h2  
v (20) = rn * h1 * (3._qp + rn) / 2._qp / h2 / h4  
v (21) = - v (1)  
v (22) = - v (3)  
v (23) = - v (5)  
v (24) = - v (7)  
!
ac (1,1) =  (a1*v(2)*(x/y)**(n+1) -a2*v(4)*(x/y)**(n-1) +a3*v(6)*(y/x)**(n) -a4*v(8)*(y/x)**(n+2))
ac (1,2) = -(a1*v(9)*(x/y)**(n+1)-a2*v(10)*(x/y)**(n-1)+a3*v(11)*(y/x)**(n)-a4*v(12)*(y/x)**(n+2))

ac (2,1) =  (b1*v(2)*(x/y)**(n+1) -b2*v(4)*(x/y)**(n-1) +b3*v(6)*(y/x)**(n) -b4*v(8)*(y/x)**(n+2))
ac (2,2) = -(b1*v(9)*(x/y)**(n+1)-b2*v(10)*(x/y)**(n-1)+b3*v(11)*(y/x)**(n)-b4*v(12)*(y/x)**(n+2))

ac (3,1) =  rhu*grx*ac(1,1)-(3._qp*rhu*za*y/(2._qp*rn+1._qp))*((x/y)**(n)-(y/x)**(n+1))
ac (3,2) =  rhu*grx*ac(1,2)
ac (3,5) =  rhu*(y/x)**(n+1)
ac (3,6) =  (rhu*y/(2._qp*rn+1._qp))*((x/y)**(n)-(y/x)**(n+1))

ac (5,1) = -(3._qp*za*y/(2._qp*rn+1._qp))*((x/y)**(n)-(y/x)**(n+1))
ac (5,5) =  (y/x)**(n+1)
ac (5,6) =  (y/(2._qp*rn+1._qp))*((x/y)**(n)-(y/x)**(n+1))

ac (6,1) =  3._qp*za*(a1*v(2)*(x/y)**(n+1)-(1+a2*v(4))*(x/y)**(n-1)+a3*v(6)*(y/x)**(n) -a4*v(8)*(y/x)**(n+2))
ac (6,2) =  3._qp*za*ac(1,2)
ac (6,6) =  (x/y)**(n-1)

ad (1,1) =  rhu*gry*( a1*v(1)*x*(x/y)**(n) -a2*v(3)*x*(x/y)**(n-2) +a3*v(5)*y*(y/x)**(n) -a4*v(7)*y*(y/x)**(n+2))
ad (1,3) =          (a1*v(13)*x*(x/y)**(n)-a2*v(14)*x*(x/y)**(n-2)+a3*v(15)*y*(y/x)**(n)-a4*v(16)*y*(y/x)**(n+2))
ad (1,4) =         -(a1*v(17)*x*(x/y)**(n)-a2*v(18)*x*(x/y)**(n-2)+a3*v(19)*y*(y/x)**(n)-a4*v(20)*y*(y/x)**(n+2))
ad (1,5) = -rhu*    (a1*v(21)*x*(x/y)**(n)-a2*v(22)*x*(x/y)**(n-2)+a3*v(23)*y*(y/x)**(n)-a4*v(24)*y*(y/x)**(n+2))

ad (2,1) =  rhu*gry*( b1*v(1)*x*(x/y)**(n) -b2*v(3)*x*(x/y)**(n-2) +b3*v(5)*y*(y/x)**(n) -b4*v(7)*y*(y/x)**(n+2))
ad (2,3) =          (b1*v(13)*x*(x/y)**(n)-b2*v(14)*x*(x/y)**(n-2)+b3*v(15)*y*(y/x)**(n)-b4*v(16)*y*(y/x)**(n+2))
ad (2,4) =         -(b1*v(17)*x*(x/y)**(n)-b2*v(18)*x*(x/y)**(n-2)+b3*v(19)*y*(y/x)**(n)-b4*v(20)*y*(y/x)**(n+2))
ad (2,5) = -rhu*    (b1*v(21)*x*(x/y)**(n)-b2*v(22)*x*(x/y)**(n-2)+b3*v(23)*y*(y/x)**(n)-b4*v(24)*y*(y/x)**(n+2))

ad (3,1) =  rhu*grx*ad(1,1)
ad (3,3) =  rhu*grx*ad(1,3)
ad (3,4) =  rhu*grx*ad(1,4)
ad (3,5) =  rhu*grx*ad(1,5)

ad (6,1) =  3._qp*za*ad(1,1)
ad (6,3) =  3._qp*za*ad(1,3)
ad (6,4) =  3._qp*za*ad(1,4)
ad (6,5) =  3._qp*za*ad(1,5)

bc (3,1) =  (c1*v(2)*(1._qp/y)*(x/y)**(n) -c2*v(4)*(1._qp/y)*(x/y)**(n-2) +c3*v(6)*(1._qp/x)*(y/x)**(n) &
            -c4*v(8)*(1._qp/x)*(y/x)**(n+2))
bc (3,2) = -(c1*v(9)*(1._qp/y)*(x/y)**(n)-c2*v(10)*(1._qp/y)*(x/y)**(n-2)+c3*v(11)*(1._qp/x)*(y/x)**(n) &
           -c4*v(12)*(1._qp/x)*(y/x)**(n+2))

bc (4,1) =  (d1*v(2)*(1._qp/y)*(x/y)**(n) -d2*v(4)*(1._qp/y)*(x/y)**(n-2) +d3*v(6)*(1._qp/x)*(y/x)**(n) &
           -d4*v(8)*(1._qp/x)*(y/x)**(n+2))
bc (4,2) = -(d1*v(9)*(1._qp/y)*(x/y)**(n)-d2*v(10)*(1._qp/y)*(x/y)**(n-2)+d3*v(11)*(1._qp/x)*(y/x)**(n) &
           -d4*v(12)*(1._qp/x)*(y/x)**(n+2))


bd (3,1) =  rhu*gry*(c1*v(1)*(x/y)**(n) -c2* v(3)*(x/y)**(n-2) +c3*v(5)*(y/x)**(n+1) -c4*v(7)*(y/x)**(n+3))
bd (3,3) =          (c1*v(13)*(x/y)**(n)-c2*v(14)*(x/y)**(n-2)+c3*v(15)*(y/x)**(n+1)-c4*v(16)*(y/x)**(n+3))
bd (3,4) =         -(c1*v(17)*(x/y)**(n)-c2*v(18)*(x/y)**(n-2)+c3*v(19)*(y/x)**(n+1)-c4*v(20)*(y/x)**(n+3))
bd (3,5) = -rhu*    (c1*v(21)*(x/y)**(n)-c2*v(22)*(x/y)**(n-2)+c3*v(23)*(y/x)**(n+1)-c4*v(24)*(y/x)**(n+3))

bd (4,1) = rhu*gry*(d1*v(1)*(x/y)**(n) -d2* v(3)*(x/y)**(n-2) +d3*v(5)*(y/x)**(n+1) -d4*v(7)*(y/x)**(n+3))
bd (4,3) =         (d1*v(13)*(x/y)**(n)-d2*v(14)*(x/y)**(n-2)+d3*v(15)*(y/x)**(n+1)-d4*v(16)*(y/x)**(n+3))
bd (4,4) =        -(d1*v(17)*(x/y)**(n)-d2*v(18)*(x/y)**(n-2)+d3*v(19)*(y/x)**(n+1)-d4*v(20)*(y/x)**(n+3))
bd (4,5) = -rhu*   (d1*v(21)*(x/y)**(n)-d2*v(22)*(x/y)**(n-2)+d3*v(23)*(y/x)**(n+1)-d4*v(24)*(y/x)**(n+3))

end subroutine matprod
!
!
!
!
!
!
!
!
!
!
!
!
!
!
subroutine prem (radsup, radinf, rigi, dens)  
!
! # This program gives the volume-averaged mean values of the
!   density,  and the lame parameters  mu and  lambda for the
!   prem-models for the reference periods at  1 s  and  200 s
!
!
! Input: 
!
!   radsup = bottom radius (km)
!   radinf = top radius  (km)
!
! Output: 
!
!   rigi = shear modulus of the layer (in units of 10**11 Pa)
!   dens = layer density (kg/m**3)
!
use strata
implicit none  
integer (i4b) :: i, idepth, itop, k, mod
real (qp) :: radsup, radinf, rigi, dens  
real (qp) :: mu (100), kappa (100), rho (100), r (100)  
real (qp) :: mean1 (0:100), mean2 (0:100), mean3 (0:100), depth, top
!
open (19, file = 'prem1.dat'  , status = 'old')  
open (29, file = 'prem200.dat', status = 'old')  
!
top = radsup  
depth = radinf  
!
!
124 format (4(e14.5))  
!
do mod = 29, 29  ! Uses the PREM at 200 s
!
    do i = 1, 94  
    read (mod, * ) r (i), kappa (i), mu (i), rho (i)  
    enddo  
    if (depth < 0._qp .or. depth > 6371._qp) then  
        print * , 'sorry, outside earth'  
        print *, top, depth  
        goto 13  
    endif  
    if (depth > top) then  
        print * , 'sorry, top is smaller than depth'  
        print *, top, depth  
        goto 13  
    endif  
    do i = 1, 94  
        if (r (i) - top >= 0._qp) then  
            itop = i  
            goto 4  
        endif  
    end do  
4   do i = 1, 94  
        if (r (i) - depth >= 0._qp) then  
            idepth = i - 1  
                if ( (r (i)  == r (i + 1) .and.r (i)  == depth) .or.r &
                   (i)  == depth) then
                   idepth = i  
                 endif  
             goto 6  
         endif  
     end do  
6     continue  
    if (mod == 19) then  
!         print*,'reference period of 1 s:'
!       print*,'------------------------'
    else  
!         print*,'reference period of 200 s:'
!       print*,'--------------------------'
    endif  
    if (r (idepth)  == r (itop) ) then  
!         print*,'homogeneous layering'
        mean1 (itop - 1) = mu (itop)  
        mean2 (itop - 1) = kappa (itop) - 2._qp / 3._qp * mu (itop)  
        mean3 (itop - 1) = rho (itop)  
        goto 11  
    endif  
!       print*,'reference top layer =',r(itop),' km'
!       print*,'reference lower layer =',r(idepth),' km'
    mean1 (idepth - 1) = 0._qp  
    mean2 (idepth - 1) = 0._qp  
    mean3 (idepth - 1) = 0._qp  
    do k = idepth, itop - 1  
        mean1 (k) = mean1 (k - 1) + (r (k + 1) **3 - r (k) **3) &
         / (r (itop) **3 - r (idepth) **3) * mu (k + 1)
        mean2 (k) = mean2 (k - 1) + (r (k + 1) **3 - r (k) **3) &
         / (r (itop) **3 - r (idepth) **3) * (kappa (k) - 2._qp / &
         3._qp * mu (k + 1) )
        mean3 (k) = mean3 (k - 1) + (r (k + 1) **3 - r (k) **3) &
         / (r (itop) **3 - r (idepth) **3) * rho (k + 1)
    end do  
!  11   print*,'mean value of mu = ',mean1(itop-1)/1q3,'x 10^11 pa'
!       print*,'mean value of lambda = ',mean2(itop-1)/1q3,'x 10^11 pa'
!       print*,'mean value of rho = ',mean3(itop-1)*1q3,'kg per m^3'
!
!
!
   11     rigi = mean1 (itop - 1) / 1.q3  
    dens = mean3 (itop - 1) * 1.q3  
!
!
!
end do  
close (19)  
close (29)  
13 continue  
end subroutine prem
!
!
!
!
!
!
!
!
!
!-------------------------------------------------------+
!
 SUBROUTINE AXIS_LOAD (Iw, TIPO, AMPLITUDE, H, MASS)
!
!-------------------------------------------------------+
!
! Important notes for the User:
! ****************************
!
!       This routine computes the coefficients SIGMA of the Legendre
! polynomials series of an axissimmetric load, in the reference system
! where the axis of symmetry coincides with the 'z' axis.
!
!       The coefficients are computed in the range [l_min:l_max] set
! by the kw Harmonic_Degrees. For IW=1, they are written on the file
! load_coeff.dat in units of kg/m**2.
!
!       The routine AXIS_LOAD can manage the axissimmetric loads with
! label C_L=10, 11, 20, 21, 30, and 40 (see the chapter concerning
! the kw Load_Geometry in the User guide).
!
!       Please note that:
!
!
! For C_L <= 21: The value of H in input is interpreted as the maximum
!                thickness of the load during the course of its time-
!                evolution (in m). For load with parabolic cross-section
!                (C_L = 20 or 21), H is the thickness at the centre of the
!                load. For any value of C_L, H must be >=0 in input.
!                The mass of the load (MASS) is computed by this routine in
!                units of kg consistently with the value of H, and made
!                available esternally. The input AMPLITUDE represents the
!                half-amplitude of the load (degrees).
!
!                     
! For C_L == 30: The inputs H and AMPLITUDE are ignored and set to zero
!                in output. The mass MASS (kg)  is taken as an input parameter
!                and must be > 0.
!
!
! For C_L == 40: The input MASS value is interpreted as the value of the
!                load at the pole of the harmonic load, in units of kg/m**2.
!                The inputs H and AMPLITUDE are ignored and set to zero
!                in output. 
! 
!
!
!
!
! Declarations ...
!
! 
!
 use COMMON; use TRIGFCN; implicit NONE
!!
!    
 real(dp)      :: AMPLITUDE    ! Half-amplitude of the load
 real(dp)      :: MASS         ! Maximum load mass
 real(dp)      :: H            ! Maximum load thickness
 real(dp)      :: PRESS        ! Mass/surface at the load pole
 real(dp)      :: aux          ! Auxiliary variable
 real(dp)      :: leg, cheb    ! Legendre & Chebichev polynomials
!
 integer (i4b) :: tipo         ! Type of harmonic load
 integer (i4b) :: l            ! Harmonic degree
 integer (i4b) :: Iw           ! Switch for writing or not on load_coeff.dat
!
!
!
!
! +--------------------------------------+
! | (Some) Warning and Error Conditions  |
! +--------------------------------------+   
!
!
  IF (Tipo<=21) THEN
!
  IF (h < 0.d0 .and. iv==0) THEN
   write(99,*)'ERROR in AXIS_LOAD. The load thickness is negative '
   write(99,*)'**** JOB ABORTED ********************************* '; Stop
  ENDIF
  IF (h < 0.d0 .and. iv==1) THEN
   write(99,*)'ERROR in AXIS_LOAD. The load thickness is negative '
   write(99,*)'**** JOB ABORTED ********************************* '
   write(*,*) 'ERROR in AXIS_LOAD. The load thickness is negative '
   write(*,*) '**** JOB ABORTED ********************************* '; Stop
  ENDIF
!
  IF (h == 0.d0 .and. iv==0) THEN
   write(99,*) 'WARNING in AXIS_LOAD. The load thickness is = 0'
  ENDIF
  IF (h == 0.d0 .and. iv==1) THEN
   write(99,*) 'WARNING in AXIS_LOAD. The load thickness is = 0'
   write(* ,*) 'WARNING in AXIS_LOAD. The load thickness is = 0'
  ENDIF
!
  IF ((amplitude <= 0.d0 .or. amplitude >= 180.d0) .and. iv==0) THEN
   write(99,*)'ERROR in AXIS_LOAD: The load half--amplitude is out of range '
   write(99,*)'********** Job aborted ************************************* '; stop
  endif 
!
  IF ((amplitude <= 0.d0 .or. amplitude >= 180.d0) .and. iv==1) THEN
   write(99,*)'ERROR in AXIS_LOAD: The load half--amplitude is out of range '
   write(99,*)'********** Job aborted ************************************* '
   write(*,*) 'ERROR in AXIS_LOAD: The load half--amplitude is out of range '
   write(*,*) '********** Job aborted ************************************* '; Stop
  endif 
!
  ENDIF   ! on tipo <= 21
!
!
!
 If (Tipo == 30) Then
!   
  if(mass < 0.d0 .and. iv==0) then
   write(99,*) 'ERROR in AXIS_LOAD. The  mass  is negative '
   write(99,*) '**** JOB ABORTED ************************* '; Stop
  endif
  if(mass < 0.d0 .and. iv==1) then
   write(99,*) 'ERROR in AXIS_LOAD. The  mass  is negative '
   write(99,*) '**** JOB ABORTED ************************* '
   write(*,*)  'ERROR in AXIS_LOAD. The  mass  is negative '
   write(*,*)  '**** JOB ABORTED ************************* '; Stop
  endif 
!
  if(mass == 0.d0 .and. iv==0) then
   write(99,*) 'WARNING in AXIS_LOAD. The  mass   is  = 0  '
  endif
  if(mass == 0.d0 .and. iv==1) then
   write(99,*) 'WARNING in AXIS_LOAD. The  mass   is  = 0  '
   write(* ,*) 'WARNING in AXIS_LOAD. The  mass   is  = 0  '
  endif  
!
  Endif    ! on tipo = 30
!		 
!
!
  If (Tipo == 40) Then
!   
  if(mass  < 0.d0 .and. iv==0) then
   write(99,*) 'ERROR in AXIS_LOAD. The <<pressure>> is negative '
   write(99,*) '**** JOB ABORTED ******************************* '; Stop
  endif
  if(mass  < 0.d0 .and. iv==1) then
   write(99,*) 'ERROR in AXIS_LOAD. The <<pressure>> is negative '
   write(99,*) '**** JOB ABORTED ******************************* '
   write (*,*) 'ERROR in AXIS_LOAD. The <<pressure>> is negative '
   write (*,*) '**** JOB ABORTED ******************************* '; Stop
  endif 
!
  if(mass == 0.d0 .and. iv==0) then
   write(99,*) 'WARNING in AXIS_LOAD. The <<pressure>>  is  = 0  '
  endif
  if(mass == 0.d0 .and. iv==1) then
   write(99,*) 'WARNING in AXIS_LOAD. The <<pressure>>  is  = 0  '
   write(* ,*) 'WARNING in AXIS_LOAD. The <<pressure>>  is  = 0  '
  endif  
!
  Endif   ! on tipo = 40
!
!
!
!
!
! # The vector sigma is initialized
!
    sigma(:) = 0.d0 
!    
!
!
!***************
!  "Disc Load" *
!***************
!       
IF(TIPO==10) THEN         
!
mass = h*2.d0*pi*rhoice*raggio*raggio*(1._dp-cosd(amplitude)) 
!
!
!
do l = lmin, lmax                 
!
 SIGMA(l) = (-0.5d0)*rhoice*h*                                &
             ( 		             			      &
 leg(l+1,0,cosd(amplitude)) - leg(l-1,0,cosd(amplitude))      &
             )                                                 
!
enddo 
!
 SIGMA(0) = rhoice*h*(1._dp - cosd(amplitude))/2._dp
!
 SIGMA(1) =  (-0.5d0)*rhoice*h*                                &
             ( 		             			      &
 leg(1+1,0,cosd(amplitude)) - leg(1-1,0,cosd(amplitude))      &
             )                                                 
!
ENDIF  
!
!
!********************************
!  "Disc load + secondary load" *
!********************************
!
IF(TIPO==11) THEN         
!
mass = h*2.d0*pi*rhoice*raggio*raggio*(1._dp-cosd(amplitude))
!
!
!
!
do l = lmin, lmax 
!
 SIGMA(l) = -(rhoice*h)*                                      &
             ( 		             			      &
 leg(l+1,0,cosd(amplitude)) - leg(l-1,0,cosd(amplitude))      &
             )/                                               &
             (1._dp + cosd(amplitude)) 
!
enddo                    
!
 SIGMA(0) = 0._dp  ! mass is conserved
!
 SIGMA(1) = -(rhoice*h)*                                      &
             ( 		             			      &
 leg(1+1,0,cosd(amplitude)) - leg(1-1,0,cosd(amplitude))      &
             )/                                               &
             (1._dp + cosd(amplitude))  
!
ENDIF  
!
!
!********************
! "Parabolic load"  *
!********************
!
IF(TIPO==20) THEN         
!
mass = h*(4.d0/3.d0)*pi*rhoice*raggio*raggio*(1._dp-cosd(amplitude))
!
!
do l = lmin, lmax 
!
aux =    (cheb(l+1,amplitude)-cheb(l+2,amplitude))/(float(l)+1.5d0) -     &
         (cheb(l-1,amplitude)-cheb(l+0,amplitude))/(float(l)-0.5d0) 
!
 SIGMA(l) = - rhoice*h*AUX/4._dp/(1._dp-cosd(amplitude))  
!
enddo                    
!
 SIGMA(0) =   rhoice*h*(1._dp - cosd(amplitude))/3._dp                    
!
aux =    (cheb(1+1,amplitude)-cheb(1+2,amplitude))/(float(1)+1.5d0) -     &
         (cheb(1-1,amplitude)-cheb(1+0,amplitude))/(float(1)-0.5d0)
!
 SIGMA(1) = - rhoice*h*AUX/4._dp/(1._dp-cosd(amplitude))
!
ENDIF  
!
!
!*************************************
! "Parabolic load" + secondary load  *
!*************************************
!
IF(TIPO==21) THEN        
!
mass = h*(4.d0/3.d0)*pi*rhoice*raggio*raggio*(1._dp-cosd(amplitude))
!
!
!
do l = lmin, lmax 
!
aux =    (cheb(l+1,amplitude)-cheb(l+2,amplitude))/(float(l)+1.5d0) -        &
         (cheb(l-1,amplitude)-cheb(l+0,amplitude))/(float(l)-0.5d0) 
!
 SIGMA(l) = aux 
! 
aux =    - (leg(l+1,0,cosd(amplitude)) - leg(l-1,0,cosd(amplitude)))       * & 
         4._dp*(1._dp - cosd(amplitude))**2/3._dp/(1._dp+cosd(amplitude))   ! Changed 'dcosd' to 'cosd'
!
 SIGMA(l) = SIGMA(l) + aux 
!
 SIGMA(l) = - rhoice*h/4._dp/(1._dp - cosd(amplitude)) * SIGMA(l) 
!
enddo                    
!
 SIGMA(0) = 0._dp  ! mass is conserved 
!
aux =    (cheb(1+1,amplitude)-cheb(1+2,amplitude))/(float(1)+1.5d0) -        &
         (cheb(1-1,amplitude)-cheb(1+0,amplitude))/(float(1)-0.5d0)
!
 SIGMA(1) = aux
! 
aux =    - (leg(1+1,0,cosd(amplitude)) - leg(1-1,0,cosd(amplitude)))       * &
         4._dp*(1._dp - cosd(amplitude))**2/3._dp/(1._dp+cosd(amplitude))  ! Changed 'dcosd' to 'cosd'
!
 SIGMA(1) = SIGMA(1) + aux
!
 SIGMA(1) = - rhoice*h/4._dp/(1._dp - cosd(amplitude)) * SIGMA(1)
!
!
ENDIF  
!
!
!*************
! Point Mass *
!*************
!
IF(TIPO==30) THEN         
!
AMPLITUDE=0.d0; H=0.d0  
!
do l = lmin, lmax 
 SIGMA(l) = mass*(2._dp*l + 1._dp)/4._dp/pi/raggio/raggio 
enddo                    
!
 SIGMA(0) = mass/4._dp/pi/raggio/raggio
!
 SIGMA(1) = mass*(2._dp*1 + 1._dp)/4._dp/pi/raggio/raggio
!
ENDIF  
!
!
!
!
!**********************
! Zonal Harmonic Load *
!**********************
!
!
IF(TIPO==40) THEN         
!
PRESS=MASS; AMPLITUDE=0.d0; H=0.d0 
! 
sigma(lmin)=PRESS 
!
    if(lmin==0) sigma(0) = PRESS 
    if(lmin/=0) sigma(0) = 0.d0 
!
ENDIF 
!
!
!
!
!
!
!
!
!
!
 If( Iw == 1) Then    ! # Writes on load_coeff.dat, if requested
!
!
 open(11,file='load_coeff.dat',status='unknown') 
!
 WRITE(11,*) '########################################'
 write(11,*) '# AXIS SYMMETRIC load of type =', tipo
 WRITE(11,*) '# degree & load (kg/m**2)      '
!
 aux=0.d0;        write(11,'(i4,2x,E14.8)') int(aux), sigma(0)
!
 do l=lmin, lmax; write(11,'(i4,2x,E14.8)') l, SIGMA(l); enddo
!
 close(11)
!
              Endif
!
!
!
!
! 
END SUBROUTINE AXIS_LOAD 
!
!
!
!
!
!
!
!
! -------------------------------------------------------------------+
!			     	        		             |
  SUBROUTINE RECT_LOAD0 (iw, FI_C, TE_C, DFI, DTE,   H, MASSA)   !
!						                     |
! -------------------------------------------------------------------+
!
!       This routine computes the harmonic coefficients FF(l,m) and
! GG(l,m) for a load characterized by a rectangular base, i.e, a
! surface mass distribution of uniform thickness H (m) placed on the
! region between two meridian and two parallels. The centroid of the
! load has longitude FI_C (deg) and colatitude TE_C (deg), and DFI, DTE
! (deg) are the widths of the load along the longitude and colatitude,
! respectively. For such a load, the compensation on a complementary
! ocean load is not possible, but that on oceans of realistic shape
! is indeed possible. The mass of the primary load MASSA is computed
! on the basis of the geometrical features of the load provided in
! input. For IW==1, the FF and GG coefficients are written on file
! load_coeff.dat
!
!
!
use COMMON; use TRIGFCN; implicit NONE
!
!!
REAL(dp),parameter ::dfi_min =  .25d0 ! Minimum allowed width in longitude (deg)
REAL(dp),parameter ::dfi_max = 90.d0 ! Maximum allowed width in longitude (deg)
REAL(dp),parameter ::dte_min =  .25d0 ! Minimum allowed width in colatitude (deg)
REAL(dp),parameter ::dte_max = 90.d0 ! Maximum allowed width in colatitude (deg)
!
integer (i4b) :: l, m; COMMON /AREA/ L,M   ! A common for the degree & order
!
real(dp) :: h, massa         ! Thickness (m) & load mass (kg)
!
real(dp) :: fi_1, fi_2, dfi  ! Top & bottom longitude bounds, and width
real(dp) :: te_1, te_2, dte  ! Left & right longitude bounds, and width
!
real(dp) :: fi_c             ! Longitude of the centroid
real(dp) :: te_c             ! Ciolatitude of the centroid
!
real(dp) :: alfa(llmin:llmax), &
            beta(llmin:llmax), &
            dom (llmin:llmax)        ! Auxiliary variables
!
real(dp) :: aux, sum       ! Auxiliary variables
real(dp) :: gamma_lm       ! Integral of P_lm
real(dp) :: common_factor  ! A constant
!
real(dp) :: funcp  ; external FUNCP     ! Function to be integrated
!
real(dp) :: GPROD 	                ! Function Gprod = (l-m)!/(l+m)!
!
!
Real(dp) :: x1, x2         ! Bounds of the integration interval
Real(dp) :: CSI(1:38)      ! Bounds of the colatitude intervals
Real(dp) :: del            ! Lenght of a single interval
Integer (i4b) :: N_sub     ! Number of subdivisions in longitude
Integer (i4b) :: Iw        ! =0/1, to print on 'load_coeff.dat' 
!
Integer (i4b) :: kk        ! A do-loop control variable
!
!
!
!
!
!
! # Defines FI_k & TE_k (k=1,2)
!
  FI_1 = FI_C - DFI/2.D0
  FI_2 = FI_C + DFI/2.D0 
!  
  TE_1 = TE_C - DTE/2.D0
  TE_2 = TE_C + DTE/2.D0 
!
! # Some tests
!
  If ((dfi < dfi_min .or. dfi > dfi_max) .and. iv==0) then
   write(99,*) 'ERROR in RECT_LOAD0: The amplitude in longitude of the '
   write(99,*) 'rectangular load is out of range *** JOB ABORTED ***** '; stop
  endif
!
  If ((dfi < dfi_min .or. dfi > dfi_max) .and. iv==1) then
   write(99,*) 'ERROR in RECT_LOAD0: The amplitude in longitude of the '
   write(99,*) 'rectangular load is out of range *** JOB ABORTED ***** '
   write (*,*) 'ERROR in RECT_LOAD0: The amplitude in longitude of the '
   write (*,*) 'rectangular load is out of range *** JOB ABORTED ***** '; stop
  endif
!
  If ((dte < dte_min .or. dte > dte_max) .and. iv==0) then
   write(99,*) 'ERROR in RECT_LOAD0: The amplitude in colatitude of the '
   write(99,*) 'rectangular load is out of range *** JOB ABORTED  ***** '; stop
  endif
!
  If ((dte < dte_min .or. dte > dte_max) .and. iv==1) then
   write(99,*) 'ERROR in RECT_LOAD0: The amplitude in colatitude of the '
   write(99,*) 'rectangular load is out of range *** JOB ABORTED  ***** '
   write(*,*)  'ERROR in RECT_LOAD0: The amplitude in colatitude of the '
   write(*,*)  'rectangular load is out of range *** JOB ABORTED  ***** '; stop
  endif
!
 IF( (FI_C < DFI/2. .OR. FI_C > 360. - DFI/2.) .and. IV ==0 ) THEN
  WRITE(99,*) 'ERROR in RECT_LOAD0: The longitude of the centroid must  '
  WRITE(99,*) 'vary in the range [dfi/2:360-dfi/2]. *** JOB ABORTED *** '; STOP
 ENDIF
!
 IF( (FI_C < DFI/2. .OR. FI_C > 360. - DFI/2.) .and. IV ==1 ) THEN
  WRITE(99,*) 'ERROR in RECT_LOAD0: The longitude of the centroid must  '
  WRITE(99,*) 'vary in the range [dfi/2:360-dfi/2]. *** JOB ABORTED *** '
  WRITE(*,*)  'ERROR in RECT_LOAD0: The longitude of the centroid must  '
  WRITE(*,*)  'vary in the range [dfi/2:360-dfi/2]. *** JOB ABORTED *** '; STOP
 ENDIF
!
 IF( (TE_C < DTE/2. .OR. TE_C > 360. - DTE/2.) .and. IV ==0 ) THEN
  WRITE(99,*) 'ERROR in RECT_LOAD0: The colatitude of the centroid must  '
  WRITE(99,*) 'vary in the range [dfi/2:360-dfi/2].  *** JOB ABORTED *** '; STOP
 ENDIF
!
 IF( (TE_C < DTE/2. .OR. TE_C > 360. - DTE/2.) .and. IV ==1 ) THEN
  WRITE(99,*) 'ERROR in RECT_LOAD0: The colatitude of the centroid must  '
  WRITE(99,*) 'vary in the range [dfi/2:360-dfi/2].  *** JOB ABORTED *** '
  WRITE(*,*)  'ERROR in RECT_LOAD0: The colatitude of the centroid must  '
  WRITE(*,*)  'vary in the range [dfi/2:360-dfi/2].  *** JOB ABORTED *** '; STOP
 ENDIF
!
!
  If (h==0.d0 .and. iv==0) then
   write(99,*) 'WARNING in RECT_LOAD0: The load Height is = 0 '
  endif
  If (h==0.d0 .and. iv==1) then
   write(99,*) 'WARNING in RECT_LOAD0: The load Height is = 0 '
   write(* ,*) 'WARNING in RECT_LOAD0: The load Height is = 0 '
  endif
!
  If (h<0.d0 .and. iv==0) Then
   Write(99,*) 'ERROR in RECT_LOAD0: The Height is < 0 '
   Write(99,*) '**** JOB ABORTED ******************** '; Stop
  Endif
  If (h<0.d0 .and. iv==1) Then
   Write(99,*) 'ERROR in RECT_LOAD0: The Height is < 0 '
   Write(99,*) '**** JOB ABORTED ******************** '
   Write(*,*)  'ERROR in RECT_LOAD0: The Height is < 0 '
   Write(*,*)  '**** JOB ABORTED ******************** '; Stop
  Endif     
!  
!
! # Coefficients depending only upon 'm'
!

        Dom (0) = 1.d0 
	Alfa(0) = (FI_2 - FI_1)*PI/180._DP 
	Beta(0) = 0.D0 

        Do m= 1, LMAX 
          DOM (m)    =     0._DP 
          ALFA(m)    =     (SIND(M*FI_2)-SIND(M*FI_1))/FLOAT(M)    !  
          BETA(m)    =     (COSD(M*FI_2)-COSD(M*FI_1))/FLOAT(M)    !   
        Enddo     
!
!
! +----------------------------------------------------------------+
!  If Dte > 5 degrees, the integration in Co-latitude is broken in 
!  several integrations on intervals with amplitude <= 5 degrees.
! +----------------------------------------------------------------+
!
        If(Dte <= 5.) Then 
!	
	            N_sub = 1
                    x1 =   Cosd(Te_1) 
		    x2 =   Cosd(Te_2)
		    del = Abs(x2 - x1) 
		    Csi(1)= x1
		    Csi(2)= x2   		    
!
		       Else 
!		  
                    x1    =    Cosd(Te_1)                 ! Left  Bound 
		    x2    =    Cosd(Te_2)                 ! Right Bound
		    N_sub =  Int( (Te_2-Te_1)/5. + 1)     ! Number of intervals
		    del   =  Abs(x2-x1)/Float(N_sub)      ! Lenght of Each interval
!		    	       
                    Csi(1) = x1 
		    Do kk=2, N_sub 
		    Csi(kk) = Csi(1) - (kk-1)*del
		    Enddo  
		    Csi(N_sub+1) = x2 
!		    
!
                       Endif 

!
!
      DO  L = Lmin, Lmax 
!
         DO  M = 0, L
!
         sum=0d0 
!
         Do kk=1, n_sub 
!
         CALL QGAUS (FUNCP, CSI(kk), CSI(kk+1), GAMMA_LM)
!
         COMMON_FACTOR = (2.D0 - DOM(m)) * ((2.D0*L+1.D0)/4.D0/PI)*    & 
    	      	          GPROD(L,M) * GAMMA_LM * RHOICE * H  
!
         Sum = Sum + Common_Factor
!
         Enddo
!
         FF(L,M) =     ALFA(m) * Sum  
         GG(L,M) =    -BETA(m) * Sum  
!
!
         ENDDO
!	  
      ENDDO 
!
!
! # Defines SIGMA(0)
!
      SIGMA(0) = RHOICE*H*((FI_2 - FI_1)*PI/180._DP)* & 
                 (COSD(TE_1)-COSD(TE_2))/4._DP/PI 
!
      FF(0,0) = SIGMA(0) 		 
      GG(0,0) = 0.D0  
!
!
! # Computes the MASS of the load
!    
      MASSA = SIGMA(0)*4.d0*RAGGIO*RAGGIO*PI  
!     
!
! # Defines the degree 1 coefficients

      DO  L = 1, 1
         DO  M = 0, L
!
         sum=0d0 
!
         Do kk=1, n_sub
           CALL QGAUS (FUNCP, CSI(kk), CSI(kk+1), GAMMA_LM)
           COMMON_FACTOR = (2.D0 - DOM(m)) * ((2.D0*L+1.D0)/4.D0/PI)*    &
    	         	    GPROD(L,M) * GAMMA_LM * RHOICE * H
!
         Sum = Sum + Common_Factor
!
         Enddo
!
         FF(L,M) =     ALFA(m) * Sum  
         GG(L,M) =    -BETA(m) * Sum  
!
         ENDDO
      ENDDO 
!
!
!
  If(iw == 1) then ! # Writes on load_coeff.dat, if requested
!
  open(11,file='load_coeff.dat',status='unknown') 
!

 WRITE(11,*) '########################################'
 write(11,*) '#     RECTANGULAR load of type 50       '
 WRITE(11,*) '#   degree,  order, & load (kg/m**2)    '
!
!
!
  aux=0.d0;         write(11,44)  int(aux), int(aux), FF(0,0), GG(0,0)
!
  do l=lmin, lmax; do m=0,l
!
                    write(11,44)  l,        m,        FF(l,m), GG(l,m)
!
  enddo; enddo
!
  close(11) 
!
  44 format(2(i4,1x),2(E14.8,1x))
!
              Endif
!
!
!
END SUBROUTINE RECT_LOAD0
!
!
!







!
!
        FUNCTION FUNCP(X)
!
! # FUNCP is defined as FUNCP = -P_lm(x)
!
        USE STRATA; IMPLICIT NONE
        REAL(DP)      :: FUNCP, LEG, X
	INTEGER (I4B) :: L,M
	COMMON /AREA/ L,M
!
        FUNCP  = - LEG (L,M,X)
!
        END FUNCTION FUNCP    
!
!
!
      SUBROUTINE QGAUS (FUNC,A,B,SS)
!
! # Computes the integral of FUNC in the
!   interval (a,b) using the GAUSS method
!
      USE STRATA; IMPLICIT NONE
      REAL(DP) :: a,b,ss,func
      EXTERNAL func
      INTEGER (I4B) j
      REAL(DP) ::  dx,xm,xr,w(5),x(5)
      SAVE w,x
      DATA w/.2955242247,.2692667193,.2190863625,.1494513491, &
      .0666713443/
      DATA x/.1488743389,.4333953941,.6794095682,.8650633666, &
      .9739065285/
      xm=0.5*(b+a)
      xr=0.5*(b-a)
      ss=0
      do 11 j=1,5
        dx=xr*x(j)
        ss=ss+w(j)*(func(xm+dx)+func(xm-dx))
11    continue
      ss=xr*ss
      return
      END SUBROUTINE QGAUS 
!  (C) Copr. 1986-92 Numerical Recipes Software W^3.
!                                                             
!
!
!
!  -------------------o
 function cheb (n,x)  !
!  -------------------o
! chebichev polynomials of the 1st kind T_n(x)
!
use COMMON; use TRIGFCN
implicit none 
real(dp) :: cheb, x
integer(i4b) :: n 
!
if (n>=0) then 
  cheb=cosd(n*x) ! x must be given in degrees
else 
write(99,*) 'ERROR in CHEB: n must be >= 0' 
write(99,*) 'The current value is     ', n
write(99,*) '*** JOB ABORTED *************'; Stop  
endif 
!
! 
end function cheb 
!
!
!
!+------------------------------+
	function leg (l,m,x)    !
!+------------------------------+ 	
        use COMMON 
        implicit none 
        real(dp) :: leg, x, pmm, somx2, fact1, pmmp1, pll
        integer (i4b) :: l,  im, i, ll, m
!
! # Computes the associated Legendre polynomial P_lm (x) in
!   double precision. The factor (-1)**m is included in P_lm.
!   l and m must be >=0, with m<=l, and |x| must be <=1.
!   (Adapted from Numerical Recipes).
!
!
        if(abs(x)>1.or.abs(m)>l) then
        write(99,*) 'ERROR in LEG: X out of range and/or |m| > l'
	write(99,*) 'x    = ', x 
	write(99,*) 'l, m = ', l, m 
        write(99,*) ' ****** JOB ABORTED ***********************'; STOP 
        endif
        im=abs(m)
        pmm=1.d0
        if(im>0) then
            somx2=sqrt((1.d0-x)*(1.d0+x))
            fact1=1.d0
            do 11 i=1,im
                pmm=-pmm*fact1*somx2
                fact1=fact1+2.d0
  11        continue
        endif
        if(l==im) then
            leg=pmm
        else
            pmmp1=x*(2*im+1)*pmm
            if(l==im+1) then
                leg=pmmp1
            else
                do 12 ll=im+2,l
                   pll=(x*(2*ll-1)*pmmp1-(ll+im-1)*pmm)/(ll-im)
                   pmm=pmmp1
                   pmmp1=pll
  12            continue
             leg=pll
            endif
         endif
         return
         end function leg    
!
!
!

!
!+-----------------------------------+
         function legdev1(l,m,x)     !
!+-----------------------------------+	 
         use COMMON
         implicit none 
         real(dp) :: leg, legdev1, x 
         integer (i4b) :: l,  m  
!
! # Computes the derivative of P_lm(cos(teta)) wrt teta.
!   See Eq. D.149 (pag. 981) of Ben-Menahem & Singh
!
        if(abs(x)>1.or.abs(m)>l) then
        write(99,*) 'ERROR in LEGDEV1: X out of range and/or |m| > l'
        write(99,*) ' ****** JOB ABORTED ***************************'; STOP 
        endif
!
        if(x==1._dp.or.x==-1._dp) then
        legdev1=0._dp 
        return
                                      endif
!
        legdev1=((l+1)*x*leg(l,m,x)-(l-m+1)*leg(l+1,m,x))/(1.d0-x)/(1.d0+x)
        legdev1 = - legdev1*sqrt(1.d0-x*x)
!
        end function legdev1 
!
!
!
!
!
!
       FUNCTION d_load_history (tim, nt)
!
! +----------------------------------------------------------+
! | For a given time TIM, this function provides the value   |
! | of the derivative of the time-history 'nt' at that time  |
! | [01.02.02] Revised on 29 August 2002                     |
! +----------------------------------------------------------+
!
       use COMMON
       implicit NONE
!
!
       Real(sp) :: d_load_history  ! Derivative of the time- history
       real(sp) :: tim             ! local time (kyr)
       REAL(sp) :: d_f0, d_f1, d_f2, &
                   d_f3, d_f4, d_f5, &
                   d_f6, d_f7, d_f8
!
!
       Integer(i4b) :: Nt          ! Label of the Th
!
!
       if (nt==0) d_load_history = d_f0 (tim)
       if (nt==1) d_load_history = d_f1 (tim)
       if (nt==2) d_load_history = d_f2 (tim)
       if (nt==3) d_load_history = d_f3 (tim)
       if (nt==4) d_load_history = d_f4 (tim)
       if (nt==5) d_load_history = d_f5 (tim)
       if (nt==6) d_load_history = d_f6 (tim)
       if (nt==7) d_load_history = d_f7 (tim)
       if (nt==8) d_load_history = d_f8 (tim)
!
       END FUNCTION d_load_history
!
!
!
!

!
!
!
       FUNCTION load_history (nt, tim)
!
! +----------------------------------------------------------+
! | For a given time TIM, this function provides the value   |
! | of the time-history 'nt' at that time [01.02.02]         |
! | Revised on 29 August 2002                                |
! +----------------------------------------------------------+
!
       use COMMON
       implicit NONE
!
       Real(sp) :: load_history  ! Value of the time- history
!
       real(sp) :: tim           ! local time variable (Kyr)
       REAL(sp) :: f0, f1, f2, &
                   f3, f4, f5, &
                   f6, f7, f8
!
!
       Integer(i4b) :: Nt        ! Label of the Th
!
!
       if (nt==0) load_history = f0 (tim)
       if (nt==1) load_history = f1 (tim)
       if (nt==2) load_history = f2 (tim)
       if (nt==3) load_history = f3 (tim)
       if (nt==4) load_history = f4 (tim)
       if (nt==5) load_history = f5 (tim)
       if (nt==6) load_history = f6 (tim)
       if (nt==7) load_history = f7 (tim)
       if (nt==8) load_history = f8 (tim)
!
       END FUNCTION load_history
!
!
!
!
!
function th_0 (t)
!
! Instantaneous loading (Heaviside)
!
 use COMMON; implicit NONE
 REAL(sp) :: th_0, f0, d_f0
 REAL(sp) :: t       ! local time
!
!
!+--------------------------------------------------------+
! Definition of the Th#0: Instantaneous loading (Heaviside)
!+--------------------------------------------------------+
!
ENTRY f0 (t)
               f0 = 0.0
IF(t >= 0.0)   f0 = 1.0
RETURN
!
!+--------------------------------------------------------+
! Derivative of the Th#0: Instantaneous loading (Heaviside)
!+--------------------------------------------------------+
!
ENTRY d_f0 (t)
               d_f0 = 0.0
RETURN
!
end function th_0
!
!
!
!
!
subroutine Convol_0_new  (t)
!
! Instantaneous loading (Heaviside)
!
 use COMMON
 implicit NONE
!
 REAL(sp) :: t          ! << local time >> variable
 REAL(sp) :: AUX        ! Auxilium
 REAL(sp) :: f0, d_f0   ! Time-history and its derivative
!
 INTEGER (i4b) :: m  ! Modes index
 INTEGER (i4b) :: L  ! Harmonic Degree
!
!
!+------------+
! Th#0 (*) ldcs
!+------------+
!
ENTRY cnv_0 (t)
!
  conv_h(:) = 0.0
  conv_l(:) = 0.0
  conv_k(:) = 0.0
!
            IF(IR==1) conv_h(:) =        h_e(:)  * f0(t)
            IF(it==1) conv_l(:) =        l_e(:)  * f0(t)
            IF(IG==1) conv_k(:) =   (1.0+k_e(:)) * f0(t)
!
!
!
  IF(Only_Elastic ==0) THEN
!
   IF(t>=0.0) THEN
!
  do l=lmin, lmax
     do m=1,nroots
         aux = exp(-t/tek(l,m))
         IF(IR==1) conv_h(l) = conv_h(l) -  r_h (l,m) * (1.0-AUX)
         IF(it==1) conv_l(l) = conv_l(l) -  r_l (l,m) * (1.0-AUX)
         IF(IG==1) conv_k(l) = conv_k(l) -  r_k (l,m) * (1.0-AUX)
     enddo
  enddo
!
  ENDIF
!
  ENDIF
!
!
RETURN
!
!
!+--------------------+
! d/dt (Th#0 (*) ldcs)
!+--------------------+
!
ENTRY dcnv_0 (t)
!
  dconv_h(:) = 0.0
  dconv_l(:) = 0.0
  dconv_k(:) = 0.0
!
!
            IF(dIR==1) dconv_h(:) =        h_e(:)  * d_f0(t)
            IF(dit==1) dconv_l(:) =        l_e(:)  * d_f0(t)
            IF(dIG==1) dconv_k(:) =   (1.0+k_e(:)) * d_f0(t)
!
!
!
  IF(Only_Elastic ==0 .and. t>= 0.0) THEN
!
  do l=lmin, lmax
     do m=1,nroots
         aux = exp(-t/tek(l,m))
         IF(dIR==1) dconv_h(l) = dconv_h(l) - r_h(l,m) * AUX /tek(l,m)
         IF(dit==1) dconv_l(l) = dconv_l(l) - r_l(l,m) * AUX /tek(l,m)
         IF(dIG==1) dconv_k(l) = dconv_k(l) - r_k(l,m) * AUX /tek(l,m)
     enddo
  enddo
!
  ENDIF
!
!
!
RETURN
!
!
RETURN
!
end subroutine  Convol_0_new
!
!
!
!
!
!
!
!
!
function th_1 (t)
!
! (Instantaneous un- loading)
!
 use COMMON; implicit NONE
 REAL(sp) :: th_1, f1, d_f1
!
 REAL(sp) :: t       ! local time
!
!
!+--------------------------------------------------+
! Definition of the Th#1: (Instantaneous un- loading)
!+--------------------------------------------------+
!
ENTRY f1 (t)
              f1 = 1.0
IF(t >= 0.0)  f1 = 0.0
RETURN
!
!+--------------------------------------------------+
! Derivative of the Th#1: (Instantaneous un- loading)
!+--------------------------------------------------+
!
ENTRY d_f1 (t)
               d_f1 = 0.0
RETURN
!
end function th_1
!
!
!
!
!
!
!
subroutine Convol_1_new (t)
!
! (Instantaneous Un- loading)
!
 use COMMON; implicit NONE
!
 REAL(sp) :: t          ! << local time >> variable
 REAL(sp) :: AUX        ! Auxilium
 REAL(sp) :: f1, d_f1   ! Time-history and its derivative
!
 INTEGER (i4b) :: m  ! Modes index
 INTEGER (i4b) :: L  ! Harmonic Degree
!
!
!+------------+
! Th#1 (*) ldcs
!+------------+
!
ENTRY cnv_1 (t)
!
  conv_h(:) = 0.0
  conv_l(:) = 0.0
  conv_k(:) = 0.0
!
!
            IF(IR==1) conv_h(:) =        h_e(:)  *f1(t)
            IF(it==1) conv_l(:) =        l_e(:)  *f1(t)
            IF(IG==1) conv_k(:) =   (1.0+k_e(:)) *f1(t)
!
!
!
  IF(only_elastic==0) THEN
!
    IF(IR==1) conv_h(:) =  conv_h(:) +   ( h_f(:) - h_e(:) )
    IF(it==1) conv_l(:) =  conv_l(:) +   ( l_f(:) - l_e(:) )
    IF(IG==1) conv_k(:) =  conv_k(:) +   ( k_f(:) - k_e(:) )

      IF( t>=0.0 ) THEN
!
    do l=lmin, lmax
     do m=1,nroots
         aux = exp(-t/tek(l,m))
         IF(IR==1) conv_h(l) = conv_h(l) + r_h(l,m) * (1.0-AUX)
         IF(it==1) conv_l(l) = conv_l(l) + r_l(l,m) * (1.0-AUX)
         IF(IG==1) conv_k(l) = conv_k(l) + r_k(l,m) * (1.0-AUX)
     enddo
    enddo
!
  ENDIF
!
  ENDIF
!
  RETURN
!
!
!
!+---------------------+
! (d/dt) Th#1 (*) ldcs
!+---------------------+
!
ENTRY dcnv_1 (t)
!
  dconv_h(:) = 0.0
  dconv_l(:) = 0.0
  dconv_k(:) = 0.0
!
!
            IF(dIR==1) dconv_h(:) =        h_e(:)  *d_f1(t)
            IF(dit==1) dconv_l(:) =        l_e(:)  *d_f1(t)
            IF(dIG==1) dconv_k(:) =   (1.0+k_e(:)) *d_f1(t)
!
!
!
  IF(only_elastic==0 .and. t>=0.0) THEN
!
    do l=lmin, lmax
     do m=1,nroots
         aux = exp(-t/tek(l,m))/tek(l,m) 
         IF(dIR==1) dconv_h(l) = dconv_h(l) + r_h (l,m) * AUX
         IF(dit==1) dconv_l(l) = dconv_l(l) + r_l (l,m) * AUX
         IF(dIG==1) dconv_k(l) = dconv_k(l) + r_k (l,m) * AUX
     enddo
    enddo
!
  ENDIF
!
!
!
  RETURN
!
!
end subroutine  Convol_1_new
!
!
!
!
!
function th_2 (t)
!
!  (loading & unloading)
!
 use COMMON; implicit NONE
 REAL(sp) :: th_2, f2, d_f2
!
 REAL(sp) :: t       ! local time
!
!
!+---------------------------------------------+
! Definition of the Th#2:  (loading & unloading)
!+---------------------------------------------+
!
!
ENTRY f2 (t)
IF(t <  -tau_lu)                 f2 = 0.0
IF(t >= -tau_lu .and. t < 0.0)   f2 = 1.0
IF(t >=  0.0 )                   f2 = 0.0
RETURN
!
!
!+---------------------------------------------+
! Derivative of the Th#2:  (loading & unloading)
!+---------------------------------------------+
!
ENTRY d_f2 (t)
               d_f2 = 0.0
RETURN
!
end function th_2
!
!
!
!
!
!
!
subroutine Convol_2_new (t)
!
! (loading & un- loading)
!
 use COMMON; implicit NONE
!
 REAL(sp) :: t          ! << local time >> variable
 REAL(sp) :: AUX        ! Auxilium
 REAL(sp) :: f2, d_f2   ! Time-history and its derivative
!
 INTEGER (i4b) :: m  ! Modes index
 INTEGER (i4b) :: L  ! Harmonic Degree
!
!
!+------------+
! Th#2 (*) ldcs
!+------------+
!
ENTRY cnv_2 (t)
!
  conv_h(:) = 0.0
  conv_l(:) = 0.0
  conv_k(:) = 0.0
!
            IF(IR==1) conv_h(:) =        h_e(:)  * f2(t)
            IF(it==1) conv_l(:) =        l_e(:)  * f2(t)
            IF(IG==1) conv_k(:) =   (1q0+k_e(:)) * f2(t)
!
!
  IF(only_elastic==0 .and. t >= -tau_lu) THEN
!
    do l=lmin, lmax
     do m=1,nroots
         aux = exp(-(t+tau_lu)/tek(l,m))
         IF(IR==1) conv_h(l) = conv_h(l) - r_h(l,m) * (1.0-AUX)
         IF(it==1) conv_l(l) = conv_l(l) - r_l(l,m) * (1.0-AUX)
         IF(IG==1) conv_k(l) = conv_k(l) - r_k(l,m) * (1.0-AUX)
     enddo
    enddo
!
  ENDIF
!
  IF(only_elastic==0 .and.  t>= 0.0 ) THEN

    do l=lmin, lmax
     do m=1,nroots
         aux = exp(-t/tek(l,m))
         IF(IR==1) conv_h(l) = conv_h(l) + r_h(l,m) *  (1.0-AUX)
         IF(it==1) conv_l(l) = conv_l(l) + r_l(l,m) *  (1.0-AUX)
         IF(IG==1) conv_k(l) = conv_k(l) + r_k(l,m) *  (1.0-AUX)
     enddo
    enddo
!
  ENDIF
!
  RETURN
!
!
!+-------------------+
! (d/dt) Th#1 (*) ldcs
!+-------------------+
!
ENTRY dcnv_2 (t)
!
  dconv_h(:) = 0.0
  dconv_l(:) = 0.0
  dconv_k(:) = 0.0
!
            IF(dIR==1) dconv_h(:) =        h_e(:)  *d_f2(t)
            IF(dit==1) dconv_l(:) =        l_e(:)  *d_f2(t)
            IF(dIG==1) dconv_k(:) =   (1q0+k_e(:)) *d_f2(t)
!
  IF(only_elastic==0 .and. t>= -tau_lu ) THEN
!
    do l=lmin, lmax
     do m=1,nroots
         aux = exp(-(t+tau_lu)/tek(l,m)) / tek(l,m)
         IF(dIR==1) dconv_h(l) = dconv_h(l) - r_h(l,m) * AUX
         IF(dit==1) dconv_l(l) = dconv_l(l) - r_l(l,m) * AUX
         IF(dIG==1) dconv_k(l) = dconv_k(l) - r_k(l,m) * AUX
     enddo
    enddo
!
  ENDIF
!
  IF(only_elastic==0 .and. t>= 0.0) THEN
!
    do l=lmin, lmax
     do m=1,nroots
         aux = exp(-t/tek(l,m)) / tek(l,m)
         IF(dIR==1) dconv_h(l) = dconv_h(l) +  r_h(l,m) * AUX
         IF(dit==1) dconv_l(l) = dconv_l(l) +  r_l(l,m) * AUX
         IF(dIG==1) dconv_k(l) = dconv_k(l) +  r_k(l,m) * AUX
     enddo
    enddo
!
  ENDIF
!
  RETURN
!
!
end subroutine  Convol_2_new
!
!
!
!
!
!
!
function th_3 (t)
!
!  (simple deglaciation)
!
 use COMMON; implicit NONE
 REAL(sp) :: th_3, f3, d_f3
!
 REAL(sp) :: t       ! local time
!
!
!
!+--------------------------------------------+
! Definition of the Th#3: (simple deglaciation)
!+--------------------------------------------+
!
!
ENTRY f3 (t)
IF(t <   0.0)                   f3 = 1.0
IF(t >=  0.0 .and. t < tau_sd)  f3 = 1.0 - t/tau_sd
IF(t >=  tau_sd )               f3 = 0.0
RETURN
!
!+--------------------------------------------+
! Derivative of the Th#3: (simple deglaciation)
!+--------------------------------------------+
!
ENTRY d_f3 (t)
IF(t <   0.0)                   d_f3 =  0.0
IF(t >=  0.0 .and. t < tau_sd)  d_f3 = -1.0/tau_sd
IF(t >=  tau_sd )               d_f3 =  0.0
RETURN
!
end function th_3
!
!
!
!
!
!
!
!
subroutine Convol_3_new (t)
!
! (simple deglaciation)
!
 use COMMON; implicit NONE
!
 REAL(sp) :: t          ! << local time >> variable
 REAL(sp) :: AUX        ! Auxilium
 REAL(sP) :: TAU        ! Auxilium
 REAL(sp) :: f3, d_f3   ! time-history and its derivative
!
 INTEGER (i4b) :: m  ! Modes index
 INTEGER (i4b) :: L  ! Harmonic Degree
!
!
!+------------+
! Th#3 (*) ldcs
!+------------+
!
ENTRY cnv_3 (t)
!
  TAU=TAU_SD
!
  conv_h(:) = 0.0
  conv_l(:) = 0.0
  conv_k(:) = 0.0
!
    IF(IR==1) conv_h(:) =        h_e(:)  * f3(t)
    IF(it==1) conv_l(:) =        l_e(:)  * f3(t)
    IF(IG==1) conv_k(:) =   (1.0+k_e(:)) * f3(t)
!
  IF(only_elastic==0) THEN
!
    IF(IR==1) conv_h(:) =  conv_h(:) +   ( h_f(:) - h_e(:) )
    IF(it==1) conv_l(:) =  conv_l(:) +   ( l_f(:) - l_e(:) )
    IF(IG==1) conv_k(:) =  conv_k(:) +   ( k_f(:) - k_e(:) )
!
  IF(t >= TAU) THEN
!
    IF(IR==1) conv_h(:) =  conv_h(:) -   ( h_f(:) - h_e(:) ) * (1.-t/tau)
    IF(it==1) conv_l(:) =  conv_l(:) -   ( l_f(:) - l_e(:) ) * (1.-t/tau)
    IF(IG==1) conv_k(:) =  conv_k(:) -   ( k_f(:) - k_e(:) ) * (1.-t/tau)
!
    do l=lmin, lmax
     do m=1,nroots
         aux =  - (1. - EXP(-(t-tau)/tek(l,m)) ) / (tau/tek(l,m))
         IF(IR==1) conv_h(l) = conv_h(l) - r_h(l,m)* aux
         IF(it==1) conv_l(l) = conv_l(l) - r_l(l,m)* aux
         IF(IG==1) conv_k(l) = conv_k(l) - r_k(l,m)* aux
     enddo
    enddo
!
  ENDIF
!
  IF(t>= 0.0 ) THEN

    do l=lmin, lmax
     do m=1,nroots
         aux =  ( t/TAU - (1.- exp(-t/tek(l,m))) / (TAU/TEK(l,m)) )
         IF(IR==1) conv_h(l) = conv_h(l) + r_h(l,m)*AUX
         IF(it==1) conv_l(l) = conv_l(l) + r_l(l,m)*AUX
         IF(IG==1) conv_k(l) = conv_k(l) + r_k(l,m)*AUX
     enddo
    enddo
!
  ENDIF
!
  ENDIF
!
  RETURN
!
!
!+-------------------+
! (d/dt) Th#3 (*) ldcs
!+-------------------+
!
ENTRY dcnv_3 (t)
!
  dconv_h(:) = 0.0
  dconv_l(:) = 0.0
  dconv_k(:) = 0.0
!
  IF(dIR==1) dconv_h(:) =       h_e(:)  *d_f3(t)
  IF(dit==1) dconv_l(:) =       l_e(:)  *d_f3(t)
  IF(dIG==1) dconv_k(:) =  (1.0+k_e(:)) *d_f3(t)
!
!
  IF(only_elastic==0) THEN
!
  IF(t>= tau) THEN
!
    IF(dIR==1) dconv_h(:) =  dconv_h(:) +   ( h_f(:) - h_e(:) )/TAU
    IF(dit==1) dconv_l(:) =  dconv_l(:) +   ( l_f(:) - l_e(:) )/TAU
    IF(dIG==1) dconv_k(:) =  dconv_k(:) +   ( k_f(:) - k_e(:) )/TAU
!
    do l=lmin, lmax
     do m=1,nroots
         AUX = EXP(-(t-tau)/tek(l,m))/tau
         IF(dIR==1) dconv_h(l) = dconv_h(l) + r_h(l,m)* AUX
         IF(dit==1) dconv_l(l) = dconv_l(l) + r_l(l,m)* AUX
         IF(dIG==1) dconv_k(l) = dconv_k(l) + r_k(l,m)* AUX
     enddo
    enddo
!
  ENDIF
!
  IF(t>= 0.0 ) THEN

    do l=lmin, lmax
     do m=1,nroots
         aux = (1. - EXP(-(t)/tek(l,m)) )/tau
         IF(dIR==1) dconv_h(l) = dconv_h(l) + r_h(l,m)*AUX
         IF(dit==1) dconv_l(l) = dconv_l(l) + r_l(l,m)*AUX
         IF(dIG==1) dconv_k(l) = dconv_k(l) + r_k(l,m)*AUX
!
     enddo
    enddo
!
  ENDIF
!
  ENDIF
!
  RETURN
!
!
  RETURN
!
!
end subroutine  Convol_3_new
!
!
!
!
!
!
!
function th_4 (t)
!
!  (saw- tooth)
!
 use COMMON; implicit NONE
 REAL(sp) :: th_4, f4, d_f4
!
 REAL(sp) :: t       ! local time
 REAL(SP) :: teta    ! tauc + dinc
 REAL(SP) ::  fn,  fn_up,  fn_do  ! Auxiliary functions
 REAL(SP) :: dfn, dfn_up, dfn_do  ! Their derivatives
!
 INTEGER(i4b) :: n
!
!
!+------------------------------------+
! Definition of the Th#4:  (saw- tooth)
!+------------------------------------+
!
ENTRY f4 (t)
!
teta=tauc+dinc
!
F4=0.0
!
DO 1 n=0,NR
!
fn   =  0.0
!
fn_up=  ( t + n*teta + tauc)/tauc
fn_do=  (-t - n*teta + dinc)/dinc
!
IF ( t >= -n*teta - tauc) fn = fn + fn_up
IF ( t >= -n*teta       ) fn = fn - fn_up
IF ( t >= -n*teta       ) fn = fn + fn_do
IF ( t >= -n*teta + dinc) fn = fn - fn_do
!
F4 = f4 + fn
!
1 CONTINUE
!
RETURN
!
!+------------------------------------+
! Derivative of the Th#4:  (saw- tooth)
!+------------------------------------+
!
ENTRY d_f4 (t)
!
teta=tauc+dinc
!
d_F4=0.0
!
DO 2 n=0,NR
!
dfn   =  0.0
!
dfn_up=  1./tauc
dfn_do= -1./dinc
!
IF ( t >= -n*teta - tauc) dfn = dfn + dfn_up
IF ( t >= -n*teta       ) dfn = dfn - dfn_up
IF ( t >= -n*teta       ) dfn = dfn + dfn_do
IF ( t >= -n*teta + dinc) dfn = dfn - dfn_do
!
d_F4 = d_f4 + dfn
!
2 CONTINUE
!
RETURN
!
end function th_4
!
!
!
!
!
!
!
subroutine Convol_4_new (t)
!
! (saw- tooth)
!
 use COMMON; implicit NONE
!
 REAL(sp) :: t          ! << local time >> variable
 REAL(sp) :: AUX        ! Auxilium
 REAL(sP) :: Teta       ! Auxilium
 REAL(sp) :: f4, d_f4   ! time-history and its derivative
!
 INTEGER (i4b) :: m  ! Modes index
 INTEGER (i4b) :: L  ! Harmonic Degree
 INTEGER (i4b) :: n  ! Index
!
!
!+------------+
! Th#4 (*) ldcs
!+------------+
!
ENTRY cnv_4 (t)
!
teta=tauc+dinc
!
  conv_h(:) = 0.0
  conv_l(:) = 0.0
  conv_k(:) = 0.0
!
    IF(IR==1) conv_h(:) =        h_e(:)  * f4(t)
    IF(it==1) conv_l(:) =        l_e(:)  * f4(t)
    IF(IG==1) conv_k(:) =   (1.0+k_e(:)) * f4(t)
!
!
!
    IF (only_elastic==0) THEN
!
    do 1 n=0, NR
!
!
       IF(t>= -n*teta) then
!
         do l=lmin, lmax
         do m=1,nroots
            aux= (1./tauc+1./dinc)* &
                   ( (t+n*teta) - (1.-EXP(-(t+n*teta)/tek(l,m)) )*tek(l,m))
            IF(IR==1) conv_h(l) = conv_h(l) + r_h(l,m)* AUX
            IF(it==1) conv_l(l) = conv_l(l) + r_l(l,m)* AUX
            IF(IG==1) conv_k(l) = conv_k(l) + r_k(l,m)* AUX
         enddo
         enddo
!
       ENDIF
!
       IF(t>= -n*teta-tauc) then
!
         do l=lmin, lmax
         do m=1,nroots
            aux= (t+n*teta+tauc)/tauc - &
                        (1.-EXP( -(t+n*teta+tauc)/tek(l,m)))/(tauc/tek(l,m))
            IF(IR==1) conv_h(l) = conv_h(l) - r_h(l,m)* AUX
            IF(it==1) conv_l(l) = conv_l(l) - r_l(l,m)* AUX
            IF(IG==1) conv_k(l) = conv_k(l) - r_k(l,m)* AUX
         enddo
         enddo
!
       ENDIF
!
       IF(t>= -n*teta+dinc) then
!
         do l=lmin, lmax
         do m=1,nroots
            aux= (t+n*teta-dinc)/dinc - &
                        (1.-EXP( -(t+n*teta-dinc)/tek(l,m)))/(dinc/tek(l,m))
            IF(IR==1) conv_h(l) = conv_h(l) - r_h(l,m)* AUX
            IF(it==1) conv_l(l) = conv_l(l) - r_l(l,m)* AUX
            IF(IG==1) conv_k(l) = conv_k(l) - r_k(l,m)* AUX
         enddo
         enddo
!
       ENDIF
!
!
 1  CONTINUE
!
    ENDIF
!
!
  RETURN
!
!
!
!
!
!+-------------------+
! (d/dt) Th#4 (*) ldcs
!+-------------------+
!
ENTRY dcnv_4 (t)
!
teta=tauc+dinc
!
  dconv_h(:) = 0.0
  dconv_l(:) = 0.0
  dconv_k(:) = 0.0
!
    IF(dIR==1) dconv_h(:) =        h_e(:)  * d_f4(t)
    IF(dit==1) dconv_l(:) =        l_e(:)  * d_f4(t)
    IF(dIG==1) dconv_k(:) =   (1.0+k_e(:)) * d_f4(t)
!
!
!
    IF (only_elastic==0) THEN
!
    do 2 n=0, NR
!
!
       IF(t>= -n*teta) then
!
         do l=lmin, lmax
         do m=1,nroots
            aux= (1./tauc+1./dinc)* &
                   ( (1+0*teta) + (0.0-EXP(-(t+n*teta)/tek(l,m)) ) )
            IF(dIR==1) dconv_h(l) = dconv_h(l) + r_h(l,m)* AUX
            IF(dit==1) dconv_l(l) = dconv_l(l) + r_l(l,m)* AUX
            IF(dIG==1) dconv_k(l) = dconv_k(l) + r_k(l,m)* AUX
         enddo
         enddo
!
       ENDIF
!
       IF(t>= -n*teta-tauc) then
!
         do l=lmin, lmax
         do m=1,nroots
            aux= 1./tauc + &
                        (0.0-EXP( -(t+n*teta+tauc)/tek(l,m)))/(tauc/1.0)
            IF(dIR==1) dconv_h(l) = dconv_h(l) - r_h(l,m)* AUX
            IF(dit==1) dconv_l(l) = dconv_l(l) - r_l(l,m)* AUX
            IF(dIG==1) dconv_k(l) = dconv_k(l) - r_k(l,m)* AUX
         enddo
         enddo
!
       ENDIF
!
       IF(t>= -n*teta+dinc) then
!
         do l=lmin, lmax
         do m=1,nroots
            aux= 1./dinc + &
                        (0.0-EXP( -(t+n*teta-dinc)/tek(l,m)))/(dinc/1.0)
            IF(dIR==1) dconv_h(l) = dconv_h(l) - r_h(l,m)* AUX
            IF(dit==1) dconv_l(l) = dconv_l(l) - r_l(l,m)* AUX
            IF(dIG==1) dconv_k(l) = dconv_k(l) - r_k(l,m)* AUX
         enddo
         enddo
!
       ENDIF
!
!
 2  CONTINUE
!
    ENDIF
!
!
  RETURN
!
!
end subroutine  Convol_4_new
!
!
!
!
!
!
function th_5 (t)
!
!  (periodic)
!
 use COMMON; implicit NONE
 REAL(sp) :: th_5, f5, d_f5
!
 REAL(sp) :: t       ! local time
!
 real(sp), parameter :: pig    = 3.1415926535897932384
!
!+----------------------------------+
! Definition of the Th#5:  (periodic)
!+----------------------------------+
!
ENTRY f5 (t)
!
F5  =  0.5*(1.0+SIN(2*pig*t/period))
!
RETURN
!
!+------------------------------------+
! Derivative of the Th#5:  (saw- tooth)
!+------------------------------------+
!
ENTRY d_f5 (t)
!
d_f5 = -0.5*(2.*pig/period)*COS(2*pig*t/period)
!
RETURN
!
end function th_5
!
!
!

!
subroutine Convol_5_new (t)
!
! (periodic)
!
 use COMMON; implicit NONE
!
 REAL(sp) :: t          ! << local time >> variable
 REAL(sp) :: AUX        ! Auxilium
 REAL(sp) :: f5, d_f5   ! time-history and its derivative
 REAL(sp) :: omega      ! 'frequency' = 2*pi/period
 real(sp), parameter :: pig    = 3.1415926535897932384
!
 INTEGER (i4b) :: m  ! Modes index
 INTEGER (i4b) :: L  ! Harmonic Degree
!
!
!+------------+
! Th#5 (*) ldcs
!+------------+
!
ENTRY cnv_5 (t)
!
  omega=2.*pig/period
!
  conv_h(:) = 0.0
  conv_l(:) = 0.0
  conv_k(:) = 0.0
!
    IF(IR==1) conv_h(:) =        h_e(:)  * f5(t)
    IF(it==1) conv_l(:) =        l_e(:)  * f5(t)
    IF(IG==1) conv_k(:) =   (1.0+k_e(:)) * f5(t)
!
  IF(only_elastic==0) THEN
!
    IF(IR==1) conv_h(:) =  conv_h(:) +  0.5* ( h_f(:) - h_e(:) )
    IF(it==1) conv_l(:) =  conv_l(:) +  0.5* ( l_f(:) - l_e(:) )
    IF(IG==1) conv_k(:) =  conv_k(:) +  0.5* ( k_f(:) - k_e(:) )
!
         do l=lmin, lmax
         do m=1,nroots
            aux= (1./(1.+omega**2*tek(l,m)**2))* &
                 ( -omega*tek(l,m)*COS(omega*t) + SIN(omega*t) )
!
            IF(IR==1) conv_h(l) = conv_h(l) - 0.5*r_h(l,m)* AUX
            IF(it==1) conv_l(l) = conv_l(l) - 0.5*r_l(l,m)* AUX
            IF(IG==1) conv_k(l) = conv_k(l) - 0.5*r_k(l,m)* AUX
         enddo
         enddo
!
  ENDIF
!
  RETURN
!
!
!
!
!
!+-------------------+
! (d/dt) Th#5 (*) ldcs
!+-------------------+
!
ENTRY dcnv_5 (t)
!
  dconv_h(:) = 0.0
  dconv_l(:) = 0.0
  dconv_k(:) = 0.0
!
    IF(dIR==1) dconv_h(:) =        h_e(:)  * d_f5(t)
    IF(dit==1) dconv_l(:) =        l_e(:)  * d_f5(t)
    IF(dIG==1) dconv_k(:) =   (1.0+k_e(:)) * d_f5(t)
!
    IF (only_elastic==0) THEN
!
         do l=lmin, lmax
         do m=1,nroots
!
            aux = (1./(1.+omega**2*tek(l,m)**2))* &
            ( +omega**2*tek(l,m)*COS(omega*t) + omega*SIN(omega*t) )
!
            IF(dIR==1) dconv_h(l) = dconv_h(l) - 0.5*r_h(l,m)* AUX
            IF(dit==1) dconv_l(l) = dconv_l(l) - 0.5*r_l(l,m)* AUX
            IF(dIG==1) dconv_k(l) = dconv_k(l) - 0.5*r_k(l,m)* AUX
!
         enddo
         enddo
!
    ENDIF
!
!
  RETURN
!
!
end subroutine  Convol_5_new
!
!
!
!
!
function th_6 (t)
!
!  (periodic)
!
 use COMMON; implicit NONE
 REAL(sp) :: th_6, f6, d_f6
!
 REAL(sp) :: t       ! local time
!
 real(sp), parameter :: pig    = 3.1415926535897932384
!
 INTEGER(i4b) :: k   ! index
!
!
!+------------------------------------------+
! Definition of the Th#6:  (piecewise linear)
!+------------------------------------------+
!
!
ENTRY f6 (t)
!
!
  F6 = 0.0
!
!
  F6 = a_hh(0)
!
  do k=0, nss
    IF(t >= t_hh(k)) F6 = F6 + (alf(k)+bet(k)*t)
  end do
!
RETURN
!
!+------------------------------------------+
! Derivative of the Th#6:  (piecewise linear)
!+------------------------------------------+
!
!
ENTRY d_f6 (t)
!
!
  D_F6 = 0.0
!
  do k=0, nss
    IF(t >= t_hh(k)) D_F6 = D_F6 + bet(k)
  end do
!
!
RETURN
!
end function th_6
!
!
!
!
!
!
!
!
!
subroutine Convol_6_new (t)
!
! (piecewise linear)
!
 use COMMON; implicit NONE
!
 REAL(sp) :: t          ! << local time >> variable
 REAL(sp) :: AUX        ! Auxilium
 REAL(sp) :: f6, d_f6   ! time-history and its derivative
 real(sp), parameter :: pig    = 3.1415926535897932384
!
 INTEGER (i4b) :: m  ! Modes index
 INTEGER (i4b) :: L  ! Harmonic Degree
 INTEGER (i4b) :: k  ! Simply an index
!
!
!
!+------------+
! Th#6 (*) ldcs
!+------------+
!
ENTRY cnv_6 (t)
!
  conv_h(:) = 0.0
  conv_l(:) = 0.0
  conv_k(:) = 0.0
!
    IF(IR==1) conv_h(:) =        h_e(:)  * f6(t)
    IF(it==1) conv_l(:) =        l_e(:)  * f6(t)
    IF(IG==1) conv_k(:) =   (1.0+k_e(:)) * f6(t)
!
  IF(only_elastic==0) THEN
!
    IF(IR==1) conv_h(:) =  conv_h(:) +  a_hh(0)* ( h_f(:) - h_e(:) )
    IF(it==1) conv_l(:) =  conv_l(:) +  a_hh(0)* ( l_f(:) - l_e(:) )
    IF(IG==1) conv_k(:) =  conv_k(:) +  a_hh(0)* ( k_f(:) - k_e(:) )
!
    do 1 k=0, NSS
!
    IF(t >= t_hh(k)) THEN
!
         do l=lmin, lmax
         do m=1,nroots
!
            AUX= (alf(k) + bet(k)*(t_hh(k) - tek(l,m)))* &
                  ( EXP(-(t-t_hh(k))/tek(l,m))  -1.0 ) - bet(k)*(t-t_hh(k))
!
            IF(IR==1) conv_h(l) = conv_h(l) + r_h(l,m)* AUX
            IF(it==1) conv_l(l) = conv_l(l) + r_l(l,m)* AUX
            IF(IG==1) conv_k(l) = conv_k(l) + r_k(l,m)* AUX
         enddo
         enddo
!
   ENDIF
!
 1 CONTINUE
!
    ENDIF
!
RETURN
!
!
!
!
!+---------------------+
! d/dt of Th#6 (*) ldcs
!+---------------------+
!
ENTRY dcnv_6 (t)
!
  dconv_h(:) = 0.0
  dconv_l(:) = 0.0
  dconv_k(:) = 0.0
!
    IF(dIR==1) dconv_h(:) =        h_e(:)  * d_f6(t)
    IF(dit==1) dconv_l(:) =        l_e(:)  * d_f6(t)
    IF(dIG==1) dconv_k(:) =   (1.0+k_e(:)) * d_f6(t)
!
  IF(only_elastic==0) THEN
!
    do 2 k=0, NSS
!
    IF(t >= t_hh(k)) THEN
!
         do l=lmin, lmax
         do m=1,nroots
!
            AUX= (alf(k) + bet(k)*(t_hh(k) - tek(l,m)))* &
                  ( - EXP(-(t-t_hh(k))/tek(l,m))/tek(l,m)  -0.0 ) - bet(k)
!
            IF(dIR==1) dconv_h(l) = dconv_h(l) + r_h(l,m)* AUX
            IF(dit==1) dconv_l(l) = dconv_l(l) + r_l(l,m)* AUX
            IF(dIG==1) dconv_k(l) = dconv_k(l) + r_k(l,m)* AUX
         enddo
         enddo
!
    ENDIF
!
 2  CONTINUE
!
    ENDIF
!
!
RETURN
!
!
end subroutine  Convol_6_new
!



!
!
function th_7 (t)
!
!  (piecewise constant)
!
 use COMMON; implicit NONE
 REAL(sp) :: th_7, f7, d_f7
!
 REAL(sp) :: t       ! local time
!
 real(sp), parameter :: pig    = 3.1415926535897932384
!
 INTEGER(i4b) :: k   ! index
!
!
!+--------------------------------------------+
! Definition of the Th#7:  (piecewise constant)
!+--------------------------------------------+
!
!
ENTRY F7 (t)
!
!
!
  F7 = a_h(0)
!
  do k=0, ns
    IF(t >= float(k)*dilta) F7 = F7 + (a_h(k+1)-a_h(k))
  end do
!
RETURN
!
!+--------------------------------------------+
! Derivative of the Th#7:  (piecewise constant)
!+--------------------------------------------+
!
!
ENTRY d_f7 (t)
!
!
  D_F7 = 0.0
!
!
RETURN
!
end function th_7
!
!

!
!
subroutine Convol_7_new (t)
!
! (piecewise constant)
!
 use COMMON; implicit NONE
!
 REAL(sp) :: t          ! << local time >> variable
 REAL(sp) :: AUX        ! Auxilium
 REAL(sp) :: f7, d_f7   ! time-history and its derivative
 real(sp), parameter :: pig    = 3.1415926535897932384
!
 INTEGER (i4b) :: m  ! Modes index
 INTEGER (i4b) :: L  ! Harmonic Degree
 INTEGER (i4b) :: k  ! Simply an index
!
!
!
!+------------+
! Th#7 (*) ldcs
!+------------+
!
ENTRY cnv_7 (t)
!
  conv_h(:) = 0.0
  conv_l(:) = 0.0
  conv_k(:) = 0.0
!
    IF(IR==1) conv_h(:) =        h_e(:)  * f7(t)
    IF(it==1) conv_l(:) =        l_e(:)  * f7(t)
    IF(IG==1) conv_k(:) =   (1.0+k_e(:)) * f7(t)
!
  IF(only_elastic==0) THEN
!
    IF(IR==1) conv_h(:) =  conv_h(:) +  a_h(0)* ( h_f(:) - h_e(:) )
    IF(it==1) conv_l(:) =  conv_l(:) +  a_h(0)* ( l_f(:) - l_e(:) )
    IF(IG==1) conv_k(:) =  conv_k(:) +  a_h(0)* ( k_f(:) - k_e(:) )
!
    do 1 k=0, NS
!
    IF(t >= float(k)*dilta) THEN
!
         do l=lmin, lmax
         do m=1,nroots
!
            AUX=  (1.0 - EXP(-(t-k*dilta)/tek(l,m))) * (a_h(k+1)-a_h(k))
!
            IF(IR==1) conv_h(l) = conv_h(l) - r_h(l,m)* AUX
            IF(it==1) conv_l(l) = conv_l(l) - r_l(l,m)* AUX
            IF(IG==1) conv_k(l) = conv_k(l) - r_k(l,m)* AUX
         enddo
         enddo
!
   ENDIF
!
 1 CONTINUE
!
    ENDIF
!
RETURN
!
!
!
!
!+---------------------+
! d/dt of Th#7 (*) ldcs
!+---------------------+
!
ENTRY dcnv_7 (t)
!
  dconv_h(:) = 0.0
  dconv_l(:) = 0.0
  dconv_k(:) = 0.0
!
    IF(dIR==1) dconv_h(:) =        h_e(:)  * d_f7(t)
    IF(dit==1) dconv_l(:) =        l_e(:)  * d_f7(t)
    IF(dIG==1) dconv_k(:) =   (1.0+k_e(:)) * d_f7(t)
!
  IF(only_elastic==0) THEN
!
    do 2 k=0, NS
!
    IF(t >= float(k)*dilta) THEN
!
         do l=lmin, lmax
         do m=1,nroots
!
            AUX=  (0.0 + EXP(-(t-k*dilta)/tek(l,m))/tek(l,m)) * (a_h(k+1)-a_h(k))
!
            IF(dIR==1) dconv_h(l) = dconv_h(l) - r_h(l,m)* AUX
            IF(dit==1) dconv_l(l) = dconv_l(l) - r_l(l,m)* AUX
            IF(dIG==1) dconv_k(l) = dconv_k(l) - r_k(l,m)* AUX
         enddo
         enddo
!
   ENDIF
!
 2 CONTINUE
!
    ENDIF
!
RETURN
!
!
!
end subroutine  Convol_7_new
!


!
!
function th_8 (t)
!
!  (piecewise constant)
!
 use COMMON; implicit NONE
 REAL(sp) :: th_8, f8, d_f8
 REAL(sp) ::       f7, d_f7
!
 REAL(sp) :: t         ! local time
!
 real(sp), parameter :: pig    = 3.1415926535897932384
 REAL(sp) :: AUX, daux ! auxilia
!
 INTEGER(i4b) :: k     ! index
!
!
!+-----------------------------------------------------------------------+
! Definition of the Th#8:  (piecewise constant with finite loading phase)
!+-----------------------------------------------------------------------+
!
!
ENTRY F8 (t)
!
  F8 = 0.0
!
!
  IF(t >=  0.0) F8 = F7(t)
!
  IF(t < -tauf) F8 = 0.0
!
  IF(-tauf <= t .AND. t < 0.0) F8 = a_h(0)*(1.0+t/tauf)
!
RETURN
!
!+----------------------------------------------------------------------+
! Derivative of the Th#8:  (piecewise constant with finite loading phase)
!+----------------------------------------------------------------------+
!
!
ENTRY d_f8 (t)
!
  D_F8 = 0.0
!
!
  IF(t >=  0.0) D_F8 = D_F7(t)
!
  IF(t < -tauf) D_F8 = 0.0
!
  IF(-tauf <= t .AND. t < 0.0) D_F8 = a_h(0)/tauf
!
RETURN
!
end function th_8
!
!
!
!
!
!
subroutine Convol_8_new (t)
!
! (piecewise constant with finite loading phase)
!
 use COMMON; implicit NONE
!
 REAL(sp) :: t          ! << local time >> variable
 REAL(sp) :: AUX        ! Auxilium
 REAL(sp) :: f8, d_f8   ! time-history and its derivative
 REAL(sp) :: f7, d_f7   ! also these are needed
 real(sp), parameter :: pig    = 3.1415926535897932384
!
 INTEGER (i4b) :: m  ! Modes index
 INTEGER (i4b) :: L  ! Harmonic Degree
 INTEGER (i4b) :: k  ! Simply an index
!
!
!
!+------------+
! Th#8 (*) ldcs
!+------------+
!
ENTRY cnv_8 (t)
!
!
    conv_h(:) =  0.0
    conv_l(:) =  0.0
    conv_k(:) =  0.0
!
!
!
    call cnv_7(t)
!
!
    IF(IR==1) conv_h(:) =   conv_h(:) +      h_e(:)  * (f8(t)-f7(t))
    IF(it==1) conv_l(:) =   conv_l(:) +      l_e(:)  * (f8(t)-f7(t))
    IF(IG==1) conv_k(:) =   conv_k(:) + (1.0+k_e(:)) * (f8(t)-f7(t))
!
  IF(only_elastic==0) THEN
!
    IF(IR==1) conv_h(:) =  conv_h(:) +  a_h(0)* ( h_e(:) - h_f(:) )
    IF(it==1) conv_l(:) =  conv_l(:) +  a_h(0)* ( l_e(:) - l_f(:) )
    IF(IG==1) conv_k(:) =  conv_k(:) +  a_h(0)* ( k_e(:) - k_f(:) )
!
!
    IF(t >= 0.0) THEN
!
         do l=lmin, lmax
         do m=1,nroots
!
            AUX=  a_h(0)*(t/tauf - (1.0 -EXP(-t/tek(l,m)) )*tek(l,m)/tauf)
!
            IF(IR==1) conv_h(l) = conv_h(l) + r_h(l,m)* AUX
            IF(it==1) conv_l(l) = conv_l(l) + r_l(l,m)* AUX
            IF(IG==1) conv_k(l) = conv_k(l) + r_k(l,m)* AUX
         enddo
         enddo
!
   ENDIF
!
    IF(t >= -tauf) THEN
!
         do l=lmin, lmax
         do m=1,nroots
!
            AUX=  a_h(0)* ((1.0 + t/tauf) - &
                          ( 1.0 -EXP(-(t+tauf)/tek(l,m)) )*tek(l,m)/tauf )
!
            IF(IR==1) conv_h(l) = conv_h(l) - r_h(l,m)* AUX
            IF(it==1) conv_l(l) = conv_l(l) - r_l(l,m)* AUX
            IF(IG==1) conv_k(l) = conv_k(l) - r_k(l,m)* AUX
         enddo
         enddo
!
   ENDIF
!
   ENDIF
!
RETURN
!
!
!
!
!+---------------------+
! d/dt of Th#8 (*) ldcs
!+---------------------+
!
ENTRY dcnv_8 (t)
!
!
    call dcnv_7(t)
!
!
    IF(dIR==1) dconv_h(:) =   dconv_h(:) +      h_e(:)  * (d_f8(t)-d_f7(t))
    IF(dit==1) dconv_l(:) =   dconv_l(:) +      l_e(:)  * (d_f8(t)-d_f7(t))
    IF(dIG==1) dconv_k(:) =   dconv_k(:) + (1.0+k_e(:)) * (d_f8(t)-d_f7(t))
!
  IF(only_elastic==0) THEN
!
!
    IF(t >= 0.0) THEN
!
         do l=lmin, lmax
         do m=1,nroots
!
            AUX=  a_h(0)*(1./tauf + (0.0 - EXP(-t/tek(l,m)) )/tauf)
!
            IF(dIR==1) dconv_h(l) = dconv_h(l) + r_h(l,m)* AUX
            IF(dit==1) dconv_l(l) = dconv_l(l) + r_l(l,m)* AUX
            IF(dIG==1) dconv_k(l) = dconv_k(l) + r_k(l,m)* AUX
         enddo
         enddo
!
   ENDIF
!
    IF(t >= -tauf) THEN
!
         do l=lmin, lmax
         do m=1,nroots
!
            AUX=  a_h(0)* ((0.0 + 1./tauf) - &
                          ( 0.0 -EXP(-(t+tauf)/tek(l,m)) )/tauf )
!
            IF(dIR==1) dconv_h(l) = dconv_h(l) - r_h(l,m)* AUX
            IF(dit==1) dconv_l(l) = dconv_l(l) - r_l(l,m)* AUX
            IF(dIG==1) dconv_k(l) = dconv_k(l) - r_k(l,m)* AUX
         enddo
         enddo
!
   ENDIF
!
   ENDIF
!
RETURN
!
!
!
end subroutine  Convol_8_new
!
!
!
!
!
!
! +----------------------------------------------+
      SUBROUTINE CONVOL_load_history   (HI, t)   !
! +----------------------------------------------+
!
 use COMMON; implicit NONE
!
 real (sp) :: t
 integer (i4b) :: hi
!
!
! Input 
!          t  =  time (kyr)
!          Hi =  Time-history label (0 <= hi <= 8)
! 
! Output   
!          conv_h,  conv_l,  conv_k  = convolutions at time t
!          dconv_h, dconv_l, dconv_k = their time-derivatives
!
!
 If ((hi<0 .or. hi>8) .and. iv==0) then
  Write(99,*) 'ERROR IN SBR. CONVOL_load_history: Unknown  time history '
  Write(99,*) 'Only histories with label in the range [0:5] are allowed '
  Write(99,*) '******** JOB ABORTED *********************************** ';STOP
 Endif
 If ((hi<0 .or. hi>8) .and. iv==1) then
  Write(99,*) 'ERROR IN SBR. CONVOL_load_history: Unknown  time history '
  Write(99,*) 'Only histories with label in the range [0:8] are allowed '
  Write(99,*) '******** JOB ABORTED *********************************** '
  Write(*,*)  'ERROR IN SBR. CONVOL_load_history: Unknown  time history '
  Write(*,*)  'Only histories with label in the range [0:8] are allowed '
  Write(*,*)  '******** JOB ABORTED *********************************** ';STOP
 Endif
!
!
 IF (HI == 0) call  cnv_0 (t)
 IF (HI == 0) call dcnv_0 (t)
!
 IF (HI == 1) call  cnv_1 (t)
 IF (HI == 1) call dcnv_1 (t)
!
 IF (HI == 2) call  cnv_2 (t)
 IF (HI == 2) call dcnv_2 (t)
!
 IF (HI == 3) call  cnv_3 (t)
 IF (HI == 3) call dcnv_3 (t)
!
 IF (HI == 4) call  cnv_4 (t)
 IF (HI == 4) call dcnv_4 (t)
!
 IF (HI == 5) call  cnv_5 (t)
 IF (HI == 5) call dcnv_5 (t)
!
 IF (HI == 6) call  cnv_6 (t)
 IF (HI == 6) call dcnv_6 (t)
!
 IF (HI == 7) call  cnv_7 (t)
 IF (HI == 7) call dcnv_7 (t)
!
 IF (HI == 8) call  cnv_8 (t)
 IF (HI == 8) call dcnv_8 (t)
!
!
END SUBROUTINE CONVOL_load_history 
!
!
!
!
!
!
!
! -----------------------------------------------------------------+
!							           |
 SUBROUTINE AXIS_DISP0 (IK, observer, load_center, VECT, D_VECT)   !
!					                           |
! -----------------------------------------------------------------+ 
!
!
 use COMMON; use TRIGFCN; implicit NONE   
!
!
! # This routine computes the components of vector VECT ==
!   (u_rad, u_teta, u_long, geoid), and its time-derivative
!   (du_rad, du_teta, du_long, e d_geoid) for axissimmetric
!   loads of type 10, 11, 20, 21, 30, o 40, in a point of
!   coordinates OBSERVER == (long., colat.), assuming that
!   the axis of simmetry of the load intersects the surface
!   of the Earth at point LOAD_CENTER = (long., colat.)
!
!
! # Fixed 'special case 5' (observer and load at the antipodes)
!   (DM Feb 20, 2017)
!

!
 real(dp), parameter :: eps = 1.D-200    ! Tolerance for antipode check
!
 real(dp) ::   vect(4)    ! == ( u_rad,  u_teta,  u_long,   geoid)
 real(dp) :: d_vect(4)    ! == (du_rad, du_teta, du_long, d_geoid)
!
 real(dp) ::   corr(4)    ! Ocean correction for vect
 real(dp) :: d_corr(4)    ! Ocean correction for d_vect
!
 real(dp) :: observer(2)       ! Long. & Colat. of the observer
 real(dp) :: long_ob, teta_ob  ! Long. & Colat. of the observer
 
 real(dp) :: load_center(2)    ! Long. & Colat. of the load center
 real(dp) :: long_ce, teta_ce  ! Long. & Colat. of the load center
!
 real(dp) ::  u,  v,  h   ! Components of the displacement along the vertical (u),
!                           colatitude (v) and geoid (h) in the reference frame
!                           of the load.
!
 real(dp) :: du, dv, dh   ! Derivatives of u, v, and h
!		           
 real(dp) :: leg, legdev1, legfun, legder     ! Legendre function & derivative
 real(dp) :: gamma, aux                       ! Auxiliary variables
!
 real(dp) :: cos_teta, sin_teta, sinx, cosx   ! Angles for the trigonometric
!                                               analysis
!
 integer (i4b) :: l, m                ! do-loops control variables
!
 integer (i4b) :: IK                  ! Ocean correction switch
!
 logical log_1, log_2, log_3, log_4   ! Auxiliary logical variables
!
!
!
! # Executable statements
!
 teta_ob   =    observer(2); long_ob  =    observer(1)  
 teta_ce   = load_center(2); long_ce  = load_center(1)  
!
 If(Long_Ob<0.D0 .or. Long_Ob>360d0) Then 
 Write(99,*) ' ERROR from Sbr. AXIS_DISP0: Longitude of the OBSERVER out of Range '
 Write(99,*) ' #### JOB ABORTED ################################################# '; Stop 
 Endif 
!
 If(Teta_Ob<0.D0 .or. Teta_Ob>180d0) Then 
 Write(99,*) ' ERROR from Sbr. AXIS_DISP0: Colatitude of the OBSERVER out of Range '
 Write(99,*) ' #### JOB ABORTED ################################################## '; Stop 
 Endif   
!
 If(Long_Ce<0.D0 .or. Long_Ce>360d0) Then 
 Write(99,*) ' ERROR from Sbr. AXIS_DISP0: Longitude of the LOAD CENTER out of Range '
 Write(99,*) ' #### JOB ABORTED #################################################### '; Stop 
 Endif 
!
 If(Teta_Ce<0.D0 .or. Teta_Ce>180d0) Then 
 Write(99,*) ' ERROR from Sbr. AXIS_DISP0: Colatitude of the LOAD CENTER out of Range '
 Write(99,*) ' #### JOB ABORTED ##################################################### '; Stop 
 Endif   
!
!
 vect(:) = 0.d0; d_vect(:) = 0.d0 
 corr(:) = 0.d0; d_corr(:) = 0.d0  
!
!
!
! # Cosine & Sine of the observer colatitude wrt the load center
!
 cos_teta = cosd(teta_ob)*cosd(teta_ce) + 			  & 
            sind(teta_ob)*sind(teta_ce)*cosd(long_ob-long_ce) 
!
 sin_teta = sqrt(1._dp - cos_teta*cos_teta) 
!
!
!
! -----------------------------------------+				         		
! #   Loop on the Legendre  polynomials    | 
! -----------------------------------------+                                         
!
  u = 0._dp;  v = 0._dp;  h = 0._dp 
 du = 0._dp; dv = 0._dp; dh = 0._dp
!
!
   DO 100 L = LMIN, LMAX     
!
!		
gamma = (3._dp/rhoea)*sigma(l)/(2._dp*l+1._dp) 
!
!
IF(IR==1.or.IG==1.or.dIR==1.or.dIG==1) legfun = gamma * leg     (l,0,cos_teta)
IF(it==1.or.dit==1)                    legder = gamma * legdev1 (l,0,cos_teta)
!
IF(IR==1)   u  = u  +  conv_h(l) * legfun 
IF(it==1)   v  = v  +  conv_l(l) * legder
IF(IG==1)   h  = h  +  conv_k(l) * legfun
!
IF(dIR==1) du = du + dconv_h(l) * legfun
IF(dit==1) dv = dv + dconv_l(l) * legder
IF(dIG==1) dh = dh + dconv_k(l) * legfun
!
!
  100 CONTINUE
!
!
!
!
! +---------------+
! | Special cases |
! +---------------+
!
!
!
! 1) Load to the pole, observer not to the pole
!
log_1 = (teta_ce  == 0._dp .or. teta_ce  == 180._dp)
log_2 =  teta_ob  /= 0._dp 
! 
if (log_1 .and. log_2) then 
  vect(1) =  u;   vect(2)=  v;   vect(3)=  0.D0;   vect(4)=  h 
d_vect(1) = du; d_vect(2)= dv; d_vect(3)=  0.D0; d_vect(4)= dh
goto 123
                      endif
!		       		       
!
! 2) Load & Observer at one of the poles.
!
log_1 = teta_ce == 0._dp   .AND. teta_ob ==   0._dp
log_2 = teta_ce == 0._dp   .AND. teta_ob == 180._dp
log_3 = teta_ce == 180._dp .AND. teta_ob ==   0._dp
log_4 = teta_ce == 180._dp .AND. teta_ob == 180._dp
!
if(log_1.or.log_2.or.log_3.or.log_4) then                                    
  vect(1) =  u ;   vect(2)=  0.D0;   vect(3)=  0.D0;   vect(4)=  h 
d_vect(1) = du ; d_vect(2)=  0.D0; d_vect(3)=  0.D0; d_vect(4)= dh 
goto 123
                                    endif 
!
! 
! 3) Load not at the pole, observer at one of the poles
!
!
! The radial displacement & geoid are computed normally,
! whereas u_teta & u_long are set to zero conventionally
! (the direction of the colatitude & longitude unit vectors
! is not defined at the poles)
!
log_1 = teta_ce /=    0._dp .and. teta_ob ==    0._dp 
log_2 = teta_ce /=  180._dp .and. teta_ob ==    0._dp 
log_3 = teta_ce /=    0._dp .and. teta_ob ==  180._dp 
log_4 = teta_ce /=  180._dp .and. teta_ob ==  180._dp 
!
if(log_1.or.log_2.or.log_3.or.log_4) then                                    
  vect(1) =  u ;   vect(2)=  0.D0;   vect(3)=  0.D0;   vect(4)=  h 
d_vect(1) = du ; d_vect(2)=  0.D0; d_vect(3)=  0.D0; d_vect(4)= dh 
goto 123
                                    endif 
!
!
! 4) Load not at the pole, observer at the centre of the load
!
log_1 =  teta_ce /=   0._dp .and. teta_ob == teta_ce   
log_2 =  teta_ce /= 180._dp .and. teta_ob == teta_ce 
log_3 =  long_ce == long_ob
!
if((log_1  .and. log_3) .or. (log_2 .and. log_3)) then 
  vect(1) =  u ;   vect(2)=  0.D0;   vect(3)=  0.D0;   vect(4)=  h 
d_vect(1) = du ; d_vect(2)=  0.D0; d_vect(3)=  0.D0; d_vect(4)= dh 
goto 123
                                                 endif
!
!
! 5) Load not at the pole, observer at its antipodes
!
!log_1 = abs(long_ce - long_ob) == 180._dp 
!log_2 =     teta_ce + teta_ob  == 180._dp
!
!if( log_1 .or. log_2 ) then  
!
log_1 = (cos_teta + 1.d0) < eps
!
if( log_1 ) then
  vect(1) =  u ;   vect(2)=  0.D0;   vect(3)=  0.D0;   vect(4)=  h 
d_vect(1) = du ; d_vect(2)=  0.D0; d_vect(3)=  0.D0; d_vect(4)= dh 
goto 123
                      endif 
!
!
!+================+
!| General case   |
!+================+
!
!
if( sin_teta /= 0._dp .and. teta_ce /= 0._dp) Then
!
 sinx =  sind(long_ob-long_ce)*sind(teta_ce)/sin_teta 
!
 cosx = (cosd(teta_ce)-cosd(teta_ob)*cos_teta)/sin_teta/sind(teta_ob) 
!
   vect(1) =  u;   vect(2) =  v*cosx;   vect(3) =  v*sinx;   vect(4) = h 
 d_vect(1) = du; d_vect(2) = dv*cosx; d_vect(3) = dv*sinx; d_vect(4) = dh   
!
                                              Endif
!
!+--------------------------------------------------+
! Correction for compensation on a realistic ocean  |
!+--------------------------------------------------+
!
 123 CONTINUE
!
     If     (IK == 0) RETURN
!
     IF     (IK == 1) Then
         Call OCEAN_CORRECTION (OBSERVER, CORR, D_CORR)
           vect(:) =   vect(:) +   corr(:)
         d_vect(:) = d_vect(:) + d_corr(:)
                     Endif
!
!
!
!
!
END SUBROUTINE AXIS_DISP0
!
!
!
!
!
!
! -------------------------------------------------------------+
!		              			    	       |
 SUBROUTINE AXIS_Pressure (observer, load_center, TI_HI, PP)   !
!					                       |
! -------------------------------------------------------------+
!
 use COMMON; use TRIGFCN; implicit NONE   
!
! # This routine computes the surface load PP (mass per unit surface)
!   for axissimmetric loads of the type 10, 11, 20, and 21 at the point
!   of coordinates OBSERVER = (long., colat.). TH is the value of the
!   time-history at current time.
!
!!
!
 real(dp) :: observer(2)       ! Observer Longitude & Colatitude
 real(dp) :: long_ob, teta_ob  ! Observer Longitude & Colatitude
 
 real(dp) :: load_center(2)    ! Long. & Colat. of the load center
 real(dp) :: long_ce, teta_ce  ! Long. & Colat. of the load center
!
 real(dp) :: PP                ! Load at point "Observer" (kg/m**2)
!
 real(dp) :: TI_HI             ! Value of the time-history
 real(dp) :: leg               ! Legendre polynomial
 real(dp) :: cos_teta          ! Cosine of the colatitude of the observer
!                              ! wrt the center pf the load.
!
 integer (i4b) :: l, m         ! do-loop controla variables
!
!
!
! +------------------------+
! | Executable statements  |
! +------------------------+
!
 teta_ob   =    observer(2); long_ob  =    observer(1)  
 teta_ce   = load_center(2); long_ce  = load_center(1)  
!
 If(Long_Ob<0.D0 .or. Long_Ob>360d0) Then 
 Write(99,*) ' ERROR from Sbr. AXIS_DISP0: Longitude of the OBSERVER out of Range '
 Write(99,*) ' #### JOB ABORTED ################################################# '; Stop 
 Endif 
!
 If(Teta_Ob<0.D0 .or. Teta_Ob>180d0) Then 
 Write(99,*) ' ERROR from Sbr. AXIS_DISP0: Colatitude of the OBSERVER out of Range '
 Write(99,*) ' #### JOB ABORTED ################################################## '; Stop 
 Endif   
!
 If(Long_Ce<0.D0 .or. Long_Ce>360d0) Then 
 Write(99,*) ' ERROR from Sbr. AXIS_DISP0: Longitude of the LOAD CENTER out of Range '
 Write(99,*) ' #### JOB ABORTED #################################################### '; Stop 
 Endif 
!
 If(Teta_Ce<0.D0 .or. Teta_Ce>180d0) Then 
 Write(99,*) ' ERROR from Sbr. AXIS_DISP0: Colatitude of the LOAD CENTER out of Range '
 Write(99,*) ' #### JOB ABORTED ##################################################### '; Stop 
 Endif   
!
!
! # Cosine of the colatitude of the observer wrt the center pf the load.
!
 cos_teta = cosd(teta_ob)*cosd(teta_ce) + 			  & 
            sind(teta_ob)*sind(teta_ce)*cosd(long_ob-long_ce) 
!
! -----------------------------------------+				         		
! #   Loop on the Legendre  polynomials    | 
! -----------------------------------------+
!
    PP = 0.D0
! 
    DO L = 0, LMAX
!
    PP = PP + sigma(l)*leg (l,0,cos_teta)
!
    Enddo 
!      
!
    PP = PP * TI_HI
!
!
END SUBROUTINE AXIS_Pressure
!
!
!
!

!
!
!----------------------------------------------------+						                
       SUBROUTINE RECT_PRESSURE (observer, th, PP)   !
!----------------------------------------------------+
!
  use COMMON; use TRIGFCN; implicit NONE   
!
!
! # This routine computes the surface load PP (mass per unit surface)
!   for non-axissimmetric loads of the type 50 at the  point of
!   coordinates OBSERVER = (long., colat.). TH is the value of the
!   time-history at current time.
!
real(dp):: th   ! Time history value at current time
real(dp):: PP   ! Load (kg/m**2) at point 'Observer'
!
real(dp) :: observer (2), long_ob, teta_ob  ! Observer coordinates
!
real(dp) :: leg ! Associated Legendre function
!
!
integer (i4b) :: l, m  ! do-loops control variables
!
!
 teta_ob  =   observer(2)
 long_ob  =   observer(1)  
!
!
 If(Long_Ob<0.D0 .or. Long_Ob>360d0) Then 
 Write(99,*) ' ERROR from Sbr. RECT_DISP0: Longitude of the OBSERVER out of Range '
 Write(99,*) ' #### JOB ABORTED ################################################# '; Stop 
 Endif 
!
 If(Teta_Ob<0.D0 .or. Teta_Ob>180d0) Then 
 Write(99,*) ' ERROR from Sbr. RECT_DISP0: Colatitude of the OBSERVER out of Range '
 Write(99,*) ' #### JOB ABORTED ################################################## '; Stop 
 Endif   
!
!
! # Computed the load (kg/m**2) at the observer position
!
!
        PP = 0.D0 
!
	DO L= 0, Lmax
          DO M=0, L 
!    	
!
           PP = PP + ( FF(L,M)* COSD(M*LONG_OB)   + & 
	               GG(L,M)* SIND(M*LONG_OB) ) * & 
	               LEG (L,M,COSD(TETA_OB))  
!
          ENDDO	         
        ENDDO         
!
!
        PP = PP * th 	
!
!
End Subroutine Rect_Pressure 
!
!
!
!
!
!--------------------------------------------------------------------------+						                
       SUBROUTINE RECT_DISP0 (INDOC, observer, VECT, D_VECT)               !   
!--------------------------------------------------------------------------+
!
!
  use COMMON; use TRIGFCN; implicit NONE   
!
! # This routine computes the components of vector  VECT =
!   (u_rad, u_teta, u_long, geoid) and its time-derivative
!   (du_rad, du_teta, du_long, d_geoid), for loads of type
!   50 at the point of coordinates OBSERVER= (long., colat.).
!!
 real(dp) ::   vect(4)    ! == ( u_rad,  u_teta,  u_long,   geoid)
 real(dp) :: d_vect(4)    ! == (du_rad, du_teta, du_long, d_geoid)
!
 real(dp) :: observer (2), long_ob, teta_ob    ! Observer
!
 real(dp) ::  u_rad    ! Radial displacement
 real(dp) ::  u_teta   ! Displacement along colatitude
 real(dp) ::  u_long   ! Displacement along colatitude
 real(dp) ::  geoid    ! Geoid height
!  
 real(dp) ::  du_rad   ! Time-derivative of u_rad
 real(dp) ::  du_teta  ! Time-derivative of u_teta
 real(dp) ::  du_long  ! Time-derivative of u_long
 real(dp) ::  d_geoid  ! Time-derivative of geoid
!
 real(dp) :: leg, legdev1   ! Associate Legendre function & its derivative
!                           ! wrt teta
!
 real(dp) :: CFA, aux;  logical log_1   ! Auxiliary variables
!
 integer (i4b) :: l, m                ! do-loops control variables
 integer (i4b) :: INDOC               ! OCEAN index (0/1)
 real(dp) :: FS(0:LLMAX), FC(0:LLMAX) ! Sines and Cosines
 real(dp) :: BRA1, BRA2, PLM, PLMD    ! Auxiliary variables
!
!
!
! +-----------------------+
! | Executable statements |
! +-----------------------+
!
 teta_ob  =   observer(2)
 long_ob  =   observer(1)  
!
 If(Long_Ob<0.D0 .or. Long_Ob>360d0) Then 
 Write(99,*) ' ERROR from Sbr. RECT_DISP0: Longitude of the OBSERVER out of Range '
 Write(99,*) ' #### JOB ABORTED ################################################# '; Stop 
 Endif 
!
 If(Teta_Ob<0.D0 .or. Teta_Ob>180d0) Then 
 Write(99,*) ' ERROR from Sbr. RECT_DISP0: Colatitude of the OBSERVER out of Range '
 Write(99,*) ' #### JOB ABORTED ################################################## '; Stop 
 Endif   
!
!
 VECT(:) = 0._DP; D_VECT(:) = 0._DP 
!
!
!
! # Trigonometric functions
!
  DO M=0, LMAX 
   FS(M)=SIND(M*LONG_OB)
   FC(M)=COSD(M*LONG_OB)
  ENDDO  
!
!
!  write(*,*) SIGMA(0)/OC(0,0)
!   read(*,*) 
! 
! # Correction for realistic oceans
!
!
        IF (INDOC == 0) THEN
            FF(:,:) =  FF(:,:)
            GG(:,:) =  GG(:,:)
        ENDIF
!
        IF (INDOC == 1) THEN 
!            FF(:,:) = FF(:,:) - OC(:,:)*SIGMA(0)/OC(0,0)   
!	     GG(:,:) = GG(:,:) - OS(:,:)*SIGMA(0)/OC(0,0) 

             FF(:,:) =  - OC(:,:)*SIGMA(0)/OC(0,0)   
 	     GG(:,:) =  - OS(:,:)*SIGMA(0)/OC(0,0) 

        ENDIF
!
! ---------------------------------------+
! # Loop on the real spherical harmonics |
! ---------------------------------------+
!
          u_rad = 0._dp;  u_teta = 0._dp;  u_long = 0._dp;   geoid = 0._dp
         du_rad = 0._dp; du_teta = 0._dp; du_long = 0._dp; d_geoid = 0._dp 
!
!
	DO L=LMIN, LMAX
!
        CFA = + (3._dp/rhoea)/(2._dp*l+1._dp) 
!
        DO M=0, L 
!
        PLM = LEG (L,M,COSD(TETA_OB)) 
!
!
        If(IR==1.or.dIR==1.or.IG==1.or.dIG==1) THEN 
!
        BRA1= (FF(L,M)*FC(M) + GG(L,M)*FS(M)) * PLM * CFA  
!
        If(IR  ==1)   u_rad  =   u_rad +  conv_h(l)*BRA1
        If(dIR==1)  du_rad  =  du_rad + dconv_h(l)*BRA1
!
        If(IG  ==1)    geoid =   geoid +  conv_k(l)*BRA1 
        If(dIG==1)  d_geoid = d_geoid + dconv_k(l)*BRA1
!
                                                         ENDIF 
!
        If(it==1.or.dit==1) THEN 
!
	PLMD= LEGDEV1(L,M,COSD(TETA_OB))
!
        BRA1= (FF(L,M)*FC(M) + GG(L,M)*FS(M)) * PLMD * CFA 
	BRA2= (GG(L,M)*FC(M) - FF(L,M)*FS(M)) *  PLM * CFA * M  
!
        If(it   == 1 .and. (sind(teta_ob) /= 0.d0))              THEN 
            u_teta =  u_teta   +  conv_l(l)*BRA1 
	    u_long =  u_long   +  conv_l(l)*BRA2/SIND(TETA_OB)
	ENDIF
!
        If(dit == 1 .and. (sind(teta_ob) /= 0.d0))              THEN 
           du_teta = du_teta + dconv_l(l)*BRA1
	   du_long = du_long + dconv_l(l)*BRA2/SIND(TETA_OB)
	ENDIF
!
                                                         ENDIF 
				 	
!
        ENDDO
!	         
        ENDDO         
!
!
!
! * Special case: Observer in one of the two poles
!
 log_1 = teta_ob == 0._dp   .OR. teta_ob ==  180._dp
!
 if(log_1) then                                    
  u_teta = 0._dp ;  u_long = 0._dp 
 du_teta = 0._dp ; du_long = 0._dp   
          endif 
!
!
!
!
!+--------------------------------------+
! Assignements to vectors VECT & D_VECT |
!+--------------------------------------+
!
   vect(1) =  u_rad;    vect(2) =  u_teta
   vect(3) =  u_long;   vect(4) =   geoid  
 d_vect(1) =  du_rad; d_vect(2) = du_teta
 d_vect(3) = du_long; d_vect(4) = d_geoid
!
!
!
end subroutine RECT_DISP0 
!
!
!
!
!
!
!
!							         
! +--------------------------------------------------------+
   SUBROUTINE OCEAN_CORRECTION  (observer, CORR, D_CORR)   !
! +--------------------------------------------------------+
!
use COMMON; use TRIGFCN; implicit NONE   
!
!  # This routine computes the correction to the fields already
!    computed via AXIS_DISP0 or RECT_DISP0 for the load distributed
!    on oceans of realistic shape. CORR = ( cu_rad, cu_teta, cu_long,
!    c_geoid) = solution vector at point OBSERVER, and  D_CORR =
!    (cud_rad, cud_teta, cud_long, cd_geoid) = time derivative of
!    CORR.
!!
!
 real(dp) :: observer(2), teta_ob, long_ob ! Observer coordinates
!
!
 real(dp) ::  u_rad    ! Radial displacement
 real(dp) ::  u_teta   ! Displacement along colatitude
 real(dp) ::  u_long   ! Displacement along colatitude
 real(dp) ::  geoid    ! Geoid height
!  
 real(dp) ::  du_rad   ! Time-derivative of u_rad
 real(dp) ::  du_teta  ! Time-derivative of u_teta
 real(dp) ::  du_long  ! Time-derivative of u_long
 real(dp) ::  d_geoid  ! Time-derivative of geoid
!
 real(dp) ::   CORR(4) ! = (  u_rad,  u_teta,  u_long,   geoid) 
 real(dp) :: D_CORR(4) ! = ( du_rad, du_teta, du_long, d_geoid)
!
 real(dp) :: leg, legdev1   ! Associate Legendre function & its derivative
!                           ! wrt teta
!
 real(dp) :: gamma, aux     ! Auxiliary variables
!
 integer (i4b) :: l, m; logical log_1   ! Auxiliary variables
!
 real(dp) :: FS(0:LLMAX), FC(0:LLMAX)   ! Sine & Cosines
!
 real(dp) :: CFA, BRA1, BRA2, PLM, PLMD ! Auxiliary variables
!
!
!
!
! +-----------------------+
! | Executable statements |
! +-----------------------+
!
 teta_ob = observer(2)
 long_ob = observer(1)  
!
 corr(:)=0.d0; d_corr(:)=0.d0
!
!
 IF(sigma(0)==0._dp) call m_3 (30662)  ! Error, if sigma(0) is zero already
!
!
!
   u_rad = 0._dp;   u_teta = 0._dp;   u_long = 0._dp;   geoid = 0._dp
  du_rad = 0._dp;  du_teta = 0._dp;  du_long = 0._dp; d_geoid = 0._dp 
!
!
! -------------------------------+
! Terms which depend only on 'm' |
! -------------------------------+
!
  do m=0, lmax 
    FC(M) = COSD (M*LONG_OB)
    FS(M) = SIND (M*LONG_OB) 
  enddo 
!
!
!
! --------------------------------------+
! # Sum on the real spherical harmonics |
! --------------------------------------+
!
	DO L=LMIN, LMAX
!
        CFA =   - ((3._dp/rhoea)/(2._dp*l+1._dp))*(SIGMA(0)/OC(0,0))
!
        DO M=0, L 
!
        If(IR==1.or.dIR==1.or.IG==1.or.dIG==1) THEN 
!
        PLM = LEG    (L,M,COSD(TETA_OB)) 
!
        BRA1= (OC(L,M)*FC(M) + OS(L,M)*FS(M)) * PLM * CFA  
!
        If(IR  ==1)   u_rad  =   u_rad +  conv_h(l)*BRA1
        If(dIR==1)  du_rad  =  du_rad + dconv_h(l)*BRA1
!
        If(IG  ==1)    geoid =   geoid +  conv_k(l)*BRA1 
        If(dIG==1)  d_geoid = d_geoid + dconv_k(l)*BRA1
!
                                               ENDIF
!
        If(it==1.or.dit==1) THEN 
!
        PLM = LEG    (L,M,COSD(TETA_OB)) 
	PLMD= LEGDEV1(L,M,COSD(TETA_OB))
!
        BRA1= (OC(L,M)*FC(M) + OS(L,M)*FS(M)) * PLMD * CFA 
	BRA2= (OS(L,M)*FC(M) - OC(L,M)*FS(M)) *  PLM * CFA * M  
!
        If(it   == 1 .and. (sind(teta_ob) /= 0.d0))              THEN 
            u_teta =  u_teta   +  conv_l(l)*BRA1 
	    u_long =  u_long   +  conv_l(l)*BRA2/SIND(TETA_OB)
	ENDIF
!
        If(dit == 1 .and. (sind(teta_ob) /= 0.d0))              THEN 
           du_teta = du_teta + dconv_l(l)*BRA1
	   du_long = du_long + dconv_l(l)*BRA2/SIND(TETA_OB)
	ENDIF
!
                                                         ENDIF 				 	
!
        ENDDO
!	         
        ENDDO         
!
!
!
!
! # Special case: Observer in one of the poles.
!
 log_1 = teta_ob == 0._dp   .OR. teta_ob ==  180._dp
!
 if(log_1) then                                    
   u_teta = 0._dp ;   u_long = 0._dp 
  du_teta = 0._dp ;  du_long = 0._dp   
          endif 
!
!
!+========================================+
! # Assignement to vectors CORR & D_CORR  |
!+========================================+
!
   corr(1) =  u_rad  ;   corr(2) =  u_teta
   corr(3) =  u_long ;   corr(4) =   geoid 
!
 d_corr(1) =  du_rad ; d_corr(2) = du_teta
 d_corr(3) = du_long ; d_corr(4) = d_geoid
!
!
!
 end subroutine OCEAN_CORRECTION  
!
!
!
!
!
!
!
! +======================+
   SUBROUTINE READ_OCEAN !
! +======================+  
! 
  use COMMON; implicit NONE   
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
! # Memorizza i coefficienti della funzione OCEANO su una base * REALE * |
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
!                              
        integer (i4b) l,m,ind_l,ind_m 
!
        open(38, file = 'oceano.128', status = 'unknown')  	                
!
	do 1 ind_l = 0, lmax
        do 2 ind_m = 0, ind_l 
!   
           read  (38, '(2(1x,i3),2(1x,e14.8))') l,m, OC(L,M),OS(L,M) 
!
  2     Continue 
  1     Continue  
!
        close(38) 
!
END SUBROUTINE READ_OCEAN 
!
!
!
!
!
!
!
!*************************************
!
  subroutine SPEC (CODE, ILM, LT)
!
!*************************************
!
use COMMON_FOR_SPECTRUM; implicit NONE 
!
!
! # La sbr SPEC helps to select the Earth model.
!
!   Input:   NV = Number of viscoelastic layers (from TASK_*.DAT)
!
!   Output:  NROOTS = Number of relaxation modes expected
!
!            R   (0:NV+1) = Radii of the interfaces
!            RHO (0:NV+1) = Density of the layers
!            RMU (0:NV+1) = Rigidity of the layers
!
!
!          where:
!
!            R(0)      = radius of the CMB
!            R(NV+1)   = radius of the Earth surface
!            RHO(0)    = density of the core
!            RHO(NV+1) = density of the lithosphere
!            RMU(0)    = rigidity of the core = 0.
!            RMU(NV+1) = rigidity of the lithosphere
!
!
!
!
!
!
! Features of the PREM model:
!
! o) We use the PREM 200s values
! o) IN PREM, the core radius is 3480 km
! o)  "   " , the earth radius is 6371 km
! o)  "   " , the "670" discontinuity is at radius R=5701 km
! o)  "   " , at the "670" discontinuity, d(rho)/rho is 12%'
! o)  "   " , there is a discontinuity at the depth 24.4 km (MOHO)
! o)  "   " , the "420" discontinuity is at radius R=5971 (it is in fact a '400')
! o)  "   " , there is a discontinuity at the depth of = 220 km (R= 6151 km)
! 
!
!
! 
!
!
real     (qp) :: LT      ! Lithospheric thickness from TASK_*.DAT (km)
real     (qp) :: rcmb    ! Radius of the CMB (fixed) 
real     (qp) :: rade    ! Radius of the Earth (fixed) 
real     (qp) :: JUNK    ! A junk 
real     (qp) :: DEL     ! Increment for the Lm 
real     (qp) :: TLM     ! Thickness of the Lower MAntle
!
parameter (rcmb = 3480._qp, rade = 6371._qp)  ! Fixed parameters 
!
integer (i4b) :: k       ! Do-loop index 
integer (i4b) :: CODE    ! A Code to label the various models, given NV
integer (i4b) :: NLM     ! Number of Lm layers 
integer (i4b) :: ILM     ! (0/1) Controls the Lm profile for NV= 7 and NV=9
INTEGER (i4b) :: idens   ! density inversions switch (0/1)
!
!
!
!
If(NV/=1.and.nv/=2.and.nv/=3.and.nv/=4.and.nv/=5.and.nv/=7.and.nv/=8.and.nv/=9.and.iv==0)Then
Write(99,*)'ERROR IN SBR. SPEC: The required model is not available '
WRITE(99,*)'in the models library. **** JOB ABORTED *************** ';STOP
Endif
If(NV/=1.and.nv/=2.and.nv/=3.and.nv/=4.and.nv/=5.and.nv/=7.and.nv/=8.and.nv/=9.and.iv==1)Then
Write (*,*)'ERROR IN SBR. SPEC: The required model is not available '
WRITE (*,*)'in the models library. **** JOB ABORTED *************** ';STOP
Endif
!
!
!
!
!
! #--------------------#
! #       NV = 1       #
! #--------------------#
!
  IF (NV == 1) THEN 
!
!                                  
! NV=1  CODE=0 ---> Averaged PREM, 30 <= LT <= 300 km.
! NV=1  CODE=1 ---> Yuen Sabadini Boschi [1982], LT= 100 km.
! NV=1  CODE=2 ---> Ondrej Cadek, LT= 70 km.
! NV=1  CODE=3 ---> Ondrej Cadek, LT= 200 km.
! NV=1  CODE=4 ---> Giunchi & Spada [2000] BUT self-g., 30 <= LT <= 300 km
! 
!
If ((CODE<0 .or. CODE>4) .and. iv==0) Then
Write(99,*) 'ERROR in SPEC: The model CODE is not available '
Write(99,*) '**** JOB ABORTED ******************************';STOP
Endif
If ((CODE<0 .or. CODE>4) .and. iv==1) Then
Write(*,*)  'ERROR in SPEC: The model CODE is not available '
Write(*,*)  '**** JOB ABORTED ******************************';STOP
Endif
!
!
!
!
!
if(CODE == 0) then    !Averaged PREM, 30 <= LT <= 300 km.                                        
!
!
nroots=4*nv
!
If((LT < 30q0 .OR. LT > 300q0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC:  The LT parameter is out of bounds'
Write(99,*) 'It  must  be  30 <= LT <= 300 km  for NV=1 and CODE=0'
Write(99,*) '******** JOB ABORTED ********************************';STOP
Endif
If((LT < 30q0 .OR. LT > 300q0) .and. iv==1) Then
Write (*,*) 'ERROR in Sbr SPEC:  The LT parameter is out of bounds'
Write (*,*) 'It  must  be  30 <= LT <= 300 km  for NV=1 and CODE=0'
Write (*,*) '******** JOB ABORTED ********************************';STOP
Endif  
!
r (0)  = rcmb; r (1) = rade-LT; r (2) = rade 
!
 call prem (r (0), 0.q0,  JUNK,      rho(0)  )    ! Core     
 call prem (r (2), r (1), rmu (2),   rho(2)  )    ! Litho
 call prem (r (1), r (0), rmu (1),   rho(1)  )    ! Mantle
 rmu(0) = 0q0 
!
endif    						    
!
!+++++++++++++++++++++++
if(CODE == 1) then   !Yuen Sabadini Boschi [1982], LT= 100 km. 				            
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10927._qp; rmu(0)=0._qp       
r (1)      = 6271._qp; rho(1)=4314._qp ; rmu(1)=1.45_qp     
r (2)      = 6371._qp; rho(2)=2689._qp ; rmu(2)=0.28_qp    
!
endif                                                       
!+++++++++++++++++++++++
!
!+++++++++++++++++++++++
if(CODE == 2) then   !Ondrej Cadek, LT= 70 km.  				            
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10927._qp; rmu(0)=0._qp       
r (1)      = 6301._qp; rho(1)=4500._qp ; rmu(1)=1.45_qp     
r (2)      = 6371._qp; rho(2)=4500._qp ; rmu(2)=1.45_qp    
!
endif     
!+++++++++++++++++++++++
!
!+++++++++++++++++++++++
if(CODE == 3) then  !Ondrej Cadek, LT= 200 km.  				            
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10927._qp; rmu(0)=0._qp       
r (1)      = 6171._qp; rho(1)=4500._qp ; rmu(1)=1.45_qp     
r (2)      = 6371._qp; rho(2)=4500._qp ; rmu(2)=1.45_qp    
!
endif   
!+++++++++++++++++++++++
!
!
!
if(CODE == 4) then  !Giunchi & Spada [2000] BUT self-g., 30 <= LT <= 300 km
nroots=4*nv
!
If((LT < 30q0 .OR. LT > 300q0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds '
Write(99,*) 'It  must  be  30 <= LT <= 300 km  for NV=1 and CODE=4'
Write(99,*) '**** JOB ABORTED ************************************';STOP
Endif
If((LT < 30q0 .OR. LT > 300q0) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds '
Write(*,*)  'It  must  be  30 <= LT <= 300 km  for NV=1 and CODE=4'
Write(*,*)  '**** JOB ABORTED ************************************';STOP
Endif  
!
r (0)      = 3480._qp     ;  rho(0)=10977._qp;  rmu(0)=0._qp
r (1)      = 6371._qp - LT;  rho(1)=4518._qp ;  rmu(1)=1.45_qp
r (2)      = 6371._qp     ;  rho(2)=3300._qp ;  rmu(2)=0.28_qp
!
endif
!
ENDIF    ! Endif on NV=1
!
!
!
! #--------------------#
! #       NV = 2       #
! #--------------------#
!
  IF (NV==2) THEN 
!
!
! NV=2  CODE=0 ---> Averaged PREM, 30 <= LT <= 300 km.
! NV=2  CODE=1 ---> Averaged PREM, BUT d(rho)/rho at 670 from PREM, 30 <= LT <= 300 km.
! NV=2  CODE=2 ---> Yuen Sabadini Boschi [1982], LT= 100 km.
! NV=2  CODE=3 ---> Bills & James [1997], LT= 100 km.
! NV=2  CODE=4 ---> Lefftz Sabadini Legros  [1994], LT= 150 km.
! NV=2  CODE=5 ---> Piersanti Postseismic, LT= 80 km.
! NV=2  CODE=6 ---> Ricard Sabadini Spada [1992], Model h, LT= 100 km.
!
!
If((CODE<0 .or. CODE>6) .and. iv==0) Then
Write(99,*)'ERROR in Sbr SPEC: The model CODE is not available '
Write(99,*)'**** JOB ABORTED **********************************';STOP
Endif
If((CODE<0 .or. CODE>6) .and. iv==1) Then
Write(*,*) 'ERROR in Sbr SPEC: The model CODE is not available '
Write(*,*) '**** JOB ABORTED **********************************';STOP
Endif
!
!
!
If(CODE == 0) then  !Averaged PREM, 30 <= LT <= 300 km.
!
          Write(99,*) 'Averaged PREM, 30 <= LT <= 300 km'
IF(iv==1) Write (*,*) 'Averaged PREM, 30 <= LT <= 300 km'
!
!
nroots=4*nv
!
If((LT < 30q0 .OR. LT > 300q0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  30 <= LT <= 300 km  for NV=2  and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
Endif
If((LT < 30q0 .OR. LT > 300q0) .and. iv==1) Then
Write(*,*) 'ERROR in Sbr SPEC: The LT parameter is out of  bounds'
Write(*,*) 'It must be  30 <= LT <= 300 km  for NV=2  and CODE=0'
Write(*,*) '**** JOB ABORTED ************************************';STOP
Endif  
!
r(3)    = rade            
r(2)    = rade-LT         
r(1)    = r(3)-670._qp   
r(0)    = rcmb           
!
call prem(r(0), 0.q0, rmu(0), rho(0) )  ! Core 
call prem(r(1), r(0), rmu(1), rho(1) )  ! Mant. Inf.
call prem(r(2), r(1), rmu(2), rho(2) )  ! Mant. Sup. 
call prem(r(3), r(2), rmu(3), rho(3) )  ! Litho
rmu (0) = 0._qp  
!
              endif
!
!
If(CODE == 1) then  !Averaged PREM, BUT d(rho)/rho at 670 from PREM, 30 <= LT <= 300 km.
!
          write(99,*) &
'Averaged PREM, BUT d(rho)/rho at 670 is from PREM, 30<=LT<=300 km'
IF(iv==1) WRITE(*,*)  &
'Averaged PREM, BUT d(rho)/rho at 670 is from PREM, 30<=LT<=300 km'
!
nroots=4*nv
!
If((LT < 30q0 .OR. LT > 300q0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  30 <= LT <= 300 km  for NV=2  and CODE=1'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
Endif
If((LT < 30q0 .OR. LT > 300q0) .and. iv==1) Then
Write(*,*) 'ERROR in Sbr SPEC:  The LT parameter is out of bounds'
Write(*,*) 'It must be  30 <= LT <= 300 km  for NV=2  and CODE=1'
Write(*,*) '**** JOB ABORTED ************************************';STOP
Endif  
!
!
r(3)    = rade           
r(2)    = rade-LT         
r(1)    = r(3)-670._qp    
r(0)    = rcmb           
!
call prem(r(0), 0.q0, rmu(0), rho(0) )  ! Core 
call prem(r(1), r(0), rmu(1), rho(1) )  ! Mant. Inf.
call prem(r(2), r(1), rmu(2), rho(2) )  ! Mant. Sup. 
call prem(r(3), r(2), rmu(3), rho(3) )  ! Lito 
rmu (0) = 0._qp  
!
! 
rho(1) = (4.41241q0+4.38071q0)/2.q0*1.q3
rho(2) = (3.99214q0+3.98399q0)/2.q0*1.q3
!
             Endif 
!
!
!
If(CODE == 2) then !Yuen Sabadini Boschi [1982], LT= 100 km.
!
                 Write(99,*) 'Yuen Sabadini Boschi [1982], LT= 100 km'
IF(IV==1)        Write(* ,*) 'Yuen Sabadini Boschi [1982], LT= 100 km'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10927._qp; rmu(0)=0._qp       ! Core 
r (1)      = 5701._qp; rho(1)=4919._qp ; rmu(1)=2.17_qp     ! Lower mantle
r (2)      = 6271._qp; rho(2)=4430._qp ; rmu(2)=0.837_qp    ! Upper mantle
r (3)      = 6371._qp; rho(3)=2689._qp ; rmu(3)=0.282_qp    ! Litosphere  
!
Endif 
!
!
!
if(CODE == 3) then  !Bills & James [1997], LT= 100 km.
!
           Write(99,*) 'Bills & James [1997], LT= 100 km'
IF(iv==1)  Write(* ,*) 'Bills & James [1997], LT= 100 km'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10925._qp; rmu(0)=0._qp       ! Core 
r (1)      = 5701._qp; rho(1)=4508._qp ; rmu(1)=1.99_qp     ! Lower mantle
r (2)      = 6271._qp; rho(2)=4120._qp ; rmu(2)=0.954_qp    ! Upper mantle
r (3)      = 6371._qp; rho(3)=2771._qp ; rmu(3)=0.315_qp    ! Litosphere  
!
endif 
!
!
!
if(CODE == 4) then  !Lefftz Sabadini Legros  [1994], LT= 150 km.
!
!
           write(99,*) 'Lefftz Sabadini Legros  [1994], LT= 150 km'
IF(iv==1)  write(* ,*) 'Lefftz Sabadini Legros  [1994], LT= 150 km'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10987._qp; rmu(0)=0._qp        ! Core 
r (1)      = 5701._qp; rho(1)=4904._qp ; rmu(1)=2.225_qp     ! Lower mantle
r (2)      = 6221._qp; rho(2)=3666._qp ; rmu(2)=0.9169_qp    ! Upper mantle
r (3)      = 6371._qp; rho(3)=3232._qp ; rmu(3)=0.6114_qp    ! Litosphere  
!
endif 
!
!
If(CODE == 5) then  !Piersanti Postseismic, LT= 80 km.
!
          Write(99,*)'Postseismic Piersanti model, LT= 80 km '
IF(iv==1) WRITE (*,*)'Postseismic Piersanti model, LT= 80 km '

nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10932._qp; rmu(0)=0._qp        ! Core 
r (1)      = 5701._qp; rho(1)=4878._qp ; rmu(1)=2.171_qp     ! Lower mantle
r (2)      = 6291._qp; rho(2)=3614._qp ; rmu(2)=0.8464_qp    ! Upper mantle
r (3)      = 6371._qp; rho(3)=3115._qp ; rmu(3)=0.5597_qp    ! Litosphere  
!
endif 
!
if(CODE == 6) then !Ricard Sabadini Spada [1992], Model h, LT= 100 km.
!
          WRITE(99,*)  'Ricard Sabadini Spada [1992], Model h, LT= 100 km'
IF(iv==1) WRITE(*, *)  'Ricard Sabadini Spada [1992], Model h, LT= 100 km'
!
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10926._qp; rmu(0)= 0._qp      ! Core 
r (1)      = 5701._qp; rho(1)=4508._qp ; rmu(1)= 1.51_qp    ! Lower mantle
r (2)      = 6271._qp; rho(2)=4120._qp ; rmu(2)= 1.38_qp    ! Upper mantle
r (3)      = 6371._qp; rho(3)=4120._qp ; rmu(3)= 1.38_qp    ! Litosphere  
!
endif 
!
!
ENDIF   ! Endif on NV=2
!
!
!
!
! #--------------------#
! #       NV = 3       #
! #--------------------#
!
!
  IF (NV == 3) THEN 
!
!
! NV=3  CODE=0 ---> Averaged PREM, 30 <= LT <= 300 km.
! NV=3  CODE=1 ---> Not completely Averaged Prem, 30 <= LT <= 300 km.
! NV=3  CODE=2 ---> Cianetti Giunchi Spada [2002], LT=  120 km.
! NV=3  CODE=3 ---> James & Morgan [1990] GRL (table #1), LT= 120 Km.
! NV=3  CODE=4 ---> Similar to W. R. Peltier [1985], LT=  120 km.
! NV=3  CODE=5 ---> Paul Johnston [benchmark, 1997], LT= 70 km.
!
! NEW -- Feb 06, 2017
! NV=3  CODE=6 ---> Vermeersen & Sabadini LT= 120 km.
! NV=3  CODE=7 ---> Model "M3-L70-V01" reference model for the "Test suite" - 2009
!
!
If( (CODE<0 .or. CODE>7) .and. iv==0) Then
Write(99,*)'ERROR in Sbr SPEC: The CODE is not available'
Write(99,*)'**** JOB ABORTED ***************************';stop
Endif
If( (CODE<0 .or. CODE>7) .and. iv==1) Then
Write(*,*) 'ERROR in Sbr SPEC: The CODE is not available'
Write(*,*) '**** JOB ABORTED ***************************';Stop
Endif
!
!
!
!
if(CODE == 0) then  ! Averaged PREM, 30 <= LT <= 300 km.
!
           Write(99,*) 'Averaged PREM, 30 <= LT <= 300 km'
IF(iv==1)  Write(* ,*) 'Averaged PREM, 30 <= LT <= 300 km'
!
nroots=4*nv
!
If((LT < 30q0 .OR. LT > 300q0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  30 <= LT <= 300 km   for NV=3 and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
Endif
If((LT < 30q0 .OR. LT > 300q0) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  30 <= LT <= 300 km   for NV=3 and CODE=0'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
Endif  
!
r(4)    = rade
r(3)    = rade-LT
r(2)    = rade-400._qp   ! Discontinuity at 400 km depth
r(1)    = rade-670._qp   ! Discontinuity at 670 km depth
r(0)    = rcmb           ! CMB radius
!
call prem(r(0), 0.q0, rmu(0), rho(0))  ! Core 
call prem(r(1), r(0), rmu(1), rho(1))  ! Lower mantle
call prem(r(2), r(1), rmu(2), rho(2))  ! Transition Zone
call prem(r(3), r(2), rmu(3), rho(3))  ! Upper mantle
call prem(r(4), r(3), rmu(4), rho(4))  ! Litho
rmu (0) = 0._qp  
!
endif 
!
!
!
If(CODE == 1) then  ! Not completely Averaged Prem, 30 <= LT <= 300 km.
!
Write(99,*) '3-layers mantle model NOT PREM averaged          '
Write(99,*) 'Lower mantle density = PREM value at depth =670- '
Write(99,*) 'Upper mantle density = PREM value at depth =400+ '
Write(99,*) 'Transition zone density = average PREM value     '
Write(99,*) 'The shear modulus is PREM-averaged               '
!
IF(iv==1) THEN
Write(*,*)  '3-layers mantle model NOT PREM averaged          '
Write(*,*)  'Lower mantle density = PREM value at depth =670- '
Write(*,*)  'Upper mantle density = PREM value at depth =400+ '
Write(*,*)  'Transition zone density = average PREM value     '
Write(*,*)  'The shear modulus is PREM-averaged               '
endif
!
nroots=4*nv
!
r(4)    = rade
r(3)    = rade-LT
r(2)    = rade-400._qp   ! Discontinuity at 400 km depth
r(1)    = rade-670._qp   ! Discontinuity at 670 km depth
r(0)    = rcmb           ! CMB radius
!
call prem(r(0), 0.q0, rmu(0), rho(0))  ! Core 
call prem(r(1), r(0), rmu(1), rho(1))  ! LM 
call prem(r(2), r(1), rmu(2), rho(2))  ! TZ 
call prem(r(3), r(2), rmu(3), rho(3))  ! UM 
call prem(r(4), r(3), rmu(4), rho(4))  ! Litho
rmu (0) = 0._qp  
!
!
rho(2) = (4.41241q0+4.38071q0)/2.q0*1.q3
rho(3) = (3.54325q0+3.51639q0)/2.q0*1.q3  
!
endif 
!
!                      
!
If(CODE == 2) then ! Cianetti Giunchi Spada [2002], LT=  120 km.
!
          Write(99,*) '3-- layer mantle model by CIANETTI et al. [2002]'
IF(iv==1) Write (*,*) '3-- layer mantle model by CIANETTI et al. [2002]'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10925._qp; rmu(0)=0._qp        ! Core 
r (1)      = 5701._qp; rho(1)=4508._qp ; rmu(1)=2.000_qp     ! Lower mantle
r (2)      = 5951._qp; rho(2)=4220._qp ; rmu(2)=1.100_qp     ! TZ  
r (3)      = 6251._qp; rho(3)=4120._qp ; rmu(3)=0.950_qp     ! Upper mantle
r (4)      = 6371._qp; rho(4)=4120._qp ; rmu(4)=0.730_qp     ! Litosphere  
!
endif 
!
!
!
if(CODE == 3) then ! James & Morgan [1990] GRL (table #1), LT= 120 Km.
!
          write(99,*) '3-- layer mantle model by James and Morgan (1990) (Tab. #1)  '
IF(iv==1) WRITE(*, *) '3-- layer mantle model by James and Morgan (1990) (Tab. #1)  '
nroots=4*nv
!
r (0)      = 3486._qp; rho(0)=11110._qp; rmu(0)=0._qp        ! Core 
r (1)      = 5701._qp; rho(1)=4900._qp ; rmu(1)=2.300_qp     ! Lower mantle
r (2)      = 5951._qp; rho(2)=3800._qp ; rmu(2)=1.450_qp     ! TZ  
r (3)      = 6251._qp; rho(3)=3550._qp ; rmu(3)=0.710_qp     ! Upper mantle
r (4)      = 6371._qp; rho(4)=2900._qp ; rmu(4)=0.400_qp     ! Litosphere  
!
endif 
!
!
!
if(CODE == 4) then ! Similar to W. R. Peltier [1985], LT=  120 km.
!
write(99,*) 'Similar to the 3-- layer mantle model by Peltier (1985) (Tab. #2)'
IF(iv==1) &
write (*,*) 'Similar to the 3-- layer mantle model by Peltier (1985) (Tab. #2)'
!
!
!
!
nroots=4*nv
!
r(4)    = rade           ! Radius of the Earth  
r(0)    = rcmb           ! Radius of the CMB 
rmu(0)  = 0.q0           
!
 call prem(r(0), 0.q0, rmu(0), rho(0))                        ! Core 
r (1) = rade-670._qp; rho(1)=4372._qp ; rmu(1)=rho(1)*6117._qp**2/1.q11 ! Lm
r (2) = rade-420._qp; rho(2)=4100._qp ; rmu(2)=rho(2)*5475._qp**2/1.q11 ! Tz
r (3) = rade-120._qp; rho(3)=3959._qp ; rmu(3)=rho(3)*5219._qp**2/1.q11 ! Um
 call prem(r(4), r(3), rmu(4), rho(4))                        ! Litho
!
endif 
!
!
!
if(CODE == 5) then   ! Paul Johnston [benchmark, 1997], LT= 70 km.
!
write(99,*) '3-- layer mantle model by Paul Johnston (benchmark, 1997)'
IF(iv==1) &
write( *,*) '3-- layer mantle model by Paul Johnston (benchmark, 1997)'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10750._qp; rmu(0)=0._qp          ! Core 
r (1)      = 5701._qp; rho(1)=4978._qp ; rmu(1)=2.2834_qp      ! Lower mantle
r (2)      = 5951._qp; rho(2)=3871._qp ; rmu(2)=1.0549_qp      ! TZ  
r (3)      = 6301._qp; rho(3)=3438._qp ; rmu(3)=0.70363_qp     ! Upper mantle
r (4)      = 6371._qp; rho(4)=3037._qp ; rmu(4)=0.50605_qp     ! Litosphere
!
endif 
!
!
if(CODE == 6) then   ! 
!
write(99,*) '3-- layer mantle model by Vermeersen and Sabadini 1996'
IF(iv==1) &
write( *,*) '3-- layer mantle model by Vermeersen and Sabadini 1996'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10932._qp; rmu(0)=0._qp        ! Core 
r (1)      = 5701._qp; rho(1)=4878._qp ; rmu(1)=2.19_qp      ! Lower mantle
r (2)      = 5951._qp; rho(2)=3857._qp ; rmu(2)=1.06_qp      ! TZ  
r (3)      = 6250._qp; rho(3)=3434._qp ; rmu(3)=0.727_qp     ! Upper mantle
r (4)      = 6371._qp; rho(4)=3184._qp ; rmu(4)=0.602_qp     ! Litosphere
!
endif 
!
!
if(CODE == 7) then   ! Test suite model M3-L70-V01 
!
write(99,*) '3-- layer mantle model (Test suite model M3-L70-V01)'
IF(iv==1) &
write( *,*) '3-- layer mantle model (Test suite model M3-L70-V01)'
!
nroots=4*nv
!
r (0)      = 3480._qp; rho(0)=10750._qp; rmu(0)=0._qp          ! Core 
r (1)      = 5701._qp; rho(1)=4978._qp ; rmu(1)=2.28340_qp     ! Lower mantle
r (2)      = 5951._qp; rho(2)=3871._qp ; rmu(2)=1.05490_qp     ! TZ  
r (3)      = 6301._qp; rho(3)=3438._qp ; rmu(3)=0.70363_qp     ! Upper mantle
r (4)      = 6371._qp; rho(4)=3037._qp ; rmu(4)=0.50605_qp     ! Litosphere
!
endif 
!
!
ENDIF  !  Endif on NV = 3.
!
!
!
!
!
! #--------------------#
! #       NV = 4       #
! #--------------------#
!
!
! 
  IF (NV == 4) THEN 
!
!                                          
! NV=4  CODE=0 ---> Averaged PREM, 40 <= LT <= 150 km.
! NV=4  CODE=1 ---> Averaged PREM, LM not evaraged in Rho, 40 <= LT <= 150 km.
! NV=4  CODE=2 ---> Averaged, LM not averaged in Rho nor Mu, 40 <= LT <= 150 km.     
!
!
if((CODE <0 .or. CODE >2) .and. iv==0) Then
Write(99,*)'ERROR in SPEC: The model CODE is not available'
Write(99,*)'**** JOB ABORTED *****************************';STOP
Endif
if((CODE <0 .or. CODE >2) .and. iv==1) Then
Write (*,*)'ERROR in SPEC: The model CODE is not available'
Write (*,*)'**** JOB ABORTED *****************************';STOP
Endif 
!
!
!
If(CODE == 0) then 
!
!
Write(99,*)'4--layers PREM- averaged mantle model with:  '
Write(99,*)'Elastic lithosphere 40 <= LT <= 150 km       '
Write(99,*)'Shallow upper mantle 1  (70 <= Thick <=180)  '
Write(99,*)'Shallow upper mantle 2  (Thick =  180 km)    '
Write(99,*)'Transition zone         (Thick =  270 km)    '
Write(99,*)'Lower mantle down to the CMB                 '
!
IF(iv==1) THEN
Write(*,*) '4--layers PREM- averaged mantle model with:  '
Write(*,*) 'Elastic lithosphere 40 <= LT <= 150 km       '
Write(*,*) 'Shallow upper mantle 1  (70 <= Thick <=180)  '
Write(*,*) 'Shallow upper mantle 2  (Thick =  180 km)    '
Write(*,*) 'Transition zone         (Thick =  270 km)    '
Write(*,*) 'Lower mantle down to the CMB                 '
ENDIF
!
nroots=4*nv
!
If((LT < 40q0 .OR. LT > 150q0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=150 km   for NV=4 and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40q0 .OR. LT > 150q0) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=150 km   for NV=4 and CODE=0'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
r(5)    = rade
r(4)    = rade-LT
r(3)    = rade-220._qp   ! 200 km depth discontinuity
r(2)    = rade-400._qp   ! 400       "        "
r(1)    = rade-670._qp   ! 670       "        "
r(0)    = rcmb           ! CMB radius
!
!
 call prem(r(0), 0.q0, rmu(0), rho(0))  ! Core 
 call prem(r(1), r(0), rmu(1), rho(1))  ! Lower mantle
 call prem(r(2), r(1), rmu(2), rho(2))  ! Transition zone
 call prem(r(3), r(2), rmu(3), rho(3))  ! Shallow upper mantle (1)
 call prem(r(4), r(3), rmu(4), rho(4))  ! Shallow upper mantle (2)
 call prem(r(5), r(4), rmu(5), rho(5))  ! Lithosphere
rmu (0) = 0._qp  
!
endif 
!
!
If(CODE == 1) Then
!
!
Write(99,*) '4-- layers mantle PREM- averaged mantle model:   '
Write(99,*) 'Elastic lithosphere with 40 <= LT <= 150 km      '
Write(99,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)   '
Write(99,*) 'Shallow upper mantle 2  (Thick   =  180 km)      '
Write(99,*) 'Transition zone         (Thick   =  270 km)      '
Write(99,*) 'Lower mantle density = that of PREM at 670-      '
Write(99,*) 'Lower mantle rigidity is PREM- averaged dal PREM '
!
IF(iv==1) THEN
Write(*,*) '4-- layers mantle PREM- averaged mantle model:   '
Write(*,*) 'Elastic lithosphere with 40 <= LT <= 150 km      '
Write(*,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)   '
Write(*,*) 'Shallow upper mantle 2  (Thick   =  180 km)      '
Write(*,*) 'Transition zone         (Thick   =  270 km)      '
WRITE(*,*) 'Lower mantle density = that of PREM at 670-      '
Write(*,*) 'Lower mantle rigidity is PREM- averaged dal PREM '
ENDIF
!
nroots=4*nv
!
!
If((LT < 40q0 .OR. LT > 150q0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=150 km   for NV=4 and CODE=1'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40q0 .OR. LT > 150q0) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=150 km   for NV=4 and CODE=1'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
!
r(5)    = rade
r(4)    = rade-LT
r(3)    = rade-220._qp   !  Discontinuita' at 220 km depth.
r(2)    = rade-400._qp   !        "           400 km "   "
r(1)    = rade-670._qp   !        "           670 km "   "
r(0)    = rcmb           !  Radius of the CMB 
!
 call prem(r(0), 0.q0,  rmu(0), rho(0))  ! Core 
 call prem(r(1), r(0),  rmu(1), rho(1))  ! LM 
 call prem(r(2), r(1),  rmu(2), rho(2))  ! TZ 
 call prem(r(3), r(2),  rmu(3), rho(3))  ! Sm(1)
 call prem(r(4), r(3),  rmu(4), rho(4))  ! Sm(2) 
 call prem(r(5), r(4),  rmu(5), rho(5))  ! Lito  
!
rmu (0) = 0._qp  
rho (1) = (4.41241+4.38071)/2.q0*1.q3   ! Lower mantle density
!
endif 
!
!
if(CODE == 2) Then
!
Write(99,*) '4-- layers mantle PREM- averaged mantle model:   '
Write(99,*) 'Elastic lithosphere with 40 <= LT <= 150 km      '
Write(99,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)   '
Write(99,*) 'Shallow upper mantle 2  (Thick   =  180 km)      '
Write(99,*) 'Transition zone         (Thick   =  270 km)      '
Write(99,*) 'Lower mantle density = that of PREM at 670-      '
Write(99,*) 'Lower mantle rigidity = that of PREM at 670-     '
!
IF(iv==1) then
Write (*,*) '4-- layers mantle PREM- averaged mantle model:   '
Write (*,*) 'Elastic lithosphere with 40 <= LT <= 150 km      '
Write (*,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)   '
Write (*,*) 'Shallow upper mantle 2  (Thick   =  180 km)      '
Write (*,*) 'Transition zone         (Thick   =  270 km)      '
Write (*,*) 'Lower mantle density = that of PREM at 670-      '
Write (*,*) 'Lower mantle rigidity = that of PREM at 670-     '
ENDIF
!
!
nroots=4*nv
!
!
If((LT < 40q0 .OR. LT > 150q0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=150 km   for NV=4 and CODE=2'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40q0 .OR. LT > 150q0) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=150 km   for NV=4 and CODE=2'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
r(5)    = rade
r(4)    = rade-LT
r(3)    = rade-220._qp   !  Discontinuita' at 220 km depth.
r(2)    = rade-400._qp   !        "           400 km "   "
r(1)    = rade-670._qp   !        "           670 km "   "
r(0)    = rcmb           !  Radius of the CMB 
!
call prem(r(0), 0.q0,  rmu(0), rho(0))  ! Core 
call prem(r(1), r(0),  rmu(1), rho(1))  ! LM 
call prem(r(2), r(1),  rmu(2), rho(2))  ! TZ 
call prem(r(3), r(2),  rmu(3), rho(3))  ! Sm(1)
call prem(r(4), r(3),  rmu(4), rho(4))  ! Sm(2) 
call prem(r(5), r(4),  rmu(5), rho(5))  ! Litho
!
rmu (0) = 0._qp  
rmu (1) = (1.584q0 + 1.639q0)/2._qp   ! Lower mantle density
rho (1) = (4.41241+4.38071)/2.q0*1.q3 ! Lower mantle rigidity
!
endif 
!
!
ENDIF  ! Endif on NV=4
!
!
!
!
!
! #--------------------#
! #       NV = 5       #
! #--------------------#
!
!
! 
  IF (NV == 5) THEN 
!
!                                          
! NV=5  CODE=0 ---> VM5a, 40 <= LT <= 100 km
!
!
if((CODE <0 .or. CODE >0) .and. iv==0) Then
Write(99,*)'ERROR in SPEC: The model CODE is not available'
Write(99,*)'**** JOB ABORTED *****************************';STOP
Endif
if((CODE <0 .or. CODE >0) .and. iv==1) Then
Write (*,*)'ERROR in SPEC: The model CODE is not available'
Write (*,*)'**** JOB ABORTED *****************************';STOP
Endif 
!
!
!
If(CODE == 0) then 
!
!
Write(99,*)'VM5a viscosity model (PREM-averaged)              '
Write(99,*)'Lithosphere                  (40 <= LT <= 100 km) '
Write(99,*)'Lower lithosphere               (Thick =   40 km) '
Write(99,*)'Upper mantle                 (280 <= Thick <=340) '
Write(99,*)'Transition zone                 (Thick =  250 km) '
Write(99,*)'Lower mantle 1                  (Thick =  590 km) '
Write(99,*)'Lower mantle 2 down to the CMB                    '
!
IF(iv==1) THEN
Write(*,*) 'VM5a viscosity model (PREM-averaged)                            '
Write(*,*) 'Lithosphere                  (40 <= LT <= 100 km) '
Write(*,*) 'Lower lithosphere               (Thick =   40 km) '
Write(*,*) 'Upper mantle                 (280 <= Thick <=340) '
Write(*,*) 'Transition zone                 (Thick =  250 km) '
Write(*,*) 'Lower mantle 1                  (Thick =  590 km) '
Write(*,*) 'Lower mantle 2 down to the CMB                    '
ENDIF
!
nroots=4*nv
!
If((LT < 40q0 .OR. LT > 100q0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=100 km   for NV=5 and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40q0 .OR. LT > 100q0) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=100 km   for NV=5 and CODE=0'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
r(6)    = rade
r(5)    = rade-LT
r(4)    = rade-LT-40._qp  
r(3)    = rade-420._qp   ! 420 km depth discontinuity
r(2)    = rade-670._qp   ! 670       "        "
r(1)    = rade-1260._qp  ! 1260      "        "
r(0)    = rcmb           ! CMB radius
!
!
 call prem(r(0), 0.q0, rmu(0), rho(0))  ! Core 
 call prem(r(1), r(0), rmu(1), rho(1))  ! Lower mantle
 call prem(r(2), r(1), rmu(2), rho(2))  ! Lower mantle
 call prem(r(3), r(2), rmu(3), rho(3))  ! Transition zone
 call prem(r(4), r(3), rmu(4), rho(4))  ! Upper mantle
 call prem(r(5), r(4), rmu(5), rho(5))  ! Lower lithosphere
 call prem(r(6), r(5), rmu(6), rho(6))  ! Lithosphere
rmu (0) = 0._qp  
!
endif 
!
!
!
ENDIF  ! Endif on NV=5
!
!
!
!
!
! #--------------------#
! #       NV = 7       #
! #--------------------#
!
!
!
  IF (NV == 7) THEN 
! 
!                                      
! NV=7  CODE=0 ---> Averaged PREM, 40 <= LT <= 150 km. 3 SM layers, 4 LM layers.
!
!
if((CODE <0 .or. CODE >0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The CODE is not available'
Write(99,*) '**** JOB ABORTED ***************************';STOP
Endif 
if((CODE <0 .or. CODE >0) .and. iv==1) Then
Write (*,*) 'ERROR in Sbr SPEC: The CODE is not available'
Write (*,*) '**** JOB ABORTED ***************************';STOP
Endif 
!
!
!
If(CODE == 0) then 
!
Write(99,*) '7-layers PREM averaged mantle model:  '
Write(99,*) 'Elastic lithosphere with thickness 40 <= LT <= 150 km  '
Write(99,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)         '
Write(99,*) 'Shallow upper mantle 2  (Thick   =  180 km)            '
Write(99,*) 'Transition zone         (Thick   =  270 km)            '
WRITE(99,*) 'Lower mantle is PREM averaged on 4 layers              '
!
IF(iv==1) THEN
Write(*,*) '7-layers PREM averaged mantle model:  '
Write(*,*) 'Elastic lithosphere with thickness 40 <= LT <= 150 km  '
Write(*,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)         '
Write(*,*) 'Shallow upper mantle 2  (Thick   =  180 km)            '
Write(*,*) 'Transition zone         (Thick   =  270 km)            '
WRITE(*,*) 'Lower mantle is PREM averaged on 4 layers              '
ENDIF
!
nroots=4*nv
!
!
If((LT < 40q0 .OR. LT > 150q0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=150 km   for NV=7 and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40q0 .OR. LT > 150q0) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=150 km   for NV=7 and CODE=0'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
r(8)    = rade                   ! Radius of the Earth 
r(7)    = rade-LT                ! Litho with thickness 40 <= LT <= 150 km   
r(6)    = rade-220._qp           ! Discontinuity at 220 km depth
r(5)    = rade-400._qp           !       "          400 km "   " 
r(4)    = rade-670._qp           !       "          670 km "   " 
r(0)    = rcmb 
!
!
NLM = 4            !  Number of v.e. layers in the Lower Mantle  (4) 
TLM = r(4)-r(0)    !  Thickness of the Lower Mantle 
!
If (ILM == 1) DEL = TLM/Float(NLM)  ! For ILM=1, the lower mantle layers
                                    ! have all the same thickness
!
If (ILM == 0) DEL = 370.q0          ! For ILM=0, the lower mantle layers
                                    ! have a thickness increasing with depth
!
Do K = 1, NLM-1 
     R(k) = R(k-1) + DEL + & 
     Float(NLM-k)*(TLM - NLM*DEL)/float(NLM)/float(NLM-1)/0.5q0 
Enddo      
!
!
Write(99,*) 'Thickness of the lower mantle layers from BOTTOM to TOP (km):'
Write(99,'(f10.4)')   r(1)-r(0)
Write(99,'(f10.4)')   r(2)-r(1)
Write(99,'(f10.4)')   r(3)-r(2)
Write(99,'(f10.4)')   r(4)-r(3)
IF(iv==1) THEN
Write(*,*)  'Thickness of the lower mantle layers from BOTTOM to TOP (km):'
Write(*,'(f10.4)')   r(1)-r(0)
Write(*,'(f10.4)')   r(2)-r(1)
Write(*,'(f10.4)')   r(3)-r(2)
Write(*,'(f10.4)')   r(4)-r(3)
ENDIF
!
!
call prem(r(0), 0.q0, rmu(0), rho(0))  ! Core 
call prem(r(1), r(0), rmu(1), rho(1))  ! LM(1) 
call prem(r(2), r(1), rmu(2), rho(2))  ! LM(2) 
call prem(r(3), r(2), rmu(3), rho(3))  ! LM(3) 
call prem(r(4), r(3), rmu(4), rho(4))  ! LM(4) 
call prem(r(5), r(4), rmu(5), rho(5))  ! TZ 
call prem(r(6), r(5), rmu(6), rho(6))  ! Sm(1)
call prem(r(7), r(6), rmu(7), rho(7))  ! Sm(2) 
call prem(r(8), r(7), rmu(8), rho(8))  ! Lito  
rmu (0) = 0._qp  
!
endif 
!
ENDIF  ! Endif on NV=7
!
!
!
!
!
! #--------------------#
! #       NV = 8       #
! #--------------------#
!
!
!
  IF (NV == 8) THEN 
! 
!                                      
! NV=7  CODE=0 ---> VM7
!
!
if((CODE <0 .or. CODE >0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The CODE is not available'
Write(99,*) '**** JOB ABORTED ***************************';STOP
Endif 
if((CODE <0 .or. CODE >0) .and. iv==1) Then
Write (*,*) 'ERROR in Sbr SPEC: The CODE is not available'
Write (*,*) '**** JOB ABORTED ***************************';STOP
Endif 
!
!
!
If(CODE == 0) then 
!
Write(99,*) '8-layers PREM averaged VM7 mantle model:  '
Write(99,*) 'Elastic lithosphere with thickness LT = 75km '
Write(99,*) 'Shallow upper mantle 1  (Thick   =   40 km)            '
Write(99,*) 'Shallow upper mantle 2  (Thick   =  305 km)            '
Write(99,*) 'Transition zone         (Thick   =  250 km)            '
Write(99,*) 'Lower mantle 1          (Thick   =  310 km)            '
Write(99,*) 'Lower mantle 2          (Thick   =  490 km)            '
Write(99,*) 'Lower mantle 3          (Thick   =  480 km)            '
Write(99,*) 'Lower mantle 4          (Thick   =  450 km)            '
Write(99,*) 'Lower mantle 5 down to the CMB                         '
!
IF(iv==1) THEN
Write(*,*)  '8-layers PREM averaged VM7 mantle model:  '
Write(*,*)  'Elastic lithosphere with thickness LT = 75km '
Write(*,*)  'Shallow upper mantle 1  (Thick   =   40 km)            '
Write(*,*)  'Shallow upper mantle 2  (Thick   =  305 km)            '
Write(*,*)  'Transition zone         (Thick   =  250 km)            '
Write(*,*)  'Lower mantle 1          (Thick   =  310 km)            '
Write(*,*)  'Lower mantle 2          (Thick   =  490 km)            '
Write(*,*)  'Lower mantle 3          (Thick   =  480 km)            '
Write(*,*)  'Lower mantle 4          (Thick   =  450 km)            '
Write(*,*)  'Lower mantle 5 down to the CMB                         '
ENDIF
!
nroots=4*nv
!
!
r(9)    = rade                   ! Radius of the Earth 
r(8)    = rade-75._qp            ! Litho with thickness LT = 75km  
r(7)    = rade-115._qp           ! Discontinuity at 115 km depth
r(6)    = rade-420._qp           !       "          420 km "   " 
r(5)    = rade-670._qp           !       "          670 km "   " 
r(4)    = rade-980._qp           !       "          980 km "   " 
r(3)    = rade-1470._qp          !       "         1470 km "   " 
r(2)    = rade-1950._qp          !       "         1950 km "   " 
r(1)    = rade-2400._qp          !       "         2400 km "   " 
r(0)    = rcmb 
!
!
!
!
call prem(r(0), 0.q0, rmu(0), rho(0))  ! Core 
call prem(r(1), r(0), rmu(1), rho(1))  ! LM(1) 
call prem(r(2), r(1), rmu(2), rho(2))  ! LM(2) 
call prem(r(3), r(2), rmu(3), rho(3))  ! LM(3) 
call prem(r(4), r(3), rmu(4), rho(4))  ! LM(4) 
call prem(r(5), r(4), rmu(5), rho(5))  ! LM(5)
call prem(r(6), r(5), rmu(6), rho(6))  ! TZ 
call prem(r(7), r(6), rmu(7), rho(7))  ! UM(1)
call prem(r(8), r(7), rmu(8), rho(8))  ! UM(2) 
call prem(r(9), r(8), rmu(9), rho(9))  ! Lito  
rmu (0) = 0._qp  
!
endif 
!
ENDIF  ! Endif on NV=8

!
!
!
!
!
!
! #--------------------#
! #       NV = 9       #
! #--------------------#
!
!
  IF (NV == 9) THEN 
! 
!                                      
! CODE =0 --> Prem mediato (sm=3, lm=6)                     
!
!
if((CODE <0 .or. CODE >0) .and. iv==0) Then
Write(99,*)'ERROR in Sbr SPEC: The  CODE is not available'
Write(99,*)'**** JOB ABORTED ****************************';STOP
Endif
if((CODE <0 .or. CODE >0) .and. iv==1) Then
Write(*,*) 'ERROR in Sbr SPEC: The  CODE is not available'
Write(*,*) '**** JOB ABORTED ****************************';STOP
Endif 
!
!
!
!
!
If(CODE == 0) then 
!
Write(99,*) '9-layers PREM averaged mantle model:  '
Write(99,*) 'Elastic lithosphere with thickness 40 <= LT <= 150 km  '
Write(99,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)         '
Write(99,*) 'Shallow upper mantle 2  (Thick   =  180 km)            '
Write(99,*) 'Transition zone         (Thick   =  270 km)            '
WRITE(99,*) 'Lower mantle is PREM averaged on 4 layers              '
!
IF(iv==1) THEN
Write(*,*) '9-layers PREM averaged mantle model:  '
Write(*,*) 'Elastic lithosphere with thickness 40 <= LT <= 150 km  '
Write(*,*) 'Shallow upper mantle 1  (70 <= Thick <=180 km)         '
Write(*,*) 'Shallow upper mantle 2  (Thick   =  180 km)            '
Write(*,*) 'Transition zone         (Thick   =  270 km)            '
WRITE(*,*) 'Lower mantle is PREM averaged on 6 layers              '
ENDIF
!
nroots=4*nv
!
!
If((LT < 40q0 .OR. LT > 150q0) .and. iv==0) Then
Write(99,*) 'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(99,*) 'It must be  40  <= LT <=150 km   for NV=9 and CODE=0'
Write(99,*) '**** JOB ABORTED ***********************************';STOP
endif
If((LT < 40q0 .OR. LT > 150q0) .and. iv==1) Then
Write(*,*)  'ERROR in Sbr SPEC: The LT parameter is out of bounds'
Write(*,*)  'It must be  40  <= LT <=150 km   for NV=9 and CODE=0'
Write(*,*)  '**** JOB ABORTED ***********************************';STOP
endif
!
!
!
!
r(10)   = rade                   ! Radius of the Earth   
r(9)    = rade-LT                ! Litho with thickness 40 <= LT <= 150 km     
r(8)    = rade-220._qp           ! Discontinuity at 220 km depth
r(7)    = rade-400._qp           !       "          400 km "   " 
r(6)    = rade-670._qp           !       "          670 km "   " 
r(0)    = rcmb                   ! Radius of the CMB 
!
!
!
NLM = 6            !  Number of v.e. layers in the Lower Mantle  (6) 
TLM = r(6)-r(0)    !  Thickness of the Lower Mantle 
!
If (ILM == 1) DEL = TLM/Float(NLM)  ! For ILM=1, the lower mantle layers
                                    ! have all the same thickness
!
If (ILM == 0) DEL = 300.q0          ! For ILM=0, the lower mantle layers
                                    ! have a thickness increasing with depth
!
!
Do K = 1, NLM-1 
     R(k) = R(k-1) + DEL + & 
     Float(NLM-k)*(TLM - NLM*DEL)/float(NLM)/float(NLM-1)/0.5q0 
Enddo      
!
!

Write(99,*) 'Lower mantle layers thickness, from BOTTOM to TOP (km):'
Write(99,'(f10.4)')   r(1)-r(0)
Write(99,'(f10.4)')   r(2)-r(1)
Write(99,'(f10.4)')   r(3)-r(2)
Write(99,'(f10.4)')   r(4)-r(3)
Write(99,'(f10.4)')   r(5)-r(4)
Write(99,'(f10.4)')   r(6)-r(5)
!
IF(iv==1) THEN
Write(99,*) 'Lower mantle layers thickness, from BOTTOM to TOP (km):'
Write(99,'(f10.4)')   r(1)-r(0)
Write(99,'(f10.4)')   r(2)-r(1)
Write(99,'(f10.4)')   r(3)-r(2)
Write(99,'(f10.4)')   r(4)-r(3)
Write(99,'(f10.4)')   r(5)-r(4)
Write(99,'(f10.4)')   r(6)-r(5)
endif

!
!
call prem(r(0),  0.q0, rmu(0),  rho(0))   ! Core 
call prem(r(1),  r(0), rmu(1),  rho(1))   ! LM 
call prem(r(2),  r(1), rmu(2),  rho(2))   ! LM 
call prem(r(3),  r(2), rmu(3),  rho(3))   ! LM 
call prem(r(4),  r(3), rmu(4),  rho(4))   ! LM 
call prem(r(5),  r(4), rmu(5),  rho(5))   ! LM 
call prem(r(6),  r(5), rmu(6),  rho(6))   ! LM 
call prem(r(7),  r(6), rmu(7),  rho(7))   ! TZ 
call prem(r(8),  r(7), rmu(8),  rho(8))   ! Sm(1)
call prem(r(9),  r(8), rmu(9),  rho(9))   ! Sm(2) 
call prem(r(10), r(9), rmu(10), rho(10))  ! Lito  
rmu (0) = 0._qp  
endif 
!
!
ENDIF  !  Endif on  NV=9
!
!
!
!
!
!
!
!
!
!
!
!
!
!  +-----------------------------------------+
!  |  Conversion of R, RMU & RHO in SI units |
!  +-----------------------------------------+
!
    r  (:) =r  (:)*1.q3        ! R in m
    rmu(:) =rmu(:)*1.q11       ! RMU in Pa 
    rho(:) =rho(:)             ! Rho in Kg/m**3
!
!
!
Write(99,*) '---------------------------------------------------------------'
Write(99,*) ' Radii, densities, shear moduli & viscosity from bottom to top '
Write(99,*) '        SI units: m, kg/m**3, Pa, Pa.s, respectively           '
Write(99,*) '---------------------------------------------------------------'
!
Do k = 0, nv+1   
  Write(99,'(i4, 1x, 5(2x, E14.8))') k, r(k), rho(k), rmu(k), vis(k) 
Enddo  
!
!
!
if(lt .ne. (r(nv+1)-r(nv))/1q3 .and. iv==0) then
!write(99,*) lt, (r(nv+1)-r(nv))/1q3  
write(99,*) ' + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  '
write(99,*) 'WARNING from sbr SPEC:   The LT value given in input is not ' 
write(99,*) 'consistent with the default for the current values of NV &  '
write(99,*) 'and CODE.  The litho thickness is set to the default value. '
write(99,*) ' + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  '
endif 
if(lt .ne. (r(nv+1)-r(nv))/1q3 .and. iv==1) then 
!write(*,*) lt, (r(nv+1)-r(nv))/1q3  
write(* ,*) ' + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  '
write(* ,*) 'WARNING from sbr SPEC:   The LT value given in input is not ' 
write(* ,*) 'consistent with the default for the current values of NV &  '
write(* ,*) 'and CODE.  The litho thickness is set to the default value. '
write(* ,*) ' + + + + + + + + + + + + + + + + + + + + + + + + + + + + +  '
endif 
!
          Write(99,'(a51,F12.4)') &
'Lithospheric thickness effectively employed (km) =', (r(nv+1)-r(nv))/1q3 
if(iv==1) Write(* ,'(a51,F12.4)') &
'Lithospheric thickness effectively employed (km) =', (r(nv+1)-r(nv))/1q3 
!
!
!
!  +-------------------------------+
!  | Test for density inversions   |
!  +-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=+
!
WRITE(99,*) 'Looking for density inversions...'
!
idens=0
Do k = 1, nv+1  
if (rho (k) - rho (k - 1)  >  0._qp)  then
idens=1
Write ( 99, * ) 'WARNING from Sbr. SPEC: a density inversion exists !'
endif
Enddo
!
IF(idens==0)             WRITE(99,*) 'No density inversions found'
IF(idens==0 .and. iv==1) WRITE(* ,*) 'No density inversions found'
!
!
!
!  +--------------------------------+
!  | Earth mass and average density |
!  +-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-+
!
IF((NV==3).AND.(CODE==7)) THEN
!
    emass = 4._qp * pi * rho (0) * r (0) **3 / 3._qp  
    do k = 1, nv+1  
        emass = emass + (4._qp / 3._qp) * pi * (r (k) **3 - r (k - 1) **3) * rho (k)
    enddo  
    rhoea = emass/(4._qp/3._qp)/pi/r(nv+1)/r(nv+1)/r(nv+1)
!
    Write(99,*) 'Using true values for the Earth mass and average density'
    Write(99,*) 'Effective mass of the model (kg)'
    Write(99,'(d14.8)'), emass
!
    Write(99,*) 'Average density of the model (kg/m^3)'
    Write(99,'(d14.8)'), rhoea
!
ELSE
    emass = ref_emass
    rhoea = ref_rhoea
!
    Write(99,*) 'Using reference values for the Earth mass and average density'
    Write(99,*) 'Mass of the model (kg)'
    Write(99,'(d14.8)'), emass
!
    Write(99,*) 'Average density of the model (kg/m^3)'
    Write(99,'(d14.8)'), rhoea
!
ENDIF
!
!
!
!
!
!
END SUBROUTINE SPEC 
!
!
!
!
!
!
  subroutine m_3 (ne)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! # A collection of Warning & Error messages for Task#3
! - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   USE for_task3
   implicit NONE
!
   INTEGER ne
!
!
   IF(ne==30010) THEN
       IF(IV/=0.and.IV/=1) THEN
       WRITE(99,*) 'ERROR in sbr task_3.dat:     The VERBOSE switch '
       WRITE(99,*) 'can only assume values 0 and 1. ** JOB ABORTED* '; STOP
       endif
   ENDIF
!
!
  IF(ne==30012) THEN
   IF  (n_armo >=2 .and. iv==0) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw Harmonic_Degrees must appear once'
   WRITE(99,*) '= JOB ABORTED ====================== JOB ABORTED ============';STOP
   ENDIF
   IF  (n_armo >=2 .and. iv==1) THEN
   WRITE(* ,*) 'ERROR in sbr TASK_3: The kw Harmonic_Degrees must appear once'
   WRITE(* ,*) '= JOB ABORTED ====================== JOB ABORTED ============'
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw Harmonic_Degrees must appear once'
   WRITE(99,*) '= JOB ABORTED ====================== JOB ABORTED ============';STOP
   ENDIF
  ENDIF
!
!
IF(ne==30020) THEN
 IF(IV==0) then
  If(lmin<llmin .or. lmax>llmax .or. lmax<lmin .or. lmin<0 .or. lmax<0) Then
  Write(99,*)'ERROR IN SBR. task_1.dat: One of the following forbidden conditions'
  Write(99,*)'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
  Write(99,*)'Check the input file task_3.dat ****** JOB ABORTED ****************';Stop
  Endif
 Endif
IF(IV==1) then
 If(lmin<llmin .or. lmax>llmax .or. lmax<lmin .or. lmin<0 .or. lmax<0) Then
  Write(99,*)'ERROR IN SBR. task_1.dat: One of the following forbidden conditions'
  Write(99,*)'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
  Write(99,*)'Check the input file task_3.dat ****** JOB ABORTED ****************'
  Write(*,*) 'ERROR IN SBR. task_1.dat: One of the following forbidden conditions'
  Write(*,*) 'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
  Write(*,*) 'Check the input file task_3.dat ****** JOB ABORTED ****************';Stop
 Endif
 Endif
ENDIF
!
!
IF(ne==30030) THEN
    IF(only_elastic==0.and.iv==1) &
           WRITE(*,*)  'Elastic AND viscous fields'
    IF(only_elastic==1.and.iv==1) &
           WRITE(*,*)  'ONLY Elastic fields'
    IF(only_elastic==0) WRITE(99,*) 'Elastic AND viscous fields'
    IF(only_elastic==1) WRITE(99,*) 'ONLY Elastic fields'
ENDIF
!
!
IF(ne==30040) THEN
    IF(DROP_MODES==1.and.iv==1) &
           WRITE(*,*)  'The ''non physical''modes will be killed'
    IF(DROP_MODES==0.and.iv==1) &
           WRITE(*,*)  'The ''non physical''modes will NOT be killed' 
!
 IF(DROP_MODES==1) WRITE(99,*)'The ''non physical''modes will be killed'
 IF(DROP_MODES==0) WRITE(99,*)'The ''non physical''modes will NOT be killed'
ENDIF
!
!
IF(ne==30050) THEN
    IF(iarmo==0.and.IV==0) then
                 WRITE(99,*) &
                 'ERROR IN SBR. TASK_3:      The KW Harmonic_Degrees '
                 WRITE(99,*) &
                 'must be activated before Make_Model can be executed';STOP
    endif
    IF(iarmo==0.and.IV==1) then
                 WRITE(99,*) &
                 'ERROR IN SBR. TASK_3:      The KW Harmonic_Degrees '
                 WRITE(99,*) &
                 'must be activated before Make_Model can be executed'
                 WRITE(*,*) &
                 'ERROR IN SBR. TASK_3:      The KW Harmonic_Degrees '
                 WRITE(*,*) &
                 'must be activated before Make_Model can be executed';STOP
    endif
ENDIF
!
!
  IF(ne==30052) THEN
   IF  (n_make >=2 .and. iv==0) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw Make_Model must appear once'
   WRITE(99,*) '====== JOB ABORTED ====================================';STOP
   ENDIF
   IF  (n_make >=2 .and. iv==1) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw Make_Model must appear once'
   WRITE(99,*) '====== JOB ABORTED ===================================='
   WRITE(* ,*) 'ERROR in sbr TASK_3: The kw Make_Model must appear once'
   WRITE(* ,*) '====== JOB ABORTED ====================================';STOP
   ENDIF
  ENDIF
!
!
IF(ne==30060) THEN
 IF(imake==1.and.IV==0) then
        WRITE(99,*)'ERROR IN SBR. TASK_3: The Kws Make_Model and External_Model'
        WRITE(99,*)'are mutually exclusive. If Make_Model is active,  External_'
        WRITE(99,*)'model must be not. **** JOB ABORTED ***********************';STOP
 endif
 IF(imake==1.and.IV==1) then
        WRITE(99,*)'ERROR IN SBR. TASK_3: The Kws Make_Model and External_Model'
        WRITE(99,*)'are mutually exclusive. If Make_Model is active,  External_'
        WRITE(99,*)'model must be not. **** JOB ABORTED ***********************'
        WRITE(*,*) 'ERROR IN SBR. TASK_3: The Kws Make_Model and External_Model'
        WRITE(*,*) 'are mutually exclusive. If Make_Model is active,  External_'
        WRITE(*,*) 'model must be not. **** JOB ABORTED ***********************';STOP
 endif
ENDIF
!
!
  IF(ne==30062) THEN
   IF  (n_exte >=2 .and. iv==0) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw External_Model must appear once'
   WRITE(99,*) '====== JOB ABORTED ======== ====== JOB ABORTED ============';STOP
   ENDIF
   IF  (n_exte >=2 .and. iv==1) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw External_Model must appear once'
   WRITE(99,*) '====== JOB ABORTED ======== ====== JOB ABORTED ============'
   WRITE(* ,*) 'ERROR in sbr TASK_3: The kw External_Model must appear once'
   WRITE(* ,*) '====== JOB ABORTED ======== ====== JOB ABORTED ============';STOP
   ENDIF
  ENDIF
!
!
IF(ne==30070) THEN
   Write(99,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '
   Write(99,*) ' The model is the one by   GIUNCHI & SPADA,  [2000] '
   Write(99,*) ' NV=1 and LT=100 km, see  GRL table (1) for details '
   Write(99,*) ' The model is Gravitating but not self -Gravitating '
   Write(99,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '
   IF(IV==1) THEN
   Write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '
   Write(*,*) ' The model is the one by   GIUNCHI & SPADA,  [2000] '
   Write(*,*) ' NV=1 and LT=100 km, see  GRL table (1) for details '
   Write(*,*) ' The model is Gravitating but not self -Gravitating '
   Write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '
   ENDIF
Endif  
!
!
IF(ne==30075) THEN
IF(iv==0 .AND. (imake==0 .AND. iexte==0)) THEN
WRITE(99,*) 'ERROR in sbr. TASK_3: One of the kws Make_Model and '
WRITE(99,*) 'External_Model must be active  before the kw Ad_Hoc '
WRITE(99,*) 'can be activated ======== JOB ABORTED ============= ';STOP
ENDIF
IF(iv==1 .AND. (imake==0 .AND. iexte==0)) THEN
WRITE(* ,*) 'ERROR in sbr. TASK_3: One of the kws Make_Model and '
WRITE(* ,*) 'External_Model must be active  before the kw Ad_Hoc '
WRITE(* ,*) 'can be activated ======== JOB ABORTED ============= '
WRITE(99,*) 'ERROR in sbr. TASK_3: One of the kws Make_Model and '
WRITE(99,*) 'External_Model must be active  before the kw Ad_Hoc '
WRITE(99,*) 'can be activated ======== JOB ABORTED ============= ';STOP
ENDIF
ENDIF
!
!
  IF(ne==30077) THEN
   IF  (n_adhoc >=2 .and. iv==0) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw Ad_Hoc must appear once'
   WRITE(99,*) '====== JOB ABORTED =========== JOB ABORTED ========';STOP
   ENDIF
   IF  (n_adhoc >=2 .and. iv==1) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw Ad_Hoc must appear once'
   WRITE(99,*) '====== JOB ABORTED =========== JOB ABORTED ========'
   WRITE(* ,*) 'ERROR in sbr TASK_3: The kw Ad_Hoc must appear once'
   WRITE(* ,*) '====== JOB ABORTED =========== JOB ABORTED ========';STOP
   ENDIF
  ENDIF
!
!
IF(ne==30080) THEN
If(Numero_di_elementi <= 0 .and. iv==0) Then
Write(99,*) 'ERROR in TASK_3: The number of elements of the aggregate  '
Write(99,*) 'declared in file  ',      filename,     '   must be >= 1 '
Write(99,*) '**** JOB ABORTED **************************************** '; Stop
Endif
If(Numero_di_elementi <= 0 .and. iv==1) Then
Write(* ,*) 'ERROR in TASK_3: The number of elements of the aggregate  '
Write(* ,*) 'declared in file  ',      filename,     '  must be >= 1 '
Write(* ,*) '**** JOB ABORTED **************************************** '
Write(99,*) 'ERROR in TASK_3: The number of elements of the aggregate  '
Write(99,*) 'declared  in file ',      filename,     '  must be >= 1 '
Write(99,*) '**** JOB ABORTED **************************************** '; Stop
Endif
If(Numero_di_elementi >= MEl .and. iv==0) Then
Write(99,*) 'ERROR in TASK_3: The number of elements of the aggregate given  '
Write(99,*) 'in file', filename, 'exceeds ', MEl, ' **** JOB ABORTED ***  ';stop
Endif
If(Numero_di_elementi >= MEl .and. iv==1) Then
Write(99,*) 'ERROR in TASK_3: The number of elements of the aggregate given  '
Write(99,*) 'in file', filename, 'exceeds ', MEl, ' **** JOB ABORTED ***  '
Write(* ,*) 'ERROR in TASK_3: The number of elements of the aggregate given  '
Write(* ,*) 'in file', filename, 'exceeds ', MEl, ' **** JOB ABORTED ***  ';stop
Endif
ENDIF
!
!
IF(ne==30090) THEN
    IF(C_L(j)/=10.and.C_L(j)/=11.and.C_L(j)/=20.and.C_L(j)/=21.and. &
                     C_L(j)/=30.and.C_L(j)/=50.and.iv==0) THEN
    WRITE(99,*) 'ERROR in sbr TASK_3: A wrong load code is given. JOB ABORTED';stop
  endif
    IF(C_L(j)/=10.and.C_L(j)/=11.and.C_L(j)/=20.and.C_L(j)/=21.and. &
                     C_L(j)/=30.and.C_L(j)/=50.and.iv==1) THEN
    WRITE(99,*) 'ERROR in sbr TASK_3: A wrong load code is given. JOB ABORTED'
    WRITE(* ,*) 'ERROR in sbr TASK_3: A wrong load code is given. JOB ABORTED';stop
    endif
ENDIF
!
!
IF(ne==30100) then
 IF ( (IOC(j) < 0 .or. IOC(j) > 1) .and. IV==0) THEN
  WRITE(99,*)'ERROR IN SBR. TASK_2: The option IOC = ', IOC, ' Is not allowed !'
   WRITE(99,*)'Choose IOC = 0 ---> to set NO compensation on real oceans  '
   WRITE(99,*)'Choose IOC = 1 ---> to set    compensation on real oceans  '
  WRITE(99,*)'******* JOB ABORTED *************************************  ';STOP
 ENDIF
 IF ( (IOC(j) < 0 .or. IOC(j) > 1) .and. IV==1) THEN
  WRITE(99,*)'ERROR IN SBR. TASK_2: The option IOC = ', IOC, ' Is not allowed !'
   WRITE(99,*)'Choose IOC = 0 ---> to set NO compensation on real oceans  '
    WRITE(99,*)'Choose IOC = 1 ---> to set    compensation on real oceans  '
     WRITE(99,*)'******* JOB ABORTED *************************************  '
     WRITE(*,*) 'ERROR IN SBR. TASK_2: The option IOC = ', IOC, ' Is not allowed !'
    WRITE(*,*) 'Choose IOC = 0 ---> to set NO compensation on real oceans  '
   WRITE(*,*) 'Choose IOC = 1 ---> to set    compensation on real oceans  '
  WRITE(*,*) '******* JOB ABORTED *************************************  ';STOP
 ENDIF
endif
!
!
IF(ne==30110) then
 IF ((C_L(j)==11.or.C_L(j)==21).and.IOC(j)==1.and.IV==0.and.Numero_di_elementi==1) THEN
 WRITE(99,*)'ERROR  IN SBR. TASK_2: If the compensation on a complementary '
  WRITE(99,*)'disc load is chosen (C_L(j) = 11 or 21),   there cannot '
  WRITE(99,*)'be compensation on a realistic ocean (IOC=1).  Try with IOC=1 '
 WRITE(99,*)'******* JOB ABORTED ***************************************** ';STOP
 ENDIF
!
 IF ((C_L(j)==11.or.C_L(j)==21).and.IOC(j)==1.and.IV==0.and.Numero_di_elementi==1) THEN
 WRITE(99,*)'ERROR  IN SBR. TASK_2: If the compensation on a complementary '
  WRITE(99,*)'disc load is chosen    (C_L = 11 or 21),   there cannot '
   WRITE(99,*)'be compensation on a realistic ocean (IOC=1).  Try with IOC=1 '
    WRITE(99,*)'******* JOB ABORTED ***************************************** '
    WRITE(*,*) 'ERROR  IN SBR. TASK_2: If the compensation on a complementary '
   WRITE(*,*) 'disc load is chosen   (C_L = 11 or 21),   there  cannot '
  WRITE(*,*) 'be compensation on a realistic ocean (IOC=1).  Try with IOC=1 '
 WRITE(*,*) '******* JOB ABORTED ***************************************** ';STOP
 ENDIF
endif
!
!
IF(ne==30120) THEN
If(Numero_di_elementi>=2.AND.iv==0.and.(C_L(j)==21.or.C_L(j)==11)) Then
Write(99,*) 'ERROR in SBR TASK_3:       For aggregates of 2 or more elements '
 Write(99,*) 'mass conservation cannot occour on complementary secondary loads'
Write(99,*) '**** JOB ABORTED ****** JOB ABORTED ****** JOB ABORTED ******** ';STOP
endif
If(Numero_di_elementi>=2.AND.iv==0.and.(C_L(j)==21.or.C_L(j)== 11)) Then
Write(99,*) 'ERROR in SBR TASK_3:       For aggregates of 2 or more elements '
 Write(99,*) 'mass conservation cannot occour on complementary secondary loads'
  Write(99,*) '**** JOB ABORTED ****** JOB ABORTED ****** JOB ABORTED ******** '
  Write(* ,*) 'ERROR in SBR TASK_3:       For aggregates of 2 or more elements '
 Write(* ,*) 'mass conservation cannot occour on complementary secondary loads'
Write(* ,*) '**** JOB ABORTED ****** JOB ABORTED ****** JOB ABORTED ******** ';STOP
endif
endif
!
!
IF(ne==30130) THEN
  If ((l_h(j)<0.or.l_h(j)>8) .and. IV==0) then
  Write(99,*) 'ERROR IN SBR. TASK_3: Only time-histories with label  '
   Write(99,*) '0,  1,  2,  ... 7, and 8 are available at the moment  '
  Write(99,*) ' ***** JOB ABORTED ********************************** ';STOP
  Endif
  If ((l_h(j)<0.or.l_h(j)>8) .and. IV==1) then
  Write(99,*) 'ERROR IN SBR. TASK_3: Only time-histories with label  '
   Write(99,*) '0,  1,  2,  ... 7, and 8 are available at the moment  '
    Write(99,*) ' ***** JOB ABORTED ********************************** '
    Write(*,*)  'ERROR IN SBR. TASK_3: Only time-histories with label  '
   Write(*,*)  '0,  1,  2,  ... 7, and 8 are available at the moment  '
  Write(*,*)  ' ***** JOB ABORTED ********************************** ';STOP
  Endif
endif
!
!
IF(ne==30140) THEN
       IF(Altezz(j) < 0.0 .and. IV==0)  THEN
       WRITE(99,*) 'ERROR IN SBR. TASK_3: The maximum thickness of the element #', J
       WRITE(99,*) 'must be >=0. ****** JOB ABORTED**** ****** JOB ABORTED *****'; STOP
        ENDIF
        IF(Altezz(j) < 0.0 .and. IV==1)  THEN
       WRITE(99,*) 'ERROR IN SBR. TASK_3: The maximum thickness of the element #', J
       WRITE(99,*) 'must be >=0. ****** JOB ABORTED**** ****** JOB ABORTED *****'
       WRITE(* ,*) 'ERROR IN SBR. TASK_3: The maximum thickness of the element #', J
       WRITE(* ,*) 'must be >=0. ****** JOB ABORTED**** ****** JOB ABORTED *****'; STOP
        ENDIF
endif
!
!
IF(ne==30145) THEN
       Write(99,*) 'Longitude of the center (deg)=',  Long_c(j)
       Write(99,*) 'Colatitude of the center (deg)=', Cola_c(j)
       Write(99,*) 'Half- amplitude (deg)=',          Amplit(j)
       Write(99,*) 'Maximum load thickness (m) =',    Altezz(j)
       IF(iv==1) then
       Write(*,*) 'Longitude of the center (deg)=',  Long_c(j)
       Write(*,*) 'Colatitude of the center (deg)=', Cola_c(j)
       Write(*,*) 'Half- amplitude (deg)=',          Amplit(j)
       Write(*,*) 'Maximum load thickness (m) =',    Altezz(j)
       endif
endif
!
!
IF(ne==30150) THEN
IF(l_h(j)==6 .OR. l_h(j)==7 .OR. l_h(j)==8) THEN
       WRITE(99,*) 'WARNING in sbr. TASK_2: The value of the maximum thickness'
       WRITE(99,*) 'given above will be superseded if LH= 6, 7, or 8'
       IF(iv==1) THEN
            WRITE(*,*) 'WARNING in sbr. TASK_2: The value of the maximum thickness'
            WRITE(*,*) 'given above will be superseded if LH= 6, 7, or 8'
       ENDIF
ENDIF
ENDIF
!
!
IF(ne==30160) THEN
        IF (Grande(j) < 0.0 .and. IV==0)  THEN
       WRITE(99,*) 'ERROR IN SBR. TASK_3: The maximum thickness of the element #', J
       WRITE(99,*) 'must be >=0. ****** JOB ABORTED**** ****** JOB ABORTED *****'; STOP
        ENDIF
        IF (Grande(j) < 0.0 .and. IV==1)  THEN
       WRITE(99,*) 'ERROR IN SBR. TASK_3: The maximum thickness of the element #', J
       WRITE(99,*) 'must be >=0. ****** JOB ABORTED**** ****** JOB ABORTED *****'
       WRITE(* ,*) 'ERROR IN SBR. TASK_3: The maximum thickness of the element #', J
       WRITE(* ,*) 'must be >=0. ****** JOB ABORTED**** ****** JOB ABORTED *****'; STOP
        ENDIF
endif
!
!
IF(ne==30165) then
      Write(99,*) 'Longitude of the point mass (deg)=',    Long_c(j)
      Write(99,*) 'Colatitude of the point mass (deg)=',    Cola_c(j)
      Write(99,*) 'Maximum mass of the point mass (kg) =',  Grande(j)
      IF(iv==1) THEN
      Write(*,*) 'Longitude of the point mass (deg)=',    Long_c(j)
      Write(*,*) 'Colatitude of the point mass (deg)=',    Cola_c(j)
      Write(*,*) 'Maximum mass of the point mass (kg) =',  Grande(j)
      ENDIF
ENDIF
!
!
IF(ne==30170) THEN
       WRITE(99,*) 'WARNING in sbr. TASK_2: The value of the maximum mass'
       WRITE(99,*) 'given above will be superseded if LH= 6, 7, or 8'
       WRITE(*, *) 'WARNING in sbr. TASK_2: The value of the maximum mass'
       WRITE(*, *) 'given above will be superseded if LH= 6, 7, or 8'
endif
!
!
IF(ne==30180) THEN
        IF(Altezz(j) < 0.0 .and. IV==0)  THEN
      WRITE(99,*) 'ERROR IN SBR. TASK_3: The maximum thickness of the element #', J
      WRITE(99,*) 'must be >=0. ****** JOB ABORTED**** ****** JOB ABORTED *****'; STOP
        ENDIF
        IF(Altezz(j) < 0.0 .and. IV==1)  THEN
      WRITE(99,*) 'ERROR IN SBR. TASK_3: The maximum thickness of the element #', J
      WRITE(99,*) 'must be >=0. ****** JOB ABORTED**** ****** JOB ABORTED *****'
      WRITE(* ,*) 'ERROR IN SBR. TASK_3: The maximum thickness of the element #', J
      WRITE(* ,*) 'must be >=0. ****** JOB ABORTED**** ****** JOB ABORTED *****'; STOP
        ENDIF
endif
!
!
IF(ne==30185) THEN
    Write(99,*) 'Longitude of the centroid (deg)=',    Long_c(j)
     Write(99,*) 'Colatitude of the centroid (deg)=',   Cola_c(j)
      Write(99,*) 'Width in longitude (deg)=',           d_long (j)
     Write(99,*) 'Width in  colatitude (deg)=',         d_cola (j)
    Write(99,*) 'Maximum load thickness (m) =',        Altezz(j)
    IF(iv==1) THEN
    Write(*,*) 'Longitude of the centroid (deg)=',    Long_c(j)
     Write(*,*) 'Colatitude of the centroid (deg)=',   Cola_c(j)
      Write(*,*) 'Width in longitude (deg)=',           d_long (j)
     Write(*,*) 'Width in  colatitude (deg)=',         d_cola (j)
    Write(*,*) 'Maximum load thickness (m) =',        Altezz(j)
    ENDIF
endif
!
!
IF(ne==30190) THEN
       WRITE(99,*) 'WARNING in sbr. TASK_2: The value of the maximum thickness'
       WRITE(99,*) 'given above will be superseded if LH= 6, 7, or 8'
       WRITE(*, *) 'WARNING in sbr. TASK_2: The value of the maximum thickness'
       WRITE(*, *) 'given above will be superseded if LH= 6, 7, or 8'
endif
!
!
IF(ne==30200) THEN
                        If (T_au(j) <=0.0 .and. IV==0) Then
   write(99,*) 'ERROR IN SBR. TASK_3: The lapse of time between loading '
   write(99,*) 'and un- loading must be > 0 kyrs **** JOB ABORTED ***** '; STOP
                        Endif
!
                        If (T_au(j) <=0.0 .and. IV==1) Then
   write(99,*) 'ERROR IN SBR. TASK_3: The lapse of time between loading '
   write(99,*) 'and un- loading must be > 0 kyrs **** JOB ABORTED ***** '
   write(* ,*) 'ERROR IN SBR. TASK_3: The lapse of time between loading '
   write(* ,*) 'and un- loading must be > 0 kyrs **** JOB ABORTED ***** '; STOP
                        Endif
!
endif
!
!
IF(ne==30210) THEN
   If (T_au(j)<=0.0 .and. IV==0) Then
   write(99,*) 'ERROR IN SBR. TASK_3:   The lenght of the melting phase '
   write(99,*) 'must be positive in time- history #3. ** JOB ABORTED *** '
  Endif
  If  (T_au(j)<=0.0 .and. IV==1) Then
   write(99,*) 'ERROR IN SBR. TASK_3:   The lenght of the melting phase '
   write(99,*) 'must be positive in time- history #3. ** JOB ABORTED *** '
   write(* ,*) 'ERROR IN SBR. TASK_3:   The lenght of the melting phase '
   write(* ,*) 'must be positive in time- history #3. ** JOB ABORTED *** '; STOP
  Endif
endif
!
!
IF(ne==30220) THEN
if((N_R(j)<0.or.N_R(j)>5) .and. IV==0) then
write(99,*) 'ERROR in sbr TASK_3:    The maximum allowed number of'
write(99,*) 'loading and unloading phases is NR=5. JOB ABORTED ** '; STOP
endif
if((N_R(j)<0.or.N_R(j)>5) .and. IV==1) then
write(99,*) 'ERROR in sbr TASK_3:    The maximum allowed number of'
write(99,*) 'loading and unloading phases is NR=5. JOB ABORTED ** '
write(* ,*) 'ERROR in sbr TASK_3:    The maximum allowed number of'
write(* ,*) 'loading and unloading phases is NR=5. JOB ABORTED ** '; STOP
endif
endif
!
!
IF(ne==30230) THEN
if(T_AUC(j)<=0.0 .and. IV==0) then
write(99,*) 'ERROR in sbr TASK_3: The length of the loading phase'
write(99,*) 'must be positive.  ***************** JOB ABORTED ** '; STOP
endif
if(T_AUC(j)<=0.0 .and. IV==1) then
write(99,*) 'ERROR in sbr TASK_3: The length of the loading phase'
write(99,*) 'must be positive.  ***************** JOB ABORTED ** '
write(* ,*) 'ERROR in sbr TASK_3: The length of the loading phase'
write(* ,*) 'must be positive.  ***************** JOB ABORTED ** '; STOP
endif
endif
!
!
IF(ne==30240) THEN
if(D_INC(j)<=0.0 .and. IV==0) then
write(99,*) 'ERROR in sbr TASK_3: The length of the un -loading phase'
write(99,*) 'must be positive.  *****************     JOB ABORTED ** '; STOP
endif
if(D_INC(j)<=0.0 .and. IV==1) then
write(99,*) 'ERROR in sbr TASK_3: The length of the un- loading phase'
write(99,*) 'must be positive.  *****************     JOB ABORTED ** '
write(* ,*) 'ERROR in sbr TASK_3: The length of the un -loading phase'
write(* ,*) 'must be positive.  *****************     JOB ABORTED ** '; STOP
endif
endif
!
!
IF(ne==30250) THEN
if( T_per(j) <= 0.0 .and. IV==0) then
write(99,*) 'ERROR in sbr TASK_3:    The period of the sinusoidal'
write(99,*) 'time-history must be positive. ***** JOB ABORTED ** '; STOP
endif
if( T_per(j) <= 0.0 .and. IV==0) then
write(99,*) 'ERROR in sbr TASK_3:    The period of the sinusoidal'
write(99,*) 'time-history must be positive. ***** JOB ABORTED ** '
write(* ,*) 'ERROR in sbr TASK_3:    The period of the sinusoidal'
write(* ,*) 'time-history must be positive. ***** JOB ABORTED ** '; STOP
endif
ENDIF
!
!
  IF(ne==30260) THEN
  IF(iv==0 .AND. n_ss(j) == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  At least two points must be given in '
  WRITE(99,*) '                    order to specify the time-history th #6 '
  WRITE(99,*) '************************* JOB ABORTED ********************* ';STOP
         ENDIF
  IF(iv==1 .AND. n_ss(j) == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  At least two points must be given in '
  WRITE(99,*) '                    order to specify the time-history th #6 '
  WRITE(99,*) '************************* JOB ABORTED ********************* '
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  At least two points must be given in '
  WRITE(* ,*) '                    order to specify the time-history th #6 '
  WRITE(* ,*) '************************* JOB ABORTED ********************* ';STOP
         ENDIF
  ENDIF
!
!
  IF(ne==30270) THEN
  IF(iv==0 .AND. n_ss(j) > n_ss_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The number of linear segments in time-'
  WRITE(99,*) 'history #6 must not exceed ', n_ss_max
  WRITE(99,*) '************************* JOB ABORTED ********************** ';STOP
         ENDIF
  IF(iv==1 .AND. n_ss(j) > n_ss_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The number of linear segments in time-'
  WRITE(99,*) 'history #6 must not exceed ', n_ss_max
  WRITE(99,*) '************************* JOB ABORTED ********************** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  The number of linear segments in time-'
  WRITE(* ,*) 'history #6 must not exceed ', nss_max
  WRITE(* ,*) '************************* JOB ABORTED ********************** ';STOP
         ENDIF
  ENDIF
!
!
 IF(ne==30280) THEN
         IF(iv==0 .AND. Tmp(j,0) /= 0.0) THEN
     WRITE(99,*) 'ERROR IN SBR. TASK_3:    The time t_0 must be =zero'
     WRITE(99,*) 'in time-history #6. ***** JOB ABORTED ************'
     endif
         IF(iv==1 .AND. Tmp(j,0) /= 0.0) THEN
     WRITE(99,*) 'ERROR IN SBR. TASK_3:    The time t_0 must be =zero'
     WRITE(99,*) 'in time-history #6. ***** JOB ABORTED ************'
     WRITE(* ,*) 'ERROR IN SBR. TASK_3:    The time t_0 must be =zero'
     WRITE(* ,*) 'in time-history #6. ***** JOB ABORTED ************';STOP
         ENDIF
 ENDIF
!
!
 IF(ne==30290) THEN
     do i=0, n_ss(j)
         IF(iv==0 .AND. Tmp(j,i) < 0.0) THEN
     WRITE(99,*) 'ERROR IN SBR. TASK_3:    The times which mark the boundaries'
     WRITE(99,*) 'between adjacent linear phases must not be negative in time-'
     WRITE(99,*) 'history #6. ***** JOB ABORTED ******************************';stop
     endif
         IF(iv==1 .AND. Tmp(j,i) < 0.0) THEN
     WRITE(99,*) 'ERROR IN SBR. TASK_3:    The times which mark the boundaries'
     WRITE(99,*) 'between adjacent linear phases must not be negative in time-'
     WRITE(99,*) 'history #6. ***** JOB ABORTED ******************************'
     WRITE(* ,*) 'ERROR IN SBR. TASK_3:    The times which mark the boundaries'
     WRITE(* ,*) 'between adjacent linear phases must not be negative in time-'
     WRITE(* ,*) 'history #6. ***** JOB ABORTED ******************************';stop
         ENDIF
     enddo
 endif
!
!
 IF(ne==30300) THEN
  do i=0, n_ss(j)
         IF(iv==0 .AND. i>=1 .and. Tmp(j,i) <= Tmp(j,i-1)) THEN
  WRITE(99,*) 'ERROR IN SBR. TASK_3:        The times which mark the boundaries'
  WRITE(99,*) 'between adjacent linear phases must be ordered from the smallest'
  WRITE(99,*) 'to the largest in time-history #6. **** JOB ABORTED ************';stop
     endif
         IF(iv==1 .AND. i>=1 .and. Tmp(j,i) <= Tmp(j,i-1)) THEN
  WRITE(99,*) 'ERROR IN SBR. TASK_3:        The times which mark the boundaries'
  WRITE(99,*) 'between adjacent linear phases must be ordered from the smallest'
  WRITE(99,*) 'to the largest in time-history #6. **** JOB ABORTED ************'
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:        The times which mark the boundaries'
  WRITE(* ,*) 'between adjacent linear phases must be ordered from the smallest'
  WRITE(* ,*) 'to the largest in time-history #6. **** JOB ABORTED ************';stop
         ENDIF
  enddo
 endif
!
!
!
 IF(ne==30320) THEN
! +++++++++++++++
  do i=0, n_ss(j)
! +++++++++++++++
  IF(iv==0 .AND. Alh(j,i) < 0.0 .AND. (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#6. *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==1 .AND. Alh(j,i)< 0.0 .AND.  (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#6. *** JOB ABORTED *************  '
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(* ,*) 'or zero at any time in th#6. *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==0 .AND. Alh(j,i) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#6. **** JOB ABORTED **************  ';STOP
         ENDIF
  IF(iv==1 .AND. Alh(j,i) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#6. **** JOB ABORTED **************  '
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#6. **** JOB ABORTED **************  ';STOP
         ENDIF
! +++++
  enddo
! +++++
ENDIF
!
!
IF(ne==30330) THEN
  IF(JUNK <= 0.0 .AND. iv==0) THEN
   WRITE(99,*) 'ERROR in sbr. TASK_3: The maximum heights (or masses) must be > 0 '
   WRITE(99,*) 'in time- history L_h==6. ******* JOB ABORTED ******************** ';stop
  endif
  IF(JUNK <= 0.0 .AND. iv==1) THEN
   WRITE(* ,*) 'ERROR in sbr. TASK_3: The maximum heights (or masses) must be > 0 '
   WRITE(* ,*) 'in time- history L_h==6. ******* JOB ABORTED ******************** '
   WRITE(99,*) 'ERROR in sbr. TASK_3: The maximum heights (or masses) must be > 0 '
   WRITE(99,*) 'in time- history L_h==6. ******* JOB ABORTED ******************** ';stop
  endif
ENDIF
!
!
IF(ne==30340) THEN
  IF( D_ilt(j) <= 0.0 .and. iv==1) THEN
  WRITE(99,*) 'ERROR in sbr. TASK_3: The length of each step during the piecewise '
  WRITE(99,*) 'constant phase must be > 0 for L_h = 7 ===== JOB ABORTED ========= ';stop
  ENDIF
  IF( D_ilt(j)<= 0.0 .and. iv==1) THEN
  WRITE(99,*) 'ERROR in sbr. TASK_3: The length of each step during the piecewise '
  WRITE(99,*) 'constant phase must be > 0 for L_h = 7 ===== JOB ABORTED ========= '
  WRITE(* ,*) 'ERROR in sbr. TASK_3: The length of each step during the piecewise '
  WRITE(* ,*) 'constant phase must be > 0 for L_h = 7 ===== JOB ABORTED ========= ';STOP
  ENDIF
ENDIF
!
!
IF(ne==30350) THEN
  IF(iv==0 .AND. n_s(j) == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  At least two points must be given in '
  WRITE(99,*) '                    order to specify the time-history th #7 '
  WRITE(99,*) '======== JOB ABORTED ============= JOB ABORTED ============ ';STOP
         ENDIF
  IF(iv==1 .AND. n_s(j) == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  At least two points must be given in '
  WRITE(99,*) '                    order to specify the time-history th #7 '
  WRITE(99,*) '======== JOB ABORTED ============= JOB ABORTED ============ '
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  At least two points must be given in '
  WRITE(* ,*) '                    order to specify the time-history th #7 '
  WRITE(* ,*) '======== JOB ABORTED ============= JOB ABORTED ============ ';STOP
         ENDIF
ENDIF
!
!
IF(ne==30360) THEN
  IF(iv==0 .AND. n_s(j) > n_s_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The number of time-steps in the piecewise'
  WRITE(99,*) 'constant phase of time-history #7 must not exceed ', n_s_max
  WRITE(99,*) '========= JOB ABORTED ============= JOB ABORTED ============ ';STOP
         ENDIF
  IF(iv==1 .AND. n_s(j) > n_s_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The number of time-steps in the piecewise'
  WRITE(99,*) 'constant phase of time-history #7 must not exceed ', n_s_max
  WRITE(99,*) '========= JOB ABORTED ============= JOB ABORTED ============ '
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  The number of time-steps in the piecewise'
  WRITE(* ,*) 'constant phase of time-history #7 must not exceed ', n_s_max
  WRITE(* ,*)  '======== JOB ABORTED ============= JOB ABORTED ============ ';STOP
         ENDIF
ENDIF
!
!
IF(ne==30370) THEN
  IF(iv==0 .AND. Alt(j,0)< 0.0 .AND. (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness before the beginning '
  WRITE(99,*) 'of the piecewise constant phase must be >=0 for time-history #7'
  WRITE(99,*) '======= JOB ABORTED ============= JOB ABORTED =================';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,0)< 0.0 .AND.  (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  The load thickness before the beginning '
  WRITE(* ,*) 'of the piecewise constant phase must be >=0 for time-history #7'
  WRITE(* ,*) '======= JOB ABORTED ============= JOB ABORTED ================='
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness before the beginning '
  WRITE(99,*) 'of the piecewise constant phase must be >=0 for time-history #7'
  WRITE(99,*) '======= JOB ABORTED ============= JOB ABORTED =================';STOP
         ENDIF
  IF(iv==0 .AND. Alt(j,0) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load mass before the beginning '
  WRITE(99,*) 'of the piecewise constant phase must be >=0 for time-history #7'
  WRITE(99,*) '======= JOB ABORTED ============= JOB ABORTED =================';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,0) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load mass before the beginning '
  WRITE(99,*) 'of the piecewise constant phase must be >=0 for time-history #7'
  WRITE(99,*) '======= JOB ABORTED ============= JOB ABORTED ================='
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load mass before the beginning '
  WRITE(99,*) 'of the piecewise constant phase must be >=0 for time-history #7'
  WRITE(99,*) '======= JOB ABORTED ============= JOB ABORTED =================';STOP
         ENDIF
ENDIF
!
!
IF(ne==30380) THEN
! ++++++++++++++
  do i=1, n_s(j)
! ++++++++++++++
!
  IF(iv==0 .AND. Alt(j,i) < 0.0 .AND. (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,i) < 0.0 .AND.  (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. *** JOB ABORTED *************  '
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(* ,*) 'or zero at any time in th#7. *** JOB ABORTED *************  ';STOP
         ENDIF
!
  IF(iv==0 .AND. Alt(j,i) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  ';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,i) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  '
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  ';STOP
         ENDIF
!
  IF(iv==0 .AND. Alt(j,i) < 0.0 .AND. (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,i)< 0.0 .AND.  (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. *** JOB ABORTED *************  '
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(* ,*) 'or zero at any time in th#7. *** JOB ABORTED *************  ';STOP
         ENDIF
!
  IF(iv==0 .AND. Alt(j,i) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  ';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,i) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  '
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  ';STOP
         ENDIF
!+++++
 enddo
!+++++
ENDIF
!
!
IF(ne==30390) then
  IF(JUNK <= 0.0 .AND. iv==0) THEN
   WRITE(99,*) 'ERROR in sbr. TASK_3: The maximum thickness (or mass) must be > 0 '
   WRITE(99,*) 'in time- history L_h==7 ======= JOB ABORTED ==================== ';stop
  endif
  IF(JUNK <= 0.0 .AND. iv==1) THEN
   WRITE(* ,*) 'ERROR in sbr. TASK_3: The maximum thickness (or mass) must be > 0 '
   WRITE(* ,*) 'in time- history L_h==7 ======= JOB ABORTED ==================== '
   WRITE(99,*) 'ERROR in sbr. TASK_3: The maximum thickness (or mass) must be > 0 '
   WRITE(99,*) 'in time- history L_h==7 ======= JOB ABORTED ==================== ';stop
  endif
endif
!
!
IF(ne==30400) then
                        If (T_au(j) <=0.0 .and. IV==0) Then
   write(99,*) 'ERROR IN SBR. TASK_3: The length of the loading phase   '
   write(99,*) 'must be > 0 kyrs in time-history #8 **** JOB ABORTED ***** '; STOP
                        Endif
!
                        If (T_au(j) <=0.0 .and. IV==1) Then
   write(99,*) 'ERROR IN SBR. TASK_3: The length of the loading phase   '
   write(99,*) 'must be > 0 kyrs in time-history #8 **** JOB ABORTED ***** '
   write(* ,*) 'ERROR IN SBR. TASK_3: The length of the loading phase   '
   write(* ,*) 'must be > 0 kyrs in time-history #8 **** JOB ABORTED ***** '; STOP
                        Endif
endif
!
!
IF(ne==30410) then
  IF(iv==0 .AND. n_s(j) == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  At least two points must be given in '
  WRITE(99,*) '                    order to specify the time-history th #8 '
  WRITE(99,*) '======== JOB ABORTED ============= JOB ABORTED ============ ';STOP
         ENDIF
  IF(iv==1 .AND. n_s(j) == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  At least two points must be given in '
  WRITE(99,*) '                    order to specify the time-history th #8 '
  WRITE(99,*) '======== JOB ABORTED ============= JOB ABORTED ============ '
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  At least two points must be given in '
  WRITE(* ,*) '                    order to specify the time-history th #8 '
  WRITE(* ,*) '======== JOB ABORTED ============= JOB ABORTED ============ ';STOP
         ENDIF
endif
!
!
IF(ne==30420) then
  IF(iv==0 .AND. n_s(j) > n_s_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The number of time-steps in the piecewise'
  WRITE(99,*) 'constant phase of time-history #8 must not exceed ', n_s_max
  WRITE(99,*) '========= JOB ABORTED ============= JOB ABORTED ============ ';STOP
         ENDIF
  IF(iv==1 .AND. n_s(j) > n_s_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The number of time-steps in the piecewise'
  WRITE(99,*) 'constant phase of time-history #8 must not exceed ', n_s_max
  WRITE(99,*) '========= JOB ABORTED ============= JOB ABORTED ============ '
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  The number of time-steps in the piecewise'
  WRITE(* ,*) 'constant phase of time-history #8 must not exceed ', n_s_max
  WRITE(* ,*)  '======== JOB ABORTED ============= JOB ABORTED ============ ';STOP
         ENDIF
endif
!
!
IF(ne==30430) THEN
  IF(iv==0 .AND. Alt(j,0)< 0.0 .AND. (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness before the beginning '
  WRITE(99,*) 'of the piecewise constant phase must be >=0 for time-history #8'
  WRITE(99,*) '======= JOB ABORTED ============= JOB ABORTED =================';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,0)< 0.0 .AND.  (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  The load thickness before the beginning '
  WRITE(* ,*) 'of the piecewise constant phase must be >=0 for time-history #8'
  WRITE(* ,*) '======= JOB ABORTED ============= JOB ABORTED ================='
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness before the beginning '
  WRITE(99,*) 'of the piecewise constant phase must be >=0 for time-history #8'
  WRITE(99,*) '======= JOB ABORTED ============= JOB ABORTED =================';STOP
         ENDIF
  IF(iv==0 .AND. Alt(j,0) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load mass before the beginning '
  WRITE(99,*) 'of the piecewise constant phase must be >=0 for time-history #8'
  WRITE(99,*) '======= JOB ABORTED ============= JOB ABORTED =================';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,0) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load mass before the beginning '
  WRITE(99,*) 'of the piecewise constant phase must be >=0 for time-history #8'
  WRITE(99,*) '======= JOB ABORTED ============= JOB ABORTED ================='
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load mass before the beginning '
  WRITE(99,*) 'of the piecewise constant phase must be >=0 for time-history #7'
  WRITE(99,*) '======= JOB ABORTED ============= JOB ABORTED =================';STOP
         ENDIF
ENDIF
!
!
IF(ne==30440) then
! ++++++++++++++
  do i=1, n_s(j)
! ++++++++++++++
!
  IF(iv==0 .AND. Alt(j,i) < 0.0 .AND. (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#8 *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,i) < 0.0 .AND.  (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#8 *** JOB ABORTED *************  '
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(* ,*) 'or zero at any time in th#8 *** JOB ABORTED *************  ';STOP
         ENDIF
!
  IF(iv==0 .AND. Alt(j,i) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#8 **** JOB ABORTED **************  ';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,i) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#8 **** JOB ABORTED **************  '
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#8 **** JOB ABORTED **************  ';STOP
         ENDIF
!
  IF(iv==0 .AND. Alt(j,i) < 0.0 .AND. (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#8 *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,i)< 0.0 .AND.  (C_L(j) <=21 .OR. C_L(j) == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#8 *** JOB ABORTED *************  '
  WRITE(* ,*) 'ERROR IN SBR. TASK_3:  The load thickness must be positive  '
  WRITE(* ,*) 'or zero at any time in th#7. *** JOB ABORTED *************  ';STOP
         ENDIF
!
  IF(iv==0 .AND. Alt(j,i) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  ';STOP
         ENDIF
  IF(iv==1 .AND. Alt(j,i) < 0.0 .AND. C_L(j) == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  '
  WRITE(99,*) 'ERROR IN SBR. TASK_3:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  ';STOP
         ENDIF
!+++++
 enddo
!+++++
ENDIF
!
!
IF(ne==30450) then
  IF(JUNK <= 0.0 .AND. iv==0) THEN
   WRITE(99,*) 'ERROR in sbr. TASK_3: The maximum thickness (or mass) must be > 0 '
   WRITE(99,*) 'in time- history L_h==8 ======= JOB ABORTED ==================== ';stop
  endif
  IF(JUNK <= 0.0 .AND. iv==1) THEN
   WRITE(* ,*) 'ERROR in sbr. TASK_3: The maximum thickness (or mass) must be > 0 '
   WRITE(* ,*) 'in time- history L_h==8 ======= JOB ABORTED ==================== '
   WRITE(99,*) 'ERROR in sbr. TASK_3: The maximum thickness (or mass) must be > 0 '
   WRITE(99,*) 'in time- history L_h==8 ======= JOB ABORTED ==================== ';stop
  endif
ENDIF
!
!
IF(ne==30456) THEN
IF(iv==0 .AND. (imake==0 .AND. iexte==0)) THEN
WRITE(99,*) 'ERROR in sbr. TASK_3: One of the kws Make_Model and '
WRITE(99,*) 'External_Model must be active  before the kw ICE_*  '
WRITE(99,*) 'can be activated ======== JOB ABORTED ============= ';STOP
ENDIF
IF(iv==1 .AND. (imake==0 .AND. iexte==0)) THEN
WRITE(* ,*) 'ERROR in sbr. TASK_3: One of the kws Make_Model and '
WRITE(* ,*) 'External_Model must be active  before the kw ICE_*  '
WRITE(* ,*) 'can be activated ======== JOB ABORTED ============= '
WRITE(99,*) 'ERROR in sbr. TASK_3: One of the kws Make_Model and '
WRITE(99,*) 'External_Model must be active  before the kw ICE_*  '
WRITE(99,*) 'can be activated ======== JOB ABORTED ============= ';STOP
ENDIF
ENDIF
!
!
IF(ne==30457) THEN
IF(iv==0 .AND. iadhoc ==1) THEN
WRITE(99,*) 'ERROR in sbr. TASK_3: Ad_Hoc and ICE_* are mutually '
WRITE(99,*) 'exclusive. Both are active ==== JOB ABORTED ======= ';STOP
ENDIF
IF(iv==1 .AND. iadhoc ==1) THEN
WRITE(* ,*) 'ERROR in sbr. TASK_3: Ad_Hoc and ICE_* are mutually '
WRITE(* ,*) 'exclusive. Both are active ==== JOB ABORTED ======= '
WRITE(99,*) 'ERROR in sbr. TASK_3: Ad_Hoc and ICE_* are mutually '
WRITE(99,*) 'exclusive. Both are active ==== JOB ABORTED ======= ';STOP
ENDIF
ENDIF
!
!
  IF(ne==30458) THEN
   IF  (n_iicest >=2 .and. iv==0) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw ICE_* must appear once'
   WRITE(99,*) '======= JOB ABORTED ========= JOB ABORTED ========';STOP
   ENDIF
   IF  (n_iicest >=2 .and. iv==1) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw ICE_* must appear once'
   WRITE(99,*) '======= JOB ABORTED ========= JOB ABORTED ========'
   WRITE(* ,*) 'ERROR in sbr TASK_3: The kw ICE_* must appear once'
   WRITE(* ,*) '======= JOB ABORTED ========= JOB ABORTED ========';STOP
   ENDIF
  ENDIF
!
!
IF(ne==30470) then
 IF ( (IICE < 0 .or. IICE > 1) .and. IV==0) THEN
  WRITE(99,*)'ERROR IN SBR. TASK_3: Expecting 0 or 1 for switch IICE'
 WRITE(99,*)'============ JOB ABORTED =============================='
 ENDIF
 IF ( (IICE  < 0 .or. IICE > 1) .and. IV==1) THEN
  WRITE(99,*)'ERROR IN SBR. TASK_3: Expecting 0 or 1 for a switch IICE'
 WRITE(99,*)'============ JOB ABORTED ================================'
   WRITE(*,*)'ERROR IN SBR. TASK_3: Expecting 0 or 1 for a switch IICE'
 WRITE(*,*)'============ JOB ABORTED =================================';STOP
 ENDIF
endif
!
!
IF(ne==30480) then
 IF ( (FRTD < 0 .or. FRTD > 1) .and. IV==0) THEN
  WRITE(99,*)'ERROR IN SBR. TASK_3: Expecting 0 or 1 for switch FRTD'
 WRITE(99,*)'============ JOB ABORTED =============================='
 ENDIF
 IF ( (FRTD  < 0 .or. FRTD > 1) .and. IV==1) THEN
  WRITE(99,*)'ERROR IN SBR. TASK_3: Expecting 0 or 1 for switch FRTD'
 WRITE(99,*)'============ JOB ABORTED =============================='
   WRITE(*,*)'ERROR IN SBR. TASK_3: Expecting 0 or 1 for switch FRTD'
 WRITE(*,*)'============ JOB ABORTED ===============================';STOP
 ENDIF
endif
!
!
IF(ne==30500) then
 IF ( (F7T8 < 0 .or. F7T8 > 1) .and. IV==0) THEN
  WRITE(99,*)'ERROR IN SBR. TASK_3: Expecting 0 or 1 for switch F7T8'
 WRITE(99,*)'============ JOB ABORTED =============================='
 ENDIF
 IF ( (F7T8 < 0 .or. F7T8 > 1) .and. IV==1) THEN
  WRITE(99,*)'ERROR IN SBR. TASK_3: Expecting 0 or 1 for switch F7T8'
 WRITE(99,*)'============ JOB ABORTED =============================='
   WRITE(*,*)'ERROR IN SBR. TASK_3: Expecting 0 or 1 for switch F7T8'
 WRITE(*,*)'============ JOB ABORTED ===============================';STOP
 ENDIF
endif
!
!
IF(ne==30510) then
                        If (Tau_8 <=0.0 .and. IV==0) Then
   write(99,*) 'ERROR IN SBR. TASK_3: The length of the loading phase   '
   write(99,*) 'must be > 0 kyrs in time-history #8 **** JOB ABORTED ***** '; STOP
                        Endif
!
                        If (Tau_8 <=0.0 .and. IV==1) Then
   write(99,*) 'ERROR IN SBR. TASK_3: The length of the loading phase   '
   write(99,*) 'must be > 0 kyrs in time-history #8 **** JOB ABORTED ***** '
   write(* ,*) 'ERROR IN SBR. TASK_3: The length of the loading phase   '
   write(* ,*) 'must be > 0 kyrs in time-history #8 **** JOB ABORTED ***** '; STOP
                        Endif
endif
!
!
IF(ne==30520) then
          If ((scale < 0.01 .OR. SCALE > 100.) .and. IV==0) Then
   write(99,*) 'ERROR IN SBR. TASK_3: The scaling factor SCALE must be > 0'
   write(99,*) '======= JOB ABORTED =======  ======= JOB ABORTED =======  '; STOP
          Endif
!
          If ((scale < 0.01 .OR. SCALE > 100.) .and. IV==0) Then
   write(99,*) 'ERROR IN SBR. TASK_3: The scaling factor SCALE must be > 0'
   write(99,*) '======= JOB ABORTED =======  ======= JOB ABORTED =======  '
   write(*,*) 'ERROR IN SBR. TASK_3: The scaling factor SCALE must be > 0'
   write(*,*) '======= JOB ABORTED =======  ======= JOB ABORTED =======  '; STOP
         Endif
endif
!
!
IF(ne==30530) then
IF(keyword(1:5)=='ICE_1') then
                        If ( (kp >= 4 .or. kp <= 0) .and. IV==0) Then
   write(99,*) 'ERROR IN SBR. TASK_3: the number of sub-aggregates of ICE1'
   write(99,*) 'must be in the range [0:3]  ======= JOB ABORTED =======  '; STOP
                        Endif
!
                        If ( (kp >= 4 .or. kp <= 0) .and. IV==1) Then
   write(99,*) 'ERROR IN SBR. TASK_3: the number of sub-aggregates of ICE1'
   write(99,*) 'must be in the range [0:3]  ======= JOB ABORTED =======  '
   write(* ,*) 'ERROR IN SBR. TASK_3: the number of sub-aggregates of ICE1'
   write(* ,*) 'must be in the range [0:3]  ======= JOB ABORTED =======  '; STOP
                        Endif
   endif
!
IF(keyword(1:5)=='ICE_2') then
                        If ( (kp >= 4 .or. kp <= 0) .and. IV==0) Then
   write(99,*) 'ERROR IN SBR. TASK_3: the number of sub-aggregates of ICE2'
   write(99,*) 'must be in the range [0:3]  ======= JOB ABORTED =======  '; STOP
                        Endif
!
                        If ( (kp >= 4 .or. kp <= 0) .and. IV==1) Then
   write(99,*) 'ERROR IN SBR. TASK_3: the number of sub-aggregates of ICE2'
   write(99,*) 'must be in the range [0:3]  ======= JOB ABORTED =======  '
   write(* ,*) 'ERROR IN SBR. TASK_3: the number of sub-aggregates of ICE2'
   write(* ,*) 'must be in the range [0:3]  ======= JOB ABORTED =======  '; STOP
                        Endif
   endif
!
IF(keyword(1:5)=='ICE_3') then
                        If ( (kp >= 10 .or. kp <= 0) .and. IV==0) Then
   write(99,*) 'ERROR IN SBR. TASK_3: the number of sub-aggregates of ICE3'
   write(99,*) 'must be in the range [0:3]  ======= JOB ABORTED =======  '; STOP
                        Endif
!
                        If ( (kp >= 10 .or. kp <= 0) .and. IV==1) Then
   write(99,*) 'ERROR IN SBR. TASK_3: the number of sub-aggregates of ICE3'
   write(99,*) 'must be in the range [0:3]  ======= JOB ABORTED =======  '
   write(* ,*) 'ERROR IN SBR. TASK_3: the number of sub-aggregates of ICE3'
   write(* ,*) 'must be in the range [0:3]  ======= JOB ABORTED =======  '; STOP
                        Endif
   endif
endif
!
!
if(ne==30540) then
   do kk=1,kp
!
       if(iv==0 .and. sub_agg(kk) /= 'ice1.dat' .and. sub_agg(kk) /= 'ice1.nam'  &
                .and. sub_agg(kk) /= 'ice1.eup' .and. sub_agg(kk) /= 'ice1.gro') then
   write(99,*)'error in sbr. task_3: the name of the sub-aggregate of ice1 must'
   write(99,*)'be one among  --->  ice1.dat, ice1.nam, ice1.eup, ice1.gro <--- '
   write(99,*)'========= job aborted =========   ========= job aborted ======= ';stop
                 endif

       if(iv==1 .and. sub_agg(kk) /= 'ice1.dat' .and. sub_agg(kk) /= 'ice1.nam'  &
                .and. sub_agg(kk) /= 'ice1.eup' .and. sub_agg(kk) /= 'ice1.gro') then
   write(99,*)'error in sbr. task_3: the name of the sub-aggregate of ice1 must'
   write(99,*)'be one among  --->  ice1.dat, ice1.nam, ice1.eup, ice1.gro <--- '
   write(99,*)'========= job aborted =========   ========= job aborted ======= '
   write(* ,*)'error in sbr. task_3: the name of the sub-aggregate of ice1 must'
   write(* ,*)'be one among  --->  ice1.dat, ice1.nam, ice1.eup, ice1.gro <--- '
   write(* ,*)'========= job aborted =========   ========= job aborted ======= ';stop
                 endif
!
   end do
endif
!
!
if(ne==30541) then
   do kk=1,kp
       if(iv==0 .and. sub_agg(kk) /= 'ice2.dat' .and. sub_agg(kk) /= 'ice2.nam'  &
                .and. sub_agg(kk) /= 'ice2.eup' .and. sub_agg(kk) /= 'ice2.gro') then
   write(99,*)'error in sbr. task_3: the name of the sub-aggregate of ice1 must'
   write(99,*)'be one among  --->  ice2.dat, ice2.nam, ice2.eup, ice2.gro <--- '
   write(99,*)'========= job aborted =========   ========= job aborted ======= ';stop
                 endif

       if(iv==1 .and. sub_agg(kk) /= 'ice2.dat' .and. sub_agg(kk) /= 'ice2.nam'  &
                .and. sub_agg(kk) /= 'ice2.eup' .and. sub_agg(kk) /= 'ice2.gro') then
   write(99,*)'error in sbr. task_3: the name of the sub-aggregate of ice2 must'
   write(99,*)'be one among  --->  ice2.dat, ice2.nam, ice2.eup, ice2.gro <--- '
   write(99,*)'========= job aborted =========   ========= job aborted ======= '
   write(* ,*)'error in sbr. task_3: the name of the sub-aggregate of ice2 must'
   write(* ,*)'be one among  --->  ice2.dat, ice2.nam, ice2.eup, ice2.gro <--- '
   write(* ,*)'========= job aborted =========   ========= job aborted ======= ';stop
                 endif
   end do
endif
!
!
IF(ne==30542) then
   do kk=1,kp
       if(iv==0 .and. sub_agg(kk) /= 'ice3.dat' .and. sub_agg(kk) /= 'ice3.and'  &
                .and. sub_agg(kk) /= 'ice3.ant' .and. sub_agg(kk) /= 'ice3.bal'  &
                .and. sub_agg(kk) /= 'ice3.bri' .and. sub_agg(kk) /= 'ice3.gro'  &
                .and. sub_agg(kk) /= 'ice3.ike' .and. sub_agg(kk) /= 'ice3.nam'  &
                .and. sub_agg(kk) /= 'ice3.pol' .and. sub_agg(kk) /= 'ice3.sib' ) then
   write(99,*)'error in sbr. task_3:  the name of the sub-aggregate of ice3 must'
   write(99,*)'be one among  ---> ice3.dat, .and., .ant, .bal, .bri, .gro, .ike,'
   write(99,*)'.nam, .pol, .sib  <--- ============== job aborted ============== ';stop
                 endif

       if(iv==1 .and. sub_agg(kk) /= 'ice3.dat' .and. sub_agg(kk) /= 'ice3.and'  &
                .and. sub_agg(kk) /= 'ice3.ant' .and. sub_agg(kk) /= 'ice3.bal'  &
                .and. sub_agg(kk) /= 'ice3.bri' .and. sub_agg(kk) /= 'ice3.gro'  &
                .and. sub_agg(kk) /= 'ice3.ike' .and. sub_agg(kk) /= 'ice3.nam'  &
                .and. sub_agg(kk) /= 'ice3.pol' .and. sub_agg(kk) /= 'ice3.sib' ) then
   write(99,*)'error in sbr. task_3:  the name of the sub-aggregate of ice3 must'
   write(99,*)'be one among  ---> ice3.dat, .and., .ant, .bal, .bri, .gro, .ike,'
   write(99,*)'.nam, .pol, .sib  <--- ============== job aborted ============== '
   write(* ,*)'error in sbr. task_3:  the name of the sub-aggregate of ice33 must'
   write(* ,*)'be one among  ---> ice3.dat, .and., .ant, .bal, .bri, .gro, .ike,'
   write(* ,*)'.nam, .pol, .sib  <--- ============== job aborted ============== ';stop
                 endif
   end do
ENDIF
!
!
IF(ne==30550) then
  IF(iv==0) then
   WRITE(99,*)'ERROR IN SBR. TASK_3: If ICE*.DAT is in the list of sub-aggregates,'
   WRITE(99,*)'no other sub-aggregate can be present ======== JOB ABORTED ========';stop
            endif
  IF(iv==1) then
   WRITE(99,*)'ERROR IN SBR. TASK_3: If ICE*.DAT is in the list of sub-aggregates,'
   WRITE(99,*)'no other sub-aggregate can be present ======== JOB ABORTED ========'
   WRITE(* ,*)'ERROR IN SBR. TASK_3: If ICE*.DAT is in the list of sub-aggregates,'
   WRITE(* ,*)'no other sub-aggregate can be present ======== JOB ABORTED ========';STOP
            endif
endif
!
!
IF(ne==30560) then
  IF(iv==0) then
   WRITE(99,*)'ERROR IN SBR. TASK_3: A sub-aggregate name must appear once in the'
   WRITE(99,*)'sub-aggregates list for ICE1 ================ JOB ABORTED ========';stop
            endif
  IF(iv==1) then
   WRITE(99,*)'ERROR IN SBR. TASK_3: A sub-aggregate name must appear once in the'
   WRITE(99,*)'sub-aggregates list for ICE1 ================ JOB ABORTED ========'
   WRITE(* ,*)'ERROR IN SBR. TASK_3: A sub-aggregate name must appear once in the'
   WRITE(* ,*)'sub-aggregates list for ICE1 ================ JOB ABORTED ========';STOP
            endif
endif
!
!
IF(ne==30570) then
          WRITE(99,*)'WARNING from SBR.  TASK_3: The conversion from Rectangles'
          WRITE(99,*)'(CL=50) to discs (CL=10) is not possible with ICE3, since'
          WRITE(99,*)'ICE3 is already in the form CL=10 ======================='
IF(iv==1) THEN
          WRITE(* ,*)'WARNING from SBR.  TASK_3: The conversion from Rectangles'
          WRITE(* ,*)'(CL=50) to discs (CL=10) is not possible with ICE3, since'
          WRITE(* ,*)'ICE3 is already in the form CL=10 ======================='
ENDIF
endif
!
!
IF(ne==30580) then
  If(iadhoc==0 .AND. iicest==0 .and. iv==0) Then
  Write(99,*)'ERROR in TASK_3: One of the kws Ad_Hoc and ICE_* must be activated'
  Write(99,*)'before Local_Study can be executed ======== JOB ABORTED ======';Stop
  Endif
  If(iadhoc==0 .AND. iicest==0 .and. iv==1) Then
  Write(99,*)'ERROR in TASK_3: One of the kws Ad_Hoc and ICE_* must be activated'
  Write(99,*)'before Local_Study can be executed ======== JOB ABORTED ======'
  Write(* ,*)'ERROR in TASK_3: One of the kws Ad_Hoc and ICE_* must be activated'
  Write(* ,*)'before Local_Study can be executed ======== JOB ABORTED ======'; STOP
  Endif
endif
!
!
IF(ne==30581) then
  If(iadhoc==0 .AND. iicest==0 .and. iv==0) Then
  Write(99,*)'ERROR in TASK_3: One of the kws Ad_Hoc and ICE_* must be activated'
  Write(99,*)'before Global_Study can be executed ======== JOB ABORTED ======';Stop
  Endif
  If(iadhoc==0 .AND. iicest==0 .and. iv==1) Then
  Write(99,*)'ERROR in TASK_3: One of the kws Ad_Hoc and ICE_* must be activated'
  Write(99,*)'before Global_Study can be executed ======== JOB ABORTED ======'
  Write(* ,*)'ERROR in TASK_3: One of the kws Ad_Hoc and ICE_* must be activated'
  Write(* ,*)'before Global_Study can be executed ======== JOB ABORTED ======'; STOP
  Endif
endif
!
!
!
  IF(ne==30582) THEN
   IF  (n_inte >=2 .and. iv==0) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw Local_Study must appear once'
   WRITE(99,*) '= ========= JOB ABORTED ========= JOB ABORTED ==============';STOP
   ENDIF
   IF  (n_inte >=2 .and. iv==1) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw Local_Study must appear once'
   WRITE(99,*) '= ========= JOB ABORTED ========= JOB ABORTED =============='
   WRITE(* ,*) 'ERROR in sbr TASK_3: The kw Local_Study must appear once'
   WRITE(* ,*) '= ========= JOB ABORTED ========= JOB ABORTED ==============';STOP
   ENDIF
  ENDIF
!
!
  IF(ne==30583) THEN
   IF  (n_este >=2 .and. iv==0) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw Global_Study must appear once'
   WRITE(99,*) '= ========= JOB ABORTED ========= JOB ABORTED ==============';STOP
   ENDIF
   IF  (n_este >=2 .and. iv==1) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: The kw Global_Study must appear once'
   WRITE(99,*) '= ========= JOB ABORTED ========= JOB ABORTED =============='
   WRITE(* ,*) 'ERROR in sbr TASK_3: The kw Global_Study must appear once'
   WRITE(* ,*) '= ========= JOB ABORTED ========= JOB ABORTED ==============';STOP
   ENDIF
  ENDIF
!
!
    IF(ne==30584) THEN
    IF(n_inte ==1 .and. IV==0) then
      WRITE(99,*)'ERROR IN SBR. TASK_3:  The KWs Global_ and Local_Study are'
      WRITE(99,*)'mutually exclusive. Both are active. ***** JOB ABORTED **********';STOP
    endif
    IF(n_inte ==1 .and. IV==1) then
      WRITE(99,*)'ERROR IN SBR. TASK_3:  The KWs Global_ and Local_Study are'
      WRITE(99,*)'mutually exclusive. Both are active. ***** JOB ABORTED **********'
      WRITE(* ,*)'ERROR IN SBR. TASK_3:  The KWs Global_ and Local_Study are'
      WRITE(* ,*)'mutually exclusive. Both are active. ***** JOB ABORTED **********';STOP
    endif
    ENDIF
!
!
IF(ne==30590) then
 If ((type_inten<=0.or.type_inten>4).and.iv==0) then
  Write(99,*) 'ERROR IN SBR. TASK_3:  the Local_Study label '
  Write(99,*) 'must vary in the range [1:4] *** JOB ABORTED *** '; stop
 Endif
 If ((type_inten<=0.or.type_inten>4).and.iv==1) then
  Write(99,*) 'ERROR IN SBR. TASK_3:  the Local_Study label '
  Write(99,*) 'must vary in the range [1:4] *** JOB ABORTED *** '
  Write(* ,*) 'ERROR IN SBR. TASK_3:  the Local_Study label '
  Write(* ,*) 'must vary in the range [1:4] *** JOB ABORTED *** '; stop
 Endif
ENDIF
!
!
IF(ne==30591) then
 If ((type_esten<=0.or.type_esten>2).and.iv==0) then
  Write(99,*) 'ERROR IN SBR. TASK_3:  the Global_Study label '
  Write(99,*) 'must vary in the range [1:2] *** JOB ABORTED *** '; stop
 Endif
 If ((type_esten<=0.or.type_esten>2).and.iv==1) then
  Write(99,*) 'ERROR IN SBR. TASK_3:  the Global_Study label '
  Write(99,*) 'must vary in the range [1:2] *** JOB ABORTED *** '
  Write(* ,*) 'ERROR IN SBR. TASK_3:  the Global_Study label '
  Write(* ,*) 'must vary in the range [1:2] *** JOB ABORTED *** '; stop
 Endif
ENDIF
!
!
IF(ne==30600) then
If (IR==0.and.dIR==0.and.it==0.and.dit==0  &
                            .and.IG==0.and.dIG==0.and.iv==0) then 
Write(99,*) 'ERROR IN SBR. TASK_3: At least one of the switches IR, DIR, etc. must '
Write(99,*) 'be =1 in order to perform an Local_Study. **** JOB ABORTED ********** ';stop
Endif          
If (IR==0.and.dIR==0.and.it==0.and.dit==0  &
                            .and.IG==0.and.dIG==0.and.iv==1) then 
Write(99,*) 'ERROR IN SBR. TASK_3: At least one of the switches IR, DIR, etc. must '
Write(99,*) 'be =1 in order to perform an Local_Study. **** JOB ABORTED ********** '
Write(* ,*) 'ERROR IN SBR. TASK_3: At least one of the switches IR, DIR, etc. must '
Write(* ,*) 'be =1 in order to perform an Local_Study. **** JOB ABORTED ********** ';stop
endif
endif
!
!
IF(ne==30610) THEN
 If (iv==0 .and. T1 > T2) then
  Write(99,*) 'ERROR IN SBR. TASK_3: In kw Local_Study T1 must be <= T2  '
  Write(99,*) '(both T1 and T2 must be expressed in kyrs Before Present, BP) '
  WRITE(99,*) 'רררררררררר JOB ABORTED רררררררררר רררר JOB ABORTED רררררררררר ';STOP
 Endif
 If (iv==1 .and. T1 > T2) then
  Write(* ,*) 'ERROR IN SBR. TASK_3: In kw Local_Study T1 must be <= T2  '
  Write(* ,*) '(both T1 and T2 must be expressed in kyrs Before Present, BP) '
  WRITE(* ,*) 'רררררררררר JOB ABORTED רררררררררר רררר JOB ABORTED רררררררררר '
  Write(99,*) 'ERROR IN SBR. TASK_3: In kw Local_Study T1 must be <= T2  '
  Write(99,*) '(both T1 and T2 must be expressed in kyrs Before Present, BP) '
  WRITE(99,*) 'רררררררררר JOB ABORTED רררררררררר רררר JOB ABORTED רררררררררר ';STOP
 Endif
ENDIF
!
!
IF(ne==30620) THEN
 If (iv==0 .and. DTI <= 0.0) then
  Write(99,*) 'ERROR IN SBR. TASK_3: The time increment must be >  0 '
  Write(99,*) 'in the Local Study #1 **** JOB ABORTED ********** ';stop
 Endif
 If (iv==1 .and. DTI <= 0.0) then
  Write(99,*) 'ERROR IN SBR. TASK_3: The time increment must be >  0 '
  Write(99,*) 'in the Local Study #1 **** JOB ABORTED ********** '
  Write(* ,*) 'ERROR IN SBR. TASK_3: The time increment must be >  0 '
  Write(* ,*) 'in the Local Study #1 **** JOB ABORTED ********** ';stop
 Endif
endif
!
!
       IF(ne==30625) THEN
       IF(iv==0 .AND. ( Long_obs(kk) < 0.d0 .or. Long_obs(kk) > 360.d0)) THEN
         WRITE(99,*) 'ERROR in sbr TASK_3: The longitude of one observer is'
         WRITE(99,*) 'out of bounds ========== JOB ABORTED ================'
       ENDIF
       IF(iv==0 .AND. ( Long_obs(kk) < 0.d0 .or. Long_obs(kk) > 360.d0)) THEN
         WRITE(* ,*) 'ERROR in sbr TASK_3: The longitude of one observer is'
         WRITE(* ,*) 'out of bounds ========== JOB ABORTED ================'
         WRITE(99,*) 'ERROR in sbr TASK_3: The longitude of one observer is'
         WRITE(99,*) 'out of bounds ========== JOB ABORTED ================';STOP
       ENDIF
       IF(iv==1 .AND. ( Cola_obs(kk) < 0.d0 .or. Cola_obs(kk) > 180.d0)) THEN
         WRITE(99,*) 'ERROR in sbr TASK_3: The colatitude of one observer is'
         WRITE(99,*) 'out of bounds ========== JOB ABORTED ================'
       ENDIF
       IF(iv==1 .AND. ( Cola_obs(kk) < 0.d0 .or. Cola_obs(kk) > 180.d0)) THEN
         WRITE(* ,*) 'ERROR in sbr TASK_3: The colatitude of one observer is'
         WRITE(* ,*) 'out of bounds ========== JOB ABORTED ================'
         WRITE(99,*) 'ERROR in sbr TASK_3: The colatitude of one observer is'
         WRITE(99,*) 'out of bounds ========== JOB ABORTED ================';STOP
       ENDIF
       ENDIF
!
!
IF(ne==30630) THEN
 IF(NOBS>NOBS_MAX .and. IV==0) THEN 
  WRITE (99,*) 'ERROR IN SBR. TASK_3: The number of observers in file',file_sparsi
  write (99,*) 'exceeds the maximum allowed', nobs_max,' *** JOB ABORTED *********'; stop 
 ENDIF 
 IF(NOBS>NOBS_MAX .and. IV==1) THEN 
  WRITE (99,*) 'ERROR IN SBR. TASK_3: The number of observers in file',file_sparsi
  write (99,*) 'exceeds the maximum allowed', nobs_max,' *** JOB ABORTED *********' 
  WRITE (* ,*) 'ERROR IN SBR. TASK_3: The number of observers in file',file_sparsi
  write (* ,*) 'exceeds the maximum allowed', nobs_max,' *** JOB ABORTED *********'; stop   
 ENDIF
endif
!
!
IF(ne==30640) THEN
IF(NOBS==0 .and. IV==0) THEN
  WRITE (99,*) 'ERROR IN SBR. TASK_3: It appears that the file',file_sparsi,'is empty'
  write (99,*) 'or that its name is not properly aligned in file task_3.dat. JOB ABORTED ****'; stop
 ENDIF 
 IF(NOBS==0 .and. IV==1) THEN 
  WRITE (99,*) 'ERROR IN SBR. TASK_3: It appears that the file',file_sparsi,'is empty'
  write (99,*) 'or that its name is not properly aligned in file task_3.dat. JOB ABORTED ****'
  WRITE (*, *) 'ERROR IN SBR. TASK_3: It appears that the file',file_sparsi,'is empty'
  write (*, *) 'or that its name is not properly aligned in file task_3.dat. JOB ABORTED ****'; stop
 ENDIF
ENDIF
!
!
IF(ne==30650) THEN
  If(iv==0.and.(L1<0.d0.OR.L1>360.d0.OR.L2<0.d0.OR.L2>360.d0.OR.L1>=L2.OR.EPL<=0.d0)) Then
   Write(99,*) 'ERROR IN SBR. TASK_3:        One of the forbidden conditions: '
   Write(99,*) 'L1<0 deg; L1>360 deg; L2<0 deg; L2>360 deg; L1 > L2; L1 = L2; '
   Write(99,*) 'EPL <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ********************************************';stop 
  Endif
  If(iv==1.and.(L1<0.d0.OR.L1>360.d0.OR.L2<0.d0.OR.L2>360.d0.OR.L1>=L2.OR.EPL<=0.d0)) Then 
   Write(99,*) 'ERROR IN SBR. TASK_3:        One of the forbidden conditions: '
   Write(99,*) 'L1<0 deg; L1>360 deg; L2<0 deg; L2>360 deg; L1 > L2; L1 = L2; '
   Write(99,*) 'EPL <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ********************************************'
   Write(* ,*) 'ERROR IN SBR. TASK_3:        One of the forbidden conditions: '
   Write(* ,*) 'L1<0 deg; L1>360 deg; L2<0 deg; L2>360 deg; L1 > L2; L1 = L2; '
   Write(* ,*) 'EPL <= 0 has been met. Please Check the input file task_2.dat '          
   Write(* ,*) '***** JOB ABORTED ********************************************';stop  
  Endif 
ENDIF
!
!
IF(ne==30660) THEN
  If(iv==0.and.(C1<0.d0.OR.C1>180.d0.OR.C2<0.d0.OR.C2>180.d0.OR.C1>=C2.OR.EPC<=0.d0)) Then
   Write(99,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(99,*) 'C1<0 deg; C1>180 deg; C2<0 deg; C2>180 deg; C1 > C2; C1 = C2; ' 
   Write(99,*) 'EPC <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ********************************************';stop 
  Endif
  If(iv==1.and.(C1<0.d0.OR.C1>180.d0.OR.C2<0.d0.OR.C2>180.d0.OR.C1>=C2.OR.EPC<=0.d0)) Then 
   Write(99,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(99,*) 'C1<0 deg; C1>180 deg; C2<0 deg; C2>180 deg; C1 > C2; C1 = C2; ' 
   Write(99,*) 'EPC <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ******************************************* ' 
   Write(* ,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(* ,*) 'C1<0 deg; C1>180 deg; C2<0 deg; C2>180 deg; C1 > C2; C1 = C2; ' 
   Write(* ,*) 'EPC <= 0 has been met. Please Check the input file task_2.dat '          
   Write(* ,*) '***** JOB ABORTED ******************************************* ';stop  
  Endif
ENDIF
!
!
 IF(ne==30661) THEN
      DO J=1, mel_eff
         IF(iv ==0 .and. I_THI==1 .and. c_l(j) /= 10 .AND. c_l(j) /=50) THEN
         WRITE(*,*) 'ERROR in sbr TASK_3: The computation of the primary load'
         WRITE(*,*) 'thickness  requires a primary load of  the type 10 or 50'
         WRITE(*,*) '*********** JOB ABORTED ********************************'
         STOP
         ENDIF
         IF(iv ==1 .and. I_THI==1 .and. c_l(j) /= 10 .AND. c_l(j) /=50) THEN
         WRITE(* ,*) 'ERROR in sbr TASK_3: The computation of the primary load'
         WRITE(* ,*) 'thickness  requires a primary load of  the type 10 or 50'
         WRITE(* ,*) '*********** JOB ABORTED ********************************'
         WRITE(99,*) 'ERROR in sbr TASK_3: The computation of the primary load'
         WRITE(99,*) 'thickness  requires a primary load of  the type 10 or 50'
         WRITE(99,*) '*********** JOB ABORTED ********************************'
         STOP
         ENDIF
      IF(iv ==0 .AND. I_THI==1 .AND. IOC(j)==1) THEN
         WRITE(*,*) 'ERROR in sbr TASK_3: The computation of the primary load'
         WRITE(*,*) 'thickness requires IOC=0 (no compensation on the oceans)'
         WRITE(*,*) '*********** JOB ABORTED ********************************'
         STOP
      ENDIF
      IF(iv ==1 .AND. I_THI==1 .AND. IOC(j)==1) THEN
         WRITE(* ,*) 'ERROR in sbr TASK_3: The computation of the primary load'
         WRITE(* ,*) 'thickness requires IOC=0 (no compensation on the oceans)'
         WRITE(* ,*) '*********** JOB ABORTED ********************************'
         WRITE(99,*) 'ERROR in sbr TASK_3: The computation of the primary load'
         WRITE(99,*) 'thickness requires IOC=0 (no compensation on the oceans)'
         WRITE(99,*) '*********** JOB ABORTED ********************************'
         STOP
      ENDIF
      ENDDO
      IF(iv ==0 .AND. I_THI==1 .and. lmin/=2) THEN
         WRITE(*,*) 'ERROR in sbr TASK_3: The computation of the primary load'
         WRITE(*,*) 'thickness requires -----> l_min=2 in kw Harmonic_Degrees'
         WRITE(*,*) '*********** JOB ABORTED ********************************'
         STOP
      ENDIF
      IF(iv ==1 .AND. I_THI==1 .and. lmin/=2) THEN
         WRITE(* ,*) 'ERROR in sbr TASK_3: The computation of the primary load'
         WRITE(* ,*) 'thickness requires -----> l_min=2 in kw Harmonic_Degrees'
         WRITE(* ,*) '*********** JOB ABORTED ********************************'
         WRITE(99,*) 'ERROR in sbr TASK_3: The computation of the primary load'
         WRITE(99,*) 'thickness requires -----> l_min=2 in kw Harmonic_Degrees'
         WRITE(99,*) '*********** JOB ABORTED ********************************'
         STOP
      ENDIF
ENDIF
!
!
 IF (NE==30662) THEN
  IF(iv==0 .AND. sigma(0)==0.d0) THEN
      WRITE(99,*) 'ERROR from sbr ocean_correction: the degree zero'
      WRITE(99,*) 'load component is already equal to zero, so that'
      WRITE(99,*) 'there is no need to perform another compensation'
      WRITE(99,*) '============== JOB ABORTED ====================='
      STOP
 ENDIF
 IF(iv==1 .AND. sigma(0)==0.d0) THEN
      WRITE(* ,*) 'ERROR from sbr ocean_correction: the degree zero'
      WRITE(* ,*) 'load component is already equal to zero, so that'
      WRITE(* ,*) 'there is no need to perform another compensation'
      WRITE(* ,*) '============== JOB ABORTED ====================='
      WRITE(99,*) 'ERROR from sbr ocean_correction: the degree zero'
      WRITE(99,*) 'load component is already equal to zero, so that'
      WRITE(99,*) 'there is no need to perform another compensation'
      WRITE(99,*) '============== JOB ABORTED ====================='
      STOP
 ENDIF
 ENDIF
!
!
IF(ne==30670) THEN
 IF(NOBS>NOBS_MAX .and. IV==0) THEN 
  WRITE (99,*) 'ERROR IN SBR. TASK_3: The number of observers exceeds the'
  write (99,*) 'maximum allowed = ', nobs_max,' *** JOB ABORTED *********'; stop
 ENDIF 
 IF(NOBS>NOBS_MAX .and. IV==1) THEN 
  WRITE (99,*) 'ERROR IN SBR. TASK_3: The number of observers exceeds the'
  write (99,*) 'maximum allowed = ', nobs_max,' *** JOB ABORTED *********'
  WRITE (* ,*) 'ERROR IN SBR. TASK_3: The number of observers exceeds the'
  write (* ,*) 'maximum allowed = ', nobs_max,' *** JOB ABORTED *********'; stop
 ENDIF
endif
!
!
IF(ne==30680) THEN
IF(IV==0) THEN
  WRITE (99,*) 'ERROR IN SBR. TASK_3: Error opening file ',file_sparsi
  write (99,*) '================= JOB ABORTED ======================== '; stop
 ENDIF 
 IF(IV==1) THEN
  WRITE (99,*) 'ERROR IN SBR. TASK_3: Error opening file ',file_sparsi
  write (99,*) '================= JOB ABORTED ======================== '
  WRITE (* ,*) 'ERROR IN SBR. TASK_3: Error opening file ',file_sparsi
  write (* ,*) '================= JOB ABORTED ======================== '; stop
 ENDIF
ENDIF
!
!
IF(ne==30690) THEN
If ( ilb/=1 .and. ilb/=2 .and. iv==0)  Then
 WRITE(99,*) 'ERROR in sbr TASK_3: The input file ', file_sparsi
 Write(99,*) 'must contain 1 or 2 in its first line. See the User guide'
 WRITE(99,*) 'JOB ABORTED ======= JOB ABORTED ======= JOB ABORTED ====='; Stop
Endif
If ( ilb/=1 .and. ilb/=2 .and. iv==1)  Then
 WRITE(99,*) 'ERROR in sbr TASK_3: The input file ', file_sparsi
 Write(99,*) 'must contain 1 or 2 in its first line. See the User guide'
 WRITE(99,*) 'JOB ABORTED ======= JOB ABORTED ======= JOB ABORTED ====='
 WRITE(* ,*) 'ERROR in sbr TASK_3: The input file ', file_sparsi
 Write(* ,*) 'must contain 1 or 2 in its first line. See the User guide'
 WRITE(* ,*) 'JOB ABORTED ======= JOB ABORTED ======= JOB ABORTED ====='; Stop
Endif
ENDIF
!
!
IF(ne==30700) THEN
  IF(iv==0 .and. &
  (Long_site_1(kk)==Long_site_2(kk).OR.Cola_site_1(kk)== Cola_site_2(kk))) THEN
            WRITE(99,*) 'ERROR in Sbr. TASK_3:  the two sites must have  '
            WRITE(99,*) 'distinct longitude and colatitide. JOB ABORTED  '; STOP
        ENDIF
  IF(iv==1 .and. &
  (Long_site_1(kk)==Long_site_2(kk).OR.Cola_site_1(kk)== Cola_site_2(kk))) THEN
            WRITE(99,*) 'ERROR in Sbr. TASK_3:  the two sites must have  '
            WRITE(99,*) 'distinct longitude and colatitide. JOB ABORTED  '
            WRITE(* ,*) 'ERROR in Sbr. TASK_3:  the two sites must have  '
            WRITE(* ,*) 'distinct longitude and colatitide. JOB ABORTED  '; STOP
        ENDIF
ENDIF
!
!
IF(ne==30710) THEN
  IF(iv==0 .AND. (long_site_1(kk) < 0.0 .OR. long_site_1(kk) > 360.0 .OR. &
                  cola_site_1(kk) < 0.0 .OR. cola_site_1(kk) > 180.0 )) THEN
  WRITE(99,*) 'ERROR from sbr TASK_3: The longitude or colatitude of site#1 '
  WRITE(99,*) 'is out of bounds ============ JOB ABORTED ================== ';STOP
            ENDIF
    IF(iv==1 .AND. (long_site_1(kk) < 0.0 .OR. long_site_1(kk) > 360.0 .OR. &
                  cola_site_1(kk) < 0.0 .OR. cola_site_1(kk) > 180.0 )) THEN
  WRITE(99,*) 'ERROR from sbr TASK_3: The longitude or colatitude of site#1 '
  WRITE(99,*) 'is out of bounds ============ JOB ABORTED ================== '
  WRITE(* ,*) 'ERROR from sbr TASK_3: The longitude of colatitude of site#1 '
  WRITE(* ,*) 'is out of bounds ============ JOB ABORTED ================== ';STOP
            ENDIF
  IF(iv==0 .AND. (long_site_2(kk) < 0.0 .OR. long_site_2(kk) > 360.0 .OR. &
                  cola_site_2(kk) < 0.0 .OR. cola_site_2(kk) > 180.0 )) THEN
  WRITE(99,*) 'ERROR from sbr TASK_3: The longitude or colatitude of site#2 '
  WRITE(99,*) 'is out of bounds ============ JOB ABORTED ================== ';STOP
            ENDIF
    IF(iv==1 .AND. (long_site_2(kk) < 0.0 .OR. long_site_2(kk) > 360.0 .OR. &
                  cola_site_2(kk) < 0.0 .OR. cola_site_2(kk) > 180.0 )) THEN
  WRITE(99,*) 'ERROR from sbr TASK_3: The longitude of colatitude of site#2 '
  WRITE(99,*) 'is out of bounds ============ JOB ABORTED ================== '
  WRITE(* ,*) 'ERROR from sbr TASK_3: The longitude of colatitude of site#2 '
  WRITE(* ,*) 'is out of bounds ============ JOB ABORTED ================== ';STOP
            ENDIF
ENDIF
!
!
IF(ne==30720) THEN
        IF ( iv == 0 .and. Name_1 (kk) == Name_2 (kk))THEN
            WRITE(99,*) 'Error in Sbr. TASK_3: the names of the two sites '
            WRITE(99,*) 'must be distinct ======= JOB ABORTED =========== '; STOP
        ENDIF

        IF ( iv == 1 .and. Name_1 (kk) == Name_2 (kk) )THEN
            WRITE(* ,*) 'Error in Sbr. TASK_3: the names of the two sites '
            WRITE(* ,*) 'must be distinct ======= JOB ABORTED =========== '
            WRITE(99,*) 'Error in Sbr. TASK_3: the names of the two sites '
            WRITE(99,*) 'must be distinct ======= JOB ABORTED =========== '; STOP
        ENDIF
ENDIF
!
!
IF(ne==30730) THEN
     IF(iv==0 .AND. NOBS == 0) THEN
     WRITE(99,*) 'ERROR in sbr TASK_3: There must be at least a couple of sites'
     WRITE(99,*) '================== JOB ABORTED ==============================';STOP
     ENDIF
     IF(iv==1 .AND. NOBS == 0) THEN
     WRITE(*, *) 'ERROR in sbr TASK_3: There must be at least a couple of sites'
     WRITE(*, *) '================== JOB ABORTED =============================='
     WRITE(99,*) 'ERROR in sbr TASK_3: There must be at least a couple of sites'
     WRITE(99,*) '================== JOB ABORTED ==============================';STOP
     ENDIF
ENDIF
!
!
IF(ne==30740) THEN
 IF(IV==0 .AND. ( JN<0 .OR. JN>3 )) THEN
   WRITE(99,*) 'ERROR in sbr TASK_3: Only 0<= JN <=3 is allowed'
   WRITE(99,*) '============= JOB ABORTED ====================='; STOP
 ENDIF
 IF(IV==1 .AND. ( JN<0 .OR. JN>3 )) THEN
   WRITE(* ,*) 'ERROR in sbr TASK_3: Only 0<= JN <=3 is allowed'
   WRITE(* ,*) '============= JOB ABORTED ====================='
   WRITE(99,*) 'ERROR in sbr TASK_3: Only 0<= JN <=3 is allowed'
   WRITE(99,*) '============= JOB ABORTED ====================='; STOP
 ENDIF
ENDIF
!
!
IF(ne==30770) THEN
IF (iv==0 .and. (LDE<LMIN.OR.LDE>LMAX.OR.MDE<0.OR.LDE<0.OR.IABS(MDE)>LDE))  THEN
Write(99,*)'ERROR IN SBR. TASK_3: One of the following forbidden conditions has been met:'
Write(99,*)'LDE < LMIN; LDE > LMAX; MDE < 0; |MDE| > LDE. Check the input file task_2.dat' 
Write(99,*)'************** JOB ABORTED **************************************************';stop   
ENDIF 
IF (iv==1 .and. (LDE<LMIN.OR.LDE>LMAX.OR.MDE<0.OR.LDE<0.OR.IABS(MDE)>LDE))  THEN 
Write(99,*)'ERROR IN SBR. TASK_3: One of the following forbidden conditions has been met:'
Write(99,*)'LDE < LMIN; LDE > LMAX; MDE < 0; |MDE| > LDE. Check the input file task_2.dat' 
Write(99,*)'************** JOB ABORTED **************************************************'  
Write(* ,*)'ERROR IN SBR. TASK_3: One of the following forbidden conditions has been met:'
Write(*, *)'LDE < LMIN; LDE > LMAX; MDE < 0; |MDE| > LDE. Check the input file task_2.dat' 
Write(* ,*)'************** JOB ABORTED **************************************************';stop 
ENDIF
ENDIF
!
!
IF(ne==30780) THEN
IF (iv==0 .and. Inorm /= 0 .and. Inorm /= 1) Then
Write(99,*)'ERROR IN SBR. TASK_3: 0 or 1 is expected for the switch Inorm. JOB ABORTED'; Stop
endif
IF (iv==1 .and. Inorm /= 0 .and. Inorm /= 1) Then 
Write(99,*)'ERROR IN SBR. TASK_3: 0 or 1 is expected for the switch Inorm. JOB ABORTED'
Write(* ,*)'ERROR IN SBR. TASK_3: 0 or 1 is expected for the switch Inorm. JOB ABORTED'; Stop
endif
ENDIF
!
!
IF(ne==30790) THEN
 If (iv==0 .and. T1 > T2) then
  Write(99,*) 'ERROR IN SBR. TASK_3: T1 must be <= T2 in KW Global_Study'
  Write(99,*) '***** JOB ABORTED ********************************************';stop    
 Endif
 If (iv==1 .and. T1 > T2) then
  Write(99,*) 'ERROR IN SBR. TASK_3: T1 must be <= T2 in KW Global_Study'
  Write(99,*) '***** JOB ABORTED ********************************************'   
  Write(* ,*) 'ERROR IN SBR. TASK_3: T1 must be <= T2 in KW Global_Study'
  Write(* ,*) '***** JOB ABORTED ********************************************';stop   
 Endif
 If (iv==0 .and. DTI <= 0.0) then
  Write(99,*) 'ERROR IN SBR. TASK_3: T1 must be <= T2 in KW Global_Study'
  Write(99,*) '***** JOB ABORTED ********************************************';stop    
 Endif
 If (iv==1 .and. DTI <= 0.0) then
  Write(99,*) 'ERROR IN SBR. TASK_3: DTI must be > 0. in KW Global_Study'
  Write(99,*) '***** JOB ABORTED ********************************************'
  Write(* ,*) 'ERROR IN SBR. TASK_3: DTI must be > 0. in KW Global_Study'
  Write(* ,*) '***** JOB ABORTED ********************************************';stop   
 Endif
ENDIF
!
!
IF(ne==30800) THEN
  If(iv==0 .and. lmin/=lmax .and. lmax/=2) Then 
  Write(99,*) 'ERROR in  TASK_3:     In order to compute the  INERTIA change'
  Write(99,*) 'and derivatives, use L_min = L_max = 2 in KW HArmonic_Degrees'    
  Write(99,*) '******* JOB ABORTED ****************************************'; Stop
  Endif 
  If(iv==1 .and. lmin/=lmax .and. lmax/=2) Then 
  Write(*, *) 'ERROR in  TASK_3:     In order to compute the  INERTIA change'
  Write(*, *) 'and derivatives, use L_min = L_max = 2 in KW HArmonic_Degrees'    
  Write(*, *) '******* JOB ABORTED ****************************************' 
  Write(99,*) 'ERROR in  TASK_3:     In order to compute the  INERTIA change'
  Write(99,*) 'and derivatives, use L_min = L_max = 2 in KW HArmonic_Degrees'    
  Write(99,*) '******* JOB ABORTED ****************************************'; Stop
  Endif 
ENDIF
!
!
IF(ne==30810) THEN
IF(n_este == 0 .AND. n_inte == 0) THEN
WRITE(99,*) 'WARNING from SBR TASK_3: The kws Local_Study & Global_Study'
WRITE(99,*) 'are both inactive.  Only  the  ldcs  are available after execution'
ENDIF
IF(iv==1 .and. n_este == 0 .AND. n_inte == 0) THEN
WRITE(* ,*) 'WARNING from SBR TASK_3: The kws Local_Study & Global_Study'
WRITE(* ,*) 'are both inactive.  Only  the  ldcs  are available after execution'
ENDIF
ENDIF
!
! - - - - - - - - - - 
  end subroutine m_3
! - - - - - - - - - - 
!
!
!
!
!
!
!
  subroutine m_2 (ne) 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! # A collection of Warning & Error messages for Task#2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
   USE for_task2
   implicit NONE
!
   INTEGER ne
!
!
IF(ne==20000) then
       IF(IV/=0.and.IV/=1) THEN
       WRITE(99,*) 'ERROR in sbr task_2.dat:     The VERBOSE switch '
       WRITE(99,*) 'can only assume values 0 and 1. ** JOB ABORTED* '; STOP
       endif
endif
!
!
IF(ne==20010) then
IF(IV==0) then
If(lmin<llmin .or. lmax>llmax .or. lmax<lmin .or. lmin<0 .or. lmax<0) Then
Write(99,*)'ERROR IN SBR. task_1.dat: One of the following forbidden conditions'
Write(99,*)'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
Write(99,*)'Check the input file task_2.dat ****** JOB ABORTED ****************';Stop
Endif
Endif
IF(IV==1) then
If(lmin<llmin .or. lmax>llmax .or. lmax<lmin .or. lmin<0 .or. lmax<0) Then
Write(99,*)'ERROR IN SBR. task_1.dat: One of the following forbidden conditions'
Write(99,*)'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
Write(99,*)'Check the input file task_2.dat ****** JOB ABORTED ****************'
Write(*,*) 'ERROR IN SBR. task_1.dat: One of the following forbidden conditions'
Write(*,*) 'lmin < llmin, lmax > llmax, lmax < lmin, lmin < 0, lmax < 0 is met '
Write(*,*) 'Check the input file task_2.dat ****** JOB ABORTED ****************';Stop
Endif
Endif
endif
!
!
IF(ne==20020) then
    IF(iarmo==0.and.IV==0) then
                 WRITE(99,*) &
                 'ERROR IN SBR. TASK_2:      The KW Harmonic_Degrees '
                 WRITE(99,*) &
                 'must be activated before Make_Model can be executed';STOP
    endif
    IF(iarmo==0.and.IV==1) then
                 WRITE(99,*) &
                 'ERROR IN SBR. TASK_2:      The KW Harmonic_Degrees '
                 WRITE(99,*) &
                 'must be activated before Make_Model can be executed'
                 WRITE(*,*) &
                 'ERROR IN SBR. TASK_2:      The KW Harmonic_Degrees '
                 WRITE(*,*) &
                 'must be activated before Make_Model can be executed';STOP
    endif
endif
!
!
IF(ne==20030) then
IF(ideja==1.and.IV==0) then
       WRITE(99,*)'ERROR IN SBR. TASK_2: The Kws Make_Model and External_Model'
       WRITE(99,*)'are mutually exclusive. If Make_Model is active,  External_'
       WRITE(99,*)'model must be not. **** JOB ABORTED ***********************';STOP
endif
IF(ideja==1.and.IV==1) then
       WRITE(99,*)'ERROR IN SBR. TASK_2: The Kws Make_Model and External_Model'
       WRITE(99,*)'are mutually exclusive. If Make_Model is active,  External_'
       WRITE(99,*)'model must be not. **** JOB ABORTED ***********************'
       WRITE(*,*) 'ERROR IN SBR. TASK_2: The Kws Make_Model and External_Model'
       WRITE(*,*) 'are mutually exclusive. If Make_Model is active,  External_'
       WRITE(*,*) 'model must be not. **** JOB ABORTED ***********************';STOP
endif
endif
!
!
IF(ne==20040) then
   Write(99,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '
   Write(99,*) ' The model is the one by   GIUNCHI & SPADA,  [2000] '
   Write(99,*) ' NV=1 and LT=100 km, see  GRL table (1) for details '
   Write(99,*) ' The model is Gravitating but not self -Gravitating '
   Write(99,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '
   IF(IV==1) THEN
   Write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '
   Write(*,*) ' The model is the one by   GIUNCHI & SPADA,  [2000] '
   Write(*,*) ' NV=1 and LT=100 km, see  GRL table (1) for details '
   Write(*,*) ' The model is Gravitating but not self -Gravitating '
   Write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++ '
   ENDIF
endif
!
!
IF(ne==20050) then
    IF(inone==0.and.IV==0) then
           WRITE(99,*)'ERROR IN SBR. TASK_2:   One of the Kws Make_Model and'
           WRITE(99,*)'External_Model MUST be activated before Load_Geometry'
           WRITE(99,*)'can be executed. **** JOB ABORTED *******************';STOP
    endif
    IF(inone==0.and.IV==1) then
           WRITE(*,*) 'ERROR IN SBR. TASK_2:   One of the Kws Make_Model and'
           WRITE(*,*) 'External_Model MUST be activated before Load_Geometry'
           WRITE(*,*) 'can be executed. **** JOB ABORTED *******************'
           WRITE(99,*)'ERROR IN SBR. TASK_2:   One of the Kws Make_Model and'
           WRITE(99,*)'External_Model MUST be activated before Load_Geometry'
           WRITE(99,*)'can be executed. **** JOB ABORTED *******************';STOP
    endif
endif
!
!
IF(ne==20060) then
    IF(n_load >=2 .and. IV==0) then
      WRITE(99,*)'ERROR IN SBR. TASK_2:  The KW Load_Geometry appears'
      WRITE(99,*) n_load, 'times in file task_2.dat ** JOB ABORTED ***********';STOP
    endif
    IF(n_load >=2 .and. IV==1) then
      WRITE(* ,*)'ERROR IN SBR. TASK_2:  The KW Load_Geometry appears'
      WRITE(* ,*) n_load, 'times in file task_2.dat ** JOB ABORTED ***********'
      WRITE(99,*)'ERROR IN SBR. TASK_2:  The KW Load_Geometry appears'
      WRITE(99,*) n_load, 'times in file task_2.dat ** JOB ABORTED ***********';STOP
    endif
endif
!
!
IF(ne==20070) then!
   IF (CL /= 10 .and. &
       CL /= 11 .and. &
       CL /= 20 .and. &
       CL /= 21 .and. &
       CL /= 30 .and. &
       CL /= 40 .and. &
       CL /= 50 .and. IV==0)       THEN
   WRITE(99,*) 'ERROR IN SBR. TASK_2: A wrong code for the load has  been '
   WRITE(99,*) 'given as parameter to Load_Geometry. The code must be one '
   WRITE(99,*) 'of the following: 10, 11, 20, 21, 30, 40, 50. JOB ABORTED ';STOP
   ENDIF
   IF (CL /= 10 .and. &
       CL /= 11 .and. &
       CL /= 20 .and. &
       CL /= 21 .and. &
       CL /= 30 .and. &
       CL /= 40 .and. &
       CL /= 50 .and. IV==1)       THEN
   WRITE(99,*) 'ERROR IN SBR. TASK_2: A wrong code for the load has  been '
   WRITE(99,*) 'given as parameter to Load_Geometry. The code must be one '
   WRITE(99,*) 'of the following: 10, 11, 20, 21, 30, 40, 50. JOB ABORTED '
   WRITE(*,*)  'ERROR IN SBR. TASK_2:  A wrong code for the load has been '
   WRITE(*,*)  'given as parameter to Load_Geometry. The code must be one '
   WRITE(*,*)  'of the following: 10, 11, 20, 21, 30, 40, 50. JOB ABORTED ';STOP
   ENDIF
endif
!
!
IF(ne==20080) then
 IF ( (IOC < 0 .or. IOC > 1) .and. IV==0) THEN
 WRITE(99,*)'ERROR IN SBR. TASK_2: The option IOC = ', IOC, ' Is not allowed !'
 WRITE(99,*)'Choose IOC = 0 ---> to set NO compensation on real oceans  '
 WRITE(99,*)'Choose IOC = 1 ---> to set    compensation on real oceans  '
 WRITE(99,*)'******* JOB ABORTED *************************************  ';STOP
 ENDIF
 IF ( (IOC < 0 .or. IOC > 1) .and. IV==1) THEN
 WRITE(99,*)'ERROR IN SBR. TASK_2: The option IOC = ', IOC, ' Is not allowed !'
 WRITE(99,*)'Choose IOC = 0 ---> to set NO compensation on real oceans  '
 WRITE(99,*)'Choose IOC = 1 ---> to set    compensation on real oceans  '
 WRITE(99,*)'******* JOB ABORTED *************************************  '
 WRITE(*,*) 'ERROR IN SBR. TASK_2: The option IOC = ', IOC, ' Is not allowed !'
 WRITE(*,*) 'Choose IOC = 0 ---> to set NO compensation on real oceans  '
 WRITE(*,*) 'Choose IOC = 1 ---> to set    compensation on real oceans  '
 WRITE(*,*) '******* JOB ABORTED *************************************  ';STOP
 ENDIF
endif
!
!
IF(ne==20090) then
 IF ( (CL==11 .or. CL==21) .and. IOC==1 .and. IV==0) THEN
 WRITE(99,*)'ERROR  IN SBR. TASK_2: If the compensation on a complementary '
 WRITE(99,*)'disc load is chosen    (CL = 11 or 21),   there cannot '
 WRITE(99,*)'be compensation on a realistic ocean (IOC=1).  Try with IOC=1 '
 WRITE(99,*)'******* JOB ABORTED ***************************************** ';STOP
 ENDIF
 IF ( (CL==11 .or. CL==21) .and. IOC==1 .and. IV==1) THEN
 WRITE(99,*)'ERROR  IN SBR. TASK_2: If the compensation on a complementary '
 WRITE(99,*)'disc load is chosen    (CL = 11 or 21),   there cannot '
 WRITE(99,*)'be compensation on a realistic ocean (IOC=1).  Try with IOC=1 '
 WRITE(99,*)'******* JOB ABORTED ***************************************** '
 WRITE(*,*) 'ERROR  IN SBR. TASK_2: If the compensation on a complementary '
 WRITE(*,*) 'disc load is chosen   (CL = 11 or 21),   there  cannot '
 WRITE(*,*) 'be compensation on a realistic ocean (IOC=1).  Try with IOC=1 '
 WRITE(*,*) '******* JOB ABORTED ***************************************** ';STOP
 ENDIF
endif
!
!
IF(ne==20100) then
 IF      (lmin/=lmax .and. IV ==0) THEN
  Write  (99,*)  'ERROR from Sub. Task_2: for an Harmonic load,'
  WRITE  (99,*)  'lmin must be = lmax.    *** JOB ABORTED **** ';STOP
 ENDIF
 IF      (lmin/=lmax .and. IV ==1) THEN
  Write  (99,*)  'ERROR from Sub. Task_2: for an Harmonic load,'
  WRITE  (99,*)  'lmin must be = lmax.    *** JOB ABORTED **** '
  Write  (*,*)   'ERROR from Sub. Task_2: for an Harmonic load,'
  WRITE  (*,*)   'lmin must be = lmax.    *** JOB ABORTED **** ';STOP
 ENDIF
endif
!
!
IF(ne==20110) then
 IF (IOC /= 0 .and. IV==0) THEN
 WRITE(99,*)'ERROR IN SBR. TASK_2: The option IOC = ', IOC, ' Is not allowed ! '
 WRITE(99,*)'Only  the choice  IOC  =  0 (NO  compensation on  real oceans) is '
 WRITE(99,*)'possible for an Harmonic Load. ******* JOB ABORTED ************** ';stop
 END IF
 IF (IOC /= 0 .and. IV==1) THEN
 WRITE(99,*)'ERROR IN SBR. TASK_2: The option IOC = ', IOC, ' Is not allowed ! '
 WRITE(99,*)'Only  the choice  IOC  =  0 (NO  compensation on  real oceans) is '
 WRITE(99,*)'possible for an Harmonic Load. ******* JOB ABORTED ************** '
 WRITE(*,*) 'ERROR IN SBR. TASK_2: The option IOC = ', IOC, ' Is not allowed ! '
 WRITE(*,*) 'Only  the choice  IOC  =  0 (NO  compensation on  real oceans) is '
 WRITE(*,*) 'possible for an Harmonic Load. ******* JOB ABORTED ************** ';stop
 END IF
endif
!
!
IF(ne==20120) then
    IF(iload==0.and.IV==0) then
           WRITE(99,*) 'ERROR IN SBR.   TASK_2: The Kw Load_Geometry'
           WRITE(99,*) 'MUST be activated before Load_History can be'
           WRITE(99,*) 'executed. **********  JOB ABORTED ********* ';STOP
    endif
    IF(iload==0.and.IV==1) then
           WRITE(99,*) 'ERROR IN SBR.   TASK_2: The Kw Load_Geometry'
           WRITE(99,*) 'MUST be activated before Load_History can be'
           WRITE(99,*) 'executed. **********  JOB ABORTED ********* '
           WRITE(* ,*) 'ERROR IN SBR.   TASK_2: The Kw Load_Geometry'
           WRITE(* ,*) 'MUST be activated before Load_History can be'
           WRITE(* ,*) 'executed. **********  JOB ABORTED ********* ';STOP
    endif
endif
!
!
IF(ne==20130) then
    IF(n_hist>=2 .and. IV==0) then
           WRITE(99,*)'ERROR IN SBR. TASK_2:  The KW Load_History  appears' 
           WRITE(99,*) n_hist, 'times. ** JOB ABORTED *************';STOP
    endif
    IF(n_hist>=2 .and. IV==1) then
           WRITE(* ,*)'ERROR IN SBR. TASK_2:  The KW Load_History  appears' 
           WRITE(* ,*) n_hist, 'times. ** JOB ABORTED *************'
           WRITE(99,*)'ERROR IN SBR. TASK_2:  The KW Load_History  appears' 
           WRITE(99,*) n_hist, 'times. ** JOB ABORTED *************';STOP
    endif
endif
!
!
IF(ne==20140) then
  If ((lh<0.or.lh>8) .and. IV==0) then
  Write(99,*) 'ERROR IN SBR. TASK_2: Only time-histories with label  '
  Write(99,*) '0,  1,  2,  ... 7, and 8 are available at the moment  '
  Write(99,*) ' ***** JOB ABORTED ********************************** ';STOP
  Endif
  If ((lh<0.or.lh>8) .and. IV==1) then
  Write(99,*) 'ERROR IN SBR. TASK_2: Only time-histories with label  '
  Write(99,*) '0,  1,  2,  ... 7, and 8 are available at the moment  '
  Write(99,*) ' ***** JOB ABORTED ********************************** '
  Write(*,*)  'ERROR IN SBR. TASK_2: Only time-histories with label  '
  Write(*,*)  '0,  1,  2,  ... 7, and 8 are available at the moment  '
  Write(*,*)  ' ***** JOB ABORTED ********************************** ';STOP
  Endif
endif
!
!
IF(ne==20150) then
        IF(CL <=21 .AND. JUNK < 0.d0 .and. IV==0)  THEN
        WRITE(99,*) 'ERROR IN SBR. TASK_2: The thickness of the load'
        WRITE(99,*) 'must be >=0. ****** JOB ABORTED ***************'; STOP
        ENDIF
        IF(CL <=21 .AND. JUNK < 0.d0 .and. IV==1)  THEN
        WRITE(99,*) 'ERROR IN SBR. TASK_2: The thickness of the load'
        WRITE(99,*) 'must be >=0. ****** JOB ABORTED ***************'
        WRITE(*,*)  'ERROR IN SBR. TASK_2: The thickness of the load'
        WRITE(*,*)  'must be >=0. ****** JOB ABORTED ***************'; STOP
        ENDIF
        IF(CL ==50 .AND. JUNK < 0.d0 .and. IV==0)  THEN
        WRITE(99,*) 'ERROR IN SBR. TASK_2: The thickness of the load'
        WRITE(99,*) 'must be >=0. ****** JOB ABORTED **********'; STOP
        ENDIF
        IF(CL ==50 .AND. JUNK < 0.d0 .and. IV==1)  THEN
        WRITE(99,*) 'ERROR IN SBR. TASK_2: The thickness of the load'
        WRITE(99,*) 'must be >=0. ****** JOB ABORTED **********'
        WRITE(*,*)  'ERROR IN SBR. TASK_2: The thickness of the load'
        WRITE(*,*)  'must be >=0. ****** JOB ABORTED **********'; STOP
        ENDIF
        IF(CL ==30 .AND. JUNK < 0.d0 .and. IV==0)  THEN
        WRITE(99,*) 'ERROR IN SBR. TASK_2: The mass of the load'
        WRITE(99,*) 'must be >=0. ****** JOB ABORTED **********'; STOP
        ENDIF
        IF(CL ==30 .AND. JUNK < 0.d0 .and. IV==1)  THEN
        WRITE(99,*) 'ERROR IN SBR. TASK_2: The mass of the load'
        WRITE(99,*) 'must be >=0. ****** JOB ABORTED **********'
        WRITE(*,*)  'ERROR IN SBR. TASK_2: The mass of the load'
        WRITE(*,*)  'must be >=0. ****** JOB ABORTED **********'; STOP
        ENDIF
        IF(CL ==40 .AND. JUNK < 0.d0 .and. IV==0)  THEN
        WRITE(99,*) 'ERROR IN SBR. TASK_2: The load pressure   '
        WRITE(99,*) 'must be >=0. ****** JOB ABORTED *******   '; STOP
        ENDIF
        IF(CL ==40 .AND. JUNK < 0.d0 .and. IV==1)  THEN
        WRITE(99,*) 'ERROR IN SBR. TASK_2: The load pressure   '
        WRITE(99,*) 'must be >=0. ****** JOB ABORTED *******   '
        WRITE(*,*)  'ERROR IN SBR. TASK_2: The load pressure   '
        WRITE(*,*)  'must be >=0. ****** JOB ABORTED *******   '; STOP
        ENDIF
endif
!
!
IF(ne==20160) THEN
  If (TAU_LU<=0.d0 .and. IV==0) Then
   write(99,*) 'ERROR IN SBR. TASK_2: The lapse of time between loading '
   write(99,*) 'and ul- loading must be > 0 kyrs **** JOB ABORTED ***** '; STOP
  Endif
  If (TAU_LU<=0.d0 .and. IV==1) Then
   write(99,*) 'ERROR IN SBR. TASK_2: The lapse of time between loading '
   write(99,*) 'and ul- loading must be > 0 kyrs **** JOB ABORTED ***** '
   write(* ,*) 'ERROR IN SBR. TASK_2: The lapse of time between loading '
   write(* ,*) 'and ul- loading must be > 0 kyrs **** JOB ABORTED ***** '; STOP
  Endif
ENDIF
!
!
IF(ne==20170) THEN
  If (TAU_SD<=0.d0 .and. IV==0) Then
   write(99,*) 'ERROR IN SBR. TASK_2:   The lenght of the melting phase '
   write(99,*) 'must be positive in time- history #3. ** JOB ABORTED *** '
  Endif
  If (TAU_SD<=0.d0 .and. IV==1) Then
   write(99,*) 'ERROR IN SBR. TASK_2:   The lenght of the melting phase '
   write(99,*) 'must be positive in time- history #3. ** JOB ABORTED *** '
   write(* ,*) 'ERROR IN SBR. TASK_2:   The lenght of the melting phase '
   write(* ,*) 'must be positive in time- history #3. ** JOB ABORTED *** '; STOP
  Endif
ENDIF
!
!
IF(ne==20180) THEN
if((NR<0.or.NR>5) .and. IV==0) then
write(99,*) 'ERROR in sbr TASK_2:    The maximum allowed number of'
write(99,*) 'loading and unloading phases is NR=5. JOB ABORTED ** '; STOP
endif
if((NR<0.or.NR>5) .and. IV==1) then
write(99,*) 'ERROR in sbr TASK_2:    The maximum allowed number of'
write(99,*) 'loading and unloading phases is NR=5. JOB ABORTED ** '
write(* ,*) 'ERROR in sbr TASK_2:    The maximum allowed number of'
write(* ,*) 'loading and unloading phases is NR=5. JOB ABORTED ** '; STOP
endif
ENDIF
!
!
!
IF(ne==20190) THEN
if(TAUC<=0.0 .and. IV==0) then
write(99,*) 'ERROR in sbr TASK_2: The length of the loading phase'
write(99,*) 'must be positive.  ***************** JOB ABORTED ** '; STOP
endif
if(TAUC<=0.0 .and. IV==1) then
write(99,*) 'ERROR in sbr TASK_2: The length of the loading phase'
write(99,*) 'must be positive.  ***************** JOB ABORTED ** '
write(* ,*) 'ERROR in sbr TASK_2: The length of the loading phase'
write(* ,*) 'must be positive.  ***************** JOB ABORTED ** '; STOP
endif
ENDIF
!
!
!
IF(ne==20200) THEN
if(DINC<=0.0 .and. IV==0) then
write(99,*) 'ERROR in sbr TASK_2: The length of the un -loading phase'
write(99,*) 'must be positive.  *****************     JOB ABORTED ** '; STOP
endif
if(DINC<=0.0 .and. IV==1) then
write(99,*) 'ERROR in sbr TASK_2: The length of the un- loading phase'
write(99,*) 'must be positive.  *****************     JOB ABORTED ** '
write(* ,*) 'ERROR in sbr TASK_2: The length of the un -loading phase'
write(* ,*) 'must be positive.  *****************     JOB ABORTED ** '; STOP
endif
ENDIF
!
!
IF(ne==20210) THEN
if(period <= 0.0 .and. IV==0) then
write(99,*) 'ERROR in sbr TASK_2:    The period of the sinusoidal'
write(99,*) 'time-history must be positive. ***** JOB ABORTED ** '; STOP
endif
if(period <= 0.0 .and. IV==0) then
write(99,*) 'ERROR in sbr TASK_2:    The period of the sinusoidal'
write(99,*) 'time-history must be positive. ***** JOB ABORTED ** '
write(* ,*) 'ERROR in sbr TASK_2:    The period of the sinusoidal'
write(* ,*) 'time-history must be positive. ***** JOB ABORTED ** '; STOP
endif
ENDIF
!
!
IF(ne==20220) THEN
         IF(iv==0 .AND. ijunk /= 6) THEN
         WRITE(99,*) 'ERROR IN SBR. TASK_2:   The file timeh_6.dat should'
         WRITE(99,*) 'have the correct time history stamp (6) on line one'
         WRITE(99,*) '******** JOB ABORTED ******************************';stop
         endif
         IF(iv==1 .AND. ijunk /= 6) THEN
         WRITE(* ,*) 'ERROR IN SBR. TASK_2:   The file timeh_6.dat should'
         WRITE(* ,*) 'have the correct time history stamp (6) on line one'
         WRITE(* ,*) '******** JOB ABORTED ******************************'
         WRITE(99,*) 'ERROR IN SBR. TASK_2:   The file timeh_6.dat should'
         WRITE(99,*) 'have the correct time history stamp (6) on line one'
         WRITE(99,*) '******** JOB ABORTED ******************************';stop
         endif
ENDIF
!
!
IF(ne==20230) THEN
         IF(iv==0 .AND. t_hh(0) /= 0.0) THEN
     WRITE(99,*) 'ERROR IN SBR. TASK_2:    The time t_0 must be zero'
     WRITE(99,*) 'in time-history #6. ***** JOB ABORTED ************'
     endif
         IF(iv==1 .AND. t_hh(0) /= 0.0) THEN
     WRITE(99,*) 'ERROR IN SBR. TASK_2:    The time t_0 must be zero'
     WRITE(99,*) 'in time-history #6. ***** JOB ABORTED ************'
     WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The time t_0 must be zero'
     WRITE(* ,*) 'in time-history #6. ***** JOB ABORTED ************';STOP
         ENDIF
ENDIF
!
!
IF(ne==20240) THEN
         IF(iv==0 .AND. t_hh(i) <= 0.0) THEN
     WRITE(99,*) 'ERROR IN SBR. TASK_2:    The times which mark the boundaries'
     WRITE(99,*) 'between adjacent linear phases must not be negative in time-'
     WRITE(99,*) 'history #6. ***** JOB ABORTED ******************************';stop
     endif
         IF(iv==1 .AND. t_hh(i) <= 0.0) THEN
     WRITE(99,*) 'ERROR IN SBR. TASK_2:    The times which mark the boundaries'
     WRITE(99,*) 'between adjacent linear phases must not be negative in time-'
     WRITE(99,*) 'history #6. ***** JOB ABORTED ******************************'
     WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The times which mark the boundaries'
     WRITE(* ,*) 'between adjacent linear phases must not be negative in time-'
     WRITE(* ,*) 'history #6. ***** JOB ABORTED ******************************';stop
         ENDIF
ENDIF
!
!
IF(ne==20250) THEN
         IF(iv==0 .AND. i>=1 .and. t_hh(i) <= t_hh(i-1)) THEN
  WRITE(99,*) 'ERROR IN SBR. TASK_2:        The times which mark the boundaries'
  WRITE(99,*) 'between adjacent linear phases must be ordered from the smallest'
  WRITE(99,*) 'to the largest in time-history #6. **** JOB ABORTED ************';stop
     endif
         IF(iv==1 .AND. i>=1 .and. t_hh(i) <= t_hh(i-1)) THEN
  WRITE(99,*) 'ERROR IN SBR. TASK_2:        The times which mark the boundaries'
  WRITE(99,*) 'between adjacent linear phases must be ordered from the smallest'
  WRITE(99,*) 'to the largest in time-history #6. **** JOB ABORTED ************'
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:        The times which mark the boundaries'
  WRITE(* ,*) 'between adjacent linear phases must be ordered from the smallest'
  WRITE(* ,*) 'to the largest in time-history #6. **** JOB ABORTED ************';stop
         ENDIF
ENDIF
!
!
IF(ne==20260) THEN
  IF(iv==0 .AND. a_hh(i) < 0.0 .AND. (CL <=21 .OR. CL == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#6. *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==1 .AND. a_hh(i) < 0.0 .AND. (CL <=21 .OR. CL == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#6. *** JOB ABORTED *************  '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:  The load thickness must be positive  '
  WRITE(* ,*) 'or zero at any time in th#6. *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==0 .AND. a_hh(i) < 0.0 .AND. CL == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#6. **** JOB ABORTED **************  ';STOP
         ENDIF
  IF(iv==1 .AND. a_hh(i) < 0.0 .AND. CL == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#6. **** JOB ABORTED **************  '
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#6. **** JOB ABORTED **************  ';STOP
         ENDIF
  IF(iv==0 .AND. a_hh(i) < 0.0 .AND. CL == 40) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load per unit surface  '
  WRITE(99,*) 'at the pole of the load must be positive or zero at any time  '
  WRITE(99,*) 'in time-history #6. ***** JOB ABORTED **********************  ';STOP
         ENDIF
  IF(iv==1 .AND. a_hh(i) < 0.0 .AND. CL == 40) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load per unit surface  '
  WRITE(99,*) 'at the pole of the load must be positive or zero at any time  '
  WRITE(99,*) 'in time-history #6. ***** JOB ABORTED **********************  '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:  The mass of the load per unit surface  '
  WRITE(* ,*) 'at the pole of the load must be positive or zero at any time  '
  WRITE(* ,*) 'in time-history #6. ***** JOB ABORTED **********************  ';STOP
         ENDIF
endif
!
!
IF(ne==20270) THEN
  IF(iv==0 .AND. nss == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  At least two points must be given in '
  WRITE(99,*) 'file timeh_6.dat in order to specify the time-history th #6 '
  WRITE(99,*) '************************* JOB ABORTED ********************* ';STOP
         ENDIF
  IF(iv==1 .AND. nss == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  At least two points must be given in '
  WRITE(99,*) 'file timeh_6.dat in order to specify the time-history th #6 '
  WRITE(99,*) '************************* JOB ABORTED ********************* '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:  At least two points must be given in '
  WRITE(* ,*) 'file timeh_6.dat in order to specify the time-history th #6 '
  WRITE(* ,*) '************************* JOB ABORTED ********************* ';STOP
         ENDIF
endif
!
!
IF(ne==20280) THEN
  IF(iv==0 .AND. nss > nss_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The number of linear segments in time-'
  WRITE(99,*) 'history #6 must not exceed ', nss_max
  WRITE(99,*) '************************* JOB ABORTED ********************** ';STOP
         ENDIF
  IF(iv==1 .AND. nss > nss_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The number of linear segments in time-'
  WRITE(99,*) 'history #6 must not exceed ', nss_max
  WRITE(99,*) '************************* JOB ABORTED ********************** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:  The number of linear segments in time-'
  WRITE(* ,*) 'history #6 must not exceed ', nss_max
  WRITE(* ,*) '************************* JOB ABORTED ********************** ';STOP
         ENDIF
endif
!
!
IF(ne==20290) THEN
  IF(iv==0 .AND. junk <= 0.0 .AND. (CL <=21 .OR. CL == 50))    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The maximum load thickness must be '
  WRITE(99,*) 'positive or zero in th#6 **************** JOB ABORTED *** ';STOP
         ENDIF
  IF(iv==1 .AND. junk <= 0.0 .AND. (CL <=21 .OR. CL == 50))    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The maximum load thickness must be '
  WRITE(99,*) 'positive or zero in th#6 **************** JOB ABORTED ***** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The maximum load thickness must be '
  WRITE(* ,*) 'positive or zero in th#6 **************** JOB ABORTED ***** ';STOP
         ENDIF
  IF(iv==0 .AND. junk <= 0.0 .AND. CL == 50)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The maximum load mass must be '
  WRITE(99,*) 'positive or zero in th#6 ********* JOB ABORTED ***** ';STOP
         ENDIF
  IF(iv==1 .AND. junk <= 0.0 .AND. CL == 50)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The maximum load mass must be '
  WRITE(99,*) 'positive or zero in th#6 ********* JOB ABORTED ******* '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The maximum load mass must be '
  WRITE(* ,*) 'positive or zero in th#6 ********* JOB ABORTED ******* ';STOP
         ENDIF
  IF(iv==0 .AND. junk <= 0.0 .AND. CL == 40)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The maximum load mass/surface must be '
  WRITE(99,*) 'positive or zero in th#6 ******************* JOB ABORTED *** ';STOP
         ENDIF
  IF(iv==1 .AND. junk <= 0.0 .AND. CL == 40)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The maximum load mass/surface must be '
  WRITE(99,*) 'positive or zero in th#6 ******************* JOB ABORTED ***** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The maximum load mass/surface must be '
  WRITE(* ,*) 'positive or zero in th#6 ******************* JOB ABORTED ***** ';STOP
         ENDIF
endif
!
!
IF(ne==20300) THEN
         IF(iv==0 .AND. ijunk /= 7) THEN
         WRITE(99,*) 'ERROR IN SBR. TASK_2:   The file timeh_7.dat should'
         WRITE(99,*) 'have the correct time history label (7) on line one'
         WRITE(99,*) '******** JOB ABORTED ******************************';stop
         endif
         IF(iv==1 .AND. ijunk /= 7) THEN
         WRITE(* ,*) 'ERROR IN SBR. TASK_2:   The file timeh_7.dat should'
         WRITE(* ,*) 'have the correct time history label (7) on line one'
         WRITE(* ,*) '******** JOB ABORTED ******************************'
         WRITE(99,*) 'ERROR IN SBR. TASK_2:   The file timeh_7.dat should'
         WRITE(99,*) 'have the correct time history label (7) on line one'
         WRITE(99,*) '******** JOB ABORTED ******************************';stop
         endif
endif
!
!
IF(ne==20310) THEN
  IF(iv==0 .AND. (dilta < dilta_min .or. dilta > dilta_max) ) THEN
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The time step in Th#7 must be in the'
  WRITE(99,*) 'range', dilta_min,'-',dilta_max,' kyrs. *** JOB ABORTED ***** ';STOP
  endif
  IF(iv==1.AND. (dilta < dilta_min .or. dilta > dilta_max) ) THEN
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The time step in Th#7 must be in the'
  WRITE(99,*) 'range', dilta_min,'-',dilta_max,' kyrs. *** JOB ABORTED ***** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The time step in Th#7 must be in the'
  WRITE(* ,*) 'range', dilta_min,'-',dilta_max,' kyrs. *** JOB ABORTED ***** ';STOP
  ENDIF
endif
!
!
IF(ne==20320) THEN
  IF(iv==0 .AND. L /= I) THEN
   WRITE(99,*) 'ERROR IN SBR. TASK_2: I expected step#', L, ',I have found step #', I
   WRITE(99,*) 'Check the input file timeh_7.dat ****** JOB ABORTED *********** ';STOP
  ENDIF
  IF(iv==1 .AND. L /= I) THEN
  WRITE(99,*) 'ERROR IN SBR. TASK_2: I expected step#', L, ',I have found step #', I
  WRITE(99,*) 'Check the input file timeh_7.dat ****** JOB ABORTED *********** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2: I expected step#', L, ',I have found step #', I
  WRITE(* ,*) 'Check the input file timeh_7.dat ****** JOB ABORTED *********** ';STOP
  ENDIF
endif
!
!
IF(ne==20330) THEN
  IF(iv==0 .AND. a_h(i) < 0.0 .AND. (CL <=21 .OR. CL == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==1 .AND. a_h(i) < 0.0 .AND. (CL <=21 .OR. CL == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. *** JOB ABORTED *************  '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:  The load thickness must be positive  '
  WRITE(* ,*) 'or zero at any time in th#7. *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==0 .AND. a_h(i) < 0.0 .AND. CL == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  ';STOP
         ENDIF
  IF(iv==1 .AND. a_h(i) < 0.0 .AND. CL == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  '
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#7. **** JOB ABORTED **************  ';STOP
         ENDIF
  IF(iv==0 .AND. a_h(i) < 0.0 .AND. CL == 40) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load per unit surface  '
  WRITE(99,*) 'at the pole of the load must be positive or zero at any time  '
  WRITE(99,*) 'in time-history #7. ***** JOB ABORTED **********************  ';STOP
         ENDIF
  IF(iv==1 .AND. a_h(i) < 0.0 .AND. CL == 40) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load per unit surface  '
  WRITE(99,*) 'at the pole of the load must be positive or zero at any time  '
  WRITE(99,*) 'in time-history #7. ***** JOB ABORTED **********************  '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:  The mass of the load per unit surface  '
  WRITE(* ,*) 'at the pole of the load must be positive or zero at any time  '
  WRITE(* ,*) 'in time-history #7. ***** JOB ABORTED **********************  ';STOP
         ENDIF
endif
!
!
IF(ne==20340) THEN
  IF(iv==0 .AND. ns == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    At least one step must be given in '
  WRITE(99,*) 'file timeh_7.dat in order to specify the time-history th #7 '
  WRITE(99,*) '************************* JOB ABORTED ********************* ';STOP
         ENDIF
  IF(iv==1 .AND. ns == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    At least one step must be given in '
  WRITE(99,*) 'file timeh_7.dat in order to specify the time-history th #7 '
  WRITE(99,*) '************************* JOB ABORTED ********************* '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    At least one step must be given in '
  WRITE(* ,*) 'file timeh_7.dat in order to specify the time-history th #7 '
  WRITE(* ,*) '************************* JOB ABORTED ********************* ';STOP
         ENDIF
endif
!
!
IF(ne==20350) THEN
  IF(iv==0 .AND. ns > ns_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The number of steps in time-'
  WRITE(99,*) 'history #7 must not exceed ', ns_max
  WRITE(99,*) '************************* JOB ABORTED ************ ';STOP
         ENDIF
  IF(iv==1 .AND. ns > ns_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The number of steps in time-'
  WRITE(99,*) 'history #7 must not exceed ', ns_max
  WRITE(99,*) '************************* JOB ABORTED ************ '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:  The number of steps in time-'
  WRITE(* ,*) 'history #7 must not exceed ', ns_max
  WRITE(* ,*) '************************* JOB ABORTED ************ ';STOP
         ENDIF
endif
!
!
IF(ne==20360) THEN
  IF(iv==0 .AND. junk <= 0.0 .AND. (CL <=21 .OR. CL == 50))    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The maximum load thickness must be '
  WRITE(99,*) 'positive or zero in th#7 **************** JOB ABORTED *** ';STOP
         ENDIF
  IF(iv==1 .AND. junk <= 0.0 .AND. (CL <=21 .OR. CL == 50))    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The maximum load thickness must be '
  WRITE(99,*) 'positive or zero in th#7 **************** JOB ABORTED ***** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The maximum load thickness must be '
  WRITE(* ,*) 'positive or zero in th#7 **************** JOB ABORTED ***** ';STOP
         ENDIF
  IF(iv==0 .AND. junk <= 0.0 .AND. CL == 50)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The maximum load mass must be '
  WRITE(99,*) 'positive or zero in th#7 ********* JOB ABORTED ***** ';STOP
         ENDIF
  IF(iv==1 .AND. junk <= 0.0 .AND. CL == 50)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The maximum load mass must be '
  WRITE(99,*) 'positive or zero in th#7 ********* JOB ABORTED ******* '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The maximum load mass must be '
  WRITE(* ,*) 'positive or zero in th#7 ********* JOB ABORTED ******* ';STOP
         ENDIF
  IF(iv==0 .AND. junk <= 0.0 .AND. CL == 40)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The maximum load mass/surface must be '
  WRITE(99,*) 'positive or zero in th#7 ******************* JOB ABORTED *** ';STOP
         ENDIF
  IF(iv==1 .AND. junk <= 0.0 .AND. CL == 40)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The maximum load mass/surface must be '
  WRITE(99,*) 'positive or zero in th#7 ******************* JOB ABORTED ***** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The maximum load mass/surface must be '
  WRITE(* ,*) 'positive or zero in th#7 ******************* JOB ABORTED ***** ';STOP
         ENDIF
endif
!
!
IF(ne==20370) THEN
         IF(iv==0 .AND. ijunk /= 8) THEN
         WRITE(99,*) 'ERROR IN SBR. TASK_2:   The file timeh_8.dat should'
         WRITE(99,*) 'have the correct time history label (8) on line one'
         WRITE(99,*) '******** JOB ABORTED ******************************';stop
         endif
         IF(iv==1 .AND. ijunk /= 8) THEN
         WRITE(* ,*) 'ERROR IN SBR. TASK_2:   The file timeh_8.dat should'
         WRITE(* ,*) 'have the correct time history label (8) on line one'
         WRITE(* ,*) '******** JOB ABORTED ******************************'
         WRITE(99,*) 'ERROR IN SBR. TASK_2:   The file timeh_8.dat should'
         WRITE(99,*) 'have the correct time history label (8) on line one'
         WRITE(99,*) '******** JOB ABORTED ******************************';stop
         endif
ENDIF
!
!
IF(ne==20380) THEN
  IF(iv==0 .AND. (dilta < dilta_min .or. dilta > dilta_max) ) THEN
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The time step in Th#8 must be in the'
  WRITE(99,*) 'range', dilta_min,'-',dilta_max,' kyrs. *** JOB ABORTED ***** ';STOP
  endif
  IF(iv==1.AND. (dilta < dilta_min .or. dilta > dilta_max) ) THEN
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The time step in Th#8 must be in the'
  WRITE(99,*) 'range', dilta_min,'-',dilta_max,' kyrs. *** JOB ABORTED ***** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The time step in Th#8 must be in the'
  WRITE(* ,*) 'range', dilta_min,'-',dilta_max,' kyrs. *** JOB ABORTED ***** ';STOP
  ENDIF
endif
!
!
IF(ne==20390) THEN
IF(iv==0 .AND. (tauf < tauf_min .or. tauf > tauf_max) ) THEN
WRITE(99,*) 'ERROR IN SBR. TASK_2: The length of the loading phase in Th#8 '
WRITE(99,*) 'must be in the range', tauf_min,'-',tauf_max,'kyrs ***JOB ABORTED**';STOP
endif
IF(iv==1 .AND. (tauf < tauf_min .or. tauf > tauf_max) ) THEN
WRITE(99,*) 'ERROR IN SBR. TASK_2: The length of the loading phase in Th#8 '
WRITE(99,*) 'must be in the range', tauf_min,'-',tauf_max,'kyrs ***JOB ABORTED**'
WRITE(* ,*) 'ERROR IN SBR. TASK_2: The length of the loading phase in Th#8 '
WRITE(* ,*) 'must be in the range', tauf_min,'-',tauf_max,'kyrs ***JOB ABORTED**';STOP
ENDIF
endif
!
!
IF(ne==20400) THEN
  IF(iv==0 .AND. L /= I) THEN
   WRITE(99,*) 'ERROR IN SBR. TASK_2: I expected step#', L, ',I have found step #', I
   WRITE(99,*) 'Check the input file timeh_8.dat ****** JOB ABORTED *********** ';STOP
  ENDIF
  IF(iv==1 .AND. L /= I) THEN
  WRITE(99,*) 'ERROR IN SBR. TASK_2: I expected step#', L, ',I have found step #', I
  WRITE(99,*) 'Check the input file timeh_8.dat ****** JOB ABORTED *********** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2: I expected step#', L, ',I have found step #', I
  WRITE(* ,*) 'Check the input file timeh_8.dat ****** JOB ABORTED *********** ';STOP
  ENDIF
endif
!
!
IF(ne==20410) THEN
  IF(iv==0 .AND. a_h(i) < 0.0 .AND. (CL <=21 .OR. CL == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#8. *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==1 .AND. a_h(i) < 0.0 .AND. (CL <=21 .OR. CL == 50)) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The load thickness must be positive  '
  WRITE(99,*) 'or zero at any time in th#8. *** JOB ABORTED *************  '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:  The load thickness must be positive  '
  WRITE(* ,*) 'or zero at any time in th#8. *** JOB ABORTED *************  ';STOP
         ENDIF
  IF(iv==0 .AND. a_h(i) < 0.0 .AND. CL == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#8. **** JOB ABORTED **************  ';STOP
         ENDIF
  IF(iv==1 .AND. a_h(i) < 0.0 .AND. CL == 30) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#8. **** JOB ABORTED **************  '
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load must be positive  '
  WRITE(99,*) 'or zero at any time in th#8. **** JOB ABORTED **************  ';STOP
         ENDIF
  IF(iv==0 .AND. a_h(i) < 0.0 .AND. CL == 40) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load per unit surface  '
  WRITE(99,*) 'at the pole of the load must be positive or zero at any time  '
  WRITE(99,*) 'in time-history #8. ***** JOB ABORTED **********************  ';STOP
         ENDIF
  IF(iv==1 .AND. a_h(i) < 0.0 .AND. CL == 40) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The mass of the load per unit surface  '
  WRITE(99,*) 'at the pole of the load must be positive or zero at any time  '
  WRITE(99,*) 'in time-history #8. ***** JOB ABORTED **********************  '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:  The mass of the load per unit surface  '
  WRITE(* ,*) 'at the pole of the load must be positive or zero at any time  '
  WRITE(* ,*) 'in time-history #8. ***** JOB ABORTED **********************  ';STOP
         ENDIF
endif
!
!
IF(ne==20420) THEN
  IF(iv==0 .AND. ns == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    At least one step must be given in '
  WRITE(99,*) 'file timeh_8.dat in order to specify the time-history th #8 '
  WRITE(99,*) '************************* JOB ABORTED ********************* ';STOP
         ENDIF
  IF(iv==1 .AND. ns == 0) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    At least one step must be given in '
  WRITE(99,*) 'file timeh_8.dat in order to specify the time-history th #8 '
  WRITE(99,*) '************************* JOB ABORTED ********************* '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    At least one step must be given in '
  WRITE(* ,*) 'file timeh_8.dat in order to specify the time-history th #8 '
  WRITE(* ,*) '************************* JOB ABORTED ********************* ';STOP
         ENDIF
endif
!
!
IF(ne==20430) THEN
  IF(iv==0 .AND. ns > ns_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The number of steps in time-'
  WRITE(99,*) 'history #8 must not exceed ', ns_max
  WRITE(99,*) '************************* JOB ABORTED ************ ';STOP
         ENDIF
  IF(iv==1 .AND. ns > ns_max) then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The number of steps in time-'
  WRITE(99,*) 'history #8 must not exceed ', ns_max
  WRITE(99,*) '************************* JOB ABORTED ************ '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:  The number of steps in time-'
  WRITE(* ,*) 'history #8 must not exceed ', ns_max
  WRITE(* ,*) '************************* JOB ABORTED ************ ';STOP
         ENDIF
endif
!
!
IF(ne==20440) THEN
  IF(iv==0 .AND. junk <= 0.0 .AND. (CL <=21 .OR. CL == 50))    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The maximum load thickness must be '
  WRITE(99,*) 'positive or zero in th#8 **************** JOB ABORTED *** ';STOP
         ENDIF
  IF(iv==1 .AND. junk <= 0.0 .AND. (CL <=21 .OR. CL == 50))    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The maximum load thickness must be '
  WRITE(99,*) 'positive or zero in th#8 **************** JOB ABORTED ***** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The maximum load thickness must be '
  WRITE(* ,*) 'positive or zero in th#8 **************** JOB ABORTED ***** ';STOP
         ENDIF
  IF(iv==0 .AND. junk <= 0.0 .AND. CL == 50)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The maximum load mass must be '
  WRITE(99,*) 'positive or zero in th#8 ********* JOB ABORTED ***** ';STOP
         ENDIF
  IF(iv==1 .AND. junk <= 0.0 .AND. CL == 50)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The maximum load mass must be '
  WRITE(99,*) 'positive or zero in th#8 ********* JOB ABORTED ******* '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The maximum load mass must be '
  WRITE(* ,*) 'positive or zero in th#8 ********* JOB ABORTED ******* ';STOP
         ENDIF
  IF(iv==0 .AND. junk <= 0.0 .AND. CL == 40)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:  The maximum load mass/surface must be '
  WRITE(99,*) 'positive or zero in th#8 ******************* JOB ABORTED *** ';STOP
         ENDIF
  IF(iv==1 .AND. junk <= 0.0 .AND. CL == 40)    then
  WRITE(99,*) 'ERROR IN SBR. TASK_2:    The maximum load mass/surface must be '
  WRITE(99,*) 'positive or zero in th#8 ******************* JOB ABORTED ***** '
  WRITE(* ,*) 'ERROR IN SBR. TASK_2:    The maximum load mass/surface must be '
  WRITE(* ,*) 'positive or zero in th#8 ******************* JOB ABORTED ***** ';STOP
         ENDIF
endif
!
!
IF(ne==20450) THEN
    IF(ihist==0.and.IV==0) then
           WRITE(99,*) 'ERROR IN SBR.   TASK_2: The Kw Load_Histoty'
           WRITE(99,*) 'MUST be activated before Local_Study can be'
           WRITE(99,*) 'executed. **********  JOB ABORTED ********* ';STOP
    endif
    IF(ihist==0.and.IV==1) then
           WRITE(99,*) 'ERROR IN SBR.   TASK_2: The Kw Load_Histoty'
           WRITE(99,*) 'MUST be activated before Local_Study can be'
           WRITE(99,*) 'executed. **********  JOB ABORTED ********* ' 
           WRITE(* ,*) 'ERROR IN SBR.   TASK_2: The Kw Load_Histoty'
           WRITE(* ,*) 'MUST be activated before Local_Study can be'
           WRITE(* ,*) 'executed. **********  JOB ABORTED ********* ';STOP
    endif
endif
!
!
IF(ne==20460) THEN
   IF(n_inte >=2 .and. IV==0) then
      WRITE(99,*)'ERROR IN SBR. TASK_2:  The KW Local_Study appears' 
      WRITE(99,*) n_inte, 'times in file task_2.dat ** JOB ABORTED ***********';STOP
    endif
    IF(n_inte >=2 .and. IV==1) then
      WRITE(* ,*)'ERROR IN SBR. TASK_2:  The KW Local_Study appears' 
      WRITE(* ,*) n_inte, 'times in file task_2.dat ** JOB ABORTED ***********'
      WRITE(99,*)'ERROR IN SBR. TASK_2:  The KW Local_Study appears' 
      WRITE(99,*) n_inte, 'times in file task_2.dat ** JOB ABORTED ***********';STOP
    endif
endif
!
!
IF(ne==20470) THEN
 If ((type_inten<=0.or.type_inten>5).and.iv==0) then 
  Write(99,*) 'ERROR IN SBR. TASK_2:  the Local_Study label ' 
  Write(99,*) 'must vary in the range [1:5] *** JOB ABORTED *** '; stop 
 Endif
 If ((type_inten<=0.or.type_inten>5).and.iv==1) then 
  Write(99,*) 'ERROR IN SBR. TASK_2:  the Local_Study label ' 
  Write(99,*) 'must vary in the range [1:5] *** JOB ABORTED *** ' 
  Write(* ,*) 'ERROR IN SBR. TASK_2:  the Local_Study label ' 
  Write(* ,*) 'must vary in the range [1:5] *** JOB ABORTED *** '; stop 
 Endif
endif
!
!
IF(ne==20480) THEN
If (IR==0.and.dIR==0.and.it==0.and.dit==0  &
                            .and.IG==0.and.dIG==0.and.iv==0) then 
Write(99,*) 'ERROR IN SBR. TASK_2: At least one of the switches IR, DIR, etc. must '
Write(99,*) 'be =1 in order to perform an Local_Study. **** JOB ABORTED ********** ';stop
Endif          
If (IR==0.and.dIR==0.and.it==0.and.dit==0  &
                            .and.IG==0.and.dIG==0.and.iv==1) then 
Write(99,*) 'ERROR IN SBR. TASK_2: At least one of the switches IR, DIR, etc. must '
Write(99,*) 'be =1 in order to perform an Local_Study. **** JOB ABORTED ********** '
Write(* ,*) 'ERROR IN SBR. TASK_2: At least one of the switches IR, DIR, etc. must '
Write(* ,*) 'be =1 in order to perform an Local_Study. **** JOB ABORTED ********** ';stop
endif
endif
!
!
IF(ne==20490) THEN
 If (iv==0 .and. tmin > tmax) then 
  Write(99,*) 'ERROR IN SBR. TASK_2: TMIN must be <= TMAX in KW Local_Study ' 
  Write(99,*) '***** JOB ABORTED ********************************************';stop    
 Endif
 If (iv==1 .and. tmin > tmax) then 
  Write(99,*) 'ERROR IN SBR. TASK_2: TMIN must be <= TMAX in KW Local_Study ' 
  Write(99,*) '***** JOB ABORTED ********************************************'   
  Write(* ,*) 'ERROR IN SBR. TASK_2: TMIN must be <= TMAX in KW Local_Study ' 
  Write(* ,*) '***** JOB ABORTED ********************************************';stop   
 Endif
endif
!
!
IF(ne==20500) THEN
 If (iv==0 .and. DELTAT < 0.0) then 
  Write(99,*) 'ERROR IN SBR. TASK_2: DELTAT must be >= 0. in KW Local_Study ' 
  Write(99,*) '***** JOB ABORTED ********************************************';stop    
 Endif
 If (iv==1 .and. DELTAT < 0.0) then 
  Write(99,*) 'ERROR IN SBR. TASK_2: DELTAT must be >= 0. in KW Local_Study ' 
  Write(99,*) '***** JOB ABORTED ********************************************'
  Write(* ,*) 'ERROR IN SBR. TASK_2: DELTAT must be >= 0. in KW Local_Study ' 
  Write(* ,*) '***** JOB ABORTED ********************************************';stop   
 Endif
endif
!
!
IF(ne==20510) THEN
 IF(NOBS>NOBS_MAX .and. IV==0) THEN 
  WRITE (99,*) 'ERROR IN SBR. TASK_2: The number of observers in file',file_sparsi 
  write (99,*) 'exceeds the maximum allowed', nobs_max,' *** JOB ABORTED *********'; stop 
 ENDIF 
 IF(NOBS>NOBS_MAX .and. IV==1) THEN 
  WRITE (99,*) 'ERROR IN SBR. TASK_2: The number of observers in file',file_sparsi 
  write (99,*) 'exceeds the maximum allowed', nobs_max,' *** JOB ABORTED *********' 
  WRITE (* ,*) 'ERROR IN SBR. TASK_2: The number of observers in file',file_sparsi 
  write (* ,*) 'exceeds the maximum allowed', nobs_max,' *** JOB ABORTED *********'; stop   
 ENDIF
endif
!
!
IF(ne==20520) THEN
IF(NOBS==0 .and. IV==0) THEN 
  WRITE (99,*) 'ERROR IN SBR. TASK_2: It appears that the file',file_sparsi,'is empty'
  write (99,*) 'or that its name is not properly aligned in file task_2.dat. JOB ABORTED ****'; stop 
 ENDIF 
 IF(NOBS==0 .and. IV==1) THEN 
  WRITE (99,*) 'ERROR IN SBR. TASK_2: It appears that the file',file_sparsi,'is empty'
  write (99,*) 'or that its name is not properly aligned in file task_2.dat. JOB ABORTED ****' 
  WRITE (*, *) 'ERROR IN SBR. TASK_2: It appears that the file',file_sparsi,'is empty'
  write (*, *) 'or that its name is not properly aligned in file task_2.dat. JOB ABORTED ****'; stop 
 ENDIF
endif
!
!
IF(ne==20530) THEN
 If (iv==0 .and. tmin > tmax) then 
  Write(99,*) 'ERROR IN SBR. TASK_2: TMIN must be <= TMAX in KW Local_Study ' 
  Write(99,*) '***** JOB ABORTED ********************************************';stop    
 Endif
 If (iv==1 .and. tmin > tmax) then 
  Write(99,*) 'ERROR IN SBR. TASK_2: TMIN must be <= TMAX in KW Local_Study ' 
  Write(99,*) '***** JOB ABORTED ********************************************'   
  Write(* ,*) 'ERROR IN SBR. TASK_2: TMIN must be <= TMAX in KW Local_Study ' 
  Write(* ,*) '***** JOB ABORTED ********************************************';stop   
 Endif
endif
!
!
IF(ne==20540) THEN
 If (iv==0 .and. DELTAT < 0.0) then 
  Write(99,*) 'ERROR IN SBR. TASK_2: DELTAT must be >= 0. in KW Local_Study ' 
  Write(99,*) '***** JOB ABORTED ********************************************';stop    
 Endif
 If (iv==1 .and. DELTAT < 0.0) then 
  Write(99,*) 'ERROR IN SBR. TASK_2: DELTAT must be >= 0. in KW Local_Study ' 
  Write(99,*) '***** JOB ABORTED ********************************************'
  Write(* ,*) 'ERROR IN SBR. TASK_2: DELTAT must be >= 0. in KW Local_Study ' 
  Write(* ,*) '***** JOB ABORTED ********************************************';stop   
 Endif
If (iv==0 .and. DELTAT == 0.0) then 
  Write(99,*) 'ERROR IN SBR. TASK_2: DELTAT is = 0. in KW Local_Study' 
  Write(99,*) '********************* JOB ABORTED ************************';stop    
 Endif
 If (iv==1 .and. DELTAT == 0.0) then 
  Write(99,*) 'ERROR IN SBR. TASK_2: DELTAT is = 0. in KW Local_Study' 
  Write(99,*) '********************* JOB ABORTED ************************'
  Write(*,*) 'ERROR IN SBR. TASK_2: DELTAT is = 0. in KW Local_Study' 
  Write(*,*) '********************* JOB ABORTED ************************';stop   
 Endif
endif
!
!
IF(ne==20550) THEN
  If (iv==0 .and. (long_obs(1) <0.d0 .OR. long_obs(1) > 360.d0)) Then 
   Write(99,*) 'ERROR IN SBR. TASK_2: The longitude of the meridian must be '
   Write(99,*) 'in the range [0:360] degrees. **** JOB ABORTED ************ '; stop  
  endif
  If (iv==1 .and. (long_obs(1) <0.d0 .OR. long_obs(1) > 360.d0)) Then 
   Write(* ,*) 'ERROR IN SBR. TASK_2: The longitude of the meridian must be '
   Write(* ,*) 'in the range [0:360] degrees. **** JOB ABORTED ************ '  
   Write(99,*) 'ERROR IN SBR. TASK_2: The longitude of the meridian must be '
   Write(99,*) 'in the range [0:360] degrees. **** JOB ABORTED ************ '; stop  
  Endif
endif
!
!
IF(ne==20560) THEN
  If(iv==0.and.(C1<0.d0.OR.C1>180.d0.OR.C2<0.d0.OR.C2>180.d0.OR.C1>=C2.OR.EPC<=0.d0)) Then 
   Write(99,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(99,*) 'C1<0 deg; C1>180 deg; C2<0 deg; C2>180 deg; C1 > C2; C1 = C2; ' 
   Write(99,*) 'EPC <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ********************************************';stop 
  Endif
  If(iv==1.and.(C1<0.d0.OR.C1>180.d0.OR.C2<0.d0.OR.C2>180.d0.OR.C1>=C2.OR.EPC<=0.d0)) Then 
   Write(99,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(99,*) 'C1<0 deg; C1>180 deg; C2<0 deg; C2>180 deg; C1 > C2; C1 = C2; ' 
   Write(99,*) 'EPC <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ******************************************* ' 
   Write(* ,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(* ,*) 'C1<0 deg; C1>180 deg; C2<0 deg; C2>180 deg; C1 > C2; C1 = C2; ' 
   Write(* ,*) 'EPC <= 0 has been met. Please Check the input file task_2.dat '          
   Write(* ,*) '***** JOB ABORTED ******************************************* ';stop  
  Endif
endif
!
!
IF(ne==20570) THEN
 IF(NOBS>NOBS_MAX .and. IV==0) THEN
  WRITE (99,*) 'ERROR IN SBR. TASK_2: The number of observers' 
  write (99,*) 'exceeds the maximum allowed', nobs_max,' *** JOB ABORTED *********'; stop 
 ENDIF 
 IF(NOBS>NOBS_MAX .and. IV==1) THEN 
  WRITE (99,*) 'ERROR IN SBR. TASK_2: The number of observers'
  write (99,*) 'exceeds the maximum allowed', nobs_max,' *** JOB ABORTED *********' 
  WRITE (* ,*) 'ERROR IN SBR. TASK_2: The number of observers'
  write (* ,*) 'exceeds the maximum allowed', nobs_max,' *** JOB ABORTED *********'; stop   
 ENDIF
endif
!
!
IF(ne==20580) THEN
  If (iv==0 .and. (cola_obs(1) <0.d0 .OR. cola_obs(1) > 180.d0)) Then 
   Write(99,*) 'ERROR IN SBR. TASK_2: The colatitude of the parallel must be '
   Write(99,*) 'in the range [0:180] degrees. **** JOB ABORTED ************* '; stop  
  endif
  If (iv==1 .and. (cola_obs(1) <0.d0 .OR. cola_obs(1) > 180.d0)) Then 
   Write(* ,*) 'ERROR IN SBR. TASK_2: The colatitude of the parallel must be '
   Write(* ,*) 'in the range [0:180] degrees. **** JOB ABORTED ************* '
   Write(99,*) 'ERROR IN SBR. TASK_2: The colatitude of the parallel must be '
   Write(99,*) 'in the range [0:180] degrees. **** JOB ABORTED ************* '; stop  
  Endif
endif
!
!
IF(ne==20590) THEN
  If(iv==0.and.(L1<0.d0.OR.L1>360.d0.OR.L2<0.d0.OR.L2>360.d0.OR.L1>=L2.OR.EPL<=0.d0)) Then 
   Write(99,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(99,*) 'L1<0 deg; L1>360 deg; L2<0 deg; L2>180 deg; L1 > L2; L1 = L2; ' 
   Write(99,*) 'EPL <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ********************************************';stop 
  Endif
  If(iv==1.and.(L1<0.d0.OR.L1>360.d0.OR.L2<0.d0.OR.L2>360.d0.OR.L1>=L2.OR.EPL<=0.d0)) Then 
   Write(99,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(99,*) 'L1<0 deg; L1>360 deg; L2<0 deg; L2>180 deg; L1 > L2; L1 = L2; ' 
   Write(99,*) 'EPL <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ********************************************'
   Write(* ,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(* ,*) 'L1<0 deg; L1>360 deg; L2<0 deg; L2>180 deg; L1 > L2; L1 = L2; ' 
   Write(* ,*) 'EPL <= 0 has been met. Please Check the input file task_2.dat '          
   Write(* ,*) '***** JOB ABORTED ********************************************';stop  
  Endif 
endif
!
!
IF(ne==20600) THEN
  If(iv==0.and.(L1<0.d0.OR.L1>360.d0.OR.L2<0.d0.OR.L2>360.d0.OR.L1>=L2.OR.EPL<=0.d0)) Then
   Write(99,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(99,*) 'L1<0 deg; L1>360 deg; L2<0 deg; L2>360 deg; L1 > L2; L1 = L2; '
   Write(99,*) 'EPL <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ********************************************';stop 
  Endif
  If(iv==1.and.(L1<0.d0.OR.L1>360.d0.OR.L2<0.d0.OR.L2>360.d0.OR.L1>=L2.OR.EPL<=0.d0)) Then 
   Write(99,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(99,*) 'L1<0 deg; L1>360 deg; L2<0 deg; L2>360 deg; L1 > L2; L1 = L2; '
   Write(99,*) 'EPL <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ********************************************'
   Write(* ,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(* ,*) 'L1<0 deg; L1>360 deg; L2<0 deg; L2>360 deg; L1 > L2; L1 = L2; ' 
   Write(* ,*) 'EPL <= 0 has been met. Please Check the input file task_2.dat '          
   Write(* ,*) '***** JOB ABORTED ********************************************';stop  
  Endif 
endif
!
!
IF(ne==20610) THEN
  If(iv==0.and.(C1<0.d0.OR.C1>180.d0.OR.C2<0.d0.OR.C2>180.d0.OR.C1>=C2.OR.EPC<=0.d0)) Then 
   Write(99,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(99,*) 'C1<0 deg; C1>180 deg; C2<0 deg; C2>180 deg; C1 > C2; C1 = C2; ' 
   Write(99,*) 'EPC <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ********************************************';stop 
  Endif
  If(iv==1.and.(C1<0.d0.OR.C1>180.d0.OR.C2<0.d0.OR.C2>180.d0.OR.C1>=C2.OR.EPC<=0.d0)) Then 
   Write(99,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(99,*) 'C1<0 deg; C1>180 deg; C2<0 deg; C2>180 deg; C1 > C2; C1 = C2; ' 
   Write(99,*) 'EPC <= 0 has been met. Please Check the input file task_2.dat '          
   Write(99,*) '***** JOB ABORTED ******************************************* ' 
   Write(* ,*) 'ERROR IN SBR. TASK_2:        One of the forbidden conditions: '
   Write(* ,*) 'C1<0 deg; C1>180 deg; C2<0 deg; C2>180 deg; C1 > C2; C1 = C2; ' 
   Write(* ,*) 'EPC <= 0 has been met. Please Check the input file task_2.dat '          
   Write(* ,*) '***** JOB ABORTED ******************************************* ';stop  
  Endif
endif
!
!
IF(ne==20620) THEN
    IF(ihist==0.and.IV==0) then
           WRITE(99,*) 'ERROR IN SBR.   TASK_2: The Kw Load_Histoty'
           WRITE(99,*) 'MUST be activated before Global_Study can be'
           WRITE(99,*) 'executed. **********  JOB ABORTED ********* ';STOP
    endif
    IF(ihist==0.and.IV==1) then
           WRITE(99,*) 'ERROR IN SBR.   TASK_2: The Kw Load_Histoty'
           WRITE(99,*) 'MUST be activated before Global_Study can be'
           WRITE(99,*) 'executed. **********  JOB ABORTED ********* ' 
           WRITE(* ,*) 'ERROR IN SBR.   TASK_2: The Kw Load_Histoty'
           WRITE(* ,*) 'MUST be activated before Global_Study can be'
           WRITE(* ,*) 'executed. **********  JOB ABORTED ********* ';STOP
    endif
endif
!
!
IF(ne==20630) THEN
    IF(n_este >=2 .and. IV==0) then
      WRITE(99,*)'ERROR IN SBR. TASK_2:  The KW Global_Study appears' 
      WRITE(99,*) n_este, 'times in file task_2.dat ** JOB ABORTED ***********';STOP
    endif
    IF(n_este >=2 .and. IV==1) then
      WRITE(* ,*)'ERROR IN SBR. TASK_2:  The KW Global_Study appears' 
      WRITE(* ,*) n_este, 'times in file task_2.dat ** JOB ABORTED ***********'
      WRITE(99,*)'ERROR IN SBR. TASK_2:  The KW Global_Study appears' 
      WRITE(99,*) n_este, 'times in file task_2.dat ** JOB ABORTED ***********';STOP
    endif
endif
!
!
IF(ne==20640) THEN
    IF(n_inte ==1 .and. IV==0) then
      WRITE(99,*)'ERROR IN SBR. TASK_2:  The KWs Global_ and Local_Study are' 
      WRITE(99,*)'mutually exclusive. Both are active. ***** JOB ABORTED **********';STOP
    endif
    IF(n_inte ==1 .and. IV==1) then
      WRITE(99,*)'ERROR IN SBR. TASK_2:  The KWs Global_ and Local_Study are' 
      WRITE(99,*)'mutually exclusive. Both are active. ***** JOB ABORTED **********'
      WRITE(* ,*)'ERROR IN SBR. TASK_2:  The KWs Global_ and Local_Study are' 
      WRITE(* ,*)'mutually exclusive. Both are active. ***** JOB ABORTED **********';STOP
    endif
endif
!
!
IF(ne==20650) THEN
 If ((type_esten <=0 .or. type_esten >3).and.iv==0) then
  Write(99,*) 'ERROR IN SBR. TASK_2:  the Global_Study label ' 
  Write(99,*) 'must vary in the range [1:3] *** JOB ABORTED *** '; stop 
 Endif
 If ((type_esten <=0 .or. type_esten >3).and.iv==1) then
  Write(99,*) 'ERROR IN SBR. TASK_2:  the Global_Study label ' 
  Write(99,*) 'must vary in the range [1:3] *** JOB ABORTED *** ' 
  Write(* ,*) 'ERROR IN SBR. TASK_2:  the Global_Study label ' 
  Write(* ,*) 'must vary in the range [1:3] *** JOB ABORTED *** '; stop 
 Endif
endif
!
!
IF(ne==20660) THEN
IF (iv==0 .and. (LDE<LMIN.OR.LDE>LMAX.OR.MDE<0.OR.LDE<0.OR.IABS(MDE)>LDE))  THEN
Write(99,*)'ERROR IN SBR. TASK_2: One of the following forbidden conditions has been met:'
Write(99,*)'LDE < LMIN; LDE > LMAX; MDE < 0; |MDE| > LDE. Check the input file task_2.dat' 
Write(99,*)'************** JOB ABORTED **************************************************';stop   
ENDIF 
IF (iv==1 .and. (LDE<LMIN.OR.LDE>LMAX.OR.MDE<0.OR.LDE<0.OR.IABS(MDE)>LDE))  THEN 
Write(99,*)'ERROR IN SBR. TASK_2: One of the following forbidden conditions has been met:'
Write(99,*)'LDE < LMIN; LDE > LMAX; MDE < 0; |MDE| > LDE. Check the input file task_2.dat' 
Write(99,*)'************** JOB ABORTED **************************************************'  
Write(* ,*)'ERROR IN SBR. TASK_2: One of the following forbidden conditions has been met:'
Write(*, *)'LDE < LMIN; LDE > LMAX; MDE < 0; |MDE| > LDE. Check the input file task_2.dat' 
Write(* ,*)'************** JOB ABORTED **************************************************';stop 
ENDIF
!
  IF (iv==0 .and. (Lde > 36 .or. Lde < 2)) Then 
  Write(99,*) 'ERROR IN SBR. TASK_2: The range of harmonics degrees for' 
  write(99,*) 'this kind of analysis is [2:36] ************ JOB ABORTED **'; Stop 
  endif
  IF (iv==1 .and. (Lde > 36 .or. Lde < 2)) Then 
  Write(99,*) 'ERROR IN SBR. TASK_2: The range of harmonics degrees for' 
  write(99,*) 'this kind of analysis is [2:36] ************ JOB ABORTED **'
  Write(* ,*) 'ERROR IN SBR. TASK_2: The range of harmonics degrees for' 
  write(* ,*) 'this kind of analysis is [2:36] ************ JOB ABORTED **'; Stop    
  endif 
endif
!
!
IF(ne==20670) THEN
IF (iv==0 .and. Inorm /= 0 .and. Inorm /= 1) Then
Write(99,*)'ERROR IN SBR. TASK_2: 0 or 1 is expected for the switch Inorm. JOB ABORTED'; Stop  
endif
IF (iv==1 .and. Inorm /= 0 .and. Inorm /= 1) Then 
Write(99,*)'ERROR IN SBR. TASK_2: 0 or 1 is expected for the switch Inorm. JOB ABORTED'  
Write(* ,*)'ERROR IN SBR. TASK_2: 0 or 1 is expected for the switch Inorm. JOB ABORTED'; Stop    
endif 
endif
!
!
IF(ne==20680) THEN
 If (iv==0 .and. tmin > tmax) then
  Write(99,*) 'ERROR IN SBR. TASK_2: TMIN must be <= TMAX in KW Global_Study' 
  Write(99,*) '***** JOB ABORTED ********************************************';stop    
 Endif
 If (iv==1 .and. tmin > tmax) then 
  Write(99,*) 'ERROR IN SBR. TASK_2: TMIN must be <= TMAX in KW Global_Study'  
  Write(99,*) '***** JOB ABORTED ********************************************'   
  Write(* ,*) 'ERROR IN SBR. TASK_2: TMIN must be <= TMAX in KW Global_Study' 
  Write(* ,*) '***** JOB ABORTED ********************************************';stop   
 Endif
 If (iv==0 .and. DELTAT <= 0.0) then
  Write(99,*) 'ERROR IN SBR. TASK_2: DELTAT must be > 0. in KW Global_Study'
  Write(99,*) '***** JOB ABORTED ********************************************';stop    
 Endif
 If (iv==1 .and. DELTAT <= 0.0) then
  Write(99,*) 'ERROR IN SBR. TASK_2: DELTAT must be > 0. in KW Global_Study'
  Write(99,*) '***** JOB ABORTED ********************************************'
  Write(* ,*) 'ERROR IN SBR. TASK_2: DELTAT must be > 0. in KW Global_Study'
  Write(* ,*) '***** JOB ABORTED ********************************************';stop   
 Endif
endif
!
!
IF(ne==20690) THEN
If(iv==0.and.(Lde1<lmin.or.Lde1>lmax.or.Lde2<lmin.or.lde2>lmax.or.Lde1>=lde2.or.Lde1<0.or.Lde2<0))THEN
Write(99,*)'ERROR IN SBR. TASK_2: One of the following forbidden conditions has been met:'
Write(99,*)'LDE1>=LDE2;  LDE1<0;  LDE2<0;  LDE1<LMIN;  LDE1>LMAX;  LDE2<LMIN;  LDE2>LMAX '                              
WRITE(99,*)'***************************** JOB ABORTED ********************************** ';Stop 
Endif
If(iv==1.and.(Lde1<lmin.or.Lde1>lmax.or.Lde2<lmin.or.lde2>lmax.or.Lde1>=lde2.or.Lde1<0.or.Lde2<0))THEN
Write(99,*)'ERROR IN SBR. TASK_2: One of the following forbidden conditions has been met:'
Write(99,*)'LDE1>LDE2;  LDE1<0;  LDE2<0;  LDE1<LMIN;  LDE1>LMAX;  LDE2<LMIN;  LDE2>LMAX '                              
WRITE(99,*)'***************************** JOB ABORTED ********************************** '
Write(* ,*)'ERROR IN SBR. TASK_2: One of the following forbidden conditions has been met:'
Write(* ,*)'LDE1>LDE2;  LDE1<0;  LDE2<0;  LDE1<LMIN;  LDE1>LMAX;  LDE2<LMIN;  LDE2>LMAX '                              
WRITE(* ,*)'***************************** JOB ABORTED ********************************** ';stop
Endif
  IF (iv==0 .and. (Lde1 > 36 .or. Lde2 > 36 .or. Lde1 < 2 .or. Lde2 < 2)) Then
  Write(99,*) 'ERROR IN SBR. TASK_2: The range of harmonics degrees for' 
  write(99,*) 'this kind of analysis is [2:36] ************ JOB ABORTED **'; Stop 
  endif
  IF (iv==1 .and. (Lde1 > 36 .or. Lde2 > 36 .or. Lde1 < 2 .or. Lde2 < 2)) Then 
  Write(99,*) 'ERROR IN SBR. TASK_2: The range of harmonics degrees for' 
  write(99,*) 'this kind of analysis is [2:36] ************ JOB ABORTED **'
  Write(* ,*) 'ERROR IN SBR. TASK_2: The range of harmonics degrees for' 
  write(* ,*) 'this kind of analysis is [2:36] ************ JOB ABORTED **'; Stop    
  endif 
endif
!
!
IF(ne==20700) THEN
  If(iv==0 .and. lmin/=lmax .and. lmax/=2) Then 
  Write(99,*) 'ERROR in  TASK_2:     In order to compute the  INERTIA change' 
  Write(99,*) 'and derivatives, use L_min = L_max = 2 in KW HArmonic_Degrees'    
  Write(99,*) '******* JOB ABORTED ****************************************'; Stop
  Endif 
  If(iv==1 .and. lmin/=lmax .and. lmax/=2) Then 
  Write(*, *) 'ERROR in  TASK_2:     In order to compute the  INERTIA change' 
  Write(*, *) 'and derivatives, use L_min = L_max = 2 in KW HArmonic_Degrees'    
  Write(*, *) '******* JOB ABORTED ****************************************' 
  Write(99,*) 'ERROR in  TASK_2:     In order to compute the  INERTIA change' 
  Write(99,*) 'and derivatives, use L_min = L_max = 2 in KW HArmonic_Degrees'    
  Write(99,*) '******* JOB ABORTED ****************************************'; Stop
  Endif 
endif
!
!
  END subroutine m_2
!
!
!
!
!
!
 subroutine printdue (iup)
! - - - - - - - - - - - - - - - - - - - - -
! Prints file headers and more, for Task#2
! - - - - - - - - - - - - - - - - - - - - -
 USE for_task2
 implicit none
!
 INTEGER IUP
!

IF(iup==106) THEN
!	
     if(k==1) then
     do ijunk=48,54,1
     call DATE_AND_TIME (date,timc)
      Write(IJUNK,*) '# done by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
     enddo 
     endif 
!
IF (IR ==1 .and. k==1) THEN
           WRITE(48, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
!          
if(cl<=30.or.cl==50) WRITE(48, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE 
if(cl==40)           WRITE(48, '(a44,1x,E14.6)') '# Mass/surface at the load pole (kg/m**2) = ', MICE     
!
           WRITE(48, '(a46)') '# -------------------------------------------- '
           WRITE(48, '(a46)') '#      lon (deg)      lat (deg)      u_rad (m) '
           WRITE(48, '(a46)') '# -------------------------------------------- '
ENDIF
!
!
IF (IT ==1 .and. k==1) THEN
           WRITE(49, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
!          
if(cl<=30.or.cl==50) WRITE(49, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE 
if(cl==40)           WRITE(49, '(a44,1x,E14.6)') '# Mass/surface at the load pole (kg/m**2) = ', MICE     
!
           WRITE(49, '(a46)') '# -------------------------------------------- '
           WRITE(49, '(a46)') '#      lon (deg)      lat (deg)      u_lon (m) '
           WRITE(49, '(a46)') '# -------------------------------------------- '
           WRITE(50, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
!          
if(cl<=30.or.cl==50) WRITE(50, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE 
if(cl==40)           WRITE(50, '(a44,1x,E14.6)') '# Mass/surface at the load pole (kg/m**2) = ', MICE     
!
           WRITE(50, '(a46)') '# -------------------------------------------- '
           WRITE(50, '(a46)') '#      lon (deg)      lat (deg)      u_lat (m) '
           WRITE(50, '(a46)') '# -------------------------------------------- '
ENDIF
!
!
IF (IG ==1 .and. k==1) THEN
           WRITE(51, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
!          
if(cl<=30.or.cl==50) WRITE(51, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE 
if(cl==40)           WRITE(51, '(a44,1x,E14.6)') '# Mass/surface at the load pole (kg/m**2) = ', MICE     
!
           WRITE(51, '(a46)') '# -------------------------------------------- '
           WRITE(51, '(a46)') '#      lon (deg)      lat (deg)      geoid (m) '
           WRITE(51, '(a46)') '# -------------------------------------------- '
ENDIF

!
!
IF (DIR ==1 .and. k==1) THEN
           WRITE(52, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
!          
if(cl<=30.or.cl==50) WRITE(52, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE 
if(cl==40)           WRITE(52, '(a44,1x,E14.6)') '# Mass/surface at the load pole (kg/m**2) = ', MICE     
!
           WRITE(52, '(a49)') '# ----------------------------------------------- '
           WRITE(52, '(a49)') '#      lon (deg)      lat (deg)     d_rad (mm/yr) '
           WRITE(52, '(a49)') '# ----------------------------------------------- '
ENDIF
!
IF (DIT ==1 .and. k==1) THEN
           WRITE(53, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
!          
if(cl<=30.or.cl==50) WRITE(53, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE 
if(cl==40)           WRITE(53, '(a44,1x,E14.6)') '# Mass/surface at the load pole (kg/m**2) = ', MICE     
!
           WRITE(53, '(a49)') '# ----------------------------------------------- '
           WRITE(53, '(a49)') '#      lon (deg)      lat (deg)     d_lon (mm/yr) '
           WRITE(53, '(a49)') '# ----------------------------------------------- '
           WRITE(54, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
!          
if(cl<=30.or.cl==50) WRITE(54, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE 
if(cl==40)           WRITE(54, '(a44,1x,E14.6)') '# Mass/surface at the load pole (kg/m**2) = ', MICE     
!
           WRITE(54, '(a49)') '# ----------------------------------------------- '
           WRITE(54, '(a49)') '#      lon (deg)      lat (deg)     d_lat (mm/yr) '
           WRITE(54, '(a49)') '# ----------------------------------------------- '
ENDIF
!
IF (DIG ==1 .and. k==1) THEN
           WRITE(55, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
!          
if(cl<=30.or.cl==50) WRITE(55, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE 
if(cl==40)           WRITE(55, '(a44,1x,E14.6)') '# Mass/surface at the load pole (kg/m**2) = ', MICE     
!
           WRITE(55, '(a49)') '# ----------------------------------------------- '
           WRITE(55, '(a49)') '#      lon (deg)      lat (deg)     d_geo (mm/yr) '
           WRITE(55, '(a49)') '# ----------------------------------------------- '
ENDIF
!
!
If(IR==1)   write (48,'(3(1x,f14.5))') OBS(1),90.Q0-OBS(2),   VECT(1)
If(it==1)   write (49,'(3(1x,f14.5))') OBS(1),90.Q0-OBS(2),   VECT(3) 
If(it==1)   write (50,'(3(1x,f14.5))') OBS(1),90.Q0-OBS(2),  -VECT(2)   
If(IG==1)   write (51,'(3(1x,f14.5))') OBS(1),90.Q0-OBS(2),   VECT(4)  
!
If(dIR==1)  write (52,'(3(1x,f14.5))') OBS(1),90.Q0-OBS(2), D_VECT(1)
If(dit==1)  write (53,'(3(1x,f14.5))') OBS(1),90.Q0-OBS(2), D_VECT(3)
If(dit==1)  write (54,'(3(1x,f14.5))') OBS(1),90.Q0-OBS(2),-D_VECT(2)
If(dIG==1)  write (55,'(3(1x,f14.5))') OBS(1),90.Q0-OBS(2), D_VECT(4)
!
!
ENDIF !(iup==106)
!
!
  IF(iup==105 .and. k==1 ) THEN
!
  call DATE_AND_TIME (date,timc)
!
  DO IJUNK=10,11
     Write(IJUNK, *) '#                '
     Write(IJUNK, *) '# done by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)  
     Write(IJUNK, *) '# Time  (kyrs) = ', TIME
  ENDDO
!  
  IF(type_inten==3) THEN
     WRITE(10, *) '# Longitude (deg) =',  obs(1)
     WRITE(11, *) '# Longitude (deg) =',  obs(1)
  ENDIF
  IF(type_inten==4) THEN
     WRITE(10, *) '# Colatitude (deg) =', obs(2)
     WRITE(11, *) '# Colatitude (deg) =', obs(2)
  ENDIF
!
  if(CL <= 30 .or. CL == 50) then 
               write(10,*) '# The mass of the load (kg) is =', mice
               write(11,*) '# The mass of the load (kg) is =', mice
  endif
!
  if(CL == 40) then 
               write(10,*) '# The mass/surface at the load pole (kg/m**2) is =', mice
	       write(11,*) '# The mass/surface at the load pole (kg/m**2) is =', mice
  endif
!  
  IF     (type_inten==3 .and. k==1) THEN
!
  WRITE(10, *) '# '
  WRITE(10, *) '#  Colat.    u_rad     u_col     u_lon     geoid    '
  WRITE(10, *) '#  (deg)      (m)       (m)       (m)       (m)     '
  WRITE(10, *) '# -------------------------------------------------'
  WRITE(11, *) '# '
  WRITE(11, *) '#  Colat.    d_rad     d_col     d_lon    d_geoid   '
  WRITE(11, *) '#  (kyrs)   (mm/yr)   (mm/yr)   (mm/yr)   (mm/yr)   '
  WRITE(11, *) '# -------------------------------------------------'
!
  ELSEIF (type_inten==4 .and. k==1) THEN
!
  WRITE(10, *) '# '
  WRITE(10, *) '#  Longi.    u_rad     u_col     u_lon     geoid    '
  WRITE(10, *) '#  (deg)      (m)       (m)       (m)       (m)     '
  WRITE(10, *) '# -------------------------------------------------'
  WRITE(11, *) '# '
  WRITE(11, *) '#  Longi.    d_rad     d_col     d_lon    d_geoid   '
  WRITE(11, *) '#  (kyrs)   (mm/yr)   (mm/yr)   (mm/yr)   (mm/yr)   '
  WRITE(11, *) '# -------------------------------------------------'
!
  ENDIF
!
  ENDIF
!
!
!
!
  IF(iup==101 .and. cl /= 40) THEN
!
  call DATE_AND_TIME (date,timc)
!
  Do JJ=10,11
  WRITE(jj, *)               '# '
  WRITE(jj, *)               '# '
      Write(JJ, *) '# done by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)    
  Write(JJ, '(a31,1x,f9.4)') '# Observer  longitude (deg) = ', OBS(1)
  Write(JJ, '(a31,1x,f9.4)') '# Observer colatitude (deg) = ', OBS(2)
  ENDDO
  WRITE(10, *) '# '
  WRITE(10, *) '#  Time      u_rad     u_col     u_lon     geoid       mass'
  WRITE(10, *) '# (kyrs)      (m)       (m)       (m)       (m)        (kg)'
  WRITE(10, *) '# --------------------------------------------------------------'
  WRITE(11, *) '# '
  WRITE(11, *) '#  Time      d_rad     d_col     d_lon    d_geoid      mass'
  WRITE(11, *) '# (kyrs)    (mm/yr)   (mm/yr)   (mm/yr)   (mm/yr)      (kg)'
  WRITE(11, *) '# --------------------------------------------------------------'
  ENDIF
!

  IF(iup==101 .and. cl == 40) THEN
!
  call DATE_AND_TIME (date,timc)
!
  Do JJ=10,11
  WRITE(jj, *)               '# '
  WRITE(jj, *)               '# '
      Write(JJ, *) '# done by Task#2 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6)    
  Write(JJ, '(a31,1x,f9.4)') '# Observer  longitude (deg) = ', OBS(1)
  Write(JJ, '(a31,1x,f9.4)') '# Observer colatitude (deg) = ', OBS(2)
  ENDDO
WRITE(10, *) '# '
WRITE(10, *) '#  Time      u_rad     u_col     u_lon     geoid     mass/surface'
WRITE(10, *) '# (kyrs)      (m)       (m)       (m)       (m)        (kg/m**2) '
WRITE(10, *) '# ---------------------------------------------------------------'
WRITE(11, *) '# '
WRITE(11, *) '#  Time      d_rad     d_col     d_lon    d_geoid    mass/surface'
WRITE(11, *) '# (kyrs)    (mm/yr)   (mm/yr)   (mm/yr)   (mm/yr)      (kg/m**2) '
WRITE(11, *) '# ---------------------------------------------------------------'
ENDIF
!
!
!
 end subroutine printdue
!
!
!
!
!
!
!
!
!
!
 subroutine printtre (iup)
!-------------------------------------
! Prints data and messages for Task#3
!-------------------------------------
 USE for_task3
 implicit none
!
 INTEGER IUP
!
 IF(iup==12) THEN
  IF(iv==1.AND.(MOD(k,10) ==0.or.k==1.or.k==NOBS).AND.  &
               (MOD(j,100)==0.or.j==1.or.j==Mel_eff))     THEN
   WRITE(*,'(a11,i3,a4,i4,a8,f9.4,a13,i3,a5,i3)')    &
   '% Observer', k, 'of ', NOBS, ' at time', time_BP,&
   ' BP, element ', J, ' of ', Mel_Eff
 ENDIF
 ENDIF
!
!
 IF(iup==3) THEN
  IF(iv==1.AND.(MOD(k,10)==0.or.k==1.or.k==NOBS).AND.  &
              (MOD(j,100)==0.or.j==1.or.j==Mel_eff))     THEN
   WRITE(*,'(a11,i6,a4,i6,a8,f9.4,a13,i3,a5,i3)')    &
   '$ Observer', k, 'of ', NOBS, ' at time', Given_Time, &
   ' BP, element ', J, ' of ', Mel_Eff
 ENDIF
 ENDIF
!
!
IF(iup==4) THEN
  IF(iv==1.AND.(MOD(j,100)==0.or.j==1.or.j==Mel_eff))  THEN
   WRITE(*,'(a10,i3,a4,i4,a8,f9.4,a13,i3,a5,i3)')    &
   '$ Couple ', k, 'of ', NOBS, ' at time', Given_Time, &
   ' BP, element ', J, ' of ', Mel_Eff
ENDIF
ENDIF
!
!
 IF(iup==5) THEN
  IF(iv==1.AND.(MOD(j,100)==0.or.j==1.or.j==Mel_eff))     THEN
   WRITE(*,'(a12,f9.4,a11,i3,a5,i3)')    &
   'Time (BP) =', time_BP, ', element ', J, ' of ', Mel_Eff
 ENDIF
 ENDIF
!
!
 IF(iup==31) THEN
!
WRITE(99,*) '*** Output Files ***'
IF(IR==1) &
WRITE(99,*) ' A_rad.dat: Long. (deg), Lat. (deg), Radial displacement (m)'
IF(IT==1) &
WRITE(99,*) ' A_lon.dat: Long. (deg), Lat. (deg), Displacement along longitude (m)'
IF(IT==1) &
WRITE(99,*) ' A_lat.dat: Long. (deg), Lat. (deg), Displacement along latitude (m)'
IF(IG==1) &
WRITE(99,*) ' A_geo.dat: Long. (deg), Lat. (deg), Geoid ondulation (m)'
IF(DIR==1) &
WRITE(99,*) 'Ad_rad.dat: Long. (deg), Lat. (deg), Velocity along radius (mm/yr)'
IF(DIT==1) &
WRITE(99,*) 'Ad_lon.dat: Long. (deg), Lat. (deg), Velocity along longitude (mm/yr)'
IF(DIT==1) &
WRITE(99,*) 'Ad_lat.dat: Long. (deg), Lat. (deg), Velocity along latitude (mm/yr)'
IF(DIG==1) &
WRITE(99,*) 'Ad_geo.dat: Long. (deg), Lat. (deg), Rate of geoid change (mm/yr)'
!
IF(i_thi==1) &
WRITE(99,*) 'load_thick.dat:  Long. (deg), Lat. (deg), Load thickness (m)'
!
IF(iv==1) THEN
WRITE(*,*) '*** Output Files ***'
IF(IR==1) &
WRITE(*,*) ' A_rad.dat: Long. (deg), Lat. (deg), Radial displacement (m)'
IF(IT==1) &
WRITE(*,*) ' A_lon.dat: Long. (deg), Lat. (deg), Displacement along longitude (m)'
IF(IT==1) &
WRITE(*,*) ' A_lat.dat: Long. (deg), Lat. (deg), Displacement along latitude (m)'
IF(IG==1) &
WRITE(*,*) ' A_geo.dat: Long. (deg), Lat. (deg), Geoid ondulation (m)'
IF(DIR==1) &
WRITE(*,*) 'Ad_rad.dat: Long. (deg), Lat. (deg), Velocity along radius (mm/yr)'
IF(DIT==1) &
WRITE(*,*) 'Ad_lon.dat: Long. (deg), Lat. (deg), Velocity along longitude (mm/yr)'
IF(DIT==1) &
WRITE(*,*) 'Ad_lat.dat: Long. (deg), Lat. (deg), Velocity along latitude (mm/yr)'
IF(DIG==1) &
WRITE(*,*) 'Ad_geo.dat: Long. (deg), Lat. (deg), Rate of geoid change (mm/yr)'
!
IF(i_thi==1) &
WRITE(*,*) 'load_thick.dat:  Long. (deg), Lat. (deg), Load thickness (m)'
ENDIF
!
 ENDIF
!
!
  IF(iup==101) THEN
!
  call DATE_AND_TIME (date,timc)
!
  if(ir==1.or.it==1.or.ig==1) then 
  WRITE(10, *)               '# '
  Write(10, *) '# done by Task#3 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  WRITE(10, *)               '# '	     
  Write(10, '(a31,1x,f9.4)') '# Observer  longitude (deg) = ', OBS(1)
  Write(10, '(a31,1x,f9.4)') '# Observer colatitude (deg) = ', OBS(2)
  WRITE(10, *) '# '
  WRITE(10, *) '# Time BP    u_rad     u_col     u_lon     geoid       mass'
  WRITE(10, *) '# (kyrs)      (m)       (m)       (m)       (m)        (kg)'
  WRITE(10, *) '# --------------------------------------------------------------'
  endif 
! 
  if(dir==1.or.dit==1.or.dig==1) then 
  WRITE(11, *)               '# '
  Write(11, *) '# done by Task#3 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  WRITE(11, *)               '# '	     
  Write(11, '(a31,1x,f9.4)') '# Observer  longitude (deg) = ', OBS(1)
  Write(11, '(a31,1x,f9.4)') '# Observer colatitude (deg) = ', OBS(2)
  WRITE(11, *) '# '
  WRITE(11, *) '# Time BP    d_rad     d_col     d_lon    d_geoid       mass'
  WRITE(11, *) '# (kyrs)    (mm/yr)   (mm/yr)   (mm/yr)   (mm/yr)       (kg)'
  WRITE(11, *) '# --------------------------------------------------------------'
  endif  
!  
  ENDIF
!
!
IF(iup==105) THEN
!
 call DATE_AND_TIME (date,timc)
  
  if(k==1) then 
  do JJ=48,56,1
  Write(JJ, *) '# done by Task#3 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
  enddo
  Endif 	     

IF (IR ==1 .and. k==1) THEN
           WRITE(48, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
           WRITE(48, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE
           WRITE(48, '(a46)') '# -------------------------------------------- '
           WRITE(48, '(a46)') '#      lon (deg)      lat (deg)      u_rad (m) '
           WRITE(48, '(a46)') '# -------------------------------------------- '
ENDIF
!
!
IF (I_THI==1 .and. k==1) THEN
           WRITE(56, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
           WRITE(56, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE
           WRITE(56, '(a46)') '# -------------------------------------------- '
           WRITE(56, '(a46)') '#      lon (deg)      lat (deg)      thick (m) '
           WRITE(56, '(a46)') '# -------------------------------------------- '
ENDIF
!
!
IF (IT ==1 .and. k==1) THEN
           WRITE(49, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
           WRITE(49, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE
           WRITE(49, '(a46)') '# -------------------------------------------- '
           WRITE(49, '(a46)') '#      lon (deg)      lat (deg)      u_lon (m) '
           WRITE(49, '(a46)') '# -------------------------------------------- '
           WRITE(50, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
           WRITE(50, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE
           WRITE(50, '(a46)') '# -------------------------------------------- '
           WRITE(50, '(a46)') '#      lon (deg)      lat (deg)      u_lat (m) '
           WRITE(50, '(a46)') '# -------------------------------------------- '
ENDIF
!
!
IF (IG ==1 .and. k==1) THEN
           WRITE(51, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
           WRITE(51, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE
           WRITE(51, '(a46)') '# -------------------------------------------- '
           WRITE(51, '(a46)') '#      lon (deg)      lat (deg)      geoid (m) '
           WRITE(51, '(a46)') '# -------------------------------------------- '
ENDIF

!
!
IF (DIR ==1 .and. k==1) THEN
           WRITE(52, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
           WRITE(52, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE
           WRITE(52, '(a49)') '# ----------------------------------------------- '
           WRITE(52, '(a49)') '#      lon (deg)      lat (deg)     d_rad (mm/yr) '
           WRITE(52, '(a49)') '# ----------------------------------------------- '
ENDIF
!
IF (DIT ==1 .and. k==1) THEN
           WRITE(53, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
           WRITE(53, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE
           WRITE(53, '(a49)') '# ----------------------------------------------- '
           WRITE(53, '(a49)') '#      lon (deg)      lat (deg)     d_lon (mm/yr) '
           WRITE(53, '(a49)') '# ----------------------------------------------- '
           WRITE(54, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
           WRITE(54, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE
           WRITE(54, '(a49)') '# ----------------------------------------------- '
           WRITE(54, '(a49)') '#      lon (deg)      lat (deg)     d_lat (mm/yr) '
           WRITE(54, '(a49)') '# ----------------------------------------------- '
ENDIF
!
IF (DIG ==1 .and. k==1) THEN
           WRITE(55, '(a19,1x, f9.4)') '# Time BP (kyrs) = ', GIVEN_TIME
           WRITE(55, '(a19,1x,E14.6)') '# Load Mass (kg) = ', MICE
           WRITE(55, '(a49)') '# ----------------------------------------------- '
           WRITE(55, '(a49)') '#      lon (deg)      lat (deg)     d_geo (mm/yr) '
           WRITE(55, '(a49)') '# ----------------------------------------------- '
ENDIF
!
!
If (IR==1)  write (48,'(3(1x,f14.5))') OBS(1),90.d0-OBS(2),     WECT(1)
If (IT==1)  write (49,'(3(1x,f14.5))') OBS(1),90.d0-OBS(2),     WECT(3)
If (IT==1)  write (50,'(3(1x,f14.5))') OBS(1),90.d0-OBS(2),    -WECT(2)
If (IG==1)  write (51,'(3(1x,f14.5))') OBS(1),90.d0-OBS(2),     WECT(4)
!
If(DIR==1)  write (52,'(3(1x,f14.5))') OBS(1),90.d0-OBS(2),   D_WECT(1)
If(DIT==1)  write (53,'(3(1x,f14.5))') OBS(1),90.d0-OBS(2),   D_WECT(3)
If(DIT==1)  write (54,'(3(1x,f14.5))') OBS(1),90.d0-OBS(2),  -D_WECT(2)
If(DIG==1)  write (55,'(3(1x,f14.5))') OBS(1),90.d0-OBS(2),   D_WECT(4)
!
If(I_thi==1)write (56,'(2(1x,f14.5),E20.8)') OBS(1),90.Q0-OBS(2), PRES/rhoice
!
ENDIF !(iup==105)
!
!
 IF(iup==129) THEN
! 
  call DATE_AND_TIME (date,timc)
!
  Write(55, *) '# done by Task#3 of TABOO on ', & 
                       date(1:4), '.', date(5:6), '.', date(7:8), & 
	     ' time=', timc(1:2), '.', timc(3:4), '.', timc(5:6) 
 
! 
  Write(55,*) '# The rates below are computed for time_BP =', given_time  
! 
  IF     (ilb==1) then
  WRITE(55,*) &
   '#-------------------------------------------------------------'
  WRITE(55,*) &
   '#Site#1  Long_1 Cola_1 Site#2  Long_2 Cola_2  _L_   _T_   _V_ '
  WRITE(55,*) &
   '#        (deg)  (deg)          (deg)  (deg)       (mm/yr)     '
  WRITE(55,*) &
   '#-------------------------------------------------------------'
  ELSEIF (ilb==2) then
  WRITE(55,*) &
   '#-------------------------------------------------------------'
  WRITE(55,*) &
   '#        Long_1 Cola_1         Long_2 Cola_2  _L_   _T_   _V_ '
  WRITE(55,*) &
   '#        (deg)  (deg)          (deg)  (deg)       (mm/yr)     '
  WRITE(55,*) &
   '#-------------------------------------------------------------'
  ENDIF
 ENDIF
!
!
 end subroutine printtre
!
!
!

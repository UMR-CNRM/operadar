!2012/12 C. Augros, Version V0
!Creation d'une table avec les S11, S22 integres sur l'ensemble des diametres,
!(à partir de la table issue de ampld.ClotildeV3.lp.f),
! pour differents rapports de mélange.
!LAM, et les coef integres sont convertis en m
!On crée donc une table pour chaque type d'hydromet (r, s, g), chaque gamme de LAM
!avec :
! LAM | Tc | ELEV | M | S11 | S22 ...
!rq : pour la glace primaire, il faudra faire une table spéciale avec M et Ci
! parce-que N depend de M et de Ci

!2013/02/20 C. Augros, Version4
!S22S11 nombre complexe => ecriture de REAL(S22S11) et AIMAG(S22S11)
!dans fichier integre
!suppression de sigbeta dans le nom de la table

!2013/10/24 C. Augros Suppression chemins pour lancement depuis Aneto
!les fichiers d'entree et l'executable doivent etre dans le meme repertoire

!2013/11/19 Version adaptée pour test fo,ction dielec et calcul des variables dpol 

!2015/11 Ajout glace primaire avec Ci=800 m-3 et Ni MésoNH

!2016/06 Mary: ajout de la fonction PSD et lecture des paramètres ICE3 dans un fichier indépendant qui s'appelle ICE3_constant.txt

!2018/01 Tony: Introduction des espèces fondantes "wg" et "ws"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Le code se découpe en 2 ou 3 étapes :
!
! Etape 0 : Réalisée seulement dans le cas où l'espèce est fondante
!           Détermination de l'équivalent en pluie M_liq du contenu de l'espèce fondante par intégration (somme de Riemman)
!           M_liq va servir à établir la distribution en taille des hydrométéores fondants (Szyrmer et al., 1999)                                
!                                                  
! Etape 1 : Calcul du contenu de l'espèce par intégration (somme de Riemman) Mint
!           La somme de Riemman étant une approximation de l'intégrale et l'intégration se faisant sur un domaine fini
!           (on n'intègre par sur la queue de la distribution), elle différe du contenu initial M.
!           Mint va donc servir à corriger la distribution en taille des hydrométéores (étape 2)
!
! Etape 2 : Calcul des variables radar
!           La distribution en taille de l'espèce est modifiée par le facteur M/Mint pour assurer que l'intégration finale
!           par somme de Riemman retombe sur le contenu initial M
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!====================
! TmatInt_201911.f90
!====================

! 2019/12 C. Augros
! Fusion of Tony and Natalia's codes from Tony's code
! includes improvement of melting layer simulation + bug corrections
! + LIMA (2-moments for rain and ice) + hail
! Fw replaced by P3: third parameter of the table (concentration No for 2-moment species
! or liquid water fraction for 1-moment melting species (snow, graupel)
! Translation of comments in english, correction of indentation
! Test mode: to check the reading of the T-matrix table => 
! line 120: set testMode=.FALSE. for normal use
! New reading of the table: the parameters (min/max/step) are read directly in the coefficient file

!2020/01 C. Augros
! min/max/step parameters including expM and expCC are writen in the first lines of the output files
! Optimization of the output files size: format modification +
! use of ";" as separator instead of spaces (2X)
! 4 first lines of TmatVarInt_LIMA_S106.2_rr02 :
!LAMmin;ELEVmin;ELEVstep;ELEVmax;Tcmin;Tcstep;Tcmax;expMmin;expMstep;expMmax;expCCmin;expCCstep;expCCmax;Fwmin;Fwstep;Fwmax
!106.20;  0.00;  4.00; 20.00;-20.00;  1.00; 40.00; -7.00;  0.01; -1.50;  1.80;  9.00;  0.10;  0.00;  0.10;  0.00
!Tc;ELEV;P3;M;zhhlg;zdrlg;rhohv;kdp
!-20.0; 0.0;0.631E+02;0.1000E-06;-.2294E+02;0.6094E-02;0.1000E+01;0.4759E-07

! 2020/01/23 - Possibility to read ICE3, or LIMA constant file:
!            Hydromet_const_file: 'ICE3_constants.txt' or 'LIMA_constants.txt'
!            CCLOUD needs to be modified to select ICE3, LIMA or LIMT
!            - Hydromet constants are consistent with MesoNH documentation 5.4 +code
!            NB: for hail, c and d parameters in the fall velocity relation (eq 7.39 in MesoNH doc)are:
!            c=201 for LIMA (value in the code+ doc: tab 7.1 p 147), 207 for ICE3 (value in the code)
!            c=207, d=0.64 correspond to Bohm (1989)
!            But these values correspond to rho_d=1.1kg/m3, whereas in Meso-NH, c and d are used
!            as if they were calculated for rho_0=1.2041kg/m3.
!            This is corrected for LIMA (=> c=201)
!            - Modification of subroutine PSD to allow the possibility of 
!            CCLOUD='LIMT' => LIMA with diagnostic mu (cf Taufour et al 2018)
!            - New step for number concentration (to reduce computation time): 
!            expCCmin=1.8,expCCmax=9.,expCCstep=0.1 !0.02

!====================
! TmatInt_202002.f90
!====================

! 2020/02/07 C. Augros
! Add new variable: Nmoments (1 or 2)
! 2-moments for ice for ICE3, LIMA and LIMT

! 2020/04/09 C. Augros
! New diameter steps with exponential range: more diameters for the lowest
! expD,expDmin,expDmax,expDstep
! new parameter: mumax for LIMT (mumax=15 in B. Vie's code version of Taufour et al 2018)

!========================
! Tmat_202004/TmatInt.f90
!========================
! 2020/04/21 C. Augros
! Correction of the integration over the PSD: use of the equivalent diameter
! (partially melted if mixed phase)
! instead of the melted diameter
! For this new integration version, all scattering coef + dpol variables 
! are named with 2 (ex S11carre2, zhh2 ...)
! The former version is also calculated but the parameters saved in the outfiles 
! are integrated with the new version  

! To do:
! - add the version suggested by Wolfensberger 2018 for the distribution of 
! wet hydrometer (linear transition between rain and graupel distribution)
! - for wet graupel: use the concentration of rain from LIMA and LIMT instead 
! of the 1-moment rain concentration

PROGRAM MakeTmatInt

use iso_fortran_env
IMPLICIT NONE

!************* Declarations ***************     

! Lookup table parameters
CHARACTER*6::canting
CHARACTER*7::DIEL
CHARACTER*4::ARfunc
CHARACTER*4::DSTYfunc
REAL:: SIGBETA,ARcnst
REAL :: LAM,FREQ 
REAL :: ELEV,ELEVmin,ELEVmax,ELEVstep !,ELEVrecinf,ELEVrecsup
REAL :: Tc,Tcmin,Tcmax,Tcstep,Tk !,Tcrecinf,Tcrecsup
REAL :: Fwmin, Fwstep, Fwmax,Fw !,Fwrecinf,Fwrecsup
REAL :: expDmin,expDmax,expDstep
REAL :: Frim ! degree of riming factor

! Sup and inf born for diameter parameters
REAL :: Drecsup,Dmrecsup,Deqrecsup,Deqmrecsup,Deqrrecsup,Deqrmrecsup
REAL :: Drecinf,Dmrecinf,Deqrecinf,Deqmrecinf,Deqrrecinf,Deqrmrecinf

! Indexes for lookup table
INTEGER :: idELEV,idTc,idD,idM,iLine,idP3,idFw !,ispecies
INTEGER :: nELEV,nTc,nD,nM,nFw,nLines
INTEGER :: nP3 !nCC,
INTEGER :: nELEVloop, nTcloop, nP3loop, nMloop

! Density and dielectric constants
REAL :: RHOLW,QMW !RHOX,Vpart,Mpart,
COMPLEX :: EPSW
REAL :: K2 !dielectric factor

! Variables for input and output files names
CHARACTER*2 :: typeh
CHARACTER*2 :: band
CHARACTER*4  :: CCLOUD ! LIMA ICE3 LIMC ICJW
CHARACTER*2 :: espece, MOMENT
CHARACTER*2 :: espece_rr ! rain
character(len=1024) :: nomfileCoef, nomfileCoef_rr, nomfileCoefInt, Hydromet_const_file
integer :: ios
character(len=256) :: exec_dir, repo

! Variables to be read in the table (Tmatrix and Rayleigh scattering)
REAL ::ELEVrec,Tcrec,Fwrec
REAL ::Drec,Dmrec,Deqrec,Deqmrec,Deqrrec,Deqrmrec
REAL :: sighh,sigvv,kdp,REdeltaco,IMdeltaco,Ah,Av
REAL :: sighhR,sigvvR,kdpR,REdeltacoR,IMdeltacoR,AhR,AvR
COMPLEX :: deltaco,deltacoR

! Arrays with all tmatrix table values  
REAL, DIMENSION(:),ALLOCATABLE ::ELEVtab,Tctab,Fwtab
REAL, DIMENSION(:),ALLOCATABLE ::Dtab,Dmtab,Deqtab,Deqmtab,Deqrtab,Deqrmtab
REAL, DIMENSION(:),ALLOCATABLE :: sighhtab,sigvvtab,REdeltacotab,IMdeltacotab
REAL, DIMENSION(:),ALLOCATABLE :: zhhtab,zdrtab,kdptab,rhohvtab,Ahtab,Avtab
REAL, DIMENSION(:),ALLOCATABLE :: sighhRtab,sigvvRtab
REAL, DIMENSION(:),ALLOCATABLE :: REdeltacoRtab,IMdeltacoRtab
REAL, DIMENSION(:),ALLOCATABLE :: zhhRtab,zdrRtab,kdpRtab,rhohvRtab
REAL, DIMENSION(:),ALLOCATABLE :: AhRtab,AvRtab


! Variables needed for integrated scattering coef and dpol variables          
INTEGER :: kTmat,Nmoments
REAL :: expMmin, expMstep, expMmax, expM, M
REAL :: expCCmin, expCCstep, expCCmax !,CC ! rain and ice concentration for LIMA
REAL :: aj,bj,nuj,alphaj,Cj,Xj,ccj,ddj ! ICE3 param
REAL :: aj_rr,bj_rr,nuj_rr,alphaj_rr,Cj_rr,Xj_rr,ccj_rr,ddj_rr ! ICE3 param (rain)
REAL :: P3min, P3step, P3max,P3,expP3
REAL*8 :: lamj           ! slope param and DSP 
REAL :: mumax,radar_cnst !max mu for LIMA with diagnostic mu (Taufour et al 2018)
REAL*8 :: N
REAL*8 :: N_rr,N_ss !,N_rr_inf,N_s2r,N_r2s,N_ss,N_ss_inf,N_d ! rain
REAL*8 :: N2, N_rr2,N_ss2
REAL :: vtr,vts,vtm,M_liq,M_ss,M_rr,Mint,phi,RRint !,KAPPA
REAL :: M_rr2, M_ss2,Mint2,RRint2

! Integrated coefficients and dual pol variables
REAL :: sighhint,sigvvint
REAL :: sighhRint,sigvvRint
COMPLEX :: deltacoint,deltacoRint
REAL :: REdeltacoint,IMdeltacoint,REdeltacoRint, IMdeltacoRint
REAL :: zhhint,zdrint,kdpint,rhohvint,Ahint,Avint!,Deltaz,Adp,Ah,Av
REAL :: zhhRint,zdrRint,kdpRint,rhohvRint,AhRint,AvRint
LOGICAL:: testMode



!**************************************************
!                    Constants
!**************************************************
REAL :: P=DACOS(-1D0) !Pi
REAL :: LIGHTSPEED = 299792458.
RHOLW=1000.


!**************************************************
!              Handle input arguments 
!**************************************************
call get_command_argument(1, exec_dir)
call get_command_argument(2, repo)
call get_command_argument(3, typeh)
call get_command_argument(4, band)
call get_command_argument(5, CCLOUD)
call get_command_argument(6, MOMENT)


!**************************************************
!                    Parameters 
!**************************************************

mumax=15.0 !6.0
testMode=.FALSE. !.TRUE.

! min/step/max values for the exponant of the hydromet content M (kg/m3):
expMmin=-7.0 !-7.0
expMstep=0.01 
expMmax=-1.5 !-1.5

! min/step/max values for the exponant of 2-moment species concentration(/m3)
expCCmin=0.0 
expCCmax=9.0
expCCstep=0.1


!*****************************************************
!               BEGINNING OF THE PROGRAM             !
!*****************************************************

WRITE(0,*) ' -----------------------------------------------------------'
WRITE (0,*) ' Conversion of the table as a function of the diameter into a table as a function of the hydrometeor content.'
WRITE (0,*) '    microphysics=',CCLOUD
WRITE (0,*) '    type=',typeh,' (',MOMENT,')'
WRITE (0,*) '    band=',trim(band)

IF (typeh=='ii') THEN
  Nmoments=2
ELSE
  IF (MOMENT=='1M') THEN 
    Nmoments=1
  ELSE IF (MOMENT=='2M') THEN 
    Nmoments=2
  ENDIF
ENDIF

!----- Reading of the hydrometeor constants for PSD calculations
Hydromet_const_file = trim(exec_dir)//'/../param/'//CCLOUD//'_constants.txt'
WRITE(0,*) ' Reading : ',trim(Hydromet_const_file)

open(unit=333,file = trim(Hydromet_const_file),action="read",iostat=ios)
read(333,*) ! Read the first line of the file
DO while(espece /= typeh) 
  read(333,*) espece,aj,bj,alphaj,nuj,Xj,Cj,ccj,ddj
ENDDO
WRITE(0,*) '    PSD constants : ',espece,aj,bj,alphaj,nuj,Xj,Cj,ccj,ddj
close(333)

!----- Reading of the rain constants 
open(unit=333,file = trim(Hydromet_const_file),action="read",iostat=ios)
read(333,*) ! Read the first line of the file
DO while(espece_rr /= 'rr') 
  read(333,*) espece_rr,aj_rr,bj_rr,alphaj_rr,nuj_rr,Xj_rr,Cj_rr,&
              ccj_rr,ddj_rr
ENDDO
close(333)


!---- Input Files
nomfileCoef  = trim(exec_dir)//'/../tables/'//trim(repo)//'/TmatCoefDiff_'//trim(band)//typeh
nomfileCoef_rr = trim(exec_dir)//'/../tables/'//trim(repo)//'/TmatCoefDiff_'//trim(band)//'rr'

WRITE(0,*) ' Reading ',trim(nomfileCoef)

!------ Reading of min/step/max parameters in the 2nd line 
! of coefficients file
OPEN (1,FILE=trim(nomfileCoef),action="read")
READ (1,*) ! Reads the first line of the file
1001 FORMAT (A6,2X,F4.1,2X,A7,2X,A4,2X,F4.2,2X,A4,2X,F6.2,2X, & !2nd line
              F5.1,2X,F6.2,2X,F6.2,2X,F6.2,2X,F6.2,2X,F6.2,2X, &
              F6.2,2X,F6.2,2X,F6.4,2X,F6.2,2X,F6.2,2X,F6.2,2X,F6.2)
READ (1,1001) canting,SIGBETA,DIEL,ARfunc,ARcnst,DSTYfunc,LAM, &
              Frim,ELEVmin,ELEVstep,ELEVmax,Tcmin,Tcstep,Tcmax, &
              expDmin,expDstep,expDmax,Fwmin,Fwstep,Fwmax

WRITE (0,*)'    --> File first line : '
WRITE (0,*)"        ","canting,SIGBETA,DIEL,ARfunc,ARcnst,DSTYfunc,LAM,", &
            "Frim,ELEVmin,ELEVstep,ELEVmax,Tcmin,Tcstep,Tcmax,", &
            "expDmin,expDstep,expDmax,Fwmin,Fwstep,Fwmax"
WRITE(0,'(A)', ADVANCE='NO') '          '
WRITE (0,1001) canting,SIGBETA,DIEL,ARfunc,ARcnst,DSTYfunc,LAM, &
                Frim,ELEVmin,ELEVstep,ELEVmax,Tcmin,Tcstep,Tcmax, &
                expDmin,expDstep,expDmax,Fwmin,Fwstep,Fwmax

WRITE(0,*) ' Reading ',trim(nomfileCoef_rr)

nTc=nint((Tcmax-Tcmin)/Tcstep)+1
nELEV=nint((ELEVmax-ELEVmin)/ELEVstep)+1
nD=nint((expDmax-expDmin)/expDstep)+1
nFw=nint((Fwmax-Fwmin)/Fwstep)+1
nLines=nTc*nELEV*nFw*nD ! number of lines of the Tmat file

WRITE(0,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
        "     Number of lines in the files : ", nLines, &
            " (nD=",nD,", nELEV=",nELEV,", nTc=",nTc,", nFw=",nFw,")"

! Allocate tables for tmatrix table values
ALLOCATE (Tctab(nLines),ELEVtab(nLines),&
          Fwtab(nLines),Dtab(nLines),Dmtab(nLines),Deqtab(nLines),&
          Deqmtab(nLines),Deqrtab(nLines),Deqrmtab(nLines),&
          sighhtab(nLines),sigvvtab(nLines),REdeltacotab(nLines),&
          IMdeltacotab(nLines),zhhtab(nLines),zdrtab(nLines),&
          kdptab(nLines),rhohvtab(nLines),Ahtab(nLines),Avtab(nLines),&
          sighhRtab(nLines),sigvvRtab(nLines),REdeltacoRtab(nLines),&
          IMdeltacoRtab(nLines),zhhRtab(nLines),zdrRtab(nLines),& 
          kdpRtab(nLines),rhohvRtab(nLines),AhRtab(nLines),AvRtab(nLines))
          
5002 FORMAT (F5.1,2X,F4.1,2X,F4.2,2X,E10.4,2X,E10.4,2X,E10.4,&
      2X,E10.4,2X,E10.4,2X,E10.4,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,&
      2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,&
      2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,&
      2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5)      

READ (1,*) ! Reads the 3d line of the file
iLine=1
DO WHILE (iLine .LE. nLines) ! Reads lines 4 to nLines
  READ( 1,5002) Tctab(iLine),ELEVtab(iLine),&
                Fwtab(iLine),Dtab(iLine),Dmtab(iLine),Deqtab(iLine),&
            Deqmtab(iLine),Deqrtab(iLine),Deqrmtab(iLine),&
            sighhtab(iLine),sigvvtab(iLine),REdeltacotab(iLine),&
            IMdeltacotab(iLine),zhhtab(iLine),zdrtab(iLine),&
            kdptab(iLine),rhohvtab(iLine),Ahtab(iLine),Avtab(iLine),& 
            sighhRtab(iLine),sigvvRtab(iLine),REdeltacoRtab(iLine),&
            IMdeltacoRtab(iLine),zhhRtab(iLine),zdrRtab(iLine),& 
            kdpRtab(iLine),rhohvRtab(iLine),AhRtab(iLine),AvRtab(iLine)
    iLine=iLine+1
  ENDDO
CLOSE(1)


! ======== TEST mode =====================
! To check the correct reading of the table
IF (testMode) THEN
  WRITE(0,*) "====== Test Mode ======"
  WRITE(0,*) " 4th line of ",nomfileCoef
  WRITE (0,*) "Tc,ELEV,Fw,D,Dm,Deq,Deqm,Deqr,Deqrm,sighh",&
        "sigvv,REAL(deltaco),AIMAG(deltaco),zhh,zdr,kdp,rhohv,Ah,Av,sighhR",&
        "sigvvR,REAL(deltacoR),AIMAG(deltacoR),zhhR,zdrR,kdpR,rhohvR,AhR,AvR"
  WRITE(0,5002) Tctab(1),ELEVtab(1),&
                Fwtab(1),Dtab(1),Dmtab(1),Deqtab(1),&
                Deqmtab(1),Deqrtab(1),Deqrmtab(1),&
                sighhtab(1),sigvvtab(1),REdeltacotab(1),&
                IMdeltacotab(1),zhhtab(1),zdrtab(1),kdptab(1),rhohvtab(1),&
                Ahtab(1),Avtab(1),& 
                sighhRtab(1),sigvvRtab(1),REdeltacoRtab(1),IMdeltacoRtab(1),& 
                zhhRtab(1),zdrRtab(1),kdpRtab(1),rhohvRtab(1),&
                AhRtab(1),AvRtab(1)
  WRITE(0,*) " Last line of Tmatrix table :",nLines
  WRITE(0,5002) Tctab(nLines),ELEVtab(nLines),&
                Fwtab(nLines),Dtab(nLines),Dmtab(nLines),Deqtab(nLines),&
                Deqmtab(nLines),Deqrtab(nLines),Deqrmtab(nLines),&
                sighhtab(nLines),sigvvtab(nLines),REdeltacotab(nLines),&
                IMdeltacotab(nLines),zhhtab(nLines),zdrtab(nLines),&
                kdptab(nLines),rhohvtab(nLines),Ahtab(nLines),Avtab(nLines),&
                sighhRtab(nLines),sigvvRtab(nLines),&
                REdeltacoRtab(nLines),IMdeltacoRtab(nLines),& 
                zhhRtab(nLines),zdrRtab(nLines),kdpRtab(nLines),& 
                rhohvRtab(nLines),AhRtab(nLines),AvRtab(nLines)
ENDIF
! ==========================================  



!=================================
!Output file

IF (.not. testMode) THEN
    nomfileCoefInt = trim(exec_dir)//'/../tables/'//trim(repo)//'/TmatCoefInt_'//CCLOUD//'_'//MOMENT//'_'//trim(band)//typeh
    PRINT *,' Creation of ',trim(nomfileCoefInt)
    
    OPEN (6,FILE=trim(nomfileCoefInt))
    
    ! ==== Description part of the output file ====
    ! Parameters values
    WRITE (6,1000)              
    1000 FORMAT ('micro;canting;SIGBETA;DIEL;ARfunc;ARcnst;DSTYfunc;LAM;' &
                  'Frim;ELEVmin;ELEVstep;ELEVmax;Tcmin;Tcstep;Tcmax;' &
                  'expMmin;expMstep;expMmax;expCCmin;expCCstep;expCCmax;' &
                  'Fwmin;Fwstep;Fwmax')
    WRITE (6,1002) CCLOUD,canting,SIGBETA,DIEL,ARfunc,ARcnst,DSTYfunc,LAM, &
                    Frim,ELEVmin,ELEVstep,ELEVmax,Tcmin,Tcstep,Tcmax, &
                    expMmin,expMstep,expMmax,expCCmin,expCCstep,expCCmax, &
                    Fwmin,Fwstep,Fwmax
    1002 FORMAT (A4,';',A6,';',F4.1,';',A7,';',A4,';',F4.2,';',A4,';',F6.2,';', &
                  F5.1,';',F6.2,';',F6.2,';',F6.2,';',F6.2,';',F6.2,';',F6.2,';', &
                  F6.2,';',F6.2,';',F6.2,';',F6.2,';',F6.2,';',F6.2,';', &
                  F6.2,';',F6.2,';',F6.2)               
        
    ! Variables names
    WRITE (6,5003)
    5003 FORMAT ('Tc;ELEV;P3;M;RR;' &
                  'sighh;sigvv;REdeltaco;IMdeltaco;' &
                  'zhh;zdr;kdp;rhohv;Ah;Av;' &
                  'sighhR;sigvvR;REdeltacoR;IMdeltacoR;' &
                  'zhhR;zdrR;kdpR;rhohvR;AhR;AvR')
ENDIF

!WRITE(0,*) '================== LAM=',LAM
LAM=LAM*1E-3

!================== Loop over temperature (°C) ====================
IF (testMode) THEN
  nTcloop=0
ELSE
  nTcloop=nTc-1
ENDIF

DO idTc=0,nTcloop
  Tc=Tcmin+idTc*Tcstep
  Tk=Tc+273.15 !T (°K) 
  
  !---------- Radar Constant     
  FREQ=(LIGHTSPEED)/LAM
  call QEPSW(Tk,FREQ,EPSW)
  QMW=SQRT(EPSW)  
  K2=((QMW**2-1.)/(QMW**2+2.))**2
  radar_cnst = (1E12*(LAM)**4)/(K2*(P**5D0)) !


  !=============== Loop over ELEV ===================
  IF (testMode) THEN
    nELEVloop=0
  ELSE
    nELEVloop=nELEV-1
  ENDIF
  DO idELEV=0,nELEVloop
    ELEV=ELEVmin+idELEV*ELEVstep

    !============= Loop over the 3d parameter of the Tmatrix table : P3 ==== !ChangeEB_Pending - I'm not sure I'm understanding this... should I say M=2 for Thom? Even if only it is 2M for rain+ice?
    !WRITE(0,*) 'Nmoments=',Nmoments 
    IF (Nmoments==2) THEN  ! if (LIMA + rain or cloud water) or ice : P3=concentration CC
      nP3=nint((expCCmax-expCCmin)/expCCstep)+1
      P3max=expCCmax 
      P3min=expCCmin 
      P3step=expCCstep 
      nFw=1
    ELSE ! 1-moment : P3 = Fw
      nP3=nint((Fwmax-Fwmin)/Fwstep)+1
      nFw=nP3
      P3max=Fwmax
      P3min=Fwmin
      P3step=Fwstep
    ENDIF
    
    nP3loop=nP3-1  
    IF (testMode) THEN
      WRITE(0,*) 'nP3loop=',nP3loop," nP3=",nP3," expCCmax=",expCCmax," expCCmin=",expCCmin
    !  nP3loop=0
    ENDIF
    
    DO idP3=0,nP3loop
      IF (Nmoments==2) THEN
        expP3=P3min+idP3*P3step
        P3=10**(expP3)
        idFw=0
      ELSE
        idFw=idP3
        P3=P3min+idP3*P3step
      ENDIF
      Fw=Fwmin+idFw*Fwstep
      
      !WRITE(0,*) '------ Fw=',Fw,'P3=',P3

      !=============== Loop over hydromet content M ============ !ChangeEB_Pending: So the hydreomet content is not read? Is it a table for all the possibilities?
                                                                 ! Should I then do the same for the air density / air denisty at the lowest level??
      nM=nint((expMmax-expMmin)/expMstep)+1
      nMloop=nM-1
      IF (testMode) THEN
        nMloop=0
      ENDIF

      DO idM=0,nMloop
        expM=expMmin+idM*expMstep
        M=10**(expM)         
        sighhint=0
        sigvvint=0
        kdpint=0
        Ahint=0
        Avint=0
        sighhRint=0
        sigvvRint=0
        kdpRint=0
        AhRint=0
        AvRint=0
        deltacoint=complex(0,0)
        deltacoRint=complex(0,0)          
        M_liq=0
        M_ss=0
        M_rr=0    
        M_ss2=0
        M_rr2=0

        ! ----------------------------------------------- !
        ! ------------------ Step 0 -------------------- !
        ! ----------------------------------------------- !

        ! Calculation of the "liquid" and "solid parts of the contents (M_rr
        ! and M_ss) with a Riemman sum (integration over diameters).
        ! => to evaluate the error (delta) resulting in doing this sum instead of 
        ! a complete integral (over the full diameter range)
        ! M_integral = M_Riemman + delta
        
        IF (testMode) THEN
          WRITE(0,*) '---- Step 0, M=',M
        ENDIF
        
                  
        Mint=0
        Mint2=0
        !============= Loop over diameters ==============
        DO idD=0,nD-1
               
          kTmat=idTc*nELEV*nFw*nD+idELEV*nFw*nD+idFw*nD+idD+1
          
          Tcrec=Tctab(kTmat)
          ELEVrec=ELEVtab(kTmat)
          Fwrec=Fwtab(kTmat)
          Drec=Dtab(kTmat)
          Dmrec=Dmtab(kTmat)
          Deqrec=Deqtab(kTmat)
          Deqmrec=Deqmtab(kTmat)
          Deqrrec=Deqrtab(kTmat)
          Deqrmrec=Deqrmtab(kTmat)
          
          sighh=sighhtab(kTmat)
          sigvv=sigvvtab(kTmat)
          kdp=kdptab(kTmat)
          Ah=Ahtab(kTmat)
          Av=Avtab(ktmat)
          REdeltaco=REdeltacotab(kTmat)
          IMdeltaco=IMdeltacotab(kTmat)
          
          sighhR=sighhRtab(kTmat)
          sigvvR=sigvvRtab(kTmat)
          kdpR=kdpRtab(kTmat)
          AhR=AhRtab(kTmat)
          AvR=AvRtab(kTmat)    
          REdeltacoR=REdeltacoRtab(kTmat)
          IMdeltacoR=IMdeltacoRtab(kTmat)                 
          
                    
          IF(idD .GT. 0) THEN 
            Drecinf=Dtab(kTmat-1)
            Dmrecinf=Dmtab(kTmat-1)
            Deqrecinf=Deqtab(kTmat-1)
            Deqmrecinf=Deqmtab(kTmat-1)
            Deqrrecinf=Deqrtab(kTmat-1)
            Deqrmrecinf=Deqrmtab(kTmat-1)
          ENDIF
          
          IF(idD .LT. nD-1) THEN
            Drecsup=Dtab(kTmat+1)
            Dmrecsup=Dmtab(kTmat+1)
            Deqrecsup=Deqtab(kTmat+1)
            Deqmrecsup=Deqmtab(kTmat+1)
            Deqrrecsup=Deqrtab(kTmat+1)
            Deqrmrecsup=Deqrmtab(kTmat+1)
          ENDIF
          
          Drec=Drec*1E-3 ! diameter in t-matrix table
          Dmrec=Dmrec*1E-3   ! max diameter
          Deqrrec=Deqrrec*1E-3 ! equivalent melted diameter  
          Deqrmrec=Deqrmrec*1E-3 ! equivalent melted max diameter 
          Deqmrec=Deqmrec*1E-3 ! 
          
          Deqmrecinf=Deqmrecinf*1E-3
          Deqmrecsup=Deqmrecsup*1E-3
          
          Dmrecinf=Dmrecinf*1E-3
          Dmrecsup=Dmrecsup*1E-3
          Deqrmrecinf=Deqrmrecinf*1E-3
          Deqrmrecsup=Deqrmrecsup*1E-3

          
          IF(idD .EQ. 0) THEN
            Dmrecinf=0 !Dmrec-Dstep
            Deqrmrecinf=0 !Deqrmrec-Dstep
            Deqmrecinf=0 
          ENDIF
          
          IF(idD .EQ. nD-1) THEN
            Dmrecsup=Dmrec+(Dmrec-Dmrecinf)
            Deqrmrecsup=Deqrmrec+(Deqrmrec-Deqrmrecinf)
            Deqmrecsup=Deqmrec+(Deqmrec-Deqmrecinf)
          ENDIF
          
          ! Distribution of the solid part of the hydrometeor
          IF (testMode) THEN
             WRITE(0,*) "----"
             WRITE(0,*) "Deqrmrec,Dmrec,Nmoments=",Deqrmrec,Dmrec,Nmoments
          ENDIF 
          IF((1-Fw)*M .GT. 0) THEN
            CALL PSD((1-Fw)*M,Dmrec,Tk,CCLOUD,P3,typeh,&
              aj,bj,nuj,alphaj,Cj,Xj,mumax,Nmoments,lamj,N_ss)
              !WRITE(0,*) "PSD solid part, N_ss=",N_ss
          ELSE
            N_ss=0
          ENDIF
          
          ! Distribution of the "liquid" part of the hydrometeor (if wet)
          ! The maximum diameter Dm is replaced by the maximum
          ! diameter of the equivalent "melted" hydrometeor Deqrm
          IF(Fw*M .GT. 0) THEN
            call PSD(Fw*M,Deqrmrec,Tk,CCLOUD,P3,'rr',&
            aj_rr,bj_rr,nuj_rr,alphaj_rr,Cj_rr,Xj_rr,mumax,Nmoments,lamj,N_rr)
            !WRITE(0,*) "PSD liquid part, N_rr=",N_ss
          ELSE
            N_rr=0
          ENDIF

          M_ss=M_ss+1000*(P/6)*(Deqrrec**3)*N_ss*(Dmrecsup-Dmrecinf)/2
          M_rr=M_rr+1000*(P/6)*(Deqrrec**3)*N_rr*(Deqrmrecsup-Deqrmrecinf)/2                       
          M_ss2=M_ss2+(aj*Dmrec**bj)*N_ss*(Dmrecsup-Dmrecinf)/2      
          M_rr2=M_rr2+(aj_rr*Deqrmrec**bj_rr)*N_rr*(Deqrmrecsup-Deqrmrecinf)/2          

          if (testMode) WRITE(0,*) "N_ss,N_rr = ",N_ss,N_rr          
                    
        ENDDO !===== End loop over diameters step 0
        
        IF (testMode) THEN
          WRITE(0,*) " ===== END STEP 0 : ================================="
          WRITE(0,*) "Fw, M = ",Fw,M
          WRITE(0,*) "M_ss,M_ss2 = ",M_ss,M_ss2           
          WRITE(0,*) "M_rr,M_rr2 = ",M_rr,M_rr2
          
          WRITE(0,*) "M, M_ss+M_rr : ",M, M_ss+M_rr
          WRITE(0,*) "M, M_ss2+M_rr2 : ", M, M_ss2+M_rr2
        ENDIF
        ! ----------------------------------------------- !
        ! ---------------- End Step 0 ------------------ !
        ! ----------------------------------------------- !

        ! ----------------------------------------------- !
        ! ------------------ Step 1 -------------------- !
        ! ----------------------------------------------- !
        ! Calculation of the liquid N_rr and solid N_ss parts of the distribution
        ! corrected by (1-Fw)*M/M_ss for N_ss and (Fw*M)/M_rr for N_rr
        ! to take into account the troncation of the distribution when doing a Riemann sum
        ! Calculation of the "combined" distribution N using Szyrmer and Zawadzki (1999) formulation
        ! (Re) Calculation of the content Mint (after) 
        
        !============= Loop over diameters step 1==============
        DO idD=0,nD-1
  
          kTmat=idTc*nELEV*nFw*nD+idELEV*nFw*nD+idFw*nD+idD+1

          Tcrec=Tctab(kTmat)
          ELEVrec=ELEVtab(kTmat)
          Fwrec=Fwtab(kTmat)
          Drec=Dtab(kTmat)
          Dmrec=Dmtab(kTmat)
          Deqrec=Deqtab(kTmat)
          Deqmrec=Deqmtab(kTmat)
          Deqrrec=Deqrtab(kTmat)
          Deqrmrec=Deqrmtab(kTmat)
          
          sighh=sighhtab(kTmat)
          sigvv=sigvvtab(kTmat)
          kdp=kdptab(kTmat)
          Ah=Ahtab(kTmat)
          Av=Avtab(ktmat)
          REdeltaco=REdeltacotab(kTmat)
          IMdeltaco=IMdeltacotab(kTmat)            
          sighhR=sighhRtab(kTmat)
          sigvvR=sigvvRtab(kTmat)
          kdpR=kdpRtab(kTmat)
          AhR=AhRtab(kTmat)
          AvR=AvRtab(ktmat)
          REdeltacoR=REdeltacoRtab(kTmat)
          IMdeltacoR=IMdeltacoRtab(kTmat)        
          
          IF(idD .GT. 0) THEN 
            Drecinf=Dtab(kTmat-1)
            Dmrecinf=Dmtab(kTmat-1)
            Deqrecinf=Deqtab(kTmat-1)
            Deqmrecinf=Deqmtab(kTmat-1)
            Deqrrecinf=Deqrtab(kTmat-1)
            Deqrmrecinf=Deqrmtab(kTmat-1)
          ENDIF
          
          IF(idD .LT. nD-1) THEN
            Drecsup=Dtab(kTmat+1)
            Dmrecsup=Dmtab(kTmat+1)
            Deqrecsup=Deqtab(kTmat+1)
            Deqmrecsup=Deqmtab(kTmat+1)
            Deqrrecsup=Deqrtab(kTmat+1)
            Deqrmrecsup=Deqrmtab(kTmat+1)
          ENDIF
  
          Drec=Drec*1E-3
          Dmrec=Dmrec*1E-3   ! Max diameter
          Deqrrec=Deqrrec*1E-3
          Deqrmrec=Deqrmrec*1E-3
          Deqmrec=Deqmrec*1E-3
          
          Dmrecinf=Dmrecinf*1E-3
          Dmrecsup=Dmrecsup*1E-3
          Deqrmrecinf=Deqrmrecinf*1E-3
          Deqrmrecsup=Deqrmrecsup*1E-3
          Deqmrecinf=Deqmrecinf*1E-3
          Deqmrecsup=Deqmrecsup*1E-3
  
          IF(idD .EQ. 0) THEN
            Dmrecinf=0 !Dmrec-Dstep
            Deqrmrecinf=0 !Deqrmrec-Dstep
            Deqmrecinf=0 
          ENDIF
          
          IF(idD .EQ. nD-1) THEN
            Dmrecsup=Dmrec+(Dmrec-Dmrecinf)
            Deqrmrecsup=Deqrmrec+(Deqrmrec-Deqrmrecinf)
            Deqmrecsup=Deqmrec+(Deqmrec-Deqmrecinf)
          ENDIF            
  
          !ChangeEB_Pending + ChangeEB_FallVelocity: Here I should add a new equation for the fall velocity in the WRF-Thompson scheme (Equation A3 from Thompson (2008))
          !In this equation apart from ccj and ddj there's a new term, f, for each hydrometeor. It is "Cj" in the configuration table (used in ICE3 MP)
          !f_cl:-; f_rr:195; f_ii:0; f_ss:100; f_gg:0
          !The equation should be: vtr = (RHOAS/RHO)**(1/2)*ccj*D**ddj*exp(Cj*D)
          !CAUTION: It needs the air density at the surface (RHOAS) and the air denisty (RHO). I still don't know how to call them - loop? 
          !CAUTION: Cloud water dos not sediement
          !Note: I have doubts about that part of the code. vtr - for rain, vts - for snow, vtm - for mixed? What happens with graupel? And cloud ice? 

          !In the Thompson MP scheme, vtx = (rho0/rho)**(1/2)*a*D**b*exp(-fD) (Thompson, 2008)
          !However, as the air density rho is not avaiable, the median of (rho0/rho)**(1/2) for a Study Case in la Cerdanya valley (Pyrinees) has been computed (=1.06778848)
          !User may want to adapt it. Maybe create a linear regression between temperature and rho?
          
          !ChangeEB_FallVelocity -  For test purposes, thom_fall experiment, use ICE3 param but with the fall velocity from thompson. 
          phi=0.246*Fw+(1-0.246)*(Fw**7)
          vtr = 1.06778848 *ccj_rr*Deqrrec**ddj*EXP(-Deqrrec*Cj)
          vts = 1.06778848 *ccj*Drec**ddj*EXP(-Drec*Cj) 
          vtm=phi*vtr+(1-phi)*vts 


!           IF (testMode) THEN
!               WRITE(0,*) " Drec,Deqrrec,Dmrec,Deqrmrec",Drec,Deqrrec,Dmrec,Deqrmrec
!           ENDIF
  
          ! Distribution of the solid part of the hydromet content
          IF((1-Fw)*M .GT. 0) THEN
            call PSD((1-Fw)*M,Dmrec,Tk,CCLOUD,P3,typeh,&
              aj,bj,nuj,alphaj,Cj,Xj,mumax,Nmoments,lamj,N_ss)
              N_ss2=N_ss
            IF(M_ss .GT. 0) THEN
              N_ss=(((1-Fw)*M)/M_ss)*N_ss
            ENDIF
            IF(M_ss2 .GT. 0) THEN
              N_ss2=(((1-Fw)*M)/M_ss2)*N_ss2
            ENDIF
          ELSE
            N_ss=0
            N_ss2=0
          ENDIF    
  
          ! Distribution of the liquid part of the hydromet content
          IF(Fw*M .GT. 0) THEN
            call PSD(Fw*M,Deqrmrec,Tk,CCLOUD,P3,'rr',&
                      aj_rr,bj_rr,nuj_rr,alphaj_rr,Cj_rr,Xj_rr,mumax,Nmoments,&
                      lamj,N_rr)
            N_rr2=N_rr
            IF(M_rr .GT. 0) THEN
              N_rr=((Fw*M)/M_rr)*N_rr
            ENDIF
            IF(M_rr2 .GT. 0) THEN
              N_rr2=((Fw*M)/M_rr2)*N_rr2
            ENDIF
          ELSE
            N_rr=0
            N_rr2=0
          ENDIF
          
          ! Combination of the liquid and solid parts of the distribution
          ! Szyrmer and Zawadzki (1999)
          ! PhD T. Le Bastard (eq 4.20 p 107)
          N = (1-Fw)*N_ss*(vts/vtm)+Fw*N_rr*(vtr/vtm)
          N2 = (1-Fw)*N_ss2*(vts/vtm)+Fw*N_rr2*(vtr/vtm)
      
          Mint=Mint+1000*(P/6)*(Deqrrec**3)*N*(Deqrmrecsup-Deqrmrecinf)/2
          Mint2=Mint2+(aj_rr*Deqrmrec**bj_rr)*N2*(Deqmrecsup-Deqmrecinf)/2
          
!           IF (testMode) THEN
!               WRITE(0,*) " N_ss,N_ss2,N_rr,N_rr2, N, N2",N_ss, N_ss2, N_rr, N_rr2,N,N2
!           ENDIF


        ENDDO !===== End loop over diameters step 1
        
        IF (testMode) THEN
          WRITE(0,*) "======= END STEP 1 ====="
          WRITE(0,*) "M, Mint, Mint2 : ",M, Mint, Mint2
        ENDIF

        ! ----------------------------------------------- !
        ! ---------------- End step 1 ------------------ !
        ! ----------------------------------------------- !


        ! ----------------------------------------------- !
        ! ------------------ Step 2 -------------------- !
        ! ----------------------------------------------- !
      
        ! The total number concentration is recomputed to apply the correcting factor M/Mint.
        RRint=0
        RRint2=0
        !============= Loop over diameters  step 2==============
        DO idD=0,nD-1
    
          !Position dans le fichier tmat des coefs de diffusion
          kTmat=idTc*nELEV*nFw*nD+idELEV*nFw*nD+idFw*nD+idD+1

          Tcrec=Tctab(kTmat)
          ELEVrec=ELEVtab(kTmat)
          Fwrec=Fwtab(kTmat)
          Drec=Dtab(kTmat)
          Dmrec=Dmtab(kTmat)
          Deqrec=Deqtab(kTmat)
          Deqmrec=Deqmtab(kTmat)
          Deqrrec=Deqrtab(kTmat)
          Deqrmrec=Deqrmtab(kTmat)
          sighh=sighhtab(kTmat)
          sigvv=sigvvtab(kTmat)
          kdp=kdptab(kTmat)
          Ah=Ahtab(kTmat)
          Av=Avtab(ktmat)
          REdeltaco=REdeltacotab(kTmat)
          IMdeltaco=IMdeltacotab(kTmat)            
          sighhR=sighhRtab(kTmat)
          sigvvR=sigvvRtab(kTmat)
          kdpR=kdpRtab(kTmat)
          AhR=AhRtab(kTmat)
          AvR=AvRtab(ktmat) 
          REdeltacoR=REdeltacoRtab(kTmat)
          IMdeltacoR=IMdeltacoRtab(kTmat)              
          
          IF(idD .GT. 0) THEN 
            Drecinf=Dtab(kTmat-1)
            Dmrecinf=Dmtab(kTmat-1)
            Deqrecinf=Deqtab(kTmat-1)
            Deqmrecinf=Deqmtab(kTmat-1)
            Deqrrecinf=Deqrtab(kTmat-1)
            Deqrmrecinf=Deqrmtab(kTmat-1)
          ENDIF
          
          IF(idD .LT. nD-1) THEN
            Drecsup=Dtab(kTmat+1)
            Dmrecsup=Dmtab(kTmat+1)
            Deqrecsup=Deqtab(kTmat+1)
            Deqmrecsup=Deqmtab(kTmat+1)
            Deqrrecsup=Deqrtab(kTmat+1)
            Deqrmrecsup=Deqrmtab(kTmat+1)
          ENDIF

          Drec=Drec*1E-3
          Dmrec=Dmrec*1E-3   ! Diamètre max de la particule pure équivalente
          Deqrrec=Deqrrec*1E-3
          Deqrmrec=Deqrmrec*1E-3
          Deqmrec=Deqmrec*1E-3
          
          Dmrecinf=Dmrecinf*1E-3
          Dmrecsup=Dmrecsup*1E-3
          Deqrmrecinf=Deqrmrecinf*1E-3
          Deqrmrecsup=Deqrmrecsup*1E-3
          Deqmrecinf=Deqmrecinf*1E-3
          Deqmrecsup=Deqmrecsup*1E-3


          IF(idD .EQ. 0) THEN
            Dmrecinf=0 !Dmrec-Dstep
            Deqrmrecinf=0 !Deqrmrec-Dstep
            Deqmrecinf=0 
          ENDIF
          
          IF(idD .EQ. nD-1) THEN
            Dmrecsup=Dmrec+(Dmrec-Dmrecinf)
            Deqrmrecsup=Deqrmrec+(Deqrmrec-Deqrmrecinf)
            Deqmrecsup=Deqmrec+(Deqmrec-Deqmrecinf)
          ENDIF
          

          IF ((typeh .EQ. 'wg') .OR. (typeh .EQ. 'ws')) THEN
            ! Vitesses de chute
            phi=0.246*Fw+(1-0.246)*(Fw**7)
            vtr=ccj_rr*(Deqrrec**ddj_rr)
            vts=ccj*(Drec**ddj)
            vtm=phi*vtr+(1-phi)*vts

            ! Distribution of the solid part of the hydromet content
            IF((1-Fw)*M .GT. 0) THEN
              call PSD((1-Fw)*M,Dmrec,Tk,CCLOUD,P3,typeh,&
                        aj,bj,nuj,alphaj,Cj,Xj,mumax,Nmoments,lamj,N_ss)
              N_ss2=N_ss
              IF(M_ss .GT. 0) THEN
                N_ss=(((1-Fw)*M)/M_ss)*N_ss
                N_ss2=(((1-Fw)*M)/M_ss2)*N_ss2
              ENDIF
            ELSE
              N_ss=0
              N_ss2=0
            ENDIF

            ! Distribution of the liquid part of the hydromet content
            IF(Fw*M .GT. 0) THEN
              call PSD(Fw*M,Deqrmrec,Tk,CCLOUD,P3,'rr',&
                        aj_rr,bj_rr,nuj_rr,alphaj_rr,Cj_rr,Xj_rr,mumax,Nmoments,lamj,N_rr)
              N_rr2=N_rr
              IF(M_rr .GT. 0) THEN
                N_rr=((Fw*M)/M_rr)*N_rr
                N_rr2=((Fw*M)/M_rr2)*N_rr2
              ENDIF
            ELSE
              N_rr=0
              N_rr2=0
            ENDIF

            ! Combination of the liquid and solid parts of the distribution
            N = (1-Fw)*N_ss*(vts/vtm)+Fw*N_rr*(vtr/vtm)
            N2 = (1-Fw)*N_ss2*(vts/vtm)+Fw*N_rr2*(vtr/vtm)
          ELSE
            !call PSD(M,(Dmrec+Dmrecinf)/2.,Tk,&
            ! aj,bj,nuj,alphaj,Cj,Xj,lamj,N)
            vtm=ccj*(Drec**ddj)
            call PSD(M,Dmrec,Tk,CCLOUD,P3,typeh,&
                      aj,bj,nuj,alphaj,Cj,Xj,mumax,Nmoments,lamj,N)
            N2=N  
          ENDIF ! ENDIF "ws" or "wg"
          
          IF(Mint .GT. 0) THEN
            N=(M/Mint)*N
            N2=(M/Mint2)*N2
          ENDIF
          
          RRint = RRint+1000*(P/6)*(Deqrrec**3)*N*vtm*3600*(Deqrmrecsup-Deqrmrecinf)/2            
          RRint2 = RRint2+1000*(P/6)*(Deqrrec**3)*N2*vtm*3600*(Deqmrecsup-Deqmrecinf)/2            
          
          sighhint=sighhint+sighh*N2*(Deqmrecsup-Deqmrecinf)/2          
          sigvvint=sigvvint+sigvv*N2*(Deqmrecsup-Deqmrecinf)/2         
          kdpint=kdpint+kdp*N2*(Deqmrecsup-Deqmrecinf)/2
          Ahint=Ahint+Ah*N2*(Deqmrecsup-Deqmrecinf)/2
          Avint=Avint+Av*N2*(Deqmrecsup-Deqmrecinf)/2
          sighhRint=sighhRint+sighhR*N2*(Deqmrecsup-Deqmrecinf)/2
          sigvvRint=sigvvRint+sigvvR*N2*(Deqmrecsup-Deqmrecinf)/2
          kdpRint=kdpRint+kdpR*N2*(Deqmrecsup-Deqmrecinf)/2
          AhRint=AhRint+AhR*N2*(Deqmrecsup-Deqmrecinf)/2
          AvRint=AvRint+AvR*N2*(Deqmrecsup-Deqmrecinf)/2
          deltaco=cmplx(REdeltaco,IMdeltaco)
          deltacoint=deltacoint+deltaco*N2*(Deqmrecsup-Deqmrecinf)/2
          deltacoR=cmplx(REdeltacoR,IMdeltacoR)
          deltacoRint=deltacoRint+deltacoR*N2*(Deqmrecsup-Deqmrecinf)/2           
  
          !WRITE(0,*) Fw,N*((Dmrecsup+Dmrec)/2-(Dmrec+Dmrecinf)/2),M,M_liq,Mint,Mint2,M/Mint,M/Mint2

        ENDDO !===== End loop over diameters step 2

        ! ----------------------------------------------- !
        ! ---------------- End step 2 ------------------ !
        ! ----------------------------------------------- !

!          !WRITE(0,*) "-- 4 --",M,Mint2,M_liq,RRint
!          IF (testMode) THEN
!            WRITE(0,*) "======= END STEP 2 ====="
!            WRITE(0,*) "M, Mint, Mint2 : ",M, Mint, Mint2
!          ENDIF

        IF (isnan(sighhint)) THEN
          sighhint=0
          sigvvint=0
          kdpint=0
          Ahint=0
          Avint=0
          sighhRint=0
          sigvvRint=0
          kdpRint=0
          AhRint=0
          AvRint=0
          REdeltacoint=0
          IMdeltacoint=0
          REdeltacoRint=0
          IMdeltacoRint=0
        ENDIF
        
        
        ! ----- Dual-pol variables with T-matrix
        CALL COMPUTE_VARPOL(sighhint,sigvvint,radar_cnst,P,deltacoint,zhhint,& 
                            zdrint,rhohvint)
                                    
        ! ----- Dual-pol variables with Rayleigh
        CALL COMPUTE_VARPOL(sighhRint,sigvvRint,radar_cnst,P,deltacoRint,& 
                            zhhRint,zdrRint,rhohvRint)

        if (.not. testMode) THEN
            !-----  Writing PSD integrated coefficients + variables           
            WRITE (6,5006) Tc,ELEV,P3,M,RRint,sighhint,sigvvint,REAL(deltacoint),&
            IMAG(deltacoint),zhhint,zdrint,kdpint,rhohvint,Ahint,Avint,& 
            sighhRint,sigvvRint,REAL(deltacoRint),IMAG(deltacoRint),& 
            zhhRint,zdrRint,kdpRint,rhohvRint,AhRint,AvRint
    
            5006 FORMAT (F5.1,';',F4.1,';',E9.3,';',E10.4,';',E10.4,';',E10.4,';',&
                          E10.4,';',E10.4,';',E10.4,';',E10.4,';',E10.4,';',&
                          E10.4,';',E10.4,';',E10.4,';',E10.4,';',E10.4,';',&  
                          E10.4,';',E10.4,';',E10.4,';',E10.4,';',E10.4,';',&
                          E10.4,';',E10.4,';',E10.4,';',E10.4)
        ENDIF                 


        IF (testMode) THEN
            WRITE (0,*) "M=",M
            WRITE (0,*) "P3=",P3
            WRITE (0,*) "zhhint=",zhhint
            WRITE (0,*) "zdrint=",zdrint
            WRITE (0,*) "kdpint=",kdpint
            WRITE (0,*) "Ahint=",Ahint
            WRITE (0,*) "Avint=",Avint
            WRITE (0,*) "RRint, RRint2=",RRint, RRint2
        ENDIF

      ENDDO !fin boucle sur M
    ENDDO !fin boucle P3 (Fw or number concentration No)
  ENDDO !fin boucle ELEV
ENDDO !fin boucle T°  

CLOSE(6)

DEALLOCATE (Tctab,ELEVtab,&
          Fwtab,Dtab,Dmtab,Deqtab,&
          Deqmtab,Deqrtab,Deqrmtab,&
          sighhtab,sigvvtab,REdeltacotab,&
          IMdeltacotab,zhhtab,zdrtab,&
          kdptab,rhohvtab,Ahtab,Avtab,sighhRtab,&
          sigvvRtab,REdeltacoRtab,&
          IMdeltacoRtab,zhhRtab,zdrRtab,& 
          kdpRtab,rhohvRtab,AhRtab,AvRtab)

WRITE(0,*) ' SAVED !'
WRITE(0,*) ' -----------------------------------------------------------'

END PROGRAM MakeTmatInt

! ========================================================================
SUBROUTINE COMPUTE_VARPOL(sighh,sigvv,radar_cnst,P,deltaco,zhh,zdr,rhohv)
! computation of Zh (dBZ), Zdr (dB) and rhohv from scattering cross sections

IMPLICIT NONE

REAL,INTENT(in) :: sighh,sigvv,radar_cnst,P
COMPLEX,INTENT(in) :: deltaco
REAL,INTENT(out) :: zhh,zdr,rhohv

REAL ::zvv


IF (sighh>0) THEN
    zhh=10*log10(radar_cnst*sighh)
ELSE
    zhh=-999.0
ENDIF
IF (sigvv>0) THEN
    zvv=10*log10(radar_cnst*sigvv)
ELSE
    zvv=-999.0
ENDIF
IF (sighh>0 .AND. sigvv>0) THEN
    zdr=zhh-zvv
    rhohv = 4*P*ABS(deltaco)/SQRT(sighh*sigvv)
ELSE
    zdr=-999.0
    rhohv=-999.0    
ENDIF

RETURN
END
! END COMPUTE_VARPOL

SUBROUTINE QEPSW(PTEMP,PFREQ,EPSW) 
! water complex dielectric function (Liebe et al., 1991)
! electromagnetic fields in exp(-i*omega*t), i.e. Im(epsw)>=0
! in  : ptemp=temperature in K
!       pfreq=frequency in Hz
! out : eps=epsilon
IMPLICIT NONE
REAL,INTENT(in) :: PFREQ,PTEMP
COMPLEX,INTENT(out) :: EPSW

REAL :: ZTHETA,ZEZ,ZEINF,ZF

ZTHETA=1-300./PTEMP
ZEZ=77.66-103.3*ZTHETA
ZEINF=0.066*ZEZ
ZF=(20.27+146.5*ZTHETA+314.*ZTHETA**2)*1.E9

EPSW=ZEINF+(ZEZ-ZEINF)/(1.-(0.,1.)*PFREQ/ZF)
RETURN
END
!END FUNCTION QEPSW


!**********************************************************************
!            SUBROUTINE PSD 
!
! Calculation of slope parameter and PSD for a given M, D
! and concentration if Nmoments=2)
! snow parameterization based on Wurtz 2021, 2023
!
! IN: M (Hydrometeor content if 1 moment (kg m-3) 
! M=rhodref * r  with rhodref = dry air density (kg m-3)
!                           r = mixing ratio (kg kg-1)
! PhD O. Caumont (2007) p 23
! D = diameter (m)
! CCLOUD (Microphysics scheme: ICE3 or LIMA)
! P3 : Liquid water fraction Fw or number concentration No=CRAINT*rhodref (m-3)
! a [kg.m^-b],b,nu,alpha,c,x : microphysics parameters
!
! OUT: lamb,N
! 1-moment and 2-moments: N (m-3)

!**********************************************************************          
SUBROUTINE PSD(M,D,Tk,CCLOUD,P3,typeh,a,b,nuconst,alpha,c,x,mumax,Nmoments,lamb,N)
IMPLICIT NONE
INTEGER,INTENT(in) :: Nmoments
REAL,INTENT(in)  :: a,b,nuconst,alpha,c,x,mumax
REAL,INTENT(in)  :: D,M,P3,Tk ! temperature in Kelvin for Wurtz parameterization of snow
CHARACTER*2,INTENT(in)  ::typeh
CHARACTER*4,INTENT(in)  ::CCLOUD
REAL*8,INTENT(out) :: lamb,N

REAL*8 :: No,mucalc,mu,Dm,nu,lambMP

!ChangeEB_PSD - Parameters for Thompson Snow PSD
REAL*8 :: M0, M2, M3, oM3, Mrat, loga_, a_, b_, tc_ndth, slam1, slam2
REAL, PARAMETER:: mu_s = 0.6357
REAL, PARAMETER:: Kap0 = 490.6
REAL, PARAMETER:: Kap1 = 17.46
REAL, PARAMETER:: Lam0 = 20.78
REAL, PARAMETER:: Lam1 = 3.29
REAL, PARAMETER:: PI = 3.1415926536
REAL, PARAMETER:: Nt_c = 100E6 !Units m-3 = 100cm-3

REAL, DIMENSION(10), PARAMETER:: &
sa = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
        0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
REAL, DIMENSION(10), PARAMETER:: &
sb = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
        0.060366,  0.000079,  0.000594, 0.0,      -0.003577/)

!Parameters for Thompson Graupel PSD
REAL*8 :: ygra1, zans1, N0_exp, lam_exp, lam_prev, ng

!ChangeEB_PSD - Parameters c*e, c*g, o*g:
REAL*8 :: cre2, crg2, org2
REAL*8 :: cie2, cie3, cig2, cig3, oig2
REAL*8 :: cse1
REAL*8 :: cce1, cce2, cce3, ccg1, ccg2, ccg3
REAL*8 :: cge1, cge2, cge3, cgg1, cgg2, cgg3, oge1, ogg1, ogg2, ogg3

! LIMA is 2-moments for cloud droplets, rain drops and pristine ice crystals
! For these 3 species, the number concentration N is prognostic
! In this code, the unit of No is m-3
! => the pronostic parameter CRAINT (P3) in MesoNH has to be multiplied by rhodref
!  before looking for the scattering coef in the T-matrix table

! LIMA Marie Taufour with diagnostic mu  
! IF (CCLOUD=='LIMT' .AND. typeh=='rr'.AND. Nmoments==2) THEN
!   Dm=(M/(a*No))**(1/b)  ! Mean volume diameter for spherical drops (b = 3)
!                         ! obtained if we consider that all drops have the same size 
!                         ! (MesoNH 5.4 Physics doc eq 7.30 p 146)   
!   mucalc=38.66*exp(-1.59*(Dm*10**3)) ! adapted from eq (7) of Taufour et al 2018 (Dm in mm!)
!   mu=max(min(mucalc,mumax),0.0) ! mu is borned: 0 <= mu <= 15
!   nu=mu+1
! ELSE
nu=nuconst
!ENDIF

! ---- Compute No and lambda -----
IF (CCLOUD=='ICE3') THEN  
  IF (typeh=='ii') THEN
    No = P3
    lamb= (a*No*GAMMA(nu+b/alpha)/(M*GAMMA(nu)))**(1./b) !PSD LIMA
  ELSE
    No   = c*(lamb**x)
    lamb = (M*GAMMA(nu)/(a*c*GAMMA(nu+b/alpha)))**(1./(x-b)) ! PSD ICE3
  ENDIF

ELSE IF (CCLOUD=='ICJW') THEN
  IF (typeh=='ss') THEN 
    IF (Tk>263.15) THEN 
      lambMP = 10**(14.554-0.0423*Tk) ! eq6 Wurtz 2023
    ELSE
      lambMP = 10**(6.226-0.0106*Tk) ! eq6 Wurtz 2023
    ENDIF
    lamb = lambMP * (1/(2*SQRT(2.0))) &
          * SQRT((GAMMA(nu+2/alpha)*GAMMA(nu+4/alpha)) &
                  /(GAMMA(nu+1/alpha)*GAMMA(nu+3/alpha))) ! eq8 Wurtz 2023
    No = (M*GAMMA(nu)*lamb**b)/(a*GAMMA(nu+b/alpha)) ! eq III.7 these J. Wurtz
  ELSE
    IF (typeh=='ii') THEN
      No = P3
      lamb= (a*No*GAMMA(nu+b/alpha)/(M*GAMMA(nu)))**(1./b)
    ELSE
      lamb = (M*GAMMA(nu)/(a*c*GAMMA(nu+b/alpha)))**(1./(x-b))
      No   = c*(lamb**x)
    ENDIF
  ENDIF

ELSE IF (CCLOUD=='LIMA') THEN
  IF (Nmoments==1) THEN
    lamb = (M*GAMMA(nu)/(a*c*GAMMA(nu+b/alpha)))**(1./(x-b))
    No   = c*(lamb**x)
  ELSE IF (Nmoments==2) THEN
    No=P3
    lamb= (a*No*GAMMA(nu+b/alpha)/(M*GAMMA(nu)))**(1./b)
  ENDIF


ELSE IF (CCLOUD=='LIMC') THEN
  IF (typeh=='ss') THEN 
    IF (Tk>263.15) THEN 
      lambMP = 10**(14.554-0.0423*Tk) ! eq6 Wurtz 2023
    ELSE
      lambMP = 10**(6.226-0.0106*Tk) ! eq6 Wurtz 2023
    ENDIF
    lamb = lambMP * (1/(2*SQRT(2.0))) &
          * SQRT((GAMMA(nu+2/alpha)*GAMMA(nu+4/alpha)) &
                  /(GAMMA(nu+1/alpha)*GAMMA(nu+3/alpha))) ! eq8 Wurtz 2023
    No = (M*GAMMA(nu)*lamb**b)/(a*GAMMA(nu+b/alpha)) ! eq III.7 these J. Wurtz
  ELSE
    IF (Nmoments==1) THEN
      IF (typeh=='ii') THEN
        No = P3
        lamb= (a*No*GAMMA(nu+b/alpha)/(M*GAMMA(nu)))**(1./b)
      ELSE
        lamb = (M*GAMMA(nu)/(a*c*GAMMA(nu+b/alpha)))**(1./(x-b))
        No   = c*(lamb**x)
      ENDIF
    ELSE IF (Nmoments==2) THEN
      No=P3
      lamb= (a*No*GAMMA(nu+b/alpha)/(M*GAMMA(nu)))**(1./b)
    ENDIF
  ENDIF

ELSE IF (CCLOUD=='THOM') THEN !ChangeEB_PSD - Add the Thompson Microphysical scheme
  IF (typeh=='rr') THEN
    cre2 = nu+1
    crg2 = GAMMA(cre2)
    org2 = 1/crg2
    lamb = (PI*1000/6 * P3/M * GAMMA(nu+3+1)/GAMMA(nu+1))**(1./b)
    No = P3*org2 * lamb**cre2

  ELSE IF (typeh=='ii') THEN
    cie2 = nu+1
    cie3 = nu+b+1
    cig2 = GAMMA(cie2)
    cig3 = GAMMA(cie3)
    oig2 = 1/cig2
    lamb = (PI*890/6 * P3/M * cig3/cig2)**(1./b)
    No = P3*oig2*lamb**cie2 

  ELSE IF (typeh=='cl') THEN
    cce1 = 1+nu
    cce2 = b+1+nu
    cce3 = nu+b+1+3
    ccg1 = GAMMA(cce1)
    ccg2 = GAMMA(cce2)
    ccg3 = GAMMA(cce3)
    lamb = (Nt_c*a*ccg2*1/ccg1/M)**(1./b)
    No = Nt_c/ccg1 * lamb**cce1
    

  ELSE IF (typeh=='ss') THEN
    lamb = 0 !In fact, snow in Thom does not havy any lambda/N0. Is lamba necessary for the code/operadar in the future?? 
    No = 0 

  ELSE IF (typeh=='gg') THEN
    cge1 = b+1
    cge2 = nu+1
    cge3 = b+nu+1
    cgg1 = GAMMA(cge1)
    cgg2 = GAMMA(cge2)
    cgg3 = GAMMA(cge3)
    oge1 = 1/cge1
    ogg1 = 1/cgg1
    ogg2 = 1/cgg2
    ogg3 = 1/cgg3

    ygra1 = alog10(MAX(1.E-9, M))
    zans1 = 3.0 + 2./7.*(ygra1+8.)
    zans1 = MAX(2., MIN(zans1, 6.))
    N0_exp = 10.**(zans1)
    lam_exp = (N0_exp*a*cgg1/M)**oge1
    lamb = lam_exp * (cgg3*ogg2*ogg1)**(1./b) !This lambda is used to predict the ng || !ChangeEB_Pending - Not sure about what lambda take... I could try with this
    !! v1 computes No from N0_exp. Is an approach similar to CR-SIM, although not exactly the same
    !No = N0_exp/(cgg2*lam_exp)*lamb**cge2 

    !!ChangeEB_Pending - v2 computes No from the diagnostic ng, 
    ng = cgg2*ogg3*M*lamb**b / a
    !!!lamb = (a*cgg3*ogg2*ng/M)**(1./b) !Recalculate lambda - Not in v0 
    No = ng*ogg2*lamb**cge2
  ENDIF
ENDIF    

  

! Compute N from N0 and lambda !ChangeEB_PSD: I commented the following line, to adapt the N(D) if it's thom_snow
!N = No*(alpha/GAMMA(nu))*(lamb**(alpha*nu))*(D**(alpha*nu-1))*EXP(-(lamb*D)**alpha)
!Extracted from Thompson MP (WRF model) - Thompson (2008)

IF (CCLOUD=='THOM') THEN
  IF (typeh=='ss') THEN
    !Add here the Thom Snow PSD. Only for Thompson and snow. 

    !Calculate the second moment (M2)
    M2 = M * 1/a 
    
    !Calculate the third moment (M3)

    cse1 = b+1
    tc_ndth = MIN(-0.1, Tk-273.15) !tc_ndths is the temperature in celsius - tc_NDTHompsonSnow
    !The WRF source code is a bit tricky there:
    loga_ = sa(1) + sa(2)*tc_ndth + sa(3)*cse1 &
        + sa(4)*tc_ndth*cse1 + sa(5)*tc_ndth*tc_ndth &
        + sa(6)*cse1*cse1 + sa(7)*tc_ndth*tc_ndth*cse1 &
        + sa(8)*tc_ndth*cse1*cse1 + sa(9)*tc_ndth*tc_ndth*tc_ndth &
        + sa(10)*cse1*cse1*cse1
    a_ = 10.0**loga_
    b_ = sb(1)+sb(2)*tc_ndth+sb(3)*cse1 + sb(4)*tc_ndth*cse1 &
        + sb(5)*tc_ndth*tc_ndth + sb(6)*cse1*cse1 &
        + sb(7)*tc_ndth*tc_ndth*cse1 + sb(8)*tc_ndth*cse1*cse1 &
        + sb(9)*tc_ndth*tc_ndth*tc_ndth+sb(10)*cse1*cse1*cse1
    M3 = a_ * M2**b_

    !Simplifications for the Snow PSD
    oM3 = 1./M3
    Mrat = M2*(M2*oM3)*(M2*oM3)*(M2*oM3)
    M0   = (M2*oM3)**mu_s
    slam1 = M2 * oM3 * Lam0
    slam2 = M2 * oM3 * Lam1

    !Line 4478 off module_mp_thompson.F of WRF source code - We also could look at line 2087...
    N = Mrat*(Kap0*DEXP(-slam1*D) + Kap1*M0*D**mu_s * DEXP(-slam2*D))
  ELSE
    N = No*D**nu*exp(-lamb*D) !Gamma function for rain, graupel, cloud and ice in WRF_Thompson
  ENDIF
ELSE
  N = No*(alpha/GAMMA(nu))*(lamb**(alpha*nu))*(D**(alpha*nu-1))*EXP(-(lamb*D)**alpha)
ENDIF

! !IF (testMode) THEN
!   WRITE (0,*) "lamb,D,alpha, nu, No,N=",lamb,D,alpha,nu,No,N
!   WRITE (0,*) "No*alpha/GAMMA(nu)",No*alpha/GAMMA(nu)
!   WRITE (0,*) "lamb**(alpha*nu)",lamb**(alpha*nu)
!   WRITE (0,*) "(D**(alpha*nu-1))",(D**(alpha*nu-1))
!   WRITE (0,*) "EXP(-(lamb*D)**alpha)",EXP(-(lamb*D)**alpha)
!   WRITE (0,*) "Ncalc=",No*(alpha/GAMMA(nu))*(lamb**(alpha*nu))*(D**(alpha*nu-1))*EXP(-(lamb*D)**alpha)
! !ENDIF      

RETURN
END
        









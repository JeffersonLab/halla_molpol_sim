c-------------------------------------------------------------------------------
c     Original Monte Carlo code of Swartz originally written for the
c     SLAC Moller polarimeter (NEWMOL); the corresponding paper:
c     M. Swartz, NIM Nucl. Instr. and Meth. A 363 (1995) 526
c
c     This is a modified version of NEWMOL for a two arm polarimeter
c     Include real momentum distributions for Fe+Co electrons
c     This version is modified for use with the CEBAF Moller Polarimeter
c     af. 10.08.94.
c
c     All subroutines for the CEBAF optics programmed by Laurens de Bever (L.B.)
c
ccc
ccc   Modifications to get the information available in ntuples for
ccc   reanalyzing with PAW, by Matthias Loppacher 12.94.
ccc   Comments and additional modifications, by Matthias Loppacher 11.96. (M.L.)
ccc   !Histo    specific for histo analysis
ccc   !Ntuple   specific for ntuple analysis
ccc
c
c-------------------------------------------------------------------------------
C
      PROGRAM TWOMOL
c
      IMPLICIT NONE
C
      character*128 argin(10)
      integer   nrarg, iargc, inunit

      real*8 quad1x, quad1y, quad3x, quad3y !offset of quads
      integer npads             !can be used for #hodoscope channels (set also in
      parameter (npads=1)       !routine Detectors)
      REAL*8 Lsol, LQ1, LQ2, LQ3, GdotL1,GdotL2, GdotL3
      REAL*8 ZDIST,VEC(6,2)
      real*8 quad1exit(2,2)
      real*8 quad2entrance(2,2),quad2exit(2,2)
      real*8 quad3entrance(2,2),quad3exit(2,2),q3mechexit(2,2)
      real*8 pipeexit(2,2),armentr(2,2),armexit(2,2)
      REAL*8 PEL,PBEAM,CST,TGEN,PHIGEN,WEIGHT,PWT(2,3)
      real*4 p_up
      real*8 pmoll(2)           !lab momentum of electrons
      real*8 mytheta,myphi
      REAL*8 XTGT,DEPTH,XDEPTH,PB
      REAL*8 WGTBM,WGTML
      real*8 WGT,ASX2,ASY2,ASZ2,wgtcx
      REAL*8 SIGNAL,SIGNL2,ASZ,ASX,ASY,DIFF,WEITOT,WEITO2
      REAL*8 siSIGNAL, siSIGNL2, DBsignal, DBsignl2
      real*8 WEITOTs(2),WEITO2s(2)
      REAL*8 BLKTL(2,3),BLKTL2(2,3)
      real*8 avgA,avgA2,davgA
      real*8 weights(2)
      real*8 TC, TM1, M1M2, M2M3,M1C, CM2, M3D(2)
      real*8 BSolLong, MinSolFld, SM1, soldx, soldy, soldxp, soldyp
      parameter (MinSolFld=0.1d0)
      real*4 Azz(npads,2), Azze(npads,2)
      real*4 DbAzz, DbAzze, Dbsig, Dbsige
      INTEGER I,J,K,L,NUMEVENTS,ITRIAL,IHD
      INTEGER solsim
      integer nplt, nplt2, nplt3, nplt4
      integer ienergy
      real*8 padtl(2,3,npads,2), padtl2(2,3,npads,2)
      INTEGER II, JJ
      LOGICAL TCKBRM,TAILS,BINDE,Q1CLIP,QFRINGE,THREEQUAD
      logical twomols, singlemol(2), scflag(3)
      CHARACTER*12 TITLE(2)
      CHARACTER*55 HEADING
      CHARACTER*80 FILENAME
      real*8 b_energy, mcur1, mcur2, r1, r2, r3, tip1, tip2, tip3
      real*8 grad1, grad2, grad3
      real*8 q1tip, q2tip, q3tip
      real*8 q1rad, rq1_1, rq1_2 !Q1 clearance, and electron 1 & 2 radii at Q1 exit
      data q1rad / 39.0/        !Q1 exit clearance (radius, in mm)

      common /collimators/ pcol,colxmin,colxmax,colymin,colymax
      real*8 pcol(7)
      real*8 colxmin(7),colxmax(7),colymin(7),colymax(7)

      real*8 col(2,7,2)
      real*8 col1z, col3z, col5z, col6z
      real*8  PPOL, PUNPOL, TGTPOL
      real*8 targetpol, twom, efficiency
      real*8  XQ1, XQ2, XQ3,XDET(2)
      real*8 ran3
      real*8 quad2offset, quad3offset
      real*8 polfact, polfacte
      real*8 dgpolfact,dgpolfacte
      logical all_coll          !function declaration
ccc
      INTEGER NWPAWC, H, iquest
      PARAMETER(NWPAWC= 5000000) ! OK for 1'100'000 events (6 effects on)
      COMMON /PAWC/ H(NWPAWC)
c     common/quest/iquest(100)          ! default 16000
      common/quest/iquest(16000) ! default 16000
      real*8 climit, crange
      common /anglerange/ climit,crange
c
      integer*4 istat,icycle    ! M.L. L.B. 94
c     real*4            ntuple(20)
c     common/cern/      ntuple
ccc
      COMMON / RANSEED / II

      REAL*8 SIGX,SIGXP,SIGY,SIGYP,X0,Y0,TX0,TY0
      common /beamparam/ SIGX,SIGY,SIGXP,SIGYP,X0,Y0,TX0,TY0

      character*1 answer
C

c     gradient parameters for quadrupoles
c     data Q1/-0.00438,-0.00578,-0.00608,-0.00572,-0.00502,-0.00408/
c     data Q2/0.001229,0.002621,0.004308,0.006221,0.008336,0.010653/

C
C     TCKBRM and TAILS switch on thick target bremsstrahlung and
C     multiple scattering with Moliere tails, respectively
C     BINDE switches on the electron binding effects
C
      DATA TITLE/' 10-UM MOLIE',' 10-UM GAUSS'/
C


      CALL HLIMIT(NWPAWC)
      iquest(10)=64000          ! 16000 -> 64000 events times 4
! has to be placed after hlimit!
      solsim = 0
      twom = 0.
C
C     this may need to be updated depending on field
      targetpol = 0.08043
c     targetpol = 1.0
C     quad magnet aperture (radius)
c     r1=47.625
      r1=47.6
      r2=127.
      r3=127.                   !3rd quad, same as 2nd
C
C     MAGNET LENGTH
C     Lsol= 125.    !This gets set later now
      LQ1 = 361.92              !250.
c     LQ2 = 965.2            ! this is 914.4+0.6*r2
      LQ2 = 990.6               ! this is 39" eff. length
      LQ3 = 990.6

C
C     THE POLARIMETER GEOMETRY: TC = TARGET to COLLIMATOR 1 DISTANCE (MM)
C     TM1 = TARGET TO 1st quad MAGNET ENTRANCE
C     SM1 = Solenoid exit to 1st Quad Entrance
C     M1M2= 1st quad MAGNET exit to 2nd quad entrance
C     M2M3= 2nd quad MAGNET exit to 3rd quad entrance
c     M3D = 3rd quad MAG exit to detectors
c     M1C = 1st quad MAG exit to collimator
c     CM2 = collimator 6 to 2nd quad MAG entrance

C     DJG Original distances before 2010

c     col1z = 2000.4            ! collimator 1 distance form target
c     col3z = 2090.9            ! collimator 3
c     col5z = 2203.6            ! collimator 5
c     col6z = 2295.7            ! collimator 6
c     TC   = 2000.4

C     DJG  New distances for Qweak and after
C     DJG  Q2 has been shifted upstream 123.604 mm, so I assume the collimator box moves with it.
C     DJG  An additional 30.48 cm was needed to make more space for the Compton, so I took it up
C     DJG  between Q1 and Q2. Not ideal, but works ok.


c     XDET =  (11.592-0.3048)*1.0E3
c     XDET =  (11.592-0.3048-0.09)*1.0E3 ! rough number from Kelly T.
c     XDET = (11190.45 - 0.06) !~11m is average targ to det. distance,
c     !measured 7/10 JAM

      XDET(1) = 11191.77        !left det   survey Fall 2011
      XDET(2) = 11188.23        !right det  survey Fall 2011

c     col1z = 2000.4-123.604-304.8
c     col3z = 2090.9-123.604-304.8
c     col5z = 2203.6-123.604-304.8
c     col6z = 2295.7-123.604-304.8
c     TC   = 2000.4-123.604-304.8

C     DJG: updated for survey results from July 2010
C     Collimator box is fiducialized relative to center of collimator 5
C     Ideal distance from target is 1757.09 mm, with 4.32 mm offset downstream

      col5z = 1757.09+4.32
      col3z = col5z-112.7
      col1z = col5z-203.2
      col6z = col5z+92.1        !note, this value of 92.1 mm is apparently 4.3 mm larger than survey
      TC=col1z

      col1z = 2000.4-123.604-304.8
      col3z = 2090.9-123.604-304.8
      col5z = 2203.6-123.604-304.8
      col6z = 2295.7-123.604-304.8
      TC   = 2000.4-123.604-304.8


C     DJG GEP/SANE position
c     XQ1 = 1.0008e3

C     2010 position
c     XQ1 = 0.8484e3
      XQ1 = 848.4 - 0.06 + 2.6  !1024.6 is nominal value from dimad. .06 is offset from targ, 2.6 from Q1 JAM

C     DJG GEP/SANE position
c     XQ2 = 3.1195e3
C     2010 position
c     XQ2 = 3.0016196e3-304.8
      XQ2 = 2696.75 - 0.06 -2.84 !z-offset for Q2 is -2.84 JAM
c     XQ3 = 4.334612e3-304.8
      XQ3 = 4029.83 - 0.06 - 0.12 !z-offset for Q3 is -0.12 JAM

      TM1  = XQ1 - LQ1/2.
C     SM1  = TM1 - LSol
      M1M2 = XQ2 - LQ2/2. - XQ1 - LQ1/2.
      M2M3 = XQ3 - LQ3/2. - XQ2 - LQ2/2.
c     M3D  = XDET - XQ3 - LQ3/2.
      M3D(1)  = XDET(1) - XQ3 - LQ3/2.
      M3D(2)  = XDET(2) - XQ3 - LQ3/2.
      M1C  = TC - Tm1 - LQ1
      CM2  = XQ2 - LQ2/2. - col6z




c--   Here we read in the parameters for this run..
c--

c-----------------
c--   If there is an input file specified as a command argument, use it.
c--   If there is no command line argumenet, try to use file mc_input.dat
c--   If that file does not exist, read input from stdin.
      nrarg = iargc()
      if (nrarg .gt. 0) then
         inunit= 31
         do j=1,nrarg
            call getarg(j,argin(j))
         enddo
         write(filename,'("infiles/",(a))') argin(1)
         i=index(filename,' ')
         write(filename(i:),'(".inp")')
         write(6,*) 'filename is', filename
         open(31,file=filename,status='old')
      else
         inunit= 32
         open(unit=32,file='infiles/mc_input.dat',status='old',err=1098)
      endif
      goto 1099
 1098 inunit=5                  !if all else failed, try to read from stdin
 1099 continue
c-----------------

C     RANDOM NUMBER GENERATOR SEED

      if(inunit.eq.5) write(6,'(a)')' Random number generator seed: '
      read(inunit,*) ii
      print*, 'Random number generator seed: ',ii


C     THE BEAM ENERGY
C

      if(inunit.eq.5) write(6,'(a)')' BEAM Energy: '
      read(inunit,*) b_energy

      PBEAM= b_energy           !beam momentem in GeV/c (E=p for e-'s)
      p_up= int(pbeam/2.0)+1.0  !nice upper limit for moller momenta histograms
C
C     THE QUADRUPOLE MAGNET CURRENTS
C
      if(inunit.eq.5) write(6,'(a)')' QUAD currents q1 & q2: '
      read(inunit,*)mcur1, mcur2

c     tip1 = q1tip(mcur1)      !pole tip fields
c     tip2 = q2tip(mcur2)
c     tip3 = q2tip(mcur2)

C     DJG Let MC handle negative currents
C     DJG note that these sign conventions are the opposite
C     DJG what is usually done in TRANSPORT
C     DJG Normally, horizontally focusing (Q1) = positive field
C     DJG           vertically focusing (Q2) = negative field

      if(mcur1.lt.0.0) then
         tip1 = q1tip(abs(mcur1)) !pole tip fields
      else
         tip1 = -q1tip(mcur1)
      endif

      if(mcur2.lt.0.0) then
         tip2 = -q2tip(abs(mcur2))
      else
         tip2 = q2tip(mcur2)
      endif

      tip3 = tip2

c     grad1=-tip1/r1           ! field gradients
      grad1=tip1/r1             ! field gradients (already built into sign of field)
      grad2=tip2/r2
      grad3=tip3/r3

      if(inunit.eq.5) print *,tip1,tip2,mcur1,mcur2,grad1,grad2,grad3
C

      if(inunit.eq.5) write(6,'(a)')' Number of events: '
      read(inunit,*)NUMEVENTS

      if(inunit.eq.5) write(6,'(a)')' scatt. Angle range (from,howmuch) [rad]:'
      read(inunit,*) climit,crange

      if(inunit.eq.5) write(6,'(a)')' BEAM OFFSET (X/Y) and ANGLE (X/Y): '
      read(inunit,*) X0,Y0,TX0,TY0


C     SIGMAS FOR THE BEAM SIZE AND DIVERGENCE AT THE M0LLER TARGET
C
C     DATA SIGX,SIGY,SIGXP,SIGYP/0.500,0.200,0.48D-6,0.55D-6/

      if(inunit.eq.5) then
         write(6,'(a)')' Beam size(X/Y) and divergence(X/Y): '
      endif

      read(inunit,*) SIGX,SIGY,SIGXP,SIGYP

      if(inunit.eq.5) then
         write(6,'(a)')' Radiative Corrections ?(y/n): '
      endif

      read(inunit,'(a)')answer
      if ((answer.eq.'y').or.(answer.eq.'Y')) then
         scflag(2) = .true.
         TCKBRM = .true.
      else
         scflag(2) = .false.
         TCKBRM = .false.
      endif

      if(inunit.eq.5) write(6,'(a)')' Multiple Scattering ?(y/n): '
      read(inunit,'(a)')answer
      if ((answer.eq.'y').or.(answer.eq.'Y')) then
         TAILS = .true.         !Multiple scattering with moliere radius
         scflag(1) = .true.     !multiple satt. before interation
         scflag(3) = .true.     !multiple satt. after interation
      else
         TAILS = .false.
         scflag(1) = .false.
         scflag(3) = .false.
      endif

      if(inunit.eq.5) write(6,'(a)')' Levchuk Effect ?(y/n): '
      read(inunit,'(a)')answer
      if ((answer.eq.'y').or.(answer.eq.'Y')) then
         BINDE = .TRUE.
      else
         BINDE = .false.
      endif



C     THE TARGET THICKNESS (MM)
      if(inunit.eq.5) write(6,'(a)')' Target Thickness [mm]: '
      read(inunit,*) xtgt
*     XTGT = 0.010
ccc   XTGT = 0.020           ! 20 UM LONGITUDINAL


      if(inunit.eq.5) write(6,'(a)')' Q2 position offset [mm]: '
      read(inunit,*) quad2offset

      quad3offset=quad2offset

C     THE COLLIMATOR POSITIONS
      if(inunit.eq.5) then
         write(6,'(a)') 'Settings of Collimators 1,2,3,4,6,7 (cm) //
     >[all .lt. 0.0]:'
      endif
      read(inunit,*) pcol(1),pcol(2),pcol(3),pcol(4),pcol(6),pcol(7)
      pcol(5)= 0.0              !fixed -- centered on beamline and irrelevant.
      do i=1,7
         colxmin(i)= -999.9
         colxmax(i)=  999.9
         colymin(i)= -999.9
         colymax(i)=  999.9
      end do
      colymax(1)= -10.0*pcol(1)
      colymin(2)=  10.0*pcol(2)
      colxmin(3)=  10.0*pcol(3)
      colxmax(4)= -10.0*pcol(4)
      colxmax(6)=  10.0*pcol(6)
      colxmin(7)= -10.0*pcol(7)

      if(inunit.eq.5) then
         write(6,'(a)')' Impose Q1 Exit Aperture Restriction ?(y/n): '
      endif
      read(inunit,'(a)')answer
      if ((answer.eq.'y').or.(answer.eq.'Y')) then
         Q1CLIP = .TRUE.
      else
         Q1CLIP = .false.
      endif

      if(inunit.eq.5) then
         write(6,'(a)')'Use 3rd order fringe eff. in quads ?(y/n): '
      endif
      read(inunit,'(a)')answer
      if ((answer.eq.'y').or.(answer.eq.'Y')) then
         QFRINGE = .TRUE.
      else
         QFRINGE = .false.
      endif

      if(inunit.eq.5) then
         write(6,'(a)')'Use 3 quad tune (high energy) ?(y/n): '
      endif
      read(inunit,'(a)')answer
      if ((answer.eq.'y').or.(answer.eq.'Y')) then
         THREEQUAD= .TRUE.
      else
         THREEQUAD = .false.
      endif

      if(.not.THREEQUAD) then
         tip2=0.0
         grad2=0.0
      endif

      if(inunit.eq.5) write(6,'(a)')'Which solenoid simulation?
     >F=full,S=simple,N=none'
      read(inunit,'(a)')answer
      if ((answer.eq.'f').or.(answer.eq.'F')) then
         SOLSIM = 1
         write(6,'(a)')'Complex field will be used'
         read(inunit,*) BSolLong
         write(6,'(a)')' Enter solenoid
     >dx,dy,dxp,dyp (mm,radians)'
         read(inunit,*) soldx, soldy, soldxp, soldyp
c     LSol=500.   !propagate further through field
         LSol=450.              !propagate further through field
      elseif ((answer.eq.'s').or.(answer.eq.'S')) then
         write(6,*) 'Simple solenoid will be used'
         SOLSIM = -1
         write(6,'(a)')' Enter Solenoid Field (T)'
         read(inunit,*) BSolLong
         write(6,'(a)')' Enter solenoid
     >dx,dy,dxp,dyp (mm,radians)'
         read(inunit,*) soldx, soldy, soldxp, soldyp
c     LSol=125. !DJG: note that Lsol is Leff/2 in this context
         LSol=155.              !DJG: Leff is 310 mm (from field map).
!But we only go through 1/2 of solenoid
      else
         write(6,'(a)') 'No solenoid effect applied'
         solsim=0
         BSolLong=0.0
         soldx=0
         soldy=0
         soldxp=0
         soldyp=0
         LSol=125.
      endif

C     DJG calculate input file-dep. distances
      if(QFRINGE) then
         TM1  = XQ1 - (LQ1-r1)/2. - r1
         M1M2 = XQ2 - (LQ2-r2)/2. - r2 - XQ1 - (LQ1-r1)/2.- r1
         M2M3 = XQ3 - (LQ3-r3)/2. - r3 - XQ2 - (LQ2-r2)/2.- r2
c     M3D  = XDET - XQ3 - (LQ3-r3)/2. - r3
         M3D(1)  = XDET(1) - XQ3 - (LQ3-r3)/2. - r3
         M3D(2)  = XDET(2) - XQ3 - (LQ3-r3)/2. - r3
         M1C  = TC - Tm1 - (LQ1+r1)
         CM2  = XQ2 - (LQ2-r2)/2. - r2 - col6z
      endif

      SM1  = TM1 - LSol


C     if(inunit.eq.5) write(6,'(a)')' Enter Solenoid Field (T)'
C     read(inunit,*) BSolLong
C
C     if(inunit.eq.5) write(6,'(a)')' Enter solenoid dx,dy,dxp,dyp (mm,radians)'
C     read(inunit,*) soldx, soldy, soldxp, soldyp


c---------------------------------------------------------------
      write(6,811) b_energy, PBEAM
 811  format(' Beam Energy: ',f6.3,' Momentum: ',f6.3)
      write(6,812) mcur1,mcur2
 812  format(' Quad Currents: ',2f8.2, ' Amps')
      write(6,813) tip1,tip2,tip3
 813  format(' Quad Tip Fields: ',3f10.4, ' Tesla')
      write(6,814) grad1,grad2,grad3
 814  format(' Quad Gradients: ',3f10.5, ' T/mm')
      GdotL1= (grad1*10.0*10.0)*(LQ1/10.0)
      GdotL2= (grad2*10.0*10.0)*(LQ2/10.0)
      write(6,8141) GdotL1,GdotL2
 8141 format(' Quad GL ',2f10.4, ' (kG/cm)*cm')
      write(6,815) NUMEVENTS
 815  format(' Generating ',i10,' events.')
      write(6,816) climit,crange
 816  format(' Scattering Angle Range: ',2f10.4)
      write(6,817) X0,Y0,TX0,TY0
 817  format(' Beam Offsets: Position x/y and Angle x/y',4f10.4)
      write(6,818) SIGX,SIGY,SIGXP,SIGYP
 818  format(' Beam size(X/Y) and divergence(X/Y):',4f10.4)
      write(6,819) scflag(2)
 819  format(' Radiative Corrections?: ',L1)
      write(6,820) scflag(1)
 820  format(' MCS Before Interaction?: ',L1)
      write(6,821) scflag(3)
 821  format(' MCS After Interaction?: ',L1)
      write(6,822) BINDE
 822  format(' Levchuk Effect??: ',L1)
      write(6,823) xtgt
 823  format(' Target Thickness: ',f10.3)
      write(6,824) pcol
 824  format(' Collimators at:',7f6.1)
      write(6,825) Q1CLIP
 825  format(' Q1 Exit Aperture Clipping?? ',L1)
      write(6,826) QFRINGE
 826  format(' Quad fringe field effects?? ',L1)
      write(6,827) THREEQUAD
 827  format(' 3 QUAD tune?? ',L1)
      write(6,828) BSolLong
 828  format(' Solenoid Field (Tesla): ',f4.2)
      write(6,829) soldx, soldy, soldxp, soldyp
 829  format(' Solenoid offsets (x,y in mm, dxp,dyp in radians):',2f6.2,2f8.5)
c     if(abs(BSolLong).le.MinSolFld) write(6,828)
c     828    format(' ... Solenoid field below cut. No solenoid effect applied.')
c---------------------------------------------------------------



      IHD=0
C
C     IHD DETERMINES THE HEADING TO BE PRINTED
C
      IF(TAILS) THEN
         IHD=IHD+1
      ELSE
         IHD=IHD+2
      ENDIF
C
C     PUT THE DETECTOR POSITION INTO THE HEADING
C
      IF(BINDE) THEN
         WRITE(HEADING,50) TITLE(IHD),Y0,X0,PBEAM
 50      FORMAT('H2 P DIST',A12,1X,'Y',F5.2,' MM',' X',
     >        F5.2,' MM',F7.3,' GEV')
      ELSE
         WRITE(HEADING,51) TITLE(IHD),Y0,X0,PBEAM
 51      FORMAT('FREE E:  ',A12,1X,'Y', F5.2,' MM', ' X',
     >        F5.2,' MM',F7.3,' GEV')
      ENDIF
C
C     INITIALIZE EVERYTHING
ccc
c     call hropen(16,'ntuple','scr2:cebaf.ntuple','nq',4096,istat)      !Ntuple
      if(nrarg.gt.0) then
         write(filename,'("hbook/",(a))') argin(1)
         i=index(filename,' ')
         write(filename(i:),'(".hbook")')
         call hropen(17,'histo',filename,'nq',4096,istat) !Histo
      else
         call hropen(17,'histo',
     &        'hbook/mollermcb.hbook','nq',4096,istat) !Histo
      endif
ccc
      CALL HBOOKM(npads,p_up)
ccc
ccc   initialize ntuple    stored in file cebaf.ntuple  M.L. 12.94
ccc
c     call hcdir('//ntuple',' ')                                        !Ntuple
c     call hbnt(10,'cebaf_ntuple',' ')                          !Ntuple
c     call hbset('bsize',4096*20,ierr)                          !Ntuple
c     call hbname(10,' ', 0, '$clear')                          !Ntuple
c     call hbname(10,'events',ntuple,'collx1,colly1,collx2
c     *  ,colly2,detx1,dety1,detx2,dety2,Axxp,Axxm,Axyp,Axym
c     *  ,Azzp,Azzm,Axxp2,Axxm2,Axyp2,Axym2,Azzp2,Azzm2')               !Ntuple
c     call hcdir('//histo',' ')                                         !Histo
ccc
      DO J=1,3
         DO I=1,2
            BLKTL(I,J)=0.
            BLKTL2(I,J)=0.
         ENDDO
      ENDDO
      avgA=0.0
      avgA2=0.0
      WEITOT=0.
      WEITO2=0.
      nplt= 0
      nplt2= 0
      nplt3= 0
      nplt4= 0
C
ccc
ccc   Define the number of event trials (need many to simulate bound e- tgt)
ccc   Code works up to 4'000'000 events (with ntuples). There was a problem
ccc   at the beginning with the CERN library that it crashed after 1'000'000
ccc   which was solved by increasing iquest(10) from 16000 to 64000 and
ccc   call hropen(16,'ntuple','scr2:cebaf.ntuple','nq',4096,istat)
ccc   instead of ...'nq',1024,istat) the default. M.L. 11.96
ccc
C
c
c________________________________________________________________________
c________________________________________________________________________
c________________________________________________________________________
C
C     LOOP THROUGH EVENTS
c________________________________________________________________________
C
      DO 900 ITRIAL=1,NUMEVENTS
ccc
ccc   Get some feed-back how many event are already done M.L. 11.96
ccc
         if(mod(itrial,100000) .eq. 0) then
            write(6,*) itrial
         endif
ccc
         twomols= .true.
         singlemol(1)= .true.
         singlemol(2)= .true.
C
C     GENERATE BEAM PHASE SPACE FOR THE INCIDENT VECTOR
C
         CALL BPHASE(VEC(1,1),PBEAM)
C
C     ACCOUNT FOR VERTICAL AND HORIZONTAL BEAM OFFSET
C
         call hfdp1(215,vec(1,1),1.0D0) !histo beamspot
         call hfdp1(216,vec(2,1),1.0D0) !histo
C
c________________________________________________________________________
         call hfill(10,1.,0.0,1.0) !count- after moller scattering
         do j=1,2
            if (singlemol(j)) call hfill(10,1.0+j,0.,1.)
         end do
c________________________________________________________________________
c________________________________________________________________________
C
C     CHOOSE A RANDOM INTERACTION DEPTH IN THE TARGET
C
         DEPTH = 0.999*RAN3(II)+0.0005
         XDEPTH=XTGT*DEPTH
C
C     MULTIPLE SCATTER AND BREM THE E- TO THE INTERACTION POINT
C
         if (scflag(1)) then
            CALL BRMSCT(VEC(1,1),XDEPTH,1,TAILS,TCKBRM,WGTBM)
            PB=SQRT(VEC(4,1)**2+VEC(5,1)**2+VEC(6,1)**2)
            if (PB .lt. 0.1*b_energy) goto 900 !low energy cutoff
         else
            call zdrift(vec(1,1),xdepth)
            wgtbm= 1.0
         endif
C
c________________________________________________________________________
         call hfill(10,4.,0.0,1.0) !count- after moller scattering
         do j=1,2
            if (singlemol(j)) call hfill(10,4.0+j,0.,1.)
         end do
c________________________________________________________________________
c________________________________________________________________________
C
C     GENERATE MOLLER SCATTERING
C
         CALL VECGEN(VEC,XTGT,CST,TGEN,PHIGEN,PEL,
     >        WEIGHT,PWT,BINDE,scflag(2), PPOL,
     >        PUNPOL, TGTPOL, targetpol)
         call hfdp1(219,cst,0.65D0) !Histo cosine of CM scat. angle
         mytheta = acos(cst)*57.2958d0
         myphi = phigen*57.2958d0
         call hfdp1(401,mytheta,weight)
         call hfdp1(4030,mytheta,1.0d0)
         call hfdp1(402,myphi,weight)
         wgtcx= weight          !save unpol'd x-section weight
         DO J=1,2               !weights for individual electrons
            weights(j) = (WGTBM*WEIGHT)/(1.D0*NUMEVENTS)
            PB=SQRT(VEC(4,j)**2+VEC(5,j)**2+VEC(6,j)**2)
            if (PB .lt. 0.1*b_energy) then !low energy cutoff
               twomols= .false. ! should be false
               singlemol(j)= .false. ! should be false
            endif
         ENDDO
c     weight for coincidence electron events
         WEIGHT = (WGTBM*WEIGHT)/(1.D0*NUMEVENTS)
C
c________________________________________________________________________
         call hfill(10,7.,0.0,1.0) !count- after moller scattering
         do j=1,2
            if (singlemol(j)) call hfill(10,7.0+j,0.,1.)
         end do
c________________________________________________________________________
C
C     MULTIPLE SCATTER AND BREM THE OUTGOING ELECTRONS
C
         XDEPTH=XTGT*(1.-DEPTH)
         WGTML=1.D0
         DO J=1,2
            if (singlemol(j)) then
               if (scflag(3)) then
                  CALL BRMSCT(VEC(1,J),XDEPTH,1,TAILS,TCKBRM,WGT)
               else
                  call zdrift(vec(1,j),xdepth)
                  wgt= 1.0
               endif
            else

            endif
            IF(WGT.LE.0.) then
               twomols= .false.
               singlemol(j)= .false.
               weights(j)= 0.0
               weight= 0.0
            else
               WGTML=WGTML*WGT
               weights(j)= weights(j)*wgt
               PB=SQRT(VEC(4,j)**2+VEC(5,j)**2+VEC(6,j)**2)
               if (PB .lt. 0.1*b_energy) then !low energy cutoff
                  twomols= .false.
                  singlemol(j)= .false.
               endif
            endif
         ENDDO
         if (twomols) then
            WEIGHT=WEIGHT*WGTML
            call hfdp1(4031,mytheta,1.0d0)
         endif

C
         call hfill(10,10.,0.0,1.0) !count- begin tracking
         do j=1,2
            if (singlemol(j)) call hfill(10,10.0+j,0.,1.)
         end do

c________________________________________________________________________
C
C     PROPAGATE FROM TARGET, AT CENTER OF SOLENOID, TO END OF SOLENOID
C
         if(solsim.lt.0) then
            do j=1,2
               call misalign(vec(1,j), soldx, soldy, soldxp, soldyp)
               call solenoid(vec(1,j),Lsol,BSolLong,10)
               call misalign(vec(1,j),-soldx,-soldy,-soldxp,-soldyp)
            enddo
         elseif(solsim.gt.0) then
            do j=1,2
               call misalign(vec(1,j),soldx,soldy,soldxp,soldyp)
               call solswim(BsolLong,vec(1,j),vec(2,j),vec(3,j),vec(4,j),
     >              vec(5,j),vec(6,j))
c     write(6,*) 'cheesy poofs',j,vec(1,j),vec(2,j),vec(3,j)
               call misalign(vec(1,j),-soldx,-soldy,-soldxp,-soldyp)
            enddo
         endif
c________________________________________________________________________
C
C     PROPAGATE TO FIRST QUAD MAGNET
C
         if(solsim.ne.0) then
            ZDIST = SM1
         else
            ZDIST = TM1
         endif
c     write(6,*) 'cheesy poofs',ZDIST,SM1
         DO J=1,2
            if (singlemol(j))  then
               if(solsim.ne.0) then
                  CALL ZDRIFT(VEC(1,J),ZDIST+(Lsol-VEC(3,J)))
               else
                  CALL ZDRIFT(VEC(1,J),ZDIST)
               endif
               CALL HFDP2(91,VEC(1,J),VEC(2,J),1.D0) !histo Q1 dist.
            endif
         ENDDO
C
c________________________________________________________________________
C
C     PROPAGATE THROUGH THE 1ST QUAD MAGNET
C
c     iii= iii +1
         DO J=1,2
            if (singlemol(j)) then
ccc
ccc   Misalignment of the quad1, done by shifting beam forth and back
ccc   Not the most elegant way but it works M.L. 9.96
ccc
ccc   vec(1,2)=vec(1,2)-2       ! M.L. 9.96 to simulate y shift of quad
ccc   vec(1,1)=vec(1,1)-2       ! M.L. 9.96 to simulate x shift of quad
               if(QFRINGE) then
                  CALL QUADRUPOLE_FRINGE(VEC(1,J),LQ1,r1,tip1)
               else
                  quad1x=0.01   !Mis-alignment of electron = -quad shift JAM
                  quad1y=0.12   !


                  CALL quad_misalign(vec(1,j),quad1x,quad1y) !Misaligns electrons JAM
                       CALL QUADRUPOLE(VEC(1,J),LQ1,grad1,10)
                       CALL quad_misalign(vec(1,j),-quad1x,-quad1y) !re-aligns electrons JAM
                    endif

ccc   vec(1,1)=vec(1,1)+2       ! M.L. 9.96 to simulate x shift of quad
ccc   vec(1,2)=vec(1,2)+2       ! M.L. 9.96 to simulate y shift of quad
                    quad1exit(1,j)= vec(1,j)
                    quad1exit(2,j)= vec(2,j)
                    CALL HFDP2(92,VEC(1,J),VEC(2,J),1.D0) !histo post-Q1 dist.      
                    if( dabs(mytheta - 90.0) .lt. 0.1) then
                       CALL HFDP2(921,VEC(1,J),VEC(2,J),1.D0) !histo Q1 exit for 90 degrees     
                    endif


                 endif
              ENDDO

c________________________________________________________________________
C
c--   Account for any clipping at Q1 exit
              if (Q1CLIP) then
                 DO J=1,2
                    if (singlemol(j)) then
                       rq1_1= sqrt(quad1exit(1,j)**2 + quad1exit(2,j)**2)
                       if( rq1_1 .ge. q1rad ) then
                          singlemol(j)=.false.
                          twomols= .false.
                       endif
                    endif
                 enddo
              endif

              if (twomols) then
                 call hfdp1(4032,mytheta,1.0d0)
              endif

c________________________________________________________________________
C
C     PROPAGATE TO THE COLLIMATOR 1
C
              ZDIST = M1C
              DO J=1,2
                 if (singlemol(j)) call zdrift(vec(1,j),zdist)
                 IF (singlemol(j)) THEN

                    CALL HFDP2(100,VEC(1,J),VEC(2,J),1.D0) !histo pre-coll dist.
                    CALL HFDP2(1001,VEC(1,J),VEC(2,J),1.D0)
                    col(1,1,j)=VEC(1,J)
                    col(2,1,j)=VEC(2,J)
                    NPLT=NPLT+1
c     ntuple(2*j-1)=vec(1,j)            !Ntuple coll x1 x2
c     ntuple(2*j)=vec(2,j)              !Ntuple coll y1 y2
                 ENDIF
              ENDDO
c________________________________________________________________________
C
C     PROPAGATE TO THE COLLIMATOR 3
C
              ZDIST = col3z-col1z
              DO J=1,2
                 if (singlemol(j)) call zdrift(vec(1,j),zdist)
                 IF (singlemol(j)) THEN
                    CALL HFDP2(100,VEC(1,J),VEC(2,J),1.D0) !histo pre-coll dist.
                    CALL HFDP2(1003,VEC(1,J),VEC(2,J),1.D0)
                    col(1,3,j)=VEC(1,J)
                    col(2,3,j)=VEC(2,J)
                    NPLT=NPLT+1
c     ntuple(2*j-1)=vec(1,j)            !Ntuple coll x1 x2
c     ntuple(2*j)=vec(2,j)              !Ntuple coll y1 y2
                 ENDIF
              ENDDO

c________________________________________________________________________
C
C     PROPAGATE TO THE COLLIMATOR 5
C
              ZDIST = col5z-col3z
              DO J=1,2
                 if (singlemol(j)) call zdrift(vec(1,j),zdist)
                 IF (singlemol(j)) THEN
                    CALL HFDP2(100,VEC(1,J),VEC(2,J),1.D0) !histo pre-coll dist.
                    CALL HFDP2(1005,VEC(1,J),VEC(2,J),1.D0)
                    col(1,5,j)=VEC(1,J)
                    col(2,5,j)=VEC(2,J)
                    NPLT=NPLT+1
c     ntuple(2*j-1)=vec(1,j)            !Ntuple coll x1 x2
c     ntuple(2*j)=vec(2,j)              !Ntuple coll y1 y2
                 ENDIF
              ENDDO

c________________________________________________________________________
C
C     PROPAGATE TO THE COLLIMATOR 6
C
              ZDIST = col6z-col5z
              DO J=1,2
                 if (singlemol(j)) call zdrift(vec(1,j),zdist)
                 IF (singlemol(j)) THEN
                    CALL HFDP2(100,VEC(1,J),VEC(2,J),1.D0) !histo pre-coll dist.
                    CALL HFDP2(1006,VEC(1,J),VEC(2,J),1.D0)
                    col(1,6,j)=VEC(1,J)
                    col(2,6,j)=VEC(2,J)
                    NPLT=NPLT+1
c     ntuple(2*j-1)=vec(1,j)            !Ntuple coll x1 x2
c     ntuple(2*j)=vec(2,j)              !Ntuple coll y1 y2
                 ENDIF
              ENDDO
c________________________________________________________________________
C
C     SIMULATE THE COLLIMATOR
C
ccc
ccc   No Collimator is used when reanalyzing with ntuples is done. We can
ccc   then simulate the collimator by cuts in the ntuple analysis M.L. 11.96
ccc
              ienergy = 1
              DO J=1,2
                 if (singlemol(j)) then
                    singlemol(j)= all_coll(col(1,1,j))
                    twomols= ( twomols .and. singlemol(j) )
                 endif
c     if (singlemol(j)) then
c     +         call collimate(vec(1,j),singlemol(j),twomols,ienergy)
                 IF ((NPLT2.LT.5000).and.(singlemol(j))) THEN
                    CALL HFDP2(101,VEC(1,J),VEC(2,J),1.D0) !histo post-coll dist.
                    NPLT2=NPLT2+1
                 ENDIF
              ENDDO
              if (twomols) then
                 call hfdp1(4033,mytheta,1.0d0)
              endif
ccc
c________________________________________________________________________
C
C     PROPAGATE TO THE 2ND QUAD MAGNET
C
              ZDIST = CM2
              DO J=1,2
                 if (singlemol(j)) call zdrift(vec(1,j),zdist)
              ENDDO
c________________________________________________________________________
C
C     PROPAGATE THROUGH THE 2ND QUAD MAGNET
C


              DO J=1,2
                 if (singlemol(j)) then
                    CALL HFDP2(1011,VEC(1,J),VEC(2,J),1.D0) !at Q2 aperture.
                 endif
              enddo
              if(dabs(VEC(1,1)).gt.100..or.dabs(vec(1,2)).gt.100.) go to 900


              DO J=1,2
                 if (singlemol(j)) then
ccc
ccc   Missalignment of the quad2, done by shifting beam forth and back
ccc   Not the most elegant way but it works M.L. 9.96
ccc
ccc   vec(1,2)=vec(1,2)-5       ! M.L. 9.96 to simulate y shift of quad

                    vec(1,j)=vec(1,j) + quad2offset ! simulate x shift of quad
                    quad2entrance(1,j)= vec(1,j)
                    quad2entrance(2,j)= vec(2,j)
c     print *,vec(1,j),j

                    if(tip2.gt.0.0) then
                       if(QFRINGE) then
                          CALL QUADRUPOLE_FRINGE(VEC(1,J),LQ2,r2,tip2)
                       else
                          CALL QUADRUPOLE(VEC(1,J),LQ2,grad2,10)
                       endif
                    else
                       if(QFRINGE) then
                          zdist=LQ2+r2
                       else
                          zdist = LQ2
                       endif
                       call zdrift(vec(1,J),zdist)
                    endif

                    vec(1,j)=vec(1,j) - quad2offset ! simulate x shift of quad
                    quad2exit(1,j)= vec(1,j)
                    quad2exit(2,j)= vec(2,j)
ccc   vec(1,2)=vec(1,2)+5       ! M.L. 9.96 to simulate y shift of quad
c     write(6,*) vec(1,j),vec(2,j),vec(3,j),vec(4,j),vec(5,j),vec(6,j),j,'l'
                 endif
              ENDDO



c________________________________________________________________________
C
C     PROPAGATE TO THE 3rd QUAD MAGNET
C
              ZDIST = M2M3
              DO J=1,2
                 if (singlemol(j)) call zdrift(vec(1,j),zdist)
              ENDDO
c________________________________________________________________________
C
C     PROPAGATE THROUGH THE 3RD QUAD MAGNET
C

c     dg        DO J=1,2
c     dg           if (singlemol(j)) then
c     dg              CALL HFDP2(1011,VEC(1,J),VEC(2,J),1.D0) !at Q2 aperture.
c     dg           endif
c     dg        enddo
c     dg        if(dabs(VEC(1,1)).gt.100..or.dabs(vec(1,2)).gt.100.) go to 900

              DO J=1,2
                 if (singlemol(j)) then
ccc
ccc   Missalignment of the quad3, done by shifting beam forth and back
ccc   Not the most elegant way but it works M.L. 9.96

                    quad3entrance(1,j)= vec(1,j)
                    quad3entrance(2,j)= vec(2,j)

                    quad3x=-0.06 !Mis-alignment of electron = -quad shift JAM
                    quad3y=-0.15 !

                    if(QFRINGE) then
                       CALL QUADRUPOLE_FRINGE(VEC(1,J),LQ3,r3,tip3)
                    else
                       CALL quad_misalign(vec(1,j),quad3x,quad3y)

                            CALL QUADRUPOLE(VEC(1,J),LQ3,grad3,10)

                            CALL quad_misalign(vec(1,j),-quad3x,-quad3y)
                         endif
                         quad3exit(1,j)= vec(1,j)
                         quad3exit(2,j)= vec(2,j)

                      endif
                   ENDDO

C     DJG This only gets us to the magnetic exit - now we need to get to the mechanical exit.
C     DJG NOTE: if you use the fringe stuff - this will not work!


c________________________________________________________________________
C
C
C     PROPAGATE TO MECHANICAL END of Q3
C
                   ZDIST = 88.9
                   DO J=1,2
                      if (singlemol(j)) call zdrift(vec(1,j),zdist)
ccc   IF ((NPLT3.LT.5000).and.(singlemol(j))) THEN
                      q3mechexit(1,j) = vec(1,j)
                      q3mechexit(2,j) = vec(2,j)
                   ENDDO


c________________________________________________________________________
C
C
C     PROPAGATE TO THE DETECTOR PLANE
C
c     ZDIST = M3D-88.9
                   DO J=1,2
                      if(vec(1,j).lt.0.0) then !left
                         ZDIST=M3D(1)-88.9
                      else      !right
                         ZDIST=M3D(2)-88.9
                      endif
                      if (singlemol(j)) call zdrift(vec(1,j),zdist)
ccc   IF ((NPLT3.LT.5000).and.(singlemol(j))) THEN
                      IF (singlemol(j)) THEN
                         CALL HFDP2(102,VEC(1,J),VEC(2,J),1.D0) !histo detector dist.
                         NPLT3=NPLT3+1
c     ntuple(2*j+3)=vec(1,j)             !Ntuple det x1 x2 M.L. 12.94
c     ntuple(2*j+4)=vec(2,j)             !Ntuple det y1 y2
                      ENDIF

                   ENDDO

C     DJG All this stuff is just for making plots are various apertures!!!

C     THIS IS UGLY!!!

C     THEN DRIFT BACKWARDS to ARM EXIT
C
                   ZDIST = -499.0
                   DO J=1,2
                      if (singlemol(j)) call zdrift(vec(1,j),zdist)
ccc   IF ((NPLT3.LT.5000).and.(singlemol(j))) THEN
                      armexit(1,j) = vec(1,j)
                      armexit(2,j) = vec(2,j)
                   ENDDO

C     THEN DRIFT BACKWARDS to ARM ENTRANCE
C
                   ZDIST = -3396.0
                   DO J=1,2
                      if (singlemol(j)) call zdrift(vec(1,j),zdist)
ccc   IF ((NPLT3.LT.5000).and.(singlemol(j))) THEN
                      armentr(1,j) = vec(1,j)
                      armentr(2,j) = vec(2,j)
                   ENDDO

C     THEN DRIFT TO ENRANCE OF "GOOFY" BOX
C
                   ZDIST = -2597.0
                   DO J=1,2
                      if (singlemol(j)) call zdrift(vec(1,j),zdist)
ccc   IF ((NPLT3.LT.5000).and.(singlemol(j))) THEN
                      pipeexit(1,j) = vec(1,j)
                      pipeexit(2,j) = vec(2,j)
                   ENDDO

C     THEN DRIFT BACK TO THE DETECTORS
C
                   ZDIST = 499.0+3396.0+2597.0
                   DO J=1,2
                      if (singlemol(j)) call zdrift(vec(1,j),zdist)
                   ENDDO


c________________________________________________________________________
C
C     SUM THE RATE OF ELECTRON PAIRS HITTING THE DETECTOR PLANE
C
                   if (twomols) then
                      WEITOT=WEITOT+WEIGHT
                      WEITO2=WEITO2+WEIGHT**2
                   endif
                   do j=1,2
                      if (singlemol(j)) then
                         WEITOTs(j)=WEITOTs(j)+WEIGHTs(j)
                         WEITO2s(j)=WEITO2s(j)+WEIGHTs(j)**2
                      endif
                   enddo
C
                   if (twomols) then
                      call hfdp1(4034,mytheta,1.0d0)
                   endif
c________________________________________________________________________
C
C     SIMULATE THE DETECTOR
C
                   call detectors(vec,padtl,padtl2,blktl,blktl2,
     +                  weight,weights,pwt,wgtcx,singlemol,twomols,avgA,avgA2)
C
                   if (twomols) then
                      call hfdp1(4035,mytheta,1.0d0)
                   endif

                   do j=1,2
                      jj= 3-j
                      IF (singlemol(j)) THEN
                         CALL HFDP2(103,VEC(1, J),VEC(2, J),1.D0) !histo detector dist.
                         CALL HFDP2(1030+J,VEC(1,J),VEC(2,J),1.D0) !histo mating  dist.
                         CALL HFDP2(1030+J,VEC(1,JJ),VEC(2,JJ),1.D0) !histo mating  dist.
                         PMOLL(j)=SQRT(VEC(4,j)**2+VEC(5,j)**2+VEC(6,j)**2)
                         CALL HFDP1(200+j,PMOLL(j),1.D0) !histo detector dist.
                         if (twomols) then
                            call hfdp1(2005+10*j,PMOLL(j),1.D0) !mom of conicidence electrons
                         end if
                      ENDIF
                   end do
ccc   IF ((NPLT4.LT.5000).and.(twomols)) THEN
                   IF (twomols) THEN
c     IF (1.eq.1) THEN
                      twom = twom + 1.
                      NPLT4=NPLT4+1
                      CALL HFDP2(104,VEC(1,1),VEC(2,1),1.D0) !histo detector coinc. dist.
                      CALL HFDP2(104,VEC(1,2),VEC(2,2),1.D0)
                      call hfdp2(105,vec(1,1),vec(2,1),1.d0)
                      call hfdp2(106,vec(1,2),vec(2,2),1.d0)
                      call hfdp2(2021,vec(1,1),pmoll(1),1.0d0)
                      call hfdp2(2022,vec(1,2),pmoll(2),1.0d0)
                      call hfdp1(403,mytheta,weight)
                      call hfdp1(404,myphi,weight)
                      call hfdp1(423,mytheta,1.0d0)
                      call hfdp1(424,myphi,1.0d0)
                      CALL HFDP2(405,col(1,1,1),col(2,1,1),1.D0)
                      CALL HFDP2(405,col(1,1,2),col(2,1,2),1.D0)
                      CALL HFDP2(406,col(1,3,1),col(2,3,1),1.D0)
                      CALL HFDP2(406,col(1,3,2),col(2,3,2),1.D0)
                      CALL HFDP2(407,col(1,5,1),col(2,5,1),1.D0)
                      CALL HFDP2(407,col(1,5,2),col(2,5,2),1.D0)
                      CALL HFDP2(408,col(1,6,1),col(2,6,1),1.D0)
                      CALL HFDP2(408,col(1,6,2),col(2,6,2),1.D0)
c--   Distribution of Moller coincidences' positions at quad...
                      CALL HFDP2(925,quad1exit(1,1),quad1exit(2,1),1.D0)
                      CALL HFDP2(925,quad1exit(1,2),quad1exit(2,2),1.D0)
                      CALL HFDP2(926,quad2entrance(1,1),quad2entrance(2,1),1.D0)
                      CALL HFDP2(926,quad2entrance(1,2),quad2entrance(2,2),1.D0)
                      CALL HFDP2(927,quad2exit(1,1),quad2exit(2,1),1.D0)
                      CALL HFDP2(927,quad2exit(1,2),quad2exit(2,2),1.D0)
                      CALL HFDP2(928,quad3entrance(1,1),quad3entrance(2,1),1.D0)
                      CALL HFDP2(928,quad3entrance(1,2),quad3entrance(2,2),1.D0)
                      CALL HFDP2(929,quad3exit(1,1),quad3exit(2,1),1.D0)
                      CALL HFDP2(929,quad3exit(1,2),quad3exit(2,2),1.D0)

                      CALL HFDP2(930,q3mechexit(1,1),pipeexit(2,1),1.D0)
                      CALL HFDP2(930,q3mechexit(1,2),pipeexit(2,2),1.D0)

                      CALL HFDP2(931,pipeexit(1,1),pipeexit(2,1),1.D0)
                      CALL HFDP2(931,pipeexit(1,2),pipeexit(2,2),1.D0)

                      CALL HFDP2(932,armentr(1,1),armentr(2,1),1.D0)
                      CALL HFDP2(932,armentr(1,2),armentr(2,2),1.D0)

                      CALL HFDP2(933,armexit(1,1),armexit(2,1),1.D0)
                      CALL HFDP2(933,armexit(1,2),armexit(2,2),1.D0)

c--   Account for any clipping at Q1 exit
                      rq1_1= sqrt(quad1exit(1,1)**2 + quad1exit(2,1)**2)
                      rq1_2= sqrt(quad1exit(1,2)**2 + quad1exit(2,2)**2)
c     if( (rq1_1 .lt. q1rad) .and. (rq1_2 .lt. q1rad) ) then
                      if( (rq1_1 .ge. q1rad) .or. (rq1_2 .ge. q1rad) ) then
c     CALL HFDP2(114,VEC(1,1),VEC(2,1),1.D0) !histo detector coinc. dist.
c     CALL HFDP2(114,VEC(1,2),VEC(2,2),1.D0)
c     call hfdp1(413,mytheta,weight)
c     call hfdp1(414,myphi,weight)
c     CALL HFDP2(415,col(1,1,1),col(2,1,1),1.D0)
c     CALL HFDP2(415,col(1,1,2),col(2,1,2),1.D0)
c     CALL HFDP2(416,col(1,3,1),col(2,3,1),1.D0)
c     CALL HFDP2(416,col(1,3,2),col(2,3,2),1.D0)
c     CALL HFDP2(417,col(1,5,1),col(2,5,1),1.D0)
c     CALL HFDP2(417,col(1,5,2),col(2,5,2),1.D0)
c     CALL HFDP2(418,col(1,6,1),col(2,6,1),1.D0)
c     CALL HFDP2(418,col(1,6,2),col(2,6,2),1.D0)

                   endif
                   if (TGTPOL.eq.1) then
                      call hfdp1(313,PPOL,1.0D0)
                   else
                      call hfdp1(314,PUNPOL,1.0D0)
                   endif
                ENDIF


 900         CONTINUE           !>>>>>>>>>>>>>>>>end of main loop<<<<<<<<<<<<<<<<<<<<
c=========================================================================
c=========================================================================

             efficiency= twom/float(NUMEVENTS)
C
C     CALCULATE THE LINESHAPE AND ASYMMETRIES
C
C
             SIGNAL=BLKTL(1,3)+BLKTL(2,3)
             avgA=avgA/signal
             avgA2=avgA2/signal
             SIGNL2=SQRT(BLKTL2(1,3)+BLKTL2(2,3))
             DIFF = BLKTL(1,3)-BLKTL(2,3)
             IF(SIGNAL.GT.0.) THEN
                ASZ=(BLKTL(2,3)-BLKTL(1,3))/SIGNAL
                ASZ2=SQRT((1.-ASZ)**2*BLKTL2(2,3)+
     +               (1.+ASZ)**2*BLKTL2(1,3))/SIGNAL
                ASX=(BLKTL(2,1)-BLKTL(1,1))/SIGNAL
                ASX2=SQRT((1.-ASX)**2*BLKTL2(2,1)+
     +               (1.+ASX)**2*BLKTL2(1,1))/SIGNAL
                ASY=(BLKTL(2,2)-BLKTL(1,2))/SIGNAL
                ASY2=SQRT((1.-ASY)**2*BLKTL2(2,2)+
     +               (1.+ASY)**2*BLKTL2(1,2))/SIGNAL
                IF(BINDE) THEN
                   ASZ=ASZ/targetpol
                   ASZ2=ASZ2/targetpol
                   ASX=ASX/targetpol
                   ASX2=ASX2/targetpol
                   ASY=ASY/targetpol
                   ASY2=ASY2/targetpol
                ENDIF
             ELSE
                ASZ=0.
                ASX=0.
                ASY=0.
                ASZ2=0.
                ASX2=0.
                ASY2=0.
             ENDIF
c
c--   ...for doubles...
c
             DBsignal= signal/2. !compensate for fact that for every beam electron
             DBsignl2= signl2/1.4142 !2 events were generated (+/- pol'n)
             Dbsig= DBsignal*6.25E+12 !convert to events/second (for 1uA beam)
             Dbsige= DBsignl2*6.25E+12
             DbAzz= asz
             DbAzze= asz2
c
c--   Calculate the polarization factor for use in analyzing the data
c
             if (asz.gt.0.0) then
                polfact= 1.0/(targetpol*asz)
                polfacte=asz2/(targetpol*asz*asz)
             else
                polfact= 0.0
                polfacte= 0.0
             endif

             davgA = sqrt(abs(avgA2-avgA**2)/twom)

             dgpolfact=1.0/avgA
             dgpolfacte=(davgA/avgA)*dgpolfact


c=========================================================================

c--   print out the run summary

             write(6,811) b_energy, PBEAM
             write(6,812) mcur1,mcur2
             write(6,813) tip1,tip2
             write(6,814) grad1,grad2
             write(6,815) NUMEVENTS
             write(6,816) climit,crange
             write(6,817) X0,Y0,TX0,TY0
             write(6,818) SIGX,SIGY,SIGXP,SIGYP
             write(6,819) scflag(2)
             write(6,820) scflag(1)
             write(6,821) scflag(3)
             write(6,822) BINDE
             write(6,823) xtgt
             write(6,*)'Efficiency: ',efficiency
             write(6,2004) polfact,polfacte
 2004        format(' Simulated polarization factor is ',f9.3,'+-',f5.3)
             write(6,2006) dgpolfact,dgpolfacte
 2006        format(' CHEESY POOFS: Gaskell polfac is  ',f9.3,'+-',f5.3)
             write(6,2005) Dbsig, Dbsige
 2005        format(' Anticipated Moller rate is ',f7.0,'+-',f4.0,
     >            ' Events/sec/microamp.')


c________________________________________________________________________
             open(23,file='effic.dat',access='append',status='unknown')
             write(23,'(2e12.4)')mcur1,mcur2,efficiency
             close(23)
c________________________________________________________________________
C
C     PRINT THE EVENT TOTALS
C
             PRINT 1000, WEITOT,SQRT(WEITO2)
 1000        FORMAT(1X,'TOTAL # OF E- ON DET = ',E11.4,'+-',E11.4)
             PRINT 1001, WEITOTs(1),SQRT(WEITO2s(1))
 1001        FORMAT(1X,'TOTAL # OF E- ON lower DET = ',E11.4,'+-',E11.4)
             PRINT 1002, WEITOTs(2),SQRT(WEITO2s(2))
 1002        FORMAT(1X,'TOTAL # OF E- ON upper DET = ',E11.4,'+-',E11.4)
c
c--------------------------------------------------
             OPEN(14,FILE='AVASYM.dat',ACCESS='APPEND',STATUS='UNKNOWN')
C
             WRITE(14,1400) HEADING
 1400        FORMAT(/,1H0,5X,A55)
             write(14,'(8E10.2)') SIGX,SIGY,SIGXP,SIGYP,X0,Y0,TX0,TY0
             write(14,'(a,2f7.1,a,f7.1)')' Quad curr. 1 & 2: ',mcur1, mcur2,
     >            ' Q2 offset = ',quad2offset
             WRITE(14,1401)
 1401        FORMAT(' Chan',8X,'Signal',11x,'ZZ-Asymmetry',4x,'XX-Asymmetry',
     >4x,'XY-Asymmetry')
             WRITE(14,1406) K,L,SIGNAL,SIGNL2,ASZ,ASZ2,ASX,ASX2,ASY,ASY2
 1406        FORMAT(2X,I1,'-',I1,1X,E10.3,'+-',E9.3,1X,F7.4,'+-',F6.4,
     >            1X,F7.4,'+-',F6.4,1X,F7.4,'+-',F6.4)
             write(14,1407) efficiency
 1407        format(' Efficiency: ',e12.4)
             write(14,1408) pcol
 1408        format(' Collimators: ',7f7.2)

             CLOSE(14)
c--------------------------------------------------

c------
c--   WRITE PERFORMANCE SUMMARY IN TERSE FORMAT TO 'summary.dat' FOR POSSIBLE
C--   EXTERNAL USE
C------
             OPEN(9,FILE='summary.dat',ACCESS='APPEND',STATUS='UNKNOWN')
             write(9,2003) mcur1, mcur2, pcol, Dbsig, Dbsige,
     $            DbAzz, DbAzze, polfact, polfacte, efficiency,
     $            X0,Y0,TX0,TY0
 2003        format(1x,2f7.1, 7f6.2, 2(1x,F9.1),
     $            2(1x,F8.5),2(1x,f5.2),e12.4,
     $            2f5.1, 2e8.2)
             close(9)
C
C-----------------------------------------------------------------------
c     now for single arm results...
             OPEN(15,FILE='AAVASYM_sing.dat',ACCESS='APPEND',STATUS='UNKNOWN')
             WRITE(15,1500) HEADING
 1500        FORMAT(/,1H0,5X,A55)
             WRITE(15,1502)
 1502        FORMAT(' Chan',6X,'Signal',11x,'ZZ-Asymmetry',4x,'XX-Asymmetry',
     >4x,'XY-Asymmetry')
             DO K=1,npads
                DO L=1,2
                   Azz(K,L)= 0.
                   Azze(K,L)= 0.
                   siSIGNAL=padtl(1,3,K,L)+padtl(2,3,K,L)
                   siSIGNL2=SQRT(padtl2(1,3,K,L)+padtl2(2,3,K,L))
                   DIFF = padtl(1,3,K,L)-padtl(2,3,K,L)
                   IF(siSIGNAL.GT.0.) THEN
                      ASZ=(padtl(2,3,K,L)-padtl(1,3,K,L))/siSIGNAL
                      ASZ2=SQRT((1.-ASZ)**2*padtl2(2,3,K,L)
     >                     +(1.+ASZ)**2*padtl2(1,3,K,L))/siSIGNAL
                      ASX=(padtl(2,1,K,L)-padtl(1,1,K,L))/siSIGNAL
                      ASX2=SQRT((1.-ASX)**2*padtl2(2,1,K,L)
     >                     +(1.+ASX)**2*padtl2(1,1,K,L))/siSIGNAL
                      ASY=(padtl(2,2,K,L)-padtl(1,2,K,L))/siSIGNAL
                      ASY2=SQRT((1.-ASY)**2*padtl2(2,2,K,L)
     >                     +(1.+ASY)**2*padtl2(1,2,K,L))/siSIGNAL
                      IF(BINDE) THEN
                         ASZ=ASZ/targetpol
                         ASZ2=ASZ2/targetpol
                         ASX=ASX/targetpol
                         ASX2=ASX2/targetpol
                         ASY=ASY/targetpol
                         ASY2=ASY2/targetpol
                      ENDIF
                   ELSE
                      ASZ=0.
                      ASX=0.
                      ASY=0.
                      ASZ2=0.
                      ASX2=0.
                      ASY2=0.
                   ENDIF
                   Azz(K,L)= ASZ
                   Azze(K,L)= ASZ2
                   WRITE(15,2001) K,L,siSIGNAL,siSIGNL2,ASZ,ASZ2,ASX,ASX2,ASY,ASY2
                ENDDO
c     write(7,2002) K, Azz(K,1), Azze(K,1), Azz(K,2), Azze(K,2)
             ENDDO
             CLOSE(15)
C
C-----------------------------------------------------------------------

             call HPAK(205,Azz(1,1)) !Histo
             call HPAKE(205,Azze(1,1)) !Histo
             call HPAK(206,Azz(1,2)) !Histo
             call HPAKE(206,Azze(1,2)) !Histo
 2001        FORMAT(1X,I2,'-',I1,1X,E10.3,'+-',E9.3,1X,F7.4,'+-',F6.4,
     >            1X,F7.4,'+-',F6.4,1X,F7.4,'+-',F6.4)
 2002        format(1x,I2,1x,4(F8.6,1x))
C

C     PRINT HISTOGRAMS
C
             call hrout(0,icycle,' ') !Histo
             call hrend('histo') !Histo

             STOP
             END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE ZDRIFT(VEC,ZDIST)
      IMPLICIT NONE
      REAL*8 VEC(6),ZDIST
      if(VEC(6).eq.0.) print *,'VEC(6)=0 !'
      VEC(1)=VEC(1)+ZDIST*VEC(4)/VEC(6)
      VEC(2)=VEC(2)+ZDIST*VEC(5)/VEC(6)
      VEC(3)=VEC(3)+ZDIST
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE BEND(VEC,BRAD,BDIR,BLEN)
      IMPLICIT NONE
      REAL*8 VEC(6),BRAD,BDIR,BLEN,THETA,PHI,ODIR,ODIR0,Z0,PZO
C
C     CALCULATE VARIOUS RELEVANT GEOMETRIC QUANTITIES
      ODIR=3-BDIR
      THETA=ATAN(VEC(3+ODIR)/VEC(6))
      Z0=VEC(3)+BRAD*SIN(THETA)
      ODIR0=VEC(ODIR)-BRAD*COS(THETA)
      PZO=SQRT(VEC(3+ODIR)**2+VEC(6)**2)
      VEC(3)=VEC(3)+BLEN
      VEC(ODIR)=ODIR0-ODIR0*SQRT(BRAD**2-(VEC(3)-Z0)**2)/ABS(ODIR0)
      PHI=ATAN((VEC(3)-Z0)/(ODIR0-VEC(ODIR)))
      VEC(BDIR)=VEC(BDIR)+BRAD*(THETA-PHI)*(VEC(3+BDIR)/PZO)
      VEC(3+ODIR)=PZO*SIN(PHI)
      VEC(6)=PZO*COS(PHI)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE HBOOKM(npads,p_up)

      implicit none

      INTEGER NWPAWC, H, npads
      real*4 p_up
      real rpads,dmx,dmxc
ccc   PARAMETER(NWPAWC= 320000)
      PARAMETER(NWPAWC= 5000000)
      COMMON /PAWC/ H(NWPAWC)

c--   Count some things
      call hbook1(10,'Scalers', 20,0.5,20.5, 0.0)
      rpads= npads
      dmxc= 60.0                !histogram limit for collimator planes
      dmx= 1000.0               !histogram limit for detector plane
      CALL HBOOK2(91,
     +     ' Q1 Entry',200,-dmxc,dmxc,100,-50.,50.,80000.)
      CALL HBOOK2(92,
     +     ' Q1 Exit',200,-dmxc,dmxc,100,-50.,50.,80000.)
      CALL HBOOK2(921,
     +     ' Q1 Exit for 90 degrees',200,-dmxc,dmxc,100,-50.,50.,80000.)
      CALL HBOOK2(925,
     +     ' Q1 Exit - coincidences',200,-dmxc,dmxc,100,-50.,50.,80000.)
      CALL HBOOK2(926,
     +     ' Q2 Entrance - coincidences',200,-160.,160.,60,-160.,160.,80000.)
      CALL HBOOK2(927,
     +     ' Q2 Exit - coincidences',200,-160.,160.,60,-160.,160.,80000.)
      CALL HBOOK2(928,
     +     ' Q3 Entrance - coincidences',200,-160.,160.,60,-160.,160.,80000.)
      CALL HBOOK2(929,
     +     ' Q3 Exit - coincidences',200,-160.,160.,60,-160.,160.,80000.)

      CALL HBOOK2(930,
     +     ' Q3 MECAHNICAL Exit - coincidences',200,-160.,160.,60,-160.,160.,80000.)
      CALL HBOOK2(931,
     +     ' Q3 PIPE Exit - coincidences',200,-160.,160.,60,-160.,160.,80000.)
      CALL HBOOK2(932,
     +     ' ARM ENTRANCE - coincidences',200,-500.,500.,60,-30.,30.,80000.)
      CALL HBOOK2(933,' ARM EXIT',500,-dmx,dmx,100,-50.,50.,80000.)
      CALL HBOOK2(100, ' COLLIMATOR PLANE - all hits',
     $     200,-dmxc,dmxc,100,-50.,50.,80000.)
      CALL HBOOK2(1001,
     +     ' At Collimator 1',200,-dmxc,dmxc,100,-50.0,50.0,80000.)
      CALL HBOOK2(1003,
     +     ' At Collimator 3',200,-dmxc,dmxc,100,-50.0,50.0,80000.)
      CALL HBOOK2(1005,
     +     ' At Collimator 5',200,-dmxc,dmxc,100,-50.0,50.0,80000.)
      CALL HBOOK2(1006,
     +     ' At Collimator 6',200,-dmxc,dmxc,100,-50.0,50.0,80000.)
      CALL HBOOK2(101,' COLLIMATOR PLANE - after collimation',
     $     200,-dmxc,dmxc,50,-50.,50.,80000.)
      CALL HBOOK2(1011,' Q2 PLANE',200,-120.,120.,60,-120.,120.,80000.)
      CALL HBOOK2(102,' DETECTOR PLANE',200,-dmx,dmx,200,-dmx,dmx,80000.)
      CALL HBOOK2(1025,' DETECTOR PLANE- ZOOM',500,-dmx,dmx,100,-50.,50.,80000.)
      CALL HBOOK2(103,' DETECTOR- at least SINGLE ',500,-dmx,dmx,100,-50.,50.,80000.)
      CALL HBOOK2(1031,' DETECTOR MATES',500,-dmx,dmx,100,-50.,50.,80000.)
      CALL HBOOK2(1032,' DETECTOR MATES',500,-dmx,dmx,100,-50.,50.,80000.)
      CALL HBOOK2(104,' DETECTOR COINC',500,-dmx,dmx,100,-50.,50.,80000.)
      CALL HBOOK2(114,' DETECTOR coinc passing Q1',500,-dmx,dmx,100,-50.,50.,80000.)
      CALL HBOOK2(105,' DETECTOR coinc LEFT',
     +     500,-dmx,-400.,100,-50.,50.,80000.)
      CALL HBOOK2(106,' DETECTOR coinc RIGHT',
     +     500,400.,dmx,100,-50.,50.,80000.)
      call hbook1(201,' Final Momentum-1',200,0.,5.,0.)
      call hbook1(202,' Final Momentum-2',200,0.,5.,0.)
      call hbook1(2015,' Accepted  Momentum-1',200,0.,5.,0.)
      call hbook1(2025,' Accepted  Momentum-2',200,0.,5.,0.)
      call hbook2(2021,' Momentum vs. x at Detector',100,-600.,-400.,100,0.0,p_up,0.)
      call hbook2(2022,' Momentum vs. x at Detector',100, 400., 600.,100,0.0,p_up,0.)
      call hbook1(203,' Lower Pad Distribution 1',npads,-4.,7.,0.)
      call hbook1(204,' Upper Pad Distribution 1',npads,-4.,7.,0.)
      call hbook1(205,' Upper Pad-1 Analyzing Powers ',npads,0.,rpads,0.)
      call hbook1(206,' Upper Pad-2 Analyzing Powers ',npads,0.,rpads,0.)
      call hbook1(215,' Beam X Distribution ',21,-10.,10.,0.)
      call hbook1(216,' Beam Y Distribution ',21,-10.,10.,0.)
      call hbook1(219,' Cosine CM Theta ', 201,-1.0,1.0,0.)
      call hbook1(310,' Polar. Momentum dist. (KeV/c)',200,-200.,200.,0.)
      call hbook1(311,' un-Polar. Momentum dist. (KeV/c)',200,-200.,200.,0.)
      call hbook1(312,' S0, Invariant Mass ',1000,.0,.05,0.)
      call hbook1(313,' COIN Polar. Momentum dist. (KeV/c)',200,-200.,200.,0.)
      call hbook1(314,' COIN un-Polar. Momentum dist. (KeV/c)',200,-200.,200.,0.)
      call hbook1(401,' all THETA', 50,75.,105.,0.)
      call hbook1(402,' all PHI', 80, 160.,200.,0.)
      call hbook1(403,' THETA acceptance', 50,75.,105.,0.)
      call hbook1(4030,' THETA acceptance test 0', 50,75.,105.,0.)
      call hbook1(4031,' THETA acceptance test 1', 50,75.,105.,0.)
      call hbook1(4032,' THETA acceptance test 2', 50,75.,105.,0.)
      call hbook1(4033,' THETA acceptance test 3', 50,75.,105.,0.)
      call hbook1(4034,' THETA acceptance test 4', 50,75.,105.,0.)
      call hbook1(4035,' THETA acceptance test 5', 50,75.,105.,0.)
      call hbook1(404,' PHI acceptance',80,160.,200.,0. )
      call hbook1(423,' THETA acceptance unweighted', 50,75.,105.,0.)
      call hbook1(424,' PHI acceptance unweighted',80,160.,200.,0. )
      CALL HBOOK2(405,' COLL 1 - accepetd electrons',200,-dmxc,dmxc,100,-50.,50.,0.)
      CALL HBOOK2(406,' COLL 3 - accepetd electrons',200,-dmxc,dmxc,100,-50.,50.,0.)
      CALL HBOOK2(407,' COLL 5 - accepetd electrons',200,-dmxc,dmxc,100,-50.,50.,0.)
      CALL HBOOK2(408,' COLL 6 - accepetd electrons',200,-dmxc,dmxc,100,-50.,50.,0.)
c     call hbook1(413,' THETA acceptance - and pass Q1', 50,75.,105.,0.)
c     call hbook1(414,' PHI acceptance - and pass Q1',80,160.,200.,0. )
c     CALL HBOOK2(415,' COLL 1 - and pass Q1',200,-dmxc,dmxc,100,-50.,50.,0.)
c     CALL HBOOK2(416,' COLL 3 - and pass Q1',200,-dmxc,dmxc,100,-50.,50.,0.)
c     CALL HBOOK2(417,' COLL 5 - and pass Q1',200,-dmxc,dmxc,100,-50.,50.,0.)
c     CALL HBOOK2(418,' COLL 6 - and pass Q1',200,-dmxc,dmxc,100,-50.,50.,0.)

      CALL HBOOK2(504,' Hodo hit pattern RvL',39,-12.,27.,39,-12.,27.,0.)
      CALL HBOOK2(505,' Hodo hits - coinc RvL',16,0.,16.,16,0.,16.,0.)
      call HBOOK1(511,' HODO diag proj (R+L)/2 for R-L"o#2',16,0.,16.,0.) 
      call HBOOK1(512,' HODO diag proj (R+L)/2 for R-L"o#3',16,0.,16.,0.) 
      call HBOOK1(513,' HODO diag proj (R-L) for R+L=15"A#2',30,-15.0,15.,0.) 
      call HBOOK1(514,' HODO diag proj (R-L) for R+L=15"A#3',30,-15.0,15.,0.) 
      call HBOOK1(515,' HODO diag proj (R-L) for R+L=15"A#7',30,-15.0,15.,0.) 
      call HBOOK1(516,' HODO diag proj (R-L) for R+L=15"A#15',30,-15.0,15.,0.) 


      call HBOOK1(600,'AZ',200,-0.3,0.3,0.0)


      print *, 'HBOOK is initialized.'
      RETURN
      END


C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE HFDP1(ID,XDP,WDP)
C
      IMPLICIT NONE
      INTEGER ID
      REAL*8  XDP, WDP
      REAL*4  X,   W
C
      X = XDP
      W = WDP
      CALL HF1(ID,X,W)          !*** MODIFIED

C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE HFDP2(ID,XDP,YDP,WDP)
C
      IMPLICIT NONE
      INTEGER ID
      REAL*8  XDP, YDP, WDP
      REAL*4  X,   Y,   W
C
      X = XDP
      Y = YDP
      W = WDP
      CALL HF2(ID,X,Y,W)        !*** MODIFIED
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE VECGEN(VEC,XTGT,COSTR,TGEN,PHIGEN,PGEN,
     >     WEIGHT,PWT,BINDE,scflag, PPOL,
     >     PUNPOL, TGTPOL, targetpol)
C
C     *************************************************************
C     GENERATE 6-VECTORS FOR THE MONTE CARLO.
C     VEC(I,1) CONTAINS THE INCIDENT BEAM COORDINATES, THE 2
C     GENERATED 6-VECTORS ARE OUTPUT AS VEC(I,1)+VEC(I,2), COSTR
C     IS THE COS(CM SCATTERING ANGLE), WEIGHT IS THE EVENT WEIGHT,
C     PWT IS EVENT WEIGHT FOR 2 POL STATES TIMES 3 BEAM-TGT POL
C     COMBINATIONS (PZ-PZ, PX-PX, PX-PZ),
C     BINDE SWITCHES ON ELECTRON BINDING EFFECTS.
C     *************************************************************
C
      IMPLICIT none
      REAL*8 VEC(6,2),PBEAM,WEIGHT,COSTR,TGEN,PGEN,PHIGEN,PWT(2,3)
      real*8 xtgt,targetpol,elmass,twoelm,em2,twopi,pi,alpha
      real*8    umin, ztgt, atgt, barn, csigml,tgtden,alum,zlum
      real*8  phirng, philim, dphspl,r,tgtmom,zzz,aaa,corfac,s0
      real*8  tumin, hbeta,rmin, r1,u1,r2,u2,x1,x2,s,phi,p1,p2
      real*8  theta1, theta2,r3,u3,r4,u4,x3,cos2st,sin2st,ds
      real*8  strfcn,delec,wtt,wzz,wxx,wyy,tx,ty,tz,tx1,ty1,arg,tz1
      real*8  tx2,ty2,tz2,costar,x4
      real*8 PPOL, PUNPOL, TGTPOL
*     REAL*4 RAN11,GAUSS1
      real*8 ran3
      real*8 climit, crange
      common /anglerange/ climit,crange
      LOGICAL BINDE,FCALL, scflag
      INTEGER II
      COMMON / RANSEED / II
      INTEGER NWPAWC, H
ccc   PARAMETER(NWPAWC= 320000)
      PARAMETER(NWPAWC= 5000000)
      COMMON /PAWC/ H(NWPAWC)

C     EXTERNAL RAN11
C
C     FCALL FLAGS THE FIRST CALL OF THE SUBROUTINE
      DATA FCALL/.TRUE./
      IF(FCALL) THEN
C
C     INITIALIZE EVERYTHING
C     THE MASS OF THE ELECTRON
         ELMASS=0.511D-3
         TWOELM=2.*ELMASS
         EM2=ELMASS**2
C
C     USEFUL CONSTANTS
         PI=3.1415927
         TWOPI=2.*PI
         ALPHA=1./137.
C
C     ZTGT,ATGT ARE THE ATOMIC NUMBER AND MASS OF THE TGT NUCLEUS
         ZTGT=26.43
         ATGT=57.14
C
C     BARN IS THE VALUE OF (HBAR*C)**2 IN GEV**2-CM**2
         BARN=0.389386E-27
C
C     THE CONSTANT OF THE M0LLER CROSS SECTION (IN THE CM FRAME)
C
         CSIGML=BARN*ALPHA**2
C
C     THE NUCLEAR LUMOSITY IS THE PRODUCT OF SEVERAL FACTORS,
C
C     THE TARGET MASS DENSITY (g/cm**3)
         TGTDEN=8.15
C
C     AND AVAGADRO'S NUMBER OVER THE MASS NUMBER (convert XTGT to cm)
         ALUM=TGTDEN*(XTGT/10.)*6.02E23/ATGT
C
C     THE LUMINOSITY OF ATOMIC ELECTRONS
         ZLUM=ZTGT*ALUM
C
C     FINALLY, IT'S NECESSARY TO DEFINE INTEGRATION LIMITS
C
C     THE LOWER LIMIT AND RANGE OF CM SCATTERING ANGLES
C
ccc   climit -0.4       crange 0.8      M.L. 12.94
ccc

C
C     THE LOWER LIMIT OF ELECTRON AZIMUTH AND RANGE
C
ccc
ccc   Has to be choosen larger than the detectorarea
ccc   e.g. +50% in all directions. M.L. 11.96
ccc
         PHIRNG=0.55            !M.L. 12.94
ccc
         PHIRNG=1.0             !
c     PHILIM=3.1416-PHIRNG/2.   !define appropriate phi space
c     PHILIM=4.7124-PHIRNG/2.
         PHILIM=2.9-PHIRNG/2.
C
C     NOW CALCULATE THE VOLUME OF PHASE SPACE ALLOWED BY THE LIMITS
C
         DPHSPL=PHIRNG*CRANGE
         FCALL=.FALSE.

      ENDIF
C
C     LET PBEAM BE THE INCIDENT ELECTRON ENERGY
C
      PBEAM=SQRT(VEC(4,1)**2+VEC(5,1)**2+VEC(6,1)**2)
      IF(BINDE) THEN
C
C     GENERATE THE TARGET ELECTRON CHARACTERISTICS
C
C     CHOOSE POLARIZED OR UNPOLARIZED, MOMENTUM, AND
C     KINEMATIC CORRECTION FACTOR
C
         R=RAN3(II)
c     tgtmom= 0.           !test
         IF(R.LE.targetpol) THEN
            CALL PEGEN(2,TGTMOM)
            zzz = RAN3(II)
            aaa = (-1.+2.*zzz)
            PPOL = tgtmom*1000000.d0*aaa
            call hfdp1(310,PPOL,1.0D0)
            CORFAC=1.+TGTMOM*aaa/ELMASS
            TGTPOL=1.
         ELSE
            CALL PEGEN(1,TGTMOM)
            zzz = RAN3(II)
            aaa = (-1.+2.*zzz)
            PUNPOL = tgtmom*1000000.d0*aaa
            call hfdp1(311,PUNPOL,1.0D0)
            CORFAC=1.+TGTMOM*aaa/ELMASS
            TGTPOL=0.
         ENDIF
      ELSE
         CORFAC=1.
         TGTPOL=1.
C     call hfdp1(310,tgtmom*1000000.,1.0D0)
C     call hfdp1(311,tgtmom*1000000.,1.0D0)
      ENDIF

C
C     GENERATE MOLLER SCATTERING (ADD VIRTUAL RAIDATION)
C
C     CHOOSE THE CM SCATTERING ANGLE
C
 1500 COSTAR=CLIMIT+CRANGE*RAN3(II)
c     if ((cmtheta .lt. 65.).or.(cmtheta .gt.110.)) cmtheta=65.
c     cmtheta= cmtheta + .1
c     costar= cosd(cmtheta)
C
C     The calculate the cm energy of the e-e- system (in GeV/c)
C
      S0=TWOELM*PBEAM*CORFAC
      call hfdp1(312,S0,1.0D0)
C
C     The correct scale for the bremsstrahlung is the minimum of T and U
C
      TUMIN=0.5*S0*(1.-ABS(COSTAR))
C
C     The constant HBETA is 1/2 of the bremsstrahlung constant beta
C
      HBETA=ALPHA/PI*(DLOG(TUMIN/EM2)-1.)
C
C     Define the minimum value of the photon energy fraction
C
      UMIN=1.d-30
C
C     Calculate the corresponding minimum random number (for U1,U2 gen)
C
      RMIN=UMIN**HBETA
C
C     Choose X1 and X2 to account for initial state bremsstrahlung
C
C     First, choose the photon energy fractions according to roughly the
C     correct distribution (U1,U2 are MUCH more important that X1,X2)
C
 1600 R1=RAN3(II)
      IF(R1.LT.RMIN) THEN
         U1=UMIN
      ELSE
         U1=R1**(1./HBETA)
      ENDIF
      R2=RAN3(II)
      IF(R2.LT.RMIN) THEN
         U2=UMIN
      ELSE
         U2=R2**(1./HBETA)
      ENDIF
C
C     Now convert them to X1,X2, and S
C
      if (scflag) then
         X1=1.-U1
         X2=1.-U2
      else
         x1= 1.0
         x2= 1.0
      endif
      S=S0*X1*X2
C
C     The cross section formally diverges at s=0, protect against
C
      IF(S.LT.(0.001**2)) GO TO 1600
C
C     Calculate the laboratory quantities from everything else
C
      PHI=PHILIM+PHIRNG*RAN3(II)


c--   If you wish to see true phi and theta symmetry, turn on the following stmt..
c     if (ran3(ii) .gt. 0.5) phi= phi - 3.1415927
c--   Otherwise, note that one detector is larger than the other so that
c--   accepted theta and phi distributions need not be symmetric.

c     phi= 4.7124
      P1=0.5*PBEAM*X1*(1.+COSTAR)
      P2=0.5*PBEAM*X1*(1.-COSTAR)
      THETA1=SQRT(CORFAC*TWOELM*X2*(1./P1-1./(X1*PBEAM)))
      THETA2=SQRT(CORFAC*TWOELM*X2*(1./P2-1./(X1*PBEAM)))
C
C     Allow for final state radiation
C
      R3=RAN3(II)
      IF(R3.LT.RMIN) THEN
         U3=UMIN
      ELSE
         U3=R3**(1./HBETA)
      ENDIF
      R4=RAN3(II)
      IF(R4.LT.RMIN) THEN
         U4=UMIN
      ELSE
         U4=R4**(1./HBETA)
      ENDIF
      if (scflag) then
         X3=1.-U3
         X4=1.-U4
      else
         x3= 1.0
         x4= 1.0
      endif
      P1=P1*X3
      P2=P2*X4
C
C     CALCULATE THE M0LLER CROSS SECTION
      COS2ST=COSTAR**2
      SIN2ST=1.-COS2ST
      DS=CSIGML/S*(3.+COS2ST)**2/SIN2ST**2
C
C     Next, include the electron structure functions (correct the
C     difference between the exact distributions and the approximate,
C     generated distributions
C
      if (scflag) then
         STRFCN=DELEC(U1,TUMIN)*U1**(1.-HBETA)/HBETA
     >        *DELEC(U2,TUMIN)*U2**(1.-HBETA)/HBETA
     >        *DELEC(U3,TUMIN)*U3**(1.-HBETA)/HBETA
     >        *DELEC(U4,TUMIN)*U4**(1.-HBETA)/HBETA
      else
         strfcn= 1.0
      endif
C
C     CALCULATE THE EVENT WEIGHT
      WEIGHT=ZLUM*DPHSPL*DS*STRFCN
C
C     CALCULATE THE EVENT WEIGHT FOR THE VARIOUS POLARIZATION DIRECTIONS
C
      WTT=TGTPOL*SIN2ST/(3.+COS2ST)**2
      WZZ=WTT*(7.+COS2ST)
      PWT(1,3)=WEIGHT*(1.-WZZ)
      PWT(2,3)=WEIGHT*(1.+WZZ)
C
C     THE TGT SPIN IS ASSUMED TO BE ALONG THE X AXIS
      WXX=WTT*SIN2ST*COS(2.*PHI)
      PWT(1,1)=WEIGHT*(1.-WXX)
      PWT(2,1)=WEIGHT*(1.+WXX)
      WYY=WTT*SIN2ST*SIN(2.*PHI)
      PWT(1,2)=WEIGHT*(1.-WYY)
      PWT(2,2)=WEIGHT*(1.+WYY)
C
C     SAVE THE CM SCATTERING ANGLE
      COSTR=COSTAR
C
C     STORE THE RESULTS OF THE EVENT GENERATION
C     (VEC(I,1) ALREADY CONTAINS THE INCIDENT BEAM DIRECTION)
C
      TX=VEC(4,1)/PBEAM
      TY=VEC(5,1)/PBEAM
      TZ=VEC(6,1)/PBEAM
C
C     ADD THE ELECTRON DIRECTION
C
      TX1=TX+THETA1*COS(PHI)
      TY1=TY+THETA1*SIN(PHI)
      ARG=1.-TX1**2-TY1**2
      IF(ARG.GT.0.) THEN
         TZ1=SQRT(ARG)
      ELSE
         TZ1=0.
      ENDIF
      TX2=TX+THETA2*COS(PHI+PI)
      TY2=TY+THETA2*SIN(PHI+PI)
      ARG=1.-TX2**2-TY2**2
      IF(ARG.GT.0.) THEN
         TZ2=SQRT(ARG)
      ELSE
         TZ2=0.
      ENDIF
      VEC(4,1)=P1*TX1
      VEC(5,1)=P1*TY1
      VEC(6,1)=P1*TZ1
      VEC(1,2)=VEC(1,1)
      VEC(2,2)=VEC(2,1)
      VEC(3,2)=VEC(3,1)
      VEC(4,2)=P2*TX2
      VEC(5,2)=P2*TY2
      VEC(6,2)=P2*TZ2
C
C     SAVE THE GENERATED LAB SCATTERING ANGLE AND MOMENTUM
C
      TGEN=THETA1
      PHIGEN=PHI
      PGEN=P1
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE BREMS(VEC,WEIGHT,THICK)
C
C     ENERGY LOSS DUE TO BREMSTRAHLUNG IN A PERMANDUR TARGET WITH
C     THICKNESS 'THICK' (RADIATION LENGTHS)
C     TSAI, REV. MOD. PHYS. 46, 815 (1974), PAGES 827-833
C
      IMPLICIT NONE
      REAL*8 VEC(6),WEIGHT,PBEAM,THICK,EXPONT,BT
      REAL*8 TX,TY,TZ,BEFF,LTK,Y,RANNUM
*     REAL*4 RAN11
*     real*8 lp, lw
      real*8  ran3
      INTEGER II
      COMMON / RANSEED / II
C     EXTERNAL RAN11
C     BEFF is Tsai's b; used equations 3.67, 3.68, 4.3 with Z = 26.43
      DATA BEFF /1.355 /
C
      PBEAM = DSQRT( VEC(4)**2 + VEC(5)**2 + VEC(6)**2 )
      TX = VEC(4) / PBEAM
      TY = VEC(5) / PBEAM
      TZ = VEC(6) / PBEAM
      LTK= THICK  / TZ
      BT=BEFF*LTK
C
C     GENERATE BREMSTRALUNG EVENTS
C
      RANNUM=RAN3(II)
      IF(RANNUM.GT.1.D-25) THEN
         EXPONT = LOG(RANNUM)/BT
         IF (EXPONT.GT.-50.) THEN
            Y=EXP(EXPONT)
         ELSE
            Y=0.
         ENDIF
      ELSE
         Y=0.
      ENDIF
      PBEAM = PBEAM*(1.-Y)
      WEIGHT=(1.-Y+0.75*Y**2)/(1./(1.+BT)+0.75*BT/(2.+BT))
C
      VEC(4) = TX*PBEAM
      VEC(5) = TY*PBEAM
      VEC(6) = TZ*PBEAM
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE BREMT(VEC,WEIGHT,THICK)
C
C     ENERGY LOSS DUE TO BREMSTRAHLUNG IN A PERMANDUR TARGET WITH
C     THICKNESS 'THICK' (RADIATION LENGTHS)
C     TSAI, REV. MOD. PHYS. 46, 815 (1974), EQUATION 3.84
C     GENERATE ELOSS THROUGH PHOTON GENERATION
C
      IMPLICIT NONE
      REAL*8 VEC(6),WEIGHT,THICK,P0
      REAL*8 RNORM,Y,YMIN
*     real*8 ZINT, rannum, expont
*     REAL*4 RAN11
      real*4 AVGGAM
      real*8  ran3
      INTEGER I,NGAM,IERROR
      LOGICAL FCALL
      INTEGER II
      COMMON / RANSEED / II
C     EXTERNAL RAN11
C     FCALL FLAGS THE FIRST CALL
      DATA FCALL/.TRUE./
      DATA YMIN/0.002D0/
C     FIRST CALL
      IF(FCALL) THEN
C
C     P0 IS THE NUMBER OF GAMMAS/RAD LENGTH
C
         P0=4./3.*(LOG(1./YMIN)-(1.-YMIN))+1./2.*(1.-YMIN**2)
C
C     RNORM IS THE NORMALIZATION FOR THE 1/Y GENERATION
C
         RNORM=LOG(1./YMIN)
         FCALL=.FALSE.
      ENDIF
C
C     INITIALIZE THE WEIGHT AND THE MEAN NUMBER OF PHOTONS TO GENERATE
C
      WEIGHT=1.0
      AVGGAM=P0*THICK
C
C     CHOOSE THE ACTUAL NUMBER OF PHOTONS TO GENERATE
C
      CALL POISSN(AVGGAM,NGAM,IERROR)
      IF(IERROR.NE.0) THEN
         PRINT 100, IERROR
 100     FORMAT(1X,'ERROR IN POISSON # GENRATION, IERROR = ',I5)
         STOP
      ENDIF
      IF(NGAM.GT.0) THEN
C
C     LOOP OVER ALL PHOTONS
C
         DO I=1,NGAM
C
C     GENERATE A PHOTON ENERGY ACCORDING TO A 1/Y DISTRIBUTION
C
            Y=EXP(-RAN3(II)*RNORM)
C
C     ADJUST THE ENERGY OF THE ELECTRON
C
            VEC(4) = (1.-Y)*VEC(4)
            VEC(5) = (1.-Y)*VEC(5)
            VEC(6) = (1.-Y)*VEC(6)
C
C     THE WEIGHT CORRECTS THE APPROXIMATE DISTRIBUTION TO THE DESIRED ONE
C
            WEIGHT=WEIGHT*RNORM/P0*(4./3.*(1.-Y)+Y**2)
         ENDDO
      ENDIF
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE BPHASE(VEC,PBEAM)
C
C     GENERATE BEAM PHASE SPACE DISTRIBUTRIONS
C
      IMPLICIT NONE
      REAL*8 VEC(6),PBEAM,RANXY,RANXYP,PHIXY,PHIXYP,PI
      REAL*8 SIGX,SIGXP,SIGY,SIGYP,TX,TY,RAD,TX0,TY0,X0,Y0
*     REAL*4 RAN11
      real*8 ran3
      INTEGER II
      COMMON / RANSEED / II

C     SIGMAS FOR THE BEAM SIZE AND DIVERGENCE AT THE M0LLER TARGET
C     and  POSITION AND ANGLE OFFSETS FOR THE BEAM
      common /beamparam/ SIGX,SIGY,SIGXP,SIGYP,X0,Y0,TX0,TY0

      DATA PI / 3.141593 /

C
      PHIXY = 2.*PI*RAN3(II)
      RANXY=1.-RAN3(II)
      IF(RANXY.GT.1.d-30) THEN
         RAD=SQRT(-2.*LOG(RANXY))
      ELSE
         RAD=0.
      ENDIF
      VEC(1)=X0+SIGX*RAD*COS(PHIXY)
      VEC(2)=Y0+SIGY*RAD*SIN(PHIXY)
      VEC(3)=0.
C
      PHIXYP=2.*PI*RAN3(II)
      RANXYP=1.-RAN3(II)
      IF(RANXYP.GT.1.d-30) THEN
         RAD=SQRT(-2.*LOG(RANXYP))
      ELSE
         RAD=0.
      ENDIF
      TX=TX0+SIGXP*RAD*COS(PHIXYP)
      TY=TY0+SIGYP*RAD*SIN(PHIXYP)
      VEC(4) = TX*PBEAM
      VEC(5) = TY*PBEAM
      VEC(6) = DSQRT( PBEAM**2 - VEC(4)**2 - VEC(5)**2 )
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE MSCAT(VEC,ZDIST,IMAT,TAILS)
C
C     MULTIPLE SCATTER THE PARTICLE OVER THE DISTANCE ZDIST
C     (IMAT=1 FOR FE, =2 FOR N2, =3 for al, =4 for He);  He added 15.04.94 (af)
C     TAILS SWITCHES ON MOLIERE SCATTERING
C
      IMPLICIT NONE
      REAL*8 VEC(6),ZDIST,PBEAM,PHIXY,PI
      REAL*8 TX,TY,TZ,THETA,SSIZE,X0(4),SIGT,RANN
*     REAL*4 RAN11
*     real*8 thet, ranxy, deltaz
      real*8 ran3
      LOGICAL TAILS
      INTEGER NSTEP,I,IMAT
      INTEGER II
      COMMON / RANSEED / II
      DATA PI / 3.141593 /
C
C     X0 stores the radiation lengths of iron, N2, AL and He
C
      DATA X0/17.6,3.04D5,89.0,5.30D6/
C
C     Check to insure that IMAT is valid
C
      IF(IMAT.LT.1.OR.IMAT.GT.4) THEN
         PRINT 5, IMAT
 5       FORMAT(/,1X,'MSCAT routine called with IMAT = ',I5)
         STOP
      ENDIF
C     C
C     C  DECIDE ON THE NUMBER OF 0.1% X0 STEPS TO TAKE
C     C
C     NSTEP=INT(1000.*ZDIST/X0(IMAT))+1
C
C     DECIDE ON THE NUMBER OF 10% X0 STEPS TO TAKE
C
      NSTEP=INT(10.*ZDIST/X0(IMAT))+1
      SSIZE=ZDIST/NSTEP
C
C     CALCULATE THE ELECTRON MOMENTUM
C
      PBEAM=SQRT(VEC(4)**2+VEC(5)**2+VEC(6)**2)
      SIGT=0.014/PBEAM*SQRT(SSIZE/X0(IMAT))
C
C     INITIALIZE THE PROJECTED X,Y ANGLES
C
      TX=VEC(4)/PBEAM
      TY=VEC(5)/PBEAM
C
C     LOOP OVER ALL STEPS
C
      DO I=1,NSTEP
         PHIXY=2.*PI*RAN3(II)
         IF(TAILS) THEN
            CALL MOLIERE(IMAT,PBEAM,SSIZE,THETA)
         ELSE
C
C     IF MOLIERE TAILS NOT SELECTED, GENERATE A GAUSSIAN DIST
C
            RANN=RAN3(II)
            IF(RANN.GT.1.d-30) THEN
               THETA=SIGT*SQRT(-2.*LOG(RANN))
            ELSE
               THETA=0.
            ENDIF
         ENDIF
         TX=TX+THETA*COS(PHIXY)
         TY=TY+THETA*SIN(PHIXY)
         TZ=SQRT(1.-TX**2-TY**2)
C
C     UPDATE THE POSITION VECTOR
C
         VEC(3)=VEC(3)+SSIZE
         VEC(1)=VEC(1)+TX/TZ*SSIZE
         VEC(2)=VEC(2)+TY/TZ*SSIZE
      ENDDO
C
C     UPDATE THE MOMENTUM VECTOR
C
      VEC(4)=TX*PBEAM
      VEC(5)=TY*PBEAM
      VEC(6)=TZ*PBEAM
C
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE BRMSCT(VEC,ZDIST,IMAT,TAILS,TCKBRM,WEIGHT)
C
C     MULTIPLE SCATTER AND BREMSSTRAHLUNG ENERGY LOSS OVER DISTANCE ZDIST
C     (IMAT=1 FOR FE, =2 FOR N2, =3 for al, =4 for He);  He added 15.04.94 (af)
C     TAILS SWITCHES ON MOLIERE SCATTERING
C     TCKBRM SELECTS THE BREMSSTRAHLUNG APPROXIMATION
C     WEIGHT IS SET TO 0 IF THE PRIMARY E- FALLS BELOW 5 GEV
C
      IMPLICIT NONE
      REAL*8 VEC(6),ZDIST,PBEAM,PHIXY,PI,WEIGHT,THICK,WGTBRM
      REAL*8 TX,TY,TZ,THETA,SSIZE,X0(4),SIGT,RANN,P0,PAVG
*     real*8 thet,ranxy,deltaz
*     REAL*4 RAN11
      real*8 ran3
      LOGICAL TAILS, TCKBRM
      INTEGER NSTEP,I,IMAT
      INTEGER II
      COMMON / RANSEED / II
      DATA PI / 3.141593 /

C
C     X0 stores the radiation lengths of iron, N2, and AL
C
      DATA X0/17.6,3.04D5,89.0,5.30D6/
C
C     Check to insure that IMAT is valid
C
      IF(IMAT.LT.1.OR.IMAT.GT.4) THEN
         PRINT 5, IMAT
 5       FORMAT(/,1X,'MSCAT routine called with IMAT = ',I5)
         STOP
      ENDIF
C
C     INITIALIZE THE EVENT WEIGHT
C
      WEIGHT=1.0
C
C     DECIDE ON THE NUMBER OF 3.3% X0 STEPS TO TAKE
C
      NSTEP=INT(30.*ZDIST/X0(IMAT))+1
      SSIZE=ZDIST/NSTEP
      THICK = SSIZE/X0(IMAT)
C
C     CALCULATE THE ELECTRON MOMENTUM
C
      PBEAM=SQRT(VEC(4)**2+VEC(5)**2+VEC(6)**2)
C
C     INITIALIZE THE PROJECTED X,Y ANGLES
C
      TX=VEC(4)/PBEAM
      TY=VEC(5)/PBEAM
C
C     LOOP OVER ALL STEPS
C
      DO I=1,NSTEP
C
C     STORE THE INITIAL MOMENTUM
C
         P0=PBEAM
C
C     BREMSTRAHLUNG
C
         IF(TCKBRM) THEN
            CALL BREMS(VEC,WGTBRM,THICK)
         ELSE
            CALL BREMT(VEC,WGTBRM,THICK)
         ENDIF
         WEIGHT=WEIGHT*WGTBRM
C
C     CALCULATE THE ELECTRON MOMENTUM AT THE END OF THE STEP
C
         PBEAM=SQRT(VEC(4)**2+VEC(5)**2+VEC(6)**2)
C
C     CALCULATE THE AVERAGE ELECTRON MOMENTUM DURING THE STEP
C
         PAVG=(P0+PBEAM)/2.
C
C     CHOOSE A RANDOM PHI COORDINATE
C
         PHIXY=2.*PI*RAN3(II)
         IF(TAILS) THEN
            CALL MOLIERE(IMAT,PAVG,SSIZE,THETA)
         ELSE
C
C     IF MOLIERE TAILS NOT SELECTED, GENERATE A GAUSSIAN DIST
C


            SIGT=0.014/PAVG*SQRT(THICK)
            RANN=RAN3(II)
            IF(RANN.GT.1.d-30) THEN
               THETA=SIGT*SQRT(-2.*LOG(RANN))
            ELSE
               THETA=0.
            ENDIF
         ENDIF
         TX=TX+THETA*COS(PHIXY)
         TY=TY+THETA*SIN(PHIXY)
         TZ=SQRT(1.-TX**2-TY**2)
C
C     UPDATE THE POSITION VECTOR
C
         VEC(3)=VEC(3)+SSIZE
         VEC(1)=VEC(1)+TX/TZ*SSIZE
         VEC(2)=VEC(2)+TY/TZ*SSIZE
      ENDDO
C
C     UPDATE THE MOMENTUM VECTOR
C
      VEC(4)=TX*PBEAM
      VEC(5)=TY*PBEAM
      VEC(6)=TZ*PBEAM
C
C

      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE MOLIERE(IMAT,PBEAM,THICK,THETA)
C
C     ********************************************************************
C     * This routine generates a scattering angle according to the       *
C     * Moliere distribution.                                            *
C     * Parameters: IMAT -- flags material as Fe(1) N2(2) AL(3) or He(4) *
C     *             PBEAM - the energy of the electron in GeV            *
C     *             THICK - the thickness of the scatterer (in mm)       *
C     *             THETA - the scattered (space angle) of the e-        *
C     ********************************************************************
C
      IMPLICIT none

      real*8 pbeam,thick,theta

C     EXTERNAL PFUNC
C
C     NETA is the number of ETA bins to be stored
C
      integer*4 neta ,imat
      PARAMETER (NETA=120)
C
C     Define a COMMON block to store the necessary information for the
C     external function PFUNC: the value of the parameter GAMMA
C
      real*8 gamma
      COMMON /PFPAR/ GAMMA
*     REAL*4 ALIM,BLIM,EPS
*     real*4 r, x0
C
C     The array ETAPRB(100,4) stores the probability distribution of the
C     variable eta for each of the four materials
C
      real*8 ETAPRB(0:NETA,4),SETA(3,NETA)
      real*8 A(4),Z(4),TOLD(4),RHO(4),B(4),T1CON(4),GAMCON(4)
*     REAL*4 RAN11
      real*8  ran3
      LOGICAL FCALL

      INTEGER II
      COMMON / RANSEED / II

      integer*4 i,ird,j
      real*8 qb,ql,del,etamax,dx,erand,eta0,p0,p1,eta
      real*8 theta1


      DATA FCALL/.TRUE./
      DATA A/55.85,14.01,26.98,4.0/,Z/26.,7.,13.,2./,
     +     RHO/7.87,1.25D-3,2.70,0.178D-3/
C
C     Check to insure that IMAT is valid
C
      IF(IMAT.LT.1.OR.IMAT.GT.4) THEN
         PRINT 5, IMAT
 5       FORMAT(/,1X,'Moliere routine called with IMAT = ',I5)
         STOP
      ENDIF
      IF(FCALL) THEN
C
C     INITIALIZE EVERYTHING
C
C     Determine the constants for the theta1's and the B values
C
         DO J=1,4
            QB=SQRT(Z(J)*(Z(J)+1.))
            QL=(Z(J)+1.)*Z(J)**(1./3.)
            T1CON(J)=0.3965D-3*QB*SQRT(RHO(J)/A(J))
            DEL=1.13
            GAMCON(J)=8.831D3*QL*RHO(J)/(A(J)*DEL)
         ENDDO
         OPEN(21,FILE='TAILFCN',STATUS='UNKNOWN')
C
C     Fill the array SETA
C

         READ(21,*,IOSTAT=IRD) SETA
         IF(IRD.NE.0) THEN
            PRINT 26
 26         FORMAT(1H0,'Read of SETA unsuccessful, generating')
C
C     READ unsuccessful, generate the data
C
            DO I=1,NETA-1
               ETAMAX=0.10*I
               CALL SUMETA(ETAMAX,SETA(1,I))
            ENDDO
            SETA(1,NETA)=1.0
            SETA(2,NETA)=0.0
            SETA(3,NETA)=0.
C
C     WRITE out data for future use
C
            REWIND 19
            WRITE(21,*) SETA
         ENDIF
         CLOSE(21)
C
C     Renormalize the distribution to have unit probability at X=10.
C
         DO J=1,4
            ETAPRB(0,J)=0.0
            ETAPRB(NETA,J)=1.0
         ENDDO
C
C     Set the old thicknesses to 0.
C
         DO J=1,4
            TOLD(J)=0.0
         ENDDO
         PRINT 40
 40      FORMAT(/,1X,' Moliere intitialization complete ')
         FCALL=.FALSE.
      ENDIF
C
C     Check to see of the routine has been initialized for this thickness
C
      DX=THICK/10.
      IF(DX.NE.TOLD(IMAT)) THEN
C
C     Calculate the value of B
C
         GAMMA=GAMCON(IMAT)*DX
         IF(GAMMA.LT.10.) GAMMA=10.
         B(IMAT)=1.153+2.583*LOG10(GAMMA)
C
C     Fill the array ETAPRB
C
         DO I=1,NETA-1
            ETAPRB(I,IMAT)=SETA(1,I)+SETA(2,I)/B(IMAT)
     >           +SETA(3,I)/B(IMAT)**2
         ENDDO
C     PRINT 45, ETAPRB
C     45      FORMAT(1X,6F11.9)
C
C     Remember the current thickness
C
         TOLD(IMAT)=DX
      ENDIF
C
C     Next, generate a random number
C
      ERAND=RAN3(II)
C
C     Next, interpolate a value of ETA such that the integral probability
C     is equal to ERAND
C
      DO I=1,NETA
         IF(ETAPRB(I,IMAT).GE.ERAND) THEN
            ETA0=0.10*(I-1)
            P0=ETAPRB(I-1,IMAT)
            P1=ETAPRB(I,IMAT)
            GO TO 100
         ENDIF
      ENDDO
      PRINT 50, ERAND
 50   FORMAT(1H0,'random ETA prob = ',E12.6,' not matched')
      THETA=0.
      RETURN
 100  ETA=ETA0+0.10*(ERAND-P0)/(P1-P0)
C     CALL HFDP1(20,X,1.)
C
C     Finally, convert it into a space angle
C
      THETA1=T1CON(IMAT)*SQRT(DX*B(IMAT))/PBEAM
      THETA=ETA*THETA1
      RETURN
      END
C
      SUBROUTINE SUMETA(ETAMAX,SETA)
C
C     ********************************************************************
C     * This function calculates the integrals of the components of the  *
C     * Moliere scattering distribution.  ETAMAX is the upper limit of   *
C     * the integrals                                                    *
C     ********************************************************************
C
      IMPLICIT none
      real*8 SETA(3), XSUM(3)
C
C     Define the number of steps
C
      real*8 etamax,xstep,xint1,xint2,fscat,xint3,coeff,x
      integer*4 nstep,nsum,i,istep

      IF(ETAMAX.LT.1.0) NSTEP=20
      IF(ETAMAX.GE.1.0) NSTEP=INT(20.*ETAMAX)
C
C     There must be an even number of steps
C
      IF(NSTEP/2*2.NE.NSTEP) NSTEP=NSTEP+1
C
C     XSTEP IS THE STEP SIZE IN X
C
      XSTEP=ETAMAX/NSTEP
C
C     XSUM IS THE INTEGRAL
C
      DO I=1,3
         XSUM(I)=0.
      ENDDO
      NSUM=NSTEP+1
      DO ISTEP=1,NSUM
         X=(ISTEP-1)*XSTEP
C
C     CALCULATE THE INTEGRAND
C
C     XINT IS THE INTEGRAND
         XINT1=X*2.*EXP(-X**2)
         XINT2=X*FSCAT(1,X)
         XINT3=X*FSCAT(2,X)
C
C     DETERMINE THE COEFFICIENT OF THE INTEGRAND (ACCORDING TO SIMPSON'S
C     RULE)
         IF(ISTEP.EQ.1.OR.ISTEP.EQ.NSUM) THEN
            COEFF=1.0
         ELSE
            COEFF=4.0
            IF((ISTEP/2)*2.NE.ISTEP) COEFF=2.0
         ENDIF
C
C     SUM THE INTEGRAND
         XSUM(1)=XSUM(1)+COEFF*XINT1
         XSUM(2)=XSUM(2)+COEFF*XINT2
         XSUM(3)=XSUM(3)+COEFF*XINT3
      ENDDO
C
C     NOW MULTIPLY BY THE NORMALIZATION AND THE BIN WIDTH
C
      DO I=1,3
         SETA(I)=XSTEP/3.*XSUM(I)
      ENDDO
      RETURN
      END
      FUNCTION FSCAT(K,ETA)
C
C     ********************************************************************
C     * This function calculates the Moliere F functions of order K.     *
C     ********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER K
      REAL*8 KFACT
C
C     Calculate K factorial
C
      KFACT=1.
      IF(K.GT.0) THEN
         DO I=1,K
            KFACT=KFACT*I
         ENDDO
      ENDIF
C
C     Define the integration limits and number of steps
C
      NSTEP=100
      XMAX=10.
C
C     XSTEP IS THE STEP SIZE IN X
C
      XSTEP=XMAX/NSTEP
C
C     XSUM IS THE INTEGRAL
C
      XSUM=0.
      NSUM=NSTEP+1
      DO ISTEP=1,NSUM
         X=(ISTEP-1)*XSTEP
C
C     CALCULATE THE INTEGRAND
C
C     XINT IS THE INTEGRAND
         Y=X**2/4.
         IF(Y.GT.0.) THEN
            Z=Y*LOG(Y)
         ELSE
            Z=0.
         ENDIF
         XINT=DBeSJ0(ETA*X)*EXP(-Y)*X*Z**K
C
C     DETERMINE THE COEFFICIENT OF THE INTEGRAND (ACCORDING TO SIMPSON'S
C     RULE)
         IF(ISTEP.EQ.1.OR.ISTEP.EQ.NSUM) THEN
            COEFF=1.0
         ELSE
            COEFF=4.0
            IF((ISTEP/2)*2.NE.ISTEP) COEFF=2.0
         ENDIF
C
C     SUM THE INTEGRAND
         XSUM=XSUM+COEFF*XINT
      ENDDO
      XSUM=XSUM/KFACT
C
C     NOW MULTIPLY BY THE NORMALIZATION AND THE BIN WIDTH
C
      FSCAT=XSTEP/3.*XSUM
      RETURN
      END
C
      FUNCTION DELEC(U,S)
C
C     **********************************************************
C     This routine calculates the electron structure function
C     as given in Alexander et al., Phys.Rev.D37, 56 (1988)
C     U is 1-X...for precision reasons.
C     **********************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 L,L1,ME,ME2
C
C     Define a bunch of physical constants
C
C
C     The electron mass and the square of the electron mass
C
      ME=0.5110034E-3
      ME2=ME**2
C
C     Alpha and pi
C
      ALPHA=1./137.036
      PI=3.141593
C
C     If X is too close to 1, reset it to a more sensible value
C
      IF(U.LT.1.e-30) U=1.e-30
C
C     Calculate several useful auxillary quantities
C
      L=DLOG(S/ME2)
      L1=L+2.*DLOG(U)
      BETA=2.*ALPHA/PI*(L-1.)
      EBEAM=SQRT(S)/2.
C
C     Calculate e-photon structure function
C
      DEPHOT=BETA/2.*U**(BETA/2.-1.)
     >     *(1.+3.*BETA/8.-BETA**2/48.*(L/3.+PI**2-47./8.))
     >     -BETA/4.*(2.-U)+BETA**2/32.
     >     *(4.*(2.-U)*DLOG(1./U)-(1.+3.*(1.-U)**2)/U*DLOG(1.-U)-6.+U)
C
C     Calculate e-e pair structure function
C
      IF(U.LE.2.*ME/EBEAM) THEN
C
C     Insufficient energy to make an electron pair
C
         DEEE=0.
      ELSE
         DEEE=(ALPHA/PI)**2*((U-2.*ME/EBEAM)**(BETA/2.)/U
     >        *(L1-5./3.)**2/12.*(1.+(1.-U)**2+BETA*(L1-5./3.)/6.)
     >        +L**2/4.*(2./3.*(1.-(1.-U)**3)/(1.-U)+U/2.+(2.-U)*DLOG(1.-U)))
      ENDIF
C
C     Sum both contributions to the electron structure function
C
      DELEC=DEPHOT+DEEE
      RETURN
      END
      SUBROUTINE PLANEX(VEC,PLANE)
C
C     **********************************************************************
C     THIS ROUTINE TRANLATES THE TRAJECTORY 6-VECTOR UNTIL IT INTERSECTS
C     A PLANE SPECIFIED BY THE SIX-VECTOR PLANE.  THE PLANE VECTOR CONTAINS
C     THE COORDINATES OF A POINT ON THE PLANE XP,YP,ZP FOLLOWED BY THE
C     DIRECTION COSINES OF A NORMAL TO THE PLANE.
C     **********************************************************************
C
      IMPLICIT none
      real*8 VEC(6),PLANE(6),DIRCOS(3)
      real*8 rnorm,ptot,a,b,zint,xint,yint
      integer*4 i
C
C     INSURE THAT THE PLANE NORMAL DIR COSINES ARE PROPERLY NORMALIZED
      RNORM=SQRT(PLANE(4)**2+PLANE(5)**2+PLANE(6)**2)
      DO 100 I=4,6
         PLANE(I)=PLANE(I)/RNORM
 100  CONTINUE
C
C     CALCULATE THE DIR COSINES FOR THE PARTICLE TRAJECTORY
      PTOT=SQRT(VEC(4)**2+VEC(5)**2+VEC(6)**2)
      DO 200 I=1,3
         DIRCOS(I)=VEC(3+I)/PTOT
 200  CONTINUE
C
C     CALCULATE SEVERAL INTERMEDIATE QUANTITIES
      A=0.
      B=0.
      DO 300 I=1,3
         A=A+PLANE(I+3)*(VEC(I)-PLANE(I))
         B=B+PLANE(I+3)*DIRCOS(I)
 300  CONTINUE
      B=B/DIRCOS(3)
C
C     CALCULATE THE Z INTERCEPT
C     ZINT=(B*VEC(3)+PLANE(6)*PLANE(3)-A)/(B+PLANE(6))
      ZINT=VEC(3)-A/B
C
C     CALCULATE THE X AND Y INTERCEPTS
      XINT=VEC(1)+DIRCOS(1)/DIRCOS(3)*(ZINT-VEC(3))
      YINT=VEC(2)+DIRCOS(2)/DIRCOS(3)*(ZINT-VEC(3))
C
C     NOW UPDATE THE POSITION COORDINATES OF VEC
      VEC(1)=XINT
      VEC(2)=YINT
      VEC(3)=ZINT
      RETURN
      END
C
C     CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE PEGEN(IDIST,PEL)
C
C     ********************************************************************
C     * This routine generates the momenta of electrons in the Fe atom.  *
C     * Parameters: IDIST - indicates unpolarized electron (1) or        *
C     *                     polarized electron (2)                       *
C     *             PEL --- the momentum of the electron in GeV/c        *
C     ********************************************************************
C
      IMPLICIT none

      integer*4 idist,npel,i,imat
      real*8 pel,pelmax,pelprb,erand,pel0,p0,p1


C
C     NPEL is the number of PEL bins to be stored
C
      PARAMETER (NPEL=150)
C
C     The array PELPRB(100,3) stores the probability distribution of the
C     variable PEL for unpolarized and polarized electrons
C
      DIMENSION PELPRB(2,0:NPEL)
*     REAL*4 RAN11
      real*8  ran3
      LOGICAL FCALL
      INTEGER II
      COMMON / RANSEED / II
      DATA FCALL/.TRUE./
C
C     Check to insure that IDIST is valid
C
      IF(IDIST.LT.1.OR.IDIST.GT.2) THEN
         PRINT 5, IMAT
 5       FORMAT(/,1X,'PEGEN routine called with IDIST = ',I5)
         STOP
      ENDIF
      IF(FCALL) THEN
C
C     INITIALIZE EVERYTHING
C
C
C     Calculate the probability distribution
C
         DO I=1,NPEL
            PELMAX=2.00*I
            CALL SUMMOM(PELMAX,PELPRB(1,I))
         ENDDO
         PELPRB(1,0)=0.
         PELPRB(2,0)=0.
         DO I=1,NPEL
            PELPRB(1,I)=PELPRB(1,I)/PELPRB(1,NPEL)
            PELPRB(2,I)=PELPRB(2,I)/PELPRB(2,NPEL)
         ENDDO
         PELPRB(1,NPEL)=1.0
         PELPRB(2,NPEL)=1.0
         DO I=0,NPEL
c     PRINT 30, 2.*I, PELPRB(1,I), PELPRB(2,I)
c     write(8,30) 2.*I, PELPRB(1,I), PELPRB(2,I)
 30         FORMAT(1X,F6.2,', ',E11.4,', ',E11.4)
         ENDDO
         PRINT 40
 40      FORMAT(/,1X,' PEGEN intitialization complete ')
         FCALL=.FALSE.
      ENDIF
C
C     Next, generate a random number
C
      ERAND=RAN3(II)
C
C     Next, interpolate a value of PEL such that the integral probability
C     is equal to ERAND
C
      DO I=1,NPEL
         IF(PELPRB(IDIST,I).GE.ERAND) THEN
            PEL0=2.00*(I-1)
            P0=PELPRB(IDIST,I-1)
            P1=PELPRB(IDIST,I)
            GO TO 100
         ENDIF
      ENDDO
      PRINT 50, ERAND
 50   FORMAT(1H0,'random PEL prob = ',E12.6,' not matched')
      PEL=0.
      RETURN
 100  PEL=PEL0+2.00*(ERAND-P0)/(P1-P0)
C
C     Renormalize to GeV/c
C
      PEL=PEL*1.D-6
      RETURN
      END
C
C
      SUBROUTINE SUMMOM(PELMAX,SUMP)
C
C     ********************************************************************
C     * This function calculates the integrals of the electron momentum  *
C     * distributions for the Fe atom.  PELMAX is the upper limit of the *
C     * integrals                                                        *
C     ********************************************************************
C
      IMPLICIT none
      real*8  SUMP(2), XSUM(2), REFMOM(8)

      integer*4 nstep,i,nsum,istep
      real*8 pelmax,xstep,x,p1,p2,p3
      real*8 p4,p5,p6,p7,p8,xtmp1,xtmp2,xint1,xint2,coeff

C
C     REFMOM is the conversion factor from dimensionless momentum to
C     KeV/c for each successive shell (assume screening)
C
C     DATA REFMOM/97.0,89.5,59.7,7.46/
      DATA REFMOM/95.1,76.5,35.4,5.60,98.8,80.2,37.3,5.60/
C
C     Define the number of steps
C
      IF(PELMAX.LT.5.0) NSTEP=20
      IF(PELMAX.GE.5.0) NSTEP=INT(4.*PELMAX)
C
C     There must be an even number of steps
C
      IF(NSTEP/2*2.NE.NSTEP) NSTEP=NSTEP+1
C
C     XSTEP IS THE STEP SIZE IN X
C
      XSTEP=PELMAX/NSTEP
C
C     XSUM IS THE INTEGRAL
C
      DO I=1,2
         XSUM(I)=0.
      ENDDO
      NSUM=NSTEP+1
      DO ISTEP=1,NSUM
         X=(ISTEP-1)*XSTEP
C
C     CALCULATE THE INTEGRAND
C
C     XINT1 IS THE UNPOLARIZED DISTRIBUTION
C
         P1=(X/REFMOM(1))**2
         P2=(X/REFMOM(2))**2
         P3=(X/REFMOM(3))**2
         P4=(X/REFMOM(4))**2
         P5=(X/REFMOM(5))**2
         P6=(X/REFMOM(6))**2
         P7=(X/REFMOM(7))**2
         P8=(X/REFMOM(8))**2
         XTMP1=2.*10.1859*P1/(1.+P1)**4/REFMOM(1)
     >        +(2.*325.949*P2*(-1.+4.*P2)**2+6.*1738.4*P2**2)
     >        /(1.+4.*P2)**6/REFMOM(2)
     >        +(2.*275.02*P3*(-1.+4.*(-1.+9.*P3)**2/(1.+9.*P3)**2)**2
     >        /(1.+9.*P3)**4+(6.*79205.7*P3**2*(-1.+9.*P3)**2
     >        +3.79*570281.*P3**3)/(1.+9.*P3)**8)/REFMOM(3)
     >        +2.*651.899*P4*(-1.+16.*P4)**2*(-4.
     >        +8.*(-1.+16.*P4)**2/(1.+16.*P4)**2)**2/(1.+16.*P4)**6
     >        /REFMOM(4)
         XTMP2=2.*10.1859*P5/(1.+P5)**4/REFMOM(5)
     >        +(2.*325.949*P6*(-1.+4.*P6)**2+6.*1738.4*P6**2)
     >        /(1.+4.*P6)**6/REFMOM(6)
     >        +(2.*275.02*P7*(-1.+4.*(-1.+9.*P7)**2/(1.+9.*P7)**2)**2
     >        /(1.+9.*P7)**4+(6.*79205.7*P7**2*(-1.+9.*P7)**2
     >        +4.705*570281.*P7**3)/(1.+9.*P7)**8)/REFMOM(7)
     >        +2.*651.899*P8*(-1.+16.*P8)**2*(-4.
     >        +8.*(-1.+16.*P8)**2/(1.+16.*P8)**2)**2/(1.+16.*P8)**6
     >        /REFMOM(8)
         XINT1=(XTMP1+XTMP2)/48.495
         XINT2=(570281.*P3**3/(1.+9.*P3)**8/REFMOM(3)
     >        +570281.*P7**3/(1.+9.*P7)**8/REFMOM(7))/2.
C
C     DETERMINE THE COEFFICIENT OF THE INTEGRAND (ACCORDING TO SIMPSON'S
C     RULE)
         IF(ISTEP.EQ.1.OR.ISTEP.EQ.NSUM) THEN
            COEFF=1.0
         ELSE
            COEFF=4.0
            IF((ISTEP/2)*2.NE.ISTEP) COEFF=2.0
         ENDIF
C
C     SUM THE INTEGRAND
         XSUM(1)=XSUM(1)+COEFF*XINT1
         XSUM(2)=XSUM(2)+COEFF*XINT2
      ENDDO
C
C     NOW MULTIPLY BY THE NORMALIZATION AND THE BIN WIDTH
C
      DO I=1,2
         SUMP(I)=XSTEP/3.*XSUM(I)
      ENDDO
      RETURN
      END

c__________________________________________________________________________
c**********************************************************************
c***  ***
c***  Subroutine Collimator                                         ***
c***  ***
c**********************************************************************
c
      subroutine Collimate(vec, singlemol, twomols, ienergy)
c     -------------------------------------------------------
c
c     Checks whether a Moller electron passes through one of the collimator
c     slits in front of the second quadrupole, at z = ~2150 mm. The slit
c     dimensions left and right are assumed to be identical but energy
c     dependent. The integer variable ienergy can take values from 1 to 6.
c
c     5-8-94
c     dB
      implicit none
c
      real*8  vec(6)
      logical singlemol, twomols
      integer ienergy
c
      real*8  xmin(6), xmax(6), ymax(6)
c
c     data xmin/40.0,32.5,29.0,26.5,24.5,23.0/
c     data xmax/48.5,41.5,37.0,33.5,31.0,29.5/
c     data ymax/25.,20.,15.,12.,11.,10./
c
      data xmin/30.,5*0./
      data xmax/52.,5*0./
      data ymax/44.,5*0./
c
c     data xmin/40.0,32.5,29.0,22.5,24.5,23.0/
c     data xmax/48.5,41.5,37.0,37.5,31.0,29.5/  !S.Danagoulian 09/00
c     data ymax/25.,20.,15.,16.,11.,10./
c

      ienergy = 1
c     if ((ienergy .lt. 1) .or. (ienergy .gt. 6)) stop  !done in main routine
c     *'Incident electron energy out of range in subroutine Collimator'
c
      if ((abs(vec(1)) .lt. xmin(ienergy)) .or.
     *     (abs(vec(1)) .gt. xmax(ienergy)) .or.
     *     (abs(vec(2)) .gt. ymax(ienergy))) then
         singlemol = .false.
         twomols = .false.
      endif
c
      return
      end



c**********************************************************************
c***  ***
c***  Function all_coll                                             ***
c***  ***
c**********************************************************************
c
      logical function all_coll(xy)
c     -------------------------------------------------------
c
c     Using the coordinates provided at the actual z positions of
c     the collimators, this function returns true (pass) or false (blocked)
c     to indicate whether the electron passed through the numbered collimator
c     or was blocked by it.
c
c     hcf 2/6/01
      implicit none
c
      real*8  xy(2,7)           !contains (x,y) pairs at each of seven collimator Z's
      integer icol

      real*8 xx,yy
      common /collimators/ pcol,colxmin,colxmax,colymin,colymax
      real*8 pcol(7)
      real*8 colxmin(7),colxmax(7),colymin(7),colymax(7)
c     data xmin / -999.0, -999.0,  -49.0, -999.0, -999.0, -999.0,   48.0 /
c     data xmax /  999.0,  999.0,  999.0,   49.0,  999.0,  -48.0,  999.0 /
c     data ymin / -999.0,  -50.0, -999.0, -999.0, -999.0, -999.0, -999.0 /
c     data ymax /   50.0,  999.0,  999.0,  999.0,  999.0,  999.0,  999.0 /

      all_coll= .true.

c--   In calling program we normally only fill 1,3,5,6, so fill others here
      xy(1,2)= xy(1,1)          !since collimators 1 and 2 are at same z
      xy(2,2)= xy(2,1)
      xy(1,4)= xy(1,3)          !since collimators 3 and 4 are at same z
      xy(2,4)= xy(2,3)
      xy(1,7)= xy(1,6)          !since collimators 6 and 7 are at same z
      xy(2,7)= xy(2,6)

      do icol=1,5
         xx= xy(1,icol)
         yy= xy(2,icol)
         if ( (xx.lt.colxmin(icol)) .or. (xx.gt.colxmax(icol)) .or.
     $        (yy.lt.colymin(icol)) .or. (yy.gt.colymax(icol)) ) all_coll=.false.
      end do
      xx= xy(1,6)
      if ( (xx.gt.colxmax(6)) .and. (xx.lt.colxmin(7)) ) all_coll=.false.


      return
      end

c**********************************************************************
c
      subroutine Collimator(vec, singlemol, twomols, ienergy)
c     -------------------------------------------------------
c
c     Checks whether a Moller electron passes through one of the collimator
c     slits in front of the second quadrupole, at z = ~2150 mm. The slit
c     dimensions left and right are assumed to be identical but energy
c     dependent. The integer variable ienergy can take values from 1 to 6.
c
c     5-8-94
c     dB
      implicit none
c
      real*8  vec(6)
      logical singlemol, twomols
      integer ienergy
c
      real*8  xmin(6), xmax(6), ymax(6)
c
      data xmin/40.0,32.5,29.0,26.5,24.5,23.0/
      data xmax/48.5,41.5,37.0,33.5,31.0,29.5/
      data ymax/25.,20.,15.,12.,11.,10./
c
      if ((ienergy .lt. 1) .or. (ienergy .gt. 6)) stop
     *     'Incident electron energy out of range in subroutine Collimator'
c
      if ((abs(vec(1)) .lt. xmin(ienergy)) .or.
     *     (abs(vec(1)) .gt. xmax(ienergy)) .or.
     *     (abs(vec(2)) .gt. ymax(ienergy))) then
         singlemol = .false.
         twomols = .false.
      endif
c     print *,' collimator'
c
      return
      end


c__________________________________________________________________________
c**********************************************************************
c***  ***
c***  Subroutine Detectors                                          ***
c***  ***
c**********************************************************************
c
      Subroutine detectors(vec,padtl,padtl2,blktl,blktl2,
     +     weight,weights,pwt,wgtcx,singlemol,twomols,avgA,avgA2)
c     --------------------------------------------
c
c     This routine simulates the detectors and properly weights which
c     dectectors were hit.                              af. 25.04.94
c
c     Checks whether a Moller electron falls into one of the two trapezoids
c     composing the detector system. The acceptances of the left and right
c     detectors need not be identical.
c     5-8-94
c     dB
      implicit none
c
      integer npads, j, K, i, ileft, iright, isum, idif,iasy
*     integer l, kblk(2)
      parameter (npads=1)
      REAL*8 VEC(6,2),padtl(2,3,npads,2),padtl2(2,3,npads,2)
      real*8 BLKTL(2,3),BLKTL2(2,3)
      real*8 WEIGHT,PWT(2,3),weights(2),wgtcx, wgt
      real*8 LOCALWT(2,3)
      real*8 AZ_local,avgA,avgA2
      logical twomols, singlemol(2)
ccc
c     real*4            ntuple(20)      ! ntuple(9..14) used for Assymetries
c     common/cern/      ntuple          ! calculation
ccc
c
      real*8 x, xleftmin, xleftmax, xrightmin, xrightmax
      real*8 xleftoffset,xrightoffset
      real*8 xhodomin, dxhodo
      real*8 xlefthodomin,xrighthodomin
      real*8  y
*     real*8 yleftmin, yleftmax, yrightmin, yrightmax
      logical left, right
c
c
c     data xleftmin/-430.0/             ! real setup M.L.Oct.95
c     data xleftmax/-550.0/            ! real setup M.L.Oct.95
c     data xleftmin/-426.5/    !3-18-2011 JAM
c     data xleftmax/-546.5/    !AVG now 486.5
c     !Note also this coordinate system appears
c     !backwards compared to survey C1305.

      data xleftmin/-425.19/
      data xleftmax/-545.19/

c     data xrightmin/418.0/   ! real setup M.L. Oct. 95
c     data xrightmax/562.0/   ! real setup M.L. Oct. 95
c     data xrightmin/410.7/    ! From survey 2010
c     data xrightmax/554.7/

      data xrightmin/413.49/    ! From survey 2010
      data xrightmax/557.49/




      data xhodomin /406.0/     !distance edge of hodo 1 to beam
      data dxhodo / 12.0/       !hodoscope segment width 12mm

      data iasy/0/

      real*8 yleft              !left detector offset 3-18-2011 JAM
      real*8 yright             !Right det. offset. Unsure of coord. system

c     yleft = -4.1    !JAM left detector offset from avg
c     yright = -2.2   !JAM right detector offset from avg

c     yleft = -0.74    !Survey Fall 2011
c     yright = -0.63   !Survey Fall 2011

C     test !!!!
C     I think the sign really is blown here..
      yleft = 0.74              !Survey Fall 2011
      yright = 0.63             !Survey Fall 2011


c     collimator/detector offsets
c     xleftoffset=-3.5 ! this is an absolute value relative to 490 mm
c     xrightoffset=-7.3

      xleftoffset=-4.81         ! this is an absolute value relative to 490 mm
      xrightoffset=-4.51

      xlefthodomin=xhodomin+xleftoffset
      xrighthodomin=xhodomin+xrightoffset

c     xlefthodomin=xhodomin
c     xrighthodomin=xhodomin



      do j=1,2

         left = .false.
         right = .false.
c     x = vec(1,j)    !JAM
c     y = vec(2,j)    !(These two lines are original code, commented out Mar 18, 2011)

         x = vec(1,j)

         if (x .GT. 0) then     !JAM this code simulates a
            y = vec(2,j) + yright !left/right detector offset
         else
            y = vec(2,j) + yleft
         endif


         if (x .le. 0.0) then
            ileft=  1+(dabs(x)-xlefthodomin)/dxhodo ! hodo paddle
            ileft= 15-ileft     !make it count like real hodo does
c     ileft= 30-ileft  !make it count like real hodo does
         else
            iright=  1+(dabs(x)-xrighthodomin)/dxhodo ! hodo paddle
         endif

c--   Make Detector Mask Horizontal Cuts...
         if ((x .ge. xleftmax) .and.
     >        (x .le. (xleftmin))) then
            left = .true.
         elseif ((x .ge. xrightmin) .and.
     >           (x .le. xrightmax)) then
            right = .true.
         endif
c
c--   ... and Detector Mask Vertical Cuts.
         if (left) then
            if (abs(y).gt.-1./12.*(x+xleftoffset)-15.83333333) then
               left = .false.
            endif
         endif
c
         if (right) then
            if (abs(y) .gt. 1./12.*(x-xrightoffset)-6.83333333) then
               right = .false.
            endif
         endif
c
         if (.not.(left .or. right)) then
            singlemol(j) = .false.
            twomols = .false.
         endif



      enddo
      CALL hfill(504,real(ileft),real(iright),1.0)
      if (twomols) then
         CALL hfill(505,real(ileft),real(iright),1.0)
         isum= (ileft + iright)
         idif= (ileft-iright)
         if(iabs(idif) .le. 2) call hfill(511,float(isum)/2.0,0.,1.0)
         if(iabs(idif) .le. 3) call hfill(512,float(isum)/2.0,0.,1.0)
c
         if(iabs(isum-15) .le. 2) call hfill(513,float(idif), 0.0, 1.0)
         if(iabs(isum-15) .le. 3) call hfill(514,float(idif), 0.0, 1.0)
         if(iabs(isum-15) .le. 7) call hfill(515,float(idif), 0.0, 1.0)
         if(iabs(isum-15) .le. 15) call hfill(516,float(idif), 0.0, 1.0)
      endif
c
C
C     SUM THE EVENT WEIGHTS FOR THE TWO POLARIZATION STATES
C     TIMES 3 ASYMMETRIES FOR EACH BLOCK PAIR
C     3 ASYMMETRIES: PZ(beam)*PZ(target)
C     PX(beam)*PX(target)
C     PX(beam)*PY(target) (should be 0 if accepted azimuth
C     is in the x or y planes)
c     do correct pad weighting
      do j=1,2
         if (singlemol(j)) then
            do k=1,3
               do i=1,2
                  WGT=PWT(I,k)*weights(j)/wgtcx
                  padtl(i,k,npads,j)=padtl(i,k,npads,j) + wgt
                  padtl2(i,k,npads,j)=padtl2(i,k,npads,j) + wgt**2
ccc
ccc   ntuple 9,10 -> Axx  11,12 -> Axy  13,14 -> Azz
ccc
c     ntuple(6+k*2+i)=WGT
c     ntuple(12+k*2+i)=WGT**2
ccc
               enddo
            enddo
         endif
      enddo
C
c     do correct block weighting... only 1 block pair for CEBAF polarimeter
      if (twomols) then         !skip only single moller events
         DO J=1,3
            DO I=1,2
               WGT=PWT(I,J)*WEIGHT/wgtcx
               BLKTL(I,J)=BLKTL(I,J)+WGT
               BLKTL2(I,J)=BLKTL2(I,J)+WGT**2
               LOCALWT(I,J)=WGT
            ENDDO
         ENDDO

         AZ_local = (LOCALWT(2,3)-LOCALWT(1,3))/
     >        (LOCALWT(1,3)+LOCALWT(2,3))
         avgA = avgA+AZ_local*(LOCALWT(1,3)+LOCALWT(2,3))
         avgA2 = avgA2+AZ_local**2*(LOCALWT(1,3)+LOCALWT(2,3))


      endif



 900  return

      end

c__________________________________________________________________________
c**********************************************************************
c***  ***
c***  Subroutine Quadrupole                                         ***
c***  ***
c**********************************************************************
c
      subroutine Quadrupole(vec, Leff, gradient, nstep)
c     -------------------------------------------------
c
c     Construct the first order matrix for a long, strong quadrupole by
c     dividing the quadrupole into nstep slices and applying the thin-lense
c     approximation for each of these slices. It then transfers a particle
c     with entrance coordinates x,y,z and momentum components px,py,pz
c     (contained in the six-dimensional vector vec) through the quadrupole.
c     Units used are:
c     x,y,z:     mm
c     px,py,pz:  GeV/c
c     Leff:      mm
c     gradient:  Tesla/mm
c
c     29-7-94
c     dB
      implicit none
c
      real*8  vec(6), Leff, gradient
      integer nstep
      real*8 rnstep
c
      real*8  unit(4,4), drft(4,4), quad(4,4), total(4,4), result(4,4)
      real*8  initial(4), final(4)
      real*8  p, Brho, c, me, pi
      real*8 rtmp, Leff_tmp
      integer i,j,k
c
      data c/2.99792458d8/
      data me/0.5110034/
      data pi/3.1415926535/
c
      rnstep = nstep
      if (nstep .le. 0) stop
     *     'Number of steps <= 0 in subroutine Quadrupole'
      p = sqrt(vec(4)**2 + vec(5)**2 + vec(6)**2)
      Brho = 1.0E12*p/c

      initial(1) = vec(1)
      initial(2) = atan(vec(4)/vec(6))
      initial(3) = vec(2)
      initial(4) = atan(vec(5)/vec(6))
c


      call Drift(unit, 0.0d0)
      call Drift(drft, 0.5d0*Leff/rnstep)


      call Thin_quad(quad, Leff*gradient/Brho/rnstep)

c
      call Matrix_product(result, quad, drft)
      do i=1,4
         do j=1,4
            quad(i,j)=result(i,j)
         enddo
      enddo

      call Matrix_product(result, drft, quad)
      do i=1,4
         do j=1,4
            quad(i,j)=result(i,j)
         enddo
      enddo
      call Matrix_product(total, unit, unit)
      do i=1,nstep
         call Matrix_product(result, quad, total)
         do j=1,4
            do k=1,4
               total(j,k)=result(j,k)
            enddo
         enddo
      enddo
      call Project(final, total, initial)

c
      vec(1) = final(1)
      vec(2) = final(3)
      vec(3) = vec(3) + Leff
      vec(4) = p*dsin(final(2))
      vec(5) = p*dsin(final(4))
      vec(6) = dsqrt(p**2 - vec(4)**2 - vec(5)**2)
c
      return
      end

c__________________________________________________________________________
c**********************************************************************
c***  ***
c***  Subroutine Quadrupole2                                         ***
c***  ***
c**********************************************************************
c
      subroutine Quadrupole2(vec, Leff, gradient, nstep)
c     -------------------------------------------------
c
c     Construct the first order matrix for a long, strong quadrupole by
c     dividing the quadrupole into nstep slices and applying the thin-lense
c     approximation for each of these slices. It then transfers a particle
c     with entrance coordinates x,y,z and momentum components px,py,pz
c     (contained in the six-dimensional vector vec) through the quadrupole.
c     Units used are:
c     x,y,z:     mm
c     px,py,pz:  GeV/c
c     Leff:      mm
c     gradient:  Tesla/mm
c
c     29-7-94
c     dB
      implicit none
c
      real*8  vec(6), Leff, gradient
      integer nstep
      real*8 rnstep
c
      real*8  unit(4,4), drft(4,4), quad(4,4), total(4,4), result(4,4)
      real*8  initial(4), final(4)
      real*8  p, Brho, c, me, pi
      real*8 rtmp, Leff_tmp
      integer i,j,k
c
      data c/2.99792458d8/
      data me/0.5110034/
      data pi/3.1415926535/
c
      rnstep = nstep
      if (nstep .le. 0) stop
     *     'Number of steps <= 0 in subroutine Quadrupole'
      p = sqrt(vec(4)**2 + vec(5)**2 + vec(6)**2)
      Brho = 1.0E12*p/c

      initial(1) = vec(1)
      initial(2) = atan(vec(4)/vec(6))
      initial(3) = vec(2)
      initial(4) = atan(vec(5)/vec(6))
c


      call Drift(unit, 0.0d0)
      call Drift(drft, 0.5d0*Leff/rnstep)


c     c      call Thin_quad(quad, Leff*gradient/Brho/rnstep)
c     c
ccc
c     c      call Matrix_product(result, quad, drft)
c     c do i=1,4
c     c   do j=1,4
c     c     quad(i,j)=result(i,j)
c     c   enddo
c     c enddo
c     c
c     c      call Matrix_product(result, drft, quad)
c     c do i=1,4
c     c   do j=1,4
c     c     quad(i,j)=result(i,j)
c     c   enddo
c     c enddo
c     c      call Matrix_product(total, unit, unit)
c     c      do i=1,nstep
c     c        call Matrix_product(result, quad, total)
c     c do j=1,4
c     c   do k=1,4
c     c     total(j,k)=result(j,k)
c     c   enddo
c     c enddo
c     c      enddo
c     c      call Project(final, total, initial)

C     Here I try to account for the reduction in effective length of Q2/Q3 when away from
C     the center

      do k=1,nstep

         rtmp=sqrt(initial(1)**2+initial(3)**2)
         rtmp=rtmp/25.4
         Leff_tmp = Leff*(1+ (0.33785E-3+0.26718*rtmp-0.50202*rtmp**2+
     >        0.29045*rtmp**3-0.057843*rtmp**4)/100.0 )
         call Thin_quad(quad, Leff_tmp*gradient/Brho/rnstep)

         call Matrix_product(result, quad, drft)
         do i=1,4
            do j=1,4
               quad(i,j)=result(i,j)
            enddo
         enddo

         call Matrix_product(result, drft, quad)
         do i=1,4
            do j=1,4
               quad(i,j)=result(i,j)
            enddo
         enddo
c     call Matrix_product(total, unit, unit)

         call Project(final, quad, initial)

         do i =1,4
            initial(i)=final(i)
         enddo
      enddo



c
      vec(1) = final(1)
      vec(2) = final(3)
      vec(3) = vec(3) + Leff
      vec(4) = p*dsin(final(2))
      vec(5) = p*dsin(final(4))
      vec(6) = dsqrt(p**2 - vec(4)**2 - vec(5)**2)
c
      return
      end





c**********************************************************************
c***  ***
c***  Subroutine Solenoid                                           ***
c***  ***
c**********************************************************************
c
      subroutine Solenoid(vec, Leff, Blong, nstep)
c     --------------------------------------------
c
c     Construct the first order matrix for a long, strong solenoid by
c     dividing the solenoid into nstep slices and applying the thin-lense
c     approximation for each of these slices. It then transfers a particle
c     with entrance coordinates x,y,z and momentum components px,py,pz
c     (contained in the six-dimensional vector vec) through the solenoid.
c     Units used are:
c     x,y,z:     mm
c     px,py,pz:  GeV/c
c     Leff:      mm
c     Blong:     Tesla
c
c     29-7-94
c     dB
      implicit none
c
      real*8  vec(6), Leff, Blong
      integer nstep
c
      real*8  unit(4,4), drft(4,4), slnd(4,4), total(4,4),result(4,4)
      real*8  initial(4), final(4)
      real*8  p, Brho, c, me, pi
      integer i,j,k
c
      data c/2.99792458E+08/
      data me/0.5110034/
      data pi/3.1415926535/
c
      if (nstep .le. 0) stop
     *     'Number of steps <= 0 in subroutine Solenoid'
      p = sqrt(vec(4)**2 + vec(5)**2 + vec(6)**2)
      Brho = 1.0E12*p/c
      initial(1) = vec(1)
      initial(2) = atan(vec(4)/vec(6))
      initial(3) = vec(2)
      initial(4) = atan(vec(5)/vec(6))
c
      call Drift(unit, 0.0d0)
      call Drift(drft, 0.5*Leff/nstep)
c     Note: calling with Leff*Blong/Brho/nstep > 0 calculates matrix for a
c     positively charge particle.
c     Need to call with strength < 0 for electron

c     call Thin_slnd(slnd, -Leff*Blong/Brho/nstep)
      call Thin_slnd2(slnd, -Leff*Blong/Brho/nstep/2., Leff)
c     call Thin_slnd(slnd, Leff*Blong/Brho/nstep)
c
      call Matrix_product(slnd, slnd, drft)
      call Matrix_product(result, drft, slnd)
      do i=1,4
         do j=1,4
            slnd(i,j)=result(i,j)
         enddo
      enddo
      call Matrix_product(total, unit, unit)
      do i=1,nstep
         call Matrix_product(result, slnd, total)
         do j=1,4
            do k=1,4
               total(j,k)=result(j,k)
            enddo
         enddo
      enddo
      call Project(final, total, initial)
c
      vec(1) = final(1)
      vec(2) = final(3)
      vec(3) = vec(3) + Leff
      vec(4) = p*sin(final(2))
      vec(5) = p*sin(final(4))
      vec(6) = sqrt(p**2 - vec(4)**2 - vec(5)**2)
c
      return
      end




c**********************************************************************
c***  ***
c***  Subroutine Drift                                              ***
c***  ***
c**********************************************************************
c
      subroutine Drift(drft, L)
c     --------------------------
c
c     Creates a 4x4 matrix for a drift space with length L
c
      implicit none
c
      real*8  drft(4,4), L
      integer i,j
c
      do i=1,4
         do j=1,4
            drft(i,j) = 0.
         enddo
         drft(i,i) = 1.
      enddo
      drft(1,2) = L
      drft(3,4) = L
c
      return
      end




c**********************************************************************
c***  ***
c***  Subroutine Thin_quad                                          ***
c***  ***
c**********************************************************************
c
c
      subroutine Thin_quad(quad, strength)
c     ------------------------------------
c
c     Creates a 4x4 matrix for a quadrupole with effective length Leff
c     and focusing strength strength = Leff * gradient / (B * rho).
c     strength > 0 means horizontally defocusing, vertically focusing
c     strength < 0 means horizontally focusing, vertically defocusing.
c
      implicit none
c
      real*8  quad(4,4), strength
      integer i,j
c
      do i=1,4
         do j=1,4
            quad(i,j) = 0.
            quad(i,i) = 1.
         enddo
      enddo
      quad(2,1) = strength
      quad(4,3) = -strength
c
      return
      end





c**********************************************************************
c***  ***
c***  Subroutine Thin_slnd                                          ***
c***  ***
c**********************************************************************
c
      subroutine Thin_slnd(slnd, strength)
c     ------------------------------------
c
c     Creates a 4x4 matrix for a solenoid with effective length Leff
c     and focusing strength strength = Leff * longitudinal field / (B * rho).
c
      implicit none
c
      real*8  slnd(4,4), strength
      integer i,j
c
      do i=1,4
         do j=1,4
            slnd(i,j) = 0.
            slnd(i,i) = 1.
         enddo
      enddo
      slnd(2,4) = strength
      slnd(4,2) = -strength
c
      return
      end

c**********************************************************************
c***  ***
c***  Subroutine Thin_slnd2                                          ***
c***  ***
c**********************************************************************
c
      subroutine Thin_slnd2(slnd, strength,Leff)
c     ------------------------------------
c
c     Creates a 4x4 matrix for a solenoid with effective length Leff
c     and focusing strength strength = Leff * longitudinal field / (B * rho).
c
      implicit none
c
      real*8  slnd(4,4), strength, Leff
      integer i,j
c
      do i=1,4
         do j=1,4
            slnd(i,j) = 0.
            slnd(i,i) = 1.
         enddo
      enddo
      slnd(2,4) = strength
      slnd(4,2) = -strength
c     New
      slnd(1,3)=strength
      slnd(3,1)=-strength
      slnd(2,1)=-strength**2/Leff
      slnd(4,3)=-strength**2/Leff


      return
      end




c**********************************************************************
c***  ***
c***  Subroutine Matrix_product                                     ***
c***  ***
c**********************************************************************
c
      subroutine Matrix_product(result, left, right)
c     ----------------------------------------------
c
c     Computes the 4x4 matrix product result = left * right
c
      implicit none
c
      real*8  result(4,4), left(4,4), right(4,4)
      real*8  tmpl(4,4), tmpr(4,4)
      integer i,j,k
c
      do i=1,4
         do j=1,4
            tmpl(i,j) = left(i,j)
            tmpr(i,j) = right(i,j)
            result(i,j) = 0.
            k=j                 !to make code work right
c     (what the hack is this for??? ML 9/96)
         enddo
      enddo
c
      do i=1,4
         do j=1,4
            do k=1,4
               result(i,j) = result(i,j) + tmpl(i,k) * tmpr(k,j)
            enddo
         enddo
      enddo
c
      return
      end




c**********************************************************************
c***  ***
c***  Subroutine Project                                            ***
c***  ***
c**********************************************************************
c
      subroutine Project(result, matrix, vector)
c     ------------------------------------------
c
c     Projects a 4x4 matrix onto a vector: result = matrix * vector
c
      implicit none
c
      real*8  matrix(4,4), result(4), vector(4)
      real*8  tmpv(4)
      integer i,j
c
      do i=1,4
         tmpv(i) = vector(i)
         result(i) = 0.
      enddo
c
      do i=1,4
         do j=1,4
            result(i) = result(i) + matrix(i,j) * tmpv(j)
         enddo
      enddo
c
      return
      end




c**********************************************************************
c***  ***
c***  Subroutine Detector                                           ***
c***  ***
c**********************************************************************
c
      subroutine Detector(vec, singlemol, twomols)
c     --------------------------------------------
c
c     Checks whether a Moller electron falls into one of the two trapezoids
c     composing the detector system. The acceptances of the left and right
c     detectors need not be identical.
c     5-8-94
c     dB
      implicit none
c
      real*8  vec(6)
      logical singlemol, twomols
c
      real*8  x, xleftmin, xleftmax, xrightmin, xrightmax
      real*8  y, yleftmin, yleftmax, yrightmin, yrightmax
      logical left, right
c
      data xleftmin/-420.0/
      data xleftmax/-560.0/
      data xrightmin/430.0/
      data xrightmax/550.0/
      data yleftmin/25.0/
      data yleftmax/35.0/
      data yrightmin/20.0/
      data yrightmax/30.0/
c
c     data xleftmin/-390.0/
c     data xleftmax/-590.0/
c     data xrightmin/390.0/
c     data xrightmax/590.0/
c     data yleftmin/40.0/
c     data yleftmax/60.0/
c     data yrightmin/40.0/
c     data yrightmax/60.0/
c
      left = .false.
      right = .false.
      x = vec(1)
      y = vec(2)
c
      if ((x .ge. xleftmax) .and. (x .le. xleftmin)) then
         left = .true.
      elseif ((x .ge. xrightmin) .and. (x .le. xrightmax)) then
         right = .true.
      endif
c
      if (left) then
         if (abs(y) .gt. yleftmin +
     *        (x-xleftmin)*(yleftmax-yleftmin)/(xleftmax-xleftmin))
     *        left = .false.
      endif
c
      if (right) then
         if (abs(y) .gt. yrightmin +
     *        (x-xrightmin)*(yrightmax-yrightmin)/(xrightmax-xrightmin))
     *        right = .false.
      endif
c
      if (.not.(left .or. right)) then
         singlemol = .false.
         twomols = .false.
      endif
c
      return
      end




c**********************************************************************
c***  ***
c***  Auxiliary subroutines                                         ***
c***  ***
c**********************************************************************
c
c
c
      real*8 function top(E,m)
      real*8 E, m
      top = sqrt(E*E-m*m)
      end
c
      real*8 function toE(p,m)
      real*8 p, m
      toE = sqrt(p*p+m*m)
      end
c
      real*8 function r2d(angle)
      real*8 angle, pi
      data pi/3.1415926535/
      r2d = 180.0*angle/pi
      end
c
      real*8 function d2r(angle)
      real*8 angle, pi
      data pi/3.1415926535/
      d2r = pi*angle/180.0
      end


      Function Rndm(idummy)
      DOUBLE PRECISION rvec, rvec2,pi
      data pi/3.14159265d0/
c     V116 Double Precision Random #.
      call rm48(rvec,1)
      call rm48(rvec2,1)

      rndm = sngl(rvec)
      return
      end


      subroutine rannor(a1,b1)

c     Normal distribution
      DOUBLE PRECISION rvec, rvec2,a,b,pi
      data pi/3.14159265d0/

      call rm48(rvec,1)
      call rm48(rvec2,1)
      a = dsin(2.d0*pi*rvec)*dsqrt(-2.d0*dlog(rvec2))

      call rm48(rvec,1)
      call rm48(rvec2,1)
      b = dsin(2.d0*pi*rvec)*dsqrt(-2.d0*dlog(rvec2))

      a1=sngl(a)
      b1=sngl(b)
      RETURN
      END

c-----------------------------------------------------------------------
      subroutine misalign(vec,dx,dy,dxp,dyp)
c--
c--   Add a displacement and slope to the trajectory vector in order to
c--   provide simulation of an offset or rotation of a component.
c--   Array 'vec(6)' is x,y,z,px,py,pz  in mm and GeV/c. It is what gets
c--   modified by this routing.
c--   The alterations are specified by
c--   dx,dy: displacements (mm, same units as vec)
c--   dxp,dyp: x and y slope offsets (dx/dz and dy/dz) which are
c--   used to modify the momentum components. They are assumed
c--   to be small.
c-----------------------------------------------------------------------
      real*8 vec(6),dx,dy,dxp,dyp
      real*8 pp2,xp,yp
c--
      vec(1)= vec(1) + dx
      vec(2)= vec(2) + dy
      pp2= ( vec(4)**2 + vec(5)**2 + vec(6)**2 )
      xp= vec(4)/vec(6)
      yp= vec(5)/vec(6)
      xp= xp + dxp
      yp= yp + dyp
      vec(6)= sqrt(pp2 / (xp**2 + yp**2 + 1.0) )
      vec(4)= xp*vec(6)
      vec(5)= yp*vec(6)
      return
      end
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
      subroutine quad_misalign(vec,dx,dy)
c--
c--   Add an x- and y-  displacement to the trajectory vector in order to
c--   simulate an offset of the quads.
c--   Since we cannot physically move the quads we must move the particle.
c--   Array 'vec(6)' is x,y,z,px,py,pz  in mm and GeV/c. Only vec() 1 and 2
c--   get modified in this routing.
c--   The alterations are specified by
c--   dx,dy: displacements (mm, same units as vec)
c--   Created by JMagee and BRislow feb, 2011.
c-----------------------------------------------------------------------
      real*8 vec(6),dx,dy
c--
      vec(1)= vec(1) + dx
      vec(2)= vec(2) + dy
      return
      end
c-----------------------------------------------------------------------

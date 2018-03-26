c      PROGRAM TEST
c      call testsnakedipole()
c     real bx,by,bz
c      call snakedipole(0.,-8400./sqrt(2.),8400./sqrt(2.),bx,by,bz,
c     1 0.3302653)
c      print*,bx,by,bz
c      END
c     Input is in mm and T, output is in T.
      subroutine snakedipole(x,y,z,bx,by,bz,b_0)
      integer in
      real dipx,dipy,dipz,fact1(5),fact2(5)
      real sx,sy,sz,x,y,z,bx,by,bz,b_0,sbx,sby,sbz
c      print*,'beginning snake dipole'
c      print*,'B_0',b_0
      sbx=0.
      sby=0.
      sbz=0.
      call g4mctosnake(x,y,z,sx,sy,sz)
c      print*,x,y,z
c      print*,sx,sy,sz
c     dipolert2 takes cm instead of mm because the function was stolen from raytrace
      dipx=sx/10.
      dipy=sy/10.
      dipz=sz/10.
      in=0
c     Determine which region you are asking for the field
c     1: before dipole entrance
c     2: in the dipole
c     3: after dipole exit
      if(z<-2000.)then
         in=0
      else if(sy>=0)then
         in=3
      else
         in=1
      endif
c      print*,in
      if(in==1.or.in==2)then
c     the coordinates y and z are flipped for dipolert2
c         print*,in
         call dipolert2(dipx, dipz, dipy, sbx, sbz, sby, b_0)
         sbx=-sbx
         sby=-sby
         sbz=-sbz
      else if(in==3)then
c         print*,in
c     You can take what SNAKE has, or set your own b_0
c         fact1(1)= 0.3302653
c         print*,b_0,"versus",0.3302653
         fact1(1)= b_0
         fact1(2)=-30. 
         fact2(1)=-1.157615215
         fact2(2)= 0.
         fact2(3)= 0.
         fact2(4)= 0.
         fact2(5)= 0.
         qrad1= 300.
         qrad2= 300.
         dl=1591.0
         call dqfringe2(sx,sy,sz,dl,sbx,sby,sbz,fact1,fact2,qrad1,qrad2)
      else if(in==0)then
c         print*,in
         sbx=0.
         sby=0.
         sbz=0.
      else
         print*,'Something has gone very wrong'
      endif
c      print*,'BSSSSS',bx,by,bz
c      if(abs(bx)>100..or.abs(by)>100..or.abs(bz)>100.)then
c         bx = 0.;
c         by = 0.;
c         bz = 0.;
c      endif
c      if ( in == 3 ) then
c         print*,sbx,sby,sbz
c      endif
      call snaketog4mc(bx,by,bz,sbx,sby,sbz)
c      print*,bx,by,bz
c      print*,'exiting snake dipole, returning to G4MC'
c      by = 0.
c      bz = 0.
      end

c     Convert from g4mc dipole coordinates to snake coordinates, both in mm
      subroutine g4mctosnake(x,y,z,sx,sy,sz)
      real xtr,ytr,ztr,xrt,yrt,zrt,r0,theta
      r0=8400.
c     Perform translation from dipole "center of cirle" to dipole exit face
      xtr=x
      ytr=y+r0*sin(45./57.29578)
      ztr=z-r0*cos(45./57.29578)
c     Perform rotation to make +z perp to dipole face, NOT direction of beam
c     Of course, those are the SNAKE coordinates
      theta=-15.
      xrt=xtr
      yrt=ytr*cos(theta/57.29578) - ztr*sin(theta/57.29578)
      zrt=ytr*sin(theta/57.29578) + ztr*cos(theta/57.29578)
      sx= zrt
      sy= yrt
      sz=-xrt
      end

c     Convert from SNAKE FIELD to G4MC FIELD
c     NOT FOR COORDINATES!!!!! FIELD ONLY
      subroutine snaketog4mc(bx,by,bz,bsx,bsy,bsz)
      real xtr,ytr,ztr,xrt,yrt,zrt,r0,theta
c      r0=8400.
c      print*,'input',sx,sy,sz
c      xrt=-sz
c      yrt=sy
c      zrt=sx
      xrt=-bsz
      yrt=-bsx
      zrt= bsy
c     print*,'change axes', xrt,yrt,zrt
      theta=15.
      bx=xrt
      bz=-yrt*cos(theta/57.29578) + zrt*sin(theta/57.29578)
      by= yrt*sin(theta/57.29578) + zrt*cos(theta/57.29578)
c      y=yrt*cos(theta/57.29578) - zrt*sin(theta/57.29578)
c      z=yrt*sin(theta/57.29578) + zrt*cos(theta/57.29578)
c      print*,'perform rotation',x,y,z
c     Perform translation from dipole "center of cirle" to dipole exit face
c      x=xtr
c      y=ytr-r0*sin(45./57.29578)
c      z=ztr+r0*cos(45./57.29578)
c     Perform rotation to make +z perp to dipole face, NOT direction of beam
c     Of course, those are the SNAKE coordinates

      end

c  dipolert2- analytic calculation of dipole magnetic fields
c            stolen from raytrace adapted for use in snake
c                                         -jjl  8/21/89
c    dipolert2 is called from snake in analyf
c                snake provides x,y,z in the c-axis system in cm
c                dipolert2 returns bx,by,bz in the same system
c
      subroutine dipolert2(x,y,z,bbx,bby,bbz,b_0)
      implicit real*8(a-h,o-z)
      real*8 dcosd,dsind
      external dcosd,dsind
      real x,y,z,bbx,bby,bbz,b_0
      real*8  ndx, k
      common  /blck10/  bx, by, bz, k, tc, dtc
      common  /blck20/  ndx,bet1,gama,delt,csc
      common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
      common  /blck22/  d, dg, s, bf, bt
      common  /blck23/  c0, c1, c2, c3, c4, c5
      common  /blck24/  rb, xc, zc
      common  /blck25/  in, mtyp
      dimension tc(2,6), dtc(6), bb(2,0:12)
c   2.     ,   4.    ,   2.    ,   8.    , 4.               --> lf1, lu1, lf2, dg, mtyp
c   0.     ,   0.    ,  25.    , 840.0   , 0.3302653        --> a  , b  , d  , rb, bf
c  45.     , -30.000 , -30.000                              --> phi, alpha, beta
c  -1.263                                                   --> ndx, bet1, gama, delt
c 130.     ,-100.    ,-100.    , 130.                       --> z11, z12, z21, z22
c   0.04725,   2.2395,   -.9768,    .7288,-.1299    ,.0222  --> c0, c1, c2, c3, c4, c5,
c   0.04725,   2.2395,   -.9768,    .7288,-.1299    ,.0222  --> c0, c1, c2, c3, c4, c5,
c   0.                                                      --> br1, br2, xcr1, xcr2, dels, dels?
c   0.                                                      --> rca, rca, scor, scor
c   2.680, -14.10, 50.0, -100.                              --> s2, s3, s4, s5, s6, s7, s8
c  -1.50, 20.26, -200.0                                     -->
      lf1  = 2.
      lu1  = 4.
      lf2  = 2.
      dg   = 8.
      mtyp = 4.
      a    = 0.
      b    = 0.
      d    = 25.
      rb   = 840.
c      bf   = 0.3302653
      bf = b_0
      phi  = 45.
      alpha= -30.
      beta = -30.
      ndx  = -1.263
c      you can turn off the index for debugging efforts
c      ndx  = 0.
      bet1 = 0.
      gama = 0.
      delt = 0.
      z11  =  130.
      z12  = -100.
      z21  = -100.
      z22  =  130.
      br1  = 0.
      br2  = 0.
      xcr1 = 0.
      xcr2 = 0.

c      write(8,*)' in dipolert2 1',x,y,z
c**** transform from c-axis to b-axis
c**** going from dex to den coordinates
      xpric=x+(2.*rb*dsind(phi/2.)*dsind((phi/2.)-beta))
      zpric=z+(2.*rb*dsind(phi/2.)*dcosd((phi/2.)-beta))
c bug check -jjl 12/14/95
c      xpric=x+(2.*rb*dsind(phi/2.)*dsind((phi/2.)-alpha))
c      zpric=z+(2.*rb*dsind(phi/2.)*dcosd((phi/2.)-alpha))
c end bug check
      cosa=dcosd(phi-alpha-beta)
      sina=dsind(phi-alpha-beta)
      tc(1,1)=(-xpric*cosa)+(zpric*sina)
      tc(1,2)=dble(y)
      tc(1,3)=(-xpric*sina)-(zpric*cosa)
c**** transform to second vfb coord system
c***
      copab =dcos( (phi-alpha-beta)/57.29578)
      sipab =dsin( (phi-alpha-beta)/57.29578)
      cospb =dcos( (phi/2.-beta)/57.29578 )
      sinpb =dsin( (phi/2.-beta)/57.29578 )
      sip2 =dsin( (phi/2.)/57.29578 )
      xt = tc(1,1)
      zt = tc(1,3)
      vxt = tc(1,4)
      vzt = tc(1,6)
      tc(2,3) = - zt  *copab +  xt  *sipab -2.*rb*sip2*cospb
      tc(2,1) = - zt  *sipab -  xt  *copab -2.*rb*sip2*sinpb
      tc(2,6) = - vzt *copab +  vxt *sipab
      tc(2,4) = - vzt *sipab -  vxt *copab
      tc(2,2)=dble(y)
c      write(8,*)' in dipolert2 1',(tc(i),i=1,3)
c      print *,' in dipolert2 1', (tc(1, i), i=1,3)
c      print *,' in dipolert2 2', (tc(2, i), i=1,3)
c
c  test for region
c
      if(tc(1,3).gt.z11)then
       bx=0.
       by=0.
       bz=0.
       go to 10
      endif
      if((tc(1,3).ge.z12).and.(tc(2,3).lt.z21)) go to 1  !entrance fringe
      if((tc(1,3).ge.z12).and.(tc(2,3).ge.z21)) go to 25 !overlapping entr&exit
      if(tc(1,3).lt.z12) go to 2                         !uniform field or exit
c      go to 15
c****
c**** in designates magnet regions for bfun
c****
  1   in = 1
c      print *,'entrance fringe'
      xc= rb*dcos( alpha/ 57.29578 )
      zc=-rb*dsin( alpha/ 57.29578 )
c****
c   2.     ,   4.    ,   2.    ,   8.    , 4.               --> lf1, lu1, lf2, dg, mtyp
c   0.     ,   0.    ,  25.    , 840.0   , 0.3302653        --> a  , b  , d  , rb, bf
c  45.     , -30.000 , -30.000                              --> phi, alpha, beta
c  -1.263                                                   --> ndx, bet1, gama, delt
c 130.     ,-100.    ,-100.    , 130.                       --> z11, z12, z21, z22
c   0.04725,   2.2395,   -.9768,    .7288,-.1299    ,.0222  --> c0, c1, c2, c3, c4, c5,
c   0.04725,   2.2395,   -.9768,    .7288,-.1299    ,.0222  --> c0, c1, c2, c3, c4, c5,
c   0.                                                      --> br1, br2, xcr1, xcr2, dels, dels?
c   0.                                                      --> rca, rca, scor, scor
c   2.680, -14.10, 50.0, -100.                              --> s2, s3, s4, s5, s6, s7, s8
c  -1.50, 20.26, -200.0                                     -->
      c0   = 0.04725
      c1   = 2.2395
      c2   = -.9768
      c3   =  .7288
      c4   = -.1299
      c5   =  .0222
      dels = 0.
      rca  = 0.
      csc = dcos( alpha/57.29578 )
      scor = 0.
      s2   =   2.680 / rb    + rca/2.d0
      s3   = -14.10  / rb**2
      s4   =  50.0   / rb**3 + rca**3/8.d0
      s5   =-100.    / rb**4
      s6   =   0.    / rb**5 + rca**5/16.d0
      s7   =   0.    / rb**6
      s8   =   0.    / rb**7 + rca**7/25.6d0
c      print*,in,s2,s3,s4,s5,s6,s7,s8
c
      call ndip(1)
      write(8,*)' return from ndip',bx,by,bz
c**** transform to c-axis system
c***
      bxdum=bx
      bzdum=bz
      bx=(-bxdum*cosa)-(bzdum*sina)
      bz=(+bxdum*sina)-(bzdum*cosa)
      go to 10
  2   continue
      if(tc(2,3).lt.z21)go to 3      !uniform field
      if(tc(2,3).gt.z22)then         !outside of fringe region
       bx=0.
       by=0.
       bz=0.
       go to 10
      endif
      if(tc(2,3).le.z22)go to 4      !exit fringe
c****
c****
c**** uniform field integration region
c****
c****
  3   in = 2
c      print *,'uniform field'
ctest      xc=-rb*dcos( beta / 57.29578 )
      xc= rb*dcos( beta / 57.29578 )
      zc=-rb*dsin( beta / 57.29578 )
c
      call ndip(1)
c      btot=btot+(by*rho/ 57.29587)
c      if(theta.eq.22.5d+00)br0=by
c      write(1,*)' ',theta,by
      go to 10
c***
c***
c**** setup for second fringe field and integration
c****
c****
c   2.     ,   4.    ,   2.    ,   8.    , 4.               --> lf1, lu1, lf2, dg, mtyp
c   0.     ,   0.    ,  25.    , 840.0   , 0.3302653        --> a  , b  , d  , rb, bf
c  45.     , -30.000 , -30.000                              --> phi, alpha, beta
c  -1.263                                                   --> ndx, bet1, gama, delt
c 130.     ,-100.    ,-100.    , 130.                       --> z11, z12, z21, z22
c   0.04725,   2.2395,   -.9768,    .7288,-.1299    ,.0222  --> c0, c1, c2, c3, c4, c5,
c   0.04725,   2.2395,   -.9768,    .7288,-.1299    ,.0222  --> c0, c1, c2, c3, c4, c5,
c   0.                                                      --> br1, br2, xcr1, xcr2, dels, dels?
c   0.                                                      --> rca, rca, scor, scor
c   2.680, -14.10, 50.0, -100.                              --> s2, s3, s4, s5, s6, s7, s8
c  -1.50, 20.26, -200.0                                     -->

  4   br   = br2
c   print *,'exit fringe'
      c0   = 0.04725
      c1   = 2.2395
      c2   = -.9768
      c3   = .7288
      c4   = -.1299
      c5   =  .0222
      dels = 0.
      rca  = 0.
      scor = 0.
      csc = dcos( beta /57.29578 )
      s2   =  -1.50 / rb    + rca/2.d0
      s3   =  20.26 / rb**2
      s4   =-200.   / rb**3 + rca**3/8.d0
      s5   =   0.   / rb**4
      s6   =   0.   / rb**5 + rca**5/16.d0
      s7   =   0.   / rb**6
      s8   =   0.   / rb**7 + rca**7/25.6d0
      in = 3
c      print*,in,s2,s3,s4,s5,s6,s7,s8
      xc=-rb*dcos( beta / 57.29578 )
      zc=-rb*dsin( beta / 57.29578 )
c
      call ndip(2)
c      btot=btot+(by*rho/ 57.29587)
c      write(1,*)'  ',theta,by
  10  continue
      bbx=bx
      bby=by
      bbz=bz
c      write(8,*)' returning ',bx,by,bz,bbx,bby,bbz
      return
c
c overlapping entrance and exit fringe fields
c
  25  in = 1                         !setup for entrance field
c      print *,'overlaping entrance/exit'
      xc= rb*dcos( alpha/ 57.29578 )
      zc=-rb*dsin( alpha/ 57.29578 )
c****
c   2.     ,   4.    ,   2.    ,   8.    , 4.               --> lf1, lu1, lf2, dg, mtyp
c   0.     ,   0.    ,  25.    , 840.0   , 0.3302653        --> a  , b  , d  , rb, bf
c  45.     , -30.000 , -30.000                              --> phi, alpha, beta
c  -1.263                                                   --> ndx, bet1, gama, delt
c 130.     ,-100.    ,-100.    , 130.                       --> z11, z12, z21, z22
c   0.04725,   2.2395,   -.9768,    .7288,-.1299    ,.0222  --> c0, c1, c2, c3, c4, c5,
c   0.04725,   2.2395,   -.9768,    .7288,-.1299    ,.0222  --> c0, c1, c2, c3, c4, c5,
c   0.                                                      --> br1, br2, xcr1, xcr2, dels, dels?
c   0.                                                      --> rca, rca, scor, scor
c   2.680, -14.10, 50.0, -100.                              --> s2, s3, s4, s5, s6, s7, s8
c  -1.50, 20.26, -200.0                                     -->

      c0   = 0.04725
      c1   = 2.2395
      c2   = -.9768
      c3   = .7288
      c4   = -.1299
      c5   =  .0222
      dels = 0.
      rca  = 0.
      csc = dcos( alpha/57.29578 )
      scor = 0.
      s2   =   2.680/ rb    + rca/2.d0
      s3   = -14.10 / rb**2
      s4   =  50.   / rb**3 + rca**3/8.d0
      s5   =-100.   / rb**4
      s6   =   0.   / rb**5 + rca**5/16.d0
      s7   =   0.   / rb**6
      s8   =   0.   / rb**7 + rca**7/25.6d0
c      print*,in,s2,s3,s4,s5,s6,s7,s8
c
      call nndip(1,bb)
c***
c***
c**** setup for exit fringe field
c****
c****
      br   = br2
      c0   = 0.04725
      c1   = 2.2395
      c2   = -.9768
      c3   = .7288
      c4   = -.1299
      c5   =  .0222
      dels = 0.
      rca  = 0.
      scor = 0.
      csc = dcos( beta /57.29578 )
      s2   =  -1.50 / rb    + rca/2.d0
      s3   =  20.26 / rb**2
      s4   =-200.   / rb**3 + rca**3/8.d0
      s5   =   0.   / rb**4
      s6   =   0.   / rb**5 + rca**5/16.d0
      s7   =   0.   / rb**6
      s8   =   0.   / rb**7 + rca**7/25.6d0
      in = 3
c      print*,in,s2,s3,s4,s5,s6,s7,s8
      xc=-rb*dcos( beta / 57.29578 )
      zc=-rb*dsin( beta / 57.29578 )
      call nndip(2,bb)
      if(y.eq.0.)then
       bx=0.
       bz=0.
       by=bb(1,0)+bb(2,0)-bf
      else
      do 26 i=0,12
       bb(1,i)=bb(1,i)+bb(2,i)-bf
  26  continue
      yg1 = y/dg
      yg2 = yg1**2
      yg3 = yg1**3
      yg4 = yg1**4
c        1         2         3         4         5         6         7
      bx = yg1 * ( (bb(1,5)-bb(1,7))*2./3. - (bb(1,6)-bb(1,8))/12. ) +
     1     yg3*( (bb(1,5)-bb(1,7))/6. - (bb(1,6)-bb(1,8))/12. -
     2     (bb(1,3) + bb(1,11) - bb(1,4) - bb(1,12)
     3      - 2.*bb(1,5) + 2.*bb(1,7) ) / 12. )
      by = bb(1,0) - yg2*( ( bb(1,1) + bb(1,9) + bb(1,5) + bb(1,7)
     1     - 4.*bb(1,0) ) *2./3. -
     2     ( bb(1,2) + bb(1,10) + bb(1,6) + bb(1,8)-4.*bb(1,0))/24.) +
     3     yg4*(-(bb(1,1)+ bb(1,9)+ bb(1,5)+ bb(1,7)- 4.*bb(1,0) )/6.+
     4     ( bb(1,2)+ bb(1,10)+ bb(1,6) + bb(1,8) - 4.*bb(1,0))/24. +
     5     ( bb(1,3)+bb(1,11)+bb(1,4) +bb(1,12)-2.*bb(1,1)-2.*bb(1,9)-
     6     2.*bb(1,5) - 2.*bb(1,7) + 4.*bb(1,0) ) / 12. )
      bz = yg1*((bb(1,1)-bb(1,9))*2./3. - (bb(1,2)-bb(1,10) ) /12. ) +
     1     yg3*( (bb(1,1)- bb(1,9))/6. - (bb(1,2) - bb(1,10) ) / 12. -
     2     ( bb(1,3) + bb(1,4) - bb(1,11) - bb(1,12) 
     3      - 2.*bb(1,1) + 2.*bb(1,9) ) / 12.  )
      bt  =dsqrt(bx*bx + by*by + bz*bz)
      endif
      go to 10
      end
c
      subroutine ndip(l)
c****
c****
c**** mtyp = 3 or 4
c**** this version of bfun is mainly for nonuniform field magnets
c**** the central field region is represented to 3'rd order on-and-
c**** off the midplane by analytic expressions. see slac no. 75
c**** fringe field regions represented by fermi type fall-off
c**** along with radial fall-off
c**** components of 'b' in fringe region evaluated by numerical methods
c****
c****
c**** the relationship between b0, ......... b12 and b(i,j) relative to
c**** axes (z,x) is given by
c****
c****
c**** b0  = b( 0, 0 )
c**** b1  = b( 1, 0 )
c**** b2  = b( 2, 0 )
c**** b3  = b( 1, 1 )
c**** b4  = b( 1,-1 )
c**** b5  = b( 0, 1 )
c**** b6  = b( 0, 2 )
c**** b7  = b( 0,-1 )
c**** b8  = b( 0,-2 )
c**** b9  = b(-1, 0 )
c**** b10 = b(-2, 0 )
c**** b11 = b(-1, 1 )
c**** b12 = b(-1,-1 )
c****
c****
      implicit real*8(a-h,o-z)
      real*8  ndx, k
      common  /blck10/  bx, by, bz, k, tc, dtc
      common  /blck20/  ndx,bet1,gama,delt,csc
      common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
      common  /blck22/  d, dg, s, bf, bt
      common  /blck23/  c0, c1, c2, c3, c4, c5
      common  /blck24/  rb, xc, zc
      common  /blck25/  in, mtyp
      dimension tc(2,6), dtc(6)
c      print*, 'in ndip'
      x = tc(l,1)
      y = tc(l,2)
      z = tc(l,3)
      dx = x - xc
      dz = z - zc
      rp =dsqrt( dx**2 + dz**2 )
      dr = rp - rb
      go to ( 1, 2, 3, 14 ), in
    7 print 8, in, mtyp
      call exit(0)
    8 format (    '0 error -go to -  in bfun   in=', i3, '   mtyp=',i4 )
    2 drr1 = dr/rb
      drr2 = drr1*drr1
      drr3 = drr2*drr1
      drr4 = drr3*drr1
c      print *, 'index was 2'
      if( y .ne. 0. )  go to 4
c****
c**** mid-plane uniform field region
c****
c      print*,'midplane is simple', dr, rb, drr1
      bx = 0.
      by = 0.
      if( mtyp .eq. 3) by=
     1     bf* ( 1. - ndx*drr1 + bet1*drr2 + gama*drr3 + delt*drr4 )
      if( mtyp .eq. 4) by= bf/ (1. + ndx*drr1 )
c      print *, by, in, mytp
      bz = 0.
      bt = by
      return
c****
c**** non mid-plane uniform field region
c****
    4 yr1 = y/rb
      yr2 = yr1*yr1
      yr3 = yr2*yr1
      yr4 = yr3*yr1
      rr1 = rb/rp
      rr2 = rr1*rr1
      rr3 = rr2*rr1
      if( mtyp .eq. 3 ) go to 11
      if( mtyp .eq. 4 ) go to 12
      go to 7
c****
c**** mtyp = 3
c****
   11 brr = bf*( ( -ndx + 2.*bet1*drr1 + 3.*gama*drr2 + 4.*delt*drr3 )
     1   *yr1 - (ndx*rr2 + 2.*bet1*rr1*(1.-rr1*drr1) +
     2   3.*gama*( 2. + 2.*rr1*drr1 - rr2*drr2 ) +
     3   4.*delt*( 6.*drr1 + 3.*rr1*drr2 - rr2*drr3 ))*yr3/6. )
      by = bf* ( 1. - ndx*drr1 + bet1*drr2 + gama*drr3 + delt*drr4 -
     1   .5*yr2*( -ndx*rr1 + 2.*bet1*( 1. + rr1*drr1) +
     2   3.*gama*drr1*( 2. + rr1*drr1) + 4.*delt*drr2*(3. + rr1*drr1) )
     3   + yr4*( -ndx*rr3 + 2.*bet1*( rr3*drr1 - rr2) +
     4   3.*gama*( 4.*rr1 - 2.*rr2*drr1 + rr3*drr2 ) +
     5   4.*delt*( 6. + 12.*rr1*drr1 - 3.*rr2*drr2 + rr3*drr3 ) )/24. )
      go to 13
c****
c**** mtyp = 4
c****
   12 dnr1 = 1. + ndx*drr1
      dnr2 = dnr1*dnr1
      dnr3 = dnr2*dnr1
      dnr4 = dnr3*dnr1
      dnr5 = dnr4*dnr1
      brr = bf*ndx*( -yr1/dnr2 + yr3*( 6.*ndx*ndx/dnr4 -
     1   2.*ndx*rr1/dnr3 - rr2/dnr2 ) /6.  )
      by = bf*( 1./dnr1 + .5*yr2*ndx*( -2.*ndx/dnr3 + rr1/dnr2) +
     2   yr4*ndx*( 24.*ndx**3 /dnr5 - 12.*ndx*ndx*rr1/dnr4 -
     3   2.*ndx*rr2/dnr3 - rr3/dnr2 ) /24.  )
c****
c****
   13 bx = brr*dx/rp
      bz = brr*dz/rp
      bt  =dsqrt(bx*bx + by*by + bz*bz)
      return
c****
c****
    1 sine = -1.
      go to 5
    3 sine = 1.
    5 if( z  .gt. 0. ) dr = x * sine*csc
c      print*,x, y, z
      call ndpp( b0, z, x, y, dr      )
      if( y  .ne. 0. )  go to 6
c****
c**** mid-plane fringing field region
c****
      bx = 0.
      by = b0
      bz = 0.
      bt   = b0
      return
c****
c**** non mid-plane fringing field region
c****
    6 if( z .gt. 0. )  go to 9
      dr1  =       (dsqrt( dx**2 + (dz+dg)**2 ) - rb )
      dr2  =       (dsqrt( dx**2 + (dz+2.*dg)**2 ) - rb )
      dr3  =       (dsqrt( (dx+dg)**2 + (dz+dg)**2 )  - rb )
      dr4  =       (dsqrt( (dx-dg)**2 + (dz+dg)**2 )  - rb )
      dr5  =       (dsqrt( (dx+dg)**2 + dz**2 ) - rb )
      dr6  =       (dsqrt( (dx+ 2.*dg)**2 + dz**2 ) - rb )
      dr7  =       (dsqrt( (dx-dg)**2 + dz**2 ) - rb )
      dr8  =       (dsqrt( (dx- 2.*dg)**2 + dz**2 ) - rb )
      dr9  =       (dsqrt( dx**2 + (dz-dg)**2 ) - rb )
      dr10 =       (dsqrt( dx**2 + (dz-2.*dg)**2 ) - rb )
      dr11 =       (dsqrt( (dx+dg)**2 + (dz-dg)**2 )  - rb )
      dr12 =       (dsqrt( (dx-dg)**2 + (dz-dg)**2 )  - rb )
      go to 10
    9 dr1  = sine* x*csc
      dr2  = dr1
      dr9  = dr1
      dr10 = dr1
      dr3  = sine* ( x + dg )*csc
      dr5  = dr3
      dr11 = dr3
      dr4  = sine*( x - dg )*csc
      dr7  = dr4
      dr12 = dr4
      dr6  = sine* ( x + 2.*dg )*csc
      dr8  = sine* ( x - 2.*dg )*csc
c****
c****
   10 call ndpp ( b1 , z + dg, x , y , dr1 )
      call ndpp ( b2 , z + 2.*dg, x , y , dr2 )
      call ndpp ( b3 , z + dg, x + dg , y , dr3 )
      call ndpp ( b4 , z + dg, x - dg , y , dr4 )
      call ndpp ( b5 , z , x + dg , y, dr5 )
      call ndpp ( b6 , z , x + 2.*dg , y , dr6 )
      call ndpp ( b7 , z , x - dg , y, dr7 )
      call ndpp ( b8 , z , x - 2.*dg , y , dr8 )
      call ndpp ( b9 , z - dg, x , y , dr9 )
      call ndpp ( b10, z - 2.*dg, x, y, dr10 )
      call ndpp ( b11, z - dg, x + dg , y , dr11 )
      call ndpp ( b12, z - dg, x - dg , y , dr12 )
      yg1 = y/dg
      yg2 = yg1**2
      yg3 = yg1**3
      yg4 = yg1**4
c
c fudge to fix dipole fringe field -4/8/02  jjl
c
c      print*,'we are fudgin'
c      print*, y,dr1
      fudge=1.30986      
      bx = fudge*yg1 * ( (b5-b7)*2./3. - (b6-b8)/12. )  +
     1     yg3*( (b5-b7)/6. - (b6-b8)/12. -
     2     (b3 + b11 - b4 - b12 - 2.*b5 + 2.*b7 ) / 12. )
      by = b0 - yg2*( ( b1 + b9 + b5 + b7 - 4.*b0 ) *2./3. -
     1     ( b2 + b10 + b6 + b8 - 4.*b0 ) / 24. ) +
     2     yg4* (-( b1 + b9 + b5 + b7 - 4.*b0 ) / 6. +
     3     ( b2 + b10 + b6 + b8 - 4.*b0 ) / 24. +
     4     ( b3 + b11 + b4 + b12 - 2.*b1 - 2.*b9 -
     5     2.*b5 - 2.*b7 + 4.*b0 ) / 12. )
      bz = fudge*yg1*( (b1 - b9 ) *2./3. - ( b2 - b10 ) /12. ) +
     1     yg3*( ( b1 - b9 ) / 6. - ( b2 - b10 ) / 12. -
     2     ( b3 + b4 - b11 - b12 - 2.*b1 + 2.*b9 ) / 12.  )
      bt  =dsqrt(bx*bx + by*by + bz*bz)
      return
   14 bx = 0.
      by = br
      bz = 0.
      bt = br
      return
      end
c
      subroutine nndip(l,bb)
c****
c****
c**** mtyp = 3 or 4
c**** this version of bfun is mainly for nonuniform field magnets
c**** the central field region is represented to 3'rd order on-and-
c**** off the midplane by analytic expressions. see slac no. 75
c**** fringe field regions represented by fermi type fall-off
c**** along with radial fall-off
c**** components of 'b' in fringe region evaluated by numerical methods
c****
c****
c**** the relationship between b0, ......... b12 and b(i,j) relative to
c**** axes (z,x) is given by
c****
c****
c**** bb(l,0)  = b( 0, 0 )
c**** bb(l,1)  = b( 1, 0 )
c**** bb(l,2)  = b( 2, 0 )
c**** bb(l,3)  = b( 1, 1 )
c**** bb(l,4)  = b( 1,-1 )
c**** bb(l,5)  = b( 0, 1 )
c**** bb(l,6)  = b( 0, 2 )
c**** bb(l,7)  = b( 0,-1 )
c**** bb(l,8)  = b( 0,-2 )
c**** bb(l,9)  = b(-1, 0 )
c**** bb(l,10) = b(-2, 0 )
c**** bb(l,11) = b(-1, 1 )
c**** bb(l,12) = b(-1,-1 )
c****
c****
c modified for use with snake in order to deal with everlapping entrance and
c exit fringing fields - jjl 5/21/91
      implicit real*8(a-h,o-z)
      real*8  ndx, k
      common  /blck10/  bx, by, bz, k, tc, dtc
      common  /blck20/  ndx,bet1,gama,delt,csc
      common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
      common  /blck22/  d, dg, s, bf, bt
      common  /blck23/  c0, c1, c2, c3, c4, c5
      common  /blck24/  rb, xc, zc
      common  /blck25/  in, mtyp
      dimension tc(2,6), dtc(6),bb(2,0:12)
      x = tc(l,1)
      y = tc(l,2)
      z = tc(l,3)
      dx = x - xc
      dz = z - zc
      rp =dsqrt( dx**2 + dz**2 )
      dr = rp - rb
      go to ( 1, 7, 3, 7 ), in
    7 print 8, in, mtyp
      call exit(0)
    8 format ('0 error -go to -  in nndip   in=', i3, '   mtyp=',i4)
c****
c****
    1 sine = -1.
      go to 5
    3 sine = 1.
    5 if( z  .gt. 0. ) dr = x * sine*csc
      call ndpp( bb(l,0), z, x, y, dr      )
      if( y  .ne. 0. )  go to 6
c****
c**** mid-plane fringing field region
c****
      return
c****
c**** non mid-plane fringing field region
c****
    6 if( z .gt. 0. )  go to 9
      dr1  =       (dsqrt( dx**2 + (dz+dg)**2 ) - rb )
      dr2  =       (dsqrt( dx**2 + (dz+2.*dg)**2 ) - rb )
      dr3  =       (dsqrt( (dx+dg)**2 + (dz+dg)**2 )  - rb )
      dr4  =       (dsqrt( (dx-dg)**2 + (dz+dg)**2 )  - rb )
      dr5  =       (dsqrt( (dx+dg)**2 + dz**2 ) - rb )
      dr6  =       (dsqrt( (dx+ 2.*dg)**2 + dz**2 ) - rb )
      dr7  =       (dsqrt( (dx-dg)**2 + dz**2 ) - rb )
      dr8  =       (dsqrt( (dx- 2.*dg)**2 + dz**2 ) - rb )
      dr9  =       (dsqrt( dx**2 + (dz-dg)**2 ) - rb )
      dr10 =       (dsqrt( dx**2 + (dz-2.*dg)**2 ) - rb )
      dr11 =       (dsqrt( (dx+dg)**2 + (dz-dg)**2 )  - rb )
      dr12 =       (dsqrt( (dx-dg)**2 + (dz-dg)**2 )  - rb )
      go to 10
    9 dr1  = sine* x*csc
      dr2  = dr1
      dr9  = dr1
      dr10 = dr1
      dr3  = sine* ( x + dg )*csc
      dr5  = dr3
      dr11 = dr3
      dr4  = sine*( x - dg )*csc
      dr7  = dr4
      dr12 = dr4
      dr6  = sine* ( x + 2.*dg )*csc
      dr8  = sine* ( x - 2.*dg )*csc
c****
c****
   10 call ndpp ( bb(l,1) , z + dg, x , y , dr1 )
      call ndpp ( bb(l,2) , z + 2.*dg, x , y , dr2 )
      call ndpp ( bb(l,3) , z + dg, x + dg , y , dr3 )
      call ndpp ( bb(l,4) , z + dg, x - dg , y , dr4 )
      call ndpp ( bb(l,5) , z , x + dg , y, dr5 )
      call ndpp ( bb(l,6) , z , x + 2.*dg , y , dr6 )
      call ndpp ( bb(l,7) , z , x - dg , y, dr7 )
      call ndpp ( bb(l,8) , z , x - 2.*dg , y , dr8 )
      call ndpp ( bb(l,9) , z - dg, x , y , dr9 )
      call ndpp ( bb(l,10), z - 2.*dg, x, y, dr10 )
      call ndpp ( bb(l,11), z - dg, x + dg , y , dr11 )
      call ndpp ( bb(l,12), z - dg, x - dg , y , dr12 )
      return
      end
c
      subroutine  ndpp ( bfld, z, x, y , dr )
c****
c****
c****
c****
      implicit real*8(a-h,o-z)
      real*8  ndx, k
      common  /blck10/  bx, by, bz, k, tc, dtc
      common  /blck20/  ndx,bet1,gama,delt,csc
      common  /blck21/  rca,dels,br,s2,s3,s4,s5,s6,s7,s8,scor
      common  /blck22/  d, dg, s, bf, bt
      common  /blck23/  c0, c1, c2, c3, c4, c5
      common  /blck24/  rb, xc, zc
      common  /blck25/  in, mtyp
      dimension tc(2,6), dtc(6)
      drr1 = dr/rb
      drr2 = drr1*drr1
      drr3 = drr2*drr1
      drr4 = drr3*drr1
c****
c**** mtyp    :                  modified iterative procedure
c****
c      print*,'beginning modified iterative precedure'
c      print*, x,y,z,dr,bfld
      xp = x
      xp2 = xp*xp
      xp3 = xp2*xp
      xp4 = xp3 * xp
      zp = -(s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
     1       s7*xp4*xp3 + s8*xp4*xp4 )
      az = (z-zp)/10.d0
      azmax = dsqrt(  x*x + z*z  )
      if( az  .gt.  azmax  )  az = azmax
      zsign = z-zp
      rinv4 = 0.
      do 11 i=1,21
      xp   = x + az*(i-11)
      xp2 = xp*xp
      xp3 = xp2*xp
      xp4 = xp3*xp
      zp = -(s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
     1       s7*xp4*xp3 + s8*xp4*xp4 )
      xxp = x-xp
      zzp = z-zp
      dd =            xxp*xxp + zzp*zzp
      if( dd  .lt.  1.d-15 )  dd = 1.d-15
      if( dd  .gt.  1.d15  )  dd = 1.d15
      rinv4 = rinv4 + 1.0d0 / (dd*dd )
   11 continue
      dp = dsqrt( 1.d0/rinv4 )
      dp = dsqrt( dp )
      s = 1.9023d0* dsign( 1.d0, zsign ) * dp/d + dels
c****
c**** first guess for closest point is
c****
c*    xp = x
c*    xp2 = xp*xp
c*    xp3 = xp2*xp
c*    xp4 = xp3*xp
c****
c**** calculate zp on curve for corresponding xp
c****
c*    zp = -( s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
c*   1   s7*xp4*xp3 + s8*xp4*xp4 )
c*    zsign = z-zp
c****
c**** slope of curve at xp, zp
c****
c*    do 4 i=1,3
c*    dzdxc = -(2.*s2*xp + 3.*s3*xp2+ 4.*s4*xp3 + 5.*s5*xp4 +
c*   1   6.*s6*xp4*xp + 7.*s7*xp4*xp2 + 8.*s8*xp4*xp3 )
c****
c**** next approximation to closest point is
c****
c*    xp = ( dzdxc*(z-zp)  +  dzdxc*dzdxc*xp + x ) / (1.+dzdxc*dzdxc)
c*    if( i  .eq.  1 )  xp = (3.*xp +  x ) / 4.
c*    xp2 = xp*xp
c*    xp3 = xp2*xp
c*    xp4 = xp3*xp
c*    zp = -( s2*xp2 + s3*xp3 + s4*xp4 + s5*xp4*xp + s6*xp4*xp2 +
c*   1   s7*xp4*xp3 + s8*xp4*xp4 )
c*  4 continue
c*    xxp = x-xp
c*    zzp = z-zp
c*    s = dsign( 1.d0,zsign) * dsqrt( xxp*xxp + zzp*zzp) / d + dels
c****
c****
c****
c****
      cs=c0+s*(c1+s*(c2+s*(c3+s*(c4+s*c5))))
      if( dabs(cs)  .gt.  70.  )  cs =dsign( 70.d0 ,cs  )
      e=dexp(cs)
      p0 = 1.0 + e
      db=bf-br
      bfld = 0.
      if( mtyp .eq. 3 ) bfld =
     1       br +( 1. - ndx*drr1 + bet1*drr2+gama*drr3+delt*drr4)*db/p0
      if( mtyp .eq. 4 ) bfld = br + ( 1./(1. +ndx*drr1) )*db/p0
c****
c**** print 100, x, y, z,  dr, s, bfld
c      print*, x, y, z,  dr, s, bfld
c*100 format( 1p6d15.4 )
c****
      return
      end

      subroutine dqfringe2(x,y,z,dl,bx,by,bz,fact1,fact2,qrad1,qrad2)
c***************************************************************************
c    dqfringe - computes overlap fringe fields for a dipole quadrupole     *
c    pair.  the dipole center exit is located at x=y=z=0. the optic axis   *
c    makes an angle th with the dipole face.  the quad is located a        *
c    distance dl along the optic axis from the dipole.                     *
c    diagrams in notebook 29 april 1987                                    *
c    modified version of dqfringe. altered to include dipolert 8/23/89 jjl *
c                                                                          *
c      fact1 describes dipole - bdip, th                                   *
c      fact2 describes quad - bquad, bhex, boct, bdec, bddec               *
c***************************************************************************
      real fact1(5),fact2(5),b(3),bpri(3)
c  pifac=pi/180
      data pifac/1.745329e-02/,ten/10./
      th=fact1(2)*pifac
c      write(7,*)' dl, th' ,dl,th
      bx=0.
      by=0.
      bz=0.
c  exit of d1 field
        xdum=x/ten
        ydum=y/ten
        zdum=z/ten
c      write(7,*)' xdum,ydum,zdum ',xdum,ydum,zdum
c        print*,' xdum,ydum,zdum ',xdum,ydum,zdum
        call dipolert2(xdum,zdum,ydum,bx,bz,by,fact1(1))
c        print*,bx,by,bz
c  get the right sign on the field
        bx=-bx
        by=-by
        bz=-bz
c      write(7,*)' dipole',bx,by,bz
c  entrance quad field
c   translate x & y axes
      dx=dl*sin(th)
      dy=dl*cos(th)
      ypri=y-dy
      xpri=x+dx
      zpri=z
c      write(7,*)' primed coord ',xpri,ypri,zpri
c      print*,dx,dy,dl,th,sin(th),cos(th)
c      print*,' primed coord ',xpri,ypri,zpri
c   rotate by th
      xdpri=(xpri*cos(th))+(ypri*sin(th))
      ydpri=(-xpri*sin(th))+(ypri*cos(th))
      zdpri=z
c      write(7,*)' dprimed coord ',xdpri,ydpri,zdpri
c      print*,' dprimed coord ',xdpri,ydpri,zdpri
      ydpri=-ydpri
      call mpoles(1,xdpri,zdpri,ydpri,fact2,qrad2,b(1),b(3),b(2))
      b(2)=-b(2)
c      write(7,*)' mpoles',(b(i),i=1,3)
c      print*,' mpoles',(b(i),i=1,3)
      bpri(1)=(b(1)*cos(-th))+(b(2)*sin(-th))
      bpri(2)=(-b(1)*sin(-th))+(b(2)*cos(-th))
      bpri(3)=b(3)
      bx=bx+bpri(1)
      by=by+bpri(2)
      bz=bz+bpri(3)
      return
      end

      subroutine mpoles(in,nx,ny,nz,ngrad,nrad,nbx,nby,nbz)
c stolen from raytrace jjl 10/8/86
c****                                                                   
c**** calculation of multipole(poles) field components                  
c****                                                                   
c****                                                                   
c****                                                                   
c**** 2 - quadrupole  (grad1)                                           
c**** 3 - hexapole    (grad2)                                           
c**** 4 - octapole    (grad3)                                           
c**** 5 - decapole    (grad4)                                           
c**** 6 - dodecapole  (grad5)                                           
c****                                                                   
c****                                                                   
      implicit real*8(a-h,o-z)                                          
      real nx,ny,nz,ngrad(5),nrad,nbx,nby,nbz
      common  /blck91/  c0, c1, c2, c3, c4, c5                          
c      data c0, c1, c2, c3, c4, c5/0.1122d0,8.5d0,-1.4982d0,3.5882d0,
c     &     -2.1209d0,1.7230d0/
c  above are not quite right -JJL 9/4/98
      data c0, c1, c2, c3, c4, c5/0.1039d0,6.27108d0,-1.51247d0,
     &     3.59946d0,-2.1323d0,1.7230d0/
      if (nrad.eq.0.)then
      write(6,*)' error in mpoles,  nrad= 0.'
      call exit(0)
      endif
      grad1 = -ngrad(1)/nrad
      grad2 =  ngrad(2)/nrad**2
      grad3 = -ngrad(3)/nrad**3
      grad4 =  ngrad(4)/nrad**4
      grad5 = -ngrad(5)/nrad**5
      rad=nrad
      d = 2. * rad                                                      
      frh  = 1.d0
      fro  = 1.d0
      frd  = 1.d0
      frdd = 1.d0
      dh  = frh *d
      do  = fro *d
      dd  = frd *d
      ddd = frdd*d
      x = nx
      y = ny
      z = nz
      x2 = x*x                                                          
      x3 = x2*x                                                         
      x4 = x3*x                                                         
      x5 = x4*x                                                         
      x6 = x5*x
      x7 = x6*x
      y2 = y*y                                                          
      y3 = y2*y                                                         
      y4 = y3*y                                                         
      y5 = y4*y                                                         
      y6 = y5*y
      y7 = y6*y
      go to ( 2, 1, 2 ) , in                                            
      print 3, in                                                       
    3 format( '  error in bpoles in= ', i5 ///)                          
      call exit(0)   
    1 continue                                                          
      b2x = grad1*y                                                     
      b2y = grad1*x                                                     
      b3x = grad2*2.*x*y                                                
      b3y = grad2*(x2-y2)                                               
      b4x = grad3*(3.*x2*y-y3)                                          
      b4y = grad3*(x3-3.*x*y2)                                          
      b5x = grad4*4.*(x3*y-x*y3)                                        
      b5y = grad4*(x4-6.*x2*y2+y4)                                      
      b6x = grad5*(5.*x4*y-10.*x2*y3+y5)                                
      b6y = grad5*(x5-10.*x3*y2+5.*x*y4)                                
      bx = b2x + b3x + b4x + b5x + b6x                                  
      by = b2y + b3y + b4y + b5y + b6y                                  
      bz = 0.                                                           
      bt =   sqrt( bx*bx + by*by )                                     
      nbx=bx
      nby=by
      nbz=bz
      return                                                            
c****
c****
c****
    2 s = z/d                                                           
      call bpls( 2, d, s, re, g1, g2, g3, g4, g5, g6 )
      b2x = grad1*( re*y - (g2/12.)*(3.*x2*y + y3) +                    
     1   (g4/384.)*(5.*x4*y + 6.*x2*y3 + y5 ) -                         
     2   (g6/23040.)*(7.*x6*y + 15.*x4*y3 + 9.*x2*y5 + y7)  )
      b2y = grad1*( re*x - (g2/12.)*(x3 + 3.*x*y2) +                    
     1   (g4/384.)*(x5 + 6.*x3*y2 + 5.*x*y4 ) -                         
     2   (g6/23040.)*(x7 + 9.*x5*y2 + 15.*x3*y4 + 7.*x*y6) )
      b2z = grad1*( g1*x*y - (g3/12.)*(x3*y + x*y3 ) +                  
     1   (g5/384.)*(x5*y +2.*x3*y3 + x*y5)  )
c****
c****
      ss = z/dh  + dsh
      call bpls( 3, dh, ss, re, g1, g2, g3, g4, g5, g6 )
      b3x = grad2*( re*2.*x*y - (g2/48.)*(12.*x3*y + 4.*x*y3 ) )        
      b3y = grad2*( re*(x2-y2) - (g2/48.)*(3.*x4 + 6.*x2*y2 - 5.*y4 ) ) 
      b3z = grad2*( g1*(x2*y - y3/3.) - (g3/48.)*(3.*x4*y+2.*x2*y3-y5)) 
c****
c****
      ss = z/do  + dso
      call bpls( 4, do, ss, re, g1, g2, g3, g4, g5, g6 )
      b4x = grad3*( re*(3.*x2*y - y3) - (g4/80.)*(20.*x4*y - 4.*y5 ) )  
      b4y = grad3*( re*(x3 - 3.*x*y2) - (g4/80.)*(4.*x5-20.*x*y4 ) )    
      b4z = grad3*g1*(x3*y - x*y3 )                                     
c****
c****
      ss = z/dd  + dsd
      call bpls( 5, dd, ss, re, g1, g2, g3, g4, g5, g6 )
      b5x = grad4*re*(4.*x3*y - 4.*x*y3)                                
      b5y = grad4*re*(x4 - 6.*x2*y2 + y4 )                              
      b5z = grad4*g1*(x4*y - 2.*x2*y3 + y5/5. )                         
c****
c****
      ss = z/ddd + dsdd
      call bpls( 6, ddd,ss, re, g1, g2, g3, g4, g5, g6 )
      b6x = grad5*re*(5.*x4*y - 10.*x2*y3 + y5 )                        
      b6y = grad5*re*(x5 - 10.*x3*y2 + 5.*x*y4 )                        
      b6z = 0.                                                          
c****
c****
      bx = b2x + b3x + b4x + b5x + b6x                                  
      by = b2y + b3y + b4y + b5y + b6y                                  
      bz = b2z + b3z + b4z + b5z + b6z                                  
      bt =   sqrt( bx*bx + by*by + bz*bz )                             
      nbx=bx
      nby=by
      nbz=bz
      return                                                            
      end                                                               
      subroutine bpls ( igp, d, s, re, g1, g2, g3, g4, g5, g6 )
c****
c****
c****
      implicit real*8 (a-h,o-z)
c****
c****
      common  /blck91/  c0, c1, c2, c3, c4, c5                          
c****
c****
      s2 = s*s
      s3 = s2*s
      s4 = s2*s2
      s5 = s4*s
      cs = c0 + c1*s + c2*s2 + c3*s3 + c4*s4 + c5*s5                    
      cp1 =(c1 + 2.*c2*s + 3.*c3*s2 + 4.*c4*s3 + 5.*c5*s4) / d          
      cp2 = (2.*c2 + 6.*c3*s + 12.*c4*s2 + 20.*c5*s3  ) / (d*d)         
      cp3 = ( 6.*c3 + 24.*c4*s + 60.*c5*s2 ) / (d**3)                   
      cp4 = ( 24.*c4 + 120.*c5*s ) / (d**4)                             
c****
      cp5 = 120.*c5/(d**5)
c****
c****
c****
      if( abs(cs) .gt. 70. )  cs = sign(70.d0, cs )                   
      e = exp(cs)                                                      
      re = 1./(1. + e)                                                  
      ere = e*re                                                        
      ere1= ere*re
      ere2= ere*ere1                                                    
      ere3= ere*ere2                                                    
      ere4= ere*ere3                                                    
c****
      ere5= ere*ere4
      ere6= ere*ere5
c****
c****
      cp12 = cp1*cp1                                                    
      cp13 = cp1*cp12                                                   
      cp14 = cp12*cp12                                                  
      cp22 = cp2*cp2                                                    
c****
      cp15 = cp12*cp13
      cp16 = cp13*cp13
      cp23 = cp2*cp22
      cp32 = cp3*cp3
c****
c****
      if( igp .eq. 6 ) return
      g1 = -cp1*ere1                                                    
c****
c****
      if( igp .eq. 5 ) return
      if( igp .eq. 4 ) go to 1
      g2 =-( cp2+cp12   )*ere1    + 2.*cp12 * ere2                      
      g3 =-(cp3 + 3.*cp1*cp2 + cp13  ) * ere1      +                    
     1   6.*(cp1*cp2 + cp13)*ere2 - 6.*cp13*ere3                        
c****
c****
      if( igp .eq. 3 ) return
1     g4 = -(cp4 + 4.*cp1*cp3 + 3.*cp22 + 6.*cp12*cp2 + cp14)*ere1  +   
     1   (8.*cp1*cp3 + 36.*cp12*cp2 + 6.*cp22 + 14.*cp14)*ere2    -     
     2   36.*(cp12*cp2 + cp14)*ere3       + 24.*cp14*ere4               
c****
c****
      if( igp .ne. 2 ) return
      g5 = (-cp5 - 5.*cp1*cp14 - 10.*cp2*cp3 - 10.*cp12*cp3 -
     1     15.*cp1*cp22 - 10.*cp13*cp2 - cp15)*ere1 +
     2     (10.*cp1*cp4 +20.*cp2*cp3 +60.*cp12*cp3 + 90.*cp1*cp22 +
     3     140.*cp13*cp2 +30.*cp15)*ere2 + (-60.*cp12*cp3 -
     4     90.*cp1*cp22 - 360.*cp13*cp2 - 150.*cp15)*ere3 +
     5     (240.*cp13*cp2 +240.*cp15)*ere4 + (-120.*cp15)*ere5
      g6 = (-6.*cp1*cp5 - 15.*cp2*cp4 - 15.*cp12*cp4 - 10.*cp32 -
     1     60.*cp1*cp2*cp3 - 20.*cp13*cp3 - 15.*cp23 - 45.*cp12*cp22 -
     2     15.*cp14*cp2 - cp16)*ere1 + (12.*cp1*cp5 + 30.*cp2*cp4 +
     3     90.*cp12*cp4 +20.*cp32 + 360.*cp1*cp2*cp3 +280.*cp13*cp3 +
     4     90.*cp23 + 630.*cp12*cp22 + 450.*cp14*cp2 + 62.*cp16)*ere2 +
     5     (-90.*cp12*cp4 - 360.*cp1*cp2*cp3 -720.*cp13*cp3 -90.*cp23 -
     6     1620.*cp12*cp22 -2250.*cp14*cp2 - 540.*cp16)*ere3 +
     7     (480.*cp13*cp3 + 1080.*cp12*cp22 + 3600.*cp14*cp2 +
     8     1560.*cp16)*ere4 + (-1800.*cp14*cp2 - 1800.*cp16)*ere5 +
     9     720.*cp16*ere6
c****
      return
      end


************************************

      FUNCTION dcosd(x)
      real*8 d2rad
      parameter (d2rad=3.141592654/180.0)
      real*8 dcosd,x

      dcosd = cos(x*d2rad)
      END

************************************ 

      FUNCTION dsind(x)
      real*8 d2rad
      parameter (d2rad=3.141592654/180.0)
      real*8 dsind,x

      dsind = sin(x*d2rad)
      END

************************************ 

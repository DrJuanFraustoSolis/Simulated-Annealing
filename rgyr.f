c**************************************************************
c
c This file contains the subroutines: rgyr
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************


      subroutine rgyr(nml, rgy, ee)

C CALCULATES THE RADIUS-OF-GYRATION AND THE END-TO-END DISTANCE
C FOR A GIVEN PROTEIN CONFORMATION
C If nml == 0, calculate the radius of gyration for all molecules
C
C     rgy  = radius-of-gyration
C     ee   = end-to-end distance
C
C REQUIREMENTS: c_alfa has to be called BEFORE call of this subroutine
C
C CALLS: NONE
C
      include 'INCL.H'
cf2py intent(in) nml
cf2py intent(out) rgy
cf2py intent(out) ee
      integer typ
      if (nml.eq.0) then
        nml1 = 1
        nml2 = ntlml
      else
        nml1 = nml
        nml2 = nml
      endif

        
      nat = iatrs2(irsml2(nml2))-iatrs1(irsml1(nml1))+1

      if (nat.le.0) then
        write (*,'(a,i4)')
     #     ' rgyr> No atoms found for molecule #',nml
        return 
      endif

      dn = dble(nat)
      dnp = 0.0
      dnh = 0.0

      dx =0.d0
      dxp = 0.0d0
      dxp = 0.0d0
      dy =0.d0
      dyp = 0.0d0
      dyh = 0.0d0
      dz =0.d0
      dzp = 0.0d0
      dzh = 0.0d0
      d2 =0.d0
      d2p = 0.0
      d2h = 0.0

      do i=iatrs1(irsml1(nml1)), iatrs1(irsml1(nml1)) + nat
         xi = xat(i)
         yi = yat(i)
         zi = zat(i)
         dx = dx + xi
         dy = dy + yi
         dz = dz + zi
         d2 = d2 + xi**2 + yi**2 + zi**2
         if(sigma(i).lt.0) then
            dxp = dxp + xi
            dyp = dyp + yi
            dzp = dzp + zi
            d2p = d2p + xi**2 + yi**2 + zi**2
            dnp = dnp + 1
         endif
         if(sigma(i).gt.0) then
            dxh = dxh + xi
            dyh = dyh + yi
            dzh = dzh + zi
            d2h = d2h + xi**2 + yi**2 + zi**2
            dnh = dnp + 1
         endif
      enddo

      dx = dx/dn
      dy = dy/dn
      dz = dz/dn
      d2 = d2/dn
      rg2 = d2 - (dx**2+dy**2+dz**2)

      dxp = dxp/dnp
      dyp = dyp/dnp
      dzp = dzp/dnp
      d2p = d2p/dnp
      rg2p = d2p - (dxp**2+dyp**2+dzp**2)

      dxh = dxh/dnh
      dyh = dyh/dnh
      dzh = dzh/dnh
      d2h = d2h/dnh

      rg2h = d2h - (dxh**2+dyh**2+dzh**2)

      if (rg2.gt.0.d0) then
        rgy  = sqrt(rg2)
      else
        rgy  = 0.d0
      endif
      
      if (rg2p.gt.0.d0) then
        rgyp  = sqrt(rg2p)
      else
        rgyp  = 0.d0
      endif

      if (rg2h.gt.0.d0) then
        rgyh  = sqrt(rg2h)
      else
        rgyh  = 0.d0
      endif

      i1=ind_alf(irsml1(nml1))  
      i2=ind_alf(irsml2(nml2))  

      ee = sqrt((xat(i2)-xat(i1))**2+(yat(i2)-yat(i1))**2
     #         +(zat(i2)-zat(i1))**2)

      return
      end


c **************************************************************
c
c This file contains the subroutines: opereg,gdtgbl,gdtreg
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      subroutine opereg(nml)

c .......................................................................
c PURPOSE: calculate regul. energy & it's partial derivatives
c          for molecule 'nml' vs. variables 'iv'
c
c  NB: if the unit axis for an internal variable coincides with a
c      global axis (i.e. for torsion or bond length variation round
c      or along 'xrfax', respectively, and bd. angle var. round
c      'zrfax'): VdW & 14 interaction partners of moving set atoms
c      should be used for calculation, instead of the mov. sets,
c      with opposite sign.
c
c      Example: By the the way the molecule-fixed system is set up,
c               changes in Phi_1 affect atomic positions BEFORE the
c               N-C^alpha bond relatively to the space-fixed system,
c               not the moving set of Phi_1.
c
c CALLS:    gdtgbl, gdtreg
c ......................................................................

      include 'INCL.H'
      include 'INCP.H'

      dimension xfat(mxat),yfat(mxat),zfat(mxat),
     #          xfrat(mxat),yfrat(mxat),zfrat(mxat),

     #          xfvr(mxvr),yfvr(mxvr),zfvr(mxvr),
     #          xfrvr(mxvr),yfrvr(mxvr),zfrvr(mxvr)

      logical   lnb


      ntlvr=nvrml(nml)
      if (ntlvr.eq.0) then
        write (*,'(a,i4)')
     #           ' opereg> No variables defined in molecule #',nml
        return
      endif

      ix2=ixrfpt(2,nml)    ! as indicator for situation noted above

      ifivr=ivrml1(nml)         ! 1st var. &
      ilavr=ifivr+ntlvr-1       ! last var. of 'nml'

c --------------------------- initializations
      do i=ifivr,ilavr
        gdeyrg(i)=0.d0
        xfvr(i)=0.d0
        yfvr(i)=0.d0
        zfvr(i)=0.d0
        xfrvr(i)=0.d0
        yfrvr(i)=0.d0
        zfrvr(i)=0.d0
      enddo

      ii=(nml-1)*6

      do i=ii+1,ii+6
        gdeygb(i) = 0.d0
      enddo

      x1=rfpt(1,nml)     ! r_1
      y1=rfpt(2,nml)
      z1=rfpt(3,nml)

      a= gbpr(4,nml)     ! alpha
      sa = sin(a)
      ca = cos(a)

      xk = yrfax(1,nml)  ! axis K
      yk = yrfax(2,nml)
      zk = yrfax(3,nml)

      eyrg = 0.d0

      do i=iatrs1(irsml1(nml)),iatrs2(irsml2(nml))

        j = ixatp(i)
        if (j.gt.0) then

          xi = xat(i) ! position of atom in internal model
          yi = yat(i)
          zi = zat(i)

          xji = xatp(j) - xi    ! x distance between internal model and PDB
          yji = yatp(j) - yi
          zji = zatp(j) - zi

          eyrg = eyrg + xji**2 + yji**2 + zji**2  ! The regularization energy is just
                                                  ! the sum over the atom distances 
                                                  ! squared.

          dx = 2.d0 * xji   ! f = - dE/dR_i
          dy = 2.d0 * yji   ! The factor of 2 comes from the derivative
          dz = 2.d0 * zji

c =============================================== global pars.

          gdeygb(ii+1) = gdeygb(ii+1) - dx   ! d(E_ij) / d(x_i)
          gdeygb(ii+2) = gdeygb(ii+2) - dy   ! d(E_ij) / d(y_i)
          gdeygb(ii+3) = gdeygb(ii+3) - dz   ! d(E_ij) / d(z_i)

c -------------------------- r = r_i - r_1
          x = xi - x1
          y = yi - y1
          z = zi - z1

          gdeygb(ii+4) = gdeygb(ii+4) +dx*y-dy*x                      !  d(E_ij) / d(a)

          gdeygb(ii+5) = gdeygb(ii+5) +z*(dy*ca-dx*sa)+dz*(x*sa-y*ca) !  d(E_ij) / d(b)

          gdeygb(ii+6) = gdeygb(ii+6) +dx*(zk*y-yk*z)+dy*(xk*z-zk*x)  !  d(E_ij) / d(g)
     #                                +dz*(yk*x-xk*y)

c =============================================== for internal vars.

          xfat(i) = dx
          yfat(i) = dy
          zfat(i) = dz

          xfrat(i) = dy*zi-dz*yi  ! g = f x r
          yfrat(i) = dz*xi-dx*zi  !
          zfrat(i) = dx*yi-dy*xi  !

        else
          xfat(i) = 0.d0
          yfat(i) = 0.d0
          zfat(i) = 0.d0
          xfrat(i) = 0.d0
          yfrat(i) = 0.d0
          zfrat(i) = 0.d0
        endif

      enddo ! atoms

      if (tesgrd) call gdtgbl(nml)

      i1s=imsml1(nml)+nmsml(nml)  ! last mov. set of 'nml' + 1
      i1a=iadml1(nml)+nadml(nml)  ! last added var. of 'nml' + 1

      do io=ilavr,ifivr,-1  ! ______ loop over vars in desc. order

        lnb = .false.  ! = true, if situation noted above takes place

        iv=iorvr(io)       ! index,
        it=ityvr(iv)       ! type,
        ia=iatvr(iv)       ! primary mov. atom,
        ib=iowat(ia)       ! "base" of current var.

        xb=xat(ib)
        yb=yat(ib)
        zb=zat(ib)

c ---------------------------------------- axis for var.

        if (it.eq.3) then      ! torsion

          ex= xtoat(ib)
          ey= ytoat(ib)
          ez= ztoat(ib)

          if (ib.eq.ix2) lnb = .true.

        elseif (it.eq.2) then  ! b.angle

          ex= xbaat(ia)
          ey= ybaat(ia)
          ez= zbaat(ia)

          if (ib.eq.ix2) lnb = .true.

        elseif (it.eq.1) then  ! b.length

          ex=xtoat(ia)
          ey=ytoat(ia)
          ez=ztoat(ia)

          if (ia.eq.ix2) lnb = .true.

        endif

        xfiv=0.0
        yfiv=0.0
        zfiv=0.0
        xfriv=0.0
        yfriv=0.0
        zfriv=0.0

        if (.not.lnb) then

          i2s=i1s-1      ! last m.s  &
          i1s=imsvr1(iv) ! 1st m.s for var. index 'iv'

          do ims=i1s,i2s  ! __ loop over moving sets

            i1=latms1(ims)   ! 1st &
            i2=latms2(ims)   ! last mov. atom in mov. set 'ims'

            do i=i1,i2  ! __ loop over atoms i ===================

              xfiv = xfiv + xfat(i)   ! f
              yfiv = yfiv + yfat(i)   !
              zfiv = zfiv + zfat(i)   !

              xfriv = xfriv + xfrat(i)  ! g
              yfriv = yfriv + yfrat(i)  !
              zfriv = zfriv + zfrat(i)  !

            enddo  ! ... atoms i
          enddo  ! ... m.s.

          i2a=i1a-1       ! last &
          i1a=iadvr1(iv)  ! 1st 'added' var. for 'iv'

          do iad=i1a,i2a  ! loop over add. var.

            lad=ladvr(iad)

            xfiv = xfiv + xfvr(lad)
            yfiv = yfiv + yfvr(lad)
            zfiv = zfiv + zfvr(lad)
            xfriv = xfriv + xfrvr(lad)
            yfriv = yfriv + yfrvr(lad)
            zfriv = zfriv + zfrvr(lad)

          enddo

        else

          do ivw=ivwat1(ia),ivwat2(ia)  ! vdW-domains of 'ia'
            do j=lvwat1(ivw),lvwat2(ivw)  ! .. their atoms

              xfiv = xfiv - xfat(j)
              yfiv = yfiv - yfat(j)
              zfiv = zfiv - zfat(j)

              xfriv = xfriv - xfrat(j)
              yfriv = yfriv - yfrat(j)
              zfriv = zfriv - zfrat(j)

            enddo
          enddo

          do i14=i14at1(ia),i14at2(ia)  ! 1-4 partn. of 'ia'
            j=l14at(i14)

            xfiv = xfiv - xfat(j)
            yfiv = yfiv - yfat(j)
            zfiv = zfiv - zfat(j)

            xfriv = xfriv - xfrat(j)
            yfriv = yfriv - yfrat(j)
            zfriv = zfriv - zfrat(j)

          enddo

        endif

        xfvr(iv) = xfiv
        yfvr(iv) = yfiv
        zfvr(iv) = zfiv
        xfrvr(iv) = xfriv
        yfrvr(iv) = yfriv
        zfrvr(iv) = zfriv

        if (it.eq.3.or.it.eq.2) then  ! torsion,b.angle

          gdeyrg(iv)= (ey*zb-ez*yb)*xfiv+(ez*xb-ex*zb)*yfiv+
     #                (ex*yb-ey*xb)*zfiv
     #               +ex*xfriv+ey*yfriv+ez*zfriv

        elseif (it.eq.1) then         ! b.length

          gdeyrg(iv)= -(ex*xfiv+ey*yfiv+ez*zfiv)

        endif

        if (tesgrd) call gdtreg(nml,iv)

      enddo  ! ... variables in desc. order

      return
      end
c **************************
      subroutine gdtgbl(nml)
C
C CALLS: bldmol,enyreg
c
c -------------------------- gradtest for 'gbpr'

      include 'INCL.H'

      parameter (del=1.d-7)


      ii=(nml-1)*6

      do i = 1,6

c ----------------------------- modify
        pro = gbpr(i,nml)
        gbpr(i,nml) = pro+del
        call bldmol(nml)

        gdn = ( enyreg(nml) - eyrg ) / del

        write (*,*) ' Gb. var #',(ii+i),': ',gdeygb(ii+i),gdn,
     #                                   abs(gdn-gdeygb(ii+i))
c ----------------------------- restore
        gbpr(i,nml) = pro
        call bldmol(nml)

      enddo  ! pars.

      return
      end
c *****************************
      subroutine gdtreg(nml,iv)

c .................................................................
c PURPOSE: calculate partial derivative of reg. energy for molecule
c          'nml' vs. variable 'iv' NUMERICALLY and compare with
c          its value obtained analytically
c
c CALLS:  setvar, enyreg
c .................................................................

      include 'INCL.H'

      parameter (del=1.d-6)

      dimension vlvrx(mxvr)

c ____________________________ get & save values of variables
      do i=1,ivrml1(ntlml)+nvrml(ntlml)-1
        it=ityvr(i)  ! type
        if (it.eq.3) then      ! torsion
          vlvrx(i)=toat(iatvr(i))
        elseif (it.eq.2) then  ! b.angle
          vlvrx(i)=baat(iatvr(i))
        elseif (it.eq.1) then  ! b.length
          vlvrx(i)=blat(iatvr(i))
        endif
      enddo

      ovr=vlvrx(iv)
      vlvrx(iv)=ovr+del       ! change variable 'iv' by 'del'
      call setvar(nml,vlvrx)

      eynw=enyreg(nml)

      gdn=(eynw-eyrg)/del    ! numerical derivative
      gda=gdeyrg(iv)         ! analytical der.

      write (*,'(1x,2a,2(e12.6,a))') nmvr(iv),': ',gda,' (',
     #       abs(gda-gdn),')'

c _________________________ restore vars
      vlvrx(iv)=ovr
      call setvar(nml,vlvrx)

      return
      end


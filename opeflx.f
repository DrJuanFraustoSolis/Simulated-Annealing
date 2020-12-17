c **************************************************************
c
c This file contains the subroutines: opeflx,gdtflx
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      subroutine opeflx(nml)

c ......................................................................
c PURPOSE: Calculate internal energy for FLEX dataset and its partial
c          derivatives vs. variables using recursive algorithm from:
c          Noguti T, Go N, J Phys Soc (Japan) v52 3685-3690 1984; Abe H,
c          Braun W, Noguti T, Go N, Comp Chem v8 239-247 1984; Mazur A K,
c          Abagyan R A, J Biomol Struct Dyn v6 815-832, which I modified
c          for atomic forces instead of simple derivatives (see Lavery R,
c          Sklenar H, Zakrzewska K, Pullman B, J Biomol Struct Dyn v3
c          989-1014 1986)
c
c CALLS:   gdtflx
c ......................................................................

      include 'INCL.H'

      dimension xfat(mxat),yfat(mxat),zfat(mxat),
     #          xtat(mxat),ytat(mxat),ztat(mxat),
     #          xfvr(mxvr),yfvr(mxvr),zfvr(mxvr),
     #          xfrvr(mxvr),yfrvr(mxvr),zfrvr(mxvr)


      eyel=0.d0
      eyvw=0.d0
      eyhb=0.d0
      eyvr=0.d0
      eysm=0.d0

      ntlvr=nvrml(nml)
      if (ntlvr.eq.0) then
        write (*,'(a,i4)')
     #           ' opeflx> No variables defined in molecule #',nml
        return
      endif

      ifivr=ivrml1(nml)
      ilavr=ifivr+ntlvr-1

      do i=ifivr,ilavr
        gdeyvr(i)=0.d0
        xfvr(i)=0.0
        yfvr(i)=0.0
        zfvr(i)=0.0
        xfrvr(i)=0.0
        yfrvr(i)=0.0
        zfrvr(i)=0.0
      enddo

      do i=iatrs1(irsml1(nml)),iatrs2(irsml2(nml))
        xfat(i)=0.d0
        yfat(i)=0.d0
        zfat(i)=0.d0

        xtat(i)=0.d0
        ytat(i)=0.d0
        ztat(i)=0.d0
      enddo

      i1s=imsml1(nml)+nmsml(nml)
      i1a=iadml1(nml)+nadml(nml)

      do io=ilavr,ifivr,-1  ! ______ loop over variables in desc. order

        iv=iorvr(io)       ! index of var.
        ia=iatvr(iv)       ! prim.mv.at

        ib=iowat(ia)       ! base
        xb=xat(ib)
        yb=yat(ib)
        zb=zat(ib)

        it=ityvr(iv)       ! type
        ic=iclvr(iv)       ! class

        fvr=0.d0

        if (it.eq.3) then      ! torsion
          ex=xtoat(ib)
          ey=ytoat(ib)
          ez=ztoat(ib)

          vr=toat(ia)

          e0=e0to(ic)
          if (e0.ne.0.d0) then
            vrn=vr*rnto(ic)
            eyvr=eyvr+e0*(1.d0+sgto(ic)*cos(vrn))
            fvr=esnto(ic)*sin(vrn)
          endif

        elseif (it.eq.2) then  ! b.angle
          ex=xbaat(ia)
          ey=ybaat(ia)
          ez=zbaat(ia)

          vr=baat(ia)

        elseif (it.eq.1) then  ! b.length
          ex=xtoat(ia)
          ey=ytoat(ia)
          ez=ztoat(ia)

          vr=blat(ia)

        endif

c ============================================ Energies & Atomic forces

        xfiv=0.d0
        yfiv=0.d0
        zfiv=0.d0
        xfriv=0.d0
        yfriv=0.d0
        zfriv=0.d0

        i2s=i1s-1      ! last m.s per 'iv'
        i1s=imsvr1(iv) ! 1st m.s

        do ims=i1s,i2s  ! __ loop over m.s

          i1=latms1(ims)
          i2=latms2(ims)

          do i=i1,i2  ! __ loop over atoms i ===================

            ity=ityat(i)
            cqi=conv*cgat(i)

            xi=xat(i)
            yi=yat(i)
            zi=zat(i)

            xfi=xfat(i)
            yfi=yfat(i)
            zfi=zfat(i)

            xti=xtat(i)
            yti=ytat(i)
            zti=ztat(i)

            do ivw=ivwat1(i),ivwat2(i)  ! loop over vdW-domains of 'i'
              do j=lvwat1(ivw),lvwat2(ivw)  ! atoms j

                jty=ityat(j)

                xj=xat(j)
                yj=yat(j)
                zj=zat(j)

                xij=xj-xi
                yij=yj-yi
                zij=zj-zi

                rij2=xij*xij+yij*yij+zij*zij
                rij6=rij2**3
                rij12=rij6*rij6
                rij=sqrt(rij2)

                cqiqj=cqi*cgat(j)

                if (epsd) then

                  sr=slp_f*rij
                  sr2=sr*sr
                  xsr=(plt-1.d0)*exp(-sr)/2.d0
                  ep=plt-(sr2+2.d0*sr+2.d0)*xsr

                  eel=cqiqj/(rij*ep)
                  eyel=eyel+eel
                  deel=eel+cqiqj*(slp_f*sr2*xsr)/(ep*ep)

                else
              
                  eel=cqiqj/rij
                  eyel=eyel+eel
                  deel=eel

                endif


                eyrp=aij(ity,jty)/rij12
                eyds=cij(ity,jty)/rij6

                c=(-12.d0*eyrp+6.d0*eyds- deel)/rij

                xij=xij/rij
                yij=yij/rij
                zij=zij/rij

                xfji=c*xij
                yfji=c*yij
                zfji=c*zij

                ijhb=ihbty(ity,jty)
                if (ijhb.ne.0.and.rij.le.cohb) then  ! HB Possible
                  if (ijhb.gt.0) then  ! i=H,j=acceptor
                    iowh=iowat(i)
                    px=xi-xat(iowh)
                    py=yi-yat(iowh)
                    pz=zi-zat(iowh)
                  else                 ! i=acceptor,j=H
                    jowh=iowat(j)
                    px=xat(jowh)-xj
                    py=yat(jowh)-yj
                    pz=zat(jowh)-zj
                  endif

                  p=sqrt(px*px+py*py+pz*pz)
                  px=px/p
                  py=py/p
                  pz=pz/p

                  cth=xij*px+yij*py+zij*pz

                  if (cth.gt.0.d0) then

                    deyrp=(ahb(ity,jty)-aij(ity,jty))/rij12
                    deyds=(chb(ity,jty)-cij(ity,jty))/rij6

                    dhb=deyrp-deyds
                    eyhb=eyhb+eyrp-eyds+cth*dhb

                    if (ijhb.gt.0) then  ! i=H
                      xti=xti +dhb * (zij*py-yij*pz)
                      yti=yti +dhb * (xij*pz-zij*px)
                      zti=zti +dhb * (yij*px-xij*py)
                    else                 ! j=H
                      xtat(j)=xtat(j) +dhb * (zij*py-yij*pz)
                      ytat(j)=ytat(j) +dhb * (xij*pz-zij*px)
                      ztat(j)=ztat(j) +dhb * (yij*px-xij*py)
                    endif

                    dhb=dhb/rij
                    hhb=cth*(7.d0*deyds-13.d0*deyrp)/rij

                    xfji=xfji+ dhb*px+ hhb*xij
                    yfji=yfji+ dhb*py+ hhb*yij
                    zfji=zfji+ dhb*pz+ hhb*zij
c __________________________________________________ No Hydrogen Bond
                  else
                    eyvw=eyvw+eyrp-eyds
                  endif
                else
                  eyvw=eyvw+eyrp-eyds
                endif

                xfi=xfi+xfji
                yfi=yfi+yfji
                zfi=zfi+zfji

                xfat(j)=xfat(j)-xfji
                yfat(j)=yfat(j)-yfji
                zfat(j)=zfat(j)-zfji

              enddo  ! ... atoms j
            enddo  ! ... vdW-domains of i

            do i14=i14at1(i),i14at2(i)  ! loop over 1-4 partn. of 'i'
              j=l14at(i14)

              jty=ityat(j)

              xj=xat(j)
              yj=yat(j)
              zj=zat(j)

              xij=xj-xi
              yij=yj-yi
              zij=zj-zi

              rij2=xij*xij+yij*yij+zij*zij
              rij6=rij2**3
              rij12=rij6*rij6
              rij=sqrt(rij2)

              cqiqj=cqi*cgat(j)

              if (epsd) then

                sr=slp_f*rij
                sr2=sr*sr
                xsr=(plt-1.d0)*exp(-sr)/2.d0
                ep=plt-(sr2+2.d0*sr+2.d0)*xsr

                eel=cqiqj/(rij*ep)
                eyel=eyel+eel
                deel=eel+cqiqj*(slp_f*sr2*xsr)/(ep*ep)

              else

                eel=cqiqj/rij
                eyel=eyel+eel
                deel=eel

              endif

              eyrp=a14(ity,jty)/rij12
              eyds=cij(ity,jty)/rij6

              c=(-12.d0*eyrp+6.d0*eyds- deel )/rij

              xij=xij/rij
              yij=yij/rij
              zij=zij/rij

              xfji=c*xij
              yfji=c*yij
              zfji=c*zij

              ijhb=ihbty(ity,jty)
              if (ijhb.ne.0.and.rij.le.cohb) then  ! HB Possible
                if (ijhb.gt.0) then  ! i=H,j=acceptor
                  iowh=iowat(i)
                  px=xi-xat(iowh)
                  py=yi-yat(iowh)
                  pz=zi-zat(iowh)
                else                 ! i=acceptor,j=H
                  jowh=iowat(j)
                  px=xat(jowh)-xj
                  py=yat(jowh)-yj
                  pz=zat(jowh)-zj
                endif

                p=sqrt(px*px+py*py+pz*pz)
                px=px/p
                py=py/p
                pz=pz/p

                cth=xij*px+yij*py+zij*pz

                if (cth.gt.0.d0) then

                  deyrp=(ahb(ity,jty)-a14(ity,jty))/rij12
                  deyds=(chb(ity,jty)-cij(ity,jty))/rij6

                  dhb=deyrp-deyds
                  eyhb=eyhb+eyrp-eyds+cth*dhb

                  if (ijhb.gt.0) then  ! i=H
                    xti=xti -dhb * (yij*pz-zij*py)
                    yti=yti -dhb * (zij*px-xij*pz)
                    zti=zti -dhb * (xij*py-yij*px)
                  else                 ! j=H
                    xtat(j)=xtat(j) +dhb * (yij*pz-zij*py)
                    ytat(j)=ytat(j) +dhb * (zij*px-xij*pz)
                    ztat(j)=ztat(j) +dhb * (xij*py-yij*px)
                  endif

                  dhb=dhb/rij
                  hhb=cth*(7.d0*deyds-13.d0*deyrp)/rij

                  xfji=xfji+ dhb*px+ hhb*xij
                  yfji=yfji+ dhb*py+ hhb*yij
                  zfji=zfji+ dhb*pz+ hhb*zij
c __________________________________________________ No Hydrogen Bond
                else
                  eyvw=eyvw+eyrp-eyds
                endif
              else
                eyvw=eyvw+eyrp-eyds
              endif

              xfi=xfi+xfji
              yfi=yfi+yfji
              zfi=zfi+zfji

              xfat(j)=xfat(j)-xfji
              yfat(j)=yfat(j)-yfji
              zfat(j)=zfat(j)-zfji

            enddo  ! ... 1-4-partners of i

            xfat(i)=xfi
            yfat(i)=yfi
            zfat(i)=zfi
            xtat(i)=xti
            ytat(i)=yti
            ztat(i)=zti

            xfiv=xfiv + xfi
            yfiv=yfiv + yfi
            zfiv=zfiv + zfi

            xfriv=xfriv + yfi*zi-zfi*yi + xti
            yfriv=yfriv + zfi*xi-xfi*zi + yti
            zfriv=zfriv + xfi*yi-yfi*xi + zti

          enddo  ! ... atoms i
        enddo  ! ... m.s.

        i2a=i1a-1       ! last 'added' var.
        i1a=iadvr1(iv)  ! 1st 'added' var.

        do iad=i1a,i2a
          lad=ladvr(iad)
          xfiv=xfiv+xfvr(lad)
          yfiv=yfiv+yfvr(lad)
          zfiv=zfiv+zfvr(lad)
          xfriv=xfriv+xfrvr(lad)
          yfriv=yfriv+yfrvr(lad)
          zfriv=zfriv+zfrvr(lad)
        enddo

        xfvr(iv)=xfiv
        yfvr(iv)=yfiv
        zfvr(iv)=zfiv
        xfrvr(iv)=xfriv
        yfrvr(iv)=yfriv
        zfrvr(iv)=zfriv

        if (it.eq.3.or.it.eq.2) then  ! torsion,b.angle

          gdeyvr(iv)= (ey*zb-ez*yb)*xfiv+(ez*xb-ex*zb)*yfiv+
     #                (ex*yb-ey*xb)*zfiv
     #               +ex*xfriv+ey*yfriv+ez*zfriv -fvr

        elseif (it.eq.1) then         ! b.length

          gdeyvr(iv)= -(ex*xfiv+ey*yfiv+ez*zfiv) -fvr

        endif

        if (tesgrd) call gdtflx(nml,iv)   !  grad.-test

      enddo  ! ... variables

      eysm= eyel+eyvw+eyhb+eyvr

      return
      end
c *****************************
      subroutine gdtflx(nml,iv)

c .....................................................................
c PURPOSE: calculate partial derivative of internal energy for molecule
c          'nml' vs. variable 'iv' NUMERICALLY and compare with
c          its value obtained analytically
c
c CALLS:  setvar, enyflx
c .....................................................................

      include 'INCL.H'

      parameter (del=1.d-7)

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
      eyol=enyflx(nml)

      vlvrx(iv)=ovr+del       ! change variable 'iv' by 'del'
      call setvar(nml,vlvrx)
      eynw=enyflx(nml)       ! new energy

      gdn=(eynw-eyol)/del    ! numerical derivative
      gda=gdeyvr(iv)         ! analytical der.

      write (*,'(1x,2a,2(e12.6,a))') nmvr(iv),': ',gda,' (',
     #       abs(gda-gdn),')'

c _________________________ restore
      vlvrx(iv)=ovr
      call setvar(nml,vlvrx)

      return
      end


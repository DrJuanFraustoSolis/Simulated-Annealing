c **************************************************************
c
c This file contains the subroutines: opeshe,gdtshe
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      subroutine opeshe(nml)

c ......................................................................
c PURPOSE: Calculate internal energy for ECEPP/3 dataset and its partial
c          derivatives vs. variables using recursive algorithm from:
c          Noguti T, Go N, J Phys Soc (Japan) v52 3685-3690 1984; Abe H,
c          Braun W, Noguti T, Go N, Comp Chem v8 239-247 1984; Mazur A K,
c          Abagyan R A, J Biomol Struct Dyn v6 815-832, which I modified
c          for atomic forces instead of simple derivatives (see Lavery R,
c          Sklenar H, Zakrzewska K, Pullman B, J Biomol Struct Dyn v3
c          989-1014 1986)
c
c CALLS:   gdtshe 
c ......................................................................

      include 'INCL.H'

      dimension xfat(mxat),yfat(mxat),zfat(mxat),
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
     #           ' opeshe> No variables defined in molecule #',nml
        return
      endif

      ifivr=ivrml1(nml)
      ilavr=ifivr+ntlvr-1

      do i=ifivr,ilavr
        gdeyvr(i)=0.d0
        xfvr(i)=0.d0
        yfvr(i)=0.d0
        zfvr(i)=0.d0
        xfrvr(i)=0.d0
        yfrvr(i)=0.d0
        zfrvr(i)=0.d0
      enddo

      do i=iatrs1(irsml1(nml)),iatrs2(irsml2(nml))
        xfat(i)=0.d0
        yfat(i)=0.d0
        zfat(i)=0.d0
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
            fvr=esnto(ic)*sin(vrn)            ! FORCE from variable
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

            do ivw=ivwat1(i),ivwat2(i)  ! loop over vdW-domains of 'i'

              do j=lvwat1(ivw),lvwat2(ivw)  ! atoms j

                jty=ityat(j)

                xij=xat(j)-xi
                yij=yat(j)-yi
                zij=zat(j)-zi

                rij2=xij*xij+yij*yij+zij*zij
                rij4=rij2*rij2
                rij6=rij4*rij2
                rij = sqrt(rij2)

                if (epsd) then

                  sr=slp*rij
                  sr2=sr*sr
                  xsr=(plt-1.d0)*exp(-sr)/2.d0
                  ep=plt-(sr2+2.d0*sr+2.d0)*xsr
                  eel=cqi*cgat(j)/(rij*ep)
                  deel=eel+cqi*cgat(j)*(slp*sr2*xsr)/(ep*ep)

                else

                  eel=cqi*cgat(j)/rij
                  deel=eel

                endif

                eyel=eyel+eel

                if (ihbty(ity,jty).ne.0) then
                  eyrp=ahb(ity,jty)/(rij6*rij6)
                  eyds=chb(ity,jty)/(rij6*rij4)
                  eyhb=eyhb+eyrp-eyds
                  c=(-12.d0*eyrp+10.d0*eyds-deel)/rij2
                else
                  eyrp=aij(ity,jty)/(rij6*rij6)
                  eyds=cij(ity,jty)/rij6
                  eyvw=eyvw+eyrp-eyds
                  c=(-12.d0*eyrp+6.d0*eyds-deel)/rij2
                endif

                xfji=c*xij
                yfji=c*yij
                zfji=c*zij

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

              xij=xat(j)-xi
              yij=yat(j)-yi
              zij=zat(j)-zi

              rij2=xij*xij+yij*yij+zij*zij
              rij4=rij2*rij2
              rij6=rij4*rij2
              rij = sqrt(rij2)

              if (epsd) then

                sr=slp*rij
                sr2=sr*sr
                xsr=(plt-1.d0)*exp(-sr)/2.d0
                ep=plt-(sr2+2.d0*sr+2.d0)*xsr
                eel=cqi*cgat(j)/(rij*ep)
                deel=eel+cqi*cgat(j)*(slp*sr2*xsr)/(ep*ep)

              else

                eel=cqi*cgat(j)/rij
                deel=eel

              end if

              eyel=eyel+eel

              if (ihbty(ity,jty).ne.0) then
                eyrp=ahb(ity,jty)/(rij6*rij6)
                eyds=chb(ity,jty)/(rij6*rij4)
                eyhb=eyhb+eyrp-eyds
                c=(-12.d0*eyrp+10.d0*eyds-deel)/rij2
              else
                eyrp=a14(ity,jty)/(rij6*rij6)
                eyds=cij(ity,jty)/rij6
                eyvw=eyvw+eyrp-eyds
                c=(-12.d0*eyrp+6.d0*eyds-deel)/rij2
              endif

              xfji=c*xij
              yfji=c*yij
              zfji=c*zij

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

            xfiv=xfiv + xfi
            yfiv=yfiv + yfi
            zfiv=zfiv + zfi

            xfriv=xfriv + yfi*zi-zfi*yi
            yfriv=yfriv + zfi*xi-xfi*zi
            zfriv=zfriv + xfi*yi-yfi*xi

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

        if (tesgrd) call gdtshe(nml,iv)

      enddo  ! ... variables

      eysm= eyel+eyvw+eyhb+eyvr

      return
      end
c *****************************
      subroutine gdtshe(nml,iv)

c .....................................................................
c PURPOSE: calculate partial derivative of internal energy for molecule
c          'nml' vs. variable 'iv' NUMERICALLY and compare with
c          its value obtained analytically
c
c CALLS:  setvar, enyshe
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
      eyol=enyshe(nml)

      vlvrx(iv)=ovr+del       ! change variable 'iv' by 'del'
      call setvar(nml,vlvrx)
      eynw=enyshe(nml)       ! new energy

      gdn=(eynw-eyol)/del    ! numerical derivative
      gda=gdeyvr(iv)         ! analytical der.

      write (*,'(1x,2a,2(e12.6,a))') nmvr(iv),': ',gda,' (',
     #       abs(gda-gdn),')'

c _________________________ restore
      vlvrx(iv)=ovr
      call setvar(nml,vlvrx)

      return
      end


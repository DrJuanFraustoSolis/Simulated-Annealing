c **************************************************************
c
c This file contains the subroutines: enyflx
c
c Copyright 2003       Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************


      real*8 function enyflx(nml)

c .......................................................................
c
c  PURPOSE: Calculate internal energy of molecule 'nml' with FLEX dataset
c
c  CALLS: none
c
c .......................................................................

      include 'INCL.H'

      ntlvr=nvrml(nml)
      if (ntlvr.eq.0) then
        write (*,'(a,i4)')
     #           ' enyflx> No variables defined in molecule #',nml
        return
      endif

      enyflx=0.0
      eyel=0.0
      eyvw=0.0
      eyhb=0.0
      eyvr=0.0

      ifivr=ivrml1(nml)
      i1s=imsml1(nml)+nmsml(nml)

      do io=ifivr+ntlvr-1,ifivr,-1  ! ______ over variables in desc. order
        iv=iorvr(io)       ! index of var.

        ia=iatvr(iv)       ! prim.mv.at
        it=ityvr(iv)       ! type
        ic=iclvr(iv)       ! class

        if (it.eq.3) then      ! torsion
          vr=toat(ia)

          e0=e0to(ic)
          if (e0.ne.0.0) eyvr=eyvr+e0*(1.0+sgto(ic)*cos(vr*rnto(ic)))

        elseif (it.eq.2) then  ! b.angle
          vr=baat(ia)
        elseif (it.eq.1) then  ! b.length
          vr=blat(ia)
        endif

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

            do ivw=ivwat1(i),ivwat2(i)  ! over vdW-domains of 'i'
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
                if(epsd) then
c --------------------------------- distance dependent dielectric constant
                sr=slp_f*rij
                ep=plt-(sr*sr+2.0*sr+2.0)*(plt-1.0)*exp(-sr)/2.0
                else
                ep=1.0d0
                end if

                eyel=eyel+cqi*cgat(j)/(rij*ep)

                evw=aij(ity,jty)/rij12-cij(ity,jty)/rij6

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

                  cth=(xij*px+yij*py+zij*pz)/(rij*
     #                 sqrt(px*px+py*py+pz*pz))

                  if (cth.gt.0.0) then
                    eyhb=eyhb+ evw + cth*(
     #                   (ahb(ity,jty)-aij(ity,jty))/rij12-
     #                   (chb(ity,jty)-cij(ity,jty))/rij6 )
                  else                             ! No Hydrogen Bond
                    eyvw=eyvw + evw
                  endif
                else
                  eyvw=eyvw + evw
                endif

              enddo  ! ... atoms j
            enddo  ! ... vdW-domains of i

            do i14=i14at1(i),i14at2(i)   !  over 1-4 partners of 'i'
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
              if(epsd) then
c --------------------------------- distance dependent dielectric constant
              sr=slp_f*rij
              ep=plt-(sr*sr+2.0*sr+2.0)*(plt-1.)*exp(-sr)/2.0
              else
              ep=1.0d0
              end if

              eyel=eyel+ cqi*cgat(j)/(rij*ep)

              evw=a14(ity,jty)/rij12-cij(ity,jty)/rij6

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

                cth=(xij*px+yij*py+zij*pz)/(rij*
     #               sqrt(px*px+py*py+pz*pz))

                if (cth.gt.0.0) then
                  eyhb=eyhb+ evw + cth*(
     #                 (ahb(ity,jty)-a14(ity,jty))/rij12-
     #                 (chb(ity,jty)-cij(ity,jty))/rij6 )
                else                             ! No Hydrogen Bond
                  eyvw=eyvw + evw
                endif
              else
                eyvw=eyvw + evw
              endif

            enddo  ! ... 1-4-partners of i

          enddo  ! ... atoms i
        enddo  ! ... m.s.

      enddo  ! ... variables

      eysm = eyel + eyvw + eyhb + eyvr
      enyflx=eysm
      return
      end


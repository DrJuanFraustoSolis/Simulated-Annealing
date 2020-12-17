c **************************************************************
c
c This file contains the subroutines: cnteny
c
c Copyright 2005       Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      subroutine cnteny(nml)

c ................................................................................
c PURPOSE: Calculate atomic contact energy of molecule 'nml' with ECEPP parameters
c
c CALLS: nursat
c ................................................................................

      include 'INCL.H'

      parameter (coeycn=2.d0)  ! min. cont. energy to display 

      dimension eyatcn(mxat),idxat(mxat)


      ieyel=0  ! =1, if count electrost. energy

      ntlvr=nvrml(nml)

      if (ntlvr.eq.0) then
        write (*,'(a,i4)')
     #           ' cnteny> No variables defined in molecule #',nml
        return
      endif

      if (ieyel.eq.1) then
        coey=5.d0*coeycn
      else
        coey=coeycn
      endif

      iat1=iatrs1(irsml1(nml))-1   ! last atom before 'nml'
      nat=iatrs2(irsml2(nml))-iat1 ! no. of atoms in 'nml'

      do i=1,nat
        eyatcn(i)=0.d0
      enddo


      ifivr=ivrml1(nml)
      i1s=imsml1(nml)+nmsml(nml)

      do io=ifivr+ntlvr-1,ifivr,-1  ! ______ over variables in desc. order
        iv=iorvr(io)       ! index of var.

        ia=iatvr(iv)       ! prim.mv.at
        it=ityvr(iv)       ! type
        ic=iclvr(iv)       ! class

        i2s=i1s-1      ! last m.s per 'iv'
        i1s=imsvr1(iv) ! 1st m.s

        do ims=i1s,i2s  ! __ loop over m.s
          i1=latms1(ims)
          i2=latms2(ims)

          do i=i1,i2  ! __ over atoms i ===================

            ii=i-iat1

            ity=ityat(i)
            cqi=conv*cgat(i)

            xi=xat(i)
            yi=yat(i)
            zi=zat(i)

            do ivw=ivwat1(i),ivwat2(i)  !  vdW-domains of 'i'
              do j=lvwat1(ivw),lvwat2(ivw)  !  atoms j

                jj=j-iat1

                jty=ityat(j)

                xij=xat(j)-xi
                yij=yat(j)-yi
                zij=zat(j)-zi

                rij2=xij*xij+yij*yij+zij*zij
                rij4=rij2*rij2
                rij6=rij4*rij2

                if (ieyel.eq.1) then

                  rij=sqrt(rij2)

                  if (epsd) then
                    sr=slp*rij
                    ep=plt-(sr*sr+2.0*sr+2.0)*(plt-1.0)*exp(-sr)/2.0
                  else
                    ep = 1.0d0
                  endif

                  ey=cqi*cgat(j)/(rij*ep)

                  eyatcn(ii)=eyatcn(ii)+.5d0*ey
                  eyatcn(jj)=eyatcn(jj)+.5d0*ey

                endif ! eyel

                if (ihbty(ity,jty).eq.0) then
                  ey=aij(ity,jty)/(rij6*rij6)-cij(ity,jty)/rij6
                else ! HB
                  ey=ahb(ity,jty)/(rij6*rij6)-chb(ity,jty)/(rij6*rij4)
                endif

                eyatcn(ii)=eyatcn(ii)+.5d0*ey
                eyatcn(jj)=eyatcn(jj)+.5d0*ey

              enddo  ! ... atoms j
            enddo  ! ... vdW-domains of i

            do i14=i14at1(i),i14at2(i)   !  over 1-4 partn. of 'i'
              j=l14at(i14)

              jj=j-iat1

              jty=ityat(j)

              xij=xat(j)-xi
              yij=yat(j)-yi
              zij=zat(j)-zi
              rij2=xij*xij+yij*yij+zij*zij
              rij4=rij2*rij2
              rij6=rij4*rij2

              if (ieyel.eq.1) then

                rij = sqrt(rij2)

                if (epsd) then
                  sr=slp*rij
                  ep=plt-(sr*sr+2.0*sr+2.0)*(plt-1.0)*exp(-sr)/2.0
                else
                  ep=1.0d0
                endif

                ey=cqi*cgat(j)/(rij*ep)

                eyatcn(ii)=eyatcn(ii)+.5d0*ey
                eyatcn(jj)=eyatcn(jj)+.5d0*ey

              endif  ! eel

              if (ihbty(ity,jty).eq.0) then
                ey=a14(ity,jty)/(rij6*rij6)-cij(ity,jty)/rij6
              else
                ey=ahb(ity,jty)/(rij6*rij6)-chb(ity,jty)/(rij6*rij4)
              endif

              eyatcn(ii)=eyatcn(ii)+.5d0*ey
              eyatcn(jj)=eyatcn(jj)+.5d0*ey

            enddo  ! ... 1-4-partners of i

          enddo  ! ... atoms i
        enddo  ! ... m.s.

      enddo  ! ... variables

      nbc=0

      do i=1,nat
        ey=eyatcn(i)
        if (ey.gt.coey) then
          nbc=nbc+1
          ir=nursat(i)
          write(*,'(1x,i4,1x,a4,1x,a4,a2,e11.4)') ir,seq(ir),nmat(i),
     #                                            ': ',ey
        endif
      enddo

      return
      end


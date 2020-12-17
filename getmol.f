c **************************
c **************************************************************
c
c This file contains the subroutines: getmol,redres
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************


      subroutine getmol(nml)

c ...................................................................
c PURPOSE:  assemble data for molecule 'nml' according to
c           its sequence using residue library 'reslib'
c
c ! Molecules must be assembled in sequential order (1 -> ntlml)
c            (or number of atoms & variables must remain the same)
c
c INPUT:    irsml1(nml),irsml2(nml),seq(irsml1()...irsml2())
c           nml>1: irsml1(nml-1),iatrs2(irsml2(nml-1))
c                  ivrrs1(irsml2(nml-1)),nvrrs(irsml2(nml-1))
c
c OUTPUT:   molecule  - ivrml1,nvrml
c           residues  - iatrs1,ixatrs,iatrs2,ivrrs1,nvrrs
c           atoms     - nmat,ityat,cgat,blat,baat,csbaat,snbaat,
c                       toat,cstoat,sntoat
c           bonds     - nbdat,iowat,iyowat,ibdat(1-mxbd,),iybdat(1-mxbd,)
c                       ! 1st atom of 'nml': iowat indicates 1st bond
c                          to a FOLLOWING atom (not previous) !
c           variables - ityvr,iclvr,iatvr,nmvr
c                                    
c CALLS:    iopfil,redres,iendst
c ...................................................................

      include 'INCL.H'

      character res*4

      if (iopfil(lunlib,reslib,'old','formatted').le.izero) then
        write (*,'(a,/,a,i3,2a)') 
     #    ' getmol> ERROR opening library of residues:',
     #    ' LUN=',lunlib,' FILE=',reslib(1:iendst(reslib))
        stop
      endif

      if (nml.eq.1) then
        ntlat=0
        ntlvr=0
      else
        i=irsml2(nml-1)            ! last res. of previous mol.
        ntlat=iatrs2(i)
        ntlvr=ivrrs1(i)+nvrrs(i)-1
      endif
      ivrml1(nml)=ntlvr + 1

      ilars=irsml2(nml)
      ifirs=irsml1(nml)

      do nrs=ifirs,ilars                 ! Residues in molecule

        res = seq(nrs)
        call tolost(res)  ! ensure lower case for residue name

        if (res(:3).eq.'nme'.and.nrs.ne.ilars) then
          write (*,'(3a)') ' getmol> residue >',res,
     #                     '< allowed at C-terminus only !'
          close(lunlib)
          stop
        elseif (res(:3).eq.'ace'.and.nrs.ne.ifirs) then
          write (*,'(3a)') ' getmol> residue >',res,
     #                     '< allowed at N-terminus only !'
          close(lunlib)
          stop
        endif

        call redres(res,nat,nxt,nvr)

        if ((nat+ntlat).gt.mxat) then
          write (*,'(a,i5)') ' getmol> number of atoms > ',mxat
          close(lunlib)
          stop
        endif
        if ((nvr+ntlvr).gt.mxvr) then
          write (*,'(a,i5)') ' getmol> number of variables > ',mxvr
          close(lunlib)
          stop
        endif

        rewind lunlib

c ___________________________________________________________ Atoms
        do i=1,nat
          n=i+ntlat
          nmat(n)=nmath(i)
          ityat(n)=ityath(i)
          cgat(n)=cgath(i)
          blat(n)=blath(i)
          ba=baath(i)
          baat(n)=ba
          csbaat(n)=cos(ba)
          snbaat(n)=sin(ba)
          to=toath(i)
          toat(n)=to
          cstoat(n)=cos(to)
          sntoat(n)=sin(to)
c ______________________________ bonds to previous & following atoms
          iow=iowath(i)
          if (iow.eq.0) then          ! 1st atom of residue
            if (nrs.eq.ifirs) then  ! 1st atom of 'nml'

              iowat(n)=ibdath(1,i)+ntlat
              iyowat(n)=iybdath(1,i)

              do j = 1,mxbd-1
                ibdath(j,i)=ibdath(j+1,i)
                iybdath(j,i)=iybdath(j+1,i)
              enddo
              ibdath(mxbd,i)=0
              iybdath(mxbd,i)=1

              nbdath(i)=nbdath(i)-1
            else                    ! connected with prev. res.
              nh=ixatrs(nrs-1)    ! atom to 'next' residue of prev.res.
              iowat(n)=nh
              iyowat(n)=1         !!! only single bonds assumed !!!

c ___________________________ correct atom to 'next' res.
              nbd=nbdat(nh)
              if (nbd.eq.mxbd) then
                write(*,'(a,i2,a,i4,2a,i4,a)') 
     #           ' getmol> need ',(mxbd+2),
     #           'th bond to connect residues ',
     #           nrs-1,seq(nrs-1),' and ',nrs,seq(nrs)
                close(lunlib)
                stop
              else  !  correct atom to 'next' res.
c _______________________________!! dihedrals for atoms bound to 'nh'
c                                   are assumed to be phase angles !!
                do j=1,nbd

                  nj=ibdat(j,nh)
                  t=toat(nj)

                  if (t.eq.0.0) then
                    write (*,'(3a,/,2a)') 
     #               ' getmol> DIHEDRAL for atom ',nmat(nj),
     #               ' should be PHASE angle with respect to atom ',
     #               nmat(n),' & therefore must be not 0.0 !!'
                    close(lunlib)
                    stop
                  endif

                  t=t+to
                  if (abs(t).gt.pi) t=t-sign(pi2,t)
                  toat(nj)=t
                  cstoat(nj)=cos(t)
                  sntoat(nj)=sin(t)

                enddo

                nbd1=nbd+1
                ibdat(nbd1,nh)=n
                iybdat(nbd1,nh)=1  ! (only single bonds !)
                nbdat(nh)=nbd1

              endif

            endif

          else                      ! connected within res.

            iowat(n)=ntlat+iow
            iyowat(n)=iyowath(i)

          endif

          nbdat(n)=nbdath(i)

          do j=1,mxbd

            ibd=ibdath(j,i)

            if (ibd.ne.0) then
              ibdat(j,n)=ibd+ntlat
              iybdat(j,n)=iybdath(j,i)
            else
              ibdat(j,n)=0
              iybdat(j,n)=1
            endif

          enddo

        enddo  ! ... atoms

c ________________________________________________________ Variables
        ivrrs1(nrs)=ntlvr+1
        mvr=0

        do i=1,nvr

          if (nrs.eq.ifirs) then

            iat=iatvrh(i)
c ____________________________________ Exclude all variables for 1st atom
c                                       & torsion for atoms bound to it
            if ( iat.eq.1.or.
     #        (iowath(iat).eq.1.and.ityvrh(i).eq.3)) goto 1

          endif

          mvr=mvr+1
          ntlvr=ntlvr+1
          ityvr(ntlvr)=ityvrh(i)
          iclvr(ntlvr)=iclvrh(i)
          iatvr(ntlvr)=iatvrh(i)+ntlat
          nmvr(ntlvr)=nmvrh(i)

    1   enddo   ! ... Variables

        nvrrs(nrs)=mvr

        iatrs1(nrs)=ntlat+1   ! first backbone atom of res.
        ixatrs(nrs)=ntlat+nxt    ! last backbone atom
        ntlat=ntlat+nat
        iatrs2(nrs)=ntlat        ! last atom of res.

      enddo             ! ... residues

      close(lunlib)

c _______________________________ Variables
      if (nml.eq.1) then
        nvrml(nml)=ntlvr
      else
        nvrml(nml)=ntlvr-ivrml1(nml) + 1
      endif

      return
      end
c **************************************
      subroutine redres(res,nat,nxt,nvrr)

c .......................................................
c PURPOSE:  read atom data for residue 'res' from library
c           (file 'lunlib' 'reslib' opened in routine calling
c            this one)
c
c OUTPUT:   nat   - number of atoms in residue
c           nxt   - atom which may bind to following residue
c           nvrr  - number of variables in residue
c           for atoms     - nmath,blath,baath(rad),toath(rad),
c                           ityath,iyowath,iowath (INSIDE residue,
c                                                  =0 if 1st atom)
c           for variables - ityvrh (1=bl/2=ba/3=to),iclvrh,iatvrh,nmvrh
c
c LIBRARY:  residue-lines:
c            '#', res, nat, nxt;  Format: a1,a4,2i4
c           atom-lines:
c           nmat,3{"fix" =' ', clvr,nmvr, blat/baat(deg)/toat(deg)},
c             cgat, ityat, iowat,ibdat1,ibdat2,ibdat3;
c           Format: a4, 3(1x,i2,a1,a3,f9.3), f7.4, i4,4i4
c
C CALLS: iendst,tolost
c
c .......................................................

      include 'INCL.H'

      dimension icl(3),ibd(mxbd)
      character blnk,fix(3),nm(3)*3,res*4,resl*4,line*132
      data blnk/' '/


      nln=0
      do i=1,mxbd
        ibd(i)=0
      enddo

      resl=res

      call tolost(resl)  ! ensure lower case for residue name

c ________________________________ find residue 'resl'
    1 line=blnk
      nln=nln+1

      read (lunlib,'(a)',end=2,err=3) line
      lg=iendst(line)

      if (lg.ge.13.and.line(1:1).eq.'#'.and.line(2:5).eq.resl) then
c _____________________________________________ read atom data for 'resl'
        read (line(6:13),'(2i4)',err=3) nat,nxt

        if (nat.gt.mxath) then
          write (*,'(a,i5)') ' redres> number of atoms > ',mxath
          close(lunlib)
          stop
        endif

        nvrr=0
        do i=1,nat

          nln=nln+1

          read (lunlib,'(a4,3(1x,i2,a1,a3,d9.3),d7.4,i4,4i4)',
     #                 end=3,err=3)
     #     nmath(i),icl(1),fix(1),nm(1),blath(i),icl(2),fix(2),nm(2),ba,
     #     icl(3),fix(3),nm(3),to,cgath(i),ity,iow,(ibd(j),j=1,mxbd)

          if (ity.le.0.or.ity.gt.mxtyat) goto 6
          ityath(i)=ity

          jow=abs(iow)

          if (res(:3).eq.'ace'.and.i.eq.1) then  ! exception from following check
            iexcp = 1
          else
            iexcp = 0
          endif

          if (iexcp.eq.0.and.i.le.jow) then
            if (i.eq.jow) then
              write (*,'(5a)') ' redres> atom ',nmath(i),' of ',
     #                          resl,' cannot preceed itself '
            else
              write (*,'(5a,i4)') ' redres> atom ',nmath(i),' of ',
     #                    resl,' should be placed AFTER atom #',jow
            endif
            goto 5
          endif

          iowath(i)=jow
          iyowath(i)=sign(1,iow)
c ____________________________________ check order & find number of bonds
c                                      (bonds closing ring must be last !)
          ib1=abs(ibd(1))
          ib2=abs(ibd(2))
          ib3=abs(ibd(3))

          if (ib1.eq.i.or.ib2.eq.i.or.ib3.eq.i) goto 4

          if (ib1.eq.0) then                  ! no bond to following
            if (ib2.ne.0.or.ib3.ne.0) goto 4
            nbdath(i)=0
          else
            if (ib1.eq.jow) goto 4
            if (ib2.eq.0) then
              if (ib3.ne.0) goto 4
              nbdath(i)=1
            else
              if ( ib2.eq.jow.or.ib2.eq.ib1.or.
     #            (ib2.gt.i.and.ib2.lt.ib1) ) goto 4
              if (ib3.eq.0) then
                nbdath(i)=2
              else
                if (ib3.eq.jow.or.ib3.eq.ib1.or.ib3.eq.ib2.or.
     #              (ib3.gt.i.and.(ib3.lt.ib1.or.ib3.lt.ib2)) ) goto 4
                nbdath(i)=3
              endif
            endif
          endif

          do j=1,mxbd
            ibdath(j,i)=abs(ibd(j))
            iybdath(j,i)=sign(1,ibd(j))
          enddo

          baath(i)=ba*cdr  ! convert angles into 'radians'
          toath(i)=to*cdr

c ______________________________ internal degrees of freedom
          do j=1,3
            if (fix(j).ne.blnk) then
              nvrr=nvrr+1

              if (nvrr.gt.mxvrh) then
                write (*,'(a,i5)') ' redres> number of variables > ',
     #                             mxvrh
                close(lunlib)
                stop
              endif

              ic=icl(j)

              if ( ic.le.0    
     #         .or.(j.eq.3.and.ic.gt.mxtyto)           ! dihedral
     #         .or.(j.eq.2.and.ic.gt.mxtyba)           ! bond angle
     #         .or.(j.eq.1.and.ic.gt.mxtybl) ) goto 7  ! b. length

              ityvrh(nvrr)=j
              iclvrh(nvrr)=ic
              iatvrh(nvrr)=i
              nmvrh(nvrr)=nm(j)

            endif

          enddo
        enddo     ! ... atoms

        return
      endif

      goto 1

c ____________________________________________________________ ERRORS
    2 write (*,'(4a)') ' redres> residue >',resl,'< NOT FOUND in ',
     #reslib(1:iendst(reslib))
      close(lunlib)
      stop

    3 write (*,'(a,i4,2a)') ' redres> ERROR reading line No. ',nln,
     #' in ',reslib(1:iendst(reslib))
      close(lunlib)
      stop

    4 write (*,'(4a)') ' redres> Incorrect order of bonds for atom ',
     #                      nmath(i),' of ',resl

    5 write (*,'(8x,2a)') '... must correct ',
     #                      reslib(1:iendst(reslib))
      close(lunlib)
      stop

    6 write (*,'(a,i2,4a)') ' redres> unknown type :',ity,
     #                   ': for atom ',nmath(i),' in residue ',resl
      close(lunlib)
      stop

    7 write (*,'(a,i2,4a)') ' redres> unknown class :',ic,
     #                   ': for variable ',nm(j),' in residue ',resl
      close(lunlib)
      stop

      end


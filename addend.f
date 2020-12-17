c **************************************************************
c
c
c This file contains the subroutines:  addend, redchg, rplgrp
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c $Id: addend.f 334 2007-08-07 09:23:59Z meinke $
c **************************************************************
      subroutine addend(nml,grpn,grpc)

c ..............................................................
c  PURPOSE:  modify terminal residues to complete bonding scheme
c            with residue 'grpn' at N- and residue 'grpc' at C-terminus
c            ! need initial co-ordinates for residues to modify
c            ! for N-terminus: may add only simple groups
c
c  CALLS:  rplgrp,tolost,redchg
c ..............................................................
 
      include 'INCL.H'

      character grpn*4,grpc*4

      character res*4,rpat*4,sbrs*4,grn*4,grc*4


      grn = grpn
      call tolost(grn)
      grc = grpc
      call tolost(grc)

      if (grn(:3).eq.'ace'.or.grc(:3).eq.'ace'
     #.or.grn(:3).eq.'nme'.or.grc(:3).eq.'nme') then

        write(*,'(2a)') ' addend> N-Acetyl (ace) or N-Methylamide (nme)'
     #   ,' should be put in SEQUENCE file, not added as end groups'

        stop
      endif

c __________________________________________ N-terminus
      ifirs=irsml1(nml)
      rpat='n   '
      res=seq(ifirs)
      call tolost(res)
      if (res(:3).ne.'ace') then
        if (grn(:3).eq.'nh2') then
          if (res(:3).eq.'hyp'.and..not.flex) then
            if(sh2) then
              sbrs='nh2+'
            else
              write (*,'(2a)') ' addend> ',
     #         ' No N-terminal Hyp possible with ECEPP/3 dataset'
              stop
            endif

          elseif (res(:3).eq.'pro') then
            sbrs='nh1 '
          elseif (res(:3).eq.'gly'.and..not.flex) then
            if (sh2) then
              sbrs='nh2 '
            else
              sbrs='nh2g'
            endif
          else
            sbrs='nh2 '
          endif

        elseif (grn.eq.'nh3+')  then
          if (res(:3).eq.'pro'.and..not.flex) then
            sbrs='nh2+'
          elseif (res(:3).eq.'gly'.and..not.flex) then
            sbrs='nh3g'
          else
            sbrs='nh3+'
          endif

        else

          write(*,'(2a)') ' addend> Can add only ',
     *     'nh2 or nh3+ to N-terminus'
          stop

        endif

        call rplgrp(nml,ifirs,rpat,sbrs)
        if (flex) call redchg(nml,ifirs,rpat,sbrs) ! Flex dataset

      else ! ace

        write(*,'(2a)') ' addend> Acetyl group',
     #     ' at N-terminus not modified'
      endif

c __________________________________________ C-terminus
      ilars=irsml2(nml)
      rpat='c   '
      res=seq(ilars)
      call tolost(res)
      
      if (res(:3).ne.'nme') then

        if (grc.eq.'cooh') then

          if (res(:3).eq.'gly'.and..not.sh2) then
            sbrs='gooh'
          else
            sbrs='cooh'
          endif

        elseif (grc.eq.'coo-') then

          if (res(:3).eq.'gly'.and..not.sh2) then
            sbrs='goo-'
          else
            sbrs='coo-'
          endif

        else

          write(*,'(2a)') ' addend> Can add only ',
     #     'cooh or coo- to C-terminus'
          stop

        endif
       
        call rplgrp(nml,ilars,rpat,sbrs)
        if (flex) call redchg(nml,ilars,rpat,sbrs) ! Flex dataset

      else  ! N'-methylamide

        write(*,'(2a)') ' addend> N-Methylamide',
     #     ' at C-terminus not modified'

      endif

c ----------------------------- net charge of molecule
      cg = 0.d0
      do i=iatrs1(irsml1(nml)),iatrs2(irsml2(nml))
        cg = cg + cgat(i)
      enddo
      if (abs(cg).gt.1.d-5) write(*,'(a,i2,a,f7.3,/)')
     #        ' addend> Net charge of molecule #'
     #        ,nml,': ',cg

      return
      end
c ****************************************
      subroutine rplgrp(nml,nrs,rpat,sbrs)

c ...............................................................
c PURPOSE: replace atom(s) rooted at atom 'rpat' in residue
c          'nrs' of molecule 'nml' by atom(s) rooted at 
c          'rpat' of residue 'sbrs' (same name of root
c          atom 'rpat' maintains bonding geometry for
c          preceeding atoms in 'nrs')
c
c          is NOT performed if 'rpat' is within mainchain,
c          except it is first/last mainchain atom of 'nml'
c
c CALLS: dihedr,iopfil,iendst,eyring,fndbrn,redres,setsys,valang
c ...............................................................

      include 'INCL.H'

      character rpat*4,sbrs*4
      logical ntbb,bb

      dimension ibd(mxbd+1),iybd(mxbd+1)


      ifirs=irsml1(nml)
      ilars=irsml2(nml)

      nxt=ixatrs(nrs)
      nfi=iatrs1(nrs)
      nla=iatrs2(nrs)
c __________________________ indices of atoms to be replaced
      do i=nfi,nla
        if (rpat.eq.nmat(i)) then
          nfirp=i
          goto 1
        endif
      enddo
      write (*,'(4a,i4,a,i4)') ' rplgrp> cannot find atom >',rpat,
     #'< to be replaced in residue ',seq(nrs),nrs,' of molecule ',nml
      stop

    1 call fndbrn(nml,nrs,nfirp,nlarp,irng1,irng2,bb)
      if (irng1.gt.0) goto 10
      ntbb=.false.

      if (bb) then  ! ..... backbone

        if (nfirp.eq.nxt.and.nrs.eq.ilars) goto 2 !  last mainchain atom

        if (nfirp.eq.nfi.and.nrs.eq.ifirs) then   !  1st MAINCHAIN ATOM !

          ntbb=.true.

          ibd(1)=iowat(nfirp)
          iybd(1)=iyowat(nfirp)

          do i=1,mxbd
            ibd(i+1)=ibdat(i,nfirp)
            iybd(i+1)=iybdat(i,nfirp)
          enddo

          ibdrg=0         ! __________________ check ring
          iybdrg=1

          do i=1,nbdat(nfirp)+1
            if (iowat(ibd(i)).ne.nfirp) then
              if (ibdrg.ne.0) then
                write (*,'(2a,i3)') 
     #             ' rplgrp> Can handle only simple ring at 1st',
     #             ' atom of molecule #',nml
                stop
              endif
              ibdrg=ibd(i)
              iybdrg=iybd(i)
            endif
          enddo
          nxtbb1=nlarp+1   ! _________________ next backbone atoms
          isgbb1=iyowat(nxtbb1)

          if (nxtbb1.eq.nxt) then
            if (nrs.lt.ilars) then
              nxtbb2=iatrs1(nrs+1)
              goto 3
            endif
          else
            do i=nxt,nxtbb1+1,-1
              if (iowat(i).eq.nxtbb1) then
                nxtbb2=i
                goto 3
              endif
            enddo
          endif
          goto 11
        else
          write (*,'(4a,i4,a,i4)') 
     #      ' rplgrp> Cannot replace BACKBONE atom ',rpat,
     #      ' of residue ',seq(nrs),nrs,' in molecule #',nml
          stop
        endif

      endif  ! N-terminus
c _________________________________ previous atoms
    2 if (nfirp.eq.nfi.and.nrs.eq.ifirs) goto 11
      nxtbb1=iowat(nfirp)
      if (nxtbb1.eq.nfi.and.nrs.eq.ifirs) goto 11
      nxtbb2=iowat(nxtbb1)
c _______________________________ get data for substituent atoms
    3 if (iopfil(lunlib,reslib,'old','formatted').le.izero) then
        write (*,'(a,/,a,i3,2a)') 
     #    ' rplgrp> ERROR opening library of residues:',
     #    ' LUN=',lunlib,' FILE=',reslib(1:iendst(reslib))
        stop
      endif
      call redres(sbrs,natsb,nxtsb,nvrsb)
      close (lunlib)
c __________________________ indices of substituent atoms
      do i=1,natsb
        if (rpat.eq.nmath(i)) then
          nfisb=i
          goto 4
        endif
      enddo
      write (*,'(4a)') ' rplgrp> Cannot find atom >',rpat,
     #'< in substituent residue ',sbrs
      stop

    4 nlasb=nfisb
      do i=1,nbdath(nfisb)
        ib=ibdath(i,nfisb)
        if (iowath(ib).lt.nfisb) goto 10
        do j=ib,natsb
          if (j.gt.ib.and.iowath(j).lt.ib) goto 5
          do k=1,nbdath(j)
            if (ibdath(k,j).lt.nfisb) goto 10
          enddo
          nlasb=j
        enddo  ! ... branch atoms
    5 enddo  ! ... branches
c _________________________________________________ local axes at 'nfirp'
      call setsys(nxtbb1,nfirp,nxtbb2,x1,x2,x3,y1,y2,y3,z1,z2,z3)

      xtoat(nfirp)=x1
      ytoat(nfirp)=x2
      ztoat(nfirp)=x3
      xbaat(nfirp)=z1
      ybaat(nfirp)=z2
      zbaat(nfirp)=z3

c _____________________ add virtual atoms
      if (ntbb) then

        ct=cstoat(nxtbb2)  ! t.angle_(+2)
        st=sntoat(nxtbb2)
        ca=csbaat(nxtbb1)  ! b.angle_(+1)
        sa=snbaat(nxtbb1)

c ------------------- Eyring
        h2=-sa*ct
        h3=-sa*st
        x1=-ca*x1+h2*y1+h3*z1
        x2=-ca*x2+h2*y2+h3*z2
        x3=-ca*x3+h2*y3+h3*z3
        dx=one/sqrt(x1*x1+x2*x2+x3*x3)
        x1=x1*dx
        x2=x2*dx
        x3=x3*dx

        xat(izero)=xat(nfirp)+x1
        yat(izero)=yat(nfirp)+x2
        zat(izero)=zat(nfirp)+x3
        z1=-st*y1+ct*z1
        z2=-st*y2+ct*z2
        z3=-st*y3+ct*z3
        dz=one/sqrt(z1*z1+z2*z2+z3*z3)
        z1=z1*dz
        z2=z2*dz
        z3=z3*dz

        ct=cstoat(nxtbb1)  ! t.angle_(+1)
        st=sntoat(nxtbb1)

c -------------------- Eyring with b.angle = 90 deg.
        xat(-ione)=xat(izero)-ct*(z2*x3-z3*x2)-st*z1
        yat(-ione)=yat(izero)-ct*(z3*x1-z1*x3)-st*z2
        zat(-ione)=zat(izero)-ct*(z1*x2-z2*x1)-st*z3

      endif
c _____________________________________________ Shift atom data
      nrp=nlarp-nfirp
      nsb=nlasb-nfisb
      if (nrp.ne.nsb) then
        nsh=nsb-nrp
        do i=nfi,nfirp-1    ! bonds to atoms after repl. group
          do j=1,mxbd
            ib=ibdat(j,i)
            if (ib.gt.nlarp.and.ib.le.nla) ibdat(j,i)=ib+nsh
          enddo
        enddo
        ilaat=iatrs2(irsml2(ntlml))
        if (nrp.gt.nsb) then   ! less atoms
          i1=nlarp+1
          i2=ilaat
          i3=1
        else  ! more atoms
          if ((ilaat+nsh).gt.mxat) then
            write (*,'(a,i5)') ' rplgrp> number of atoms > ',mxat
            stop
          endif

          i1=ilaat
          i2=nlarp+1
          i3=-1
        endif

        do i=i1,i2,i3
          ii=i+nsh
          nbdat(ii)=nbdat(i)

          ibd(1)=iowat(i)
          iybd(1)=iyowat(i)
          do j = 1,mxbd
            ibd(j+1)=ibdat(j,i)
            iybd(j+1)=iybdat(j,i)
          enddo

          do j=1,mxbd+1
            if (ibd(j).gt.nfirp) ibd(j)=ibd(j)+nsh
          enddo

          iowat(ii)=ibd(1)
          iyowat(ii)=iybd(1)
          do j = 1,mxbd
            ibdat(j,ii)=ibd(j+1)
            iybdat(j,ii)=iybd(j+1)
          enddo

          ityat(ii)=ityat(i)
          nmat(ii)=nmat(i)
          cgat(ii)=cgat(i)
          blat(ii)=blat(i)
          baat(ii)=baat(i)
          snbaat(ii)=snbaat(i)
          csbaat(ii)=csbaat(i)
          toat(ii)=toat(i)
          sntoat(ii)=sntoat(i)
          cstoat(ii)=cstoat(i)
          xat(ii)=xat(i)
          yat(ii)=yat(i)
          zat(ii)=zat(i)
          xtoat(ii)=xtoat(i)
          ytoat(ii)=ytoat(i)
          ztoat(ii)=ztoat(i)
          xbaat(ii)=xbaat(i)
          ybaat(ii)=ybaat(i)
          zbaat(ii)=zbaat(i)

        enddo
c ____________________________________________ Shift residue data
        do i=nrs+1,irsml2(ntlml)
          iatrs1(i)=iatrs1(i)+nsh
          iatrs2(i)=iatrs2(i)+nsh
          ixatrs(i)=ixatrs(i)+nsh
        enddo
        iatrs2(nrs)=nla+nsh
        if (nxt.gt.nlarp) ixatrs(nrs)=nxt+nsh
      else
        nsh=0
      endif
c _________________________________________ Correct data of 'nfirp'
      ish=nfirp-nfisb
      ityat(nfirp)=ityath(nfisb)

      if(.not.sh2) cgat(nfirp)=cgath(nfisb)  ! NOT for 'ECEPP/2'

      nb=nbdath(nfisb)           ! _______________ Bonds

      do i=1,mxbd
        ib=ibdath(i,nfisb)
        if (ib.ne.0) then
          ibd(i)=ib+ish
          iybd(i)=iybdath(i,nfisb)
        else
          ibd(i)=0
          iybd(i)=1
        endif

      enddo

      if (ntbb) then

        ibd(mxbd+1)=0
        iybd(mxbd+1)=1

        ibd(nb+1)=nxtbb1+nsh   ! bond to next backbone atom
        iybd(nb+1)=isgbb1

        if (ibdrg.ne.0) then
          nb=nb+1
          if (nb.gt.mxbd) then
            write (*,'(6a,/,2a,3(i4,a))')
     #      ' rplgrp> Cannot add atoms following ',rpat,
     #      ' from group ',sbrs,' to atom ',rpat,
     #      ' of residue ',seq(nrs),nrs,' in molecule #',nml,
     #      ' because need >',(mxbd+1),' bonds'
            stop
          endif
          ibd(nb+1)=ibdrg+nsh
          iybd(nb+1)=iybdrg
        endif
        iowat(nfirp)=ibd(1)
        iyowat(nfirp)=iybd(1)
        do i=1,mxbd
          ibdat(i,nfirp)=ibd(i+1)
          iybdat(i,nfirp)=iybd(i+1)
        enddo
      else  ! not 'ntbb'

        do i=1,mxbd
          ibdat(i,nfirp)=ibd(i)
          iybdat(i,nfirp)=iybd(i)
        enddo

      endif
      nbdat(nfirp)=nb
c _________________________________________ Add data for substituent
      ii=nfirp
      do i=nfisb+1,nlasb
        ii=ii+1
        nbdat(ii)=nbdath(i)

        iow=iowath(i)+ish
        iowat(ii)=iow
        iyowat(ii)=iyowath(i)

        do j=1,mxbd
          ib=ibdath(j,i)
          if (ib.ge.nfisb) then
            ibdat(j,ii)=ib+ish
          else
            ibdat(j,ii)=ib
          endif
          iybdat(j,ii)=iybdath(j,i)
        enddo

        ityat(ii)=ityath(i)
        nmat(ii)=nmath(i)
        cgat(ii)=cgath(i)
        blat(ii)=blath(i)
        ba=baath(i)
        baat(ii)=ba
        csbaat(ii)=cos(ba)
        snbaat(ii)=sin(ba)
        to=toath(i)
        toat(ii)=to
        cstoat(ii)=cos(to)
        sntoat(ii)=sin(to)
        call eyring(ii,iow)
        if (ntbb) then        ! reset some internal coordinates
          if (iow.eq.nfirp) then
            ba=valang(izero,nfirp,ii)
            baat(ii)=ba
            csbaat(ii)=cos(ba)
            snbaat(ii)=sin(ba)
            to=dihedr(-ione,izero,nfirp,ii)
            toat(ii)=to
            cstoat(ii)=cos(to)
            sntoat(ii)=sin(to)
          elseif (iowat(iow).eq.nfirp) then
            to=dihedr(izero,nfirp,iow,ii)
            toat(ii)=to
            cstoat(ii)=cos(to)
            sntoat(ii)=sin(to)
          endif

        endif  ! ntbb

      enddo  ! substituent atoms
c ___________________________________________________ Take care of Variables
c (assume variables of replaced group/substituent to be stored CONSECUTIVELY)

      ilavr=ivrml1(ntlml)+nvrml(ntlml)-1
      ifivr=ivrrs1(nrs)   ! variables to be replaced (#=ivrrp,last=lvrrp)
      ivrrp=0
      do i=ifivr,ifivr+nvrrs(nrs)-1
        iat=iatvr(i)
        if (nfirp.lt.iat.and.iat.le.nlarp) then
          ivrrp=ivrrp+1
          lvrrp=i
        endif
      enddo
      if (ivrrp.eq.0) then   ! No variables to replace
        do i=ifivr,ilavr
          if (iatvr(i).gt.nlarp) then
            lvrrp=i-1
            goto 6
          endif
        enddo
        lvrrp=ilavr
      endif
    6 ivrsb=0   ! variables from substituent (#=ivrsb)
      do i=1,nvrsb
        iat=iatvrh(i)
        if (nfisb.lt.iat.and.iat.le.nlasb) then
          ity=ityvrh(i)
          if (ntbb.and.iowath(iat).eq.nfisb.and.ity.gt.2) goto 7
          ivrsb=ivrsb+1
          nmvrh(ivrsb)=nmvrh(i)
          ityvrh(ivrsb)=ityvrh(i)
          iclvrh(ivrsb)=iclvrh(i)
          iatvrh(ivrsb)=iat
        endif
    7 enddo
      if (nsh.ne.0) then   ! if # of atoms changed
        do i=lvrrp+1,ilavr
          iatvr(i)=iatvr(i)+nsh
        enddo
      endif

      if (ivrrp.ne.ivrsb) then   ! shift data for variables
        jsh=ivrsb-ivrrp
        if (ivrrp.gt.ivrsb) then
          i1=lvrrp+1
          i2=ilavr
          i3=1
        else
          if ((ilavr+jsh).gt.mxvr) then
            write (*,'(a,i5)') ' rplgrp> number of variables > ',mxvr
            stop
          endif
          i1=ilavr
          i2=lvrrp+1
          i3=-1
        endif
        do i=i1,i2,i3
          ii=i+jsh
          ityvr(ii)=ityvr(i)
          iclvr(ii)=iclvr(i)
          nmvr(ii)=nmvr(i)
          iatvr(ii)=iatvr(i)
        enddo

        do i=nrs+1,irsml2(ntlml)
          ivrrs1(i)=ivrrs1(i)+jsh
        enddo
        nvrrs(nrs)=nvrrs(nrs)+jsh
        nvrml(nml)=nvrml(nml)+jsh
        do i=nml+1,ntlml
          ivrml1(i)=ivrml1(i)+jsh
        enddo
      endif
      ii=lvrrp-ivrrp     ! Add variables for substitutent
      do i=1,ivrsb
        ii=ii+1
        nmvr(ii)=nmvrh(i)
        ityvr(ii)=ityvrh(i)
        iclvr(ii)=iclvrh(i)
        iatvr(ii)=iatvrh(i)+ish
      enddo

      return
c __________________________________________ Errors
   10 write (*,'(3a,/,2a,i4,a,i4,/,2a)') 
     #   ' rplgrp> Cannot replace atom(s) following ',rpat,
     #   ' from INSIDE a ring','    in residue: ',seq(nrs),nrs,
     #   ' in molecule #',nml,' or in substitute: ',sbrs
      stop
   11 write (*,'(4a,i4,a,i4,/,a)') 
     #   ' rplgrp> Cannot replace atom(s) following ',rpat,
     #   ' of residue ',seq(nrs),nrs,' in molecule #',nml,
     #   ' since necessary 2 previous atoms are not available' 
      stop

      end
c ****************************************
      subroutine redchg(nml,nrs,rpat,sbrs)

c .........................................................
c PURPOSE: read and place atomic point charges from residue
c          'sbrs' to residue 'nrs' of molecule 'nml'
c          from library 'chglib' with LUN=lunchg, if ilib=1
c                       'reslib' with LUN=lunlib, if ilib=2
c
c CALLS: iopfil,iendst,tolost
c ........................................................

      include 'INCL.H'

      character rpat*4,sbrs*4,atnm*4,res*4,cgty*5,line*100

      ifirs=irsml1(nml)
      ilars=irsml2(nml)
      res=seq(nrs)
      call tolost(res)

      if (nrs.eq.ifirs.or.nrs.eq.ilars) then

        if (rpat.eq.'n   '.or.rpat.eq.'c   ') then
          if (nrs.eq.ifirs.and.nrs.eq.ilars) goto 10 ! Dont have this yet

          if (rpat.eq.'n   '.and.nrs.eq.ifirs) then
            if (res(1:3).eq.'pro'.or.res(1:3).eq.'hyp') then
              if (sbrs.eq.'nh1 ') cgty='n'//res
              if (sbrs.eq.'nh2+') cgty='+'//res
            else
              if (sbrs.eq.'nh2 ') cgty='n'//res
              if (sbrs.eq.'nh3+') cgty='+'//res
            endif
          elseif (rpat.eq.'c   '.and.nrs.eq.ilars) then
            if (sbrs.eq.'cooh') cgty='c'//res
            if (sbrs.eq.'coo-') cgty='-'//res
          else
            goto 10
          endif
        else
          write (*,*) ' redchg> dont know which end goup is present'
          stop
        endif
        ilib=1
      else
        cgty(1:)=sbrs
        ilib=2
      endif

      if (ilib.eq.1) then

        if (iopfil(lunchg,chgfil,'old','formatted').le.izero) then
          write (*,'(a,/,a,i3,2a)') 
     #      ' redchg> ERROR opening library of charges:',
     #      ' LUN=',lunchg,' FILE=',chgfil(1:iendst(chgfil))
          stop
        endif

    1   line=' '
        read (lunchg,'(a)',end=3) line
        l=iendst(line)

        if (l.ge.10.and.line(1:1).eq.'#'.and.line(2:6).eq.cgty) then
          read (line(7:10),'(i4)') nchg
          if (nchg.le.mxath) then
            read (lunchg,'(10(2x,a4,1x))') (nmath(i),i=1,nchg)
            read (lunchg,'(10(f6.4,1x))') (cgath(i),i=1,nchg)
            close(lunchg)
            do i=iatrs1(nrs),iatrs2(nrs)
              atnm=nmat(i)
              do j=1,nchg
                if (nmath(j).eq.atnm) then
                  cgat(i)=cgath(j)
                  goto 2
                endif
              enddo
              write (*,'(6a)') ' redchg> Cannot find atom: ',atnm,
     #                        ' for entry: ',cgty,' in library: ',
     #                        chgfil(1:iendst(chgfil))
              stop
    2       enddo
            return
          else
            write (*,'(4a)')
     #       ' redchg> must increase MXATH to read data for entry: ',
     #        cgty,' in library: ',chgfil(1:iendst(chgfil))
            close(lunchg)
            stop
          endif
        endif
        goto 1
    3   write (*,'(4a)')
     #   ' redchg> Cannot find entry: ',cgty,' in library: ',
     #      chgfil(1:iendst(chgfil))
        close(lunchg)
        stop

      elseif (ilib.eq.2) then

        if (iopfil(lunlib,reslib,'old','formatted').le.izero) then
          write (*,'(a,/,a,i3,2a)') 
     #      ' redchg> ERROR opening library of residues:',
     #      ' LUN=',lunlib,' FILE=',reslib(1:iendst(reslib))
          stop
        endif

    4   line=' '
        read (lunlib,'(a)',end=6) line
        l=iendst(line)

        if (l.ge.9.and.line(1:1).eq.'#'.and.line(2:5).eq.cgty(1:4)) then
          read (line(6:9),'(i4)') nchg
          if (nchg.le.mxath) then
            read (lunlib,'(a4,42x,d7.4)') (nmath(i),cgath(i),i=1,nchg)
            close(lunlib)
            do i=iatrs1(nrs),iatrs2(nrs)
              atnm=nmat(i)
              do j=1,nchg
                if (nmath(j).eq.atnm) then
                  cgat(i)=cgath(j)
                  goto 5
                endif
              enddo
              write (*,'(6a)') ' redchg> Cannot find atom: ',atnm,
     #                        ' for entry: ',cgty,' in library: ',
     #                        reslib(1:iendst(reslib))
              stop
    5       enddo
            return
          else
            write (*,'(4a)')
     #       ' redchg> must increase MXATH to read data for entry: ',
     #        cgty,' in library: ',reslib(1:iendst(reslib))
            close(lunchg)
            stop
          endif
        endif
        goto 4
    6   write (*,'(4a)')
     #   ' redchg> Cannot find entry: ',cgty,' in library: ',
     #      reslib(1:iendst(reslib))
        close(lunchg)
        stop

      endif

   10 write (*,'(4a)')
     #    ' redchg> Do not have charges for N/C-terminal residue ',
     #    res,' modified with group :',sbrs
      stop

      end


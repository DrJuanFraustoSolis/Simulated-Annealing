c**************************************************************
c
c This file contains the subroutines: setmvs,fndbrn
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************


      subroutine setmvs(nml)

c ......................................................................
c PURPOSE: 1. ORDER variables according to rules:
c             variables with same base: 1st comes TORSION (can be only
c               one with this base, since PHASE a. assumed to be FIXED),
c               after this, for atoms branching from this base:
c               for a b.angle & b.length with common primary moving
c               atom=branch atom - b.angle comes 1st
c
c             iorvr(i), i=i_fivr_ml,i_lavr_ml -> indices of ordered var.
c
c          2. define NON-OVERLAPPING moving sets of atoms in molecule
c             'nml' related to local variables
c
c             nmsml(i_ml) - number of moving sets per molecule
c             imsvr1(i_vr),imsvr2() - indices of 1st/last m.s for var. 'i_vr'
c                                     in 'latms1' & 'latms2'
c             latms1(i_ms),latms2() - range of atoms of i-th m.s
c
c          3. define indices of next-following variables for each var.,
c             which complete its physical moving set ('added' variables)
c          
c             nadml(i_ml) - number of 'added' var.s per molecule
c             iadvr1(i_vr),iadvr2() - indices of 1st/last 'added' var. for
c                                     var. 'i_vr' in 'ladvr'
c             ladvr() - indices of 'added' variables
c
c          4. define index of corresponding variable for each atom
c
c ! routine must be called successively for molecules 1 -> ntlml
c
c CALLS: fndbrn, nursvr
c ......................................................................

      include 'INCL.H'

      logical bb

      parameter (mxh=10)
      dimension lvw1h(mxh),lvw2h(mxh),l1h(mxh),l2h(mxh)


      ntlvr=nvrml(nml)

      if (nml.eq.1) then
        imsml1(1)=1
        nms=0
        iadml1(1)=1
        nad=0
      else
        imsml1(nml)=imsml1(nml-1)+nmsml(nml-1)
        nms=imsml1(nml)-1
        iadml1(nml)=iadml1(nml-1)+nadml(nml-1)
        nad=iadml1(nml)-1
      endif

      if (ntlvr.eq.0) then
        write (*,'(a,i4)')
     #           ' setmvs> No variables defined in molecule #',nml
        nmsml(nml)=0
        nadml(nml)=0
        return
      endif
c _________________ Take index of primary atom for each variable 
c                   (i.e. index of atom moved by variable) to
c                   sort variables, handling variables with same base:
c                   modify indices to obtain appropriate order

      ifirs=irsml1(nml)
      ilars=irsml2(nml)

      ifivr=ivrml1(nml)
      ilavr=ifivr+ntlvr-1
      ifiat=iatrs1(ifirs)
      ilaat=iatrs2(ilars)

      do n=ifirs,ilars   ! ______________________ Residues
        ib=ivrrs1(n)
        do i=ib,ib+nvrrs(n)-1  ! _________________ Variables
          ia=iatvr(i)
          io=iowat(ia)     ! ('ia' cannot be 1st atom of 'nml')
          it=ityvr(i)
          if (it.eq.3) then       ! torsion
            do j=1,nbdat(io)
              ii=ibdat(j,io)
              if (iowat(ii).eq.io) ia=min(ia,ii)
            enddo
            iadvr1(i)=ia*10
          elseif (it.eq.2) then   ! bond angle
            iadvr1(i)=ia*10+1
          elseif (it.eq.1) then   ! bond length
            iadvr1(i)=ia*10+2
          endif
          iorvr(i)=i  ! (initialize for sorting)
        enddo  ! ... Variables
      enddo  ! ... Residues
c ___________________________________ Sort variables in ascending order
c                        (i.e. from start of molecule/base of branches)
c array 'iorvr' gives indices of (1st,2nd, ... ,n-th) variables;
c as can be found in arrays for variables (example: ityvr(iorvr())
      k=ilavr
      l=ifivr+ntlvr/2
      ii=ifivr-1
    1 if (l.gt.ifivr) then
        l=l-1
        io=iorvr(l)
        n=iadvr1(io)
      else
        io=iorvr(k)
        n=iadvr1(io)
        iorvr(k)=iorvr(ifivr)
        k=k-1
        if (k.eq.ifivr) then
          iorvr(k)=io
          goto 2
        endif
      endif
      i=l
      j=l+l-ii
      do while (j.le.k)
        if (j.lt.k.and.iadvr1(iorvr(j)).lt.iadvr1(iorvr(j+1))) j=j+1
        if (n.lt.iadvr1(iorvr(j))) then
          iorvr(i)=iorvr(j)
          i=j
          j=j+j-ii
        else
          j=k+1
        endif
      enddo
      iorvr(i)=io
      goto 1
c ______________________________ Find non-overlapping ranges of atoms (moving
c                                sets) for each variable
   2  nms=imsml1(nml)-1

      do io=ifivr,ilavr  ! _____ Loop over variables in 'ascendent' order
        iv=iorvr(io)
        ir=nursvr(iv)     ! residue for variable 'iv'
        ia=iatvr(iv)      ! primary mov. atom
        ib=iowat(ia)      ! base
c __________________________ First, determine complete mov. set for 'iv'
        it=ityvr(iv)
        if (it.eq.3) then       ! torsion
          i1=0
          do i=1,nbdat(ib)
            j=ibdat(i,ib)
            if (iowat(j).eq.ib) then  ! excl. ring
              call fndbrn(nml,ir,j,k,irg1,irg2,bb)
              if (bb) k=ilaat
              if (i1.ne.0) then  ! combine ranges
                if (j.gt.(i2+1).or.k.lt.(i1-1)) then
                  write (*,'(3a,/,2a,i4,a,i3)') 
     #             ' setmvs> Cannot combine disjunct ranges of atom',
     #             ' indices for torsion ',nmvr(iv),' in residue ',
     #             seq(ir),ir,' of molecule # ',nml
                  stop
                else
                  if (j.lt.i1) i1=j
                  if (k.gt.i2) i2=k
                endif
              else
                i1=j
                i2=k
              endif
            endif
          enddo
        elseif (it.eq.2.or.it.eq.1) then   ! b. angle, b. length
          i1=ia
          call fndbrn(nml,ir,i1,i2,irg1,irg2,bb)
          if (bb) i2=ilaat
        endif

        if ((nms+1).gt.mxms) then
          write (*,'(a,i4,a,i5)') ' setmvs> Molecule # ',nml,
     #    ': Number of moving sets > ',mxms
          stop
        endif

        imsvr1(iv)=nms+1  ! index of 1st
        imsvr2(iv)=nms+1  !        & last m.s for var. 'iv'

c ______________ Next, exclude overlaps between mov. set for 'iv' and the
c                m.s. for 'previous' variables by reducing/splitting those

        do jo=ifivr,io-1  ! prev. variables ...
          jv=iorvr(jo)

          j1s=imsvr1(jv)               ! index of 1st m.s. for 'jv'
          jns=imsvr2(jv)-j1s+1         ! # of m.s. for 'jv'

          j=j1s
          do while (j.lt.(j1s+jns))  ! while there are m.s. for 'jv'
            j1=latms1(j)  ! 1st &
            j2=latms2(j)  ! last atom of m.s. 'j'
            if (i1.le.j2.and.i2.ge.j1) then  ! Overlap
              ja=0
              if (i1.gt.j1) then
                if (i2.gt.j2) goto 6
                ja=1
                latms2(j)=i1-1
              endif
              if (i2.lt.j2) then
                if (i1.lt.j1) goto 6
                if (ja.gt.0) then  ! +1 moving set
                  nms=nms+1
                  if (nms.gt.mxms) then
                    write (*,'(a,i4,a,i5)') ' setmvs> Molecule # ',
     #               nml,': Number of moving sets > ',mxms
                     stop
                  endif
                  jns=jns+1
                  do k=nms,j+2,-1  ! shift ranges of m.s.
                    latms1(k)=latms1(k-1)
                    latms2(k)=latms2(k-1)
                  enddo
                  do ko=jo+1,io
                    k=iorvr(ko)
                    imsvr1(k)=imsvr1(k)+1
                    imsvr2(k)=imsvr2(k)+1
                  enddo
                  latms2(j+1)=j2
                endif
                latms1(j+ja)=i2+1
                ja=ja+1
              endif
              if (ja.eq.0) then  ! -1 moving set
                nms=nms-1
                jns=jns-1
                do k=j,nms
                  latms1(k)=latms1(k+1)
                  latms2(k)=latms2(k+1)
                enddo
                do ko=jo+1,io
                  k=iorvr(ko)
                  imsvr1(k)=imsvr1(k)-1
                  imsvr2(k)=imsvr2(k)-1
                enddo
              else
                j=j+ja
              endif
            else  ! No overlap
              j=j+1
            endif 
          enddo  ! mov. sets for 'jv'
          imsvr2(jv)=j1s+jns-1

        enddo  ! prev. variables
c _______________________________ Finally, add moving set for 'iv'
        nms=nms+1
        latms1(nms)=i1
        latms2(nms)=i2
      enddo  ! variables
      nmsml(nml)=nms-imsml1(nml)+1
c _____________________________ Determine index of moving set for each atom
      do ia=ifiat,ilaat
        ixmsat(ia)=0
      enddo
      do is=imsml1(nml),nms
        do ia=latms1(is),latms2(is)
          ixmsat(ia)=is
        enddo
      enddo
c _____________________________ Determine indices of variables which moving
c                               set sets have to be added (=are related) to
c                               those of a given variable

      i=iorvr(ifivr)  ! initialize index of CURRENT var.
      ii=imsvr1(i)    !    -"-     index of its 1st m.s

      do io=ifivr,ilavr-1  ! ________ loop over variables

        ic=i               ! save index of CURRENT var.
        ia=iatvr(i)        ! ist primar.mv.atom
        ib=iowat(ia)       ! its base
        it=ityvr(i)        ! its type
        is=ii              ! index of its 1st m.s

        n=nad+1
        iadvr1(i)=n         ! # of its 1st 'added' var.

        i=iorvr(io+1)       ! index of next-in-order var.
        ii=imsvr1(i)        ! index of its 1st m.s

        do jo=io+1,ilavr  ! ______ over following-in-order var.
          j=iorvr(jo)        ! index of var.
          ja=iatvr(j)        ! its prim.mv.at
          jb=iowat(ja)       ! its base

c _______________ current var. is torsion & shares base with var. 'j'
          if (it.eq.3.and.jb.eq.ib) then
            do k=n,nad  !  ? has this branch been registered before ?
              if (iatvr(ladvr(k)).eq.ja) goto 3
            enddo
            nad=nad+1
            if (nad.gt.mxvr) then
              write (*,'(a,i4,a,i5)') ' setmvs> Molecule # ',nml,
     #                         ': Number of added variables > ',mxvr
              stop
            endif
            ladvr(nad)=j   ! save index of 'added' variable
          endif

    3     if (is.lt.ii) then   ! _____ current var. has any m.s:
            do k=is,ii-1       ! ? base of var. 'j' within m.s ?
              if (latms1(k).le.jb.and.jb.le.latms2(k)) then
                do l=n,nad
                  if (iatvr(ladvr(l)).eq.ja) goto 4
                enddo
                nad=nad+1
                if (nad.gt.mxvr) then
                  write (*,'(a,i4,a,i5)') ' setmvs> Molecule # ',nml,
     #                         ': Number of added variables > ',mxvr
                  stop
                endif
                ladvr(nad)=j
              endif
    4       enddo
          else                 ! _____ current var. has no m.s:
            if (ja.eq.ia) then ! ? share prim.mv.at with var. 'j' ?
              do k=n,nad
                if (iatvr(ladvr(k)).eq.ja) goto 5
              enddo
              nad=nad+1
              if (nad.gt.mxvr) then
                write (*,'(a,i4,a,i5)') ' setmvs> Molecule # ',nml,
     #                       ': Number of added variables > ',mxvr
                stop
              endif
              ladvr(nad)=j
            endif
          endif
    5   enddo  ! ... following-in-order variables
        iadvr2(ic)=nad       ! last 'added' var. for current var.
      enddo  ! ... variables

      iadvr1(i)=nad+1  ! don't forget last variable
      iadvr2(i)=nad

      nadml(nml)=nad-iadml1(nml)+1
c _____________________________________ Summary
c      do io=ilavr,ifivr,-1
c        iv=iorvr(io)
c        ib=iowat(iatvr(iv))
c        i1s=imsvr1(iv)
c        i2s=imsvr2(iv)
c        if (i1s.le.i2s) then
c          do i=i1s,i2s
c            i1=latms1(i)
c            i2=latms2(i)
c            if (i.eq.i1s) then
c              write (*,'(a,i3,7a,i4,3a,i4,a)') 'res # ',nursvr(iv),
c     #        ' var: ',nmvr(iv),' base:',nmat(ib),'    atoms= ',
c     #        nmat(i1),'(',i1,') - ',nmat(i2),'(',i2,')'
c            else
c              write (*,'(39x,2a,i4,3a,i4,a)')
c     #        nmat(i1),'(',i1,') - ',nmat(i2),'(',i2,')'
c            endif
c          enddo
c        else
c          write (*,'(a,i3,5a)') 'res # ',nursvr(iv),
c     #    ' var: ',nmvr(iv),' base:',nmat(ib),'  No atoms'
c        endif
c        i1a=iadvr1(iv)
c        i2a=iadvr2(iv)
c        if (i1a.le.i2a) then
c          write (*,'(a,30(1x,a))') ' Depending variables:',
c     #                    (nmvr(ladvr(i)),i=i1a,i2a)
c        else
c          write (*,'(a)') ' No dep. variables'
c        endif
c      enddo
c _____________________________________ Summary - End
 
      return

    6 write (*,'(a,i4,/,2(a,i5),a)') 
     # ' setmvs> Error in atom numbering of molecule # ',nml,
     # ': atom ranges for variables # ',iv,' and # ',jv,
     # ' overlap only PARTLY'
      stop

      end
c *******************************************************
      subroutine fndbrn(nml,nrs,ifirg,ilarg,irg1,irg2,bb)

c .........................................................
c PURPOSE: determine range [ifirg,ilarg] of atom indices
c          for branch starting from atom 'ifirg' of residue
c          'nrs' in molecule 'nml'
c OUTPUT:  BB          - .t. if 'ifirg' is a backbone atom
c          IRG1 & IRG2 - atom indices of ring-closing bond,
c                        if 'ifirg' is INSIDE a ring, but NOT
c                        its 1st atom ( in 'multiple' rings
c                        only LAST closing bond is given !)
c
c CALLS: none
c
c .........................................................

      include 'INCL.H'

      logical bb
      dimension ibd(4)

      ilarg=ifirg

      bb=.false.
      irg1=0

      ifi=iatrs1(nrs)
      ila=iatrs2(nrs)
      ixt=ixatrs(nrs)

      if (ifirg.eq.ifi) then ! = 1st mainchain atom
        bb=.true.
        if (nrs.ne.irsml1(nml)) then
          ilarg=ila
        else  ! 1st residue of 'nml'

          ibd(1)=iowat(ifirg)
          ibd(2)=ibdat(1,ifirg)
          ibd(3)=ibdat(2,ifirg)
          ibd(4)=ibdat(3,ifirg)

          il=0
          do i=1,nbdat(ifirg)+1
            ib=ibd(i)
            if (ib.gt.il.and.iowat(ib).eq.ifirg) il=ib
          enddo
          if (il.gt.0) ilarg=il-1
        endif
      else
        if (ifirg.eq.ixt) bb=.true.
        do i=1,nbdat(ifirg)      ! ______________ check bonds
          ib=ibdat(i,ifirg)
          if (iowat(ib).eq.ifirg) then  ! branch
            do j=ib,ila
              if (j.gt.ib.and.iowat(j).lt.ib) goto 1
              if (j.eq.ixt) bb=.true.
              do k=1,nbdat(j)
                jb=ibdat(k,j)
                if (jb.lt.ifirg) then  ! ring
                  irg1=j
                  irg2=jb
                endif
              enddo
              ilarg=j
            enddo  ! ... branch atoms
          elseif (ib.lt.ifirg) then  ! ring
            irg1=ifirg
            irg2=ib
          endif
    1   enddo  ! ... bonds
      endif

      return
      end


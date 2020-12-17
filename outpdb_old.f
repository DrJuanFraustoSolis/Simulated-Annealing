c **************************************************************
c
c This file contains the subroutines: outpdb
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      subroutine outpdb(nml,fileName)

c ..............................................
c  PURPOSE:  write coordinates of molecule 'nml'
c            in PDB-format (with specialities for hydrogens)
c
c  INPUT:    nml - number of molecule
c
c            npdb - unit of output-file
c
c  CALLS:    toupst,iendst
c ..............................................

      include 'INCL.H'

      dimension ibd(4)
      character chid,cdin,res*3,atnm*5,linty*6,linout*80
      character*(*) fileName


      cdin=' '       !!! residue insert code
      occ=one        !!! occupancy
      bva=zero       !!! B-value

      i0 = ichar('0')
      i9 = ichar('9')


      if (nml.lt.0.or.nml.gt.ntlml) then
        write(*,*) ' outpdb>  No such molecule #',nml,' !'
        return
      elseif (nml.gt.0) then
        im1 = nml
        im2 = nml
      else
        im1 = 1
        im2 = ntlml
      endif
	
      if (ibegst(fileName).gt.0) then
        iout = 99
	write(*,*) 'aqui se quedo el programa'
        open(iout, file=fileName)
      else
        iout = 6
      endif

      iat=0

	
      do iml = im1,im2

        if (ntlml.eq.1) then
          chid=' '
        else
          chid = char(64 + iml)
        endif

        irs=0
        ifirs=irsml1(iml)
        ifiat=iatrs1(ifirs)

        do nrs=ifirs,irsml2(iml)

          irs=irs+1
          res(1:)=seq(nrs)(1:3)

          if (res.ne.'ace'.and.res.ne.'nme') then
            linty = 'ATOM  '
          else
            linty = 'HETATM'
          endif

          do i=iatrs1(nrs),iatrs2(nrs)
            iat=iat+1

            atnm=' '
            atnm(2:5)=nmat(i)

            if (atnm(2:2).eq.'h') then ! hydrogens by PDB convention

              j = iendst(atnm)
              if (ichar(atnm(j:j)).ge.i0.and.
     #            ichar(atnm(j:j)).le.i9) then
                atnm(1:1)=atnm(j:j)
                atnm(j:j)=' '
              endif

            endif

            call toupst(atnm)
            call toupst(res)
 
            linout = ' '
            write (linout,1) linty,iat,atnm,res(1:3),chid,irs,cdin,
     #                       xat(i),yat(i),zat(i),occ,bva
            write(iout,'(a80)') linout

          enddo  ! atoms
        enddo  ! residues

        iat = iat + 1
        linout = ' '
        write(linout,2) 'TER   ',iat,res(1:3),chid
        write(iout,'(a80)') linout

      enddo  ! molecules

c ______________________________________ connectivity
c                                        ( only bonds i-j with i<j)

      do iml = im1,im2

        ifirs=irsml1(iml)

        if (nml.gt.0) then
          iat = iatrs1(ifirs) - 1
        else
          iat = 1 - iml
        endif

        do nrs=ifirs,irsml2(iml)
          nfi=iatrs1(nrs)

          do i=nfi,iatrs2(nrs)

            if (nbdat(i).gt.0) then
              if (nrs.eq.ifirs.and.i.eq.nfi) then
                ibd(1)=iowat(i)
                ibd(2)=ibdat(1,i)
                ibd(3)=ibdat(2,i)
                ibd(4)=ibdat(3,i)
                jj=4
              else
                ibd(1)=ibdat(1,i)
                ibd(2)=ibdat(2,i)
                ibd(3)=ibdat(3,i)
                jj=3
              endif
              nbd=0
              do j=1,jj
                if (ibd(j).gt.i) then
                  nbd=nbd+1
                  ibd(nbd)=ibd(j)
                endif
              enddo

              if (nbd.gt.0) then
                linout = ' '
                write (linout,3) 'CONECT',(i-iat),((ibd(j)-iat),j=1,nbd)
                write(iout,'(a80)') linout
              endif

            endif  ! bonds

          enddo  ! atoms
        enddo  ! residues

      enddo  ! molecules

      linout = ' '
      write (linout,'(a3)') 'END'
      write(iout,'(a80)') linout
      
      close(iout)

    1 format (a6,i5,1x,a5,a3,1x,a1,i4,a1,3x,3f8.3,2(1x,f5.2))
    2 format (a6,i5,6x,a3,1x,a1)
    3 format (a6,5i5)

      return
      end

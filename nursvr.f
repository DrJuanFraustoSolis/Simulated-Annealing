c**************************************************************
c
c This file contains the subroutines: nursvr, nursat
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************
      integer*4 function nursvr(ivr)

c ...........................................................
c  PURPOSE: defines index of residue for given variable 'ivr'
c
c  CALLS: none
c
c ...........................................................
      include 'INCL.H'

      do i=ntlml,1,-1
        ifirs=irsml1(i)
        if (ivr.ge.ivrrs1(ifirs).and.nvrml(i).gt.0) then
          do j=irsml2(i),ifirs,-1
            if (ivr.ge.ivrrs1(j).and.nvrrs(j).gt.0) then
              nursvr=j
              return
            endif
          enddo
        endif
      enddo

      write (*,'(a,i5)') ' nursvr > Cannot find variable # ',ivr
      stop

      end

c **********************************
      integer*4 function nursat(iat)

c .......................................................
c  PURPOSE: defines index of residue for given atom 'iat'
c .......................................................

      include 'INCL.H'

      do i=1,ntlml

        ifirs=irsml1(i)
        ilars=irsml2(i)

        if (iat.ge.iatrs1(ifirs).and.iat.le.iatrs2(ilars)) then

          do j=ifirs,ilars

            if (iat.ge.iatrs1(j).and.iat.le.iatrs2(j)) then

              nursat=j

              return
            endif

          enddo

        endif
      enddo

      write (*,'(a,i5)') ' nursat > Cannot find atom # ',iat
      stop

      end

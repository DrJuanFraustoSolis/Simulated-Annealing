C **************************************************************
c
c This file contains the subroutines: gradient
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************


      subroutine gradient()

c -------------------------------------------
c PURPOSE: calculate energy & gradients
c
c CALLS:   opeflx,opereg,opeshe,opesol,setvar
c -------------------------------------------

      include 'INCL.H'


      esm = 0.d0

      do i = 1,ntlml  ! molecules

        call setvar(i,vlvr)  ! set variables & rebuild

        if (flex) then
          call opeflx(i)
        else
          call opeshe(i)
        endif

        esm = esm + eysm

        if (itysol.lt.0) then

          call opesol(i)
          esm = esm + eysl

          ivr1=ivrml1(i)
          ivr2=ivr1+nvrml(i)-1

          do j=ivr1,ivr2
            gdeyvr(j) = gdeyvr(j)+gdeysl(j)
          enddo

        else if (itysol.eq.0) then

          eysl = 0.d0
        else

          write(*,*)  'gradient> Set itysol < 0'
          stop
        endif

        if (ireg.eq.1)  call opereg(i)

      enddo

      eysm = esm

      return
      end

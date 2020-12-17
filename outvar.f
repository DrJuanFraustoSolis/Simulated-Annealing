! **************************************************************
!
! This file contains the subroutines: outvar
!
! Copyright 2005       Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************
 
      subroutine outvar(nml,fileName)
  
!--------------------------------------------------------------------
!       Output of variables of the current protein conformation 
!- 
!       nml != 0       :     molecule index
!       nml  = 0       :     all molecules
!       fileName  = '' :     write to standard output
!       fileName != '' :     write into file with name
!
! CALLS: iendst,iopfil,ibegst,nursvr
!-------------------------------------------------------------------

      include 'INCL.H'

      character*(*) fileName

      if (nml.lt.0.or.nml.gt.ntlml) then
        write(*,*) ' outvar>  No such molecule #',nml,' !'
        return
      endif

      if (ibegst(fileName).gt.0) then
        iout = 98
        open(iout, file=fileName, status='unknown')
      else
        iout = 6
      endif

      call outvbs(nml, iout)

      if (iout.ne.6) close(iout)

      return
      end

      subroutine outvbs(nml, iout)

      include 'INCL.H'

      character mlfd*7,strg(6)*17

      if (nml.lt.0.or.nml.gt.ntlml) then
        write(*,*) ' outvbs>  No such molecule #',nml,' !'
        return
      elseif (nml.gt.0) then
        im1 = nml
        im2 = nml
      else
        im1 = 1
        im2 = ntlml
      endif

      do iml = im1,im2

        write(mlfd,'(i4,a3)') iml,' : '

! --------------------------------- global pars.

        if ( gbpr(1,iml).ne.zero
     #   .or.gbpr(2,iml).ne.zero
     #   .or.gbpr(3,iml).ne.zero
     #   .or.gbpr(4,iml).ne.zero
     #   .or.gbpr(5,iml).ne.zero
     #   .or.gbpr(6,iml).ne.zero ) then

          do i = 1,3
            write(strg(i),'(f17.6)') gbpr(i,iml)
          enddo
          do i = 4,6
            write(strg(i),'(f17.6)') (gbpr(i,iml)*crd)
          enddo

          write(iout,'(1x,2a,1x,12a)')
     #    '@ ',mlfd,(strg(i)(ibegst(strg(i)):),',',i=1,5),
     #    strg(6)(ibegst(strg(6)):)

        endif

        ifivr = ivrml1(iml)
        is = irsml1(iml)-1

        do i = ifivr,ifivr+nvrml(iml)-1

          if (.not.fxvr(i)) then

            write(iout,'(3x,a,i3,1x,a,1x,a,1x,a,1x,f10.3)')
     #        mlfd,(nursvr(i)-is),':',nmvr(i),':',vlvr(i)*crd
          else

            it=ityvr(i)

            if (it.eq.3) then

              write(iout,'(3x,a,i3,1x,a,1x,a,1x,a,1x,f10.3,1x,a)')
     #          mlfd,(nursvr(i)-is),':',nmvr(i),':',vlvr(i)*crd
     #          ,' &'
            endif

          endif

        enddo  ! internal vars

      enddo  ! molecules

      end subroutine outvbs

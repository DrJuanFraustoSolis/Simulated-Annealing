c**************************************************************
c
c This file contains the subroutines: regul
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      
      subroutine regul(nml, iter, nsteps, acc)

c ----------------------------------------------------------
c PURPOSE: regularization of PDB-structure into SMMP geometry
c
c          @param nml molecule to be regularized
c          @param iter number of iterations during regularization
c          @param nsteps maximum number of steps in minimization
c          @param acc acceptance criterium for minimization
c
c CALLS:   minim, cnteny, outvar,rmsdopt
c ----------------------------------------------------------

      include 'INCL.H'
      include 'INCP.H'
      
cf2py intent(in) nml
cf2py intent(in) iter
cf2py intent(in) nsteps
cf2py intent(in) acc

      dimension rm(3,3),av1(3),av2(3)
      logical ishy(mxvr),fxvro(mxvr)



      wtrg = 1.d0
      wtey = 0.d0

      write(*,'(/,a,2(a,f4.2),/)')
     #  ' ====================== Regularization only',
     #  '   Wt(energy) = ',wtey,'  Wt(regul.) = ',wtrg

      call minim(1, nsteps, acc)

      write(*,*) ' '
      write(*,*) ' ---------- contacts after 1st regularization'
      write(*,*) ' '
      call cnteny(nml)
      write(*,*) ' '

      nrs = irsml2(nml)-irsml1(nml)+1
      call rmsdopt(nml,1,nrs,ixatp,xatp,yatp,zatp,0,rm,av1,av2,rmsd)

      write(*,*) ' RMSD = ',rmsd

c --------------------------------------- fix vars. defined in PDB


      do i = ivrml1(nml),nvrml(nml) 
        fxvro(i) = fxvr(i)  ! save
        if (isrfvr(i)) fxvr(i) = .true.  ! fix vars. defined in ref.str.
      enddo  ! vars.
      ireg = 0

      write(*,'(/,a,2(a,f4.2),/)')
     #  ' ====================== Internal Energy for Hydrogens only',
     #  '   Wt(energy) = ',wtey,'  Wt(regul.) = ',wtrg

      call minim(1, nsteps, acc)

      write(*,*) ' '
      write(*,*) ' ---------- contacs after Emin. for hydrogens'
      write(*,*) ' '
      call cnteny(nml)

      do i = ivrml1(nml),nvrml(nml) 
        fxvr(i) = fxvro(i)  ! restore
      enddo  ! vars.
      ireg = 1


      wtrg = 1.d0
      wtey = 0.d0

      n=iter            
      dn=1.d0/dble(n)

      do it = 1,n

        wtrg = 1.d0 - dn*dble(it)
        wtey = 1.d0 - wtrg

        write(*,'(/,a,i2,2(a,e11.3),/)')
     #    ' ================ Minimization #',it,
     #        '   Wt(energy) = ',wtey,'  Wt(regul.) = ',wtrg

        call minim(1, nsteps, acc)

        nrs = irsml2(nml)-irsml1(nml)+1
        call rmsdopt(nml,1,nrs,ixatp,xatp,yatp,zatp,0,rm,av1,av2,rmsd)

        write(*,*) ' '
        write(*,*) ' RMSD = ',rmsd

      enddo

      write(*,*) ' '
      write(*,*) ' ---------- contacts after full regularization'
      write(*,*) ' '
      call cnteny(nml)

c      call outpdb(nml,12)

c Output of dihedral angles of the regularized structure
      write(*,*) 'Dihedral angles of the regularized structure;'
      call outvar(nml, 'regd.var')

      ireg = 0

      return
      end


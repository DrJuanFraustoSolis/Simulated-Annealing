c **************************************************************
c
c This file contains the subroutines: canon,can_weight 
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      subroutine  canon(nequi, nswp, nmes, temp, lrand)
C -----------------------------------------------------------------
C PURPOSE: CANONICAL SIMULATION OF PROTEINS USING METROPOLIS UPDATES
C
C CALLS:  addang,energy,metropolis,hbond,helix,outvar,outpdb,rgyr
C
C-------------------------------------------------------------------
      include 'INCL.H'

cf2py intent(in) nequi
cf2py intent(in) nswp
cf2py intent(in) nmes
cf2py intent(in) temp
cf2py logical optional, intent(in):: lrand = 1

c     external rand
      external can_weight
     
      logical lrand
c      parameter(lrand=.false.)
c      parameter(nequi=10, nswp=1000,nmes=10)
c      parameter(temp=300.0)
C     lrand=.true.: creates random start configuration 
C     nequi: Number of sweeps for equilibrisation of system
      integer nequi
C     nswp:  Number of sweeps for simulation run
      integer nswp
c     nmes:  Number of sweeps between measurments
      integer nmes
C     temp:  Temperature of simulation
      double precision temp
C
!      common/bet/beta

      character*80 file

c     Define files for output:
      open(13,file='time.d')

 
      beta=1.0/ ( temp * 1.98773d-3 )

c _________________________________ random start
      if(lrand) then
       do i=1,nvr
        iv=idvr(i)  ! provides index of non-fixed variable
        dv=axvr(iv)*(grnd()-0.5)
        vr=addang(pi,dv)
        vlvr(iv)=vr
       enddo
      end if

      eol = energy()
      write (*,'(a,e12.5,/)')  'energy of start configuration:',eol

C Write start configuration in pdb-format into file
      call outpdb(0,'start.pdb')

c =====================Equilibration by  Metropolis
      acz = 0.0d0
      do nsw=1,nequi
         call metropolis(eol,acz,can_weight)
      end do
      write(*,*) 'Energy after equilibration:',eol

C======================Simulation in canonical ensemble
      acz = 0.0d0
      do nsw=0,nswp
        call metropolis(eol,acz,can_weight)
c
        if(mod(nsw,nmes).eq.0) then
C Measure radius of gyration and end-to-end distance
C rgy: radius of gyration
C ee:  end-to-end distance
         call rgyr(1,rgy,ee)
C Measure helicity 
C nhel: number of helical residues
c mhel: number of helical segments
c nbet: number of sheet-like residues
c mbet: number of sheet-like segments
         call helix(nhel,mhel,nbet,mbet) 
C Measure number of hydrogen bonds (mhb)
        do i=1,ntlml
         call hbond(i,mhb,0)
        end do
C Write down information on actual conformation
         write(13,'(i5,2f12.3,5i7)')  nsw,  eol, rgy,
     &                              nhel,mhel,nbet,mbet,mhb
        end if
C
      end do

      acz = acz/dble(nsw*nvr)
      write(*,*) 'acceptance rate:',acz
      write(*,*)
c ------------ Output Dihedreals of final configuration
      write(*,*) 'last energy',eol
      call outvar(0,'lastconf.var')
C     Output final conformation as pdb-file
      call outpdb(0,'final.pdb')

      close(11)
      close(12)
      close(13)
c =====================


       end

c ********************************************************
      real*8 function can_weight(x)
c
c CALLS: none
c

      implicit real*8 (a-h,o-z) 

      common/bet/beta

      can_weight = beta*x

      return

      end

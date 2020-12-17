c **************************************************************
c
c This file contains the subroutines: mulcan_sim,muca_weight2
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      subroutine  mulcan_sim
C
C PURPOSE: PERFORM A MULTICANONICAL SIMULATION
C REQUIRES AS INPUT THE MULTICANONICAL PARAMETER AS CALCULATED
C BY THE SUBROUTINE mulcan_par
C
c CALLS: addang, contacts,energy,metropolis
c
      include 'INCL.H'

c     external rand
      external muca_weight2

      logical restart

      parameter(restart=.false.)
      parameter(kmin=-12,kmax=20,ebin=1.0d0)
      parameter(nsweep=100000,nequi=100)
      Parameter(nsave=1000,nmes=10)
C
C     restart: .true. =  restart of simulation
C              .false. = start of simulation with random configuration
C     kmin,kmax: Range of multicanonical parameter
C     ebin:      bin size for multicanonical parameter
C     nequi: Number of sweeps for equilibrisation
C     nsweep:  Number of sweeps for simulation run
C     nsave:  Number of sweeps after which actual configuration is saved
C             for re-starts
C     nmes: Number of sweeps between measurments
C

      dimension xhist(kmin:kmax),ihist(kmin:kmax)
      common/muca2/b(kmin:kmax),alpha(kmin:kmax)

 
C FILE with last conformation (for re-starts)
      open(8,file='EXAMPLES/start.d')
C File with contact map of reference configuration
      open(9,file='EXAMPLES/enkefa.ref')
C File with multicanonical parameter
      open(10,file='EXAMPLES/muca.d')
C Result file: Time series of certain quantities
      open(11, file='EXAMPLES/time.d')


      do j=kmin,kmax
       ihist(j)= 0
      end do

      nresi=irsml2(1)-irsml1(1) + 1
c     nresi:  Number of residues

C READ  REFERENCE CONTACT MAP
      nci = 0
      do i=1,nresi
       read(9,*) (iref(i,j) , j=1,nresi)
      end do
      do i=1,nresi
       do j=nresi,i+3,-1
        if(iref(i,j).eq.1) nci = nci + 1
       end do
      end do
      write(*,*) 'Number of contacts in reference conformation:',nci

C READ IN FIELDS WITH MULTICANONICAL PARAMETER
      Do j=kmin,kmax
       read(10,*) i,b(i),alpha(i)
      end do
C

      if(restart) then
       read(8,*) nswm, eol_old
       read(8,*) (xhist(j), j=kmin,kmax)
       do i=1,nvr
        read(8,*) j, x
        iv = idvr(j)
        vlvr(iv) = x
       end do
       write(*,*) 'Last iteration, energy:',nswm,eol_old
      else
c _________________________________ random start
       do i=1,nvr
        iv=idvr(i)  ! provides index of non-fixed variable
        dv=axvr(iv)*(grnd()-0.5)
        vr=addang(pi,dv)
        vlvr(iv)=vr
       enddo
      end if
c
      eol = energy()
      write (*,'(e12.5,/)')  eol
      call contacts(nhy,nhx,dham)
      write(*,*) 'Number of contacts in start configuration:',nhy
      write(*,*) 'Number of native contacts in start configuration:',
     &            nhx
      do i=1,nresi
       write(*,'(62I1)') (ijcont(i,j), j=1,nresi)
      end do
      write(*,*)
C

      
      if(.not.restart) then
c =====================Equilibrization by  Metropolis
       do nsw=1,nequi
        call metropolis(eol,acz,muca_weight2)
       end do
       do i=kmin,kmax
        ihist(i) = 0
       end do 
       nswm = 1
      end if

C======================Simulation
      acz = 0.0d0
C LOOP OVER SWEEPS
      do nsw=nswm,nsweep
C
C METROPOLIS UPDATE
       call metropolis(eol,acz,muca_weight2)
       muold = min(kmax,max(kmin,int(eol/ebin+sign(0.5d0,eol))))
       ihist(muold) = ihist(muold) + 1
C
C  SAVE ACTUAL CONFORMATIONS FOR RE-STARTS:
       if(mod(nsw,nsave).eq.0) then
        rewind 8
        write(8,*) nswm, eol
        write(8,*) (xhist(j), j=kmin,kmax)
        do i=1,nvr
         iv = idvr(i)
         write(8,*) i,vlvr(iv)
        end do
       end if
C Measurements after NMES sweeps
       if(mod(nsw,nmes).eq.0) then
C Take a histogram of energy
        do i=kmin,kmax
         xhist(i) = xhist(i) + ihist(i)
         ihist(i) = 0
        end do
c Calculate contacts in actual configuartion and compare with reference
C configuration
c       call contacts(nhx,nhy,dham)
C nhx : Number of contcats in actual conformation
C nhy : Number of contacts which are identical in actual and reference 
C       configuration
C dham:  Hamming distance between actual and reference configuration      
C
        write(11,'(i7,f12.2,2i8,f12.4)')  nsw,eol,nhx,nhy,dham
       end if
      end do
C END OF SIMULATION

C FINAL OUTPUT:
      acz = acz/dble(nsw*nvr)
      write(*,*) 'last energy',eol
      write(*,*) 'aczeptance rate:',acz

C WRITE DOWN (UN-REWEIGHTED) HISTOGRAM OF MULTICANONICAL SIMULATION
      do i=kmin,kmax
       if(xhist(i).gt.0.0d0) then
        write(*,*) i,xhist(i)
       end if
      end do
c =====================
      close(8)
      close(9)
      close(10)
      close(11)

      return
      end

c ************************************************************
      real*8 function muca_weight2(x)

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      Parameter(kmin=-12,kmax=20,ebin=1.0d0)
 
      common/muca2/b(kmin:kmax),alpha(kmin:kmax)

      muold = min(kmax,max(kmin,int(x/ebin+sign(0.5d0,x))))
      muca_weight2 = b(muold)*x + alpha(muold)

      return

      end



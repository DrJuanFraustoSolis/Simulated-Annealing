! **************************************************************
!
! This file contains the subroutines: mulcan_par, muca_weight
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************

      subroutine  mulcan_par
!
! PURPOSE: CALCULATION OF ESTIMATORS FOR MULTICANONICAL WEIGHTS
!          USING A METHOD DESCRIBED FIRST IN: Bernd Berg, 
!          {\it J.~Stat.~Phys.} {\bf 82}, 331~(1996). THE METHOD
!          IS STABLE, BUT MAY NOT ALWAYS BE THE FASTEST WAY OF
!          GETTING PARAMETERS
!
! CALLS: addang,energy,metropolis
!

      include 'INCL.H'

!     external rand
      external muca_weight

      logical  l_iter

      parameter(l_iter=.false.)
      parameter(kmin=-12,kmax=20,ebin=1.0d0)
      parameter(xmin = kmin,xmax=kmax)
      parameter(nsweep=100000,nup=5000)
      parameter(temp=1000.0)
!
!     l_iter: .true. =  for iteration of multicanonical parameters
!             .false. = starts new calculation of multicanonica parameters
!     kmin,kmax: Range of multicanonical parameter
!     ebin:      bin size for multicanonical parameter
!     nup: Number of sweeps  between updates of multicanonical parameters
!     nsweep:  Number of sweeps for simulation run
!     temp: temperature for initial canonical run 
!

      dimension xhist(kmin:kmax),g1(kmin:kmax),ent(kmin:kmax)
      common/muca/b(kmin:kmax),alpha(kmin:kmax),beta
      common/jhist/ihist(kmin:kmax)


! Files with multicanonical parameter
      open(11,file='EXAMPLES/mpar_full.d')
      open(12,file='EXAMPLES/muca.d')


      beta=1.0/ ( temp * 1.98773d-3 )
      do j=kmin,kmax
       ihist(j)= 0
      end do


      if(l_iter) then
! READ IN FIELDS WITH MULTICANONICAL PARAMETER
       Do j=kmin,kmax
        read(11,*) i,xhist(i),g1(i),b(i),alpha(i),ent(i)
       end do
      else
       do j=kmin,kmax
        b(j) = 0.0d0
        alpha(j) = 0.0d0
        g1(j) = 0.0d0
        xhist(j) = 0.0d0
        ent(j) = 0.0d0
       end do
      end if
!
! _________________________________ random start
       do i=1,nvr
        iv=idvr(i)  ! provides index of non-fixed variable
        dv=axvr(iv)*(grnd()-0.5)
        vr=addang(pi,dv)
        vlvr(iv)=vr
       enddo
!
      eol = energy()
      write (*,'(a,e12.5,/)')  'Energy of start configuration: ',eol
      write(*,*)

      call outpdb(1, 'start.pdb')
!

      
!======================Simulation
      acz = 0.0d0
! LOOP OVER SWEEPS
      do nsw=1,nsweep
!
! METROPOLIS UPDATE
       call metropolis(eol,acz,muca_weight)

       muold = int(min(xmax,max(xmin,eol/ebin+sign(0.5d0,eol))))
       ihist(muold) = ihist(muold) + 1
       write (*,*) nsw, eol, acz
!
! ITERATE MULTICANONICAL WEIGHTS EVERY NUP SWEEPS 
       if(mod(nsw,nup).eq.0) then
        xhist(kmax) = xhist(kmax) + ihist(kmax)
        xhist(kmax-1) = xhist(kmax-1) + ihist(kmax-1)
        rewind 11
        do i=kmax-2,kmin,-1
         xhist(i) = xhist(i) + ihist(i)
         if((ihist(i).gt.1).and.(ihist(i+1).gt.1)) then
          g0 = float(ihist(i+1))*float(ihist(i))
          g0 = g0/float(ihist(i)+ihist(i+1))
          g1(i) = g1(i) + g0
          g0 = g0/g1(i)
          bi_1 = log(float(ihist(i+1)))
          bi = log(float(ihist(i)))
          b(i) = b(i) + g0*(bi_1 - bi)/ebin
         else
          if(xhist(i).lt.2.0) then
           b(i) = max(0.0d0,b(i+1))
          end if
         end if
        end do
        b(kmax) = b(kmax-2) 
        b(kmax-1) = b(kmax-2) 
        alpha(kmax) = 0
        ent(kmax) = (beta+b(kmax))*float(kmax)*ebin + alpha(kmax)
        write(11,'(i5,5f15.4)') kmax,xhist(kmax),g1(kmax),
     &                          b(kmax),alpha(kmax),ent(kmax)
        ihist(kmax) = 0
        do i=kmax-1,kmin,-1
         ihist(i) = 0
         alpha(i) = alpha(i+1) + (b(i+1) - b(i))*(i+0.5d0)*ebin
         ent(i) = (beta+b(i))*float(i)*ebin + alpha(i)
         write(11,'(i5,5f15.4)') i,xhist(i),g1(i),b(i),alpha(i),ent(i)
        end do
       end if 
!
      end do
! END OF SIMULATION

! FINAL OUTPUT:
      acz = acz/dble(nsw*nvr)
      write(*,*) 'last energy',eol
      write(*,*) 'aczeptance rate:',acz
!
! OUTPUT OF FINAL HISTOGRAM
!
      write(*,*) 'Histogram:'
      do i=kmin,kmax
       if(xhist(i).gt.0.0d0) then
        write(*,*) i,xhist(i)
       end if
      end do
!
! WRITE DOWN MULTICANONICAL PARAMETER
      alpha(kmax)=0.0d0
      b(kmax) = b(kmax)+beta
      do i= kmax-1,kmin,-1
       b(i) = b(i) + beta
       alpha(i) = alpha(i+1) + (b(i+1) - b(i))*(i+0.5d0)*ebin
      end do
      rewind 12
      do i=kmin,kmax
       write(12,*) i,b(i),alpha(i)
      end do
  
      close(11)
      close(12)
!
! =====================
!
      return
      end

!**********************************************************************

      real*8 function muca_weight(x)

      implicit real*8 (a-h,o-z)
      implicit integer*8 (i-n)

      parameter(kmin=-12,kmax=20,ebin=1.0d0)

      common/muca/b(kmin:kmax),alpha(kmin:kmax),beta
    
      xmin = kmin
      xmax = kmax 
      muold = int(min(xmax,max(xmin,x/ebin+sign(0.5d0,x))))
      muca_weight = (beta+b(muold))*x + alpha(muold)

      return

      end
 
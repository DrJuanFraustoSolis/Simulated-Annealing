! **************************************************************
!
! This file contains the module multicanonical with the 
! subroutines: mulcan_par, muca_weight, mulcan_sim, and 
! mulcan_weight2
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku Hu
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! Changed to module by Jan H. Meinke
! 
! **************************************************************

module multicanonical

      real*8, private, allocatable :: b(:), alpha(:)
      real*8, private :: ebin, beta
      real*8, private :: xmin, xmax

      contains
      subroutine  mulcan_par(nsweep, nup, temp, kmin, kmax, binWidth, l_iter)
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
!
!     l_iter: .true. =  for iteration of multicanonical parameters
!             .false. = starts new calculation of multicanonica parameters
!     kmin,kmax: Range of multicanonical parameter
!     ebin:      bin size for multicanonical parameter
!     nup: Number of sweeps  between updates of multicanonical parameters
!     nsweep:  Number of sweeps for simulation run
!     temp: temperature for initial canonical run 
!
      integer, intent(in) :: kmin, kmax, nsweep, nup
      real*8, intent(in) :: binWidth, temp
      logical,  intent(in) :: l_iter

      integer :: i,j
      integer :: iv,nsw,muold
      real*8 :: bi,bi_1,vr,g0,dv
      real*8 :: acz

      real*8, allocatable :: xhist(:),g1(:),ent(:)
      integer, allocatable :: ihist(:)

      xmin = kmin
      xmax = kmax
      ebin = binWidth 

! Allocate the arrays
      allocate(b(kmin:kmax), alpha(kmin:kmax))
      allocate(xhist(kmin:kmax), g1(kmin:kmax), ent(kmin:kmax) )
      allocate(ihist(kmin:kmax))

      if (.not.(allocated(b).and.allocated(alpha).and.allocated(xhist) &
                .and.allocated(g1).and.allocated(ent).and.allocated(ihist))) then
         stop "Unable to allocate memory in mulcan_par. Exiting."
      end if

! Files with multicanonical parameter
      open(11,file='mpar_full.d')
      open(12,file='muca.d')

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
        iv=idvr(i)   ! provides index of non-fixed variable
        dv=axvr(iv)*(grnd()-0.5)
        vr=addang(pi,dv)
        vlvr(iv)=vr
       enddo

      eol = energy()
      write (*,'(a,e12.5,/)')  'Energy of start configuration: ',eol
      write(*,*)

      call outpdb(1, 'start.pdb')

      
!======================Simulation
      acz = 0.0d0
! Loop over sweeps
      do nsw=1,nsweep

       call metropolis(eol,acz,muca_weight)

       muold = int(min(xmax, max(xmin, eol/ebin+sign(0.5d0,eol))))
       ihist(muold) = ihist(muold) + 1
!
! Iterate multicanonical weights every nup sweeps 
       if(mod(nsw,nup).eq.0) then
        xhist(kmax) = xhist(kmax) + ihist(kmax)
        xhist(kmax-1) = xhist(kmax-1) + ihist(kmax-1)
        rewind 11
        do i=kmax-2,kmin,-1
         xhist(i) = xhist(i) + ihist(i)
         if((ihist(i).gt.1).and.(ihist(i+1).gt.1)) then
          g0 = Float(ihist(i+1))*float(ihist(i))
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
        write(11,'(i5,5f15.4)') kmax,xhist(kmax),g1(kmax),              &
     &                          b(kmax),alpha(kmax),ent(kmax)
        ihist(kmax) = 0
        do i=kmax-1,kmin,-1
         ihist(i) = 0
         alpha(i) = alpha(i+1) + (b(i+1) - b(i))*(i+0.5d0)*ebin
         ent(i) = (beta+b(i))*float(i)*ebin + alpha(i)
         write(11,'(i5,5f15.4)') i,xhist(i),g1(i),b(i),alpha(i),ent(i)
        end do
       end if ! parameter update
      end do ! loop over sweeps

! Final output
      acz = acz/dble(nsw*nvr)
      write(*,*) 'last energy',eol
      write(*,*) 'acceptance rate:',acz

      write(*,*) 'Histogram:'
      do i=kmin,kmax
       if(xhist(i).gt.0.0d0) then
        write(*,*) i,xhist(i)
       end if
      end do

! Save multicanonical parameters to file
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

      deallocate(b, alpha)
      deallocate(xhist, g1, ent )
      deallocate(ihist)

      return
      end subroutine mulcan_par


      subroutine  mulcan_sim(nequi, nsweep, nmes, nsave, kmin, kmax, binWidth, restart)
!
! PURPOSE: PERFORM A MULTICANONICAL SIMULATION
! REQUIRES AS INPUT THE MULTICANONICAL PARAMETER AS CALCULATED
! BY THE SUBROUTINE mulcan_par
!
! CALLS: addang, contacts,energy,metropolis
!
      include 'INCL.H'

!     restart: .true. =  restart of simulation
!              .false. = start of simulation with random configuration
!     kmin,kmax: Range of multisvn+ssh://nicole.nic.kfa-juelich.de/home/meinke/repositories/smmp/smmp/trunk/SMMPcanonical parameter
!     ebin:      bin size for multicanonical parameter
!     nequi: Number of sweeps for equilibrisation
!     nsweep:  Number of sweeps for simulation run
!     nsave:  Number of sweeps after which actual configuration is saved
!             for re-starts
!     nmes: Number of sweeps between measurments
      integer, intent(in) :: kmin, kmax, nsweep, nmes, nequi, nsave
      real*8, intent(in) :: binWidth
      logical,  intent(in) :: restart

      real*8, allocatable :: xhist(:),g1(:),ent(:)
      integer, allocatable :: ihist(:)

      integer :: i,j
      integer :: iv,nsw,muold
      real*8 :: bi,bi_1,vr,g0,dv
      real*8 :: acz
      
      xmin = dble(kmin)
      xmax = dble(kmax)
      ebin = binWidth 

      allocate(b(kmin:kmax), alpha(kmin:kmax))
      allocate(xhist(kmin:kmax), g1(kmin:kmax), ent(kmin:kmax) )
      allocate(ihist(kmin:kmax))

      if (.not.(allocated(b).and.allocated(alpha).and.allocated(xhist) &
                .and.allocated(g1).and.allocated(ent).and.allocated(ihist))) then
         stop "Unable to allocate memory in mulcan_sim. Exiting."
      end if

 

! FILE with last conformation (for re-starts)
      open(8,file='start.d')
! File with contact map of reference configuration
! FIXME: This must go. Reference structure needs to be read in main()
      open(9,file='enkefa.ref')
! File with multicanonical parameter
      open(10,file='muca.d')
! Result file: Time series of certain quantities
      open(11, file='time.d')

      do j=kmin,kmax
       ihist(j)= 0
       xhist(j) = 0
      end do

      nresi=irsml2(1)-irsml1(1) + 1 ! nresi:  Number of residues
! Read  reference contact map
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

! Read in fields with multicanonical parameter
      Do j=kmin,kmax
       read(10,*) i,b(i),alpha(i)
      end do


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
! _________________________________ random start
       do i=1,nvr
        iv=idvr(i)  ! provides index of non-fixed variable
        dv=axvr(iv)*(grnd()-0.5)
        vr=addang(pi,dv)
        vlvr(iv)=vr
       enddo
      end if

      eol = energy()
      write (*,'(e12.5,/)')  eol
      call contacts(nhy,nhx,dham)
      write(*,*) 'Number of contacts in start configuration:',nhy
      write(*,*) 'Number of native contacts in start configuration:',   &
     &            nhx
      do i=1,nresi
       write(*,'(62I1)') (ijcont(i,j), j=1,nresi)
      end do
      write(*,*)


      
      if(.not.restart) then
! _________________Equilibrization by  Metropolis
       do nsw=1,nequi
        call metropolis(eol,acz,muca_weight2)
       end do
       do i=kmin,kmax
        ihist(i) = 0
       end do 
       nswm = 1
      end if

!======================Simulation
      acz = 0.0d0
! Loop over sweeps
      do nsw=nswm,nsweep
       call metropolis(eol,acz,muca_weight2)
       muold = min(kmax,max(kmin,int(eol/ebin+sign(0.5d0,eol))))
       ihist(muold) = ihist(muold) + 1

!  SAVE ACTUAL CONFORMATIONS FOR RE-STARTS:
       if(mod(nsw,nsave).eq.0) then
        rewind 8
        write(8,*) nswm, eol
        write(8,*) (xhist(j), j=kmin,kmax)
        do i=1,nvr
         iv = idvr(i)
         write(8,*) i,vlvr(iv)
        end do
       end if
! Measurements after NMES sweeps
       if(mod(nsw,nmes).eq.0) then
! Take a histogram of energy
        do i=kmin,kmax
         xhist(i) = xhist(i) + ihist(i)
         ihist(i) = 0
        end do
! Calculate contacts in actual configuartion and compare with reference
! configuration
!       call contacts(nhx,nhy,dham)
! nhx : Number of contcats in actual conformation
! nhy : Number of contacts which are identical in actual and reference 
!       configuration
! dham:  Hamming distance between actual and reference configuration      
!
        write(11,'(i7,f12.2,2i8,f12.4)')  nsw,eol,nhx,nhy,dham
       end if
      end do ! End of simulation

      acz = acz/dble(nsw*nvr)
      write(*,*) 'last energy',eol
      write(*,*) 'acceptance rate:',acz

! WRITE DOWN (UN-REWEIGHTED) HISTOGRAM OF MULTICANONICAL SIMULATION
      do i=kmin,kmax
       if(xhist(i).gt.0.0d0) then
        write(*,*) i,xhist(i)
       end if
      end do

      close(8)
      close(9)
      close(10)
      close(11)

      deallocate(b, alpha)
      deallocate(xhist, g1, ent )
      deallocate(ihist)

      return
      end subroutine mulcan_sim

! ************************************************************
      real*8 function muca_weight(x)
      
      real*8, intent(in) :: x
      
      real*8 :: beta
      common /bet/beta

      integer :: muold


      muold = int(min(xmax,max(xmin,x/ebin+sign(0.5d0,x))))
      muca_weight = (beta+b(muold))*x + alpha(muold)

      return

      end function muca_weight

! ************************************************************
      real*8 function muca_weight2(x)

      real*8, intent(in) :: x
      
      real*8 :: beta
      common /bet/beta
      
      integer :: muold

      muold = int(min(xmax,max(xmin,x/ebin+sign(0.5d0,x))))
      muca_weight2 = b(muold)*x + alpha(muold)

      return

      end function muca_weight2

end module


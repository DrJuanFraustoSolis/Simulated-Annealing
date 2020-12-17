! **************************************************************
!
! This file contains the subroutines: partem_s
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
!
!
! **************************************************************
      
      subroutine partem_s(num_rep, nequi, nswp, nmes, nsave, newsta,
     &                     switch)
!
! PURPOSE: SIMULATION OF PROTEINS BY PARALLEL TEMPERING ALGORITHM
!          ON SINGLE-PROCESSOR MACHINE.
!
! CALLS: addang,energy,metropolis,rgyr,setvar,
!
      include 'INCL.H'
      external can_weight
! TODO Store global coordinates in pgbpr
      logical :: newsta
      integer :: num_rep, nequi, nswp, nmes, nsave, switch
      real*8 :: ttemp
      character*80 :: filebase, fileNameMP

      parameter(gasc=0.00198773d0)
!     newsta: .true. for re-start; .false. for random start configurations
!     no:     Number of replicas
!     nvrmax: Maximal number of dihdral angles
!     nequi:  Number of sweeps for equilibrisation
!     nswp:   Number of sweeps for simulation
!     nmes:   Number of sweeps between measurements
!     gasc:   Gas constant
!
      real*8, allocatable :: coor_G(:, :),pbe(:),
     & eol(:),acc(:),temp(:), pgbpr(:, :, :)
      integer, allocatable :: ipoi(:)
!     temp: temperatures for each replica
!     coor_G: dihedral angles for ALL replicas
!     pbe:    inverse temperatures for each replica
!     ipoi:   Points to replica ipoi(k) which is currently
!             at inverse temperature pbe(k)
!     eol:    energy of each replica 
!     acc:    accepatance rate of each replica
!
!     beta:   inverse temperature of single replica
!
!
!     Allocate the arrays
      allocate(coor_G(nvr, num_rep))
      allocate(pgbpr(6, ntlml, num_rep))
      allocate(temp(num_rep), pbe(num_rep), eol(num_rep))
      allocate(acc(num_rep), ipoi(num_rep))
      
      if (.not.(allocated(coor_G).and.allocated(temp)
     &          .and. allocated(pbe).and. allocated(eol)
     &          .and. allocated(acc).and. allocated(ipoi))) then
        stop "Unable to allocate memory in partem_s. Exiting."
      end if
!     Initialize arrays
      do i = 1, num_rep
         do j = 1, nvr
            coor_G(j, i) = 0.0
         end do
         do nml = 1, ntlml
            do j = 1, 6
               pgbpr(j, nml, i) = 0.0
            end do
         end do
         temp(i) = 0.0
         pbe(i) = 0.0
         eol(i) = 0.0
         acc(i) = 0.0
         ipoi(i) = 0
      end do
        

!     READ IN TEMPERATURES
      open(11, file='temperatures', status='old')
      do i=1,num_rep
         read(11,*) j,ttemp
         temp(j) = ttemp
         pbe(j) = 1.0d0/(ttemp * gasc)
         ipoi(j) = j
      end do
      close(11)

      open(13,file='time.d',status='unknown')
      
      if(.not.newsta) then
! READ START Values
         open(15,file='par_R.in',form='unformatted')
         do k = 1, num_rep
            read(15) nstart
            read(15) eol(k)
            lunvar = 14
            filebase = 'conf_0000.var'
            varfil = fileNameMP(filebase, 6, 9, k)
            call redvar()
         end do
         close(15)
      else
         nstart=1
         do i=1,num_rep
            if (switch.ne.0) then
               do j=1,nvr
                  iv = idvr(j)
                  if (switch.gt.0) then
                     dv=axvr(iv)*(grnd()-0.5)
                     vr=addang(pi,dv)
                  else
                     vr = pi
                  end if
                  coor_G(j,i)=vr
               end do
            end if
            do nml = 1, ntlml
               do j = 1, 6
                  pgbpr(j, nml, i) = gbpr(j, nml)
               end do
            end do            
         end do
      
!      Equilibrization for each replica (No replica exchange move)
         do k=1,num_rep
            beta=pbe(k)
            do i=1,nvr
               iv =idvr(i)
               vlvr(iv) = coor_G(i,k)
            end do
            do nml = 1, ntlml
               do j = 1, 6
                  gbpr(j, nml) = pgbpr(j, nml, k)
               end do
            end do
            energ = energy()
            do nsw=1,nequi
               CALL METROPOLIS(energ,acz,can_weight)
            end do
            write(*,*) 'Start energy after equilibration for replica:', 
     &                 k, energ
            do i=1,nvr
               iv = idvr(i)
               coor_G(i,k) = vlvr(iv)
            end do
            do nml = 1, ntlml
               do j = 1, 6
                  pgbpr(j, nml, k) = gbpr(j, nml)
               end do
            end do
            eol(k) = energ
         end do
      end if
            
! Now begins the simulation with Multiple Markov Chains
      iswitch = 1
      do nsw=nstart,nswp
!-------------------Ordinary  MC sweep for each replica
         do k1=1,num_rep
            k = ipoi(k1)
            beta = pbe(k1)
            energ = eol(k)
            acz= acc(k)
            do i=1,nvr
               iv = idvr(i)
               vlvr(iv) = coor_G(i,k)
            end do
            do i=1,ntlml
               do j = 1, 6
                  gbpr(j, nml) = pgbpr(j, i, k)
               end do
               call setvar(i,vlvr)
            end do
            CALL METROPOLIS(energ,acz,can_weight)
! 
            if(mod(nsw,nmes).eq.0) then
! Measure and store here all quantities you want to analyse later
! ee: end-to-end distance
! rgy: radius of gyration
               call rgyr(0,ee,rgy)
               temp0=1.0d0/beta/gasc
               write(13,'(i8,i3,4f12.3)')  nsw,k1,temp0,energ,ee,rgy
            end if
!
            do i=1,nvr
               iv = idvr(i)
               coor_G(i,k) = vlvr(iv)
            end do
            do nml = 1, ntlml
               do j = 1, 6
                  pgbpr(j, nml, k) = gbpr(j, nml)
               end do
            end do
        
            eol(k) = energ
            acc(k) = acz
            acz = 0.0d0
         end do
!------------------Exchange of Replicas
         if (.not.num_rep.gt.1) cycle ! Skip exchange if we only have a single replica.
         if(iswitch.eq.1) then
            num_rep1 = num_rep-1
            nu  = 1
            iswitch = 2
         else
            num_rep1 = num_rep
            nu  = 2
            iswitch = 1
         end if
         do i=nu,num_rep1,2     ! labels (inverse) temperatures
            j=i+1
            if(i.eq.num_rep) j=1
            in=ipoi(i)      
            jn=ipoi(j)
            delta=-pbe(i)*eol(jn)-pbe(j)*eol(in)
     &         +pbe(i)*eol(in)+pbe(j)*eol(jn)
            ra = grnd()
            if(ra.le.exp(delta)) then ! Metropolis to switch temperatures
               ipoi(i) = jn             ! between replicas
               ipoi(j) = in
            end if
         end do
!-----------End exchange of Markov chains
!
      end do
!
! Write down final conformations for re-starts
      open(15,file='par_R.in',form='unformatted')
      do k=1,num_rep
         write(15) nswp
         write(15) eol
         do i=1,nvr
            iv = idvr(i)
            vlvr(iv) = coor_G(i,k)
         end do
         do i=1,ntlml
            do j = 1, 6
               gbpr(j, nml) = pgbpr(j, i, k)
            end do
            call setvar(i,vlvr)
         end do
         filebase = "conf_0000.var"
         call outvar(0, fileNameMP(filebase, 6, 9, k))
      end do
      close(15)
      close(13)
      return
      end


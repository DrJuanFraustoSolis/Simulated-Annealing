c**************************************************************
c     
c This file contains the subroutines: partem_p 
C Compared to the version in the main distribution, this
C routine doesn't write the rmsd nor native contacts to the time
C series.
c     
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku Hu
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c     
c     **************************************************************

      subroutine  partem_p(num_rep, nequi, nswp, nmes, nsave, newsta,
     &                     switch, rep_id, partem_comm)
C     
C     PURPOSE: SIMULATION OF PROTEINS BY PARALLEL TEMPERING ALGORITHM 
C     ON PARALLEL COMPUTERS USING MPI
C     
C     switch: Choses the starting configuration:
C     -1 - stretched configuration
C     0 - don't change anything
C     1 - random start configuration
C     
c     CALLS:  addang,contacts,energy,hbond,helix,iendst,metropolis,
c     outvar,(rand),rgyr
C     
      include '../INCL.H'
      include '../INCP.H'
      include '../incl_lund.h'
      include 'mpif.h'

      logical newsta
      integer switch, partem_comm, rep_id, nsave
      
      external can_weight

C     nequi:  number of Monte Carlo sweeps for thermalization
C     nswp:   number of Monte Carlo sweeps
C     nmes:   number of Monte Carlo sweeps between measurments
C     newsta: .true. for new simulations, .false. for re-start

      dimension  eavm(MAX_PROC),sph(MAX_PROC),intem(MAX_PROC),
     &     inode(MAX_PROC), geavm(MAX_PROC), gsph(MAX_PROC)
      double precision    pbe(MAX_PROC),yol(MAX_PROC),acy(MAX_PROC),
     &     acy1(MAX_PROC),acx1(MAX_PROC),
     &     rgyrp(MAX_PROC), eol0,rgyp,acz0
      
      double precision    e_min, e_minp(MAX_PROC), e_minpt(MAX_PROC)
      integer   h_max, h_maxp(MAX_PROC)
c     Order of replica exchange
      integer   odd
!     Counter to keep random number generators in sync
      integer randomCount
      
c     Collect partial energies. Only the root writes to disk. We have to
c     collect the information from the different replicas and provide 
c     arrays to store them.
c     eyslr    storage array for solvent energy
c     eyelp     -      "        - coulomb energy
c     eyvwp     -      "        - van-der-Waals energy
c     eyhbp     -      "        - hydrogen bonding energy
c     eysmi    -      "        - intermolecular interaction energy
      double precision eyslr(MAX_PROC)
      double precision eyelp(MAX_PROC),eyvwp(MAX_PROC),eyhbp(MAX_PROC), 
     &     eyvrp(MAX_PROC),eysmip(MAX_PROC)
c     Collect information about accessible surface and van-der-Waals volume
c     asap      storage array for solvent accessible surface
c     vdvolp     storage array for van-der-Waals volume
      double precision asap(MAX_PROC), vdvolp(MAX_PROC)

      integer nhelp(MAX_PROC),nbetp(MAX_PROC), mhbp(MAX_PROC),
     &     nctotp(MAX_PROC)
      integer imhbp(MAX_PROC)
      character*80 filebase, fileNameMP, tbase0,tbase1
c     frame     frame number for writing configurations
c     trackID   configuration that should be tracked and written out
c     dir          direction in random walk
c     -1 - visited highest temperature last
c     1 - visited lowest temperature last
c     0 - haven't visited the boundaries yet.
c     dirp      storage array for directions.
      integer frame, trackID, dir
      integer dirp(MAX_PROC)

      frame = ifrrm
      trackID = 1
      odd = 1
      write (*,*) 'Starting parallel tempering.'
      write (*,*) 'parameters, ',switch,newsta,nmes,nswp,nmes,
     &            rep_id, num_rep, partem_comm, myrank
      call flush(6)
C     
c     
C     File with temperatures 
      open(11,file='temperatures_abeta',status='old')
      
      tbase0='trj_00000'
      open(18,file=fileNameMP(tbase0,5,9,rep_id),status='unknown')
      if (rep_id.eq.0.and.myrank.eq.0) then
c     File with time series of simulation
         open(14,file='ts.d',status='unknown')
      endif
      
C     READ IN TEMPERATURES
      do i=1,num_rep
         read(11,*) j,temp
         pbe(j) = 1.0d0/( temp * 1.98773d-3 )
      end do
      close(11)

c     nresi:  number of residues
      nresi=irsml2(1)-irsml1(1)+1
C     
C     Initialize variables
      do i=1,num_rep      
         acx1(i) = 0.0d0
         acy(i) = 0.0d0
         eavm(i) = 0.0d0
         sph(i) = 0.0d0
         geavm(i) =0.0d0
         gsph(i) = 0.0d0
         e_minp(i) = 1.0d15
         h_maxp(i) = 0 
         dirp(i) = 0
      end do
      dirp(1) = 1
      dirp(num_rep) = -1
      e_min = 1.0d15
      h_max = 0
      dir = dirp(rep_id + 1)

c     _________________________________ Initialize Variables
      if(newsta) then
         iold=0
         do i=1,num_rep
            inode(i) = i
            intem(i) = i
         end do
c     _________________________________ initialize starting configuration
         if (switch.ne.0) then
            do i=1,nvr
               iv=idvr(i)       ! provides index of non-fixed variable
               if (switch.gt.0) then
                  dv=axvr(i)*(grnd()-0.5)
                  vr=addang(pi,dv)
               else
                  vr = pi 
               endif 
               vlvr(iv)=vr
            enddo
         endif
      else
         if(rep_id.eq.0.and.myrank.eq.0) then
            open(13,file='par_R.in', status='unknown')
            read(13,*) iold
            do i=1,num_rep
               read(13,*) j,inode(i),intem(i),yol(i),e_minp(i),h_maxp(i)
            end do
            jold=(iold/nmes)*num_rep
            rewind 14
            do i=1,jold
               read(14,*) idum1,idum2,idum3,dummy
     &              ,dummy, dummy, dummy
     &              ,dummy, idum1, idum2, idum3
               call flush(6)
            end do
            close(13)
         end if
         CALL MPI_BCAST(IOLD,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(INTEM,num_rep,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(INODE,num_rep,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(YOL,num_rep,MPI_DOUBLE_PRECISION,0,
     #        MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(E_MINP, num_rep, MPI_DOUBLE_PRECISION, 0, 
     #        MPI_COMM_WORLD, IERR)
         CALL MPI_BCAST(h_maxp,num_rep,MPI_INTEGER,0,MPI_COMM_WORLD,
     $        IERR)
      end if
      
      BETA = pbe(inode(rep_id+1))
      e_min = e_minp(inode(rep_id+1))
      h_max = h_maxp(inode(rep_id+1))
      write (*,*) "E_min=",e_min," for ", intem(rep_id + 1) 
      eol=energy()
      if(.not.newsta.and.abs(yol(rep_id + 1) - eol).gt.0.1) then
         write(*,*) rep_id, ' Warning: yol(rep_id).ne.eol:'
         write(*,*) rep_id, yol(rep_id + 1), eol
      endif
C     Start of simulation
      write (*,*) '[',rep_id, myrank, beta, partem_comm,
     &            '] Energy before equilibration:', eol
c     =====================Equilibration by canonical Metropolis
      do nsw=1,nequi
         call metropolis(eol,acz,can_weight)
      end do
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
      write (*,*) '[',rep_id,'] Energy after equilibration:', eol
      call flush(6)
C     
C======================Multiple Markov Chains
      acz = 0
      do nsw=1,nswp
c------------First ordinary Metropolis 
         call metropolis(eol,acz,can_weight)
         iold = iold + 1	
         eol0 = eol
         if(mod(iold,nmes).eq.0) then
            if ((rep_id + 1).eq.trackID.and.myrank.eq.0) then
               frame = iold /nmes
               filebase = "frame_00000.pdb"
               call outpdb(0, fileNameMP(filebase, 7, 11, frame))
            endif
            acz0 = acz
C     Measure global radius of gyration
            call rgyr(0,rgy,ee)  
            rgyp = rgy
C     Measure Helicity and Sheetness
            call helix(nhel,mhel,nbet,mbet)
C     Measure Number of hydrogen bonds
            mhb = 0
            do i = 1, ntlml
               call hbond(i,tmhb,-1) 
               mhb = mhb + 1
            enddo
            call interhbond(imhb) 
C     Measure total number of contacts (NCTOT) and number of
C     native contacts (NCNAT)
            call contacts(nctot,ncnat,dham)
c     Add tracking of lowest energy configuration
            if (eol.lt.e_min) then
c     Write out configuration
               i=rep_id+1
               j=inode(i)
               e_min = eol
               filebase = "c_emin_0000.pdb"
               call outpdb(0, fileNameMP(filebase, 8, 11, j))
               filebase = "c_emin_0000.var"
               call outvar(0, fileNameMP(filebase, 8, 11, j))
               filebase = "c_emin_0000.dat"
               open(15, file=fileNameMP(filebase, 8, 11, j),
     &              status="unknown")
!     write(15,'(i8,2i4,f6.2,2f8.2,5i8)') iold,i,j,pbe(i), 
               write(15,*) iold,i,j,pbe(i), 
     &              eol, eysl, eyel, eyvw, eyhb, eyvr, eysmi,asa,
     &              vdvol, rgy, nhel, nbet, mhb, imhb, nctot
               close(15)
            endif
c     Add tracking of configuration with larges hydrogen contents.
            if ((mhb + imhb).gt.h_max) then
c     Write out configuration
               i = rep_id + 1
               j = inode(i)
               h_max = mhb + imhb
               filebase = "c_hmax_0000.pdb"
               call outpdb(0,fileNameMP(filebase,8,11,j))
               filebase = "c_hmax_0000.var"
               call outvar(0,fileNameMP(filebase,8,11,j))
               filebase = "c_hmax_0000.dat"
               open(15, file=fileNameMP(filebase, 8, 11, j),
     &              status="unknown")
!     write(15,'(i8,2i4,f6.2,2f8.2,5i8)') iold,i,j,pbe(i), 
               write(15,*) iold,i,j,pbe(i), 
     &              eol, eysl, eyel, eyvw, eyhb, eyvr, eysmi,asa,
     &              vdvol, rgy, nhel, nbet, mhb, imhb, nctot
               close(15)
            endif

C     
C--------------------Gather measurement data
! I only use the master node of each replica for data collection. The
! variable partem_comm provides the appropriate communicator.
            if (partem_comm.ne.MPI_COMM_NULL) then
               CALL MPI_GATHER(RGYP,1,MPI_DOUBLE_PRECISION,RGYRP,1,
     &              MPI_DOUBLE_PRECISION, 0,partem_comm,IERR)
               CALL MPI_GATHER(NHEL,1,MPI_INTEGER,NHELP,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(NBET,1,MPI_INTEGER,NBETP,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(MHB,1,MPI_INTEGER,MHBP,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(iMHB,1,MPI_INTEGER,iMHBP,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(NCTOT,1,MPI_INTEGER,NCTOTP,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(dir,1,MPI_INTEGER,dirp,1,MPI_INTEGER,
     &              0,partem_comm,IERR)
               CALL MPI_GATHER(acz0,1,MPI_DOUBLE_PRECISION,acy1,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(e_min,1,MPI_DOUBLE_PRECISION,e_minp,1,
     &              MPI_DOUBLE_PRECISION,0, partem_comm,IERR)
               CALL MPI_GATHER(EOL0,1,MPI_DOUBLE_PRECISION,YOL,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eysl,1,MPI_DOUBLE_PRECISION,eyslr,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eyel,1,MPI_DOUBLE_PRECISION,eyelp,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eyvw,1,MPI_DOUBLE_PRECISION,eyvwp,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eyhb,1,MPI_DOUBLE_PRECISION,eyhbp,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eyvr,1,MPI_DOUBLE_PRECISION,eyvrp,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(eysmi,1,MPI_DOUBLE_PRECISION,eysmip,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(asa,1,MPI_DOUBLE_PRECISION,asap,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
               CALL MPI_GATHER(vdvol,1,MPI_DOUBLE_PRECISION,vdvolp,1,
     &              MPI_DOUBLE_PRECISION,0,partem_comm,IERR)

c     Write trajectory
               write (18,*) '@@@',iold,inode(rep_id+1)
               call outvbs(0,18)
               write (18,*) '###'
!                call flush(18)
c     Write current configuration
               if ((mod(iold, nsave).eq.0)) then
                  filebase = "conf_0000.var"
                  call outvar(0, fileNameMP(filebase, 6, 9, rep_id+1))
               endif
            endif

            if(rep_id.eq.0.and.myrank.eq.0) then
               randomCount = 0
c  Update acceptance, temperature wise average of E and E^2 used to calculate
c  specific heat. 
               do i=1,num_rep
                  j=intem(i)
                  acy(i)=0.0
c  Above: contents of acy1 are added to acy(i) a few lines down. 
c  acy1(intem(i)) contains information received from the node at temperature
c  i, on how many updates have been accepted in node intem(i). Since acz
c  is not reset to 0 every cycle, acy(i) must be set to 0 here. Else, there 
c  will be serious double counting and the values of acceptance printed 
c  will be simply wrong.
                  e_minpt(i)=e_minp(intem(i))
               end do
               do i=1, num_rep
                  e_minp(i) = e_minpt(i)
                  j=intem(i)
                  acy(i)=acy(i)+acy1(j)
                  eavm(i)= eavm(i)+yol(j)
                  sph(i) = sph(i)+yol(j)*yol(j)
               enddo


C     Write measurements to the time series file ts.d
               do i=1,num_rep
                  j=intem(i)
                     write(14,*) iold, i, j, pbe(i),
     &                 yol(j), eyslr(j), eysmip(j),
     &                 rgyrp(j), nhelp(j), nbetp(j), imhbp(j)
                      
               end do
c     Write the current parallel tempering information into par_R.in
               if ((mod(iold, nsave).eq.0))
     &         then
                  open(13,file='par_R.in', status='unknown')
                  write(13,*) iold
                  do i=1,num_rep
                     write(13,*) i,inode(i),intem(i),yol(i),e_minp(i),
     &                    h_maxp(i)
                  end do
C     -------------------------- Various statistics of current run
c               swp=nswp-nequi
                  swp=nsw
                  write(13,*) 'Acceptance rate for change of chains:'
                  do k1=1,num_rep
                     temp=1.0d0/pbe(k1)/0.00198773
                     write(13,*) temp, acx1(k1)*2.0d0*nmes/swp 
c  Above: it's the acceptance rate of exchange of replicas. Since a 
c  replica exchange is attempted only once every nmes sweeps, the 
c  rate should be normalized with (nmes/swp).
                  end do 
                  write(13,*)
                  do k1=1,num_rep
                     k = intem(k1)   
                     temp=1.0d0/pbe(k1)/0.00198773
                     beta = pbe(k1)
                     geavm(k1) = nmes*eavm(k1)/swp
                     gsph(k1)  = (nmes*sph(k1)/swp-geavm(k1)**2)
     #                    *beta*beta/nresi
                     write(13,'(a,2f9.2,i4,f12.3)') 
     &                    'Temperature, Node,local acceptance rate:',
     &                    beta,temp,k,acy(k1)/dble(nsw*nvr)
c  Above: Changed (nswp-nequi) in the denominator of acceptance as 
c  acceptance values are initialized to 0 after equilibration cycles are
c  finished. Note also that since this is being written in the middle of
c  the simulation, it is normalized to nsw instead of nswp.
                     write(13,'(a,3f12.2)') 
     &                    'Last Energy, Average Energy, Spec. Heat:', 
     &                    yol(k),geavm(k1),gsph(k1) 
                     write(13,*) 
                  end do
                  close(13)
! Finally, flush the time series file to ensure that we can do a proper 
! restart.
                  call flush(14)
                  call flush(18)
               end if   

C--------------------Parallel Tempering  update
c     Swap with right neighbor (odd, even)            
               if(odd.eq.1) then
                  nu=1
                  no1 = num_rep-1
c     Swap with left neighbor (even, odd)
               else
                  nu = 2
                  no1 = num_rep
               end if
               do i=nu,no1,2
                  j=i+1
c     Periodic bc for swaps
                  if(i.eq.num_rep) j=1
                  iidx=intem(i)
                  jn=intem(j)
                  wij=exp(-pbe(i)*yol(jn)-pbe(j)*yol(iidx)
     &                 +pbe(i)*yol(iidx)+pbe(j)*yol(jn))
                  rd=grnd()
                  randomCount = randomCount + 1
                  if(wij.ge.rd) then
                     acx1(i) = acx1(i)+1
                     intem(i) = jn
                     intem(j) = iidx
                     inode(in)= j
                     inode(jn)= i
                  end if
               end do
c     ---------------- End Loop over nodes which creates a new temperature
c     map for all nodes, at the node with rank 0. 
c     
               odd = 1 - odd
            end if
c     End of "if (myrank.eq.0) ...". The block above includes PT update and 
c     writing of observables into the time series file etc. 
           
c     Below: Communicate new temperature-node map to all nodes
            CALL MPI_BCAST(INTEM,num_rep,MPI_INTEGER,0,MPI_COMM_WORLD,
     &           IERR)
            CALL MPI_BCAST(INODE,num_rep,MPI_INTEGER,0,MPI_COMM_WORLD,
     &           IERR)
         CALL MPI_BCAST(E_MINP,num_rep,MPI_DOUBLE_PRECISION,0,
     &           MPI_COMM_WORLD,IERR)
            CALL MPI_BCAST(H_MAXP,num_rep,MPI_INTEGER,0,MPI_COMM_WORLD,
     &           IERR)
c Synchronize random number generators for replica 0
            if (rep_id.eq.0) then
               CALL MPI_BCAST(randomCount,1,MPI_INTEGER,0,my_mpi_comm,
     &                IERR)
               if (myrank.ne.0) then
!                  write (*,*) '[', myrank,'] Missed', randomCount, 
!     &                            'random numbers.'
                  do i = 1, randomCount
                     rd = grnd()
!                     write (*,*) '[', myrank,'] rd=', rd
                  enddo
               endif
            endif

            BETA=PBE(INODE(rep_id+1))
            e_min = e_minp(inode(rep_id+1))
            h_max = h_maxp(inode(rep_id+1))
            if (INODE(rep_id + 1).eq.1) dir = 1
            if (INODE(rep_id + 1).eq.num_rep) dir = -1

         endif
c        End of "if (mod(iold,nmes).eq.0) ..."
      end do
c-----------End Loop over sweeps
c     
C     OUTPUT:
C--------------------For Re-starts:
      nu = rep_id + 1
      filebase = "conf_0000.var"
      call outvar(0, fileNameMP(filebase, 6, 9, nu))
      e_final=energy()
      if (partem_comm.ne.MPI_COMM_NULL) then
         write (*,*) rep_id, ' E_final', e_final
      endif
      eol0 = eol
      acz0 = acz
      if (partem_comm.ne.MPI_COMM_NULL) then
         CALL MPI_GATHER(EOL0,1,MPI_DOUBLE_PRECISION,YOL,1,
     #        MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
         CALL MPI_GATHER(acz0,1,MPI_DOUBLE_PRECISION,acy1,1,
     #     MPI_DOUBLE_PRECISION,0,partem_comm,IERR)
      endif
      
      if(rep_id.eq.0.and.myrank.eq.0) then
         close(14)
         open(13,file='par_R.in', status='unknown')
         write(13,*) iold
         do i=1,num_rep
            write(13,*) i,inode(i),intem(i),yol(i),e_minp(i),h_maxp(i)
         end do
C     -------------------------- Various statistics of current run
         swp=nswp
         write(13,*) 'Acceptance rate for change of chains:'
         do k1=1,num_rep
            temp=1.0d0/pbe(k1)/0.00198773
            write(13,*) temp, acx1(k1)*2.0d0*nmes/swp
         end do 
         write(13,*)
         do k1=1,num_rep
            k = intem(k1)   
            temp=1.0d0/pbe(k1)/0.00198773
            beta = pbe(k1)
            geavm(k1) = nmes*eavm(k1)/swp
            gsph(k1)  = (nmes*sph(k1)/swp-geavm(k1)**2)*beta*beta/nresi
            write(13,'(a,2f9.2,i4,f12.3)') 
     &           'Temperature, Node,local acceptance rate:',
     &           beta,temp,k,acy(k1)/dble((nswp)*nvr)
            write(13,'(a,3f12.2)') 
     &           'Last Energy, Average Energy, Spec. Heat:', 
     &           yol(k),geavm(k1),gsph(k1) 
            write(13,*) 
         end do
         close(13)
c         close(16)
      end if
      close(18)

c     =====================
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

      return

      end

!**************************************************************
! Parallel tempering of a system with two Abeta 16--22 peptides
! using 8 replicas.
!
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku Hu
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!     
! CALLS init_energy,init_molecule,partem_p
!     
! **************************************************************
      program pmain

      include '../INCL.H'
      include '../INCP.H'
      include '../incl_lund.h'
      include 'mpif.h'

      character*80 libdir
      character*80 in_fil,ou_fil,filebase, varfile
      character*80 fileNameMP

      character grpn*4,grpc*4
      logical newsta

!c    Number of replicas 
      integer num_replica
!c    Number of processors per replica
      integer num_ppr
!c    Range of processor for crating communicators
      integer proc_range(3)
!c    Array of MPI groups
      integer group(MAX_REPLICA), group_partem
!c    Array of MPI communicators
      integer comm(MAX_REPLICA), partem_comm
!c    Array of nodes acting as masters for the energy calculation.
      integer ranks(MAX_REPLICA)
!c    Configuration switch
      integer switch
      integer rep_id
!     set number of replicas
      double precision eols(MAX_REPLICA)


      common/updstats/ncalls(5),nacalls(5)


!     MPI stuff, and random number generator initialisation

      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myrank,ierr)
      call mpi_comm_size(mpi_comm_world,num_proc,ierr)

!       call VTSetup()
      enysolct = 0
      seed = 8368
      call sgrnd(seed)          ! Initialize the random number generator

!     =================================================== Energy setup
      libdir='../SMMP/'      
!     Directory for SMMP libraries

!     The switch in the following line is now not used.
      flex=.false.              ! .true. for Flex  / .false. for ECEPP

!     Choose energy type with the following switch instead ...
      ientyp = 0
!     0  => ECEPP2 or ECEPP3 depending on the value of sh2
!     1  => FLEX 
!     2  => Lund force field
!     3  => ECEPP with Abagyan corrections
!     

      sh2=.false.               ! .true. for ECEPP/2; .false. for ECEPP3
      epsd=.false.              ! .true. for  distance-dependent epsilon

      itysol= 1                 !  0: vacuum
                                ! >0: numerical solvent energy
                                ! <0: analytical solvent energy & gradients
      isolscl=.false.
      tesgrd=.false.            ! .true. to check analytical gradients

      call init_energy(libdir)

!     calculate CPU time using MPI_Wtime()
      startwtime = MPI_Wtime()


!     ================================================= Structure setup
      grpn = 'nh2'              ! N-terminal group
      grpc = 'cooh'             ! C-terminal group

      iabin = 1                 ! =0: read from PDB-file
                                ! =1: ab Initio from sequence (& variables)

      in_fil='abeta.seq'        ! Sequence file
      varfile = ' '

      newsta=.true.
      boxsize = 1000.0d0    ! Only relevant for multi-molecule systems
      num_replica = 8      ! Number of independent replicas. The file
                            !   temperatures must have at least as many
                            !   entries
      nequi=100             ! Number of MC sweeps before measurements 
                            !   and replica exchanges are started 
      nswp=10000           ! Number of sweeps
      nmes=10               ! Interval for measurements and replica exchange
      nsave=1000            ! Not used at the moment
                    
      switch = -1           ! How should the configuration be   
                            !   initialized?
                            ! -1 stretched chain 
                            !  0 don't do anything
                            !  1 initialize each angle to a random value
      
      ifrm=0
      ntlml = 0

! Decide if and when to use BGS, and initialize Lund data structures 
      bgsprob=0.75    ! Prob for BGS, given that it is possible
! upchswitch= 0 => No BGS 1 => BGS with probability bgsprob 
! 2 => temperature dependent choice 
      upchswitch=1
      rndord=.true.
      
      if (ientyp.eq.2) call init_lundff
!     =================================================================
!     Distribute nodes to parallel tempering tasks
!     I assume that the number of nodes available is an integer 
!     multiple n of the number of replicas. Each replica then gets n
!     processors to do its energy calculation.
      num_ppr = num_proc / num_replica

      call mpi_comm_group(mpi_comm_world,  group_world, error)

!     The current version doesn't require a separate variable j. I
!     could just use i * num_ppr but this way it's more flexible.
      j = 0
      do i = 1, num_replica 
         ranks(i) = j 
         proc_range(1) = j
         proc_range(2) = j + num_ppr - 1
         proc_range(3) = 1
         call mpi_group_range_incl(group_world, 1, proc_range, group(i)
     &                              ,error)
         write (*,*) "Assigning rank ", j, proc_range, 
     &               "to group", group(i)
         call flush(6)
         j = j + num_ppr
      enddo

      do i = 1, num_replica
         call mpi_comm_create(mpi_comm_world, group(i), comm(i),error)
         if (comm(i).ne.MPI_COMM_NULL) then
             my_mpi_comm = comm(i)
             rep_id = i - 1
             write (*,*) rep_id, "has comm", my_mpi_comm
             call flush(6)
         endif
      enddo

!     Setup the communicator used for parallel tempering
      write (*,*) "PTGroup=", ranks(:num_replica)
      call flush(6)
      call mpi_group_incl(group_world, num_replica, ranks, group_partem,
     &                    error)
      call mpi_comm_create(mpi_comm_world, group_partem, partem_comm, 
     &                     error)

      if (partem_comm.ne.MPI_COMM_NULL) then
         write (*,*) partem_comm,myrank, "is master for ", rep_id, "."
      endif

      call mpi_comm_rank(my_mpi_comm,myrank,ierr)
      call mpi_comm_size(my_mpi_comm,no,ierr)

      write (*,*) "My new rank is ", myrank, "of", no
      call flush(6)
! = Done setting up communicators =====================================
      num_mol = 4

      if (newsta) then
         varfile = ' '
         do i = 1, num_mol - 1
            call init_molecule(iabin, grpn, grpc,in_fil,varfile)
         end do 
         varfile = 'abeta.var'
         call init_molecule(iabin, grpn, grpc,in_fil,varfile)
      else 
         filebase = "conf_0000.var"
         varfile = ' '
         do i = 1, num_mol - 1
            call init_molecule(iabin, grpn, grpc,in_fil,varfile)
         end do 
         call init_molecule(iabin, grpn, grpc,in_fil,
     &        fileNameMP(filebase, 6, 9, rep_id + 1))
      endif
      if (ientyp.eq.3) call init_abgn

!     ========================================  start of parallel tempering run
      write (*,*) "There are ", no,
     &            " processors available for ",rep_id
      call flush(6)
      do nml = 1, num_mol
         call distributeWorkLoad(no, nml)
      end do
      
!     READ  REFERENCE CONTACT MAP. 
!     Contacts needs a reference map, and we would like to get the
!     number of existing contacts. We just use a dummy map, since 
!     Abeta doesn't have a native structure.
      open(12, file = 'abeta.ref', status ="old")
      nresi=irsml2(nml)-irsml1(nml)+1
      write (*,*) 'Number of residues:', nresi
      do i=1,nresi
         read(12,*) (iref(i,j), j=1,nresi)
         write (*,*) iref(i,i)
!          do j = 1, nresi
!             read(12, *) iref(i,j)
!             write (*,*) i,j,iref(i,j)
!          end do
      end do
      nci = 0
      do i=1,nresi
         do j=nresi,i+3,-1
            if(iref(i,j).eq.1) nci = nci + 1
         end do
      end do
      
      call setvar(1,vlvr)

      do nml = 1, num_mol
         call setvar(nml,vlvr)
         e1 = enyshe(nml)
         write (*,*) 'Energy of molecule',nml, ' is ', e1,
     &         '(', eyvw, eyel, eyhb, eyvr, ')'
      end do
      
      ei = eninteract()
      write (*,*) ' Interaction energy is', ei
      et = energy()
      write (*,*) 'Solvent energy is', eysl
      write (*,*) 'Total energy is ', et
      
      call partem_p(num_replica, nequi, nswp, nmes, nsave, newsta,
     &              switch, rep_id, partem_comm)
!     ========================================  end of parallel tempering run
!     calculate CPU time using MPI_Wtime()
      endwtime = MPI_Wtime()
      if(myrank.eq.0.and.rep_id.eq.0) then

         write(*,*) "time for simulation using ", num_proc,
     &        " processors =",  endwtime - startwtime, " seconds"
         call flush(6)
      endif

!     ========================================  End of main
      CALL mpi_finalize(ierr)

      end


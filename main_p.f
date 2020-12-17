c     **************************************************************
c     
c     This file contains the   main (PARALLEL TEMPERING  JOBS ONLY,
C     FOR SINGULAR PROCESSOR JOBS USE main)
C     
C     This file contains also the subroutine: p_init_molecule
c     
c     Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c     Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c     
C     CALLS init_energy,p_init_molecule,partem_p
C     
c     **************************************************************
      program pmain

      include 'INCL.H'
      include 'INCP.H'
      include 'incl_lund.h'
      include 'mpif.h'

      character*80 libdir
      character*80 in_fil,ou_fil,filebase, varfile
      character*80 fileNameMP

      character grpn*4,grpc*4
      logical newsta

cc    Number of replicas 
      integer num_replica
cc    Number of processors per replica
      integer num_ppr
cc    Range of processor for crating communicators
      integer proc_range(3)
cc    Array of MPI groups
      integer group(MAX_REPLICA), group_partem
cc    Array of MPI communicators
      integer comm(MAX_REPLICA), partem_comm
cc    Array of nodes acting as masters for the energy calculation.
      integer ranks(MAX_REPLICA)
cc    Configuration switch
      integer switch
      integer rep_id
c     set number of replicas
      double precision eols(MAX_REPLICA)


      common/updstats/ncalls(5),nacalls(5)


c     MPI stuff, and random number generator initialisation

      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,myrank,ierr)
      call mpi_comm_size(mpi_comm_world,num_proc,ierr)

!       call VTSetup()
      enysolct = 0
      seed = 8368
      call sgrnd(seed)          ! Initialize the random number generator

c     =================================================== Energy setup
      libdir='SMMP/'      
c     Directory for SMMP libraries

c     The switch in the following line is now not used.
      flex=.false.              ! .true. for Flex  / .false. for ECEPP

c     Choose energy type with the following switch instead ...
      ientyp = 0
c     0  => ECEPP2 or ECEPP3 depending on the value of sh2
c     1  => FLEX 
c     2  => Lund force field
c     3  => ECEPP with Abagyan corrections
c     

      sh2=.false.               ! .true. for ECEPP/2; .false. for ECEPP3
      epsd=.false.              ! .true. for  distance-dependent epsilon

      itysol= 1                 !  0: vacuum
                                ! >0: numerical solvent energy
                                ! <0: analytical solvent energy & gradients
      isolscl=.false.
      tesgrd=.false.            ! .true. to check analytical gradients

      call init_energy(libdir)

c     calculate CPU time using MPI_Wtime()
      startwtime = MPI_Wtime()


c     ================================================= Structure setup
      grpn = 'nh2'              ! N-terminal group
      grpc = 'cooh'             ! C-terminal group

      iabin = 1                 ! =0: read from PDB-file
                                ! =1: ab Initio from sequence (& variables)

      in_fil='EXAMPLES/1bdd.seq'        ! Sequence file
      varfile = ' '

      newsta=.true.
      boxsize = 1000.0d0    ! Only relevant for multi-molecule systems
      num_replica = 1      ! Number of independent replicas. The file
                            !   temperatures must have at least as many
                            !   entries
      nequi=10              ! Number of MC sweeps before measurements 
                            !   and replica exchanges are started 
      nswp=500000           ! Number of sweeps
      nmes=10               ! Interval for measurements and replica exchange
      nsave=1000            ! Not used at the moment
                    
      switch = -1           ! How should the configuration be   
                            !   initialized?
                            ! -1 stretched chain 
                            !  0 don't do anything
                            !  1 initialize each angle to a random value
      
      ifrm=0
      ntlml = 0

c Decide if and when to use BGS, and initialize Lund data structures 
      bgsprob=0.6    ! Prob for BGS, given that it is possible
c upchswitch= 0 => No BGS 1 => BGS with probability bgsprob 
c 2 => temperature dependent choice 
      upchswitch=1
      rndord=.true.
      if (ientyp.eq.2) call init_lundff
c     =================================================================
c     Distribute nodes to parallel tempering tasks
c     I assume that the number of nodes available is an integer 
c     multiple n of the number of replicas. Each replica then gets n
c     processors to do its energy calculation.
      num_ppr = num_proc / num_replica

      call mpi_comm_group(mpi_comm_world,  group_world, error)

c     The current version doesn't require a separate variable j. I
c     could just use i * num_ppr but this way it's more flexible.
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

c     Setup the communicator used for parallel tempering
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

      if (newsta) then
         varfile = 'EXAMPLES/1bdd.var'
         call init_molecule(iabin, grpn, grpc,in_fil,varfile)
      else 
         filebase = "conf_0000.var"
         call init_molecule(iabin, grpn, grpc,in_fil,
     &        fileNameMP(filebase, 6, 9, rep_id + 1))
      endif
      if (ientyp.eq.3) call init_abgn

      nml = 1

c     RRRRRRRRRRMMMMMMMMMMMMSSSSSSSSSSDDDDDDDDDDDDD
      call rmsinit(nml,'EXAMPLES/1bdd.pdb')
c     RRRRRRRRRRMMMMMMMMMMMMSSSSSSSSSSDDDDDDDDDDDDD

!     READ  REFERENCE CONTACT MAP
      open(12, file = 'EXAMPLES/1bdd.ref', status ="old")
      nresi=irsml2(nml)-irsml1(nml)+1
      do i=1,nresi
         read(12,*) (iref(i,j), j=1,nresi)
      end do
      nci = 0
      do i=1,nresi
         do j=nresi,i+3,-1
            if(iref(i,j).eq.1) nci = nci + 1
         end do
      end do

c     ========================================  start of parallel tempering run
      write (*,*) "There are ", no,
     &            " processors available for ",rep_id
      call flush(6)
      nml = 1
      call distributeWorkLoad(no, nml)
      
      call partem_p(num_replica, nequi, nswp, nmes, nsave, newsta,
     &              switch, rep_id, partem_comm)
c     ========================================  end of parallel tempering run
c     calculate CPU time using MPI_Wtime()
      endwtime = MPI_Wtime()


      if(my_pt_rank.eq.0) then
         write(*,*) "time for simulation using ", num_proc,
     &        " processors =",  endwtime - startwtime, " seconds"
         call flush(6)
      endif

      print *,'update type, num calls, accepted calls '
      do i=1,5
         print *,i,ncalls(i),nacalls(i)
      enddo

c     ========================================  End of main
      CALL mpi_finalize(ierr)

      end


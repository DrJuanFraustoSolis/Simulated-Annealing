!**************************************************************     
! Parallel tempering of Met-Enaphalin on a single CPU
!
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku Hu
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************
      
      program main

      include '../INCL.H'
      include '../INCP.H'
      common/updstats/ncalls(5),nacalls(5)
      character*80 libdir, seqfile, varfile
      character grpn*4,grpc*4
      logical bgsposs,newsta
      integer switch

c =================================================== Energy setup

c            Directory for SMMP libraries
c     Change the following directory path to where you want to put SMMP
c     libraries of residues. 
      libdir='../SMMP/'

c      The switch in the following line is now not used.
      flex=.false.        ! .true. for Flex  / .false. for ECEPP

c     Choose energy type with the following switch instead ...
      ientyp = 0
c        0  => ECEPP2 or ECEPP3 depending on the value of sh2
c        1  => FLEX 
c        2  => Lund force field
c        3  => ECEPP with Abagyan corrections
c

      sh2=.false.         ! .true. for ECEPP/2; .false. for ECEPP3
      epsd=.false.        ! .true. for  distance-dependent  dielectric
                          !  permittivity

      itysol= 1    !  0: vacuum
                   ! >0: numerical solvent energy
                   ! <0: analytical solvent energy & gradients

      call init_energy(libdir)

c ================================================= Structure setup

      grpn = 'nh2' ! N-terminal group
      grpc = 'cooh'! C-terminal group

      iabin = 1  ! =0: read from PDB-file
                 ! =1: ab Initio from sequence (& variables)
      seqfile='enkefa.seq'
      varfile='enkefa.var'
      ntlml = 0
      write (*,*) 'Solvent: ', itysol
c     Initialize random number generator.
      call sgrnd(31433)
      
      if (itysol.eq.0.and.ientyp.eq.3) then
         print *,'Can not use Abagyan entropic corrections without '
         print *,'solvent term. '
         stop
      endif

      call init_molecule(iabin,grpn,grpc,seqfile,varfile)
c Decide if and when to use BGS, and initialize Lund data structures 
      bgsprob=0.75   ! Prob for BGS, given that it is possible
c upchswitch= 0 => No BGS 1 => BGS with probability bgsprob 
c 2 => temperature dependent choice 
      upchswitch=1
      rndord=.true.
      call init_lund
      if (ientyp.eq.2) call init_lundff
      if (ientyp.eq.3) call init_abgn
      

c ========================================  Add your task down here
      num_rep = 5
      nequi = 100
      nsweep = 10000
      nmes = 10
      newsta = .true.
      switch = 1
c     parallel tempering on a single CPU
      eol = energy()
      write (*,*) "Energy before randomization:", eol
      call partem_s(num_rep, nequi, nsweep, nmes, nmes, newsta, switch)
      eol = energy()
      write (*,*) "Final energy:", eol

c ========================================  End of main      
       end

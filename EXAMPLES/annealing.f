!**************************************************************     
! Minimize the vacuum energy of Met-Enkaphalin using SMMP
!
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku Hu
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************

      program main
      
      include "../INCL.H"
      
!! Storage space for force-field library paths, seuence-file name, and variable-file name
      character*80 libdir, seqfile, varfile
!! Storage space for end groups
      character grpn*4,grpc*4
      integer nequi, nsweeps, nmes
      double precision tmax, tmin
      logical lrand
      
! =================================================== Energy setup

!     Directory for SMMP libraries
!     Change the following directory path to where you want to put SMMP
!     libraries of residues. 
      libdir='../SMMP/'

!!     Choose energy type with the following switch
!        0  => ECEPP2 or ECEPP3 depending on the value of sh2
!        1  => FLEX 
!        2  => Lund force field
!        3  => ECEPP with Abagyan corrections
!
      ientyp = 0

      sh2=.true.         ! .true. for ECEPP/2; .false. for ECEPP3
      epsd=.false.        ! .true. for  distance-dependent  dielectric
                          !  permittivity

      itysol= 0    !  0: vacuum
                   ! >0: numerical solvent energy
                   ! <0: analytical solvent energy & gradients

      call init_energy(libdir)

! ================================================= Structure setup

      grpn = 'nh2' ! N-terminal group
      grpc = 'cooh'! C-terminal group

      iabin = 1  ! =0: read from PDB-file
                 ! =1: ab Initio from sequence (& variables)
      seqfile='enkefa.seq'
!      varfile='EXAMPLES/enkefa.var'
      varfile = 'enkefa.ann'
      ntlml = 0
      
      call init_molecule(iabin,grpn,grpc,seqfile,varfile)

      seed = 81236
      call sgrnd(seed)
      upchswitch= 0 !Turn off biased Gaussian step
      nequi=100
      nsweeps=100000
      nmes=1000
      tmax=1000
      tmin=100
      lrand=.true. 
      call anneal(nequi, nsweeps, nmes, tmax, tmin, lrand)

      end program main
      

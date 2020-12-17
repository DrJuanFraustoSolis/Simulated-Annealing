c Simulated annealing (SA)
c Copyright (C) 2020  Dr. Juan Paulo Sánchez Hernández, and Dr. Juan Frausto Solis
c Copyright (C) 2005 Frank Eisenmenger, U.H.E. Hansmann, Shura Hayryan, Chin-Ku Hu

c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or (at
c your option) any later version.
c 
c This program is distributed in the hope that it will be useful, but
c WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c General Public License for more details.
c 
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
c USA.

c **************************************************************
c This file contains the:  main (SINGLE PROCESSOR JOBS ONLY,
C                                FOR PARALLEL JOBS USE pmain)
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
C CALLS: init_energy,init_molecule
C CALLS TASK SUBROUTINE: anneal,canon,elp,minim,mulcan_par,
c                        mulcan_sim,partem_s, or regul
C CAN ALSO CALL MEASUREMENT ROUTINES: cnteny,contacts,helix,hbond,
C                                    outpdb,outvar,rgyr,
C                                    rmsinit and rsmdfun,zimmer
c $Id: main.f 334 2007-08-07 09:23:59Z meinke $
c **************************************************************
      
      program main

      include 'INCL.H'
      include 'INCP.H'

      common/updstats/ncalls(5),nacalls(5)
      character*80 libdir, seqfile, varfile
      character grpn*4,grpc*4
      character(8) x1
      logical lrand,bgsposs
	  real initial,final,namefile1
 

c      character*10 b(3)
c      integer date_time(8)
c      character*50 filename
c      real tiempo
c	  integer time

c =================================================== Energy setup

c            Directory for SMMP libraries
c     Change the following directory path to where you want to put SMMP
c     libraries of residues. 
      libdir='./SMMP/'

c      The switch in the following line is now not used.
      flex=.false.        ! .true. for Flex  / .false. for ECEPP

c     Choose energy type with the following switch instead ...
      ientyp = 0
c        0  => ECEPP2 or ECEPP3 depending on the value of sh2
c        1  => FLEX 
c        2  => Lund force field
c        3  => ECEPP with Abagyan corrections
c

      sh2=.true.         ! .true. for ECEPP/2; .false. for ECEPP3
      epsd=.false.        ! .true. for  distance-dependent  dielectric
                          !  permittivity

      itysol= 0    !  0: vacuum
                   ! >0: numerical solvent energy
                   ! <0: analytical solvent energy & gradients

      call init_energy(libdir)

c ================================================= Structure setup

      grpn = 'nh2' ! N-terminal group
      grpc = 'cooh'! C-terminal group

      iabin = 1  ! =0: read from PDB-file
                 ! =1: ab Initio from sequence (& variables)
      seqfile='EXAMPLES/1in3.seq'
      varfile='EXAMPLES/1in3.var'
!       varfile = ' '
      
      ntlml = 0
      write (*,*) 'Solvent: ', itysol

c     Initialize random number generator.
	CALL SYSTEM_CLOCK(seed)
c	call sgrnd(seed)
c	call srand(seed)

c     Initialize random number generator.
c      call sgrnd(31433)
      
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
 
		call cpu_time(initial)
			call sa(g_energy)
		call cpu_time(final)
      			total = final - initial
	
	namefile1=g_energy*10000
	write(x1,'(I8)')  int(namefile1)
    
      	write(*,*) 'ENERGY FOUND',g_energy,total
	open(16, file='./RESULTS/ENERGY'//x1//'.txt',status='new')
	write(16,*) 'ENERGY FOUND',g_energy,total
	close (16)

c ========================================  End of main      
       end

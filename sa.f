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


c***********************SA************************************
c
c -This file contains a simulated annealing (SA)
c -The SA has two fundamental parts: 1) The cooling scheme with a start and final temperature, 
c  and 2) The convergence criteria
c 
c -The firts version was developed by Dr. Juan Paulo Sánchez Hernández, 
c and Dr. Juan Frausto Solis
c  2015-2020
c***************************************************************


      subroutine  sa(g_energy)

C --------------------------------------------------------------
C PURPOSE: SIMULATED ANNEALING SEARCH OF LOWEST-POTENTIAL-ENERGY
C          CONFORMATIONS OF PROTEINS
C
C CALLS: addang,energy,metropolis,outvar,outpdb,rgyr,setvar,zimmer
C
C ---------------------------------------------------------------
C ---------------------------------------------------------------

C------------------------------------------
C---------------------SA-----------------
C------------------------------------------

      include 'INCL.H'

      character(8) x2
      logical lrand
      parameter(lrand=.true.)

        dimension vlvrm(mxvr)
	dimension idvrm(mxvr)
	dimension axvrm(mxvr)

        dimension vlvr_aux(mxvr)
	dimension idvr_aux(mxvr)
	dimension axvr_aux(mxvr)

	real*8 e,bangl
c	real inicio,final,total
	real initT,namefile
c	integer seed1,p,THREADS,ID,MaxThreads
 	
	s_window = 90.0d0

	s_slope = 0.000000000001d0
	slope = 0
	C = 0
	D = 0
	cnt = 1
	cnt2 = 0

	initT = secnds(0.0)
       nresi=irsml2(1)-irsml1(1) + 1     
c ==================
c == random start ==    !!!random start when no template  
c ==================

c       if(lrand) then
c        	do i=1,nvr
c         	iv=idvr(i)  ! provides index of non-fixed variable
c        	 e = rand(int(initT))
c        	 dv=axvr(iv)*e
c	         vr=addang(pi,dv)
c          	vlvr(iv)=vr
c        	enddo
c
c  	  end if	

c ====================================================	
        eol = energy()
	ymin = eol
 	epsilon = eol

c =======================================================
c == Write start configuration in pdb-format into file ==
c =======================================================
c      do i=1,ntlml
c        call outpdb(i,11)
c      end do

c==========================
c==	Simulated Annealing  ==
c==========================
      write(*,*) 'SIMULATED ANNEALING'
	call cpu_time(start)

c ====================================================
c == Initial parameters for Simulated Annealing ==
c ====================================================
	
c ====================================================
c ========== Parameters of 1in3 ==============
c ====================================================	          

c====== This parameters were tunning by J. Frausto-Solís, H. Sanvicente-Sánchez, and F. Imperial-Valenzuela, "ANDYMARK: An Analytical Method to Establish Dynamically the Length of the Markov Chain in Simulated Annealing for the Satisfiability Problem," in Simulated Evolution and Learning, Springer Berlin Heidelberg, 2006, pp. 269–276.
c============================================================================================


	currtem = 44778080443048300000000000000000000000.0d0 !Initial Temperature 1in3
	stop_temp = 0.0000000200604823484395 !Final Temperature 1in3  

        alpha = 0.95
	blmax = 360.0d0
	bbeta = 1.00075046918563 !
	
c =================================
c == Start the Simulated Annealing == 
c =================================
      do while(currtem.GE.stop_temp)
        propon = 0.0d0
        acepta = 0.0d0
        ycurr = 0.0d0
c =====================================
c == Start Metropolis ==
c == Crescent Markov Chain ==
c =====================================
        do while(propon.LE.NINT(blmax))
          propon = propon + 1.0d0
          ycurr = eol

          call metropolis(eol,currtem,acepta)
c =================================================
c == Store the local lowest-energy conformation  ==
c =================================================
          if (eol.LT.ycurr) then
            ycurr = eol
          end if
c =================================================
c == Store the global lowest-energy conformation ==
c =================================================
          if (eol.LT.ymin) then
            ymin = eol
            write(*,*) 'MINIMA:', currtem, ymin
c ==================================
c == backup of Minimum point   ==
c == idvrm, axvrm, vlvrm  ==
c ==================================
       	do i=1,nvr
		  idvrm(i)  = idvr(i)
		  iv        = idvr(i)
		  axvrm(iv) = axvr(iv)
		  vlvrm(iv) = vlvr(iv)
       	enddo
          end if
c==================================================
c   The convergence criteria for dynamic equilibrium
c==================================================

			if ((currtem.LT.s_window).AND. (cnt.LT.3)) then
				D = D + ymin
				C = C + cnt * ymin
				cnt = cnt + 1
			end if

			if (cnt.EQ.3) then
		        slope = (((12 * C)-(6*(cnt-1)*D))/(cnt**3 - cnt))	
			cnt = 1	

			 if (slope.LT.s_slope) then
c----------------------------Overheating strategy-------------------------------------------
				currtem = currtem*0.618
				bangl= bangl+1
		         endif		
			
			  if (bangl.EQ.2) then
				currtem = stop_temp
			  endif

			endif

        enddo
c ========================
c == End of Metropolis ==
c ========================

c ==============================
c == Lower the Temperature ==
c ==============================
        currtem = alpha * currtem
	write(*,*) currtem,ymin
c ==============================
c == Crescent Markov Chain ==
c ==============================
        blmax = bbeta * blmax

      enddo
	g_energy = ymin
c ====================================
c == End of the Simulated Annealing ==
c ====================================
	  namefile=g_energy*10000
	  write(x2,'(I8)') int(namefile) 	
          open(16, file='./RESULTS/Estruc'//x2//'.pdb',status='new')	  
c =======================================================
c == Write start configuration in pdb-format into file ==
c =======================================================
       do i=1,ntlml
         call outpdb(i,16)
       end do

       close(16)

      end

c **************************************************************
c
c This file contains the subroutines: contacts,c_alfa,c_cont
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************
      
       subroutine contacts(ncn,nham2,dham)

c ..............................................................
c
c CALCULATES NUMBER OF CONTACTS IN GIVEN CONFORMATION, NUMBER OF
c CONTACTS WHICH ARE THE SAME IN GIVEN AND REFERENCE ONFORMATION,
c AND THE HAMMING DISTANCE BETWEEN GIVEN  CONFORMATION AND THE 
c REFERENCE CONFORMATIONa
c
c CALLS: c_cont
c ..............................................................

      include 'INCL.H'

!f2py integer, intent(out) :: ncn, nham2 
!f2py double precistion, intent(out) :: dham

      call c_cont(1,ncode)

      ncn=0
      nham=0
      nham2=0
      nresi=irsml2(1)-irsml1(1)+1
      do i=1,nresi
        do j=nresi,i+3,-1

        if (ijcont(i,j).eq.1) then
          ncn=ncn+1
          if (iref(i,j).eq.1) nham2 = nham2+1
        end if

        nham = nham + abs(ijcont(i,j)-iref(i,j))
       end do
      end do

      if (ncn.ne.0.and.nci.ne.0) then
        dham = float(nham)/float(ncn)/float(nci)
      else
        dham = 1.0
      end if

      return
      end 


c ********************************* 
      subroutine c_alfa(nmol,ncode)

c ......................................................
c    Calculates the indices of C-alpha atoms and 
c    stores in the array ind_alf(mxrs)
c                        
c    Usage: call c_alfa(nmol,ncode)
c
c           nmol - index of the molecule
c           ncode ---> not in use in the current version
c
c    OUTPUT:  ind_alf(mxrs)
c
c CALLS: none
c ......................................................

      include 'INCL.H'
              
      do n_res=irsml1(nmol),irsml2(nmol) ! Over res. 
        do ia=iatrs1(n_res),iatrs2(n_res) ! Over the atoms of res. 

c     Check for C_alpha atoms

         if (nmat(ia)(1:2).eq.'ca') then
           ind_alf(n_res)=ia
         endif

       enddo ! Over the atoms of res. 
      enddo ! Over the res. 
 
      return
      end

c **********************************
      subroutine c_cont (nmol,ncode)

c..............................................................
c  Calculates the matrix of contacts between aminoacid residues 
c  of the molecule "nmol" according to  L.Mirny and E.Domany, 
c  PROTEINS:Structure, Function, and Genetics 26:391-410 (1996)
c                
c  Two residues are in contact if their C_alpha atoms are
c  closer than 8.5 Angstrem
c
c  Usage: call c_cont(nmol,ncode)
c
c       Where nmol is the index of the molecule (always 1, in the 
c       current version of SMM)
c       ncode ---> not in use in the current version
c
c  IMPORTANT: Before the first call of this subroutine  "c_alfa"
c          must be called to calculate the inices of C_alpha atoms.
c          (ONLY ONCE)
c
c   OUTPUT: The output of this routine is the contact matrix
c          ijcont(mxrs,mxrs) 
c
c              ijcont(i,j)=0---> residues i and j are not in contact
c              ijcont(i,j)=1---> ---------''----- are in contact
c              ijcont(i,j)=2---> residues i and j are adjacent
c
c    NOTE:  Adjacent residues are always in contact (and therefore not
c           counted)
c
c         Here "mxrs" is the maximum number of residues for SMM
c         Obviously, this subroutine calculates only NxN part
c         of that matrix, N -is the number of res. in "nmol"
c 
c CALLS:  none
c..............................................................

       include 'INCL.H'

       rcut=8.5   ! Domany
              
              
       do nr_i=irsml1(nmol),irsml2(nmol) ! Over res. i

          ijcont(nr_i,nr_i)=2
          if(nr_i+1.le.irsml2(nmol)) then
              ijcont(nr_i,nr_i+1)=2
              ijcont(nr_i+1,nr_i)=2
              if(nr_i+2.le.irsml2(nmol)) then
                 ijcont(nr_i,nr_i+2)=2
                 ijcont(nr_i+2,nr_i)=2
               end if
           end if

           do nr_j=nr_i+3,irsml2(nmol) ! Over res. j 

c             write(*,'(2i3)'),nr_i,nr_j

              ic=0

              ialf=ind_alf(nr_i)
              jalf=ind_alf(nr_j)

              rij2=(xat(ialf)-xat(jalf))**2
     #              +(yat(ialf)-yat(jalf))**2
     #                   + (zat(ialf)-zat(jalf))**2
              if(sqrt(rij2).lt.rcut) ic=1

c             write(*,'(2i3)'),nr_i,nr_j

              ijcont(nr_i,nr_j)=ic
              ijcont(nr_j,nr_i)=ic ! The matrix is symmetrical

           end do ! Over res. j
       end do ! Over res. i

       return
       end


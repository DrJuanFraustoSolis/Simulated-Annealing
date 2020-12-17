c **************************************************************
c
c This file contains the subroutines: init_energy,setpar
C This file contains a BLOCK DATA statement
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************


      subroutine init_energy(libdir)

c ----------------------------------------------
c PURPOSE: initialize energy parameters
c        0  => ECEPP2 or ECEPP3 depending on the value of sh2
c        1  => FLEX 
c        2  => Lund force field
c        3  => ECEPP with Abagyan corrections
c
c
c CALLS:   setpar, tessel,iendst
c
c contains: BLOCK DATA
c ----------------------------------------------

      include 'INCL.H'

      character libdir*(*),tesfil*80
      
      if (ientyp.eq.1) then
          flex = .true.
      else 
          flex = .false.
      end if
      
      lunlib=10
      ll=iendst(libdir)

      if (flex) then                    !  Flex

        reslib=libdir(1:ll)//'lib.flx'
        lunchg=12
        chgfil=libdir(1:ll)//'charges'

      else

        if (sh2) then                   !  Scheraga (ECEPP/2)
          reslib=libdir(1:ll)//'lib.sh2'
        else                            !  Scheraga (ECEPP/3)
          reslib=libdir(1:ll)//'lib.sh3'
        endif

      endif

      call setpar()  ! Initialize force field parameters


C----Initialize solvation part if necessary
      write (*,*) 'init_energy: itysol = ',itysol
      write(*,*) 'init_energy: esol_scaling = ',isolscl

      its = iabs(itysol)

      if (its.gt.0.and.its.le.mxtysol) then

        ll=iendst(libdir)
        tesfil = libdir(1:ll)//'tes.dat'

        open(unit=20,file=tesfil,status='old',err=10)

        call tessel()

        close(20)

      else
 
        if (itysol.ne.0) then
          write(*,'(a)') ' init_energy>  undefined solvent type !'
          stop
        endif

      endif

c ___________________________ initialise COMMON 'con_r'
      idloa=ichar('a')
      idloz=ichar('z')
      idupa=ichar('A')
      idupz=ichar('Z')

      return

 10   write (*, '(a)') 'Cannot open the file with surface points'
      stop

      end

c *********************
      subroutine setpar

c __________________________________________________________
c PURPOSE: initialize parameter set for empirical potentials
c          depending on variable 'flex'
c
c CALLS:   None
c __________________________________________________________

      include 'INCL.H'

      dimension hbc(mxhbdo,mxhbac),hba(mxhbdo,mxhbac)
      logical do(mxtyat),ac(mxtyat)


      tesgrd = .false.  ! numerical check of analytical gradients

c ______________________________________ Lennard-Jones parameters
      if (flex) then

        plt=plt_f
        slp=slp_f
        cohb=cohb_f

        do i=1,mxtyat
          do j=1,i
            cij(i,j)=c_f(i,j)
            cij(j,i)=cij(i,j)
            aij(i,j)=a_f(i,j)*1000
            aij(j,i)=aij(i,j)
          enddo
          do j=1,mxtyat
            a14(i,j)=aij(i,j)
          enddo

          do(i)=do_f(i)
          ac(i)=ac_f(i)
        enddo

        do i=1,mxhbdo
          do j=1,mxhbac
            hbc(i,j)=chb_f(i,j)
            hba(i,j)=ahb_f(i,j)
          enddo
        enddo

        do i=1,mxtyto
          e0to(i)=e0to_f(i)/2.
          sgto(i)=sgto_f(i)
          rnto(i)=rnto_f(i)
        enddo

      else  ! ---------------------- Scheraga:

        conv=conv/eps_s

        do i=1,mxtyat
          atpl(i)=atpl(i)/100
          efel(i)=efel(i)/100
          emin(i)=emin(i)/1000
          rmin(i)=rmin(i)/100
        enddo

        do i=1,mxtyat
          ri=rmin(i)
          ai=atpl(i)
          aei=sqrt(ai/efel(i))
cc          aic=ai/ehm             !!  ICM
cc          do j=i,mxtyat          !!  -"-
          aic=ai*ehm                 !!  comment for ICM:
          cij(i,i)=aic*ai/(aei+aei)  !!        -"-
          a=.5*cij(i,i)*ri**6
          aij(i,i)=a                 !!
          a14(i,i)=.5*a              !!
          do j=i+1,mxtyat            !!
            aj=atpl(j)
c _______ Constant for 6-12 attractive term (Slater-Kirkwood formula)
            c=aic*aj/(aei+sqrt(aj/efel(j)))
            cij(i,j)=c
            cij(j,i)=c
c ____________________________ repulsive term (form. 3 & 6 of ref 2)
            rij=.5*(ri+rmin(j))
            a=.5*c*rij**6
            aij(i,j)=a
            aij(j,i)=a
            a=.5*a
            a14(i,j)=a
            a14(j,i)=a
          enddo
          do(i)=do_s(i)
          ac(i)=ac_s(i)
        enddo

c +++++++++++++++++++++++++++++++++
        cij(1,1)=45.5d0
        aij(1,1)=14090.0d0
        cij(2,2)=45.5d0
        aij(2,2)=14380.0d0
        cij(3,3)=45.5d0
        aij(3,3)=8420.0d0
        cij(4,4)=45.5d0
        aij(4,4)=8420.0d0
        cij(5,5)=45.5d0
        aij(5,5)=11680.0d0
        cij(6,6)=45.5d0
        aij(6,6)=14380.0d0
        cij(7,7)=370.5d0
        aij(7,7)=906100.0d0
        cij(8,8)=766.6d0
        aij(8,8)=1049000.0d0
        cij(9,9)=509.5d0
        aij(9,9)=653600.0d0
        cij(10,10)=217.2d0
        aij(10,10)=125600.0d0
        cij(11,11)=369.0d0
        aij(11,11)=170200.0d0
        cij(12,12)=217.2d0
        aij(12,12)=125600.0d0
        cij(13,13)=401.3d0
        aij(13,13)=375200.0d0
        cij(14,14)=401.3d0
        aij(14,14)=375200.0d0
        cij(15,15)=401.3d0
        aij(15,15)=375200.0d0
        cij(16,16)=2274.4d0
        aij(16,16)=5809000.0d0
        cij(17,17)=45.5d0
        aij(17,17)=5340.0d0
        cij(18,18)=370.5d0
        aij(18,18)=909000.0d0
c +++++++++++++++++++++++++++++++++
        do i=1,mxtyat
          a14(i,i)=.5*aij(i,i)
        enddo
c +++++++++++++++++++++++++++++++++

        do i=1,mxtyat
c	  write( *, '(18f14.6)' )  ( a14(i,j), j = 1, mxtyat )
        enddo

        do i=1,mxhbdo
          do j=1,mxhbac
            hbc(i,j)=chb_s(i,j)
            hba(i,j)=ahb_s(i,j)
          enddo
        enddo

        do i=1,mxtyto
          e0to(i)=e0to_s(i)/2.
          sgto(i)=sgto_s(i)
          rnto(i)=rnto_s(i)
        enddo


      endif
c -------------------------------------------- Hydrogen Bond Parameters
      do i=1,mxtyat
        do j=1,mxtyat
          ihbty(i,j)=0
          chb(i,j)=0.
          ahb(i,j)=0.
        enddo
      enddo

      iac=0
      ido=0
      do i=1,mxtyat
        if (do(i)) then
          ido=ido+1
          jac=0
          do j=1,i-1
            if (ac(j)) then
              jac=jac+1
              ihbty(i,j)=1
              ihbty(j,i)=-1
              chb(i,j)=hbc(ido,jac)
              chb(j,i)=chb(i,j)
              ahb(i,j)=hba(ido,jac)
              ahb(j,i)=ahb(i,j)
            endif
          enddo
        elseif (ac(i)) then
          iac=iac+1
          jdo=0
          do j=1,i-1
            if (do(j)) then
              jdo=jdo+1
              ihbty(i,j)=-1
              ihbty(j,i)=1
              chb(i,j)=hbc(jdo,iac)
              chb(j,i)=chb(i,j)
              ahb(i,j)=hba(jdo,iac)
              ahb(j,i)=ahb(i,j)
            endif
          enddo
        endif
      enddo

      do i=1,mxtyto  ! eases calculation of torsional derivatives
        esnto(i)=e0to(i)*sgto(i)*rnto(i)
      enddo

      return
      end
c **************
      BLOCK DATA

      include 'INCL.H'

c  Atom types ------------------------------------------------------------
c                                    Original types  -Scheraga:  -Flex:
c  H  1 - with aliphatic carbon                               1      12
c     2 - with aromatic carbon                                3      13
c     3 - with non-sp3 types of nitrogen                      2       1
c     4 - with sp3-hybr. nitrogen                             2       2
c     5 - with oxygen                                         4       1
c     6 - with sulfur                                         3(was 5)1
c  C  7 - sp3-hybr. carbon                                    6,9     3
c     8 - sp2-carbon (carbonyl,carboxyl,carboxylate)          7,11    4
c     9 - aromatic carbon                                     8,10    4
c  O 10 - hydroxyl, ester oxygen (inc. water)                 18,19   8
c    11 - carbonyl oxygen                                     17      9
c    12 - carboxylate oxygen                                  18,19  10
c  N 13 - aliph. nitrogen with 0/1 hydrogen & charged N       13-15   6
c    14 - nitrogen with two hydrogens                         13-15   5
c    15 - all other nitrogens (+ sp2-hybrid. in heteroc.)     13-15   7
c  S 16 - any sulfur                                          20,21  18,19
c  H 17 - H-delta of Pro, Hyp of ECEPP/3 dataset               5(new) -
c  C 18 - C-delta of Pro, Hyp of ECEPP/3 dataset              12(new) -

c  Classes for torsional potential ---------------------------------------
c
c   1 : 'Omega' = C'(pept.)-N(pept.)  [Cpept-Npept]
c   2 : 'Phi'   = N(pept.)-C(sp3)     [C4-Npept]
c   3 : 'Psi'   = C(sp3)-C'(pept.)    [C4-Cpept]
c   4 : 'Chi1'  = C(sp3)-C(sp3)       [C4-C4]
c   5 : C(sp3)-OH (Hydroxyl)          [C4-OH]
c   6 : C(sp3)-NH2                    [C4-NH2]
c   7 : C(sp3)-NH3+                   [C4-NH3+]
c   8 : C(sp3)-NH-(guanidyl)          [C4-NHX]
c   9 : C(sp3)-COOH(carboxyl)         [C4-COO]
c  10 : C(sp3)-COO-(carboxylate)      [C4-COO]
c  11 : C(sp3)-CO(sp2 of amide)       [C4-Cpept]
c  12 : C(sp3)-C(aromatic ring)       [C4-C3]
c  13 : C(sp3)-S                      [C4-SC4]
c  14 : C(sp3)-SH                     [C4-SH]
c  15 : C(aromatic ring)-OH           [C3-OH]
c ________________________________________________ "rigid" torsions:
c  16 : C(carboxyl)-OH                [C3-OH]
c  17 : -NH-C(sp2 of guanidyl)        [C3-NHX]
c  18 : -C(sp3)-NH2 (guanidyl)        [not in Flex]
c  19 : -C(sp3)-NH2 (amide)           [Cpept-Npept]

      data conv/332.d0/  ! to convert electrost. energy into [kcal/mole]

c ------------------------- ECEPP/3 potential --------------------------------
c 1) Momany F.A McGuire R.F Burgess A.W Scheraga H.A J Phys Chem v79 2361-2381
c    1975
c 2) Nemethy G Pottle M.S Scheraga H.A, J Phys Chem v87 1883-1887 1983
c 3) Sippl M.J Nemethy G Scheraga H.A J Phys Chem v88 6231-6233 1984
c 4) Nemethy G Gibson K.D Palmer K.A Yoon C.N Paterlini G Zagari A Rumsey S
c    Scheraga H.A J Phys Chem v96 6472-6484 1992
c ----------------------------------------------------------------------------

      data eps_s/2.d0/  ! Distance-INdependent diel. constant
c     data eps_s/6.d0/  ! Distance-INdependent diel. constant
      data plt/78.d0/,  slp/0.3d0/   ! Parameters for Epsilon(R)

      data ehm /362.55d0/  !  Angstrom**2/3 * kcal / mol  ! from KONF90
cc      data ehm /362.09561409d0/  !  Angstrom**2/3 * kcal / mol
c From:
c   1.5
c * elementary charge       = 4.80325   *e+2  Angstrom**3/2 * g**1/2 * s**(-1)
c * Planck's constant/2*Pi  = 1.0545887 *e-34 Joule * s
c * Avogadro's number       = 6.022045  *e+23 mol**(-1)
c / sqrt (mass of electron) = sqrt (9.109534  *e-28 g )
c / thermal equivalent      = 4.1868    *e+3  Joule * kcal**(-1)
cc      data ehm /362.36d0/         ! calculated using Tab II in ref. 2
cc      data 1/ehm /2.757670d-3/    ! 3*sqrt(m)/(2*e*h) taken from ICM

c ---------------------- atomic polarizabilties (*100,[Angstrom**3])
c                1   2   3   4   5   6   7    8    9  10  11  12
      data atpl/42.,42.,42.,42.,42.,42.,93.,151.,115.,59.,84.,59.,
c               13  14  15   16  17  18
     #          93.,93.,93.,220.,42.,93./
c ---------------------- effective numbers of electrons (*100,ref. 2)
c                1   2   3   4   5   6    7    8    9   10   11   12
      data efel/85.,85.,85.,85.,85.,85.,520.,520.,520.,700.,700.,700.,
c                13   14   15    16  17   18
     #          610.,610.,610.,1480.,85.,520./
c ------------------------- min. pairwise 6-12 energy (*1000,[kcal/mol])
c                1   2   3   4   5   6   7    8   9  10   11  12
      data emin/37.,36.,61.,61.,44.,36.,38.,140.,99.,94.,200.,94.,
c                13   14   15   16  17  18
     #          107.,107.,107.,223.,99.,38./
c ---------------------------- opt. pairwise distance (*100,[Angstrom])
c                 1    2    3    4    5    6    7    8    9   10   11
      data rmin/292.,293.,268.,268.,283.,293.,412.,374.,370.,324.,312.,
c                12   13   14   15   16   17   18
     #          324.,351.,351.,351.,415.,248.,412./
c ---------------------------------------------- Hydrogen-bond donors
c                  1       2       3      4      5      6
      data do_s/.false.,.false.,.true.,.true.,.true.,.false.,
c                  7       8       9      10      11      12
     #          .false.,.false.,.false.,.false.,.false.,.false.,
c                 13      14      15      16      17      18
     #          .false.,.false.,.false.,.false.,.false.,.false./
c -------------------------------------------- Hydrogen-bond acceptors
c                  1       2       3       4       5       6
      data ac_s/.false.,.false.,.false.,.false.,.false.,.false.,
c                  7       8       9      10     11     12
     #          .false.,.false.,.false.,.true.,.true.,.true.,
c                 13     14     15     16      17      18
     #          .true.,.true.,.true.,.false.,.false.,.false./
cc     #        .false.,.true.,.true.,.false.,.false.,.false./ !! ICM
c --------------------------------- HB-parameters (/1000,attraction)
      data chb_s/2624.,2624.,4610.,.0,  ! given as:
     #           4014.,4014.,5783.,.0,  !  (ac_typ x do_typ)
     #           2624.,2624.,4610.,.0,  ! to be used:
     #           8244.,8244.,8244.,.0,  !  (DO_typ x AC_typ)
     #           8244.,8244.,8244.,.0,  ! i.e.:
     #           8244.,8244.,8244.,.0/  !  ( 3-5 x 10-15 )
c --------------------------------- HB-parameters (/1000,repulsion)
      data ahb_s/ 5890., 5890.,11220.,.0,
     #           12040.,12040.,16583.,.0,  ! 13344 -> 16583 = Ref. 3
     #            5890., 5890.,11220.,.0,
     #           32897.,32897.,32897.,.0,
     #           32897.,32897.,32897.,.0,
     #           32897.,32897.,32897.,.0/

c                   1  2  3   4   5   6   7  8  9 10 11 12 13  14  15
      data e0to_s /20.,0.,0.,2.7,.6,1.8,1.8,0.,0.,0.,0.,0.,2.,1.5,3.5
c                   16 17  18  19
     #             ,8.,18.,20.,15./
      data sgto_s /-1.,0.,0., 1.,1., 1., 1.,0.,0.,0.,0.,0.,1., 1.,-1.
     #            ,-1.,-1.,-1.,-1./
      data rnto_s / 2.,0.,0., 3.,3., 3., 3.,0.,0.,0.,0.,0.,3., 3., 2.
     #             ,2., 2., 2., 2./

c ---------------------------- Flex potential ----------------------------------
c Lavery R Sklenar H Zakrzewska K Pullman B J Biomol Struct Dyn v3 989-1014 1986
c VdW-parameters from: Zhurkin V.B Poltiev V.I Florent'ev V.L Molekulyarnaya
c                      Biologiya v14 116 1980
c ------------------------------------------------------------------------------

      data plt_f/78.d0/,  slp_f/0.16d0/   ! Parameters for Epsilon(R)
      data cohb_f/6.d0/  ! Cut-off distance betw. H- & acceptor atom for HB
      data c_f/  ! ----------- Lennard-Jones C6-parameters (attraction)
     #40.,40.,40.,40.,40.,40.,100.,126.,126.,86.,121.,105.,126.,105.,
     #146.,213.1,3*0.,40.,40.,40.,40.,40.,100.,126.,126.,86.,121.,105.,
     #126.,105.,146.,213.1,4*0.,40.,40.,40.,40.,100.,126.,126.,86.,121.,
     #105.,126.,105.,146.,213.1,5*0.,40.,40.,40.,100.,126.,126.,86.,121.
     #,105.,126.,105.,146.,213.1,6*0.,40.,40.,100.,126.,126.,86.,121.,
     #105.,126.,105.,146.,213.1,7*0.,40.,100.,126.,126.,86.,121.,105.,
     #126.,105.,146.,213.1,8*0.,250.,316.,316.,217.,305.,264.,316.,264.,
     #367.,489.,9*0.,400.,400.,274.,385.,334.,400.,334.,464.,537.4,10*0.
     #,400.,274.,385.,334.,400.,334.,464.,537.4,11*0.,200.,283.,245.,
     #278.,233.,330.,424.,12*0.,400.,347.,391.,327.,465.,583.,13*0.,300.
     #,339.,284.,403.,530.,14*0.,400.,334.,467.,556.5,15*0.,280.,391.,
     #484.,16*0.,550.,673.4,17*0.,246.,38*0./
      data a_f/  ! ---- Lennard-Jones A12-parameters (/1000,repulsion)
     #7.74,7.74,7.74,7.74,7.74,7.74,70.6,81.6,81.6,31.3,42.2,36.6,71.4,
     #62.,78.3,189.3,3*0.,7.74,7.74,7.74,7.74,7.74,70.6,61.7,61.7,31.3,
     #17.8,15.4,53.7,62.,58.7,189.3,4*0.,7.74,7.74,7.74,7.74,70.6,81.6,
     #81.6,31.3,42.2,36.6,71.4,62.,78.3,189.3,5*0.,7.74,7.74,7.74,70.6,
     #61.7,61.7,31.3,17.8,15.4,53.7,62.,58.7,189.3,6*0.,7.74,7.74,70.6,
     #81.6,81.6,31.3,42.2,36.6,71.4,62.,78.3,189.3,7*0.,7.74,70.6,81.6,
     #81.6,31.3,42.2,36.6,71.4,62.,78.3,189.3,8*0.,512.,601.,601.,256.,
     #349.,302.,538.,464.,598.,1196.,9*0.,704.,704.,298.,406.,351.,630.,
     #544.,699.,1203.5,10*0.,704.,298.,406.,351.,630.,544.,699.,1203.5,
     #11*0.,129.,176.,153.,269.,233.,303.,561.8,12*0.,240.,208.,366.,
     #317.,413.,772.5,13*0.,180.,317.,274.,358.,702.3,14*0.,565.,488.,
     #629.,1105.8,15*0.,421.,544.,976.8,16*0.,705.,1259.5,17*0.,503.3,
     #38*0./
c ---------------------------------------------- Hydrogen-bond donors
c                  1       2       3      4      5      6
      data do_f/.false.,.false.,.true.,.true.,.true.,.true.,
c                  7       8       9      10      11      12
     #          .false.,.false.,.false.,.false.,.false.,.false.,
c                 13      14      15      16      17      18
     #          .false.,.false.,.false.,.false.,.false.,.false./
c -------------------------------------------- Hydrogen-bond acceptors
c                  1       2       3       4       5       6
      data ac_f/.false.,.false.,.false.,.false.,.false.,.false.,
c                  7       8       9      10     11     12
     #          .false.,.false.,.false.,.true.,.true.,.true.,
c                 13      14      15     16     17      18
     #          .false.,.false.,.true.,.true.,.false.,.false./
c --------------------------------- HB-parameters (/1000,attraction)
      data chb_f/180.,180.,160. ,226.8,  ! given as  (ac_typ x do_typ)
     #           175.,175.,150. ,305.8,  ! to be used as:
     #           100.,100., 85. , 85. ,  !   (Do_typ x Ac_typ)
     #           165.,165.,150. ,305.8,  ! i.e.:
     #           720.,720.,643.1,845.6,  !   ( 3-6 x 10-12,15,16 )
     #             0.,  0.,  0. ,  0. /
c --------------------------------- HB-parameters (/1000,repulsion)
      data ahb_f/  6600.,  6600.,  6200., 12855.,
     #             6400.,  6400.,  7100., 29216.,
     #             2400.,  2400.,  2400.,  2000.,
     #             6200.,  6200.,  7100., 29216.,
     #           185000.,185000.,172301.,308235.,
     #                0.,      0.,    0.,     0./

c                   1  2   3   4  5   6   7  8   9 10 11  12  13 14  15
      data e0to_f /20.,2.,.8,2.6,.6,2.1,1.8,3.2,.5,.5,.8,1.2,1.3,1.,6.2
c                  16  17 18 19
     #           ,6.2, 8.,0.,20./
      data sgto_f /-1.,1.,-1., 1.,1., 1., 1.,1.,1.,1.,-1.,1., 1.,1.,-1.
     #           ,-1.,-1.,0.,-1./
      data rnto_f / 2.,3., 3., 3.,3., 6., 3.,3.,6.,6., 3.,6., 3.,3., 2.
     #            ,2., 2.,0.,2./

      data (rsnmcd(i),i=1,nrsty)/  ! Names for all amino acid residue types
     # 'ala ','arg ','arg+','asn ','asp ','asp-','cys ','cyss','gln ',
     # 'glu ','glu-','gly ','his ','hise','hisd','his+','hyp ','hypu',
     # 'ile ','leu ','lys ','lys+','met ','phe ','cpro','pro ','cpru',
     # 'prou','pron','pro+','ser ','thr ','trp ','tyr ','val ' /
  
      data (onltcd(i),i=1,nrsty)/  ! One-letter codes for amino acid types
     # 'A',   'R',   'R',   'N',   'D',   'D',   'C',   'C',   'Q',
     # 'E',   'E',   'G',   'H',   'H',   'H',   'H',   'P',   'P',
     # 'I',   'L',   'K',   'K',   'M',   'F',   'P',   'P',   'P',
     # 'P',   'P',   'P',   'S',   'T',   'W',   'Y',   'V' /

c     The vdW radii (in Angstr.) for the atomic groups and
c      coefficients for their solvation free energy (kcal/molxA**2)

c  Method: 

c  itysol=1 : OONS --> T.Ooi, et al, 
c                      Proc. Natl. Acad. Sci. USA 8 (1987) 3086-3090.
C  itysol=2 : JRF  --> J.Vila, et al, 
c                      PROTEINS: Struct Funct Genet 10(1991) 199-218.
C  itysol=3 : WE92 --> L.Wesson, D.Eisenberg, 
c                      Protein Science 1 (1992) 227-235.
C  itysol=4 : SCH1 --> D.Eisenberg, et al,
c                      Chem Scrip 29A (1989) 217-221.
C  itysol=5 : SCH2 --> A.H.Juffer, et al,
c                      Proteine Science 4 (1995) 2499-2509.
C  itysol=6 : SCH3 --> L.Wesson, D.Eisenberg, 
c                      Protein Science 1 (1992) 227-235.
C  itysol=7 : SCH4 --> C.A. Schiffer, et al,
c                      Mol. Simul. 10(1993) 121-149.
C  itysol=8 : EM86 --> D.Eisenberg, A.D. Mclachlan, 
c                      Nature 319 (1986) 199-203.
C  itysol=9 : BM   --> B. Freyberg, et al,
c                      J. Mol. Biol. 233 (1993) 275-292.

c ATOM
c TYPE OONS    JRF     WE92    SCH1    SCH2    SCH3    SCH4   EM86   BM

c 1   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.000 0.000
c 2   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.000 0.000
c 3   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.000 0.000
c 4   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.000 0.000
c 5   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.000 0.000
c 6   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.000 0.000
c 7   0.0080  0.2160  0.0120  0.0180  0.0130  0.0040  0.0325  0.016 1.000
c 8   0.4270 -0.7320  0.0120  0.0180  0.0130  0.0040  0.0325  0.016 1.000
c 9  -0.0080 -0.6780  0.0120  0.0180  0.0130  0.0040  0.0325  0.016 1.000
c10  -0.1720 -0.9100 -0.1160 -0.0090 -0.0070 -0.1130 -0.0175 -0.006 0.000
c11  -0.0380 -0.2620 -0.1750 -0.0090 -0.0070 -0.1660 -0.2800 -0.006 0.000
c12  -0.0380 -0.9100 -0.1750 -0.0370 -0.1120 -0.1660 -0.2800 -0.024 0.000
c13  -0.1320 -0.3120 -0.1860 -0.0380 -0.0870 -0.1690 -0.2175 -0.05  0.000
c14  -0.1320 -0.3120 -0.1160 -0.0090 -0.0070 -0.1130 -0.0175 -0.006 0.000
c15  -0.1320 -0.3120 -0.1160 -0.0090 -0.0070 -0.1130 -0.0175 -0.006 0.000
c16  -0.0210 -0.2810 -0.0180  0.0050 -0.0036 -0.0170 -0.0090  0.021 0.000
c17   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.000 0.000
c18   0.0080  0.2160  0.0120  0.0180  0.0130  0.0040  0.0325  0.016 1.000

       data coef_sl/54*0.0,
     # 0.008, 0.216, 0.012, 0.018, 0.013, 0.004, 0.0325, 0.016,1.000,
     # 0.427,-0.732, 0.012, 0.018, 0.013, 0.004, 0.0325, 0.016,1.000,
     # -.008,-0.678, 0.012, 0.018, 0.013, 0.004, 0.0325, 0.016,1.000,
     # -.172,-0.910,-0.116,-0.009,-0.007,-0.113,-0.0175,-0.006,0.000,
     # -.038,-0.262,-0.116,-0.009,-0.007,-0.113,-0.0175,-0.006,0.000,
     # -.038,-0.910,-0.175,-0.037,-0.112,-0.166,-0.2800,-0.024,0.000,
     # -.132,-0.312,-0.186,-0.038,-0.087,-0.169,-0.2175,-0.05, 0.000,
     # -.132,-0.312,-0.116,-0.009,-0.007,-0.113,-0.0175,-0.006,0.000,
     # -.132,-0.312,-0.116,-0.009,-0.007,-0.113,-0.0175,-0.006,0.000,
     # -.021,-0.281,-0.018, 0.005,-.0036,-0.017,-0.0090, 0.021,0.000,
     #  9*0.0,
     # 0.008, 0.216, 0.012, 0.018, 0.013, 0.004, 0.0325, 0.016,1.000/

       data rad_vdw/54*0.,
     #              2*2.0,6*1.9,2.0,
     #              2*1.55,6*1.9,1.5,
     #              2*1.75,6*1.9,1.85,
     #              9*1.4,
     #              9*1.4,
     #              9*1.4,
     #              2*1.55,6*1.7,1.5,
     #              2*1.55,6*1.7,1.5,
     #              2*1.55,6*1.7,1.5,
     #              2*2.0,6*1.8,1.85,
     #              9*0.,
     #              2*2.0,6*1.9,2.0/

      end


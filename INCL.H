      implicit integer*4 (i-n)
      implicit real*8 (a-h,o-z)
      character*255 version

      integer MAX_REPLICA, MAX_PROC, mxml, mxrs, mxat, mxbd, mxath
      integer mxvr, mxms, mxvw, mx14, mxvrh, mxtyat, mxhbdo
      integer mxhbac, mxtybl, mxtyba, mxtyto, mxrt, mxvrrt, mxrtsu
      integer nrsty
      double precision pi, pi2, pi4, crd, cdr, zero, one
      integer izero, ione

      parameter(MAX_PROC=2048)
      parameter (MAX_REPLICA = 2048)                                       ! max. number of replicas
      parameter (mxml=10)                                                ! max. number of molecules
      
      !*****************************************
      !INCREASING THE NUMBER OF RESIDUES (mxrs) AND THE
      !  NUMBER OF ATOMS (mxat) to 300 and 5600 
      ! for the protein INOL 
      ! this number could be increased or reduced
      !******************************************

      parameter(mxrs=300)                                               ! max. total number of residues
      parameter(mxat=5600)                                             ! max. total number of atoms
      parameter(mxbd=3)                                                 ! max. number of bonds to following atoms
      parameter(mxath=100)                                              ! max. number of atoms in help-arrays
      parameter(mxvr=mxrs*8)                                            ! max. number of local variables
      parameter(mxms=mxvr*3)                                            ! max. total number of moving sets
      parameter(mxvw=mxat*4)                                            ! max. number of vwd domains
      parameter(mx14=mxat*4)                                            ! max. number of '1-4' partners          
      parameter(mxvrh=mxath)                                            ! max. number of variables in help
      parameter(mxtyat=18)                                              ! max. # of energetic atom-types
      parameter(mxhbdo=4)                                                ! max. types of Hydrogens as donors in HB
     
      parameter (mxhbac=6)                                               ! max. types of atoms as acceptors in HB
      parameter(mxtybl=0)                                               ! max. # of types of bond length potentials
      parameter(mxtyba=0)                                               ! - " -   b. angle pot.
      parameter(mxtyto=19)                                              ! - " -   torsional pot.
      parameter(mxrt=500)                                               ! max. # of restraint types
      parameter(mxvrrt=mxrt*2)                                          ! max. # of restraint variables
      parameter(mxrtsu=mxvr*6)                                          ! max. # of pointers subsets -> restraints
      parameter(nrsty=35)                                                ! max. # of residue types

      parameter (pi=3.141592653589793d0)
      parameter(pi2=2.d0*pi)
      parameter(pi4=4.d0*pi)                               
      parameter(crd=57.2957795130823d0)                                
      parameter(cdr=1.745329251994330d-2)

      parameter (izero=0)
      parameter(ione=1)                                       
      parameter(zero=0.d0) 
      parameter(one=1.d0)
      parameter (version=                                               &
     &           '$Id: INCL.H 335 2007-08-09 09:32:21Z meinke $')

      integer seed
      integer enysolct
      common /counter/enysolct
      common /commonrandom/ seed
      common /bounds/ boxsize

      integer myrank, no, my_mpi_comm
      integer workPerProcessor(0:mxml, 0:MAX_PROC)
      common/paral/myrank,no,my_mpi_comm,workPerProcessor
      common /bet/beta
!==   For solvent============================
      parameter (mxpoint=2000)                                          !
      parameter(mxbox=80000)                                          !
      parameter(mxinbox=80000)                                        !
      parameter(mxtysol=9)                                               !Max number of solvent types
      parameter (rwater=1.4) !The Wdv radius of the water molecule

!==   Solvent types 
      logical isolscl
      common /isolty/ itysol, isolscl
      common /tes/ spoint(mxpoint,3)
      common /pnt/ npnt
      common /sol/ coef_sl(mxtysol,mxtyat),rad_vdw(mxtysol,mxtyat)
      common/ressurf/surfres(mxrs)
      dimension rsstrg(mxrs)
      common /abgncor/rsstrg

!----  The vdW radii of the atoms and solvation coefficients
       common /solvent/ rvdw(mxat),sigma(mxat)

!     vdvol - van der Waals volume of the molecule
!     asa   - Surface of molecule
      common /v_sol/ vdvol,asa,asap,asah

!     gradan - analytical gradients of solv. energy against atomic coordinates
      common /deriv/ gradan(mxat,3)


!-----For calculation of hydrogen bonds
      parameter(mxtyhb=100) !  Maximum number of hydrogen bond types
      common /h_bond/ ntyhb,ixtyhb(mxtyhb),nutyhb(mxtyhb),              &
     &                mmhb(mxml, mxml)
!-----For molecule parameter
      common/mol_par/ nvr

!-----Hydrophobic/polar radius of gyration
      common /rgy/ rgyh, rgyp
!=======================================  
       
      character nmrt(mxrt)*4,nmvrrt(mxvrrt)*3
      character rsnmcd(nrsty)*4,onltcd(nrsty)
      character nmath(mxath)*4,nmvrh(mxvrh)*3
      logical flex, tesgrd,epsd,sh2
! -----Switches for Bias Gaussian Steps and random order in variable  updates
! upchswitch=0 => No BGS
! upchswitch=1 => For backbone updates, use BGS with probability bgsprob
! upchswitch=2 => More sophisticated, temperature dependent choice between
! simple pivot update and BGS
      logical rndord
      integer upchswitch
      double precision bgsprob
      common /updchois/rndord,upchswitch,bgsprob
      common /con_i/ idloa,idloz,idupa,idupz      
!---- The contact matrix and the indices of C_alpha atoms
      common/cont_i_j/ ijcont(mxrs,mxrs),ind_alf(mxrs)
      common /jref/ iref(mxrs,mxrs),nci

!-----For calculation of Zimmermann code of conformation:
      character zimm*(mxrs)
      common/zimme/zimm 
      
!  lunlib - LUN of residue library
!  lunseq - LUN of sequence file
!  lunchg - LUN of charges for terminal residues (=0 if not necessary)
!  lunvar - LUN of file with commands to set variables
!  lunrst - LUN for variable restraints
!  reslib - name of residue library
!  seqfil - name of sequence file
!  chgfil - name of file with charges for terminal residues
!  varfil - name of file with commands to set variables
!  rstfil - name of file with variable restraints

      common /fil_i/ lunlib,lunseq,lunchg,lunvar,lunrst
      character reslib*80,seqfil*80,chgfil*80,varfil*80,rstfil*80
      common /fil_c/ reslib,seqfil,chgfil,varfil,rstfil

!  ntlml - total number of molecules
!  irsml1,irsml2 - indices of 1st/last residue of each molecule
!  ivrml1,nvrml - index of 1st/number of loc. variables in molecule
!  imsml1,nmsml - index of 1st/number of moving sets in molecule
!  iadml1,nadml - index of 1st/number of 'added varaibles' in molecule
!  nmml  - name of molecule

! -------------------------------------------------------------- new
!  ixrfpt - indices of 3 reference points (first 3 backbone atoms)
!           per molecule
!  ixatp - pointer of each SMMP atom to corresponding atom in a ref.
!          structure, (=0 if no corr. atom in ref. str.)
!  isrfvr() = true, if i-th variable is defined in ref. structure
!  rfpt - coordinates of a reference point (first bb atom) in molecule
!  xrfax,yrfax,zrfax - direction cosines of 3 ref. axes (axes J,K,L
!                      according to global angles alpha,beta,gamma)
!  gbpr - 6 global parameters per molecule (3 shifts, 3 angles)
! --------------------------------------------------------------------

      common /mol_i/ ntlml,irsml1(mxml),irsml2(mxml),ivrml1(mxml),      &
     &nvrml(mxml),imsml1(mxml),nmsml(mxml),iadml1(mxml),nadml(mxml)     &
     &,ixrfpt(3,mxml)                                                   &
     &,ixatp(mxat)                                                      &
     &,ireg
! -------------------

      logical isrfvr(mxvr)
      common /refstr/ isrfvr

      character nmml(mxml)*80
      common /mol_c/ nmml

! ----------------------------------------- new
      common /mol_r/ rfpt(3,mxml),gbpr(6,mxml),                         &
     &xrfax(3,mxml),yrfax(3,mxml),zrfax(3,mxml)
! ---------------------------------------------

!  iatrs1,iatrs2  - indices of 1st (backbone)/last atom of residue
!  ixatrs  - index of last backbone atom (may be bound to next res.)
!  ivrrs1,nvrrs  - index of 1st/number of local variable for residue
!  seq    - names of residues

      common /res_i/ iatrs1(mxrs),iatrs2(mxrs),ixatrs(mxrs),ivrrs1(mxrs)&
     &              ,nvrrs(mxrs)
      character seq(mxrs)*4
      common /res_c/ seq

!  ityat  - atom class
!  iowat  - index of preceeding atom (=0 if none)
!  iyowat - type of bond with prec. atom (=1:single bond,=-1:not single)
!  nbdat  - number of bonds with following atoms
!  ibdat  - indices of following covalently bonded atoms
!  iybdat - types of covalently bonds to following atoms (see 'iyowat')
!  ixmsat - index of the moving set for each atom
!  cgat   - point charge at atom
!  xbaat  - |
!  ybaat  - | axis to alter valence-angle for i-th atom
!  zbaat  - | (placed at 'iowat(i)')
!  baat   - valence angle 'iowat(iowat(i)), iowat(i), i'  [rad.]
!  snbaat - Sin(baat)
!  csbaat - Cos(baat)
!  xtoat  - |
!  ytoat  - | axis between 'iowat(i)' and i-th atom
!  ztoat  - |
!  toat   - dihedral 'iowat(iowat(iowat(i))), iowat(iowat(i)),
!                     iowat(i), i'  [rad.]
!  sntoat - Sin(toat)
!  cstoat - Cos(toat)
!  blat   - valence-length 'iowat(i)-i'
!  xat    - |
!  yat    - | cartesian coordinates of atom
!  zat    - |
!  nmat   - atom name

      common /atm_i/ ityat(mxat),iowat(mxat),iyowat(mxat),nbdat(mxat),  &
     &               ibdat(mxbd,mxat),iybdat(mxbd,mxat),ixmsat(mxat)    
      common /atm_r/ cgat(mxat),                                        &
     &xbaat(mxat),ybaat(mxat),zbaat(mxat),baat(mxat),snbaat(mxat),      &
     &csbaat(mxat),                                                     &
     &xtoat(mxat),ytoat(mxat),ztoat(mxat),toat(mxat),sntoat(mxat),      & 
     &cstoat(mxat),                                                     &
     &blat(mxat),                                                       &
     &xat(-1:mxat),yat(-1:mxat),zat(-1:mxat)                            
      character nmat(mxat)*4
      common /atm_c/ nmat

!  olvlvr - stores initial values for internal variables
!  nmvr  - name of variable
!  iatvr - index of primary moving atom for loc. variable
!  ityvr - type of variable ( 1 = valence length,
!                             2 = valence angle,
!                             3 = dihedral angle)
!  iclvr - indicates class for corresponding variable potential
!  iorvr - indicates appropriate order of variables
!  imsvr - index of 1st moving set per variable as stored in 
!  'latms1,latms2'
!  imsvr1,imsvr2 - indices of 1st & last moving set per variable 
!  'latms1,latms2'
!  latms1,latms2 - lists of indices for 1st and last atom of each moving set
!  iadvr1,iadvr2 - indices of 1st & last related variable per var. in 'ladvr'
!  ladvr - list of variable indices which are related to other variables
!  fxvr - .true. if variable will be fixed throughout calculations

      common /var_r/ olvlvr(mxvr),vlvr(mxvr),axvr(mxvr)
      common /var_i/ iatvr(mxvr),ityvr(mxvr),iclvr(mxvr),iorvr(mxvr),   &
     &               imsvr1(mxvr),imsvr2(mxvr),latms1(mxms),latms2(mxms)&
     &              ,iadvr1(mxvr),iadvr2(mxvr),ladvr(mxvr),idvr(mxvr)
      character nmvr(mxvr)*3
      common /var_c/ nmvr
      logical fxvr(mxvr)
      common /var_l/ fxvr



!  Energy parametrization for SCHERAGA set
!
!  ehm   - 1.5*(elem. charge)/(Planck\'s const./2Pi * electr. mass)
!  atpl  - atomic polarizabilities
!  efel  - empirical numbers of effective electrons per atom
!  emin  - minima of pairwise Lennard-Jones terms of like atoms
!  rmin  - corresponding pair-distance of energy minima
!  chb_s,ahb_s - parameter of attraction/repulsion terms of HB-potential
!  e0to_s,sgto_s,rnto_s - barrier[kcal/mole],sign,multiplicity of tors.potential
!  eps_s - distance-INdependent diel. constant
!  do_s,ac_s - flag to define atom-types serving as donors/acceptors in HB

      common /sche_r/ ehm,atpl(mxtyat),efel(mxtyat),emin(mxtyat),       &
     &    rmin(mxtyat),chb_s(mxhbdo,mxhbac),ahb_s(mxhbdo,mxhbac),       &
     &            e0to_s(mxtyto),sgto_s(mxtyto),rnto_s(mxtyto),         &
     &            eps_s
      logical do_s(mxtyat),ac_s(mxtyat)
      common /sche_l/ do_s,ac_s

!  Energy parametrization for FLEX set
!
!  c_f,a_f - parameter of attraction/repulsion terms of LJ-potential
!  chb_f,ahb_f - parameter of attraction/repulsion terms of HB-potential
!  e0to_s,sgto_s,rnto_s - barrier[kcal/mole],sign,multiplicity of tors.potential
!  plt_f,slp_f - parameters for sigmoidal distance-dependent epsilon
!  cohb_f - cut-off distance between H- & acceptor atom in Hydrogen-Bonds
!  do_f,ac_f - flags atom-types serving as donors/acceptors in HB

      common /flex_r/ c_f(mxtyat,mxtyat),a_f(mxtyat,mxtyat),            &
     &            chb_f(mxhbdo,mxhbac),ahb_f(mxhbdo,mxhbac),            &
     &            e0to_f(mxtyto),sgto_f(mxtyto),rnto_f(mxtyto),         &
     &            plt_f,slp_f,cohb_f
      logical do_f(mxtyat),ac_f(mxtyat)
      common /flex_l/ do_f,ac_f

!  Currently chosen energy parametrization
!
!  ihbty =  1, i - HB donor hydrogen & j - HB acceptor
!        = -1, i - HB acceptor & j - HB donor hydrogen
!        =  0, otherwise
!  cij,aij - parameters of attraction/repulsion terms of LJ-potential
!  a14 - special value for repulsion in "1-4"-interactions
!  chb,ahb - parameters of attraction/repulsion terms of HB-potential
!  e0to,sgto,rnto - barrier*0.5 !!,sign,multiplicity of torsional potential
!           esnto = e0/2 * sign * multipl.           
!  conv - parameter to convert electrostat. energy into [kcal/mole]
!  plt,slp,cohb - see FLEX
!  flex - .true. if FLEX-parametrization is to be used (obsolete)
!  ientyp : 0 => ECEPP, 1 => FLEX, 2 => Lund, 3 => ECEPP+Abagyan
!  tesgrd - .true. if GRADTEST to be included

      common /epar_i/ ihbty(mxtyat,mxtyat)
      common /epar_r/ cij(mxtyat,mxtyat),aij(mxtyat,mxtyat),            &
     & a14(mxtyat,mxtyat),chb(mxtyat,mxtyat),ahb(mxtyat,mxtyat),        &
     & e0to(mxtyto),sgto(mxtyto),rnto(mxtyto),esnto(mxtyto),            &
     & conv, plt,slp, cohb
      
      common /epar_l/ flex, tesgrd,sh2,epsd,ientyp

!  ENERGY & -MINIMIZATION

!  nvwml,ivwml1 - number/index of 1st vdW-domain per molecule
!  ivwat1,ivwat2 - index of 1st/last vdW-domain per atom
!  lvwat1,lvwat2 - lists of indices of 1st and last atom in each vdw domain
!  n14ml,i14ml1 - number/index of 1st '1-4' interaction partner per molecule
!  i14at1,i14at2 - index of 1st/last '1-4' interact. partner per atom
!  ixatvw - index of respective atom for each vdW-domain
!  ixatvw - index of respective atom for each '14'-interaction partner

! ---------------------------------------------------------------------------- new
!  eyel,eyvw,eyhb,eyvr,eysm,eysl,eyrg - Coulomb,vdW,HB,variable-potential & total 
!  eyslh, eyslp - hydrophilic (h)/phobic (p) conributions to the solvent energy
!   potential energy, solvation, regularization energy
!  gdeyvr - partial derivative of total energy vs. variable (Gradient)
!  gdeyrg - partial derivative of regularization energy vs. variable (Gradient)
!  gdeysl - partial derivative of solvatation energy vs. variable (Gradient)
!  wtey - weight to energy & energy gradient
!  wtrg - weight to reg. target & its gradient
      common /eny_i/ nvwml(mxml),ivwml1(mxml),n14ml(mxml),i14ml1(mxml), &
     &               ivwat1(mxat),ivwat2(mxat),i14at1(mxat),i14at2(mxat)&
     &              ,lvwat1(mxvw),lvwat2(mxvw),l14at(mx14)              &
     &              ,ixatvw(mxvw),ixat14(mx14)
      common /eny_r/ eyel,eyvw,eyhb,eyvr,eysm,eysl,eyrg,eyslh,eyslp,eyab&
     &               ,eysmi,eyhbi,eyeli,eyvwi                           &
     &               ,gdeyvr(mxvr),gdeyrg(mxvr),gdeysl(mxvr)            &
     &               ,gdeygb(6*mxml),wtey,wtrg

!  ntlrt  - total # of restraint types
!  ntlsu  - total # of variable subsets
!  nmrt   - name of r.t
!  pbrt   - probability for r.t
!  i1vrrt,i2vrrt - 1st,last restraint variable for a r.t
!  nmvrrt - name of restraint variable
!  vuvrrt - average value for r.v
!  dvvrrt - standard deviation for r.v
!  i1rtsu - point to 1st restraint type for a variable subset
!  nprtsu - # of restraint types for each subset
!  iprtsu - pointers from subsets --> restr. types
!  i1vrsu,i2vrsu - 1st,last pointer in 'ipvrsu' per var. subset
!  ipvrsu - pointers var. subset -> variables

      common /rstr_i/ ntlrt,ntlsu,i1vrrt(mxrt),i2vrrt(mxrt),            &
     &                iprtsu(mxrtsu),i1rtsu(mxvr),nrtsu(mxvr),          &
     &                i1vrsu(mxvr),i2vrsu(mxvr),ipvrsu(mxvr)

      common /rstr_c/ nmrt,nmvrrt
      common /rstr_r/ pbrt(mxrt),vuvrrt(mxvrrt),dvvrrt(mxvrrt)

!  rsnmcd - 4-letter codes for NRSTY residue types
!  onltcd - respective 1-letter codes for residue types

      common /rsnm_c/ rsnmcd,onltcd

!   arrays for help

      common /help_i/ ityath(mxath),iowath(mxath),iyowath(mxath),       &
     &                nbdath(mxath),ibdath(mxbd,mxath),                 &
     &                iybdath(mxbd,mxath),                              &
     &                ityvrh(mxvrh),iclvrh(mxvrh),iatvrh(mxvrh)
      common /help_r/ blath(mxath),baath(mxath),toath(mxath),           &
     &                cgath(mxath)
      common /help_c/ nmath,nmvrh


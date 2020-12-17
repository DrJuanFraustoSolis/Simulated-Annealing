c *******************************************************************
c SMMP version of Anders Irback's force field, to be called the Lund
c force field. This file contains the function enylun, which in turn 
c calls all the terms in the energy function. The terms Bias (ebias), 
c Hydrogen bonds (ehbmm and ehbms), Hydrophobicity (ehp) and the 
c Excluded volume (eexvol and eloexv) are also implemented in this 
c file.
c
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
      subroutine init_lundff
      include 'INCL.H'
      include 'incl_lund.h'
      character mynm*4
      logical prlvr

      print *,'initializing Lund forcefield'
c     Some parameters in the Lund force field.
c     The correspondence between internal energy scale and kcal/mol
      eunit=1.3315
c     Bias
      kbias=100.0*eunit
c      print *,'Bias'
c     Hydrogen bonds
      epshb1=3.1*eunit
      epshb2=2.0*eunit
      sighb=2.0
      cthb=4.5
      cthb2=cthb*cthb
      powa=0.5
      powb=0.5
      blhb=-30.0*(((sighb/cthb)**10-(sighb/cthb)**12))/cthb2
      alhb=-(5*((sighb/cthb)**12)-6*((sighb/cthb)**10))-blhb*cthb2
      sighb2=sighb*sighb
      cdon=1.0
      cacc=(1.0/1.23)**powb
      csacc=(1.0/1.25)**powb
c      print *,'Hydrogen bonds'
c     Hydrophobicity
c      print *,'Hydrophobicity with nhptyp = ',nhptyp

      hpstrg(1)=0.0*eunit
      hpstrg(2)=0.1*eunit
      hpstrg(3)=0.1*eunit
      hpstrg(4)=0.1*eunit
      hpstrg(5)=0.9*eunit
      hpstrg(6)=2.8*eunit
      hpstrg(7)=0.1*eunit
      hpstrg(8)=2.8*eunit
      hpstrg(9)=3.2*eunit

      do i=1,mxrs
         do j=1,6
            ihpat(i,j)=0
         enddo
         nhpat(i)=0
      enddo
      do i=irsml1(1),irsml2(ntlml)
         mynm=seq(i)
         call tolost(mynm)
         if ((mynm.eq.'pro').or.(mynm.eq.'cpro')
     #           .or.(mynm.eq.'cpru').or.(mynm.eq.'prou')
     #           .or.(mynm.eq.'pron').or.(mynm.eq.'pro+')) then
            prlvr=.true.        ! residue i is a proline variant
         else 
            prlvr=.false.
         endif 

         if (mynm.eq.'ala') then 
            nhpat(i)=1
            ihpat(i,1)=iCa(i)+2
         else if (mynm.eq.'val') then
            nhpat(i)=3
            ihpat(i,1)=iCa(i)+2
            ihpat(i,2)=iCa(i)+4
            ihpat(i,3)=iCa(i)+8
         else if (prlvr) then
            nhpat(i)=3
            ihpat(i,1)=iCa(i)+2
            ihpat(i,2)=iCa(i)+5
            ihpat(i,3)=iCa(i)+8
         else if (mynm.eq.'leu') then
            nhpat(i)=4
            ihpat(i,1)=iCa(i)+2
            ihpat(i,2)=iCa(i)+5
            ihpat(i,3)=iCa(i)+7
            ihpat(i,3)=iCa(i)+11
         else if (mynm.eq.'ile') then
            nhpat(i)=4
            ihpat(i,1)=iCa(i)+2
            ihpat(i,2)=iCa(i)+4
            ihpat(i,3)=iCa(i)+8
            ihpat(i,3)=iCa(i)+11
         else if (mynm.eq.'met') then
            nhpat(i)=4
            ihpat(i,1)=iCa(i)+2
            ihpat(i,2)=iCa(i)+5
            ihpat(i,3)=iCa(i)+8
            ihpat(i,3)=iCa(i)+9
         else if (mynm.eq.'phe') then
            nhpat(i)=6
            ihpat(i,1)=iCa(i)+5
            ihpat(i,2)=iCa(i)+6
            ihpat(i,3)=iCa(i)+8
            ihpat(i,3)=iCa(i)+10
            ihpat(i,3)=iCa(i)+12
            ihpat(i,3)=iCa(i)+14
         else if (mynm.eq.'tyr') then
            nhpat(i)=6
            ihpat(i,1)=iCa(i)+5
            ihpat(i,2)=iCa(i)+6
            ihpat(i,3)=iCa(i)+8
            ihpat(i,3)=iCa(i)+10
            ihpat(i,3)=iCa(i)+13
            ihpat(i,3)=iCa(i)+15
         else if (mynm.eq.'trp') then
            nhpat(i)=6
            ihpat(i,1)=iCa(i)+10
            ihpat(i,2)=iCa(i)+11
            ihpat(i,3)=iCa(i)+13
            ihpat(i,3)=iCa(i)+15
            ihpat(i,3)=iCa(i)+17
            ihpat(i,3)=iCa(i)+19
         endif
      enddo
c      print *,'Hydrophobicity'

c     Excluded volume and local pair excluded volume terms
      exvk=0.1*eunit
      exvcut=4.3
      exvcut2=exvcut*exvcut

      sigsa(1)=1.0  ! hydrogen
      sigsa(2)=1.0  ! hydrogen
      sigsa(3)=1.0  ! hydrogen
      sigsa(4)=1.0  ! hydrogen
      sigsa(5)=1.0  ! hydrogen
      sigsa(6)=1.0  ! hydrogen
      sigsa(7)=1.75 ! carbon
      sigsa(8)=1.75 ! carbon
      sigsa(9)=1.75 ! carbon
      sigsa(10)=1.42 ! oxygen
      sigsa(11)=1.42 ! oxygen
      sigsa(12)=1.42 ! oxygen
      sigsa(13)=1.55 ! nitrogen
      sigsa(14)=1.55 ! nitrogen
      sigsa(15)=1.55 ! nitrogen
      sigsa(16)=1.77 ! sulfur
      sigsa(17)=1.0  ! hydrogen
      sigsa(18)=1.75 ! carbon

      do i=1,mxtyat
         do j=1,mxtyat
            sig2lcp(i,j)=(sigsa(i)+sigsa(j))**2
            asalcp(i,j)=-7*((sig2lcp(i,j)/exvcut2)**6.0)
            bsalcp(i,j)=6*((sig2lcp(i,j)/exvcut2)**6.0)/exvcut2
         enddo
      enddo
c      print *,'Local pair excluded volume constants'

      exvlam=0.75
      exvcutg=exvcut*exvlam
      exvcutg2=exvcutg*exvcutg

      do i=1,mxtyat
         do j=1,mxtyat
            sig2exv(i,j)=(exvlam*(sigsa(i)+sigsa(j)))**2
            asaexv(i,j)=-7*((sig2exv(i,j)/exvcutg2)**6.0)
            bsaexv(i,j)=6*((sig2exv(i,j)/exvcutg2)**6.0)/exvcutg2
         enddo
      enddo
c      print *,'General excluded volume constants'

c     Initialization of the connections matrix matcon(i,j). The index
c     i runs from -mxconr to +mxconr, and j from 1 to mxat. 
c     matcon(i2-i1,i1) = 0, if the distance between atoms i1 and i2 is fixed
c                      = 2, if atoms i1 and i2 are separated by 3 covalent
c                           bonds and their distance can change
c                      = 1, for all other pairs
c     if abs(i2-i1) > mxconr, the atoms are assumed to be separated by 
c     many bonds, and with no restriction on their distances. On a protein 
c     molecule made of natural amino acids, atoms with indices separated
c     by more than 35 can not be connected by three covalent bonds.

      do i=1,mxat
         do j=-mxconr,mxconr
            matcon(j,i)=1
         enddo
         matcon(0,i)=0
      enddo
c     continued...
      do iml=1,ntlml
         do iat1=iatrs1(irsml1(iml)),iatrs2(irsml2(iml))
            do j=1,nbdat(iat1)
               iat2=ibdat(j,iat1)
               matcon(iat2-iat1,iat1)=0 ! directly bonded atoms
               matcon(iat1-iat2,iat2)=0 !
               do k=1,nbdat(iat2)
                  iat3=ibdat(k,iat2)
                  matcon(iat3-iat1,iat1)=0 !
                  matcon(iat1-iat3,iat3)=0 ! 2 covalent bonds
                  do l=1,nbdat(iat3)
                     iat4=ibdat(l,iat3)
                     if (iat4.ne.iat2) then
                        matcon(iat4-iat1,iat1)=2 !
                        matcon(iat1-iat4,iat4)=2 ! 3 covalent bonds
                     endif
                     do m=1,nbdat(iat2)
                        iat3p=ibdat(m,iat2)
                        if (iat3p.ne.iat3) then 
                           matcon(iat4-iat3p,iat3p)=2 ! 3 covalent bonds
                           matcon(iat3p-iat4,iat4)=2 !
                        endif
                     enddo
                  enddo
               enddo
               do k=1,nbdat(iat1)
                  iat3=ibdat(k,iat1)
                  matcon(iat3-iat2,iat2)=0 ! also 2 bonds iat2-iat1-iat3
                  matcon(iat2-iat3,iat3)=0 !
               enddo
            enddo
         enddo

c         print *,'going to initialize connections for first residue'
c         print *,'iN,iCa,iC =',iN(irsml1(iml)),
c     #        iCa(irsml1(iml)),iC(irsml1(iml))
         do iat1=iN(irsml1(iml))+1,iCa(irsml1(iml))-1
c            print *,'connections for iat1 = ',iat1
            matcon(iat1-iN(irsml1(iml)),iN(irsml1(iml)))=0
            matcon(iN(irsml1(iml))-iat1,iat1)=0
            matcon(iat1-iCa(irsml1(iml)),iCa(irsml1(iml)))=0
            matcon(iCa(irsml1(iml))-iat1,iat1)=0

            matcon(iat1-iCa(irsml1(iml))-1,iCa(irsml1(iml))+1)=2
            matcon(iCa(irsml1(iml))+1-iat1,iat1)=2
            matcon(iat1-iCa(irsml1(iml))-2,iCa(irsml1(iml))+2)=2
            matcon(iCa(irsml1(iml))+2-iat1,iat1)=2
            matcon(iat1-iC(irsml1(iml)),iC(irsml1(iml)))=2
            matcon(iC(irsml1(iml))-iat1,iat1)=2
         enddo

c     Below: for certain residues, some atoms separated by 3 or more bonds
c     do not change distance. So, the connection matrix term for such pairs
c     should be zero. 

         do irs=irsml1(iml),irsml2(iml)
            if (irs.eq.irsml1(iml)) then
               iatoff=1
            else 
               iatoff=0
            endif
            if (irs.eq.irsml2(iml)) then
               iatmrg=2
            else
               iatmrg=0
            endif
            mynm=seq(irs)
            call tolost(mynm)
            if ((mynm.eq.'pro').or.(mynm.eq.'cpro')
     #              .or.(mynm.eq.'cpru').or.(mynm.eq.'prou')
     #              .or.(mynm.eq.'pron').or.(mynm.eq.'pro+')) then
               prlvr=.true.     ! residue i is a proline variant
            else 
               prlvr=.false.
            endif 

            if ((mynm.eq.'arg').or.(mynm.eq.'arg+')) then
               do iat1=iatoff+iatrs1(irs)+13,iatrs2(irs)-2-iatmrg
                  do iat2=iat1+1,iatrs2(irs)-2-iatmrg
                     matcon(iat2-iat1,iat1)=0
                     matcon(iat1-iat2,iat2)=0
                  enddo
               enddo
            else if ((mynm.eq.'his').or.(mynm.eq.'hise')
     #              .or.(mynm.eq.'hisd').or.(mynm.eq.'his+')) then
               do iat1=iatoff+iatrs1(irs)+7,iatrs2(irs)-2-iatmrg
                  do iat2=iat1+1,iatrs2(irs)-2-iatmrg
                     matcon(iat2-iat1,iat1)=0
                     matcon(iat1-iat2,iat2)=0
                  enddo
               enddo
            else if (mynm.eq.'phe') then 
               do iat1=iatoff+iatrs1(irs)+7,iatrs2(irs)-2-iatmrg
                  do iat2=iat1+1,iatrs2(irs)-2-iatmrg
                     matcon(iat2-iat1,iat1)=0
                     matcon(iat1-iat2,iat2)=0
                  enddo
               enddo
            else if (mynm.eq.'tyr') then
               do iat1=iatoff+iatrs1(irs)+7,iatrs2(irs)-2-iatmrg
                  do iat2=iat1+1,iatrs2(irs)-2-iatmrg
                     if (iat2.ne.iatrs1(irs)+15) then
                        matcon(iat2-iat1,iat1)=0
                        matcon(iat1-iat2,iat2)=0
                     endif
                  enddo
               enddo
            else if (mynm.eq.'trp') then
               do iat1=iatoff+iatrs1(irs)+7,iatrs2(irs)-2-iatmrg
                  do iat2=iat1+1,iatrs2(irs)-2-iatmrg
                     matcon(iat2-iat1,iat1)=0
                     matcon(iat1-iat2,iat2)=0
                  enddo
               enddo
            else if (prlvr) then
c           Proline. Many more distances are fixed because of the fixed 
c           phi angle
               do iat1=iatoff+iatrs1(irs),iatrs2(irs)-2-iatmrg
                  do iat2=iat1+1,iatrs2(irs)-2-iatmrg
                     matcon(iat2-iat1,iat1)=0
                     matcon(iat1-iat2,iat2)=0
                  enddo
               enddo
c           distances to the C' atom of the previous residue are also fixed
               if (irs.ne.irsml1(iml)) then
                  iat1=iowat(iatrs1(irs))
                  do iat2=iatrs1(irs),iatrs2(irs)-1-iatmrg
                     matcon(iat2-iat1,iat1)=0
                     matcon(iat1-iat2,iat2)=0
                  enddo
               endif
            endif
         enddo
      enddo
c     finished initializing matrix conmat
c      print *,'Connections matrix'

c     Local pair excluded volume
      do i=1,mxml
         ilpst(i)=1
         ilpnd(i)=0
      enddo
      do i=1,50*mxrs
         lcp1(i)=0
         lcp2(i)=0
      enddo
      ilp=0
      do iml=1,ntlml
         do iat1=iatrs1(irsml1(iml)),iatrs2(irsml2(iml))
            do iat2=iat1+1,iatrs2(irsml2(iml))
               if ((iat2-iat1.le.mxconr).and.
     #                 matcon(iat2-iat1,iat1).eq.2) then
                  ilp=ilp+1
                  lcp1(ilp)=iat1
                  lcp2(ilp)=iat2
               endif
            enddo
         enddo

         ilpnd(iml)=ilp
         if (iml.lt.ntlml) then 
            ilpst(iml+1)=ilp+1
         endif
c         print *,'molecule ',iml,' lc pair range ',ilpst(iml),ilpnd(iml)
c         print *,'local pair list'
         do lci=ilpst(iml),ilpnd(iml)
            iat1=lcp1(lci)
            iat2=lcp2(lci)
c            print *,lci,iat1,iat2,matcon(iat2-iat1,iat1)
         enddo
      enddo

      print *,'finished initializing Lund force field'
      end

      integer function ihptype(iaa)
      include 'INCL.H'
      character mynm*4
      mynm=seq(iaa)
      ityp=-1
      call tolost(mynm)
      if (mynm.eq.'ala') then
         ityp=1
      else if ((mynm.eq.'val').or.(mynm.eq.'leu').or.(mynm.eq.'ile')
     #        .or.(mynm.eq.'met').or.(mynm.eq.'pro').or.(mynm.eq.'cpro')
     #        .or.(mynm.eq.'cpru').or.(mynm.eq.'prou')
     #        .or.(mynm.eq.'pron').or.(mynm.eq.'pro+')) then
         ityp=2
      else if ((mynm.eq.'phe').or.(mynm.eq.'tyr').or.(mynm.eq.'trp')) 
     #        then
         ityp=3
      endif
      ihptype=ityp
      return 
      end

      real*8 function ebiasrs(irsd)
      include 'INCL.H'
      include 'incl_lund.h'
      dimension q1(2),q2(2)
      data q1/-0.2,0.2/
      data q2/0.42,-0.42/ 

      et=0.0
      do i=0,1
         iat1=iN(irsd)+i
         do j=0,1
            iat2=iC(irsd)+j
            xij=xat(iat1)-xat(iat2)
            yij=yat(iat1)-yat(iat2)
            zij=zat(iat1)-zat(iat2)
            et=et+q1(i+1)*q2(j+1)/sqrt(xij*xij+yij*yij+zij*zij)
         enddo
      enddo
      ebiasrs=kbias*et
      return
      end
c     Evaluates backbone backbone hydrogen bond strength for residues 
c     i and j, taking the donor from residue i and acceptor from residue j
      real*8 function ehbmmrs(i,j)
      include 'INCL.H'
      include 'incl_lund.h'
      double precision r2,r4,r6,dx,dy,dz,ca,cb
      integer d1,d2,a1,a2
      d1=iN(i)
      d2=d1+1
      a1=iC(j)+1   ! dipoles are numbered from -ve to +ve charge
      a2=iC(j)     ! so, acceptor 1 is the Oxygen, and a2, the carbon
      dx=xat(a1)-xat(d2)
      dy=yat(a1)-yat(d2)
      dz=zat(a1)-zat(d2)
      r2=dx*dx+dy*dy+dz*dz
      if (r2.gt.cthb2) then 
c         print *,'hbmm = 0 ',cthb2,r2,a1,a2,d1,d2
c         print *,'a1,a2,d1,d2,r2 = ',a1,a2,d1,d2,r2,sighb2,cthb
         ehbmmrs=0
         return
      endif
      ca=(xat(d2)-xat(d1))*dx+(yat(d2)-yat(d1))*dy+(zat(d2)-zat(d1))*dz
      cb=(xat(a2)-xat(a1))*dx+(yat(a2)-yat(a1))*dy+(zat(a2)-zat(a1))*dz
      if (powa.gt.0.and.ca.le.0) then
c         print *,'hbmm, returning 0 because of angle a'
         ehbmmrs=0
         return
      endif
      if (powb.gt.0.and.cb.le.0) then
c         print *,'hbmm, returning 0 because of angle b'
         ehbmmrs=0
         return
      endif
      r6=sighb2/r2
      r4=r6*r6
      r6=r6*r4
      evlu=((ca*ca/r2)**(0.5*powa))*((cb*cb/r2)**(0.5*powb))
      evlu=evlu*(r6*(5*r6-6*r4)+alhb+blhb*r2)
c      print *,'found hbmm contribution ',evlu
      ehbmmrs=epshb1*evlu
      return
      end
      real*8 function enylun(nml)
c     nml = 1 .. ntlml. No provision exists to handle out of range values
c     for nml inside this function.
      include 'INCL.H'
      include 'incl_lund.h'
      character mynm*4
      logical prlvr
      eyhb=0.0   ! backbone-backbone and sidechain-backbone HB
      eyel=0.0   ! Bias, or local electrostatic term
      eysl=0.0   ! Hydrophobicity, replaces the solvent term of SMMP
      eyvr=0.0   ! Local pair excluded volume, in a sense a variable potential
      eyvw=0.0   ! atom-atom repulsion, excluded volume
c     atom-atom repulsion is calculated on a system wide basis, instead of
c     molecule by molecule for efficiency. Look into function exvlun.

      istres=irsml1(nml)
      indres=irsml2(nml)

c     First, all terms that can be calculated on a residue by residue basis
      do i=istres,indres
         mynm=seq(i)
         call tolost(mynm)
         if ((mynm.eq.'pro').or.(mynm.eq.'cpro')
     #           .or.(mynm.eq.'cpru').or.(mynm.eq.'prou')
     #           .or.(mynm.eq.'pron').or.(mynm.eq.'pro+')) then
            prlvr=.true.        ! residue i is a proline variant
         else 
            prlvr=.false.
         endif 

c     Bias, or local electrostatic term. Excluded from the list are
c     residues at the ends of the chain, glycine and all proline variants
         if ((i.ne.istres).and.(i.ne.indres).and.
     #           .not.prlvr.and.mynm.ne.'gly') then
            eyel=eyel+ebiasrs(i)
         endif
c     Backbone--backbone hydrogen bonds
         shbm1=1.0
         shbm2=1.0
         if ((i.eq.istres).or.(i.eq.indres)) shbm1=0.5
c     Residue i contributes the donor, and j, the acceptor, so both i and
c     j run over the whole set of amino acids.
c     No terms for residue i, if it is a proline variant.
         if (.not.prlvr) then  
            do j=istres,indres
               if ((j.eq.istres).or.(j.eq.indres)) shbm2=0.5
               etmp=ehbmmrs(i,j)
               eyhb=eyhb+shbm1*shbm2*etmp
            enddo
         endif
c     Hydrophobicity, only if residue i is hydrophobic to start with
         ihpi=ihptype(i)
         if (ihpi.ge.0) then 
c        Unlike hydrogen bonds, the hydrophobicity potential is symmetric
c        in i and j. So, the loop for j runs from i+1 to the end.

            do j=i+1,indres
               ihpj=ihptype(j)
               if (ihpj.ge.0) then
                  etmp=ehp(i,j,ihpi,ihpj)
                  if (j.eq.(i+1)) etmp=0
                  if (j.eq.(i+2)) etmp=0.5*etmp
                  eysl=eysl+etmp
               endif 
            enddo
         endif
      enddo

c     Terms that are not calculated residue by residue ...

c     Local pair or third-neighbour excluded volume
c     Numerically this is normally the term with largest positive
c     contribution to the energy in an equilibrated stystem. 

      i1=ilpst(nml)
      i2=ilpnd(nml)
      etmp=0.0
      do i=i1,i2
         etmp1=0.0
         iat1=lcp1(i)
         iat2=lcp2(i)
         iatt1=ityat(iat1)
         iatt2=ityat(iat2)
         xij=xat(iat1)-xat(iat2)
         yij=yat(iat1)-yat(iat2)
         zij=zat(iat1)-zat(iat2)
         r2=(xij*xij + yij*yij + zij*zij)
         if (r2.le.exvcut2) then
            r6=sig2lcp(iatt1,iatt2)/r2 
            r6=r6*r6*r6
            etmp1=(r6*r6+asalcp(iatt1,iatt2)+bsalcp(iatt1,iatt2)*r2)
            if (etmp1.ge.500) then
               print *,'local pair contribution ',iat1,iat2,etmp1
            endif
            etmp=etmp+etmp1
         endif
c         print *,'pair : ',iat1,iat2,' contribution ',etmp1 
c         print *,exvcut2,r2
      enddo
      eyvr=exvk*etmp
      eysm=eyel+eyhb+eyvr+eysl
      enylun=eysm
      end
      real*8 function eninlun()
      include 'INCL.H'
        eysmi=0.0
        eyeli=0.0
        eyvwi=0.0
        eyhbi=0.0
        eninlun=0.0
        return
      end
      real*8 function ehp(i1,i2,ihp1,ihp2)
      include 'INCL.H'
      include 'incl_lund.h'
      dimension r2min(12)

      b2=20.25
      a2=12.25
c      ihp1=ihptype(i1)
c      ihp2=ihptype(i2)
      if ((ihp1.le.0).or.(ihp2.le.0)) then
         ehp=0.0
         return
      endif
      ni=nhpat(i1)
      nj=nhpat(i2)
      do i=1,ni+nj
         r2min(i)=100000000  ! quite far
      enddo
      do i=1,ni
         k=ihpat(i1,i)
         do j=1,nj
            l=ihpat(i2,j)
            xij=xat(k)-xat(l)
            yij=yat(k)-yat(l)
            zij=zat(k)-zat(l)
            dtmp=(xij*xij + yij*yij + zij*zij)
            if (dtmp.le.r2min(i)) r2min(i)=dtmp
            if (dtmp.le.r2min(ni+j)) r2min(ni+j)=dtmp
         enddo
      enddo
      sum=0
      do i=1,ni+nj
         if (r2min(i).le.b2) then 
            if (r2min(i).lt.a2) then 
               sum=sum+1
            else 
               sum=sum+(b2-r2min(i))/(b2-a2)
            endif
         endif
      enddo
      ehp=-hpstrg((ihp1-1)*nhptyp+ihp2)*sum/(ni+nj)
      return
      end

      real*8 function exvlun(nml)
      include 'INCL.H'
      include 'incl_lund.h'
c     For multi-chain systems it makes little sense to split the calculation
c     of this term into an 'interaction part' and a contribution from 
c     individual molecules. So, normally this should always be called with
c     argument nml=0. Only for diagnostic reasons, you might want to find
c     the contribution from one molecule in a multi-chain system assuming 
c     there was no other molecule.
      dimension isort(mxat),ngbr(mxat),locccl(mxat),incell(mxcell)
      dimension icell(mxat)
      if (nml.eq.0) then 
         istat=iatrs1(irsml1(1))
         indat=iatrs2(irsml2(ntlml))
      else 
         istat=iatrs1(irsml1(nml))
         indat=iatrs2(irsml2(nml))
      endif

      eyvw=0.0
c     The beginning part of this implementation is very similar to the 
c     assignment of cells to the atoms during calculation of solvent 
c     accessible surface area. So, much of that part is similar. But 
c     unlike the accessible surface calculations, this term is symmetric
c     in any two participating atoms. So, the part after the assignment
c     of cells differs even in the structure of the code.

      do i=1,mxcell
         incell(i)=0
      enddo
c      print *,'evaluating general excluded volume :',istat,',',indat
c     Find minimal containing box
      xmin=xat(istat)
      ymin=yat(istat)
      zmin=zat(istat)
      xmax=xmin
      ymax=ymin
      zmax=zmin

      do i=istat,indat
         if (xat(i).le.xmin) then
            xmin=xat(i)
         else if (xat(i).ge.xmax) then
            xmax=xat(i)
         endif
         if (yat(i).le.ymin) then
            ymin=yat(i)
         else if (yat(i).ge.ymax) then
            ymax=yat(i)
         endif
         if (zat(i).le.zmin) then
            zmin=zat(i)
         else if (zat(i).ge.zmax) then
            zmax=zat(i)
         endif
      enddo

      sizex=xmax-xmin
      sizey=ymax-ymin
      sizez=zmax-zmin
c     Number of cells along each directions that fit into the box.
      ndx=int(sizex/exvcutg)+1
      ndy=int(sizey/exvcutg)+1
      ndz=int(sizez/exvcutg)+1

      nxy=ndx*ndy
      ncell=nxy*ndz
c      print *,'Number of cells along x,y,z = ',ndx,',',ndy,',',ndz
      if (ncell.ge.mxcell) then
         print *,'exvlun> required number of cells',ncell,
     #        ' exceeded the limit ',mxcell
         print *,'recompile with a higher mxcell.'
         stop
      endif
c     Expand box to contain an integral number of cells along each direction
      shiftx=(dble(ndx)*exvcutg-sizex)/2.0
      shifty=(dble(ndy)*exvcutg-sizey)/2.0
      shiftz=(dble(ndz)*exvcutg-sizez)/2.0
      xmin=xmin-shiftx
      ymin=ymin-shifty
      zmin=zmin-shiftz
      xmax=xmax+shiftx
      ymax=ymax+shifty
      zmax=zmax+shiftz

c     Set occupied cells to zero. Note that the maximum number of occupied
c     cells is the number of atoms in the system.
      nocccl=0
      do i=1,mxat
         locccl(i)=0
      enddo

c     Put atoms in cells
      do j=istat,indat
         mx=min(int(max((xat(j)-xmin)/exvcutg,0.0d0)),ndx-1)
         my=min(int(max((yat(j)-ymin)/exvcutg,0.0d0)),ndy-1)
         mz=min(int(max((zat(j)-zmin)/exvcutg,0.0d0)),ndz-1)
         icellj=mx+my*ndx+mz*nxy+1
         icell(j)=icellj
         if (icellj.gt.mxcell) then 
            print *,'exvlun> bad cell index ',icellj,' for atom ',j
            stop
         else 
            if (incell(icellj).eq.0) then 
c           previously unoccupied cell
               nocccl=nocccl+1
               locccl(nocccl)=icellj
            endif
            incell(icellj)=incell(icellj)+1
         endif
      enddo
c      print *,'finished assigning cells. nocccl = ',nocccl
c     Cummulative occupancy of i'th cell
      do i=1,ncell
         incell(i+1)=incell(i+1)+incell(i)
      enddo
c      print *,'finished making cumulative cell sums'
c     Sorting atoms by their cell index
      do i=istat,indat
         j=icell(i)
         jj=incell(j)
         isort(jj)=i
         incell(j)=jj-1
      enddo
c      print *,'sorted atoms by cell index'
      etmp=0.0
      do icl=1,nocccl
c     loop through occupied cells
         lcell=locccl(icl)
         ix=mod(lcell-1,ndx)
         iy=(mod(lcell-1,nxy)-ix)/ndx
         iz=(lcell-1-ix-ndx*iy)/nxy
c         print *,'icl=',icl,'absolute index of cell = ',lcell
c         print *,'iz,iy,ix = ',iz,iy,ix
c     find all atoms in current cell and all its forward-going neighbours
         nex=min(ix+1,ndx-1)
         ney=min(iy+1,ndy-1)
         nez=min(iz+1,ndz-1)
         nsame=0
         nngbr=0
         do jx=ix,nex
            do jy=iy,ney
               do jz=iz,nez
                  jcl=jx+ndx*jy+nxy*jz+1
                  do ii=incell(jcl)+1,incell(jcl+1)
c                    count the total number of neighbours
                     nngbr=nngbr+1
                     if (jx.eq.ix.and.jy.eq.iy.and.jz.eq.iz) then
c                    count how many neighbours are from the same cell
                        nsame=nsame+1
                     endif
                     ngbr(nngbr)=isort(ii)
                  enddo
               enddo
            enddo
         enddo
c     A few more cells need to be searched, so that we cover 13 of the 26
c     neighbouring cells. 
c        1
         jx=ix+1
         jy=iy
         jz=iz-1
         jcl=jx+ndx*jy+nxy*jz+1
         do ii=incell(jcl)+1,incell(jcl+1)
            nngbr=nngbr+1
            ngbr(nngbr)=isort(ii)
         enddo
c        2
         jx=ix
         jy=iy-1
         jz=iz+1
         jcl=jx+ndx*jy+nxy*jz+1
         do ii=incell(jcl)+1,incell(jcl+1)
            nngbr=nngbr+1
            ngbr(nngbr)=isort(ii)
         enddo
c        3
         jx=ix-1
         jy=iy+1
         jz=iz
         jcl=jx+ndx*jy+nxy*jz+1
         do ii=incell(jcl)+1,incell(jcl+1)
            nngbr=nngbr+1
            ngbr(nngbr)=isort(ii)
         enddo
c        4
         jx=ix+1
         jy=iy+1
         jz=iz-1
         jcl=jx+ndx*jy+nxy*jz+1
         do ii=incell(jcl)+1,incell(jcl+1)
            nngbr=nngbr+1
            ngbr(nngbr)=isort(ii)
         enddo
c        5
         jx=ix+1
         jy=iy-1
         jz=iz+1
         jcl=jx+ndx*jy+nxy*jz+1
         do ii=incell(jcl)+1,incell(jcl+1)
            nngbr=nngbr+1
            ngbr(nngbr)=isort(ii)
         enddo
c        6
         jx=ix+1
         jy=iy-1
         jz=iz-1
         jcl=jx+ndx*jy+nxy*jz+1
         do ii=incell(jcl)+1,incell(jcl+1)
            nngbr=nngbr+1
            ngbr(nngbr)=isort(ii)
         enddo

c         print *,'atoms in same cell ',nsame
c         print *,'atoms in neighbouring cells ',nngbr
         do i1=1,nsame
c        Over all atoms from the original cell
            iat1=ngbr(i1)
            do i2=i1,nngbr
c           Over all atoms in the original+neighbouring cells
               iat2=ngbr(i2)
               xij=xat(iat1)-xat(iat2)
               yij=yat(iat1)-yat(iat2)
               zij=zat(iat1)-zat(iat2)
               r2=(xij*xij+yij*yij+zij*zij)

               if (r2.le.exvcutg2) then 
                  if (abs(iat2-iat1).gt.mxconr.or.
     #                 matcon(iat2-iat1,iat1).eq.1) then
                     iatt1=ityat(iat1)
                     iatt2=ityat(iat2)
                     r6=sig2exv(iatt1,iatt2)/r2
                     r6=r6*r6*r6
                     etmp1=r6*r6+asaexv(iatt1,iatt2)
     #                    +bsaexv(iatt1,iatt2)*r2
                     etmp=etmp+etmp1
c                     if (etmp1.ge.2000) then
c                        print *,'contribution ',iat1,iat2,etmp1
c                        call outpdb(1,'EXAMPLES/clash.pdb')
c                        stop
c                     endif
                  endif
               endif
            enddo
         enddo
      enddo
c      irs=1
c      do iat=iatrs1(irs),iatrs2(irs)
c         do j=-mxconr,mxconr
c            print *,iat,j,':',matcon(j,iat)
c         enddo
c      enddo
c      irs=irsml2(1)
c      do iat=iatrs1(irs),iatrs2(irs)
c         do j=-mxconr,mxconr
c            print *,iat,j,':',matcon(j,iat)
c         enddo
c      enddo

      eyvw=exvk*etmp
      exvlun=eyvw
      return
      end

      real*8 function exvbrfc()
c     Brute force excluded volume evaluation
      include 'INCL.H'
      include 'incl_lund.h'

      etmp=0.0
      etmp1=0.0
      print *,'max connection radius is ',mxconr
      do iat1=iatrs1(irsml1(1)),iatrs2(irsml2(ntlml))
         do iat2=iat1+1,iatrs2(irsml2(ntlml))
            xij=xat(iat1)-xat(iat2)
            yij=yat(iat1)-yat(iat2)
            zij=zat(iat1)-zat(iat2)
            r2=(xij*xij+yij*yij+zij*zij)

            if (r2.le.exvcutg2) then 
               if (abs(iat2-iat1).gt.mxconr.or.
     #                 matcon(iat2-iat1,iat1).eq.1) then
                  iatt1=ityat(iat1)
                  iatt2=ityat(iat2)
                  r6=sig2exv(iatt1,iatt2)/r2
                  r6=r6*r6*r6
                  etmp1=r6*r6+asaexv(iatt1,iatt2)
     #                 +bsaexv(iatt1,iatt2)*r2
                  etmp=etmp+etmp1
                  if (iat1.eq.43.and.iat2.eq.785) then
                     print *,'contribution ',iat1,iat2,etmp1
                     print *,'r2 = ',r2
                     print *,'atom types : ',iatt1,iatt2
                     print *,'sig2exv = ',sig2exv(iatt1,iatt2)
                     print *,'asa = ',asaexv(iatt1,iatt2)
                     print *,'bsa = ',bsaexv(iatt1,iatt2)

c     call outpdb(1,'EXAMPLES/clash.pdb')
c     stop
                  endif
               else 
c                  print *,'atoms ', iat1,' and ',iat2,' were close',
c     #                 'but matcon is ',matcon(iat2-iat1,iat1)
               endif
            endif
         enddo
      enddo
      exvbrfc=etmp*exvk
      return
      end

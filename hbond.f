c**************************************************************
c
c This file contains the subroutines: hbond,chhb,ishybd,
c                                     ishybdo,nursat,interhbond
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************


      subroutine hbond(nml,mhb,ipr)
c .................................................................
c PURPOSE: find hydrogen bonds in molecule 'nml'
c
c          prints HBonds, if ipr > 0
c
c OUTPUT: mhb - number of hyd.bds. of type i->i+4
c  
c   to INCL.H:
c
c         ntyhb  - number of different types of hyd. bds. found
c         nutyhb - number of hyd.bds. found for each type
c         ixtyhb - index for each type of hyd. bd. composed as
c                  (atom idx. of H) * 1000 + atm.idx. of acceptor
c
c CALLS: chhb,ishybd  (ishybdo),nursat
C
c................................................................

      include 'INCL.H'

cf2py intent(out) mhb  
      parameter (atbase=mxat)      
      logical ishb


      do i=1,mxtyhb
        nutyhb(i) = 0
        ixtyhb(i) = 0
      enddo
      ntyhb=0
      if (nml.eq.0) then
        ntlvr = nvr
        ifivr = ivrml1(1)
      else
        ntlvr=nvrml(nml)
        if (ntlvr.eq.0) then
          write (*,'(a,i4)')
     #           ' hbond> No variables defined in molecule #',nml
          return
        endif

        ifivr=ivrml1(nml)
c Index of last moving set
        i1s=imsml1(nml)+nmsml(nml)
      endif 
c Loop over all variables
      do io=ifivr+ntlvr-1,ifivr,-1  
c Get index of variable
        iv=iorvr(io)       
c Index of next to last moving set
        i2s=i1s-1      
c Index of moving set belonging to iv
        i1s=imsvr1(iv) 
c Loop over all moving sets between the one belonging to iv and the 
c next to last one
        do ims=i1s,i2s  
c First atom in moving set
          i1=latms1(ims)
c Last atom in moving set
          i2=latms2(ims)
c Loop over all atoms in moving set.
          do i=i1,i2  
c Loop over van der Waals domains of atom i
            do ivw=ivwat1(i),ivwat2(i)
c Loop over atoms in van der Waals domain.  
              do j=lvwat1(ivw),lvwat2(ivw)  

                call ishybd(i,j,ishb,ih,ia)   ! Thornton criteria

                if (ishb) then

                  ixhb=ih*atbase+ia

                  do k=1,ntyhb
                    if (ixhb.eq.ixtyhb(k)) then
                      nutyhb(k)=nutyhb(k)+1
                      goto 1
                    endif
                  enddo

                  if (ntyhb.lt.mxtyhb) then
                    ntyhb=ntyhb+1
                    nutyhb(ntyhb)=1
                    ixtyhb(ntyhb)=ixhb
                  else
                    write(*,*) ' hbond> increase parameter MXTYHB'
                    stop
                  endif

                endif  ! have h.b.

    1         enddo  ! ... atoms j
            enddo  ! ... vdW-domains of i

            do i14=i14at1(i),i14at2(i)   !  over 1-4 partners of 'i'
              j=l14at(i14)

              call ishybd(i,j,ishb,ih,ia)   ! Thornton criteria
c              call ishybdo(i,j,ishb,ih,ia)  

              if (ishb) then

                ixhb=ih*atbase+ia

                do k=1,ntyhb
                  if (ixhb.eq.ixtyhb(k)) then
                    nutyhb(k)=nutyhb(k)+1
                    goto 2
                  endif
                enddo

                if (ntyhb.lt.mxtyhb) then
                  ntyhb=ntyhb+1
                  nutyhb(ntyhb)=1
                  ixtyhb(ntyhb)=ixhb
                else
                  write(*,*) ' hbond> increase parameter MXTYHB'
                  stop
                endif

              endif ! have h.b.

    2       enddo  ! ... 1-4-partners of i

          enddo  ! ... atoms i
        enddo  ! ... m.s.

      enddo  ! ... variables

      mhb=0

c     do inhb=1,ntyhb
c      mhb = mhb+nutyhb(inhb)
c     enddo

      if (ipr.gt.0) write(*,'(1x,a,/)') ' hbond>  Hydrogen Bonds:'

      if (ntyhb.gt.0) then

        ii = 0
        do i=1,ntyhb
          jj = nutyhb(i)
          do j = 1,jj
            ii = ii + 1
            ix = ixtyhb(ii)

            id =ix / atbase         ! donor atom
            nd = nursat(id)       !  & residue

            ia = ix - id * atbase   ! acceptor atom
            na = nursat(ia)

            n = nd - na

            if (n.gt.4) mhb = mhb+1  ! only count these

            if (ipr.gt.0) then  

              if (n.gt.0) then
                write(*,'(1x,i3,a2,a4,a3,i3,1x,a4,a7,a4,a3,i3,1x,a4,a9,
     #                    i2)')
     #           ii,') ',nmat(ia),' ( ',na,seq(na),' ) <-- ',nmat(id),
     #           ' ( ', nd,seq(nd),' ) = i +',n
              else
                write(*,'(1x,i3,a2,a4,a3,i3,1x,a4,a7,a4,a3,i3,1x,a4,a9,
     #                    i2)')
     #           ii,') ',nmat(ia),' ( ',na,seq(na),' ) <-- ',nmat(id),
     #           ' ( ', nd,seq(nd),' ) = i -',abs(n)
              endif

              call chhb(ia,id)
            endif
 
          enddo
         enddo
      endif

      return
      end
c .....................................................................
c Calculates hydrogen bonds between different chains.
c 
c @return number of intermolecular hydrogen bonds. Returns 0 if only 
c         one molecule is present. The value is returned in the 
c         variable mhb.
c
c @author Jan H. Meinke <j.meinke@fz-juelich.de>
c                                            
c .....................................................................
      subroutine interhbond(mhb)

      include 'INCL.H'
      
cf2py intent(out) mhb     
      
      logical ishb
      
      integer*4 mhb
      mhb = 0
      do iml = 1, ntlml
        do jml = iml + 1, ntlml
          mmhb(iml, jml) = 0
          do ires= irsml1(iml), irsml2(iml)
            do jres= irsml1(jml), irsml2(jml)
              do iat = iatrs1(ires), iatrs2(ires)
                do jat = iatrs1(jres), iatrs2(jres)
                  call ishybd(iat,jat,ishb,ih,ia)
                  if (ishb) then
                    mhb = mhb + 1
                    mmhb(iml, jml) = mmhb(iml, jml) + 1
                  endif
                enddo ! jat    
              enddo ! iat
            enddo ! jres
          enddo ! ires
          mmhb(jml, iml) = mmhb(iml, jml)
        enddo ! jml
      enddo ! iml
      
      end ! subroutine interhbond
c ************************
      subroutine chhb(i,j)

      include 'INCL.H'

      ihb = ihbty(ityat(i),ityat(j))

      if (ihb.gt.0) then
        ih=i
        ia=j
      else
        ih=j
        ia=i
      endif

      dah=sqrt((xat(ih)-xat(ia))**2+(yat(ih)-yat(ia))**2+
     #          (zat(ih)-zat(ia))**2)

      id=iowat(ih)

      dad=sqrt((xat(id)-xat(ia))**2+(yat(id)-yat(ia))**2+
     #          (zat(id)-zat(ia))**2)
      adha=valang(id,ih,ia)*crd

      ib=iowat(ia)

      ahab=valang(ih,ia,ib)*crd
      adab=valang(id,ia,ib)*crd

      write(*,*) '  '
      write(*,*) ' Dah: ',dah,' Dad: ',dad
      write(*,*) ' Adha: ',adha,' Ahab: ',adha,' Adab: ',adab
      write(*,*) '  '

      return
      end
c *************************************
      subroutine ishybd(i,j,ishb,ih,ia)
      

c ..........................................................
c  PURPOSE: checks for hydrogen bond between atoms 'i' & 'j'
c           according to geometric criteria
c 
c  OUTPUT:  logical 'ishb' - true, if have Hydrogen bond
c           ih - index of Hydrogen atom
c           ia - index of Acceptor atom
c
c  [I.K.McDonald,J.M.Thornton,Satisfying hydrogen bond
c   potential in proteins.J.Mol.Biol.238(5),777-793 (1994)]
c
c  D: hydrogen(=H) donor, A: acceptor, B: atom bound to A
c
c  Dis_HA <= 2.5 & Dis_DA <= 3.9 & Angle(D-H-A) > 90 &
c  Angle(H-A-B) > 90 & Angle(D-A-B) > 90
c ..........................................................

      include 'INCL.H'

      parameter (cdad=3.9d0,
     #           cdah=2.5d0,
     #           cang=110.d0)
c     #           cang=90.d0)

      logical ishb

      cahb = cang * cdr
      ishb = .false.


      if (i.le.0.or.j.le.0) return

      ihb = ihbty(ityat(i),ityat(j))

      if (ihb.eq.0) then
        return
      elseif (ihb.gt.0) then
        ih=i
        ia=j
      else
        ih=j
        ia=i
      endif

      if (sqrt((xat(ih)-xat(ia))**2+(yat(ih)-yat(ia))**2+
     #         (zat(ih)-zat(ia))**2).gt.cdah) return

      id=iowat(ih)

      if (id.le.0.or.sqrt((xat(id)-xat(ia))**2+(yat(id)-yat(ia))**2+
     #                    (zat(id)-zat(ia))**2).gt.cdad
     #           .or.valang(id,ih,ia).lt.cahb) return

      ib=iowat(ia)

      if (ib.gt.0.and.valang(ih,ia,ib).ge.cahb
     #           .and.valang(id,ia,ib).ge.cahb) ishb=.true.

      return
      end
c **************************************

      subroutine ishybdo(i,j,ishb,ih,ia)

c ..........................................................
c  PURPOSE: checks for hydrogen bond between atoms 'i' & 'j'
c           according to geometric criteria
c 
c  OUTPUT:  logical 'ishb' - true, if have Hydrogen bond
c           ih - index of Hydrogen atom
c           ia - index of Acceptor atom
c
c  D: hydrogen(=H) donor, A: acceptor
c
c    Dis_AH <= 2.5 & Angle(D-H-A) >= 160
c ...........................................................

      include 'INCL.H'

      parameter (cdah=2.5d0,
     #           cang=140.d0)

      logical ishb

      cahb = cang * cdr


      ishb = .false.

      if (i.le.0.or.j.le.0) return

      ihb = ihbty(ityat(i),ityat(j))

      if (ihb.eq.0) then
        return
      elseif (ihb.gt.0) then
        ih=i
        ia=j
      else
        ih=j
        ia=i
      endif

      if (sqrt((xat(ih)-xat(ia))**2+(yat(ih)-yat(ia))**2+
     #         (zat(ih)-zat(ia))**2).gt.cdah) return

      id=iowat(ih)

      if (id.gt.0.and.valang(id,ih,ia).ge.cahb)  ishb=.true.

      return
      end


! ****************************************************************
! Trial version implementing the semi-local conformational update
! BGS (Biased Gaussian Steps). This file presently contains the
! functions initlund, bgsposs and bgs. 
!
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
! ****************************************************************

! Subroutine initlund: Initializes data structures used frequently 
! in connection with Biased Gaussian Steps and the energy functions
! from Anders Irback's protein folding model. Calls: none.
!
      subroutine init_lund
      include 'INCL.H'
      include 'incl_lund.h'
      logical bgsposs
      do i=1,mxrs
         iN(i)=-1
         iCa(i)=-1
         iC(i)=-1
         iphi(i)=-34
         ipsi(i)=-35
      enddo
!      print *,'total number of variables = ',nvr

      do i=1,ntlml
         npprs=1
         do j=ivrml1(i),ivrml1(i)+nvrml(i)-1
            mlvr(j)=i
            if (nmvr(j).eq.'phi') then 
               iphi(npprs)=j
! Now if the residue is a proline, there is no phi angle in the variable
! list in SMMP, and so iphi(j) will remain at the initial value. 
! So, make sure you never use iphi(i) for proline. 
            endif
            if (nmvr(j).eq.'psi'.or.nmvr(j).eq.'pst') then 
               ipsi(npprs)=j
               npprs=npprs+1
            endif
         enddo
         do j=irsml1(i),irsml2(i)
            iN(j)=iatrs1(j) 
            do k=iatrs1(j),iatrs2(j)
               if (nmat(k)(1:2).eq.'ca') then 
                  iCa(j)=k 
               endif
               if (nmat(k)(1:1).eq.'c') then 
                  iC(j)=k 
               endif 
            enddo
!            print *,'determined phi,psi serial for residue',j,' as '
!     #           ,iphi(j),ipsi(j)
         enddo
      enddo
      abgs=300.0
      bbgs=10.0
      bgsnvar=0
      do i=1,nvr
         if (bgsposs(i)) then
            bgsnvar=bgsnvar+1
            bgsvar(bgsnvar)=i
         endif
      enddo
      end

! Checks if it is possible to perform a BGS update starting at the 
! variable indexed ipos. Calls: none.
      logical function bgsposs(ips)
      include 'INCL.H'
      include 'incl_lund.h'
      logical ians

      jv=idvr(ips)
      iaa=nursvr(jv)
      ians=.true.
!      print *,'evaluating bgs possibility for ',ips,nmvr(jv)
      if (nmvr(jv).ne.'phi') then 
!         print *,'bgs not possible because variable name is ',nmvr(jv)
         ians=.false.
      else if (iaa.gt.(irsml2(mlvr(jv))-3)) then 
!         print *,'bgs impossible, residue too close to end'
!         print *,'iaa = ',iaa,' end = ',irsml2(mlvr(jv))
         ians=.false.
      else 
         nnonfx=0 
         do i=iaa,iaa+3
            if (iphi(i).gt.0) then 
               if (.not.fxvr(iphi(i))) then
                  nnonfx=nnonfx+1
               endif
            endif
            if (.not.fxvr(ipsi(i))) then 
               nnonfx=nnonfx+1
            endif
         enddo
         if (nnonfx.lt.6) then
!            print *,iaa,'bgs impossible because ndof = ',nnonfx
            ians=.false.
         endif
      endif
      if (ians) then
!         print *,'bgs is possible for angle ',ips,jv,nmvr(jv)
      endif
      bgsposs=ians
      return
      end

! Biased Gaussian Steps. Implements a semi-local conformational update
! which modifies the protein backbone locally in a certain range of 
! amino acids. The 'down-stream' parts of the molecule outside the 
! region of update get small rigid body translations and rotations.
!
! Use the update sparingly. It is rather local, and is not of great 
! value if we need big changes in the conformation. It is recommended
! that this update be used to refine a structure around a low energy 
! candidate structure. Even at low energies, if you always
! perform BGS, the chances of coming out of that minimum are small. 
! So, there is a probability bgsprob, which decides whether BGS or the 
! normal single angle update is used. 
!
! Calls: energy, dummy (function provided as argument), addang, (rand)
!

      integer function bgs(eol1,dummy)
      include 'INCL.H'
      include 'incl_lund.h'
      external dummy
      dimension xiv(8,3),bv(8,3),rv(3,3),dv(3,8,3)
      dimension ab(8), A(8,8),p(8),ppsi(8)
      double precision ovr(mxvr)
! Initialize
!      print *,'using BGS on angle ',nmvr(idvr(ivar))
      if (bgsnvar.eq.0) then 
         bgs=0
         goto 171
      endif
      ivar=1+grnd()*bgsnvar
      do i=1,8
         iph(i)=-50000
         dph(i)=0
      enddo
      nph=0
      jv=idvr(ivar)
      ia=nursvr(jv)
! Get BGS matrices based on coordinates of atoms in 4 amino acids
      do i=1,4
         icurraa=ia+i-1
         if (iphi(icurraa).gt.0.and..not.fxvr(iphi(icurraa))) then 
            nph=nph+1
            xiv(nph,1)=xat(iCa(icurraa))
            xiv(nph,2)=yat(iCa(icurraa))
            xiv(nph,3)=zat(iCa(icurraa))
            bv(nph,1)=xiv(nph,1)-xat(iN(icurraa))
            bv(nph,2)=xiv(nph,2)-yat(iN(icurraa))
            bv(nph,3)=xiv(nph,3)-zat(iN(icurraa))
            ab(nph)=bv(nph,1)*bv(nph,1)+bv(nph,2)*bv(nph,2)
     #           +bv(nph,3)*bv(nph,3)
            iph(nph)=iphi(icurraa)
         endif 
         if (.not.fxvr(ipsi(icurraa))) then
            nph=nph+1
            xiv(nph,1)=xat(iC(icurraa))
            xiv(nph,2)=yat(iC(icurraa))
            xiv(nph,3)=zat(iC(icurraa))
            bv(nph,1)=xiv(nph,1)-xat(iCa(icurraa))
            bv(nph,2)=xiv(nph,2)-yat(iCa(icurraa))
            bv(nph,3)=xiv(nph,3)-zat(iCa(icurraa))
            ab(nph)=bv(nph,1)*bv(nph,1)+bv(nph,2)*bv(nph,2)
     #           +bv(nph,3)*bv(nph,3)
            iph(nph)=ipsi(icurraa)
         endif
      enddo
      rv(1,1)=xat(iCa(ia+3))
      rv(1,2)=yat(iCa(ia+3))
      rv(1,3)=zat(iCa(ia+3))
      rv(2,1)=xat(iC(ia+3))
      rv(2,2)=yat(iC(ia+3))
      rv(2,3)=zat(iC(ia+3))
      rv(3,1)=xat(iC(ia+3)+1)
      rv(3,2)=yat(iC(ia+3)+1)
      rv(3,3)=zat(iC(ia+3)+1)

      do i=1,3
         do j=1,nph
            dv(i,j,1)=(1.0/ab(j))*(bv(j,2)*(rv(i,3)-xiv(j,3))-
     c           bv(j,3)*(rv(i,2)-xiv(j,2)))
            dv(i,j,2)=(-1.0/ab(j))*(bv(j,1)*(rv(i,3)-xiv(j,3))-
     c           bv(j,3)*(rv(i,1)-xiv(j,1)))
            dv(i,j,3)=(1.0/ab(j))*(bv(j,1)*(rv(i,2)-xiv(j,2))-
     c           bv(j,2)*(rv(i,1)-xiv(j,1)))
         enddo
      enddo
      do i=1,nph 
         do j=i,nph
            A(i,j)=0
            do k=1,3
               do l=1,3
                  A(i,j)=A(i,j)+dv(k,i,l)*(dv(k,j,l))
               enddo
            enddo
            A(i,j)=bbgs*A(i,j)
            if (i.eq.j) then
               A(i,j)=A(i,j)+1
            endif
            A(i,j)=0.5*abgs*A(i,j)
         enddo
      enddo
      do i=1,nph
         do j=i,nph
            sum=A(i,j)
            do k=i-1,1,-1
               sum=sum-A(i,k)*A(j,k)
            enddo
            if (i.eq.j) then 
               p(i)=sqrt(sum) 
            else 
               A(j,i)=sum/p(i)
            endif
         enddo
      enddo
! Generate 8 Gaussian distributed small random angles
      do i=1,8,2
         r1=grnd()
!        In the rare event that this random number is 0, just take the next
         if (r1.le.0) r1=grnd()
         r1=sqrt(-log(r1))
         r2=grnd()
         ppsi(i)=r1*cos(pi2*r2)
         ppsi(i+1)=r1*sin(pi2*r2)
      enddo
      do i=1,nph
         dph(i)=0
      enddo
! Solve lower triangular matrix to get dphi proposals 
      do i=nph,1,-1
         sum=ppsi(i)
         do k=i+1,nph
            sum=sum-A(k,i)*dph(k)
         enddo
         dph(i)=sum/p(i)
      enddo
! Calculate intrinsic (non-Boltzmann) weight for forward process
!      print *,'calculating intrinsic weight for forward process'
      sum=0
      do i=1,nph
         sum=sum+ppsi(i)*ppsi(i)
      enddo 
      wfw=exp(-sum)
      do i=1,nph
         wfw=wfw*p(i)
      enddo

! Reconstruct chain and calculate new energy
!      print *,'going to assign changes to the chain'
      ovr = vlvr
      do i=1,nph
         vlvr(iph(i))=addang(vlvr(iph(i)),dph(i))
      enddo
      enw = energy()
! Calculate weight for reverse process for detail balance
!      print *,'proceeding to calculate weight for the reverse process'
      nph=0
      do i=1,4
         icurraa=ia+i-1
         if (iphi(icurraa).gt.0.and..not.fxvr(iphi(icurraa))) then 
            nph=nph+1
            xiv(nph,1)=xat(iCa(icurraa))
            xiv(nph,2)=yat(iCa(icurraa))
            xiv(nph,3)=zat(iCa(icurraa))
            bv(nph,1)=xiv(nph,1)-xat(iN(icurraa))
            bv(nph,2)=xiv(nph,2)-yat(iN(icurraa))
            bv(nph,3)=xiv(nph,3)-zat(iN(icurraa))
            ab(nph)=bv(nph,1)*bv(nph,1)+bv(nph,2)*bv(nph,2)
     #           +bv(nph,3)*bv(nph,3)
            iph(nph)=iphi(icurraa)
         endif 
         if (.not.fxvr(ipsi(icurraa))) then 
            nph=nph+1
            xiv(nph,1)=xat(iC(icurraa))
            xiv(nph,2)=yat(iC(icurraa))
            xiv(nph,3)=zat(iC(icurraa))
            bv(nph,1)=xiv(nph,1)-xat(iCa(icurraa))
            bv(nph,2)=xiv(nph,2)-yat(iCa(icurraa))
            bv(nph,3)=xiv(nph,3)-zat(iCa(icurraa))
            ab(nph)=bv(nph,1)*bv(nph,1)+bv(nph,2)*bv(nph,2)
     #           +bv(nph,3)*bv(nph,3)
            iph(nph)=ipsi(icurraa)
         endif
      enddo
      rv(1,1)=xat(iCa(ia+3))
      rv(1,2)=yat(iCa(ia+3))
      rv(1,3)=zat(iCa(ia+3))
      rv(2,1)=xat(iC(ia+3))
      rv(2,2)=yat(iC(ia+3))
      rv(2,3)=zat(iC(ia+3))
      rv(3,1)=xat(iC(ia+3)+1)
      rv(3,2)=yat(iC(ia+3)+1)
      rv(3,3)=zat(iC(ia+3)+1)

      do i=1,3
         do j=1,nph
            dv(i,j,1)=(1.0/ab(j))*(bv(j,2)*(rv(i,3)-xiv(j,3))-
     c           bv(j,3)*(rv(i,2)-xiv(j,2)))
            dv(i,j,2)=(-1.0/ab(j))*(bv(j,1)*(rv(i,3)-xiv(j,3))-
     c           bv(j,3)*(rv(i,1)-xiv(j,1)))
            dv(i,j,3)=(1.0/ab(j))*(bv(j,1)*(rv(i,2)-xiv(j,2))-
     c           bv(j,2)*(rv(i,1)-xiv(j,1)))
         enddo
      enddo
      do i=1,nph 
         do j=i,nph
            A(i,j)=0
            do k=1,3
               do l=1,3
                  A(i,j)=A(i,j)+dv(k,i,l)*(dv(k,j,l))
               enddo
            enddo
            A(i,j)=bbgs*A(i,j)
            if (i.eq.j) then 
               A(i,j)=A(i,j)+1
            endif
            A(i,j)=0.5*abgs*A(i,j)
         enddo
      enddo
      do i=1,nph
         do j=i,nph
            sum=A(i,j)
            do k=i-1,1,-1
               sum=sum-A(i,k)*A(j,k)
            enddo
            if (i.eq.j) then 
               p(i)=sqrt(sum) 
            else 
               A(j,i)=sum/p(i)
            endif
         enddo
      enddo
      do i=1,nph
         ppsi(i)=p(i)*dph(i)
         do j=i+1,nph
            ppsi(i)=ppsi(i)+A(j,i)*dph(j)
         enddo  
      enddo
      sum=0
      do i=1,nph
         sum=sum+ppsi(i)*ppsi(i)
      enddo
      wbw=exp(-sum) 
      do i=1,nph
         wbw=wbw*p(i)
      enddo
!     Acceptance condition (includes additional weight)
      rd=grnd()
!      print *,'generated  selection random number ',rd
!      print *,'wfw/wbw = ',wfw/wbw
      rd=-log(rd*wfw/wbw)
!      print *,'modified rd = ',rd
!      print *,'before calculating energy change'
      delta = dummy(enw)-dummy(eol1)
!      print *,'delta = ',delta
!       call outpdb(0,'after.pdb')
!      print *,'after outpdb for after.pdb'
!      do i=1,nph 
!         print *,'BGS>',i,iph(i),vlvr(iph(i)),dph(i)
!      enddo
      if (rd.ge.delta) then 
!     accept
         eol1=enw
         bgs=1
!         print *,'BGS move accepted'
      else 
!     reject
         vlvr = ovr
!          enw=energy()
!          if (abs(enw-eol1).gt.0.000001) then
!             write(*,*) 'rejected bgs move: energy change :',eol1,enw
!          endif
!         write(*,*) 'rejected bgs move: ',eol1,enw,wfw,wbw,rd
         bgs=0
      endif
 171  continue
      return
      end


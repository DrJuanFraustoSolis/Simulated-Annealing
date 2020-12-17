c**************************************************************
c
c This file contains the subroutines: zimmer
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c CALLS: none
c
c **************************************************************

      
         subroutine zimmer(nresi) 

C Calculates the Zimmerman-code of a configuration (Zimmerman et. al.
C Macromolecules, vol. 10 (1977) 1-9.)
C
C Note the difference in Notations:
C    SMMP:                      Zimmerman, et.al.:
C         A                                       A
C         B                                       B
C         C                                       C
C         D                                       D
C         E                                       E
C         F                                       F
C         G                                       G
C         H                                       H
C         a                                       A*
C         b                                       B*
C         c                                       C*
C         d                                       D*
C         e                                       E*
C         f                                       F*
C         g                                       G*
C         h                                       H*
C
      include 'INCL.H'
cf2py intent(in) nresi
      character*1 zim

     
      do j=1,nresi
         write(zimm(j:j),'(a1)') '-'
      end do

      do i=1,nvr
        iv=idvr(i)
        if(nmvr(iv).eq.'phi') then 
         xphi = vlvr(iv)*crd
         xpsi = vlvr(idvr(i+1))*crd

         if(xphi.le.-110.0d0) then
           if(xpsi.le.-140.0d0.or.xpsi.gt.110.0d0)  zim='E'
           if(xpsi.le.-90.0d0.and.xpsi.gt.-140.0d0) zim='H'
           if(xpsi.le.-40.0d0.and.xpsi.gt.-90.0d0)  zim='G'
           if(xpsi.le.20.0d0.and.xpsi.gt.-40.0d0)   zim='B'
           if(xpsi.le.110.0d0.and.xpsi.gt.20.0d0)   zim='D'
         else if(xphi.le.-40.0d0) then
           if(xpsi.le.-140.0d0.or.xpsi.gt.130.0d0)  zim='F'
           if(xpsi.le.-90.0d0.and.xpsi.gt.-140.0d0) zim='H'
           if(xpsi.le.-10.0d0.and.xpsi.gt.-90.0d0)  zim='A'
           if(xpsi.le.50.0d0.and.xpsi.gt.-10.0d0)   zim='B'
           if(xpsi.le.130.0d0.and.xpsi.gt.50.0d0)   zim='C'
         else if(xphi.le.0.0d0) then
                                                    zim='H'
         else if(xphi.le.40.0d0) then
                                                    zim='h'
         else if(xphi.le.110.0d0) then
           if(xpsi.le.-130.0d0.or.xpsi.gt.140.0d0)  zim='f'
           if(xpsi.le.140.0d0.and.xpsi.gt.90.0d0)   zim='h'
           if(xpsi.le.90.0d0.and.xpsi.gt.10.0d0)    zim='a'
           if(xpsi.le.10.0d0.and.xpsi.gt.-50.0d0)   zim='b'
           if(xpsi.le.-50.0d0.and.xpsi.gt.-130.0d0) zim='c'
         else
           if(xpsi.le.-110.0d0.or.xpsi.gt.140.0d0)  zim='e'
           if(xpsi.le.140.0d0.and.xpsi.gt.90.0d0)   zim='h'
           if(xpsi.le.90.0d0.and.xpsi.gt.40.0d0)    zim='g'
           if(xpsi.le.40.0d0.and.xpsi.gt.-20.0d0)   zim='b'
           if(xpsi.le.-20.0d0.and.xpsi.gt.-110.0d0) zim='d'
         end if
         nres=nursvr(iv)
         write(zimm(nres:nres),'(a1)') zim
        end if
      end do

      return
      end


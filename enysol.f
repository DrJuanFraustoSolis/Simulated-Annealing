c **************************************************************
c
c This file contains the subroutines: enysol,tessel
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      
      real*8 function enysol(nmol)

      include 'INCL.H'
c --------------------------------------------------------------
c
c     Double Cubic Lattice algorithm for calculating the
c     solvation energy of proteins using 
c     solvent accessible area method.
c
c     if nmol == 0 do solvation energy over all residues.
c CALLS: nursat
c
c -------------------------------------------------------------
c TODO: Check the solvent energy for multiple molecules     
      dimension numbox(mxat),inbox(mxbox+1),indsort(mxat),look(mxat)
      dimension xyz(mxinbox,3),radb(mxinbox),radb2(mxinbox)
      logical surfc(mxpoint)

!      common/ressurf/surfres(mxrs)
      
      eyslh = 0.0
      eyslp = 0.0
      if (nmol.eq.0) then
         nrslow=irsml1(1)
         nrshi=irsml2(ntlml)
      else
         nrslow = irsml1(nmol)
         nrshi  = irsml2(nmol)
      endif
      nlow = iatrs1(nrslow)
      nup = iatrs2(nrshi)
      do i=nrslow,nrshi
       surfres(i) = 0.0d0
      end do 

      numat= nup - nlow + 1

      do i=1,mxbox+1
         inbox(i)=0
      end do
     
      asa=0.0d0
      vdvol=0.0d0
      eysl=0.0d0

      avr_x=xat(nlow)
      avr_y=yat(nlow)
      avr_z=zat(nlow)
      xmin=xat(nlow)
      ymin=yat(nlow)
      zmin=zat(nlow)
      xmax=xmin
      ymax=ymin
      zmax=zmin

      rmax=rvdw(nlow)

      do j=nlow+1,nup
         if(xat(j).le.xmin) then
            xmin=xat(j)
         else if(xat(j).ge.xmax) then
            xmax=xat(j)
         end if
         avr_x=avr_x+xat(j)
         if(yat(j).le.ymin) then
            ymin=yat(j)
         else if(yat(j).ge.ymax) then
            ymax=yat(j)
         end if
         avr_y=avr_y+yat(j)
         if(zat(j).le.zmin) then
           zmin=zat(j)
         else if(zat(j).ge.zmax) then
           zmax=zat(j)
         end if
         avr_z=avr_z+zat(j)
         if(rvdw(j).ge.rmax) rmax=rvdw(j)
      end do

      avr_x=avr_x/dble(numat)
      avr_y=avr_y/dble(numat)
      avr_z=avr_z/dble(numat)
      diamax=2.d0*rmax

c  The sizes of the big box

      sizex=xmax-xmin
      sizey=ymax-ymin
      sizez=zmax-zmin

c  How many maximal diameters in each size ?

      ndx=sizex/diamax + 1
      ndy=sizey/diamax + 1
      ndz=sizez/diamax + 1

c We may need the number of quadratic boxes in (X,Y) plane

      nqxy=ndx*ndy

c   The number of cubic boxes of the size "diamax"

      ncbox=nqxy*ndz
      if(ncbox.ge.mxbox) then
       print *,'enysol> bad ncbox',ncbox
       stop
      end if
        
c Let us shift the borders to home all boxes

      shiftx=(dble(ndx)*diamax-sizex)/2.d0
      shifty=(dble(ndy)*diamax-sizey)/2.d0
      shiftz=(dble(ndz)*diamax-sizez)/2.d0
      xmin=xmin-shiftx
      ymin=ymin-shifty
      zmin=zmin-shiftz
      xmax=xmax+shiftx
      ymax=ymax+shifty
      zmax=zmax+shiftz

c Finding the box of each atom

      do j=nlow,nup
        mx=min(int(max((xat(j)-xmin)/diamax,0.0d0)),ndx-1)
        my=min(int(max((yat(j)-ymin)/diamax,0.0d0)),ndy-1)
        mz=min(int(max((zat(j)-zmin)/diamax,0.0d0)),ndz-1)
        nboxj=mx+my*ndx+mz*nqxy+1
        numbox(j)=nboxj
        if (nboxj.gt.mxbox) then
         write(*,'(a)') 'enysol> bad mxboxe-2'
         write(*,*) 'diagnostics ...'
         write(*,*) 'atom ',j
         write(*,*) 'position ',xat(j),yat(j),zat(j)
         write(*,*) 'box indices ',mx,my,mz
         write(*,*) 'resulting boxindex and limit ',nboxj,mxbox

         stop
        else
         inbox(nboxj)=inbox(nboxj)+1
        end if
      end do

c  Summation over the boxes

      do i=1,ncbox
        inbox(i+1)=inbox(i+1)+inbox(i)
      end do
         
        
c   Sorting the atoms by the their box numbers

      do i=nlow,nup
         j=numbox(i)
         jj=inbox(j)
         indsort(jj)=i
         inbox(j)=jj-1
      end do   
         
c Getting started

      do iz=0,ndz-1 ! Over the boxes along Z-axis
       do iy=0,ndy-1 !Over the boxes along Y-axis
        do ix=0,ndx-1 !Over the boxes along X-axis

           ibox=ix+iy*ndx+iz*nqxy + 1

c  Does this box contain atoms ?

           lbn=inbox(ibox+1)-inbox(ibox)

           if(lbn.gt.0) then ! There are some atoms
             nsx=max(ix-1,0)
             nsy=max(iy-1,0)
             nsz=max(iz-1,0)
             nex=min(ix+1,ndx-1)
             ney=min(iy+1,ndy-1)
             nez=min(iz+1,ndz-1)
                     
c  Atoms in the boxes around

             jcnt=1
             do  jz=nsz,nez
              do  jy=nsy,ney
               do  jx=nsx,nex
                   jbox=jx+jy*ndx+jz*nqxy+1
                   do  ii=inbox(jbox)+1,inbox(jbox+1)
                    if(rvdw(indsort(ii)).gt.0.0d0) then
                       look(jcnt)=indsort(ii)
                       jcnt=jcnt+1
                     end if
                    end do
               end do
              end do
             end do     
                             
             do  ia=inbox(ibox)+1,inbox(ibox+1)
               jbi=indsort(ia)
               trad=rvdw(jbi)
               if(trad.gt.0.0) then
                 nnei=0
                 do  ib=1,jcnt-1
                   jtk=look(ib)
                   if(jtk.ne.jbi)then
                     dx=(xat(jbi)-xat(jtk))/trad
                     dy=(yat(jbi)-yat(jtk))/trad
                     dz=(zat(jbi)-zat(jtk))/trad
                     dd=dx*dx+dy*dy+dz*dz
                     akrad=rvdw(jtk)/trad
                     dr=1.0d0+akrad
                     dr=dr*dr
cc if contact
                     if(dd.le.dr) then
                       nnei=nnei+1
                       xyz(nnei,1)=dx
                       xyz(nnei,2)=dy
                       xyz(nnei,3)=dz
                       radb(nnei)=akrad
                       radb2(nnei)=akrad*akrad
                     end if
                   end if
                 end do
cc
                  do il=1,npnt
                     surfc(il)=.false.
                  end do

c  Check overlap

                  lst=1
                  do  il=1,npnt
                   sdd = 0.0d0
                   do ilk=1,3
                     sdd = sdd +(xyz(lst,ilk)+spoint(il,ilk))**2
                   end do
                   if(sdd.gt.radb2(lst)) then
                     do  ik=1,nnei
                      sdd =0.0d0
                      do ilk=1,3
                        sdd = sdd +(xyz(ik,ilk)+spoint(il,ilk))**2
                      end do
                      if(sdd.le.radb2(ik)) then
                         lst=ik
                         go to 99
                      end if
                     end do
 99                  continue
    
                     if(ik.gt.nnei)then
                       surfc(il)=.true.
                     end if
                   end if
                  end do

                 icount=0
                 dx=0.0d0
                 dy=0.0d0
                 dz=0.0d0
                 do il=1,npnt
                  if(surfc(il)) then
                    icount=icount+1
                    dx=dx+spoint(il,1)
                    dy=dy+spoint(il,2)
                    dz=dz+spoint(il,3)
                  end if
                 end do
                 sdr=4.d0*pi*trad*trad/dble(npnt)
                 area   = sdr*dble(icount)
                 volume = sdr/3.0d0*(trad*dble(icount)
     #                               +(xat(jbi)-avr_x)*dx
     #                               +(yat(jbi)-avr_y)*dy
     #                               +(zat(jbi)-avr_z)*dz)

                 asa=asa+area
                 vdvol=vdvol+volume
                 eysl=eysl+area*sigma(jbi)
c Separate hydrophilic (h) and hyrdophobic (p) contributions to eysl
                 if (sigma(jbi).lt.0) then
                    eyslp = eyslp + area * sigma(jbi)
                    asap = asap + area
                 endif
                 if (sigma(jbi).gt.0) then
                    eyslh = eyslh + area * sigma(jbi)
                    asah = asah + area
                 endif
c Measure how much a residue is solvent accessible:
                 jres = nursat(jbi)
                 surfres(jres) = surfres(jres) + area
               end if
             end do 
           end if
          end do
         end do
       end do

c
       if (isolscl) then
          nhx=0
          mhx=0
          nbt=0
          mbt=0
          call helix(nhx,mhx,nbt,mbt)
          eysl=((nhx+nbt)*eysl)/(irsml2(ntlml)-irsml1(1))
       endif

       enysol = eysl

       return
       end
c *********************
      subroutine tessel
      include 'INCL.H'
      character lin*80

c    Skipping comment lines, which begin with '!'  

      read(20,'(a)') lin
      do while(lin(1:1).eq.'!')
        read (20,'(a)') lin
      end do

c   The first non-comment line is the number of the surface points

      read(lin(1:5),'(i5)') npnt
c        write(*,'(a,i5)') 'the number of points---->',npnt

c    Read the surface points   

      do i=1,npnt
         read(20,'(3f20.10)') spoint(i,1),spoint(i,2),spoint(i,3)
         
c        write(31,'(3f20.10)') spoint(i,1),spoint(i,2),spoint(i,3)
      end do
 
      return
 
      end


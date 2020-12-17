c **************************************************************
c
c This file contains the subroutines: dihedr,valang
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      real*8 function dihedr(i1,i2,i3,i4)

c .............................................
c  PURPOSE: return dihedral angle (i1,i2,i3,i4)
c           [in rad.]
c
c  INPUT:   i1,i2,i3,i4 - indices of four atoms
c
c  CALLS:   none
c .............................................

      include 'INCL.H'

      x1=xat(i2)-xat(i1)
      y1=yat(i2)-yat(i1)
      z1=zat(i2)-zat(i1)
      x2=xat(i3)-xat(i2)
      y2=yat(i3)-yat(i2)
      z2=zat(i3)-zat(i2)
      ux1=y1*z2-z1*y2
      uy1=z1*x2-x1*z2
      uz1=x1*y2-y1*x2
      x1=xat(i4)-xat(i3)
      y1=yat(i4)-yat(i3)
      z1=zat(i4)-zat(i3)
      ux2=z1*y2-y1*z2
      uy2=x1*z2-z1*x2
      uz2=y1*x2-x1*y2

      u1=ux1*ux1+uy1*uy1+uz1*uz1
      u2=ux2*ux2+uy2*uy2+uz2*uz2
      u=u1*u2

      if (u.ne.zero) then
        a=(ux1*ux2+uy1*uy2+uz1*uz2)/sqrt(u)
        a=max(a,-one)
        a=min(a,one)
        dihedr=acos(a)
        if (ux1*(uy2*z2-uz2*y2)+uy1*(uz2*x2-ux2*z2)+
     #      uz1*(ux2*y2-uy2*x2).lt.zero) dihedr =-dihedr
        return
      else
        write (*,'(a,4i5)')' dihedr> Error in coordinates of atoms #: '
     #                     ,i1,i2,i3,i4

        write (*,*) 'stored coordinates are xvals :',
     #       xat(i1),xat(i2),xat(i3),xat(i4) 
        write (*,*) 'yvals:', yat(i1),yat(i2),yat(i3),yat(i4) 
        write (*,*) 'zvals:', zat(i1),zat(i2),zat(i3),zat(i4) 
        call outvar(0,'crash.var')        
        stop
      endif

      end
c ************************************
      real*8 function valang(i1,i2,i3)

c .........................................
c  PURPOSE: return valence angle (i1,i2,i3)
c           [in rad.] with 'i2' as vertex
c
c  INPUT:   i1,i2,i3 - indices of 3 atoms
c
c  CALLS:   none
c .............................................

      include 'INCL.H'
      h1=xat(i2)
      h2=yat(i2)
      h3=zat(i2)
      x1=xat(i1)-h1
      x2=yat(i1)-h2
      x3=zat(i1)-h3
      y1=xat(i3)-h1
      y2=yat(i3)-h2
      y3=zat(i3)-h3

      x=x1*x1+x2*x2+x3*x3
      y=y1*y1+y2*y2+y3*y3
      u=x*y

      if (u.ne.zero) then

        a=(x1*y1+x2*y2+x3*y3)/sqrt(u)
        a=max(a,-one)
        a=min(a,one)
        valang=acos(a)
        return

      else
        write (*,'(a,3i5)')' valang> Error in coordinates of atoms #: '
     #                     ,i1,i2,i3
        write (*,*) 'stored coordinates are xvals :',
     #       xat(i1),xat(i2),xat(i3) 
        write (*,*) 'yvals:', yat(i1),yat(i2),yat(i3) 
        write (*,*) 'zvals:', zat(i1),zat(i2),zat(i3) 
        call outvar(0,'crash.var') 
        stop
      endif

      end


c **************************************************************
c
c This file contains the subroutines: difang,addang
c
c Copyright 2003       Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      real*8 function difang(a1,a2)

c ......................................................
c  PURPOSE:  difang = a2 - a1  with:  -pi < difang <= pi
c            
c  INPUT:    a1,a2-two angles [rad.]
c
c  CALLS: none
c
c ......................................................

      implicit real*8 (a-h,o-z)

      parameter (pi=3.141592653589793d0,
     #           pi2=2.d0*pi)

      d=mod((a2-a1),pi2)
      if (abs(d).le.pi) then
        difang=d
      else
        difang=d-sign(pi2,d)
      endif

      return
      end
c *********************************
      real*8 function addang(a1,a2)

c ......................................................
c  PURPOSE:  addang = a1 + a2  with:  -pi < addang <= pi
c            
c  INPUT:    a1,a2-two angles [rad.]
c
c  CALLS: none
c
c ......................................................

      implicit real*8 (a-h,o-z)

      parameter (pi=3.141592653589793d0,
     #           pi2=2.d0*pi)

      d=mod((a1+a2),pi2)
      if (abs(d).le.pi) then
        addang=d
      else
        addang=d-sign(pi2,d)
      endif

      return
      end


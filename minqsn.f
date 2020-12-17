c **************************************************************
c
c This file contains the subroutines: minqsn,mc11a,mc11e
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      subroutine minqsn(n,mxn,x,f,g,scal,acur,h,d,w,xa,ga,xb,gb,maxfun,
     #                  nfun)

c .............................................................
c  PURPOSE: Quasi-Newton minimizer
c
c           Unconstrained local minimization of function FUNC
c           vs. N variables by quasi-Newton method using BFGS-
c           formula to update hessian matrix; approximate line
c           searches performed using cubic extra-/interpolation
c           [see Gill P.E., Murray W., Wright M.H., Practical
c            Optimization, Ch. 2.2.5.7, 4.3.2.1 ff.,4.4.2.2.,
c            4.5.2.1]
c
c  INPUT:   X,F,G - variables, value of FUNC, gradient at START
c           SCAL  - factors to reduce(increase) initial step &
c                   its lower bound for line searches, diagonal
c                   elements of initial hessian matrix
c           MXN - maximal overall number of function calls
c
c  OUTPUT:  X,F,G - variables, value of FUNC, gradient at MINIMUM
c           NFUN  - overall number of function calls used
c  
c  ARRAYS:  H - approximate hessian matrix in symmetric storage
c               (dimension N(N+1)/2)
c           W,D,XA,XB,GA,GB - dimension N
c
c  CALLS:   MOVE - external to calculate function for current X
c                  and its gradients
c           MC11E- solve system H*D=-G for search direction D, where
c                  H is given in Cholesky-factorization
c           MC11A- update H using BFGS formula, factorizise new H
c                  according to Cholesky (modified to maintain its
c                  positive definiteness)
c
c  PARAMETERS:
c  
c  EPS1 - checks reduction of FUNC during line search
c         ( 0.0001 <= EPS1 < 0.5 )
c  EPS2 - controls accuracy of line search (reduce to increase
c         accuracy; EPS1 < EPS2 <= 0.9 )
c  ACUR - fractional precision for determination of variables
c         (should not be smaller than sqrt of machine accuracy)
c  TINY - prevent division by zero during cubic extrapolation
c .............................................................

      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      parameter ( eps1=0.1d0,
     #            eps2=0.7d0,
     #            tiny=1.d-32,

     #            zero=0.d0,
     #            izero=0,
     #            ione=1 )

      dimension x(mxn),g(mxn),scal(mxn),h(mxn*(mxn+1)/2),d(mxn),w(mxn),
     #          xa(mxn),ga(mxn),xb(mxn),gb(mxn)

      nfun=0
      itr=0
      dff=0.
c _______________ hessian to a diagonal matrix depending on scale
      c=0.
      do i=1,n
        c=max(c,abs(g(i)*scal(i)))
      enddo
      if (c.le.0.) c=1.

      n1=n+1
      i1=(n*n1)/2
      do i=1,i1
        h(i)=0.
      enddo

      j=1
      do i=1,n
        h(j)=.01*c/scal(i)**2
        j=j+n1-i
      enddo

    1 isfv=1    !  Re-start Search from Best Point so far

      fa=f
      do i=1,n
        xa(i)=x(i)
        ga(i)=g(i)
      enddo

    2 itr=itr+1    ! Start New Line-search from A

c ______________ search direction of the iteration
      do i=1,n
        d(i)=-ga(i)
      enddo
      call MC11E (h,n,mxn,d,w,n)

      c=0.
      dga=0.
      do i=1,n
        di=d(i)
        c=max(c,abs(di/scal(i)))
        dga=dga+ga(i)*di           ! directional derivative
      enddo
      c=1./c

      if (dga.ge.0.) goto 5      ! search is uphill

      fmin=fa
      gmin=dga
      stmin=0.

      stepub=0.        ! initial upper and
      steplb=acur*c       ! lower bound on step

c ________________________ initial step of the line search
      if (dff.gt.0.) then
        step=min(1.d0,(dff+dff)/(-dga))
      else
        step=min(1.d0,c)
      endif

    3 if (nfun.ge.maxfun) then
cc        write (*,*) ' minfor> exceeded max. number of function calls'
        return
      endif

      c=stmin+step    ! Step along Search direction A->B
      do i=1,n
        xb(i)=xa(i)+c*d(i)
      enddo

      nfun=nfun+1
      call MOVE(nfun,n,fb,xb,gb)

      isfv=min(2,isfv)

      if (fb.le.f) then
        if (fb.eq.f) then
          gl1=0.
          gl2=0.
          do i=1,n
            si=scal(i)**2
            gl1=gl1+si*g(i)**2
            gl2=gl2+si*gb(i)**2
          enddo
          if (gl2.ge.gl1) goto 4
        endif
c ______________ store function value if it is smallest so far
        f=fb
        do i=1,n
          x(i)=xb(i)
          g(i)=gb(i)
        enddo

        isfv=3
      endif

    4 dgb=0.
      do i=1,n
        dgb=dgb+gb(i)*d(i)  ! directional derivative at B
      enddo

      if (fb-fa.le.eps1*c*dga) then  ! sufficient reduction of F

        stmin=c
        fmin=fb
        gmin=dgb

        stepub=stepub-step      ! new upper bound on step

c _______________________________ next step by extrapolation
        if (stepub.gt.0.) then
          step=.5*stepub
        else
          step=9.*stmin
        endif
        c=dga+3.*dgb-4.*(fb-fa)/stmin
        if (c.gt.0.) step=min(step,stmin*max(1.d0,-dgb/c))

        if (dgb.lt.eps2*dga) goto 3  ! line minimization still not
                                     ! accurate enough -> further step
        isfv=4-isfv

        if (stmin+step.gt.steplb) then  ! line minim. complete ->
                                        ! update Hessian
          do i=1,n
            xa(i)=xb(i)
            xb(i)=ga(i)
            d(i)=gb(i)-ga(i)
            ga(i)=gb(i)
          enddo

          ir=-n
          sig=1./dga
          call MC11A (h,n,mxn,xb,sig,w,ir,ione,zero)
          ir=-ir
          sig=1./(stmin*(dgb-dga))
          call MC11A (h,n,mxn,d,sig,d,ir,izero,zero)

          if (ir.eq.n) then  ! Start new line search
            dff=fa-fb
            fa=fb
            goto 2
          else               ! rank of new matrix is deficient
            write (*,*) ' minfor> rank of hessian < number of variables'
            return
          endif

        endif

      else if (step.gt.steplb) then  ! insufficient reduction of F:
                                     ! new step by cubic interpolation
        stepub=step ! new upper bound

        c=gmin+dgb-3.*(fb-fmin)/step
        c=c+gmin-sqrt(c*c-gmin*dgb)     !! may be sqrt ( <0 )
        if (abs(c).lt.tiny) then  ! prevent division by zero
          step=.1*step
        else
          step=step*max(.1d0,gmin/c)
        endif
        goto 3     ! -> reduced step along search direction

      endif

    5 if (isfv.ge.2) goto 1   ! -> Restart from best point so far

      nfun=nfun+1
      call MOVE(nfun,n,f,x,g)

      return
      end
c ***********************************************
      subroutine mc11a(a,n,mxn,z,sig,w,ir,mk,eps)
c
c CALLS: none
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      dimension a(mxn*(mxn+1)/2),z(mxn),w(mxn)

      if (n.gt.1) then
        if (sig.eq.0..or.ir.eq.0) return
        np=n+1
        if (sig.lt.0.) then
          ti=1./sig
          ij=1
          if (mk.eq.0) then
            do i=1,n
              w(i)=z(i)
            enddo
            do i=1,n
              if (a(ij).gt.0.) then
                v=w(i)
                ti=ti+v**2/a(ij)
                if (i.lt.n) then
                  do j=i+1,n
                    ij=ij+1
                    w(j)=w(j)-v*a(ij)
                  enddo
                endif
                ij=ij+1
              else
                w(i)=0.
                ij=ij+np-i
              endif
            enddo
          else
            do i=1,n
              if (a(ij).ne.0.) ti=ti+w(i)**2/a(ij)
              ij=ij+np-i
            enddo
          endif
          if (ir.le.0) then
            ti=0.
            ir=-ir-1
          else if (ti.gt.0.) then
            if (eps.ne.0.) then
              ti=eps/sig
            else
              ir=ir-1
              ti=0.
            endif
          else if (mk.le.1) then
            goto 1
          endif
          mm=1
          tim=ti
          do i=1,n
            j=np-i
            ij=ij-i
            if (a(ij).ne.0.) tim=ti-w(j)**2/a(ij)
            w(j)=ti
            ti=tim
          enddo
          goto 2
        endif
    1   mm=0
        tim=1./sig
    2   ij=1
        do i=1,n
          v=z(i)
          if (a(ij).gt.0.) then
            al=v/a(ij)
            if (mm.le.0) then
              ti=tim+v*al
            else
              ti=w(i)
            endif
            r=ti/tim
            a(ij)=a(ij)*r
            if (r.eq.0..or.i.eq.n) goto 3
            b=al/ti
            if (r.gt.4.) then
              gm=tim/ti
              do j=i+1,n
                ij=ij+1
                y=a(ij)
                a(ij)=b*z(j)+y*gm
                z(j)=z(j)-v*y
              enddo
            else
              do j=i+1,n
                ij=ij+1
                z(j)=z(j)-v*a(ij)
                a(ij)=a(ij)+b*z(j)
              enddo
            endif
            tim=ti
            ij=ij+1
          else if (ir.gt.0.or.sig.lt.0..or.v.eq.0.) then
            ti=tim
            ij=ij+np-i
          else
            ir=1-ir
            a(ij)=v**2/tim
            if (i.eq.n) return
            do j=i+1,n
              ij=ij+1
              a(ij)=z(j)/v
            enddo
            return
          endif
        enddo
    3   if (ir.lt.0) ir=-ir
      else
        a(1)=a(1)+sig*z(1)**2
        ir=1
        if (a(1).gt.0.) return
        a(1)=0.
        ir=0
      endif

      return
      end
c ************************************
      subroutine mc11e(a,n,mxn,z,w,ir)
c
c CALLS: none
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)

      dimension a(mxn*(mxn+1)/2),z(mxn),w(mxn)

      if (ir.lt.n) return  ! rank of matrix deficient

      w(1)=z(1)
      if (n.gt.1) then
        do i=2,n
          ij=i
          i1=i-1
          v=z(i)
          do j=1,i1
            v=v-a(ij)*z(j)
            ij=ij+n-j
          enddo
          z(i)=v
          w(i)=v
        enddo
        z(n)=z(n)/a(ij)
        np=n+1
        do nip=2,n
          i=np-nip
          ii=ij-nip
          v=z(i)/a(ii)
          ip=i+1
          ij=ii
          do j=ip,n
            ii=ii+1
            v=v-a(ii)*z(j)
          enddo
          z(i)=v
        enddo
      else
        z(1)=z(1)/a(1)
      endif

      return
      end


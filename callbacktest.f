      subroutine metropolis(eol,enw,dummy)
cf2py real*8 intent(in,out) eol
cf2py real*8 intent(in,out) enw
        external dummy
        delta =  dummy(enw) - dummy(eol)
        write (*,*) delta
      end
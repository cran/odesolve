c -------- mymod.f -> mymod.so (or dll or whatever) ------
c compile with R CMD SHLIB mymod.f or (on Windows) Rcmd SHLIB mymod.f

c copied directly from the lsoda documentation,with common block added

      subroutine mymod(odeparms)
      external odeparms
      integer N
      double precision parms(3)
      common /myparms/parms

      N = 3
      call odeparms(N, parms)
      return
      end

      subroutine myderivs (neq, t, y, ydot)
      double precision t, y, ydot, parms(3)
      integer neq
      dimension y(3), ydot(3)
      common /myparms/parms

      ydot(1) = -parms(1)*y(1) + parms(2)*y(2)*y(3)
      ydot(3) = parms(3)*y(2)*y(2)
      ydot(2) = -ydot(1) - ydot(3)
      return
      end

      subroutine myjac (neq, t, y, ml, mu, pd, nrowpd)
      integer neq, ml, mu, nrowpd
      double precision y(*), pd(nrowpd,*), t, parms(3)
      common /myparms/parms

      pd(1,1) = -parms(1)
      pd(2,1) =  parms(1)
      pd(3,1) = 0.0
      pd(1,2) = parms(2)*y(3)
      pd(2,2) = -parms(2)*y(3) - 2*parms(3)*y(2)
      pd(3,2) = 2*parms(3)*y(2)
      pd(1,3) = parms(2)*y(2)
      pd(2,3) = -parms(2)*y(2)
      pd(3,3) = 0.0

      return
      end

 

c
c - rhs function
c
      subroutine  fex (neq, t, y, ydot, rpar, ipar)
      implicit none
      integer  neq
      real*8   t, y(*), ydot(*)
      integer  ipar(*)
      real*8   rpar(*)
      integer  tcgetsrc, ierr
      external tcgetsrc

c get src term from tchem library
      ierr = tcgetsrc( y, neq, ydot )

      return
      end
c
c - Jacobian
c
      subroutine  jex(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
      implicit none
      integer  neq, ml, mu, nrpd
      real*8   t, y(neq), pd(nrpd,*)
      integer  ipar(*)
      real*8   rpar(*)
      integer  tcgetjactynanl, ierr
      external tcgetjactynanl

      ierr = tcgetjactynanl( y, neq-1, pd )

      return 
      end

      SUBROUTINE GENPRM (NX, NSIZE)
      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)

      dimension nx(nsize), ny(nsize)

      external rndstart
      external rndend
      external rndunif

      call rndstart()

      n = nsize

      do 10 i=1,nsize
        ny(i) = nx(i)
 10   continue

      do 20 i=1,nsize
        j = n * rndunif() + 1
        nx(i) = ny(j)
        ny(j) = ny(n)
        n = n - 1
 20   continue

      call rndend()
      return
      end

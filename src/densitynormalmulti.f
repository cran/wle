      SUBROUTINE DENNORMU (DATI,NSIZE,NVAR,
     & DINV,den)

C NON PARAMETRIC DENSITY ESTIMATOR
C Claudio Agostinelli 
C DAIS, Universita' di Venezia, Italy
C E-mail: claudio@unive.it
C April, 24, 2015
C Version: 0.1
C
C PARAMETER:
C NAME:  I/O:    TYPE:  DIMENSIONS:   DESCRIPTIONS:
C DATI   input    D      NSIZE*NVAR   matrix of the data
C NSIZE  input    I      1            length of the data 
C NVAR   input    I      1            number of variables
C DINV   input    D      1            inverse of the var/cov bandwidth
C den    output   D      NSIZE        the kernel density

      implicit double precision(d)
      implicit integer (n,i,j,k)

      parameter(dzero=0.0d00)
      parameter(ddue=2.0d00)

      dimension den(nsize)
      dimension dati(nsize,nvar)
      dimension dinv(nvar,nvar)

      dsize=nsize
      do 100 i=1,nsize
        den(i)=dzero
        do 200 j=1,nsize
          dtemp=dzero
          do 300 k=1,nvar
            do 400 kk=1,nvar
              dtemp = dtemp + (dati(i,k)-dati(j,k))
     &           *(dati(i,k)-dati(j,kk))*dinv(k,kk)
 400        continue
 300      continue
          den(i)=den(i)+dexp(-dtemp/ddue)
 200    continue
 100  continue
      return
      end

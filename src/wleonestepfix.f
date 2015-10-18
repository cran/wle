      SUBROUTINE WLEONESTEPFIX (YDATA,XDATA,INTER,NSIZE,NCOL,
     & NVAR,DPARAM,DVAR,IRAF,RK,
     & NSTEP,rparam,varia,resid,totpesi,dpesi,
     & dden,dmod,ddelta)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     One Step WLE in the 
C     normal regression linear model
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Padova
C             35121 Padova
C             ITALIA
C
C     E-mail: claudio@stat.unipd.it
C
C     August, 2 2001
C
C     Version: 0.3
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Copyright (C) 2001 Claudio Agostinelli
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; version 2 of the License.
C 
C
C   This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PARAMETER:
C     NAME:     I/O:    TYPE:  DIMENSIONS:   DESCRIPTIONS:
C     YDATA     input    D      NSIZE        vector of the dependent variable
C     XDATA     input    D      NSIZE*NVAR   matrix of the explanatory variables
C     INTER     input    I      1            if 1 then intercept is insert in the explanatory variables
C                                            otherwise intercept must set to be 0
C     NSIZE     input    I      1            length of the data
C     NCOL      input    I      1            number of explanatory variables in the matrix XDATA
C     NVAR      input    I      1            number of first explanatory variables to be used in the calculation
C     DPARAM     input    D      NVAR+INTER   initial estimates of the parameter
C     DVAR      input    D      1            initial estimate of the residuals variance (if set to zero the variance from the residuals of the initial parameters is used)  
C     IRAF      input    I      1            type of RAF
C                                            1: Hellinger distance 
C                                            2: Negative Exponential disparity 
C                                            3: Chi squared disparity
C     RK        input    D      1            smoothing parameter
C     NSTEP      input    I      1           number of steps
C
C     rparam    output   D      (NCOL+INTER) the WLE parameters
C     varia     output   D      1            the WLE variance of the residulas
C     resid     output   D      NSIZE        the WLE residuals
C     totpesi   output   D      1            the total sum of the weights
C     dpesi      output   D      NSIZE        the weights
C     dden       output   D      NSIZE        the kernel density
C     dmod       output   D      NSIZE        the smoothed model
C     ddelta     output   D      NSIZE        the Pearson residuals

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The code use external subroutines dues to:
C
C
C                      BLAS level 1, 2, 3
C
C               Basic Linear Algebra Subprograms
C
C                      Dongarra, J.J.
C                      Du Croz, J.
C                      Hammarling, S.
C                      Hanson, R. J.
C
C                      Argonne National Laboratory
C                      9700 South Cass Avenue
C                      Argonne, Illnois 60439
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The code use external subroutines dues to:
C
C***BEGIN PROLOGUE  DQRSL
C***PURPOSE  Apply the output of DQRDC to compute coordinate transfor-
C            mations, projections, and least squares solutions.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D9, D2A1
C***TYPE      DOUBLE PRECISION (SQRSL-S, DQRSL-D, CQRSL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, ORTHOGONAL TRIANGULAR,
C             SOLVE
C***AUTHOR  Stewart, G. W., (U. of Maryland)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The code use external subroutines dues to:
C
C***BEGIN PROLOGUE  DQRDC
C***PURPOSE  Use Householder transformations to compute the QR
C            factorization of an N by P matrix.  Column pivoting is a
C            users option.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D5
C***TYPE      DOUBLE PRECISION (SQRDC-S, DQRDC-D, CQRDC-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, ORTHOGONAL TRIANGULAR,
C             QR DECOMPOSITION
C***AUTHOR  Stewart, G. W., (U. of Maryland)

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dqu=2.0d00)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     rerr: The smallest double precision number can be treated
C     on the machine as denominator in a division 
C     this value should work in any machines
C
      parameter(rerr=1.0d-40)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      dimension d(nsize),rm(nsize),delta(nsize),adelta(nsize)
      dimension ds(nsize),rw(nsize)
      dimension ydata(nsize),xdata(nsize,ncol) 
      dimension rdata(nsize)
      dimension xidata(nsize,nvar+inter)
      dimension rparam(ncol+inter) 
      dimension resid(nsize)
      dimension dpesi(nsize)
      dimension dden(nsize)
      dimension dmod(nsize)
      dimension ddelta(nsize)

      dimension dparam(nvar+inter)

      dimension xparam(nsize)
      dimension rwmat(nsize,nsize)
      dimension wxidata(nsize,nvar+inter),wydata(nsize)
      dimension dparold(nvar+inter)

C      dimension storep(nvar+inter+2)
C      dimension storew(nsize)
C      dimension storer(nsize)

      dimension work(nvar+inter),jpvt(nvar+inter),qraux(nvar+inter)
      dimension qy(nsize), qty(nsize)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C     the blas: dgemm subroutine give matrix-matrix products
C     the blas: dgemv subroutine give matrix-vector products 
      external dgemm, dgemv
C     the slatec: dqrsl subroutine give least square parameters
      external dqrsl
C      the slatec: dqrdc subroutine give the QR decomposition nedeed by the dqrsl
      external dqrdc
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      dsize=nsize
      tot=dzero

C
C     Check the presence of the intercept
C

      npre=nvar+inter
      dpre=npre

      do 10 i=1, nsize
      do 20 j=1, nvar
         xidata(i,j)=xdata(i,j)
 20   continue
      if (inter.eq.1) then
         xidata(i,npre)=duno
      endif
      do 15 ik=1, nsize
         rwmat(i,ik)=dzero
 15   continue
 10   continue

C
C     Residuals of the whole dataset from the starting values
C
      call dgemv('N',nsize,npre,duno,xidata,nsize,dparam,1,
     & dzero,xparam,1)

      do 115 i=1, nsize
         rdata(i)=ydata(i)-xparam(i)
 115  continue

C
C     Variance of the residuals
C
      if(dvar.eq.dzero) then
      var=dzero
      do 110 i=1,nsize
         var=var+(rdata(i)**ddue)         
 110  continue
      var=var/(dsize-dpre)
      else
      var=dvar
      endif

      rnows=var
      rnowh=rk*rnows
      rnowhs=(duno+rk)*rnows

C     Iteration Steps

      do 9999 istep=1,nstep

      do 125 i=1,npre
         dparold(i)=dparam(i)
 125  continue

      rolds=rnows

      do 130 i=1,nsize
         d(i)=dzero
         rm(i)=dzero

         do 120 ik=1,nsize
              ds(ik)=dexp(-((rdata(i)-rdata(ik))**dqu)/
     &               (dqu*rnowh))/dsqrt(rnowh)
              d(i)=d(i)+ds(ik)
 120     continue

         d(i)=d(i)/dsize      
         rm(i)=dexp(-((rdata(i))**dqu)/
     &         (dqu*rnowhs))/dsqrt(rnowhs)

 130  continue

      do 140 i=1,nsize
         if(rm(i).gt.rerr) then
             delta(i)=d(i)/rm(i)-duno

C Type of RAF
C IRAF = 1 Hellinger distance
C IRAF = 2 Negative Exponential disparity
C IRAF = 3 Chi Squared disparity

          if(iraf.eq.1) then   
             adelta(i)=dqu*(dsqrt(delta(i)+duno)-duno)
          endif
          if(iraf.eq.2) then
             adelta(i)=dqu - (dqu+delta(i))*dexp(-delta(i))
          endif   
          if(iraf.ne.3) then
             rw(i)=(adelta(i)+duno)/(delta(i)+duno)
          else
             rw(i)=duno - (delta(i) / (delta(i) + ddue))**ddue             
          endif

             if(rw(i).lt.dzero) then 
                rw(i)=dzero
             endif
             if(rw(i).gt.duno) then
                rw(i)=duno
             endif
             
         else
             rw(i)=dzero
             delta(i)=duno/rerr
         endif         
 140  continue

      tot=dzero

      do 150 i=1,nsize
         rwmat(i,i)=dsqrt(rw(i))
         tot=tot+rw(i)
C         write(*,*) rw(i)
 150  continue   

      tot=tot/dsize

C
C     Weighted Least square parameters from the whole sample
C
C
C     Parameters
C

      call dgemv('N',nsize,nsize,duno,rwmat,nsize,ydata,1,
     & dzero,wydata,1)

      call dgemm('N','N',nsize,npre,nsize,duno,rwmat,nsize,
     & xidata,nsize,dzero,wxidata,nsize) 

      call dqrdc (wxidata,nsize,nsize,npre,qraux,jpvt,work,0)

      call dqrsl (wxidata,nsize,nsize,npre,qraux,wydata,qy,qty,     
     & dparam,rdata,xparam,00111,info)

      call dgemv('N',nsize,npre,duno,xidata,nsize,dparam,1,
     & dzero,xparam,1)

      if(info.ne.0) then 
         go to 8888
      endif   

C         write(*,*) info

C      do 160 i=1,npre
C         write(*,*) dparam(i)
C 160  continue

C
C     Residuals in the whole dataset
C
      do 165 i=1, nsize
         rdata(i)=ydata(i)-xparam(i)
 165  continue
      
C
C     Variance of the residuals
C
      var=dzero

      do 167 i=1,nsize
         var=var+(rw(i)*(rdata(i)**ddue))
 167  continue

      var=var/((tot*dsize)-dpre)

      rnows=var
      rnowh=rk * rnows
      rnowhs=(duno+rk)*rnows

 9999 continue

C Nstep done

 8888 continue

C write down the results and return

         do 205 i=1,npre
            rparam(i)=dparam(i)
 205     continue
         do 207 i=npre+1,ncol
            rparam(i)=dzero
 207     continue
         varia=var
         totpesi=tot
CCCC 210  continue

         do 230 i=1,nsize
            dpesi(i)=rw(i)
            dden(i)=d(i)
            dmod(i)=rm(i)
            ddelta(i)=delta(i)
            resid(i)=rdata(i)
 230     continue

      return
      end

      SUBROUTINE WLS (YDATA,XDATA,DPESI,NSIZE,NPRE,NCOL,NTYPE,
     & dparam,iinfo)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Weighted Least Square
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Padova
C             35121 Padova
C             ITALIA
C
C     E-mail: claudio@stat.unipd.it
C
C     October, 10 1999
C
C     Version: 0.2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Copyright (C) 1999 Claudio Agostinelli
C
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
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
C     YDATA     input    D     NSIZE         vector of the dependent variable
C     XDATA     input    D     NSIZE*NPRE    matrix of the explanatory variables
C     DPESI     input    D     NSIZE         vector of the weights
C     NSIZE     input    I     1             number of observations
C     NPRE      input    I     1             number of parameter
C     NCOL      input    I     1             number of the first explanatory variables to be used
C     NTYPE     input    I     1             if 1 the dpesi are used.
C     dparam    output   D     NPRE          vector of the parameter 
C     iinfo     output   I     1
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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

      dimension ydata(nsize)
      dimension wydata(nsize)
      dimension xdata(nsize,npre)
      dimension wxdata(nsize,ncol)
      dimension dparam(npre)
      dimension tparam(ncol) 
      dimension dtemp(nsize)
      dimension dfit(nsize)
      dimension rwmat(nsize,nsize)
      dimension dpesi(nsize)

      dimension work(npre),jpvt(npre),qraux(npre)
      dimension qy(nsize), qty(nsize)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C     the slatec: dqrsl subroutine give least square parameters
      external dqrsl
C      the slatec: dqrdc subroutine give the QR decomposition nedeed by the dqrsl
      external dqrdc
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC



      do 610 i=1,npre
         dparam(i)=dzero
 610  continue

      do 630 j=1,nsize
         wydata(j)=ydata(j)
         do 620 i=1,ncol
            wxdata(j,i)=xdata(j,i)
 620     continue
 630  continue

      if(ntype.eq.1) then

      do 750 i=1,nsize
         do 751 j=1,nsize
            if(i.eq.j) then
               rwmat(i,j)=dsqrt(dpesi(i))
            else
               rwmat(i,j)=dzero
            endif
 751     continue
 750  continue

C
C
C     The weighted observations
C
C
      call dgemv('N',nsize,nsize,duno,rwmat,nsize,ydata,1,
     & dzero,wydata,1)

      call dgemm('N','N',nsize,ncol,nsize,duno,rwmat,nsize,
     & xdata,nsize,dzero,wxdata,nsize)

      endif

      call dqrdc (wxdata,nsize,nsize,ncol,qraux,jpvt,work,0)

      call dqrsl (wxdata,nsize,nsize,ncol,qraux,wydata,qy,qty,
     & tparam,dtemp,dfit,00111,iinfo)

      do 640 i=1,ncol
         dparam(i)=tparam(i)
 640  continue      

      return
      end

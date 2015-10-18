      SUBROUTINE MLECP (YDATA,XDATA,INTER,NSIZE,NVAR,
     & NREP,
     & VAR,
     & cp,param,varia,resid,iinfo)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Mallows Cp based on 
C     Least Squares in the normal regression linear model
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
C     YDATA     input    D      NSIZE        vector of the dependent variable
C     XDATA     input    D      NSIZE*NVAR   matrix of the explanatory variables
C     INTER     input    I      1            if 1 then intercept is insert in the explanatory variables
C                                            otherwise intercept must set to be 0
C     NSIZE     input    I      1            length of the data
C     NVAR      input    I      1            number of explanatory variables
C     NREP      input    I      1            number of model to be reported
C     VAR       input    D      1            variance to be used in the evaluation of the wcp, 
C                                            if zero a full model residual variance is used 
C
C     cp        output   D      NREP*(NVAR+INTER+1) 
C     param     output   D      NREP*(NVAR+INTER) the LS parameters
C     varia     output   D      NREP the LS variance of the residuals
C     resid     output   D      NREP*NSIZE the LS residuals for each model
C
C     iinfo      output   I      1
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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

      dimension ydata(nsize),xdata(nsize,nvar) 
      dimension xidata(nsize,nvar+inter), xiidata(nsize,nvar+inter)
      dimension param(nrep,nvar+inter), xparam(nsize)
      dimension wparam(nvar+inter)
      dimension varia(nrep) 
      dimension resid(nrep,nsize),wresid(nsize) 
      dimension nmodel(nvar+inter)
      dimension cp(nrep,nvar+inter+1)

      dimension work(nvar+inter),jpvt(nvar+inter),qraux(nvar+inter)
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

      dsize=nsize
      iinfo=0
      dvar = 0.0D+00
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
 10   continue
      
C From the biggest model to the smallest one.

C      nmaxmod=(2**npre)-1
       nmaxmod=nrep

      do 900 imodel=nmaxmod,1,-1

         call dectobin (imodel,npre,nmodel)
         
C         write(*,*) nmodel

         ncol=0

         do 30 i=1,npre
            if (nmodel(i).eq.1) then
            ncol=ncol+1
               do 40 j=1,nsize
                  xiidata(j,ncol)=xidata(j,i)
 40            continue
            endif
 30      continue


      call dqrdc (xiidata,nsize,nsize,ncol,qraux,jpvt,work,0)

      call dqrsl (xiidata,nsize,nsize,ncol,qraux,ydata,qy,qty,     
     & wparam,wresid,xparam,00111,info)

      if(info.ne.0) then
         iinfo=1
         return
      endif   

      wvaria=dzero
      do 35 i=1,nsize
         wvaria=wvaria+wresid(i)**ddue
 35   continue

      dcol=ncol

      wvaria=wvaria/(dsize-dcol)

       if (ncol.eq.npre) then
          dvar=wvaria
      endif    

      if (var.ne.dzero) then
         dvar=var
      endif   

         ipos=0
      do 60 i=1,npre
         cp(imodel,i)=nmodel(i)
C         write(*,*) cp(imodel,i)

         if (nmodel(i).eq.1) then
           ipos=ipos+1
           param(imodel,i)=wparam(ipos)
         else
           param(imodel,i)=dzero 
         endif

 60   continue

      do 100 i=1,nsize
         resid(imodel,i)=wresid(i)
C         write(*,*) imodel,i
 100  continue

      wsum=dzero
      do 80 isize=1,nsize
         wsum=wsum+(wresid(isize)**ddue)
 80   continue

      cp(imodel,npre+1)=
     & wsum/dvar - dsize + (2*ncol) 
C      write(*,*) cp(imodel,npre+1)
      varia(imodel)=wvaria
CCC 70   continue NOT USED 20151017

 900  continue
      return
      end

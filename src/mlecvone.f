      SUBROUTINE MLECVONE (YDATA,XDATA,INTER,NSIZE,NVAR,
     &  NMCCV,IMODEL,NSPLIT,ISEED,
     &  cv,info)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Cross-Validation 
C     in the normal regression linear model
C     for one model
C
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Venezia
C             30125 VENEZIA
C             ITALIA
C
C     E-mail: claudio@unive.it
C
C     December, 04 2003
C
C     Version: 0.3
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Copyright (C) 2003 Claudio Agostinelli
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
C     NMCCV     input    I      1            number of MonteCarlo replication for the Cross-Validation
C     NSPLIT    input    I      1            dimension of the split for the costruction-validation sample size
C
C     cv        output   D      nvar+inter+1 
C     info      output   I      1
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C The code use an external subroutine dues to:
C
C
C
C                                     RANLIBF
C
C            Library of Fortran Routines for Random Number Generation
C
C                            Compiled and Written by:
C
C                                 Barry W. Brown
C                                  James Lovato                             C
C
C                     Department of Biomathematics, Box 237
C                     The University of Texas, M.D. Anderson Cancer Center
C                     1515 Holcombe Boulevard
C                     Houston, TX      77030
C
C
C This work was supported by grant CA-16672 from the National Cancer Institute.
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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

      parameter(rerr=1.0d-20)

      dimension ydata(nsize),xdata(nsize,nvar) 
      dimension xidata(nsize,nvar+inter) 
      dimension nmodel(nvar+inter)
      dimension cv(nvar+inter+1)

      dimension dparam(nvar+inter)
      dimension tparam(nvar+inter)

      dimension nstart(nsize)
      dimension wycsub(nsplit), wxcsub(nsplit,nvar+inter)
      dimension wyvsub(nsize-nsplit) 
      dimension wxvsub(nsize-nsplit,nvar+inter) 
      dimension xparam(nsize-nsplit)
      dimension wxdata(nsplit,nvar+inter)
      dimension ddd(nsplit)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C     the ranlib: genprm subroutine generate a 
      external genprm
C     the blas: dgemm subroutine give matrix-matrix products
C     the blas: dgemv subroutine give matrix-vector products 
      external dgemm, dgemv
C     the slatec: dqrsl subroutine give least square parameters
      external dqrsl
C      the slatec: dqrdc subroutine give the QR decomposition nedeed by the dqrsl
      external dqrdc
C
      external setall
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call setall(iseed,iseed)      

      dsize=nsize
      info=0
      iinfo=0
      dcv=dzero
      nvalid=nsize-nsplit
      dvalid=nvalid
      dconv=dzero
      dmccv=nmccv

      do 5 i=1,nsize
         nstart(i)=i
 5    continue

C
C     Check the presence of the intercept
C
      npre=nvar+inter
      dpre=npre

      call dectobin (imodel,npre,nmodel)

      do 10 i=1, nsize
      do 20 j=1, nvar   
         xidata(i,j)=xdata(i,j)
 20   continue
      if (inter.eq.1) then
         xidata(i,npre)=duno
      endif
 10   continue

      do 9999 imc=1,nmccv
C
C     Costruction-Validation Sample
C

      call genprm(nstart,nsize)

C      write(*,*) nstart


      do 110 i=1,nsplit
         wycsub(i)=ydata(nstart(i))
      do 120 j=1,npre   
         wxcsub(i,j)=xidata(nstart(i),j)
 120  continue
 110  continue

      do 130 i=nsplit+1,nsize
         wyvsub(i-nsplit)=ydata(nstart(i))
      do 140 j=1,npre   
         wxvsub(i-nsplit,j)=xidata(nstart(i),j)
 140     continue
 130  continue

CCC         write(*,*) nmodel

         ncol=0
         do 30 i=1,npre
            if (nmodel(i).eq.1) then
            ncol=ncol+1
               do 40 j=1,nsplit
                  wxdata(j,ncol)=wxcsub(j,i)
 40            continue
            endif
 30      continue

      call wls(wycsub,wxdata,ddd,nsplit,npre,ncol,0,tparam,iinfo)

      if (iinfo.ne.0) then
         dconv=dconv+duno
         info=3
      else   
      ipos=0
      do 200 i=1,npre
         if (nmodel(i).eq.1) then
           ipos=ipos+1
           dparam(i)=tparam(ipos)
         else
           dparam(i)=dzero 
         endif
 200  continue

      call dgemv('N',nvalid,npre,duno,wxvsub,nvalid,dparam,1,
     & dzero,xparam,1)
           
      dfit=dzero
      do 100 i=1,nvalid
         if(abs(wyvsub(i)-xparam(i)).gt.rerr) then
            dfit=dfit+(wyvsub(i)-xparam(i))**ddue
         endif
 100  continue
         dfit=dfit/nvalid   
         dcv=dcv+dfit
      endif

C     End of the Monte Carlo replications 
 9999 continue

      dcv=dcv/(dmccv-dconv)  
      do 898 i=1,npre
         cv(i)=nmodel(i)
 898  continue
      cv(npre+1)=dcv

      return
      end


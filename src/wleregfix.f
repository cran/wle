      SUBROUTINE WLEREGFIX (YDATA,XDATA,INTER,NSIZE,NCOL,
     & NVAR,NBOOT,NGRP,NREP,IRAF,RK,
     & RPREC,REQUAL,IMAX,rparam,varia,resid,totpesi,dpesi,
     & dden,dmod,ddelta,nsame,nsol,nconv)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Bootstrap roots search for WLE in the 
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
C     August, 2, 2001
C
C     Version: 0.4
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
C     NBOOT     input    I      1            number of bootstrap replication
C     NGRP      input    I      1            dimension of the subsample 
C     NREP      input    I      1            number of solution be reported
C     IRAF      input    I      1            type of RAF
C                                            1: Hellinger distance 
C                                            2: Negative Exponential disparity 
C                                            3: Chi squared disparity
C     RK        input    D      1            smoothing parameter
C     RPREC     input    D      1            precision of the convergence 
C                                            in absolute value
C     REQUAL    input    D      1            when two roots are said equal,
C                                            must be less than RPREC 
C     IMAX      input    I      1            maximum number of iterations for each starting points
C
C     rparam    output   D      NREP*(NCOL+INTER) the WLE parameters
C     varia     output   D      NREP         the WLE variance of the residulas
C     resid     output   D      NREP*NSIZE   the WLE residuals for each roots
C     totpesi   output   D      NREP         the total sum of the weights
C     dpesi     output   D      NREP*NSIZE   the weights
C     dden      output   D      NREP*NSIZE   the residuals density
C     dmod      output   D      NREP*NSIZE   the smoothed model
C     ddelta    output   D      NREP*NSIZE   the Pearson residuals
C     nsame     output   I      NREP         frequencies of each root
C     nsol      output   I      1            the total number of solutions
C     nconv     output   I      1            number of boostrap sampling that 
C                                            does not converge
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
      dimension rdata(nsize),rsub(ngrp)
      dimension xidata(nsize,nvar+inter)
      dimension rparam(nrep,ncol+inter), varia(nrep) 
      dimension resid(nrep,nsize), totpesi(nrep)
      dimension dpesi(nrep,nsize), nsame(nrep)
      dimension dden(nrep,nsize), dmod(nrep,nsize)
      dimension ddelta(nrep,nsize) 
      dimension dparam(nvar+inter)

      dimension xparam(nsize),xsubpar(ngrp)
      dimension rwmat(nsize,nsize)
      dimension wxidata(nsize,nvar+inter),wydata(nsize)
      dimension dparold(nvar+inter)

      dimension ysub(ngrp), xsub(ngrp,nvar+inter) 
      dimension nstart(nsize)
      dimension nrand(nboot,ngrp)

      dimension work(nvar+inter),jpvt(nvar+inter)
      dimension qraux(nvar+inter)
      dimension qy(nsize), qty(nsize)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C     the ranlib: genprm subroutine generate a permuation of an array 
      external genprm
C     the blas: dgemm subroutine give matrix-matrix products
C     the blas: dgemv subroutine give matrix-vector products 
      external dgemm, dgemv
C     the slatec: dqrsl subroutine give least square parameters
      external dqrsl
C      the slatec: dqrdc subroutine give the QR decomposition nedeed by the dqrsl
      external dqrdc
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      dgrp=ngrp
      dsize=nsize
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

      do 30 i=1, nrep
         nsame(i)=0
 30   continue

      do 50 i=1,nsize
         nstart(i)=i
 50   continue

      nconv=0 
      nsol=0

C Start bootstrapping

      do 900 iboot=1,nboot

C Initial values 

7777  call genprm(nstart,nsize)

      nagain=0

      do 60 icheck=2,iboot
      do 70 istart=1,ngrp
      do 80 iistart=1,ngrp
   	 if(nstart(istart).eq.nrand(icheck,iistart)) then
             nagain=nagain+1
         endif
 80   continue   
 70   continue
         if(nagain.eq.ngrp) then
             goto 7777
         else
            nagain=0
         endif
 60   continue

      do 90 i=1,ngrp
         ysub(i)=ydata(nstart(i))
      do 95 j=1,npre   
         xsub(i,j)=xidata(nstart(i),j)
 95   continue
         nrand(iboot,i)=nstart(i)
 90   continue

C
C     Least square parameters from the subsample
C
C
C     Parameters from the starting values
C

      call dqrdc (xsub,ngrp,ngrp,npre,qraux,jpvt,work,0)

      call dqrsl (xsub,ngrp,ngrp,npre,qraux,ysub,qy,qty,     
     & dparam,rsub,xsubpar,00111,info)

      if (info.ne.0) then
         nconv=nconv+1
         go to 7777
      endif   

C      do 100 j=1, npre
C         write(*,*) dparam(j)
C 100  continue

C
C     Variance of the residuals
C

      var=dzero

      do 110 i=1,ngrp
         var=var+(rsub(i)**ddue)         
 110  continue

      var=var/(dgrp-dpre)

C
C     Residuals of the whole dataset from the starting values
C
      call dgemv('N',nsize,npre,duno,xidata,nsize,dparam,1,
     & dzero,xparam,1)

      dmedia=dzero
      do 115 i=1, nsize
         rdata(i)=ydata(i)-xparam(i)
         dmedia=dmedia+rdata(i)
 115  continue
 
      dmedia=dmedia/dsize

      rnows=var
      rnowh=rk*rnows
      rnowhs=(duno+rk)*rnows
      iter=0

C     Iteration Steps until convergence achieved

 9999 continue

      iter=iter+1

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
	 rm(i)=dexp(-((rdata(i)-dmedia)**dqu)/
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
         nconv=nconv+1
         go to 8888
      endif   

C         write(*,*) info

C      do 160 i=1,npre
C         write(*,*) dparam(i)
C 160  continue

C
C     Residuals in the whole dataset
C
      dmedia=dzero
      do 165 i=1, nsize
         rdata(i)=ydata(i)-xparam(i)
         dmedia=dmedia+rw(i)*rdata(i)
 165  continue
      dmedia=dmedia/(tot*dsize)

C
C     Variance of the residuals
C

      var=dzero

      do 167 i=1,nsize
         var=var+(rw(i)*((rdata(i)-dmedia)**ddue))
 167  continue

      var=var/((tot*dsize)-dpre)

      rnows=var
      rnowh=rk * rnows
      rnowhs=(duno+rk)*rnows

      diffp=dzero

      do 170 i=1,npre
         diff=abs(dparam(i)-dparold(i))
         if (diff.gt.diffp) then
            diffp=diff
         endif   
 170  continue

      diffs=abs(rnows-rolds)

      if(iter.gt.imax) then 
         nconv=nconv+1
         goto 8888 
      endif  

      if(diffp.gt.rprec.or.diffs.gt.rprec) then
         goto 9999
      endif   

C Convergence achieved

C
C Is this a new root, then store it
C
      if(nsol.eq.0) then
         nsol=nsol+1
         nsame(1)=1 
         do 175 i=1,npre 
            rparam(nsol,i)=dparam(i)
 175     continue
         varia(nsol)=rnows
         totpesi(nsol)=tot
         do 180 i=1,nsize
            dpesi(nsol,i)=rw(i)
            dden(nsol,i)=d(i)
            dmod(nsol,i)=rm(i)
            ddelta(nsol,i)=delta(i)
            resid(nsol,i)=rdata(i)
 180     continue
      else
         do 185 isol=1,nsol
            diffp=dzero
            do 187 i=1,npre
               diff=abs(dparam(i)-rparam(isol,i))
               if (diff.gt.diffp) then
                  diffp=diff
               endif   
 187        continue
            diffs=abs(rnows-varia(isol))
            if(diffp.lt.requal.and.diffs.lt.requal
     &          ) then
                   nsame(isol)=nsame(isol)+1
                   goto 190
            endif
 185     continue
                   nsol=nsol+1
		   nsame(nsol)=1
                   do 195 i=1,npre 
                      rparam(nsol,i)=dparam(i)
 195               continue
                   varia(nsol)=rnows
                   totpesi(nsol)=tot
                   do 200 i=1,nsize
                      dpesi(nsol,i)=rw(i)
                      dden(nsol,i)=d(i)
                      dmod(nsol,i)=rm(i)
                      ddelta(nsol,i)=delta(i)
                      resid(nsol,i)=rdata(i)
 200               continue
 190               continue
      endif 

      if(nsol.eq.nrep) goto 6666

8888  continue 
C       
C end bootstrapping
C

 900  continue

 6666 continue

      return
      end




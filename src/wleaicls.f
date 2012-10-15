CCCCC NON TERMINATA !!!!!

      SUBROUTINE WLEAICLS (YDATA,XDATA,INTER,NSIZE,NVAR,
     & NBOOT,NGRP,NREP,IRAF,RK,
     & RPREC,REQUAL,IMAX,VAR,NMAXSOL,DMINPESI,ALPHA,
     & cp,param,varia,resid,totpesi,pesi,nsame,info)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Weighted Akaike Information Criterion based on 
C     WLE in the normal regression linear model
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Padova
C             35121 Padova
C             ITALIA
C
C     E-mail: claudio@unive.it
C
C     September, 20 2010
C
C     Version: 0.1
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Copyright (C) 2010 Claudio Agostinelli
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
C     NVAR      input    I      1            number of explanatory variables
C     NBOOT     input    I      1            number of bootstrap replication
C     NGRP      input    I      1            dimension of the subsample 
C     NREP      input    I      1            number of model to be reported
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
C     VAR       input    D      1            variance to be used in the evaluation of the waic, 
C                                            if zero a full model residual variance is used 
C     NMAXSOL   input    I      1            Maximum number of solutions to be considered for the full model 
C     DMINPESI  input    D      1            Minimum percentage of the weights for which a roots is considered
C     ALPHA     input    D      1            the coefficient of the number of the parameters.
C
C
C     cp        output   D      NREP*(NVAR+INTER+1) the aic
C     param     output   D      NREP*(NVAR+INTER) the WLE parameters
C     varia     output   D      NREP the WLE variance of the residuals
C     resid     output   D      NREP*NSIZE the WLE residuals for each model
C     totpesi   output   D      NREP the total sum of the weights
C     pesi      output   D      NMAXSOL*NSIZE the weights
C     nsame     output   I      NMAXSOL frequencies of each root
C     info      output   I      1
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     nsol               I      NMAXSOL         the total number of solutions
C     nconv              I      NMAXSOL         number of boostrap sampling that 
C                                            does not converg

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dqu=2.0d00)
      parameter(dpi=3.141592653d00)

      dimension ydata(nsize),xdata(nsize,nvar) 
      dimension xidata(nsize,nvar+inter), xiidata(nsize,nvar+inter)
      dimension param(nrep,nvar+inter)
      dimension dparam(nvar+inter), tparam(nvar+inter)
      dimension xparam(nsize)
      dimension rwmat(nsize,nsize)
      dimension wxidata(nsize,nvar+inter),wydata(nsize)
      dimension wparam(nmaxsol,nvar+inter),wvaria(nmaxsol)
      dimension wtotpesi(nmaxsol), wresid(nmaxsol,nsize) 
      dimension wpesi(nmaxsol,nsize)
      dimension dden(nmaxsol,nsize), dmod(nmaxsol,nsize)
      dimension ddelta(nmaxsol,nsize), nwsame(nmaxsol) 
      dimension varia(nrep)
      dimension resid(nrep,nsize)
      dimension nmodel(nvar+inter)
      dimension cp(nrep,nvar+inter+1)
      dimension dpesi(nsize)
      dimension ddd(nsize)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C     wleregfix: subroutine for evaluating a wle estimators
C                for the normal regression model 
      external wleregfix
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
      info=0
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

      do 15 imodel=nmaxmod,1,-1

         do 25 i=1,npre
           cp(imodel,i)=nmodel(i)
 25      continue
         cp(imodel,npre+1)=dzero
 15   continue

      call wleregfix(ydata,xidata,0,nsize,npre,
     &  npre,nboot,ngrp,nmaxsol,iraf,rk,rprec,requal,imax,
     &  wparam,wvaria,wresid,wtotpesi,wpesi,
     &  dden,dmod,ddelta,nwsame,nsol,iconv)

C      write(*,*) 'sono arrivato qui 1'

      if(iconv.eq.nboot) then
         info=1
         return
      endif

      indice=0
      dvar=wvaria(1)+duno
      do 30 i=1,nsol
         if(dvar.gt.wvaria(i).and.dminpesi
     &      .lt.wtotpesi(i)) then
            dvar=wvaria(i)
            indice=i
         endif
 30   continue
      if(indice.eq.0) then
         info=2
         return
      endif
      dtotpesi=wtotpesi(indice)
      do 40 isize=1,nsize
         dpesi(isize)=wpesi(indice,isize)
 40   continue

      do 50 i=1,nsize
         rwmat(i,i)=dsqrt(dpesi(i))
 50   continue

C      write(*,*) 'sono arrivato qui 2'

C
C     The weighted observations
C

      call dgemv('N',nsize,nsize,duno,rwmat,nsize,ydata,1,
     & dzero,wydata,1)

      call dgemm('N','N',nsize,npre,nsize,duno,rwmat,nsize,
     & xidata,nsize,dzero,wxidata,nsize)

C From the biggest model to the smallest one.

C      nmaxmod=(2**npre)-1
       nmaxmod=nrep

      do 900 imodel=nmaxmod,1,-1

         call dectobin (imodel,npre,nmodel)
         
C         write(*,*) nmodel

         ncol=0

         do 60 i=1,npre
            if (nmodel(i).eq.1) then
            ncol=ncol+1
               do 70 j=1,nsize
                  xiidata(j,ncol)=wxidata(j,i)
 70            continue
            endif
 60      continue


      call wls(wydata,xiidata,ddd,nsize,npre,ncol,0,tparam,iinfo)

      ipos=0
      do 80 i=1,npre
         if (nmodel(i).eq.1) then
           ipos=ipos+1
           dparam(i)=tparam(ipos)
         else
           dparam(i)=dzero 
         endif
 80   continue

C      write(*,*) 'sono arrivato qui 9'

      call dgemv('N',nsize,npre,duno,wxidata,nsize,dparam,1,
     & dzero,xparam,1)
           
      do 90 i=1,nsize
        resid(imodel,i) = wydata(i)-xparam(i)
 90   continue


      wsum=dzero
      do 120 isize=1,nsize
         wsum=wsum+dpesi(isize)*
     &   (resid(imodel,isize)**ddue)
 120      continue

         cp(imodel,npre+1)=
     &   wsum/dvar + (dtotpesi*dsize)*
     &   dlog(ddue*dpi*dvar)+  
     &   (alpha*ncol)
C      write(*,*) cp(nmaxsol*(imodel-1)+isol,npre+1)

C      varia(nmaxsol*(imodel-1)+isol)=wvaria(isol)
C      nsame(nmaxsol*(imodel-1)+isol)=nwsame(isol)
C      totpesi(nmaxsol*(imodel-1)+isol)=wtotpesi(isol)
 
 900  continue
      return
      end
























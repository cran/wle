      SUBROUTINE WLEAIC (YDATA,XDATA,INTER,NSIZE,NVAR,
     & NBOOT,NGRP,NREP,IRAF,RK,
     & RPREC,REQUAL,IMAX,VAR,NMAXSOL,DMINPESI,IPESI,ALPHA,
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
C     VAR       input    D      1            variance to be used in the evaluation of the wcp, 
C                                            if zero a full model residual variance is used 
C     NMAXSOL   input    I      1            Maximum number of solutions to be considered for each model 
C     DMINPESI  input    D      1            Minimum percentage of the weights for which a roots is considered
C     IPESI     input    I      1            weights by full model or by reduced model  
C     ALPHA     input    D      1            the coefficient of the number of the parameters.
C
C
C     cp        output   D      (NREP*NMAXSOL)*(NVAR+INTER+1) 
C     param     output   D      (NREP*NMAXSOL)*(NVAR+INTER) the WLE parameters
C     varia     output   D      (NREP*NMAXSOL) the WLE variance of the residuals
C     resid     output   D      (NREP*NMAXSOL)*NSIZE the WLE residuals for each model
C     totpesi   output   D      NREP*NMAXSOL the total sum of the weights
C     pesi      output   D      (NREP*NMAXSOL)*NSIZE the weights
C     nsame     output   I      (NREP*NMAXSOL) frequencies of each root
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
      dimension param(nrep*nmaxsol,nvar+inter)
      dimension wparam(nmaxsol,nvar+inter)
      dimension varia(nrep*nmaxsol),wvaria(nmaxsol) 
      dimension resid(nrep*nmaxsol,nsize),wresid(nmaxsol,nsize) 
      dimension totpesi(nrep*nmaxsol),wtotpesi(nmaxsol)
      dimension pesi(nrep*nmaxsol,nsize),wpesi(nmaxsol,nsize)
      dimension dden(nmaxsol,nsize), dmod(nmaxsol,nsize)
      dimension ddelta(nmaxsol,nsize) 
      dimension nsame(nrep*nmaxsol), nwsame(nmaxsol)
      dimension nmodel(nvar+inter)
      dimension cp(nrep*nmaxsol,nvar+inter+1)
      dimension dpesi(nsize)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C     wleregfix: subroutine for evaluating a wle estimators
C                for the normal location model 
      external wleregfix
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

       call wleregfix(ydata,xiidata,0,nsize,npre,
     & ncol,nboot,ngrp,nmaxsol,iraf,rk,rprec,requal,imax,
     & wparam,wvaria,wresid,wtotpesi,wpesi,
     & dden,dmod,ddelta,nwsame,nsol,nconv)

       if(nconv.eq.nboot) then
          info=1
       endif

       if (ncol.eq.npre) then
          indice=0
          dvar=wvaria(1)+duno
C         write(*,*) nsol
          do 50 i=1,nsol
             if(dvar.gt.wvaria(i).and.dminpesi
     &            .lt.wtotpesi(i)) then
                dvar=wvaria(i)
                indice=i
C                write(*,*) wtotpesi
             endif
 50       continue
          if(indice.eq.0) then
             info=2
             return
          endif
          dtotpesi=wtotpesi(indice)
          do 55 isize=1,nsize
             dpesi(isize)=wpesi(indice,isize)
 55       continue
      endif    

      if (var.ne.dzero) then
         dvar=var
      endif   

      if (nsol.gt.0) then 

      do 70 isol=1,nsol
         ipos=0
      do 60 i=1,npre
         cp(nmaxsol*(imodel-1)+isol,i)=nmodel(i)
C         write(*,*) cp(nmaxsol*(imodel-1)+isol,i)

         if (nmodel(i).eq.1) then
           ipos=ipos+1
           param(nmaxsol*(imodel-1)+isol,i)=wparam(isol,ipos)
         else
           param(nmaxsol*(imodel-1)+isol,i)=dzero 
         endif

 60   continue

      do 100 i=1,nsize
         resid(nmaxsol*(imodel-1)+isol,i)=wresid(isol,i)
C         write(*,*) nmaxsol*(imodel-1)+isol,i,nmaxsol
         pesi(nmaxsol*(imodel-1)+isol,i)=wpesi(isol,i)
 100  continue

      if(ipesi.eq.0) then
         wsum=dzero
         do 80 isize=1,nsize
            wsum=wsum+dpesi(isize)*
     &           (wresid(isol,isize)**ddue)
 80      continue

         cp(nmaxsol*(imodel-1)+isol,npre+1)=
     &        wsum/dvar + (dtotpesi*dsize)*
     &        dlog(ddue*dpi*dvar)+  
     &       (alpha*ncol)
C      write(*,*) cp(nmaxsol*(imodel-1)+isol,npre+1)
      else
         wsum=dzero
         do 85 isize=1,nsize
            wsum=wsum+wpesi(isol,isize)*
     &           (wresid(isol,isize)**ddue)
 85      continue

         cp(nmaxsol*(imodel-1)+isol,npre+1)=
     &   wsum/dvar + (wtotpesi(isol)*dsize)*
     &        dlog(ddue*dpi*dvar)+  
     &       (alpha*ncol)
C      write(*,*) cp(nmaxsol*(imodel-1)+isol,npre+1)
      endif
      varia(nmaxsol*(imodel-1)+isol)=wvaria(isol)
      nsame(nmaxsol*(imodel-1)+isol)=nwsame(isol)
      totpesi(nmaxsol*(imodel-1)+isol)=wtotpesi(isol)
 70   continue

      endif

 900  continue
      return
      end
























      SUBROUTINE WLEAICFAST (YDATA,XDATA,INTER,NSIZE,NVAR,
     & NBOOT,NGRP,NREP,IRAF,RK,
     & RPREC,REQUAL,IMAX,VAR,NMAXSOL,DMINPESI,ALPHA,
     & cp,wparam,wvaria,wresid,wtotpesi,wpesi,nwsame,info)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Weighted Akaike Information Criterion based on 
C     WLE in the normal regression linear model
C     FAST VERSION USING WLS INSTEAD OF WLE IN THE SUBMODELS
C     Author: Claudio Agostinelli 
C             DAIS
C             Universita' di Venezia
C             30121 Padova
C             ITALIA
C
C     E-mail: claudio@unive.it
C
C     July, 7 2011
C
C     Version: 0.1
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Copyright (C) 2011 Claudio Agostinelli
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
C     ALPHA     input    D      1            the coefficient of the number of the parameters.
C
C
C     cp        output   D      NREP*(NVAR+INTER+1) 
C     wparam     output   D      NREP*(NVAR+INTER) the WLS parameters
C     wvaria     output   D      NREP the WLS variance of the residuals
C     wresid     output   D      NREP*NSIZE the WLS residuals for each model
C     wtotpesi   output   D      NREP the total sum of the weights
C     wpesi      output   D      NSIZE the weights of the full model
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     nsol               I      NMAXSOL         the total number of solutions for teh full model
C     nconv              I      NMAXSOL         number of boostrap sampling that 
C                                            does not converge for the full model

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dqu=2.0d00)
      parameter(dpi=3.141592653d00)

      dimension ydata(nsize),xdata(nsize,nvar) 
      dimension xidata(nsize,nvar+inter)
      dimension wydata(nsize)
      dimension wxiidata(nsize,nvar+inter),wxidata(nsize,nvar+inter)
      dimension wparam(nmaxsol,nvar+inter)
      dimension wvaria(nmaxsol) 
      dimension wresid(nmaxsol,nsize) 
      dimension wtotpesi(nmaxsol)
      dimension wpesi(nmaxsol,nsize)
      dimension dden(nmaxsol,nsize), dmod(nmaxsol,nsize)
      dimension ddelta(nmaxsol,nsize) 
      dimension nwsame(nmaxsol)
      dimension nmodel(nvar+inter)
      dimension cp(nrep,nvar+inter+1)
      dimension dpesi(nsize)
      dimension ddd(nsize)
      dimension rwmat(nsize,nsize)
      dimension dparam(nvar+inter), rparam(nvar+inter)
      dimension dresid(nsize)
      dimension xparam(nsize)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C     wleregfix: subroutine for evaluating a wle estimators
C                for the normal regression model 
      external wleregfix
C     the blas: dgemm subroutine give matrix-matrix products
C     the blas: dgemv subroutine give matrix-vector products 
      external dgemm, dgemv
C
      external wls
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      dsize=nsize
      info=0
C     Check the presence of the intercept
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

C      nmaxmod=(2**npre)-1
      nmaxmod=nrep
CCCCCC FOR THE FULL MODEL
      call wleregfix(ydata,xidata,0,nsize,npre,
     & npre,nboot,ngrp,nmaxsol,iraf,rk,rprec,requal,imax,
     & wparam,wvaria,wresid,wtotpesi,wpesi,
     & dden,dmod,ddelta,nwsame,nsol,nconv)

      if (nconv.eq.nboot) then
        info=1
      endif

      if (nsol.gt.0) then 
        indice=0
        dvar=wvaria(1)+duno
CCCC         write(*,*) nsol
        do 50 i=1,nsol
          if (dvar.gt.wvaria(i).and.dminpesi
     &           .lt.wtotpesi(i)) then
            dvar=wvaria(i)
            indice=i
C                write(*,*) wtotpesi
          endif
 50     continue
        if (indice.eq.0) then
          info=2
          return
        endif
        dtotpesi=wtotpesi(indice)
        do 55 isize=1,nsize
          dpesi(isize)=wpesi(indice,isize)
 55     continue
        if (var.ne.dzero) then
          dvar=var
        endif

        do 140 j=1,nsize
          do 150 i=1,nsize
            if (i.eq.j) then
              rwmat(i,j)=dsqrt(dpesi(i))
            else
              rwmat(i,j)=dzero
            endif
 150      continue
 140    continue
C      write(*,*) 'sono arrivato qui 2'
C        write(*,*) rwmat

C     The weighted observations
        call dgemv('N',nsize,nsize,duno,rwmat,nsize,ydata,1,
     & dzero,wydata,1)

        call dgemm('N','N',nsize,npre,nsize,duno,rwmat,nsize,
     & xidata,nsize,dzero,wxidata,nsize)

C      write(*,*) 'sono arrivato qui 3'

C        write(*,*) wydata
C        write(*,*) wxidata

C From the biggest model to the smallest one.
        do 900 imodel=nmaxmod,1,-1
        call dectobin (imodel,npre,nmodel)
C         write(*,*) nmodel
         ncol=0
         do 300 i=1,npre
           if (nmodel(i).eq.1) then
             ncol=ncol+1
             do 400 j=1,nsize
               wxiidata(j,ncol)=wxidata(j,i)
 400         continue
           endif
 300     continue
         call wls(wydata,wxiidata,ddd,nsize,npre,ncol,0,dparam,iinfo)
         ipos=0
         do 60 i=1,npre
           cp(imodel,i)=nmodel(i)
C           write(*,*) cp(imodel,i)
           if (nmodel(i).eq.1) then
             ipos=ipos+1
             rparam(i)=dparam(ipos)
           else
             rparam(i)=dzero
           endif
 60      continue

CCCCCCCCCCCCCCCCC Fitted values 
         call dgemv('N',nsize,npre,duno,xidata,nsize,rparam,1,
     &     dzero,xparam,1)
CCCCCCCCCCCCCCCCC Residuals           
         do 100 i=1,nsize
           dresid(i) = (ydata(i)-xparam(i))**ddue
 100     continue
         wsum=dzero
         do 80 isize=1,nsize
           wsum=wsum+dpesi(isize)*
     &       (dresid(isize)**ddue)
 80      continue

         cp(imodel,npre+1)=
     &     wsum/dvar + (dtotpesi*dsize)*
     &     dlog(ddue*dpi*dvar)+  
     &     (alpha*ncol)
C      write(*,*) cp(imodel,npre+1)
 900  continue
      endif
      return
      end
























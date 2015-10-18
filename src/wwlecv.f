       SUBROUTINE WWLECV (YDATA,XDATA,INTER,NSIZE,NVAR,
     & NBOOT,NGRP,NREP,NMCCV,NSPLIT,IRAF,RK,
     & RPREC,REQUAL,IMAX,NMAXSOL,DMINPESI,
     & cv,wparam,wvaria,wresid,wtotpesi,wpesi,nwsame,indice,
     & info)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Weighted Cross-Validation based on 
C     WLE in the normal regression linear model
C     In questo caso usiamo stima di verosimiglianza pesata anche nei
C     sottomodelli per le stime dei residui. Pesi fissi al modello piu'
C     grande come in WLECV. Per confronto con WBIC
C     Author: Claudio Agostinelli 
C             DAIS
C             Universita' Ca' Foscari
C             30121 Venezia
C             ITALIA
C
C     E-mail: claudio@.unive.it
C
C     October, 13 2012
C
C     Version: 0.2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Copyright (C) 2012 Claudio Agostinelli
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
C     NMCCV     input    I      1            number of MonteCarlo replication for the Cross-Validation
C     NSPLIT    input    I      1            dimension of the split for the costruction-validation sample size
C
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
C     NMAXSOL   input    I      1            Maximum number of solutions to be considered for each model 
C     DMINPESI  input    D      1            Minimum percentage of the weights for which a roots is considered
C
C     cv        output   D      NREP*(NVAR+INTER+1) 
C     wparam    output   D      NMAXSOL*(NVAR+INTER) the WLE parameters
C     wvaria    output   D      NMAXSOL the WLE variance of the residuals
C     wresid    output   D      NMAXSOL*NSIZE the WLE residuals for each model
C     wtotpesi  output   D      NMAXSOL the total sum of the weights
C     wpesi     output   D      NMAXSOL*NSIZE the weights
C     nwsame    output   I      NMAXSOL frequencies of each root
C     indice    output   I      1       position of the root used for the weights
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

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dqu=2.0d00)

      parameter(rerr=1.0d-20)

      dimension ydata(nsize),xdata(nsize,nvar) 
      dimension xidata(nsize,nvar+inter)
      dimension rwmat(nsize,nsize)
      dimension wxidata(nsize,nvar+inter),wydata(nsize)
      dimension wparam(nmaxsol,nvar+inter)
      dimension wvaria(nmaxsol)
      dimension wresid(nmaxsol,nsize) 
      dimension wtotpesi(nmaxsol)
      dimension wpesi(nmaxsol,nsize)
      dimension dden(nmaxsol,nsize), dmod(nmaxsol,nsize)
      dimension ddelta(nmaxsol,nsize) 
      dimension nwsame(nmaxsol)
      dimension nmodel(nvar+inter)
      dimension dconv(nrep)
      dimension cv(nrep,nvar+inter+1)
      dimension dpesi(nsize)

      dimension wcresid(nsplit)
      dimension wctotpesi(nsplit)
      dimension wcpesi(nsplit)
      dimension dcden(nsplit)
      dimension dcmod(nsplit)
      dimension dcdelta(nsplit)

      dimension dparam(nvar+inter)
      dimension tparam(nvar+inter)

      dimension nstart(nsize)
      dimension nrand(nmccv,nsplit)
      dimension ycsub(nsplit) 
      dimension xcsub(nsplit,nvar+inter)
      dimension yvsub(nsize-nsplit) 
      dimension xvsub(nsize-nsplit,nvar+inter) 
      dimension wyvsub(nsize-nsplit) 
      dimension wxvsub(nsize-nsplit,nvar+inter) 
      dimension xparam(nsize-nsplit)
      dimension xcdata(nsplit,nvar+inter)

C     wleregfix: subroutine for evaluating a wle estimators
C                for the normal location model 
      external wleregfix
      external genprm
C     the blas: dgemm subroutine give matrix-matrix products
C     the blas: dgemv subroutine give matrix-vector products 
      external dgemm, dgemv
      external dectobin

      dsize=nsize
      info=0
      iinfo=0
      nvalid=nsize-nsplit
      dvalid=nvalid
      dmccv=nmccv

      do 5 i=1,nsize
         nstart(i)=i
 5    continue

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
C      nmaxmod=(2**npre)-1
C
      nmaxmod=nrep
 
C      write(*,*) nmaxmod
C      write(*,*) nrep

      do 65 imodel=nmaxmod,1,-1

C      write(*,*) imodel
 
         call dectobin (imodel,npre,nmodel)
      do 60 i=1,npre
         cv(imodel,i)=nmodel(i)
 60   continue
         cv(imodel,npre+1)=dzero
         dconv(imodel)=dzero
 65   continue

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
          do 50 i=1,nsol
             if(dvar.gt.wvaria(i).and.dminpesi
     &            .lt.wtotpesi(i)) then
                dvar=wvaria(i)
                indice=i
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

      do 150 i=1,nsize
         rwmat(i,i)=dsqrt(dpesi(i))
 150  continue

C
C     The weighted observations
C

      call dgemv('N',nsize,nsize,duno,rwmat,nsize,ydata,1,
     & dzero,wydata,1)

      call dgemm('N','N',nsize,npre,nsize,duno,rwmat,nsize,
     & xidata,nsize,dzero,wxidata,nsize)

C      write(*,*) 'sono arrivato qui 2'

C
C     Begin of the Monte Carlo replications  
C
      do 9999 imc=1,nmccv
C
C     Costruction-Validation Sample
C
7777  call genprm(nstart,nsize)

      nagain=0

      do 90 icheck=2,imc
      do 70 istart=1,nsplit
      do 80 iistart=1,nsplit
         if(nstart(istart).eq.nrand(icheck,iistart)) then
             nagain=nagain+1
         endif
 80   continue   
 70   continue
         if(nagain.eq.nsplit) then
             goto 7777
         else
            nagain=0
         endif
 90   continue

C      write(*,*) 'sono arrivato qui 3'

      do 110 i=1,nsplit
         ycsub(i)=ydata(nstart(i))
      do 120 j=1,npre   
         xcsub(i,j)=xidata(nstart(i),j)
 120  continue
         nrand(imc,i)=nstart(i)
 110  continue

C      write(*,*) 'sono arrivato qui 4'

      dtotval=dzero
      do 130 i=nsplit+1,nsize
        yvsub(i-nsplit)=ydata(nstart(i))
        wyvsub(i-nsplit)=wydata(nstart(i))        
        dtotval=dtotval+dpesi(nstart(i))
        do 140 j=1,npre
          xvsub(i-nsplit,j)=xidata(nstart(i),j)
          wxvsub(i-nsplit,j)=wxidata(nstart(i),j)
 140    continue
 130  continue

C      write(*,*) 'sono arrivato qui 5'

C
C From the biggest model to the smallest one.
C

      do 900 imodel=nmaxmod,1,-1

         call dectobin (imodel,npre,nmodel)
         
C         write(*,*) nmodel

         ncol=0
         do 30 i=1,npre
            if (nmodel(i).eq.1) then
            ncol=ncol+1
               do 40 j=1,nsplit
                  xcdata(j,ncol)=xcsub(j,i)
 40            continue
            endif
 30      continue

C      write(*,*) 'sono arrivato qui 6'

       call wleregfix(ycsub,xcdata,0,nsplit,npre,
     &  ncol,nboot,ngrp,1,iraf,rk,rprec,requal,imax,
     &  tparam,wvaria,wcresid,wctotpesi,wcpesi,
     &  dcden,dcmod,dcdelta,nwsame,nsol,iconv)

C      write(*,*) 'sono arrivato qui 7'

      if (iconv.eq.nboot) then
         dconv(imodel)=dconv(imodel)+duno
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

C      write(*,*) 'sono arrivato qui 9'

      call dgemv('N',nvalid,npre,duno,wxvsub,nvalid,dparam,1,
     & dzero,xparam,1)
           
      dfit=dzero
      do 100 i=1,nvalid
         if(abs(wyvsub(i)-xparam(i)).gt.rerr) then
            dfit=dfit+(wyvsub(i)-xparam(i))**ddue
         endif
 100  continue
         dfit=dfit/dtotval   

      cv(imodel,npre+1)=cv(imodel,npre+1)+dfit

      endif

 900  continue

C      write(*,*) 'sono arrivato qui 10'

C     End of the Monte Carlo replications 
 9999 continue

C      write(*,*) 'sono arrivato qui 11'

      do 210 i=1,nmaxmod
C         write(*,*) dconv(i)
         cv(i,npre+1)=cv(i,npre+1)/(dmccv-dconv(i))
 210  continue
      return
      end

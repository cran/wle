      SUBROUTINE WSTEP (YDATA,XDATA,INTER,NSIZE,NVAR,
     & NBOOT,NGRP,NREP,IRAF,RK,NTYPE,
     & RPREC,REQUAL,IMAX,NMAXSOL,DMINPESI,FIN,FOUT,NMETHOD,
     & dstep,dfparam,dfvar,dfresid,dftot,dfpesi,nfsame,
     & indice,info,imodel,nfsol)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Weighted Stepwise, Forward and Backward regression based on 
C     WLE in the normal regression linear model
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Venezia
C             30121 Venezia
C             ITALIA
C
C     E-mail: claudio@unive.it
C
C     April, 10, 2005
C
C     Version: 0.5
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Copyright (C) 2005 Claudio Agostinelli
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
C     NBOOT     input    I      1            number of bootstrap replication
C     NGRP      input    I      1            dimension of the subsample 
C     NREP      input    I      1            number of model to be reported
C     IRAF      input    I      1            type of RAF
C                                            1: Hellinger distance 
C                                            2: Negative Exponential disparity 
C                                            3: Chi squared disparity
C     RK        input    D      1            smoothing parameter
C     NTYPE     input    I      1            type of stepwise               C                                             1: Forward
C                                             2: Backward
C                                             3: Stepwise
C
C     RPREC     input    D      1            precision of the convergence 
C                                            in absolute value
C     REQUAL    input    D      1            when two roots are said equal,
C                                            must be less than RPREC 
C     IMAX      input    I      1            maximum number of iterations for each starting points                        
C     NMAXSOL   input    I      1            Maximum number of solutions to be considered for each model 
C     DMINPESI  input    D      1            Minimum percentage of the weights for which a roots is considered
C     FIN       input    D      1
C     FOUT      input    D      1
C     NMETHOD   input    I      1
C
C     dstep     output   D      NREP*(NVAR+INTER+1) 
C     dfparam   output   D      NMAXSOL*(NVAR+INTER) the WLE parameters
C     dfvaria   output   D      NMAXSOL the WLE variance of the residuals
C     dfresid   output   D      NMAXSOL*NSIZE the WLE residuals for each model
C     dftot     output   D      NMAXSOL the total sum of the weights
C     dfpesi    output   D      NMAXSOL*NSIZE the weights
C     nfsame    output   I      NMAXSOL frequencies of each root
C     indice    output   I      1       position of the root used for the weights
C     info      output   I      1
C     imodel    output   I      1       number of model scanned
C     nfsol     output   I      1       number of roots found
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

      dimension ydata(nsize),xdata(nsize,nvar) 
      dimension xidata(nsize,nvar+inter)
      dimension xisub(nsize,nvar+inter)

C      dimension wxidata(nsize,nvar+inter)

      dimension wparam(nmaxsol,nvar+inter)
      dimension wvaria(nmaxsol)
      dimension avaria(nmaxsol)
      dimension aavaria(nvar+inter,nmaxsol)
 
      dimension wresid(nmaxsol,nsize) 
      dimension wtotpesi(nmaxsol)
      dimension wpesi(nmaxsol,nsize)
      dimension dden(nmaxsol,nsize), dmod(nmaxsol,nsize)
      dimension ddelta(nmaxsol,nsize) 
      dimension nwsame(nmaxsol)
      dimension nmodel(nvar+inter)
      dimension dstep(nrep,nvar+inter+1)
      dimension dpesi(nsize)

      dimension dfparam(nmaxsol,nvar+inter)
      dimension dfvar(nmaxsol)
      dimension dfpesi(nmaxsol,nsize)
      dimension dfden(nmaxsol,nsize), dfmod(nmaxsol,nsize)
      dimension dfdelta(nmaxsol,nsize) 

      dimension dfresid(nmaxsol,nsize)
      dimension dftot(nmaxsol)

C      dimension dfsame(nmaxsol)
C      dimension nstart(nsize)
C      dimension wydata(nsize) 
C      dimension wxdata(nsize,nvar+inter) 

      dimension dparam(nvar+inter)
      dimension xparam(nsize)

      dimension nvarin(nvar+inter)
      dimension ntemp(nvar+inter)

      dimension ftest(nvar+inter,nmaxsol)
      dimension ftunico(nvar+inter)
      
      dimension indsol(nvar+inter)


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C     wleregfix: subroutine for evaluating a wle estimators
C                for the normal location model 
      external wleregfix
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

      dsize=nsize
      info=0
      iinfo=0
      nstep=0
      istep=1

      if(ntype.eq.3) then
         nstep=1
         ntype=1
      endif

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

      imodel=0
      iconv=0

      call wleregfix(ydata,xidata,0,nsize,npre,
     &  npre,nboot,ngrp,nmaxsol,iraf,rk,rprec,requal,imax,
     &  dfparam,dfvar,dfresid,dftot,dfpesi,
     &  dfden,dfmod,dfdelta,nfsame,nfsol,iconv)

C       write(*,*) 'nsize: ',nsize
C       write(*,*) 'npre: ',npre
C       write(*,*) 'ngrp: ',ngrp
C       write(*,*) 'nmaxsol: ',nmaxsol
C       write(*,*) 'iraf: ',iraf
C       write(*,*) 'rprec: ',rprec
C       write(*,*) 'requal: ',requal
C       write(*,*) 'imax: ',imax
C       write(*,*) 'dfpesi: ',dfpesi
C       write(*,*) 'dfparam',dfparam


C       write(*,*) 'iconv: ', iconv , ' nboot: ', nboot


       if(iconv.eq.nboot) then
          info=1
          return
       endif

          indice=0
          dvar=dfvar(1)+duno
          do 50 i=1,nfsol
             if(dvar.gt.dfvar(i).and.dminpesi
     &            .lt.dftot(i)) then
                dvar=dfvar(i)
                indice=i
             endif
 50       continue
          if(indice.eq.0) then
             info=2
             return
          endif
          dtotpesi=dftot(indice)
          do 55 isize=1,nsize
             dpesi(isize)=dfpesi(indice,isize)
 55       continue

C          write(*,*) 'dpesi: ',dpesi

          if(ntype.eq.1) then

             ntotin=0

             do 1050 i=1,nvar
                nmodel(i)=0
 1050        continue
             nmodel(npre)=1

             dmediay=0.0d00
             do 3200 i=1,nsize
                dmediay=dmediay+dpesi(i)*ydata(i)
 3200        continue
             dmediay=dmediay/dtotpesi

             dvarmod=0.0d00
             do 3100 i=1,nsize
                dvarmod=dvarmod+dpesi(i)*(ydata(i)-dmediay)**ddue
 3100        continue

          elseif(ntype.eq.2) then

             ntotin=npre

             do 1055 i=1,npre
                nmodel(i)=1
 1055        continue

             dvarmod=dvar*((dtotpesi*dsize)-npre)
                
          endif   

 9999 continue

C      write(*,*) 'nmodel: ', nmodel
C      write(*,*) 'ntype: ', ntype
      

      ipos=1
      ntotin=0

      if(ntype.eq.1) then

         do 1070 i=1,npre
            ntotin=ntotin+nmodel(i)
            nvarin(i)=0
            if(nmodel(i).eq.0) then
               nvarin(ipos)=i
               ipos=ipos+1
            endif   
 1070    continue

         ntotin=npre-ntotin

C         write(*,*) 'ntotin: ',ntotin

      elseif(ntype.eq.2) then   

         do 1071 i=1,npre
            ntotin=ntotin+nmodel(i)
            nvarin(i)=0
            if(nmodel(i).eq.1) then
               nvarin(ipos)=i
               ipos=ipos+1
            endif   
 1071    continue

C         write(*,*) 'ntotin: ',ntotin

      endif   

      do 1060 ipos=1,ntotin
         do 1065 i=1,npre
            ntemp(i)=nmodel(i)
 1065    continue


        if(ntype.eq.1) then
           ntemp(nvarin(ipos))=1
        elseif(ntype.eq.2) then   
           ntemp(nvarin(ipos))=0
        endif   

         jpos=1

      do 1025 j=1,npre
         if(ntemp(j).eq.1) then
            do 1020 i=1,nsize
               xisub(i,jpos)=xidata(i,j)       
 1020       continue
            jpos=jpos+1
         endif   
 1025 continue

C      write(*,*) 'nvarin: ', nvarin      
C      write(*,*) 'ntemp: ', ntemp

      jpos=jpos-1

      if(nmethod.eq.1) then

C         write(*,*) 'sono qui'

         call wls (ydata,xisub,dpesi,nsize,npre,jpos,1,
     &        dparam,iiinfo)
      
C         write(*,*) dparam

         call dgemv('N',nsize,npre,duno,xisub,nsize,dparam,1,
     &        dzero,xparam,1)


C         write(*,*) xparam

         do 116 ij=1,nsize
            wresid(1,ij)=ydata(ij)-xparam(ij)
 116     continue
         nsol=1

      else

         call wleregfix(ydata,xisub,0,nsize,npre,
     &        jpos,nboot,ngrp,nmaxsol,iraf,rk,rprec,requal,imax,
     &        wparam,wvaria,wresid,wtotpesi,wpesi,
     &        dden,dmod,ddelta,nwsame,nsol,iconv)
      
      endif
 
         do 1080 isol=1,nsol
            avaria(isol)=dzero
            do 1090 i=1,nsize
               avaria(isol)=avaria(isol)
     &              +dpesi(i)*wresid(isol,i)**ddue
 1090       continue
 
C         write(*,*) 'avaria ',avaria
C         write(*,*) 'dvar ',dvar
C         write(*,*) 'dvarmod ',dvarmod

         djpos=jpos
         if (ntype.eq.1) then
            ftest(ipos,isol) = (dvarmod - avaria(isol))
     &           / (avaria(isol)/((dtotpesi*dsize)-djpos))
         elseif(ntype.eq.2) then
            ftest(ipos,isol) = (avaria(isol) - dvarmod)
     &           / (dvarmod/((dtotpesi*dsize)-djpos-duno))            
         endif   


            aavaria(ipos,isol)=avaria(isol)
 1080    continue

         ftunico(ipos)=ftest(ipos,1)-duno

         do 2000 isol=1,nsol
            if(ftunico(ipos).lt.ftest(ipos,isol)) then
               ftunico(ipos)=ftest(ipos,isol)
               indsol(ipos)=isol
            endif   
 2000    continue

 1060  continue

C       write(*,*) 'ftunico: ', ftunico

       if(ntype.eq.1) then

          ftmax=ftunico(1)-duno
          indftmax=0
          do 2010 iselect=1,ntotin    
             if(ftunico(iselect).gt.ftmax) then
                ftmax=ftunico(iselect)
                indftmax=iselect
             endif
 2010     continue

          if(ftmax.gt.fin) then
             nmodel(nvarin(indftmax))=1
             dvarmod=aavaria(indftmax,indsol(indftmax))

             imodel=imodel+1

C             write(*,*) 'imodel ',imodel

             do 60 i=1,npre
                dstep(imodel,i)=nmodel(i)
 60          continue
             dstep(imodel,npre+1)=ftmax

             ntot=0
             do 3000 i=1,npre
                ntot=ntot+nmodel(i)
 3000        continue

             if(nstep.eq.1) then
                if(istep.eq.1) then
                   ntype=2
                   istep=0
                else
                   ntype=1
                   istep=1
                endif
             endif  

             if(ntot.eq.npre) then
                goto 8888
             endif   

             goto 9999
          endif   

       elseif(ntype.eq.2) then

          if(nstep.eq.1) then
             if(istep.eq.1) then
                ntype=2
                istep=0
             else
                ntype=1
                istep=1
             endif
          endif

          ftmax=ftunico(1)+duno
          indftmax=0
          do 2011 iselect=1,ntotin    
             if(ftunico(iselect).lt.ftmax) then
                ftmax=ftunico(iselect)
                indftmax=iselect
             endif
 2011     continue

C          write(*,*) 'ftmax ',ftmax

          if(ftmax.lt.fout) then
             nmodel(nvarin(indftmax))=0
             dvarmod=aavaria(indftmax,indsol(indftmax))

             imodel=imodel+1

C             write(*,*) 'imodel ',imodel

             do 61 i=1,npre
                dstep(imodel,i)=nmodel(i)
 61          continue
             dstep(imodel,npre+1)=ftmax

             ntot=0
             do 3001 i=1,npre
                ntot=ntot+nmodel(i)
 3001        continue

C             write(*,*) 'ntot ',ntot

             if(ntot.eq.0) then
                goto 8888
             endif   

             goto 9999
          endif

          if(nstep.eq.1) then
             goto 9999
          endif   
       endif

 8888  continue
       return
       end

















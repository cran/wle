      SUBROUTINE STEP (YDATA,XDATA,INTER,NSIZE,NVAR,
     & NREP,NTYPE,
     & FIN,FOUT,
     & dstep,info,imodel)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Stepwise, Forward and Backward regression
C      in the normal regression linear model
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Padova
C             35121 Padova
C             ITALIA
C
C     E-mail: claudio@stat.unipd.it
C
C     September, 6, 2001
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
C     NVAR      input    I      1            number of explanatory variables
C     NREP      input    I      1            number of model to be reported
C     NTYPE     input    I      1            type of stepwise               C                                             1: Forward
C                                             2: Backward
C                                             3: Stepwise
C
C     FIN       input    D      1
C     FOUT      input    D      1
C
C     dstep     output   D      NREP*(NVAR+INTER+1) 
C     info      output   I      1
C     imodel    output   I      1       number of model scanned
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

      dimension aavaria(nvar+inter)
 
      dimension wresid(nsize) 
      dimension nmodel(nvar+inter)
      dimension dstep(nrep,nvar+inter+1)
      dimension dpesi(nsize)

      dimension dparam(nvar+inter)
      dimension xparam(nsize)

      dimension nvarin(nvar+inter)
      dimension ntemp(nvar+inter)

      dimension ftest(nvar+inter)
      dimension ftunico(nvar+inter)
      
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

      dsize=nsize
      info=0
      nstep=0
      istep=1

      do 1 i=1, nsize
         dpesi(i)=duno
 1    continue

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

          if(ntype.eq.1) then

             ntotin=0

             do 1050 i=1,nvar
                nmodel(i)=0
 1050        continue
             nmodel(npre)=1

            mediay=0.0d00
             do 3200 i=1,nsize
                mediay=mediay+ydata(i)
 3200        continue
             mediay=mediay/dsize

             dvarmod=0.0d00
             do 3100 i=1,nsize
                dvarmod=dvarmod+(ydata(i)-mediay)**ddue
 3100        continue

          elseif(ntype.eq.2) then

             ntotin=npre

             do 1055 i=1,npre
                nmodel(i)=1
 1055        continue


        call wls (ydata,xidata,dpesi,nsize,npre,npre,1,
     &        dparam,info)
      
C         write(*,*) dparam

         call dgemv('N',nsize,npre,duno,xidata,nsize,dparam,1,
     &        dzero,xparam,1)

C         write(*,*) xparam

         do 111 ij=1,nsize
            wresid(ij)=ydata(ij)-xparam(ij)
 111        continue

            dvarmod=dzero
            do 112 i=1,nsize
               dvarmod=dvarmod
     &              +wresid(i)**ddue
 112        continue
                
          endif  

9999  continue
C      write(*,*) ''
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

C         write(*,*) 'sono qui'

         call wls (ydata,xisub,dpesi,nsize,npre,jpos,1,
     &        dparam,info)
      
C         write(*,*) dparam

         call dgemv('N',nsize,npre,duno,xisub,nsize,dparam,1,
     &        dzero,xparam,1)


C         write(*,*) xparam

         do 116 ij=1,nsize
            wresid(ij)=ydata(ij)-xparam(ij)
 116     continue


            avaria=dzero
            do 1090 i=1,nsize
               avaria=avaria
     &              +wresid(i)**ddue
 1090       continue

C         write(*,*) 'avaria ',avaria
C         write(*,*) 'dvarmod ',dvarmod

         if (ntype.eq.1) then
            ftest(ipos) = (dvarmod - avaria)
     &           / (avaria/(dsize-jpos-1))
         elseif(ntype.eq.2) then
            ftest(ipos) = (avaria - dvarmod)
     &           / (dvarmod/(dsize-jpos))            
         endif   


            aavaria(ipos)=avaria
 1080    continue


               ftunico(ipos)=ftest(ipos)

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
             dvarmod=aavaria(indftmax)

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
             dvarmod=aavaria(indftmax)

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

















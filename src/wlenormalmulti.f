      SUBROUTINE WLENORMMULTI (DATA,NSIZE,NVAR,NBOOT,
     & NGRP,NREP,IRAF,RK,
     & RPREC,REQUAL,IMAX,dmedia,varia,totpesi,pesi,nsame,
     & nsol,nconv)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Bootstrap roots search for WLE in the 
C     normal location and scale problem
C     in the multivarite context
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
C     DATA      input    D      NSIZE*NVAR   vector of the data
C     NSIZE     input    I      1            length of the data 
C     NVAR      input    I      1            number of variables
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
C     IMAX      input    I      1            maximum number of iteration for each starting points.
C
C     dmedia     output   D      NREP*NVAR    the WLE mean
C     varia     output   D      NREP*NVAR*NVAR the WLE variance and covariance matrix
C     totpesi   output   D      NREP         the total sum of the weights
C     pesi      output   D      NREP*NSIZE   the weights
C     nsame     output   I      1            frequencies of each root
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
      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dqu=2.0d00)
      parameter(dzero=0.0d00)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     rerr: The smallest double precision number can be treated
C     on the machine as denominator in a division 
C     this value should work in any machines
C
      parameter(rerr=1.0d-65)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dimension d(nsize),rm(nsize),delta(nsize),adelta(nsize)
      dimension ds(nsize),rw(nsize)
      dimension data(nsize,nvar), datat(nvar,nsize)
      dimension dmedia(nrep,nvar),varia(nrep,nvar,nvar) 
      dimension totpesi(nrep)
      dimension pesi(nrep,nsize),nsame(nboot)
      dimension rloc(nvar),dvsuno(nsize),dvguno(ngrp)
      dimension rsca(nvar,nvar)
      dimension rlocold(nvar), dmw(nsize,nsize), rnowsold(nvar,nvar)
      dimension rnows(nvar,nvar), rnowh(nvar,nvar)
      dimension rnowhs(nvar,nvar),dmvuno(nvar,nvar)

      dimension sub(ngrp,nvar), nstart(nsize)
      dimension storep(nboot,nvar)
      dimension storev(nboot,nvar,nvar)
      dimension storew(nboot,nsize)
      dimension storet(nboot)
      dimension nrand(nboot,ngrp)

      dimension ipvt(nvar)
      dimension ztemp(nvar), dwork(nvar)
      dimension datai(nvar), dataik(nvar)

      dimension ddeth(2), ddeths(2)
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C     the ranlib: genprm subroutine generate a permuation of an array 
      external genprm
C     the blas: dgemm subroutine give matrix-matrix products
C     the blas: dgemv subroutine give matrix-vector products 
      external dgemm, dgemv
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dsize=nsize
      dgrp=ngrp


C      write(*,*) ngrp
C      write(*,*) nvar

      do 30 i=1,nboot
         nsame(i)=0
      do 40 j=1,nvar
         storep(i,j)=dzero
      do 45 jj=1,nvar
         storev(i,j,jj)=dzero
 45   continue
 40   continue
 30   continue
      
      do 50 i=1,nsize
         nstart(i)=i
         dvsuno(i)=duno
 50   continue

      do 55 j=1,ngrp
         dvguno(j)=duno
 55   continue

      do 66 j=1,nvar
         do 67 jj=1,nvar
            if(j.eq.jj) then
               dmvuno(j,jj)=duno
            else
               dmvuno(j,jj)=dzero
            endif   
 67      continue
 66      continue

C      nconv=0 
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
         do 95 j=1,nvar
            sub(i,j)=data(nstart(i),j)
 95      continue
         nrand(iboot,i)=nstart(i)
 90   continue

      call dgemv('T',ngrp,nvar,duno,sub,ngrp,dvguno,1,
     & dzero,rloc,1)

      do 110 j=1,nvar
         rloc(j)=rloc(j)/dgrp
 110  continue

C      write(*,*) 'rloc :', rloc

      call dgemm('T','N',nvar,nvar,ngrp,duno,sub,ngrp,
     & sub,ngrp,dzero,rsca,nvar) 

C      write(*,*) 'rsca :', rsca

      do 115 j=1,nvar
         do 117 jj=1,nvar
            rnows(j,jj)=rsca(j,jj)/dgrp-rloc(j)*rloc(jj)
 117     continue
 115  continue

C      write(*,*) 'rnows :', rnows

      iter=0

C     Iteration Steps until convergence achieved

 9999 continue

      iter=iter+1

C      write(*,*) 'iter ', iter

      do 121 i=1,nvar
         rlocold(i)=rloc(i)
         do 122 j=1,nvar
            rnowsold(i,j)=rnows(i,j)
 122     continue
 121  continue

      call dgemm('N','N',nvar,nvar,nvar,rk,rnows,nvar,
     & dmvuno,nvar,dzero,rnowh,nvar) 

      call dgemm('N','N',nvar,nvar,nvar,(rk+duno),rnows,nvar,
     & dmvuno,nvar,dzero,rnowhs,nvar) 

C      write(*,*) 'rk :', rk      
C      write(*,*) 'rnowh :', rnowh
C      write(*,*) 'rnowhs :', rnowhs
  
      call dgeco(rnowh,nvar,nvar,ipvt,rcond,ztemp)
      call dgedi(rnowh,nvar,nvar,ipvt,ddeth,dwork,11)	

      deth=ddeth(1)*(10.0d00**ddeth(2))

C      write(*,*) 'Inversa rnowh: ', rnowh
C      write(*,*) 'Determinante rnowh: ', deth
C      write(*,*) 'ddeth: ', ddeth

 
      call dgeco(rnowhs,nvar,nvar,ipvt,rcond,ztemp)
      call dgedi(rnowhs,nvar,nvar,ipvt,ddeths,dwork,11)

      deths=ddeths(1)*(10.0d00**ddeths(2))

C      write(*,*) 'Inversa rnowhs: ', rnowhs
C      write(*,*) 'Determinante rnowhs: ', deths
C      write(*,*) 'ddeths: ', ddeths

      do 130 i=1,nsize
         d(i)=dzero
         rm(i)=dzero

         do 120 ik=1,nsize
 
            do 125 k=1,nvar
               dataik(k)=data(i,k)-data(ik,k)
 125        continue

            call dgemv('N',nvar,nvar,duno,rnowh,nvar,dataik,1,
     &           dzero,ztemp,1)

            dtemp=ddot(nvar,dataik,1,ztemp,1)

C            write(*,*) 'ztemp: ', ztemp
C            write(*,*) 'dtemp: ', dtemp
C            write(*,*) 'deth : ', deth
            
             ds(ik)=dexp(- dtemp / ddue)/dsqrt(deth)

C             write(*,*) 'ds  ', ds(ik)  

             d(i)=d(i)+ds(ik)

 120      continue

          d(i)=d(i)/nsize

C          write(*,*) 'd ', d(i)

          do 126 k=1,nvar
             datai(k)=data(i,k)-rloc(k)
 126      continue

          call dgemv('N',nvar,nvar,duno,rnowhs,nvar,datai,1,
     &           dzero,ztemp,1)

          dtemp=ddot(nvar,datai,1,ztemp,1)

   	  rm(i)=dexp(- dtemp / ddue)/dsqrt(deths)

C          write(*,*) 'rm ', rm(i)

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
             rw(i)=1-((delta(i)**ddue) / ((delta(i)**ddue) +
     &       ddue))
          endif

             if(rw(i).lt.0.0) then 
                rw(i)=dzero
             endif
             if(rw(i).gt.1.0) then
                rw(i)=duno
             endif
             
         else
             rw(i)=dzero
         endif         
 140   continue

      tot=dzero

      do 150 i=1,nsize
         tot=tot+rw(i)
 150   continue   

C       write(*,*) 'rw ', rw
C       write(*,*) 'tot ', tot      

      do 166 j=1,nsize
         do 167 jj=1,nsize
            if(j.eq.jj) then
               dmw(j,jj)=rw(j)
            else
               dmw(j,jj)=dzero
            endif   
167      continue
166      continue

C         write(*,*) 'dmw ', dmw

      call dgemv('T',nsize,nvar,duno,data,nsize,rw,1,
     & dzero,rloc,1)

      do 155 j=1,nvar
         rloc(j)=rloc(j)/tot
 155  continue

C      write(*,*) 'rloc :', rloc

      call dgemm('T','N',nvar,nsize,nsize,duno,data,nsize,
     & dmw,nsize,dzero,datat,nvar) 

C      write(*,*) 'datat ',datat

      call dgemm('N','N',nvar,nvar,nsize,duno,datat,nvar,
     & data,nsize,dzero,rsca,nvar) 

C      write(*,*) 'rsca :', rsca

      do 175 j=1,nvar
         do 177 jj=1,nvar
            rnows(j,jj)=rsca(j,jj)/tot-rloc(j)*rloc(jj)
 177     continue
 175  continue

C      write(*,*) 'rnows :', rnows

      diffm=dzero
      diffs=dzero

      do 181 i=1,nvar
         diffm=max(diffm,abs(rlocold(i)-rloc(i)))
         do 182 j=1,nvar
            diffs=max(diffs,abs(rnowsold(i,j)-rnows(i,j)))
 182     continue
 181  continue

      if(iter.gt.imax) then 
         nconv=nconv+1
         goto 8888 
      endif  

      if(diffm.gt.rprec.or.diffs.gt.rprec) then
         goto 9999
      endif   

C Convergence achieved

C      write(*,*) 'Converge ', iboot

C
C Is this a new root, then store it
C
      if(nsol.eq.0) then
         nsol=nsol+1
         nsame(1)=1 
         do 191 i=1,nvar
            storep(nsol,i)=rloc(i)
            do 192 j=1,nvar
               storev(nsol,i,j)=rnows(i,j)
 192        continue
 191     continue

         storet(nsol)=tot

         do 193 i=1,nsize
            storew(nsol,i)=rw(i)
 193        continue

      else

         do 300 isol=1,nsol

            diffm=dzero
            diffs=dzero

            do 310 i=1,nvar
               diffm=max(diffm,abs(rloc(i)-storep(isol,i)))
               do 320 j=1,nvar
                  diffs=max(diffs,abs(rnows(i,j)-storev(isol,i,j)))
 320           continue
 310        continue
          
            if(diffm.lt.requal.and.diffs.lt.requal
     &          ) then
                   nsame(isol)=nsame(isol)+1
                   goto 290
            endif
 300     continue
                   nsol=nsol+1
		   nsame(nsol)=1
                   do 330 i=1,nvar
                      storep(nsol,i)=rloc(i)
                      do 340 j=1,nvar
                         storev(nsol,i,j)=rnows(i,j)
 340                  continue
 330               continue

                   storet(nsol)=tot

                   do 200 i=1,nsize
                      storew(nsol,i)=rw(i)
 200               continue
 290               continue
      endif 

      if(nsol.eq.nrep) goto 6666

8888  continue 

C       
C end bootstrapping
C

 900  continue

C write down the results and return

 6666 continue

      do 410 isol=1,nsol 
         do 411 i=1,nvar
            dmedia(isol,i)=storep(isol,i)
            do 412 j=1,nvar
               varia(isol,i,j)=storev(isol,i,j)
 412        continue
 411     continue

         totpesi(isol)=storet(isol)
  410  continue

      do 420 isol=1,nsol 
         do 430 i=1,nsize
            pesi(isol,i)=storew(isol,i)
 430     continue 
 420  continue

      return
      end





































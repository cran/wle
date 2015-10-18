      SUBROUTINE WLENORM (DDATI,NSIZE,NLENGTH,NBOOT,NGRP,NREP,
     & IRAF,RK,
     & RPREC,REQUAL,IMAX,dmedia,varia,totpesi,pesi,
     & dden,dmod,ddelta,nsame,nsol,nconv)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Bootstrap roots search for WLE in the 
C     normal location and scale problem
C     
C     Author: Claudio Agostinelli 
C             DAIS
C             Ca' Foscari University
C             30121 Venice
C             ITALIA
C
C     E-mail: claudio@unive.it
C
C     October, 13, 2012
C
C     Version: 0.6
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
C     DDATI     input    D      NLENGTH      vector of the data
C     NSIZE     input    I      1            length of the data 
C     NLENGTH   input    I      1            
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
C     dmedia     output   D      NREP         the WLE mean
C     varia     output   D      NREP         the WLE variance
C     totpesi   output   D      NREP         the total sum of the weights
C     pesi      output   D      NREP*NSIZE   the weights
C     dden      output   D      NREP*NSIZE   the kernel density
C     dmod      output   D      NREP*NSIZE   the smoothed model
C     delta     output   D      NREP*NSIZE   the Pearson residuals
C     nsame     output   I      NREP         frequencies of each root
C     nsol      output   I      1            the total number of solutions
C     nconv     output   I      1            number of boostrap sampling that 
C                                            does not converge
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     rerr: The smallest double precision number can be treated
C     on the machine as denominator in a division 
C     this value should work in any machines
C
      parameter(rerr=1.0d-65)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dimension d(nsize),rm(nsize),delta(nsize),adelta(nsize)
      dimension ds(nsize),rw(nsize)
      dimension ddati(nlength)
      dimension dati(nsize)
      dimension dmedia(nrep),varia(nrep), totpesi(nrep)
      dimension pesi(nrep,nsize),nsame(nrep)
      dimension dden(nrep,nsize),dmod(nrep,nsize) 
      dimension ddelta(nrep,nsize)

      dimension sub(ngrp), nstart(nsize)
      dimension nrand(nboot,ngrp)

      external genprm

      do 1 i=1,nsize
         dati(i)=ddati(i)
 1    continue

      do 30 i=1,nrep
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
         sub(i)=dati(nstart(i))
         nrand(iboot,i)=nstart(i)
 90   continue

      totloc=0.0d00
      totvar=0.0d00

      do 100 i=1,ngrp

         totloc=totloc+sub(i)

 100  continue   

      rloc=totloc/ngrp

      do 110 i=1,ngrp

         totvar=totvar+((sub(i)-rloc)**2)

 110  continue

      rsca=totvar/ngrp

      rnowm=rloc
      rnows=rsca
      rnowh=rk*rnows
      rnowhs=(duno+rk)*rnows
      iter=0
      
C     Iteration Steps until convergence achieved

 9999 continue

      iter=iter+1

      do 130 i=1,nsize
         d(i)=0.0d00
         rm(i)=0.0d00

         do 120 ik=1,nsize
      ds(ik)=dexp(-((dati(i)-dati(ik))**ddue)/(ddue*rnowh))/
     &       dsqrt(rnowh)
            d(i)=d(i)+ds(ik)
 120      continue

         d(i)=d(i)/nsize      
         rm(i)=dexp(-((dati(i)-rnowm)**ddue)/(ddue*rnowhs))/
     &        dsqrt(rnowhs)

 130  continue

      do 140 i=1,nsize
         if(rm(i).gt.rerr) then
             delta(i)=d(i)/rm(i)-duno


C Type of RAF
C IRAF = 1 Hellinger distance
C IRAF = 2 Negative Exponential disparity
C IRAF = 3 Chi Squared disparity

          if(iraf.eq.1) then   
             adelta(i)=ddue*(dsqrt(delta(i)+duno)-duno)
          endif
          if(iraf.eq.2) then
             adelta(i)=ddue - (ddue+delta(i))*dexp(-delta(i))
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
 140   continue

      tot=dzero
      totloc=dzero
      totvar=dzero

      do 150 i=1,nsize
         tot=tot+rw(i)
         totloc=totloc+rw(i)*dati(i)
 150   continue   

      roldm=rnowm
      rnowm=totloc/tot

      do 160 i=1,nsize
         totvar=totvar+rw(i)*((dati(i)-rnowm)**2)
 160   continue

      rolds=rnows
      rnows=totvar/tot
      rnowh=rk * rnows
      rnowhs=(duno+rk)*rnows

      diffm=abs(rnowm-roldm)
      diffs=abs(rnows-rolds)

      if(iter.gt.imax) then 
         nconv=nconv+1
         goto 8888 
      endif  

      if(diffm.gt.rprec.or.diffs.gt.rprec) then
         goto 9999
      endif   

C Convergence achieved

C
C Is this a new root, then store it
C
      if(nsol.eq.0) then
         nsol=nsol+1
         nsame(1)=1 
         dmedia(nsol)=rnowm
         varia(nsol)=rnows
         totpesi(nsol)=tot
         do 170 i=1,nsize
            pesi(nsol,i)=rw(i)
            dden(nsol,i)=d(i)            
            dmod(nsol,i)=rm(i)
            ddelta(nsol,i)=delta(i)

 170     continue
      else
         do 180 isol=1,nsol
            diffm=abs(rnowm-dmedia(isol))
            diffs=abs(rnows-varia(isol))
            if(diffm.lt.requal.and.diffs.lt.requal
     &          ) then
                   nsame(isol)=nsame(isol)+1
                   goto 190
            endif
 180     continue
                   nsol=nsol+1
                   nsame(nsol)=1
      		   dmedia(nsol)=rnowm
      		   varia(nsol)=rnows
      		   totpesi(nsol)=tot
                   do 200 i=1,nsize
                      pesi(nsol,i)=rw(i)
                      dden(nsol,i)=d(i)            
                      dmod(nsol,i)=rm(i)
                      ddelta(nsol,i)=delta(i)
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


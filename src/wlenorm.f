      SUBROUTINE WLENORM (DATA,NSIZE,NBOOT,NGRP,NREP,IRAF,RK,
     & RPREC,REQUAL,IMAX,dmedia,varia,totpesi,pesi,nsame,nsol,
     & nconv)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Bootstrap roots search for WLE in the 
C     normal location and scale problem
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
C     DATA      input    D      NSIZE        vector of the data
C     NSIZE     input    I      1            length of the data 
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

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dqu=2.0d00)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     rerr: The smallest double precision number can be treated
C     on the machine as denominator in a division 
C     this value should work in any machines
C
      parameter(rerr=1.0d-65)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dimension d(nsize),rm(nsize),delta(nsize),adelta(nsize)
      dimension ds(nsize),rw(nsize)
      dimension data(nsize)
      dimension dmedia(nrep),varia(nrep), totpesi(nrep)
      dimension pesi(nrep,nsize),nsame(nboot)
  
      dimension sub(ngrp), nstart(nsize)
      dimension storep(nboot,3)
      dimension storew(nboot,nsize)
      dimension nrand(nboot,ngrp)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C     the genprm subroutine generate a permuation of an array 
      external genprm
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do 30 i=1,nboot
         nsame(i)=0
      do 40 j=1,3
         storep(i,j)=0.0d00

 40   continue
 30   continue
      
      do 50 i=1,nsize
         nstart(i)=i
 50   continue

C      nconv=0 
C      nsol=0


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
         sub(i)=data(nstart(i))
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
      ds(ik)=dexp(-((data(i)-data(ik))**dqu)/(dqu*rnowh))/dsqrt(rnowh)
            d(i)=d(i)+ds(ik)
 120      continue

         d(i)=d(i)/nsize      
	rm(i)=dexp(-((data(i)-rnowm)**dqu)/(dqu*rnowhs))/dsqrt(rnowhs)

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
                rw(i)=0.0d00
             endif
             if(rw(i).gt.1.0) then
                rw(i)=1.0d00
             endif
             
         else
             rw(i)=0.0d00
         endif         
 140   continue

      tot=0.0d00
      totloc=0.0d00
      totvar=0.0d00

      do 150 i=1,nsize
         tot=tot+rw(i)
         totloc=totloc+rw(i)*data(i)
 150   continue   

      roldm=rnowm
      rnowm=totloc/tot

      do 160 i=1,nsize
         totvar=totvar+rw(i)*((data(i)-rnowm)**2)
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
         storep(nsol,1)=rnowm
         storep(nsol,2)=rnows
         storep(nsol,3)=tot
         do 170 i=1,nsize
            storew(nsol,i)=rw(i)
 170     continue
      else
         do 180 isol=1,nsol
            diffm=abs(rnowm-storep(isol,1))
            diffs=abs(rnows-storep(isol,2))
            if(diffm.lt.requal.and.diffs.lt.requal
     &          ) then
                   nsame(isol)=nsame(isol)+1
                   goto 190
            endif
 180     continue
                   nsol=nsol+1
		   nsame(nsol)=1
      		   storep(nsol,1)=rnowm
      		   storep(nsol,2)=rnows
      		   storep(nsol,3)=tot
                   do 200 i=1,nsize
                      storew(nsol,i)=rw(i)
 200               continue
 190               continue
      endif 

      if(nsol.eq.nrep) goto 6666

8888  continue 
C       
C end bootstrapping
C

 900  continue

C write down the results and return

 6666 continue

      do 210 isol=1,nsol 
         dmedia(isol)=storep(isol,1)
         varia(isol)=storep(isol,2)
         totpesi(isol)=storep(isol,3)
  210  continue

      do 220 isol=1,nsol 
         do 230 i=1,nsize
            pesi(isol,i)=storew(isol,i)
 230     continue 
 220  continue

      return
      end






































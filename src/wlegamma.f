      SUBROUTINE WLEGAMMA(DATI,NSIZE,IRAF,RK,IUSE,DSUP,
     & RPREC,RPRECINT,RLAMBDA,ROMEGA, rw,d,rm)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Weights for the 
C     gamma model  
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Padova
C             35121 Padova
C             ITALIA
C
C     E-mail: claudio@stat.unipd.it
C
C     February, 19,   2001
C
C     Version: 0.1
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Copyright (C) 2001 Claudio Agostinelli
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
C     DATI      input    D      NSIZE        vector of the data
C     NSIZE     input    I      1            length of the data 
C     IRAF      input    I      1            type of RAF
C                                            1: Hellinger distance 
C                                            2: Negative Exponential disparity 
C                                            3: Chi squared disparity
C     RK        input    D      1            smoothing parameter
C     rw        output   D      NSIZE        the weights
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)


C     Area comune
      common/comune/ dx, dh, dl, doo, dgam 
C      save dx, dh, dl, doo, dgam 


      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dpi=3.141592654d00)
      parameter(nlimit=10000)
      parameter(nlenw=20000)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     rerr: The smallest double precision number can be treated
C     on the machine as denominator in a division 
C     this value should work in any machines
C
      parameter(rerr=1.0d-65)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dimension d(nsize),rm(nsize),delta(nsize),adelta(nsize)
      dimension ds(nsize),rw(nsize)
      dimension dati(nsize)

      dimension iwork(nlimit+1)
      dimension work(nlenw)
      dimension points(2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C
      external dmod
C
C      external dqagi
      external dqagp
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do 130 i=1,nsize
         d(i)=0.0d00
         rm(i)=0.0d00

         do 120 ik=1,nsize
            ds(ik)=dexp(-((dati(i)-dati(ik))**ddue)/(ddue*rk)) 
     &             + exp(-((dati(i)+dati(ik))**ddue)/(ddue*rk))
            d(i)=d(i)+ds(ik)
 120      continue

         d(i)=d(i)/nsize      

         dx=dati(i)
         dl=rlambda
         doo=romega
         dh=rk
         dgam=dgamma(romega)


         if (iuse.eq.1) then      
C         rprecint=rprec*1.00d-4

         points(1)=0.0d00
         points(2)=dsup

         call dqagp(dmod,0.0d00,dsup,2,points,0.0d00,
     &    rprecint,dtemp,abserr,neval,ier,nlimit,
     &    nlenw,nlast,iwork,work)

C         call dqagi(dmod,0.0d00,1,0.00d00,rprecint,dtemp,
C     &     abserr,neval,ier,nlimit,nlenw,nlast,iwork,work)

         else
              dtemp=(dl**doo) * (dx**(doo-duno)) * 
     &      (dexp(-dl*dx)) /dgam
              dtemp=dtemp*dsqrt(ddue*dpi*dh)
         endif


         rm(i)=dtemp

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

      return
      end

      double precision function dmod(t)
      implicit double precision (a-h,o-z)
      common/comune/ dx, dh, dl, doo, dgam
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)

      dmod=(dexp(-((dx-t)**ddue)/(ddue*dh))+
     &      dexp(-((dx+t)**ddue)/(ddue*dh)))
     &      * (dl**doo) * (t**(doo-duno)) * 
     &      (dexp(-dl*t)) /dgam

      return
      end








































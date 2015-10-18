      SUBROUTINE wlew2(DATI,NSIZE,ZDATA,NZSIZE,IRAF,RK,
     & DVAR,d,rm,totpesi,rw)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     
C     
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Padova
C             35121 Padova
C             ITALIA
C
C     E-mail: claudio@stat.unipd.it
C
C     September, 27, 2001
C
C     Version: 0.3
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PARAMETER:
C     NAME:     I/O:    TYPE:  DIMENSIONS:   DESCRIPTIONS:
C     DATI      input    D      NSIZE        vector of the data we want the weights
C     NSIZE     input    I      1            length of the vector dati 
C     ZDATA     input    D      NSIZE        vector of the data used for the empirical cumulative distribution function
C     NZSIZE    input    I      1            lenght of the vector zdata
C     IRAF      input    I      1            type of RAF
C                                            1: Hellinger distance 
C                                            2: Negative Exponential disparity 
C                                            3: Chi squared disparity
C     RK        input    D      1            smoothing parameter
C     DVAR      input    D      1         initial estimate for the variance
C     totpesi   output   D      1         the total sum of the weights
C     rw      output   D     NSIZE     the weights
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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

      dimension d(nsize), rm(nsize), delta(nsize), adelta(nsize)
      dimension rw(nsize), dati(nsize), zdata(nzsize)
 
      rnowh=rk * dvar
      rnowhs=(duno+rk)*dvar

      do 130 i=1,nsize
         d(i)=0.0d00
         rm(i)=0.0d00

         do 120 ik=1,nzsize
            d(i)=d(i)+dexp(-((dati(i)-zdata(ik))**ddue)/
     &       (ddue*rnowh))/dsqrt(rnowh)
 120      continue

         d(i)=d(i)/nzsize      
         rm(i)=dexp(-((dati(i))**ddue)/(ddue*rnowhs))/dsqrt(rnowhs)

 130  continue

      do 140 i=1,nsize
         if(rm(i).gt.rerr) then
             delta(i)=d(i)/rm(i)-duno

C Type of RAF
C IRAF = 1 Hellinger distance
C IRAF = 2 Negative Exponential disparity
C IRAF = 3 Symmetric Chi Squared disparity

          if(iraf.eq.1) then   
             adelta(i)=ddue*(dsqrt(delta(i)+duno)-duno)
          elseif(iraf.eq.2) then
             adelta(i)=ddue - (ddue+delta(i))*dexp(-delta(i))
          endif   
          if(iraf.ne.3) then
             rw(i)=(adelta(i)+duno)/(delta(i)+duno)
          else
             rw(i)=1-(delta(i) / (delta(i) +
     &       ddue))**ddue
          endif

             if(rw(i).lt.0.0) then 
                rw(i)=dzero
             elseif(rw(i).gt.1.0) then
                rw(i)=duno
             endif
             
         else
             rw(i)=dzero
         endif         

 140   continue

      totpesi=dzero

      do 150 i=1,nsize
         totpesi=totpesi+rw(i)
 150   continue   

      return
      end






      SUBROUTINE WLEMAH(XDAT,NXSIZE,NXVAR,DSUMW,IRAF,DTAU,
     & RK,rw,d,rm)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Weights for the 
C     Squared Mahalanobis Distance using exact distribution  
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Scienze Ambientali, Informatica e Statistica
C             Universita' Ca' Foscari
C             30121 VENEZIA
C             ITALIA
C
C     E-mail: claudio@unive.it
C
C     August, 21, 2013
C
C     Version: 0.1
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Copyright (C) 2013 Claudio Agostinelli
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
C     XDAT      input    D      NXSIZE       vector of the data 
C                                            (squared Mahalanobis distances)
C     NXSIZE    input    I      1            rows of X
C     NXVAR     input    I      1            columns of X 
C     IRAF      input    I      1            type of RAF
C                                            1: Hellinger distance 
C                                            2: Negative Exponential disparity 
C                                            3: Chi squared disparity
C                                            4: PDM
C                                            5: GKL
C     DTAU      input    D      1            parameter used in PDM and GKL
C     RK        input    D      1            smoothing parameter
C     rw        output   D      NXSIZE        the weights
C     d         output   D      NXSIZE        
C     rm        output   D      NXSIZE        
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      parameter(dzero=0.0d00)
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)
      parameter(dpi=3.141592654d00)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     rerr: The smallest double precision number can be treated
C     on the machine as denominator in a division 
C     this value should work in any machines
C
      parameter(rerr=1.0d-65)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dimension d(nxsize),rm(nxsize),delta(nxsize),adelta(nxsize)
      dimension ds(nxsize),rw(nxsize)
      dimension xdat(nxsize)

      external ddbeta

      dxsize=nxsize
      dxvar=nxvar
      dh=rk
      do 110 i=1,nxsize
        d(i)=dzero
        rm(i)=dzero
        do 100 ik=1,nxsize
          ds(ik)=dexp(-((xdat(i)-xdat(ik))**ddue)/(ddue*rk)) 
     &          + exp(-((xdat(i)+xdat(ik))**ddue)/(ddue*rk))
          d(i)=d(i)+ds(ik)
 100    continue
        d(i)=d(i)/(dxsize*dsqrt(ddue*dpi*rk))*(dsumw-1)**ddue/dsumw
        dx=xdat(i)*dsumw/(dsumw-1)**ddue
      call ddbeta(dx, dxvar/ddue, (dsumw - dxvar - duno)/ddue, 0, rm(i))
 
C        write(*,*) 'd', d(i)
C        write(*,*) 'rm', rm(i)
        if(rm(i).gt.rerr) then
          delta(i)=d(i)/rm(i)-duno
C        write(*,*) 'delta', delta(i)

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
          if(iraf.eq.4) then 
            adelta(i)=dtau*((delta(i) + duno)**(duno/dtau) - duno)
          endif
          if(iraf.eq.5) then 
            adelta(i)=dlog(dtau*delta(i)+duno)/dtau
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
 110  continue
      return
      end

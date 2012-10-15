      SUBROUTINE WLECHISQ(XDAT,YDAT,NXSIZE,NYSIZE,IRAF,DTAU,
     & RK,IUSE,DSUP,RPRECINT,ROMEGA,rw,d,rm)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Weights for the 
C     Chi Square model  
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Scienze Ambientali, Informatica e Statistica
C             Universita' Ca' Foscari
C             30121 VENEZIA
C             ITALIA
C
C     E-mail: claudio@unive.it
C
C     February, 21, 2011
C
C     Version: 0.2
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    Copyright (C) 2011 Claudio Agostinelli
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
C     YDAT        input    D      NYSIZE       vector where to evaluate weights 
C     NXSIZE    input    I      1            length of XDAT
C     NYSIZE    input    I      1            lenght of YDAT 
C     IRAF      input    I      1            type of RAF
C                                            1: Hellinger distance 
C                                            2: Negative Exponential disparity 
C                                            3: Chi squared disparity
C                                            4: PDM
C                                            5: GKL
C     DTAU      input    D      1            parameter used in PDM and GKL
C     RK        input    D      1            smoothing parameter
C     IUSE      input    I      1            use the smoothing model or not    
C     DSUP      input    D      1 
C     RPRECINT  input    D      1
C     ROMEGA    input    D      1            omega parameter of the gamma dist shape
C     rw        output   D      NXSIZE        the weights
C     d         output   D      NXSIZE        
C     rm        output   D      NXSIZE        
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

C     Area comune
      common/comchi/ dx, dh, doo
C      save dx, dh, doo

      parameter(dzero=0.0d00)
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
      dimension d(nysize),rm(nysize),delta(nysize),adelta(nysize)
      dimension ds(nysize),rw(nysize)
      dimension xdat(nxsize),ydat(nysize)

      dimension iwork(nlimit+1)
      dimension work(nlenw)
      dimension points(2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C
      external rmod
      external dqagp
      external dgammac
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dxsize=nxsize
      doo=romega
      dh=rk
      ddmax = dzero
      do 100 i=1,nxsize
        ddmax = max(xdat(i), ddmax)
 100  continue
      dsup = ddmax+dsup
      do 130 i=1,nysize
        d(i)=dzero
        rm(i)=dzero

        do 120 ik=1,nxsize
          ds(ik)=dexp(-((ydat(i)-xdat(ik))**ddue)/(ddue*rk)) 
     &          + exp(-((ydat(i)+xdat(ik))**ddue)/(ddue*rk))
          d(i)=d(i)+ds(ik)
 120    continue
        d(i)=d(i)/(dxsize*dsqrt(ddue*dpi*rk))

        dx=ydat(i)
        if (iuse.eq.1) then
          points(1)=dzero
          points(2)=12.0d00*rk*romega + dsup
          call dqagp(rmod,dzero,dsup,2,points,dzero,
     &      rprecint,dtemp,abserr,neval,ier,nlimit,
     &      nlenw,nlast,iwork,work)
          rm(i)=dtemp/dsqrt(ddue*dpi*rk)
        else
          rm(i)=dgammac(dx,doo,ddue,0)
        endif
 130  continue

      do 140 i=1,nysize
C        write(*,*) 'd', d(i)
C        write(*,*) 'rm', rm(i)
        if(rm(i).gt.rerr) then
          delta(i)=d(i)/rm(i)-duno
C          write(*,*) 'delta', delta(i)

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

            if(rw(i).lt.0.0) then 
              rw(i)=dzero
            endif
            if(rw(i).gt.1.0) then
              rw(i)=duno
            endif
             
          else
            rw(i)=dzero
            delta(i)=duno/rerr
          endif         
 140    continue
      return
      end

      double precision function rmod(t)
      implicit double precision (a-h,o-z)
      common/comchi/ dx, dh, doo
      parameter(dzero=0.0d00)     
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)

      rmod=dgammac(t,doo,ddue,0)
      rmod=rmod*(dexp(-((dx-t)**ddue)/(ddue*dh))+
     &      dexp(-((dx+t)**ddue)/(ddue*dh)))
      return
      end

      SUBROUTINE WLEGAMMA(DATI,DY,NXSIZE,NYSIZE,IRAF,DTAU,
     & RK,IUSE,DSUP,RPREC,RPRECINT,RLAMBDA,ROMEGA,rw,d,rm)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Weights for the 
C     gamma model  
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' Ca' Foscari
C             30121 VENEZIA
C             ITALIA
C
C     E-mail: claudio@unive.it
C
C     July, 29, 2011
C
C     Version: 0.6-1
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
C     DATI      input    D      NXSIZE       vector of the data
C     DY        input    D      NYSIZE       vector where to evaluate weights 
C     NXSIZE    input    I      1            length of DATI
C     NYSIZE    input    I      1            lenght of DY 
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
C     RPREC     input    D      1
C     RPRECINT  input    D      1
C     RLAMBDA   input    D      1            lambda parameter of the gamma dist
C     ROMEGA    input    D      1            omega parameter of the gamma dist
C     rw        output   D      NXSIZE        the weights
C     d         output   D      NXSIZE        
C     rm        output   D      NXSIZE        
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)


C     Area comune
      common/comune/ dx, dh, dl, doo
C Tolto dgam dalla common area
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
      dimension d(nysize),rm(nysize),delta(nysize),adelta(nysize)
      dimension ds(nysize),rw(nysize)
      dimension dati(nxsize),dy(nysize)

      dimension iwork(nlimit+1)
      dimension work(nlenw)
      dimension points(2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     EXTERNAL SUBROUTINE
C
      external dmodel
C
      external dqagp
      external dgammac
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      dxsize=nxsize
      do 130 i=1,nysize
         d(i)=0.0d00
         rm(i)=0.0d00

         do 120 ik=1,nxsize
            ds(ik)=dexp(-((dy(i)-dati(ik))**ddue)/(ddue*rk)) 
     &             + dexp(-((dy(i)+dati(ik))**ddue)/(ddue*rk))
            d(i)=d(i)+ds(ik)
 120      continue
         d(i)=d(i)/(dxsize*dsqrt(ddue*dpi*rk))

         dx=dy(i)
         dl=rlambda
         doo=romega
         dh=rk
         if (iuse.eq.1) then
C         rprecint=rprec*1.00d-4
           points(1)=dzero
           points(2)=dsup

           call dqagp(dmodel,dzero,dsup,2,points,dzero,
     &       rprecint,dtemp,abserr,neval,ier,nlimit,
     &       nlenw,nlast,iwork,work)

           dtemp=dtemp/dsqrt(ddue*dpi*rk)
         else
            dtemp=dgammac(dx,doo,duno/dl,0)
         endif
         rm(i)=dtemp

 130  continue

      do 140 i=1,nysize
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
                rw(i)=0.0d00
             endif
             if(rw(i).gt.1.0) then
                rw(i)=1.0d00
             endif
             
         else
             rw(i)=0.0d00
             delta(i)=1.0d00/rerr
         endif         
 140   continue

CCCCC ADDED 20151017
      rprec = 0.0D+00
      return
      end

      double precision function dmodel(t)
      implicit double precision (a-h,o-z)
      common/comune/ dx, dh, dl, doo
C Tolto dgam dalla common area
      parameter(duno=1.0d00)
      parameter(ddue=2.0d00)

      dmodel=dgammac(t,doo,duno/dl,0)
      dmodel=dmodel*(dexp(-((dx-t)**ddue)/(ddue*dh))+
     &      dexp(-((dx+t)**ddue)/(ddue*dh)))
      return
      end

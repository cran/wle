      SUBROUTINE WLENMIX (DATI,NSIZE,IP,RLOC,RSCA,RRPR, 
     & NBOOT,NREP,IRAF,RK,RPREC,REQUAL,IMAX,IALL,DMIN,
     & dmedia,varia,
     & rprop,totpesi,pesi,dden,dmod,ddelta,nsame,nsol,nconv)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Bootstrap roots search for WLE in the 
C     normal location and scale problem in mixture models
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Padova
C             35121 Padova
C             ITALIA
C
C     E-mail: claudio@stat.unipd.it
C
C     August, 12, 2001
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
C     DATI      input    D      NSIZE        vector of the data
C     NSIZE     input    I      1            length of the data 
C     NBOOT     input    I      1            number of bootstrap replication
C     IP        input    I      1            number of components
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
C     IALL      input    I      1            look for all the components
C     DMIN      input    D      1            minimum weights
C
C     dmedia    output   D      NREP         the WLE mean
C     varia     output   D      NREP         the WLE variance
C     totpesi   output   D      NREP         the total sum of the weights
C     pesi      output   D      NREP*NSIZE   the weights
C     dden      output   D      NREP*NSIZE   the kernel density estimator
C     dmod      output   D      NREP*NSIZE   the smoothed model
C     ddelta    output   D      NREP*NSIZE   the Pearson residuals
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
      dimension dati(nsize)
      dimension dmedia(nrep,ip),varia(nrep,ip), totpesi(nrep)
      dimension rprop(nrep,ip), rpr(ip), roldr(ip)
      dimension rrpr(nboot,ip), rloc(nboot,ip), rsca(nboot,ip)
      dimension pesi(nrep,nsize),nsame(nrep)
      dimension dden(nrep,nsize),dmod(nrep,nsize)
      dimension ddelta(nrep,nsize)
      dimension totloc(ip), totvar(ip), tot(ip)
      dimension rnowm(ip), rnows(ip), rnowhs(ip), dtau(nsize,ip)
      dimension rrnowm(ip), rrnows(ip), rrrpr(ip)
      dimension totw(ip), roldm(ip), rolds(ip), norder(ip)
      dimension diffm(ip), diffs(ip), diffp(ip)
      dimension ddati(nsize)

      do 30 i=1,nrep
         nsame(i)=0
 30   continue

      nconv=0 
      nsol=0

C Start bootstrapping

      do 900 iboot=1,nboot

      rnowh=dzero
      do 115 i=1,ip
         rnowm(i)=rloc(iboot,i)
         rnows(i)=rsca(iboot,i)
         rpr(i)=rrpr(iboot,i)
         rnowh=rnowh+rpr(i)*rnows(i)
 115  continue
         rnowh=rk*rnowh
      do 117 i=1,ip
         rnowhs(i)=rnowh+rnows(i)
 117  continue

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
       do 125 j=1,ip
       rm(i)=rm(i)+rpr(j)*dexp(-((dati(i)-rnowm(j))**ddue)/
     &    (ddue*rnowhs(j)))/dsqrt(rnowhs(j))
 125   continue

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
             rw(i)= duno - (delta(i) / (delta(i)+ddue))**ddue
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

      do 141 j=1,ip
         tot(j)=dzero
         totw(j)=dzero
         totloc(j)=dzero
         totvar(j)=dzero
 141  continue

      do 142 i=1,nsize
         dtemp=dzero
         do 145 j=1,ip
           dtau(i,j)=rpr(j)*dexp(-((dati(i)-rnowm(j))**ddue)/
     &     (ddue*rnows(j)))/dsqrt(rnows(j))
           dtemp=dtemp+dtau(i,j)
 145     continue
         do 147 j=1,ip
           dtau(i,j)=dtau(i,j)/dtemp
 147     continue
 142  continue

      do 155 j=1,ip
         do 150 i=1,nsize
            totw(j)=totw(j)+rw(j)
            tot(j)=tot(j)+rw(i)*dtau(i,j)
            totloc(j)=totloc(j)+rw(i)*dtau(i,j)*dati(i)
 150     continue
            roldm(j)=rnowm(j)
            rnowm(j)=totloc(j)/tot(j)
 155  continue

      do 165 j=1,ip
         do 160 i=1,nsize
            totvar(j)=totvar(j)+rw(i)*dtau(i,j)*
     &                ((dati(i)-rnowm(j))**2)
 160     continue
         rolds(j)=rnows(j)
         rnows(j)=totvar(j)/tot(j)
 165  continue

      totrpr=dzero
      do 1160 j=1,ip
         roldr(j)=rpr(j)
         rpr(j)=tot(j)/totw(j)
         if (rpr(j).lt.dzero) then
             rpr(j)=rprec/ddue
         endif
         totrpr=totrpr+rpr(j)
 1160 continue

      do 1170 j=1,ip
         rpr(j)=rpr(j)/totrpr
 1170 continue

      rnowh=dzero
      do 167 j=1,ip
         rnowh=rnowh+rpr(j)*rnows(j)
 167  continue
         rnowh=rk*rnowh
      do 169 j=1,ip
         rnowhs(j)=rnowh+rnows(j)
 169  continue

      do 1180 j=1,ip
         diffm(j)=abs(rnowm(j)-roldm(j))
         diffs(j)=abs(rnows(j)-rolds(j))
         diffp(j)=abs(rpr(j)-roldr(j))
 1180 continue

         call wlemax(diffm,ip,ddiffm)
         call wlemax(diffs,ip,ddiffs)
         call wlemax(diffp,ip,ddiffp)

      if(iter.gt.imax) then 
         nconv=nconv+1
         goto 8888 
      endif  

C      write(*,*) ddiffm
C      write(*,*) ddiffs
C      write(*,*) ddiffp


      if(ddiffm.gt.rprec.or.ddiffs.gt.rprec.or.
     &   ddiffp.gt.rprec) then
         goto 9999
      endif   

C Convergence achieved

C   
C Did you find all the components?
C

      if (iall.eq.1) then

          icomp=0
          ipos=0

C          write(*,*) rnowm
C          write(*,*) rnows
C          write(*,*) rpr

          do 1100 i=1,(ip-1)
             do 1110 j=(i+1),ip
                if (abs(rnowm(i)-rnowm(j)).lt.requal) then
                    if (abs(rnows(i)-rnows(j)).lt.requal) then
                        rpr(i)=rpr(i)+rpr(j)
                        rpr(j)=dzero
                        icomp=icomp+1
                        ipos=j
                    endif
                endif
 1110        continue
 1100     continue   


C          write(*,*) rnowm
C          write(*,*) rnows
C          write(*,*) rpr

          if (icomp.gt.0) then

C          write(*,*) 'proviamo a trovare tutte le componenti'

              j=0
              do 1120 i=1,nsize
                 if (rw(i).lt.dmin) then
                     j=j+1
                     ddati(j)=dati(i)
                 endif
 1120         continue
                 if (j.gt.2) then
                     call wlenorm(ddati,j,nsize,3,2,1,
     &                            iraf,rk,rprec,requal,imax,
     &                            tm,tv,ttot,tpesi,
     &                            tden,tmod,tdelta,ntsame,ntsol,
     &                            ntconv)
                    if (ntsol.eq.1) then
                        rnowm(ipos)=tm
                        rnows(ipos)=tv
                        rpr(ipos)=j/nsize
                        dtot=duno+rpr(ipos)
                        do 1130 k=1,ip
                           rpr(k)=rpr(k)/dtot
 1130                   continue

C          write(*,*) 'aggiunta una componente'
C          write(*,*) ipos
C          write(*,*) rnowm
C          write(*,*) rnows
C          write(*,*) rpr

                        goto 9999
                    endif
C fine ntsol.eq.1

                 endif
C fine j.gt.2

          endif
C fine icomp.gt.0


      endif

C
C Is this a new root, then store it
C

      call wleord1(rpr,ip,rrrpr,norder)

C      write(*,*) rpr
C      write(*,*) norder

      call wleord2(rnowm,ip,rrnowm,norder)
      call wleord2(rnows,ip,rrnows,norder)

C      write(*,*) rrnowm
C      write(*,*) rrnows

      if(nsol.eq.0) then
         nsol=nsol+1
         nsame(1)=1 
         do 1190 j=1,ip
         dmedia(nsol,j)=rrnowm(j)
         varia(nsol,j)=rrnows(j)
         rprop(nsol,j)=rrrpr(j)
 1190    continue
         totpesi(nsol)=totw(1)
         do 170 i=1,nsize
            pesi(nsol,i)=rw(i)
            dden(nsol,i)=d(i)
            dmod(nsol,i)=rm(i)
            ddelta(nsol,i)=delta(i)
 170     continue
      else
         do 180 isol=1,nsol

            do 1200 j=1,ip
               diffm(j)=abs(rrnowm(j)-dmedia(isol,j))
               diffs(j)=abs(rrnows(j)-varia(isol,j))
               diffp(j)=abs(rrrpr(j)-rprop(isol,j))
 1200        continue

             call wlemax(diffm,ip,ddiffm)
             call wlemax(diffs,ip,ddiffs)
             call wlemax(diffp,ip,ddiffp)

C      write(*,*) 'nuova radice'
C      write(*,*)  rrnowm
C      write(*,*)  rrnows
C      write(*,*)  rrrpr
 
C      write(*,*) 'valore confrontato'
C      write(*,*) dmedia(isol,1),dmedia(isol,2),dmedia(isol,3)
C      write(*,*) varia(isol,1),varia(isol,2),varia(isol,3)
C      write(*,*) rprop(isol,1),rprop(isol,2),rprop(isol,3)

C      write(*,*) 'differenze'
C      write(*,*) diffm
C      write(*,*) diffs
C      write(*,*) diffp

C      write(*,*) 'Sono radici uguali'
C      write(*,*) ddiffm
C      write(*,*) ddiffs
C      write(*,*) ddiffp

             if (ddiffm.lt.requal.and.ddiffs.lt.requal
     &          .and.ddiffp.lt.requal) then
                   nsame(isol)=nsame(isol)+1
                   goto 190
             endif
 180     continue
            nsol=nsol+1
            nsame(nsol)=1
            do 1210 j=1,ip
               dmedia(nsol,j)=rrnowm(j)
               varia(nsol,j)=rrnows(j)
               rprop(nsol,j)=rrrpr(j)
 1210        continue
             totpesi(nsol)=totw(1)
             do 200 i=1,nsize
                pesi(nsol,i)=rw(i)
                dden(nsol,i)=d(i)
                dmod(nsol,i)=rm(i)
                ddelta(nsol,i)=delta(i)
 200         continue
 190      continue
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE wleord1 (d,n,temp,norder)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)
 
      dimension d(n)
      dimension dd(n)
      dimension temp(n)
      dimension norder(n)

      do 30 k=1,n
            dd(k)=d(k)
 30   continue

      call wlemax(dd,n,ddmax)

      do 10 j=1,n

      call wlemin(dd,n,ddmin)
      temp(j)=ddmin

      do 20 i=1,n
          if (ddmin.eq.dd(i)) then
              norder(j)=i
              dd(i)=ddmax+10.0d00
              goto 10
          endif
 20   continue
 10   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE wleord2 (d,n,temp,norder)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)
 
      dimension d(n)
      dimension temp(n)
      dimension norder(n)

      do 10 j=1,n
      temp(j)=d(norder(j))
 10   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      SUBROUTINE wlemax (d,n,ddmax)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)
 
      dimension d(n)

      ddmax=d(1)
      do 10 i=2,n
          if (ddmax.lt.d(i)) then 
              ddmax=d(i)
          endif
 10   continue

      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
     
      SUBROUTINE wlemin (d,n,ddmin)

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)
 
      dimension d(n)

      ddmin=d(1)
      do 10 i=2,n
          if (ddmin.gt.d(i)) then 
              ddmin=d(i)
          endif
 10   continue

      return
      end



      SUBROUTINE DECTOBIN (IDEC,IDIM,ibin)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Decimal number to binary number
C     
C     Author: Claudio Agostinelli 
C             Dipartimento di Statistica
C             Universita' di Padova
C             35121 Padova
C             ITALIA
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
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PARAMETER:
C     NAME:     I/O:    TYPE:  DIMENSIONS:   DESCRIPTIONS:
C     IDEC      input    I     1             integer in decimal format
C     IDIM      input    I     1             max number of digit of the binary format
C     ibin      output   I     IDIM          vector represent the integer number in binary format

      implicit double precision(a-h,o-z)
      implicit integer (n,i,j)

      dimension ibin(idim)
 
      dresult=idec

      do 10 i=idim,1,-1         
         ibin(i)=int(dresult/dble(2**(i-1)))
         dresult=dresult-dble((2**(i-1))*ibin(i))
 10   continue
      return
      end

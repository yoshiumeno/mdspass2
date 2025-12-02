!C**********************************************************************
!C                        INITIAL-NANO-ARMCHAIR
!C**********************************************************************
      IMPLICIT DOUBLE PRECISION(A,B,D-H,O-Z),INTEGER(I-N),COMPLEX*16(C)
      parameter(na=4,man=20000)
      dimension :: x(na),y(na),z(na), &
           xnt(man),ynt(man),znt(man), &
           a1(3),a2(3),a3(3), &
           cr(na),nw(30) 
      character (3) :: strm,strn

      natom=1
      iatom=0
      ao=1.42d0
      print *, "enter nwall"
      read *, nwall

      do 2000 iwall=1,nwall
      print *, "enter m"
      read *, n1
      print *, "enter n"
      read *, n2
      print *, "enter nlayer"
      read *, nlayer

      print '(1x,"TUBE_TYPE ",5x,"(",i3,",",i3,")")',n1,n2

!c      katom=na

      pi=4.0D0*DATAN(1.0D0)  !±ß¼þÎ¨¦Ð


      a1(1)=ao*dsqrt(3.0d0)
      a1(2)=0.0d0
      a1(3)=0.0d0

      a2(1)=0.0d0
      a2(2)=ao*3.0d0
      a2(3)=0.0d0

      a3(1)=0.0d0
      a3(2)=0.0d0
      a3(3)=6.696d0/1.42d0*ao
                                 !c=6.696 ¢ò¡Ê£úÊý¸þ¡Ë
                                 !c=1.42  ¢ò¡Ê£ø¡¢£ùÊý¸þ¡Ë

!c     ¸¶»Ò¤Î½é´üÇÛÃÖ
      x(1)=0.0d0
      y(1)=0.0d0
      z(1)=0.0d0

      x(2)=0.0d0
      y(2)=ao
      z(2)=0.0d0

      x(3)=ao*dsqrt(3.0d0)/2.0d0
      y(3)=-ao/2.0d0
      z(3)=0.0d0

      x(4)=ao*dsqrt(3.0d0)/2.0d0
      y(4)=ao*3.0d0/2.0d0
      z(4)=0.0d0


      at11=a1(1)
      at12=0.0d0
      at21=a1(1)/2.0d0
      at22=-ao*3.0d0/2.0d0

      atube1=dble(at11*n1)+dble(at21*n2)
      atube2=dble(at12*n1)+dble(at22*n2)
      rtubex=dsqrt(atube1**2.0d0+atube2**2.0d0)

      tubecos=atube1/rtubex
      tubesin=atube2/rtubex

      ctube=dcmplx(tubecos,tubesin)
      do i=1,na
         cr(i)=dcmplx(x(i),y(i))*dconjg(ctube)
         x(i)=dreal(cr(i))
         y(i)=dimag(cr(i))
         write(90,*) x(i),y(i),z(i)
      enddo

      nn1=2*n1+n2
      nn2=2*n2+n1
      call lcm(nn1,nn2,nc)
!c      print *, n1,n2,nc

      as1=dble(at11*nc/nn1)-dble(at21*nc/nn2)
      as2=dble(at12*nc/nn1)-dble(at22*nc/nn2)
      rtubey=dsqrt(as1**2.0d0+as2**2.0d0)
!c      print *, nc/nn1,nc/nn2

      ctrans=dcmplx(a1(1),a1(2))
      ctrans=ctrans*dconjg(ctube)
      a1(1)=dreal(ctrans)
      a1(2)=dimag(ctrans)
      ctrans=cmplx(a2(1),a2(2))
      ctrans=ctrans/ctube
      a2(1)=dreal(ctrans)
      a2(2)=dimag(ctrans)




      do i=1,na
         do ix=-100,100
            do iy=-100,100
               xx=x(i)+ix*a1(1)+iy*a2(1)
               yy=y(i)+ix*a1(2)+iy*a2(2)
               if(xx.ge.0.0d0-1.0d-2) then
                  if(xx.lt.rtubex-1.0d-2) then
                     if(yy.ge.0.0d0-1.0d-2) then
                        if(yy.lt.rtubey-1.0d-2) then
                           iatom=iatom+1
                           xnt(iatom)=xx
                           ynt(iatom)=yy
                           znt(iatom)=0.0d0
                        endif
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
 

!c      open(10,file='cood-gra.dat')
!c      do i=1,iatom
!c         write(10,*) xnt(i),ynt(i),znt(i)
!c      enddo
!c      write(10,*) "  "
!c      write(10,*) 0.0d0,0.0d0
!c      write(10,*) rtubex,0.0d0
!c      write(10,*) rtubex,rtubey
!c      write(10,*) 0.0d0,rtubey
!c      write(10,*) 0.0d0,0.0d0

      do i=natom,iatom
         znt(i)=ynt(i)
!c         znt(i)=xnt(i)
         xx=xnt(i)
!c         xx=ynt(i)

!c         xnt(i)=rtubex/(2.0d0*pi)*dcos(2.0d0*pi*xx/rtubex)
!c         ynt(i)=rtubex/(2.0d0*pi)*dsin(2.0d0*pi*xx/rtubex)

         xnt(i)=rtubex/(2.0d0*pi)*dcos(2.0d0*pi*xx/rtubex+2.0d0*pi/96.0d0)
         ynt(i)=rtubex/(2.0d0*pi)*dsin(2.0d0*pi*xx/rtubex+2.0d0*pi/96.0d0)
      enddo
      print '(1x,"DIAMETER",d16.8,1x,"nm")', rtubex/pi/10.0d0
      nw(iwall)=iatom*nlayer
      natom=iatom+1
 2000 continue
      katom=iatom
      do i=1,nlayer-1
         do j=1,katom
            iatom=iatom+1
            xnt(iatom)=xnt(j)
            ynt(iatom)=ynt(j)
            znt(iatom)=znt(j)+rtubey*i
         enddo
      enddo



      open(50,file='r-nano.dat')
      open(40,file='trans-v-nano.dat')

      print '(1x,"NUMBER_OF_ATOMS",i9)', iatom

      a1(1)=rtubex*2.0d0
      a1(2)=0.0d0
      a1(3)=0.0d0
      a2(1)=0.0d0
      a2(2)=rtubex*2.0d0
      a2(3)=0.0d0
      a3(1)=0.0d0
      a3(2)=0.0d0
      a3(3)=rtubey*nlayer
      
      write(40,222) a1(1),a1(2),a1(3)
      write(40,222) a2(1),a2(2),a2(3)
      write(40,222) a3(1),a3(2),a3(3)

      do i=1,iatom
         write(50,222) xnt(i),ynt(i),znt(i)
      enddo

      open(60,file='tube.xyz')

      call conv_num_char(n1,strm)
      call conv_num_char(n2,strn)
      write(60,'(i10)') iatom
      write(60,'(a20)') 'Nanotube'//strm//strn
      do i=1,iatom

!c         xnt2=xnt(i)
!c         ynt2=ynt(i)
!c         znt2=znt(i)
!c         if(abs(xnt(i)).lt.1.0d-5) xnt2=0.0d0
!c         if(abs(ynt(i)).lt.1.0d-5) ynt2=0.0d0
!c         if(abs(znt(i)).lt.1.0d-5) znt2=0.0d0
!c         write(60,*) "c",xnt2,ynt2,znt2

         write(60,111) xnt(i),ynt(i),znt(i)
         
      enddo

 111  format("c",e16.8,e16.8,e16.8)
!c 112  format("o",d16.8,d16.8,d16.8)
!c 111  format("12",d26.16,d26.16,d26.16)
 222  format(1x,d26.16,d26.16,d26.16)

      end




      subroutine lcm(n1,n2,nc)
      IMPLICIT DOUBLE PRECISION(A,B,D-H,O-Z),INTEGER(I-N),COMPLEX*16(C)

      dif=n1-n2
      if(dif.lt.0) then
         m1=n1
         m2=n2
      else
         m1=n2
         m2=n1
      endif
      mm1=m1
      nc=1
!c      print *, m1,m2
 1000 continue
      do i=2,mm1
         if(mod(m1,i).eq.0) then
            if(mod(m2,i).eq.0) then
!c               print *, i,m1,m2
               m1=m1/i
               m2=m2/i
               nc=nc*i
!c               print *, m1,m2
               go to 1000
            endif
         endif
      enddo
      nc=nc*m1*m2

      return
      end
               
      subroutine conv_num_char(ll,str)
      integer :: l
      character (3) :: str
      l=abs(ll)
      ibase=ichar('0')
      iii=int(l/100)
      ii=int((l-iii*100)/10)
      i=int((l-iii*100-ii*10)/1)
      if (ll<0) stop 'conv_num_char error'
      str=char(ibase+iii)//char(ibase+ii)//char(ibase+i)
      end subroutine conv_num_char

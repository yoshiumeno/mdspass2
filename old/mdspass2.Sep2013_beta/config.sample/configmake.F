! This program creates 'CONFIG' file.
!
! For your information:
!  Al fcc a=4.05 ang
!  Cu fcc a=3.61
!  Ni fcc a=3.52
!  Fe bcc a=2.87
!  Au fcc a=4.08
!  Si dia a=5.43
!
      program configmake
      implicit double precision (a-h,o-z)

!      integer :: iix = 90
      integer :: iix = 5
      integer :: iiy = 5
      integer :: iiz = 5
      character (2) :: asp = 'Al' ! Atom species
      character (15) :: potential_func = 'EAM Mishin' ! Potential type
      character (3) :: structure = 'fcc' ! Crystal structure
!      character (3) :: fmt = 'abs' ! Format type
      character (3) :: fmt = 'fra' ! Format type
      double precision :: a = 4.05 ! Lattice constant (ang)

      double precision :: ang= 1.0d-10 ! 1ang=10^-10 m
      double precision :: hmat(3,3)
      double precision, allocatable :: rx(:),ry(:),rz(:) ! Atom positions

      a=a*ang
      if (structure=='fcc') then
         n=4*iix*iiy*iiz
      else if (structure=='bcc') then
         n=2*iix*iiy*iiz
      else if (structure=='fc1') then
         n=3*iix*iiy*iiz
      else if (structure=='dia') then
         n=8*iix*iiy*iiz
      endif

      allocate (rx(n),ry(n),rz(n))

      if (structure=='fcc') then
         rx(1)=0.0d0
         ry(1)=0.0d0
         rz(1)=0.0d0
         rx(2)=0.0d0
         ry(2)=a/2.0d0
         rz(2)=a/2.0d0
         rx(3)=a/2.0d0
         ry(3)=0.0d0
         rz(3)=a/2.0d0
         rx(4)=a/2.0d0
         ry(4)=a/2.0d0
         rz(4)=0.0d0
         ncell=4
         hmat=0.0d0
         hmat(1,1)=a
         hmat(2,2)=a
         hmat(3,3)=a
      else if (structure=='dia') then
         rx(1)=0.0d0
         ry(1)=0.0d0
         rz(1)=0.0d0
         rx(2)=0.0d0
         ry(2)=a/2.0d0
         rz(2)=a/2.0d0
         rx(3)=a/2.0d0
         ry(3)=0.0d0
         rz(3)=a/2.0d0
         rx(4)=a/2.0d0
         ry(4)=a/2.0d0
         rz(4)=0.0d0
         do i=5,8
            rx(i)=rx(i-4)+a/4.0d0
            ry(i)=ry(i-4)+a/4.0d0
            rz(i)=rz(i-4)+a/4.0d0
         enddo
         ncell=8
         hmat=0.0d0
         hmat(1,1)=a
         hmat(2,2)=a
         hmat(3,3)=a
      else if (structure=='bcc') then
         rx(1)=0.0d0
         ry(1)=0.0d0
         rz(1)=0.0d0
         rx(2)=a/2.0d0
         ry(2)=a/2.0d0
         rz(2)=a/2.0d0
         ncell=2
         hmat=0.0d0
         hmat(1,1)=a
         hmat(2,2)=a
         hmat(3,3)=a
      else if (structure=='fc1') then
         rx(1)=0.0d0
         ry(1)=0.0d0
         rz(1)=0.0d0
         rx(2)=a/sqrt(2.0)/2.0d0
         ry(2)=a*sqrt(3.0/2.0)/2.0/3.0
         rz(2)=a*sqrt(3.0)/3.0d0
         rx(3)=a/sqrt(2.0)
         ry(3)=a*sqrt(3.0/2.0)/2.0/3.0*2.0
         rz(3)=a*sqrt(3.0)/3.0d0*2.0
         ncell=3
         hmat=0.0d0
         hmat(1,1)=a/sqrt(2.0)    ; hmat(2,1)=0.0; hmat(3,1)=0.0
         hmat(1,2)=a/sqrt(2.0)/2.0; hmat(2,2)=a*sqrt(3.0/2.0)/2.0; hmat(3,2)=0.0
         hmat(1,3)=0.0; hmat(2,3)=0.0; hmat(3,3)=a*sqrt(3.0)
      endif
         

      iii=ncell
      jjj=iii
      if (iix.ge.2) then
         do k=2,iix
            do j=1,jjj
               iii=iii+1
               rx(iii)=rx(iii-jjj)+hmat(1,1)
               ry(iii)=ry(iii-jjj)+hmat(2,1)
               rz(iii)=rz(iii-jjj)+hmat(3,1)
            enddo
         enddo
         jjj=iii
      endif
      if (iiy.ge.2) then
         do k=2,iiy
            do j=1,jjj
               iii=iii+1
               rx(iii)=rx(iii-jjj)+hmat(1,2)
               ry(iii)=ry(iii-jjj)+hmat(2,2)
               rz(iii)=rz(iii-jjj)+hmat(3,2)
            enddo
         enddo
         jjj=iii
      endif
      if (iiz.ge.2) then
         do k=2,iiz
            do j=1,jjj
               iii=iii+1
               rx(iii)=rx(iii-jjj)+hmat(1,3)
               ry(iii)=ry(iii-jjj)+hmat(2,3)
               rz(iii)=rz(iii-jjj)+hmat(3,3)
            enddo
         enddo
      endif

      hmat(1:3,1)=hmat(1:3,1)*dble(iix)
      hmat(1:3,2)=hmat(1:3,2)*dble(iiy)
      hmat(1:3,3)=hmat(1:3,3)*dble(iiz)

      open (30,FILE='CONFIG')
      write (30,*) asp
      write (30,*) potential_func
      write (30,*) n
      write (30,*) a/ang
      do j=1,3
         xx=hmat(1,j)/a; yy=hmat(2,j)/a; zz=hmat(3,j)/a
         write (30,900) xx,yy,zz
      enddo
      write (30,*) fmt
      do i=1,n
!         xx=rx(i)/ang; yy=ry(i)/ang; zz=rz(i)/ang
         xx=rx(i)/hmat(1,1); yy=ry(i)/hmat(2,2); zz=rz(i)/hmat(3,3)

! default
         write (30,900) xx,yy,zz

! example 1
!         if (zz<0.5) then
!            write (30,920) xx,yy,zz,'Si','F','F','F'
!         else
!            write (30,910) xx,yy,zz,'Si'
!         endif

! example 2
!         if (mod(i-1,8)<4) then
!            write (30,910) xx,yy,zz,'Si'
!         else
!            write (30,910) xx,yy,zz,'C'
!         endif
      enddo
      close (30)
 900  format (3f20.15)
 910  format (3f20.15,3x,a2)
 920  format (3f20.15,3x,4a3)

      end program

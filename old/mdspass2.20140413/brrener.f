*************************************************************
***                              ****************************
***   (Brennerポテンシャル,初期温度298K)*****
***   (粒子登録法、領域分割法、リストベクトル)***************
***   (速度スケーリング法  温度298K )           **************
*************************************************************

c     原子の個数

      SUBROUTINE brrener(N,X,Y,Z,XN,YN,ZN,XO,YO,ZO
     &            ,a1,a2,a3,irx,iry,irz
     &            ,vol,a,rs,nkousi
     &            ,istep,sumt,enkin,tt,enpot
     &            ,fx,fy,fz,lvi,lvj
     &            ,num,num1,num2,numa,dist,epsilon,ep,kind
     &            )

C      SUBROUTINE brrener(N,X,Y,Z,XN,YN,ZN,XO,YO,ZO
C     &            ,a1,a2,a3,gamma,xstr,ystr,zstr
C     &            ,vol,a,rs,nkousi
C     &            ,istep,sumt,enkin,tt,enpot
C     &            ,fx,fy,fz,num
C     &            ,lf,dilf,djlf,nc,dnc,h1,h2,hx,hy
C     &            ,fsij,fsik,fsjk,fskl,flik,nci,ncj,nconj
C     &            ,ncon,sumbij,sumbji
C     &            ,sp,disp,djsp
C     &            ,dbijx,dbjix,dbijy,dbjiy,dbijz,dbjiz
C     &            ,lvi,lvj,XX,YY,c,enp,dis,la,lb,lm,alpha
C     &            )

      implicit double precision (a-h,o-z)
      double precision X(N),Y(N),Z(N)
      double precision XN(N),YN(N),ZN(N),XO(N),YO(N),ZO(N)
      double precision a1(3),a2(3),a3(3)
      double precision dtx(n),dty(n),dtz(n)
      double precision vx(n),vy(n),vz(n)
      double precision px(n),py(n),pz(n)
      double precision fx(n),fy(n),fz(n)

      double precision lf(0:10,0:10,0:10)
      double precision dilf(0:10,0:10,0:10)
      double precision djlf(0:10,0:10,0:10),nc(n)
     &                ,dncx(n),dncy(n),dncz(n)
      double precision h1(0:10),h2(0:10)
      double precision hx(0:10),hy(0:10)
      double precision fsij,fsik,fsjk,fskl,flik
      double precision nci,ncj,nconj,ncon(n,n)
      double precision sumbij(n,n,-irx:irx,-iry:iry,-irz:irz)
      double precision sumbji(n,n,-irx:irx,-iry:iry,-irz:irz)
      double precision sp(n,n,-irx:irx,-iry:iry,-irz:irz)
     &                ,disp(n,n,-irx:irx,-iry:iry,-irz:irz)
     &                ,djsp(n,n,-irx:irx,-iry:iry,-irz:irz)
      double precision dbijx(n,n,-irx:irx,-iry:iry,-irz:irz)
     &                ,dbjix(n,n,-irx:irx,-iry:iry,-irz:irz)
      double precision dbijy(n,n,-irx:irx,-iry:iry,-irz:irz)
     &                ,dbjiy(n,n,-irx:irx,-iry:iry,-irz:irz)
      double precision dbijz(n,n,-irx:irx,-iry:iry,-irz:irz)
     &                ,dbjiz(n,n,-irx:irx,-iry:iry,-irz:irz)


      dimension lvi(n**2),lvj(n**2)
      double precision XX(0:10),YY(0:10),c(0:10,3),enp(2),dis(n**2)
      double precision la(0:10),lb(0:10),lm(0:10),alpha(0:10)

      integer cx,cy,cz,cx1,cy1,cz1

      character*1 num(0:9)
      character*4 kind

      clx=a1(1)
      cly=a2(2)
      clz=a3(3)

      nl=n**2
      nf1=10
      nf2=10

      num(0)='0'
      num(1)='1'
      num(2)='2'
      num(3)='3'
      num(4)='4'
      num(5)='5'
      num(6)='6'
      num(7)='7'
      num(8)='8'
      num(9)='9'


c     ファイルオープン

c-----変数の設定-----------------------------------------------

      dt=1.0d-15           !時間ステップ

c     ポテンシャルのパラメータ
      do  i=0,nf1
         do  j=0,nf1
            do  k=0,nf2
               lf(i,j,k)=0.0d0
               dilf(i,j,k)=0.0d0
               djlf(i,j,k)=0.0d0
            enddo
         enddo
      enddo
c----------------------------パラメータ1
C      re=1.315d0
C      de=6.325d0
C      be=1.5d0
C      s=1.29d0
C      th=0.80469d0
C      r1=1.7d0
C      r2=2.0d0
C      a0=0.011304d0
C      c2=(19.0d0)**2
C      d2=(2.5d0)**2
C
C      lf(1,1,1)=0.1511d0
C      lf(2,2,1)=0.075d0
C      lf(2,3,1)=-0.0465d0
C      lf(2,3,2)=-0.0465d0
C      lf(1,2,1)=0.0126d0
C      lf(1,2,2)=-0.0355d0
C      lf(1,3,1)=-0.1130d0
C      lf(0,3,1)=-0.1220d0
C      lf(0,2,1)=0.0320d0
C      lf(0,1,1)=0.1100d0
C      lf(1,3,2)=-0.1120d0
C      lf(0,3,2)=-0.1220d0
C      lf(0,2,2)=-0.0445d0
C      lf(1,1,2)=0.0074d0
C
C      dilf(3,1,1)=-0.1160d0
C      dilf(3,2,1)=-0.13205d0
C      dilf(3,1,2)=-0.0610d0
C      dilf(2,3,2)=0.02225d0
C      dilf(2,4,2)=-0.03775d0
C      dilf(3,4,2)=0.0565d0
C      dilf(3,4,1)=0.0565d0
C      dilf(3,2,2)=-0.1065d0
Cc----------------------------パラメータ2
      re=1.39d0
      de=6.0d0
      be=2.1d0
      s=1.22d0
      th=0.5d0
      r1=1.7d0
      r2=2.0d0
      a0=0.00020813d0
      c2=(330.0d0)**2
      d2=(3.5d0)**2

      lf(1,1,1)=0.1264d0
      lf(2,2,1)=0.0605d0
      lf(2,3,1)=-0.0363d0
      lf(2,3,2)=-0.0363d0
      lf(1,2,1)=0.0120d0
      lf(1,2,2)=-0.0243d0
      lf(1,3,1)=-0.0903d0
      lf(0,3,1)=-0.0904d0
      lf(0,2,1)=0.0427d0
      lf(0,1,1)=0.0996d0
      lf(1,3,2)=-0.0903d0
      lf(0,3,2)=-0.0904d0
      lf(0,2,2)=-0.0269d0
      lf(1,1,2)=0.0108d0

      dilf(3,1,1)=-0.0950d0
      dilf(3,2,1)=-0.10835d0
      dilf(3,1,2)=-0.0452d0
      dilf(2,3,2)=0.01345d0
      dilf(2,4,2)=-0.02705d0
      dilf(3,4,2)=0.04515d0
      dilf(3,4,1)=0.04515d0
      dilf(3,2,2)=-0.08760d0

c------------------------------------------
      do  i=nf1,1,-1
         do  j=i-1,0,-1
            do  k=0,nf2
               lf(i,j,k)=lf(j,i,k)
               djlf(j,i,k)=dilf(i,j,k)
            enddo
         enddo
      enddo

      do  i=1,nf1
         do  j=1,nf1
            do  k=0,nf2
               djlf(j,i,k)=dilf(i,j,k)
            enddo
         enddo
      enddo

      do  i=0,nf1
         do  j=0,nf1
            do  k=3,nf2
               lf(i,j,k)=lf(j,i,2)
            enddo
         enddo
      enddo
c      do 5000 nnk=0,1000
c      rc=8.0d-9
      rc=r2   !切断距離
      rc2=rc**2       ! 二乗を作っておく(計算時間短縮のため)

      lstep=50        !粒子登録法のサイクル

      enk=0.0d0       !系の運動エネルギーの平均値

c-----設定温度
C      t=1.0d2
      t=0.0d0
      if(istep.eq.1)sumt=0.0d0
c-------------
      bl=1.380662d-23      !ボルツマン定数
      pi=3.14259265d0      !円周率

      am=1.6726485d-27     !陽子・中性子の質量
      bm=1.6749543d-27
      wm=am*6.0d0+bm*6.0d0 !原子の質量

c-----リストベクトルの登録距離
      rcl=rc+dsqrt((3.0d0*bl*t)/wm)*dt*2.0d0
      rcl2=rcl**2   ! 二乗を作っておく(計算時間短縮のため)


c         print *,istep
c-----------------------------------continue
C         if (istep.eq.1)then
C         call continue(N,X,Y,Z,xn,yn,zn,xo,yo,zo,istep
C     &                ,a1,a2,a3,num,num1,num2,numa,dist
C     &                ,enpot,epsilon,ep,kind,dt,sumt)
C         endif

c-------------------------------------------continue
      if (istep.eq.1) then
      do  k=1,n
         vx(k)=0.0d0
         vy(k)=0.0d0
         vz(k)=0.0d0
      enddo

      call velocity(n,vx,vy,vz,wm,t)

      do   k=1,n
         xo(k)=x(k)-vx(k)*dt
         yo(k)=y(k)-vy(k)*dt
         zo(k)=z(k)-vz(k)*dt
      enddo
      endif
Cc-----------------------------------
         if(istep.eq.1) go to 311
C         do k=1,n
C       if(mod(istep,2000).eq.1) write(7,*) istep,k
C     &                                    ,x(k),y(k),z(k)
C     &                                    ,xo(k),yo(k),zo(k)
C       enddo
         CALL RELAX(N,X,Y,Z,xn,yn,zn,xo,yo,zo,FX,FY,FZ
     &             ,dt,wm,bl,sumt,enkin,t,tt,istep
     &             ,a1,a2,a3,px,py,pz,dtx,dty,dtz)
                               !力から変位を求める
c-----基本セルから飛び出した原子を修正
Cc--------ナノチューブ
C         do k=1,n
C
C           if (xn(k).le.-a1(1)/2.0d0)then
C            xn(k)=xn(k)+a1(1)
C            x(k)=x(k)+a1(1)
C           endif
C           if (xn(k).gt.a1(1)/2.0d0)then
C            xn(k)=xn(k)-a1(1)
C            x(k)=x(k)-a1(1)
C           endif
C
C           if (yn(k).le.-a2(2)/2.0d0)then
C            yn(k)=yn(k)+a2(2)
C            y(k)=y(k)+a2(2)
C           endif
C           if (yn(k).gt.a2(2)/2.0d0)then
C            yn(k)=yn(k)-a2(2)
C            y(k)=y(k)-a2(2)
C           endif
C
C           if (zn(k).le.0.0d0)then
C            zn(k)=zn(k)+a3(3)
C            z(k)=z(k)+a3(3)
C           endif
C           if (zn(k).gt.a3(3))then
C            zn(k)=zn(k)-a3(3)
C            z(k)=z(k)-a3(3)
C           endif
C         enddo
c---------グラファイト
         do k=1,n

           if (xn(k).le.-a1(1)/2.0d0)then
            xn(k)=xn(k)+a1(1)
            x(k)=x(k)+a1(1)
           endif
           if (xn(k).gt.a1(1)/2.0d0)then
            xn(k)=xn(k)-a1(1)
            x(k)=x(k)-a1(1)
           endif

           if (yn(k).le.-a2(2)/2.0d0)then
            yn(k)=yn(k)+a2(2)
            y(k)=y(k)+a2(2)
           endif
           if (yn(k).gt.a2(2)/2.0d0)then
            yn(k)=yn(k)-a2(2)
            y(k)=y(k)-a2(2)
           endif

           if (zn(k).le.-a3(3)/2.0d0)then
            zn(k)=zn(k)+a3(3)
            z(k)=z(k)+a3(3)
           endif
           if (zn(k).gt.a3(3)/2.0d0)then
            zn(k)=zn(k)-a3(3)
            z(k)=z(k)-a3(3)
           endif
         enddo

c-----原子配置のデータを一つずらす
      do k=1,n
        xo(k)=x(k)
        yo(k)=y(k)
        zo(k)=z(k)
      enddo

      do k=1,n
        x(k)=xn(k)
        y(k)=yn(k)
        z(k)=zn(k)
      enddo
C      do k=1,n
C       if(mod(istep,2000).eq.0) write(7,*) istep,k
C     &                                    ,x(k),y(k),z(k)
C     &                                    ,xo(k),yo(k),zo(k)
C      enddo

C 33      if (istep.ge.2000)then
C         if(mod(istep,2000).eq.1) then
C         call tensile(N,X,Y,Z,xn,yn,zn,xo,yo,zo,istep
C     &                ,a1,a2,a3,num,num1,num2,numa,dist
C     &                ,enpot,enkin,epsilon,ep,kind,sumt)
C         endif
C         endif

c--------------------
 311  irangex=irx
      irangey=iry
      irangez=irz

      if(mod(istep,lstep).eq.1)then
         call bookkeep(n,nl,lvi,lvj,X,Y,Z,XN,YN,ZN,XO,YO,ZO,rcl,
     &     rcl2,lvn,clx,cly,clz,irangex,irangey,irangez)
      endif

C      do  ln=1,lvn
C            print *,lvi(ln),lvj(ln)
C       enddo



c-----原子に働く力の計算
      do  i=1,n   !  FORCEのリセット
         fx(i)=0.0d0
         fy(i)=0.0d0
         fz(i)=0.0d0
         nc(i)=0.0d0
         dncx(i)=0.0d0
         dncy(i)=0.0d0
         dncz(i)=0.0d0

         do j=1,n
               do cx=-irangex,irangex
               do cy=-irangey,irangey
               do cz=-irangez,irangez
             sumbij(i,j,cx,cy,cz)=0.0d0          !<
             sumbji(i,j,cx,cy,cz)=0.0d0          !<
             dbijx(i,j,cx,cy,cz)=0.0d0
             dbjix(i,j,cx,cy,cz)=0.0d0
             dbijy(i,j,cx,cy,cz)=0.0d0
             dbjiy(i,j,cx,cy,cz)=0.0d0
             dbijz(i,j,cx,cy,cz)=0.0d0
             dbjiz(i,j,cx,cy,cz)=0.0d0
             enddo
             enddo
             enddo
         enddo
         sumdis=0.0d0
      enddo

      enpot=0.0d0               ! ポテンシャルエネルギーのリセット


      do 4700 i=1,n
         do 4710 j=1,n
            if (j.ne.i) then
               do cx=-irangex,irangex
               do cy=-irangey,irangey
               do cz=-irangez,irangez
               nna=nna+1                              ! iとjとの距離
                  xxr=x(i)-(x(j)+dble(cx)*clx)
                  yyr=y(i)-(y(j)+dble(cy)*cly)
                  zzr=z(i)-(z(j)+dble(cz)*clz)
                  if (((abs(xxr).le.rc).and.(abs(yyr).le.rc))
     &                    .and.(abs(zzr).le.rc)) then
                     r2=xxr**2+yyr**2+zzr**2 ! 距離の自乗
                     if (r2.lt.rc2) then
                        r=dsqrt(r2)
                        rr=r ! Åに直す

                        call fij (rr,fsij,dfsij)


                        nc(i)=nc(i)+fsij
                        dncx(i)=dncx(i)+dfsij*xxr/rr
                        dncy(i)=dncy(i)+dfsij*yyr/rr
                        dncz(i)=dncz(i)+dfsij*zzr/rr

                    endif
                  endif
               enddo
               enddo
               enddo
            endif
 4710    continue
 4700 continue





      do 5200 ln=1,lvn          ! リストベクトルによるループ

       i=lvi(ln)
       j=lvj(ln)

       do cx=-irangex,irangex
       do cy=-irangey,irangey
       do cz=-irangez,irangez

c     Ｖr(r）とＶa(r）を計算

                                       ! iとjとの距離
        xxr=x(i)-(x(j)+dble(cx)*clx)
        yyr=y(i)-(y(j)+dble(cy)*cly)
        zzr=z(i)-(z(j)+dble(cz)*clz)

        if (((abs(xxr).le.rc).and.(abs(yyr).le.rc))
     &              .and.(abs(zzr).le.rc)) then
         r2=xxr**2+yyr**2+zzr**2 ! 距離の自乗

         if (r2.lt.rc2) then
          r=dsqrt(r2)
          rr=r          ! Åに直す

          call fij (rr,fsij,dfsij)

          vr=fsij*de*exp(-dsqrt(2.0d0*s)*be*(rr-re))/(s-1.0d0)
          va=fsij*de*s*exp(-dsqrt(2.0d0/s)*be*(rr-re))/(s-1.0d0)

          sumni=0.0d0
          sumnj=0.0d0

          do 5300 k=1,n
            do cx1=-irangex,irangex
            do cy1=-irangey,irangey
            do cz1=-irangez,irangez
            if((k.eq.i).and.(cx1.eq.0)
     &                 .and.(cy1.eq.0).and.(cz1.eq.0))go to 7000
            if((k.eq.j).and.(cx1.eq.cx)
     &                 .and.(cy1.eq.cy).and.(cz1.eq.cz))go to 7000

                                        ! iとkとの距離
             xxi=x(i)-(x(k)+dble(cx1)*clx)
             yyi=y(i)-(y(k)+dble(cy1)*cly)
             zzi=z(i)-(z(k)+dble(cz1)*clz)

             rik2=xxi**2+yyi**2+zzi**2 ! 距離の自乗

             rik=dsqrt(rik2)
             rrik=rik  ! Åに直す
                                       ! jとkとの距離
             xxj=(x(j)+dble(cx)*clx)-(x(k)+dble(cx1)*clx)
             yyj=(y(j)+dble(cy)*cly)-(y(k)+dble(cy1)*cly)
             zzj=(z(j)+dble(cz)*clz)-(z(k)+dble(cz1)*clz)

             rjk2=xxj**2+yyj**2+zzj**2 ! 距離の自乗

             rjk=dsqrt(rjk2)
             rrjk=rjk    ! Åに直す

             if( (rik2.le.rc2).or.(rjk2.le.rc2)) then
C              print *,rrik,rrjk
              call fij (rrik,fsik,dfsik)
              call fij (rrjk,fsjk,dfsjk)

              cosijk=(rr**2+rrik**2-rrjk**2)
     &                  /(2.0d0*rr*rrik)
              cosjik=(rr**2+rrjk**2-rrik**2)
     &                  /(2.0d0*rr*rrjk)



              gi=a0*(1.0d0+c2/d2-c2/(d2+(1.0d0+cosijk)**2))
              gj=a0*(1.0d0+c2/d2-c2/(d2+(1.0d0+cosjik)**2))

              dgi=a0*c2*2.0d0*(1+cosijk)
     &                 /(d2+(1.0d0+cosijk)**2)**2
              dgj=a0*c2*2.0d0*(1+cosjik)
     &                 /(d2+(1.0d0+cosjik)**2)**2

              sumbij(i,j,cx,cy,cz)=sumbij(i,j,cx,cy,cz)+gi*fsik
              sumbji(i,j,cx,cy,cz)=sumbji(i,j,cx,cy,cz)+gj*fsjk
              dth1=(rr**2-rrik**2+rrjk**2)
     &                      /(2.0d0*(rr**2)*rrik)
              dth2=(rrik**2-rr**2+rrjk**2)
     &                      /(2.0d0*rr*(rrik**2))

              dth3=(rr**2-rrjk**2+rrik**2)
     &                      /(2.0d0*(rr**2)*rrjk)
              dth4=-rrik/(rr*rrjk)

              dbijx(i,j,cx,cy,cz)=dbijx(i,j,cx,cy,cz)
     &     +(dgi*fsik*dth1*xxr/r
     &     +(dgi*fsik*dth2+gi*dfsik)*xxi/rik)
              dbijy(i,j,cx,cy,cz)=dbijy(i,j,cx,cy,cz)
     &     +(dgi*fsik*dth1*yyr/r
     &     +(dgi*fsik*dth2+gi*dfsik)*yyi/rik)
              dbijz(i,j,cx,cy,cz)=dbijz(i,j,cx,cy,cz)
     &     +(dgi*fsik*dth1*zzr/r
     &     +(dgi*fsik*dth2+gi*dfsik)*zzi/rik)

              dbjix(i,j,cx,cy,cz)=dbjix(i,j,cx,cy,cz)
     &     +dgj*fsjk*(dth3*xxr/r+dth4*xxi/rik)
              dbjiy(i,j,cx,cy,cz)=dbjiy(i,j,cx,cy,cz)
     &     +dgj*fsjk*(dth3*yyr/r+dth4*yyi/rik)
              dbjiz(i,j,cx,cy,cz)=dbjiz(i,j,cx,cy,cz)
     &     +dgj*fsjk*(dth3*zzr/r+dth4*zzi/rik)

              xik=nc(k)-fsik    !Xikの計算
              xjk=nc(k)-fsjk    !Xjkの計算

              call flij (xik,flik,dflik)
              call flij (xjk,fljk,dfljk)
C              print *,i,j,cx,cy,cz,k,cx1,cy1,cz1,fsik,flik,fsjk,fljk
              sumni=sumni+fsik*flik
              sumnj=sumnj+fsjk*fljk

           endif
7000       continue
           enddo
            enddo
            enddo
 5300     continue

          ncon(i,j)=1.0d0+sumni+sumnj

c     F(i,j,k)の計算
          nci=nc(i)
          ncj=nc(j)
          nconj=ncon(i,j)

C          print *,i,j,cx,cy,cz,nci,ncj,nconj

          call sub1 (n,nf1,nf2,lf,nci,ncj,nconj,hx,hy,
     &           ff,c,la,lb,lm,alpha)
          sp(i,j,cx,cy,cz)=ff
C          print *,ff
          call sub1 (n,nf1,nf2,dilf,nci,ncj,nconj,hx,hy,
     &                           ff,c,la,lb,lm,alpha)
          disp(i,j,cx,cy,cz)=ff

          call sub1 (n,nf1,nf2,djlf,nci,ncj,nconj,hx,hy,
     &                           ff,c,la,lb,lm,alpha)
          djsp(i,j,cx,cy,cz)=ff

c     Bij,Bjiを計算
          bij=(1.0d0+sumbij(i,j,cx,cy,cz))**(-th)
          bji=(1.0d0+sumbji(i,j,cx,cy,cz))**(-th)

          bbij=(bij+bji)/2.0d0+sp(i,j,cx,cy,cz)

c     iがjから感じるポテンシャルエネルギー
          enpoteach=(vr-bbij*va)

          enpot=enpot+enpoteach/2.0d0 ! ポテンシャルエネルギーの累計
C          print *,enpot,enpoteach
         endif
        endif
       enddo
       enddo
       enddo
 5200 continue




      do 4200 ln=1,lvn          ! リストベクトルによるループ

       i=lvi(ln)
       j=lvj(ln)
C      print *,i,j
       do cx=-irangex,irangex
       do cy=-irangey,irangey
       do cz=-irangez,irangez
                                       ! iとjとの距離
        xxr=x(i)-(x(j)+dble(cx)*clx)
        yyr=y(i)-(y(j)+dble(cy)*cly)
        zzr=z(i)-(z(j)+dble(cz)*clz)

        if (((abs(xxr).le.rc).and.(abs(yyr).le.rc))
     &              .and.(abs(zzr).le.rc)) then
         r2=xxr**2+yyr**2+zzr**2 ! 距離の自乗

         if (r2.lt.rc2) then
          r=dsqrt(r2)
          rr=r         ! Åに直す

          dis(ln)=rr
          sumdis=sumdis+rr
          call fij (rr,fsij,dfsij)


c     角度を計算


          sumdnix=0.0d0
          sumdniy=0.0d0
          sumdniz=0.0d0


          do 4300 k=1,n

            do cx1=-irangex,irangex
            do cy1=-irangey,irangey
            do cz1=-irangez,irangez
            if((k.eq.i).and.(cx1.eq.0)
     &                 .and.(cy1.eq.0).and.(cz1.eq.0))go to 7001
            if((k.eq.j).and.(cx1.eq.cx)
     &                 .and.(cy1.eq.cy).and.(cz1.eq.cz))go to 7001                                        ! iとkとの距離
             xxi=x(i)-(x(k)+dble(cx1)*clx)
             yyi=y(i)-(y(k)+dble(cy1)*cly)
             zzi=z(i)-(z(k)+dble(cz1)*clz)

             rik2=xxi**2+yyi**2+zzi**2 ! 距離の自乗

             rik=dsqrt(rik2)
             rrik=rik  ! Åに直す
                                        ! jとkとの距離
             xxj=(x(j)+dble(cx)*clx)-(x(k)+dble(cx1)*clx)
             yyj=(y(j)+dble(cy)*cly)-(y(k)+dble(cy1)*cly)
             zzj=(z(j)+dble(cz)*clz)-(z(k)+dble(cz1)*clz)

             rjk2=xxj**2+yyj**2+zzj**2 ! 距離の自乗

             rjk=dsqrt(rjk2)
             rrjk=rjk ! Åに直す



              if(rjk2.le.rc2) then

               call fij (rrik,fsik,dfsik)
               call fij (rrjk,fsjk,dfsjk)

               sumdnix=sumdnix+dfsik*disp(i,j,cx,cy,cz)*xxi/rik
               sumdniy=sumdniy+dfsik*disp(i,j,cx,cy,cz)*yyi/rik
               sumdniz=sumdniz+dfsik*disp(i,j,cx,cy,cz)*zzi/rik
c-----------------------------------

              va=fsjk*de*s
     &           *exp(-dsqrt(2.0d0/s)*be*(rrjk-re))/(s-1.0d0)



              cosjki=(rr**2+rrjk**2-rrik**2)
     &                  /(2.0d0*rr*rrjk)
              coskji=(rrik**2+rrjk**2-rr**2)
     &                  /(2.0d0*rrik*rrjk)


              gj=a0*(1.0d0+c2/d2-c2/(d2+(1.0d0+cosjki)**2))
              gk=a0*(1.0d0+c2/d2-c2/(d2+(1.0d0+coskji)**2))

              dgj=a0*c2*2.0d0*(1+cosjki)
     &                 /(d2+(1.0d0+cosjki)**2)**2
              dgk=a0*c2*2.0d0*(1+coskji)
     &                 /(d2+(1.0d0+coskji)**2)**2



              dth5=-rr/(rrik*rrjk)
              dth6=(rrik**2-rrjk**2+rr**2)
     &                      /(2.0d0*(rrik**2)*rrjk)

              dbjkrij=(dgj*fsij*dth3+gj*dfsij)
              dbjkrik=dgj*fsij*dth4
              dbkjrij=dgk*fsik*dth5
              dbkjrik=(dgk*fsik*dth6+gk*dfsik)

              kx=cx1-cx
              ky=cy1-cy
              kz=cz1-cz

              if ((abs(kx).gt.irx).or.
     &            (abs(ky).gt.iry).or.(abs(kz).gt.irz)) go to 7002

              fx(i)=fx(i)+(-va/4.0d0*th*(
     &          ((1.0d0+sumbij(j,k,kx,ky,kz))**(-th-1.0d0))
     &          *(dbjkrij*xxr/r+dbjkrik*xxi/rik)
     &          +((1.0d0+sumbji(j,k,kx,ky,kz))**(-th-1.0d0))
     &          *(dbkjrij*xxr/r+dbkjrik*xxi/rik))
     &          +va/2.0d0*(dfsij*disp(j,k,kx,ky,kz)*xxr/r
     &          +dfsik*djsp(j,k,kx,ky,kz)*xxi/rik)
     &          )*1.60219d-9

              fy(i)=fy(i)+(-va/4.0d0*th*(
     &           ((1.0d0+sumbij(j,k,kx,ky,kz))**(-th-1.0d0))
     &           *(dbjkrij*yyr/r+dbjkrik*yyi/rik)
     &           +((1.0d0+sumbji(j,k,kx,ky,kz))**(-th-1.0d0))
     &           *(dbkjrij*yyr/r+dbkjrik*yyi/rik))
     &           +va/2.0d0*(dfsij*disp(j,k,kx,ky,kz)*yyr/r
     &           +dfsik*djsp(j,k,kx,ky,kz)*yyi/rik)
     &           )*1.60219d-9

              fz(i)=fz(i)+(-va/4.0d0*th*(
     &           ((1.0d0+sumbij(j,k,kx,ky,kz))**(-th-1.0d0))
     &           *(dbjkrij*zzr/r+dbjkrik*zzi/rik)
     &           +((1.0d0+sumbji(j,k,kx,ky,kz))**(-th-1.0d0))
     &           *(dbkjrij*zzr/r+dbkjrik*zzi/rik))
     &           +va/2.0d0*(dfsij*disp(j,k,kx,ky,kz)*zzr/r
     &           +dfsik*djsp(j,k,kx,ky,kz)*zzi/rik)
     &           )*1.60219d-9


7002         continue
             endif
7001         enddo
            enddo
            enddo
 4300     continue
C          if (i.eq.1)then
C          write (6,902),i,j,fx(1),fy(1),fz(1)
C          endif




c     Bij,Bjiとそれぞれの微分を計算
          bij=(1.0d0+sumbij(i,j,cx,cy,cz))**(-th)
          bji=(1.0d0+sumbji(i,j,cx,cy,cz))**(-th)

          bbij=(bij+bji)/2.0d0+sp(i,j,cx,cy,cz)

c     Ｖr(r）とＶa(r）を計算
          va=fsij*de*s*exp(-dsqrt(2.0d0/s)*be*(rr-re))/(s-1.0d0)
          dvr=de*exp(-dsqrt(2.0d0*s)*be*(rr-re))
     &           *(dfsij-dsqrt(2.0d0*s)*be*fsij)/(s-1.0d0)
          dva=de*s*exp(-dsqrt(2.0d0/s)*be*(rr-re))
     &           *(dfsij-dsqrt(2.0d0/s)*be*fsij)/(s-1.0d0)

c     Fの計算

                                ! iが受ける力を累計
          fx(i)=fx(i)-va/2.0d0*th*(
     &       ((1.0d0+sumbij(i,j,cx,cy,cz))**(-th-1.0d0))
     &       *dbijx(i,j,cx,cy,cz)
     &      +((1.0d0+sumbji(i,j,cx,cy,cz))**(-th-1.0d0))
     &       *dbjix(i,j,cx,cy,cz)
     &       )*1.60219d-9
          fy(i)=fy(i)-va/2.0d0*th*(
     &       ((1.0d0+sumbij(i,j,cx,cy,cz))**(-th-1.0d0))
     &        *dbijy(i,j,cx,cy,cz)
     &      +((1.0d0+sumbji(i,j,cx,cy,cz))**(-th-1.0d0))
     &       *dbjiy(i,j,cx,cy,cz)
     &       )*1.60219d-9
          fz(i)=fz(i)-va/2.0d0*th*(
     &       ((1.0d0+sumbij(i,j,cx,cy,cz))**(-th-1.0d0))
     &        *dbijz(i,j,cx,cy,cz)
     &      +((1.0d0+sumbji(i,j,cx,cy,cz))**(-th-1.0d0))
     &        *dbjiz(i,j,cx,cy,cz)
     &       )*1.60219d-9

C          if (i.eq.1)then
C          write (6,903),i,j,-va/2.0d0*th*(
C     &       ((1.0d0+sumbij(i,j,cx,cy,cz))**(-th-1.0d0))
C     &       *dbijx(i,j,cx,cy,cz)
C     &      +((1.0d0+sumbji(i,j,cx,cy,cz))**(-th-1.0d0))
C     &       *dbjix(i,j,cx,cy,cz))
C     &       ,-va/2.0d0*th*(
C     &       ((1.0d0+sumbij(i,j,cx,cy,cz))**(-th-1.0d0))
C     &        *dbijy(i,j,cx,cy,cz)
C     &      +((1.0d0+sumbji(i,j,cx,cy,cz))**(-th-1.0d0))
C     &       *dbjiy(i,j,cx,cy,cz))
C     &        ,-va/2.0d0*th*(
C     &       ((1.0d0+sumbij(i,j,cx,cy,cz))**(-th-1.0d0))
C     &        *dbijz(i,j,cx,cy,cz)
C     &      +((1.0d0+sumbji(i,j,cx,cy,cz))**(-th-1.0d0))
C     &        *dbjiz(i,j,cx,cy,cz))
C            endif

          frij=-(dvr-bbij*dva)*1.60219d-9


                                ! iが受ける力を累計
          fx(i)=fx(i)+frij/r*xxr
     &         +va*(dncx(i)*disp(i,j,cx,cy,cz)
     &             +dfsij*xxr/r*djsp(i,j,cx,cy,cz))
     &       *1.60219d-9
          fy(i)=fy(i)+frij/r*yyr
     &         +va*(dncy(i)*disp(i,j,cx,cy,cz)
     &             +dfsij*yyr/r*djsp(i,j,cx,cy,cz))
     &       *1.60219d-9
          fz(i)=fz(i)+frij/r*zzr
     &         +va*(dncz(i)*disp(i,j,cx,cy,cz)
     &             +dfsij*zzr/r*djsp(i,j,cx,cy,cz))
     &       *1.60219d-9
C          if (i.eq.1)then
C          write (6,905)i,j,frij/r*xxr,frij/r*yyr,frij/r*zzr
C          write (6,905)i,j,-va*dncx(i)*disp(i,j,cx,cy,cz)
C     &         ,-va*dncy(i)*disp(i,j,cx,cy,cz)
C     &         ,-va*dncz(i)*disp(i,j,cx,cy,cz)
C          write (6,905)i,j,-va*xxr/r*dfsij*djsp(i,j,cx,cy,cz)
C     &         ,-va*dfsij*yyr/r*djsp(i,j,cx,cy,cz)
C     &         ,-va*dfsij*zzr/r*djsp(i,j,cx,cy,cz)
C  905 FORMAT(2x,2HC ,2I3,3D24.16)
C
C          write (6,904),i,j,frij/r*xxr
C     &         -va*(dncx(i)*disp(i,j,cx,cy,cz)+dfsij*djsp(i,j,cx,cy,cz))
C     &         ,frij/r*yyr
C     &         -va*(dncy(i)*disp(i,j,cx,cy,cz)+dfsij*djsp(i,j,cx,cy,cz))
C     &         ,frij/r*zzr
C     &         -va*(dncz(i)*disp(i,j,cx,cy,cz)+dfsij*djsp(i,j,cx,cy,cz))
C  902 FORMAT(2x,2HA ,2I3,3D24.16)
C  903 FORMAT(2x,2HB ,2I3,3D24.16)
C  904 FORMAT(2x,2HC ,2I3,3D24.16)
C          endif


         endif
        endif
       enddo
       enddo
       enddo
 4200 continue

C      do i=1,n
C      write (6,901) I,FX(I),FY(I),FZ(I)
C      enddo
C  901 FORMAT(2x,8Hforce   ,I6,3D24.16)


      RETURN
      END

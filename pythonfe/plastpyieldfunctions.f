! calculates value of different yield functions
      subroutine plastpyieldfunctions(idfy, d, Man, fy, dfydM, dfydkap, 
     & kap, dfydmandman, devproj, dmandfp, dfydfp, hardforce, dfydf,
     & dMdf, fp)

      implicit none


      integer idfy, i, j, k, l, o

      real (kind=8) d(*), Man(3,3), fy, dfydM(3,3), kap, dfydkap
      real (kind=8) devMan(3,3), trMan, id(3,3), HH, man0, ndevman
      real (kind=8) frac23,dfydmandman(3,3,3,3), devproj(3,3,3,3)
      real (kind=8) dmandfp(3,3,3,3), dfydfp(3,3), hardforce
      real (kind=8) dfydf(3,3), dMdf(3,3,3,3), frac32, alpha
      real (kind=8) SVT,SVC,RR,RT,Heavyside,fc,ft,f1,XX,TT,T0
      real (kind=8) dXdkap, dTdkap, dfcdkap, dftdkap, df1dkap
      real (kind=8) df1dM(3,3), dfcdM(3,3), dftdM(3,3)
      real (kind=8) dfcdMdM(3,3,3,3), dftdMdM(3,3,3,3)
      real (kind=8) ff,gg,hhh,ll,mm,nn,oo,pp,qq
      real (kind=8) N1(3),N2(3),N3(3),NI1(3),NI2(3),NI3(3)
      real (kind=8) M11,M12,M13,M21,M22,M23,M31,M32,M33, fp(3,3)
      real (kind=8) NI1S(3),NI2S(3),NI3S(3),SNI1(3),SNI2(3),SNI3(3)
      real (kind=8) dfydNI1(3),dfydNI2(3),dfydNI3(3)
      real (kind=8) dfydMdNI1(3,3,3),dfydMdNI2(3,3,3),dfydMdNI3(3,3,3)
      real (kind=8) dNI1dfp(3,3,3),dNI2dfp(3,3,3),dNI3dfp(3,3,3)
      real (kind=8) dfydMdfp(3,3,3,3),NNI1,NNI2,NNI3

      ! remember the hardforce is the derivative of the yield function
      ! with respect to the q with negative sign

      ! the q is the partial derivative of the free Helmholtz energy 
      ! with respect to the hardening variable

      data id/ 1.d0 , 0.d0 , 0.d0 ,
     *         0.d0 , 1.d0 , 0.d0 ,
     *         0.d0 , 0.d0 , 1.d0 /

      do i = 1,3
            do j = 1,3
                  do k = 1,3
                        do l = 1,3
      devproj(i,j,k,l) = id(i,k) * id(j,l) - id(i,j) * id(k,l) / 3.d0 
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      if(idfy==1) then

      ! here boring von Mises with linear Hardening
      trMan = Man(1,1) + Man(2,2) + Man(3,3)
      devMan(:,:) = Man(:,:) - trMan / 3.d0 * id(:,:)
      frac23 = dsqrt(2.d0/3.d0)

      ndevman = 0.d0
      do i = 1,3
            do j = 1,3
            ndevman = ndevman + devMan(i,j) * devMan(i,j)
            enddo ! j
      enddo ! i
      ndevman = dsqrt(ndevman)

      man0 = d(1)
      HH = d(2)
      fy = ndevman - frac23 * (man0 + HH * kap)

      dfydkap = -frac23 * HH

      if(ndevman.gt.1.d-10) then
            dfydM = devMan / ndevman ! in theory multiplied with the
            do i = 1,3
                  do j = 1,3
                        do k = 1,3
                              do l = 1,3
      dfydmandman(i,j,k,l) = (devproj(i,j,k,l) - dfydM(i,j) 
     &                        * dfydM(k,l)) / ndevman
                              enddo ! l
                        enddo ! k
                  enddo ! j
            enddo ! i
      else ! ndevman
            dfydM(:,:) = 0.d0
            dfydmandman(:,:,:,:) = 0.d0
      endif ! ndevman

      ! deviatoric projection tensor
      ! Nevertheless, that would let the dev stay dev
      

      ! here the hardening term in free Helmholtz is
      ! 0.5d0 * HH * kap^2
      ! => q = HH * kap
      ! dfy/dq = -2/3
      ! hardforce = 2/3
      hardforce = frac23

      elseif(idfy==2) then

      man0 = d(1)
      HH = d(2)
      alpha = d(3)

      ! here Drucker-Prager with linear hardening
      trMan = Man(1,1) + Man(2,2) + Man(3,3)
      devMan(:,:) = Man(:,:) - trMan / 3.d0 * id(:,:)
      frac32 = dsqrt(3.d0/2.d0)

      ndevman = 0.d0
      do i = 1,3
            do j = 1,3
            ndevman = ndevman + devMan(i,j) * devMan(i,j)
            enddo ! j
      enddo ! i
      ndevman = dsqrt(ndevman)

      fy = frac32 * ndevman + alpha * trMan - (man0 + HH * kap)

      dfydkap = - HH

      dfydM(:,:) = alpha * id(:,:)
      dfydmandman = 0.d0

      if(ndevman.gt.1.d-10) then

      do i = 1,3
            do j = 1,3
            dfydM(i,j) = dfydM(i,j) + frac32 * devMan(i,j) /ndevman
                  do k = 1,3
                        do l = 1,3
      dfydmandman(i,j,k,l) = frac32 * (devproj(i,j,k,l) - devMan(i,j) *
     & devMan(k,l) / (ndevman**2.d0))/ndevman
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      endif ! ndevman

      hardforce = 1.d0

      elseif(idfy==3) then
      
      ! Drucker Prager with caps

      man0 = d(1)
      HH = d(2)
      alpha = d(3)
      SVT = d(4)
      RT = d(5)
      T0 = d(6)
      SVC = d(7)
      RR = d(8)

      frac32 = 3.d0/2.d0

      trMan = Man(1,1) + Man(2,2) + Man(3,3)
      devMan(:,:) = Man(:,:) - trMan / 3.d0 * id(:,:)

      ! get all the functions 
      ! linear hardening
      f1 = man0 - alpha * trMan + HH * kap

      ! compression cap
      XX = RR * (man0 - alpha * SVC + HH * kap)
      fc = 1.d0 - Heavyside(SVC - trMan) * ((trMan - SVC)**2.d0) / 
     & (XX**2.d0)

      ! tension cap
      TT = T0 + RT * HH * kap
      ft = 1.d0 - Heavyside(trMan - SVT) * ((trMan - SVT)**2.d0) / 
     & ((TT - SVT)**2.d0)

      ndevman = 0.d0
      do i = 1,3
            do j = 1,3
            ndevman = ndevman + devMan(i,j) * devMan(i,j)
            enddo ! j
      enddo ! i
      ndevman = dsqrt(ndevman)

      fy = frac32 * ndevman**2.d0 - f1**2.d0 * fc * ft

      ! derivatives with respect to hardening
      df1dkap = HH
      dXdkap = RR * HH
      dTdkap = RT * HH

      dfcdkap = 2.d0 * Heavyside(SVC - trMan) * ((trMan - SVC)**2.d0) /
     & (XX**3.d0) * dXdkap

      dftdkap = 2.d0 * Heavyside(trMan - SVT) * ((trMan - SVT)**2.d0) /
     & ((TT - SVT)**3.d0) * dTdkap

      dfydkap = -2.d0 * f1 * df1dkap * fc * ft - f1**2.d0 * dfcdkap * ft
     & - f1**2.d0 * fc * dftdkap

      ! derivatives with respect to stress
      df1dM = -alpha * id
      dfcdM = -2.d0 * Heavyside(SVC - trMan) *(trMan - SVC) / (XX**2.d0)
     & * id
      dftdM = -2.d0 * Heavyside(trMan - SVT) *(trMan - SVT) / 
     & ((TT - SVT)**2.d0) * id

      dfydM = 3.d0 * devMan - 2.d0 * f1 * df1dM * fc * ft - f1**2.d0 *
     & dfcdM * ft - f1**2.d0 * fc * dftdM

      ! second derivatives
      do i = 1,3
            do j = 1,3
                  do k = 1,3
                        do l = 1,3
      dfcdMdM(i,j,k,l) = -2.d0 * Heavyside(SVC - trMan) / (XX**2.d0) *
     & id(i,j) * id(k,l)

      dftdMdM(i,j,k,l) = -2.d0 * Heavyside(trMan - SVT) / 
     & ((TT - SVT)**2.d0) * id(i,j) * id(k,l)
                        enddo ! L
                  enddo ! k
            enddo ! j
      enddo ! i

      do i = 1,3
            do j = 1,3
                  do k = 1,3
                        do l = 1,3
      dfydmandman(i,j,k,l) = 3.d0 * devproj(i,j,k,l) 
     & -2.d0 * df1dM(i,j) * df1dM(k,l) * fc * ft
     & -2.d0 * f1 * df1dM(i,j) * dfcdM(k,l) * ft
     & -2.d0 * f1 * df1dM(i,j) * fc * dftdM(k,l)
     & -2.d0 * f1 * df1dM(k,l) * dfcdM(i,j) * ft
     & - f1**2.d0 * dfcdMdM(i,j,k,l) * ft
     & - f1**2.d0 * dfcdM(i,j) * dftdM(k,l)
     & -2.d0 * f1 * df1dM(k,l) * fc * dftdM(i,j)
     & - f1**2.d0 * dfcdM(k,l) * dftdM(i,j)
     & - f1**2.d0 * fc * dftdMdM(i,j,k,l)
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      hardforce = 1.d0

      elseif(idfy==5) then

      ! quadratic Hill without rotatable vectors

      ff = d(1)
      gg = d(2)
      hhh = d(3)
      ll = d(4)
      mm = d(5)
      nn = d(6)
      oo = d(7)
      pp = d(8)
      qq = d(9)
      hh = d(10)

      fy = ff * (Man(2,2) - Man(3,3))**(2.d0) 
     &   + gg * (Man(3,3) - Man(1,1))**(2.d0)
     &   + hhh *(Man(1,1) - Man(2,2))**(2.d0)
     &   + ll * Man(2,3)**(2.d0)
     &   + mm * Man(3,1)**(2.d0)
     &   + nn * Man(1,2)**(2.d0)
     &   + oo * Man(3,2)**(2.d0)
     &   + pp * Man(1,3)**(2.d0)
     &   + qq * Man(2,1)**(2.d0)
     &   - 1.d0 - hh * kap

      dfydkap = -hh
      hardforce = 1.d0

      dfydM(1,1) = gg * 2.d0 * (Man(1,1) - Man(3,3)) 
     &           + hhh * 2.d0 * (Man(1,1) - Man(2,2))
      dfydM(2,2) = ff * 2.d0 * (Man(2,2) - Man(3,3))
     &           + hhh * 2.d0 * (Man(2,2) - Man(1,1))
      dfydM(3,3) = ff * 2.d0 * (Man(3,3) - Man(2,2))
     &           + gg * 2.d0 * (Man(3,3) - Man(1,1))
      dfydM(1,2) = 2.d0 * nn * Man(1,2)
      dfydM(2,3) = 2.d0 * ll * Man(2,3)
      dfydM(1,3) = 2.d0 * pp * Man(1,3)
      dfydM(2,1) = 2.d0 * qq * Man(2,1)
      dfydM(3,2) = 2.d0 * oo * Man(3,2)
      dfydM(3,1) = 2.d0 * mm * Man(3,1)

      dfydmandman = 0.d0  ! since quite some terms are zero
      dfydmandman(1,1,1,1) = 2.d0 * gg + 2.d0 * hhh
      dfydmandman(1,1,2,2) = - 2.d0 * hhh
      dfydmandman(1,1,3,3) = - 2.d0 * gg
      dfydmandman(2,2,1,1) = -2.d0 * hhh
      dfydmandman(2,2,2,2) = 2.d0 * ff + 2.d0 * hhh
      dfydmandman(2,2,3,3) = -2.d0 * ff
      dfydmandman(3,3,1,1) = -2.d0 * gg 
      dfydmandman(3,3,2,2) = - 2.d0 * ff
      dfydmandman(3,3,3,3) = 2.d0 * ff + 2.d0 * gg
      dfydmandman(1,2,1,2) = 2.d0 * nn
      dfydmandman(2,3,2,3) = 2.d0 * ll
      dfydmandman(1,3,1,3) = 2.d0 * pp
      dfydmandman(2,1,2,1) = 2.d0 * qq
      dfydmandman(3,2,3,2) = 2.d0 * oo
      dfydmandman(3,1,3,1) = 2.d0 * mm

      elseif(idfy==6) then
      
      ! quadratic Hill in finite strains
      
      ff = d(1)
      gg = d(2)
      hhh = d(3)
      ll = d(4)
      mm = d(5)
      nn = d(6)
      oo = d(7)
      pp = d(8)
      qq = d(9)
      hh = d(10)
      N1 = d(11:13)
      N2 = d(14:16)
      N3 = d(17:19)

      do i = 1,3
        NI1(i) = 0.d0
        NI2(i) = 0.d0
        NI3(i) = 0.d0
        do j = 1,3
          NI1(i) = NI1(i) + fp(i,j) * N1(j)
          NI2(i) = NI2(i) + fp(i,j) * N2(j)
          NI3(i) = NI3(i) + fp(i,j) * N3(j)
        enddo ! j
      enddo ! i

      NNI1 = dsqrt(NI1(1)**2.d0+NI1(2)**2.d0+NI1(3)**2.d0)
      NNI2 = dsqrt(NI2(1)**2.d0+NI2(2)**2.d0+NI2(3)**2.d0)
      NNI3 = dsqrt(NI3(1)**2.d0+NI3(2)**2.d0+NI3(3)**2.d0)

      ! prevent scaling of some things
      NI1(1:3) = NI1(1:3) / NNI1
      NI2(1:3) = NI2(1:3) / NNI2
      NI3(1:3) = NI3(1:3) / NNI3

      M11 = 0.d0
      M12 = 0.d0
      M13 = 0.d0
      M21 = 0.d0
      M22 = 0.d0
      M23 = 0.d0
      M31 = 0.d0
      M32 = 0.d0
      M33 = 0.d0
      do i = 1,3
            NI1S(i) = 0.d0
            NI2S(i) = 0.d0
            NI3S(i) = 0.d0
            SNI1(i) = 0.d0
            SNI2(i) = 0.d0
            SNI3(i) = 0.d0
            do j = 1,3
      ! scalar variables 
      M11 = M11 + NI1(i) * NI1(j) * Man(i,j)
      M12 = M12 + NI1(i) * NI2(j) * Man(i,j)
      M13 = M13 + NI1(i) * NI3(j) * Man(i,j)
      M21 = M21 + NI2(i) * NI1(j) * Man(i,j)
      M22 = M22 + NI2(i) * NI2(j) * Man(i,j)
      M23 = M23 + NI2(i) * NI3(j) * Man(i,j)
      M31 = M31 + NI3(i) * NI1(j) * Man(i,j)
      M32 = M32 + NI3(i) * NI2(j) * Man(i,j)
      M33 = M33 + NI3(i) * NI3(j) * Man(i,j)

      ! vector variables
      NI1S(i) = NI1S(i) + NI1(j) * Man(j,i)
      NI2S(i) = NI2S(i) + NI2(j) * Man(j,i)
      NI3S(i) = NI3S(i) + NI3(j) * Man(j,i)
      SNI1(i) = SNI1(i) + Man(i,j) * NI1(j)
      SNI2(i) = SNI2(i) + Man(i,j) * NI2(j)
      SNI3(i) = SNI3(i) + Man(i,j) * NI3(j)

                  ! third order tensor variables
                  do k = 1,3
                        ! old, used for without normalizing
                        !dNI1dfp(i,j,k) = id(i,j) * N1(k)
                        !dNI2dfp(i,j,k) = id(i,j) * N2(k)
                        !dNI3dfp(i,j,k) = id(i,j) * N3(k)
      dNI1dfp(i,j,k) = 1/NNI1 * (id(i,j) - NI1(i) * NI1(j)) * N1(k)
      dNI2dfp(i,j,k) = 1/NNI2 * (id(i,j) - NI2(i) * NI2(j)) * N2(k)
      dNI3dfp(i,j,k) = 1/NNI3 * (id(i,j) - NI3(i) * NI3(j)) * N3(k)
                  enddo ! k
            enddo ! j
      enddo ! i

      fy = ff * (M22 - M33)**(2.d0) 
     &   + gg * (M33 - M11)**(2.d0)
     &   + hhh *(M11 - M22)**(2.d0)
     &   + ll * M23**(2.d0)
     &   + mm * M31**(2.d0)
     &   + nn * M12**(2.d0)
     &   + oo * M32**(2.d0)
     &   + pp * M13**(2.d0)
     &   + qq * M21**(2.d0)
     &   - 1.d0 - hh * kap
     
      ! derivatives of yield function
      do i = 1,3
            dfydNI1(i) = 
     &-2.d0 * gg * (M33 - M11) * (NI1S(i) + SNI1(i))
     &+2.d0* hhh * (M11 - M22) * (NI1S(i) + SNI1(i))
     &+2.d0 * mm * M31 * NI3S(i)
     &+2.d0 * nn * M12 * SNI2(i)
     &+2.d0 * pp * M13 * SNI3(i)
     &+2.d0 * qq * M21 * NI2S(i)

            dfydNI2(i) = 
     & 2.d0 * ff * (M22 - M33) * (NI2S(i) + SNI2(i))
     &-2.d0* hhh * (M11 - M22) * (NI2S(i) + SNI2(i))
     &+2.d0 * ll * M23 * SNI3(i)
     &+2.d0 * nn * M12 * NI1S(i)
     &+2.d0 * oo * M32 * NI3S(i) 
     &+2.d0 * qq * M21 * SNI1(i)

            dfydNI3(i) = 
     &-2.d0 * ff * (M22 - M33) * (NI3S(i) + SNI3(i))
     &+2.d0 * gg * (M33 - M11) * (NI3S(i) + SNI3(i))
     &+2.d0 * ll * M23 * NI2S(i)
     &+2.d0 * mm * M31 * SNI1(i)
     &+2.d0 * oo * M32 * SNI2(i)
     &+2.d0 * pp * M13 * NI1S(i)

            do j = 1,3
      dfydM(i,j) = 
     & 2.d0 * ff * (M22 - M33) * (NI2(i) * NI2(j) - NI3(i) * NI3(j))
     &+2.d0 * gg * (M33 - M11) * (NI3(i) * NI3(j) - NI1(i) * NI1(j))
     &+2.d0* hhh * (M11 - M22) * (NI1(i) * NI1(j) - NI2(i) * NI2(j))
     &+2.d0 * ll * M23 * (NI2(i) * NI3(j))
     &+2.d0 * mm * M31 * (NI3(i) * NI1(j))
     &+2.d0 * nn * M12 * (NI1(i) * NI2(j))
     &+2.d0 * oo * M32 * (NI3(i) * NI2(j))
     &+2.d0 * pp * M13 * (NI1(i) * NI3(j))
     &+2.d0 * qq * M21 * (NI2(i) * NI1(j))

                  do k = 1,3

      dfydMdNI1(i,j,k) = 
     &-2.d0 * gg * (NI3(i) * NI3(j) - NI1(i) * NI1(j))
     &           * (NI1S(k) + SNI1(k))
     &-2.d0 * gg * (M33 - M11) * (id(i,k) * NI1(j) + NI1(i) * id(j,k))
     &+2.d0* hhh * (NI1(i) * NI1(j) - NI2(i) * NI2(j))
     &           * (NI1S(k) + SNI1(k))
     &+2.d0* hhh * (M11 - M22) * (id(i,k) * NI1(j) + NI1(i) * id(j,k))
     &+2.d0 * mm * NI3(i) * NI1(j) * NI3S(k)
     &+2.d0 * mm * M31 * NI3(i) * id(j,k)
     &+2.d0 * nn * NI1(i) * NI2(j) * SNI2(k)
     &+2.d0 * nn * M12 * id(i,k) * NI2(j)
     &+2.d0 * pp * NI1(i) * NI3(j) * SNI3(k)
     &+2.d0 * pp * M13 * id(i,k) * NI3(j)
     &+2.d0 * qq * NI2(i) * NI1(j) * NI2S(k) 
     &+2.d0 * qq * M21 * NI2(i) * id(j,k)

      dfydMdNI2(i,j,k) = 
     & 2.d0 * ff * (NI2(i) * NI2(j) - NI3(i) * NI3(j))
     &           * (SNI2(k) + NI2S(k))
     &+2.d0 * ff * (M22 - M33) * (id(i,k) * NI2(j) + NI2(i) * id(j,k))
     &-2.d0* hhh * (NI1(i) * NI1(j) - NI2(i) * NI2(j))
     &           * (SNI2(k) + NI2S(k))
     &-2.d0* hhh * (M11 - M22) * (id(i,k) * NI2(j) + NI2(i) * id(j,k))
     &+2.d0 * ll * NI2(i) * NI3(j) * SNI3(k) 
     &+2.d0 * ll * M23 * id(i,k) * NI3(j)
     &+2.d0 * nn * NI1(i) * NI2(j) * NI1S(k) 
     &+2.d0 * nn * M12 * NI1(i) * id(j,k)
     &+2.d0 * oo * NI3(i) * NI2(j) * NI3S(k)
     &+2.d0 * oo * M32 * NI3(i) * id(j,k)
     &+2.d0 * qq * NI2(i) * NI1(j) * SNI1(k)
     &+2.d0 * qq * M21 * id(i,k) * NI1(j)

      dfydMdNI3(i,j,k) = 
     &-2.d0 * ff * (NI2(i) * NI2(j) - NI3(i) * NI3(j))
     &           * (SNI3(k) + NI3S(k))
     &-2.d0 * ff * (M22 - M33) * (id(i,k) * NI3(j) + NI3(i) * id(j,k))
     &+2.d0 * gg * (NI3(i) * NI3(j) - NI1(i) * NI1(j))
     &           * (SNI3(k) + NI3S(k))
     &+2.d0 * gg * (M33 - M11) * (id(i,k) * NI3(j) + NI3(i) * id(j,k))
     &+2.d0 * ll * NI2(i) * NI3(j) * NI2S(k) 
     &+2.d0 * ll * M23 * NI2(i) * id(j,k)
     &+2.d0 * mm * NI3(i) * NI1(j) * SNI1(k) 
     &+2.d0 * mm * M31 * id(i,k) * NI1(j)
     &+2.d0 * oo * NI3(i) * NI2(j) * SNI2(k) 
     &+2.d0 * oo * M32 * id(i,k) * NI2(j)
     &+2.d0 * pp * NI1(i) * NI3(j) * NI1S(k) 
     &+2.d0 * pp * M13 * NI1(i) * id(j,k)

                        do l = 1,3
      dfydmandman(i,j,k,l) = 
     & 2.d0 * ff * (NI2(i) * NI2(j) - NI3(i) * NI3(j)) 
     &           * (NI2(k) * NI2(l) - NI3(k) * NI3(l))
     &+2.d0 * gg * (NI3(i) * NI3(j) - NI1(i) * NI1(j)) 
     &           * (NI3(k) * NI3(l) - NI1(k) * NI1(l))
     &+2.d0* hhh * (NI1(i) * NI1(j) - NI2(i) * NI2(j)) 
     &           * (NI1(k) * NI1(l) - NI2(k) * NI2(l))
     &+2.d0 * ll * NI2(i) * NI3(j) * NI2(k) * NI3(l)
     &+2.d0 * mm * NI3(i) * NI1(j) * NI3(k) * NI1(l)
     &+2.d0 * nn * NI1(i) * NI2(j) * NI1(k) * NI2(l)
     &+2.d0 * oo * NI3(i) * NI2(j) * NI3(k) * NI2(l)
     &+2.d0 * pp * NI1(i) * NI3(j) * NI1(k) * NI3(l)
     &+2.d0 * qq * NI2(i) * NI1(j) * NI2(k) * NI1(l)
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      ! already calculate something for dfydfp
      do i = 1,3
            do j = 1,3
                  dfydfp(i,j) = 0.d0
                  do k = 1,3
      dfydfp(i,j) = dfydfp(i,j) 
     & + dfydNI1(k) * dNI1dfp(k,i,j)
     & + dfydNI2(k) * dNI2dfp(k,i,j)
     & + dfydNI3(k) * dNI3dfp(k,i,j)
                  enddo ! k
            enddo ! j
      enddo ! i

      ! add new term
      do i = 1,3
            do j = 1,3
                  do k = 1,3
                        do l = 1,3
                        dfydMdfp(i,j,k,l) = 0.d0
                              do o = 1,3
      dfydMdfp(i,j,k,l) = dfydMdfp(i,j,k,l) 
     & + dfydMdNI1(i,j,o) * dNI1dfp(o,k,l)
     & + dfydMdNI2(i,j,o) * dNI2dfp(o,k,l)
     & + dfydMdNI3(i,j,o) * dNI3dfp(o,k,l)
                              enddo ! o
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      dfydkap = -hh
      hardforce = 1.d0

      else ! idfy
            write(*,*) 'wrong yield function specified'
      endif ! idfy

      do i = 1,3
            do j = 1,3
                  dfydfp(i,j) = 0.d0
                  dfydf(i,j) = 0.d0
                  do k = 1,3
                        do l = 1,3
      dfydfp(i,j) = dfydfp(i,j) + dfydM(k,l) * dmandfp(k,l,i,j)
      dfydf(i,j)  = dfydf(i,j)  + dfydM(k,l) * dMdf(k,l,i,j)
                        enddo ! l
                  enddo ! k
            enddo ! j
      enddo ! i

      end subroutine
      

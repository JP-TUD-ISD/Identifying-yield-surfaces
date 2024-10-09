! material subroutine, which calculates elastic response with Neo-Hooke and value of yield function
      subroutine plastwo(d, f, detf, fyflag, Piola, tang, fyo)

      implicit none

      logical fyflag
      real (kind=8) d(*), f(3,3,2), detf(2), Piola(3,3)
      real (kind=8) tang(3,3,3,3)

      ! stuff, which is needed locally
      integer i, j, k, idfy, dstart, dende, idkap 
      real (kind=8) fpn(3,3), fp(3,3), hardn, hard, fpinv(3,3), jp
      real (kind=8) fe(3,3), je, young, nu, kappa, mu
      real (kind=8) dpsivoldfe(3,3), dpsiisodfe(3,3)
      real (kind=8) dpsivoldfedfe(3,3,3,3), dpsiisodfedfe(3,3,3,3)
      real (kind=8) dpsidfe(3,3), dpsidfedfe(3,3,3,3), Mandel(3,3)
      real (kind=8) dPioladfp(3,3,3,3),dfedfp(3,3,3,3),dfedf(3,3,3,3)
      real (kind=8) dfedfdfp(3,3,3,3,3,3),dMdfp(3,3,3,3),dPdf(3,3,3,3)
      real (kind=8) dMdf(3,3,3,3), fy, dfydM(3,3), dfydkap
      real (kind=8) dfydmandman(3,3,3,3), devproj(3,3,3,3), hardforce
      real (kind=8) dfydf(3,3), totFP(3,3,3,3)
      real (kind=8) dfydfp(3,3)
      real (kind=8) id(3,3), fyo, pi

      data id/ 1.d0 , 0.d0 , 0.d0 ,
     *         0.d0 , 1.d0 , 0.d0 ,
     *         0.d0 , 0.d0 , 1.d0 /

      data pi/3.14159265359d0/

      ! ahhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh
      ! material parameters

      young = d(7)
      nu = d(8)
      idfy = int(d(9))
      if(idfy==1) then
        ! here just state the from to values of d() that need
        ! to be passed to the yieldfunction sub
        dstart = 10 
        dende = 11
      elseif(idfy==2) then
        dstart = 10
        dende = 12
      elseif(idfy==3) then
        dstart = 10
        dende = 17
      elseif(idfy==6) then
        dstart = 19
        dende = 40
            
        ! calculate norm of directional vectors
        d(29:31) = d(29:31) / dsqrt(d(29)**2.d0+d(30)**2.d0+d(31)**2.d0)
        d(32:34) = d(32:34) / dsqrt(d(32)**2.d0+d(33)**2.d0+d(34)**2.d0)
        
        ! convert degrees into radians
        d(38:40) = d(38:40) * pi / 180.d0
        
        ! rotate directional vectors
        call rotateAB(d(38:40), d(29:31), d(32:34))
        
        ! cross product of one and two is third direction
        d(35) = d(30) * d(34) - d(31) * d(33)
        d(36) = d(31) * d(32) - d(29) * d(34)
        d(37) = d(29) * d(33) - d(30) * d(32)
        ! cross product of normed vectors is normed
        
        ! from the yield stresses in the different directions,
        ! derive F,G,H,L,M,N,O,P,Q in yield function
            d(19) = 0.5d0 * (1.d0/(d(14)**2.d0) + 1.d0/(d(18)**2.d0) 
     &      - 1.d0/(d(10)**2.d0)) ! ff
            d(20) = 0.5d0 * (1.d0/(d(18)**2.d0) + 1.d0/(d(10)**2.d0) 
     &      - 1.d0/(d(14)**2.d0)) ! gg
            d(21) = 0.5d0 * (1.d0/(d(10)**2.d0) + 1.d0/(d(14)**2.d0) 
     &      - 1.d0/(d(18)**2.d0)) ! hhh
            d(22) = 1.d0 / (d(15)**2.d0) ! ll
            d(23) = 1.d0 / (d(16)**2.d0) ! mm
            d(24) = 1.d0 / (d(11)**2.d0) ! nn
            d(25) = 1.d0 / (d(17)**2.d0) ! oo
            d(26) = 1.d0 / (d(12)**2.d0) ! pp
            d(27) = 1.d0 / (d(13)**2.d0) ! qq
            
      endif

      idkap = int(d(dende + 1))
      kappa = young / (3.d0 * (1.d0 - 2.d0 * nu))
      mu    = young / (2.d0 * (1.d0 + nu))
      fpn   = id
      hardn = 0.d0

      ! elastic trial
      fp(:,:) = fpn(:,:)
      hard    = hardn

      ! invert fp
      fpinv(:,:) = fp(:,:)
      call invert3new(fpinv(:,:), jp) ! my version, det(F^p) = jp is calculated

      do i = 1,3
        do j = 1,3
          fe(i,j) = 0.d0
          do k = 1,3
            fe(i,j) = fe(i,j) + f(i,k,1) * fpinv(k,j)
          end do ! k
        end do ! j
      end do ! i

      je = detf(1) / jp  ! jp should be 1.d0 ... but be safe

      ! now for the modular approach the material model needs to be called
      ! fe needs to be passed as f, for my purely elastic materials
      ! the transition from dpsidfe to dpsidf takes place in a later sub

      call stressandtangentisotropicneoP(
     &        fe(:,:), je, kappa, mu, 
     &        dpsivoldfe(:,:), dpsiisodfe(:,:), dpsivoldfedfe(:,:,:,:),
     &        dpsiisodfedfe(:,:,:,:) )

      dpsidfe(:,:)        = dpsivoldfe(:,:) + dpsiisodfe(:,:)
      dpsidfedfe(:,:,:,:) = dpsivoldfedfe(:,:,:,:) + 
     &                      dpsiisodfedfe(:,:,:,:)

      call calcPandMandel(dpsidfe(:,:), dpsidfedfe(:,:,:,:), fe(:,:), 
     &                    fp(:,:), fpinv(:,:), Piola(:,:), Mandel(:,:),
     &                    dPioladfp(:,:,:,:),
     &                    dfedfp(:,:,:,:), dfedf(:,:,:,:),
     &                    dfedfdfp(:,:,:,:,:,:), 
     &                    dMdfp(:,:,:,:), dPdf(:,:,:,:), dMdf(:,:,:,:))

      if(fyflag)then
        call plastpyieldfunctions(
     &          idfy, d(dstart:dende), Mandel(:,:), fy, 
     &          dfydM(:,:), dfydkap, hard, dfydmandman(:,:,:,:),
     &          devproj(:,:,:,:), dMdfp(:,:,:,:), dfydfp(:,:),
     &          hardforce, dfydf(:,:), 
     &          dMdf(:,:,:,:), id )

      if(fy.gt.fyo) fyo = fy


      endif ! fyflag

      totFP(:,:,:,:) = 0.d0


      call assemblealgotangPplast(dPdf(:,:,:,:), dPioladfp(:,:,:,:),
     & totFP(:,:,:,:), tang(:,:,:,:))

2000  format(1(1p,1e16.8))

      end subroutine

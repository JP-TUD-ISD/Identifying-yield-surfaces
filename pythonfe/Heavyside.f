! Heavy side function as fortran function
      function Heavyside(x) result(H)

            implicit none

            real*8 x,H
            ! Heavyside function defined like Imad Diss

            ! I guess this case is safe to just capture by an if condition
            if(x.eq.0.d0) then
                  H = 0.5d0
                  return
            endif

            H = 0.5d0 * ( 1.d0 + sign(1.d0, x) )
            ! the sign function of fortran returns the value of the first argument with the sign of the second argument
            ! therefore, if x is positive a 1 = 0.5 * (1 + 1) is returned
            ! for x = 0 a 1/2 = 0.5d0 * (1 + 0) -> I am not sure though this will happen, so the if condition
            ! for x negative a 0 = 0.5d0 * (1 + (-1) ) is returned

      end function

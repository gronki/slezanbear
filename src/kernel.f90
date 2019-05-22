module kernel

  use globals
  implicit none

  type modelpar
    real(dp) :: alpha  = 5e-6
    real(dp) :: relabs = 8.0
    real(dp) :: hscale = 8000
    real(dp) :: beta   = 0.1
    real(dp) :: skybg  = 3e-7
  end type

contains

  elemental real(dp) function th(x,w)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: w
    if (present(w)) then
      th = 1 / ( 1 + exp(-4*x/w) )
    else
      th = 1 / ( 1 + exp(-4*x) )
    end if
  end function

  !-----------------------------------------------------------------------------

  elemental subroutine kern(par, A0, R, D, H, Hobs, J1, J2, J3)
    type(modelpar), intent(in) :: par
    real(dp), intent(in) :: A0, R, D, H, Hobs
    real(dp), intent(out) :: J1, J2, J3
    real(dp) :: L, costh, cosg, tau1, tau2, chi

    ! distance from the source to the scattering point
    L = sqrt(H**2 + (1 + H / R) * D**2)

    ! emission and scattering angles
    associate (Q => merge(D**2 / (2 * R * L), 0.0_dp, L .gt. 0))
      costh = H / L + Q
      cosg  = H / L - Q * (1 + H / R)
    end associate

    chi = (par % alpha) * (par % relabs)

    tau1 = chi * L * FF(-H / (par % Hscale))
    tau2 = chi * (par % Hscale) * exp(-Hobs / (par % Hscale)) &
         &     * (1 - exp(-(H - Hobs) / (par % Hscale)))

    J1 = (par % alpha) * exp(-H / (par % Hscale)) &
    &   * (A0 / 2) / (2 * pi * L**2 + A0) * merge(1, 0, cosg > 0)
    J2 = E(cosg) * g(costh)
    J3 = exp(-tau1) * exp(-tau2)

  contains

    elemental real(dp) function ff(x) result(y)
      real(dp), intent(in) :: x
      if ( abs(x) > 1e-2 ) then
        y = (exp(x) - 1) / x
      else
        y = 1 + x * (1 + x * (1 + x * (1 + x / 5) / 4) / 3) / 2
      end if
    end function

    elemental real(dp) function g(costh)
      real(dp), intent(in) :: costh
      g = 3 * (1 + costh**2) / 4
    end function

    elemental real(dp) function E0(cosg) result(E)
      real(dp), intent(in) :: cosg
      E = 1.5 * cosg
    end function

    elemental real(dp) function E1(cosg) result(E)
      real(dp), intent(in) :: cosg
      E = 6 * cosg * (1 - cosg)
    end function

    elemental real(dp) function E2(cosg) result(E)
      real(dp), intent(in) :: cosg
      E = 15 * cosg  * (1 - cosg)**2
    end function

    elemental real(dp) function E3(cosg) result(E)
      real(dp), intent(in) :: cosg
      E = 30 * cosg  * (1 - cosg)**3
    end function

    elemental real(dp) function E4(cosg) result(E)
      real(dp), intent(in) :: cosg
      E = 52.5_dp * cosg * (1 - cosg)**4
    end function

    elemental real(dp) function E5(cosg) result(E)
      real(dp), intent(in) :: cosg
      E = 84 * cosg * (1 - cosg)**5
    end function

    elemental real(dp) function E6(cosg) result(E)
      real(dp), intent(in) :: cosg
      E = 126 * cosg * (1 - cosg)**6
    end function

    elemental real(dp) function E7(cosg) result(E)
      real(dp), intent(in) :: cosg
      E = 180 * cosg * (1 - cosg)**7
    end function

    elemental real(dp) function E(cosg)
      real(dp), intent(in) :: cosg
      E = 1.5 * cosg * (1 - par % beta) + 3 * (1 - cosg) * (par % beta)
    end function

  end subroutine

  !-----------------------------------------------------------------------------

  elemental real(dp) function ity2mag(ity) result(mag)
    real(dp), intent(in) :: ity
    mag = 5.5 - 2.5 * log10(ity)
  end function

  elemental real(dp) function mag2ity(mag) result(ity)
    real(dp), intent(in) :: mag
    ity = 10**((5.5 - mag) / 2.5)
  end function

end module

module kernel

  use globals
  implicit none

  type modelpar
    real(fp) :: skybg  = 3e-7
    real(fp) :: alpha  = 0.2 / 1.086
    real(fp) :: relabs = 0.75
    real(fp) :: hscale = 8000
    real(fp) :: beta   = 0.3
  end type

contains

  elemental real(fp) function th(x,w)
    real(fp), intent(in) :: x
    real(fp), intent(in), optional :: w
    if (present(w)) then
      th = 1 / ( 1 + exp(-4*x/w) )
    else
      th = 1 / ( 1 + exp(-4*x) )
    end if
  end function

  !-----------------------------------------------------------------------------

  elemental subroutine kern(par, A0, R, D, H, Hobs, J1, J2, J3)
    type(modelpar), intent(in) :: par
    real(fp), intent(in) :: A0, R, D, H, Hobs
    real(fp), intent(out) :: J1, J2, J3
    real(fp) :: L, Q, costh, cosg, tau1, tau2, sct, chi

    ! distance from the source to the scattering point
    L = sqrt(H**2 + (1 + H / R) * D**2)

    ! emission and scattering angles
    Q = merge(D**2 / (2 * R * L), 0.0_fp, L .ne. 0)
    costh = H / L + Q
    cosg  = H / L - Q * (1 + H / R)

    sct = (par % alpha) / (par % Hscale)
    chi = sct * (par % relabs)

    tau1 = chi * L * FF(-H / (par % Hscale))
    tau2 = chi * (par % Hscale) * exp(-Hobs / (par % Hscale)) &
         &     * (1 - exp(-(H - Hobs) / (par % Hscale)))

    J1 = sct * exp(-H / (par % Hscale)) * (A0 / 2) / (2 * pi * L**2 + A0)
    J2 = E(cosg) * g(costh)
    J3 = exp(-tau1) * exp(-tau2)

  contains

    elemental real(fp) function ff(x) result(y)
      real(fp), intent(in) :: x
      if ( abs(x) > 1e-2 ) then
        y = (exp(x) - 1) / x
      else
        y = 1 + x * (1 + x * (1 + x * (1 + x / 5) / 4) / 3) / 2
      end if
    end function

    elemental real(fp) function g(costh)
      real(fp), intent(in) :: costh
      g = 0.75_fp * (1 + costh**2)
    end function

    elemental real(fp) function E0(cosg) result(E)
      real(fp), intent(in) :: cosg
      E = merge(1, 0, cosg > 0)
    end function

    elemental real(fp) function E1(cosg) result(E)
      real(fp), intent(in) :: cosg
      E = 1.5 * cosg * merge(1, 0, cosg > 0)
    end function

    elemental real(fp) function E2(cosg) result(E)
      real(fp), intent(in) :: cosg
      E = 6 * cosg * (1 - cosg) * merge(1, 0, cosg > 0)
    end function

    elemental real(fp) function E3(cosg) result(E)
      real(fp), intent(in) :: cosg
      E = 15  * cosg * (1 - cosg)**2 * merge(1, 0, cosg > 0)
    end function

    elemental real(fp) function E4(cosg) result(E)
      real(fp), intent(in) :: cosg
      E = 30  * cosg * (1 - cosg)**3 * merge(1, 0, cosg > 0)
    end function

    elemental real(fp) function E(cosg)
      real(fp), intent(in) :: cosg
      E = E0(cosg) * (1 - par % beta) + E4(cosg) * (par % beta)
    end function

  end subroutine

  !-----------------------------------------------------------------------------

  elemental real(fp) function ity2mag(ity) result(mag)
    real(fp), intent(in) :: ity
    mag = 5.5 - 2.5 * log10(ity)
  end function

  elemental real(fp) function mag2ity(mag) result(ity)
    real(fp), intent(in) :: mag
    ity = 10**((5.5 - mag) / 2.5)
  end function

end module

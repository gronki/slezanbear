module kernel

  use globals
  implicit none

  type modelpar
    real :: skybg = 0
    real :: alpha = 0.4
    real :: relabs = 0.75
    real :: hscale = 5000.0
  end type

contains

  elemental function th(x,w)
    real, intent(in) :: x
    real, intent(in), optional :: w
    real :: th
    if (present(w)) then
      th = 1 / ( 1 + exp(-4*x/w) )
    else
      th = 1 / ( 1 + exp(-4*x) )
    end if
  end function

  !-----------------------------------------------------------------------------

  elemental subroutine kern(par, A0, R, D, H, Hobs, J0, J1)
    type(modelpar), intent(in) :: par
    real, intent(in) :: A0, R, D, H, Hobs
    real, intent(out) :: J0, J1
    real :: L, Q, costh, cosg, tau1, tau2, sct, chi

    ! distance from the source to the scattering point
    L = sqrt(H**2 + (1 + H / R) * D**2)

    ! emission and scattering angles
    Q = D**2 / (2*R*H)
    costh = H / L * (1 + Q)
    cosg  = H / L * (1 - Q * (1 + H / R))

    sct = (par % alpha) / (par % Hscale)
    chi = sct * (par % relabs)

    tau1 = chi * L * FF(-H / (par % Hscale))
    tau2 = chi * (par % Hscale) * exp(-Hobs / (par % Hscale)) &
         &     * (1 - exp(-(H - Hobs) / (par % Hscale)))

    J0 = E(cosg) * A0 / (2 * pi * L**2 + A0) * exp(-tau1) / 2
    J1 = sct * exp(-H / (par % Hscale)) * g(costh) * exp(-tau2)

  contains

    elemental function ff(x) result(y)
      real, intent(in) :: x
      real :: y
      if ( abs(x) > 1e-2 ) then
        y = (exp(x) - 1) / x
      else
        y = 1 + x * (1 + x * (1 + x * (1 + x / 5) / 4) / 3) / 2
      end if
    end function

    elemental function g(costh)
      real, intent(in) :: costh
      real :: g
      g = 1
    end function

    elemental function E0(cosg) result(E)
      real, intent(in) :: cosg
      real :: E
      E = merge(1, 0, cosg > 0)
    end function

    elemental function E1(cosg) result(E)
      real, intent(in) :: cosg
      real :: E
      E = 1.5 * cosg * merge(1, 0, cosg > 0)
    end function

    elemental function E2(cosg) result(E)
      real, intent(in) :: cosg
      real :: E
      E = 6 * cosg * (1 - cosg) * merge(1, 0, cosg > 0)
    end function

    elemental function E3(cosg) result(E)
      real, intent(in) :: cosg
      real :: E
      E = 15  * cosg * (1 - cosg)**2 * merge(1, 0, cosg > 0)
    end function

    elemental function E4(cosg) result(E)
      real, intent(in) :: cosg
      real :: E
      E = 30  * cosg * (1 - cosg)**3 * merge(1, 0, cosg > 0)
    end function

    elemental function E(cosg)
      real, intent(in) :: cosg
      real :: E
      E = E0(cosg)
    end function

  end subroutine


end module

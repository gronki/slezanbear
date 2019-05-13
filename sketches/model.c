#include <math.h>

const double Re = 6371;
const double deg_in_rad = 180 / M_PI;

typedef struct {
  float lat;
  float lng;
} latlng;

typedef struct {
  double alpha;
  double krel;
  double Hscale;
  double albedo;
  double cosgmax;
} model_params;

//------------------------------------------------------------------------------

double FF(double x) {
  if ( fabs(x) > 1e-2 ) {
    return (exp(x) - 1) / x;
  } else {
    return 1 + x * (1 + x * (1 + x * (1 + x / 5) / 4) / 3) / 2;
  }
}

//------------------------------------------------------------------------------

inline double angdist(double lat1, double lng1, double lat2, double lng2) {
  return acos(cos(lat1) * cos(lat2) * cos(lng1 - lng2) + sin(lat1) * sin(lat2));
}

inline void geo2xyz(double lat, double lng, double *x, double *y, double *z) {
  *x = cos(lat / deg_in_rad) * cos(lng / deg_in_rad);
  *y = cos(lat / deg_in_rad) * sin(lng / deg_in_rad);
  *z = sin(lat / deg_in_rad);
}

inline void xyz2geo(double x, double y, double z, double *lat, double *lng) {
  *lat = atan2(z, sqrt(x*x + y*y)) * deg_in_rad;
  *lng = atan2(y, x) * deg_in_rad;
}

inline void arr2geo(float gt[6], float i, float j, float *lat, float *lng) {
  *lng = gt[0] + gt[1] * i + gt[2] * j;
  *lat = gt[3] + gt[4] * i + gt[5] * j;
}

inline void geo2arr(float gt[6], float lat, float lng, float *i, float *j) {
  double det = gt[1] * gt[5] - gt[2] * gt[4];
  *i = (   gt[5] * (lng - gt[0]) - gt[2] * (lat - gt[3]) ) / det + 1;
  *j = ( - gt[4] * (lng - gt[0]) + gt[1] * (lat - gt[3]) ) / det + 1;
}

//------------------------------------------------------------------------------

double g(double costh) { return 1; }
double E(double cosg) {
  return 1;
}

double kernel(model_params *par, double A0, double R, double D, double H, double Hobs) {

  double L, Q, costh, cosg, kernel, tau1, tau2, J0, J1, J2, sct, chi;

  L = sqrt(H*H + (1 + H / R) * D*D);

  // emission and scattering angles
  Q = (D * D) / (2 * R * H);
  costh = H / L * (1 + Q);
  cosg  = H / L * (1 - Q * (1 + H / R));

  sct = (par -> alpha) / (par -> Hscale);
  chi = (par -> alpha) / (par -> Hscale) * (par -> krel);

  tau1 = chi * L * FF(-H / (par -> Hscale));
  tau2 = chi * (par -> Hscale) * exp(-Hobs / (par -> Hscale))
    * (1 - exp(-(H - Hobs) / (par -> Hscale)));

  J0 = E(cosg) / 2 * A0 / (2 * M_PI * L*L + A0) * exp(-tau1);
  J1 = sct * exp(- H / (par -> Hscale)) * g(costh);
  J2 = exp(-tau2);

  return J0 * J1 * J2;
}


float interpolmap( int nx, int ny, float map[nx][ny], float gt[6],
    float lat, float lng) {
  int i, ki, kj;
  float ri,rj,xi,xj;

  geo2arr(gt, lat, lng, &xi, &xj);

  ri = xi - floor(xi);
  rj = xj - floor(xj);

  ki = round(xi - ri);
  kj = round(xj - rj);

  if ( ki >= 0 && ki < nx - 1 && kj >= 0 && kj < ny - 1 ) {
    return    map[ki  ][kj  ] * (1 - ri)  * (1 - rj)
            + map[ki+1][kj  ] * ri        * (1 - rj)
            + map[ki  ][kj+1] * (1 - ri)  * rj
            + map[ki+1][kj+1] * ri        * rj;
  }
  return 0;
}


void a() {
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      if (mapi[i][j] == 0) continue;
      par1.alpha = par.alpha * exp(-hsrc[i][j] / (par.hscale));

      for (int k = 0; k < nh; k++) {

        checkline(maph, gth,
          lat,          lng,          hsct[k],
          latsrc[i][j], lngsrc[i][j], hsrc[i][j],
          &Jatt[k]);

        if ( Jatt[k] < 1e-2 ) {
          J[k] = 0;
          continue;
        }

        kernel(&par1, asrc[i][j], Re + hsrc[i][j],
          Re * angdist(lat, lng, latsrc[i][j], lngsrc[i][j]),
          hsct[k] - hsrc[i][j], hobs - hsrc[i][j], &J0[k], &J1[k]);

        J[k] = J0[k] * J1[k] * Jatt[k] * mapi[i][j] * 1e-9;
      }
    }
  }
}


  // pure real function integrate(y, x) result(integral)
  //   real, intent(in), dimension(:) :: x, y
  //   integer i
  //   integral = 0
  //   do i = 1, size(y) - 1
  //     integral = integral + (y(i) + y(i+1)) * (x(i+1) - x(i)) / 2
  //   end do
  // end function
  //
  // !-----------------------------------------------------------------------------
  //
  // subroutine checkline(lat1, lng1, h1, lat2, lng2, h2, attn)
  //   integer, parameter :: nsectmax = 1000
  //   integer :: nsect,i
  //   real, parameter :: smooth_meters = 1.5
  //   type(geopoint) :: p1, p2
  //   real, intent(out) :: attn
  //   real, dimension(nsectmax) :: t, latx, lngx, hx, hterr
  //   real :: dist
  //
  //   nsect = max(min(nsectmax, int(radearth * angdist(lat1, lng2, lat2, lng2) / 0.01)), 5)
  //
  //   forall (i = 1:nsect) t(i) = (i - 1) / real(nsect - 1)
  //
  //   call georay(lat1, lng1, h1 / radearth, lat2, lng2, h2 / radearth, t(1:nsect), latx(1:nsect), lngx(1:nsect), hx(1:nsect))
  //
  //   do i = 1, nsect
  //     call interpolmap(maph, gth, latx(i), lngx(i), hterr(i))
  //   end do
  //
  //   attn = th(minval(hx(1:nsect) * radearth - hterr(1:nsect) * 1e-3), smooth_meters * 1e-3)
  //
  // end subroutine checkline
  //
  // !-----------------------------------------------------------------------------
  //
  // subroutine sbm_map(lat0, lng0, mapi, gti, maph, gth, isky, extn)
  //
  //   real :: hobs, A0, tauz
  //   real, dimension(size(mapi,1),size(mapi,2)) :: hsrc, asrc, lat, lng
  //   integer i,j
  //   integer, parameter :: nh = 500
  //   real, dimension(6) :: gti, gth, gto
  //   real J0(nh), J1(nh), JJ(nh)
  //
  //
  //   do concurrent (i = 1:size(mapi,1), j = 1:size(mapi,2))
  //     call arr2geo(gti, i, j, lat(i,j), lng(i,j))
  //     asrc(i,j) = abs(gti(2)*gti(6)) / deg_in_rad**2 * radearth**2 * cos(lat(i,j) / deg_in_rad)
  //     call interpolmap(maph, gth, lat(i,j), lng(i,j), hsrc(i,j))
  //   end do
  //
  //   call interpolmap(maph, gth, lat0, lng0, hobs)
  //
  //   Hobs = Hobs * 1e-3
  //   Hsrc(:,:) = Hsrc * 1e-3
  //
  //   tauz = (alpha / Hscale * relabs) * Hscale * exp(-Hobs / Hscale)
  //   isky = I0 * exp(-tauz)
  //   do concurrent (i = 1:size(mapi,1), j = 1:size(mapi,2))
  //     call kernel(alpha * exp(-hsrc(i,j) / Hscale), relabs, Hscale, albedo, cosgmax, A0, radearth + hsrc(i,j), radearth * angdist(lat0,lng0,lat,lng), h - hsrc(i,j), hobs - hsrc(i,j), J0, J1)
  //     JJ(:) = J0(:) * J1(:)
  //     isky = isky + integrate(JJ,h) * mapi(i,j)
  //   end do
  //
  //
  // end subroutine

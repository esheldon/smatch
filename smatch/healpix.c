#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "healpix.h"
#include "defs.h"
#include "vector.h"



int64 hpix_npix(int64 nside) {
    return 12*nside*nside;
}

double hpix_area(int64 nside) {
    int64 npix = hpix_npix(nside);
    return 4.0*M_PI/npix;
}


struct healpix* hpix_new(int64 nside) {
    struct healpix* hpix;

    if (nside < 1 || nside > NS_MAX) {
        printf("nside out of range [%d, %d]\n", 1, NS_MAX);
        exit(EXIT_FAILURE);
    }

    hpix = malloc(sizeof(struct healpix));
    if (hpix == NULL) {
        printf("Could not allocate struct healpix\n");
        exit(EXIT_FAILURE);
    }

    hpix->nside = nside;
    hpix->npix = hpix_npix(nside);
    hpix->area = hpix_area(nside);
    hpix->ncap = 2*nside*(nside-1); // number of pixels in the north polar cap

    return hpix;
}

// usage:  hpix=hpix_delete(hpix);
struct healpix* hpix_delete(struct healpix* hpix) {
    free(hpix);
    return NULL;
}

// ra,dec in degrees
// radius in radians
void hpix_disc_intersect(
        const struct healpix* hpix,
        double x, double y, double z, double radius, 
        lvector* listpix) {

    // this is from the f90 code
    // this number is acos(2/3)
    //double fudge = 0.84106867056793033/hpix->nside; // 1.071* half pixel size

    // this is from the c++ code
    double fudge = 1.362*M_PI/(4*hpix->nside);

    radius += fudge;
    hpix_disc_contains(hpix, x, y, z, radius, listpix);
}

void hpix_disc_contains(
        const struct healpix* hpix,
        double x0, double y0, double z0, double radius, 
        lvector* listpix) {

    //double vector0[3];
    int64 nside=hpix->nside;
    double cosang = cos(radius);

    // this does not alter the storage
    vector_resize(listpix, 0);

    double dth1 = 1. / (3.0*nside*nside);
    double dth2 = 2. / (3.0*nside);

    double phi0=0.0;
    if ((x0 != 0.) || (y0 != 0.)) {
        // in (-Pi, Pi]
        phi0 = atan2(y0, x0);
    }
    double cosphi0 = cos(phi0);
    double a = x0*x0 + y0*y0;

    //     --- coordinate z of highest and lowest points in the disc ---
    double rlat0  = asin(z0);    // latitude in RAD of the center
    double rlat1  = rlat0 + radius;
    double rlat2  = rlat0 - radius;
    double zmax;
    if (rlat1 >=  M_PI_2) {
        zmax =  1.0;
    } else {
        zmax = sin(rlat1);
    }
    int64 irmin = hpix_ring_num(hpix, zmax);
    irmin = i64max(1, irmin-1); // start from a higher point, to be safe

    double zmin;
    if (rlat2 <= -M_PI_2) {
        zmin = -1.;
    } else {
        zmin = sin(rlat2);
    }
    int64 irmax = hpix_ring_num(hpix, zmin);
    irmax = i64min(4*nside-1, irmax + 1); // go down to a lower point

    //double z, tmp=0;
    double tmp=0;
    int64 iz=0;
    for (iz=irmin; iz<= irmax; iz++) {

        double z;
        if (iz <= nside-1) { // north polar cap
              z = 1.  - iz*iz*dth1;
        } else if (iz <= 3*nside) { // tropical band + equat.
            z = (2*nside-iz) * dth2;
        } else {
            tmp = 4*nside-iz;
            z = - 1. + tmp*tmp*dth1;
        }
        double b = cosang - z*z0;
        double c = 1. - z*z;
        //double x = (cosang-z*z0)/sqrt((1-z0)*(1+z0));

        double dphi;
        if ((x0==0.) && (y0==0.)) {
            dphi=M_PI;
            if (b > 0.) {
                goto SKIP2; // out of the disc, 2008-03-30
            }
            goto SKIP1;
        } 

        double cosdphi = b / sqrt(a*c);
        if (fabs(cosdphi) <= 1.) {
              dphi = acos(cosdphi); // in [0,Pi]
        } else {
            if (cosphi0 < cosdphi) {
                goto SKIP2; // out of the disc
            }
            dphi = M_PI; // all the pixels at this elevation are in the disc
        }
SKIP1:
        hpix_in_ring(hpix, iz, phi0, dphi, listpix);

SKIP2:
        // we have to put something here
        continue;

    }


}

int64 i64max(int64 v1, int64 v2) {
    return v1 > v2 ? v1 : v2;
}
int64 i64min(int64 v1, int64 v2) {
    return v1 < v2 ? v1 : v2;
}

void hpix_in_ring(
        const struct healpix* hpix, 
        int64 iz, 
        double phi0, 
        double dphi, 
        lvector* plist) {

    int64 nr, ir, ipix1,i;
    double shift=0.5;
    int64 nside = hpix->nside;

    if (iz<nside) {
        // north pole
        ir = iz;
        nr = ir*4;
        ipix1 = 2*ir*(ir-1);        //    lowest pixel number in the ring
    } else if (iz>(3*nside)) {
        // south pole
        ir = 4*nside - iz;
        nr = ir*4;
        ipix1 = hpix->npix - 2*ir*(ir+1); // lowest pixel number in the ring
    } else {
        // equatorial region
        ir = iz - nside + 1;           //    within {1, 2*nside + 1}
        nr = nside*4;
        if ((ir&1)==0) shift = 0;
        ipix1 = hpix->ncap + (ir-1)*nr; // lowest pixel number in the ring
    }

    int64 ipix2 = ipix1 + nr - 1;  //    highest pixel number in the ring
 

    if (dphi > (M_PI-1e-7)) {
        for (i=ipix1; i<=ipix2; ++i) {
            vector_push(plist, i);
        }
    } else {

        // M_1_PI is 1/pi
        int64 ip_lo = (int64)( floor(nr*.5*M_1_PI*(phi0-dphi) - shift) )+1;
        int64 ip_hi = (int64)( floor(nr*.5*M_1_PI*(phi0+dphi) - shift) );
        int64 pixnum = ip_lo+ipix1;
        if (pixnum<ipix1) {
            pixnum += nr;
        }
        for (i=ip_lo; i<=ip_hi; ++i, ++pixnum) {
            if (pixnum>ipix2) {
                pixnum -= nr;
            }
            vector_push(plist, pixnum);
        }
    }

}


int64 hpix_ring_num(const struct healpix* hpix, double z) {
    int64 nside=hpix->nside;

    // rounds double to nearest long long int
    int64 iring = llrintl( nside*(2.-1.5*z) );

    // north cap
    if (z > M_TWOTHIRD) {
        iring = llrintl( nside* sqrt(3.*(1.-z)) );
        if (iring == 0) {
            iring = 1;
        }
    } else if (z < -M_TWOTHIRD) {
        iring = llrintl( nside* sqrt(3.*(1.+z)) );

        if (iring == 0) {
            iring = 1;
        }
        iring = 4*nside - iring;
    }

    return iring;
}

int64 hpix_eq2pix(const struct healpix* hpix, double ra, double dec) {
    int64 nside=hpix->nside;
    int64 ipix=0;
    double theta=0, phi=0;

    hpix_radec_degrees_to_thetaphi_radians(ra, dec, &theta, &phi);

    double z = cos(theta);
    double za = fabs(z);

    // in [0,4)
    double tt = fmod(phi, M_TWO_PI)/M_PI_2;

    if (za <= M_TWOTHIRD) {
        double temp1 = nside*(.5 + tt);
        double temp2 = nside*.75*z;

        int64 jp = (int64)(temp1-temp2); // index of  ascending edge line
        int64 jm = (int64)(temp1+temp2); // index of descending edge line
        int64 ir = nside + 1 + jp - jm;  // in {1,2n+1} (ring number counted from z=2/3)
        int64 kshift = 1 - (ir % 2);      // kshift=1 if ir even, 0 otherwise
        
        int64 nl4 = 4*nside;
        int64 ip = (int64)( ( jp+jm - nside + kshift + 1 ) / 2); // in {0,4n-1}

        ip = ip % nl4;

        ipix = hpix->ncap + nl4*(ir-1) + ip;

    } else { 
        // North & South polar caps
        double tp = tt - (int64)(tt);   // MODULO(tt,1.0_dp)

        double tmp = nside * sqrt( 3.0*(1.0 - za) );
        int64 jp = (int64)(tp*tmp);              // increasing edge line index
        int64 jm = (int64)((1.0 - tp) * tmp); // decreasing edge line index

        int64 ir = jp + jm + 1;        // ring number counted from the closest pole
        int64 ip = (int64)( tt * ir);     // in {0,4*ir-1}

        if (ip >= 4*ir) {
            ip = ip - 4*ir;
        }
        if (z>0.) {
            ipix = 2*ir*(ir-1) + ip;
        } else {
            ipix = hpix->npix - 2*ir*(ir+1) + ip;
        }

    }

    return ipix;
}


void hpix_eq2xyz(double ra, double dec, double* x, double* y, double* z) {

    double theta=0, phi=0;

    hpix_radec_degrees_to_thetaphi_radians(ra, dec, &theta, &phi);

    double sintheta = sin(theta);
    *x = sintheta * cos(phi);
    *y = sintheta * sin(phi);
    *z = cos(theta);

}

void hpix_radec_degrees_to_thetaphi_radians(double ra, double dec, double* theta, double* phi) {

    if (ra < 0.0 || ra > 360.) {
        printf("ra out of range [0,360.]\n");
        exit(EXIT_FAILURE);
    }
    if (dec < -90. || dec > 90.) {
        printf("dec out of range [-90.,90.]\n");
        exit(EXIT_FAILURE);
    }

    *phi = ra*D2R;
    *theta = -dec*D2R + M_PI_2;
}

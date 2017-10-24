/****
 * This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * In jurisdictions that recognize copyright laws, the author or authors
 * of this software dedicate any and all copyright interest in the
 * software to the public domain. We make this dedication for the benefit
 * of the public at large and to the detriment of our heirs and
 * successors. We intend this dedication to be an overt act of
 * relinquishment in perpetuity of all present and future rights to this
 * software under copyright law.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 * For more information, please refer to <https://unlicense.org>
 */

/*
 * Flatten:
 * Classes for applying map projections to geographic coordinates, and vice versa.
 *
 * Variable names heavily influenced by the IOGP's geomatics guidance note 7, part 2.
 * Also based on Snyder's "Map Projections--A Working Manual" (USGS Professional Paper 1395).
 *
 */

#ifndef __H_FLATTEN
#define __H_FLATTEN

#include <cmath>


namespace Flatten {

    const double pi = 3.141592653589793;

    class MapProjection {
    public:
        virtual void geo2proj(double lat, double lon, double& x, double& y) = 0;
        virtual void proj2geo(double x, double y, double& lat, double& lon) = 0;

        virtual unsigned char enviProjectionType() = 0;

        static double deg2rad(double deg) {
            return (deg / 180.0) * pi;
        }
        static double rad2deg(double rad) {
            return (rad / pi) * 180.0;
        }
    };

    class EquirectangularProjection : public MapProjection {
    public:
        /* 
         * Equations originally from http://mathworld.wolfram.com/CylindricalEquidistantProjection.html
         * and its citations 
         */

        // - phi_1: latitude of natural origin (degrees)
        // - lambda_o: longitude of natural origin (degrees)
        // - radius: radius (local) (m)
        // - FE: false easting (x direction) (m)
        // - FN: false northing (y direction) (m)
        EquirectangularProjection(double phi_1, double lambda_o, double radius, double FE, double FN) :
            phi_1(deg2rad(phi_1)),
            lambda_o(deg2rad(lambda_o)),
            r(radius),
            FE(FE),
            FN(FN)
        {
        }

        virtual void geo2proj(double lat, double lon, double& x, double& y) {
            x = FE + r * (lon - lambda_o) * std::cos(phi_1);
            y = FN + r * lat;
        }

        virtual void proj2geo(double x, double y, double& lat, double& lon) {
            lat = (y - FN) / r;
            lon = lambda_o + (x - FE) / (r * std::cos(phi_1));
        }

        virtual unsigned char enviProjectionType() {
            return 17;
        }

    private:
        double
            phi_1, // lat of natural origin, rad
            lambda_o, // lon of natural origin, rad
            r, // radius (local), m
            FE, // false easting, m
            FN; // false northing, m
    };

    class PolarStereographicProjection : public MapProjection {
    public:

        // - phi_o: latitude of origin (must be +/- 90) (degrees)
        // - lambda_o: longitude of origin (degrees)
        // - rad_pol: polar radius (m)
        // - rad_eq: equatorial radius (m)
        // - k_o: scale factor at natural origin
        // - FE: false easting (x direction) (m)
        // - FN: false northing (y direction) (m)
        PolarStereographicProjection(double phi_o, double lambda_o, double rad_pol, double rad_eq, double k_o, double FE, double FN) :
            phi_o(deg2rad(phi_o)),
            lambda_o(deg2rad(lambda_o)),
            rad_eq(rad_eq),
            k_o(k_o),
            FE(FE),
            FN(FN)
        {
            double radius_ratio = rad_pol / rad_eq;
            this->ecc = std::sqrt(1 - radius_ratio*radius_ratio); // eccentricity
        }

        virtual void geo2proj(double lat, double lon, double& x, double& y) {

            double t;
            if (phi_o < 0) { // South Pole centered
                t = std::tan(pi / 4 + lat / 2) / std::pow((1 + ecc*sin(lat)) / (1 - ecc*std::sin(lat)), ecc / 2);
            }
            else { // North Pole centered
                t = std::tan(pi / 4 - lat / 2) * std::pow((1 + ecc*sin(lat)) / (1 - ecc*std::sin(lat)), ecc / 2);
            }

            double rho = (2 * rad_eq*k_o*t) / std::sqrt(std::pow(1 + ecc, 1 + ecc) * std::pow(1 - ecc, 1 - ecc));

            // In the IOGP docs, theta here is known as omega in the N. Pole case
            double theta = lon - lambda_o;

            double dE = rho*std::sin(theta);
            double dN = rho*std::cos(theta);

            if (phi_o >= 0) {
                dN = -dN;
            }

            x = dE + FE;
            y = dN + FN;

        }

        virtual void proj2geo(double x, double y, double& lat, double& lon) {
            //rho_prime = std::sqrt((x-false_easting)*(x-false_easting) + (y-false_northing)*(y-false_northing));
            double rho_prime = std::hypot(x - FE, y - FN);
            double t_prime = rho_prime * std::sqrt(std::pow(1 + ecc, 1 + ecc) * std::pow(1 - ecc, 1 - ecc)) / (2 * rad_eq * k_o);

            double lat_conformal;
            if (phi_o < 0) { // South Pole centered
                lat_conformal = 2 * atan(t_prime) - pi / 2;
            }
            else { // North Pole centered
                lat_conformal = pi / 2 - 2 * atan(t_prime);
            }

            lat = lat_conformal +
                (std::pow(ecc, 2) / 2 + 5 * std::pow(ecc, 4) / 24 + std::pow(ecc, 6) / 12 + 13 * std::pow(ecc, 8) / 360) * std::sin(2 * lat_conformal) +
                (7 * std::pow(ecc, 4) / 48 + 29 * std::pow(ecc, 6) / 240 + 811 * std::pow(ecc, 8) / 11520) * std::sin(4 * lat_conformal) +
                (7 * std::pow(ecc, 6) / 120 + 81 * std::pow(ecc, 8) / 1120) * std::sin(6 * lat_conformal) +
                (4279 * std::pow(ecc, 8) / 161280) * std::sin(8 * lat_conformal);

            if (phi_o < 0) {
                lon = lambda_o + atan2(x - FE, y - FN);
            }
            else {
                lon = lambda_o + atan2(x - FE, FN - y);
            }
        }

        virtual unsigned char enviProjectionType() {
            return 31;
        }

    private:
        double
            phi_o, // lat of origin, rad
            lambda_o, // lon of origin, rad
            rad_eq, // equatorial radius (semimajor axis), m
            ecc, // eccentricity, unitless
            k_o, // scale factor at natural origin, unitless
            FE, // false easting, m
            FN; // false northing, m
    };
};

#endif /* __H_FLATTEN */

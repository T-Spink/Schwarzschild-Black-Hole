#ifndef GLOBALS_H
#define GLOBALS_H

#include "integrator.h"
#include "constants.h"
#include <Eigen/Dense>

struct globals
{
    // ========================================= // 
    // inputs

    // numerical params
    const double delta = 1e-1;
    const IntegrationMethod int_method = RK45;
    const double rk45_tol = 1e-6;

    // physical params
    const double Rs    = 1;
    const double disk0 = 3*Rs;
    const double diskf = 6.;

    // graphics params/ initial conditions
    const double R0 = 9.f;
    const double theta0 = M_PI_2-0.1;
    const double phi0 = 0.f;
    const int n_pixels = 300;
    const int n_window = 500;
    const double fov_degs = 90.;
    const bool show_Rs = false;
    

    // ========================================= // 
    // variables

    unsigned char* pixels = nullptr;
    std::atomic<int> render_count=0;

    int N_pixels;
    double pixel_sf;
    Eigen::Vector<double, 3> position;
    Eigen::Vector<double,3> u,v,n;
    Eigen::Matrix<double,4,4> tetrad;
    Eigen::Matrix<double,3,3> cart_to_spherical;


    // ========================================= // 
    // methods

    globals(void)
    {
        N_pixels = n_pixels*n_pixels;
        pixel_sf = R0*tan(fov_degs*rads_per_degs/2.f) / (n_pixels/2.);

        pixels = new unsigned char[n_pixels * n_pixels * 3];
        memset(pixels, 0, 3*N_pixels*sizeof(unsigned char));

        // position of camera
        position << R0, theta0, phi0;

        // determine normal to and basis vectors of the image plane 
        n << 
            R0*sin(theta0)*cos(phi0), 
            R0*sin(theta0)*sin(phi0),
            R0*cos(theta0)
        ;
        u << -n[2]/n.norm(), 0, n[0]/n.norm();
        v = n.cross(u).normalized();
    
        // initialise tetrad matrix
        const double f = (1.f-Rs/R0);
        tetrad = Eigen::Matrix<double,4,4>::Zero();
        tetrad(0,0) = 1.f/std::sqrt(f);
        tetrad(1,1) = std::sqrt(f);
        tetrad(2,2) = 1.f/R0;
        tetrad(3,3) = 1.f/(R0*sin(theta0));

        // calculate tranformation matrix
        cart_to_spherical <<  
            (sin(theta0)*cos(phi0)),
            (sin(theta0)*sin(phi0)),
            (cos(theta0)          ),
            (cos(theta0)*cos(phi0)),
            (cos(theta0)*sin(phi0)),
            (-sin(theta0)         ),
            (-sin(phi0)           ),
            (cos(phi0)            ),
            (0.                   )
        ;
    };
    
    ~globals(void)
    {
        delete[] pixels;
    };


    // ========================================= // 
};


#endif
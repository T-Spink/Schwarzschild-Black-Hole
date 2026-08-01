#ifndef GEODESICS_H
#define GEODESICS_H

#include "integrator.h"
#include <cstring>

enum class EndDomain{
    Error,
    BlackHole,
    Disk,
    Space
};

class geodesics
{

public:

    geodesics(void);
    void init_state(const double (&u)[8]);
    void init_Rs(const double Rs_);
    void init_annulus(const double annulus_inner_radius_, const double annulus_outer_radius_);
    bool init_solver(const double delta, const IntegrationMethod int_method, const double tol);
    void reset_state(void);
    ~geodesics(void);

    void traverse_geodesic(const double (&IC)[8], EndDomain &end);
    void get_RGB(EndDomain end, unsigned char (&rgb)[3]);

private:

    double state[8];  // t,r,theta,phi,dt,dR,dtheta,dphi
    double g[4], ginv[4];
    double dg[4][4][4];
    double delta = 1e-1;
    const int max_integration_steps = 1e6;

    double initial_state[8];
    double annulus_inner_radius = 3.;
    double annulus_outer_radius = 6.;
    const double annulus_width_by_2 = 0.05;
    const unsigned char annulus_inner_rgb[3] = {255, 230, 200}; // colour at disk0
    const unsigned char annulus_outer_rgb[3] = {80, 30, 10};  // colour at diskf
    const double R0_prop_for_ray_escaped = 1.05;
    const double Rs_prop_for_in_black_hole = 1.05;
    double R_ray_in_black_hole = 1.05;

    integrator* solver = nullptr;

};

#endif
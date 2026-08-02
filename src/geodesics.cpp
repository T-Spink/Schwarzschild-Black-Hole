#include "geodesics.h"
#include "globals.h"
#include "integrator.h"
#include <cstring>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <algorithm>
#include <numbers>

static void calc_metric(const double (&state)[globals::n_8_vector], double (&g)[globals::n_4_vector]);

geodesics::geodesics(integrator* const solver_ , const double delta_) : solver(solver_), delta(delta_)
{
    R_ray_in_black_hole = Rs_prop_for_in_black_hole * globals::Rs;
}

void geodesics::init_Rs(const double Rs_)
{
    Rs = Rs_;
    R_ray_in_black_hole = Rs_prop_for_in_black_hole * Rs;
}

void geodesics::init_annulus(const double annulus_inner_radius_, const double annulus_outer_radius_)
{
    annulus_inner_radius = annulus_inner_radius_;
    annulus_outer_radius = annulus_outer_radius_;
}

void geodesics::traverse_geodesic(const double (&initial_state)[globals::n_8_vector], EndDomain &end, unsigned char (&rgb)[globals::rgb_length])
{
    // --------------------------------------------------------- //
    // traverse the geodesic until an end domain is reached

    double state[globals::n_8_vector];
    memcpy(&state[0], &initial_state[0], globals::n_8_vector * sizeof(double));
    const double R_ray_escaped = R0_prop_for_ray_escaped * state[1];
    int cc = 0;
    while (++cc)
    {
        // integrate the geodesic equation
        if (!solver->integrate(delta, state)) {
            end = EndDomain::Error;
            break;
        }

        // check end conditions
        double R = state[1];
        double theta = state[2];
        double z = R * std::cos(theta);

        // check if in black hole
        if (R <= R_ray_in_black_hole){
            end = EndDomain::BlackHole;
            break;
        }

        // check if hit disc
        if ((fabs(z) <= annulus_width_by_2) && 
            (annulus_inner_radius < R) && 
            (R < annulus_outer_radius)){
            end = EndDomain::Disk;
            break;
        }

        // check if disapeared
        if (R_ray_escaped <= R) {
            end = EndDomain::Space;
            break;
        }

        // check max iterations
        if (cc >= max_integration_steps){
            end = EndDomain::Error;
            break;
        }
    
    }


    // --------------------------------------------------------- //
    // return an rgb based on end domain

    switch (end)
    {
        case EndDomain::Disk:
        {
            const double R = state[1];

            // estimate base rgb
            const double t = (R - annulus_inner_radius) / (annulus_outer_radius - annulus_inner_radius);
            const double theta = std::sin(t * std::numbers::pi / 2.);

            // linear interpolate between the colours
            for (int ii = 0; ii < globals::rgb_length; ++ii)
                rgb[ii] = annulus_inner_rgb[ii] + (annulus_outer_rgb[ii] - annulus_inner_rgb[ii]) * theta;

            // calculate relativistic effects
            const double tt = 1. / std::sqrt(1. - 3. / (2. * R));
            const Eigen::Vector<double,4> disk_4vel = {tt, 0, 0, tt * std::sqrt(1. / (2. * std::pow(R, 3)))};
            const Eigen::Vector<double,4> cam_4vel = {1. / (1. - 1. / R), 0., 0., 0.};
            double g[globals::n_4_vector];
            calc_metric(state, g);
            double disk_contr = 0, cam_contr = 0;
            for (int kk = 0; kk < globals::n_4_vector; ++kk){
                disk_contr += disk_4vel[kk] * state[kk + globals::n_4_vector] * g[kk]; // on disk
                cam_contr += cam_4vel[kk] * initial_state[kk + globals::n_4_vector] * g[kk]; // at camera
            }
            double intensity_sf = std::pow(cam_contr / disk_contr, 3);

            // apply intensity
            for (int kk = 0; kk < globals::rgb_length; ++kk)
                rgb[kk] = std::min((double)globals::rgb_element_max, intensity_sf * rgb[kk]);

            break;
        }
        case EndDomain::Space:
        case EndDomain::BlackHole:
            memcpy(&rgb[0], &globals::black_rgb[0], globals::rgb_length * sizeof(unsigned char));
            break;
        case EndDomain::Error:
        default:
            memcpy(&rgb[0], &globals::green_rgb[0], globals::rgb_length * sizeof(unsigned char));
            break;
    }

    // --------------------------------------------------------- //

}

void f_geodesic(const double* const state, double* const f)
{
    std::fill_n(f, globals::n_8_vector, 0.);
    double g[globals::n_4_vector], 
    ginv[globals::n_4_vector] = {}, 
    dg[globals::n_4_vector][globals::n_4_vector][globals::n_4_vector] = {}; // g is set to 0 in the calc_metric

    double u[globals::n_8_vector];
    memcpy(&u, state, globals::n_8_vector * sizeof(double));

    /*
    // check null condition satisfied
    double ds=0.;
    for (int ii=0; ii<4; ii++)
        ds+=g(ii,ii)*u[ii+4]*u[ii+4];
    
    // check for numerical drifting
    if (fabs(ds)>1e-2){
        ***
    }
    */

    // calculate values
    calc_metric(u, g);
    
    // ----
    // calculate ginv
    
    for (int ii = 0; ii < globals::n_3_vector; ++ii)
        ginv[ii] = 1. / g[ii];
    
    // handle possible divide by zero at theta=0,pi
    ginv[3] = (g[3] > 1e-6) ? (1. / g[3]) : (1e6);


    // ----
    // calculate dg

    const double R = u[1];
    const double theta = u[2];

    // dgdR
    const double f_ = 1. - globals::Rs / R;
    const double st = std::sin(theta);
    dg[1][0][0] = -globals::Rs / (R * R);
    dg[1][1][1] = -globals::Rs / (R * R * f_ * f_);
    dg[1][2][2] = 2. * R;
    dg[1][3][3] = 2. * R * st * st;
       
    // dgdtheta
    dg[2][3][3] = R * R * std::sin(2. * theta);


    // ----
    // calculate duds=f

    for (int mu = 0; mu < globals::n_4_vector; ++mu)
        f[mu] = u[mu + globals::n_4_vector];

    for (int mu = 0; mu < globals::n_4_vector; ++mu){

        for (int alpha = 0; alpha < globals::n_4_vector; ++alpha){
            for (int beta = 0; beta < globals::n_4_vector; ++beta){

                double lambda = dg[beta][mu][alpha] + dg[alpha][mu][beta] - dg[mu][alpha][beta]; // simplified since tensors are diagonal
                f[mu + globals::n_4_vector] += lambda * u[alpha + globals::n_4_vector] * u[beta + globals::n_4_vector];

            }
        }   
        f[mu + globals::n_4_vector] *= -0.5 * ginv[mu];

    }


    // ----
}

static void calc_metric(const double (&state)[globals::n_8_vector], double (&g)[globals::n_4_vector])
{
    memset(&g, 0, globals::n_4_vector * sizeof(double));

    const double R = state[1];
    const double theta = state[2];
    const double f = (1. - globals::Rs / R);

    g[0] = -f;
    g[1] = 1. / f;
    g[2] = R * R;
    g[3] = R * R * std::sin(theta) * std::sin(theta);  // where theta=0,pi this is zero!

}

#include "geodesics.h"
#include "integrator.h"
#include <cstring>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <algorithm>

static double Rs = 1; // this is currently set by all threads without locking etc ... not safe

static void calc_metric(const double (&state)[8], double (&g)[4]);
void f_geodesic(const double* u_ptr, double* f);

geodesics::geodesics()
{
    R_ray_in_black_hole = Rs_prop_for_in_black_hole * Rs;
    reset_state();
}

bool geodesics::init_solver(const double delta_, const IntegrationMethod int_method, const double tol)
{
    if (solver){
        std::cerr << "cannot init_solver more than once." << std::endl;
        return false;
    }
    solver = new integrator(&f_geodesic, 8, int_method);
    delta = delta_;
    const int max_its = 1e6;
    solver->init_adaptive(tol, max_its);
    return true;
}

void geodesics::init_Rs(const double Rs_)
{
    Rs = Rs_;
    R_ray_in_black_hole = Rs_prop_for_in_black_hole * Rs;
}

void geodesics::init_annulus(const double R0, const double Rf)
{
    annulus_inner_radius = R0;
    annulus_outer_radius = Rf;
}

geodesics::~geodesics(void)
{
    if (solver) 
        delete solver;
}

void geodesics::reset_state(void)
{
    std::fill_n(state, 8, 0.);
    std::fill_n(initial_state, 8, 0.);
    for (int ii = 0; ii < 4; ++ii){
        for (int jj = 0; jj < 4; ++jj)
            std::fill_n(dg[ii][jj], 4, 0.);
    }
    std::fill_n(g, 4, 0.);
    std::fill_n(ginv, 4, 0.);
    delta = 1e-1;
}

void geodesics::get_RGB(EndDomain end, unsigned char (&rgb)[3])
{
    // get colour
    switch (end)
    {
        case EndDomain::Disk:
        {
            const double R = state[1];

            // estimate base rgb
            const double t = (R - annulus_inner_radius) / (annulus_outer_radius - annulus_inner_radius);
            const double theta = std::sin(t * M_PI_2);

            // linear interpolate between the colours
            for (int ii = 0; ii < 3; ++ii)
                rgb[ii] = annulus_inner_rgb[ii] + (annulus_outer_rgb[ii] - annulus_inner_rgb[ii]) * theta;

            // calculate relativistic effects
            const double tt = 1. / std::sqrt(1. - 3. / (2. * R));
            const Eigen::Vector<double,4> disk_4vel = {tt, 0, 0, tt * std::sqrt(1. / (2. * std::pow(R, 3)))};
            const Eigen::Vector<double,4> cam_4vel = {1. / (1. - 1. / R), 0., 0., 0.};
            double g[4] = {};
            calc_metric(state, g);
            double disk_contr = 0, cam_contr = 0;
            for (int kk = 0; kk < 4; ++kk){
                disk_contr += disk_4vel[kk] * state[kk + 4] * g[kk]; // on disk
                cam_contr += cam_4vel[kk] * initial_state[kk + 4] * g[kk]; // at camera
            }
            double intensity_sf = std::pow(cam_contr / disk_contr, 3);

            // apply intensity
            for (int kk = 0; kk < 3; ++kk)
                rgb[kk] = std::min(255., intensity_sf * rgb[kk]);

            break;
        }
        case EndDomain::Space:
        case EndDomain::BlackHole:
            memset(rgb, 0, 3 * sizeof(unsigned char));
            break;
        case EndDomain::Error:
        default:
            rgb[0] = 0;
            rgb[1] = 255;
            rgb[2] = 0;
            break;
    }

}

void geodesics::traverse_geodesic(const double (&IC)[8], EndDomain &end)
{
    if (!solver){
        std::cerr << "solver must be initialised in geodesics" << std::endl;
        end = EndDomain::Error;
        return;
    }

    memcpy(state        , IC, 8 * sizeof(double));
    memcpy(initial_state, IC, 8 * sizeof(double));
    const double R_ray_escaped = R0_prop_for_ray_escaped * state[1];

    for (int cc = 0; cc < max_integration_steps; ++cc)
    {
        // integrate the geodesic equation
        if (!solver->integrate(delta, state)) {
            end = EndDomain::Error;
            return;
        }

        // check end conditions
        double R = state[1];
        double theta = state[2];
        double z = R * std::cos(theta);

        // check if in black hole
        if (R <= R_ray_in_black_hole){
            end = EndDomain::BlackHole;
            return;
        }

        // check if hit disc
        if ((fabs(z) < annulus_width_by_2) && 
            (annulus_inner_radius < R) && 
            (R < annulus_outer_radius)){
                end = EndDomain::Disk;
                return;
            }

        // check if disapeared
        if (R_ray_escaped < R) {
            end = EndDomain::Space;
            return;
        }
    
    }

    // error
    end = EndDomain::Error;
    return;

}

void f_geodesic(const double* const state, double* const f)
{
    std::fill_n(f, 8, 0.);
    double g[4], ginv[4] = {}, dg[4][4][4] = {}; // g is set to 0 in the calc_metric

    double u[8];
    memcpy(&u, state, 8 * sizeof(double));

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
    
    for (int ii = 0; ii < 3; ++ii)
        ginv[ii] = 1. / g[ii];
    
    // handle possible divide by zero at theta=0,pi
    ginv[3] = (g[3] > 1e-6) ? (1. / g[3]) : (1e6);


    // ----
    // calculate dg

    const double R = u[1];
    const double theta = u[2];

    // dgdR
    const double f_ = 1. - Rs / R;
    const double st = std::sin(theta);
    dg[1][0][0] = -Rs / (R * R);
    dg[1][1][1] = -Rs / (R * R * f_ * f_);
    dg[1][2][2] = 2. * R;
    dg[1][3][3] = 2. * R * st * st;
       
    // dgdtheta
    dg[2][3][3] = R * R * std::sin(2. * theta);


    // ----
    // calculate duds=f

    for (int mu = 0; mu < 4; ++mu)
        f[mu] = u[mu + 4];

    for (int mu = 0; mu < 4; ++mu){

        for (int alpha = 0; alpha < 4; ++alpha){
            for (int beta = 0; beta < 4; ++beta){

                double lambda = dg[beta][mu][alpha] + dg[alpha][mu][beta] - dg[mu][alpha][beta]; // simplified since tensors are diagonal
                f[mu + 4] += lambda * u[alpha + 4] * u[beta + 4];

            }
        }   
        f[mu + 4] *= -0.5 * ginv[mu];

    }


    // ----
}

static void calc_metric(const double (&state)[8], double (&g)[4])
{
    memset(&g, 0, 4 * sizeof(double));

    const double R = state[1];
    const double theta = state[2];
    const double f = (1. - Rs / R);

    g[0] = -f;
    g[1] = 1. / f;
    g[2] = R * R;
    g[3] = R * R * std::sin(theta) * std::sin(theta);  // where theta=0,pi this is zero!

}

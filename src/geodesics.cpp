#include "geodesics.h"
#include "integrator.h"
#include <cstring>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <algorithm>

geodesics::geodesics()
{
    reset();
    Rs = 1;
    disk_R0 = 3;
    disk_Rf = 6;
    delta=1e-1;
    solver = nullptr;
    end = Error;
}

void geodesics::init_solver(const double delta_, const IntegrationMethod int_method, const double tol)
{
    if (solver){
        std::cerr<<"cannot init_solver twice."<<std::endl;
        exit(1);
    }
    solver = new integrator(this, &geodesics::rhs, 8);
    delta = delta_;
    solver->init_integration(int_method);
    solver->init_rk45_tol(tol);
}

geodesics::~geodesics(void)
{
    if (solver) delete solver;
}

void geodesics::init_disk(const double R0, const double Rf)
{
    disk_R0 = R0;
    disk_Rf = Rf;
}

void geodesics::reset(void)
{
    std::fill_n(state, 8, 0.);
    std::fill_n(initial_state, 8, 0.);
    for (int ii=0; ii<4; ii++){
        for (int jj=0; jj<4; jj++)
            std::fill_n(dg[ii][jj], 4, 0.);
    }
    std::fill_n(g, 4, 0.);
    std::fill_n(ginv, 4, 0.);
    end = Error;
}

void geodesics::GetRGB(const double (&IC)[8], unsigned char (&rgb)[3])
{
    end = FindEndDomain(IC);    

    // get colour
    memset(rgb, 0, 3*sizeof(unsigned char));
    switch (end)
    {
        case Space:
        case BlackHole:
            //memset(&rgb,0,3*sizeof(unsigned char)); 
            return;

        case Disk:
        {
            const double R = state[1];

            // estimate base rgb
            const double t = (R-disk_R0)/(disk_Rf-disk_R0);
            const double theta = sin(t*M_PI_2);

            // linear interpolate between the colours
            for (int ii=0; ii<3; ii++)
                rgb[ii] = rgb0[ii] + (rgbf[ii]-rgb0[ii])*theta;

            // calculate relativistic effects
            const double tt = 1./sqrt(1.-3./(2.*R));
            const Eigen::Vector<double,4> disk_4vel = {tt, 0, 0, tt*sqrt(1/(2.*pow(R,3)))};
            const Eigen::Vector<double,4> cam_4vel = {1./(1.-1./R), 0., 0., 0.};
            calc_metric(state);
            double disk_contr = 0, cam_contr = 0;
            for (int kk=0; kk<4; kk++){
                disk_contr += disk_4vel[kk]*        state[kk+4]*g[kk]; // on disk
                cam_contr += cam_4vel[kk]*initial_state[kk+4]*g[kk]; // at camera
            }
            double intensity_sf = pow(cam_contr/disk_contr,3);

            // apply intensity
            for (int kk=0; kk<3; kk++)
                rgb[kk] = std::min(255., intensity_sf*rgb[kk]);

            return;
        }

        case Error:
        default:
            rgb[0] = 0;
            rgb[1] = 255;
            rgb[2] = 0;
            return;
    }

}

EndDomain geodesics::FindEndDomain(const double (&IC)[8])
{
    if (!solver){
        std::cerr<<"solver must be initialised in geodesics"<<std::endl;
        exit(1);
    }

    memcpy(state        , IC, 8*sizeof(double));
    memcpy(initial_state, IC, 8*sizeof(double));
    double Rmax = 1.05*state[1];

    const int max_steps = 1e6;
    int cc=0;
    for (;;)
    {
        // integrate the geodesic equation
        if (!solver->integrate(delta, state)) 
            return EndDomain::Error;

        // check end conditions
        double R = state[1];
        double theta = state[2];
        double z = R*cos(theta);

        // check if in black hole
        if (R<=1.05*Rs)
            return EndDomain::BlackHole;

        // check if hit disc
        if ((fabs(z)<0.05) && (disk_R0<R) && (R<disk_Rf))
            return EndDomain::Disk;

        // check if disapeared
        if (Rmax<R) 
            return EndDomain::Space;

        // error
        if (cc++>max_steps)
            return EndDomain::Error;
    
    }

}

void geodesics::rhs(const double* u, double* f)
{
    std::fill_n(f, 8, 0.);

    double u_[8];
    memcpy(&u_, u, 8*sizeof(double));

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

    // calcaulate values
    calc_metric(u_);
    calc_metric_inv(u_);
    calc_dmetric(u_);


    // calculate duds=f

    for (int mu=0; mu<4; mu++)
        f[mu] = u[mu+4];

    for (int mu=0; mu<4; mu++){

        for (int alpha=0; alpha<4; alpha++){
            for (int beta=0; beta<4; beta++){

                double lambda = dg[beta][mu][alpha] + dg[alpha][mu][beta] - dg[mu][alpha][beta]; // simplified since tensors are diagonal
                f[mu+4] += lambda * u[alpha+4] * u[beta+4];

            }
        }   
        f[mu+4] *= -0.5 * ginv[mu];

    }
}

void geodesics::calc_metric(const double (&state)[8])
{
    const double R = state[1];
    const double theta = state[2];
    const double f = (1-Rs/R);

    g[0] = -f;
    g[1] = 1/f;
    g[2] = R*R;
    g[3] = R*R*sin(theta)*sin(theta);  // where theta=0,pi this is zero!

    return;
}

void geodesics::calc_metric_inv(const double (&state)[8])
{
    for (int ii=0; ii<3; ii++)
        ginv[ii] = 1./g[ii];
    
    // handle possible divide by zero at theta=0,pi
    ginv[3] = (g[3] > 1e-6)?(1./g[3]):(1e6);

    return;
}

void geodesics::calc_dmetric(const double (&state)[8])
{
   const double R = state[1];
   const double theta = state[2];

    // dgdR
        double f = 1-Rs/R;
        dg[1][0][0] = -Rs/(R*R);
        dg[1][1][1] = -Rs / (R*R*f*f);
        dg[1][2][2] = 2.*R;
        dg[1][3][3]  = 2*R*sin(theta)*sin(theta);
       
    // dgdtheta
        dg[2][3][3] = R*R*sin(2.*theta);

    return;
}
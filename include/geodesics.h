#ifndef GEODESICS_H
#define GEODESICS_H

#include <cstring>

class integrator;

void f_geodesic(const double* u_ptr, double* f);

enum class EndDomain{
    NotSet,
    BlackHole,
    Disk,
    Space,
    Error
};

class geodesics
{

public:

    geodesics(integrator* const solver_, const double delta_);
    void init_Rs(const double Rs_);
    void init_annulus(const double annulus_inner_radius_, const double annulus_outer_radius_);
    void set_initial_conditions(const double delta_, const double (&u)[8]);
    void reset_state(void);
    ~geodesics(void){};

    void traverse_geodesic(const double (&initial_condition)[8], EndDomain &end, unsigned char (&rgb)[3]);


private:

    const double delta;

    const int max_integration_steps = 1e6;
    double annulus_inner_radius = 3.;
    double annulus_outer_radius = 6.;
    const double annulus_width_by_2 = 0.05;
    const unsigned char annulus_inner_rgb[3] = {255, 230, 200}; // colour at disk0
    const unsigned char annulus_outer_rgb[3] = {80, 30, 10};  // colour at diskf
    const double R0_prop_for_ray_escaped = 1.05;
    const double Rs_prop_for_in_black_hole = 1.05;
    double R_ray_in_black_hole = 1.05;
    double Rs = 1.;

    integrator* const solver = nullptr;


};

#endif
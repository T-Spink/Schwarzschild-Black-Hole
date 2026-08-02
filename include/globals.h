#ifndef GLOBALS_H
#define GLOBALS_H

#include "integrator.h"
#include <Eigen/Dense>
#include <numbers>
#include <atomic>

namespace globals
{
    // ======================================================================================= //
    // INPUTS

    // numerical integration
    constexpr double delta = 1e-1;
    constexpr IntegrationMethod int_method = IntegrationMethod::RK2;
    constexpr double rk45_tolerance = 1e-6;
    constexpr int rk45_max_steps = 1e6;

    // boundary
    constexpr double Rs = 1.;
    constexpr double annulus_inner_radius = 3. * Rs;
    constexpr double annulus_outer_radius = 6. * Rs;

    // observer
    static constexpr double observer_R = 9.;
    static constexpr double observer_theta = std::numbers::pi / 2. - 0.1;
    static constexpr double observer_phi = 0.;
    constexpr double fov_degs = 90.;
    
    // grahics
    constexpr int n_pixels = 300;
    constexpr int n_window = 500;    
    constexpr bool show_Rs = false;
    constexpr int N_pixels = n_pixels * n_pixels;


    // ======================================================================================= //

    // constants
    constexpr double rads_per_degs = std::numbers::pi / 180.;
    constexpr int n_3_vector = 3;  // [R, theta, phi]
    constexpr int n_4_vector = 4;  // [t, R, theta, phi]
    constexpr int n_8_vector = 8;  // [t, R, theta, phi, td, Rd, thetad, phid]
    constexpr int rgb_length = 3;  // [red, green, blue]
    constexpr int rgb_element_max = 255;
    constexpr unsigned char red_rgb[globals::rgb_length] = {rgb_element_max, 0, 0};
    constexpr unsigned char green_rgb[globals::rgb_length] = {0, rgb_element_max, 0};
    constexpr unsigned char black_rgb[globals::rgb_length] = {};

    // global values calculated in globals.cpp
    extern integrator* const solver;
    extern const double pixel_sf;
    extern const Eigen::Vector<double, n_3_vector> observer_position;
    extern const Eigen::Vector<double, n_3_vector> u;
    extern const Eigen::Vector<double, n_3_vector> v;
    extern const Eigen::Vector<double, n_3_vector> n;
    extern const Eigen::Matrix<double,globals::n_4_vector,globals::n_4_vector> tetrad;
    extern const Eigen::Matrix<double, n_3_vector, n_3_vector> cart_to_spherical;
    extern unsigned char* const pixels;
    extern std::atomic<int> render_count;


    // ======================================================================================= //
}

#endif
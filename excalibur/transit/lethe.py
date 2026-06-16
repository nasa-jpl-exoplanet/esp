'''grid generation for LETHE method'''

import numpy as np
from scipy import integrate
from scipy.interpolate import RectBivariateSpline
import pickle


def grid_generation(step_z_1=0.01, step_z_2=0.001, step_rprs=0.002):
    '''
    step_z_1 : Step for grid in separation from 0.0 to 0.8
    step_z_2 : Step for grid in separation from 0.8 to 1.2
    step_rprs : Step for grid in radius from 0.0 to 0.2
    '''

    # Parameters generation
    z_1 = np.arange(0.0, 0.8, step_z_1)
    z_2 = np.arange(0.8, 1.2 + step_z_2, step_z_2)
    z_grid = np.concatenate((z_1, z_2))
    rprs_grid = np.arange(0.0, 0.2 + step_rprs, step_rprs)

    # Parameters saving
    path = "/proj/sdp/data/LETHE"
    np.save(path + "/z_grid.npy", z_grid)
    np.save(path + "/rprs_grid.npy", rprs_grid)

    # First set of integrals
    def f1(r, theta, alpha, z):
        return (
            np.pow(
                np.max([0.0, 1.0 - z**2 - 2.0 * z * r * np.cos(theta) - r**2]),
                alpha,
            )
            * r
        )

    # Second set of integrals
    def f2(r, theta, alpha, z):
        x = 1.0 - z**2 - 2.0 * z * r * np.cos(theta) - r**2
        return f1(r, theta, alpha, z) * np.log(np.sign(x) * x) / 2.0

    # powers of the LD model
    alpha_test_1 = np.array([1.0 / 4.0, 2.0 / 4.0, 3.0 / 4.0, 4.0 / 4.0])
    alpha_test = np.concatenate((alpha_test_1, 2.0 * alpha_test_1))

    def compute_integral_1(z, rprs, alpha):
        res, err = integrate.dblquad(
            f1,
            0.0,
            2.0 * np.pi,  # bounds for theta
            0.0,
            rprs,  # bounds for r
            args=(alpha, z),
        )
        return res

    def compute_integral_2(z, rprs, alpha):
        res, err = integrate.dblquad(
            f2,
            0.0,
            2.0 * np.pi,  # bounds for theta
            0.0,
            rprs,  # bounds for r
            args=(alpha, z),
        )
        return res

    interp_names = [
        '/interpolator_G/f_G_0_25.pkl',
        '/interpolator_G/f_G_0_50.pkl',
        '/interpolator_G/f_G_0_75.pkl',
        '/interpolator_G/f_G_1_00.pkl',
        '/interpolator_F/f_F_0_50.pkl',
        '/interpolator_F/f_F_1_00.pkl',
        '/interpolator_F/f_F_1_50.pkl',
        '/interpolator_F/f_F_2_00.pkl',
    ]

    z, rprs = np.meshgrid(z_grid, rprs_grid, indexing='ij')
    # Computation of the values
    for index, name in enumerate(interp_names):
        if index <= 3:
            # Integrals computation
            vec_integral_1 = np.vectorize(
                lambda z_grid, rprs_grid: compute_integral_1(
                    z_grid, rprs_grid, alpha_test[index]
                )
            )
            integral = vec_integral_1(z, rprs)
            np.save(
                path + "/grid_G/" + str(alpha_test[index]) + ".npy", integral
            )

            # Interpolation computation
            interpolator = RectBivariateSpline(z_grid, rprs_grid, integral)
            with open(path + name, "wb") as file:
                pickle.dump(interpolator, file)
        else:
            # Integrals computation
            vec_integral_2 = np.vectorize(
                lambda z_grid, rprs_grid: compute_integral_2(
                    z_grid, rprs_grid, alpha_test[index]
                )
            )
            integral = vec_integral_2(z, rprs)
            np.save(
                path + "/grid_F/" + str(alpha_test[index]) + ".npy", integral
            )

            # Interpolation computation
            interpolator = RectBivariateSpline(z_grid, rprs_grid, integral)
            with open(path + name, "wb") as file:
                pickle.dump(interpolator, file)

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
    z_grid = np.concatenate(
        (
            np.arange(0.0, 0.8, step_z_1),
            np.arange(0.8, 1.2 + step_z_2, step_z_2),
        )
    )
    rprs_grid = np.arange(0.0, 0.2 + step_rprs, step_rprs)

    # powers of the LD model
    alpha_grid = np.array(
        [
            1.0 / 4.0,
            2.0 / 4.0,
            3.0 / 4.0,
            4.0 / 4.0,
            2.0 / 4.0,
            4.0 / 4.0,
            6.0 / 4.0,
            8.0 / 4.0,
        ]
    )

    # Parameters saving
    np.save("/proj/sdp/data/LETHE" + "/parameters/z_grid.npy", z_grid)
    np.save("/proj/sdp/data/LETHE" + "/parameters/rprs_grid.npy", rprs_grid)

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

    def compute_integral_1(z, rprs, alpha):
        res = integrate.dblquad(
            f1,
            0.0,
            2.0 * np.pi,  # bounds for theta
            0.0,
            rprs,  # bounds for r
            args=(alpha, z),
        )[0]
        return res

    vectorize_integral_1 = np.vectorize(compute_integral_1)

    def compute_integral_2(z, rprs, alpha):
        res = integrate.dblquad(
            f2,
            0.0,
            2.0 * np.pi,  # bounds for theta
            0.0,
            rprs,  # bounds for r
            args=(alpha, z),
        )[0]
        return res

    vectorize_integral_2 = np.vectorize(compute_integral_2)

    # integrals computation
    integral_1 = vectorize_integral_1(
        z_grid[:, None, None],
        rprs_grid[None, :, None],
        alpha_grid[None, None, :4],
    )

    integral_2 = vectorize_integral_2(
        z_grid[:, None, None],
        rprs_grid[None, :, None],
        alpha_grid[None, None, 4:],
    )

    integrals_list = [
        integral_1[:, :, 0],
        integral_1[:, :, 1],
        integral_1[:, :, 2],
        integral_1[:, :, 3],
        integral_2[:, :, 0],
        integral_2[:, :, 1],
        integral_2[:, :, 2],
        integral_2[:, :, 3],
    ]

    interpolation_names = [
        '/interpolator_G/f_G_0_25.pkl',
        '/interpolator_G/f_G_0_50.pkl',
        '/interpolator_G/f_G_0_75.pkl',
        '/interpolator_G/f_G_1_00.pkl',
        '/interpolator_F/f_F_0_50.pkl',
        '/interpolator_F/f_F_1_00.pkl',
        '/interpolator_F/f_F_1_50.pkl',
        '/interpolator_F/f_F_2_00.pkl',
    ]

    interpolators_list = []
    for integral in integrals_list:
        interpolator = RectBivariateSpline(z_grid, rprs_grid, integral)
        interpolators_list.append(interpolator)

    for name, interpolator in zip(interpolation_names, interpolators_list):
        with open("/proj/sdp/data/LETHE" + name, "wb") as file:
            pickle.dump(interpolator, file)

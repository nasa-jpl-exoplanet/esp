'''grid generation for LETHE method'''

import numpy as np
from scipy import integrate
from scipy.interpolate import RectBivariateSpline

PATH = "/proj/sdp/data/LETHE"
grid_names = [
    '/grid_G/0.25.npy',
    '/grid_G/0.5.npy',
    '/grid_G/0.75.npy',
    '/grid_G/1.0.npy',
    '/grid_F/0.5.npy',
    '/grid_F/1.0.npy',
    '/grid_F/1.5.npy',
    '/grid_F/2.0.npy',
]


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


# First integrals computation
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


# Second integrals computation
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
    np.save(PATH + "/parameters/z_grid.npy", z_grid)
    np.save(PATH + "/parameters/rprs_grid.npy", rprs_grid)

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

    for name, integral in zip(grid_names, integrals_list):
        np.save(PATH + name, integral)

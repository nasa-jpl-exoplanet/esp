import numpy as np 
import matplotlib.pyplot as plt
from scipy import integrate
import time

def grid_generation(z_samples=10000, rprs_samples=10000):
    
    debut = time.time()

    # Parameters generation
    zstep_2 = 0.39/(z_samples-1)
    zstep_1 = 10.0*zstep_2
    z_1 = np.arange(0.,0.85,zstep_1)
    z_2 = np.arange(0.85,1.15,zstep_2)
    z_3 = np.arange(1.15,1.2 + zstep_1, zstep_1)
    z_grid = np.concatenate((z_1, z_2, z_3))
    rprs_grid = np.linspace(0.,0.2,rprs_samples)
    
    # Parameters saving
    path = "/proj/sdp/data/LETHE/" + str(z_samples) + "_" + str(rprs_samples)
    np.save(path + "/z_grid.npy",z_grid)
    np.save(path + "/rprs_grid.npy",rprs_grid)
    
    # First set of integrals
    def f1(r, theta, alpha, z):
        return np.pow(np.max([0., 1.0 - z**2 - 2.0*z*r*np.cos(theta) - r**2]),alpha) * r
    
    # Second set of integrals
    def f2(r, theta, alpha, z):
        x = 1.0 - z**2 - 2.0*z*r*np.cos(theta) - r**2
        return f1(r, theta, alpha, z) * np.log(np.sign(x) * x) / 2.0
    
    # powers of the LD model
    alpha_test_1 = np.array([1.0/4.0, 2.0/4.0, 3.0/4.0, 4.0/4.0])
    alpha_test_2 = 2.0*alpha_test_1
    
    def compute_integral_1(z, rprs, alpha):
        res, err = integrate.dblquad(
            f1,
            0., 2.0*np.pi,    # bounds for theta
            0., rprs,           # bounds for r
            args=(alpha, z)
        )
        return res
    
    def compute_integral_2(z, rprs, alpha):
        res, err = integrate.dblquad(
            f2,
            0., 2.0*np.pi,    # bounds for theta
            0., rprs,           # bounds for r
            args=(alpha, z)
        )
        return res
    
    int = [1,2]
    
    for i in int :
        if (i == 1):
            for alpha in alpha_test_1:
                XR, RP = np.meshgrid(z_grid, rprs_grid, indexing='ij')
                vec_integral_1 = np.vectorize(lambda z_grid, rprs_grid: compute_integral_1(z_grid, rprs_grid, alpha))
                I = vec_integral_1(XR, RP)
                np.save(path + "/grid_G/" + str(alpha) + ".npy", I)
        else:
            for alpha in alpha_test_2:
                XR, RP = np.meshgrid(z_grid, rprs_grid, indexing='ij')
                vec_integral_2 = np.vectorize(lambda z_grid, rprs_grid: compute_integral_2(z_grid, rprs_grid, alpha))
                I = vec_integral_2(XR, RP)
                np.save(path + "/grid_F/" + str(alpha) + ".npy", I)
    
    fin = time.time()
    
    print(f"Temps d'exécution : {fin - debut:.6f} secondes")
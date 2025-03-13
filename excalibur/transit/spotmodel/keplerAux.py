import numpy as np


def keplerfunc(Marr, eccarr):
    """
    Solve Kepler's Equation using a Householder method.
    Marr, eccarr can be float or array, but we broadcast them so that
    indexing works properly.
    """
    # Convert Marr to a float array:
    Marr = np.array(Marr, dtype=float, ndmin=1)

    # If eccarr is just a single float or length-1 array,
    # broadcast it to the shape of Marr
    eccarr = np.array(eccarr, dtype=float, ndmin=1)
    if eccarr.size == 1 and Marr.size > 1:
        # broadcast a single e across all M
        eccarr = np.full_like(Marr, eccarr[0])

    conv = 1.0e-12  # convergence criterion
    k = 0.85

    # First guess
    Earr = Marr + np.sign(np.sin(Marr)) * k * eccarr
    fiarr = Earr - eccarr * np.sin(Earr) - Marr

    convd = np.where(np.abs(fiarr) > conv)[0]  # indices not converged
    nd = len(convd)
    count = 0

    # Householder iteration
    while nd > 0:
        count += 1

        #  Geoff: THIS IS NOT USED
        # M = Marr[convd]
        ecc = eccarr[convd]
        E = Earr[convd]
        fi = fiarr[convd]

        fip = 1.0 - ecc * np.cos(E)
        fipp = ecc * np.sin(E)
        fippp = 1.0 - fip

        # 3rd-order correction
        d1 = -fi / fip
        d2 = -fi / (fip + 0.5 * d1 * fipp)
        d3 = -fi / (fip + 0.5 * d2 * fipp + (d2**2) * fippp / 6.0)
        E = E + d3

        Earr[convd] = E
        fiarr = Earr - eccarr * np.sin(Earr) - Marr

        convd = np.where(np.abs(fiarr) > conv)[0]
        nd = len(convd)

        # Optional safety check to avoid infinite loops
        # if count > 5000:
        #     print("Warning: Kepler iteration not converging.")
        #     break

    # Return a scalar if there's only 1 element:
    if Earr.size == 1:
        return Earr[0]
    return Earr

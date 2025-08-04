import os
import numpy as np
from sympy import Symbol, solve
from excalibur.cerberus.tea_code import readconf as rc
from excalibur.cerberus.tea_code import format as form

def balance(a, b, verb=0, loc_out=None):
    """
    Generates an initial guess for TEA by satisfying the mass balance equation.

    Parameters
    ----------
    a : ndarray
        Stoichiometric coefficients of species (shape: [nspec, natom])
    b : ndarray
        Elemental abundances (shape: [natom])
    verb : int
        Verbosity level (0-2)
    loc_out : str or None
        Output location if saving (not used currently)

    Returns
    -------
    y : ndarray
        Initial guess for species abundances
    y_bar : float
        Sum of all species abundances
    """
    
    nspec, natom = a.shape
    verb = 0

    if verb > 1:
        print("b values:", b)
    # Find a block of species that covers all elements
    for n in range(nspec - natom + 1):
        a_chunk = a[n:n + natom, :]
        if np.all(np.sum(a_chunk, axis=0) != 0):
            free_id = list(range(n, n + natom))
            if verb > 1:
                print("Free variables selected:", free_id)
            break
    else:
        raise ValueError("Could not find a suitable block for free variables.")
    scale = 0.1
    nofit = True

    while nofit:
        pre_free = np.full(free_id[0], scale)
        post_free = np.full(nspec - free_id[-1] - 1, scale)

        free = [Symbol(f'y_unknown_{i}') for i in range(natom)]
        y_init = np.concatenate((pre_free, free, post_free))

        eqs = [sum(a[i, j] * y_init[i] for i in range(nspec)) - b[j] for j in range(natom)]
        result = solve(eqs, free, dict=True)
        
        if not result:
            scale /= 10
            if verb > 1:
                print(f"Trying smaller scale: {scale}")
            continue

        result = result[0]
        values = np.array([result[f] for f in free], dtype=float)

        if np.any(values < 0):
            scale /= 10
            if verb > 1:
                print(f"Negative values found. Trying smaller scale: {scale}")
        else:
            nofit = False
            if verb > 1:
                print(f"Viable scale found: {scale}")

    for i, idx in enumerate(free_id):
        y_init[idx] = values[i]

    if verb > 1:
        print('\nChecks:')
        for j in range(natom):
            lhs = round(np.sum(a[:, j] * y_init), 2)
            rhs = round(b[j], 2)
            if lhs == rhs:
                print(f'Equation {j + 1} satisfied.')
            else:
                print(f'Equation {j + 1} NOT satisfied.')

    y = np.array(y_init, dtype=np.float64)
    y_bar = y.sum()

    return y, y_bar
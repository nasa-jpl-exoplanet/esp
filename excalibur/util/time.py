'''function extracted from data and transit'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-locals,too-many-positional-arguments

import numpy as np


def time2z(time,
           ipct,
           tknot,
           sma,
           orbperiod,
           ecc,
           tperi=None,
           epsilon=1e-10,
           marray=True,
           ):
    '''
    G. ROUDIER: Time samples in [Days] to separation in [R*]
    '''
    if tperi is not None:
        ft0 = (tperi - tknot) % orbperiod
        ft0 /= orbperiod
        if ft0 > 0.5:
            ft0 += -1e0
            pass
        M0 = 2e0 * np.pi * ft0
        E0 = (solveme(np.array([M0]), ecc, epsilon) if marray
              else solveme(M0, ecc, epsilon)
              )
        realf = np.sqrt(1e0 - ecc) * np.cos(float(E0) / 2e0)
        imagf = np.sqrt(1e0 + ecc) * np.sin(float(E0) / 2e0)
        w = np.angle(np.complex(realf, imagf))
        if abs(ft0) < epsilon:
            w = np.pi / 2e0
            tperi = tknot
            pass
        pass
    else:
        w = np.pi / 2e0
        tperi = tknot
        pass
    ft = (time - tperi) % orbperiod
    ft /= orbperiod
    sft = ft
    sft[(sft > 0.5)] += -1e0
    M = 2e0 * np.pi * ft
    E = solveme(M, ecc, epsilon)
    realf = np.sqrt(1.0 - ecc) * np.cos(E / 2e0)
    imagf = np.sqrt(1.0 + ecc) * np.sin(E / 2e0)
    f = []
    for r, i in zip(realf, imagf):
        cn = complex(r, i)
        f.append(2e0 * np.angle(cn))
        pass
    f = np.array(f)
    r = sma * (1e0 - ecc**2) / (1e0 + ecc * np.cos(f))
    z = r * np.sqrt(
        1e0**2 - (np.sin(w + f) ** 2) * (np.sin(ipct * np.pi / 180e0)) ** 2
    )
    z[sft < 0] *= -1e0
    return z, sft


# --------------- ----------------------------------------------------
# -- TRUE ANOMALY NEWTON RAPHSON SOLVER -- ---------------------------
def solveme(M, e, eps):
    '''
    G. ROUDIER: Newton Raphson solver for true anomaly
    M is a numpy array
    '''
    E = M.copy()
    for i in np.arange(E.shape[0]):
        while abs(E[i] - e * np.sin(E[i]) - M[i]) > eps:
            num = E[i] - e * np.sin(E[i]) - M[i]
            den = 1.0 - e * np.cos(E[i])
            E[i] = E[i] - num / den
            pass
        pass
    return E

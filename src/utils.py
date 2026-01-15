import numpy as np

binary_comps = {
    "AlAs": {"Z": 23, "a": 5.6611, "Eg": 3.099, "dSO": 0.28},
    "InAs": {"Z": 41, "a": 6.0584, "Eg": 0.417, "dSO": 0.39},
    "AlP": {"Z": 14, "a": 5.4672, "Eg": 3.63, "dSO": 0.07},
    "InP": {"Z": 32, "a": 5.8697, "Eg": 1.4263, "dSO": 0.108},
    "GaSb": {"Z": 41, "a": 6.0959, "Eg": 0.812, "dSO": 0.76},
    "AlSb": {"Z": 32, "a": 6.1355, "Eg": 2.386, "dSO": 0.676},
    "InSb": {"Z": 50, "a": 6.4794, "Eg": 0.235, "dSO": 0.81},

    "GaAs": {"Z": 32, "a": 5.65325, "Eg": 1.519, "dSO": 0.341, 
             "m*e": 0.067, "alpha": 0.5405, "beta": 204,
             "gamma1": 6.98, "gamma2": 2.06, "gamma3": 2.93, 
             "c11": 1221, "c12": 566, "ac": -7.17, "av": 1.16,
             "VBO": -.8},
            
    "GaP": {"Z": 23, "a": 5.4505, "Eg": 2.886+0.1081*(1-1/np.tanh(164/300)), "dSO": 0.08, 
            "m*e": 0.13, "alpha": 0.5771, "beta": 372,
            "gamma1": 4.05, "gamma2": 0.49, "gamma3": 2.93, 
            "c11": 1405, "c12": 639, "ac": -8.2, "av": -1.7,
            "VBO": -1.27}
}

def linear_int(data, el1, el2, el3, x, bow=0, params=binary_comps):
        return x * params[f'{el1+el2}'][data] + (1-x) * params[f'{el1+el3}'][data] + x*(1-x)*bow
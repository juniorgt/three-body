from bodySytem import BodySystem

init_setup = {
    "name": "Figure-8",
    "ODESolver": "rk",
    "method": "euler",
    "G": 1,
    "M": [1, 1, 1],
    "y1": [-0.97000436, 0.4662036850, 0.24208753, 0.4323657300],
    "y2": [0.0, -0.933240737, 0.0, -0.86473146],
    "y3": [0.97000436, 0.4662036850, -0.24208753, 0.4323657300],
    "T": 6.3259,
    "h": 0.5,
}


BS = BodySystem(init_setup)
BS.run_simulation()
BS.plot()

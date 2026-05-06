import numpy as np
from pyswe import Settings, SWE

param = {"path_lib": "build/libwannier.dylib",
	 "path_tb": "hmcase0_tb.dat",
	 "nr1": 400,
     "nr2": 400,
     "nr3": 1,
     "tmax": 90.0,
     "dt": 21.97e-3,
     "intensity": 5e10,
     "lambda": 3000.0,
     "tmax_field": 80.0,
     "pol_vec": np.array([0.0, 1.0, 0.0]),
     "phi_vec": np.array([0.0, 0.0, 0.0])}

settings = Settings(param)

swe = SWE(settings)
swe.run_simulation()
t, jx, jy, jz = swe.get_current()
swe.delete()

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../pySWE")
from pyswe import Settings, SWE

param = {"path_lib": "../build/libwannier.so",
	 "path_tb": "../hmcase0_tb.dat",
	 "nr1": 400,
     "nr2": 400,
     "nr3": 1,
     "tmax": 90.0,
     "dt": 21.97e-3,
     "intensity": 1e12,
     "lambda": 3000.0,
     "tmax_field": 80.0,
     "pol_vec": np.array([0.0, 1.0, 0.0]),
     "phi_vec": np.array([0.0, 0.0, 0.0])}

settings = Settings(param)

swe = SWE(settings)
swe.run_simulation()
t, jx, jy, jz = swe.get_current()
print(jx)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(np.abs(np.fft.rfft(jx.real))**2 + np.abs(np.fft.rfft(jy.real))**2)
ax.set_yscale("log")
ax.set_ylim(1e-7, 1e3)
ax.set_xlim(-1, 500)
plt.savefig("test_swe.png")
swe.delete()

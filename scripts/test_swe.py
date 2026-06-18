import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../pySWE")
from pyswe import Settings, SWE

param = {"path_lib": "../build/libwannier.so",
	 "path_tb": "../hmcase0_tb.dat",
	 "nr1": 600,
     "nr2": 600,
     "nr3": 1,
     "tmax": 150.0,
     "dt": 21.97e-3/2.0,
     "intensity": 5e10,
     "lambda": 3000.0,
     "tmax_field": 80.0,
     "pol_vec": np.array([0.0, 1.0, 0.0]),
     "phi_vec": np.array([0.0, 0.0, 0.0])}

settings = Settings(param)

swe = SWE(settings)
swe.run_simulation()
t, jx, jy, jz = swe.get_current()
dt = t[1] - t[0]
print(jx)

tstart=param["tmax_field"]/param["dt"]*0.90
mask_j = np.ones(t.shape[0], dtype=np.float64)
mask_j = np.array([1 if ti<tstart else (np.exp(-(ti-tstart)**2/200**2))**1 for ti in t])
acc_x = np.gradient(mask_j*jx, dt)
acc_y = np.gradient(mask_j*jy, dt)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(np.abs(np.fft.rfft(acc_x.real))**2 + np.abs(np.fft.rfft(acc_y.real))**2)
ax.set_yscale("log")
#ax.set_ylim(1e-7, 1e3)
ax.set_xlim(-1,300)
plt.savefig("test_swe.png")
swe.delete()

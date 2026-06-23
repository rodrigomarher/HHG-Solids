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
     "tmax": 120.0,
     "nt": 8192,
     "intensity": 1e12,
     "lambda": 3000.0,
     "tmax_field": 120.0,
     "pol_vec": np.array([1.0, 0.0, 0.0]),
     "phi_vec": np.array([0.0, 0.0, 0.0])}

settings = Settings(param)

swe = SWE(settings)
swe.run_simulation()
t, jx, jy, jz = swe.get_current()
dt = t[1] - t[0]
print(jx)

tstart=param["tmax_field"]*param["nt"]/param["tmax"]*0.90
mask_j = np.ones(t.shape[0], dtype=np.float64)
mask_j = np.array([1 if ti<tstart else (np.exp(-(ti-tstart)**2/200**2))**1 for ti in range(t.shape[0])])
acc_x = np.gradient(mask_j*jx, dt)
acc_y = np.gradient(mask_j*jy, dt)

fig =plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(t, mask_j, color='k', alpha=0.7)
ax2 = ax.twinx()
ax.plot(t, jx)
ax.plot(t, jy)
plt.savefig("test_swe_mask.png")

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(np.abs(np.fft.rfft(acc_x.real))**2 + np.abs(np.fft.rfft(acc_y.real))**2)
ax.set_yscale("log")
ax.set_ylim(1e-8, 1e1)
ax.set_xlim(-1,400)
plt.savefig("test_swe.png")
swe.delete()

import numpy as np
import sys
sys.path.append("../pySWE")
from pyswe import Settings, SWE
import pickle as pckl
from mpi4py import MPI

param = {"path_lib": "../build/libwannier.so",
			 "path_tb": "../hmcase4_tb.dat",
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

mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()

def calculate(pol_angle):
	pol_vec = np.array([np.sin(pol_angle), np.cos(pol_angle), 0.0])
	param["pol_vec"] = pol_vec
	dt = param["dt"]
	print("python ", param["pol_vec"], flush=True)

	settings = Settings(param)
	swe = SWE(settings)
	swe.run_simulation()
	t, jx, jy, jz = swe.get_current()
	swe.delete()

	acc_x = np.gradient(jx, dt)
	acc_y = np.gradient(jy, dt)
	acc_z = np.gradient(jz, dt)
	t = np.arange(jy.shape[0])

	tstart=param["tmax_field"]/dt*0.95
	mask_j = np.ones(t.shape[0])
	mask_j = np.array([1 if ti<tstart else (np.exp(-(ti-tstart)**2/200**2))**1 for ti in t])

	hhg_x = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(mask_j*acc_x)))
	hhg_y = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(mask_j*acc_y)))
	hhg_z = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(mask_j*acc_z)))

	result  = {"acc_x": acc_x, "acc_y": acc_y, "acc_z": acc_z, "hhg_x": hhg_x, "hhg_y": hhg_y, "hhg_z": hhg_z}	
	return result

def main():
	if mpi_rank == 0:
		print("Starting pol_scan")	
	n_angles = 360
	pol_angles = np.linspace(0,2*np.pi, n_angles)
	nw = int(param["tmax"]/param["dt"])

	dt = param["dt"]
	dw = 2.0*np.pi/(param["tmax"])
	wmax = nw*dw/2.0
	freq = np.linspace(-wmax, wmax, nw)
	c = 300
	wl = param["lambda"]
	w0 = 2.0*np.pi*c/wl

	pol_idx = 0

	# Load balancing between mpi ranks
	batch_size = int(n_angles//mpi_size)
	batch_remainder = n_angles % mpi_size
	if mpi_rank < batch_remainder:
		idx_start = int(mpi_rank*(batch_size+1))
		idx_end = idx_start + batch_size + 1
	else:
		idx_start = batch_remainder * (batch_size + 1) + (mpi_rank - batch_remainder) * batch_size
		idx_end = idx_start + batch_size

	print(f"Rank: {mpi_rank}, {idx_start} - {idx_end}", flush=True)	

	local_data = []
	for idx in range(idx_start, idx_end):
		print(f"Rank: {mpi_rank}, idx: {idx}, pol_angle: {pol_angles[idx]}", flush=True)
		result = calculate(pol_angles[idx])	
		hhg_x = result["hhg_x"]
		hhg_y = result["hhg_y"]
		hhg_z = result["hhg_z"]
		acc_x = result["acc_x"]
		acc_y = result["acc_y"]
		acc_z = result["acc_z"]

		local_data.append([pol_angles[idx], hhg_x, hhg_y, hhg_z, acc_x, acc_y, acc_z])

	gathered = mpi_comm.gather(local_data, root=0)

	if mpi_rank == 0:
		print("Saving data")
		data_sim_hhg = np.zeros((n_angles, 3, nw), dtype=np.complex128)
		data_sim_acc = np.zeros((n_angles, 3, nw), dtype=np.complex128)
		pol_angles_gathered = []
		idx = 0
		for data in gathered:
			for row in data:
				pol_angles_gathered.append(row[0])
				data_sim_hhg[idx, 0, :] = row[1]
				data_sim_hhg[idx, 1, :] = row[2]
				data_sim_hhg[idx, 2, :] = row[3]
				data_sim_acc[idx, 0, :] = row[4]
				data_sim_acc[idx, 1, :] = row[5]
				data_sim_acc[idx, 2, :] = row[6]
				idx += 1

		data = {"pol_angles": pol_angles_gathered,
				"freq": freq/w0,
				"omega0": w0,
				"data_sim": data_sim_hhg,
				"data_acc": data_sim_acc}

		with open("hmcase4_I1e12_nk400.pckl", "wb") as f:
			pckl.dump(data, f)
if __name__=="__main__":
	main()

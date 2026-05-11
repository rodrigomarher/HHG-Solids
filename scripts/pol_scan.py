import numpy as np
from pyswe import Settings, SWE
import pickle as pckl
from mpi4py import MPI

param = {"path_lib": "build/libwannier.dylib",
			 "path_tb": "hmcase0_tb.dat",
			 "nr1": 200,
			 "nr2": 200,
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

	settings = Settings(param)
	swe = SWE(settings)
	swe.run_simulation()
	t, jx, jy, jz = swe.get_current()
	swe.delete()

	acc_x = np.gradient(jx, dt)
	acc_y = np.gradient(jy, dt)
	acc_z = np.gradient(jz, dt)

	hhg_x = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(acc_x)))
	hhg_y = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(acc_y)))
	hhg_z = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(acc_z)))
	
	return [hhg_x, hhg_y, hhg_z]

def main():
	if mpi_rank == 0:
		print("Starting pol_scan")	
	n_angles = 15
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

	print(f"Rank: {mpi_rank}, {idx_start} - {idx_end}")	

	local_data = []
	for idx in range(idx_start, idx_end):
		print(f"Rank: {mpi_rank}, idx: {idx}, pol_angle: {pol_angles[idx]}")
		hhg_x, hhg_y, hhg_z = calculate(pol_angles[idx])	
		
		local_data.append([pol_angles[idx], hhg_x, hhg_y, hhg_z])

	gathered = mpi_comm.gather(local_data, root=0)

	if mpi_rank == 0:
		print("Saving data")
		data_sim = np.zeros((n_angles, 3, nw), dtype=np.complex128)
		pol_angles_gathered = []
		idx = 0
		for data in gathered:
			for row in data:
				pol_angles_gathered.append(row[0])
				data_sim[idx, 0, :] = hhg_x
				data_sim[idx, 1, :] = hhg_y
				data_sim[idx, 2, :] = hhg_z
				idx += 1

		data = {"pol_angles": pol_angles_gathered,
				"freq": freq/w0,
				"omega0": w0,
				"data_sim": data_sim}

		with open("test.pckl", "wb") as f:
			pckl.dump(data, f)
if __name__=="__main__":
	main()

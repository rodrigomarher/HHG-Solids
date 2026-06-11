#!/usr/bin/env python3
import numpy as np

#filepath = '../SimData/GrSymTB_nk700_PolScan.tar'
filepath = 'hmcase0_I1e12_nk400.pckl'
figname = "hmcase0_carlos_"
#q = np.arange(3,29,1)
q_max = 70
P = 1
import numpy as np
import matplotlib.pyplot as plt
import tarfile
#from numba import njit, jit, prange
from scipy.interpolate import splrep, interp1d
from prop_q.prop_q_ctypes import prop_q
import pickle as pckl


# Custom wrap of scipy's splrep
def custom_splrep(x, y, k=3):
    from scipy import interpolate
    """
    Custom wrap of scipy's splrep for calculating spline coefficients, 
    which also check if the data is equispaced.
    
    """
    
    # Check if x is equispaced
    x_diff = np.diff(x)
    equi_spaced = all(np.round(x_diff,5) == np.round(x_diff[0],5))
    dx = x_diff[0]
    
    # Calculate knots & coefficients (cubic spline by default)
    t,c,k = interpolate.splrep(x,y, k=k) 
    
    return (t,c,k,equi_spaced,dx) 

# Numba accelerated implementation of scipy's splev
#@njit(cache=True)
def numba_splev(x, coeff): # spline is extrapolated from the end spans for points not in the support.

  t,c,k = coeff

  n = t.size
  m = x.size

  k1  = k+1
  k2  = k1+1
  nk1 = n - k1

  l  = k1
  l1 = l+1

  y = np.zeros(m)

  for i in range(m):

    # fetch a new x-value arg
    arg = x[i]

    # search for knot interval t[l] <= arg <= t[l+1]
    while not ((arg >= t[l-1]) or (l1 == k2)):
      l1 = l
      l  = l-1
    while not ((arg < t[l1-1]) or (l == nk1)):
      l = l1
      l1 = l+1
    
    # evaluate the non-zero b-splines at arg.    
    h  = np.zeros(20)
    hh = np.zeros(19)
    
    h[0] = 1.0
    
    for j in range(k):

      for ll in range(j+1):

        hh[ll] = h[ll]
      h[0] = 0.0
      for ll in range(j+1):
        li = l + ll 
        lj = li - j - 1
        if(t[li] != t[lj]):
          f = hh[ll]/(t[li]-t[lj])
          h[ll] += f*(t[li]-arg)
          h[ll+1] = f*(arg-t[lj])
        else:
          h[ll+1] = 0.0
          break
    sp = 0.0
    ll = l - 1 - k1
    for j in range(k1):
      ll += 1
      sp += c[ll]*h[j]
    y[i] = sp
  return y


#numba_splev(np.array([4.0], dtype=np.complex128), hhg_q_rcp_coefs)
#@njit(parallel=True)
#def prop_q(coefs_interp, angles, freq,q, theta=None, omega=None):
#    phi = np.linspace(0,2*np.pi, 128)
#    #if theta == None:
#    #    theta = np.linspace(0/1000, 40/1000, 128)
#    #if omega == None:
#    #    omega = np.linspace(0,2*np.pi, 256)
#    rho_array = np.arange(22,35,0.5)
#    k = q*2*np.pi/3
#    far_field_q = np.zeros((theta.shape[0], omega.shape[0]), dtype=np.complex128)
#    
#    for i in range(theta.shape[0]):
#        print(i,theta[i])
#        for j in range(omega.shape[0]):
#            for phii in phi:
#                for rho_i in rho_array:
#                    tmp_real = numba_splev(np.array([phii]), coefs_interp[0])[0]
#                    tmp_imag = numba_splev(np.array([phii]), coefs_interp[1])[0]
#                    tmp = tmp_real + 1j*tmp_imag
#                    far_field_q[i,j] += rho_i*tmp*np.exp(-1j*k*theta[i]*rho_i*np.cos(omega[j]-phii)) 
#    return far_field_q

def plot_data(path, q):
    with open(path, "rb") as f:
        simdata = pckl.load(f)

    

    delta_q = 0.75
    pol_angles = np.array(simdata["pol_angles"])
    freq = np.array(simdata["freq"])
    
    hhg_x = np.array([simdata["data_sim"][i,0,:] for i in range(pol_angles.shape[0])])
    hhg_y = np.array([simdata["data_sim"][i,1,:] for i in range(pol_angles.shape[0])])
    
    mask = np.zeros(freq.shape)
    idx_q_min = np.argmin(np.abs(freq - q + delta_q))
    idx_q_max = np.argmin(np.abs(freq - q - delta_q))
    idx_q = np.argmin(np.abs(freq - q))
    mask[idx_q_min:idx_q_max] =1.0
    
    #hhg_q_x = mask.reshape(1,-1)*hhg_x
    #hhg_q_y = mask.reshape(1,-1)*hhg_y
    hhg_q_x = hhg_x[:,idx_q]
    hhg_q_y = hhg_y[:,idx_q]
    
    hhg_q_rcp = 1/np.sqrt(2)*(hhg_q_x + 1j*hhg_q_y)
    hhg_q_lcp = 1/np.sqrt(2)*(hhg_q_x - 1j*hhg_q_y)
    
    #hhg_q_rcp = np.sum(hhg_q_rcp, axis=1)
    #hhg_q_lcp = np.sum(hhg_q_lcp, axis=1)
    
    
    #hhg_q_x = np.sum(hhg_q_x, axis=1)
    #hhg_q_y = np.sum(hhg_q_y, axis=1)
    
    
    
    hhg_q_rcp_interp = interp1d(pol_angles, hhg_q_rcp, kind='cubic')
    hhg_q_lcp_interp = interp1d(pol_angles, hhg_q_lcp,  kind='cubic')
    
    hhg_q_x_interp = interp1d(pol_angles, hhg_q_x, kind='cubic')
    hhg_q_y_interp = interp1d(pol_angles, hhg_q_y,  kind='cubic')
    
    
    
    new_angles = np.linspace(0,2*np.pi, 1024)
    
    hhg_rcp = hhg_q_rcp_interp(P*new_angles)#*np.exp(-1j*1*new_angles)
    hhg_rcp /= np.abs(hhg_rcp).max()
    hhg_lcp = hhg_q_lcp_interp(P*new_angles)#*np.exp(1j*1*new_angles)
    hhg_lcp /= np.abs(hhg_lcp).max()
    
    l_arr = np.arange(-15, 15, 1)
    dnew_angles = new_angles[1] - new_angles[0]
    cl_rcp = np.array([np.sum((hhg_rcp)*np.exp(-1j*(i)*new_angles)*dnew_angles) for i in l_arr])/(2*np.pi)
    cl_lcp = np.array([np.sum((hhg_lcp)*np.exp(-1j*(i)*new_angles)*dnew_angles) for i in l_arr])/(2*np.pi)
    
    idx_rcp = np.argmin(np.abs(l_arr+1))
    idx_lcp = np.argmin(np.abs(l_arr-1))
    # cl_rcp /= np.abs(cl_rcp[idx_rcp])
    # cl_lcp /= np.abs(cl_lcp[idx_lcp])
    
    idx_rcp_7 = np.argmin(np.abs(l_arr+7))
    idx_rcp_0 = np.argmin(np.abs(l_arr+1))
    idx_rcp_5 = np.argmin(np.abs(l_arr-5))
    
    idx_lcp_7 = np.argmin(np.abs(l_arr-7))
    idx_lcp_0 = np.argmin(np.abs(l_arr-1))
    idx_lcp_5 = np.argmin(np.abs(l_arr+5))
    
    c_rcp_neg = cl_rcp[idx_rcp_5]/cl_rcp[idx_rcp_0]
    c_rcp_plus = cl_rcp[idx_rcp_7]/cl_rcp[idx_rcp_0]
    
    c_lcp_neg = cl_lcp[idx_lcp_7]/cl_lcp[idx_lcp_0]
    c_lcp_plus = cl_lcp[idx_lcp_5]/cl_lcp[idx_lcp_0]
    
    
    theta = np.linspace(0/1000, 30/1000, 256)
    omega = np.linspace(0,2*np.pi, 512)
    Theta, Omega = np.meshgrid(theta, omega)
    X = Theta*np.cos(Omega)
    Y = Theta*np.sin(Omega)
    
    
    hhg_q_rcp_coefs_real = splrep(pol_angles, hhg_q_rcp.real, k=3)
    hhg_q_lcp_coefs_real = splrep(pol_angles, hhg_q_lcp.real, k=3)
    hhg_q_rcp_coefs_imag = splrep(pol_angles, hhg_q_rcp.imag, k=3)
    hhg_q_lcp_coefs_imag = splrep(pol_angles, hhg_q_lcp.imag, k=3)
    
    #hhg_q_rcp_coefs = (np.deg2rad(pol_angles),hhg_q_rcp_coefs._spline.c,3,True, np.deg2rad(pol_angles)[1]-np.deg2rad(pol_angles)[0])
    #hhg_q_lcp_coefs = (np.deg2rad(pol_angles),hhg_q_lcp_coefs._spline.c,3,True, np.deg2rad(pol_angles)[1]-np.deg2rad(pol_angles)[0])
    
    #far_field_q_rcp = prop_q([hhg_q_rcp_coefs_real, hhg_q_rcp_coefs_imag], pol_angles, freq, q, theta=theta, omega=omega)
    #far_field_q_rcp /= np.abs(far_field_q_rcp).max()
    coefs_interp = (hhg_q_rcp_coefs_real, hhg_q_rcp_coefs_imag)
    far_field_q_rcp = prop_q(coefs_interp, pol_angles, freq, q, theta=theta, omega=omega)
    #far_field_q_rcp /= np.abs(far_field_q_rcp).max()

    coefs_interp = (hhg_q_lcp_coefs_real, hhg_q_lcp_coefs_imag)
    far_field_q_lcp = prop_q(coefs_interp, pol_angles, freq, q, theta=theta, omega=omega)
    #far_field_q_lcp /= np.abs(far_field_q_lcp).max()
    #
    #far_field_q_lcp = prop_q([hhg_q_lcp_coefs_real, hhg_q_lcp_coefs_imag], pol_angles, freq, q, theta=theta, omega=omega)
    #far_field_q_lcp /= np.abs(far_field_q_lcp).max()
    
    print(f"q: {q}")
    print("\t c_rcp_neg: ",np.abs(c_rcp_neg), np.angle(c_rcp_neg)/np.pi, "$\\pi$")
    print("\t c_rcp_pos: ",np.abs(c_rcp_plus), np.angle(c_rcp_plus)/np.pi, "$\\pi$")
    print("\t c_lcp_neg: ",np.abs(c_lcp_neg), np.angle(c_rcp_plus)/np.pi,"$\\pi$")
    print("\t c_rcp_pos: ",np.abs(c_lcp_plus), np.angle(c_lcp_plus)/np.pi,"$\\pi$")

    print("E_max (lcp): ",np.abs(far_field_q_lcp).max(), " E_max (rcp): ", np.abs(far_field_q_rcp).max())
    max_value = np.abs(far_field_q_lcp).max()**2 if np.abs(far_field_q_lcp).max() > np.abs(far_field_q_rcp).max() else np.abs(far_field_q_rcp
).max()**2
    
    fig = plt.figure(figsize=(12.5,4),constrained_layout=True)
    ax = fig.add_subplot(1,3,1)
    ax.pcolormesh(X*1e3,Y*1e3, np.abs(far_field_q_rcp.T)**2, vmax = max_value*1.03, cmap='turbo', shading='gouraud')
    ax.set_xlabel("X divergence (mrad)", fontsize=15)
    ax.set_ylabel("Y divergence (mrad)", fontsize=15)
    
    # Move left and bottom spines outward by 10 points
    ax.spines.left.set_position(('outward', 10))
    ax.spines.bottom.set_position(('outward', 10))
    # Hide the right and top spines
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(labelsize=15)
    ax.set_xticks([-20,0,20])
    ax.set_yticks([-20,0,20])
    ax.set_title("Right polarization", fontsize=15)
    
    axins = ax.inset_axes(
        [0.7, 0.00, 0.35, 0.35],
        xticks=[], yticks=[])
    
    mask_intensity = np.ones(far_field_q_rcp.shape)*0.6
    mask_intensity[(np.abs(far_field_q_rcp)/np.abs(far_field_q_rcp).max())**2>0.05] = 0.
    
    axins.pcolormesh(X*1e3,Y*1e3, np.angle(far_field_q_rcp.T)**1, vmin=-np.pi, vmax=np.pi,cmap='hsv', shading='gouraud', )
    axins.pcolormesh(X*1e3,Y*1e3, np.angle(far_field_q_rcp.T)**1, cmap='grey', shading='gouraud', alpha=mask_intensity.T)
    #axins.set_xlim(-20,20)
    #axins.set_ylim(-20,20)
    #axins.set_xlabel("Time", labelpad=7, fontsize=12)
    #axins.set_ylabel("X axis", labelpad=3, fontsize=12)
    #axins.set_title("Driving field", fontsize=12,)
    axins.axis('off')
    
    
    ax = fig.add_subplot(1,3,2)
    ax.pcolormesh(X*1e3,Y*1e3, np.abs(far_field_q_lcp.T)**2, vmax = max_value*1.03, cmap='turbo', shading='gouraud')
    ax.set_xlabel("X divergence (mrad)", fontsize=15)
    ax.set_ylabel("Y divergence (mrad)", fontsize=15)
    
    # Move left and bottom spines outward by 10 points
    ax.spines.left.set_position(('outward', 10))
    ax.spines.bottom.set_position(('outward', 10))
    # Hide the right and top spines
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(labelsize=15)
    ax.set_xticks([-20,0,20])
    ax.set_yticks([-20,0,20])
    ax.set_title("Left polarization", fontsize=15)
    
    axins = ax.inset_axes(
        [0.7, 0.00, 0.35, 0.35],
        xticks=[], yticks=[])
    
    mask_intensity = np.ones(far_field_q_lcp.shape)*0.6
    mask_intensity[(np.abs(far_field_q_lcp)/np.abs(far_field_q_lcp).max())**2>0.05] = 0.
    
    axins.pcolormesh(X*1e3,Y*1e3, np.angle(far_field_q_lcp.T)**1, vmin=-np.pi, vmax=np.pi,cmap='hsv', shading='gouraud', )
    axins.pcolormesh(X*1e3,Y*1e3, np.angle(far_field_q_lcp.T)**1, cmap='grey', shading='gouraud', alpha=mask_intensity.T)
    #axins.set_xlim(-20,20)
    #axins.set_ylim(-20,20)
    #axins.set_xlabel("Time", labelpad=7, fontsize=12)
    #axins.set_ylabel("X axis", labelpad=3, fontsize=12)
    #axins.set_title("Driving field", fontsize=12,)
    axins.axis('off')
    
    ax = fig.add_subplot(1,3,3)
    
    ax.scatter(l_arr,np.abs(cl_rcp)**1, label='RCP', s=80)
    ax.scatter(l_arr, np.abs(cl_lcp)**1, label='LCP', s=80)
    ax.set_ylabel("$|c_{\\ell}^q|$", fontsize=15)
    ax.set_xlabel("$\\ell$ value", fontsize=15)
    ax.tick_params(labelsize=15)
    ax.legend(loc=5, fontsize=15)
    
    str_title = "$q=" + str(int(q))+"$"
    fig.text(0.01,0.955, str_title, fontsize=15)
    
    name_fig = figname +"_"+f"q{int(q)}_new.png"
    print(name_fig)
    ax.set_rasterized(True)
    ax.set_xlim(-10,10)
    plt.savefig(name_fig, dpi=600)
 
    plt.clf()
    plt.close(fig)
    
def main():
    [plot_data(filepath, qi) for qi in q]

if __name__=='__main__':
    main()

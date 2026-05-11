import ctypes as ct
import numpy as np

ND_POINTER_1_int32 = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C") 
ND_POINTER_1_double = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C") 
ND_POINTER_1_float = np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags="C") 
ND_POINTER_1_complex = np.ctypeslib.ndpointer(dtype=np.complex128, ndim=1, flags="C")

class SWE:
    def __init__(self, settings):
        self.settings = settings
        self.swe_cpp = ct.CDLL(settings.path_lib)
        self.swe_cpp.SWESim_new.argtypes = []
        self.swe_cpp.SWESim_new.restype = ct.c_void_p
        self.swe_cpp.SWESim_run_simulation.argtypes = [ct.c_void_p]
        self.swe_cpp.SWESim_run_simulation.restype = ct.c_void_p
        self.swe_cpp.SWESim_test_files.argtypes = [ct.c_void_p]
        self.swe_cpp.SWESim_test_files.restypes = ct.c_void_p
        self.swe_cpp.SWESim_set_path_tb.argtypes = [ct.c_void_p, ct.c_char_p]
        self.swe_cpp.SWESim_set_path_tb.restype = ct.c_void_p
        self.swe_cpp.SWESim_set_settings.argtypes = [ct.c_void_p, ct.c_void_p]
        self.swe_cpp.SWESim_set_settings.restype = ct.c_void_p
        self.swe_cpp.SWESim_init.argtypes = [ct.c_void_p]
        self.swe_cpp.SWESim_init.restype = ct.c_void_p
        self.swe_cpp.SWESim_restart.argtypes = [ct.c_void_p]
        self.swe_cpp.SWESim_restart.restype = ct.c_void_p
        self.swe_cpp.SWESim_get_current.argtypes = [ct.c_void_p, ND_POINTER_1_double, ND_POINTER_1_complex, ND_POINTER_1_complex, ND_POINTER_1_complex]
        self.swe_cpp.SWESim_get_current.restype = ct.c_void_p
        self.swe_cpp.SWESim_delete.argtypes = [ct.c_void_p]
        self.swe_cpp.SWESim_delete.restype = ct.c_void_p
        
        self.ptr = self.swe_cpp.SWESim_new()
        self.set_path_tb(settings.path_tb.encode("UTF-8"))
        self.set_settings(settings.ptr)
        self.init()
 
    def run_simulation(self):
        self.swe_cpp.SWESim_run_simulation(self.ptr)
    def test_files(self):
        self.swe_cpp.SWESim_test_files(self.ptr)
    def set_path_tb(self, val):
        self.swe_cpp.SWESim_set_path_tb(self.ptr, val)
    def set_settings(self, val):
        self.swe_cpp.SWESim_set_settings(self.ptr, val)
    def init(self):
        self.swe_cpp.SWESim_init(self.ptr)
    def restart(self):
        self.swe_cpp.SWESim_restart(self.ptr)
    def delete(self):
        self.swe_cpp.SWESim_delete(self.ptr)
    def get_current(self):
        t  = np.zeros(self.settings.nt, dtype=np.float64)
        jx = np.zeros(self.settings.nt, dtype=np.complex128)
        jy = np.zeros(self.settings.nt, dtype=np.complex128)
        jz = np.zeros(self.settings.nt, dtype=np.complex128)
        self.swe_cpp.SWESim_get_current(self.ptr, t, jx, jy, jz) 
        return np.array([t, jx, jy, jz])

class Settings:
    def __init__(self, param):
        self.path_lib = param["path_lib"]
        self.path_tb = param["path_tb"]
        self.nt = int(param["tmax"]/param["dt"])
        self.settings_swe = ct.CDLL(self.path_lib)
        self.settings_swe.Settings_new.argtypes = []
        self.settings_swe.Settings_new.restype = ct.c_void_p
        self.settings_swe.Settings_set_nr1.argtypes = [ct.c_void_p, ct.c_int]
        self.settings_swe.Settings_set_nr1.restype = ct.c_void_p
        self.settings_swe.Settings_set_nr2.argtypes = [ct.c_void_p, ct.c_int]
        self.settings_swe.Settings_set_nr2.restype = ct.c_void_p
        self.settings_swe.Settings_set_nr3.argtypes = [ct.c_void_p, ct.c_int]
        self.settings_swe.Settings_set_nr3.restype = ct.c_void_p
        self.settings_swe.Settings_set_tmax.argtypes = [ct.c_void_p, ct.c_double]
        self.settings_swe.Settings_set_tmax.restype = ct.c_void_p
        self.settings_swe.Settings_set_dt.argtypes = [ct.c_void_p, ct.c_double]
        self.settings_swe.Settings_set_dt.restype = ct.c_void_p
        self.settings_swe.Settings_set_intensity.argtypes = [ct.c_void_p, ct.c_double]
        self.settings_swe.Settings_set_intensity.restype = ct.c_void_p
        self.settings_swe.Settings_set_lambda.argtypes = [ct.c_void_p, ct.c_double]
        self.settings_swe.Settings_set_lambda.restype = ct.c_void_p
        self.settings_swe.Settings_set_tmax_field.argtypes = [ct.c_void_p, ct.c_double]
        self.settings_swe.Settings_set_tmax_field.restype = ct.c_void_p
        self.settings_swe.Settings_set_pol_vec.argtypes = [ct.c_void_p, ND_POINTER_1_double]
        self.settings_swe.Settings_set_pol_vec.restype = ct.c_void_p
        self.settings_swe.Settings_set_phi_vec.argtypes = [ct.c_void_p, ND_POINTER_1_double]
        self.settings_swe.Settings_set_phi_vec.restype = ct.c_void_p

        self.ptr = self.settings_swe.Settings_new()    
        self.settings_swe.Settings_set_nr1(self.ptr, param["nr1"])
        self.settings_swe.Settings_set_nr2(self.ptr, param["nr2"])
        self.settings_swe.Settings_set_nr3(self.ptr, param["nr3"])
        self.settings_swe.Settings_set_tmax(self.ptr, param["tmax"])
        self.settings_swe.Settings_set_dt(self.ptr, param["dt"])
        self.settings_swe.Settings_set_intensity(self.ptr, param["intensity"])
        self.settings_swe.Settings_set_lambda(self.ptr, param["lambda"])
        self.settings_swe.Settings_set_tmax_field(self.ptr, param["tmax_field"])
        self.settings_swe.Settings_set_pol_vec(self.ptr, param["pol_vec"])
        self.settings_swe.Settings_set_phi_vec(self.ptr, param["phi_vec"])

import numpy as np
import pandas as pd

class SetLayer:
    def __init__(self, L, dx, k, c, rho):
        self.L = L  # Length [m]
        self.dx = dx  # dx [m]
        self.div = round(L/dx)  # Number of division [-] 
        self.k = k  # Thermal conductivity [W/mK]
        self.C = c*rho  # Volumetric heat capacity [J/m^3K]
        self.R = self.dx/k  # Thermal resistance [m^2K/W]
        self.K = k/self.dx  # Thermal conductance # [W/m^2K]
        self.alpha = k/self.C  # Thermal diffusivity [m^2/s]

class SetConstruction:
    def __init__(self, *layers):
        self.div_lst = [Lidx.div for Lidx in layers]  # Division list
        self.layer_lst = list(layers)  # Layer list
        self.lN = len(self.layer_lst)  # Number of layers [-]
        self.N = sum(self.div_lst)  # Number of nodes [-]

    def dx(self):
        dx = np.array([self.layer_lst[Lidx].dx for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # dx [m]
        return dx
    
    def dx_L(self):
        dx = np.array([self.layer_lst[Lidx].dx for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # dx [m]
        dx_L = np.array([dx[0]/2] + [(dx[i-1] + dx[i])/2 for i in range(1, self.N)])  
        return dx_L
    
    def dx_R(self):
        dx = np.array([self.layer_lst[Lidx].dx for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # dx [m]
        dx_R = np.array([(dx[i] + dx[i+1])/2 for i in range(self.N-1)] + [dx[self.N-1]/2]) 
        return dx_R

    def K(self):
        R = np.array([self.layer_lst[Lidx].R for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # Thermal resistance [m^2K/W]
        K = 1/R  # Thermal conductance [W/m^2K]  
        return K
    
    def K_L(self):
        R = np.array([self.layer_lst[Lidx].R for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # Thermal resistance [m^2K/W]
        K_L = 1/np.array([R[0]/2] + [(R[i-1] + R[i])/2 for i in range(1, self.N)])  # Left interface thermal conductance [W/m^2K]
        return K_L
    
    def K_R(self):
        R = np.array([self.layer_lst[Lidx].R for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # Thermal resistance [m^2K/W]
        K_R = 1/np.array([(R[i] + R[i+1])/2 for i in range(self.N-1)] + [R[self.N-1]/2])  # Right interface thermal conductance [W/m^2K]
        return K_R
    
    def R(self):
        R = np.array([self.layer_lst[Lidx].R for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # Thermal resistance [m^2K/W]
        return R
    
    def R_L(self):
        R = np.array([self.layer_lst[Lidx].R for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # Thermal resistance [m^2K/W]
        R_L = np.array([R[0]/2] + [(R[i-1] + R[i])/2 for i in range(1, self.N)])  # Left interface thermal resistance [m^2K/W]
        return R_L
    
    def R_R(self):
        R = np.array([self.layer_lst[Lidx].R for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # Thermal resistance [m^2K/W]
        R_R = np.array([(R[i] + R[i+1])/2 for i in range(self.N-1)] + [R[self.N-1]/2])  # Right interface thermal resistance [m^2K/W]
        return R_R
    
    def K_tot(self):
        R = np.array([self.layer_lst[Lidx].R for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # Thermal resistance [m^2K/W]
        R_tot = np.sum(R)  # Total thermal resistance [m^2K/W]
        K_tot = 1/R_tot
        return K_tot
    
    def R_tot(self):
        R = np.array([self.layer_lst[Lidx].R for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # Thermal resistance [m^2K/W]
        R_tot = np.sum(R)  # Total thermal resistance [m^2K/W]
        return R_tot
    
    def C(self):
        C = np.array([self.layer_lst[Lidx].C for Lidx in range(self.lN) for _ in range(self.div_lst[Lidx])])  # Volumetric heat capacity [J/m^3K]
        return C

def D2K(Degree):  # Degree to Kelvin
    return Degree + 273.15

def K2D(Kelvin):  # Kelvin to Degree
    return Kelvin - 273.15

def arr2df(arr):  # Array to DataFrame
    return pd.DataFrame(arr)

def df2arr(df):
    return df.values

def cm2in(cm):  # cm to inch
    return cm / 2.54

def half_time_vector(vec):
    return (vec[:-1] + vec[1:]) / 2

def half_time_matrix(mat, axis):
    if axis == 0:
        return (mat[:-1, :] + mat[1:, :]) / 2
    elif axis == 1:
        return (mat[:, :-1] + mat[:, 1:]) / 2

def eliminate_pre_simulation_data(data, pre_simulation_time_step):
    return data[pre_simulation_time_step:, :]

def TDMA(Construction, T, T_L, T_R, t):
    N = Construction.N
    dx, = Construction.dx()
    K_L, K_R = Construction.K_L(), Construction.K_R()
    C = Construction.C()
    a_list = [-t*K_L[i] for i in range(N)]
    b_list = [2*dx[i]*C[i]+t*K_L[i]+t*K_R[i] for i in range(N)]
    c_list = [-t*K_R[i] for i in range(N)]

    # Set A matrix
    A = np.zeros((N, N))
    np.fill_diagonal(A[1:], a_list[1:])
    np.fill_diagonal(A[:], b_list[:])
    np.fill_diagonal(A[:, 1:], c_list[:-1])
    A_inv = np.linalg.inv(A)

    # Set B matrix
    B = []
    for i in range(1, N-1):
        B.append([t*K_L[i]*T[i-1]+(2*dx[i]*C[i]-t*K_L[i]-t*K_R[i])*T[i]+t*K_R[i]*T[i+1]])
    B.insert(0, [2*t*K_L[0]*T_L[0]+(2*dx[0]*C[0]-t*K_L[0]-t*K_R[0])*T[0]+t*K_R[0]*T[1]])
    B.append([t*K_L[N-1]*T[N-2]+(2*dx[N-1]*C[N-1]-t*K_L[N-1]-t*K_R[N-1])*T[N-1]+2*t*K_R[N-1]*T_R[N-1]])
    
    # Calculate next time step temperature
    T_new = [np.dot(A_inv, B)[i,0] for i in range(N)]
    return np.array(T_new)

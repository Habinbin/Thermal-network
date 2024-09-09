import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from tqdm import tqdm
from construction import SetLayer, SetConstruction
import dartwork_mpl as dm
import thermal_network as tn

dm.use_dmpl_style()

### Constant & variable ---------------------------------------------------------------

## convective variable
h_co = 8 # convective heat transfer coefficient [W/m2K]
R_co = 1/h_co # convective thermal resistance [m2K/W]

## Time variable
dt = 20 #[s] 
PST = 120 #PS:Pre Simulation[h]
TST = 144 #TS:Total Simulation[h]
MST = TST - PST #MS:Main Simulation[h]

tN = int(TST*tn.h2s/dt+2) # Number of time steps +2 (120.0 - 144.0 + dt [h]) 
ts_s = np.array([dt*i for i in range(tN)]) # time step list [t1,t2,..., tN]
ts_m = ts_s*tn.s2m; #minute
ts_h = ts_s*tn.s2h; #hour

ts_PST = int(PST*tn.h2s/dt) # Number of pre-simulation time steps
ts_TST = int(TST*tn.h2s/dt)+1 # Number of total simulation time steps
 


### Set wall ---------------------------------------------------------------

## Set layer
CON = tn.SetLayer(L=0.01, dx=0.01, k=1.4, c= 1000, rho=2200)
INS = tn.SetLayer(L=0.01, dx=0.01, k=0.03,c= 1400, rho=  25)

## Set wall
ConThick = 15
InsThick = 5
TotalThick = ConThick + InsThick
pushcm = 5
pushnum = int(1+ConThick/pushcm)
WallList = [[CON]*(ConThick-(i*pushcm)) + [INS]*InsThick + [CON]*(i*pushcm) for i in range(pushnum)] # 단열재를 5 cm씩 밖으로 밀면서 새로운 벽 생성
WallIndexList = [[(ConThick-(i*pushcm))] + [InsThick] + [(i*pushcm)] for i in range(pushnum)] # [ExteriorConThick, InsulationThick, InteriorConThick] index list



### Name list ---------------------------------------------------------------
WallNameList = ["W"+f"{i+1}" for i in range(pushnum)] # "W1", "W2", "W3", ...
x_node_BC = ["0 cm"] + [f"{i+0.5} cm" for i in range(0,TotalThick)] + [f"{TotalThick} cm"] 
x_node = [f"{i+0.5} cm" for i in range(0,TotalThick)]
x_interface = [f"{i} cm" for i in range(0,TotalThick+1)]
time_index = [f"{round(ts_h[i],3)} h" for i in range(ts_PST,ts_TST)]


### Calculation condition ---------------------------------------------------------------
T_init = tn.D2K(20) 
Toa = np.array([tn.D2K(20 + 10*(math.sin(2*math.pi*n/tn.d2h))) for n in ts_h])
T0 = Toa.copy()

Tis = np.array([tn.D2K(20) for _ in ts_h])
q_rad = np.array([max(0, 200*(math.sin(2*math.pi*n/tn.d2h))) for n in ts_h])

### Run ---------------------------------------------------------------
# for Widx in tqdm(WallList):

Widx = WallList[0]

## Set wall & matrix, vector
Wall = tn.SetConstruction(*Widx) # *Widx: unpacking list  

# Define system parameters
N = Wall.N
K, K_L, K_R = Wall.K()
R, R_L, R_R = Wall.R()
dx, dx_L, dx_R = Wall.dx()
C = Wall.C()
K_tot = Wall.K_tot()
R_tot = Wall.R_tot()


# Base vector and matrix 
Vector = np.array([None for _ in range(tN)])  # (tN,) vector
Matrix = np.full((tN, N), None).astype(float) # (tN x N) matrix

T_L, T, T_R = Matrix.copy(), Matrix.copy(), Matrix.copy() #Temperature (Fully unsteady state)
q_in, q, q_out, Carnot_eff_L, Carnot_eff_R = Matrix.copy(), Matrix.copy(), Matrix.copy(), Matrix.copy(), Matrix.copy() #Heatflux (Fully unsteady state)

# Initial Temperature, boundary condition
T_init = tn.D2K(20) 
T, T_L, T_R = Matrix.copy(), Matrix.copy(), Matrix.copy()
T[0,:], T_L[0,:], T_R[0,:] = T_init, T_init, T_init
T_R[:,N-1] = Tis[:] 
q_in[0,:], q[0,:], q_out[0,:] = 0, 0, 0
Carnot_eff_L[0,:] = 1 - (T0[0])/(T_L[0,:]) # Carnot efficiency
Carnot_eff_R[0,:] = 1 - (T0[0])/(T_R[0,:]) # Carnot efficiency


# broadcasting
T0 = T0.reshape(len(T0),1)

## TDMA
for n in range(tN-1):

    # Temperature
    T[n+1,:] = tn.TDMA(Wall, T[n,:], T_L[n,:], T_R[n,:], dt) # Calculate next time step temperature T[n+1,:]
    # T_L[n+1,0] = (1 - R_co/(R_co + R_L[0]))*(T0[n+1,0] - T[n+1,0]) + T[n+1,0] # Tos
    T_L[n+1,0] = (T[n+1,0]/R_L[0] + Toa[n+1]/R_co + q_rad[n+1])/(1/R_L[0] + 1/R_co) # Tos
    T_L[n+1,1:] = T[n+1,:-1] + (T[n+1,1:] - T[n+1,:-1]) * (dx[:-1] / (dx[:-1] + dx[1:])) # Linear interpolation
    T_R[n+1,:-1] = T[n+1,:-1] + (T[n+1,1:] - T[n+1,:-1]) * (dx[:-1] / (dx[:-1] + dx[1:])) # Linear interpolation

    # Heatflux
    q_in[n+1,1:] = K_L[1:]*(T[n+1,:-1] - T[n+1,1:])
    q_in[n+1,0] = K_L[0]*(T_L[n+1,0] - T[n+1,0])
    q_out[n+1,:-1] = K_R[:-1]*(T[n+1,:-1] - T[n+1,1:])
    q_out[n+1,-1] = K_R[-1]*(T[n+1,-1] - T_R[n+1,-1])
    q[n+1,:] = (q_in[n+1,:] + q_out[n+1,:])/2

    # Carnot efficiency
    Carnot_eff_L[n+1,:] = 1 - T0[n+1]/T_L[n+1,:] # Carnot efficiency
    Carnot_eff_R[n+1,:] = 1 - T0[n+1]/T_R[n+1,:] # Carnot efficiency

## eliminate pre-simulation data
T0 = tn.eliminate_pre_simulation_data(T0, ts_PST)
T_L = tn.eliminate_pre_simulation_data(T_L, ts_PST)
T = tn.eliminate_pre_simulation_data(T, ts_PST)
T_R = tn.eliminate_pre_simulation_data(T_R, ts_PST)
q_in = tn.eliminate_pre_simulation_data(q_in, ts_PST)
q = tn.eliminate_pre_simulation_data(q, ts_PST)
q_out = tn.eliminate_pre_simulation_data(q_out, ts_PST)
Carnot_eff_L = tn.eliminate_pre_simulation_data(Carnot_eff_L, ts_PST)
Carnot_eff_R = tn.eliminate_pre_simulation_data(Carnot_eff_R, ts_PST)

# half time step (time step - 1)
T0_hf = tn.half_time_matrix(T0, axis=0)
T_L_hf = tn.half_time_matrix(T_L, axis=0)
T_hf = tn.half_time_matrix(T, axis=0)
T_R_hf = tn.half_time_matrix(T_R, axis=0)
q_in_hf = tn.half_time_matrix(q_in, axis=0)
q_hf = tn.half_time_matrix(q, axis=0)
q_out_hf = tn.half_time_matrix(q_out, axis=0)
Carnot_eff_L_hf = tn.half_time_matrix(Carnot_eff_L, axis=0)
Carnot_eff_R_hf = tn.half_time_matrix(Carnot_eff_R, axis=0)

# reshape for broadcasting
Tos = T_L[:,0].reshape(len(T_L),1) 
Tis = T_R[:,-1].reshape(len(T_R),1) 
Tos_hf = T_L_hf[:,0].reshape(len(T_L_hf),1) 
Tis_hf = T_R_hf[:,-1].reshape(len(T_R_hf),1)
qis_hf = q_in_hf[:,0].reshape(len(q_in_hf),1) 
qos_hf = q_out_hf[:,-1].reshape(len(q_out_hf),1)
Carnot_eff_R_hf = Carnot_eff_R_hf[:,-1].reshape(len(Carnot_eff_R_hf),1)


## (Fully unsteady state)------------------------------------------------------------

# data organization
T = np.concatenate((Tos_hf,T_hf[:,:],Tis_hf),axis=1) 
q_node = q_hf
q_intf = np.concatenate((q_in_hf ,qos_hf),axis=1)
CXcR= (1/K)*(T0_hf*(q_hf/T_hf)**2) # Exergy consumption rate [W/m2]
Carnot_eff_hf = np.concatenate((Carnot_eff_L_hf, Carnot_eff_R_hf),axis=1) # Carnot efficiency

# dataframe
T_df = pd.DataFrame(tn.K2D(T), columns=x_node_BC)
q_node_df = pd.DataFrame(q_node, columns=x_node)
q_intf_df = pd.DataFrame(q_intf, columns=x_interface)
CXcR_df = pd.DataFrame(CXcR, columns=x_node)
Carnot_eff_df = pd.DataFrame(Carnot_eff_hf, columns=x_interface)

# indexing
T_df.index = time_index
q_node_df.index = time_index
q_intf_df.index = time_index
CXcR_df.index = time_index

plt.plot(tn.K2D(Tos))
    # # Save
    # filepath = "../output/"
    # filename = f"{WallNameList[WallList.index(Widx)]}.xlsx"
    # with pd.ExcelWriter(filepath + filename) as writer:
    #     T_df.to_excel(writer, sheet_name='T', index=False)
    #     q_node_df.to_excel(writer, sheet_name='q_node', index=False)
    #     q_intf_df.to_excel(writer, sheet_name='q_intf', index=False)
    #     CXcR_df.to_excel(writer, sheet_name='XcR', index=False)
    #     Carnot_eff_df.to_excel(writer, sheet_name='Carnot_eff', index=False)
from fenics import *
import numpy as np
import sympy as sym
from NC5_NonDim import NC5_NonDim
from pdgf import pdgf
from sympy import symbols
from sympy import atan2,Abs
import os

#Upload the mesh
mesh = Mesh('Mesh.xml')  #imports the converted mesh
###############################################################

# Optimization options for the form compiler
parameters['form_compiler']['cpp_optimize'] = True
ffc_options = {'optimize': True}
###############################################################

# define function space
V = VectorFunctionSpace(mesh, 'P', 1)
R = VectorFunctionSpace(mesh, 'R', 0)
SS = FunctionSpace(mesh,'P',1)
TSpace = TensorFunctionSpace(mesh,'P',1)
N = V.dim()
d = mesh.geometry().dim()
###############################################################

# Facet functions
Volume = MeshFunction('size_t' , mesh , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
bnd_mesh = MeshFunction('size_t', mesh , 'Mesh_facet_region.xml')  #saves the boundary info of the mesh
###############################################################


#Changing the directory to the lustre directory
os.chdir("HPC lustre directory")
######################################################################

#Plotting the volume and facet
file=File('eccentric_domain_2D_NC_2_plots_main/boundary.pvd')
file<<bnd_mesh
file=File('eccentric_domain_2D_NC_2_plots_main/volume.pvd')
file<<Volume
###############################################################

# Construct integration measure using these markers
ds = Measure('ds', subdomain_data=bnd_mesh)
dx = Measure('dx', subdomain_data=Volume)
###############################################################


# Define functions
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
u  = Function(V)                 # Displacement from previous iteration
###############################################################


# Elasticity parameter for the Neo-Hookean cylinder (intima)
mu_i, nu_i, beta_i, eta_i, rho_i, phi_i = 27.9, 0.49, 170.88, 263.66, 0.51, pi/3

# Elasticity parameter for the Hopzafel cylinder (Media)
mu_m, nu_m, beta_m, eta_m, rho_m, phi_m = 1.27, 0.49, 8.21, 21.60, 0.25, pi/8

## Elasticity parameter for the Hopzafel cylinder (adventitia)
mu_a, nu_a, beta_a, eta_a, rho_a, phi_a = 7.56, 0.49, 85.03, 38.57, 0.5, pi/2.5
###############################################################




#fibers for intitma
#in case of 3-D you should add x[2] for z
x, y, z= sym.symbols('x[0], x[1], x[2]')
theta = sym.atan2(y,x)
t = theta + pi
r_1 =  (0.714858)+(0.060788)*sym.cos(t)+(0.014752)*sym.cos(2*t)+(-0.006091)*sym.cos(3*t)+(0.004795)*sym.cos(4*t)+(0.000545)*sym.cos(5*t)+(-0.001939)*sym.cos(6*t)+(0.000596)*sym.cos(7*t)+(-0.002890)*sym.cos(8*t)+(0.079346)*sym.sin(t)+(-0.015798)*sym.sin(2*t)+(-0.007070)*sym.sin(3*t)+(-0.013298)*sym.sin(4*t)+(0.001152)*sym.sin(5*t)+(0.003392)*sym.sin(6*t)+(-0.001299)*sym.sin(7*t)+(0.001143)*sym.sin(8*t)
r_2 = (1.034536)+(-0.030658)*sym.cos(t)+(-0.066222)*sym.cos(2*t)+(0.010444)*sym.cos(3*t)+(-0.004373)*sym.cos(4*t)+(0.002543)*sym.cos(5*t)+(-0.000342)*sym.cos(6*t)+(0.001576)*sym.cos(7*t)+(-0.000734)*sym.cos(8*t)+(-0.163745)*sym.sin(t)+(0.042981)*sym.sin(2*t)+(-0.006253)*sym.sin(3*t)+(0.000991)*sym.sin(4*t)+(0.005113)*sym.sin(5*t)+(-0.001116)*sym.sin(6*t)+(-0.004498)*sym.sin(7*t)+(0.002001)*sym.sin(8*t)
f_1_x = (r_1)*sym.cos(t)
f_1_y = (r_1)*sym.sin(t)
f_2_x = (r_2)*sym.cos(t)
f_2_y = (r_2)*sym.sin(t)
w_A_x = sym.diff(f_1_x,theta)
w_A_y = sym.diff(f_1_y,theta)
w_C_x = sym.diff(f_2_x,theta)
w_C_y = sym.diff(f_2_y,theta)
AE = ((x-f_1_x)**2+(y-f_1_y)**2)**0.5
EC = ((x-f_2_x)**2+(y-f_2_y)**2)**0.5
u_s = (AE*w_A_x+EC*w_C_x)/(AE+EC)
v_s = (AE*w_A_y+EC*w_C_y)/(AE+EC)
w_s = (sym.tan(phi_i)*(u_s**2+v_s**2)**0.5)
M_0 = 1/((u_s**2+v_s**2)*(1+sym.tan(phi_i)**2))**0.5
U_s = M_0*u_s
V_s = M_0*v_s
W_s = M_0*w_s
THETA_s = sym.atan2(V_s,U_s)
SinT_0 = sym.sin(THETA_s)
nSinT_0= -sym.sin(THETA_s)
CosT_0 = sym.cos(THETA_s)
u_s_code= sym.printing.ccode(U_s)
v_s_code= sym.printing.ccode(V_s)
w_s_code= sym.printing.ccode(W_s)
SinT_0_code=sym.printing.ccode(SinT_0)
nSinT_0_code=sym.printing.ccode(nSinT_0)
CosT_0_code = sym.printing.ccode(CosT_0)
n_i=Expression((u_s_code,v_s_code), degree=0)
T_i = Expression(((CosT_0_code,SinT_0_code),(nSinT_0_code,CosT_0_code)),degree=0)
b3_i = Expression(w_s_code,degree=0)
###############################################################

#fibers for hopzafel (media)
r_3 = (1.234536)+(-0.030658)*sym.cos(t)+(-0.066222)*sym.cos(2*t)+(0.010444)*sym.cos(3*t)+(-0.004373)*sym.cos(4*t)+(0.002543)*sym.cos(5*t)+(-0.000342)*sym.cos(6*t)+(0.001576)*sym.cos(7*t)+(-0.000734)*sym.cos(8*t)+(-0.163745)*sym.sin(t)+(0.042981)*sym.sin(2*t)+(-0.006253)*sym.sin(3*t)+(0.000991)*sym.sin(4*t)+(0.005113)*sym.sin(5*t)+(-0.001116)*sym.sin(6*t)+(-0.004498)*sym.sin(7*t)+(0.002001)*sym.sin(8*t)
f_2_x = (r_2)*sym.cos(t)
f_2_y = (r_2)*sym.sin(t)
f_3_x = (r_3)*sym.cos(t)
f_3_y = (r_3)*sym.sin(t)
v_A_x = sym.diff(f_2_x,theta)
v_A_y = sym.diff(f_2_y,theta)
v_C_x = sym.diff(f_3_x,theta)
v_C_y = sym.diff(f_3_y,theta)
AB = ((x-f_2_x)**2+(y-f_2_y)**2)**0.5
BC = ((x-f_3_x)**2+(y-f_3_y)**2)**0.5
u_q = (AB*v_A_x+BC*v_C_x)/(AB+BC)
v_q = (AB*v_A_y+BC*v_C_y)/(AB+BC)
w_q = (sym.tan(phi_m)*(u_q**2+v_q**2)**0.5)
M_1 = 1/((u_q**2+v_q**2)*(1+sym.tan(phi_m)**2))**0.5
U_q = M_1*u_q
V_q = M_1*v_q
W_q = M_1*w_q
THETA_q = sym.atan2(V_q,U_q)
SinT_1 = sym.sin(THETA_q)
nSinT_1= -sym.sin(THETA_q)
CosT_1 = sym.cos(THETA_q)
u_q_code= sym.printing.ccode(U_q)
v_q_code= sym.printing.ccode(V_q)
w_q_code= sym.printing.ccode(W_q)
SinT_1_code=sym.printing.ccode(SinT_1)
nSinT_1_code=sym.printing.ccode(nSinT_1)
CosT_1_code = sym.printing.ccode(CosT_1)
n_m=Expression((u_q_code,v_q_code), degree=0)
T_m = Expression(((CosT_1_code,SinT_1_code),(nSinT_1_code,CosT_1_code)),degree=0)
b3_m = Expression(w_q_code,degree=0)
###############################################################

#fibers for hopzafel (adventitia)
r_4 = (1.334536)+(-0.030658)*sym.cos(t)+(-0.066222)*sym.cos(2*t)+(0.010444)*sym.cos(3*t)+(-0.004373)*sym.cos(4*t)+(0.002543)*sym.cos(5*t)+(-0.000342)*sym.cos(6*t)+(0.001576)*sym.cos(7*t)+(-0.000734)*sym.cos(8*t)+(-0.163745)*sym.sin(t)+(0.042981)*sym.sin(2*t)+(-0.006253)*sym.sin(3*t)+(0.000991)*sym.sin(4*t)+(0.005113)*sym.sin(5*t)+(-0.001116)*sym.sin(6*t)+(-0.004498)*sym.sin(7*t)+(0.002001)*sym.sin(8*t)
f_4_x = (r_4)*sym.cos(t)
f_4_y = (r_4)*sym.sin(t)
u_A_x = sym.diff(f_3_x,theta)
u_A_y = sym.diff(f_3_y,theta)
u_C_x = sym.diff(f_4_x,theta)
u_C_y = sym.diff(f_4_y,theta)
AD = ((x-f_3_x)**2+(y-f_3_y)**2)**0.5
DC = ((x-f_4_x)**2+(y-f_4_y)**2)**0.5
u_p = (AD*u_A_x+DC*u_C_x)/(AD+DC)
v_p = (AD*u_A_y+DC*u_C_y)/(AD+DC)
w_p = (sym.tan(phi_a)*(u_p**2+v_p**2)**0.5)
M_2 = 1/((u_p**2+v_p**2)*(1+sym.tan(phi_a)**2))**0.5
U_p = M_2*u_p
V_p = M_2*v_p
W_p = M_2*w_p
THETA_p = sym.atan2(V_p,U_p)
SinT_2 = sym.sin(THETA_p)
nSinT_2= -sym.sin(THETA_p)
CosT_2 = sym.cos(THETA_p)
u_p_code= sym.printing.ccode(U_p)
v_p_code= sym.printing.ccode(V_p)
w_p_code= sym.printing.ccode(W_p)
SinT_2_code=sym.printing.ccode(SinT_2)
nSinT_2_code=sym.printing.ccode(nSinT_2)
CosT_2_code = sym.printing.ccode(CosT_2)
n_a = Expression((u_p_code,v_p_code), degree=0)
T_a = Expression(((CosT_2_code ,SinT_2_code),(nSinT_2_code,CosT_2_code )),degree=0)
b3_a = Expression(w_p_code,degree=0)
###############################################################

# blood pressure direction
n = FacetNormal(mesh)
###############################################################

#growth for the media
g_inv_m = Expression('1/g_m', degree=0, g_m=1)
ginv_m = (T_m.T)*(T_m)
G_m = det(ginv_m)
###############################################################

#growth for the adventitia
ginv_a = Expression('1/g_a', degree=0, g_a=1)
ginv_a = (T_a.T)*(T_a)
G_a = det(ginv_a)
###############################################################

#animation file for NC4
vtkfile2 = File('D_2D_1/LDL.pvd')
vtkfile3 = File('D_2D_1/Macrophage.pvd')
vtkfile4 = File('D_2D_1/Deadcells.pvd')
vtkfile5 = File('D_2D_1/oxygen.pvd')
vtkfile6 = File('D_2D_1/MCP.pvd')
vtkfile7 = File('D_2D_1/PDGF1.pvd')
vtkfile8 = File('D_2D_1/PDGF2.pvd')
volumefile1 = File('D_2D_1/Volume_intima.pvd')
volumefile2 = File('D_2D_1/Volume_artery.pvd')
test1 = File('D_2D_1/g_i_test.pvd')
###############################################################



#Parameters for the PDGF induced Growth
g_alpha = Constant(1)
g_beta = Constant(1)
PDGF_0 = 0.2
m = 4
GAMMA_Max = 6.5
kappa = 5
alpha=1
beta = 0.25
###############################################################

#loop parameters
sigma, t = 0.0,0.0
dsig, dt = 2.0, 0.001
sigmax, t_max= 16.0, 3.0
k = int(0)
Energy = []
###############################################################

#Loop
while t<=t_max:

     #creating the growth tensor as a function of PDGF
     if k==0:
         G=Expression((('1','0'),('0','1')),degree=1)
     else:
         G = Expression((('g_t','0'),('0','g_r')),degree=1,g_t = g_beta, g_r = g_alpha)
     GG = (T_i.T)*(G)*(T_i)
     ginv_i = inv(GG)
     G_i = det(GG)

     if k>0:
         test_g = project(Expression('g_r', degree=1, g_r=g_alpha),SS)
         test_g.rename('test_g','test_g')
         test1<<(test_g,t) #plot in case you want to see the Growth
     ###############################################################


     # Kinematics for initma
     II = Identity(d)            # Identity tensor
     F =II + grad(u)            # Deformation gradient
     Fe_i= variable(F*ginv_i)
     C_i = F.T*F                   # Right Cauchy-Green tensor
     Ce_i= Fe_i.T*Fe_i
     I4_i = dot(n_i, Ce_i*n_i)+b3_i**2
     J = det(F)
     ###############################################################

     # Kinematics for media
     Fe_m= variable(F*ginv_m)                #Technically does not do anything because ginv_m is 1
     C_m = F.T*F                   # Right Cauchy-Green tensor
     Ce_m= Fe_m.T*Fe_m
     I4_m = dot(n_m, Ce_m*n_m)+b3_m**2
     ###############################################################

     # Kinematics for adventitia
     Fe_a= variable(F*ginv_a)                 #Technically does not do anything because ginv_a is 1
     C_a = F.T*F                   # Right Cauchy-Green tensor
     Ce_a= Fe_a.T*Fe_a
     I4_a = dot(n_a, Ce_a*n_a)+b3_a**2
     ###############################################################

     # Invariants of deformation tensors (intima)
     I_i = tr(C_i)
     Ie_i=tr(Ce_i)+ Constant(1)
     J_i = det(Fe_i)   #Jacobian of F or Fe?
     ###############################################################

     # Invariants of deformation tensors (media)
     I_m = tr(C_m)
     Ie_m=tr(Ce_m)+ Constant(1)
     J_m = det(Fe_m)   #Jacobian of F or Fe?
     ###############################################################

      # Invariants of deformation tensors (Adventitia)
     I_a = tr(C_a)
     Ie_a=tr(Ce_a)+ Constant(1)
     J_a = det(Fe_a)
     ###############################################################

     #deformed normal vector
     NansonOp = J*inv(F).T
     deformed_NN = dot(NansonOp,n)
     NormN = sqrt(dot(deformed_NN,deformed_NN))
     deformed_N = as_vector([deformed_NN[0],deformed_NN[1]])
     ######################################################################

     # Stored strain energy density (neo-Hookean Hopzafel model)
     psi_i = (mu_i/2)*(Ie_i - 3)+(eta_i/beta_i)*(exp(beta_i*(rho_i*(I4_i-1)**2+(1-rho_i)*(Ie_i-3)**2))-1)+((mu_i*nu_i)/(1-2*nu_i))*(J_i-1)**2- mu_i*ln(J_i)
     psi_m = (mu_m/2)*(Ie_m - 3)+(eta_m/beta_m)*(exp(beta_m*(rho_m*(I4_m-1)**2+(1-rho_m)*(Ie_m-3)**2))-1)+((mu_m*nu_m)/(1-2*nu_m))*(J_m-1)**2- mu_m*ln(J_m)
     psi_a = (mu_a/2)*(Ie_a - 3)+(eta_a/beta_a)*(exp(beta_a*(rho_a*(I4_a-1)**2+(1-rho_a)*(Ie_a-3)**2))-1)+((mu_a*nu_a)/(1-2*nu_a))*(J_a-1)**2- mu_a*ln(J_a)
     ###############################################################

     #Stress variables for the weak form
     TT_i = G_i*diff(psi_i,Fe_i)*(ginv_i.T)
     TT_m = diff(psi_m,Fe_m)*(ginv_m.T)
     TT_a = diff(psi_a,Fe_a)*(ginv_a.T)
     ######################################################################

     # Total strain energy
     Pi_E = psi_i*dx(5)+psi_m*dx(6)+psi_a*dx(7)
     ######################################################################

     #The weak form
     try:
         Pi = inner(TT_i,grad(v))*dx(5)+inner(TT_m,grad(v))*dx(6)+inner(TT_a,grad(v))*dx(7)+(sigma)*dot(deformed_N,v)*ds(1)
     except:
         Pi = inner(TT_i,grad(v))*dx(5)+inner(TT_m,grad(v))*dx(6)+inner(TT_a,grad(v))*dx(7)
     ######################################################################

     Energy.append(assemble(Pi_E))

     # Compute first variation of Pi (directional derivative about u in the direction of v)
     gradPi = derivative(Pi, u, v)
     ###############################################################

     # Compute Jacobian of gradPi
     hPi = derivative(gradPi, u, du)
     ###############################################################

     #Solve
     solve(Pi == 0, u ,
      solver_parameters={"newton_solver": {"relative_tolerance": 5e-9,
      "absolute_tolerance": 5e-9,"maximum_iterations": 10}})
     ###############################################################


     #If statements for pressure plotting exclusion and saving memory by saving every 10 steps
     if sigma<=sigmax:
          sigma+=dsig
     if sigma>sigmax:
          PDGF, kappa =pdgf(k,t,u,vtkfile7,vtkfile8,kappa)
          temp = project(PDGF,SS)
          MAX = temp.vector().max()
          print('Max PDGF is', MAX)
          if k%10==0:
              NC5_NonDim(t,u,vtkfile2,vtkfile3,vtkfile4,vtkfile5,vtkfile6,volumefile1,volumefile2)
          t+=dt
          k+=1
     if sigma==sigmax+dsig:
          sigma-=dsig

     ###############################################################

     #Printing the loop parmeters and updating them
     print("sigma is", sigma,flush=True)
     print("t is", t,flush=True)
     print("k is", k,flush=True)
     ###############################################################

     #Updating the growth using the PDGF answer
     if k!=0:
         g_alpha =project(g_alpha*exp(alpha*Constant(0.3)*PDGF**m/(PDGF_0**m+PDGF**m)*GAMMA_Max*Constant(dt)),SS)
         g_beta = project(g_beta*exp(beta*Constant(0.3)*PDGF**m/(PDGF_0**m+PDGF**m)*GAMMA_Max*Constant(dt)),SS)
         print('g max is', g_alpha.vector().max())
     ###############################################################


#saving the energy in a csv file
import csv
with open('Energy.csv', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    wr.writerow(Energy)
######################################################################

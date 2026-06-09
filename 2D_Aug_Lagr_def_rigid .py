# -*- coding: utf-8 -*-

# ------------- A CHANGER ------------- #
# Problem considered: 2D
# Version considering Augmented Lagrangian method developed
#         in Getfem (displacement)
# ------------- A CHANGER ------------- #

import numpy as np
import getfem as gf
import sys,os

def linsolve(M, B): # Call Superlu to solve a sparse linear system
    return (((gf.linsolve_superlu(M, B))[0]).T)[0]
    #return (((gf.linsolve_mumps(M, B))[:]).T)[0]
#-----
exemple = 2
#-----
os.system('rm ./results/*.vtu')
np.set_printoptions(threshold=sys.maxsize)
fic=open('energy.txt','w')

gf.util('trace level', 0)

md = gf.Model("real")
#-----
# Données physique du problème
#-----
#E =   1.e2                    # Young's modulus [Pa]
#nu =  0.35                     # Poisson's ratio [-]
#rho = 500.0                  # density
rho = 1.e-2
#-----
# Pour le frottement
mu = 0.0
#--
#C_nu  = 1.e+04
#C_tau = 1.e+03
#eps = 1.e-10 # Pour éviter de diviser par zero... on ajoute epsilon au denominateur s'il est petit
pen_factor = 1.e+12 # Pour imposer les CL de type Dirichlet
#-----
a_plane = -0.0 # pente du socle rigide
y_plane = -0.0 # Hauteur du socle rigide 
#-----------------------------------------------------------------
#Time related parameters
T_max = 1.0 # 10.e-3 # 400.e-4 #0.25 # 1.0 #  0.4 # 
dt =1.e-4
#steps = 10001 #25 #  1001 #
#dt = T_max/steps
steps = T_max/dt
#-----------------------------------------------------------------
# Visu parameters
visu = True
delta_visu = 200 # Nombre d'itérations entre 2 images
#-----------------------------------------------------------------
# Test d'arret
Newton_ite_max = 20
eps_Newton = 1.e-8
#-----------------------------------------------------------------

#-----
#Elements finis
#-----
fem_order = 1
Dim = 2                            # Dimension du problème (toujours 2 ici)
GT    = 'GT_PK'                    # Type d'Eléments Finis : triangle
FEM   = 'FEM_PK'
Integ = 'IM_TRIANGLE(7)'           # Méthode d'intégration sur le triangle
### If we use quadrangle
#GT    = 'GT_QK'                    # Type d'Eléments Finis : : quadrangle
#FEM   = 'FEM_QK'
#Integ = 'IM_QUAD(9)'               # Méthode d'intégration sur le quadrangle

# traction_force = 10 #Surfacic load on the Neumann boundary of the body in the vertical direction
# surfaceLoad=[0, 0, -traction_force]
# volumeLoad=[0, 0, 0]
#-----
#m=gf.Mesh('regular simplices', np.arange(0,1+meshSize,meshSize), 
#        np.arange(0,4+meshSize,meshSize),np.arange(0,1+meshSize,meshSize))
# Integration method
#mim=gf.MeshIm(m, gf.Integ('IM_TETRAHEDRON(8)'))
#-----
if (exemple == 1):
   nbElems = 5.
   meshSize = 1./nbElems
   m = gf.Mesh('import','structured','GT="'+GT+'('+str(Dim)+','+str(fem_order)+')"; \
               SIZES=[%f,%f];NOISED=0;NSUBDIV=[%d,%d];' % (4.,1.,2*nbElems, nbElems));
#-----
if (exemple == 2): # Ring
   center = np.array([-0.0 , +1.001]);
   rad_int = 0.7
   rad_ext = 1.0
#   m = gf.Mesh('import','gmsh','./anneau/anneau2D_1.msh')
#   m.translate(center)
   mo1 = gf.MesherObject('ball', center, rad_ext)
   mo2 = gf.MesherObject('ball', center, rad_int)
   mo3 = gf.MesherObject('set minus', mo1, mo2)
   print('Meshes generation')
   m = gf.Mesh('generate', mo3, (rad_ext-rad_int)/10., 2) #, [4,20])
#-----
if (exemple == 3): # Ball
#   center = np.array([0.0 , +0.5]);
   center = np.array([0.2 , +0.46]);
   rad_ext = 0.50
   mo1 = gf.MesherObject('ball', center, rad_ext)
   print('Meshes generation')
   m = gf.Mesh('generate', mo1, rad_ext/12., 2) #, [4,20])
#-----
#m.export_to_pos('Gonzalez.pos')
#
mf   = gf.MeshFem(m, Dim)
mf.set_fem(gf.Fem(FEM+'('+str(Dim)+','+str(fem_order)+')'))
mf1D = gf.MeshFem(m, 1) # MeshFem 1D for normal and tangential components of the stress  tensor
mf1D.set_fem(gf.Fem(FEM+'('+str(Dim)+','+str(fem_order)+')'))
## For multipliers
mf1 = gf.MeshFem(m, 1)
mf1.set_fem(gf.Fem(FEM+'('+str(Dim)+','+str(fem_order)+')'))
#mf1.set_classical_discontinuous_fem(0) 
#mf1Dnc = gf.MeshFem(m, 1)
#mf1Dnc.set_classical_discontinuous_fem(0) 

#  Integration method on each element
mim = gf.MeshIm(m, gf.Integ(Integ)) 
#-----------------------------------------------------------------
#Dirichlet and Neumann Regions
#Boundaries creation
#-----------------------------------------------------------------
DIRICHLET_BOUNDARY = 1
NEUMANN_BOUNDARY = 2
CONTACT_BOUNDARY = 3

if (exemple == 1):
#   m.set_region(DIRICHLET_BOUNDARY,
#                m.outer_faces_with_direction([-1.,  0.], 0.001)) # Left side
   m.set_region(NEUMANN_BOUNDARY,
                m.outer_faces_with_direction([ 0.,  1.], 0.001)) # Top side
   m.set_region(CONTACT_BOUNDARY,
                m.outer_faces_with_direction([ 0., -1.], 0.001)) # Bottom side
#---
if (exemple == 2):
    fb1 = m.outer_faces_in_box([center[0]-rad_int-0.01, center[1]-rad_int-0.01], [center[0]+rad_int+0.01, center[1]+rad_int+0.01])  # Boundary of the hole
    fb2 = m.outer_faces_with_direction([0., -1.], np.pi/2.0) # Contact boundary of the wheel
    HOLE_BOUND=4;
    m.set_region(HOLE_BOUND, fb1)
    m.set_region(CONTACT_BOUNDARY, fb2)
    m.region_subtract(CONTACT_BOUNDARY, HOLE_BOUND)
#---
if (exemple == 3):
    fb = m.outer_faces_in_box([center[0]-rad_ext-0.01,center[1]-rad_ext-0.01],[center[0]+rad_ext+0.01,center[1]+rad_ext+0.01])
    m.set_region(CONTACT_BOUNDARY, fb)
#---
#    m.export_to_pos('Gonzalez.pos')
#-----------------------------------------------------------------
md.add_fem_variable("u", mf)    # displacements
# Initialisation
u0 = md.interpolation("[0.,0.]", mf)
md.add_initialized_fem_data("u_old", mf, u0)
md.set_variable('u', u0)
#
U = md.variable('u')
U_old = md.variable('u_old')  #np.copy(U)
NTDL = np.size(U)
#
#Mass_matrix = -(2.*rho/(dt*dt))*gf.asm_mass_matrix(mim,mf)
#
#-----------------------------------------------------------------
# boundary conditions
#-----------------------------------------------------------------
if (exemple==1): # Ici, on impose une force constante au cours du temps
#    #md.add_Dirichlet_condition_with_multipliers(mim, 'u', 0, DIRICHLET_BOUNDARY)
#    u_d = md.interpolation("[0., 0.]",mf)
#    md.add_initialized_fem_data("DirichletData", mf, u_d)
#    #md.add_Dirichlet_condition_with_multipliers(mim, 'u', mf,  DIRICHLET_BOUNDARY, "DirichletData")
#    md.add_Dirichlet_condition_with_simplification("u", DIRICHLET_BOUNDARY, "DirichletData")
#--
    ff=mf.eval('[0.,-100.*'+str(dt)+']')
#-------
if (exemple==2): # Ici, la force imposee est nulle
    ff=mf.eval('[0.,0.*'+str(dt)+']')
#-------
if (exemple==3): # Ici, la force imposee (pesanteur)
    ff=mf.eval('[-1000.,-2000.00]')
#-----------------------------------------------------------------
#    
md.add_initialized_fem_data('f', mf, ff)
#    md.add_linear_term(mim, "-f.Test_u", NEUMANN_BOUNDARY)
md.add_source_term_brick(mim, "u", "f", NEUMANN_BOUNDARY)
if (exemple==3):
    md.add_source_term_brick(mim, "u", "f") # pesanteur (force interne)
#md.assembly()
#sec_mem = md.rhs()
##print('Second membre : '+str(sec_mem))
#-------
if (exemple==1): # Dans cet exemple, la vitesse initiale est nulle
#   v0 = mf.eval("[0.,-x*x/16.]")
   v0 = mf.eval("[0.,0]")
if (exemple==2): # Dans cet exemple, la vitesse initiale est imposée
   v0 = mf.eval("[0.5,-1]")
#   v0 = mf.eval("[-x-0.2,0]") # Exemple pour faire de la deformation independemment du contact
#v0 = md.interpolation("[0.,0.,0.]", mf)
if (exemple==3): # Dans cet exemple, la vitesse initiale est parallele au support
#   v0 = mf.eval("[0.1,-0.050]")
   v0 = mf.eval("[0.3,-0.150]")
md.add_initialized_fem_data("v_old", mf, v0)
V = md.variable('v_old') # Pas terrible si vitesse initiale non nulle 
V_old = V
md.add_initialized_fem_data("v", mf, v0)    # velocities
md.set_time_step(dt)
#md.add_initialized_data("K", E/(3.*(1.-nu)))
#md.add_initialized_data("G", E/(2.*(1.+nu)))
#
md.add_macro('F',  '(Id(meshdim)+Grad_u)')
md.add_macro('C','Right_Cauchy_Green(F)')
#md.add_macro('C',"(Id(meshdim)+Grad_u+Grad_u')") # Version linéarisée
#
md.add_macro('F_old', '(Id(meshdim)+Grad_u_old)')
md.add_macro('C_old','Right_Cauchy_Green(F_old)')
#md.add_macro('C_old',"(Id(meshdim)+Grad_u_old+Grad_u_old')") # Version linéarisée
#
md.add_macro('F_avg', '0.5*(F+F_old)')
md.add_macro('C_avg', '0.5*(C+C_old)')
#
# Loi de Ogden/Ciarlet-Geymonat
#clambda = E*nu/((1.+nu)*(1.-2.*nu)) # First Lame coefficient (N/cm^2)
#cmu = E/(2.*(1.+nu))               # Second Lame coefficient (N/cm^2)
#clambda = E*nu/((1.+nu)*(1.-nu))
#a = 0.5*(max(0.,0.5*cmu-0.25*clambda)+0.5*cmu)
#c1 = 0.5*cmu-a
#c2 = 0.25*clambda-0.5*cmu+a
#c3 = 0.5*cmu+0.25*clambda
a = 0.35   #0.5*(max(0.,0.5*cmu-0.25*clambda)+0.5*cmu)
c1 = 0.5   #0.5*cmu-a
c2 = 0.5e-2   #0.25*clambda-0.5*cmu+a
c3 = 0.35     #0.5*cmu+0.25*clambda
para = [a,c1,c2,c3]
print(para)
md.add_initialized_data("c_", para)
#
# Energie elastique et contrainte Ciarlet-Geymonat en 2D
md.add_macro("W(CC)", "c_(1)*(Trace(CC)-2.)+c_(2)*(Trace(CC)+Det(CC)-3.)+c_(3)*(Det(CC)-1)-(c_(1)+2.*c_(2)+c_(3))*log(Det(CC))")
md.add_macro("S(CC)", "2.*( (c_(1)+c_(2))*Id(meshdim) + (c_(2)+c_(3))*Det(CC)*Inv(CC) - (c_(1)+2.*c_(2)+c_(3))*Inv(CC) )")
#
#
# Gonzalez:
md.add_initialized_data("eps", 1.e-8)
md.add_macro("SS", "S(C_avg) + (2.*(W(C)-W(C_old))-S(C_avg):(C-C_old))/(Norm_sqr(C-C_old)+eps)*(C-C_old)")
#md.add_macro("SS", "S(C_avg) + (2.*W(C)-2.*W(C_old)-S(C_avg):(C-C_old))/(Norm_sqr(C-C_old)+eps)*(C-C_old)")
# Implicit Euler:
#md.add_macro("SS", "S(C)")
# Explicit Euler
#md.add_macro("SS", "S(C_old)")
 
# elasticity + inertia
md.add_initialized_data("rho", rho)
md.add_initialized_data("dt", dt)
# Gonzalez:
#Mass = md.add_linear_term(mim,"(2.*rho/(dt*dt))*u.Test_u")
md.add_linear_term(mim, "2.*rho/dt*((u-u_old)/dt-v_old).Test_u")
# Implicit Euler: -
#md.add_linear_term(mim, "rho/dt*((u-u_old)/dt-v_old).Test_u")
#
#md.add_nonlinear_term(mim, "(SS):Grad_Test_u")
md.add_nonlinear_term(mim, "(F_avg.SS):Grad_Test_u")
#md.add_linear_term(mim, "(F_old.S(C_old)):Grad_Test_u") # To try explicit
# Linear
#md.add_nonlinear_term(mim, "S(C):Grad_Test_u")
#
#Mass_matrix = md.matrix_term(Mass, 0)
# 
mf.export_to_vtu('./results/Displacement_%05i.vtu' %0, mf, U, 'Displacement',V,'Velocity')
#--------
# Initial energy
energy = gf.asm_generic(mim, 0, "0.5*rho*Norm_sqr(v_old) + W(C)", -1, md)
print("==> mech energy = ", energy)
fic.write(f"Time step: {0.0: 9.5f}, Energy: { energy:.17f}\n")
#
#--------
#  Frictional contact
#--------
md.add_initialized_data('mu',[mu])
#---
# Nitsche method:
#---
#Points = mf1D.basic_dof_nodes() #: Points of the mesh
##Points = Points[:,Contact_S_dof]
#Dist = Points[1,:]  # - y_plane :Distance signée : ici, le socle est à y=0, donc distance = y_point
#print('===> Dist'+str(Dist))
#md.add_initialized_fem_data('Dist', mf1D, Dist)
#theta_contact = 0  # Also possible : +1 or -1
#md.add_initialized_data('theta_contact',[theta_contact])
#gamma0 = 1.e+4
#md.add_initialized_data('gamma0',[gamma0])
#Neumannterm = md.Neumann_term('u', CONTACT_BOUNDARY)
#md.add_Nitsche_contact_with_rigid_obstacle_brick(mim, 'u', Neumannterm, "Dist", 
#                               "gamma0", CONTACT_BOUNDARY, theta_contact, 'mu')
#---
# Augmented Lagrangian (obstacle is 'y<0')
r = 1.e+6
md.add_initialized_data('r',[r])
md.add_filtered_fem_variable('lambda_n', mf1, CONTACT_BOUNDARY)
md.add_filtered_fem_variable('lambda_t', mf1, CONTACT_BOUNDARY)
md.add_nodal_contact_with_rigid_obstacle_brick(mim, 'u', 'lambda_n', 'lambda_t', 
                                               'r','mu',CONTACT_BOUNDARY,'y',2)
# Time loop
for timeStepIndex,timeStep in enumerate(np.arange(0.,T_max+dt,dt)):
    print('Time step: %9.5f' % timeStep)
    dU = dt*V #0.*U #   Une initialisation pour pouvoir calculer la contrainte à la première itération : correct ???
    U = U + dU
    md.set_variable("u",U)
    dU_old = 0.*U # Variable pour faire un test qui évite les cycles sur 2 itérations
#=====
# Contact (Nitche method)
#=====
#    Dist = Points[1,:]+U[1::2] #- y_plane # Mise a jour de la distance
#    print('===> Min Dist : '+str(min(Dist)))
#    md.set_variable('Dist',Dist)
#    # debug
#    if timeStepIndex % 100 == 0:  # in mỗi 100 step
#        print(f"t={timeStep:.5f}, Dist_min={Dist.min():.6f}")
#=====
# Resolution (contact frottant et hyperelasticité)
# ===== Solve nonlinear system =====
#=====
    md.solve('lsolver','mumps','max_res',1E-06,'max_iter',20,'very_noisy')
## ---
    U = md.variable("u")
    ###print(U)
#=====
# Updated of the velocity
#   Gonzalez:
    V = 2.*(U-U_old)/dt-V_old
#Post processing of the solution
    if ((visu==True)&(timeStepIndex%delta_visu==0)):
        mf.export_to_vtu('./results/Displacement_%05i.vtu' %(timeStepIndex+1), mf, U, 'Displacement',V,'Velocity')
#===== 
    # Update u_old and v_old
    U_old = U
    md.set_variable('u_old', U_old)
    V_old = V
    md.set_variable('v_old', V_old)
#===== 
    # Energy
    energy = gf.asm_generic(mim, 0, "0.5*rho*Norm_sqr(v_old) + W(C)", -1, md)
    print("==> mech energy = ", energy)
    #fic.write(f"Time step: {timeStep: 9.5f}, Energy: { energy:.17f}\n")
    fic.write(f"{timeStep: 9.5f} { energy:.17f}\n")  #easier to plot
#    print(f"Nb active contact points : {nb_active: 5d}")
#    print(f"Nb active slip points    : {nb_slip_active: 5d}")
#    print(f"Nb active stick points   : {nb_stick_active: 5d}")
    
#
#    md.shift_variables_for_time_integration()
    md.next_iter()
#===== 
# Fin de la boucle de pas de temps
#=====   

fic.close()

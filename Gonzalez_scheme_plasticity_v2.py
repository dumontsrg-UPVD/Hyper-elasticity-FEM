# -*- coding: utf-8 -*-

# ------------- A CHANGER ------------- #
# Problem considered: 2D !!!
# to treat problems with frictional contact and plasticity
# Not very good for the moment !!
# ------------- A CHANGER ------------- #
import sys,os
import numpy as np
from scipy.linalg import expm
import getfem as gf

def linsolve(matr, b): # Call Superlu to solve a sparse linear system
#    return (((gf.linsolve_superlu(matr, b))[0]).T)[0]
    return (((gf.linsolve_mumps(matr, b))[:]).T)[0]
#
def compute_dev(matr):
    m_dev = matr - 0.5*np.trace(matr)*Identity
    # It is 1./3. in 3D !!!
    # Rmq : tensordot fait la somme des Aij*Bij
    norm_h_m_dev = np.sqrt(0.5)*np.sqrt(np.tensordot(m_dev,m_dev.T,axes=2))
    # It is 1.5 in 3D !!!
    return (m_dev,norm_h_m_dev)
#-----
exemple = 4
#-----
os.system('rm ./results/*.vtu')
np.set_printoptions(threshold=sys.maxsize)
fic = open('energy.txt','w',encoding="locale")

gf.util('trace level', 0)

md = gf.Model("real")
#-----
# Données physique du problème
#-----
E =   1.e6                    # Young's modulus [Pa]
nu =  0.3                     # Poisson's ratio [-]
#rho = 500.0                  # density
rho = 1000.0
#-----
# Pour le frottement
mu = 0.5
#--
C_nu  = 1.e+06
C_tau = 1.e+03
eps = 1.e-10 # Pour éviter de diviser par zero... on ajoute epsilon au denominateur s'il est petit
pen_factor = 1.e+12 # Pour imposer les CL de type Dirichlet
#-----
a_plane = -0.0 # pente du socle rigide
y_plane = -0.0 # Hauteur du socle rigide
#-----------------------------------------------------------------
#Time related parameters
T_max = 1.e+0 #2.0 # 10.e-3 # 400.e-4 #0.25 # 1.0 #  0.4 #
dt =1.e-2
#steps = 10001 #25 #  1001 #
#dt = T_max/steps
steps = T_max/dt
#-----------------------------------------------------------------
# Visu parameters
visu = True
delta_visu = 100 # Nombre d'itérations entre 2 images
#-----------------------------------------------------------------
# Test d'arret
Newton_ite_max = 5
eps_Newton = 1.e-8
#-----------------------------------------------------------------
#-----
#Elements finis
#-----
fem_order = 1
Dim = 2                            # Dimension du problème (toujours 2 ici)
GT    = 'GT_PK'                    # Type d'Eléments Finis : triangle
FEM   = 'FEM_PK'
Integ = 'IM_TRIANGLE(1)'           # Méthode d'intégration sur le triangle
### If we use quadrangle
#GT    = 'GT_QK'                    # Type d'Eléments Finis : : quadrangle
#FEM   = 'FEM_QK'
#Integ = 'IM_QUAD(9)'               # Méthode d'intégration sur le quadrangle
#
# traction_force = 10 #Surfacic load on the Neumann boundary of the body in the vertical direction
# surfaceLoad=[0, 0, -traction_force]
# volumeLoad=[0, 0, 0]
#-----
#m=gf.Mesh('regular simplices', np.arange(0,1+meshSize,meshSize),
#        np.arange(0,4+meshSize,meshSize),np.arange(0,1+meshSize,meshSize))
# Integration method
#mim=gf.MeshIm(m, gf.Integ('IM_TETRAHEDRON(8)'))
#-----
if exemple == 1 :
    nbElems = 1.
    meshSize = 1./nbElems
    m = gf.Mesh('import','structured','GT="'+GT+'('+str(Dim)+','+str(fem_order)+')"; \
               SIZES=[%f,%f];NOISED=0;NSUBDIV=[%d,%d];' % (4.,1.,2*nbElems, nbElems))
#-----
if exemple == 2 : # Ring
    center = np.array([-0.2 , +0.5])
    rad_int = 0.35
    rad_ext = 0.50
#   m = gf.Mesh('import','gmsh','./anneau/anneau2D_1.msh')
#   m.translate(center)
    mo1 = gf.MesherObject('ball', center, rad_ext)
    mo2 = gf.MesherObject('ball', center, rad_int)
    mo3 = gf.MesherObject('set minus', mo1, mo2)
    print('Meshes generation')
    m = gf.Mesh('generate', mo3, (rad_ext-rad_int)/10., 2) #, [4,20])
#-----
if exemple == 3 : # Ball
#   center = np.array([0.0 , +0.5]);
    center = np.array([0.2 , +0.46])
    rad_ext = 0.50
    mo1 = gf.MesherObject('ball', center, rad_ext)
    print('Meshes generation')
    m = gf.Mesh('generate', mo1, rad_ext/12., 2) #, [4,20])
#-----
if exemple == 4 :
    nbElems = 16.
    meshSize = 1./nbElems
    m = gf.Mesh('import','structured','GT="'+GT+'('+str(Dim)+','+str(fem_order)+')"; \
                SIZES=[%f,%f];NOISED=0;NSUBDIV=[%d,%d];' % (1.,1.,nbElems, nbElems))
#-----
#m.export_to_pos('Gonzalez.pos')
#
mf   = gf.MeshFem(m, Dim)
mf.set_fem(gf.Fem(FEM+'('+str(Dim)+','+str(fem_order)+')'))
mf1D = gf.MeshFem(m, 1) # MeshFem 1D for normal and tangential components of the stress  tensor
mf1D.set_fem(gf.Fem(FEM+'('+str(Dim)+','+str(fem_order)+')'))
#mf1Dnc = gf.MeshFem(m, 1)
#mf1Dnc.set_classical_discontinuous_fem(0)

#  Integration method on each element
mim = gf.MeshIm(m, gf.Integ(Integ))
mimd1 = gf.MeshImData(mim, -1)
mimd4 = gf.MeshImData(mim, -1, [2,2]) # ou 4 ???
#-----------------------------------------------------------------
# Initialisation of plasticity
#-----------------------------------------------------------------
#convex_index = mim.convex_index()
#Nc = np.size(convex_index)
#print(Nc)
#d_Gamma = np.zeros(Nc)
#F_tr = np.eye(Dim*Nc,Dim*Nc) # Plastic part of the Deformation
#M_k = np.eye(Dim*Nc,Dim*Nc) # Mandel tensor
Y1_tr = 0.3
Y2_tr = 0.6
c_tr = 5.e8
md.add_initialized_data("c_tr", c_tr)
#md.add_internal_im_variable("dksi", mimd1)   # plastic multiplier increase in the current step
Identity = np.eye(Dim,Dim,dtype = 'float64')
#md.add_im_data("dGammak", mimd1)        # plastic multiplier value at previous time step
md.add_im_data("dGamma", mimd1)          # plastic multiplier value at current time step
#md.add_im_data("Fp0", mimd4)            # Plastic part of F
md.add_im_data("invF_tr", mimd4)         # Inverse of the Plastic part of F
md.add_im_data("Mk", mimd4)              # Mandel tensor at the previous step
md.add_im_data("M_trial", mimd4)         # Mandel tensor at the current step
M_dev = np.zeros((Dim,Dim)) # Deviator of the Mandel matrix at a Gauss point
#md.set_variable("F_tr_n",               # Plastic part of F at
#                 np.tile(Identity, # previous step.
#                 mimd4.nbpts()))        #
md.set_variable("invF_tr",               # Inverse of the Plastic part of F at
                 np.tile(Identity,  # previous step.
                 mimd4.nbpts()))         #
#md.set_variable("Mk",                   # Mandel Tensor at the
#                np.tile(Identity,  # previous step.
#                mimd4.nbpts()))         #
#md.set_variable("M_trial",                # Mandel Tensor at the
#                np.tile(Identity,  # current step.
#                mimd4.nbpts()))         #
md.set_variable("dGamma",np.tile([0],mimd1.nbpts()))
dGamma = np.zeros(mimd1.nbpts())
coeff_expm = np.zeros(mimd1.nbpts())
#F_tr_n  = np.zeros((mimd1.nbpts(),2,2),dtype = 'f') # Plastic part of F at time n
F_tr_n = np.zeros((mimd1.nbpts(),Dim,Dim),dtype = 'float64')
for i in range(mimd1.nbpts()):
    F_tr_n[i,:,:] = Identity
F_tr_n1 = np.zeros((mimd1.nbpts(),Dim,Dim),dtype = 'float64') # Plastic part of F at time n+1
invF_tr = np.copy(F_tr_n)                  # Inverse of the plastic part of F at time n+1
#md.add_macro("invCp0", "[[[1,0,0],[0,0,0],[0,0,0]],"+\
#                       " [[0,0,0],[0,1,0],[0,0,0]],"+\
#                       " [[0,0,0],[0,0,0],[0,0,1]],"+\
#                       " [[0,1,0],[1,0,0],[0,0,0]]].invCp0vec") #Vec4ToMat3x3
#md.add_macro("devlogbetr", "Deviator(Logm(F3d*invCp0*F3d'))")
#-----------------------------------------------------------------
#Dirichlet and Neumann Regions
#Boundaries creation
#-----------------------------------------------------------------
DIRICHLET_BOUNDARY = 1
NEUMANN_BOUNDARY = 2
CONTACT_BOUNDARY = 3

if exemple == 1 :
#   m.set_region(DIRICHLET_BOUNDARY,
#                m.outer_faces_with_direction([-1.,  0.], 0.001)) # Left side
    m.set_region(NEUMANN_BOUNDARY,
                 m.outer_faces_with_direction([ 0.,  1.], 0.001)) # Top side
    m.set_region(CONTACT_BOUNDARY,
                 m.outer_faces_with_direction([ 0., -1.], 0.001)) # Bottom side
#---
if exemple == 2 :
    fb1 = m.outer_faces_in_box([center[0]-rad_int-0.01, center[1]-rad_int-0.01], 
                               [center[0]+rad_int+0.01, center[1]+rad_int+0.01])  # Boundary of the hole
    fb2 = m.outer_faces_with_direction([0., -1.], np.pi/4.5) # Contact boundary of the wheel
    HOLE_BOUND=4
    m.set_region(HOLE_BOUND, fb1)
    m.set_region(CONTACT_BOUNDARY, fb2)
    m.region_subtract(CONTACT_BOUNDARY, HOLE_BOUND)
#---
if exemple == 3 :
    fb = m.outer_faces_in_box([center[0]-rad_ext-0.01,center[1]-rad_ext-0.01],
                              [center[0]+rad_ext+0.01,center[1]+rad_ext+0.01])
    m.set_region(CONTACT_BOUNDARY, fb)
#---
if exemple == 4 :
    m.set_region(NEUMANN_BOUNDARY,
                 m.outer_faces_with_direction([ 1., 0.], 0.001)) # Right side
    m.set_region(DIRICHLET_BOUNDARY,
                 m.outer_faces_with_direction([ 0., 1.], 0.001)) # Top side
    DIRICHLET_BOUNDARY2 = 4
    m.set_region(DIRICHLET_BOUNDARY2,
                 m.outer_faces_with_direction([ 0.,-1.], 0.001)) # Bottom side
    DIRICHLET_BOUNDARY3 = 5
    m.set_region(DIRICHLET_BOUNDARY3,
                 m.outer_faces_with_direction([ -1., 0.], 0.001)) # Left side
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
if exemple==1 : # Ici, on impose une force constante au cours du temps
#    #md.add_Dirichlet_condition_with_multipliers(mim, 'u', 0, DIRICHLET_BOUNDARY)
#    u_d = md.interpolation("[0., 0.]",mf)
#    md.add_initialized_fem_data("DirichletData", mf, u_d)
#    #md.add_Dirichlet_condition_with_multipliers(mim, 'u', mf,  DIRICHLET_BOUNDARY, "DirichletData")
#    md.add_Dirichlet_condition_with_simplification("u", DIRICHLET_BOUNDARY, "DirichletData")
#--
    ff=mf.eval('[0.,-100.*'+str(dt)+']')
#-------
if exemple==2 : # Ici, la force imposee est nulle
    ff=mf.eval('[0.,0.*'+str(dt)+']')
#-------
if exemple==3 : # Ici, la force imposee (pesanteur)
    ff=mf.eval('[-1000.,-2000.00]')
#-------
if exemple==4 :
    ff=mf.eval('[-20.,0.00]')  # Neumann condition
    u_d = md.interpolation("[0., 0.]",mf)
    md.add_initialized_fem_data("DirichletData", mf, u_d)
#    md.add_Dirichlet_condition_with_simplification("u", DIRICHLET_BOUNDARY, "DirichletData")
#    md.add_Dirichlet_condition_with_simplification("u", DIRICHLET_BOUNDARY2, "DirichletData")
#    md.add_Dirichlet_condition_with_simplification("u", DIRICHLET_BOUNDARY3, "DirichletData")
    H1 = np.array([[0,0],[0,1]],dtype="float64")
    md.add_initialized_data("H1",H1)
    md.add_generalized_Dirichlet_condition_with_penalization(mim, "u", pen_factor,
                                    DIRICHLET_BOUNDARY, "DirichletData", "H1")
    md.add_generalized_Dirichlet_condition_with_penalization(mim, "u", pen_factor,
                                    DIRICHLET_BOUNDARY2, "DirichletData", "H1")
    
    H2 = np.array([[1,0],[0,0]],dtype="float64")
    md.add_initialized_data("H2",H2)
    md.add_generalized_Dirichlet_condition_with_penalization(mim, "u", pen_factor,
                                    DIRICHLET_BOUNDARY3, "DirichletData", "H2")
#-----------------------------------------------------------------
#
md.add_initialized_fem_data('f', mf, ff)
#md.add_linear_term(mim, "-f.Test_u", NEUMANN_BOUNDARY)
md.add_source_term_brick(mim, "u", "f", NEUMANN_BOUNDARY)
if exemple==3 :
    md.add_source_term_brick(mim, "u", "f") # pesanteur (force interne)
#-------
if exemple==1 : # Dans cet exemple, la vitesse initiale est nulle
#   v0 = mf.eval("[0.,-x*x/16.]")
    v0 = mf.eval("[0.,0]")
if exemple==2 : # Dans cet exemple, la vitesse initiale est imposée
    v0 = mf.eval("[0.5,-1]")
#   v0 = mf.eval("[-x-0.2,0]") # Exemple pour faire de la deformation independemment du contact
if exemple==3 : # Dans cet exemple, la vitesse initiale est parallele au support
#   v0 = mf.eval("[0.1,-0.050]")
    v0 = mf.eval("[0.3,-0.150]")
if exemple==4 : # Dans cet exemple, la vitesse initiale est nulle
    v0 = mf.eval("[0.,0]")
md.add_initialized_fem_data("v_old", mf, v0)
V = md.variable('v_old') # Pas terrible si vitesse initiale non nulle
V_old = V
md.add_initialized_fem_data("v", mf, v0)    # velocities
#-------
md.add_initialized_data("K", E/(3.*(1.-nu)))
md.add_initialized_data("G", E/(2.*(1.+nu)))
#
md.add_macro('F','(Id(meshdim)+Grad_u).invF_tr') # Elastic part of F
md.add_macro('C','Right_Cauchy_Green(F)')        # Elastic Right Cauchy Green tensor
#md.add_macro('C',"(Id(meshdim)+Grad_u+Grad_u')") # Version linéarisée
#
md.add_macro('F_old', '(Id(meshdim)+Grad_u_old).invF_tr')
md.add_macro('C_old','Right_Cauchy_Green(F_old)')
#md.add_macro('C_old',"(Id(meshdim)+Grad_u_old+Grad_u_old')") # Version linéarisée
#
md.add_macro('F_avg', '0.5*(F+F_old)')
md.add_macro('C_avg', '0.5*(C+C_old)')
#
# Loi de ???
#md.add_macro("SS", "S(C_avg) + ( (2*W(C)-2*W(C_old))/(Norm(C-C_old)+eps) - S(C_avg):Normalized(C-C_old) )*Normalized(C-C_old)")
#md.add_macro('S(CC)', 'K/2.*log(Det(CC))*Inv(CC)+G*pow(Det(CC),-4./3.)*(Id(meshdim)-1./3.*Trace(CC)*Inv(CC))')
#md.add_macro('W(CC)', 'K/8.*sqr(log(Det(CC)))+G/2.*(pow(Det(CC),-4./3.)*(Trace(CC))-3.)')
# Loi de Kirchoff-Saint Venant
#md.add_macro("E(CC)", "0.5*(CC-Id(meshdim))")
#md.add_macro("W(CC)", "0.5*K*pow(Trace(E(CC)),2) + G*(E(CC):E(CC))") #"G*Matrix_i2(E(CC))")
#md.add_macro("S(CC)", "K*Trace(E(CC))*Id(meshdim)+2*G*E(CC)")
#
# Loi de Ogden/Ciarlet-Geymonat
clambda = E*nu/((1.+nu)*(1.-2.*nu)) # First Lame coefficient (N/cm^2)
cmu = E/(2.*(1.+nu))               # Second Lame coefficient (N/cm^2)
clambda = E*nu/((1.+nu)*(1.-nu))
a = 0.5*(max(0.,0.5*cmu-0.25*clambda)+0.5*cmu)
c1 = 0.5*cmu-a
c2 = 0.25*clambda-0.5*cmu+a
c3 = 0.5*cmu+0.25*clambda
if exemple==4 :
    c1 = 0.5e6
    c2 = 0.5e4
    c3 = 0.35e6
    a = 0.
para = [a,c1,c2,c3]
md.add_initialized_data("c_", para)
# Pour les rendre variable sur le maillage (pour plus tard, inutile pour le moment)
#mimd4 = gf.MeshImData(mim, -1, [4])
#md.add_im_data('c_', mimd4)
#md.set_variable('c_',np.tile(para, mimd4.nbpts()))             # Parameters
#
# Energie elastique et contrainte Ciarlet-Geymonat en 2D
md.add_macro("W(CC)", "c_(1)*(Trace(CC)-2.)+c_(2)*(Trace(CC)+Det(CC)-3.)+c_(3)*(Det(CC)-1)-(c_(1)+2.*c_(2)+c_(3))*log(Det(CC))")
md.add_macro("S(CC)", "2.*( (c_(1)+c_(2))*Id(meshdim) + (c_(2)+c_(3))*Det(CC)*Inv(CC) - (c_(1)+2.*c_(2)+c_(3))*Inv(CC) )")
#
# Gonzalez:
md.add_initialized_data("eps", 1.e-8)
md.add_macro("SS", "S(C_avg) + (2.*(W(C)-W(C_old))-S(C_avg):(C-C_old))/max(Norm_sqr(C-C_old),eps)*(C-C_old)")
# Implicit Euler:
#md.add_macro("SS", "S(C)")
# Explicit Euler
#md.add_macro("SS", "S(C_old)")
#
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

#md.add_nodal_contact_with_rigid_obstacle_brick(mim, varname_u, multname_n, multname_t=None, *args)
#
mf.export_to_vtu('./results/Displacement_%05i.vtu' %0, mf, U, 'Displacement',
                 V,'Velocity')
#--------
# Initial energy
energy = gf.asm_generic(mim, 0, "0.5*rho*Norm_sqr(v_old) + W(C)", -1, md)
print("==> mech energy = ", energy)
fic.write(f"Time step: {0.0: 9.5f}, Energy: { energy:.17f}\n")
#
#Contact_length = np.size(mf1D.dof_on_region(CONTACT_BOUNDARY))
Contact_U_dof = mf.dof_on_region(CONTACT_BOUNDARY)
Contact_S_dof = mf1D.dof_on_region(CONTACT_BOUNDARY)
Contact_length = np.size(Contact_S_dof)
print('Nombre de points en contact potentiel : '+str(Contact_length))
print('Points potentiellement en contact : '+str(Contact_S_dof))
Contact_force_dof = mf.dof_on_region(NEUMANN_BOUNDARY)
print('DOF Force : '+str(Contact_force_dof))

#--------
# For frictional contact
#--------
Points = mf1D.basic_dof_nodes() #: il ne faudrait avoir que ceux sur la frontière....
Points = Points[:,Contact_S_dof]
#
#Normale = np.array([ 0.,-1.], dtype='float') # Outward normal! Pourrait dépendre du point facilement
norm_n = np.sqrt(1+a_plane*a_plane)
Normale = np.array([ a_plane,-1.], dtype='float')/norm_n # Outward normal
y_vect = np.array([0,y_plane], dtype='float')
Matr_normale = Normale.reshape(-1,1)@Normale.reshape(1,-1)
Proj_tangent = np.eye(Dim)-Matr_normale
Matr_normale = pen_factor*Matr_normale
Matr_tangent = pen_factor*Proj_tangent
#--------

# Time loop
for timeStepIndex,timeStep in enumerate(np.arange(0.,T_max+dt,dt)):
    print('Time step: %9.5f' % timeStep)
    dU = dt*V #0.*U #   Une initialisation pour pouvoir calculer la contrainte à la première itération : correct ???
    U = U + dU
    md.set_variable("u",U)
#=====
#  For the plasticity
#=====
    dGamma[:] = 0.
#=======
# Frictional contact:
#=======
#    dF_fric = np.zeros(NTDL)
    F_fric  = np.zeros(NTDL) # Nodal Friction forces
    G       = np.zeros(NTDL) # Splipping Contact Forces
#    Stress_points = np.zeros(NTDL) # En fait, ce serait plutot la force nodale de contact => c'est F_fric
    if timeStepIndex > 0 :
        force = md.add_source_term_brick(mim, "u", "f", NEUMANN_BOUNDARY)
    dU_old = 0.*U # Variable pour faire un test qui évite les cycles sur 2 itérations
#=====
# Debut de la boucle de Newton
#=====
    for Newton_loop in range(Newton_ite_max): # Ameliorer cette boucle : on estime qu'il faut avoir convergé en 10 itérations max !!!
#
#=====
#  For the plasticity
#=====
# Mandel Tensor : M = C_avg*S(C_avg)
        md.set_variable("Mk",md.interpolation("C_avg.SS",mimd4,-1)) #C_avg.S(C_avg) ???   C*S(C) ???
        M = np.reshape(md.variable("Mk"),(mimd4.nbpts(),Dim,Dim))
# Mandel trial (Computed with GetFem, possible to compute it by myself, using lines with #!!)
#  !! It is sqrt(1.5) in 3D !!!
        md.set_variable("M_trial",md.interpolation("Mk+c_tr*dGamma*Deviator(Mk)/"+
                                                   "max(Norm(Deviator(Mk)),eps)/sqrt(0.5)",mimd4,-1)) #
        M_trial = np.reshape(md.variable("M_trial"),(mimd4.nbpts(),Dim,Dim))
#!!        M_trial = np.copy(M)
        if Newton_loop == 0:
            Ma_trial = M # To compute error at the first iteration
#==
# Von Mises test for plasticity (2D !!!)
        for ig in range(mimd1.nbpts()):  # Loop on Gauss points
            # Deviator of M => M_dev
            M_dev,Norm_H_M_dev = compute_dev(M[ig,:,:])
#!!            M_trial[ig,:,:] =  M_trial[ig,:,:] +c_tr*dGamma[ig]*M_dev[:,:]/max(Norm_H_M_dev,eps)
            M_trial_dev,Norm_H_M_trial = compute_dev(M_trial[ig,:,:])
#            print(Norm_H_M_trial)
            if Norm_H_M_trial > Y1_tr :
                Y_tr = Y1_tr
                if Norm_H_M_trial > Y2_tr :
                    Y_tr = Y2_tr
                    #print("Plasticity 2")
            #    else:
            #        print("Plasticity 1")
                dGamma[ig] = (Y_tr-Norm_H_M_dev)/c_tr
                M_trial[ig,:,:] = M[ig,:,:] + (Y_tr-Norm_H_M_dev)*M_dev[:,:]/Norm_H_M_dev
                M_trial_dev,Norm_H_M_trial = compute_dev(M_trial[ig,:,:])
#
                F_tr_n1[ig,:,:] = np.matmul(F_tr_n[ig,:,:],expm(dGamma[ig]/Norm_H_M_trial*M_trial_dev))
            else:
            #    print('Elasticity')
                dGamma[ig] = 0.
                coeff_expm[ig] = 0.
                F_tr_n1[ig,:,:] = F_tr_n[ig,:,:]                
            # Computation of the inverse of F
            invF_tr[ig,:,:] = np.linalg.inv(F_tr_n1[ig,:,:])
            # Pour test  : annule tout ce qui précède en plasticité !!!
            #invF_tr[ig,:,:] = Identity
            #print(invF_tr[ig,:,:])
#==
#        md.set_variable("M_trial",np.reshape(M_trial,Dim*Dim*mimd4.nbpts()))
        md.set_variable("invF_tr",np.reshape(invF_tr,Dim*Dim*mimd4.nbpts()))
        md.set_variable("dGamma",dGamma)
#        Error_tr = gf.compute(mf, md.interpolation("M_trial-Mk",mimd4), 'L2 norm', mim)
        Error_tr = np.linalg.norm(np.reshape(M_trial-Ma_trial,Dim*Dim*mimd4.nbpts())) # A ameliorer !!!
        Error_tr = Error_tr / np.linalg.norm(np.reshape(M_trial,Dim*Dim*mimd4.nbpts()))
        Ma_trial = np.copy(M_trial)
#        Error_tr = Error_tr/np.linalg.norm(np.reshape(M,Dim*Dim*mimd4.nbpts()))
#        print('Erreur :'+str(Error_tr)+'   '+str(np.linalg.norm(np.reshape(M,Dim*Dim*mimd4.nbpts()))))
#        Error_tr = 0.000
#=======
# Frictional contact:
#=======
        nb_active = 0
        nb_slip_active = 0
        nb_stick_active = 0
        ind_active_points=np.array([], dtype = 'int')
        ind_stick_active_points=np.array([], dtype = 'int')
        ind_slip_active_points=np.array([], dtype = 'int')
#
        md.assembly()
        K_tangent = md.tangent_matrix()
#
        F  = md.rhs()
#--- Pas terrrible, mais voir comment améliorer cela... Sert à calculer les forces de frottements à la fin.
#        K1 = md.tangent_matrix() # Voir comment on pourrait faire une copie de l'objet K_tangent...
        K1 = gf.Spmat('copy' , K_tangent) ## A essayer...
#        K1 = gf.Spmat('add' , K_tangent , Mass_matrix)
        F1 = np.copy(F)  #md.rhs()
#
        G = np.zeros(NTDL)
#=====
        for i in range(Contact_length):
            ind_u = Contact_U_dof[Dim*i] # DOF du noeud i en contact
#=== Contact
            # Gap au milieux du pas de temps (schéma de Gonzalez)
            Gap = Normale@(y_vect[0:Dim] - Points[0:Dim,i]
                           - 0.5*(U_old[ind_u:ind_u+Dim]+U[ind_u:ind_u+Dim]))
#
            if Gap <= 1.e-9 :
#                print('indice : '+str(ind_u))
###                Normal_stress_points = Normale@Stress_points[ind_u:ind_u+2]
                Normal_stress_points = Normale@F_fric[ind_u:ind_u+Dim]
#                Active = ( Normal_stress_points + C_nu*Gap <= 0.) # Ce n'est pas Gap, mais la vitesse, à améliorer...
                V_normal = Normale@V[ind_u:ind_u+Dim]
                Active = ( Normal_stress_points - C_nu*V_normal <= 0.e-9)
#                print('Test Contact : '+str(Normal_stress_points - C_nu*V_normal)+'  +++  '+str(Normal_stress_points)+'  +++  '+str(C_nu*V_normal) +'  +++  ')
                if Active:
#                    print('----- '+str(ind_u+1))
#                    print('***  '+str(-Normal_stress_points)+'  ***  '+str(- C_nu*Gap))
                    ind_active_points = np.append(ind_active_points,ind_u+1)  # voir comment modifer avec la normale...
                    # Ici, on impose la vitesse à la fin du pas de temps égale à zéro : pas bon...
                    #delta_impose_deplac = -Normale@(U_old[ind_u:ind_u+2]-U[ind_u:ind_u+2]+0.5*dt*V_old[ind_u:ind_u+2])
                    # Ici, on impose la vitesse au milieu du pas de temps nulle
                    delta_impose_deplac = Normale@(U_old[ind_u:ind_u+Dim]-U[ind_u:ind_u+Dim])
#
                    K_tangent.add([range(ind_u,ind_u+Dim)],[range(ind_u,ind_u+Dim)],Matr_normale)
 ###                   F[ind_u+1] += -pen_factor*delta_impose_deplac
                    F[ind_u:ind_u+Dim] += pen_factor*delta_impose_deplac*Normale
                    nb_active += 1
#=== Friction
###                    Tangent_stress_points = Tangent@Stress_points[ind_u:ind_u+2]
#-#-                    Tangent_stress_points = Tangent@F_fric[ind_u:ind_u+2]
                    Tangent_stress_points = F_fric[ind_u:ind_u+Dim] - Normal_stress_points*Normale
#                    print(Tangent_stress_points)
#-#-                    V_tau = Tangent@V[ind_u:ind_u+2]
                    V_tau = V[ind_u:ind_u+Dim] - V_normal*Normale
#                    print(V_tau)
                    Gliss_dir = Tangent_stress_points+C_tau*V_tau
#-#-                    Active_fric = (abs(Gliss_dir)+mu*Normal_stress_points>0)
                    Active_fric = (np.linalg.norm(Gliss_dir,2)+mu*Normal_stress_points>0)
                    if Active_fric:  # On impose une force tangentielle
#                        print(f"Glisse : {ind_u:5d}")
                        ind_slip_active_points = np.append(ind_slip_active_points,ind_u)
                        G[ind_u:ind_u+Dim] = mu*Normal_stress_points*Gliss_dir/max(np.linalg.norm(Gliss_dir,2),eps)
#-#-!!!                        G[ind_u:ind_u+2] = mu*Normal_stress_points*Gliss_dir/max(abs(Gliss_dir),eps)
                        nb_slip_active += 1
                    else: # On impose une vitesse tangentielle nulle (formule approchée)
#                        print(f'Adhere : {ind_u:5d}')
                        ind_stick_active_points = np.append(ind_stick_active_points,ind_u)
                        # En vitesse, a la fin du pas de temps
                        #delta_impose_deplac = Tangent@(U_old[ind_u:ind_u+2]-U[ind_u:ind_u+2]+0.5*dt*V_old[ind_u:ind_u+2])
                        # En vitesse au milieux du pas de temps
#-#-                        delta_impose_deplac = Tangent@(U_old[ind_u:ind_u+2]-U[ind_u:ind_u+2])
                        delta_impose_deplac = Matr_tangent@(U_old[ind_u:ind_u+Dim]-U[ind_u:ind_u+Dim])
                        # En déplacement (deplacement tangentiel nul à l'étape suivante):
#                        delta_impose_deplac = -Tangent@U[ind_u:ind_u+2]
                        F[ind_u:ind_u+Dim] += delta_impose_deplac  #-#-pen_factor*delta_impose_deplac #-#-*Tangent
##===
##                        K_tangent.add([ind_u],[ind_u],pen_factor) # Voir comment adapter ceci à une tangente quelconque
                        K_tangent.add([range(ind_u,ind_u+Dim)],[range(ind_u,ind_u+Dim)],Matr_tangent)
###                        F[ind_u] += pen_factor*delta_impose_deplac
                        nb_stick_active += 1
#             if (ind_u==2) : print('***  '+str(-Normal_stress_points)+'  ***  '+str(- C_nu*Gap))
#=== End Friction
#=====
# Resolution (contact frottant et plasticité)
#=====
#        print('second membre : '+str(bb))
#        md.solve('lsolver','mumps','max_res',1E-14,'max_iter',md.nbdof(),'very_noisy')
#        dU = md.variable('u') - U
#        U = md.variable('u')
        dU = linsolve(K_tangent, F-G)
        U = U + dU
##        print('U = '+str(U))
# ---
        md.set_variable("u",U)
# ---
# Calcul de la force de friction pour l'iteration suivante : on pourrait ne pas mettre -G, pour ne pas remettre +G au debut de la boucle...
###        F_fric = F1 - G - K1.mult(dU)
        F_fric = F1 - K1.mult(dU)
#        print(F_fric[ind_active_points])
#=====
# Updated of the velocity
#   Gonzalez:
        V = 2.*(U-U_old)/dt-V_old
###        F_fric_i = F_fric # update the previous friction force
#        md.set_variable("v",V)
##        print('Vitesse :   '+str(V[ind_active_points])) # Doit etre nul
##        print('deplacement a imposer :  '+str(str(U[ind_active_points])+'  '+str(U_old[ind_active_points]+0.5*dt*V_old[ind_active_points])))
#   Implicit Euler;
#        V = (U-U_old)/dt
#        V = md.variable('v')
#        U = U+dt*V
#=====
# Test d'arret sur le contact :
#=====
        if (nb_active>0) or (Newton_loop > 0):
#           print('Active set : '+str(ind_active_points))#+'   '+str(bb)+'   '+str(L_fric))
#           print('Stick      : '+str(ind_stick_active_points))
#           print('Slip       : '+str(ind_slip_active_points))
#           print(f"Energy: { energy:.7f}\n")
            L2error = gf.compute(mf, dU, 'L2 norm', mim)
            if Newton_loop > 1 :
                L2error = min(L2error,gf.compute(mf, dU_old+ dU, 'L2 norm', mim)) # pour eviter les cycles sur 2 itérations
            print('-----> Newton iteration: '+str(Newton_loop)+
                  ', L2 error: '+str(L2error+Error_tr))
            if  L2error+Error_tr < eps_Newton : break
#=====
        if  nb_active ==0 and Newton_loop == 0 :
            print('-----> Newton iteration: '+str(Newton_loop)+
                  ', No contact')
            L2error = gf.compute(mf, dU, 'L2 norm', mim)
            print('-----> L2 error: '+str(L2error+Error_tr))
            if L2error+Error_tr < eps_Newton : break
        dU_old = dU
#=====
# Fin de la boucle de Newton
#=====
    F_tr_n = np.copy(F_tr_n1) # update de Plastic part of F
#=====
    print('Active set : '+str(ind_active_points))#+'   '+str(bb)+'   '+str(L_fric))
    print('Stick      : '+str(ind_stick_active_points))
    #Post processing of the solution
    if visu & timeStepIndex%delta_visu==0 :
        mf.export_to_vtu('./results/Displacement_%05i.vtu' %(timeStepIndex+1),
                         mf, U, 'Displacement',V,'Velocity')
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
    fic.write(f"Time step: {timeStep: 9.5f}, Energy: { energy:.17f}\n")
#    print(f"Nb active contact points : {nb_active: 5d}")
#    print(f"Nb active slip points    : {nb_slip_active: 5d}")
#    print(f"Nb active stick points   : {nb_stick_active: 5d}")
#
#    md.shift_variables_for_time_integration()
#=====
# Fin de la boucle de pas de temps
#=====
fic.close()

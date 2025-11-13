# -*- coding: utf-8 -*-

# ------------- A CHANGER ------------- #
# Problem considered:
# div(\sigma(u)) + f_volumic = \rho d(du/dt)/dt in \Omega
# u=u_d on \Gamma_D
# \sigma(u).n=f_surface on \Gamma_N
# u(t=0)=0
# du/dt(t=0)=0
# ------------- A CHANGER ------------- #

import numpy as np 
import getfem as gf
import sys,os

def linsolve(M, B): # Call Superlu to solve a sparse linear system
    return (((gf.linsolve_superlu(M, B))[0]).T)[0]
#-----
os.system('rm *.vtu')
np.set_printoptions(threshold=sys.maxsize)
fic=open('energy.txt','w')

gf.util('trace level', 1)

md = gf.Model("real")

E =   1.e6                    # Young's modulus [Pa]
nu =  0.2                     # Poisson's ratio [-]
#rho = 500.0                  # density
rho = 1000.1
#-----
# Pour le contact
mu = 0.5
#--
C_nu  = 1.e+06
C_tau = 1.e+03
eps = 1.e-10 # Pour éviter de diviser par zero... on ajoute epsilon au denominateur s'il est petit
pen_factor = 1.e+12 # Pour imposer les CL de type Dirichlet
#-----
z_plane = -0.1 # Hauteur du socle rigide 
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#Time related parameters
T_max = 2. # 10.e-3 # 400.e-4 #0.25 # 1.0 #  0.4 # 
dt =1.e-4
steps = T_max/dt
#-----------------------------------------------------------------
# Visu parameters
visu = True
delta_visu = 100 # Nombre d'itérations entre 2 images
#-----------------------------------------------------------------
# Test d'arret
eps_Newton = 1.e-7
#-----------------------------------------------------------------

#-----
#Elements finis
#-----
fem_order = 1
Dim = 3                            # Dimension du problème (toujours 3 ici)
FEM   = 'FEM_PK'

nbElems = 4.
meshSize = 1./nbElems

# traction_force = 10 #Surfacic load on the Neumann boundary of the body in the vertical direction
# surfaceLoad=[0, 0, -traction_force]
# volumeLoad=[0, 0, 0]
#-----
m=gf.Mesh('regular simplices', np.arange(0,4+meshSize,meshSize), 
        np.arange(0,1+meshSize,meshSize),np.arange(0,1+meshSize,meshSize))
# Integration method
mim=gf.MeshIm(m, gf.Integ('IM_TETRAHEDRON(8)'))
#-----
#m = gf.Mesh('cartesian', np.arange(0,1+meshSize,meshSize), 
#        np.arange(0,4+meshSize,meshSize),np.arange(0,1+meshSize,meshSize))
# The mesh_im stores the integration methods for each tetrahedron
#mim = gf.MeshIm(m, gf.Integ('IM_GAUSS_PARALLELEPIPED(3,4)'))
#mf.set_fem(gf.Fem('FEM_QK(3,1)'))
#-----
# Create a MeshFem for u field of dimension 3 (i.e. a vector field)
mf = gf.MeshFem(m,Dim) # Displacement
# Choice of the finite element method
mf.set_classical_fem(fem_order)
mf2D = gf.MeshFem(m, 1) # MeshFem 2D for normal and tangential components of the stress  tensor
mf2D.set_fem(gf.Fem(FEM+'('+str(Dim)+','+str(fem_order)+')'))

#Dirichlet and Neumann conditions
#Boundaries creation
DIRICHLET_BOUNDARY = 1
NEUMANN_BOUNDARY = 2
CONTACT_BOUNDARY = 3

m.set_region(DIRICHLET_BOUNDARY,
             m.outer_faces_with_direction([-1., 0., 0.], 0.001)) # Left side
m.set_region(NEUMANN_BOUNDARY,
             m.outer_faces_with_direction([1.,  0., 0.], 0.001)) # Right side
m.set_region(CONTACT_BOUNDARY,
             m.outer_faces_with_direction([0.,  0.,-1.], 0.001)) # Bottom side

#-----------------------------------------------------------------
md.add_fem_variable("u", mf)    # displacements

# Initialisation
u0 = md.interpolation("[0.,0.,0.]", mf)
md.add_initialized_fem_data("u_old", mf, u0)
md.set_variable('u', u0)
#
U = md.variable('u')
U_old = md.variable('u_old')  #np.copy(U)
#
v0 = mf.eval("[0.,0.,-x*x/16.]")
#v0 = md.interpolation("[0.,0.,0.]", mf)
md.add_initialized_fem_data("v_old", mf, v0)
V = md.variable('v_old') # Pas terrible si vitesse initiale non nulle 
V_old = V
md.add_initialized_fem_data("v", mf, v0)    # velocities
#
NTDL = np.size(U)
#-----------------------------------------------------------------
# boundary conditions
#-----------------------------------------------------------------
ff=mf.eval('[0.,0.,0.*'+str(dt)+']')
md.add_initialized_fem_data('f', mf, ff)

u_d = md.interpolation("[0., 0., 0.]",mf)
md.add_initialized_fem_data("DirichletData", mf, u_d)
#md.add_Dirichlet_condition_with_multipliers(mim, 'u', mf,  DIRICHLET_BOUNDARY, "DirichletData")
#md.add_Dirichlet_condition_with_simplification("u", DIRICHLET_BOUNDARY, "DirichletData")
md.add_Dirichlet_condition_with_penalization(mim,"u", 1.e12,DIRICHLET_BOUNDARY, "DirichletData")
#-----------------------------------------------------------------
md.add_initialized_data("K", E/(3.*(1.-2.*nu)))
md.add_initialized_data("G", E/(2.*(1.+nu)))

md.add_macro('F',  '(Id(meshdim)+Grad_u)')
md.add_macro('C','Right_Cauchy_Green(F)')
#
md.add_macro('F_old', '(Id(meshdim)+Grad_u_old)')
md.add_macro('C_old','Right_Cauchy_Green(F_old)')
#
md.add_macro('F_avg', '0.5*(F+F_old)')
md.add_macro('C_avg', '0.5*(C+C_old)')

# Loi de ???
#md.add_macro('S(CC)', 'K/2.*log(Det(CC))*Inv(CC)+G*pow(Det(CC),-1./3.)*(Id(meshdim)-1./3.*Trace(CC)*Inv(CC))')
#md.add_macro('W(CC)', 'K/8.*sqr(log(Det(CC)))+G/2.*(pow(Det(CC),-1./3.)*(Trace(CC))-3.)')
# Loi de Kirchoff-Saint Venant
#md.add_macro("E(CC)", "0.5*(CC-Id(meshdim))")
#md.add_macro("W(CC)", "K/2*pow(Trace(E(CC)),2) + G*(E(CC):E(CC))")
#md.add_macro("S(CC)", "K*Trace(E(CC))*Id(meshdim)+2*G*E(CC)")
# Loi de Ogden/Ciarlet-Geymonat
clambda = E*nu/((1+nu)*(1-2*nu)) # First Lame coefficient (N/cm^2)
cmu = E/(2*(1+nu))               # Second Lame coefficient (N/cm^2)
a = 0.5*(max(0.,0.5*cmu-0.25*clambda)+0.5*cmu)
c1 = 0.5*cmu-a
c2 = 0.25*clambda-0.5*cmu+a
c3 = 0.5*cmu+0.25*clambda
para = [a,c1,c2,c3]
md.add_initialized_data("c_", para)
# Pour les rendre variable sur le maillage
#mimd4 = gf.MeshImData(mim, -1, [4])
#md.add_im_data('c_', mimd4)
#md.set_variable('c_',np.tile(para, mimd4.nbpts()))             # Parameters
#
#md.add_macro("W(CC)", "c_(1)*(Trace(CC)-3)+c_(2)*(Matrix_i2(CC)-3)+c_(3)*(Det(CC)-1)-(c_(1)+2*c_(2)+c_(3))*log(Det(CC))")
#md.add_macro("W(CC)", "c_(1)*Trace(CC)+c_(2)*Matrix_i2(CC)+c_(3)*Det(CC)-(c_(1)+2*c_(2)+c_(3))*log(Det(CC))")
#md.add_macro('S(CC)', "2*((c_(1)+c_(2)*Trace(CC))*CC - c_(2)*(CC.CC) +(c_(3)*Det(CC)- (c_(1)+2*c_(2)+c_(3)))*Id(meshdim))")
# Energie elastique et contrainte Ciarlet-Geymonat en 3D
md.add_macro("W(CC)", "c_(1)*(Trace(CC)-3.)+c_(2)*(Matrix_i2(CC)-3.)+c_(3)*(Det(CC)-1)-(c_(1)+2.*c_(2)+c_(3))*log(Det(CC))")
md.add_macro("S(CC)", "2.*((c_(1)+c_(2)*Trace(CC))*Id(meshdim) - c_(2)*CC +(c_(3)*Det(CC)- (c_(1)+2*c_(2)+c_(3)))*Inv(CC))")
#-----------------------------------------------------------------
#
# Gonzalez:
md.add_initialized_data("eps", 1.e-15)
#md.add_macro("SS", "S(C_avg) + (2.*(W(C)-W(C_old))-S(C_avg):(C-C_old))/max(Norm_sqr(C-C_old),eps)*(C-C_old)")
#md.add_macro("SS", "S(C_avg) + (2.*(W(C)-W(C_old))-S(C_avg):(C-C_old))/(Norm_sqr(C-C_old)+eps)*(C-C_old)")
md.add_macro("SS", "S(C_avg) + ( 2.*(W(C)-W(C_old))/max(Norm(C-C_old),eps) - S(C_avg):Normalized(C-C_old) )*Normalized(C-C_old)")
# Implicit Euler:
#md.add_macro("SS", "S(C)")
# Explicit Euler
#md.add_macro("SS", "S(C_old)")
 
# elasticity + inertia
md.add_initialized_data("rho", rho)
md.add_initialized_data("dt", dt)
###md.add_linear_term(mim, "(u-u_old-0.5*dt*(v+v_old)).Test_u")
###md.add_linear_term(mim, "rho/dt*(v-v_old).Test_v")
###md.add_nonlinear_term(mim, "(F_old*SS):Grad_Test_v")
# Gonzalez:
md.add_linear_term(mim, "2.*rho/dt*((u-u_old)/dt-v_old).Test_u")
#md.add_linear_term(mim, "(2.*rho*u/dt/dt).Test_u")
# Implicit Euler: -
#md.add_linear_term(mim, "rho/dt*((u-u_old)/dt-v_old).Test_u")
#
md.add_nonlinear_term(mim, "(F_avg.SS):Grad_Test_u")
#md.add_linear_term(mim, "(F_old.S(C_old)):Grad_Test_u") # To try explicit
#-------
#md.add_initialized_fem_data('f', mf, md.interpolation('[0., 0.,-1.]', mf))
#md.add_linear_term(mim, "-f.Test_u", NEUMANN_BOUNDARY)
#md.add_source_term_brick(mim, "u", "f", NEUMANN_BOUNDARY)
#-------
mf.export_to_vtu('./results/Displacement_%05i.vtu' %0, mf, U, 'Displacement',V,'Velocity')
#-----------------------------------------------------------------
# Initial energy
#-----------------------------------------------------------------
energy = gf.asm_generic(mim, 0, "0.5*rho*Norm_sqr(v_old) + W(C)", -1, md)
print("==> mech energy = ", energy)
fic.write(f"Time step: {0.0: 9.5f}, Energy: { energy:.17f}\n")
#
Contact_U_dof = mf.dof_on_region(CONTACT_BOUNDARY)
Contact_S_dof = mf2D.dof_on_region(CONTACT_BOUNDARY)
Contact_length = np.size(Contact_S_dof)
print('Nombre de points en contact potentiel : '+str(Contact_length))
print('Points potentiellement en contact : '+str(Contact_S_dof))
Contact_force_dof = mf.dof_on_region(NEUMANN_BOUNDARY)
print('DOF Force : '+str(Contact_force_dof))
print('DOF Deplacement'+str())
#--------
Points = mf2D.basic_dof_nodes() #: il ne faudrait avoir que ceux sur la frontière....
Points = Points[:,Contact_S_dof]
#
Normale = np.array([ 0.,0.,-1.], dtype='float') # Outward normal! Pourrait dépendre du point facilement
#Tangent = np.array([+1., 0.], dtype='float')
#n1 = Normale[0] ; n2 = Normale[1] ;
#Matr_normale = np.array([[n1*n1,n1*n2],[n1*n2,n2*n2]], dtype='float') # A Ameliorer, à mettre sous forme vectorielle, c'est n.n^t
Matr_normale = Normale.reshape(-1,1)@Normale.reshape(1,-1)
Proj_tangent = np.eye(Dim)-Matr_normale
Matr_normale = pen_factor*Matr_normale
Matr_tangent = pen_factor*Proj_tangent
#--------

#-----------------------------------------------------------------
# Time loop
#-----------------------------------------------------------------
for timeStepIndex,timeStep in enumerate(np.arange(0.,T_max+dt,dt)):
    print('Time step: %9.5f' % timeStep)
    dU = dt*V #0.*U #   Une initialisation pour pouvoir calculer la contrainte à la première itération : correct ???
    U = U + dU
    md.set_variable("u",U)
#    dF_fric = np.zeros(NTDL)
    F_fric  = np.zeros(NTDL) # Nodal Friction forces
    G       = np.zeros(NTDL) # Splipping Contact Forces
#    Stress_points = np.zeros(NTDL) # En fait, ce serait plutot la force nodale de contact => c'est F_fric
    if (timeStepIndex>0) :
       force = md.add_source_term_brick(mim, "u", "f", NEUMANN_BOUNDARY)
    dU_old = 0.*U # Variable pour faire un test qui évite les cycles sur 2 itérations
#===== 
# Debut de la boucle de Newton
#=====    
    for contact_Newton_boucle in range(10): # Ameliorer cette boucle : on estime qu'il faut avoir convergé en 10 itérations max !!!
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
#        Stress_points[0:NTDL] = F_fric[0:NTDL] + G[0:NTDL] # Forces imposées par tout les statuts de contact : 
                                                           # on pourrait améliorer avec ind_active, ind_slip_active et ind_stick_active
                                                           # et ne pas ajouter G ici et ne pas l'enlever à la fin de la boucle
        G = np.zeros(NTDL)
#=====
        for i in range(Contact_length):
            ind_u = Contact_U_dof[Dim*i] # DOF du noeud i en contact
#=== Contact
            # Gap au milieux du pas de temps (schéma de Gonzalez)
            Gap = Normale@(z_plane - Points[0:Dim,i] - 0.5*(U_old[ind_u:ind_u+Dim]+U[ind_u:ind_u+Dim])) 
            if (Gap<=1.e-9):
                Normal_stress_points = Normale@F_fric[ind_u:ind_u+Dim]
                V_normal = Normale@V[ind_u:ind_u+Dim]
                Active = ( Normal_stress_points - C_nu*V_normal <= 0.e-9)
#                print('Test Contact : '+str(Normal_stress_points - C_nu*V_normal)+'  +++  '+str(Normal_stress_points)+'  +++  '+str(C_nu*V_normal) +'  +++  ')
                if Active:
                    # Ici, on impose la vitesse au milieu du pas de temps nulle                    
                    delta_impose_deplac = Normale@(U_old[ind_u:ind_u+Dim]-U[ind_u:ind_u+Dim])
                    K_tangent.add([range(ind_u,ind_u+Dim)],[range(ind_u,ind_u+Dim)],Matr_normale)
                    F[ind_u:ind_u+Dim] += pen_factor*delta_impose_deplac*Normale
                    nb_active += 1
#=== Friction
                    Tangent_stress_points = F_fric[ind_u:ind_u+Dim] - Normal_stress_points*Normale
                    V_tau = V[ind_u:ind_u+Dim] - V_normal*Normale
                    Gliss_dir = Tangent_stress_points+C_tau*V_tau
                    Active_fric = (np.linalg.norm(Gliss_dir,2)+mu*Normal_stress_points>0)
                    if Active_fric:  # On impose une force tangentielle
#                        print(f"Glisse : {ind_u:5d}")
                        ind_slip_active_points = np.append(ind_slip_active_points,ind_u) 
                        G[ind_u:ind_u+Dim] = mu*Normal_stress_points*Gliss_dir/max(np.linalg.norm(Gliss_dir,2),eps)
                        nb_slip_active += 1
                    else: # On impose une vitesse tangentielle nulle (formule approchée)
#                        print(f'Adhere : {ind_u:5d}')
                        ind_stick_active_points = np.append(ind_stick_active_points,ind_u)
                        # En vitesse au milieux du pas de temps
                        delta_impose_deplac = Matr_tangent@(U_old[ind_u:ind_u+Dim]-U[ind_u:ind_u+Dim])
                        # En déplacement (deplacement tangentiel nul à l'étape suivante):
#                        delta_impose_deplac = -Tangent@U[ind_u:ind_u+2] 
                        F[ind_u:ind_u+Dim] += delta_impose_deplac 
##===                    
                        K_tangent.add([range(ind_u,ind_u+Dim)],[range(ind_u,ind_u+Dim)],Matr_tangent)
                        nb_stick_active += 1
#=== End Friction   
#=====
# Resolution (contact frottant et plasticité)
#=====
#        print('second membre : '+str(bb))
        dU = linsolve(K_tangent, F-G)
        U = U + dU
# ---
        md.set_variable("u",U)
# ---
# Calcul de la force de friction pour l'iteration suivante : on pourrait ne pas mettre -G, pour ne pas remettre +G au debut de la boucle...
        F_fric = F1 - K1.mult(dU)
#=====
# Updated of the velocity
#   Gonzalez:
        V = 2.*(U-U_old)/dt-V_old
#===== 
# Test d'arret sur le contact :
#=====   
        if (nb_active>0) or (contact_Newton_boucle>0):
           print('Active set : '+str(ind_active_points))#+'   '+str(bb)+'   '+str(L_fric))
#           print('Stick      : '+str(ind_stick_active_points))
#           print('Slip       : '+str(ind_slip_active_points))
#           print(f"Energy: { energy:.7f}\n")
           L2error = gf.compute(mf, dU, 'L2 norm', mim)
           L2error = min(L2error,gf.compute(mf, dU_old+ dU, 'L2 norm', mim)) # pour eviter les cycles sur 2 itérations
           print('-----> Contact iteration: '+str(contact_Newton_boucle)+', L2 error: '+str(L2error))
           if (L2error<eps_Newton): break
#=====
        if  (nb_active ==0) and (contact_Newton_boucle==0):
            print('-----> Contact iteration: '+str(contact_Newton_boucle)+', No contact')
            L2error = gf.compute(mf, dU, 'L2 norm', mim)
            if (contact_Newton_boucle>1):
               L2error = min(L2error,gf.compute(mf, dU_old+ dU, 'L2 norm', mim)) # pour eviter les cycles sur 2 itérations
            print('-----> L2 error: '+str(L2error))            
            if (L2error<eps_Newton): break
        dU_old = dU
#===== 
# Fin de la boucle de Newton
#===== 
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
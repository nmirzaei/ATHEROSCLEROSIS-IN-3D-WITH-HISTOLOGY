from fenics import *
import numpy as np
import os
######################################################################

#This function calculates the concentration of chemicals in the deformed region
def NC5_NonDim(kk,u,vtkfile2, vtkfile3, vtkfile4, vtkfile5,vtkfile6,volumefile1,volumefile2) :

        #Changing the directory to the home directory
        os.chdir("HPC home directory")
        ######################################################################
        kkk=kk*100
        #Uploading mesh
        Domain=Mesh('Mesh.xml')
        Bulk = MeshFunction('size_t' , Domain , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
        Boundary = MeshFunction('size_t', Domain , 'Mesh_facet_region.xml')
        ######################################################################

        #Changing the directory to the lustre directory
        os.chdir("HPC lustre directory")
        ######################################################################

        #Creating function spaces
        R_1 = VectorFunctionSpace(Domain, "R", 0)
        V_1 = VectorFunctionSpace(Domain,'P',1)
        W_1 = FunctionSpace(Domain,'P',1)
        ######################################################################

        #setting the extrapolation True for u just in case
        u.set_allow_extrapolation(True)
        ######################################################################

        #Defining measures
        ds = Measure('ds', subdomain_data=Boundary)
        dx = Measure('dx', subdomain_data= Bulk)
        ######################################################################

        #Finding Center of mass for the reference domain
        position_1 = Function(V_1)
        position_1.assign(Expression(["x[0]", "x[1]","x[2]"], element=V_1.ufl_element()))
        c_1 = TestFunction(R_1)
        volume_1 = assemble(Constant(1.0)*dx(domain=Domain))
        centroid_1 = assemble(dot(c_1, position_1)*dx)
        f_1 = centroid_1 / volume_1
        f_np_1 = f_1.get_local() # numpy array
        ######################################################################


        #We define the static values before the mesh movement if we want them to move with the mesh movements
        l = Expression('L_0*exp(-s*(x[0])*(x[0])-s*(x[1]+1.3)*(x[1]+1.3))',degree=0, s=0.05, L_0 = 2 )
        LL = project(l,W_1)
        #####################################################


        #Defining material points for tracking the lesion
        aa = Expression('x[0]',degree=0)
        bb = Expression('x[1]',degree=0)

        aaa = project(aa,W_1)
        bbb = project(bb,W_1)
        ######################################################################


        #Moving the mesh and recentering it
        ALE.move(Domain,u)
        new_V_2 = VectorFunctionSpace(Domain, 'P', 1)
        new_R_1 = VectorFunctionSpace(Domain, "R", 0)
        new_position_1 = Function(new_V_2)
        new_position_1.assign(Expression(["x[0]", "x[1]","x[2]"], element=new_V_2.ufl_element()))
        new_c_1 = TestFunction(new_R_1)
        new_volume_1 = assemble(Constant(1.0)*dx(domain=Domain))
        new_centroid_1 = assemble(dot(new_c_1, new_position_1)*dx(domain=Domain))
        new_f_1 = new_centroid_1 / new_volume_1
        new_f_np_1 = new_f_1.get_local() # numpy array
        dev_1 = Constant(f_np_1 - new_f_np_1)
        #recentering the intima and the artery
        ALE.move(Domain , dev_1)
        ######################################################################


        #This stuff here is to extract a subdomain from the original domain and mark everything according to the original Domain
        Sub_Mesh=SubMesh(Domain,Bulk,7)
        surface_marker = MeshFunction("size_t", Sub_Mesh, Sub_Mesh.topology().dim() - 1, 0)
        ncells = MeshFunction("size_t", Sub_Mesh, Sub_Mesh.topology().dim())
        volume_marker = MeshFunction("size_t", Sub_Mesh, Sub_Mesh.topology().dim())
        vmap = Sub_Mesh.data().array('parent_vertex_indices', 0)
        cmap = Sub_Mesh.data().array('parent_cell_indices', Sub_Mesh.topology().dim())

        n = 0
        for c in cells(Sub_Mesh):
             parent_cell = Cell(Domain, cmap[c.index()])
             volume_marker.array()[c.index()] = Bulk.array()[parent_cell.index()]
             for f in facets(parent_cell):
                  for g in facets(c):
                      g_vertices = vmap[g.entities(0)]
                      if set(f.entities(0)) == set(g_vertices):
                           surface_marker.array()[g.index()] = Boundary.array()[f.index()]
                  n=n+1
        ######################################################################

        #Defining measures for the deformed region
        ds = Measure('ds', subdomain_data=surface_marker)
        dx = Measure('dx', subdomain_data= volume_marker)
        ######################################################################



        #parameters
        mu = 68+kkk/30
        nu_1 = 10.68
        kappa_1 = 2.17e3
        kappa_2 = 200
        nu_2 = 0.02068
        nu_3 = 1.068e3
        nu_4 = 500
        beta_2 = 1.2e5-275*kkk
        beta_3 = 1-kkk/400
        beta_4 = 100
        beta_5 = 10
        gamma_min = 0.03
        gamma_max = 12
        C_crit = 0.3 # - kk/50#This determines the lack of oxygen
        m = 4
        ######################################################################



        #Defining mixed elements for a system of PDEs
        P2 = FiniteElement('P', tetrahedron,2)
        P1 = FiniteElement('P', tetrahedron,1)
        element = MixedElement([P2, P2, P2, P2])
        ######################################################################

        #Defining a FunctionSpace on the deformed Sub_mesh and a function that corresponds to the severity of Necrosis
        SS = FunctionSpace(Sub_Mesh,'P',1)
        Pi = Function(SS)
        ######################################################################


        #refer to the line 40 for explainantion
        a = interpolate(aaa,SS)
        b = interpolate(bbb,SS)
        ######################################

        #Defining Mixed Function Space
        Mixed_Space = FunctionSpace(Sub_Mesh, element)
        ######################################################################

        #Defining a Function on the mixed Space and test functions
        U = Function(Mixed_Space)
        v1, v2, v3, v4 =TestFunctions(Mixed_Space)

        #Spliting the mixed function into the corresponding variables
        M, N, C, Q= split(U)
        ######################################################################


        #The static values that we defined before mesh movement and layer separation need to be interpolated in the moved separated layer
        L=interpolate(LL,SS)
        ######################################################################


        #Dirichlet BC for oxygen
        bcs = [DirichletBC(Mixed_Space.sub(0), 1.15+kk/30000, surface_marker, 3),DirichletBC(Mixed_Space.sub(2), 1, surface_marker, 3)]
        ######################################################################


        #Defining our functional for our system of PDEs
        F= -nu_1*dot(grad(M),grad(v1))*dx-kappa_2*M*v1*ds(4)+mu*M*dot(grad(L+Q),grad(v1))*dx-(gamma_min+(gamma_max-gamma_min)*(C_crit**m/(C_crit**m+C**m)))*M*v1*dx\
           -nu_2*dot(grad(N),grad(v2))*dx+(gamma_min+(gamma_max-gamma_min)*(C_crit**m/(C_crit**m+C**m)))*M*v2*dx-N*v2*dx\
           -nu_3*dot(grad(C),grad(v3))*dx-beta_2*C*v3*dx-beta_3*C*M*v3*dx\
           -nu_4*dot(grad(Q),grad(v4))*dx-beta_4*Q*v4*dx+beta_5*L*M*v4*dx\
        ######################################################################

        #Solve for the mixed function
        solve(F == 0, U,bcs, solver_parameters={"newton_solver": {"relative_tolerance": 5e-10,
           "absolute_tolerance": 5e-10,"maximum_iterations": 500}})
        ######################################################################

        #Extract the variavles from the mixed solution
        M_, N_, C_, Q_= U.split()
        ######################################################################

        #Renaming for the purpose of animation
        L.rename('L','L')
        M_.rename('M','M')
        N_.rename('N','N')
        C_.rename('C','C')
        Q_.rename('Q','Q')
        Bulk.rename('Bulk','Bulk')
        ######################################################################


        #Saving the frames
        vtkfile2<<(L,kk)
        vtkfile3<<(M_,kk)
        vtkfile4<<(N_,kk)
        vtkfile5<<(C_,kk)
        vtkfile6<<(Q_,kk)
        volumefile1<< (volume_marker,kk)
        volumefile2<< (Bulk,kk)
        ######################################################################

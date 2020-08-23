from fenics import *
import numpy as np
import os
##################################################

#Thisfunction calculates the concentration of PDGFs
def pdgf(counter,kk,u,vtkfile7,vtkfile8,cc) :

        #Changing the directory to the lustre directory
        os.chdir("HPC home directory")
        ######################################################################
        #Uploading two meshes one for deformed one for reference region
        Domain =Mesh('Mesh.xml')
        Bulk = MeshFunction('size_t' , Domain , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
        Boundary = MeshFunction('size_t', Domain , 'Mesh_facet_region.xml')
        Domain1 = Mesh('Mesh.xml')
        Bulk1 = MeshFunction('size_t' , Domain1 , 'Mesh_physical_region.xml' )  #saves the interior info of the mesh
        Boundary1 = MeshFunction('size_t', Domain1 , 'Mesh_facet_region.xml')
        #######################################################################

        #Changing the directory to the lustre directory
        os.chdir("HPC lustre directory")
        ######################################################################


        #Defining Function Spaces for deformed and non-deformed regions
        R_1 = VectorFunctionSpace(Domain, "R", 0)
        V_1 = VectorFunctionSpace(Domain,'P',1)
        W_1 = FunctionSpace(Domain,'P',1)
        W_2 = FunctionSpace(Domain1,'P',1)
        #######################################################################

        #Allowing extrapolation for u
        u.set_allow_extrapolation(True)
        #######################################################################

        #Defining dof to vertex and vertex to dof maps for deformed and reference regions
        dof2v_W_1 = np.array(dof_to_vertex_map(W_1), dtype=int)
        dof2v_W_2 = np.array(dof_to_vertex_map(W_2), dtype=int)
        v2dof_W_1 = np.array(vertex_to_dof_map(W_1), dtype=int)
        v2dof_W_2 = np.array(vertex_to_dof_map(W_2), dtype=int)
        #######################################################################

        #Defining Measures
        ds = Measure('ds', subdomain_data=Boundary)
        dx = Measure('dx', subdomain_data= Bulk)
        #######################################################################


        #Finding Center of mass for the reference domain
        position_1 = Function(V_1)
        position_1.assign(Expression(["x[0]", "x[1]"], element=V_1.ufl_element()))
        c_1 = TestFunction(R_1)
        volume_1 = assemble(Constant(1.0)*dx(domain=Domain))
        centroid_1 = assemble(dot(c_1, position_1)*dx)
        f_1 = centroid_1 / volume_1
        f_np_1 = f_1.get_local() # numpy array
        #######################################################################

        #Defining material points for tracking the lesion
        aa1 = Expression('sqrt(pow(x[0]-c,2))',degree=0,c=0.2)
        bb1 = Expression('sqrt(pow(x[1]-c,2))',degree=0,c=-0.55)

        aaa1 = project(aa1,W_1)
        bbb1 = project(bb1,W_1)
        #######################################################################


        #Moving the mesh and recentering
        ALE.move(Domain,u)
        new_V_2 = VectorFunctionSpace(Domain, 'P', 1)
        new_R_1 = VectorFunctionSpace(Domain, "R", 0)
        new_position_1 = Function(new_V_2)
        new_position_1.assign(Expression(["x[0]", "x[1]"], element=new_V_2.ufl_element()))
        new_c_1 = TestFunction(new_R_1)
        new_volume_1 = assemble(Constant(1.0)*dx(domain=Domain))
        new_centroid_1 = assemble(dot(new_c_1, new_position_1)*dx(domain=Domain))
        new_f_1 = new_centroid_1 / new_volume_1
        new_f_np_1 = new_f_1.get_local() # numpy array
        dev_1 = Constant(f_np_1 - new_f_np_1)
        #recentering the intima and the artery
        ALE.move(Domain , dev_1)
        #######################################################################

        #PDGF Parameters
        #kappa = Expression('20/(a+(tt/3))',degree=0,a=1,tt=kk)
        kappa = cc
        nu = 100
        beta = 20
        #######################################################################

        #Function space on the deformed domain
        SS= FunctionSpace(Domain,'P',1)
        #######################################################################


        #refer to the line 51 for explainantion
        a1 = interpolate(aaa1,SS)
        b1 = interpolate(bbb1,SS)
        #######################################################################

        #Defining Function P for PDGF on the deformed domain and a test function
        P = Function(SS)
        v = TestFunction(SS)
        #######################################################################
        P0 = Expression('0.1*exp(-s*(x[0]-c1)*(x[0]-c1)-s*(x[1]-c2)*(x[1]-c2))',degree=0,s=0.3,c1 = -0.42,c2=-0.64)
        DBC = DirichletBC(SS, P0, Boundary , 1)

        #Defining the functional
        F= -nu*dot(grad(P),grad(v))*dx-beta*P*v*dx#+kappa*conditional(lt(a1,Constant(1)),1,0)*conditional(lt(b1,Constant(0.25)),1,0)*(1-P)*v*ds(1)
        #######################################################################

        #Solve for P in the deformed domain
        solve(F== 0, P ,DBC, solver_parameters={"newton_solver": {"relative_tolerance": 5e-10,
                    "absolute_tolerance": 5e-10,"maximum_iterations": 500}})
        #######################################################################

        #Define the function PP on the reference domain and put its values in an auxillary vector Array
        PP = Function(W_2)
        Array = PP.vector().get_local()
        #######################################################################

        #Replace the values of the auxillary vector with the values of P from the deformed domain
        #We do so by setting the values of PP on the dofs from the reference domain to be the values of P for the same dofs in the deformed domain
        Array[v2dof_W_2[dof2v_W_1]] = P.vector().get_local()
        #######################################################################

        #Let PP be the auxillary Array
        PP.vector()[:]=Array
        #######################################################################

        #Renaming for animation plots
        P.rename('P','P')
        PP.rename('PP','PP')
        #######################################################################

        #We plot only every 10 steps to save memory
        #Also we track PDGF in the deformed and also its representation in the undeformed region
        if counter%10==0:
            vtkfile7<<(P,kk)
            vtkfile8<<(PP,kk)
        #######################################################################

        #We use the values of PDGF in the undeformed region to create the growth tensor because it applies to the reference domain at each step
        return PP, kappa
        #######################################################################

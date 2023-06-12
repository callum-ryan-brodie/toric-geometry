
import copy
import numpy as np


# Some useful fans
p1xp1_fan = Fan2d(rays=[[0,-1],[1,0],[0,1],[-1,0]])
p2_fan = Fan2d(rays=[[-1,0],[0,-1],[1,1]])
dp3_fan = Fan2d(rays=[[1,0],[1,1],[0,1],[-1,0],[-1,-1],[0,-1]])


#
# This computes the hodge number h^{k,k} of a toric variety associated to a fan f
# input: f,K
# output: h^{k,k}

def hodge_toric(f,k):
    n = f.dim()
    h = 0      
    for i in range(k,n+1):
        h = h + ((-1)**(i-k)) * binomial(i,k) * len(f(n-i))   
        
    return h


#
# Checks if a CY hypersurface (from the anti-canonical bundle) in an ambient toric variety tv is smooth
# input: tv
# output: true/false

def antican_hypsur_misses_ambsings(tv):

    f = tv.fan()
    
    CR = tv.cohomology_ring()
    K = CR(tv.K())
    
    dim = tv.dimension()

    for i in range(0,dim+1):
        for c in f(i):
            if not c.is_smooth():
                if K*CR(c) != CR(0):
                    return False
    
    return True


#
# Checks that a hypersurface defined by a divisor D in an ambient toric variety tv doesn't intersect singularities of tv
# input: tv, D
# output: boolean

def hypsur_misses_ambsings(tv,D):

    f = tv.fan()
    
    CR = tv.cohomology_ring()

    for dim in range(0,tv.dimension()+1):

        for i in range(0,dim+1):
            for c in f(i):
                if not c.is_smooth():
                    if D*CR(c) != CR(0):
                        return False
    
    return True


#
# Returns the Polyhedron P corresponding to a toric variety TV
# input: TV
# output: P

def poly_of(TV):

    raylist = []
    for cone in TV.fan().cones(1):
        raylist.append(cone.rays()[0].list())

    return Polyhedron(vertices=raylist)


#
# Returns the LatticePolytope LP corresponding to a toric variety TV
# input: TV
# output: LP

def latpoly_of(TV):

    raylist = []
    for cone in TV.fan().cones(1):
        raylist.append(cone.rays()[0].list())

    return LatticePolytope(raylist)


#
# Returns the number N of non-polynomial (complex-structure) deformations of a CY 
# defined as the anti-canonical hypersurface in a lattice polytope P of dimension >=4 
# input: P
# output: N

def nonpolys(P):

    PD = P.polar()
    n = P.dim()        
    
    nonpolynum = 0 
        
    for faceP in P.faces(1):
    
        for facePD in PD.faces(n-2):
    
            intersections = [vector(dualpt)*vector(pt) for dualpt in list(facePD.as_polyhedron().vertices_matrix().transpose()) for pt in list(faceP.as_polyhedron().vertices_matrix().transpose())]
    
            #print intersections 
            if all(int == -1 for int in intersections):
                
                facePintpoints = LatticePolytope(faceP.as_polyhedron().vertices_matrix().transpose()).interior_points()
                facePDintpoints = LatticePolytope(facePD.as_polyhedron().vertices_matrix().transpose()).interior_points()
        
                nonpolynum = nonpolynum + len(facePintpoints)*len(facePDintpoints)
                
                break

    return nonpolynum


#
# Returns the Hodge number h11 of the CY hypersurface in a toric variety with a
# lattice polytope P of dimension >=4 
# input: P
# output: h11

def myh11(P):
    
    PD = P.polar()
    n = P.dim()

    h11 = PD.npoints() - n - 1
    
    for face in PD.faces(n-1):
       h11 = h11 - len(face.interior_points())
       
    h11 = h11 + nonpolys(P.polyhedron())
    
    return h11


#
# Returns the Hodge number h21 of the CY hypersurface in a toric variety with a
# lattice polytope P of dimension >=4 
# input: P
# output: h21

def myh21(P):
    
    PD = P.polar()
    n = P.dim()
    
    h21 = P.npoints() - n - 1
    
    for face in P.faces(n-1):
       h21 = h21 - len(face.interior_points())
       
    h21 = h21 + nonpolys(PD.polyhedron())
    
    return h21


#
# Given a fan and a vector, returns the fan of the subvariety corresponding to setting the coordinate of that vector equal to zero
# input: fan, vec1
# output: new_fan

def restr_to_subvar_onevec(fan,vec1):

    cones_to_project= []
    rays_involved = []
    
    for cone in fan(fan.dim()): #check for which cones contain vec1 (we will project these)
        new_cone = []
        mat = cone.rays().matrix()
        vec1_cont = False
        for i in range(0,mat.nrows()):
            if mat[i].list() == vec1:
                vec1_cont = True
        if vec1_cont: #if vector is in this cone
            for row in mat: #record this cone for our list
                if (not row.list() == vec1): #but we don't need the rays that equal vec1
                        new_cone.append(row.list())
                        if not row.list() in rays_involved:
                            rays_involved.append(row.list())
            cones_to_project.append(new_cone) #add the cone to the list
    
    vec_mat = matrix([vec1])
    ker_bas = vec_mat.right_kernel().basis()
    basis_for_orth_space = [ker_bas[i].list() for i in range (0,len(ker_bas))]
    
    final_cones = []
    
    print "Under the projection the rays involved projected as follows:"
    for i in range(0,len(rays_involved)):
        print str(rays_involved[i]) + " went to " + str([np.dot(rays_involved[i],basis_for_orth_space[j]) for j in range (0,len(basis_for_orth_space))])
    
    for cone in cones_to_project:
        new_cone = []
        for ray in cone: #project onto basis of orthogonal space
            ray = [np.dot(ray,basis_for_orth_space[i]) for i in range (0,len(basis_for_orth_space))]
            new_cone.append(ray)
        final_cones.append(new_cone)

    if final_cones == []:
	print "Error in restr_to_subvar_onevec: projected to the empty set."
	return
    
    return Fan([Cone(i) for i in final_cones],)#lattice=ToricLattice(len(final_cones[0]), "N"))


#
# Given a fan and two vectors, returns the fan of the subvariety corresponding to setting the coordinates of those vector equal to zero
# input: fan, vec1, vec2
# output: new_fan

def restr_to_subvar_twovec(fan,vec1,vec2):

    cones_to_project= []
    rays_involved = []
    
    for cone in fan(fan.dim()): #check for which cones contain vec1 and vec2 (we will project these)
        new_cone = []
        mat = cone.rays().matrix()
        vec1_cont = False
        vec2_cont = False
        for i in range(0,mat.nrows()):
            if mat[i].list() == vec1:
                vec1_cont = True
            if mat[i].list() == vec2:
                vec2_cont = True
        if vec1_cont and vec2_cont: #if both vectors are in this cone
            for row in mat: #record this cone for our list
                if (not row.list() == vec1) and (not row.list() == vec2): #but we don't need the rays that equals vec1 or vec2
                        new_cone.append(row.list())
                        if not row.list() in rays_involved:
                            rays_involved.append(row.list())
            cones_to_project.append(new_cone) #add the cone to the list
    
    vec_mat = matrix([vec1,vec2])
    ker_bas = vec_mat.right_kernel().basis()
    basis_for_orth_space = [ker_bas[i].list() for i in range (0,len(ker_bas))]
    
    final_cones = []
    
    print "Under the projection the rays involved projected as follows:"
    for i in range(0,len(rays_involved)):
        print str(rays_involved[i]) + " went to " + str([np.dot(rays_involved[i],basis_for_orth_space[j]) for j in range (0,len(basis_for_orth_space))])
    
    for cone in cones_to_project:
        new_cone = []
        for ray in cone: #project onto basis of orthogonal space
            ray = [np.dot(ray,basis_for_orth_space[i]) for i in range (0,len(basis_for_orth_space))]
            new_cone.append(ray)
        final_cones.append(new_cone)
    
    return Fan([Cone(i) for i in final_cones])


#
# Checks whether all the dot products between vectors in two lists are non-negative
# (For use in the fan_is_refinement function)
def all_dot_nonneg(veclist1,veclist2):
    for vec1 in veclist1:
        for vec2 in veclist2:
            if np.dot(vec1,vec2) < 0:
                return False
    return True


#
# Checks whether one fan f1 is a refinement of another f2
# input: f1, f2
# output: True or False

def fan_is_refinement(f1,f2):
    for cone1 in f1(f1.dim()):
        cone1_in_cone2 = False
        for cone2 in f2(f2.dim()):
            #put cone generators into a more sensible form
            cone1_gens = []
            cone2_dual_gens = []
            for vec1 in cone1.rays():
                cone1_gens.append(vec1.list())
            for vec2 in cone2.dual().rays():
                cone2_dual_gens.append(vec2.list())
            if all_dot_nonneg(cone1_gens,cone2_dual_gens):
                cone1_in_cone2 = True
                break
        if cone1_in_cone2 == False:
	    #print cone1.rays()
            return False
    return True


#
# Takes a fan and takes the product with a P1
# input: fan of B2
# output: fan of B2xP1

def product_p1_with(old_fan):
    
    dim = old_fan.dim()

    # take the (n)d cones and embed them in (n+1)d (also makes into list of list of vectors)
    old_cones_embedded = []
    for cone in old_fan.cones(dim):
        embedded_cone = []
        for row in cone.rays().matrix():
            embedded_cone.append(row.list()+[0])
        old_cones_embedded.append(embedded_cone)

    new_cones = []
    for embedded_cone in old_cones_embedded:
    	new_cones.append(embedded_cone+[[0]*dim+[-1]])
    	new_cones.append(embedded_cone+[[0]*dim+[ 1]])
    return Fan([Cone(i) for i in new_cones])


#
# Takes an n-dimensional fan fan1 and an m-dimensional fan fan2 and takes the product to give new_fan
# input: old_fan
# output: new_fan

def product_fans(fan1,fan2):
    
    dim1 = fan1.dim()
    dim2 = fan2.dim()

    # take the (dim1)d cones in fan1 and embed them in (dim1+dim2)d (also makes into list of list of vectors)
    old_fan1_cones_embedded = []
    for cone in fan1.cones(dim1):
        embedded_cone = []
        for row in cone.rays().matrix():
            embedded_cone.append(row.list()+[0]*dim2)
        old_fan1_cones_embedded.append(embedded_cone)

    # take the (dim1)d cones in fan2 and embed them in (dim1+dim2)d (also makes into list of list of vectors)
    old_fan2_cones_embedded = []
    for cone in fan2.cones(dim2):
        embedded_cone = []
        for row in cone.rays().matrix():
            embedded_cone.append([0]*dim1+row.list())
        old_fan2_cones_embedded.append(embedded_cone)

    new_cones = []
    for embedded_cone1 in old_fan1_cones_embedded:
	for embedded_cone2 in old_fan2_cones_embedded:
    	    new_cones.append(embedded_cone1+embedded_cone2)
    return Fan([Cone(i) for i in new_cones])


#
# Fibres a P123 over a base B to give the ambient space A of an elliptically fibered CY
# input: fan of B
# output: fan of A

def fibre_p123_over(basefan):

    base_dim = basefan.dim()
    zeroes = [0] * base_dim

    fibre_conelist = [[zeroes+[2,3],zeroes+[-1,0]],[zeroes+[2,3],zeroes+[0,-1]],[zeroes+[-1,0],zeroes+[0,-1]]]
    
    base_conelist = []    
    for cone in basefan.cones(base_dim):
        base_cone = []
        for row in cone.rays().matrix():
            base_cone.append(row.list()+[2,3])
        base_conelist.append(base_cone)
    
    totalspace_conelist = []
    for base_cone in base_conelist:
        for fibre_cone in fibre_conelist:
            cone = base_cone + fibre_cone
            totalspace_conelist.append(cone)
	
    return Fan([Cone(i) for i in totalspace_conelist])


#
# Fibres a P112 over a base B to give the ambient space A of an elliptically fibered Calabi-Yau
# input: fan of B
# output: fan of A

def fibre_p112_over(basefan):

    base_dim = basefan.dim()
    zeroes = [0] * base_dim

    fibre_conelist = [[zeroes+[1,2],zeroes+[-1,0]],[zeroes+[1,2],zeroes+[0,-1]],[zeroes+[-1,0],zeroes+[0,-1]]]
    
    base_conelist = []    
    for cone in basefan.cones(base_dim):
        base_cone = []
        for row in cone.rays().matrix():
            base_cone.append(row.list()+[1,2])
        base_conelist.append(base_cone)
    
    fivefold_conelist = []
    for base_cone in base_conelist:
        for fibre_cone in fibre_conelist:
            cone = base_cone + fibre_cone
            fivefold_conelist.append(cone)
	
    return Fan([Cone(i) for i in fivefold_conelist])


#
# Gives the lattice polytope of a P123 fibered over a base B which forms the ambient space A of an elliptically fibered Calabi-Yau
# input: fan of B
# output: lattice polytope of A

def latpoly_of_fibered_p123_over(fan):

    base_dim = fan.dim()

    ray_list = []
    for ray in fan(1):
        ray_list.append(ray.rays().matrix().list()+[2,3])
    ray_list += [[0]*base_dim+[-1,0],[0]*base_dim+[0,-1],[0]*base_dim+[2,3]]

    return(LatticePolytope(ray_list))


#
# Gives the lattice polytope of a P112 fibered over a base B which forms the ambient space A of an elliptically fibered Calabi-Yau
# input: fan of B
# output: lattice polytope of A


def latpoly_of_fibered_p112_over(fan):

    base_dim = fan.dim()

    ray_list = []
    for ray in fan(1):
        ray_list.append(ray.rays().matrix().list()+[1,2])
    ray_list += [[0]*base_dim+[-1,0],[0]*base_dim+[0,-1],[0]*base_dim+[1,2]]

    return(LatticePolytope(ray_list))


#
# Generates a tower of rays in 3d over initialvec, from z = n_min to z = n_max
# input: initialvec, n_min, n_max
# output: tower of rays

def gen_res_rays(initialvec,n_min,n_max):

    resolutionray_list = []
    for i in range (n_min,n_max+1):
        resolutionray = copy.copy(initialvec)
        resolutionray[2] = initialvec[2] + i
        resolutionray_list.append(resolutionray)
    return resolutionray_list


#
# Takes a ray and returns a list of rays, each with z added to the original nth element where z runs from z = m_min to z = m_max
# input: initialvec, n, m_min, m_max
# output: tower of rays

def gen_ray_tower(initialvec,n,m_min,m_max):

    resolutionray_list = []
    for i in range (m_min,m_max+1):
        resolutionray = copy.copy(initialvec)
        resolutionray[n] = initialvec[n] + i
        resolutionray_list.append(resolutionray)
    return resolutionray_list


#
# For a CY threefold which is a hypersurface in a toric variety with polytope P, computes the Hodge numbers 
# h0i, i = 0,1,2 of a divisor which descends from a toric divisor of the ambient space corresponding to a lattice point pt
# input: P, pt
# output: (h00,h10,h20,h30) 

def h0i(P,pt):
    
    # first we need to know in which face pt is sitting and if it is in the polytope at all

    if vector(pt) not in P.integral_points():
        print 'not contained in polytope'  

        return 0

    else:

        FF = FaceFan(P)
        d = P.dim()
    
        for i in range(d+1):
            for cone in FF.cones(i):
                if cone.relative_interior_contains(pt):
                    Pfaceverts = cone.rays()
                    k = i-1
        #print Pfaceverts
        
        
    # compute which are the vertices of the dual face on the M-lattice polytope    
        
        PD=P.polar()
        for face in PD.faces(d-k-1):

            intersections = [vector(dualpt)*vector(pt) for dualpt in list(face.as_polyhedron().vertices_matrix().transpose())]
            # print intersections 
            if all(int == -1 for int in intersections):
                Pdualfacevert = list(face.as_polyhedron().vertices_matrix().transpose())
                #print Pdualfacevert
                break   
                
    # now we need to know how many interior points this face has, we do this case by case
    
        dualfacep = Polyhedron(vertices=Pdualfacevert)     
        PCdF = PointConfiguration(dualfacep.integral_points())
        numintpts = len(PCdF.face_interior(dim=d-k-1))
        
        if d-k == d:
            h00=1
            h10=0
            h20=0
	    h30=numintpts 
            return (h00,h10,h20,h30)  
        if d-k == d-1:
            h00=1
	    h10=0
            h20=numintpts
            h30=0 
            return (h00,h10,h20,h30)  
	if d-k == d-2:
            h00=1
	    h10=numintpts
            h20=0
            h30=0 
            return (h00,h10,h20,h30) 
        if d-k == d-3:
            h00= numintpts + 1
            h10=0
            h20=0  
	    h30=0 
            return (h00,h10,h20,h30)     
        else:
            #print 'higher codimension than 3'
            return (0,0,0,0)


#
# Checks whether the assumption in Bertinis theorem holds for O(div) for some divisor div
# (Bertini's theorem: if there's nowhere that all sections vanish, then generic section in compact complex manifold is smooth)
# input: fan of X, and divisor div
# output: Boolean result

def bertinis_theorem_holds_for_div(fan,div):
    
    tv = ToricVariety(fan)
    
    num_coords = len(fan(1))
    
    zero_coord_set_list = []
    for top_dim_cone in fan(fan.dim()): #collect all the largest sets of coordinates that can vanish together
    
        ray_list = [ray.list() for ray in top_dim_cone.rays()]
        
        fan_coord_list = []
        i=0
        for c in fan(1):
            if c.rays().matrix().list() in ray_list:
                fan_coord_list.append(i)
            i += 1
            
        zero_coord_set_list = zero_coord_set_list + [fan_coord_list]
    
    monoms = div.sections_monomials()
    monoms_degrees = [monom.degrees() for monom in monoms] #now collect all monomials

    for monom_degrees in monoms_degrees:
        if len(monom_degrees) != num_coords:
            print "Error: Degrees aren't being handled correctly"
            return None
    
    for zero_coord_set in zero_coord_set_list: #go through all the sets of vanishing coordinates
        monoms_all_zero = True
        for monom_degrees in monoms_degrees: #for each go through all sets of monomials
            monom_zero = False
            for i in zero_coord_set:
                if monom_degrees[i] != 0: #check if a coordinate being set to zero appears in this monomomial
                    monom_zero = True
            if monom_zero == False:
                monoms_all_zero = False
                break
        if monoms_all_zero == True:
            #print "Found somewhere where all monomials vanish"
            return False
            break
            
    return True


#
# Checks whether the assumption in Bertinis theorem holds for the anti-canonical bundle of a space X
# (Bertini's theorem: if there's nowhere that all sections vanish, then generic section in compact complex manifold is smooth)
# input: fan of X
# output: Boolean result

def bertinis_theorem_holds_for_antican(fan):
    tv = ToricVariety(fan)
    
    num_coords = len(fan(1))
    
    zero_coord_set_list = []
    for top_dim_cone in fan(fan.dim()): #collect all the largest sets of coordinates that can vanish together
    
        ray_list = [ray.list() for ray in top_dim_cone.rays()]
        
        fan_coord_list = []
        i=0
        for c in fan(1):
            if c.rays().matrix().list() in ray_list:
                fan_coord_list.append(i)
            i += 1
            
        zero_coord_set_list = zero_coord_set_list + [fan_coord_list]
    
    tv_K = -tv.K()
    tv_K_monoms = tv_K.sections_monomials()
    tv_K_monoms_degrees = [monom.degrees() for monom in tv_K_monoms] #now collect all monomials

    for monom_degrees in tv_K_monoms_degrees:
        if len(monom_degrees) != num_coords:
            print "Error: Degrees aren't being handled correctly"
            return None
    
    for zero_coord_set in zero_coord_set_list: #go through all the sets of vanishing coordinates
        monoms_all_zero = True
        for monom_degrees in tv_K_monoms_degrees: #for each go through all sets of monomials
            monom_zero = False
            for i in zero_coord_set:
                if monom_degrees[i] != 0: #check if a coordinate being set to zero appears in this monomomial
                    monom_zero = True
            if monom_zero == False:
                monoms_all_zero = False
                break
        if monoms_all_zero == True:
            #print "Found somewhere where all monomials vanish"
            return False
            break
            
    return True


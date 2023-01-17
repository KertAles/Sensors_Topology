import gudhi
from utils import read_data, distance, numberOfComponents, powerset, collapse
import numpy as np
from gudhi import AlphaComplex
from itertools import chain, combinations
from mogutda import SimplicialComplex
import numpy as np

from export_triangulation_to_ply import export_ply


import math
import miniball
import matplotlib.pyplot as plt

# Get Čech complex
def cech2(S, epsilon) :
    # Get all edges using VR
    vertices, edges = VR_edges_only(S, epsilon*2)

    sx_dict = {}
    sx_dict[0] = vertices
    sx_dict[1] = edges
    unchecked_vert = []

    # Make adjacency list
    G = {}
    for vert in vertices :
        G[vert[0]] = []
        unchecked_vert.append(S[vert[0]])
    for edge in edges :
        G[edge[0]].append(S[edge[1]])
        G[edge[1]].append(S[edge[0]])

    # Go through all vertices
    for vertex in vertices :
        unchecked_vert.remove(S[vertex[0]])
        relevant_vert = list(set(unchecked_vert) & set(G[vertex[0]]))
        # Get all subsets of neighbouring vertices(that haven't been checked yet)
        combinations = powerset(relevant_vert)
        for j, comb in enumerate(combinations) :
            if len(comb) >= 3 :
                break
            # Check subsets with more than 1 vertex - when we add the current vertex, it makes an edge, which is boring
            if len(comb) > 1:
                # Make a set of points that include the current vertex
                comb = list(comb)
                comb.append(S[vertex[0]])
                comb = np.array(comb)
                # Get the miniball
                _, r2 = miniball.get_bounding_ball(comb)

                # miniball returns the square of the radius
                if math.sqrt(r2) < epsilon :
                    dim = len(comb) - 1
                    if dim not in sx_dict:
                        sx_dict[dim] = []

                    # Get the indices of the used points - coordinates were used for miniball calculation
                    nu_sx = []
                    for v in comb :
                        v_t = tuple(v)
                        nu_sx.append(S.index(v_t))
                    nu_sx = tuple(sorted(nu_sx))
                    if nu_sx not in sx_dict[dim]:
                        sx_dict[dim].append(nu_sx)

    return sx_dict

# Generate a VR complex with max_dim=1
def VR_edges_only(S, epsilon) :
    edges = []
    vertices = []

    for i, v1 in enumerate(S):
        vertices.append((i,))
        for j, v2 in enumerate(S[i+1:]):
            if distance(v1, v2) < epsilon:
                edges.append((i, j+i+1))

    return (vertices, edges)


# Find the optimal parameter r, using VR complex that only generates edges.
# The function returns the optimal r, so that the VR complex it generates is connected
def find_optimal_r(points, start, stop, init_step, eps=1e-6, factor=2, verbose=True) :
    step = init_step
    curr_r = -1

    bottom_r = start
    top_r = stop
    while step > eps:
        if verbose :
            print('Searching interval from ' + str(bottom_r) + ' to ' + str(top_r) + ' with step ' + str(step))
        for r in np.arange(bottom_r, top_r, step):
            rips_complex = VR_edges_only(points, r)

            components = numberOfComponents(rips_complex[0], rips_complex[1])

            if components == 1 :
                curr_r = r
                top_r = r
                break
            else:
                bottom_r = r
        step /= factor

    if curr_r > 0:
        print('Found connected VR complex at r = ' + str(curr_r))
        return curr_r, rips_complex
    else:
        print('Found no sphere on given interval.')
        return -1, None



"""
    ac = AlphaComplex(points=points)

    stree = ac.create_simplex_tree()
    filtered_limits = []
    for filtered_value in stree.get_filtration():
        filtered_limits.append(math.sqrt(filtered_value[1]))

    for R in filtered_limits :
        if R < start :
            continue
        if R > stop :
            break

        vertices = 0
        edges = 0
        faces = 0

        for filtered_value in stree.get_filtration():
            if filtered_value[1] <= R :
                dim = len(filtered_value[0])
                if dim == 1 :
                    vertices += 1
                elif dim == 2 :
                    edges += 1
                elif dim == 3 :
                    faces += 1
            else :
                break

        euler = vertices - edges + faces
        #print(str(R) + ' --- ' + str(euler))
"""

def construct_complex(cech) :
    complex = cech[2]

    for sx1 in cech[1] :
        is_included = False

        for sx2 in complex :
            if set(sx1).issubset(set(sx2)) :
                is_included = True
                break

        if not is_included :
            complex.append(sx1)

    for sx1 in cech[0]:
        is_included = False

        for sx2 in complex:
            if set(sx1).issubset(set(sx2)):
                is_included = True
                break

        if not is_included:
            complex.append(sx1)

    return complex

# Find the optimal parameter R, using the Čech complex.
def find_optimal_R(points, start, stop, init_step, eps=1e-6, factor=2, verbose=True) :
    step = init_step
    curr_R = -1

    bottom_R = start
    top_R = stop
    while step > eps :
        if verbose :
            print('Searching interval from ' + str(bottom_R) + ' to ' + str(top_R) + ' with step ' + str(step))
        for R in np.arange(bottom_R, top_R, step) :

            #cech_cx = cech(points, R)
            ac = AlphaComplex(points)
            stree = ac.create_simplex_tree(max_alpha_square=R**2)

            cech_sxes = []
            for simplex in stree.get_simplices():
                cech_sxes.append(tuple(simplex[0]))

            sphere = SimplicialComplex(simplices=cech_sxes)
            euler = sphere.euler_characteristics()

            betti0 = sphere.betti_number(0)
            betti1 = sphere.betti_number(1)
            betti2 = sphere.betti_number(2)

            """
            if 2 in cech_cx :
                cech_sxes = construct_complex(cech_cx)
                sphere = SimplicialComplex(simplices=cech_sxes)
                betti0 = sphere.betti_number(0)
                betti1 = sphere.betti_number(1)
                betti2 = sphere.betti_number(2)
                euler = sphere.euler_characteristics()
                #print(betti0)
                #print(betti1)
                #print(betti2)
            """
            if verbose :
                print(str(R) + ' --- b0 : ' + str(betti0) + ' b1: ' + str(betti1) + ' b2: ' + str(betti2) + ' Euler: ' + str(euler))
            if betti0 == 1 and betti1 == 0 and euler >= 2:
                curr_R = R
                top_R = R
                break
            else :
                bottom_R = R

        step /= factor

    if curr_R > 0:
        print('Found complex with Euler characteristic of a sphere at R = ' + str(curr_R))
        return curr_R, cech_sxes
    else :
        print('Found no sphere on given interval.')
        return -1, cech_sxes

def checkObsolescence(points, r, R, verbose=True) :
    obsolete_points = []
    vital_points = []
    points_copy = points.copy()

    for p in points :
        if p in obsolete_points :
            continue

        points_copy.remove(p)
        #print(len(points_copy))

        #cech_cx = cech(points_copy, R)
        ac = AlphaComplex(points_copy)
        stree = ac.create_simplex_tree(max_alpha_square=R ** 2)

        cech_sxes = []
        for simplex in stree.get_simplices():
            cech_sxes.append(tuple(simplex[0]))

        sphere = SimplicialComplex(simplices=cech_sxes)
        euler = sphere.euler_characteristics()

        betti0 = sphere.betti_number(0)
        betti1 = sphere.betti_number(1)
        betti2 = sphere.betti_number(2)

        rips_complex = VR_edges_only(points_copy, r)
        vr_components = numberOfComponents(rips_complex[0], rips_complex[1])

        if verbose :
            print(str(R) + ' --- b0 : ' + str(betti0) + ' b1: ' + str(betti1) + ' b2: ' + str(betti2) + ' Euler: ' + str(euler))

        if betti0 == 1 and vr_components == 1 and euler == 2:
            if verbose :
                print('Found obsolete point ' + str(p))
            obsolete_points.append(p)
        else :
            if verbose :
                print(str(p) + ' not obsolete')
            vital_points.append(p)
            points_copy.append(p)

    return vital_points, obsolete_points


if __name__ == '__main__':
    points = read_data('./data/sensors.txt')
    verbose = False

    r, vr_cx = find_optimal_r(points, 0.25, 0.7, 0.1, verbose=verbose)
    R, cech_cx = find_optimal_R(points, r/2, r*2, 0.1, eps=1e-6, verbose=verbose)
    
    #cech_sxes = cech_cx
    #sphere = SimplicialComplex(simplices=cech_sxes)
    #print(f'Optimal sphere :  Euler : {sphere.euler_characteristics()} -- b0 : {sphere.betti_number(0)} -- b1 : {sphere.betti_number(1)} -- b2 : {sphere.betti_number(2)}')

    ac = AlphaComplex(points)
    stree = ac.create_simplex_tree(max_alpha_square=R ** 2)
    cech_sxes = []
    for simplex in stree.get_simplices() :
        cech_sxes.append(tuple(simplex[0]))

    #cech_cx = construct_complex(cech(points, R))
    #betti0 = -1
    #betti1 = -1
    #cech_sxes = cech_cx
    sphere = SimplicialComplex(simplices=cech_sxes)
    print(f'Optimal sphere :  Euler : {sphere.euler_characteristics()} -- b0 : {sphere.betti_number(0)} -- b1 : {sphere.betti_number(1)} -- b2 : {sphere.betti_number(2)}')
    data = export_ply(cech_sxes, points)

    f = open("sphere.ply", "w")
    f.write(data)
    f.close()

    vital_points, obsoletePoints = checkObsolescence(points, r, R, verbose=verbose)
    print('Found vital points : ')
    print(vital_points)
    print('Found obsolete points : ')
    print(obsoletePoints)


    ac = AlphaComplex(vital_points)
    stree = ac.create_simplex_tree(max_alpha_square=R ** 2)
    cech_sxes = []
    for simplex in stree.get_simplices() :
        cech_sxes.append(tuple(simplex[0]))

    sphere = SimplicialComplex(simplices=cech_sxes)
    print(f'Trimmed sphere :  Euler : {sphere.euler_characteristics()} -- b0 : {sphere.betti_number(0)} -- b1 : {sphere.betti_number(1)} -- b2 : {sphere.betti_number(2)}')

    data = export_ply(cech_sxes, points)

    f = open("sphere_no_obsolete_points.ply", "w")
    f.write(data)
    f.close()

import gudhi
from utils import read_data, distance, numberOfComponents, powerset
import numpy as np

import math
import miniball

# Get Čech complex
def cech(S, epsilon) :
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

        for comb in combinations :
            # Check subsets with more than 1 vertex - when we add the current vertex, it makes an edge, which is boring
            if len(comb) > 1 :
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
def find_optimal_r(points, start, stop, delta=1e-3) :

    for r in np.arange(start, stop, delta) :
        rips_complex = VR_edges_only(points, r)

        components = numberOfComponents(rips_complex[0], rips_complex[1])
        if components  == 1:
            print('Found connected VR complex at r = ' + str(r))
            return r, rips_complex

    print('Found no connected VR complex on given interval.')
    return -1, None


# Find the optimal parameter R, using the Čech complex.
# The function returns the optimal parameter R, so that the Čech complex it generates
# has the same Euler characteristic as a sphere. Checking if connected isn't necessary, as the R
# is always larger than r/2 which makes a connected complex in VR( and Čech with R is a subset of VR with R/2).
def find_optimal_R(points, start, stop, delta=1e-3) :

    for R in np.arange(start, stop, delta) :
        print(R)
        cech_cx = cech(points, R)
        vertices = len(cech_cx[0])
        edges = len(cech_cx[1])
        faces = 0
        if 2 in cech_cx :
            faces = len(cech_cx[2])

        euler = vertices - edges + faces
        print(str(R) + ' --- ' + str(euler))
        if euler == 2 :
            print('Found complex with Euler characteristic of a sphere at R = ' + str(R))
            return R, cech_cx

    print('Found no sphere on given interval.')
    return -1, None



if __name__ == '__main__':
    points = read_data('./data/sensors.txt')

    r, vr_cx = find_optimal_r(points, 0.25, 0.7)
    R, cech_cx = find_optimal_R(points, r, 2.0, 0.5)
import gudhi
from utils import read_data
import numpy as np


def findComponents(V, E):
    components = []

    G = {}
    for vert in V:
        G[vert] = []
    for edge in E:
        G[edge[0]].append(edge[1])
        G[edge[1]].append(edge[0])

    V_search = V

    while len(V_search) > 0:
        vert = V_search[0]
        nu_comp = expandComponent(vert, G)
        components.append(nu_comp)

        for c_vert in nu_comp:
            V_search.remove(c_vert)

    return components

def expandComponent(vert, G):
    component = []
    component.append(vert)

    to_search = []
    curr_vert = vert

    for neighbour in G[curr_vert]:
        to_search.append(neighbour)

    while len(to_search) > 0:
        curr_vert = to_search.pop(0)
        component.append(curr_vert)

        for neighbour in G[curr_vert]:
            if not (neighbour in component) and not (neighbour in to_search):
                to_search.append(neighbour)

    return component
def find_optimal_r(points, start, stop, delta=1e-3) :
    #diff = stop - start

    #while diff < eps :
    for r in np.arange(start, stop, delta) :
        rips_complex = gudhi.RipsComplex(points=points, max_edge_length=r)
        simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
        #result_str = 'Rips complex is of dimension ' + repr(simplex_tree.dimension()) + ' - ' + \
        #             repr(simplex_tree.num_simplices()) + ' simplices - ' + \
        #             repr(simplex_tree.num_vertices()) + ' vertices.'
        #print(result_str)
        vertices = simplex_tree.num_vertices()
        edges = 0
        faces = 0
        edge_list = []
        vert_list = []
        for filtered_value in simplex_tree.get_filtration():
            dim = len(filtered_value[0])
            if dim == 1 :
                vert_list.append(filtered_value[0][0])
            elif dim == 2 :
                edge_list.append(filtered_value[0])
                edges += 1
            elif dim == 3 :
                faces += 1

        euler = vertices - edges + faces
        if euler == 2 :
            components = len(findComponents(vert_list, edge_list))
            if components  == 1:
                return r, rips_complex

    print('Found no sphere on given interval.')
    return -1, None


if __name__ == '__main__':
    points = read_data('./data/sensors.txt')

    find_optimal_r(points, 0.6, 0.7)
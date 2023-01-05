import math
from itertools import chain, combinations

def read_data(path) :
    points = []

    f = open(path, "r")
    data_string = f.read()

    for point in data_string.split('{') :
        if len(point) > 0 :
            point_split = point.split(',')
            point = (float(point_split[0]), float(point_split[1]), float(point_split[2].split('}')[0]))

            points.append(point)

    return points

def distance(v1, v2) :
    return math.sqrt((v1[0] - v2[0])**2 + (v1[1] - v2[1])**2)

# generates all subsets
def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def numberOfComponents(V, E) :
    return len(findComponents(V, E))

# Function used to find all components
# If only one component found - connected graph
def findComponents(V, E):
    components = []
    G = {}
    for vert in V:
        G[vert[0]] = []
    for edge in E:
        G[edge[0]].append(edge[1])
        G[edge[1]].append(edge[0])

    V_search = V
    while len(V_search) > 0:
        vert = V_search[0]
        nu_comp = expandComponent(vert[0], G)
        components.append(nu_comp)

        for c_vert in nu_comp:
            V_search.remove((c_vert,))

    return components

# Subfunction for findComponents
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




if __name__ == '__main__':
    path = './data/sensors.txt'

    data = read_data(path)
    print(len(data))
    print(data)
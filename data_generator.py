import math
import matplotlib.pyplot as plt
import numpy as np
from gudhi.alpha_complex import AlphaComplex
from mogutda import SimplicialComplex

from export_triangulation_to_ply import export_ply
from main import find_optimal_r, find_optimal_R


def plot_sphere_with_points(r, points):
    plt.rcParams["figure.figsize"] = [15.00, 10.0]
    plt.rcParams["figure.autolayout"] = True
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    # r = 0.05
    u, v = np.mgrid[0:2 * np.pi:30j, 0:np.pi:20j]
    x = r * np.cos(u) * np.sin(v)
    y = r * np.sin(u) * np.sin(v)
    z = r * np.cos(v)
    # ax.plot_surface(x, y, z, cmap=plt.cm.YlGnBu_r)
    ax.plot_surface(x, y, z, color='w', alpha=0.6)

    for point in points:
        ax.scatter(point[0], point[1], point[2], color='r')

    plt.show()


def sphere2cartesian(theta, phi, r):
    return [r * math.sin(theta) * math.cos(phi), r * math.sin(theta) * math.sin(phi), r * math.cos(theta)]


def generate_equidistributed_points_on_sphere(N, r, epsilon=0):
    # N = 56
    # r = 1
    # epsilon = 0.0001
    N_count = 0
    a = (4 * math.pi * r) / N
    d = math.sqrt(a)

    M_theta = round(math.pi / d)
    d_theta = math.pi / M_theta
    d_phi = a / d_theta
    generated_points = []
    for m in range(M_theta):
        theta = math.pi * (m + 0.5) / M_theta
        M_phi = round((2 * math.pi * math.sin(theta)) / d_phi)
        for n in range(M_phi - 1):
            phi = (2 * math.pi * n) / M_phi
            generated_points.append(sphere2cartesian(theta, phi, r + epsilon))
            N_count += 1
    return generated_points


# print(N_count)
# plot_sphere_with_points(1, generated_points)
points = generate_equidistributed_points_on_sphere(56, 1)
optimal_r, vr_cx = find_optimal_r(points, 0, 1, 0.5)
print(optimal_r)
optimal_R, alpha_cx = find_optimal_R(points, 0, 1, 0.5)
print(optimal_R)

ac = AlphaComplex(points)
stree = ac.create_simplex_tree(max_alpha_square=optimal_R ** 2)
cech_sxes = []
for simplex in stree.get_simplices() :
    cech_sxes.append(tuple(simplex[0]))

sphere = SimplicialComplex(simplices=cech_sxes)
betti0 = sphere.betti_number(0)
betti1 = sphere.betti_number(1)
betti2 = sphere.betti_number(2)
euler = sphere.euler_characteristics()
print(betti0)
print(betti1)
print(betti2)
print(euler)

data = export_ply(cech_sxes, points)

f = open("sphere_generated.ply", "w")
f.write(data)
f.close()

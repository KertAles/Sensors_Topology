import math
import matplotlib.pyplot as plt
import numpy as np
from main import find_optimal_r, find_optimal_R
from export_triangulation_to_ply import export_ply
from gudhi import AlphaComplex


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
#points = generate_equidistributed_points_on_sphere(56, 1)
#optimal_r = find_optimal_r(points, 0, 1)
#print(optimal_r)

# x = []
# y = []
# z = []
#
# for point in points:
#     x.append(point[0])
#     y.append(point[1])
#     z.append(point[2])
# print(x)
# print(y)
# print(z)

def sphericalcoordinate(x, y) :
    return [math.cos(x) * math.cos(y), math.sin(x) * math.cos(y), math.sin(y)]
def NX(n, x):
    pts = []

    start = (-1. + 1. / (n - 1.))
    increment = (2. - 2. / ( n - 1. ) ) / (n - 1.)
    for j in range(0 , n):
        s = start + j * increment
        pts.append(
        sphericalcoordinate(s * x, math.pi / 2. * math.copysign(1, s) * (1. - math.sqrt(1. - abs(s)))))
    return pts


def generate_points(n):
    return NX(n, 0.1 + 1.2 * n)


points = generate_points(50)
print(points)
#plot_sphere_with_points(1,points)

r, vr_cx = find_optimal_r(points, 0.25, 0.7, 0.1)
R, cech_cx = find_optimal_R(points, r/2, r*2, 0.1, eps=1e-6)

cech_sxes = cech_cx
'''sphere = SimplicialComplex(simplices=cech_sxes)
print(sphere.euler_characteristics())
print(sphere.betti_number(0))
print(sphere.betti_number(1))
print(sphere.betti_number(2))'''

ac = AlphaComplex(points)
stree = ac.create_simplex_tree(max_alpha_square=R ** 2)

cech_sxes = []
for simplex in stree.get_simplices() :
    cech_sxes.append(tuple(simplex[0]))

data = export_ply(cech_sxes, points)

f = open("fifty-points.ply", "w")
f.write(data)
f.close()
print(r)

import numpy as np
import random
import math
from statistics import mean
#from data_generator import plot_sphere_with_points
from main import find_optimal_r, find_optimal_R

def generate_random_points_on_sphere(n):
    points = []
    for i in range(n):
        theta = 2 * math.pi * random.uniform(0,1)
        phi = math.acos(1 - 2 * random.uniform(0,1))
        x = math.sin(phi) * math.cos(theta)
        y = math.sin(phi) * math.sin(theta)
        z = math.cos(phi)
        points.append((x,y,z))
    return points

rs = []
Rs = []

for i in range(100):
    points = generate_random_points_on_sphere(50)
    r, vr_cx = find_optimal_r(points, 0.25, 0.7, 0.1)
    R, cech_cx = find_optimal_R(points, r/2, r*2, 0.1, eps=1e-6)
    rs.append(r)
    Rs.append(R)

print("Average r: ", mean(rs))
print("Average R: ", mean(Rs))

#plot_sphere_with_points(1,points)
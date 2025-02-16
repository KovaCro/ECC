import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # For 3D plotting
import matplotlib.patches as patches

def plot_with_elliptic_curve_2d(a, b, range=[-5, 5]):
    y, x = np.ogrid[range[0]:range[1]:100j, range[0]:range[1]:100j]
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.contour(x.ravel(), y.ravel(), pow(y, 2) - pow(x, 3) - x * a - b, [0], colors='red')
    ax.grid()
    return plt

def plot_with_finite_elliptic_2d(a, b, p):
    pointsx = []
    pointsy = []
    for x in range(p):
        y_squared = (x**3 + a * x + b) % p
        for y in range(p):
            if (y**2) % p == y_squared:
                pointsx.append(x)
                pointsy.append(y)
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.grid()
    ax.scatter(pointsx, pointsy)
    ax.set_xlim(-p, p)
    ax.set_ylim(-p, p)
    return plt

def plot_elliptic_curve_3d(a, b, filename, bbox=(-2.5,2.5)):
    def elliptic_curve_3d(x, y, z, a, b):
      return z * y**2 - x**3 - a * z**2 * x - b * z**3
    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    A = np.linspace(xmin, xmax, 300)
    B = np.linspace(xmin, xmax, 15)
    A1,A2 = np.meshgrid(A,A)

    for z in B:
        X,Y = A1,A2
        Z = elliptic_curve_3d(X,Y,z, a, b)
        cset = ax.contour(X, Y, Z+z, [z], zdir='z', colors='red')

    for y in B:
        X,Z = A1,A2
        Y = elliptic_curve_3d(X,y,Z, a, b)
        cset = ax.contour(X, Y+y, Z, [y], zdir='y', colors='red')

    for x in B:
        Y,Z = A1,A2
        X = elliptic_curve_3d(x,Y,Z, a, b)
        cset = ax.contour(X+x, Y, Z, [x], zdir='x', colors='red')
    ax.set_zlim3d(zmin,zmax)
    ax.set_xlim3d(xmin,xmax)
    ax.set_ylim3d(ymin,ymax)

    plt.savefig(filename, dpi=300)


plot_with_elliptic_curve_2d(-1, 5, [-5, 5]).savefig("1.png")
plot_with_elliptic_curve_2d(-3, 0, [-4, 4]).savefig("2.png")
plot_elliptic_curve_3d(-3, 0, "3.png")

substraction_plot = plot_with_elliptic_curve_2d(-1, 5, [-5, 5])
substraction_plot.plot([0, 0], [2.236, -2.236], ls='--', c='black')
substraction_plot.plot([0, 0], [2.236, -2.236], 'o')
substraction_plot.text(0, 2.336, 'P', ha='center', va='bottom', fontsize='large')
substraction_plot.text(0, -2.436, '-P', ha='center', va='top', fontsize='large')
substraction_plot.savefig("4.png")

addition_plot_1 = plot_with_elliptic_curve_2d(-1, 5, [-5, 5])
addition_plot_1.axline([2.374, -4], [-1.796, -1], ls='-', c='black')
addition_plot_1.plot([2.374, -1.796, -0.06, -0.06], [-4, -1, -2.25, 2.25], 'o')
addition_plot_1.text(2.374, -4.3, 'P', ha='center', va='top', fontsize='large')
addition_plot_1.text(-1.796, -1.3, 'Q', ha='center', va='top', fontsize='large')
addition_plot_1.text(-0.06, -2.5, 'R', ha='center', va='top', fontsize='large')
addition_plot_1.plot([-0.06, -0.05], [-2.25, 2.25], ls='--', c='black')
addition_plot_1.text(-0.06, 2.5, 'P + Q', ha='center', va='bottom', fontsize='large')
addition_plot_1.savefig("5.png")

addition_plot_2 = plot_with_elliptic_curve_2d(-3, 0, [-5, 5])
addition_plot_2.axline([-1, 1.414], [2, 1.414], ls='-', c='black')
addition_plot_2.plot([-1, 2, 2], [1.414, 1.414, -1.414], 'o')
addition_plot_2.text(-1, 1.514, 'P', ha='center', va='bottom', fontsize='large')
addition_plot_2.text(2, 1.514, 'R', ha='center', va='bottom', fontsize='large')
addition_plot_2.plot([2, 2], [1.414, -1.414], ls='--', c='black')
addition_plot_2.text(1.9, -1.614, '2P', ha='center', va='top', fontsize='large')
addition_plot_2.savefig("6.png")

addition_finite_plot = plot_with_finite_elliptic_2d(1, 1, 5)
addition_finite_plot.text(3, 0.8, 'P', ha='center', va='top', fontsize='large')
addition_finite_plot.text(2, 3.8, 'Q', ha='center', va='top', fontsize='large')
addition_finite_plot.text(4, 1.8, 'P+Q', ha='center', va='top', fontsize='large')
addition_finite_plot.plot([-2.5, 2.5], [-5, 5], ls='--', c='orange')
addition_finite_plot.plot([2.5, 5], [-5, 0], ls='--', c='orange')
addition_finite_plot.plot([-5, 0], [-5, 5], ls='--', c='orange')
addition_finite_plot.plot([0, 5], [-5, 5], ls='--', c='orange')
addition_finite_plot.plot([4, 4], [3, 2], ls='--', c='red')
addition_finite_plot.savefig("7.png")
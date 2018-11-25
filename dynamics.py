import numpy as np

pi = np.pi
sin = np.sin
cos = np.cos

m = 3.5
J1 = 1.4
J2 = 1.6
l = 0.4
t = 2
q0 = 0.1 * pi
qf = pi / 6
g = 9.8


def getM(q):
    return np.array([
        [J1 + J2 + m * (l ** 2) + m * (q[1] ** 2) - m * 2 * l**2 / (sin(q[0])**2) * (1 - cos(q[0])), 0],
        [0, -m*l + m * 2 * l * (1 - cos(q[0])) / (sin(q[0]) ** 2)]
    ])



def getC(q, dq):
    return np.array([
        [2 * q[1] * dq[1] * m * dq[0] - m * (l**2) * 2 / (sin(q[0]) ** 3) * (dq[0]**2) * (1 - cos(q[0]))**2],
        [m * 2 * l / (sin(q[0])**3) * dq[0]**2 * (1-cos(q[0]))**2]
    ])


def getG(q):
    return np.array([
        [-m*g * (l/2 * cos(q[0]) + l * cos(q[0]) + q[1] * sin(q[0]))],
        [m * g * cos(q[0])]
    ])
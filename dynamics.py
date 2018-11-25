import numpy as np

pi = np.pi
sin = np.sin
cos = np.cos

m = 3.5
J1 = 1.4
J2 = 1.6
l = 0.4
tau = 2
q0 = 0.1 * pi
qf = pi / 6
g = 9.8

dt = 0.01

# only 1st joint angle
def getq(t):
    return q0 + (qf-q0) * (t/tau - 1 / (2*pi) * sin (2*pi*t/tau))


def getqdot(t):
    return (qf-q0) * (1 - cos(2*pi*t/tau)) / tau


def getqdotdot(t):
    return 2 * pi * (qf-q0) * sin (2*pi*t/tau) / tau


def getM(q):
    return np.array([
        [J1 + J2 + m * (l ** 2) + m * (q[1] ** 2) - m * 2 * l**2 / (sin(q[0])**2) * (1 - cos(q[0])), 0],
        [0, -m*l + m * 2 * l * (1 - cos(q[0])) / (sin(q[0]) ** 2)]
    ])



def getC(q, dq):
    return np.array([
        2 * q[1] * dq[1] * m * dq[0] - m * (l**2) * 2 / (sin(q[0]) ** 3) * (dq[0]**2) * (1 - cos(q[0]))**2,
        m * 2 * l / (sin(q[0])**3) * dq[0]**2 * (1-cos(q[0]))**2
    ])


def getG(q):
    return np.array([
        -m*g * (l/2 * cos(q[0]) + l * cos(q[0]) + q[1] * sin(q[0])),
        m * g * cos(q[0])
    ])


# offset of 2nd prismatic joint depending on angle of 1st joint (linear trajectory)
def getS(q):
    return 2 * l * (1 - cos(q)) / sin(q)

def getSdot(q, qdot):
    return 2 * l * qdot * (1 - cos(q)) / (sin(q)**2)

def getS2dot(q, qdot, q2dot):
    return 2 * l / (sin(q)**3) * (q2dot * sin(q) * (1-cos(q)) + (qdot**2) * (1-cos(q))**2)


import matplotlib.pyplot as plt

time = np.arange(0, tau + 0.005, dt)  # add 0.005 to create array with correct size
q1 = getq(time)
q1dot = getqdot(time)
q1dotdot = getqdotdot(time)

q = np.array([
    q1,
    getS(q1)
])

qdot = np.array([
    q1dot,
    getSdot(q1, q1dot)
])

q2dot = np.array([
    q1dotdot,
    getS2dot(q1, q1dot, q1dotdot)
])

forces = np.zeros(shape=(2, time.shape[0]))

for (i,), cur_time in np.ndenumerate(time):
    forces[:, i] = getM(q[:, i]).dot(q2dot[:, i]) + getC(q[:, i], qdot[:, i]) # + getG(q[:, i])

plt.xlabel('time, ms')
plt.title('force M')
plt.plot(time, forces[0, :])
plt.show()

plt.xlabel('time, ms')
plt.title('force P')
plt.plot(time, forces[1, :])
plt.show()

## simulation

q_real = np.zeros(shape=forces.shape)
qdot_real = np.zeros(shape=forces.shape)
q2dot_real = np.zeros(shape=forces.shape)

q_real[:, 0] = np.array([q0, getS(q0)])

for (i,), cur_time in np.ndenumerate(time):
    if i == 0:
        continue
    q_current = q_real[:, i-1]
    c_current = getC(q_current, qdot_real[:, i-1])
    g_current = 0 #getG(q_current)
    u = forces[:, i-1]
    M_inv = np.linalg.inv(getM(q_current))
    q2dot_real[:, i-1] = M_inv.dot(u - c_current - g_current)
    qdot_real[:, i] = qdot_real[:, i-1] + q2dot_real[:, i-1] * dt
    q_real[:, i] = q_real[:, i-1] + (qdot_real[:, i] + qdot_real[:, i-1]) / 2 * dt

plt.xlabel('time, ms')
plt.title('q1')
plt.plot(time, q[0, :], 'g')
plt.plot(time, q_real[0, :], 'r')
plt.show()

plt.xlabel('time, ms')
plt.title('q2')
plt.plot(time, q[1, :], 'g')
plt.plot(time, q_real[1, :], 'r')
plt.show()

plt.xlabel('time, ms')
plt.title('q2 velocity')
plt.plot(time, qdot[1, :], 'g')
plt.show()

plt.xlabel('time, ms')
plt.title('q2 acceleration')
plt.plot(time, q2dot[1, :], 'g')
plt.show()



from numpy import array, sin, cos, matmul, sqrt, amin, amax
# ---------- Load trajectory ----------

t = []
x = []
y = []
z = []
u_1 = []
u_2 = []
u_3 = []
u_4 = []
phi = []
theta = []
psi = []

trajectories = []

with open("output.txt", "r") as f:
    for line in f.readlines():
        if "nan" in line or "error" in line:
            print(line)
            exit()
        line = line.strip().split(",")
        if len(line) > 2:
            x.append(float(line[5]))
            y.append(float(line[6]))
            z.append(float(line[7]))
            phi.append(float(line[11]))
            theta.append(float(line[12]))
            psi.append(float(line[13]))
            t.append(float(line[0]))
            u_1.append(float(line[1]))
            u_2.append(float(line[2]))
            u_3.append(float(line[3]))
            u_4.append(float(line[4]))

with open("tmp.txt", "r") as f:
    for line in f.readlines():
        traj = []
        line = line.strip().split(",")
        for i in range(int((len(line) - 1)/3)):
            traj.append([float(line[3*i]),
                         float(line[3*i+1]),
                         float(line[3*i+2])])
        trajectories.append(traj)

dt = t[1] - t[0]

u_1 = array(u_1)
u_2 = array(u_2)
u_3 = array(u_3)
u_4 = array(u_4)

u_min = amin([amin(u_1), amin(u_2), amin(u_3), amin(u_4)])
u_max = amax([amax(u_1), amax(u_2), amax(u_3), amax(u_4)])

# ---------- Animate trajectory ----------

from vpython import *

scene = canvas(title='Quadcopter MPC',
                width=1200, height=900,
                center=vector(0,0,0), background=color.black)
scene.caption = "ADD CAPTION HERE LATER"

side = 10.0
thk = side / 100
s2 = 2*side - thk
s3 = 2*side + thk
floor = box(pos=vector(0, 0, -side), size=vector(s2, s2, thk), color = color.gray(0.7))

start = sphere(pos=vector(-2,-2,0), radius=0.05)
target = sphere(pos=vector(2,2,0), radius=0.05)

obstacles = []
obstacles.append( sphere(pos=vector(0,0,0), radius=0.2, color = color.red) )

box1 = box(color = color.red)
box2 = box(color = color.red)

thrust1 = arrow(shaftwidth=0.01)
thrust2 = arrow(shaftwidth=0.01)
thrust3 = arrow(shaftwidth=0.01)
thrust4 = arrow(shaftwidth=0.01)

scene.forward = vector(1.0, 1.0, -sqrt(2)/2)
scene.up = vector(0.0, 0.0, 1.0)
scene.camera.follow(box1)
scene.range = 2

time_scale = 1.0;

c = curve()

m = 0
while(True):
    n = m % len(t)

    center_mass = vector(x[n],y[n],z[n])

    rate(time_scale * 1.0 / dt)

    s_phi = sin(phi[n])
    s_theta = sin(theta[n])
    s_psi = sin(psi[n])

    c_phi = cos(phi[n])
    c_theta = cos(theta[n])
    c_psi = cos(psi[n])

    R = array([[c_phi*c_theta, c_phi*s_theta*s_psi - c_psi*s_phi, s_phi*s_psi + c_phi*c_psi*s_theta],
               [c_theta*s_phi, c_phi*c_psi + s_phi*s_theta*s_psi, c_psi*s_phi*s_theta - c_phi*s_psi],
               [-s_theta, c_theta*s_psi, c_theta*c_psi]])

    i_hat = matmul(R, array([1, 0, 0]))
    j_hat = matmul(R, array([0, 1, 0]))
    k_hat = matmul(R, array([0, 0, 1]))

    box1.pos = center_mass
    box1.axis = vector(i_hat[0], i_hat[1], i_hat[2])
    box1.up = vector(k_hat[0], k_hat[1], k_hat[2])
    box1.length = 0.45
    box1.width = 0.025
    box1.height = 0.025

    box2.pos = center_mass
    box2.axis = vector(j_hat[0], j_hat[1], j_hat[2])
    box2.up = vector(k_hat[0], k_hat[1], k_hat[2])
    box2.length = 0.45
    box2.width = 0.025
    box2.height = 0.025

    motor_dir = matmul(R, array([0.0, 0.0, 1.0]))
    motor_dir = vector(motor_dir[0], motor_dir[1], motor_dir[2])

    motor1_pos = matmul(R, array([0.225, 0.0, 0.0]))
    motor1_pos = center_mass + vector(motor1_pos[0], motor1_pos[1], motor1_pos[2])
    motor2_pos = matmul(R, array([0.0, 0.225, 0.0]))
    motor2_pos = center_mass + vector(motor2_pos[0], motor2_pos[1], motor2_pos[2])
    motor3_pos = matmul(R, array([-0.225, 0.0, 0.0]))
    motor3_pos = center_mass + vector(motor3_pos[0], motor3_pos[1], motor3_pos[2])
    motor4_pos = matmul(R, array([0.0, -0.225, 0.0]))
    motor4_pos = center_mass + vector(motor4_pos[0], motor4_pos[1], motor4_pos[2])

    thrust1.pos = motor1_pos
    thrust1.axis = motor_dir
    thrust1.length = u_1[n] / 10
    thrust2.pos = motor2_pos
    thrust2.axis = motor_dir
    thrust2.length = u_2[n] / 10
    thrust3.pos = motor3_pos
    thrust3.axis = motor_dir
    thrust3.length = u_3[n] / 10
    thrust4.pos = motor4_pos
    thrust4.axis = motor_dir
    thrust4.length = u_4[n] / 10

    c.clear()
    for i in range(len(trajectories[0])):
        c.append(vector(trajectories[n][i][0],
                        trajectories[n][i][1],
                        trajectories[n][i][2]))

    m += 1
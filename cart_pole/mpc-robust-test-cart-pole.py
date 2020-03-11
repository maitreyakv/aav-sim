import os, sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
from tqdm import tqdm

np.random.seed(0)

num_trials = 2000

err_mass_cart_min = -2.0
err_mass_cart_max = 2.0

err_mass_pole_min = -2.0
err_mass_pole_max = 4.0

err_length_pole_min = -0.2
err_length_pole_max = 0.6

err_mass_cart = np.random.uniform(err_mass_cart_min, err_mass_cart_max, num_trials)
err_mass_pole = np.random.uniform(err_mass_pole_min, err_mass_pole_max, num_trials)
err_length_pole = np.random.uniform(err_length_pole_min, err_length_pole_max, num_trials)

z_thresh = 0.1
theta_thres = np.deg2rad(1.0)

is_success = np.zeros(num_trials)

if "-c" in sys.argv:
    for i in tqdm(range(num_trials)):
        with open("input-robust-test.ini.in") as f:
            new_text = f.read()
            new_text = new_text.replace("@error_mass_cart", str(err_mass_cart[i]))
            new_text = new_text.replace("@error_mass_pole", str(err_mass_pole[i]))
            new_text = new_text.replace("@error_length_pole", str(err_length_pole[i]))

        with open("input-robust-test.ini", "w") as f:
            f.write(new_text)

        os.system("./cmake-build-debug/cartpole-sim input-robust-test.ini > /dev/null 2>&1")

        try:
            t, u, z, z_dot, theta, theta_dot = np.loadtxt("output.txt", delimiter=",", unpack=True)

            if (np.abs(z[-1] - 0.0) < z_thresh) and (np.abs(theta[-1] - np.pi) < theta_thres):
                is_success[i] = 1
            else:
                is_success[i] = 0
        except Exception as e:
            is_success[i] = -1

    plot_data = {"err_mass_cart": err_mass_cart,
                 "err_mass_pole": err_mass_pole,
                 "err_length_pole": err_length_pole,
                 "is_success": is_success}
    pickle.dump(plot_data, open("mpc-robust-plot-data.pkl", "wb"))
else:
    plot_data = pickle.load(open("mpc-robust-plot-data.pkl", "rb"))
    err_mass_cart = plot_data["err_mass_cart"]
    err_mass_pole = plot_data["err_mass_pole"]
    err_length_pole = plot_data["err_length_pole"]
    is_success = plot_data["is_success"]

marker_size = 20

plt.figure(figsize = (5.5,5))
plt.scatter(np.ma.masked_where(is_success != -1, err_mass_cart),
            np.ma.masked_where(is_success != -1, err_mass_pole),
            s=marker_size, marker="x", c="red")
plt.scatter(np.ma.masked_where(is_success != 0, err_mass_cart),
            np.ma.masked_where(is_success != 0, err_mass_pole),
            s=marker_size, marker="+", c="orange")
plt.scatter(np.ma.masked_where(is_success != 1, err_mass_cart),
            np.ma.masked_where(is_success != 1, err_mass_pole),
            s=marker_size, marker="o", c="green")

plt.xlim(err_mass_cart_min, err_mass_cart_max)
plt.ylim(err_mass_pole_min, err_mass_pole_max)
plt.xlabel("Error in cart mass [kg]", fontsize=16)
plt.ylabel("Error in pole mass [kg]", fontsize=16)
plt.xticks(size=14)
plt.yticks(size=14)
plt.tick_params(axis='both', which='major', direction='in', length=0, width=1)
plt.tick_params(axis='both', which='minor', direction='in', length=0, width=1)
plt.gcf().subplots_adjust(left=0.2)
plt.savefig("mpc-robust-cart-pole-mass-cart-mass-pole.png", dpi=500)

plt.figure(figsize = (5.5,5))
plt.scatter(np.ma.masked_where(is_success != -1, err_mass_cart),
            np.ma.masked_where(is_success != -1, err_length_pole),
            s=marker_size, marker="x", c="red")
plt.scatter(np.ma.masked_where(is_success != 0, err_mass_cart),
            np.ma.masked_where(is_success != 0, err_length_pole),
            s=marker_size, marker="+", c="orange")
plt.scatter(np.ma.masked_where(is_success != 1, err_mass_cart),
            np.ma.masked_where(is_success != 1, err_length_pole),
            s=marker_size, marker="o", c="green")
plt.xlim(err_mass_cart_min, err_mass_cart_max)
plt.ylim(err_length_pole_min, err_length_pole_max)
plt.xlabel("Error in cart mass [kg]", fontsize=16)
plt.ylabel("Error in pole length [m]", fontsize=16)
plt.xticks(size=14)
plt.yticks(size=14)
plt.tick_params(axis='both', which='major', direction='in', length=0, width=1)
plt.tick_params(axis='both', which='minor', direction='in', length=0, width=1)
plt.gcf().subplots_adjust(left=0.2)
plt.savefig("mpc-robust-cart-pole-mass-cart-length-pole.png", dpi=500)

plt.figure(figsize = (5.5,5))
plt.scatter(np.ma.masked_where(is_success != -1, err_mass_pole),
            np.ma.masked_where(is_success != -1, err_length_pole),
            s=marker_size, marker="x", c="red")
plt.scatter(np.ma.masked_where(is_success != 0, err_mass_pole),
            np.ma.masked_where(is_success != 0, err_length_pole),
            s=marker_size, marker="+", c="orange")
plt.scatter(np.ma.masked_where(is_success != 1, err_mass_pole),
            np.ma.masked_where(is_success != 1, err_length_pole),
            s=marker_size, marker="o", c="green")
plt.xlim(err_mass_pole_min, err_mass_pole_max)
plt.ylim(err_length_pole_min, err_length_pole_max)
plt.xlabel("Error in pole mass [kg]", fontsize=16)
plt.ylabel("Error in pole length [m]", fontsize=16)
plt.xticks(size=14)
plt.yticks(size=14)
plt.tick_params(axis='both', which='major', direction='in', length=0, width=1)
plt.tick_params(axis='both', which='minor', direction='in', length=0, width=1)
plt.gcf().subplots_adjust(left=0.2)
plt.savefig("mpc-robust-cart-pole-mass-pole-length-pole.png", dpi=500)

plt.show()

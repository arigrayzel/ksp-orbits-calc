import numpy as np
import sys

kerbin = {'a':1.360e10,'r':6e5,'T':9.203e6,'mu':3.5316e12, 'name':"Kerbin"}
eve = {'a':9.8327e9,'r':7e5,'T':5.658e6,'mu':8.172e12, 'name':"Eve"}
duna = {'a':2.0726e10,'r':3.2e5,'T':1.7315e7,'mu':3.0136e11, 'name':"Duna"}

bodies = {'kerbin':kerbin, 'eve':eve, 'duna':duna}
mu = 1.1723e18

def hohmann_transfer(origin, dest, r_p):
    try:
        origin = bodies[origin]
        dest = bodies[dest]
        r_p = float(r_p)*1000 + origin['r'] #parking orbit ought to be passed as km
    except:
        print("Unrecognized planet name(s). Format should be: \"trajectories hohmann <origin> <destination> <parking orbit> (km)\"")
        return

    if dest['a']>origin['a']:
        sign = -1
    else:
        sign = 1

    a_transfer = (dest['a']+origin['a'])/2
    t_transfer = np.pi/np.sqrt(mu/a_transfer**3)

    omega_o = 360/origin['T'] #degrees/s
    omega_d = 360/dest['T']

    omega_diff = omega_d-omega_o
    transfer_period = np.abs(360/omega_diff) #s

    delta_theta_transfer = 180 #Hohmann transfer orbit requires 180 degrees angular travel
    delta_theta_d = t_transfer/dest['T']*360
    delta_theta_launch = delta_theta_transfer - delta_theta_d
    print("To perform a Hohmann transfer from {} to {}, burn when the planets have a relative angle of {:.2f} degrees".format(origin["name"], dest["name"], delta_theta_launch))

############################## Delta-V calculations

    v_park = np.sqrt(origin['mu']/r_p)
    if sign==1:
        v_inf = np.sqrt(mu/origin['a'])*(1-np.sqrt(2*dest['a']/(origin['a']+dest['a'])))
    else:
        v_inf = np.sqrt(mu/origin['a'])*(np.sqrt(2*dest['a']/(dest['a']+origin['a']))-1)
    v_H = np.sqrt(v_inf**2 + 2*origin['mu']/r_p)
    delta_v_burn = v_H-v_park
    print("V_inf:{}".format(v_H))
    print("Other term:{}".format(v_park))
    print("To transfer intercept {}, apply a Delta-v of {:.2f} m/s".format(dest["name"], delta_v_burn))
if __name__ == "__main__":

    hohmann_names = ['h', 'hohmann', 'Hohmann']
    if sys.argv[1] in hohmann_names:
        hohmann_transfer(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        print("did not recognize command {}".format(sys.argv[2]))

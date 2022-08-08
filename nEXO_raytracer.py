'''
Simulates a number of Rn222 decays in the LXe cell

author: Nick Sanders

date: 8-8-2022

To run:
python nEXO_raytracer <n_decays> <file_id>
'''

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys
plt.rcParams['figure.figsize'] = [10, 10]
np.seterr(invalid='ignore')

n_sim = int(sys.argv[1])
run_id = int(sys.argv[2])

#Cell Dimensions
real_r = 2.915 * .0254  #in meters
real_h = 9.09 * .0254   #in meters
real_dim = [real_r, real_h]
PMT_dim = 0.5 * .0254   #Half side length of the PMT

#PTFE Properties
p_ref = 0.72  #Probability photon reflects off wall of chamber (in readme)

#LXe Properties
attenuation = 0.364  #Attenuation length in LXe
density = 3057  #kg/m3
LXe_mass = 8  #mass in kg in the cell
LXe_vol = LXe_mass / density
LXe_height = LXe_vol / (np.pi * real_r**2)
n_LXe = 1.69

#Misc.
p_det = 0.3  #Probability of photon creating photoelecron
'''
We use the decay of Rn222 as a test case. This can be adjusted to any desired decay just by scaling by the energy
'''
E_decay = 5590  #keV
QperE = 72 * 0.99  #quanta/keV * fraction in photons for alpha decay
Qperdecay = E_decay * QperE

def randLocation(dim):
    '''
    Picks a random location for a decay to occur (cylindrical coordinates)
    '''
    location = [dim[0]*np.sqrt(np.random.random()), 2*np.pi*np.random.random(), dim[1]*np.random.random()]
    return location

def randDirection():
    '''
    Picks a random direction for a photon to travel
    '''
    direction = [np.arccos(1 - 2*np.random.random()), 2*np.pi*np.random.random()]
    return direction

def nextSurface(x, y, dx, dy, radius):
    '''
    Finds the time when the ray next intersects the cylinder's rectangular surface (side of the cell)
    Geometry worked out here: https://www.cl.cam.ac.uk/teaching/1999/AGraphHCI/SMAG/node2.html
    '''
    a = dx**2 + dy**2
    b = (2 * x * dx) + (2 * y * dy)
    c = x**2 + y**2 - radius**2
    disc = np.sqrt(b**2 - (4 * a * c))
    t_int1 = (-b + disc) / (2 * a)  #'time' when ray intersects cylinder
    t_int2 = (-b - disc) / (2 * a)
    t_int = max(t_int1, t_int2)  #choose time that happens in the future (both never >0 as start would be outside cylinder)
    return t_int

def nextEdge(z, dz, height):
    '''
    if dz > 0:
    Finds the time when the ray next intersects the top/bottom of the cell
    '''
        target = height
    else:
        target = 0
    t_int = (target - z) / dz
    return t_int

def nextLG(z, dz, LXe_height):
    '''
    Finds the time when the ray next intersects the liquid/gas barrier. Returns infinity if moving away
    '''
    t_int = (LXe_height - z) / dz
    if t_int <= 0:  #photon is not moving toward LG barrier
        return np.inf
    return t_int

def updatePos(x, y, z, dx, dy, dz, travelled, t):
    '''
    Updates the current position of the photon and travelled distance
    '''
    xdist = dx * t
    ydist = dy * t
    zdist = dz * t
    travelled += np.sqrt(xdist**2 + ydist**2 + zdist**2)
    x += xdist
    y += ydist
    z += zdist
    return x, y, z, travelled

def hitCyl(x, y):
    '''
    Diffusely reflects photon off the wall of the cylinder
    '''
    phi = np.arctan(y/x)
    if x < 0:
        phi += np.pi  #Correction due to range of arctan being [-pi/2, pi/2]
    rand_theta = np.arccos(1 - 2*np.random.random())
    rand_phi = phi + np.pi/2 + np.pi*np.random.random()
    dx = np.cos(rand_phi) * np.sin(rand_theta)
    dy = np.sin(rand_phi) * np.sin(rand_theta)
    dz = np.cos(rand_theta)
    return dx, dy, dz

def reflectEdge(dz):
    '''
    Diffusely reflects photon off the top/bottom of LXe cell
    '''
    (rand_theta, rand_phi) = randDirection()
    dz_new = np.cos(rand_theta)
    if np.sign(dz_new) != np.sign(dz):  #Properly going the other way
        dx = np.cos(rand_phi) * np.sin(rand_theta)
        dy = np.sin(rand_phi) * np.sin(rand_theta)
        return dx, dy, dz_new
    else:
        dz = -dz_new
        dx = np.cos(rand_phi) * np.sin(rand_theta)
        dy = np.sin(rand_phi) * np.sin(rand_theta)
        return dx, dy, dz

def hitEdge(x, y, dx, dy, dz, PMT_dim):
    '''
    Checks if photon hits PMT, then reflects if not
    '''
    if abs(x) < PMT_dim and abs(y) < PMT_dim:  #ray has hit the PMT
        return dx, dy, dz, 1
    dx, dy, dz = reflectEdge(dz)
    return dx, dy, dz, 0

def refract(dx, dy, dz, in_L):
    '''
    Refracts photon through the liquid/gas boundary
    '''
    if in_L:  #moving from liquid to gas
        n1, n2 = n_LXe, 1
    else:  #moving from gas to liquid
        n1, n2 = 1, n_LXe
    prefactor = (n1 / n2) * np.sqrt(dx**2 + dy**2)
    if prefactor >= 1:  #total internal reflection
        dz = -dz
        return dx, dy, dz, in_L
    else:  #refracts
        atan = np.arctan(dy/dx)
        sx, sy, sz = np.sign(dx), np.sign(dy), np.sign(dz)
        dx = prefactor * np.cos(atan)
        dy = prefactor * np.sin(atan)
        dz = np.cos(np.arcsin(prefactor))
        #sign checks due to range of arctan
        if np.sign(dx) != sx:
            dx = -dx
        if np.sign(dy) != sy:
            dy = -dy
        if np.sign(dz) != sz:
            dz = -dz
        return dx, dy, dz, not in_L
    
def rayTrace(location, direction, dimensions, PMT_dim, p_ref):
    '''
    Traces the path of a single photon.
    Returns 0 if not detected
    Returns 1 if detected by bottom PMT
    Returns 2 if detected by top PMT
    '''
    (radius, height) = dimensions
    (r, phi, z) = location
    (polar, azimuth) = direction
    travelled = 0
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    dx = np.cos(azimuth) * np.sin(polar)
    dy = np.sin(azimuth) * np.sin(polar)
    dz = np.cos(polar)
    in_L = True  #Decays always happen in the liquid
    
    while True:
        t_cyl = nextSurface(x, y, dx, dy, radius)
        t_edg = nextEdge(z, dz, height)
        t_LG = nextLG(z, dz, LXe_height)
        if t_cyl < t_edg and t_cyl < t_LG:
            #REFLECT OFF CYLINDER
            if np.random.random() > p_ref:  #photon is absorbed by the wall
                return 0
            x, y, z, travelled = updatePos(x, y, z, dx, dy, dz, travelled, t_cyl)
            if travelled > attenuation:  #photon is absorbed by LXe
                return 0
            dx, dy, dz = hitCyl(x, y)
        elif t_edg < t_cyl and t_edg < t_LG:
            #CHECK IF IN PMT, REFLECT IF NOT
            x, y, z, travelled = updatePos(x, y, z, dx, dy, dz, travelled, t_edg)
            if travelled > attenuation:  #photon is absorbed by LXe
                return 0
            dx, dy, dz, detected = hitEdge(x, y, dx, dy, dz, PMT_dim)
            if detected == 1:
                if dz < 0:
                    return 1
                elif dz > 0:
                    return 2
                else:
                    raise Exception('Failed to find PMT')
            if np.random.random() > p_ref:  #photon is absorbed by the wall
                return 0
        else:
            #REFRACT OFF LG BOUNDARY
            x, y, z, travelled = updatePos(x, y, z, dx, dy, dz, travelled, t_LG)
            if travelled > attenuation:  #photon is absorbed by LXe
                return 0
            dx, dy, dz, in_L = refract(dx, dy, dz, in_L)

def loop(N_decay, name):
    '''
    Simulates N_decay number of Rn222 decays in the LXe cell
    Writes results to 3 files in directory ./data
    det1_counts{name}.npy  counts in bottom PMT
    det2_counts{name}.npy  counts in top PMT
    locs{name}.npy  cylindrical coordinates of the decay's location (r, phi, theta)
    '''
    decay_dim = [real_r, LXe_height]
    det1_counts = []
    det2_counts = []
    locs = []
    for i in tqdm(range(N_decay)):  #for each decay event
        location = randLocation(decay_dim)
        N_photons = int(Qperdecay)
        det1_count = 0
        det2_count = 0
        for j in range(N_photons):  #for each photon in the decay, choose a random direction and trace its path
            direction = randDirection()
            detected = rayTrace(location, direction, real_dim, PMT_dim, p_ref)
            if detected == 1:
                det1_count += 1
            elif detected == 2:
                det2_count += 1
        N_det1 = int(p_det * det1_count)
        N_det2 = int(p_det * det2_count)
        det1_counts.append(N_det1)
        det2_counts.append(N_det2)
        locs.append(location)
    np.save(f'./data/det1_counts{name}.npy', det1_counts)
    np.save(f'./data/det2_counts{name}.npy', det2_counts)
    np.save(f'./data/locs{name}.npy', locs)
 
if __name__ =='__main__':
    loop(n_sim, run_id)

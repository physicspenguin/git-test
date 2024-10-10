import numpy as np
import matplotlib
import matplotlib.pyplot as plt

class PowerClass():
  def __init__(self):
    pass

BSQ = 8
beam_div = 0.3  # mrad, 1/e2 points (full angle)
receiver_fov = beam_div
w = beam_div * 1  # mm (at 1 m range)
w_0 = 0.03  # mm
R = 1 # Range [m]
tot_power = 4  # W
eta_atm = 1  # efficiency of the atmosphere
eta_sys = 1  # efficiency of the system
D_r = 1  # [m] diameter of the receiver aperture

tot_power = 2*tot_power / (np.pi* w_0**2)  # reversing/undoing Eq. 2.3 in Carlsson

xs, ys = np.meshgrid(np.linspace(-receiver_fov/2, receiver_fov/2, 100),
                     np.linspace(-receiver_fov/2, receiver_fov/2, 100))

r2s = np.square(xs) + np.square(ys)
Is = (w_0/w) ** 2 * np.exp(-2*r2s/w**2)

def polar2cart(r, phi):
     return (r*np.sin(phi), r*np.cos(phi))


prev_radius_m = 0
prev_radius = 0
plt.figure()

# plt.imshow(Is, extent=(-receiver_fov/2, receiver_fov/2, -receiver_fov/2, receiver_fov/2),
#            cmap='viridis', clim=(0, np.max(Is)))
# cmap = matplotlib.cm.get_cmap('Accent')

power_subrays = []
n_subrays = []
area_subrays = []
rec_power_subrays = []

for i_circle in range(BSQ):
    n_sub = int(i_circle * np.pi * 2) if i_circle > 0 else 1  # special case for central subray
    angle_rho = receiver_fov/2 * (i_circle/(BSQ-0.5))  # subtract 0.5 so that the distance between
    ## the outermost circle and the boundary is the same as 1/2 the distance between the circles.
    for i_sub in range(n_sub):
        angle_phi = i_sub/n_sub * 2 * np.pi
        x,y = polar2cart(angle_rho, angle_phi)
        plt.scatter(x, y, color=cmap(i_circle/(BSQ+1)))

        if i_circle > 0:
            line_end = polar2cart(angle_rho + receiver_fov/4/(BSQ-0.5), angle_phi - 1/n_sub*np.pi)
            line_start = polar2cart(angle_rho - receiver_fov/4/(BSQ-0.5), angle_phi - 1/n_sub*np.pi)
            plt.plot([line_start[0], line_end[0]], [line_start[1], line_end[1]], 'k')

    radius = (angle_rho + receiver_fov/2 * (0.5/(BSQ-0.5)))
    circle = plt.Circle((0, 0), radius, fill=False)
    plt.gca().add_artist(circle)

    power_subray = -1 * np.pi * w_0 * w_0 / (2*n_sub) * (
        np.exp(-2/(w*w) * (radius)**2) -
        np.exp(-2/(w*w) * (prev_radius)**2)
    ) * tot_power

    radius_m = radius * R
    print(f"inner radius: {prev_radius_m:.2f}, outer radius: {radius_m:.2f}")
    area_subray = (np.pi * (radius_m**2 - prev_radius_m**2)  # outer circle - inner circle
                    / n_sub)  # divided by n
    print(f"area_subray: {area_subray:.4e} (x{n_sub:3d} = {area_subray*n_sub:.4e} ) ")
    prev_radius_m = radius_m  # save radius for the next circle
    prev_radius = radius  # save radius for the next circle
    sigma = area_subray  # here do something regarding incidence angle if appropriate

    # Laser/radar eq.
    beam_div_from_area = np.sqrt(area_subray) / np.pi / R * 2
    rec_power_subray = (power_subray * D_r ** 2) / (4 * np.pi * R ** 4 * beam_div_from_area ** 2) * eta_sys * eta_atm  * sigma

    power_subrays.append(power_subray)
    area_subrays.append(area_subray)
    n_subrays.append(n_sub)
    rec_power_subrays.append(rec_power_subray)

power_subrays = np.array(power_subrays)
rec_power_subrays = np.array(rec_power_subrays)
n_subrays = np.array(n_subrays)

# receiver fov
circle = plt.Circle((0, 0), receiver_fov / 2, fill=False, color='red', linestyle=':', linewidth=3)
plt.gca().add_artist(circle)
# beam div
circle = plt.Circle((0, 0), beam_div / 2, fill=False, color='blue', linestyle=':', linewidth=3)
plt.gca().add_artist(circle)

plt.axis('equal')
plt.xlim([-receiver_fov/2, receiver_fov/2])
plt.ylim([-receiver_fov/2, receiver_fov/2])
plt.title('Subrays and irradiance')
plt.tight_layout()
plt.show()


plt.figure()
plt.bar(range(len(area_subrays)), area_subrays*n_subrays, color='grey', label='Sum over subrays at circle')
plt.bar(range(len(area_subrays)), area_subrays, color='blue', label='Single subray at circle')
plt.title(f'Total area: {np.sum(area_subrays*n_subrays):.3f} [m²]\n'
           'Area per subray at different concentric levels')
plt.xlabel('Circle index')
plt.ylabel('Subray area')
plt.legend()
plt.tight_layout()
plt.show()

plt.figure()
plt.bar(range(len(power_subrays)), power_subrays*n_subrays, color='grey', label='Sum over areas at circle')
plt.bar(range(len(power_subrays)), power_subrays, color='blue', label='Single area at circle')
plt.title('Power per subray at different concentric levels')
plt.xlabel('Circle index')
plt.ylabel('Subray power')
plt.legend()
plt.tight_layout()
plt.show()

plt.figure()
plt.pie(power_subrays * n_subrays, labels=range(len(power_subrays)))
plt.suptitle('Power contribution from subray levels')
plt.title(f'Total power: {np.sum(power_subrays*n_subrays):.3f} [W]\n'
          f'(Receiver FoV multiplier: {receiver_fov/beam_div:.2f})')
plt.tight_layout()
plt.show()

print(power_subrays)


plt.figure()
plt.bar(range(len(rec_power_subrays)), rec_power_subrays*n_subrays, color='grey', label='Sum over subrays at circle')
plt.bar(range(len(rec_power_subrays)), rec_power_subrays, color='blue', label='Single subray at circle')
plt.title(f'Total received power: {np.sum(rec_power_subrays*n_subrays):.3e} [W]\n'
           'Received power per subray at different concentric levels')
plt.xlabel('Circle index')
plt.ylabel('Subray received power')
plt.legend()
plt.tight_layout()
plt.show()

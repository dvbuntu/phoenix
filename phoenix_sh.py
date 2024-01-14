import shapely as sh
# conda activate tf # this works, I guess
from matplotlib.collections import PatchCollection
from matplotlib.patches import Circle, Polygon, Wedge
import matplotlib.pyplot as plt
plt.ion()
import numpy as np
from shapely.plotting import plot_line, plot_points, plot_polygon
from skspatial.objects import LineSegment, Line

# should read this from a config?  But like, another python config?
a = 1
gap = a/4
mid_x = 0
mid_y = 5.5*a
tip_x = -4*a
tip_y = 4*mid_y/3
c_x = 2*a # right polyline here
c_y = 7*a/4
phi = np.pi/6
alph = np.pi/3
sig = np.pi/18
beta = np.pi-np.arctan2(tip_y-mid_y, tip_x-mid_x)

bottom = Line.from_points([-10*a,0], [10*a,0])
wing_tip = Line([tip_x,tip_y], [np.sin(sig), -np.cos(sig)])

def next_line(prev_line, bottom, wing_tip, space):
    next_line = prev_line.offset_curve(space, join_style=2, mitre_limit=gap)

    next_line_start = Line.from_points(next_line.coords[1], next_line.coords[0])
    next_line_int_bot = np.round(next_line_start.intersect_line(bottom),4)

    next_line_end = Line.from_points(next_line.coords[-2], next_line.coords[-1])
    next_line_int_tip = np.round(next_line_end.intersect_line(wing_tip),4)
    next_line =sh.LineString([tuple(next_line_int_bot)] + next_line.coords[1:-1] + [tuple(next_line_int_tip)]) 
    return next_line

patches = []

# White stripe
w1 = sh.LineString([(a,0), (c_x, c_y), (c_x, c_y+a/2),
                    (mid_x, mid_y), (tip_x, tip_y)])
w2 = next_line(w1, bottom, wing_tip, a)

W = Polygon(w1.coords[:] + list(reversed(w2.coords[:])),
        alpha=1.,
        color='#000000', # Black
        fill=False
        )

patches.append(W)

# pink 1
p11 = next_line(w2, bottom, wing_tip, gap)
p12 = next_line(p11, bottom, wing_tip, a)

P1 = Polygon(p11.coords[:] + list(reversed(p12.coords[:])),
        alpha=1.,
        color='#F5A9B8', # pink
        fill=True
        )

patches.append(P1)

# blue 1
blue_bot = Line.from_points([-10*a,a], [10*a,a])
b11 = next_line(p12, blue_bot, wing_tip, gap)
b12 = next_line(b11, blue_bot, wing_tip, a)
## adjust strip bottom chunk
b12_ext = Line.from_points(b12.coords[2], b12.coords[1])
b11_piece = Line.from_points(b11.coords[0], b11.coords[1])
b1_int = np.round(b12_ext.intersect_line(b11_piece),4)

B1 = Polygon(b11.coords[:] + list(reversed(b12.coords[1:])) + [tuple(b1_int)],
        alpha=1.,
        color='#5BCEFA', # blue
        fill=True
        )

patches.append(B1)

# head
head1 = next_line(w1, bottom, wing_tip, -gap)
head_base = head1.coords[5:7]
#head_base = sh.LineString(

p = PatchCollection(patches, match_original=True)
fig, ax  = plt.subplots()
ax.add_collection(p)
ax.set_xlim(-8*a,8*a)
ax.set_ylim(0,12*a)
plt.gca().set_aspect('equal')


sh.plotting.plot_line(head1)


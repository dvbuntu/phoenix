import shapely as sh
import shapely.ops
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
tip_x = -3.6*a
tip_y = 5*mid_y/4
c_x = 2.15*a # right polyline here
c_y = 8*a/4
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
## redo the lines to hit at this intersection
blue_bot = Line.from_points([-10*a,b1_int[1]], [10*a,b1_int[1]])
b11 = next_line(p12, blue_bot, wing_tip, gap)
b12 = next_line(b11, blue_bot, wing_tip, a)

B1 = Polygon(b11.coords[:] + list(reversed(b12.coords[1:])) + [tuple(b1_int)],
        alpha=1.,
        color='#5BCEFA', # blue
        fill=True
        )

patches.append(B1)



# head
#head_bot = Line.from_points([-10*a,c_y+a/4], [10*a,c_y+a/4])
#head_turn = (head_bot.point[1] + mid_y)/2
#head_top_cutoff = Line.from_points([-10*a,head_turn], [10*a,head_turn])
w1_piece = sh.LineString(w1.coords[2:4])
h1 = w1_piece.offset_curve(-gap, join_style=2, mitre_limit=gap)
h1s = sh.ops.substring(h1, 0, 0.5, normalized=True)

h1sl = sh.LineString([h1s.coords[1], [h1s.coords[1][0]-gap/2, h1s.coords[1][1]]])
h1sd = sh.LineString([h1s.coords[0], [h1s.coords[0][0], h1s.coords[0][1]-gap/2]])

## get normal direction
d0 = Line.from_points(*w1_piece.coords[:]).direction
d = np.array([-d0[1],d0[0]])

h1sn1 = Line(h1sl.coords[1], d/np.sqrt(np.sum(d**2)))
h1sn1 = sh.LineString([h1sn1.point, h1sn1.to_point(-a)])

h1sn2 = Line(h1sd.coords[1], d/np.sqrt(np.sum(d**2)))
h1sn2 = sh.LineString([h1sn2.point, h1sn2.to_point(-a/2)])

h1sb = Line(h1sn2.coords[1], -d0/np.sqrt(np.sum(d0**2)))
h1sb = sh.LineString([h1sb.point, h1sb.to_point(-a/2)])

h1st = sh.LineString([h1sn1.coords[1], [h1sb.coords[1][0], h1sn1.coords[1][1]]])
h1sr = sh.LineString([[h1sb.coords[1][0], h1sb.coords[1][1]],
            [h1sb.coords[1][0], h1sn1.coords[1][1]]])

if False:
	# debugging...
	colors = plt.color_sequences['tab10']
	for i,hh in enumerate([h1s,
					h1sl,
					h1sn1,
					h1sd,
					h1st,
					h1sr,
					h1sb,
					h1sn2]):
		sh.plotting.plot_line(hh, color=colors[i])


h1s = h1s.union(h1sl)
h1s = h1s.union(h1sn1)
h1s = h1s.union(h1sd)
h1s = h1s.union(h1st)
h1s = h1s.union(h1sr)
h1s = h1s.union(h1sb)
h1s = h1s.union(h1sn2)
h1s = sh.ops.linemerge(h1s.normalize())

P2 = Polygon(sh.get_coordinates(h1s),
        alpha=1.,
        color='#F5A9B8', # pink
        fill=True
        )
patches.append(P2)

# talon
p3t = [(c_x+a/2, c_y+a/4)]
p3t.append((p3t[-1][0]+a/2, p3t[-1][1]-2*a/3))
p3t.append((p3t[-1][0],     p3t[-1][1]-a))
p3t.append((p3t[-1][0]-a/4, 0))
p3t.append((p3t[0][0]-a/3, p3t[-2][1]))
p3t.append((p3t[-1][0], p3t[-4][1]))
p3s = sh.LineString(p3t)

P3 = Polygon(p3s.coords[:],
        alpha=1.,
        color='#F5A9B8', # pink
        fill=True
        )
patches.append(P3)

        
# beak
b2_pt = [(h1st.coords[-1][0]+gap, h1st.coords[-1][1]+gap)]

ddir = d/np.sqrt(np.sum(d**2))
d0dir = d0/np.sqrt(np.sum(d0**2))
b2t1 = Line(b2_pt[-1], ddir)
b2t1 = sh.LineString([b2t1.point, b2t1.to_point(-a/2)])

b2t1h = Line(b2t1.coords[-1], np.array([1,0]))
b2t1h = sh.LineString([b2t1h.point, b2t1h.to_point(a/2)])

b2t1d = Line(b2t1h.coords[-1], d0dir)
b2t1d = sh.LineString([b2t1d.point, b2t1d.to_point(a/2)])


#h1s = h1s.union(h1sl)
#h1s = h1s.union(h1sn1)
#h1s = h1s.union(h1sd)

b2_pt.append(b2t1.coords[-1])
p3t.append((b2_pt[-1][0]+a/2, b2_pt[-1][1]-2*a/3))




b2s = sh.LineString(b2_pt)



# plot the patches
p = PatchCollection(patches, match_original=True)
fig, ax  = plt.subplots()
ax.add_collection(p)
ax.set_xlim(-8*a,8*a)
ax.set_ylim(0,12*a)
plt.gca().set_aspect('equal')


#sh.plotting.plot_line(h1)
#sh.plotting.plot_line(p3s,color='red')
sh.plotting.plot_line(b2s,color='red')

plt.savefig('latest.png')



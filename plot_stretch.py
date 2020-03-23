import matplotlib.pyplot as plt
import numpy as np

#Colors
c_red = '#BC0000'
c_green = '#10A100'
c_blue = '#1000FF'
c_purple = '#C200CD'
c_cyan = '#00BBEF'
c_yellow = '#D6D600'
label_s = 7 #Label size. For ticks
legend_s = 5 #legend size
line_w = 1 #linewidth
mark_s = 2 #markersize
path = "elastic_plot.out"
dat = np.loadtxt(path)

############################################################################
#plot stretch
macolor = c_blue #Our colour

plt.plot( dat[:,0], dat[:,1], color=macolor) 
plt.fill_between( dat[:,0], dat[:,1] - dat[:,2],  dat[:,1] + dat[:,2],
alpha=0.5, edgecolor=macolor, facecolor=macolor, linewidth=0)


##############################################################################

plt.grid(True)
plt.xlabel('sub-length (bp)')
plt.ylabel('Stretch (pN)')
#plt.legend()

plt.show()


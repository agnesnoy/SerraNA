import matplotlib.pyplot as plt
import numpy as np

width=4
height=6

#Colors
c_red = '#BC0000'
c_green = '#10A100'
c_blue = '#1000FF'
c_purple = '#C200CD'
c_cyan = '#00BBEF'
c_yellow = '#D6D600'

label_s = 9  #Label size. For ticks
legend_s = 6 #legend size
title_s = 9  #Title size
axis_s = 7   #Axes sizes
line_w = 1   #linewidth
mark_s = 3    #markersize

#Array with colours to help iterate.
colours=[c_red, c_green, c_blue, c_purple, c_cyan, c_yellow]

#PLOTING SECTION
############################################################################
fig, axs = plt.subplots(3, figsize=(width,height), dpi=300, tight_layout=True)


#Plot bending angles at different lengths
#----------------------------------------------------------------------------------
column=17  #This is the column for the bending angle
filedat="structural_"

l=2 #Initial
step=10
for s in range(0,3):
    newpath=filedat+str(l)+"mer.out"
    macolor = colours[s] #Our colour

    malabel = str(l)+"mer"
    dat = np.loadtxt(newpath)
    axs[0].plot( dat[:,0], dat[:,column], "-o", markersize=mark_s, lw=line_w, color=macolor,label=malabel) 
    l = l + step

#Plot twist for 4 different regions [a,b]
#----------------------------------------------------------------------------------
column=3  #This is the column for the twist elastic constant
filedat="elastic_"

l1=6 #Length
l2=5 #diference
l=-l2
for s in range(1,5):
    l = l + l1
    newpath=filedat+"[" +str(l)+":"+str(l+l2)+"].out"
    macolor = colours[s] #Our colour

    malabel = "[" +str(l)+":"+str(l+l2)+"]"
    dat = np.loadtxt(newpath)
    axs[1].plot( dat[:,0], dat[:,column], "k", lw=line_w, color=macolor,label=malabel)
    axs[1].fill_between( dat[:,0], dat[:,column] - dat[:,column+1],  dat[:,column] + dat[:,column+1],
    alpha=0.25, edgecolor=macolor, facecolor=macolor, linewidth=0)


#Plot stretch modulus of the whole fragment
#----------------------------------------------------------------------------------
column=1  #This is the column for the stretch modulus
filedat="elastic_plot.out"

macolor = c_blue #Our colour

malabel = "B"
dat = np.loadtxt(filedat)
axs[2].plot( dat[:,0], dat[:,column], "-o", markersize=mark_s, lw=line_w, color=macolor,label=malabel)
axs[2].fill_between( dat[:,0], dat[:,column] - dat[:,column+1],  dat[:,column] + dat[:,column+1],
alpha=0.25, edgecolor=macolor, facecolor=macolor, linewidth=0)

##############################################################################

#Let's sort the labels
for s in range(0,3):
    if s == 0:
       labely="Bending angle (deg)"
       labelx='Position along the DNA'
    elif s == 1:
       labely="Twist (nm)"
       labelx='mer (bp)'
    elif s == 2:
       labely="Stretch (pN)"
       labelx='mer (bp)'

    axs[s].set_xlabel(labelx,fontsize=label_s)
    axs[s].set_ylabel(labely,fontsize=label_s)
    axs[s].tick_params(axis='both', which='both', labelsize=axis_s)
    axs[s].legend(loc='best',fontsize=legend_s)

#plt.show()
plt.savefig("my_first_result.pdf",dpi=300)

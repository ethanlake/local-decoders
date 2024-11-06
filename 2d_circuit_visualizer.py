import numpy as np
import matplotlib.pyplot as plt
import h5py 
import argparse
import sys 
from matplotlib.widgets import Slider
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation

parser = argparse.ArgumentParser() 
parser.add_argument('-hist',default='no') # spacetime history of spins to be drawn (optional)
parser.add_argument('-gates',default='no') # draws the gates given by a particular input h5 file (optional)
parser.add_argument('-save_animation',default='no') # saves animation for max_animation_time (defined below in the code) if not equal to "no" 
args = parser.parse_args() 

def loaddata(fin):
    with h5py.File(fin, "r") as f:
        key_list = [key for key in f.keys() if not key.startswith("_")] # saving dictionaries to jld2 files can sometimes cause the creation of additional keys describing the data structure that we don't care about; this gets rid of them 
        data_dict = dict.fromkeys(key_list)
        for key in key_list:
            data_dict[key] = f[key][()]  
    return data_dict 

### asthetics ###
gr = np.array([94, 94, 94])/255 # gray 
lgr = np.array([161, 161, 161])/255
mycol = np.array([33, 61, 204])/255 # dark blue 
mxcol = np.array([109, 222, 247])/255 # light blue 
tcol = np.array([255, 56, 56])/255 # red 
rcol = np.array([196, 162, 51])/255 # orange
idcol = np.array([51, 184, 87])/255 # green 
stringcol = tcol 
anyoncol = gr
anyonsize = 200
stringlw = 6
marg = 0.2
bmarg = marg 

wallcol = gr
awall = 1.
lwwall = .5
ms = 1.5
lw = 5
buff = .2 
alpha = .15 # of the drawn gates 

if args.gates != 'no': 
    data = loaddata(args.gates)
    n = data["n"]
    l = data["l"]
    L = data["L"]
    gate = str(data["gate"])
    gates = data["gates"].T 
    maxtime = np.max(gates[:,-1])+1
    ngates = np.shape(gates)[0]

    height = np.max(gates[:,-1]) # depth equal to maximum local time
    width = np.max(gates[:,1])

if args.hist != 'no': 
    hist_data = loaddata(args.hist)
    if args.gates == 'no':
        maxtime = np.shape(hist_data["err_hist"].T)[0]
        L = hist_data["L"]
    noise_hist = hist_data["noise_hist"].T 
    err_hist = hist_data["err_hist"].T  
    synd_hist = hist_data["synd_hist"].T
    periods = hist_data["periods"]
    if args.gates != 'no':
        full_gates = np.copy(gates)
        if periods > 1: # if there are more periods in the history data than in the gate data, repeat the gate data so that it has the same length as the history data 
            for per in range(1,periods): 
                time_shifted_gates = np.copy(gates)
                time_shifted_gates[:,-1] += per*(maxtime-1)
                full_gates = np.vstack((full_gates,time_shifted_gates))
        gates = full_gates
        maxtime = periods * maxtime 

sz = 6.5
fig, ax = plt.subplots(figsize=(sz,sz))

if args.save_animation == "no":
    plt.subplots_adjust(bottom=0.25) 
    slider_ax = plt.axes([0.2, 0.1, 0.65, 0.03])  # [left, bottom, width, height]
else: 
    slider_ax = plt.axes([-.5, -.5, 0, 0])  # [left, bottom, width, height]
plot_slider = Slider(slider_ax,r'',0,maxtime-1,valinit=0, valstep=1)

# Function to update plot
def update(val):
    # Clear current shapes
    ax.clear()
    ax.axis('off') 
    ax.set_xlim(-marg,L-marg)
    ax.set_ylim(-marg,L-marg)
    draw_grid = True #args.gates != "no"
    if draw_grid: 
        for i in range(L): 
            lw = 4 if i%3 == 0 else 1.5
            lw = 5 if i%9 == 0 else (2.5 if i%3 == 0 else 1.)
            ax.axvline(i,alpha=1,color=lgr,zorder=-1,lw=lw)
            ax.axhline(i,alpha=1,color=lgr,zorder=-1,lw=lw)
    showwalls = args.hist == 'no'
    if showwalls: 
        ax.plot([-buff,L+buff],[-buff,-buff],lw=lw,color=wallcol)
        ax.plot([-buff,L+buff],[L+buff,L+buff],lw=lw,color=wallcol)
        ax.plot([-buff,-buff],[-buff,L+buff],lw=lw,color=wallcol)
        ax.plot([L+buff,L+buff],[-buff,L+buff],lw=lw,color=wallcol)

        ax.plot([L/3,L/3],[-buff,L+buff],lw=lw/2,color=wallcol)
        ax.plot([2*L/3,2*L/3],[-buff,L+buff],lw=lw/2,color=wallcol)
        ax.plot([-buff,L+buff],[L/3,L/3],lw=lw/2,color=wallcol)
        ax.plot([-buff,L+buff],[2*L/3,2*L/3],lw=lw/2,color=wallcol)

    time = int(plot_slider.val)+1

    if args.gates != 'no': 
        thisgates = gates[gates[:,-1] == time,:]  # gates at chosen time slice 
        nthisgates = len(thisgates)

        for i in range(nthisgates):
            thisgate = thisgates[i]
            gtype = thisgate[0]
            go = thisgate[1]
            gx = thisgate[2]-1
            gy = thisgate[3]-1

            if gtype == 1: # identity gate 
                dx = 1
                dy = 1 
                col = idcol 
                tup = (gx,gy)

            if gtype == 2: # swap gate 
                dx = (2 if go == 1 else 1) 
                dy = (2 if go == 2 else 1)
                col = tcol 
                tup = (gx+buff/2,gy+buff/2-.5) if go == 1 else (gx+buff/2-.5,gy+buff/2)
            if gtype == 3: # m gate 
                dx = (3 if go == 1 else 1) 
                dy = (3 if go == 2 else 1)
                col = mycol 
                tup = ((gx+buff/2,gy+buff/2-.5) if go == 1 else (gx+buff/2-.5,gy+buff/2))
            if gtype == 4: # R0 
                dx = 3 
                dy = 3
                col = rcol 
                tup = (gx,gy)
            
            if gtype != 1: 
                rect = patches.Rectangle(tup, dx-buff, dy-buff, linewidth=2, edgecolor='none', facecolor=col,alpha=alpha)
                ax.add_patch(rect)

    if args.hist != 'no': 
        print("synd sums = ",np.sum(synd_hist[-1,:,:]))
        print("synd sums (end) = ",np.sum(synd_hist[-2,:,:]))
        print("size(synd_hist) = ",np.shape(synd_hist))
        thisxerrs = err_hist[time-1,:,:,0]
        for errind in np.argwhere(thisxerrs == 1): 
            if errind[0] == L-1: # on the boundary 
                ax.plot([errind[0],errind[0]+1-bmarg],[errind[1],errind[1]],lw=stringlw,color=stringcol)
                ax.plot([-bmarg,0],[errind[1],errind[1]],lw=stringlw,color=stringcol)
            else:
                ax.plot([errind[0],errind[0]+1],[errind[1],errind[1]],lw=stringlw,color=stringcol)

        thisyerrs = err_hist[time-1,:,:,1]
        for errind in np.argwhere(thisyerrs == 1): 
            if errind[1] == L-1: # on the boundary 
                ax.plot([errind[0],errind[0]],[errind[1],errind[1]+1-bmarg],lw=stringlw,color=stringcol)
                ax.plot([errind[0],errind[0]],[-bmarg,0],lw=stringlw,color=stringcol)
            else:
                ax.plot([errind[0],errind[0]],[errind[1],errind[1]+1],lw=stringlw,color=stringcol)

        thissynds = synd_hist[time-1,:,:]
        for syndind in np.argwhere(thissynds == 1): 
            ax.scatter(syndind[0],syndind[1],s=anyonsize,color=anyoncol,alpha=1,zorder=100)

        print("t, sum(anyons) = ",time-1,np.sum(synd_hist[time-1,:,:]))
    fig.canvas.draw_idle()

if args.save_animation != 'no':
    fps = 14
    dpi = 125 

    def animate(frame):
        plot_slider.set_val(frame)  
        update(frame)
    
    # max animation time (in units of automaton time steps) needs to be controlled by hand:
    max_animation_time = 330
    animation = FuncAnimation(fig, animate, frames=np.arange(0, max_animation_time))
    
    animation.save(args.save_animation, writer='ffmpeg', dpi=dpi,fps=fps)

else: 
    # allow the slider to be controlled by arrow keys 
    def on_key(event):
        if event.key == 'left':
            plot_slider.set_val(max(plot_slider.val - 1, plot_slider.valmin))
        elif event.key == 'right':
            plot_slider.set_val(min(plot_slider.val + 1, plot_slider.valmax))
    fig.canvas.mpl_connect('key_press_event', on_key)

    # connect the slider to the update function
    plot_slider.on_changed(update)

    plt.show()


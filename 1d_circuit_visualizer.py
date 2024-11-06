import numpy as np
import matplotlib.pyplot as plt
import h5py 
import argparse
import matplotlib.cm as cm 
import sys 
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman'] + plt.rcParams['font.serif']
pi = np.pi 

cmap = cm.coolwarm

parser = argparse.ArgumentParser() 
parser.add_argument('-fin',default="data/gates.h5") # a list of gates comprising the circuit to be plotted 
parser.add_argument('-history',default="no") # (optional) spacetime spin history to be overlaid on circuit diagram 
parser.add_argument('-save',default='no') # (optional) name of image file to be saved  
parser.add_argument('-nocircuit',action='store_true') # if true, doesn't draw the circuit architecture 
parser.add_argument('-symbolic',action='store_true') # if true, plots circuit elements as symbolic blocks. if false, plots the "microscopic" circuit (dots and lines)

args = parser.parse_args() 

def loaddata(fin):
    with h5py.File(fin, "r") as f:
        key_list = list(f.keys())
        data_dict = dict.fromkeys(key_list)
        for key in key_list:
            data_dict[key] = f[key][()]  

    return data_dict 

### asthetics ### 
# define some colors 
gray = np.array([94, 94, 94])/255  
black = np.array([0,0,0])/255  
dark_red = np.array([230, 35, 21])/255 
light_green = np.array([81, 204, 33])/255 
dark_blue = np.array([33, 61, 204])/255  
light_blue = np.array([109, 222, 247])/255 
red = np.array([255, 56, 56])/255 

# colors of various gates 
xcol = black # identity 
ycol = dark_blue # swap 
zcol = dark_red # EC 
wallcol = gray # dividers between gates 

# thicknesses and opacities 
awall = 1.
lwwall = .5
lw = 1.
ms = 1.5
buff = .2 # filler space between walls and shapes 

### load data about spacetime spin history ### 
if args.history != "no":
    hist_data = loaddata(args.history)
    try:
        history = hist_data["history"].T 
    except: 
        history = hist_data["average_history"].T
    print(history)

### load data about circuit to be plotted ### 
data = loaddata(args.fin)
L = data["L"]; l = data["l"]; n = data["n"]

gates = data["gates"].T 
ngates = np.shape(gates)[0]

height = np.max(gates[:,-1]) # depth equal to maximum local time
width = np.max(gates[:,1])
rat = width/height 
sz = 50 if l == 3 else (8 if l == 2 else 5) # works well for n = 3 
sz = 10

fig,ax = plt.subplots(figsize=(2,8))
ax.axis('off') 

### draw the circuit ### 
if not args.nocircuit: 
    local_times = np.zeros(L)
    time = 0

    # draw dividing lines between gates if desired 
    if args.symbolic: 
        ax.plot([-.5,-.5],[0,height],lw=lwwall,alpha=awall,c=wallcol)
        for i in range(height+1):
            ax.plot([-.5,L-.5],[i,i],lw=lwwall,alpha=awall,c=wallcol)

    # draw the gates 
    for i in range(ngates):
        gate = gates[i]
        tloc = gate[-1]-1
        if tloc < height+1: 
            gloc = gate[1]-1
            if args.symbolic:  
                if gate[0] == 1: # identity gate  
                    ax.plot([gloc,gloc],[tloc+buff,tloc+1-buff],c=xcol,ms=ms,lw=lw)
                    ax.plot([gloc+.5,gloc+.5],[tloc,tloc+1],c=wallcol,ms=ms,lw=lwwall,alpha=awall)
                if gate[0] == 2: # swap gate 
                    # looks like an X: 
                    ax.plot([gloc+buff,gloc+1-buff],[tloc+buff,tloc+1-buff],c=ycol,ms=ms,lw=lw)
                    ax.plot([gloc+1-buff,gloc+buff],[tloc+buff,tloc+1-buff],c=ycol,ms=ms,lw=lw)
                    ax.plot([gloc+1+.5,gloc+1+.5],[tloc,tloc+1],c=wallcol,lw=lwwall,alpha=awall)

                    # looks like a solid colored block: 
                    # ax.plot([gloc+2-1+.5,gloc+2-1+.5],[tloc,tloc+1],c=wallcol,ms=ms,lw=lwwall,alpha=awall)
                    # rect = patches.Rectangle((gloc-.5+buff,buff+tloc), 2-2*buff, 1-2*buff, linewidth=2, edgecolor='none', facecolor=ycol)
                    # ax.add_patch(rect)
                if gate[0] == 3: # EC gate 
                    ax.plot([gloc+n-1+.5,gloc+n-1+.5],[tloc,tloc+1],c=wallcol,ms=ms,lw=lwwall,alpha=awall)
                    rect = patches.Rectangle((gloc-.5+buff,buff+tloc), n-2*buff, 1-2*buff, linewidth=2, edgecolor='none', facecolor=zcol)
                    ax.add_patch(rect)

            else: # "wireframe" drawings 
                if gate[0] == 1: # identity  
                    ax.plot([gloc,gloc],[tloc,tloc+1],marker='o',c=xcol,ms=ms,lw=lw)
                if gate[0] == 2: # swap  
                    ax.plot([gloc,gloc+1],[tloc,tloc+1],marker='o',c=ycol,ms=ms,lw=lw)
                    ax.plot([gloc+1,gloc],[tloc,tloc+1],marker='o',c=ycol,ms=ms,lw=lw)
                if gate[0] == 3: # EC  
                    for j in range(n): 
                        for k in range(n): 
                            if j == k:
                                a = 1. 
                            else: 
                                a = .6
                            ax.plot([gloc+j,gloc+k],[tloc,tloc+1],marker='o',c=zcol,ms=ms,lw=lw,alpha=a)

if args.history != "no":
    # optionally show the places at which noise occurs: 
    # noise_hist = hist_data["noise_hist"].T 
    # print("plotting noise: ",noise_hist)
    # noisecmap = ListedColormap(['none', 'yellow'])  # 'none' for transparent, 'red' for red color
    # ax.imshow(noise_hist,cmap=noisecmap,origin='lower')

    historycmap = ListedColormap(['none', 'black'])  
    ax.imshow(history,cmap=historycmap,alpha=.45,origin='lower')

if args.save != 'no':
    print("saving...")
    plt.savefig(args.save, dpi=None,orientation='portrait', bbox_inches='tight', pad_inches=0.1,facecolor='w',edgecolor='w')
else: 
    plt.show()

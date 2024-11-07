import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm 
from matplotlib.widgets import Slider

# these need to first be defined outide of any function since we are treating them as global later 
nu = 1; beta = .125; gamma = 1.75

def scaling_plotter(data,plot,xc,nu0=1,gamma0=1.75,beta0=.125,raw=False,d=2,title=""): 
    """
    makes interactive scaling plots for various critical exponents. 
    data: a dictionary with the following entries 
        * "Ls": system sizes 
        * "xs": raw "temperatures" (not reduced temperatures) for each system size (#Ls x #(x data points) matrix)
        * "chis": susceptibilities 
        * "mags": magnetizations 
        * "binds": binder cumulants 
    plot: quantitiy to plot, one of ["binds" "mags" "chis"]
    xc: value of critical point 
    nu,gamma,beta: estimates for critical exponents 
    raw: if true, just plots the raw data without scaling
    d: dimension; used only for checking hyperscaling relations 
    title: plot title 
    critical exponents tuned using key presses. left/right tunes nu, right/left tunes gamma, and ,/. tunes beta 
    """

    ### asthetics for plots ## 
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Computer Modern Roman'] + plt.rcParams['font.serif']
    mew = 1.1
    ms = 7.5
    lw = 2.2 
    cmap = cm.coolwarm 

    fig,ax = plt.subplots(figsize=(2.5,2.5))
    ax.minorticks_off()

    Ls = data["Ls"]; chis = data["chis"]; mags = data["mags"]; binds = data["binds"]

    ### plot raw unscaled data ### 
    if raw: 
        for l in range(len(Ls)): 
            col = cmap((l+1) / len(Ls))
            L = Ls[l]

            xs = data["xs"][l]; ts = (xs - xc)/xc 
            ys = data[plot][l] # stuff to plot at this system size 

            ax.plot(xs,ys,c=col,label=r'$%d$'%Ls[l],marker='o',mfc='w',mew=mew,ms=ms,lw=lw)
        
        ax.set_xlabel(r'$t$')
        if plot == "mags": 
            ax.set_ylabel(r'$M$')
        elif plot == "chis": 
            ax.set_ylabel(r'$\chi$')
        elif plot == "binds": 
            ax.set_ylabel(r'$B$')

    ### plot scaled data (in an interactive way) ### 
    else:  
        lines0 = [] # will hold the plot lines; dynamically updated 

        # draw lines with guessed values of critical exponents 
        for l in range(len(Ls)): 
            col = cmap((l+1) / len(Ls))
            L = Ls[l]

            xs = data["xs"][l]
            ts = (xs-xc)/xc # reduced temperatures 
            ys = data[plot][l] # stuff to plot at this system size 

            if plot == "binds": 
                line0, = ax.plot(ts * Ls[l]**(1/nu0),ys,c=(0,0,0,0),label=r'$%d$'%Ls[l],marker='o',mfc=col,mew=mew,ms=ms,lw=lw)
            elif plot == "chis": 
                line0, = ax.plot(ts * Ls[l]**(1/nu0),ys * L**(-gamma0 / nu0),c=(0,0,0,0),label=r'$%d$'%Ls[l],marker='o',mfc=col,mew=mew,ms=ms,lw=lw)
            elif plot == "mags": 
                line0, = ax.plot(ts * Ls[l]**(1/nu0),ys * L**(beta0 / nu0),c=(0,0,0,0),label=r'$%d$'%Ls[l],marker='o',mfc=col,mew=mew,ms=ms,lw=lw)

            lines0.append(line0)

        def update(nu,beta,gamma):
            print("nu       = ",nu)
            print("gamma    = ",gamma)
            print("beta     = ",beta)
            print("beta(hs) = ",(d*nu-gamma)/2)
            print("eta (hs) = ",2-gamma/nu)
            print("Delta_ep = ",d-1/nu) # dimension of lightest Z2 neutral field
            print(" ")

            # will need to dynamically update the plot ranges 
            xmin = 1e10 
            xmax = -1e10
            ymin = 1e10 
            ymax = -1e10

            for l in range(len(Ls)): 
                L = Ls[l]
                xdat = (data["xs"][l]-xc)/xc * L**(1/nu)

                if plot == "chis":
                    ydat = data["chis"][l] * L**(-gamma/nu)
                elif plot == "mags": 
                    ydat = data["mags"][l] * L**(beta/nu)
                elif plot == "binds": 
                    ydat = data["binds"][l] 
                xmin = min(xmin,np.min(xdat))
                xmax = max(xmax,np.max(xdat))
                ymin = min(ymin,np.min(ydat))
                ymax = max(ymax,np.max(ydat))

                # update the line data 
                lines0[l].set_xdata(xdat) 
                lines0[l].set_ydata(ydat)
            
            ax.set_xlim(xmin*1.05,xmax*1.05)
            ax.set_ylim(ymin*.95,ymax*1.05)
            fig.canvas.draw_idle()

        # define a function that updates values of critical exponents based on arrow key presses 
        def on_key(event):
            global nu,beta,gamma
            d = .01 # incriment by which exponents are tuned 
            numin = 0; numax = 4
            gammamin = 0; gammamax = 4
            betamin = 0; betamax = 4  
            if event.key == 'left':
                nu = max(numin,nu - d) 
            elif event.key == 'right':
                nu = min(numax,nu + d) 
            elif event.key == 'up':
                gamma = min(gammamax,gamma+d)
            elif event.key == 'down':
                gamma = max(gammamin,gamma-d)
            elif event.key == '.': 
                beta = max(betamin,beta-d)
            elif event.key == ',': 
                beta = min(betamax,beta+d)

            update(nu,beta,gamma)
            
        fig.canvas.mpl_connect('key_press_event', lambda event : on_key(event))

        ax.legend(title=r'$L$')
        if plot == "binds": 
            ylab = r'$B$'
        elif plot == "chis": 
            ylab = r'$L^{-\gamma/\nu}\chi$'
        elif plot == "mags": 
            ylab = r'$L^{\beta/\nu} m$'
        ax.set(xlabel=r'$tL^{1/\nu}$',ylabel = ylab)
        ax.set_title(title)
    plt.show()
                



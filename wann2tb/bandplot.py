import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import lineStyles
red=(241.0/255,64.0/255,64.0/255)
blue=(26.0/255,111.0/255,223.0/255)
orange = (247.0/255,144.0/255,61.0/255)
green = (89.0/255,169.0/255,90.0/255)
purple = (103.0/255,80.0/255,131.0/255)
font = {'family' : 'Times New Roman',  
        'weight' : 'light',  
        'size'   : 25,  
        }
def band_plot(bandlist,ed=-2,eu=2,legend=True,filename=""):
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
    })
    for band in bandlist:
        color='k'
        label=''
        de=0
        fatband=[]
        bandid=[]
        if len(band)>1:
            color=band[1]
        if len(band)>2:
            label=band[2]
        if len(band)>3:
            fatband=band[3]
        if len(band)>4:
            bandid=band[4]
        if len(band)>5:
            de=band[5]
        kline=band[0][0]
        solution=band[0][1]
        kticks=band[0][2]
        klabels=band[0][3]
        wl=band[0][4]
        if fatband==[]:
            if bandid==[]:
                bandid=range(len(solution))
            lbool=True
            for i in bandid:
                if lbool:
                    plt.plot(kline,np.real(solution[i])+de, color=color, linewidth=2.0,label=label,zorder=2)
                    lbool=False
                else:
                    plt.plot(kline,np.real(solution[i])+de, color=color, linewidth=2.0,zorder=2)
        else:
            cl=[]
            if bandid==[]:
                bandid=range(len(solution))
            for i in bandid:
                cl.append([])
                plt.plot(kline,np.real(solution[i])+de, color='gray', linewidth=1.0, alpha=1.0,zorder=1)
            for i in bandid:
                for v in wl[i]:
                    cl[i].append(abs(np.matmul(np.matmul(v.transpose(),np.diag(fatband)),np.conj(v)))/abs(np.dot(np.conj(v.transpose()),v)))
            sc_k=[]
            sc_e=[]
            sc_c=[]
            for i in bandid:
                sc_k+=kline
                etemp=np.real(solution[i])+de
                sc_e+=etemp.tolist()
                sc_c+=cl[i]
            plt.scatter(sc_k,sc_e,s=(2*np.array(sc_c))**3,marker='o',color=color,zorder=2)
            plt.scatter(0,-100,s=100.0,marker='o',color=color,label=label,zorder=2)
            #plt.colorbar()
            plt.clim(0,1.0)
    plt.xticks(kticks,klabels)
    plt.tick_params(labelsize=25)
    plt.axhline(y=0, xmin=0, xmax=1,linestyle= '--',linewidth=0.5,color='0.5')
    for i in kticks[1:-1]:
        plt.axvline(x=i, ymin=0, ymax=1,linestyle= '--',linewidth=0.5,color='0.5')
    plt.ylim(ed,eu)
    plt.xlim(kticks[0],kticks[-1])
    if legend:
        plt.legend(loc='upper left', prop=font)
    plt.ylabel('E (eV)', fontdict=font)
    plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
    fig = plt.gcf()
    fig.set_size_inches(8,6)
    if filename != None:
        plt.savefig(filename,dpi=500)
    plt.show()
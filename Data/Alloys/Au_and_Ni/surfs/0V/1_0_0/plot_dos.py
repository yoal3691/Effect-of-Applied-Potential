# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 12:38:54 2021

@author: Yousef Alsunni

Onriginal Created on Fri May 15 13:42:30 2020

@author: zaba1157
"""

import json
from scipy.ndimage.filters import gaussian_filter1d
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import matplotlib
import numpy as np
from scipy.signal import argrelextrema
from scipy.signal import find_peaks



def get_data(filename):
    with open(filename,'r') as f:
        data = json.load(f)
    return data

def get_site_plot_dict(site_data,site_plot_args):
    site_plot_dict = site_plot_args
    for label in site_plot_args:
        if label == 'Total':
            for spin in ['spinup','spindn']:
                relEs = [x for x in site_data[spin]['E']] 
                sumdos = site_data[spin]['dos']['Total']
                site_plot_dict[label][spin] = {'relE':relEs,'sumdos':sumdos}
        else:       
            sites = [str(x) for x in site_plot_args[label]['sites']]
            els = [x for x in site_plot_args[label]['orbitals'].keys()]
            
            for spin in ['spinup','spindn']:
                relEs = [x for x in site_data[spin]['E']]
                sumdos = []
                for site in site_data[spin]['dos'].keys():
                    if site in sites:
                        el = list(site_data[spin]['dos'][site].keys())[0]
                        if el in els:
                            for orb in site_plot_args[label]['orbitals'][el]:
                                dos = site_data[spin]['dos'][site][el][orb]
                                if len(sumdos) == 0:
                                    sumdos = dos
                                else:
                                    sumdos = [x+dos[i] for i,x in enumerate(sumdos)]
                                    
                site_plot_dict[label][spin] = {'relE':relEs,'sumdos':sumdos}
            
    return site_plot_dict
                                          
def plot_site_dos(site_plot_dict,dos_dict,sigma,xlims=[-10,6]):
    colors =['rebeccapurple','red','cyan','silver','yellow','olive','magenta','royalblue', 'hotpink','lawngreen','darkcyan','saddlebrown']
    linestyles = ['solid','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid','solid']
    labels = []
    for index, label in enumerate(site_plot_dict.keys()):
        labels.append(label)
        for spin in ['spinup','spindn']:
            relEs = site_plot_dict[label][spin]['relE']
            sumdos = site_plot_dict[label][spin]['sumdos']
            
            if spin == 'spinup':                
                scale = max([x for i,x in enumerate(sumdos)
                            if relEs[i] < max(xlims)
                            and relEs[i] > min(xlims)])
#                scale = 1
                scaled_dos = [x/scale for x in sumdos]
            elif spin == 'spindn':
                scale = min([x for i,x in enumerate(sumdos)
                            if relEs[i] < max(xlims)
                            and relEs[i] > min(xlims)])
#                scale = -1
                
                scaled_dos = [-x/scale for x in sumdos]
            
            diff = [relEs[i + 1] - relEs[i]
                    for i in range(len(relEs) - 1)]
            avgdiff = sum(diff) / len(diff)
            smeared_dos = gaussian_filter1d(scaled_dos, sigma / avgdiff)
            
            #Do the plotting
            nE = [relEs[i] for i in range(len(relEs)) if relEs[i]<=0]
            nEs = [0]*len(relEs)
            nEs[:len(nE)] = nE
            ax.plot(relEs, smeared_dos,colors[index],alpha=1,linewidth=0.6)
            ax.fill_between(nEs, 0, smeared_dos, color=colors[index], alpha=0.1)
            
        
    
    lines = [Line2D([0], [0], color=c, linewidth=2, linestyle=linestyles[ii]) for ii, c in enumerate(colors)]
    plt.legend(lines, labels,loc='upper right', bbox_to_anchor=(1.4, 0.9))
    plt.ylim(-1.2, 1.2)
    plt.xlim(min(xlims), max(xlims))
    plt.yticks([])
    plt.axhline(y=0,color='k',linestyle='solid',linewidth=1)
    plt.axvline(x=0,color='k',linestyle='solid',linewidth=1)
    plt.setp(ax.spines.values(), linewidth=2)
    plt.xlabel('${E-{E_{fermi}}}$ (eV)')
    plt.ylabel('Density of States (a.u.)')
    plt.tight_layout()
    plt.show()
#    sumdosup = site_plot_dict[label]['spinup']['sumdos']    
#    peaks, _ = find_peaks(sumdosup,distance=50)
#    print("Peak Positions:")
#    for i in peaks:
#        print(relEs[i])

if __name__ == '__main__':
    
    pre = 'DOS'
    
    
    site_dos_file = pre+'_all_sites.json'
    summed_dos_file = pre+'.json'
    sigma = 0.01
    xlimits = [-6,3]
        
    # Write center data?
    write_file = True
    moments_write_file = pre+'_site_moments_dos.json'
    
    #plotting args
    figsize = (6,5)
    dpi = 200
    
    site_data = get_data(site_dos_file)
    summed_dos_data = get_data(summed_dos_file)
    
    matplotlib.rc('font',size=18)
    matplotlib.rc('axes',titlesize = 16)
    matplotlib.rc('figure',titlesize = 14)
    matplotlib.rc('xtick',labelsize=14)
    matplotlib.rc('ytick',labelsize=14)
    matplotlib.rc('legend',fontsize = 12)
    #rcParams['axes.linewidth'] = 1

    
    
        #plot site projected dos         
    """
    Example: COOH on Fe2P surface
    Fe: #4
    O*: #79
    Other O: #78
    
    COOH sites: [77,78,79,80]
    
    Example:
            
    ads_sites = [77,78,79,80]
        
        surf_sites = [x for x in site_data['spinup']['dos'].keys() 
        if x not in ads_sites and x != 'Total' ]
        
       # site plot args only plots what is provided in the key value pairs
        
       site_plot_args = {'P surf':{'sites':surf_sites,
                                'orbitals':{'P':['s','p']}},
            
                    'Fe surf':{'sites':surf_sites,
                                 'orbitals':{'Fe':['s','p','d']}},
                                 
                    'COOH*':{'sites':ads_sites,
                           'orbitals':{'O':['s','p'],'C':['s','p'],'H':['s']}},
                           
                    'Fe-d ads':{'sites':[4],
                           'orbitals':{'Fe':['d']}},
                                
                    'O*-p ads':{'sites':[79],
                           'orbitals':{'O':['p']}},                        
                    }        
        
        """
        
    
    ads_sites = []        
    surf_sites = [x for x in site_data['spinup']['dos'].keys() 
    if x not in ads_sites and x != 'Total' ]
    all_sites = [x for x in site_data['spinup']['dos'].keys() 
    if  x != 'Total' ]
    # site plot args only plots what is provided in the key value pairs        
    site_plot_args = {
#                     '$C_{(s)}$':{'sites':[65],
#                             'orbitals':{'C':['s']}},      
#
#                     '$Ag_{(total)}$':{'sites':[1],
#                             'orbitals':{'Ag':['s','p','d']}},      

#                     '$Ag_{(s)}$':{'sites':[1],
#                             'orbitals':{'Ag':['s']}},      
##
#                     '$Ag_{(p)}$':{'sites':[1],
#                             'orbitals':{'Ag':['p']}},      
#
                     '$Ni_{(s)}$':{'sites':[5],
                             'orbitals':{'Ni':['s']}},      

#                     '$C_{(px)}$':{'sites':[65],
#                             'orbitals':{'C':['px']}},      
#
#                     '$O_{(px)}$':{'sites':[66,67],
#                             'orbitals':{'O':['px']}},      
#
#                     '$C_{(py)}$':{'sites':[65],
#                             'orbitals':{'C':['py']}},      
#
#                     '$O_{(py)}$':{'sites':[66,67],
#                             'orbitals':{'O':['py']}},      
#
#                     '$C_{(pz)}$':{'sites':[65],
#                             'orbitals':{'C':['pz']}},      
#
#                     '$O_{(pz)}$':{'sites':[67,67],
#                             'orbitals':{'O':['pz']}},      
#
#                    '$Ag_{(s)}$':{'sites':[19],
#                           'orbitals':{'Ag':['s']}},                                 
##
##                    '$Ag_{(p)}$':{'sites':[18],
##                           'orbitals':{'Ag':['p']}},                                 
##
#                    '$Ag_{(d)}$':{'sites':[19],
#                           'orbitals':{'Ag':['d']}},                                 
#
#                    '$Ag_{(dxy)}$':{'sites':[18],
#                           'orbitals':{'Ag':['dxy']}},                                 
#
#                    '$Ag_{(dyz)}$':{'sites':[18],
#                           'orbitals':{'Ag':['dyz']}},                                 
#
#                    '$Ag_{(dz2)}$':{'sites':[18],
#                           'orbitals':{'Ag':['dz2']}},                                 
#
#                    '$Ag_{(dxz)}$':{'sites':[18],
#                           'orbitals':{'Ag':['dxz']}},                                 
#
#                    '$Ag_{(dx2-y2)}$':{'sites':[18],
#                           'orbitals':{'Ag':['dx2-y2']}},                                 
#
                            }    
        
        
    site_plot_dict =  get_site_plot_dict(site_data,site_plot_args)
        
    fig = plt.figure(figsize=figsize,dpi=dpi)
    ax = plt.axes()        
    plot_site_dos(site_plot_dict,summed_dos_data,sigma,xlims=xlimits) 

        
    
y_values= site_plot_dict['$Ni_{(s)}$']['spinup']['sumdos']

x_values= site_plot_dict['$Ni_{(s)}$']['spinup']['relE']

new_x_values = np.array([x for x in x_values if x <= 0])

new_y_values = np.array(y_values[:len(new_x_values)])
#
##new_x_values = np.array(new_x_values[15840:])
##new_y_values = np.array(new_y_values[15840:])
#    
area_under_the_curve = 2*np.trapz(new_y_values, x=new_x_values)     

print('Area Under Curve = '+str(area_under_the_curve))    
    
    
    
    
    
    
    
    
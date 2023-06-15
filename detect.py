import numpy as np
import matplotlib.pyplot as plt
import argparse
from const import elnames, ionstages, molnames, ifacod
from normalize import spec_normalize
import pandas as pd
import os

# initialize arrays for the data
wl = []
wll = []   
wlm = []   
wlj = []   
fl = []
bb = []
iext = []
icnt = [] 
idl = []
iml = []  
idm = []  
idj = []  
imm = []  
imj = [] 

# temp fl
fl_temp= []
fl_norm = []

# maxpoints allowed in the file
maxpoint = 2000000
maxmol = 326

# storing the data in a dictionary key -> element / molecule, value -> tuple(list[x-cordinate for mapping the text], list[text string to be mapped])
ele_pos = {}

# get the data from the file and store it in the arrays
def read_file(file):
    global fl, bb, fl_temp, wl, wll, wlm, wlj, iext, icnt, idl, iml, idm, idj, imm, imj

    filename, file_extension = os.path.splitext(file)

    if file_extension == '.xz':
        data = pd.read_table(file, compression='xz', delim_whitespace=True, header=None)
    else:
        data = pd.read_table(file, delim_whitespace=True, header=None)

    if len(data) > maxpoint:
        print('Maximum number of points reached: ', maxpoint)
    else:
        wl = data.iloc[:, 0].astype(float).tolist()  # wldum
        wll = data.iloc[:, 5].astype(float).tolist()  # wlldum
        wlm = data.iloc[:, 8].astype(float).tolist()  # wlmdum
        wlj = data.iloc[:, 11].astype(float).tolist()  # wljdum
        fl = data.iloc[:, 1].astype(float).tolist()  # fldum
        bb = data.iloc[:, 2].astype(float).tolist()  # bbdum
        iext = data.iloc[:, 3].astype(float).tolist()  # iextdum
        icnt = data.iloc[:, 4].astype(float).tolist()  # icntdum
        idl = data.iloc[:, 6].astype(float).tolist()  # idldum
        iml = data.iloc[:, 7].astype(float).tolist()  # imldum
        idm = data.iloc[:, 9].astype(float).tolist()  # idmdum
        idj = data.iloc[:, 12].astype(float).tolist()  # idjdum
        imm = data.iloc[:, 10].astype(float).tolist()  # immdum
        imj = data.iloc[:, 13].astype(float).tolist()  # imjdum

    print('Done reading file: ', file)
    print('Number of lines read: ', len(wl))

    # manipulation in fl and bb
    fl = 10**(np.array(fl) - 40)
    bb = 10**(np.array(bb) - 40)

    fl_temp = fl.copy()


def plot(ax, xmin, xmax, c='black', lw = '1'):
    global wl
    wl = np.array(wl)
    mask = np.logical_and(wl>=xmin, wl<=xmax)
    ax.plot(wl[mask], fl[mask], c=c, lw=lw)


def plotid(xmin,xmax,flag,molflg,delatom,delmol,excode=0,exlin=0):
    
    nlines = len(wl)
    wllo = 0.0
    idlo = 0
    wlmo = 0.0
    idmo = 0

    def helper_plotid(i):
        nonlocal wllo, idlo, wlmo, idmo
    
        if molflg == 0: return
        el = int(idm[i])
        if el == exlin: return
        if ifacod[el] == excode: return
        if imm[i] > iext[i] + delmol or imm[i] > iml[i] or wlm[i] < xmin or wlm[i] > xmax: return
        if wlmo == wlm[i] and idmo == idm[i]: return
        if flag == 1:
            if el <= maxmol:
                try:
                    if min([el,maxmol,int(ifacod[el])+1]) <= 0: ###ifacod was present instead of ifacod[el]
                        print("bad IDM", wlm[i], wlmo, el, idm[1], ifacod)
                except:
                    if min([el,maxmol]) <= 0: ### this cond. is not present in initial code
                        print("bad IDM", wlm[i], wlmo, el, idm[1])
                    
                sid = molnames[el].strip()
                swl = format(wlm[i],'.1f')
                labstr = sid + ' ' + swl
            else:
                sid = str(int(idm[i]))
                swl = format(wlm[i],'.1f')
                labstr = sid + ' ' + swl
            
        if flag == 0:
            sid = str(el)
            swl = format(wlm[i],'.1f')
            labstr = sid + ' ' + swl
        
        if flag == 2:
            sid = str(ifacod[el])
            swl = format(wlm[i],'.1f')
            labstr = sid + ' ' + swl
        
        if sid not in ele_pos:
            ele_pos[sid] = ([wlm[i]], [labstr])
        else:
            ele_pos[sid][0].append(wlm[i])
            ele_pos[sid][1].append(labstr)

        wlmo = wlm[i]
        idmo = idm[i]



    for i in range(nlines):
        if iml[i] > icnt[i] + delatom or iml[i] > imm[i] or wll[i] < xmin or wll[i] > xmax:
            helper_plotid(i)
            continue
        if wllo == wll[i] and idlo == idl[i]:
            helper_plotid(i)
            continue
        
        if flag == 1:
            el = int(idl[i]/100)
            ion = int(idl[i]) - el*100
            if el <= 112:
                sid = elnames[min([el,112])]
                sion = ionstages[ion]
                swl = format(wll[i],'.1f')
                labstr = sid + "[" + sion + "]" + " " + swl
            else:
                sid = str(int(idl[i]))
                swl = format(wll[i],'.1f')
                labstr = sid + ' ' + swl

        if flag == 0 or flag == 2:
            sid = str(int(idl[i]))
            swl = format(wll[i],'.1f') 
            labstr = sid + ' ' + swl

        if sid not in ele_pos:
            ele_pos[sid] = ([wll[i]], [labstr])
        else:
            ele_pos[sid][0].append(wll[i])
            ele_pos[sid][1].append(labstr)

        wllo = wll[i]
        idlo = idl[i]


last_file = ''
def driver_gui(file,color,norm,ele,xmin,xmax,flag,molflg,delatom,delmol,excode=0,exlin=0):
    global ele_pos, fl, fl_norm,last_file

    if len(wl) == 0 or file != last_file:
        read_file(file)
        last_file = file
        fl_norm = spec_normalize(fl)

    # normalize fl
    if norm:
        print("converting to normalized data")
        fl = fl_norm.copy()
    else:
        fl = fl_temp.copy()
    

    fig, ax = plt.subplots()
    plot(ax, xmin, xmax, color, 0.5)
    pre_ymax = ax.get_ylim()[1]

    plotid(xmin,xmax,flag,molflg,delatom,delmol,excode,exlin)

    ystart = 0
    ytext = pre_ymax   

    def draw_line():
        for i in range(len(ele_pos[e][0])):
            ax.plot([ele_pos[e][0][i], ele_pos[e][0][i]], [ystart, ytext], color='blue', linestyle= (0, (5, 10)), linewidth=0.5)
            ax.text(ele_pos[e][0][i], ytext, " " + ele_pos[e][1][i], rotation='vertical', horizontalalignment='center', verticalalignment='bottom', color='blue', fontsize = 'xx-small') 
   
    
    if type(ele) == list:
        
        if ele[0].strip() == 'all-atoms':
            for e in ele_pos.keys():
                if flag == 1:
                    draw_line()
                if flag == 0:
                    draw_line()

        elif ele[0].strip() == 'all-molecules':
            for e in ele_pos.keys():
                if flag == 1:
                    if e in molnames:
                        draw_line()
                if flag == 0:
                    if int(e) >= 0 and int(e) <= maxmol:
                        draw_line()
        else:
            for e in ele:
                if e in ele_pos.keys():
                    draw_line()
                else:
                    print("element " + e + " not found")
    else:
        print("ele must be a string or a list of strings")

    ele_pos = {}
    return fig, ax


def driver(file,ele,xmin,xmax,flag,molflg,delatom,delmol,excode=0,exlin=0):

    read_file(file)

    fig, ax = plt.subplots(figsize=(10,5))
    plot(ax, xmin, xmax)
    pre_ymax = ax.get_ylim()[1]
    plt.ylim(ax.get_ylim()[0]*0.5, pre_ymax*1.7)
    plotid(xmin,xmax,flag,molflg,delatom,delmol,excode,exlin)

    ystart = 0
    ytext = pre_ymax + 0.05*pre_ymax        

    if type(ele) == list:
        
        if ele[0].lower() == 'all':
            for e in ele_pos.keys():
                for i in range(len(ele_pos[e][0])):
                    ax.plot([ele_pos[e][0][i], ele_pos[e][0][i]], [ystart, ytext], 'r', linestyle= (0, (5, 10)), linewidth=0.5)
                    ax.text(ele_pos[e][0][i], ytext, " " + ele_pos[e][1][i], rotation='vertical', horizontalalignment='center', verticalalignment='bottom', color='r') 
        else:
            for e in ele:
                if e in ele_pos.keys():
                    for i in range(len(ele_pos[e][0])):
                        ax.plot([ele_pos[e][0][i], ele_pos[e][0][i]], [ystart, ytext], 'r', linestyle= (0, (5, 10)), linewidth=0.5)
                        ax.text(ele_pos[e][0][i], ytext, " " + ele_pos[e][1][i], rotation='vertical', horizontalalignment='center', verticalalignment='bottom', color='r')
                else:
                    print("element " + e + " not found")
    else:
        print("ele must be a string or a list of strings")

    return fig, ax


if __name__ == '__main__':
   
    # parse command line arguments
    parser = argparse.ArgumentParser(description='plot spectrum with identified lines')
    parser.add_argument('file', type=str, help='data file')
    parser.add_argument('ele', type=str, nargs='+', help='element or list of elements to plot')
    parser.add_argument('xmin', type=float, help='minimum wavelength to plot')
    parser.add_argument('xmax', type=float, help='maximum wavelength to plot')
    parser.add_argument('flag', type=int, help='0: print numeric identification, 1: print text identification, 2: print text form element number for atomic lines and ifacod for moleciular lines')
    parser.add_argument('molflg', type=int, help='0: do not plot molecular lines, 1: plot molecular lines')
    parser.add_argument('delatom', type=float, help='parameter for atom print threshold, gives difference between iml and icnt for lines to plot')
    parser.add_argument('delmol', type=float, help='parameter for molecular print threshold, gives difference between imm and iext for lines to plot')
    parser.add_argument('--excode', type=int, default=0, help='molecular codes (FA) to exclude from plotting ')
    parser.add_argument('--exlin', type=int, default=0, help='molecular line sets (uselin codes) to exclude from plotting ')
    parser.add_argument('-p', action='store_true', help='print all of the elements in the spectrum')

    args = parser.parse_args()
    driver(args.file,args.ele, args.xmin, args.xmax, args.flag, args.molflg, args.delatom, args.delmol, args.excode, args.exlin)
    
    if args.p == True:
        print(*list(ele_pos.keys()))

    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('Flux (erg/s/cm^2/Angstrom)')
    plt.title('Spectrum')

    plt.show()
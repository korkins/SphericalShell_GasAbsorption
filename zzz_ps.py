# -*- coding: utf-8 -*-
#
import numpy as np
import matplotlib.pyplot as plt

nvz = 18
nsz = 10
naz = 5
color = 'rgbkm'
az_legend = ['Az=0', 'Az=45', 'Az=90', 'Az=135', 'Az=180']

nlines = nvz*nsz*naz
ncols = 4

fname = "RERR_SS_PS.txt"

fdat = np.loadtxt(fname)
dat = np.zeros((nlines, ncols))
dat = fdat

x = np.zeros(nvz)
x[0:nvz] = dat[0:nvz, 2]
print(x)

plt.figure()
for isz in range(0, nsz):
    sza = dat[isz*naz*nvz, 0]
    plt.subplot(2,5,isz+1)
    y = np.zeros(nvz)
    icolor = -1
    for iaz in range(0, naz):
        icolor = icolor+1
        for ivz in range(0, nvz):
            ix = isz*naz*nvz + iaz*nvz + ivz
            y[ivz] = dat[ix, 3]
        plt.plot(x, y, color[icolor])
    plt.grid(True)
    plt.title('sza='+str(sza))
    if isz == 0:
        plt.ylabel('Relative Error, %')
    if isz == 5:
        plt.ylabel('Relative Error, %')
    if isz > 4:
        plt.xlabel('VZA, degrees')
    if isz == 9:
        plt.legend(az_legend, loc = 0)
#plt.tight_layout()
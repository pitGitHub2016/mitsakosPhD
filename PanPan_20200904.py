# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 10:14:14 2020

@author: DTryfonopoulos
"""
#%% Load The data

import get_MRD, numpy as np
filename = '2D_Rad_FINAL_1_TE5.MRD'

data = get_MRD.get_mrd_3d(filename) 
dim = data[1]
para = data[2]
data = np.squeeze(data[0])

#%% Calculatet he Shift
nrows = np.array(np.size(data,0))
ncols = np.array(np.size(data,1))

target =16384; 
newcentre = target/2;
pad2 = (target-ncols)/2
pad=np.zeros((1, np.int(pad2)));
pad = np.squeeze(pad)
fout = data*0;
xshifts = np.zeros(nrows,)

for i in range(0,nrows):
    d0 = data[i,:];
    d1 = np.fft.fftshift(np.fft.fft(d0));
    d1 =np.concatenate((pad, d1, pad))
    d1 = np.fft.ifft(d1);
    
    Y= np.max(np.abs(d1));
    xcentre = np.argmax(np.abs(d1))
    
    xshift= (newcentre -xcentre)/target * ncols; 
    xshift = xshift #Shift value in pixels 
    
    H = np.fft.fftshift(np.fft.fft(d0));
    
    xF = np.arange(-ncols/2,ncols/2,1)
    
    #H = H.* np.e(-1j * 2 * np.pi.*(xF * xshift/ncols));
    H = np.multiply(H, np.exp((np.multiply((-1j * 2 * np.pi),(xF * xshift/ncols)))))
    
    d1 = np.fft.ifft(np.fft.fftshift(H));
    fout[i,:] = d1;
    
    xshifts[i] = xshift;


plt.plot(xshifts, '*'); plt.title('Shift per Spoke')



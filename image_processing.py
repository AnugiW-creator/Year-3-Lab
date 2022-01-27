# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 09:36:26 2022

@author: anugi
"""


from astropy.table import Table
from astropy.io import fits
import astropy.cosmology as cosmo
from astropy.units import quantity
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

hdulist = fits.open(r'C:\Users\anugi\OneDrive\Documents\Physics\Year 3\Labs\Astronomical Image Processing\A1_mosaic\A1_mosaic.fits')
#%%
pixel = hdulist[0].data
print(pixel[0])

columns = len(pixel[0])
rows = len(pixel)
all_pixels = []
for i in range(rows):
    pix_arr = pixel[i]
    for j in range(columns):
        all_pixels.append(pix_arr[j])
print(all_pixels)

plt.hist(all_pixels, bins=900, color='deeppink')
plt.show()
#%%
print('max pixel value:', max(all_pixels))
print('min pixel value:', min(all_pixels))

#%%
def Gaussian(x, sigma, mu, norm):
    gaus = norm*(1/(sigma*(sp.sqrt(2*sp.pi))))*(sp.exp(-(1/2)*((x-mu)/sigma)**2))
    return gaus
#%%
anarray = sp.linspace(3350,3550,num=400,endpoint=True)

newarray = Gaussian(anarray, 50,3450,200*1000000)

plt.plot(anarray, newarray)
plt.show()
#%%
########### IDENTIFYING BACKGROUND ############
background = []
for m in range(len(all_pixels)):
    value = all_pixels[m]
    if value>=3300 and value<=3600:
        background.append(value)
        
values,bins,patches= plt.hist(background, bins=300, color='deeppink')

#Test parameters
test_sigma = 50 #Standard deviation
test_mu = 3450 #Mean
test_norm = 100*1000000 #Normalisation
p0 = [test_sigma, test_mu, test_norm]

fit, fit_cov = curve_fit(Gaussian, xdata = bins[0:300], ydata = values, p0=p0, maxfev=5000)
print(fit)
data_fit = Gaussian(bins[0:300], *fit)
plt.plot(bins[0:300], data_fit, color = 'blueviolet')

print(" sigma = " + str(fit[0])+ '+/-' +str (np.sqrt(fit_cov[0,0])))
print(" mu = " + str(fit[1])+ '+/-' +str (np.sqrt(fit_cov[1,1])))
print(" norm = " + str(fit[2])+ '+/-' +str (np.sqrt(fit_cov[2,2])))
plt.show()
#%%
twodarray = np.array(([2,3,4,5], [3,4,5,6], [4,5,6,7])) # Create test 2d array
maskarray = np.zeros((3,4)) # Zeros array with same dimensions as 2d array

for m in range(3):
    for n in range(4):
        element = twodarray[m][n]
        bloop = maskarray[m][n]
        if element <= 4: # Pick arbitrary value to exclude below
            #bloop = 1
            maskarray[m][n] = 1 # Change corresponding mask element
print(maskarray)
#%%
mask = np.ones((4611,2570)) # Create mask of 1s to indicate valid objects and background

for i in range(4611):
    for j in range(2570): # For loop through whole image
        # girl what the fuck
        impix = pixel[i][j] # Assign variable to pixel
        if impix <= 3454: # Exclude values below chosen background
            mask[i][j] = 0 # Change 1 to 0
print(mask)

#%%
zoom = []
for p in range(len(all_pixels)):
    if all_pixels[p] >=4000: #and pix[p] <= 4000:
        zoom.append(all_pixels[p])

plt.hist(zoom, bins = 100)
#%%
mask_ones = np.ones((4611,2570))
for i in range (len(pixel)):
    for j in range (len(pixel[0])):
        if pixel[i,j] <= 3455:
            mask_ones[i,j] = 0
#%%
################ IDENTIFYING STARS/NON-GALAXIES ################
objects = []
for m in range(len(all_pixels)):
    value = all_pixels[m]
    if value>=3600:# and value<=3600:
        objects.append(value)
        
values,bins,patches= plt.hist(objects, bins=10000, color='deeppink')

# =============================================================================
# #Test parameters
# test_sigma = 50 #Standard deviation
# test_mu = 3450 #Mean
# test_norm = 100*1000000 #Normalisation
# p0 = [test_sigma, test_mu, test_norm]
# 
# fit, fit_cov = curve_fit(Gaussian, xdata = bins[0:300], ydata = values, p0=p0, maxfev=5000)
# print(fit)
# data_fit = Gaussian(bins[0:300], *fit)
# plt.plot(bins[0:300], data_fit, color = 'blueviolet')
# 
# print(" sigma = " + str(fit[0])+ '+/-' +str (np.sqrt(fit_cov[0,0])))
# print(" mu = " + str(fit[1])+ '+/-' +str (np.sqrt(fit_cov[1,1])))
# print(" norm = " + str(fit[2])+ '+/-' +str (np.sqrt(fit_cov[2,2])))
# =============================================================================

plt.show()
#%%
string = 'box(7.9,6.8)'
splitstring = string.split()
print(splitstring)
#%%
############# LOAD REGIONS TO BE MASKED FOR LEFT SIDE OF IMAGE ################
x, y, width, height, z = sp.loadtxt('C:/Users/anugi/OneDrive/Documents/Physics/Year 3/Labs/Astronomical Image Processing/left_regions.txt', dtype=str, unpack=True)
print(x)
print(width)
print(height)

########### DO LEFT SIDE ##############
mask = np.ones((4611,2570)) # Create mask of 1s to indicate valid objects and background

for k in range(len(x)):
    xcentrefloat = float(x[k])
    ycentrefloat = float(y[k])
    widthfloat = float(width[k])
    heightfloat = float(height[k])
    #widthfloat = width[k]
    #heightfloat = height[k]
    xcentre = int(xcentrefloat) - 1
    ycentre = int(ycentrefloat) - 1
    width1 = (int(widthfloat) + 1)/2
    height1 = (int(heightfloat) + 1)/2    
    widthedge1 = xcentre - width1
    widthedge2 = xcentre + width1
    heightedge1 = ycentre - height1
    heightedge2 = ycentre + height1
    for i in range (4611):
        for j in range (2570):
            if i >= heightedge1 and i <=heightedge2:
                if j >= widthedge1 and j <= widthedge2:
                    mask[i][j] = 0
print(mask)
#yote = mask[0]

####################### LOAD REGIONS FOR RIGHT SIDE ####################
x, y, width, height, z = sp.loadtxt('C:/Users/anugi/OneDrive/Documents/Physics/Year 3/Labs/Astronomical Image Processing/mask_1.txt', dtype=str, unpack=True)
print(x)
print(width)
print(height)

############## DO RIGHT SIDE ###############
for k in range(len(x)):
    xcentrefloat = float(x[k])
    ycentrefloat = float(y[k])
    widthfloat = float(width[k])
    heightfloat = float(height[k])
    #widthfloat = width[k]
    #heightfloat = height[k]
    xcentre = int(xcentrefloat) - 1
    ycentre = int(ycentrefloat) - 1
    width1 = (int(widthfloat) + 1)/2
    height1 = (int(heightfloat) + 1)/2    
    widthedge1 = xcentre - width1
    widthedge2 = xcentre + width1
    heightedge1 = ycentre - height1
    heightedge2 = ycentre + height1
    for i in range (4611):
        for j in range (2570):
            if i >= heightedge1 and i <=heightedge2:
                if j >= widthedge1 and j <= widthedge2:
                    mask[i][j] = 0
print(mask)

############# LOAD EDGES #######################
x, y, width, height, z = sp.loadtxt('C:/Users/anugi/OneDrive/Documents/Physics/Year 3/Labs/Astronomical Image Processing/mask_3.txt', dtype=str, unpack=True)
print(x)
print(width)
print(height)

############## EDGES ###############
for k in range(len(x)):
    xcentrefloat = float(x[k])
    ycentrefloat = float(y[k])
    widthfloat = float(width[k])
    heightfloat = float(height[k])
    #widthfloat = width[k]
    #heightfloat = height[k]
    xcentre = int(xcentrefloat) - 1
    ycentre = int(ycentrefloat) - 1
    width1 = (int(widthfloat) + 1)/2
    height1 = (int(heightfloat) + 1)/2    
    widthedge1 = xcentre - width1
    widthedge2 = xcentre + width1
    heightedge1 = ycentre - height1
    heightedge2 = ycentre + height1
    for i in range (4611):
        for j in range (2570):
            if i >= heightedge1 and i <=heightedge2:
                if j >= widthedge1 and j <= widthedge2:
                    mask[i][j] = 0
print(mask)
#%%
############### REPLOT GAUSSIAN FOR BACKGROUND ####################
background = []
for m in range(4611):
    for n in range(2570):
        impix = pixel[m][n]
        maskpos = mask[m][n]
        if maskpos == 1:
            #background.append(impix)
            if impix>=3300 and impix<=3600:
                background.append(impix)

def Gaussian(x, sigma, mu, norm):
    gaus = norm*(1/(sigma*(sp.sqrt(2*sp.pi))))*(sp.exp(-(1/2)*((x-mu)/sigma)**2))
    return gaus

bins=np.linspace(3299.5, 3600.5, num=301, endpoint=False)
print(bins)

values,b,patches= plt.hist(background, bins=bins, color='deeppink')
print(b)
#Test parameters
test_sigma = 50 #Standard deviation
test_mu = 3450 #Mean
test_norm = 100*1000000 #Normalisation
p0 = [test_sigma, test_mu, test_norm]

fit, fit_cov = curve_fit(Gaussian, xdata = b[0:300], ydata = values, p0=p0, maxfev=5000)
print(fit)
data_fit = Gaussian(b[0:300], *fit)
plt.plot(b[0:300], data_fit, color = 'blueviolet')

print(" sigma = " + str(fit[0])+ '+/-' +str (np.sqrt(fit_cov[0,0])))
print(" mu = " + str(fit[1])+ '+/-' +str (np.sqrt(fit_cov[1,1])))
print(" norm = " + str(fit[2])+ '+/-' +str (np.sqrt(fit_cov[2,2])))
plt.show()
#%%
sigma = fit[0]
mu = fit[1]

cut_background = mu + 4.5 * sigma
print(cut_background)
#%%

###### Choose to cut off background from 3500 #######
################## CUT OFF BACKGROUND #######################

for i in range(4611):
    for j in range(2570): # For loop through whole image
        # girl what the fuck
        impix = pixel[i][j] # Assign variable to pixel
        if impix <= 3500: # Exclude values below chosen background
            mask[i][j] = 0 # Change 1 to 0
print(mask)

# =============================================================================
# pot_sources = [] # lmao
# pot_position = []
# for i in range(4611):
#     for j in range(2570): # For loop through whole image
#         # girl what the fucK
#         impix = pixel[i][j]
#         maskpix = mask[i][j]
#         if maskpix == 1:
#             pot_sources.append(impix)
# print(len(pot_sources))
# =============================================================================

######################### FINDING THE SOURCES #####################
pix_value = []
pix_pos = []
for i in range(4611):
    for j in range(2570):
        impix = pixel[i][j]
        maskpix = mask[i][j]
        if maskpix == 1:
            pix_value.append(impix)
            pix_pos.append(([i,j]))
copy_pixvalue = pix_value
copy_pixpos = pix_pos
copy_mask = mask
#%%
################### CHECK INDEX OF MAXIMUM PIXEL VALUE ################
# =============================================================================
# max_value = max(pix_value)
# max_index = pix_value.index(max_value) # Find index 
# sepindex = pix_pos[max_index]
# print(max_value)
# print(max_index)
# print(pix_pos[max_index])
# print(sepindex[0])
# =============================================================================
#%%

'''
NOTE: Indices are reversed:
    Takes form of [y, x]
'''
pix_value = copy_pixvalue
pix_pos = copy_pixpos
############# DEFINE CIRCLE AND ELLIPSE ############
def Circle(x, y, x0, y0):
    r = ((x - x0)**2 + (y - y0)**2)**(1/2)
    return r

def Ellipse(x, y, x0, y0, a, b):
    one = ((x - x0)**2)/(a**2) + ((y - y0)**2)/(b**2)
    return one

################## FIND SMALLEST SOURCES ###############
sources = []

for i in range(4611):
    for j in range(2570):
        impix = pixel[i][j]
        maskpix = mask[i][j]
        index = [i, j]
        if maskpix == 1:
            #print(index)
        #if j == 700 and i == 350:
            #break
            max_value = max(pix_value)
            max_index = pix_value.index(max_value) # Find index 
            sepindex = pix_pos[max_index] # make sure index of value in loop matches brightest pixel
            if sepindex == index:
                x0 = sepindex[1] # Separate x and y coordinates
                y0 = sepindex[0]
                r1 = 21 # Choose aperture radius from region file containing radius of brightest source
                r2 = 23 # Choose bg aperture size
                bg_aperture = []
                aperture = []
                bg_indices = []
                indices = []
                print('pixvalue:', len(pix_value))
                for k in range(len(pix_value)):
                    print(k)
                    if k == len(pix_value):
                        break
                    pv = pix_value[k]
                    pi = pix_pos[k]
                    x = pi[1]
                    y = pi[0]
                    radius = Circle(x, y, x0, y0)
                    #print(k)
                    print('pixvalue before:', len(pix_value), 'and pi:', pi)
                    if radius >= r1 and radius <= r2: # Check if pixel is within background aperture
                        bg_aperture.append(pv)
                        bg_indices.append(pi)
                        pix_value.pop(k) # Remove pixel and its indices from list
                        pix_pos.pop(k)
                        print('new pixvalue, bg:', len(pix_value))
                    elif radius <= r1: # Check if pixel is in main aperture (part of the source)
                        aperture.append(pv)
                        indices.append(pi)
                        pix_value.pop(k) # Remove pixel and its indices from list
                        pix_pos.pop(k)
                        print('new pixvalue:', len(pix_value), 'and pi:', pi)
                    
                #print('new pixvalue:', len(pix_value), 'and max value:', max_value)
                if len(bg_aperture) == 0:
                    fileindex = str(k)
                    #comb = np.vstack(aperture) # turn to 2d array
                    f = open('Sources/source'+fileindex, 'x') # Create file to put data in
                    np.savetxt(f, aperture)
                    f.close()
                    file = open('Index/index'+fileindex, 'x')
                    np.savetxt(file, indices)
                    for ind in range(len(indices)):
                        pixind = pix_pos[ind]
                        xpos = pixind[1]
                        ypos = pixind[0]
                        mask[ypos][xpos] = 0 # Change value in pixel position in mask to 0
                else:
                    continue
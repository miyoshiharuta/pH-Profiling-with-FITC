'''
Created on Feb 19, 2017

@author: Surf32
'''
import tifffile
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib._cm import _CMRmap_data
import matplotlib
from builtins import range
from bokeh.models.ranges import Range
matplotlib.rcParams.update({'font.size': 12})
import pandas as pd
import skimage.morphology
from  czifile import CziFile
import osDB
import pyqtgraph
import scipy
import skimage

pH55 = r"/media/daniel/Windows/Users/dnabu/Desktop/ResearchYogaWindows/Miyoshi/7List/phImages/pH5.5standard.tif"
pH60 = r"/media/daniel/Windows/Users/dnabu/Desktop/ResearchYogaWindows/Miyoshi/7List/phImages/pH6.0standard.tif"
pH65 = r"/media/daniel/Windows/Users/dnabu/Desktop/ResearchYogaWindows/Miyoshi/7List/phImages/pH6.5standard.tif"

#get image data
#create a pandas dataframe to hold data
index = ['pH55', 'pH60', 'pH65']
pHframe = pd.DataFrame(index=index)
pHframe['Files'] = ''
for row in index:
    pHframe['Files'].loc[row] = eval(row)

#get image data
pHframe['Image']=''
for row in pHframe.index:
    pHframe['Image'].loc[row] = tifffile.imread(pHframe['Files'].loc[row])
    
#shape?
print(pHframe['Image'].iloc[0].shape)

#get the mean inensity for each rgb channel
pHframe['Mean_intensity'] = ''
rgb = [0, 1, 2]
for row in pHframe.index:
    mean1 = []
    for c in rgb:
        mean1.append(np.mean(pHframe['Image'].loc[row][:, :, c]))
    pHframe['Mean_intensity'].loc[row] = mean1
    
#get green channel
pHframe['Mean_intensity'] = ''
rgb = [0, 1, 2]
for row in pHframe.index:
    pHframe['Mean_intensity'].loc[row] = np.mean(pHframe['Image'].loc[row][:, :, 1])

pHframe['pH'] = [5.5, 6.0, 6.5]
plt.plot(pHframe['pH'].values, pHframe['Mean_intensity'].values, '.')
par= np.polyfit(pHframe['pH'].values, pHframe['Mean_intensity'].values, 1, full=True)
slope = par[0][0]
intercept=par[0][1]
xl = [min(pHframe['pH'].values), max(pHframe['pH'].values)]
yl = [slope*xx + intercept  for xx in xl]

plt.plot(xl, yl)
savefile = r'/media/daniel/Windows/Users/dnabu/Desktop/ResearchYogaWindows/Miyoshi/7List/Fig2/RegressionPlot.svg'
plt.savefig(savefile)


#normal intensity values within image
def normalizeImage(img):    
    img = img - np.min(img)
    return img / np.max(img)

#determines the range
def setrange(x, xrange, limits):
    #x = number in index
    #range = the + and - of the given index
    #limits = [xin xmax] in possible range
    if x-xrange < limits[0]:
        x1 = limits[0]
    else:
        x1 = x-xrange
    if x + xrange > limits[1]:
        x2 = limits[1]
    else:
        x2 = x +xrange
    return np.arange(x1, x2)

#### now look at images
img11 = r"/media/daniel/Windows/Users/dnabu/Desktop/ResearchYogaWindows/Miyoshi/7List/rootImagesIsoCropped/11pHeffectFITC copy.tif"
img18 = r"/media/daniel/Windows/Users/dnabu/Desktop/ResearchYogaWindows/Miyoshi/7List/rootImagesIsoCropped/18pH5.510mMMES1mMKCl1mMCaCl2wroot copy.tif"

imgFrame = pd.DataFrame(index=['img11', 'img18'])
imgFrame["Files"] = ''
for row in imgFrame.index:
    imgFrame['Files'].loc[row] = eval(row)

#get image data
imgFrame['Image'] = ''
for row in imgFrame.index:
    imgFrame['Image'].loc[row] = tifffile.imread(imgFrame['Files'].loc[row])

#create mask of image
imgFrame['Mask']=''

for row in imgFrame.index:
    img = imgFrame['Image'].loc[row][:, :, 1]
    plt.imshow(img)
    reverse_index = np.where(img.flatten() <75)
    imgmask = np.zeros(img.flatten().shape, dtype=np.bool)
    imgmask[reverse_index] = 1
    imgmask = imgmask.reshape(img.shape)
    fig1 =plt.figure()
    ax1 = fig1.add_axes([0.1,0.1,0.8,0.8])
    ax1.imshow(imgmask)
    imgmask2 = scipy.ndimage.binary_fill_holes(imgmask)
    #remove single points
    struct = skimage.morphology.disk(20, dtype = np.bool)
    imgmask3 = skimage.morphology.binary_erosion(imgmask2, struct)
    plt.imshow(imgmask2)
    #dilate imgmask2
    #struct = scipy.ndimage.generate_binary_structure(5,5).astype(imgmask.dtype)
    struct = skimage.morphology.disk(50, dtype = np.bool)
    imgmask4 = skimage.morphology.binary_dilation(imgmask3, struct)
    indexTrue = np.where(imgmask2.flatten())
    indexTrueDilated = np.where(imgmask4.flatten())
    indexUnique = np.setdiff1d(indexTrueDilated, indexTrue)
    mask5 = np.zeros(img.flatten().shape, dtype = np.bool)
    mask5[indexUnique] = 1
    mask5 = mask5.reshape(img.shape)
    imgFrame['Mask'].loc[row] = mask5
    
import matplotlib.colors as colors
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap
cmap = plt.get_cmap('CMRmap')
new_cmap = truncate_colormap(cmap, 0, 0.8)
#creating graph for pH from masked region
xrange = 11
#get min and max limits for images
imgFrame['Limits']=''
for row in imgFrame.index:
    img = imgFrame['Image'].loc[row][:, :, 1]
    min1 = np.min(img)
    max1 = np.max(img)
    imgFrame['Limits'].loc[row] = [min1, max1]
limits = imgFrame['Limits'].values
limits = np.vstack(limits)
climits = [np.min(limits[:, 0]), np.max(limits[:,1])]

matplotlib.rcParams.update({'font.size': 12})
savefolder = r"/media/daniel/Windows/Users/dnabu/Desktop/ResearchYogaWindows/Miyoshi/7List/Fig2"
for row in imgFrame.index:
    fig1 = plt.figure( dpi=300)
    img = imgFrame['Image'].loc[row][:, :, 1].copy()
    limits = [0, img.shape[0]]
    maskIndex = imgFrame['Mask'].loc[row]
    mask1 = np.zeros(img.shape, dtype=np.bool)
    mask1[maskIndex] = 1
    ax11 = fig1.add_axes([0.1, 0.2, 0.4, 0.83])
    ax11.axis('off')
    imgax1 = ax11.imshow(img, clim = climits)
    imgax1.set_cmap(new_cmap)
    ticks1 = np.linspace(climits[0], climits[1], 4)
    cb = fig1.colorbar(imgax1, location = 'top', shrink =0.8, ticks = ticks1)#, shrink) #, orientation = 'horizental')
    ax11.set_visible('off')
    #replace tick marks on colorbar with pH
    ticks1 = cb.ax.get_xticks()
    ticks = cb.ax.get_xticklabels()
    pHticks =[]
    for i in ticks:
        pHticks.append("{0:.2f}".format((float(i.get_text()) - intercept)/slope))
    #cb.set_ticks(ticks1.tolist())
    cb.set_ticklabels(pHticks)
    cb.ax.tick_params(labelsize = 12)
    
    ax1 = fig1.add_axes([0.1, 0.1, 0.4, 0.8])
    ax1.axis('off')
    #mask1 = np.flipud(mask1)
    imgax1 = ax1.imshow(img, clim = climits)
    imgax1.set_cmap(new_cmap)
    contours = skimage.measure.find_contours(mask1, 0.8)
    for n, contour in enumerate(contours):
        ax1.plot(contour[:, 1], contour[:, 0], linewidth = 1, color = 'w')
    #cbaxes = fig1.add_axes([0.1, 0.95, 0.4, 0.05])
    #cb = fig1.colorbar(imgax1, cax=cbaxes) #, orientation = 'horizental')
    #cb = fig1.colorbar(imgax1, location = 'top', aspect=2.5, shrink =0.8)#, shrink) #, orientation = 'horizental')
    
    #plot mean pH along y axis
    ax2 = fig1.add_axes([0.45, 0.1, 0.4, 0.8])
    #get mean values
    imgNan = img.copy()
    mask1I = np.invert(mask1)
    imgNan[mask1I] = 0
    mean_data = []
    for y in range(0, img.shape[0]):
        data = imgNan[setrange(y, xrange, limits), :]
        data =data[data>0]
        mean_data.append(np.mean(data))
    ax2.plot(mean_data, range(0, img.shape[0]))
    ax2.set_xlim(climits)
    ax2.set_ylim([0, img.shape[0]])
    ticks = ax2.get_xticks()
    pHticks =[]
    for i in ticks:
        pHticks.append("{0:.2f}".format((float(i) - intercept)/slope))
    ax2.set_xticklabels(pHticks)
    ax2.set_xlabel('pH')
    ax2.get_yaxis().set_visible(False)
    savefile=os.path.join(savefolder, row + 'singleanalysis.jpeg')
    fig1.savefig(savefile)
    
#combine images
plt.close('all')
fig1 = plt.figure(figsize = [4,4], dpi=300)
row = 'img18'
img = imgFrame['Image'].loc[row][:, :, 1].copy()
limits = [0, img.shape[0]]
maskIndex = imgFrame['Mask'].loc[row]
mask1 = np.zeros(img.shape, dtype=np.bool)
mask1[maskIndex] = 1
ax11 = fig1.add_axes([0.1, 0.3, 0.2, 0.54])
ax11.axis('off')
imgax1 = ax11.imshow(img, clim = climits)
imgax1.set_cmap(new_cmap)
ticks1 = np.linspace(climits[0], climits[1], 3)
cb = fig1.colorbar(imgax1, location = 'top', shrink =0.8, ticks = ticks1)#, shrink) #, orientation = 'horizental')
ax11.set_visible('off')
#replace tick marks on colorbar with pH
ticks1 = cb.ax.get_xticks()
ticks = cb.ax.get_xticklabels()
pHticks =[]
for i in ticks:
    pHticks.append("{0:.2f}".format((float(i.get_text()) - intercept)/slope))
#cb.set_ticks(ticks1.tolist())
cb.set_ticklabels(pHticks)
cb.ax.tick_params(labelsize = 8)

ax1 = fig1.add_axes([0.1, 0.1, 0.2, 0.8])
ax1.axis('off')
#mask1 = np.flipud(mask1)
imgax1 = ax1.imshow(img, clim = climits)
imgax1.set_cmap(new_cmap)
contours = skimage.measure.find_contours(mask1, 0.8)
for n, contour in enumerate(contours):
    ax1.plot(contour[:, 1], contour[:, 0], linewidth = 1, color = 'w')

imgL=imgFrame['Image'].loc['img18'][:, :, 1].copy()
row = 'img11'
img = imgFrame['Image'].loc[row][-imgL.shape[0]:, -imgL.shape[1]:, 1].copy()
limits = [0, img.shape[0]]
maskIndex = imgFrame['Mask'].loc[row]
mask1 = np.zeros(img.shape, dtype=np.bool)
mask1[maskIndex] = 1
ax1 = fig1.add_axes([0.3, 0.1, 0.2, 0.8])
ax1.axis('off')
#mask1 = np.flipud(mask1)
imgax1 = ax1.imshow(img, clim = climits)
imgax1.set_cmap(new_cmap)
contours = skimage.measure.find_contours(mask1, 0.8)
for n, contour in enumerate(contours):
    ax1.plot(contour[:, 1], contour[:, 0], linewidth = 1, color = 'w')


#plot mean pH along y axis
ax2 = fig1.add_axes([0.5, 0.1, 0.2, 0.8])
#get mean values
for row in imgFrame.index:
    imgNan = imgFrame['Image'].loc[row][:, :, 1].copy()
    limits = [0, img.shape[0]]
    maskIndex = imgFrame['Mask'].loc[row]
    mask1 = np.zeros(img.shape, dtype=np.bool)
    mask1[maskIndex] = 1
    mask1I = np.invert(mask1)
    imgNan[mask1I] = 0
    mean_data = []
    for y in range(0, img.shape[0]):
        data = imgNan[setrange(y, xrange, limits), :]
        data =data[data>0]
        mean_data.append(np.mean(data))
    ax2.plot(mean_data, range(0, img.shape[0]))
# Calculate the proper aspect for the second axes
'''
aspect0 = ax1.get_aspect()
if aspect0 == 'equal':
    aspect0 = 1.0
dy = np.abs(np.diff(ax2.get_ylim()))
dx = np.abs(np.diff(ax2.get_xlim()))

#aspect = aspect0 / (float(dy) / dx)
aspect = aspect0 / (float(dy) / dx)
ax2.set_aspect(aspect)
'''
ax2.set_aspect(0.53)
ax2.set_xlim(climits)
ax2.set_ylim([0, img.shape[0]])
ticks = ax2.get_xticks()
pHticks =[]
for i in ticks:
    pHticks.append("{0:.2f}".format((float(i) - intercept)/slope))
ax2.set_xticklabels(pHticks)
ax2.set_xlabel('pH')
ax2.get_yaxis().set_visible(False)
savefile=os.path.join(savefolder, 'combineddata.svg')
fig1.savefig(savefile)
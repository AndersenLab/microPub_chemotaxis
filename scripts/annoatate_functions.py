import os
import numpy as np
from skimage.transform import hough_circle
from skimage.feature import peak_local_max, canny
from skimage.feature import blob_dog, blob_log, blob_doh
from skimage.draw import circle_perimeter
from skimage.morphology import binary_closing, binary_dilation, binary_erosion
from skimage.util import img_as_bool, img_as_int, img_as_float
from scipy import ndimage as ndi
from numpy import invert
from skimage.draw import circle
from skimage.io import imread
from clint.textui import puts_err, indent, colored
import matplotlib.pyplot as plt
import warnings
import hashlib
import pickle
from skimage.measure import regionprops
from matplotlib import colors

#########################################################
# define ci function to use below
def ci(v):
    return ((v[0] + v[3]) - (v[1] + v[2])) / float(v[6]) # takes input (v) from ret = [tl, tr, bl, br, n, total_q, total]
###########################################################

img = "/Users/tim/repos/microPub_chemotaxis/data/20180321_crop/1_A_0001.jpg"
radii_range = 900, 940
fp_sigma = 1
# define find plate funciton
def find_plate(img, radii_range, fp_sigma):
    # Read image, and convert to floating point
    img = img_as_bool(imread(img, mode="F", flatten=True))

    # Detect edges
    edges = canny(img, sigma=fp_sigma)

    # Find circles
    hough_radii = np.arange(radii_range[0], radii_range[1], 2)
    hough_res = hough_circle(edges, hough_radii)

    centers, accums, radii = [], [], []

    for radius, h in zip(hough_radii, hough_res):
        # For each radius, extract two circles
        num_peaks = 1
        peaks = peak_local_max(h, num_peaks=num_peaks)
        centers.extend(peaks)
        accums.extend(h[peaks[:, 0], peaks[:, 1]])
        radii.extend([radius] * num_peaks)

    center = np.mean(centers, axis=0)
    radius = (sum(radii) * 1.0) / len(radii)
    return center, radius

# test define plate funtion
# set parameters radii_range
radii_range = 900, 940
img = "/Users/tim/repos/microPub_chemotaxis/data/20180321_crop/1_A_0001.jpg"
small = 100
large = 1200
fp_sigma = 1
extra_crop = 20
#find_plate(radii_range=radii_range, img="/Users/tim/repos/microPub_chemotaxis/data/20180321_crop/1_A_0001.jpg", fp_sigma = 0.5)

# define the crop and filter function. The small and large objects are hard coded.
def crop_and_filter_plate(img, radii_range, extra_crop, small = 100, large = 1200, debug=True):
    # get file name
    fname = os.path.splitext(os.path.basename(img))[0]
    # find the center and radius of the plate in the image using find_plate function
    center, radius = find_plate(img, radii_range, fp_sigma)
    if debug:
        with indent(4):
            puts_err(colored.blue("Center: " + str(center[0]) + "," + str(center[1])))
            puts_err(colored.blue("Radius: " + str(radius)))
            puts_err(colored.blue("Cropping plate"))
    # define the x and y coords of center from fing_plate
    y, x = center

    # Filter background THIS IS NOT FILTERING BACKGROUND IT IS CROPPING THE IMAGE
    img = imread(img, flatten = True) # load the image
    t_crop = int(y - radius) # cropping top
    b_crop = int(y + radius)
    l_crop = int(x - radius)
    r_crop = int(x + radius)
    img = img[t_crop:b_crop]
    img = img[:,l_crop:r_crop] # finalize the image cropping by defining pixel extent for img as circle radius

    if debug:
        with indent(4):
            puts_err(colored.blue("Circular crop"))
        plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/" + fname + ".05_crop.png", img) # save the cropped image

    # Redefine x,y,radius; Generate circle mask.
    mask = np.zeros(img.shape, dtype=np.uint8) # make an array of zeros that fits cropped image pixel dimesnions
    x, y, radius = [img.shape[0]/2] * 3 # set x, y, and radius to 1/2 img dimension in pixels. x, y is the origin of the circle mask with radius. 
    
    # apply the circle mask to image
    rr, cc = circle(y, x, radius-extra_crop) # set pixel coords of circle for mask and crop the radius with extra-crop value
    mask[rr, cc] = 1 # fill mask array with 1s where the plate is
    img[mask == 0] = False # set non-1 pixels in img to false

    if debug:
        with indent(4):
            puts_err(colored.blue("Performing edge detection"))
        plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/" + fname + ".06_mask.png", img)

    # Apply a canny filter
    img = canny(img, sigma=1.5, mask = mask == 1, low_threshold = 0.20, high_threshold = 0.30) # find the edges of objects within the plate

    # Set mask  on canny edge detection img
    mask = np.zeros(img.shape, dtype=np.uint8) # make an array of zeros the size of the img
    rr, cc = circle(y, x, radius-extra_crop) # set pixel coords of circle for mask and crop the radius with extra-crop value
    mask[rr, cc] = 1
    img[mask == 0] = False

    if debug:
        with indent(4):
            puts_err(colored.blue("Binary  Dilation"))
        plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/" + fname + ".07_edges.png", img, cmap='copper')

    # Dilate
    img = binary_dilation(binary_closing(img))

    if debug:
        with indent(4):
            puts_err(colored.blue("Binary Fill"))
        plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/" + fname + ".08_dilation.png", img, cmap='copper')

    # Fill edges
    img = ndi.binary_fill_holes(img)

    if debug:
        with indent(4):
            puts_err(colored.blue("Apply filters"))
        plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/" + fname + ".09_fill.png", img, cmap='copper')

    # Remove small particles
    label_objects, nb_labels = ndi.label(img) # create labels for continous objects identified in img binary (ideally worms)
    sizes = np.bincount(label_objects.ravel()) # count the number of pixels in each bin and return an array with a size for each including the background

    # Label by mask
    reg_props = regionprops(label_objects) # describe regions for all objects see help(regionprops) for details on the properties. Note: background ignored b/c value = 0
    axis_length = np.array([x.minor_axis_length for x in reg_props]) # grab minor axis length
    ecc = np.array([x.eccentricity for x in reg_props]) # grab eccentricity
    solidity = np.array([x.solidity for x in reg_props]) # solidarity
    filters = np.zeros(len(reg_props)+1, dtype='int32') # setup a filter array with 0s with length eaual to reg_groups+1
    filters[filters == 0] = 4 # set all to 4
    filters[sizes < small*5] = 5 # set sizes that are less than 5X small pixel count to 5, this is a worm or worm clump: NEED TO TEST AND MEANS OUR RESOLTION IS HARDCODED!
    filters[sizes < small] = 1 # set small objects to 1, this is debris
    filters[sizes > large] = 2 # set large objects to 2, this is a larger object than expected even accounting for worm clumps.
    filters[0] = 0 # set the backround to 0, i.e., the first number in array (filters[0]).

    if debug:
        if not os.path.exists("debug/" + fname + "/"):
            os.makedirs("debug/" + fname + "/")
        for reg in reg_props:
            plt.imsave("debug/" + fname + "/" + str(reg.label) + ".png", reg.image, cmap='copper') # save each object and name with object label

    filter_img = label_objects.copy()
    for k,v in enumerate(filters):
        filter_img[filter_img == k] = v

    if debug:
        plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/" + fname + ".10_filters.png", filter_img, cmap='copper')

    filters[filters < 4] = 0 # set small (1), large (2), background (0) to zero in filters
    filters[0] = 0 # force background to 0
    img = filters[label_objects] # insert filter values to label_objects and assign to img
    # Rescale and weight
    img[img == 4] = 1 # assign large worm clumps (4) the same value as worms 1
    img[img == 5] = 1 # assign worms or worm clumps (5) to 1

    if debug:
        plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/" + fname + ".11_filtered.png", img, cmap='copper')

    return img # return the result

# set variables for annotation
img = img
n_radius_divisor = 5
extra_crop = 20
img_raw = "/Users/tim/repos/microPub_chemotaxis/data/20180321_crop/1_A_0001.jpg"

###################################################
def pixel_counts(img, n_radius_divisor):
    r = img.shape[0]/2 # set a radius

    mask = np.zeros(img.shape, dtype=np.uint8) # make an empty array of zeros matching img shape
    x, y, radius = [img.shape[0]/2] * 3 # assign x, y, and radius to same value
    radius = radius / n_radius_divisor # divide radius by divisor THE DEFAULT IS SET TO 5
    rr, cc = circle(y, x, radius)
    mask[rr, cc] = 1

    n, q = img.copy(), img.copy()
    q[mask == 1] = 0
    n[mask == 0] = 0
    tl = sum(q[0:r,0:r].flatten()) # add all 1s in quadrant 1 after masked origin in the center is set to 0
    tr = sum(q[0:r,r:].flatten()) # same for other quads
    bl = sum(q[r:,0:r].flatten())
    br = sum(q[r:,r:].flatten())
    plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/q1.png", q[0:r,0:r], cmap = 'copper') # save quadrant
        # debug start
        img_raw = "/Users/tim/repos/microPub_chemotaxis/data/20180321_crop/1_A_0001.jpg"
        img_raw = imread(img_raw, flatten = True) # load the image
        img_raw2 = img_raw[t_crop:b_crop]
        img_raw3 = img_raw2[:,l_crop:r_crop]
        img_raw4 = img_raw3[mask == 1] = 0
        plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/test_mask_on_raw.png", img_raw3) # ah, the inner circle can be plotted like this but we need coords from find_plate()
        # end debug
    plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/q2.png", q[0:r,r:], cmap = 'copper')
    plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/q3.png", q[r:,0:r], cmap = 'copper')
    plt.imsave("/Users/tim/repos/microPub_chemotaxis/data/temp/debug/q4.png", q[r:,r:], cmap = 'copper')
    n = sum(n.flatten()) # n is the pixelcount for the origin. 
    total_q = sum(q.flatten()) # total pixel count in quadrants not the origin
    total = sum(img.flatten()) # total pixel count across entire plate
    ret = [tl, tr, bl, br, n, total_q, total]
    ci_val = ci(ret) # use ci function from above to calculate chemotaxis index
    return ret + [ci_val]

return ((v[0] + v[3]) - (v[1] + v[2])) / float(v[6])

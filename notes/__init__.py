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

_program = "ct"
__version__ = "0.0.1"

def make_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)


def memoize(func):
    def decorated(*args, **kwargs):
        if not os.path.exists('_cache'):
            os.makedirs('_cache')
        digest = hashlib.md5(pickle.dumps([args, kwargs])).hexdigest()
        cache_fname = '_cache/' + digest + ".pkl"
        if os.path.exists(cache_fname):
            with open(cache_fname) as f:
                cache = pickle.load(f)
        else:
            cache = func(*args)
            # update the cache file
            with open(cache_fname, 'wb') as f:
                pickle.dump(cache, f)
        return cache
    return decorated


def suppress_warning(f):
    def decorated(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return f(*args, **kwargs)
    return decorated


def ci(v):
    """
        Calculate the chemotaxis index
    """
    return ((v[0] + v[3]) - (v[1] + v[2])) / float(v[6])


@memoize
@suppress_warning
def find_plate(img, radii_range):
    """
        Identifies the location of the plate
    """
    # Read image, and convert to floating point
    img = img_as_bool(imread(img, mode="F", flatten=True))

    # Detect edges
    edges = canny(img, sigma=2)

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

@suppress_warning
def crop_and_filter_plate(img, radii_range, extra_crop, small = 100, large = 1200, debug= False):
    fname = os.path.splitext(os.path.basename(img))[0]
    center, radius = find_plate(img, radii_range)
    if debug:
        with indent(4):
            puts_err(colored.blue("Center: " + str(center[0]) + "," + str(center[1])))
            puts_err(colored.blue("Radius: " + str(radius)))
            puts_err(colored.blue("Cropping plate"))

    y, x = center

    # Filter background
    img = imread(img, flatten = True)
    t_crop = int(y - radius + extra_crop)
    b_crop = int(t_crop + radius*2 - extra_crop*2)
    l_crop = int(x - radius + extra_crop)
    r_crop = int(x + radius - extra_crop)
    img = img[t_crop:b_crop]
    img = img[:,l_crop:r_crop]

    if debug:
        with indent(4):
            puts_err(colored.blue("Circular crop"))
        plt.imsave("debug/" + fname + ".05_crop.png", img)

    # Redefine x,y,radius; Generate circle mask.
    mask = np.zeros(img.shape, dtype=np.uint8)
    x, y, radius = [img.shape[0]/2] * 3

    rr, cc = circle(y, x, radius)
    mask[rr, cc] = 1
    img[mask == 0] = False

    if debug:
        with indent(4):
            puts_err(colored.blue("Performing edge detection"))
        plt.imsave("debug/" + fname + ".06_mask.png", img)

    # Apply a canny filter
    img = canny(img, sigma=1.5, mask = mask == 1, low_threshold = 0.20, high_threshold = 0.30)

    # Remove the edge
    mask = np.zeros(img.shape, dtype=np.uint8)
    rr, cc = circle(y, x, radius-3)
    mask[rr, cc] = 1
    img[mask == 0] = False

    if debug:
        with indent(4):
            puts_err(colored.blue("Binary  Dilation"))
        plt.imsave("debug/" + fname + ".07_edges.png", img, cmap='copper')

    # Dilate
    img = binary_dilation(binary_closing(img))

    if debug:
        with indent(4):
            puts_err(colored.blue("Binary Fill"))
        plt.imsave("debug/" + fname + ".08_dilation.png", img, cmap='copper')

    # Fill edges
    img = ndi.binary_fill_holes(img)

    if debug:
        with indent(4):
            puts_err(colored.blue("Apply filters"))
        plt.imsave("debug/" + fname + ".09_fill.png", img, cmap='copper')

    # Remove small particles
    label_objects, nb_labels = ndi.label(img)
    sizes = np.bincount(label_objects.ravel())

    # Label by mask
    reg_props = regionprops(label_objects)
    axis_length = np.array([x.minor_axis_length for x in reg_props])
    ecc = np.array([x.eccentricity for x in reg_props])
    solidity = np.array([x.solidity for x in reg_props])
    filters = np.zeros(len(reg_props)+1, dtype='int32')
    filters[filters == 0] = 4
    filters[sizes < small*5] = 5
    filters[sizes < small] = 1
    filters[sizes > large] = 2
    filters[0] = 0

    if debug:
        if not os.path.exists("debug/" + fname + "/"):
            os.makedirs("debug/" + fname + "/")
        for reg in reg_props:
            plt.imsave("debug/" + fname + "/" + str(reg.label) + ".png", reg.image, cmap='copper')

    filter_img = label_objects.copy()
    for k,v in enumerate(filters):
        filter_img[filter_img == k] = v

    if debug:
        plt.imsave("debug/" + fname + ".10_filters.png", filter_img, cmap='copper')

    filters[filters < 4] = 0
    filters[0] = 0
    img = filters[label_objects]

    # Rescale and weight
    img[img == 4] = 1
    img[img == 5] = 1

    if debug:
        plt.imsave("debug/" + fname + ".11_filtered.png", img, cmap='copper')

    return img


def pixel_counts(img, n_radius_divisor):
    r = img.shape[0]/2

    mask = np.zeros(img.shape, dtype=np.uint8)
    x, y, radius = [img.shape[0]/2] * 3
    radius = radius / n_radius_divisor
    rr, cc = circle(y, x, radius)
    mask[rr, cc] = 1

    n, q = img.copy(), img.copy()
    q[mask == 1] = 0
    n[mask == 0] = 0
    tl = sum(q[0:r,0:r].flatten())
    tr = sum(q[0:r,r:].flatten())
    bl = sum(q[r:,0:r].flatten())
    br = sum(q[r:,r:].flatten())
    plt.imsave("debug/q1.png", q[0:r,0:r], cmap = 'copper')
    plt.imsave("debug/q2.png", q[0:r,r:], cmap = 'copper')
    plt.imsave("debug/q3.png", q[r:,0:r], cmap = 'copper')
    plt.imsave("debug/q4.png", q[r:,r:], cmap = 'copper')
    n = sum(n.flatten())
    total_q = sum(q.flatten())
    total = sum(img.flatten())
    ret = [tl, tr, bl, br, n, total_q, total]
    ci_val = ci(ret)
    return ret + [ci_val]



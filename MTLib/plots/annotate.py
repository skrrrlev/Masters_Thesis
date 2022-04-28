
import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial

def annotate_scatter_plot(strings: "list[str]", xs: np.ndarray, ys:np.ndarray, k_nearest: int=3, distance_scale_factor:float=2, pixel_offset:float=12, min_x_offset:float=0, min_y_offset:float=0, verbose:bool=False, **kwargs) -> None:
    '''
    Annotate a scatter plot in such a way, that the annotations are placed away from nearby points.
    - <strings>: list of annotations to add
    - <xs>: numpy array of the x-values
    - <ys>: numpy array of the y-values
    - <k_nearest>: number of nearest points to take into account
    - <distance_scale_factor>: high values will put a larger emphasis on closer points, while low values will put a larger emphasis on distant points.
    - <pixel_offset>: The amount of pixels the centre of the annotation will be from the centre of the data point.
    - <min_x_offset>: Minimum x-distance in pixels from centre of point to centre of annotation. If you have x-uncertainties on your points, variable is very useful.
    - <min_y_offset>: Minimum y-distance in pixels from centre of point to centre of annotation. If you have y-uncertainties on your points, variable is very useful.
    - <verbose>: If True, will print the position and pixel offset associated with each string.

    **kwargs are added to the plt.annotate(...) call: 
    You could for example include: zorder = X, which is nice to make sure that the annotations are on top of the scattered points.
    Or you could set the text-size with size=10.
    '''
    n = len(strings)
    for i in range(n):

        '''Extract the i'th data point'''
        y = ys[i] 
        x = xs[i]

        '''Create a tree-structure with all the other data-points and extact the <k_nearest> nearest'''
        tree = spatial.KDTree(list(zip(np.delete(xs,i), np.delete(ys,i)))) 
        k_nearest_dist, k_nearest_idx = tree.query(np.array([x,y]),k_nearest)

        '''Calculate the weights of each data point based on the distance and the <distance_scale_factor>'''
        inverted_dist = 1/k_nearest_dist**distance_scale_factor
        inverted_dist_sum = np.sum(inverted_dist)
        weights = inverted_dist/inverted_dist_sum
        
        '''Get the weighted average x and y value of the <k_nearest> nearest neighbours'''
        weighted_x, weighted_y = 0,0
        for j in range(k_nearest):
            weighted_x += weights[j]*tree.data[k_nearest_idx[j]][0]
            weighted_y += weights[j]*tree.data[k_nearest_idx[j]][1]

        '''Calculate the distance from the point to the weighted mean'''
        diff_x = x-weighted_x
        diff_y = y-weighted_y
        dist = np.sqrt(diff_x**2+diff_y**2)

        '''Get the pixel offset in x and y'''
        pixel_offset_x = pixel_offset*(diff_x/dist)
        if pixel_offset_x >= 0:
            pixel_offset_x = max((pixel_offset_x,min_x_offset))
        else:
            pixel_offset_x = min((pixel_offset_x,-min_x_offset))
        pixel_offset_y = pixel_offset*(diff_y/dist)
        if pixel_offset_y >= 0:
            pixel_offset_y = max((pixel_offset_y,min_y_offset))
        else:
            pixel_offset_y = min((pixel_offset_y,-min_y_offset))

        if verbose:
            print(f'{i+1}) string: "{strings[i]}"')
            print(f'\tposition      | x: {str(np.round(x,15)):20}  y: {str(np.round(y,15)):20}')
            print(f'\tpixel-offsets | x: {str(np.round(pixel_offset_x,15)):20}  y: {str(np.round(pixel_offset_y,15)):20}')

        '''Add the annotations'''
        plt.annotate(strings[i], (x,y), xytext=(pixel_offset_x,pixel_offset_y), textcoords='offset pixels', ha='center', va='center', **kwargs)

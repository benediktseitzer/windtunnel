import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import csv
from tsp_solver.greedy import solve_tsp

class building():
    """
    Class used to describe the properties of a building. Used to easier plot the different buildings in a flow or
    concentration experiment.
    """
    def __init__(self, name, type, x, y, x_extent, y_extent, z_extent):
        self.name = name
        self.type = type
        self.x_pos = x
        self.y_pos = y
        self.x_extent = x_extent
        self.y_extent = y_extent
        self.z_extent = z_extent # will be probably stay unused for a while (unless side views are needed)

        self.boundaries=[]

        if type == 'rect':
            self.patch = patches.Rectangle((self.x_pos, self.y_pos), self.x_extent, self.y_extent,
                                            edgecolor='none',linewidth=1.5, fill=1,facecolor=[0,0,0,0.5],alpha=0.4)

        self.calc_boundaries()

    def calc_boundaries(self):
        self.boundaries.append([[self.x_pos,self.y_pos],
                                [self.x_pos+self.x_extent,self.y_pos]]) #lower
        self.boundaries.append([[self.x_pos+self.x_extent,self.y_pos],
                                [self.x_pos + self.x_extent, self.y_pos + self.y_extent]]) # right
        self.boundaries.append([[self.x_pos + self.x_extent, self.y_pos + self.y_extent],
                                [self.x_pos, self.y_pos + self.y_extent]])  # upper
        self.boundaries.append([[self.x_pos, self.y_pos + self.y_extent],
                                [self.x_pos, self.y_pos]])  # left

class configuration():
    '''
    Takes positions, extents and types of buildings from a .csv file and generates a configuration which can be used for
    plotting and grid generation. To initialize an object of the configuration class only a .csv file of the building
    parameters is needed. The following structure is required:
    Name    type(rect)  x_pos   y_pos   x_extent    y_extent    z_extent
    '''
    def __init__(self,path):

        self.buildings = []
        with open(path,'r') as file:
            for line in csv.reader(file,
                                   delimiter="\t"):  # You can also use delimiter="\t" rather than giving a dialect.
                self.buildings.append(building(line[0], line[1], float(line[2]), float(line[3]), float(line[4]),
                                                  float(line[5]), float(line[6])))

        self.domain_extents = np.array([np.min([struct.x_pos for struct in self.buildings]) - 100,
                     np.max([struct.x_pos + struct.x_extent for struct in self.buildings]) + 100,
                                      np.min([struct.y_pos for struct in self.buildings]) - 100,
                                      np.max([struct.y_pos + struct.y_extent for struct in self.buildings]) + 100])

    def plot_configuration(self,figure_size):
        fig, ax = plt.subplots(1,figsize=figure_size)
        ax.grid(linewidth=0.5)
        for struct in self.buildings:
            ax.add_patch(struct.patch)

        ax.set_xlim(self.domain_extents[0:2])

        ax.set_ylim(self.domain_extents[2:])
        return fig,ax

    def add_configuration(self,ax):
        ax.grid(linewidth=0.5)
        for struct in self.buildings:
            ax.add_patch(struct.patch)

        ax.set_xlim(self.domain_extents[0:2])

        ax.set_ylim(self.domain_extents[2:])
        return ax

    def filter_points(self,x,y,tolerance=0,scale = 1):
        x = x.astype(float)
        y = y.astype(float)
        for building in self.buildings:
            mask = (x + tolerance > building.x_pos) & (x - tolerance < building.x_pos+building.x_extent) & \
                   (y + tolerance > building.y_pos) & (y - tolerance < building.y_pos+building.y_extent)
            x[mask] = np.nan
            y[mask] = np.nan

        return x,y

    def get_building_boundaries(self):
        boundaries=[]
        for building in self.buildings:
            boundaries.append(building.boundaries)

        return boundaries




def intersects(s0,s1):
    # assumes line segments are stored in the format [(x0,y0),(x1,y1)]
    dx0 = s0[1][0]-s0[0][0]
    dx1 = s1[1][0]-s1[0][0]
    dy0 = s0[1][1]-s0[0][1]
    dy1 = s1[1][1]-s1[0][1]
    p0 = dy1*(s1[1][0]-s0[0][0]) - dx1*(s1[1][1]-s0[0][1])
    p1 = dy1*(s1[1][0]-s0[1][0]) - dx1*(s1[1][1]-s0[1][1])
    p2 = dy0*(s0[1][0]-s1[0][0]) - dx0*(s0[1][1]-s1[0][1])
    p3 = dy0*(s0[1][0]-s1[1][0]) - dx0*(s0[1][1]-s1[1][1])
    return (p0*p1<=0) & (p2*p3<=0)

def cost_func(angles,angle_cost=4):
    '''
    arbitrary cost function which punishes small angles: for visuals follow
    plt.figure()
    plt.plot(np.arange(91),(np.cos(np.deg2rad(np.arange(91)*4))+0.2)*((np.cos(np.deg2rad(np.arange(91)*4))+0.2)>0)*4)
    :param angles:
    :return:cost
    '''
    cost = (np.cos(np.deg2rad(angles*4))+0.2)*((np.cos(np.deg2rad(angles*4))+0.2)>0)*angle_cost
    cost[cost<1] = 1
    cost[angles==0] = 1
    cost[angles==90]= 1

    return cost

def get_metangle(x, y):
    '''Get meteorological angle of input vector.
    Args:
        x (numpy.ndarray): X-components of input vectors.
        y (numpy.ndarray): Y-components of input vectors.
    Returns:
        (numpy.ma.core.MaskedArray): Meteorological angles.
    '''
    mask = np.logical_and(x == 0, y == 0)  # Vectors (0, 0) not valid.
    met_ang = ma.masked_array((90 - np.degrees(np.arctan2(x, y)) + 360) % 360, mask=mask,
                                  fill_value=np.nan)
    return met_ang

def optimize_grid(points,configuration,avoid_buildings=True,angle_cost=10):

    if points.shape[1]==3:
        points2d = points[:, :2] # throwing away third dimension (z)


    obstacle = np.zeros([len(points2d), len(points2d)])

    for i, point in enumerate(points2d):
        angles = np.rad2deg(np.arcsin(np.abs(point - points2d)[:, 0] /
                                         (np.sqrt((point - points2d)[:, 0] ** 2 + (point - points2d)[:, 1] ** 2))))
        angles = np.nan_to_num(angles)

        for j, node in enumerate(points2d):
            obstacle[i, j] = np.any([intersects([point, node], side) for building in configuration.get_building_boundaries()
                                     for side in building])

        if i == 0:
            dist_all = np.hypot((point - points2d)[:, 0], (point - points2d)[:, 1]) * cost_func(angles, angle_cost)
        else:
            dist_all = np.vstack(
                [dist_all, np.hypot((point - points2d)[:, 0], (point - points2d)[:, 1]) * cost_func(angles, angle_cost)])
    if avoid_buildings:
        dist_all[obstacle==1] = 1000000

    path = solve_tsp(dist_all, endpoints=(0, None))

    return points[path]
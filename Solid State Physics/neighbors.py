# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 12:47:30 2020

@author: Yuchen Wang
"""

import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import cv2

class Lattice():
    '''
    A class to represent a Bravais Lattice.
    '''
    #general functions to generate new points for a Bravais Lattice
    def bcc_generator(amin, amax):
        '''
        Generate new points for BCC inside a cube with side length amin
        and outside a cube with side length amax.
        '''
        assert amin == 0, 'Not implemented for amin > 0'
        assert amin % 2 == 0 and amax % 2 == 0 and amax >= 4, 'Size should be even numbers that >= 4.'
        xs = np.arange(-amax, amax, 2)
#        xs = np.arange(amin, amax, 2)
        xs = np.unique(np.concatenate((-xs[::-1], xs)))
        x1, x2, x3 = np.meshgrid(xs, xs, xs)
        x1, x2, x3 = x1.reshape((-1,1)), x2.reshape((-1,1)), x3.reshape((-1,1))
        gpoints = np.concatenate((x1, x2, x3), axis=1)
#        points = points[((x1>=amin) | (x2>=amin) | (x3>=amin)).flatten()]
        ngpoints = gpoints+np.array((1,1,1))
        points = np.concatenate((gpoints, ngpoints), axis=0)
        return points, gpoints, ngpoints
        
    def fcc_generator(amin, amax):
        '''
        Generate new points for FCC inside a cube with side length amin
        and outside a cube with side length amax.
        '''
        assert amin == 0, 'Not implemented for amin > 0'
        assert amin % 2 == 0 and amax % 2 == 0 and amax >= 4, 'Size should be even numbers that >= 4.'
#        xs = np.arange(amin, amax, 1)
        xs = np.arange(-amax, amax, 1)
        xs = np.unique(np.concatenate((-xs[::-1], xs)))
        x1, x2, x3 = np.meshgrid(xs, xs, xs)
        x1, x2, x3 = x1.reshape((-1,1)), x2.reshape((-1,1)), x3.reshape((-1,1))
        points = np.concatenate((x1, x2, x3), axis=1)
#        points = points[((x1>=amin) | (x2>=amin) | (x3>=amin)).flatten()]
        is_point = np.sum(points, axis=1) % 2 == 0
        points = points[is_point]
        gpoints = points[np.all(np.mod(points, 2)==0, axis=1)]
        ngpoints = points[np.any(np.mod(points, 2)==1, axis=1)]
        return points, gpoints, ngpoints
        
    lat_generators = {'bcc': bcc_generator,
                     'fcc': fcc_generator}
    
    def __init__(self, size, lat_class='bcc'):
        self.lat_class = lat_class
        self.generator = Lattice.lat_generators[lat_class]
        self.size = size
        self.array, self.gpoints, self.ngpoints = self.generator(0, size)
        self.ds = np.sqrt(np.sum(self.array**2, axis=1))
        
    def enlarge(self, l):
#        new_points = self.generator(self.size, self.size+l)
        new_points, new_gpoints, new_ngpoints = self.generator(0, self.size+l)
        new_ds = np.sqrt(np.sum(new_points**2, axis=1))
        self.size += l
#        self.array = np.concatenate((self.array, new_points), axis=0)
#        self.ds = np.concatenate((self.ds, new_ds))
        self.array = new_points
        self.gpoints, self.ngpoints = new_gpoints, new_ngpoints
        self.ds = new_ds
        
    def find_neighbors(self, n_shells):
        self.n_shells = n_shells
        self.shells = []
        self.distances = []
        self.numbers = []
        self.cul_numbers = []
        self.points = []
        n = 0
        
        self.unique_ds = np.unique(np.sort(self.ds))[1:]
        assert self.unique_ds[0] > 0
        self.valid_ds = np.sum(self.unique_ds<self.size-2)
        
        while self.valid_ds < n_shells:
            self.enlarge(4)
            self.unique_ds = np.unique(np.sort(self.ds))[1:]
            assert self.unique_ds[0] > 0
            self.valid_ds = np.sum(self.unique_ds<self.size-2)
        
        for n_shell, distance in zip(range(n_shells), self.unique_ds[:n_shells]):
            self.shells.append(n_shell+1)
            self.distances.append(distance)
            points = self.array[self.ds == distance]
            self.points.append(points)
            self.numbers.append(points.shape[0])
            n += points.shape[0]
            self.cul_numbers.append(n)
#        self.points = np.array(self.points).astype('int')
        
        self.result = pd.DataFrame({'shell': self.shells,
                                    'distance': self.distances,
                                    'number': self.numbers,
                                    'total_number': self.cul_numbers, 
                                    'coordinates': self.points}) #columns=['shell', 'distance', 'number', 'coordinates'])
        return self.result
    
    def visualize(self, n_shells=None, save_each_frame=False, fps=1/3, size=(256, 144), save=False, fname='temp'):
        if n_shells==None:
            n_shells = self.n_shells
        fig = plt.figure()
        ax = Axes3D(fig)
#        useful_points = self.points
#        fps = 10
#        size = (1280, 720)
        all_points = []
        if save:
            videowriter = cv2.VideoWriter(fname, cv2.VideoWriter_fourcc('M', 'J', 'P', 'G'), fps, size)
        ax.scatter([0],[0],[0])
        for shell, points, dis in zip(self.shells[:n_shells], self.points[:n_shells], self.distances[:n_shells]):
            size = np.floor(dis/2)*2
            useful_points = self.array[np.all(np.abs(self.array)<size, axis=1)]
            for line in ax.lines:
                if line.get_color() == 'r':
                    line.set_alpha(0)
#            useful_points = points
            ax.scatter(useful_points[:,0], useful_points[:,1], useful_points[:,2], color='black')
            ax.scatter(points[:,0], points[:,1], points[:,2], color='black')
            self.plot_grid(ax, size=size, gpoints = True if self.lat_class == 'fcc' else False)
#            all_points = points #.append(points)
#            grid = all_points[np.all(np.mod(all_points, 2) == 1, axis=1)]
#            ax.scatter(grid[:,0], grid[:,1], grid[:,2], color='blue')
#            not_grid = all_points[np.any(np.mod(all_points, 2) != 1, axis=1)]
#            ax.scatter(not_grid[:,0], not_grid[:,1], not_grid[:,2], color='black')
            
            for point in points:
                plt.plot([0, point[0]], [0, point[1]], [0, point[2]], color='r')
            plt.title('{}, {} nearest neighbors: {} points, distance = {:.3f}'.format(fname, shell, points.shape[0], dis))
            if save_each_frame:
                plt.savefig('{}_{} nearest neighbors- {} points, distance = {:.3f}.png'.format(fname, shell, points.shape[0], dis))
            plt.pause(3)
    
    def plot_grid(self, ax, size=np.inf, gpoints=True, color='green', linewidth=1):
        p = self.gpoints if gpoints else self.ngpoints
        p = p[np.all(np.abs(p)<=size, axis=1)]
        xs, ys, zs = np.unique(p[:,0]), np.unique(p[:,1]), np.unique(p[:,2])
        for y in ys:
            for z in zs:
                ax.plot([np.min(xs), np.max(xs)], [y,y],[z,z], color=color, linewidth=linewidth)
        for z in zs:
            for x in xs:
                ax.plot([x, x], [np.min(ys),np.max(ys)],[z,z], color=color, linewidth=linewidth)
        for x in xs:
            for y in ys:
                ax.plot([x, x], [y,y],[np.min(zs),np.max(zs)], color=color, linewidth=linewidth)
                


if __name__ == '__main__':
    bcc=Lattice(4, 'bcc')
    bcc_result = bcc.find_neighbors(20)
    bcc_result.to_csv('bcc_data.csv')
    bcc.visualize(n_shells=3, fname='bcc', save_each_frame=True)
#    print(bcc.result)
    fcc=Lattice(4, 'fcc')
    fcc_result = fcc.find_neighbors(20)
    fcc_result.to_csv('fcc_data.csv')
    fcc.visualize(n_shells=3, fname='fcc', save_each_frame=True)

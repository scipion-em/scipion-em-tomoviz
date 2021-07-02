# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


import pyvistaqt as pvqt
import pyvista as pv
import pymeshfix as pm

import numpy as np
import matplotlib

from scipy.ndimage import zoom, gaussian_filter
# from scipy.spatial.distance import pdist

from skimage import measure
from skimage.morphology import binary_dilation, binary_erosion, ball

import pyworkflow.utils as pwutils

from pwem.emlib.image import ImageHandler


class MrcPlot(object):
    '''
    Class to visualize MRC files
    Input paramters:
         - tomo_mrc (Path (Str) - Optional): File containing a Volume (in MRC format)
         - mask_mrc (Path (Str) - Optional): File containing a Mask (in MRC format)
         - points (Path (Str) - Optional): File containing a SetOfCoordinates3D (in TEXT format)
         - normals (Path (Str) - Optional): File containing a Set of Normals (in TEXT format)
         - binning (Float - Optional):  Binning factor to be applied to tomo_mrc and mask_mrc (Very useful to save time)
         - sigma (Float - Optional):  Gaussian Filter width
         - triangulation (Bool - Optional): Tells if the representation of the Tomomgra is Voxel based on Triangle based
    Usage:
         import MRCPlot
         plt = MRCPlot(tomo_mrc=tomo_mrc, mask_mrc=mask_mrc, binning=2)
         plt.initializePlot()
    '''

    def __init__(self, tomo_mrc=None, mask_mrc=None, points=None, normals=None,
                 binning=None, sigma=1., triangulation=False):
        if binning is None:
            if tomo_mrc is not None:
                self.binning = self.getBinning(tomo_mrc)
            elif mask_mrc is not None:
                self.binning = self.getBinning(mask_mrc)
            else:
                self.binning = 0
        self.tomo = self.readMRC(tomo_mrc, order=5, binning=self.binning) if tomo_mrc is not None else None
        self.mask = self.readMRC(mask_mrc, binning=self.binning) if mask_mrc is not None else None
        self.points = np.loadtxt(points, delimiter=' ') if points is not None else None
        self.normals = np.loadtxt(normals, delimiter=' ') if normals is not None else None
        self.save_basename = pwutils.removeBaseExt(tomo_mrc) if tomo_mrc is not None and points is not None else None

        # Get Pyvista Objects
        if isinstance(self.tomo, np.ndarray):
            self.pv_tomo_slice = self.gridFromMRC(self.tomo)
            self.pv_tomo, self.opacities = self.isovolumes(self.tomo, triangulation=triangulation, sigma=sigma)
        if isinstance(self.mask, np.ndarray):
            labels = np.unique(self.mask)[1:]
            self.pv_masks = [self.surfaceFromMRC(self.mask, label=label) for label in labels]
        if isinstance(self.points, np.ndarray):
            self.points_ids = self.points[:, 3]
            self.points = np.column_stack([self.points[:, 1], self.points[:, 0], self.points[:, 2]])
            self.points /= 2 ** self.binning  # Binning Scaling
            self.pv_points = pv.PolyData(self.points)
        if isinstance(self.normals, np.ndarray):
            self.normals = np.column_stack([self.normals[:, 1], self.normals[:, 0], self.normals[:, 2]])
            # vecLength = np.amax(pdist(self.pv_points.points))
            self.normals /= np.linalg.norm(self.normals, axis=1)[:, np.newaxis]
            # self.normals *= vecLength
            self.pv_normals = pv.pyvista_ndarray(self.normals)

        self.tomo_actor = []
        self.tomo_slice_actor = None
        self.mask_actors = []
        self.points_actor = None
        self.normals_actor = None

        self.plt = pvqt.BackgroundPlotter()
        self.plt.main_menu.clear()

        pos = 0.

        if isinstance(self.tomo, np.ndarray):
            pos += 45.
            self.buttonTomo = self.plt.add_checkbox_button_widget(callback=self.plotTomo, position=(pos, 10.))
            self.plt.add_text('Tomogram', position=(pos, 65.), font_size=12)
            pos += 170.
            self.buttonSliceTomo = self.plt.add_checkbox_button_widget(callback=self.toogleSlice, position=(pos, 10.))
            self.plt.add_text('Slice Mode', position=(pos, 65.), font_size=12)

        if isinstance(self.mask, np.ndarray):
            pos += 170. if pos != 0 else 45.
            self.plt.add_text('Vesicles', position=(pos, 65.), font_size=12)
            self.buttonGraph = self.plt.add_checkbox_button_widget(callback=self.plotMasks, position=(pos, 10.))

        if isinstance(self.points, np.ndarray):
            pos += 170. if pos != 0 else 45.
            self.plt.add_text('Coords', position=(pos, 65.), font_size=12)
            self.buttonPoints = self.plt.add_checkbox_button_widget(callback=self.plotPoints, position=(pos, 10.))

            # Picking Callbacks
            def removeSelection(selection):
                self.pv_points.remove_cells(selection.active_scalars, inplace=True)
                self.points_ids = np.delete(self.points_ids, selection.active_scalars)
                if self.normals is not None:
                    self.pv_normals = np.delete(self.pv_normals, selection.active_scalars, 0)
                    self.plt.remove_actor(self.normals_actor)
                    if self.buttonNormals.GetRepresentation().GetState():
                        self.normals_actor = self.plt.add_arrows(self.pv_points.cell_centers().points, self.pv_normals,
                                                                 mag=10, color='red', reset_camera=False)

            def enableRemoveSelection():
                self.plt.enable_cell_picking(mesh=self.pv_points,
                                             callback=removeSelection,
                                             font_size=12, opacity=0)

            # Picking Controls
            self.plt.main_menu.addAction('Point Removal', enableRemoveSelection)

        if isinstance(self.normals, np.ndarray):
            pos += 170. if pos != 0 else 45.
            self.plt.add_text('Directions', position=(pos, 65.), font_size=12)
            self.buttonNormals = self.plt.add_checkbox_button_widget(callback=self.plotNormals, position=(pos, 10.))

    def getBinning(self, file):
        dim = ImageHandler().read(file + ':mrc').getDimensions()
        return int(np.floor(max(dim) / 400))

    def readMRC(self, file, binning=1, order=0, swapaxes=True):
        image = ImageHandler().read(file + ':mrc')
        data = np.squeeze(image.getData())
        data = zoom(data, 1 / (2 ** binning), order=order, prefilter=False)
        data = np.swapaxes(np.swapaxes(data, 0, 2), 1, 0) if swapaxes else data
        return data

    def gridFromMRC(self, data):
        '''Function to convert an MRC file into an Structure Grid in VTK'''

        # Create a Grid from the data
        x = np.arange(0, data.shape[0] + 1)
        y = np.arange(0, data.shape[1] + 1)
        z = np.arange(0, data.shape[2] + 1)
        x, y, z = np.meshgrid(x, y, z)

        # Copy the Grid to an Structure Grid VTK object (Pyvista)
        grid = pv.StructuredGrid(x, y, z)

        # Set the cell values. Previous reordering of the axis was needed to flatten properly the array
        grid.cell_arrays["values"] = data.flatten(order="K")

        return grid

    def surfaceFromMRC(self, data, label=1):
        '''Function to convert an MRC file into an Structure Surface in VTK'''

        # Get Only mesh corresponding to a given label (smooth the result to fill holes in the mask)
        data = data == label
        data = binary_erosion(binary_dilation(data, selem=ball(4)), selem=ball(1))

        # Triangulate coordinates using marching cubes algorithm
        grid = self.marchingCubes(data)

        # Fix the mesh
        mfix = pm._meshfix.PyTMesh(False)
        mfix.load_array(grid.points, grid.faces.reshape((-1, 4))[:, 1:])
        mfix.join_closest_components()
        mfix.fill_small_boundaries(refine=True)
        vert, faces = mfix.return_arrays()

        # Reconstruct fixed mesh
        triangles = np.empty((faces.shape[0], 4))
        triangles[:, -3:] = faces
        triangles[:, 0] = 3
        grid = pv.PolyData(vert, triangles.astype(int))

        # Final smoothing using pyvista
        grid = grid.smooth(n_iter=1000).clean()

        return grid

    def marchingCubes(self, volume, level=None, triangulation=True):
        vertices, faces = measure.marching_cubes(volume, spacing=(1, 1, 1), level=level, allow_degenerate=False)[:2]
        faces = np.column_stack((3 * np.ones((len(faces), 1), dtype=np.int), faces)).flatten()
        grid = pv.PolyData(vertices.astype(int), faces) if triangulation else pv.PolyData(vertices.astype(int))
        return grid

    def histogram(self, volume):
        hist, edges = np.histogram(volume, bins=100)
        hist = hist / np.sum(hist)
        bin_centers = np.mean(np.vstack([edges[0:-1], edges[1:]]), axis=0)
        return hist, bin_centers

    def contours(self, hist, bin_centers):
        opacities = 1 - hist[np.where((hist > np.std(hist)) * (hist < np.amax(hist)))]
        contour_values = bin_centers[np.where((hist > np.std(hist)) * (hist < np.amax(hist)))]
        return contour_values, opacities

    def isovolumes(self, volume, range=0.01, sigma=None, triangulation=True):
        volume = volume if sigma is None else gaussian_filter(volume, sigma=sigma)
        hist, bin_centers = self.histogram(volume)
        contour_values, opacities = self.contours(hist, bin_centers)
        opacities = (range / (np.amax(opacities) - np.amin(opacities))) * (opacities - np.amin(opacities))
        return [self.marchingCubes(volume, level, triangulation) for level in contour_values], opacities

    def downsamplingPC(self, coords, voxel_size):
        non_empty_voxel_keys, inverse, nb_pts_per_voxel = np.unique(((coords - np.min(coords, axis=0))
                                                                     // voxel_size).astype(int), axis=0,
                                                                    return_inverse=True,
                                                                    return_counts=True)
        idx_pts_vox_sorted = np.argsort(inverse)
        voxel_grid = {}
        grid_barycenter, grid_candidate_center = [], []
        last_seen = 0
        for idx, vox in enumerate(non_empty_voxel_keys):
            voxel_grid[tuple(vox)] = coords[idx_pts_vox_sorted[last_seen:last_seen + nb_pts_per_voxel[idx]]]
            grid_barycenter.append(np.mean(voxel_grid[tuple(vox)], axis=0))
            last_seen += nb_pts_per_voxel[idx]
        return np.asarray(grid_barycenter)

    def plotTomo(self, value):
        if value:
            cmap = matplotlib.cm.get_cmap('bone')  # Bone also looks nice
            cmap_ids = np.linspace(0, 1, len(self.pv_tomo))
            self.tomo_actor = [self.plt.add_mesh(actor, show_scalar_bar=False, opacity=op, color=cmap(cid),
                                                 render_points_as_spheres=True)
                               for actor, op, cid in zip(self.pv_tomo, self.opacities, cmap_ids)]
            self.plt.reset_camera()
        else:
            for actor in self.tomo_actor:
                self.plt.remove_actor(actor)
            self.tomo_actor = []

    def toogleSlice(self, value):
        if value:
            self.tomo_slice_actor = self.plt.add_mesh_slice(self.pv_tomo_slice, normal='z', cmap="gray", show_scalar_bar=False,
                                                            outline_translation=False, origin_translation=False)
        else:
            self.plt.remove_actor(self.tomo_slice_actor)
            self.plt.clear_plane_widgets()
            self.tomo_slice_actor = None

    def plotMasks(self, value):
        if value:
            cmap = matplotlib.cm.get_cmap('Set3')
            cmap_ids = np.linspace(0, 1, len(self.pv_masks))
            self.mask_actors = [self.plt.add_mesh(mask, show_scalar_bar=False, color=cmap(cmap_id), smooth_shading=True)
                                for mask, cmap_id in zip(self.pv_masks, cmap_ids)]
        else:
            for actor in self.mask_actors:
                self.plt.remove_actor(actor)
            self.graph_actor = []

    def plotPoints(self, value):
        if value:
            self.points_actor = self.plt.add_mesh(self.pv_points, show_scalar_bar=False, color='orange',
                                                  render_points_as_spheres=True, reset_camera=False)
        else:
            self.plt.remove_actor(self.points_actor)
            self.points_actor = None

    def plotNormals(self, value):
        if value:
            self.normals_actor = self.plt.add_arrows(self.pv_points.cell_centers().points, self.pv_normals,
                                                     mag=10, color='red', reset_camera=False)
        else:
            self.plt.remove_actor(self.normals_actor)
            self.normals_actor = None

    def initializePlot(self):
        self.plt.app.exec_()

        # Save Points and Normals
        if self.save_basename is not None:
            np.savetxt(self.save_basename + '_indices.txt', self.points_ids, delimiter=' ')


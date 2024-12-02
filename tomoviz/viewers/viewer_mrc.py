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


import os
import pyvistaqt as pvqt
import pyvista as pv
from PyQt5.QtGui import QMovie
from pyvista.utilities import generate_plane
import pymeshfix as pm
import vtk
from PyQt5 import QtGui
from PyQt5.QtWidgets import QWidget, QLabel, QDesktopWidget
from PyQt5.QtCore import Qt, QObject, QThread, pyqtSignal, QRectF

import numpy as np
import matplotlib

from scipy.ndimage import zoom, gaussian_filter
# from scipy.spatial.distance import pdist

from skimage import measure
from skimage.morphology import binary_dilation, binary_erosion, ball

import pyworkflow.utils as pwutils

from pwem.emlib.image import ImageHandler

import tomoviz

COLOR_MAP = 'gist_rainbow_r' #'terrain'


class MrcPlot(object):
    '''
    Class to visualize MRC files
    Input paramters:
         - tomo_mrc (Path (Str) - Optional): File containing a Volume (in MRC format)
         - mask_mrc (Path (Str) - Optional): File containing a Mask (in MRC format)
         - points (Path (Str) - Optional): File containing a SetOfCoordinates3D (in TEXT format)
         - boxSize (Double - Optional): Box size associated to points
         - normals (Path (Str) - Optional): File containing a Set of Normals (in TEXT format)
         - binning (Float - Optional):  Binning factor to be applied to tomo_mrc and mask_mrc (Very useful to save time)
         - sigma (Float - Optional):  Gaussian Filter width
         - triangulation (Bool - Optional): Tells if the representation of the Tomogram is Voxel based on Triangle based
    Usage:
         import MRCPlot
         plt = MRCPlot(tomo_mrc=tomo_mrc, mask_mrc=mask_mrc, binning=2)
         plt.initializePlot()
    '''

    def __init__(self, tomo_mrc=None, mask_mrc=None, points=None, boxSize=None, normals=None,
                 binning=None, sigma=1., triangulation=False, discrete=False):
        if binning is None:
            if tomo_mrc is not None:
                self.binning = self.getBinning(tomo_mrc)
            elif mask_mrc is not None:
                self.binning = self.getBinning(mask_mrc)
            else:
                self.binning = 0
        else:
            self.binning = binning
        self.tomo = (tomo_mrc, triangulation, sigma)
        self.mask = mask_mrc
        self.points = np.loadtxt(points, delimiter=' ') if points is not None else None
        self.normals = np.loadtxt(normals, delimiter=' ') if normals is not None else None
        self.boxSize = boxSize / 2 ** self.binning if boxSize is not None else None
        self.save_basename = pwutils.removeBaseExt(tomo_mrc) if tomo_mrc is not None and points is not None else None
        self.discrete = discrete

        # Get Pyvista Objects
        if isinstance(self.points, np.ndarray):
            self.points_ids = self.points[:, 3]
            self.group_ids = self.points[:, 4]
            self.points = np.column_stack([self.points[:, 1], self.points[:, 0], self.points[:, 2]])
            self.points /= 2 ** self.binning  # Binning Scaling
            self.pv_points = pv.PolyData(self.points)
            scalar_colors = np.zeros([self.points.shape[0]])
            unique_ids = np.unique(self.group_ids)
            for group_id in unique_ids:
                idp = np.where(self.group_ids == group_id)
                scalar_colors[idp] = group_id
            self.pv_points["colors"] = scalar_colors
        if isinstance(self.normals, np.ndarray):
            self.normals = np.column_stack([self.normals[:, 1], self.normals[:, 0], self.normals[:, 2]])
            # vecLength = np.amax(pdist(self.pv_points.points))
            self.normals /= np.linalg.norm(self.normals, axis=1)[:, np.newaxis]
            # self.normals *= vecLength
            self.pv_normals = pv.pyvista_ndarray(self.normals)

        self.tomo_actor = []
        self.tomo_slice_actor = None
        self.mask_actors = []
        self.points_actor = []
        self.normals_actor = None
        self.box_actor = {}

        self.first_reset = True

        # Theme
        from pyvista.themes import DocumentTheme
        customTheme = DocumentTheme()
        customTheme.cmap = COLOR_MAP
        customTheme.color_cycler = 'default'
        pv.set_plot_theme(customTheme)

        self.plt = pvqt.BackgroundPlotter(title='Scipion tomoviz viewer')
        self.plt.main_menu.clear()
        plugin_path = os.path.dirname(tomoviz.__file__)
        self.plt.app.setWindowIcon(QtGui.QIcon(os.path.join(plugin_path, "icon_square.png")))
        self.loading_screen = LoadingScreen()

        pos = 0.

        if tomo_mrc is not None:
            pos += 45.
            self.buttonTomo = self.plt.add_checkbox_button_widget(callback=self.plotTomo, position=(pos, 10.))
            self.plt.add_text('Tomogram', position=(pos, 65.), font_size=12)
            pos += 170.
            self.buttonSliceTomo = self.plt.add_checkbox_button_widget(callback=self.toogleSlice, position=(pos, 10.))
            self.plt.add_text('Slice mode', position=(pos, 65.), font_size=12)

            self.plt.clear_events_for_key("Up")
            # Callback to move the slice with the up arrows in the keyboard
            def sliceUp():
                if self.tomo_slice_actor is not None:
                    alg_tomo = vtk.vtkCutter()
                    current_plane = vtk.vtkPlane()
                    alg_tomo.SetInputDataObject(self.pv_tomo_slice)
                    plane_sliced_tomo = self.plt.plane_sliced_meshes[-1]
                    plane_widget = self.plt.plane_widgets[-1]
                    plane_widget.GetPlane(current_plane)
                    origin = np.asarray(current_plane.GetOrigin())
                    normal = np.asarray(current_plane.GetNormal())
                    # Normalize normal
                    normal /= np.linalg.norm(normal)
                    # move plane one unit in the direction of the normal
                    origin += normal
                    # Create a new plane to update the position and perform the update
                    plane = generate_plane(normal, origin)
                    alg_tomo.SetCutFunction(plane)
                    alg_tomo.Update()
                    plane_sliced_tomo.shallow_copy(alg_tomo.GetOutput())
                    plane_widget.SetOrigin(origin)
                    plane_widget.SetNormal(normal)
                    plane_widget.UpdatePlacement()
            self.plt.add_key_event("Up", sliceUp)

            self.plt.clear_events_for_key("Down")
            # Callback to move the slice with the down arrows in the keyboard
            def sliceDown():
                if self.tomo_slice_actor is not None:
                    alg_tomo = vtk.vtkCutter()
                    current_plane = vtk.vtkPlane()
                    alg_tomo.SetInputDataObject(self.pv_tomo_slice)
                    plane_sliced_tomo = self.plt.plane_sliced_meshes[-1]
                    plane_widget = self.plt.plane_widgets[-1]
                    plane_widget.GetPlane(current_plane)
                    origin = np.asarray(current_plane.GetOrigin())
                    normal = np.asarray(current_plane.GetNormal())
                    # Normalize normal
                    normal /= np.linalg.norm(normal)
                    # move plane one unit in the direction of the normal
                    origin -= normal
                    # Create a new plane to update the position and perform the update
                    plane = generate_plane(normal, origin)
                    alg_tomo.SetCutFunction(plane)
                    alg_tomo.Update()
                    plane_sliced_tomo.shallow_copy(alg_tomo.GetOutput())
                    plane_widget.SetOrigin(origin)
                    plane_widget.SetNormal(normal)
                    plane_widget.UpdatePlacement()
            self.plt.add_key_event("Down", sliceDown)

            # Callback to place plane normal along X
            def sliceX():
                if self.tomo_slice_actor is not None:
                    alg_tomo = vtk.vtkCutter()
                    current_plane = vtk.vtkPlane()
                    alg_tomo.SetInputDataObject(self.pv_tomo_slice)
                    plane_sliced_tomo = self.plt.plane_sliced_meshes[-1]
                    plane_widget = self.plt.plane_widgets[-1]
                    plane_widget.GetPlane(current_plane)
                    origin = np.asarray(current_plane.GetOrigin())
                    normal = np.asarray([1, 0, 0])
                    # Create a new plane to update the position and perform the update
                    plane = generate_plane(normal, origin)
                    alg_tomo.SetCutFunction(plane)
                    alg_tomo.Update()
                    plane_sliced_tomo.shallow_copy(alg_tomo.GetOutput())
                    plane_widget.SetOrigin(origin)
                    plane_widget.SetNormal(normal)
                    plane_widget.UpdatePlacement()
            self.plt.add_key_event("x", sliceX)

            # Callback to place plane normal along Y
            def sliceY():
                if self.tomo_slice_actor is not None:
                    alg_tomo = vtk.vtkCutter()
                    current_plane = vtk.vtkPlane()
                    alg_tomo.SetInputDataObject(self.pv_tomo_slice)
                    plane_sliced_tomo = self.plt.plane_sliced_meshes[-1]
                    plane_widget = self.plt.plane_widgets[-1]
                    plane_widget.GetPlane(current_plane)
                    origin = np.asarray(current_plane.GetOrigin())
                    normal = np.asarray([0, 1, 0])
                    # Create a new plane to update the position and perform the update
                    plane = generate_plane(normal, origin)
                    alg_tomo.SetCutFunction(plane)
                    alg_tomo.Update()
                    plane_sliced_tomo.shallow_copy(alg_tomo.GetOutput())
                    plane_widget.SetOrigin(origin)
                    plane_widget.SetNormal(normal)
                    plane_widget.UpdatePlacement()
            self.plt.add_key_event("y", sliceY)

            # Callback to place plane normal along Z
            def sliceZ():
                if self.tomo_slice_actor is not None:
                    alg_tomo = vtk.vtkCutter()
                    current_plane = vtk.vtkPlane()
                    alg_tomo.SetInputDataObject(self.pv_tomo_slice)
                    plane_sliced_tomo = self.plt.plane_sliced_meshes[-1]
                    plane_widget = self.plt.plane_widgets[-1]
                    plane_widget.GetPlane(current_plane)
                    origin = np.asarray(current_plane.GetOrigin())
                    normal = np.asarray([0, 0, 1])
                    # Create a new plane to update the position and perform the update
                    plane = generate_plane(normal, origin)
                    alg_tomo.SetCutFunction(plane)
                    alg_tomo.Update()
                    plane_sliced_tomo.shallow_copy(alg_tomo.GetOutput())
                    plane_widget.SetOrigin(origin)
                    plane_widget.SetNormal(normal)
                    plane_widget.UpdatePlacement()
            self.plt.add_key_event("z", sliceZ)

            # Callback to reset plane origin
            def reserOrigin():
                if self.tomo_slice_actor is not None:
                    alg_tomo = vtk.vtkCutter()
                    current_plane = vtk.vtkPlane()
                    alg_tomo.SetInputDataObject(self.pv_tomo_slice)
                    plane_sliced_tomo = self.plt.plane_sliced_meshes[-1]
                    plane_widget = self.plt.plane_widgets[-1]
                    plane_widget.GetPlane(current_plane)
                    origin = np.asarray(self.pv_tomo_slice.center)
                    normal = np.asarray(np.asarray(current_plane.GetNormal()))
                    # Create a new plane to update the position and perform the update
                    plane = generate_plane(normal, origin)
                    alg_tomo.SetCutFunction(plane)
                    alg_tomo.Update()
                    plane_sliced_tomo.shallow_copy(alg_tomo.GetOutput())
                    plane_widget.SetOrigin(origin)
                    plane_widget.SetNormal(normal)
                    plane_widget.UpdatePlacement()
            self.plt.add_key_event("o", reserOrigin)

        if mask_mrc is not None:
            pos += 170. if pos != 0 else 45.
            self.plt.add_text('Vesicles', position=(pos, 65.), font_size=12)
            self.buttonGraph = self.plt.add_checkbox_button_widget(callback=self.plotMasks, position=(pos, 10.))

        if isinstance(self.points, np.ndarray):
            pos += 170. if pos != 0 else 45.
            self.plt.add_text('Coords', position=(pos, 65.), font_size=12)
            self.buttonPoints = self.plt.add_checkbox_button_widget(callback=self.plotPoints, position=(pos, 10.))

            if self.boxSize is not None:
                pos += 170.
                self.plt.add_text('Boxes', position=(pos, 65.), font_size=12)
                self.buttonBoxes = self.plt.add_checkbox_button_widget(callback=self.plotBoxes, position=(pos, 10.))

            # Picking Callbacks
            def removeSelection(selection):
                selected = selection['orig_extract_id']
                self.pv_points.remove_cells(selected, inplace=True)
                ids_removed = self.points_ids[selected]
                self.points_ids = np.delete(self.points_ids, selected)
                if self.normals is not None:
                    self.pv_normals = np.delete(self.pv_normals, selected, 0)
                    self.plt.remove_actor(self.normals_actor)
                    if self.buttonNormals.GetRepresentation().GetState():
                        self.normals_actor = self.plt.add_arrows(self.pv_points.cell_centers().points, self.pv_normals,
                                                                 mag=10, color='red', reset_camera=False)
                if self.buttonBoxes.GetRepresentation().GetState():
                    for point_id in ids_removed:
                        self.plt.remove_actor(self.box_actor[int(point_id - 1)])

            def enableRemoveSelection():
                self.plt.enable_cell_picking(mesh=self.pv_points,
                                             callback=removeSelection,
                                             font_size=12, opacity=0)

            # Picking Controls
            self.plt.main_menu.addAction('Point removal', enableRemoveSelection)

        if isinstance(self.normals, np.ndarray):
            pos += 170. if pos != 0 else 45.
            self.plt.add_text('Directions', position=(pos, 65.), font_size=12)
            self.buttonNormals = self.plt.add_checkbox_button_widget(callback=self.plotNormals, position=(pos, 10.))

    def getBinning(self, file):
        dim = ImageHandler().read(file).getDimensions()
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
        grid.cell_data["values"] = data.flatten(order="K")

        return grid

    def surfaceFromMRC(self, data, binarize=True, label=1):
        '''Function to convert an MRC file into an Structure Surface in VTK'''

        # Get Only mesh corresponding to a given label (smooth the result to fill holes in the mask)
        if binarize:
            data = data == label
            data = binary_erosion(binary_dilation(data, footprint=ball(4)), footprint=ball(1))

            # Triangulate coordinates using marching cubes algorithm
            grid = self.marchingCubes(data)

        else:
            grid = self.marchingCubes(data, level=label)

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
        faces = np.column_stack((3 * np.ones((len(faces), 1), dtype=int), faces)).flatten()
        grid = pv.PolyData(vertices.astype(int), faces) if triangulation else pv.PolyData(vertices.astype(int))
        return grid

    def histogram(self, volume):
        hist, edges = np.histogram(volume, bins=100)
        hist = hist / np.sum(hist)
        bin_centers = np.mean(np.vstack([edges[0:-1], edges[1:]]), axis=0)
        return hist, bin_centers

    def contours(self, hist, bin_centers):
        if self.discrete:
            logic_slicing = np.where(hist != 0)
            sliced_hist = hist[logic_slicing]
            logic_slicing_2 = np.where(sliced_hist != 0)
            contour_values = bin_centers[logic_slicing]
            contour_values = contour_values[logic_slicing_2]
            opacities = np.ones(contour_values.shape) * 0.3
        else:
            logic_slicing = np.where((hist > np.std(hist)) * (hist < np.amax(hist)))
            sliced_hist = hist[logic_slicing]
            logic_slicing_2 = np.where(sliced_hist < np.mean(sliced_hist))
            opacities = 1 - sliced_hist
            contour_values = bin_centers[logic_slicing]
            # print(1 - opacities)
            # print(np.mean(1-opacities))
            # print(contour_values)
            contour_values = contour_values[logic_slicing_2]
            # print(contour_values)
            opacities = opacities[logic_slicing_2]
        return contour_values, opacities

    def isovolumes(self, volume, range=0.01, sigma=None, triangulation=True):
        volume = volume if sigma is None else gaussian_filter(volume, sigma=sigma)
        hist, bin_centers = self.histogram(volume)
        contour_values, opacities = self.contours(hist, bin_centers)
        if not self.discrete:
            opacities = (range / (np.amax(opacities) - np.amin(opacities))) * (opacities - np.amin(opacities))
            return [self.marchingCubes(volume, level, triangulation) for level in contour_values], opacities
        else:
            return [self.surfaceFromMRC(volume, binarize=False, label=level) for level in contour_values], opacities

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
            # First check if tomogram is already loaded in memory, otherwise load it
            if isinstance(self.tomo, tuple):
                self.loading_screen.startAnimation()
                self.loadInMemory(source='tomo', sliceMode=False)
            else:
                self.showTomogram(load=False)
        else:
            for actor in self.tomo_actor:
                self.plt.remove_actor(actor)
            self.tomo_actor = []

    def toogleSlice(self, value):
        if value:
            # First check if tomogram is already loaded in memory, otherwise load it
            if isinstance(self.tomo, tuple):
                self.loading_screen.startAnimation()
                self.loadInMemory(source='tomo', sliceMode=True)
            else:
                self.showSlices(load=False)
        else:
            self.plt.remove_actor(self.tomo_slice_actor)
            self.plt.clear_plane_widgets()
            self.tomo_slice_actor = None

    def plotMasks(self, value):
        if value:
            # First check if mask is already loaded in memory, otherwise load it
            if isinstance(self.mask, str):
                self.loading_screen.startAnimation()
                self.loadInMemory(source='mask')
            else:
                self.showMask(load=False)
        else:
            for actor in self.mask_actors:
                self.plt.remove_actor(actor)
            self.graph_actor = []

    def plotPoints(self, value):
        if value:
            reset_camera = False
            if self.first_reset:
                self.first_reset = False
                reset_camera=True

            self.points_actor.append(self.plt.add_mesh(self.pv_points, show_scalar_bar=False, scalars="colors", categories=True,
                                           render_points_as_spheres=True, reset_camera=reset_camera))
        else:
            for actor in self.points_actor:
                self.plt.remove_actor(actor)
            self.points_actor = []

    def plotBoxes(self, value):
        if value:
            for point_id in self.points_ids:
                idp = int(point_id - 1)
                cube = pv.Cube(self.points[idp],
                               x_length=self.boxSize, y_length=self.boxSize, z_length=self.boxSize)
                self.box_actor[idp] = self.plt.add_mesh(cube, show_scalar_bar=False, style='wireframe',
                                                        color='blue')
        else:
            for actor in self.box_actor.values():
                self.plt.remove_actor(actor)
            self.box_actor = {}

    def plotNormals(self, value):
        if value:
            reset_camera = False
            if self.first_reset:
                self.first_reset = False
                reset_camera = True

            self.normals_actor = self.plt.add_arrows(self.pv_points.cell_centers().points, self.pv_normals,
                                                         mag=5, color='red', reset_camera=reset_camera)
        else:
            self.plt.remove_actor(self.normals_actor)
            self.normals_actor = None

    def initializePlot(self):
        # By default coordinates button is toggled
        if self.points is not None:
            self.buttonPoints.GetRepresentation().SetState(True)
            self.plotPoints(True)

        self.plt.app.exec_()

        # Save Points and Normals
        if self.save_basename is not None:
            np.savetxt(self.save_basename + '_indices.txt', self.points_ids, delimiter=' ')

    def loadInMemory(self, source, sliceMode=False):
        self.thread = QThread()
        self.worker = Worker(self)
        self.worker.moveToThread(self.thread)
        if source == 'tomo':
            self.thread.started.connect(self.worker.runTomoLoading)
            if sliceMode:
                self.thread.finished.connect(lambda: self.showSlices(load=True))
            else:
                self.thread.finished.connect(lambda: self.showTomogram(load=True))
        if source == 'mask':
            self.thread.started.connect(self.worker.runMaskLoading)
            self.thread.finished.connect(lambda: self.showMask(load=True))
        self.worker.finished.connect(self.thread.quit)
        # self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.start()

    def showTomogram(self, load=True):
        if load:
            self.tomo, self.pv_tomo, self.opacities, self.pv_tomo_slice = self.worker.output
            self.loading_screen.stopAnimation()
        if self.pv_tomo:
            if self.discrete:
                cmap = matplotlib.cm.get_cmap('Set3')
                cmap_ids = np.linspace(0, 1, len(self.pv_tomo))
                self.tomo_actor = [self.plt.add_mesh(actor, show_scalar_bar=False, opacity=op, color=cmap(cid),
                                                     render_points_as_spheres=True)
                                   for actor, op, cid in zip(self.pv_tomo, self.opacities, cmap_ids)]
            else:
                cmap = matplotlib.cm.get_cmap('bone')  # Greys also looks nice
                cmap_ids = np.linspace(0, 1, len(self.pv_tomo))
                self.tomo_actor = [self.plt.add_mesh(actor, show_scalar_bar=False, opacity=3 * op, color=cmap(cid),
                                                     render_points_as_spheres=True)
                                   for actor, op, cid in zip(self.pv_tomo, self.opacities, cmap_ids)]
            self.plt.reset_camera()
        else:
            self.buttonTomo.GetRepresentation().SetState(False)

    def showSlices(self, load=True):
        if load:
            self.tomo, self.pv_tomo, self.opacities, self.pv_tomo_slice = self.worker.output
            self.loading_screen.stopAnimation()
        self.tomo_slice_actor = self.plt.add_mesh_slice(self.pv_tomo_slice, normal='z', cmap="gray",
                                                        show_scalar_bar=False,
                                                        outline_translation=False, origin_translation=False)

    def showMask(self, load=True):
        if load:
            self.mask, self.pv_masks = self.worker.output
            self.loading_screen.stopAnimation()
        cmap = matplotlib.cm.get_cmap('Set3')
        cmap_ids = np.linspace(0, 1, len(self.pv_masks))
        self.mask_actors = [self.plt.add_mesh(mask, show_scalar_bar=False, color=cmap(cmap_id), smooth_shading=True)
                            for mask, cmap_id in zip(self.pv_masks, cmap_ids)]


class LoadingScreen(QWidget):
    def __init__(self):
        super().__init__()
        self.setFixedSize(200, 200)
        self.setWindowFlags(Qt.WindowStaysOnTopHint | Qt.CustomizeWindowHint)
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())

        radius = 20.0
        path = QtGui.QPainterPath()
        self.resize(440, 220)
        path.addRoundedRect(QRectF(self.rect()), radius, radius)
        mask = QtGui.QRegion(path.toFillPolygon().toPolygon())
        self.setMask(mask)

        self.label_animation = QLabel(self)
        plugin_path = os.path.dirname(tomoviz.__file__)
        self.movie = QMovie(os.path.join(plugin_path, "loading.gif"))
        self.label_animation.setMovie(self.movie)

    def startAnimation(self):
        self.movie.start()
        self.show()

    def stopAnimation(self):
        self.movie.stop()
        self.close()

class Worker(QObject):
    finished = pyqtSignal()

    def __init__(self, viewer):
        super().__init__()
        self.viewer = viewer

    def runTomoLoading(self):
        triangulation = self.viewer.tomo[1]
        sigma = self.viewer.tomo[2]
        tomo = self.viewer.readMRC(self.viewer.tomo[0], order=5, binning=self.viewer.binning)
        if isinstance(tomo, np.ndarray):
            pv_tomo_slice = self.viewer.gridFromMRC(tomo)
            try:
                pv_tomo, opacities = self.viewer.isovolumes(tomo, triangulation=triangulation, sigma=sigma)
            except:
                print("3D Tomogram view could not be computed, only slice mode will be available")
                pv_tomo = None
                opacities = None
        self.output = (tomo, pv_tomo, opacities, pv_tomo_slice)
        self.finished.emit()

    def runMaskLoading(self):
        mask = self.viewer.readMRC(self.viewer.mask, binning=self.viewer.binning)
        if isinstance(mask, np.ndarray):
            labels = np.unique(mask)[1:]
            pv_masks = [self.viewer.surfaceFromMRC(mask, label=label) for label in labels]
        self.output = (mask, pv_masks)
        self.finished.emit()

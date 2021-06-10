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
from pyvista.utilities import generate_plane
import pymeshfix as pm
import vtk
import numpy as np
import matplotlib

from pwem.emlib.image import ImageHandler

class MrcPlot(object):
    '''
    Class to visualize MRC files
    Input paramters:
         - tomo_mrc (Optional): File containing a Volume (in MRC format)
         - mask_mrc (Optional): File containing a Mask (in MRC format)
    Usage:
         import MRCPlot
         plt = MRCPlot(vti_file=path_vti, graph_file=path_graph, net_file=path_net, peaks_file=path_peaks)
         plt.initializePlot()
    '''

    def __init__(self, tomo_mrc=None, mask_mrc=None):
        self.tomo = self.readMRC(tomo_mrc) if tomo_mrc is not None else None
        self.mask = self.readMRC(mask_mrc) if mask_mrc is not None else None

        # Get Pyvista Objects
        if isinstance(self.tomo, np.ndarray):
            self.pv_tomo = self.gridFromMRC(self.tomo)
        if isinstance(self.mask, np.ndarray):
            labels = np.unique(self.mask)[1:]
            self.pv_masks = [self.surfaceFromMRC(self.mask, label=label) for label in labels]

        self.tomo_actor = None
        self.mask_actors = []

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

    def readMRC(self, file):
        return np.squeeze(ImageHandler().read(file).getData())

    def gridFromMRC(self, data):
        '''Function to convert an MRC file into an Structure Grid in VTK'''

        # Reorder Axis (This is needed as Numpy doesn't follow the Z,Y,X convention from Scipion but the X,Y,Z
        # axis convention
        data = np.swapaxes(np.swapaxes(data, 0, 2), 1, 0)

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

        # Reorder Axis (This is needed as Numpy doesn't follow the Z,Y,X convention from Scipion but the X,Y,Z
        # axis convention
        data = np.swapaxes(np.swapaxes(data, 0, 2), 1, 0)

        # Get Only mesh corresponding to a given label
        data = data * (data == label).astype(int)

        # Create a Grid from the data
        x = np.arange(0, data.shape[0])
        y = np.arange(0, data.shape[1])
        z = np.arange(0, data.shape[2])
        x, y, z = np.meshgrid(x, y, z)

        # Copy Data to Poly Data VTK object (Pyvista)
        grid = pv.StructuredGrid(x, y, z)

        # Set the cell values. Previous reordering of the axis was needed to flatten properly the array
        grid.point_arrays["values"] = data.flatten(order="K").astype(int)

        # Triangulate coordinates
        grid = grid.contour()
        grid = pv.PolyData(grid.points).delaunay_3d(alpha=1).triangulate().extract_geometry()

        # Fix the mesh
        mfix = pm._meshfix.PyTMesh(False)
        mfix.load_array(grid.points, grid.faces.reshape((-1, 4))[:, 1:])
        mfix.fill_small_boundaries(nbe=100, refine=True)
        vert, faces = mfix.return_arrays()

        # Reconstruct fixed mesh
        triangles = np.empty((faces.shape[0], 4))
        triangles[:, -3:] = faces
        triangles[:, 0] = 3
        grid = pv.PolyData(vert, triangles.astype(int))

        # Final fixing using pyvista and smoothing
        grid = grid.fill_holes(10)
        grid = grid.smooth(n_iter=500).clean()

        return grid

    def plotTomo(self, value):
        if value:
            self.tomo_actor = self.plt.add_mesh_slice(self.pv_tomo, normal='z', cmap="bone", show_scalar_bar=False,
                                                      outline_translation=False, origin_translation=False)
            self.buttonSliceTomo.GetRepresentation().SetState(True)
            self.plt.reset_camera()
        else:
            self.plt.remove_actor(self.tomo_actor)
            self.buttonSliceTomo.GetRepresentation().SetState(False)
            self.tomo_actor = None

    def toogleSlice(self, value):
        if value:
            if self.buttonTomo.GetRepresentation().GetState():
                plane_sliced_mesh = self.plt.plane_sliced_meshes[0]
                alg = vtk.vtkCutter()
                alg.SetInputDataObject(self.pv_tomo)

                def callback(normal, origin):
                    plane = generate_plane(normal, origin)
                    alg.SetCutFunction(plane)
                    alg.Update()
                    plane_sliced_mesh.shallow_copy(alg.GetOutput())

                self.plt.add_plane_widget(callback=callback, bounds=self.pv_tomo.bounds,
                                          factor=1.25, normal='z',
                                          origin_translation=False,
                                          outline_translation=False,
                                          origin=self.pv_tomo.center)
            else:
                self.buttonSliceTomo.GetRepresentation().SetState(False)
        else:
            self.plt.clear_plane_widgets()

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

    def initializePlot(self):
        self.plt.app.exec_()


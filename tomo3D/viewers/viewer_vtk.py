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
import vtk

class VtkPlot(object):
    '''
    Class to visualize VTK files
    Input paramters:
         - vti_file (Optional): File containing a Volume (Volume VTK Object)
         - graph_file (Optional): File containing a Graph (PolyData VTK Object)
         - net_file (Optional): File containing a Net (PolyData VTK Object)
         - peaks_file (Optional): File containing Peaks - Coordinates (PolyData VTK Object)
    Usage:
         import VtkPlot
         plt = VtkPlot(vti_file=path_vti, graph_file=path_graph, net_file=path_net)
         plt.initializePlot()
    '''

    def __init__(self, vti_file=None, graph_file=None, net_file=None, peaks_file=None):
        self.vti = pv.read(vti_file) if vti_file is not None else None
        self.graph = pv.read(graph_file) if graph_file is not None else None
        self.net = pv.read(net_file) if net_file is not None else None
        self.peaks = pv.read(peaks_file) if peaks_file is not None else None

        self.vti_actor = None
        self.graph_actor = None
        self.net_actor = None
        self.peaks_actor = None
        self.vectors_actor = None

        self.plt = pvqt.BackgroundPlotter()
        self.plt.main_menu.clear()

        pos = 0.

        def function_builder(vtk_obj, vtk_actor, vtk_button, vtk_actor_name, prop):
            '''
            Function to create in a automatic manner the contents of the different custom menus in
            pyvistaqt
            '''
            def function():
                self.plt.remove_actor(vtk_actor)
                setattr(self, vtk_actor_name,
                        self.plt.add_mesh(vtk_obj, show_scalar_bar=False, colormap="cool", scalars=prop))
                vtk_button.GetRepresentation().SetState(True)

            return function

        if self.vti:
            pos += 45.
            self.buttonVti = self.plt.add_checkbox_button_widget(callback=self.plotVti, position=(pos, 10.))
            self.plt.add_text('Tomogram', position=(pos, 65.), font_size=12)
            pos += 170.
            self.buttonSliceVti = self.plt.add_checkbox_button_widget(callback=self.toogleSlice, position=(pos, 10.))
            self.plt.add_text('Slice Mode', position=(pos, 65.), font_size=12)

        if self.graph:
            pos += 170. if pos != 0 else 45.
            self.plt.add_text('Graph', position=(pos, 65.), font_size=12)
            self.buttonGraph = self.plt.add_checkbox_button_widget(callback=self.plotGraph, position=(pos, 10.))

            # Color By Menu
            self.callbacks_graph = {}
            graph_properties = self.graph.array_names
            graph_menu = self.plt.main_menu.addMenu('Color Graph By')
            for prop in graph_properties:
                if 'normal' not in prop:
                    self.callbacks_graph[prop] = function_builder(self.graph, self.graph_actor,
                                                                  self.buttonGraph, 'graph_actor', prop)
                    graph_menu.addAction(prop, self.callbacks_graph[prop])

        if self.net:
            pos += 170. if pos != 0 else 45.
            self.plt.add_text('Net', position=(pos, 65.), font_size=12)
            self.buttonNet = self.plt.add_checkbox_button_widget(callback=self.plotNet, position=(pos, 10.))

            # Color By Menu
            self.callbacks_net = {}
            net_properties = self.net.array_names
            net_menu = self.plt.main_menu.addMenu('Color Net By')
            for prop in net_properties:
                if 'normal' not in prop:
                    self.callbacks_net[prop] = function_builder(self.net, self.net_actor,
                                                                self.buttonNet, 'net_actor', prop)
                    net_menu.addAction(prop, self.callbacks_net[prop])

        if self.peaks:
            self.peaks.set_active_vectors('smb_normal')
            pos += 170. if pos != 0 else 45.
            self.plt.add_text('Peaks', position=(pos, 65.), font_size=12)
            self.buttonPeaks = self.plt.add_checkbox_button_widget(callback=self.plotPeaks, position=(pos, 10.))
            pos += 170.
            self.buttonVectors = self.plt.add_checkbox_button_widget(callback=self.plotVectors, position=(pos, 10.))
            self.plt.add_text('Directions', position=(pos, 65.), font_size=12)

    def plotVti(self, value):
        if value:
            self.vti_actor = self.plt.add_mesh_slice(self.vti, normal='z', cmap="bone", show_scalar_bar=False,
                                                     outline_translation=False, origin_translation=False)
            self.buttonSliceVti.GetRepresentation().SetState(True)
            self.plt.reset_camera()
        else:
            self.plt.remove_actor(self.vti_actor)
            self.buttonSliceVti.GetRepresentation().SetState(False)
            self.vti_actor = None

    def toogleSlice(self, value):
        if value:
            if self.buttonVti.GetRepresentation().GetState():
                plane_sliced_mesh = self.plt.plane_sliced_meshes[0]
                alg = vtk.vtkCutter()
                alg.SetInputDataObject(self.vti)

                def callback(normal, origin):
                    plane = generate_plane(normal, origin)
                    alg.SetCutFunction(plane)
                    alg.Update()
                    plane_sliced_mesh.shallow_copy(alg.GetOutput())

                self.plt.add_plane_widget(callback=callback, bounds=self.vti.bounds,
                                          factor=1.25, normal='z',
                                          origin_translation=False,
                                          outline_translation=False,
                                          origin=self.vti.center)
            else:
                self.buttonSliceVti.GetRepresentation().SetState(False)
        else:
            self.plt.clear_plane_widgets()

    def plotGraph(self, value):
        if value:
            self.graph_actor = self.plt.add_mesh(self.graph, show_scalar_bar=False, colormap="cool")
        else:
            self.plt.remove_actor(self.graph_actor)
            self.graph_actor = None

    def plotNet(self, value):
        if value:
            self.net_actor = self.plt.add_mesh(self.net, show_scalar_bar=False, colormap="cool")
        else:
            self.plt.remove_actor(self.net_actor)
            self.net_actor = None

    def plotPeaks(self, value):
        mag = 0.02 * self.vti.dimensions[0]
        if value:
            self.peaks_actor = self.plt.add_points(self.peaks, color="orange", point_size=mag,
                                                   render_points_as_spheres=True, show_scalar_bar=False)
        else:
            self.plt.remove_actor(self.peaks_actor)
            self.peaks_actor = None

    def plotVectors(self, value):
        mag = 0.05 * self.vti.dimensions[0]
        if value:
            self.vectors_actor = self.plt.add_arrows(self.peaks.points, self.peaks.active_vectors, color="red",
                                                     mag=mag)
        else:
            self.plt.remove_actor(self.vectors_actor)
            self.vectors_actor = None

    def initializePlot(self):
        self.plt.show()

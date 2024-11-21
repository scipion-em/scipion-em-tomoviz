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
import numpy as np

class VtkPlot(object):
    '''
    Class to visualize VTK files
    Input paramters:
         - vti_file (Path (Str) - Optional): File containing a Volume (Volume VTK Object)
         - graph_file (Path (Str) - Optional): File containing a Graph (PolyData VTK Object)
         - net_file (Path (Str) - Optional): File containing a Net (PolyData VTK Object)
         - peaks_file (Path (Str) - Optional): File containing Peaks - Coordinates (PolyData VTK Object)
         - surf_file (Path (Str) - Optional): File containing a Surface (PolyData VTK Object)
    Usage:
         import VtkPlot
         plt = VtkPlot(vti_file=path_vti, graph_file=path_graph, net_file=path_net, peaks_file=path_peaks)
         plt.initializePlot()
    '''

    def __init__(self, vti_file=None, graph_file=None, net_file=None, peaks_file=None, surf_file=None):
        self.vti = pv.read(vti_file) if vti_file is not None else None
        self.graph = pv.read(graph_file) if graph_file is not None else None
        self.net = pv.read(net_file) if net_file is not None else None
        self.peaks = pv.read(peaks_file) if peaks_file is not None else None
        self.surf = pv.read(surf_file) if surf_file is not None else None

        self.vti_actor = None
        self.graph_actor = None
        self.net_actor = None
        self.net_slice_actor = None
        self.peaks_actor = None
        self.vectors_actor = None
        self.surf_actor = None

        self.plane_widgets = []
        self.cut_origin = None
        self.cut_normal = None

        self.plt = pvqt.BackgroundPlotter(title='Scipion tomoviz viewer', window_size=(1200, 800))
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
            self.buttonSliceVti = self.plt.add_checkbox_button_widget(callback=self.toogleTomoSlice, position=(pos, 10.))
            self.plt.add_text('Tomo Slice', position=(pos, 65.), font_size=12)

        if self.graph:
            pos += 170. if pos != 0 else 45.
            self.plt.add_text('Graph', position=(pos, 65.), font_size=12)
            self.buttonGraph = self.plt.add_checkbox_button_widget(callback=self.plotGraph, position=(pos, 10.))

            # Color By Menu
            self.callbacks_graph = {}
            graph_properties = np.sort(self.graph.array_names)
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
            pos += 170.
            self.buttonSliceNet = self.plt.add_checkbox_button_widget(callback=self.toogleNetSlice, position=(pos, 10.))
            self.plt.add_text('Net Slice', position=(pos, 65.), font_size=12)

            # Color By Menu
            self.callbacks_net = {}
            net_properties = np.sort(self.net.array_names)
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

        if self.surf:
            pos += 170. if pos != 0 else 45.
            self.plt.add_text('Surface', position=(pos, 65.), font_size=12)
            self.buttonPeaks = self.plt.add_checkbox_button_widget(callback=self.plotSurface, position=(pos, 10.))

    def plotVti(self, value):
        if value:
            self.vti_actor = self.plt.add_mesh_slice(self.vti, normal='z', cmap="bone", show_scalar_bar=False,
                                                     outline_translation=False, origin_translation=False)

            # In order to set our own callback to the widget, we need to turn off the current widget
            # And turning it on again with the custom function we want to use
            self.plt.plane_widgets[-1].Off()
            del self.plt.plane_widgets[-1]
            self.buttonSliceVti.GetRepresentation().SetState(True)
            self.toogleTomoSlice(True)

            self.plt.reset_camera()
        else:
            idx = self.plane_widgets.index(self.vti_actor)
            self.plt.plane_widgets[idx].Off()
            del self.plt.plane_widgets[idx]
            del self.plane_widgets[idx]
            self.plt.remove_actor(self.vti_actor)
            self.buttonSliceVti.GetRepresentation().SetState(False)
            self.vti_actor = None

            # We need to check if we have a Net Slice and remove if it exists
            if self.net and self.buttonSliceNet.GetRepresentation().GetState():
                self.plt.remove_actor(self.net_slice_actor)
                self.net_slice_actor = None
                self.buttonSliceNet.GetRepresentation().SetState(False)

            # We remove all the slice meshes and restore the plane widget variable of the viewer
            del self.plt.plane_sliced_meshes
            self.plane_widgets = []

    def toogleTomoSlice(self, value):
        if value:
            if self.buttonVti.GetRepresentation().GetState():
                if self.vti_actor in self.plane_widgets:
                    idx = self.plane_widgets.index(self.vti_actor)
                    self.plt.plane_widgets[idx].On()
                else:
                    # Since Tomo button will toggle this automatically, it is safe to append here the tomo
                    # actor to the currently active plane widgets
                    self.plane_widgets.append(self.vti_actor)

                    # Create all the vtkCutter objects that will be updated
                    alg_tomo, alg_net = vtk.vtkCutter(), vtk.vtkCutter()
                    alg_tomo.SetInputDataObject(self.vti)
                    alg_net.SetInputDataObject(self.net)

                    # This callback is called everytime the Tomo Slice (Cut) is updated
                    # by the user (move, rotation) or because it is newly created
                    def callback(normal, origin):
                        # Update cut normal and origin stored in the widget
                        self.cut_normal = normal
                        self.cut_origin = origin
                        # Create a new Plane based on new normal/origin
                        plane = generate_plane(normal, origin)
                        # Update Tomo Slice (Cut) position based on new plane
                        alg_tomo.SetCutFunction(plane)
                        alg_tomo.Update()
                        idx = self.plane_widgets.index(self.vti_actor)
                        plane_sliced_tomo = self.plt.plane_sliced_meshes[idx]
                        plane_sliced_tomo.shallow_copy(alg_tomo.GetOutput())
                        # This callback also takes into account the Net Plane
                        # In case it exists, it is update so it matches the current orientation of the
                        # Tomogram Slice
                        if self.net and self.buttonSliceNet.GetRepresentation().GetState():
                            alg_net.SetCutFunction(plane)
                            alg_net.Update()
                            idx = self.plane_widgets.index(self.net_slice_actor)
                            plane_sliced_net = self.plt.plane_sliced_meshes[idx]
                            plane_sliced_net.shallow_copy(alg_net.GetOutput())

                    self.plt.add_plane_widget(callback=callback, bounds=self.vti.bounds,
                                              factor=1.25, normal='z',
                                              origin_translation=False,
                                              outline_translation=False,
                                              origin=self.vti.center)
            else:
                self.buttonSliceVti.GetRepresentation().SetState(False)
        else:
            idx = self.plane_widgets.index(self.vti_actor)
            self.plt.plane_widgets[idx].Off()

    def plotGraph(self, value):
        if value:
            self.graph_actor = self.plt.add_mesh(self.graph, show_scalar_bar=False, colormap="cool")
        else:
            self.plt.remove_actor(self.graph_actor)
            self.graph_actor = None

    def plotNet(self, value):
        if value:
            self.net_actor = self.plt.add_mesh(self.net, show_scalar_bar=False, colormap="cool")
            if self.buttonSliceNet.GetRepresentation().GetState():
                self.buttonSliceNet.GetRepresentation().SetState(False)
                self.plt.remove_actor(self.net_slice_actor)
                self.net_slice_actor = None
        else:
            self.plt.remove_actor(self.net_actor)
            self.net_actor = None

    def toogleNetSlice(self, value):
        if value:
            if self.buttonVti.GetRepresentation().GetState():
                # We create a new Plane to perform a cut to Net
                alg = vtk.vtkCutter()
                alg.SetInputDataObject(self.net)
                plane_sliced_mesh = pv.wrap(alg.GetOutput())
                plane = generate_plane(self.cut_normal, self.cut_origin)
                alg.SetCutFunction(plane)
                alg.Update()
                plane_sliced_mesh.shallow_copy(alg.GetOutput())

                # Create the actor and add the plane to Pyvista plane_slice_meshes and plane_widgets to update
                # it in the future
                self.net_slice_actor = self.plt.add_mesh(plane_sliced_mesh, show_scalar_bar=False, colormap="cool")
                self.plt.plane_sliced_meshes.append(plane_sliced_mesh)
                self.plane_widgets.append(self.net_slice_actor)

                if self.buttonNet.GetRepresentation().GetState():
                    self.buttonNet.GetRepresentation().SetState(False)
                    self.plt.remove_actor(self.net_actor)
                    self.net_actor = None
            else:
                self.buttonSliceNet.GetRepresentation().SetState(False)
        else:
            self.plt.remove_actor(self.net_slice_actor)
            self.net_slice_actor = None
            # Since we only add the tomo before the Net Slice, we know that the last appended element to
            # the sliced meshes is going to be the Net Slice Mesh
            del self.plt.plane_sliced_meshes[-1]
            del self.plane_widgets[-1]


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

    def plotSurface(self, value):
        if value:
            self.surf_actor = self.plt.add_mesh(self.surf, opacity=0.6, color="white", show_scalar_bar=False)
        else:
            self.plt.remove_actor(self.surf_actor)
            self.surf_actor = None

    def initializePlot(self):
        self.plt.app.exec_()

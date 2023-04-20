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
import logging
logger = logging.getLogger(__name__)
import pyvista as pv
import numpy as np
from scipy.spatial.distance import pdist
from multiprocessing import Process

from PyQt5.QtWidgets import QApplication, QMainWindow, QDockWidget
from pyvistaqt.plotting import QtInteractor
from PyQt5.QtCore import Qt

class TriangulationPlot(object):
    '''
    Class to visualize triangulation(s) and/or point cloud(s) with their associated normals
    Input paramters:
         - meshes (Mandatory):    List containing the pyvista triangulations to be shown by the viewer
         - clouds (Optional):     List containing the original point clouds associated to the triangulations.
                                  This parameter is useful to compare the coordinates adjusted to the mesh and
                                  the original ones
         -extNormals (Optional):  List containing the original normals associated to each vertex in the
                                  Delaunay triangulation. In general, this parameter is used to compare the
                                  normals computed inside Scipion (Pyvista normals stored in the meshes objects)
                                  and those normals obtained with other software (PySeg, Dynamo...)
    Usage:
         import TriangulationPlot
         plt = TriangulationPlot([mesh], clouds=[cloud], extNormals=[extNormal])
         plt.initializePlot()
    '''


    def __init__(self, meshes, clouds=None, extNormals_List=None, extNormals_coords=None):
        self.clouds = clouds
        self.meshes = meshes
        self.extNormals_List = extNormals_List
        self.extNormals_coords = extNormals_coords
        self.initialize_PyQt()
        self.actor_meshes = []
        self.actor_normals = []
        self.actor_extNormals = []
        self.actor_points = []
        self.actor_ori_points = []
        self.buttonMesh = self.p.add_checkbox_button_widget(callback=self.plotMesh, position=(65., 10.),
                                                            color_on='white')
        self.p.add_text('Mesh', position=(10.+48, 65.), font_size=12)
        self.p.add_text('Normals', position=(200.+44, 65.), font_size=12)
        self.buttonNormals = self.p.add_checkbox_button_widget(callback=self.plotNormals, position=(265., 10.))
        self.p.add_text('Cloud', position=(420.+33, 65.), font_size=12)
        self.buttonPoints = self.p.add_checkbox_button_widget(callback=self.plotPoints, position=(465., 10.),
                                                              color_on='black')
        if self.clouds is not None:
            self.p.add_text('Original Cloud', position=(598., 65.), font_size=12)
            self.buttonOriPoints = self.p.add_checkbox_button_widget(callback=self.plotOriPoints, position=(665., 10.),
                                                                     color_on='white')
        if self.extNormals_List is not None:
            self.p.add_text('External Normals', position=(598.+198, 65.), font_size=12)
            self.buttonExtNormals = self.p.add_checkbox_button_widget(callback=self.plotExtNormals, position=(865., 10.),
                                                                      color_on='red')

        self.mainWindow.show()

    def initialize_PyQt(self):
        self.app = QApplication([])
        self.mainWindow = QMainWindow()
        self.mainWindow.default_style_sheet = self.mainWindow.styleSheet()
        self.mainWindow.save_commands = True  # stores commands internally when True
        self.mainWindow.hold = False
        self.mainWindow.load_dialog = None
        self.mainWindow.app = self.app
        self.mainWindow.resize(1000, 600)
        self.p = QtInteractor(self.mainWindow)
        self.p.app_window = self.mainWindow
        self.mainWindow.dock_vtk = QDockWidget('Mesh Viewer', self.mainWindow)
        self.mainWindow.dock_vtk.setWidget(self.p.interactor)
        self.mainWindow.dock_vtk.setMinimumSize(700, 600)
        self.mainWindow.addDockWidget(Qt.LeftDockWidgetArea, self.mainWindow.dock_vtk)
        self.p.add_toolbars()

    def plotMesh(self, value):
        if value:
            for mesh in self.meshes:
                self.actor_meshes.append(self.p.add_mesh(mesh, show_edges=True))
            self.buttonPoints.GetRepresentation().SetState(False)
            if self.buttonOriPoints.GetRepresentation().GetState:
                for actor in self.actor_ori_points:
                    self.p.remove_actor(actor)
                self.buttonOriPoints.GetRepresentation().SetState(False)
        else:
            for actor in self.actor_meshes:
                self.p.remove_actor(actor)

    def plotNormals(self, value):
        if value:
            for mesh in self.meshes:
                vecLength = np.amax(pdist(mesh.points))
                normals = vecLength * mesh.point_normals
                self.actor_normals.append(self.p.add_arrows(mesh.points, normals,
                                          mag=0.1, color='blue'))
        else:
            for actor in self.actor_normals:
                self.p.remove_actor(actor)

    def plotExtNormals(self, value):
        if value:
            for idn, normals in enumerate(self.extNormals_List):
                normals = pv.pyvista_ndarray(normals)
                vecLength = np.amax(pdist(self.meshes[idn].points))
                normals = vecLength * (normals / np.linalg.norm(normals, axis=1)[:, np.newaxis])
                if self.extNormals_coords is None:
                    areZero = np.where((self.meshes[idn].point_normals == (0, 0, 0)).all(axis=1))
                    normals[areZero] = np.array((0, 0, 0))
                    self.actor_extNormals.append(self.p.add_arrows(self.meshes[idn].points, normals,
                                                                   mag=0.1, color='red'))
                else:
                    extNormals_coords = pv.pyvista_ndarray(self.extNormals_coords)
                    self.actor_extNormals.append(self.p.add_arrows(extNormals_coords, normals,
                                                                   mag=0.1, color='red'))
        else:
            for actor in self.actor_extNormals:
                self.p.remove_actor(actor)

    def plotPoints(self, value):
        if value:
            for mesh in self.meshes:
                self.actor_points.append(self.p.add_points(mesh, color='black', point_size=10,
                                                           render_points_as_spheres=True))
            self.buttonMesh.GetRepresentation().SetState(False)
        else:
            for actor in self.actor_points:
                self.p.remove_actor(actor)

    def plotOriPoints(self, value):
        if value:
            for cloud in self.clouds:
                self.actor_ori_points.append(self.p.add_points(cloud, color='white', point_size=10,
                                                               render_points_as_spheres=True))
            if self.buttonMesh.GetRepresentation().GetState():
                for actor in self.actor_meshes:
                    self.p.remove_actor(actor)
                self.buttonMesh.GetRepresentation().SetState(False)
        else:
            for actor in self.actor_ori_points:
                self.p.remove_actor(actor)

    def initializePlot(self):
        self.app.exec_()
        
# The following functions are used to create automatically a "gui thread" to avoid the 'exec_' loop 
# from freezing the main thread

def guiThread(classObj, methodName, *args, **kwargs):
    '''
    Create a new process to prevent the exec_ loop of the GUI from blocking the main thread. In order to work,
    the GUI classes must be instantiated withing the process.
    '''
    proc = Process(target=instantiateClass, args=(classObj, methodName, *args,), kwargs=kwargs)
    proc.start()
    
def instantiateClass(classObj, methodName, *args, **kwargs):
    '''
    Create an instance of any class and call a visualization method (or any other method).
        - classObj (class): Class to be instantiated
        - methodName (string): Method from the class to be called after the instantiation
        - *args (list): extra argument needed to call the class method
        - **kwargs (dict): arguments to be passed to the contructor of the class
    '''
    try:
        instance = classObj(**kwargs)
        runMethod(instance, methodName, *args)
    except Exception as e:
        logger.info('Cannot create instance of class %s' % classObj, exc_info=e)
    
def runMethod(instance, methodName, *args):
    '''
    Execute a method from an instantiated class:
        - instance (obj): object instantiated from a given class
        - methodName (string): method belonging to instance to be called
        - *args (list): extra arguments needed by the method
    '''
    try:
        method = getattr(instance, methodName)
        method(*args)
    except AttributeError as e:
        logger.error(methodName + ' is not a member the class', exc_info=e)

# **************************************************************************
# *
# * Authors:     David Herreros (dherreros@cnb.csic.es)
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
import numpy as np
import time
import os
from pyworkflow.protocol import Protocol
from pyworkflow import utils as pwutils
from pyworkflow.gui.dialog import ToolbarListDialog
from pyworkflow.gui.tree import TreeProvider
from tomo.objects import SetOfTiltSeriesCoordinates

from .viewer_triangulations import TriangulationPlot, guiThread
from .viewer_mrc import MrcPlot
from ..utils import delaunayTriangulation

from tomo.utils import extractVesicles, initDictVesicles, normalFromMatrix
import tomo.constants as const

class Tomo3DTreeProvider(TreeProvider):
    """ Populate Tree from SetOfTomograms or SetOfTiltSeries"""

    def __init__(self, tomoList):
        TreeProvider.__init__(self)
        self.tomoList = tomoList

    def getColumns(self):
        return [('Item', 300), ("# coords", 100)]

    def getObjectInfo(self, tomo):

        tomogramName = tomo.getTsId()

        return {'key': tomogramName, 'parent': None,
                'text': tomogramName,
                'values':(tomo.count),
                'tags': ("done")}

    def getObjectPreview(self, obj):
        return (None, None)

    def getObjectActions(self, obj):
        return []

    def _getObjectList(self):
        """Retrieve the object list"""
        return self.tomoList

    def getObjects(self):
        objList = self._getObjectList()
        return objList

    def configureTags(self, tree):
        tree.tag_configure("done", foreground="black")


class Tomo3DDialog(ToolbarListDialog):
    """
    This class extend from ListDialog to allow calling
    a pyvista viewer subprocess from a list of Tomograms.
    """

    def __init__(self, parent, coordinates, extnormals, **kwargs):
        self.dictVesicles, _ = initDictVesicles(coordinates)
        self.coordinates = coordinates
        self.extnormals = extnormals
        self.provider = kwargs.get("provider", None)
        ToolbarListDialog.__init__(self, parent,
                                   "Tomogram List",
                                   allowsEmptySelection=False,
                                   itemDoubleClick=self.doubleClickOnTomogram,
                                   lockGui=False,
                                   cancelButton=True,
                                   **kwargs)

    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        # tomoName = pwutils.removeBaseExt(self.tomo.getFileName())
        tomoName = pwutils.removeBaseExt(self.tomo.get())
        self.dictVesicles = extractVesicles(self.coordinates, self.dictVesicles, tomoName)
        self.createViewer()

    def createViewer(self):
        # tomoName = pwutils.removeBaseExt(self.tomo.getFileName())
        tomoName = pwutils.removeBaseExt(self.tomo.get())
        vesicles = self.dictVesicles[tomoName]['vesicles']
        shells = []
        for vesicle in vesicles:
            shells.append(delaunayTriangulation(vesicle))

        if self.extnormals is not None:
            normals = []
            extcoords = []
            normals.append([])
            for subtomo in self.extnormals:  # separate normals by tomoName/volId
                if pwutils.removeBaseExt(subtomo.getVolName()) == tomoName:
                    normal = normalFromMatrix(subtomo.getTransform().getMatrix())
                    coord = subtomo.getCoordinate3D()
                    coordSubtomo = [coord.getX(const.SCIPION),
                                    coord.getY(const.SCIPION),
                                    coord.getZ(const.SCIPION)]
                    normals[0].append(normal)
                    extcoords.append(coordSubtomo)

        else:
            normals = self.dictVesicles[tomoName]['normals']
            extcoords = None

        classArgs = {'meshes': shells, 'clouds': vesicles, 'extNormals_List': normals, 'extNormals_coords': extcoords}
        guiThread(TriangulationPlot, 'initializePlot', **classArgs)

class ViewerMRCDialog(ToolbarListDialog):
    """
    This class extend from ListDialog to allow calling
    a pyvista viewer subprocess from a list of Tomograms.
    """

    def __init__(self, parent, coords, protocol:Protocol, **kwargs):
        self.coords = coords
        self.prot = protocol
        self.provider = kwargs.get("provider", None)
        ToolbarListDialog.__init__(self, parent,
                                   "Tomogram list",
                                   allowsEmptySelection=False,
                                   itemDoubleClick=self.doubleClickOnTomogram,
                                   lockGui=False,
                                   cancelButton=True,
                                   **kwargs)

    def doubleClickOnTomogram(self, tomo=None):
        tomo_path = tomo.getFileName()
        coord_list = []
        direction_list = []
        boxSize = 32

        def fromCoordinated3D():
            for coord in self.coords.iterCoordinates(volume=tomo):
                direction = normalFromMatrix(coord.getMatrix())
                position = [coord.getX(const.BOTTOM_LEFT_CORNER),
                            coord.getY(const.BOTTOM_LEFT_CORNER),
                            coord.getZ(const.BOTTOM_LEFT_CORNER),
                            coord.getObjId()]
                position.append(coord.getGroupId()) if coord.getGroupId() is not None else position.append(0)
                coord_list.append(position)
                direction_list.append(direction)
            boxSize = self.coords.getBoxSize()

        def fromTiltSeriesCoordinates():
            for coord in self.coords.iterItems(where="_tsId='%s'" % tomo.getTsId()):
                direction = normalFromMatrix(np.eye(4))
                position = [coord.getX(),
                            coord.getY(),
                            coord.getZ(),
                            coord.getObjId(),
                            0] # Group Id
                coord_list.append(position)
                direction_list.append(direction)

        if not isinstance(self.coords, SetOfTiltSeriesCoordinates):
            fromCoordinated3D()
        else:
            fromTiltSeriesCoordinates()

        posFile = self.prot.getPath('positions.txt')
        directionsFile= self.prot.getPath('directions.txt')
        np.savetxt(posFile, np.asarray(coord_list))
        np.savetxt(directionsFile, np.asarray(direction_list))
        viewer_args = {'tomo_mrc': tomo_path, 'points': posFile, 'normals': directionsFile,
                       'boxSize': boxSize}

        guiThread(MrcPlot, 'initializePlot', **viewer_args)

    def haveCoordinatesChanged(self, coordinateFiles):
        hasChanged = False
        if len(coordinateFiles):
            for coordFile, oldCreationTime in coordinateFiles.items():
                if os.path.isfile(coordFile):
                    newCreationTime = time.ctime(os.path.getctime(coordFile))
                    if newCreationTime != oldCreationTime:
                        hasChanged = True
                        break
        return hasChanged
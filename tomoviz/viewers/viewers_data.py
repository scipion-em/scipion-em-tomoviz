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
import os.path
import time
import numpy as np
import pyworkflow.viewer as pwviewer
from pyworkflow.gui.dialog import askYesNo, showInfo
from pyworkflow.utils.properties import Message
import pyworkflow.utils as pwutils

from pwem.protocols import EMProtocol
import pwem.viewers.views as vi
from .views_tkinter_tree import Tomo3DTreeProvider
from .views_tkinter_tree import ViewerMRCDialog

import tomo.objects
from ..protocols import XmippProtFilterbyNormal


class TomoVizDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    using pyvista
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        tomo.objects.SetOfCoordinates3D,
        XmippProtFilterbyNormal,
        tomo.objects.SetOfTiltSeriesCoordinates,
        tomo.objects.SetOfSubTomograms
    ]
    _name = 'Open with TomoViz'

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return vi.ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        itemField = tomo.objects.Tomogram.TS_ID_FIELD

        if issubclass(cls, tomo.objects.SetOfCoordinates3D):
            outputCoords = obj
            tomos = outputCoords.getPrecedents()
            groupAttribute = tomo.objects.Coordinate3D.TOMO_ID_ATTR

        elif issubclass(cls, tomo.objects.SetOfSubTomograms):
            outputCoords = obj.getCoordinates3D()
            if outputCoords is None:
                showInfo("3D coordinates missing",
                         "This set of subtomograms is not associated with a set of coordinates 3d and cannot be used in this coordinate viewer.",
                         parent=self.getParent().root)
                return []
            tomos = outputCoords.getPrecedents()

            # Keep subtomos as the main collection, since Subtomograms subsets do not affect the coordinates set
            outputCoords = obj
            groupAttribute = tomo.objects.SubTomogram.VOL_NAME_FIELD

        elif issubclass(cls, EMProtocol):
            outputCoords = obj.inputMeshes.get()
            tomos = outputCoords.getPrecedents()
            groupAttribute = "_volId"

        elif issubclass(cls, tomo.objects.SetOfTiltSeriesCoordinates):
            outputCoords = obj
            tomos = outputCoords.getSetOfTiltSeries()
            groupAttribute = "_tsId"
            itemField = groupAttribute

        volIds = outputCoords.aggregate(["MAX", "COUNT"], groupAttribute, [groupAttribute])
        volIds = [(d[groupAttribute], d["COUNT"]) for d in volIds]

        tomoList = []
        coordinateFiles = dict()

        for objId in volIds:
            tomogram = tomos[{itemField:objId[0]}].clone()
            tomogram.count = objId[1]
            tomoList.append(tomogram)
            # Check if _indices.txt files exists for the tomogram and get the
            # creation time
            basename = pwutils.removeBaseExt(tomogram.getFileName())
            indices_file = basename + '_indices.txt'
            if os.path.isfile(indices_file):
                coordinateFiles[indices_file] = time.ctime(
                    os.path.getctime(indices_file))

                # tomoList = [String(tomos[objId].getFileName()) for objId in volIds]
        tomoProvider = Tomo3DTreeProvider(tomoList)
        viewer = ViewerMRCDialog(self._tkRoot, outputCoords, self.protocol,
                                 provider=tomoProvider,
                                 allowSelect=False)

        # For now we only generate 3d coordinates
        if viewer.haveCoordinatesChanged(coordinateFiles) and issubclass(cls, tomo.objects.SetOfCoordinates3D):
            import tkinter as tk
            frame = tk.Frame()
            if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, frame):
                protocol = self.protocol
                suffix = protocol._getOutputSuffix(tomo.objects.SetOfCoordinates3D)
                updated_set = protocol._createSetOfCoordinates3D(tomos, suffix)
                updated_set.setName("Selected Coordinates")
                updated_set.setPrecedents(tomos)
                updated_set.setSamplingRate(tomos.getSamplingRate())
                updated_set.setBoxSize(outputCoords.getBoxSize())
                for item in tomoList:
                    basename = pwutils.removeBaseExt(item.getFileName())
                    indices_file = basename + '_indices.txt'
                    if os.path.isfile(indices_file):
                        indices = np.loadtxt(indices_file, delimiter=' ')
                        for index in indices:
                            updated_set.append(outputCoords[index].clone())
                        pwutils.cleanPath(indices_file)
                name = protocol.OUTPUT_PREFIX + suffix
                args = {}
                args[name] = updated_set
                protocol._defineOutputs(**args)
                protocol._defineSourceRelation(tomos, updated_set)
                protocol._updateOutputSet(name, updated_set, state=updated_set.STREAM_CLOSED)

        # TODO: This should return a list of views and not already launching the dialog.
        return []

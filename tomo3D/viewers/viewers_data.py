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


import pyworkflow.viewer as pwviewer
from pyworkflow.object import String

from pwem.protocols import EMProtocol
import pwem.viewers.views as vi
from .views_tkinter_tree import Tomo3DTreeProvider
from .views_tkinter_tree import Tomo3DDialog

import tomo.objects
from ..protocols import XmippProtFilterbyNormal


class Tomo3DDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    using pyvista
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        tomo.objects.SetOfCoordinates3D,
        XmippProtFilterbyNormal
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return vi.ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        if issubclass(cls, tomo.objects.SetOfCoordinates3D):
            outputCoords = obj
        elif issubclass(cls, EMProtocol):
            outputCoords = obj.meshCoords

        tomos = outputCoords.getPrecedents()

        volIds = outputCoords.aggregate(["MAX"], "_volId", ["_volId"])
        volIds = [d['_volId'] for d in volIds]

        # tomoList = [tomos[objId].clone() for objId in volIds]
        tomoList = [String(tomos[objId].getFileName()) for objId in volIds]
        tomoProvider = Tomo3DTreeProvider(tomoList)

        if issubclass(cls, tomo.objects.SetOfCoordinates3D):
            Tomo3DDialog(self._tkRoot, outputCoords, None, provider=tomoProvider)
        elif issubclass(cls, EMProtocol):
            Tomo3DDialog(self._tkRoot, outputCoords, obj.outputset, provider=tomoProvider)

        return views

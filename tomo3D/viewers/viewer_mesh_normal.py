# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
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
from pyworkflow.protocol.params import EnumParam
import pwem.viewers.views as vi
import pwem.viewers as viewers
from .views_tkinter_tree import Tomo3DTreeProvider
from .views_tkinter_tree import Tomo3DDialog
from ..protocols import XmippProtFilterbyNormal

import tomo.objects

VOLUME_SLICES = 1
VOLUME_TOMO3D = 0

class ProtMeshNormal3DDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    using pyvista
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [XmippProtFilterbyNormal]

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('displayVol', EnumParam,
                      choices=['tomo3D', 'slices'], default=VOLUME_TOMO3D,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Display volume with',
                      help='*tomo3D*: display subtomograms with pyvista.\n*slices*: display subtomograms as 2D slices '
                           'along z axis.\n')

    def _getVisualizeDict(self):
        return {
            'displayVol': self._showVolumes,
        }


    # =========================================================================
    # Show Volumes
    # =========================================================================

    def _showVolumes(self, paramName=None):
        if self.displayVol == VOLUME_TOMO3D:
            return self._showVolumesTomo3D()
        elif self.displayVol == VOLUME_SLICES:
            return self._showVolumesXmipp()

    def _showVolumesTomo3D(self):
        outputSubtomos = self.protocol.outputset.get()
        tomos = self.protocol.inputTomos.get()

        volIds = outputSubtomos.aggregate(["MAX"], "_volId", ["_volId"])
        volIds = [d['_volId'] for d in volIds]

        tomoList = [tomos[objId].clone() for objId in volIds]
        tomoProvider = Tomo3DTreeProvider(tomoList)

        return Tomo3DDialog(self._tkRoot, outputSubtomos, provider=tomoProvider)

    def _showVolumesXmipp(self):
        # Write an sqlite with all tomograms selected for visualization.
        sampling, setOfObjects = self._getSetAndSampling()

        viewParams= {viewers.showj.MODE: viewers.showj.MODE_MD}
        view = viewers.views.ObjectView(self._project, setOfObjects.strId(),
                                        setOfObjects.getFileName(), viewParams=viewParams)
        view.setMemory(viewers.showj.getJvmMaxMemory() + 2)

        return [view]

    def _getSetAndSampling(self):
        setOfObjects = self._getOutput()
        sampling = setOfObjects.getSamplingRate()
        return sampling, setOfObjects


    # def __init__(self, **kwargs):
    #     pwviewer.Viewer.__init__(self, **kwargs)
    #     self._views = []
    #
    # def _getObjView(self, obj, fn, viewParams={}):
    #     return vi.ObjectView(
    #         self._project, obj.strId(), fn, viewParams=viewParams)
    #
    # def _visualize(self, obj, **kwargs):
    #     views = []
    #     cls = type(obj)
    #
    #     if issubclass(cls, tomo.objects.SetOfSubTomograms):
    #         outputSubtomos = obj
    #
    #         # tomos = outputSubtomos.getPrecedents()
    #         tomos = self.protocol.inputTomos
    #
    #         volIds = outputSubtomos.aggregate(["MAX"], "_volId", ["_volId"])
    #         volIds = [d['_volId'] for d in volIds]
    #
    #         tomoList = [tomos[objId].clone() for objId in volIds]
    #         tomoProvider = Tomo3DTreeProvider(tomoList)
    #
    #         Tomo3DDialog(self._tkRoot, outputSubtomos, provider=tomoProvider)
    #
    #     return views
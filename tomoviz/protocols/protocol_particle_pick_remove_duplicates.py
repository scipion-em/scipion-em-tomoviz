# **************************************************************************
# *
# * Authors:    Carlos Oscar Sorzano (coss@cnb.csic.es)
# *             Tomas Majtner (tmajtner@cnb.csic.es)  -- streaming version
# *             David Herreros Calero (dherreros@cnb.csic.es) -- Tomo version
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
# *  e-mail address 'coss@cnb.csic.es'
# *
# **************************************************************************
"""
Consensus picking protocol
"""

import numpy as np

from pyworkflow import BETA
import pyworkflow.protocol.params as params
from pyworkflow.utils import getFiles, removeBaseExt, moveFile, removeExt
from pwem.objects import Set

from .protocol_particle_pick_consensus import (ProtTomoConsensusPicking,
                                               consensusWorker, getReadyTomos)

import tomo.constants as const
from tomo.objects import SetOfCoordinates3D, Coordinate3D


class ProtTomoPickingRemoveDuplicates(ProtTomoConsensusPicking):
    """
    This protocol removes coordinates that are closer than a given threshold.
    The remaining coordinate is the average of the previous ones.
    """

    _label = 'remove duplicates'
    _devStatus = BETA
    outputName = 'outputCoordinates'
    FN_PREFIX = 'purgedCoords_'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam,
                      pointerClass='SetOfCoordinates3D',
                      label="Input coordinates", important=True,
                      help='Select the set of 3D coordinates to compare')
        form.addParam('consensusRadius', params.IntParam, default=10,
                      label="Radius",
                      help="All coordinates within this radius (in pixels) "
                           "are presumed to correspond to the same particle")

        # FIXME: It's not using more than one since
        #         self.stepsExecutionMode = STEPS_SERIAL
        # form.addParallelSection(threads=4, mpi=0)

#--------------------------- INSERT steps functions ----------------------------
    # --------------------------- INSERT steps functions ---------------------------

    def _processStep(self, tomoName):
        self.removeDuplicatesStep(tomoName)

    def getMainInput(self):
        return self.inputCoordinates.get()

    def defineRelations(self, outputSet):
        self._defineTransformRelation(self.getMainInput(), outputSet)

    def removeDuplicatesStep(self, tomoName):
        tomogram = self.tomograms[tomoName]
        tomoId = tomogram.getTsId(),
        print("Removing duplicates for tomogram %s: '%s'"
              % (tomoId, tomoName))

        coords = []
        vesicles = []
        transformations = []
        coordArray = np.asarray([x.getPosition(const.SCIPION) for x in
                                 self.getMainInput().iterCoordinates(tomogram)],
                                 dtype=int)
        vIds = np.asarray([x.getGroupId() for x in
                           self.getMainInput().iterCoordinates(tomogram)],
                           dtype=int)
        trMat = np.asarray([x.getMatrix() for x in
                            self.getMainInput().iterCoordinates(tomogram)])

        coords.append(np.asarray(coordArray, dtype=int))
        vesicles.append(np.asarray(vIds, dtype=int))
        transformations.append(trMat)

        fnTmp = self._getTmpPath('%s%s.txt' % (self.FN_PREFIX, tomoName))
        consensusWorker(coords, vesicles, transformations, 1,
                                         self.consensusRadius.get(),  fnTmp)

        self.generateOutput(tomogram, fnTmp)

    def _validate(self):
        errors = []
        return errors

    def _summary(self):
        return ["Radius = %d" % self.consensusRadius]

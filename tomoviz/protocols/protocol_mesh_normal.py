# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
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

from os import path
import numpy as np
from pyworkflow import BETA
from pyworkflow.protocol.params import PointerParam, FloatParam, BooleanParam, IntParam
import pwem.convert.transformations as tfs
from pwem.protocols import EMProtocol
from tomo.objects import SubTomogram, Coordinate3D
from tomo.protocols import ProtTomoBase
from tomo.utils import normalFromMatrix
import tomo.constants as const
from ..utils import delaunayTriangulation, computeNormals


class XmippProtFilterbyNormal(EMProtocol, ProtTomoBase):
    """ This protocol takes surfaces or ROIs (SetOfMeshes) and a SetOfSubtomograms or SetOfCoordinates3D with
    transformation matrix and filters them by different criteria related with the normal direction."""

    _label = 'filter by normal'
    _devStatus = BETA

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input', PointerParam, pointerClass="SetOfSubTomograms, SetOfCoordinates3D",
                      label='Subtomograms/Coordinates', help='SetOfSubtomograms to filter.')
        form.addParam('inputMeshes', PointerParam, label="Vesicles", pointerClass='SetOfMeshes',
                      help='Select the vesicles in which the subtomograms/coordinates are.')
        form.addParam('tilt', BooleanParam, default=False,
                      label='Filter by tilt angle',
                      help='Remove items depending on their tilt angle.')
        form.addParam('maxtilt', IntParam, default=150, label='Maximum allowed tilt', condition='tilt',
                      help='Remove the items that have a tilt angle bigger than the one specified in here, '
                           'considering tilt angle between 0 and 180 degrees.')
        form.addParam('mintilt', IntParam, default=30, label='Minimum allowed tilt', condition='tilt',
                      help='Remove the items that have a tilt angle smaller than the one specified in here, '
                           'considering tilt angle between 0 and 180 degrees.')
        form.addParam('normalDir', BooleanParam, default=True,
                      label='Filter by normal',
                      help='Remove the items that have a normal direction not equal to the normal direction of '
                           'the vesicle in the coordinate of the particle.')
        form.addParam('tol', FloatParam, default=5, label='Tolerance in degrees',
                      condition='normalDir',
                      help='Tolerance (in degrees) when comparing between subtomogram/coordinate and mesh normal '
                           'directions.')

        # form.addParam('topBottom', BooleanParam, default=True,
        #               label='Remove subtomograms in the top and bottom of the vesicle',
        #               help='Remove the subtomograms that have been picked from the top and bottom parts because they '
        #                    'had a different view.')
        # form.addParam('mwDir', BooleanParam, default=True,
        #               label='Remove subtomograms in the missing wedge direction',
        #               help='Remove the subtomograms that are in the missing wedge direction because they are highly '
        #                    'affected by the missing wedge.')
        # form.addParam('mwDir', EnumParam, default=True, label='Missing wedge direction',
        #               help='Missing wedge direction of the tomograms.')

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeNormalStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -------------------------------
    def computeNormalStep(self):
        inSet = self.input.get()
        tiltBool = self.tilt.get()
        normalBool = self.normalDir.get()
        inMeshes = self.inputMeshes.get()
        if normalBool:
            tol = self.tol.get() * np.pi / 180

            groupIdList = []
            tomoNameList = []
            for meshPoint in inMeshes:
                groupId = meshPoint.getGroupId()
                if groupId not in groupIdList:
                    groupIdList.append(groupId)
                    tomoNameList.append(str(meshPoint.getVolumeName()))

            meshDict = {key: {'tomoName': [tomoName], 'points': [], 'normals': []}
                        for key, tomoName in zip(groupIdList, tomoNameList)}

            for meshPoint in inMeshes.iterCoordinates():
                if meshDict[meshPoint.getGroupId()]["tomoName"][0] == meshPoint.getVolumeName():
                    meshDict[meshPoint.getGroupId()]["points"].append(meshPoint.getPosition(const.SCIPION))

            for i in meshDict:
                meshDict[i]["normals"] = self._getNormalVesicleList(np.asarray(meshDict[i]["points"]))

        if self._getInputisSubtomo(inSet.getFirstItem()):
            self.outSet = self._createSetOfSubTomograms()
        else:
            self.outSet = self._createSetOfCoordinates3D(inMeshes.getPrecedents())
        self.outSet.copyInfo(inSet)

        if tiltBool:
            if self._getInputisSubtomo(inSet.getFirstItem()):
                for item in inSet:
                    tilt = self._getTilt(item)
                    if self.maxtilt.get() > tilt > self.mintilt.get():
                        if normalBool:
                            self._filterByNormal(item, tol, meshDict)
                        else:
                            self.outSet.append(item)
            else:
                for item in inSet.iterCoordinates(volume=None):
                    tilt = self._getTilt(item)
                    if self.maxtilt.get() > tilt > self.mintilt.get():
                        if normalBool:
                            self._filterByNormal(item, tol, meshDict)
                        else:
                            self.outSet.append(item)

        if normalBool and not tiltBool:
            if self._getInputisSubtomo(inSet.getFirstItem()):
                for item in inSet:
                    self._filterByNormal(item, tol, meshDict)
            else:
                for item in inSet.iterCoordinates(volume=None):
                    self._filterByNormal(item, tol, meshDict)

    def createOutputStep(self):
        self._defineOutputs(outputset=self.outSet)
        self._defineSourceRelation(self.input.get(), self.outSet)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        if not self.normalDir.get() and not self.tilt.get():
            validateMsgs.append('Some filter should be switched to "Yes"')
        if not self.input.get().getFirstItem().hasTransform():
            validateMsgs.append('Input coordinates/subtomograms should have transform matrix.')
        return validateMsgs

    def _summary(self):
        summary = []
        if not self.isFinished():
            summary.append("Output subtomograms/coordinates not ready yet.")
        else:
            if self.tilt:
                summary.append("Remove input items by tilt angle (max allowed tilt: %d, min allowed tilt: %d)" %
                               (self.maxtilt.get(), self.mintilt.get()))
            if self.normalDir:
                summary.append("Remove input items by normal direction (tolerance of %0.2f degrees)" %
                               self.tol.get())
        return summary

    def _methods(self):
        methods = []
        if not self.isFinished():
            methods.append("Output subtomograms/coordinates not ready yet.")
        else:
            if self.tilt:
                methods.append("Input items with tilt angle bigger than %d or smaller than %d have been removed." %
                               (self.maxtilt.get(), self.mintilt.get()))
            if self.normalDir:
                methods.append("Input items that are not perpendicular to the membrane have been removed.")

            # if self.topBottom:
            #     methods.append("Particles in the top and bottom parts of the vesicles have been removed.")
            # if self.mwDir:
            #     methods.append("Particles in the missing wedge direction have been removed.")
        return methods

    # --------------------------- UTILS functions --------------------------------------------
    def _getInputisSubtomo(self, item):
        if isinstance(item, SubTomogram):
            return True
        else:
            return False

    def _getVesicleId(self, item):
        if isinstance(item, SubTomogram):
            c = item.getCoordinate3D()
        else:
            c = item
        if c.hasGroupId():
            vesicleId = c.getGroupId()
        else:  # For now it works with several vesicles in the same tomo just for input items with groupId
            vesicleId = 1
        return vesicleId

    def _getNormalVesicleList(self, points):
        triangulation = delaunayTriangulation(points)
        normalsList = computeNormals(triangulation, associateCoords=True)
        return normalsList

    def _getNormalVesicle(self, normalsList, item):
        if self._getInputisSubtomo(item):
            normSubtomo = normalFromMatrix(item.getTransform().getMatrix())
            coord = item.getCoordinate3D()
        else:
            normSubtomo = normalFromMatrix(item.getMatrix())
            coord = item
        coors = np.asarray([coord.getX(const.SCIPION),
                            coord.getY(const.SCIPION),
                            coord.getZ(const.SCIPION)])
        points, normals = zip(*normalsList)
        points = np.asarray(points)
        idx = np.argmin(np.sum((points - coors) ** 2, axis=1))
        return normSubtomo, normals[idx]

    def _getTilt(self, item):
        if self._getInputisSubtomo(item):
            _, tilt, _ = tfs.euler_from_matrix(item.getTransform().getMatrix(), axes='szyz')
        else:
            _, tilt, _ = tfs.euler_from_matrix(item.getMatrix(), axes='szyz')
        tilt = -np.rad2deg(tilt)
        return tilt

    def _filterByNormal(self, item, tol, meshDict):
        meshfromDict = meshDict[self._getVesicleId(item)]
        if meshfromDict["tomoName"][0] == path.basename(item.getVolName()):
            normSubtomo, normVesicle = self._getNormalVesicle(meshfromDict["normals"], item)
            if abs(normSubtomo[0] - normVesicle[0]) < tol and abs(normSubtomo[1] - normVesicle[1]) < tol and \
                    abs(normSubtomo[2] - normVesicle[2]) < tol:
                self.outSet.append(item)

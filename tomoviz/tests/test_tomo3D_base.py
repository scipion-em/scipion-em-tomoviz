# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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

import pyworkflow.utils as pwutils
import pwem.protocols as emprot
from pyworkflow.tests import BaseTest, setupTestProject

from . import DataSet
import tomoviz.protocols
import tomo.protocols


# class TestTomoPickingConsenus(BaseTest):
#     """This class check if the consensus from a series of SetOfCoordinates3D works properly."""
#
#     @classmethod
#     def setUpClass(cls):
#         setupTestProject(cls)
#         cls.dataset = DataSet.getDataSet('reliontomo')
#         cls.tomogram_1 = cls.dataset.getFile('tomo1')
#         cls.tomogram_2 = cls.dataset.getFile('tomo2')
#
#     def _importSetOfCoordinates(self, tomoFile, pattern):
#         protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
#                                               filesPath=tomoFile,
#                                               samplingRate=1,
#                                               objLabel=pwutils.removeBaseExt(tomoFile))
#
#         self.launchProtocol(protImportTomogram)
#
#         output = getattr(protImportTomogram, 'outputTomograms', None)
#
#         self.assertIsNotNone(output,
#                              "There was a problem with tomogram output")
#
#         protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
#                                   objLabel='Coordinates ' + pwutils.removeBaseExt(tomoFile),
#                                   auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_AUTO,
#                                   filesPath=self.dataset.getPath(),
#                                   importTomograms=protImportTomogram.outputTomograms,
#                                   filesPattern=pattern, boxSize=32,
#                                   samplingRate=1)
#         self.launchProtocol(protImportCoordinates3d)
#
#         return getattr(protImportCoordinates3d, 'outputCoordinates', None)
#
#     def _runTomoPickingConsensus(self, coordinates, consensus, mode, radius=100):
#         choices_mode = ['>=', '=']
#         choices_consensus = {-1: 'AND', 1: 'OR'}
#         protPickingConsensus = self.newProtocol(tomoviz.protocols.ProtTomoConsensusPicking,
#                                                 inputCoordinates=coordinates,
#                                                 consensus=consensus,
#                                                 mode=mode,
#                                                 consensusRadius=radius,
#                                                 objLabel='Consenus - ' + choices_consensus[consensus]
#                                                           + ' - ' + choices_mode[mode])
#         self.launchProtocol(protPickingConsensus)
#         return protPickingConsensus
#
#     def test_picking_consenus(self):
#         coords = []
#         coords.append(self._importSetOfCoordinates(self.tomogram_1, '*.tbl'))
#         coords.append(self._importSetOfCoordinates(self.tomogram_2, '*.tbl'))
#
#         # AND + >=
#         protConsensus = self._runTomoPickingConsensus(coords, -1, 0)
#         output = getattr(protConsensus, 'consensusCoordinates', None)
#         self.assertTrue(output,
#                              "There was a problem with consenus output")
#         self.assertTrue(output.getSize() == 20)
#         self.assertTrue(output.getBoxSize() == 32)
#         self.assertTrue(output.getSamplingRate() == 1)
#
#         # AND + =
#         protConsensus = self._runTomoPickingConsensus(coords, -1, 1, radius=50)
#         output = getattr(protConsensus, 'consensusCoordinates', None)
#         self.assertTrue(output,
#                              "There was a problem with consenus output")
#         self.assertTrue(output.getSize() == 5)
#         self.assertTrue(output.getBoxSize() == 32)
#         self.assertTrue(output.getSamplingRate() == 1)
#
#         # OR + >=
#         protConsensus = self._runTomoPickingConsensus(coords, 1, 0)
#         output = getattr(protConsensus, 'consensusCoordinates', None)
#         self.assertTrue(output,
#                              "There was a problem with consenus output")
#         self.assertTrue(output.getSize() == 83)
#         self.assertTrue(output.getBoxSize() == 32)
#         self.assertTrue(output.getSamplingRate() == 1)
#
#         # OR + =
#         protConsensus = self._runTomoPickingConsensus(coords, 1, 1)
#         output = getattr(protConsensus, 'consensusCoordinates', None)
#         self.assertTrue(output,
#                              "There was a problem with consenus output")
#         self.assertTrue(output.getSize() == 63)
#         self.assertTrue(output.getBoxSize() == 32)
#         self.assertTrue(output.getSamplingRate() == 1)
#
#         return output
#
# class TestTomoPickingRemoveDuplicates(BaseTest):
#     """This class check if the removal of duplicates from a series
#      of SetOfCoordinates3D works properly."""
#
#     @classmethod
#     def setUpClass(cls):
#         setupTestProject(cls)
#         cls.dataset = DataSet.getDataSet('reliontomo')
#         cls.tomogram_1 = cls.dataset.getFile('tomo1')
#         cls.tomogram_2 = cls.dataset.getFile('tomo2')
#
#     def _importSetOfCoordinates(self, tomoFile, pattern):
#         protImportTomogram = self.newProtocol(tomo.protocols.ProtImportTomograms,
#                                               filesPath=tomoFile,
#                                               samplingRate=5,
#                                               objLabel=pwutils.removeBaseExt(tomoFile))
#
#         self.launchProtocol(protImportTomogram)
#
#         output = getattr(protImportTomogram, 'outputTomograms', None)
#
#         self.assertIsNotNone(output,
#                              "There was a problem with tomogram output")
#
#         protImportCoordinates3d = self.newProtocol(tomo.protocols.ProtImportCoordinates3D,
#                                   objLabel='Coordinates ' + pwutils.removeBaseExt(tomoFile),
#                                   auto=tomo.protocols.ProtImportCoordinates3D.IMPORT_FROM_AUTO,
#                                   filesPath=self.dataset.getPath(),
#                                   importTomograms=protImportTomogram.outputTomograms,
#                                   filesPattern=pattern, boxSize=32,
#                                   samplingRate=5)
#         self.launchProtocol(protImportCoordinates3d)
#
#         return getattr(protImportCoordinates3d, 'outputCoordinates', None)
#
#     def _joinSets(self, sets):
#         protJoinSets = self.newProtocol(emprot.ProtUnionSet,
#                                         inputSets=sets,
#                                         objLabel='Duplicated Coordinates 3D')
#         self.launchProtocol(protJoinSets)
#         return getattr(protJoinSets, 'outputSet', None)
#
#     def _runTomoPickingRemoveDuplicates(self, coordinates):
#         protPickingConsensus = self.newProtocol(tomoviz.protocols.ProtTomoPickingRemoveDuplicates,
#                                                 inputCoordinates=coordinates,
#                                                 objLabel='Remove Duplicates')
#         self.launchProtocol(protPickingConsensus)
#         return protPickingConsensus
#
#     def test_picking_remove_duplicates(self):
#         coords = []
#         coords.append(self._importSetOfCoordinates(self.tomogram_1, '*.tbl'))
#         coords.append(self._importSetOfCoordinates(self.tomogram_1, '*.tbl'))
#         all_coords = self._joinSets(coords)
#         protRemoveDuplicates = self._runTomoPickingRemoveDuplicates(all_coords)
#         output = getattr(protRemoveDuplicates, 'outputCoordinates', None)
#         self.assertTrue(output,
#                              "There was a problem with consenus output")
#         self.assertTrue(output.getSize() == 63)
#         self.assertTrue(output.getBoxSize() == 32)
#         self.assertTrue(output.getSamplingRate() == 5)
#
#         return output


if __name__ == 'main':
    pass

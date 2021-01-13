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

import sys
import pyvista as pv
import numpy as np
from scipy.spatial import cKDTree


def rotation_matrix_from_vectors(vec1, vec2):
    # a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    a, b = vec1, vec2
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    if s != 0:
        tr = np.zeros([4, 4])
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
        tr[0:3, 0:3] = rotation_matrix
        tr[-1, -1] = 1
        return tr
    else:
        return np.eye(4)

def inverse_Transformation(transformation):
    translation = np.eye(4)
    rotation = np.eye(4)
    translation[0:-1, -1] = -transformation[0:-1, -1]
    rotation[:3, :3] = transformation[:3, :3]
    return np.linalg.inv(rotation).dot(translation)

def delaunayTriangulation(cloud, adjustCloud=True):
    cloud_pv = pv.PolyData(cloud)
    cloud_pv.GlobalWarningDisplayOff()
    mesh = cloud_pv.delaunay_3d()
    shell = mesh.extract_geometry().triangulate()

    # If points lie outside the shell, they are readjusted to compute their normals
    if adjustCloud:
        shell.subdivide(4, inplace=True)
        areZero = np.where((shell.point_normals == (0, 0, 0)).all(axis=1))
        points = np.asarray(shell.points)
        # points[areZero] = np.array((0, 0, 0))
        points = np.delete(points, areZero, axis=0)
        tree = cKDTree(points, leafsize=100)
        _, idc = tree.query(cloud)
        newCoords = points[idc]
        # newCoords = np.asarray([points[np.argmin(np.linalg.norm(points - point, axis=1))]
        #                         for point in cloud])
        shell = delaunayTriangulation(newCoords, adjustCloud=False)
    return shell

def computeNormals(triangulation, associateCoords=False):
    triangulation.compute_normals(inplace=True)
    normals = triangulation.point_normals
    points = triangulation.points

    # Check if coordinates are repeated (due to neighbour search) and copy the normal
    # so it is not (0,0,0)
    _, unique_indices, unique_inverse = np.unique(points, return_index=True, return_inverse=True, axis=0)
    unique_normals = normals[unique_indices]
    normals = unique_normals[unique_inverse]

    # NOT USED
    # for i in range(len(points)):  # generate pairs
    #     for j in range(i + 1, len(points)):
    #         if np.array_equal(points[i], points[j]):  # compare rows
    #             normals[j] = normals[i]
    #         else:
    #             pass

    # Sometimes, points may be redundant to the mesh
    # Assign the closest normal to them
    # Poosible cKDTree?
    areZero = np.where((normals == (0, 0, 0)).all(axis=1))
    if areZero[0].size != 0:
        redundant = points[areZero]
        points = np.delete(points, areZero, axis=0)
        ngNormals = np.asarray([np.argmin(np.linalg.norm(points - point, axis=1))
                                for point in redundant])
        normals[areZero] = normals[ngNormals]

    if associateCoords:
        return [(np.asarray(point), np.asarray(normal)) for point, normal in zip(points, normals)]
    else:
        return np.asarray(normals)


# --------------- Quaternion utilities ---------------

def rotation_to_quaternion(rot, type="+"):
    q = np.zeros((1, 4))
    q[0, 0] = 0.5 * np.sqrt(1 + rot[0, 0] + rot[1, 1] + rot[2, 2])
    q[0, 1] = 0.5 * ((rot[2, 1] - rot[1, 2]) / np.abs(rot[2, 1] - rot[1, 2])) * np.sqrt(1 + rot[0, 0] - rot[1, 1] - rot[2, 2])
    q[0, 2] = 0.5 * ((rot[0, 2] - rot[2, 0]) / np.abs(rot[0, 2] - rot[2, 0])) * np.sqrt(1 - rot[0, 0] + rot[1, 1] - rot[2, 2])
    q[0, 3] = 0.5 * ((rot[1, 0] - rot[0, 1]) / np.abs(rot[1, 0] - rot[0, 1])) * np.sqrt(1 - rot[0, 0] - rot[1, 1] + rot[2, 2])
    if type == "+":
        return q
    elif type == "-":
        return -q

def quaternion_to_rotation(q):
    A = np.zeros((3, 3))
    I = np.eye(3)
    A[0, 1] = -q[0, 3]
    A[0, 2] = q[0, 2]
    A[1, 0] = q[0, 3]
    A[1, 2] = -q[0, 1]
    A[2, 0] = -q[0, 2]
    A[2, 1] = q[0, 1]
    rot = I + 2 * q[0, 0] * A + 2 * A.dot(A)
    return rot

def quaternion_distance(q1, q2):
    q = q1 - q2
    module2_q = q[0, 0] ** 2 + np.linalg.norm(q[0, 1:]) ** 2
    dist = 2 * np.arccos(1 - module2_q / 2)
    return dist

def quaternion_inverse(q):
    q[1:] = -q[1:]
    q /= np.linalg.norm(q) ** 2
    return q

def weighted_tensor(tuple_q):
    T = np.zeros((4, 4))
    w_i = 1 / len(tuple_q)
    for q in tuple_q:
        T += w_i * q.T.dot(q)
    return T

def normalized_Principal_Eigenvector(T, q_bar):
    new_q_bar = T.dot(q_bar.T)
    new_q_bar /= np.linalg.norm(new_q_bar)
    return new_q_bar.T

def mean_quaternion(T, q_bar_0=None, tol=1e-3):
    if q_bar_0 == None:
        q_bar_0 = np.random.rand(1, 4)
    cost = sys.maxsize
    while (cost > tol):
        q_bar_f = normalized_Principal_Eigenvector(T, q_bar_0)
        cost = quaternion_distance(q_bar_f, q_bar_0)
        q_bar_0 = q_bar_f
    return q_bar_f, cost


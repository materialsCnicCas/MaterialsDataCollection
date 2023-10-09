#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project    : DataCollectionSoftware 
@File       : cifParser.py
@IDE        : PyCharm 
@Author     : zychen@cnic.cn
@Date       : 2023/9/15 19:00 
@Description: 
"""

import math
import re
import time

import numpy as np
import spglib


import pymongo
import os

from utils.tools.periodicTable import PTable
work_time = time.strftime("%Y-%m-%d %H:%M:%S")

class CifParse:
    """
    解析cif文件
    """
    def __init__(self, filename):
        self.filename = filename
        self.cellParameters = dict()
        self.density = 0
        self.structuretype = None
        with open(self.filename, encoding='utf-8') as f:
            self.lines = f.readlines()
        for i, line in enumerate(self.lines):
            line = line.strip('\n').strip(' ')
            if line.find('_database_code_ICSD') != -1:
                line2 = line.split()
                self.collectioncode = int(line2[1])
            elif line.find('_chemical_formula_structural') != -1:
                line2 = line.split()
                self.SimplestFormula = ''.join(line2[1:])
                if "'" in self.SimplestFormula:
                    self.SimplestFormula = self.SimplestFormula.replace("'", "")
            elif line.find('_chemical_name_structure_type') != -1:
                line2 = line.split()
                self.structuretype = line2[1]
            elif line.find('_cell_length_a') != -1:
                line2 = line.split()
                if "'" in line2[1]:
                    line2[1] = line2[1][1:-1]
                self.cellParameters['a'] = float(line2[1].split('(')[0])
            elif line.find('_cell_length_b') != -1:
                line2 = line.split()
                if "'" in line2[1]:
                    line2[1] = line2[1][1:-1]
                self.cellParameters['b'] = float(line2[1].split('(')[0])
            elif line.find('_cell_length_c') != -1:
                line2 = line.split()
                if "'" in line2[1]:
                    line2[1] = line2[1][1:-1]
                self.cellParameters['c'] = float(line2[1].split('(')[0])
            elif line.find('_cell_angle_alpha') != -1:
                line2 = line.split()
                if "'" in line2[1]:
                    line2[1] = line2[1][1:-1]
                self.cellParameters['alpha'] = float(line2[1].split('(')[0])
            elif line.find('_cell_angle_beta') != -1:
                line2 = line.split()
                if "'" in line2[1]:
                    line2[1] = line2[1][1:-1]
                self.cellParameters['beta'] = float(line2[1].split('(')[0])
            elif line.find('_cell_angle_gamma') != -1:
                line2 = line.split()
                if "'" in line2[1]:
                    line2[1] = line2[1][1:-1]
                self.cellParameters['gamma'] = float(line2[1].split('(')[0])
            elif line.find('_cell_volume') != -1:
                line2 = line.split()
                if "'" in line2[1]:
                    line2[1] = line2[1][1:-1]
                self.cellVolume = float(line2[1])
            elif line.find('_cell_formula_units_Z') != -1:
                line2 = line.split()
                if "'" in line2[1]:
                    line2[1] = line2[1][1:-1]
                self.z = float(line2[1])
            elif line.find('_exptl_crystal_density_diffrn') != -1:
                line2 = line.split()
                if "'" in line2[1]:
                    line2[1] = line2[1][1:-1]
                self.density = float(line2[1])
            elif line.find('_space_group_IT_number') != -1:
                line2 = line.split()
                if "'" in line2[1]:
                    line2[1] = line2[1][1:-1]
                self.spacegroupNumber = int(line2[1])
            elif line.find('_space_group_name_H-M_alt') != -1:
                line2 = line.split()
                self.fullsymbol = ''.join(line2[1:])
                if "'" in self.fullsymbol:
                    self.fullsymbol = self.fullsymbol[1:-1]
            elif line.find('_space_group_symop_operation_xyz') !=-1:
                self.symopstart = i+1
                for j,a in enumerate(self.lines[i+1:]):
                    if a.find('loop_') != -1:
                        self.symopend = self.symopstart + j
                        break
            elif line.find('_atom_site_occupancy') != -1:
                self.sitesstart = i + 1
                for j, a in enumerate(self.lines[i+1:]):
                    if a.find('loop_') != -1 or a.find('#End') != -1:
                        self.sitesend = self.sitesstart + j
                        break

    def count(self, a):
        length = len(a)
        n = 0
        sum = 0
        while n < length:
            if a[n] == '(':
                i, m = self.count(a[n + 1:])
                n = n + 1 + i
                sum += m
            elif a[n] == ')':
                b = ''
                n = n + 1
                while n < length:
                    if '0' <= a[n] <= '9' or a[n] == '.':
                        b += a[n]
                        n += 1
                    else:
                        break
                if b == '':
                    b = 1
                sum = sum * float(b)
                return n, sum
            else:
                x = ''
                y = ''
                while n < length:
                    if 'A' <= a[n] <= 'Z':
                        if x != '':
                            break
                        x += a[n]
                    elif '0' <= a[n] <= '9' or a[n] == '.':
                        y += a[n]
                    else:
                        if a[n] != '(' and a[n] != ')':
                            x += a[n]
                        else:
                            break
                    n += 1
                if y == '':
                    y = 1
                sum += PTable().detailed[x]["Atomic mass"] * float(y)
        return n, sum

    def getDensity(self):
        """
        获取密度
        :return:
        """
        if not self.density:
            formula = self.SimplestFormula
            _, sum = self.count(formula)
            rou = sum * self.z * 10**24 / (self.cellVolume * 6.02 * 10**23)
            return rou
        else:
            return self.density

    def getElementCount(self):
        """
        获取元素种类数
        :return:
        """
        formula = self.SimplestFormula
        visited = []
        s = ''
        for x in formula:
            if 'A' <= x <= 'Z':
                if s != '' and s not in visited:
                    visited.append(s)
                    s = ''
                s += x
            elif 'a' <= x <= 'z':
                s += x
            else:
                if s != '' and s not in visited:
                    visited.append(s)
                    s = ''
        if s != '' and s not in visited:
            visited.append(s)
        return len(visited)

    def getReducedSites(self):
        sites = []
        for id in range(self.sitesstart, self.sitesend):
            site = dict()
            data = self.lines[id]
            if "'" in data:
                data = data.replace("'", "")
            data = data.split()
            a = ''
            b = ''
            for j in data[0]:
                if '0' <= j <= '9' or j == '.':
                    b += j
                else:
                    a += j
            site['AtomicSymbol'] = a
            site['AtomicLabel'] = float(b)
            ox = data[1][-1] + data[1][len(a):-1]
            site['OxState'] = float(ox)
            site['Multiplicity'] = int(data[2])
            site['WyckoffSymbol'] = data[3]
            position = []
            for j in data[4:7]:
                x = j
                if '(' in x:
                    x = x.split('(')[0]
                position.append(float(x))
            site['Position'] = position
            site['Occupancy'] = float(data[-1].split('(')[0])
            sites.append(site)
        return sites

    def getComposition(self):
        composition = []
        symbol_visited = dict()
        sites = self.getReducedSites()
        for x in sites:
            spec = x['AtomicSymbol']
            if spec not in symbol_visited:
                atom = dict()
                length = len(symbol_visited)
                symbol_visited[spec] = length
                atom['AtomicSymbol'] = spec
                atom['AtomicNumber'] = PTable().detailed[spec]["Atomic no"]
                atom['AtomicMass'] = PTable().detailed[spec]["Atomic mass"]
                atom['Amount'] = x['Multiplicity'] * x['Occupancy']
                composition.append(atom)
            else:
                id = symbol_visited[spec]
                composition[id]['Amount'] += x['Multiplicity'] * x['Occupancy']
        return composition

    def getSpaceGroup(self):
        spacegroup = dict()
        spacegroup['SpacegroupNumber'] = self.spacegroupNumber
        spacegroup['FullSpacegroupsymbol'] = self.fullsymbol
        return spacegroup

    def getSymOperations(self):
        symops = []
        for id in range(self.symopstart, self.symopend):
            data = self.lines[id]
            data = data.split("'")
            xyz_string = data[1]
            rot_matrix = np.zeros((3, 3))
            trans = np.zeros(3)
            toks = xyz_string.strip().replace(" ", "").lower().split(",")
            re_rot = re.compile(r"([+-]?)([\d\.]*)/?([\d\.]*)([x-z])")
            re_trans = re.compile(r"([+-]?)([\d\.]+)/?([\d\.]*)(?![x-z])")
            for i, tok in enumerate(toks):
                # build the rotation matrix
                for m in re_rot.finditer(tok):
                    factor = -1.0 if m.group(1) == "-" else 1.0
                    if m.group(2) != "":
                        factor *= float(m.group(2)) / float(m.group(3)) if m.group(3) != "" else float(m.group(2))
                    j = ord(m.group(4)) - 120
                    rot_matrix[i, j] = factor
                # build the translation vector
                for m in re_trans.finditer(tok):
                    factor = -1 if m.group(1) == "-" else 1
                    num = float(m.group(2)) / float(m.group(3)) if m.group(3) != "" else float(m.group(2))
                    trans[i] = num * factor
            #print(rot_matrix,trans)
            rotation_matrix = np.array(rot_matrix)
            translation_vec = np.array(trans)
            if rotation_matrix.shape != (3, 3):
                raise ValueError("Rotation Matrix must be a 3x3 numpy array.")
            if translation_vec.shape != (3,):
                raise ValueError("Translation vector must be a rank 1 numpy array with 3 elements.")
            affine_matrix = np.eye(4)
            affine_matrix[0:3][:, 0:3] = rotation_matrix
            affine_matrix[0:3][:, 3] = translation_vec
            symops.append(affine_matrix)
        return symops

    def getOperationsFromHallNumber(self):
        spacegroupNumber =  self.spacegroupNumber
        operations = spglib.get_symmetry_from_database(spacegroupNumber)
        return operations

    def getCrystalSystem(self):
        n = self.spacegroupNumber
        if n in [1, 2]:
            return 'Triclinic'
        elif n in range(3, 16):
            return 'Monoclinic'
        elif n in range(16, 75):
            return 'Orthorhombic'
        elif n in range(75, 143):
            return 'Tetragonal'
        elif n in range(143, 146) or n == 147 or n in range(149, 155) or n in range(156, 160) or n in range(162, 166):
            return 'Trigonal'
        elif n in [146, 148, 155, 160, 161, 166, 167]:
            return 'Rhombohedral'
        elif n in range(168, 195):
            return 'Hexagonal'
        elif n in range(195, 231):
            return 'Cubic'

    def getPointGroup(self):
        n = self.spacegroupNumber
        if n == 1:
            pointgroupSymbol = '1'
            SchoenfliesNotation = 'C1'
        elif n == 2:
            pointgroupSymbol = '-1'
            SchoenfliesNotation = 'Ci'
        elif n in range(3, 6):
            pointgroupSymbol = '2'
            SchoenfliesNotation = 'C2'
        elif n in range(6, 10):
            pointgroupSymbol = 'm'
            SchoenfliesNotation = 'Cs'
        elif n in range(10, 16):
            pointgroupSymbol = '2/m'
            SchoenfliesNotation = 'C2h'
        elif n in range(16, 25):
            pointgroupSymbol = '222'
            SchoenfliesNotation = 'D2'
        elif n in range(25, 47):
            pointgroupSymbol = 'mm2'
            SchoenfliesNotation = 'C2v'
        elif n in range(47, 75):
            pointgroupSymbol = 'mmm'
            SchoenfliesNotation = 'D2h'
        elif n in range(75, 81):
            pointgroupSymbol = '4'
            SchoenfliesNotation = 'C4'
        elif n in range(81, 83):
            pointgroupSymbol = '-4'
            SchoenfliesNotation = 'S4'
        elif n in range(83, 89):
            pointgroupSymbol = '4/m'
            SchoenfliesNotation = 'C4h'
        elif n in range(89, 99):
            pointgroupSymbol = '422'
            SchoenfliesNotation = 'D4'
        elif n in range(99, 111):
            pointgroupSymbol = '4mm'
            SchoenfliesNotation = 'C4v'
        elif n in range(111, 123):
            pointgroupSymbol = '-42m'
            SchoenfliesNotation = 'D2d'
        elif n in range(123, 143):
            pointgroupSymbol = '4/mmm'
            SchoenfliesNotation = 'D4h'
        elif n in range(143, 147):
            pointgroupSymbol = '3'
            SchoenfliesNotation = 'C3'
        elif n in range(147, 149):
            pointgroupSymbol = '-3'
            SchoenfliesNotation = 'S6'
        elif n in range(149, 156):
            pointgroupSymbol = '32'
            SchoenfliesNotation = 'D3'
        elif n in range(156, 162):
            pointgroupSymbol = '3m'
            SchoenfliesNotation = 'C3v'
        elif n in range(162, 168):
            pointgroupSymbol = '-3m'
            SchoenfliesNotation = 'D3d'
        elif n in range(168, 174):
            pointgroupSymbol = '6'
            SchoenfliesNotation = 'C6'
        elif n == 174:
            pointgroupSymbol = '-6'
            SchoenfliesNotation = 'C3h'
        elif n in range(175, 177):
            pointgroupSymbol = '6/m'
            SchoenfliesNotation = 'C6h'
        elif n in range(177, 183):
            pointgroupSymbol = '622'
            SchoenfliesNotation = 'D6'
        elif n in range(183, 187):
            pointgroupSymbol = '6mm'
            SchoenfliesNotation = 'C6v'
        elif n in range(187, 191):
            pointgroupSymbol = '-6m2'
            SchoenfliesNotation = 'D3h'
        elif n in range(191, 195):
            pointgroupSymbol = '6/mmm'
            SchoenfliesNotation = 'D6h'
        elif n in range(195, 200):
            pointgroupSymbol = '23'
            SchoenfliesNotation = 'T'
        elif n in range(200, 207):
            pointgroupSymbol = 'm-3'
            SchoenfliesNotation = 'Th'
        elif n in range(207, 215):
            pointgroupSymbol = '432'
            SchoenfliesNotation = 'O'
        elif n in range(215, 220):
            pointgroupSymbol = '-43m'
            SchoenfliesNotation = 'Td'
        else:
            pointgroupSymbol = 'm-33m'
            SchoenfliesNotation = 'Oh'
        doc = {
            'PointgroupSymbol': pointgroupSymbol,
            'SchoenfliesNotation': SchoenfliesNotation
        }
        return doc

    def getBravaisLattice(self):
        latticesystem = self.getCrystalSystem()
        spacesymbol = self.getSpaceGroup()["FullSpacegroupsymbol"]
        if latticesystem == "Triclinic":
            return "PrimitiveTriclinic"
        elif latticesystem == "Monoclinic":
            if "P" in spacesymbol:
                return "PrimitiveMonoclinic"
            elif "C" in spacesymbol:
                return "BaseCenteredMonoclinic"
        elif latticesystem == "Orthorhombic":
            if "P" in spacesymbol:
                return "PrimitiveOrthorhombic"
            elif "C" in spacesymbol or "A" in spacesymbol:
                return "BaseCenteredOrthorhombic"
            elif "I" in spacesymbol:
                return "BodyCentredOrthorhombic"
            elif "F" in spacesymbol:
                return "FaceCenteredOrthorhombic"
        elif latticesystem == "Tetragonal":
            if "P" in spacesymbol:
                return "PrimitiveTetragonal"
            elif "I" in spacesymbol:
                return "BodyCenteredTetragonal"
        elif latticesystem == "Rhombohedral":
            return "PrimitiveRhombohedral"
        elif latticesystem == "Hexagonal" or latticesystem == 'Trigonal':
            return "PrimitiveHexagonal"
        elif latticesystem == "Cubic":
            if "P" in spacesymbol:
                return "PrimitiveCubic"
            elif "I" in spacesymbol:
                return "BodyCentredCubic"
            elif "F" in spacesymbol:
                return "FaceCenteredCubic"

    def getLatticeVectors(self):
        # assume a, b, c in Angstroms,
        # alpha, beta, gamma in degrees
        cellParameters = self.cellParameters
        #{'a': 4.3983, 'b': 4.3983, 'c': 2.873, 'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0}
        a = float(cellParameters['a'])
        b = float(cellParameters['b'])
        c = float(cellParameters['c'])
        alpha = float(cellParameters['alpha'])
        beta = float(cellParameters['beta'])
        gamma = float(cellParameters['gamma'])
       # a, b, c, alpha, beta, gamma = [float(x) for x in cellParameters]
        cfac = (2 * math.pi) / 360
        alpha = alpha * cfac
        beta = beta * cfac
        gamma = gamma * cfac
        s = 0.5 * (alpha + beta + gamma)
        vol = 2 * a * b * c * (math.sin(s) * math.sin(s - alpha) *
                               math.sin(s - beta) * math.sin(s - gamma)) ** 0.5
        lattice_vectors = [
            ["%.15f" % (vol / (math.sin(alpha) * b * c)),
             "%.15f" % (a * math.cos(gamma) *
                        math.sin(alpha)),
             "%.15f" % (a * math.cos(beta))],
            ["0.000000000000000",
             "%.15f" % (b * math.sin(alpha)),
             "%.15f" % (b * math.cos(alpha))],
            ["0.000000000000000", "0.000000000000000", "%.15f" % c]
        ]
        return lattice_vectors

    def getSitesFromReducedSites(self):
        """
        获取点位坐标信息
        :return:
        """
        Newsites = []
        Newsites_temp = []
        reducedsites = self.getReducedSites()
        operations = self.getSymOperations()
        for site in reducedsites:
            position=site['Position']
            position_point = np.array([position[0], position[1], position[2], 1])
            for op in operations:
                Newsite = dict()
                NewsitePos = np.dot(op, position_point)[0:3]
                Newsite['AtomicSymbol'] = site['AtomicSymbol']
                Newsite['AtomicLabel'] = site['AtomicLabel']
                Newsite['OxState'] = site['OxState']
                Newsite['WyckoffSymbol'] = site['WyckoffSymbol']
                Newsite['Position'] = NewsitePos
                Newsite['Occupancy'] = site['Occupancy']
                Newsites_temp.append(Newsite)
        Newsites_tempTxt = [str(line) for line in Newsites_temp]
        Newsites_tempTxtSet = set(Newsites_tempTxt)
        for item in Newsites_tempTxtSet:
            orignal_index = Newsites_tempTxt.index(item)
            Unique_site = Newsites_temp[orignal_index]
            Newsites.append(Unique_site)
        return Newsites

    def getCentering(self):
        spacesymbol = self.getSpaceGroup()["FullSpacegroupsymbol"]
        bravaislattice = self.getBravaisLattice()
        if bravaislattice in ['PrimitiveTriclinic', 'PrimitiveMonoclinic', 'PrimitiveOrthorhombic', 'PrimitiveTetragonal',
                              'PrimitiveHexagonal', 'PrimitiveCubic']:
            return 'Primitive'
        elif bravaislattice in ['BodyCentredOrthorhombic', 'BodyCenteredTetragonal', 'BodyCentredCubic']:
            return 'BodyCentered'
        elif bravaislattice in ['FaceCenteredOrthorhombic', 'FaceCenteredCubic']:
            return 'FaceCentered'
        elif bravaislattice == 'BaseCenteredOrthorhombic' and 'A' in spacesymbol:
            return 'A-faceCentered'
        elif bravaislattice == 'BaseCenteredMonoclinic' or (bravaislattice == 'BaseCenteredOrthorhombic' and 'C' in spacesymbol):
            return 'C-faceCentered'
        elif bravaislattice == 'PrimitiveRhombohedral':
            return 'Rhombohedral'


def cifToBson(filename, user_id, group_id):
    cif = CifParse(filename)
    doc = {
        'SimplestFormula': cif.SimplestFormula,
        'CellParameters': cif.cellParameters,
        'CellVolume': cif.cellVolume,
        'Z': cif.z,
        'Density': cif.getDensity(),
        'ElementCount': cif.getElementCount(),
        'ReducedSites': cif.getReducedSites(),
        'Composition': cif.getComposition(),
        'SpaceGroup': cif.getSpaceGroup(),
        'CrystalSystem': cif.getCrystalSystem(),
        'PointGroup': cif.getPointGroup(),
        'BravaisLattice': cif.getBravaisLattice(),
        'Centering': cif.getCentering(),
        'StructureType': cif.structuretype,
        'CollectionCode': cif.collectioncode,
        'Files': [filename],
        'User': {
            'user': user_id,
            'group': group_id
        },
        'DateCreated': work_time
    }
    return doc


# def to_mongo(doc: dict, database: str, collection: str, host, port):
#     """
#     将文档数据存入MongoDB中
#     :param collection:
#     :param database:
#     :param doc:
#     :param host:
#     :param port:
#     :return:
#        mongo中structure文档的id
#     """
#     conn = pymongo.MongoClient(host=host, port=port)
#     db = conn[database]
#     collection = db[collection]
#     object_id = collection.insert_one(doc).inserted_id
#     conn.close()
#     return object_id

#
# def searchfile(rootpath):
#     for filepath, dirnames, filenames in os.walk(rootpath):
#         for filename in filenames:
#             path1 = os.path.join(filepath, filename)
#             try:
#                 doc = toBson(path1)
#             except:
#                 continue
#             to_mongo(doc, 'ICSD', 'Structure', 'localhost', 27017)


# filename = 'D:/icsd_raw/YourCustomFileName_CollCode99461.cif'
# doc = tobson(filename)
# database = 'MatDB'
# collection = 'cif'
# id = to_mongo(doc, database, collection, 'localhost', 27017)
# rootpath = 'D:/icsd_raw'
# searchfile(rootpath)
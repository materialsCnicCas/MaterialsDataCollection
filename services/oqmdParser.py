#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project    : DataCollectionSoftware 
@File       : oqmdParser.py
@IDE        : PyCharm 
@Author     : zychen@cnic.cn
@Date       : 2023/9/4 15:55 
@Description: 
"""
import copy
import time
from math import sqrt, pi

import pymongo
import spglib as spg
from qmpy import Entry, Structure, Calculation, Atom, Prototype, Element, Site

from pymatgen.core.structure import Structure as PyStr
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga

from dao.toMongo import many_to_mongo
from utils.tools.calType import CalType
from utils.tools.highSymmetryKPath import *


def getSimplestFormula(composition):
    """
    获取最简化学式
    :return:
    """
    formula = ''
    nums = []
    common = 1
    for atom in composition:
        nums.append(int(atom["Amount"]))
    for i in range(min(nums)):
        flag = True
        for num in nums:
            if num % (i + 1) != 0:
                flag = False
                break
        if flag:
            common = i + 1
    for atom in composition:
        if atom["Amount"] / common != 1:
            formula = formula + atom["AtomicSymbol"] + str(int(atom["Amount"] / common))
        else:
            formula = formula + atom["AtomicSymbol"]
    return formula, common


def getBravaisLattice(latticesystem, spacenumber, spacesymbol):
    """
    获取bravaislattice
    :return:
    """
    if spacenumber in [146, 148, 155, 160, 161, 166, 167]:
        latticesystem = "Rhombohedral"
    elif latticesystem == "Trigonal":
        latticesystem = "Hexagonal"
    if latticesystem == "Triclinic":
        return "SimpleTriclinic"
    elif latticesystem == "Monoclinic":
        if "P" in spacesymbol:
            return "SimpleMonoclinic"
        elif "C" in spacesymbol:
            return "BaseCenteredMonoclinic"
    elif latticesystem == "Orthorhombic":
        if "P" in spacesymbol:
            return "SimpleOrthorhombic"
        elif "C" in spacesymbol or "A" in spacesymbol:
            return "BaseCenteredOrthorhombic"
        elif "I" in spacesymbol:
            return "BodyCentredOrthorhombic"
        elif "F" in spacesymbol:
            return "FaceCenteredOrthorhombic"
    elif latticesystem == "Tetragonal":
        if "P" in spacesymbol:
            return "SimpleTetragonal"
        elif "I" in spacesymbol:
            return "BodyCenteredTetragonal"
    elif latticesystem == "Rhombohedral":
        return "SimpleRhombohedral"
    elif latticesystem == "Hexagonal":
        return "SimpleRhombohedral"
    elif latticesystem == "Cubic":
        if "P" in spacesymbol:
            return "SimpleCubic"
        elif "I" in spacesymbol:
            return "BodyCentredCubic"
        elif "F" in spacesymbol:
            return "FaceCenteredCubic"


def getHighSymmetryKPath(bravaislattice, lattice, reclattice, prilattice):
    """
    生成布里渊区高对称点路径的lattice类型
    :return:
    """
    # TODO: a,b,c参数的来源需再校核
    a = lattice[0]
    b = lattice[1]
    c = lattice[2]
    a0 = sqrt(np.dot(a, a))
    b0 = sqrt(np.dot(b, b))
    c0 = sqrt(np.dot(c, c))
    alpha = np.arccos(np.dot(a, c) / (a0 * c0))

    a = reclattice[0]
    b = reclattice[1]
    c = reclattice[2]
    a1 = sqrt(np.dot(a, a))
    b1 = sqrt(np.dot(b, b))
    c1 = sqrt(np.dot(c, c))
    kalpha = np.arccos(np.dot(a, c) / (a1 * c1))
    kbeta = np.arccos(np.dot(b, c) / (b1 * c1))
    kgamma = np.arccos(np.dot(a, b) / (a1 * b1))

    a = prilattice[0]
    c = prilattice[2]
    a2 = sqrt(np.dot(a, a))
    c2 = sqrt(np.dot(c, c))
    alpha2 = np.arccos(np.dot(a, c) / (a2 * c2))
    kpath = {}
    if bravaislattice == 'SimpleCubic':
        kpath = cubic()
    elif bravaislattice == 'FaceCenteredCubic':
        kpath = fcc()
    elif bravaislattice == 'BodyCenteredCubic':
        kpath = bcc()
    elif bravaislattice == 'SimpleTetragonal':
        kpath = tet()
    elif bravaislattice == 'BodyCenteredTetragonal':
        if a0 > c0:
            kpath = bctet1(a0, c0)
        else:
            kpath = bctet2(a0, c0)
    elif bravaislattice == 'SimpleOrthorhombic':
        kpath = orc()
    elif bravaislattice == 'FaceCenteredOrthorhombic':
        if 1 / a0 ** 2 > 1 / b0 ** 2 + 1 / c0 ** 2:
            kpath = orcf1(a0, b0, c0)
        elif 1 / a0 ** 2 < 1 / b0 ** 2 + 1 / c0 ** 2:
            kpath = orcf2(a0, b0, c0)
        else:
            kpath = orcf3(a0, b0, c0)
    elif bravaislattice == 'BodyCenteredOrthorhombic':
        kpath = orci(a0, b0, c0)
    elif bravaislattice == 'BaseCenteredOrthorhombic':
        kpath = orcc(a0, b0)
    elif bravaislattice == 'SimpleHexagonal':
        kpath = hex()
    elif bravaislattice == 'SimpleRhombohedral':
        if alpha2 < 90 * pi / 180:
            kpath = rhl1(alpha2)
        else:
            kpath = rhl2(alpha2)
    elif bravaislattice == 'SimpleMonoclinic':
        kpath = mcl(b0, c0, alpha)
    elif bravaislattice == 'BaseCenteredMonoclinic':
        if kgamma > 90 * pi / 180:
            kpath = mclc1(a0, b0, c0, alpha)
        elif kgamma == 90 * pi / 180:
            kpath = mclc2(a0, b0, c0, alpha)
        elif kgamma < 90 * pi / 180 and b0 * cos(alpha) / c0 + b0 ** 2 * sin(alpha) ** 2 / a0 ** 2 < 1:
            kpath = mclc3(a0, b0, c0, alpha)
        elif kgamma < 90 * pi / 180 and b0 * cos(alpha) / c0 + b0 ** 2 * sin(alpha) ** 2 / a0 ** 2 == 1:
            kpath = mclc4(a0, b0, c0, alpha)
        elif kgamma < 90 * pi / 180 and b0 * cos(alpha) / c0 + b0 ** 2 * sin(alpha) ** 2 / a0 ** 2 > 1:
            kpath = mclc5(a0, b0, c0, alpha)
    elif bravaislattice == 'SimpleTriclinic':
        if kalpha > 90 * pi / 180 and kbeta > 90 * pi / 180 and kgamma >= 90 * pi / 180:
            kpath = tria()
        elif kalpha < 90 * pi / 180 and kbeta < 90 * pi / 180 and kgamma <= 90 * pi / 180:
            kpath = trib()
    if kpath != {}:
        for key in kpath['kpoints'].keys():
            kpath['kpoints'][key] = kpath['kpoints'][key].tolist()
    return kpath


def getStructureDoc(structure: Structure) -> dict:
    formula = structure.name
    generalFormula = structure.composition.generic
    numberOfElments = structure.ntypes
    numberofSites = structure.nsites
    # composition
    compositions = []
    comp = structure.comp
    for atomic, num in comp.items():
        element = Element.objects.get(symbol=atomic)
        doc = {
            "AtomicSymbol": atomic,
            "AtomicNumber": element.z,
            "AtomicMass": element.mass,
            "Amount": num
        }
        compositions.append(doc)
    simplestFormula, common = getSimplestFormula(compositions)  # 最大公约数 common
    # stoichiometry
    stoichiometry = copy.deepcopy(compositions)
    for sto in stoichiometry:
        sto["AtomicPercent"] = sto["Amount"] / numberofSites
    # Lattice
    lattice = [[structure.x1, structure.x2, structure.x3],
               [structure.y1, structure.y2, structure.y3],
               [structure.z1, structure.z2, structure.z3]]
    # latticeParameters
    lattice_martrix = np.array(lattice, dtype=np.float64).reshape((3, 3))
    a = lattice_martrix[0]
    b = lattice_martrix[1]
    c = lattice_martrix[2]
    a0 = sqrt(np.dot(a, a))
    b0 = sqrt(np.dot(b, b))
    c0 = sqrt(np.dot(c, c))
    alpha = np.arccos(np.dot(a, c) / (a0 * c0))
    beta = np.arccos(np.dot(b, c) / (b0 * c0))
    gamma = np.arccos(np.dot(a, b) / (a0 * b0))
    latticeParameters = {
        'a': a0,
        'b': b0,
        'c': c0,
        'alpha': alpha,
        'beta': beta,
        'gamma': gamma
    }
    volume = structure.volume

    # latticeSystem
    latticeSystem = structure.spacegroup.lattice_system
    # reciprocalLattice
    reciprocalLattice = np.linalg.inv(lattice_martrix)
    # Sites
    sites = []
    atoms = Atom.objects.filter(structure_id=structure.id)
    sites_object = Site.objects.filter(structure_id=structure.id)
    sites_dict = {}
    for site in sites_object:
        if site.wyckoff is not None:
            sites_dict[str(site.id)] = (site.wyckoff.symbol, site.wyckoff.multiplicity)
        else:
            sites_dict[str(site.id)] = (None, None)
    for atom in atoms:
        # position force atom occuancy magmom
        position = [atom.x, atom.y, atom.z]
        force = [atom.fx, atom.fy, atom.fz]
        magmom = atom.magmom
        occupancy = atom.occupancy
        atom_doc = {
            "AtomicSymbol": atom.element.symbol,
            "AtomicNumber": atom.element.z,
            "AtomicMass": atom.element.mass
        }
        site = {
            "Position": position,
            "Force": force,
            "Atom": atom_doc,
            "Occupancy": occupancy,
            "Magmom": magmom,
            "OxState": atom.ox,
            "Multiplicity": sites_dict[str(atom.site_id)][1],
            "WyckoffSymbol": sites_dict[str(atom.site_id)][0]
        }
        sites.append(site)
    # spaceGroup、 pointgroup
    pos = [i["Position"] for i in sites]
    numbers = [i["Atom"]["AtomicNumber"] for i in sites]
    cell = (lattice_martrix, pos, numbers)
    dataset = spg.get_symmetry_dataset(cell, symprec=0.01, angle_tolerance=0.005, hall_number=0)
    if dataset is not None:
        pointgroup = dataset['pointgroup']
        spacegroup_type = spg.get_spacegroup_type(dataset["hall_number"])
        spacegroup = {
            "SpaceGroupSymbol": spacegroup_type['international_short'],
            "SpaceGroupInternationalSymbol": spacegroup_type['international_full'],
            "SchoenfliesSymbol": spacegroup_type['schoenflies'],
            "SpaceGroupNumber": spacegroup_type['number'],
            "hallNumber": dataset['hall_number'],
            "hallSymbol": spacegroup_type['hall_symbol']
        }
    else:
        pointgroup = None
        spacegroup = {
            "SpaceGroupSymbol": structure.spacegroup.hm,
            "SpaceGroupInternationalSymbol": structure.spacegroup.hall,
            "SchoenfliesSymbol": structure.spacegroup.schoenflies,
            "SpaceGroupNumber": structure.spacegroup_id,
            "hallNumber": None,
            "hallSymbol": None
        }
    # primitiveLattice、 conventionalLattice
    primitiveLattice, _, _ = spg.find_primitive(cell, symprec=0.001, angle_tolerance=1)
    struct = PyStr(lattice_martrix, numbers, pos)
    conventional_struct = sga(struct).get_conventional_standard_structure()
    conventionalLattice = conventional_struct.lattice.matrix

    # bravaisLattice
    spacenumber = spacegroup["spacegroupNumber"]
    spacesymbol = spacegroup["spacegroupSymbol"]
    bravaisLattice = getBravaisLattice(latticeSystem, spacenumber, spacesymbol)
    # highSymmetryKpath
    highSymmetryKpath = getHighSymmetryKPath(bravaisLattice, conventionalLattice, reciprocalLattice, primitiveLattice)
    # protoType
    protoType = Prototype.objects.filter(structure_id=structure.id)
    prototype = None
    for pto in protoType:
        prototype = pto.name
    # CellStress
    sxx = structure.sxx
    sxy = structure.sxy
    syy = structure.syy
    syz = structure.syz
    szx = structure.szx
    szz = structure.szz
    cellStress = [[sxx, sxy, szx], [sxy, syy, syz], [szx, syz, szz]]
    structure_doc = {
        "Formula": formula,
        "GeneralFormula": generalFormula,
        "SimplestFormula": simplestFormula,
        "NumberOfElements": numberOfElments,
        "NumberOfSites": numberofSites,
        "Composition": compositions,
        "Stoichiometry": stoichiometry,
        "LatticeParameters": latticeParameters,
        "Lattice": lattice,
        "Volume": volume,
        "PrimitiveLattice": primitiveLattice.tolist(),
        "ConventionalLattice": conventionalLattice.tolist(),
        "SpaceGroup": spacegroup,
        "PointGroup": pointgroup,
        "LatticeSystem": latticeSystem,
        "BravaisLattice": bravaisLattice,
        "HighSymmetryKpath": highSymmetryKpath,
        "ReciprocalLattice": reciprocalLattice.tolist(),
        "Sites": sites,
        "Prototype": prototype,
        "CellStress": cellStress
    }

    return structure_doc, common


def oqmdMain(database, collections, host, port, userId, groupId, outFile):
    work_time = time.strftime("%Y-%m-%d %H:%M:%S")
    queryset = Entry.objects.all()
    entry_ids = [entry.id for entry in queryset]
    print('Number of entries in OQMD: %d' % len(entry_ids), file=outFile)
    error_ids = []
    geoMongoDoc = []
    staticMongoDoc = []
    count = 0
    for entry_id in entry_ids:
        # print("Extracting Entry_id: ", entry_id)
        if count == 1000:
            break
        count += 1
        try:
            # if True:
            structures_object = Structure.objects.filter(entry_id=entry_id)
            structures = {}
            for struct in structures_object:
                structures[str(struct.id)] = getStructureDoc(struct)
            calculations = Calculation.objects.filter(entry_id=entry_id)
            for calculation in calculations:
                input_id = calculation.input_id
                output_id = calculation.output_id
                inStrucuture_doc, common_1 = structures[str(input_id)] if input_id is not None else (None, 1)
                outStrucuture_doc, common_2 = structures[str(output_id)] if output_id is not None else (None, 1)
                nsteps = calculation.nsteps  # of ionic steps
                converged = True if calculation.converged == 1 else False  # 收敛性
                runtime = calculation.runtime
                parameters = calculation.settings
                configuration = calculation.configuration
                energy = calculation.energy * max(common_1, common_2) if calculation.energy is not None else None
                energy_pa = calculation.energy_pa
                band_gap = calculation.band_gap
                linearMagneticMoment = calculation.magmom
                if configuration is not None and 'static' in configuration and CalType.StaticCalculation in collections:
                    # StaticCalculation
                    static_doc = {
                        "InputStructure": inStrucuture_doc,
                        "OutputStructure": outStrucuture_doc,
                        "Energy": energy,
                        "Parameter": parameters,
                        "EnergyPerAtom": energy_pa,
                        "EleConvergence": converged,
                        "LinearMagneticMoment": linearMagneticMoment,
                        "GapFromStatic": band_gap,
                        "User": {
                            "user": userId,
                            "group": groupId
                        },
                        'Files': [calculation.id],
                        'DateCreated': work_time,
                        'CalculationType': CalType.StaticCalculation
                    }
                    staticMongoDoc.append(static_doc)
                    if len(staticMongoDoc) >= 1000:
                        many_to_mongo(staticMongoDoc, database, CalType.StaticCalculation, host, port)
                        staticMongoDoc = []

                elif CalType.GeometryOptimization in collections:
                    # GeometryOptimization
                    geometry_doc = {
                        "InStructure": inStrucuture_doc,
                        "OutStructure": outStrucuture_doc,
                        "Parameter": parameters,
                        "NSteps": nsteps,
                        "IonConvergence": converged,
                        "Energy": energy,
                        "EnergyPerAtom": energy_pa,
                        "ElapsedTime": runtime,
                        "GapFromGeo": band_gap,
                        "User": {
                            "user": userId,
                            "group": groupId
                        },
                        'Files': [calculation.id],
                        'DateCreated': work_time,
                        'CalculationType': CalType.GeometryOptimization
                    }
                    geoMongoDoc.append(geometry_doc)
                    if len(geoMongoDoc) >= 1000:
                        many_to_mongo(geoMongoDoc, database, CalType.GeometryOptimization, host, port)
                        geoMongoDoc = []


        except Exception as e:
            print("Error Entry_id: ", entry_id, '   error: ', e, file=outFile)
            error_ids.append(entry_id)
            continue

    if len(geoMongoDoc) != 0:
        many_to_mongo(geoMongoDoc, database, CalType.GeometryOptimization, host, port)
    if len(staticMongoDoc) != 0:
        many_to_mongo(staticMongoDoc, database, CalType.StaticCalculation, host, port)
    return error_ids

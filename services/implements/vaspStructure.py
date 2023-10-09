from typing import List

from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.core.lattice import Lattice as PyLattice

from services.implements.vaspParser import VaspParser
from services.structure import Structure
from utils.tools.lattice import Lattice
from utils.tools.sites import Site


class VaspStructure(Structure):
    def __init__(self, vaspParser: VaspParser = None, lattice: Lattice = None, composition: List = None,
                 sites: List[Site] = None,
                 isinit: bool = False):
        self.vaspParser = vaspParser
        if self.vaspParser is not None and not self.vaspParser.hasSetup:
            self.vaspParser = None
        self.isinit = isinit
        if vaspParser is not None:
            composition = vaspParser.getComposition()
            if isinit:
                lattice = vaspParser.lattice_s
                sites = vaspParser.sites_s
            else:
                lattice = vaspParser.lattice_e
                sites = vaspParser.sites_e

        super().__init__(lattice, composition, sites)

    def setup(self):
        if self.lattice is None or self.sites is None or self.composition is None:
            raise ValueError('Lattice Sites Composition 都不能为空！')
        self.formula = self.getFormula()
        self.simplestFormula = self.getSimplestFormula()
        self.numberOfElements = self.getNElements()
        self.numberOfSites = self.getNSites()
        self.stoichiometry = self.getStoichiometry()
        self.spaceGroup, self.pointGroup = self.getSpacePointGroup()
        if self.vaspParser is not None:
            if self.isinit:
                self.volume = self.vaspParser.volume_s
                self.latticeParameters = self.vaspParser.latticeParameters_s
            else:
                self.volume = self.vaspParser.volume_e
                self.latticeParameters = self.vaspParser.latticeParameters_e
        else:
            self.volume = self.getVolume()
            self.latticeParameters = self.getLatticeParameters()

        self.crystalSystem = self.getCrystalSystem()
        self.bravaisLattice = self.getBravaisLattice()
        self.reciprocalLattice = self.getReciprocalLattice()
        self.primitiveLattice = self.getPrimitiveLattice()
        self.conventionalLattice = self.getConventionalLattice()
        self.highSymmetryKPath = self.getHighSymmetryKPath()
        self.sitesDoc = self.getSitesList()
        self.hashId = self.getHashValue()
        # 判断空间群是否匹配
        spg = SpaceGroup(self.spaceGroup['spacegroupSymbol'])
        lattice = PyLattice(self.conventionalLattice)
        if not spg.is_compatible(lattice, tol=0.01):
            raise ValueError('space group not compatible with lattice')

    def toBson(self):
        doc = {
            "Formula": self.formula,
            "SimplestFormula": self.simplestFormula,
            "GeneralFormula": self.generalFormula,  # 新增
            "NumberOfElements": self.numberOfElements,
            "NumberOfSites": self.numberOfSites,
            "Composition": self.composition,
            "Stoichiometry": self.stoichiometry,
            "LatticeParameters": self.latticeParameters,
            "Lattice": self.lattice.matrix.tolist(),
            "Volume": self.volume,
            "SpaceGroup": self.spaceGroup,
            "PointGroup": self.pointGroup,
            "CrystalSystem": self.crystalSystem,
            "BravaisLattice": self.bravaisLattice,
            "ReciprocalLattice": self.reciprocalLattice.tolist(),
            "PrimitiveLattice": self.primitiveLattice.tolist(),
            "ConventionalLattice": self.conventionalLattice.tolist(),
            "HighSymmetryKpath": self.highSymmetryKPath,
            "Sites": self.sitesDoc,  # 有变化
            "Prototype": self.prototype,  # 新增
            "CellStress": self.cellStress,  # 新增
            "HashValue": self.hashId
        }
        return doc

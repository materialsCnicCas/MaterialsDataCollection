import linecache
import math
import os.path
import time
import xml.etree.cElementTree as ET
import numpy as np

from services.fileParser import FileParser
from utils.assist.incar import Incar
from utils.assist.outcar import OutCar
from utils.assist.oszicar import OsziCar
from utils.tools.calType import CalType
from utils.tools.lattice import Lattice
from utils.tools.periodicTable import PTable
from utils.tools.sites import Site, Atom

from utils.assist.helper import parseVarray
from utils.tools.Electronic import Spin


class VaspParser(FileParser):

    def __init__(self, vaspPath, incarPath, outcarPath, kPointPath, oszicarPath):
        super().__init__(vaspPath)
        self.vaspPath = vaspPath
        self.incarPath = incarPath
        self.outcarPath = outcarPath
        self.kPointPath = kPointPath
        self.oszicarPath = oszicarPath

        tree = ET.parse(vaspPath)
        if tree is None:
            raise ValueError('File content error, not parse!')
        self.root = tree.getroot()
        self.incar = Incar(self.incarPath) if os.path.exists(self.incarPath) else None
        self.outcar = OutCar(self.outcarPath) if os.path.exists(self.outcarPath) else None
        self.oszicar = OsziCar(self.oszicarPath) if os.path.exists(self.oszicarPath) else None
        self.calculationType = None
        self.simplestFormula = None
        self.generalFormula = None
        self.numberOfElements = None
        self.spaceGroup = None
        self.pointGroup = None
        self.crystalSystem = None
        self.resourceUsage = None
        self.cellStress_s = None
        self.cellStress_e = None
        self.electronicProperties = None
        self.magneticProperties = None
        self.atomicCharge = None
        self.atomicMagnetization = None
        self.totalDos = None
        self.partialDos = None
        self.kPoints = None
        self.eigenValues = None
        self.gapFromBand = None
        self.gapFromDos = None
        self.projectedEigen = None
        self.dielectricData = None
        self.opticalProperties = None
        self.linearMagneticMoment = None
        self.elasticProperties = None

    def setup(self):
        self.calculationType = self.getCalculationType()
        self.lattice_s, self.latticeParameters_s = self.getLatticeParameters(isinit=True)
        self.lattice_e, self.latticeParameters_e = self.getLatticeParameters(isinit=False)
        self.composition = self.getComposition()
        self.parameters = self.getParameters()
        self.sites_s = self.getSites(isinit=True)
        self.sites_e = self.getSites(isinit=False)
        self.volume_s = self.getVolume(isinit=True)
        self.volume_e = self.getVolume(isinit=False)
        self.numberOfSites = self.getNumberOfSites()
        self.startTime = self.getStartTime()
        self.software = self.getSoftware()
        self.resourceUsage = self.getResourceUsage()
        self.kPoints = self.getKPoints()

        if self.calculationType == CalType.StaticCalculation:
            self.electronicSteps = self.getElectronicSteps()
            self.thermoDynamicProperties = self.getThermoDynamicProperties()
            if self.outcar is not None:
                self.atomicCharge, self.atomicMagnetization = self.outcar.getAtomicChargeAndAtomicMagnetization()
            if self.oszicar is not None:
                self.linearMagneticMoment = self.oszicar.getLinearMagneticMoment()

        if self.calculationType == CalType.GeometryOptimization:
            self.electronicSteps = self.getElectronicSteps()
            self.ionicSteps = self.getIonicSteps()
            self.gapFromBand = self.getGapFromBand()

        if self.calculationType == CalType.BandStructure:
            self.electronicSteps = self.getElectronicSteps()
            self.thermoDynamicProperties = self.getThermoDynamicProperties()
            if self.outcar is not None:
                self.atomicCharge, self.atomicMagnetization = self.outcar.getAtomicChargeAndAtomicMagnetization()
            if self.oszicar is not None:
                self.linearMagneticMoment = self.oszicar.getLinearMagneticMoment()
            self.eigenValues = self.getEigenValues()
            self.projectedEigen = self.getProjectedEigenvalOnIonOrbitals()
            self.gapFromBand = self.getGapFromBand()

        if self.calculationType == CalType.DensityOfStates:
            self.electronicSteps = self.getElectronicSteps()
            self.thermoDynamicProperties = self.getThermoDynamicProperties()
            if self.outcar is not None:
                self.atomicCharge, self.atomicMagnetization = self.outcar.getAtomicChargeAndAtomicMagnetization()
            if self.oszicar is not None:
                self.linearMagneticMoment = self.oszicar.getLinearMagneticMoment()
            self.totalDos = self.getTotalDos()
            self.partialDos = self.getPartialDos()
            self.gapFromDos = self.getGapFromDos()

        if self.calculationType == CalType.ElasticProperties:
            self.ionicSteps = self.getIonicSteps()
            self.electronicSteps = self.getElectronicSteps()
            if self.outcar is not None:
                self.elasticProperties = self.outcar.getElasticProperties()

        if self.calculationType == CalType.DielectricProperties:
            self.electronicSteps = self.getElectronicSteps()
            self.thermoDynamicProperties = self.getThermoDynamicProperties()
            if self.outcar is not None:
                self.atomicCharge, self.atomicMagnetization = self.outcar.getAtomicChargeAndAtomicMagnetization()
            if self.oszicar is not None:
                self.linearMagneticMoment = self.oszicar.getLinearMagneticMoment()
            self.opticalProperties = self.getOpticalProperties()

        if self.calculationType == CalType.MagneticProperties:
            self.electronicSteps = self.getElectronicSteps()
            self.thermoDynamicProperties = self.getThermoDynamicProperties()
            if self.outcar is not None:
                self.atomicCharge, self.atomicMagnetization = self.outcar.getAtomicChargeAndAtomicMagnetization()
            if self.oszicar is not None:
                self.linearMagneticMoment = self.oszicar.getLinearMagneticMoment()

        self.hasSetup = True

    def getCalculationType(self):
        # getParameters
        para_list = ['IBRION', 'LORBIT', 'LOPTICS', 'LEPSILON', 'LCALCEPS', 'ISIF', 'MAGMOM']
        parameters = self.findPara(para_list)
        if parameters['IBRION'] == 1 or parameters['IBRION'] == 2 or parameters['IBRION'] == 3:
            return CalType.GeometryOptimization
        if parameters['IBRION'] == -1 and parameters['LORBIT']:
            return CalType.DensityOfStates
        if parameters['IBRION'] == -1 and linecache.getline(self.kPointPath,
                                                            3).lower() == 'l':  # kpoints generation para
            return CalType.BandStructure
        if parameters['IBRION'] == -1:
            return CalType.StaticCalculation
        if (parameters['IBRION'] == 5 or parameters['IBRION'] == 6) and parameters['ISIF'] >= 3:
            return CalType.ElasticProperties
        if ('LOPTICS' in parameters.keys() and parameters['LOPTICS']) or parameters['IBRION'] == 7 or parameters[
            'IBRION'] == 8 \
                or (parameters['IBRION'] == 5 and (('LEPSILON' in parameters.keys() and parameters['LEPSILON']) or (
                'LCALCEPS' in parameters.keys() and parameters['LCALCEPS']))) \
                or (parameters['IBRION'] == 6 and (('LEPSILON' in parameters.keys() and parameters['LEPSILON']) or (
                'LCALCEPS' in parameters.keys() and parameters['LCALCEPS']))):
            return CalType.DielectricProperties
        if 'MAGMOM' in parameters.keys():
            return CalType.MagneticProperties
        raise ValueError('无法判断提取类型，无法提取')

    def getComposition(self):
        composition = {}
        composition_full = []
        child = self.root.find("./atominfo/array[@name='atoms']/set")
        for child2 in child:
            count = 1
            if child2[0].text.strip(' ') in composition.keys():
                count += composition[child2[0].text.strip(' ')]
            composition.update({child2[0].text.strip(' '): count})
        for spec in composition.keys():
            atom = {
                "AtomicSymbol": spec,
                "AtomicNumber": PTable().detailed[spec]["Atomic no"],
                "AtomicMass": PTable().detailed[spec]["Atomic mass"],
                "Amount": composition[spec]
            }
            composition_full.append(atom)
        return composition_full

    def getNumberOfSites(self):
        child = self.root.find("./atominfo/atoms")
        return int(child.text)

    def getLatticeParameters(self, isinit: bool = False):
        if isinit:
            child = self.root.find("./structure[@name='initialpos']/crystal/varray[@name='basis']")
        else:
            child = self.root.find("./structure[@name='finalpos']/crystal/varray[@name='basis']")
        a = child[0].text.split()
        a = [float(i) for i in a]
        b = child[1].text.split()
        b = [float(i) for i in b]
        c = child[2].text.split()
        c = [float(i) for i in c]
        lattice = Lattice([a, b, c])

        a0 = math.sqrt(np.dot(a, a))
        b0 = math.sqrt(np.dot(b, b))
        c0 = math.sqrt(np.dot(c, c))
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
        return lattice, latticeParameters

    def getVolume(self, isinit: bool = False):
        if isinit:
            child = self.root.find("./structure[@name='initialpos']/crystal/i[@name='volume']")
        else:
            child = self.root.find("./structure[@name='finalpos']/crystal/i[@name='volume']")
        return float(child.text)

    def getSites(self, isinit: bool = False):
        specs = []
        sites = []
        forces = []
        paras = self.parameters
        child = self.root.find("./atominfo/array[@name='atoms']/set")
        for child2 in child:
            specs.append(child2[0].text.strip(' '))
        if isinit:
            child = self.root.find("./structure[@name='initialpos']/varray")
            forces = [None for _ in range(len(specs))]
        else:
            force_child = self.root.find("./calculation[last()]/varry[@name='forces']")
            for force in force_child:
                forces.append([float(x) for x in force.text.split()])

            child = self.root.find("./structure[@name='finalpos']/varray")

        for i in range(len(specs)):
            p = child[i].text.split()
            p = [float(x) for x in p]
            sites.append(Site(p, Atom(specs[i]), forces[i]))

        if 'MAGMOM' in paras.keys() and len(paras['MAGMOM']) == len(sites):
            magmom = paras['MAGMOM']
        else:
            magmom = [0. for _ in range(len(sites))]

        for i in range(len(sites)):
            sites[i].set_magmom(magmom[i])

        return sites

    def getParameters(self):
        calType = self.calculationType
        kpoints = {}
        pseudopotential = []
        parameters_dict = {}
        geometry_optim = ['PREC', 'ISPIN', 'ISTART', 'ICHARG', 'ENCUT', 'NELM', 'ALGO', 'EDIFF', 'EDIFFG',
                          'ISMEAR', 'SIGMA', 'GGA', 'NBANDS', 'NEDOS', 'IBRION', 'ISIF', 'NSW', 'LORBIT',
                          'ISYM', 'LDAU', 'LHFCALC', 'LWAVE', 'LCHARG']
        StaticDosBand_calculation = ['PREC', 'ISPIN', 'ISTART', 'ICHARG', 'ENCUT', 'NELM', 'ALGO', 'EDIFF', 'EDIFFG',
                                     'ISMEAR', 'SIGMA', 'GGA', 'NBANDS', 'NEDOS', 'IBRION', 'ISIF', 'LORBIT',
                                     'ISYM', 'LDAU', 'LHFCALC', 'LWAVE', 'LCHARG', 'LELF', 'LOPTICS', 'CSHIFT',
                                     'MAGMOM', 'LAECHG']
        dielectric_calculation = ['PREC', 'ISPIN', 'ISTART', 'ICHARG', 'ENCUT', 'NELM', 'ALGO', 'EDIFF', 'EDIFFG',
                                  'ISMEAR', 'SIGMA', 'GGA', 'NBANDS', 'NEDOS', 'IBRION', 'ISIF', 'LORBIT',
                                  'ISYM', 'LWAVE', 'LCHARG', 'LOPTICS', 'CSHIFT', 'MAGMOM']
        timepredict = ['PREC', 'ISPIN', 'IBRION', 'EDIFF', 'EDIFFG', 'ISIF', 'NSW', 'ALGO', 'ENCUT', 'NELECT',
                       'LOPTICS', 'LHFCALC']
        if calType == CalType.GeometryOptimization:
            parameters_dict = self.findPara(geometry_optim)
        elif calType == CalType.StaticCalculation or CalType.DensityOfStates or CalType.BandStructure or CalType.ElasticProperties or CalType.MagneticProperties:
            parameters_dict = self.findPara(StaticDosBand_calculation)
            if parameters_dict['LDAU']:
                parameters_dict['LDAUTYPE'] = self.findPara(['LDAUTYPE'])['LDAUTYPE']
                parameters_dict['LDAUL'] = self.findPara(['LDAUL'])['LDAUL']
                parameters_dict['LDAUU'] = self.findPara(['LDAUU'])['LDAUU']
                parameters_dict['LDAUJ'] = self.findPara(['LDAUJ'])['LDAUJ']
                parameters_dict['LDAUPRINT'] = self.findPara(['LDAUPRINT'])['LDAUPRINT']
            if parameters_dict['LHFCALC']:
                parameters_dict['AEXX'] = self.findPara(['AEXX'])['AEXX']
                parameters_dict['HFSCREEN'] = self.findPara(['HFSCREEN'])['HFSCREEN']
                parameters_dict['LMAXFOCK'] = self.findPara(['LMAXFOCK'])['LMAXFOCK']
        elif calType == CalType.DielectricProperties:
            parameters_dict = self.findPara(dielectric_calculation)
        else:
            parameters_dict = self.findPara(timepredict)
        child = self.root.find("./kpoints/generation")
        # Kpoint信息： band 提取与其他不同，分情况
        if calType == CalType.BandStructure:
            if child.attrib['param'] == 'listgenerated':
                for child2 in child:
                    if child2.attrib == {'type': 'int', 'name': 'divisions'}:
                        kpoints['KgridDivision'] = child2.text.strip()
                HighsymPoints = []
                for v in child.findall('v'):
                    HighsymPoints.append([float(i) for i in v.text.split()])
                kpoints['HighsymPoints'] = HighsymPoints
            parameters_dict['KpointsPath'] = kpoints
        else:
            if child.attrib['param'] == 'Gamma':
                kpoints['GammaCentered'] = True
            else:
                kpoints['GammaCentered'] = False
            for child2 in child:
                if child2.attrib == {'type': 'int', 'name': 'divisions'}:
                    kpoints['KgridDivision'] = [int(i) for i in child2.text.split()]
                if child2.attrib == {'name': 'shift'}:
                    kpoints['Meshshift'] = [float(i) for i in child2.text.split()]
            child = self.root.find("./kpoints/varray[@name='kpointlist']")
            count = 0
            for _ in child:
                count += 1
            kpoints['totalnums'] = count
            parameters_dict['Kpoints'] = kpoints
        # 赝势信息
        child = self.root.find("./atominfo/array[@name='atomtypes']/set")
        for child2 in child:
            pseudopotential.append(child2[4].text.strip(' '))
        parameters_dict['PseudoPotential'] = pseudopotential
        if self.incar is not None:
            parameters_dict = self.incar.fill_parameters(parameters_dict)
        return parameters_dict

    def getSoftware(self):
        software = {}
        for child in self.root:
            if child.tag == 'generator':
                for child2 in child:
                    if child2.attrib == {'name': "program", 'type': "string"}:
                        software['SoftwareName'] = child2.text.strip(' ')
                    elif child2.attrib == {'name': "version", 'type': "string"}:
                        software['SoftwareVersion'] = child2.text.strip(' ')
                    elif child2.attrib == {'name': "subversion", 'type': "string"}:
                        software['Subversion'] = child2.text.strip(' ').replace("    ", " ")
                    elif child2.attrib == {'name': "platform", 'type': "string"}:
                        software['Platform'] = child2.text.strip(' ')
                return software

    def getStartTime(self):
        sdata = ''
        stime = ''
        for child in self.root:
            if child.tag == 'generator':
                for child2 in child:
                    if child2.attrib == {'name': "date", 'type': "string"}:
                        sdata = child2.text
                    elif child2.attrib == {'name': "time", 'type': "string"}:
                        stime = child2.text
                timearray = time.strptime(sdata + stime, "%Y %m %d %H:%M:%S ")
                timestyle = time.strftime("%Y-%m-%d %H:%M:%S", timearray)
                return timestyle

    def getResourceUsage(self):
        resourceUsage = self.outcar.getResourceUsage()
        return resourceUsage

    def getIonicSteps(self):
        ionicsteps = {}
        sites = self.sites_e
        sites_new = []
        lattice = None
        structures = []
        all_forces = []
        stresses = []
        energy_next = 0.0
        energys = []
        cputimes = []
        totalenergydiffs = []
        for child in self.root:
            if child.tag == 'calculation':
                for child2 in child:
                    if child2.tag == 'structure':
                        for child3 in child2:
                            if child3.tag == 'crystal':
                                a = child3[0][0].text.split()
                                a = [float(i) for i in a]
                                b = child3[0][1].text.split()
                                b = [float(i) for i in b]
                                c = child3[0][2].text.split()
                                c = [float(i) for i in c]
                                lattice = Lattice([a, b, c])
                            if child3.tag == 'varray':
                                i = 0
                                sites_new = []
                                for site in sites:
                                    coords = [float(i) for i in child3[i].text.split()]
                                    site.coords = coords
                                    sites_new.append(site)
                                    i += 1
                        composition = self.composition
                        from services.implements.vaspStructure import VaspStructure
                        structure = VaspStructure(lattice=lattice, composition=composition, sites=sites_new)
                        structure.setup()
                        structure_doc = {
                            "Formula": structure.formula,
                            "SimplestFormula": structure.simplestFormula,
                            "GeneralFormula": structure.generalFormula,  # 新增
                            "NumberOfElements": structure.numberOfElements,
                            "NumberOfSites": structure.numberOfSites,
                            "Composition": structure.composition,
                            "Stoichiometry": structure.stoichiometry,
                            "LatticeParameters": structure.latticeParameters,
                            "Lattice": structure.lattice.matrix.tolist(),
                            "Volume": structure.volume,
                            "SpaceGroup": structure.spaceGroup,
                            "PointGroup": structure.pointGroup,
                            "CrystalSystem": structure.crystalSystem,
                            "BravaisLattice": structure.bravaisLattice,
                            "ReciprocalLattice": structure.reciprocalLattice.tolist(),
                            "PrimitiveLattice": structure.primitiveLattice.tolist(),
                            "ConventionalLattice": structure.conventionalLattice.tolist(),
                            "HighSymmetryKpath": structure.highSymmetryKPath,
                            "Sites": structure.sitesDoc,
                            "Prototype": structure.prototype,
                            "CellStress": structure.cellStress,
                            'HashValue': structure.hashId
                        }
                        structures.append(structure_doc)
                    if child2.attrib == {'name': 'forces'}:
                        i = 0
                        forces = []
                        for child3 in child2:
                            force = {'Force': [float(x) for x in child3.text.split()],
                                     'Atom': {
                                         'AtomicSymbol': sites[i].atom.atomicsymbol,
                                         'AtomicNumber': sites[i].atom.atomicnumber,
                                         'AtomicMass': sites[i].atom.atomicmass
                                     }
                                     }
                            forces.append(force)
                            i += 1
                        all_forces.append(forces)
                    if child2.attrib == {'name': 'stress'}:
                        a = [float(i) for i in child2[0].text.split()]
                        b = [float(i) for i in child2[1].text.split()]
                        c = [float(i) for i in child2[2].text.split()]
                        stress = [a, b, c]
                        stresses.append(stress)
                    if child2.tag == 'energy':
                        for child3 in child2:
                            if child3.attrib == {'name': 'e_fr_energy'}:
                                energy_now = float(child3.text)
                                energys.append(energy_now)
                    if child2.tag == 'time':
                        cputime = float(child2.text.split()[1])
                        cputimes.append(cputime)
        ionicsteps['TotalEnergy'] = energys
        ionicsteps['IonStepCpuTime'] = cputimes
        ionicsteps['StepStructure'] = structures
        ionicsteps['StepForces'] = all_forces
        ionicsteps['StepStress'] = stresses

        energy_pre = 0.0
        for energy in energys:
            totalenergydiffs.append(energy - energy_pre)
            energy_pre = energy

        ionicsteps['TotalEnergyDiff'] = totalenergydiffs
        para = self.parameters['EDIFFG']
        if para >= 0:
            if totalenergydiffs[-1] <= para:
                ionicsteps['IonConvergency'] = True
            else:
                ionicsteps['IonConvergency'] = False
        else:
            ionicsteps['IonConvergency'] = True
            for line in all_forces[-1]:
                for ele in line['Force']:
                    if abs(ele) > abs(para):
                        ionicsteps['IonConvergency'] = False
        return ionicsteps

    def getElectronicSteps(self):
        electronicSteps = []
        for child in self.root:
            if child.tag == 'calculation':
                energys = []
                cputimes = []
                totalenergydiffs = []
                for child2 in child:
                    if child2.tag == 'scstep':
                        for child3 in child2:
                            if child3.attrib == {'name': 'total'}:
                                times = child3.text.split()
                                if len(times) == 1:
                                    times = times[0]
                                    times = [times[: times.find('.') + 3], times[times.find('.') + 3:]]
                                if '*' in times[1]:
                                    cputime = 'NAN'
                                else:
                                    cputime = float(times[1])
                                cputimes.append(cputime)
                            if child3.tag == 'energy':
                                for child4 in child3:
                                    if child4.attrib == {'name': 'e_fr_energy'}:
                                        energy = float(child4.text)
                                        energys.append(energy)

                energy_pre = 0.0
                for energy in energys:
                    totalenergydiffs.append(energy - energy_pre)
                    energy_pre = energy

                para = self.parameters['EDIFF']
                if para >= totalenergydiffs[-1]:
                    eleconvergency = True
                else:
                    eleconvergency = False
                doc = {
                    'TotalEnergy': energys,
                    'EleStepCpuTime': cputimes,
                    'TotalEnergyDiff': totalenergydiffs,
                    'EleConvergency': eleconvergency
                }
                electronicSteps.append(doc)
        return electronicSteps

    def getThermoDynamicProperties(self):
        child = self.root.find("./calculation[last()]/energy/i[@name='e_fr_energy']")
        totalenergy = float(child.text)
        child = self.root.find("./calculation[last()]/dos/i[@name='efermi']")
        numberofatoms = int(self.root.find("./atominfo/atoms").text)
        fermienergy = float(child.text)
        energyPerAtom = totalenergy / numberofatoms
        formation_energy = 0.0
        composition = self.composition
        energy_atoms = 0.
        for comp in composition:
            if comp['AtomicSymbol'] in PTable().atom_energy:
                energy_atoms += PTable().atom_energy[comp['AtomicSymbol']] * comp['Amount']
            else:
                energy_atoms = 0.
                break
        if energy_atoms != 0.:
            formation_energy = (totalenergy - energy_atoms) / self.numberOfSites

        doc = {
            'TotalEnergy': totalenergy,
            'FermiEnergy': fermienergy,
            'EnergyPerAtom': energyPerAtom,
            'FormationEnergy': formation_energy
        }
        return doc

    def getReciprocalLattice(self, isinit: bool = False):
        """
                获取reciprocallattice
                :param isinit:
                :return:
                """
        if isinit:
            child = self.root.find("./structure[@name='initialpos']/crystal/varray[@name='rec_basis']")
        else:
            child = self.root.find("./structure[@name='finalpos']/crystal/varray[@name='rec_basis']")
        a = child[0].text.split()
        a = [float(i) for i in a]
        b = child[1].text.split()
        b = [float(i) for i in b]
        c = child[2].text.split()
        c = [float(i) for i in c]
        return np.array([a, b, c]).reshape((3, 3))

    def findPara(self, parameters: list):
        """
        获取指定的参数
        :return:
        """
        parameters_dict = {}
        xpath1 = "./incar/i"
        xpath2 = "./incar/v"
        xpath3 = "./parameters//i"
        xpath4 = "./parameters//v"
        for para in parameters:
            flag = False
            if para == 'ALGO':
                parameters_dict[para] = 'Normal'
            for xpath in [xpath1, xpath2, xpath3, xpath4]:
                child = self.root.findall(xpath)
                for child2 in child:
                    if child2.attrib['name'] == para:
                        flag = True
                        if 'type' in child2.attrib.keys():
                            if child2.attrib['type'] == 'int':
                                parameters_dict[para] = int(child2.text.strip())
                            elif child2.attrib['type'] == 'string':
                                parameters_dict[para] = child2.text.strip(' ')
                            elif child2.attrib['type'] == 'logical':
                                if 'F' in child2.text:
                                    parameters_dict[para] = False
                                else:
                                    parameters_dict[para] = True
                        else:
                            if child2.tag == 'v':
                                parameters_dict[para] = [float(t) for t in child2.text.split()]
                            elif child2.tag == 'i':
                                parameters_dict[para] = float(child2.text.strip())
                        break
                if flag:
                    break
        if 'ENCUT' in parameters and 'ENCUT' not in parameters_dict.keys():
            if parameters_dict['PREC'] == 'normal':
                enmax = []
                composition = self.composition
                for doc in composition:
                    symbol = doc['AtomicSymbol']
                    em = PTable().enmax[symbol]
                    enmax.append(em)
                parameters_dict['ENCUT'] = max(enmax)
        return parameters_dict

    def getTotalDos(self):
        """
        提取总电子态密度数据
        @return: dict or None if no total
        """
        child = self.root.find("./calculation[last()]/dos/total/array/set")
        if child is None:
            return None
        IsSpinPolarized = False
        NumberOfGridPoints = 0
        Energies = list()
        TdosData = {}
        for s in child.findall("set"):
            spin = Spin.up if s.attrib["comment"] == "spin 1" else Spin.down
            data = np.array(parseVarray(s))
            Energies = list(data[:, 0])
            TdosData[spin] = list(data[:, 1])
            NumberOfGridPoints = len(data[:, 0])
            if spin == Spin.down:
                IsSpinPolarized = True
        return {
            "IsSpinPolarized": IsSpinPolarized,
            "NumberOfGridPoints": NumberOfGridPoints,
            "Energies": Energies,
            "TdosData": TdosData
        }

    def getPartialDos(self):
        """
        分原子轨道态密度
        :return:
        """
        child = self.root.find("./calculation[last()]/dos/partial/array")
        if child is None:
            return None
        IsSpinPolarized = False
        NumberOfGridPoints = 301
        NumberOfIons = 0  #
        DecomposedLength = 0  # 分波态密度的投影数
        IsLmDecomposed = False  # 是否计算了轨道投影
        Energies = {}
        PartialDosData = []
        LDecomposed = []
        for filed in child.findall("field"):
            if "energy" not in filed.text.strip():
                LDecomposed.append(filed.text.strip())
        DecomposedLength = len(LDecomposed)
        IsLmDecomposed = True if DecomposedLength == 9 or DecomposedLength == 16 else False
        irons = child.find("set").findall("set")
        NumberOfIons = len(irons)

        for iron in irons:
            energiesOfIron = []
            dosOfIron = {}
            for s in iron.findall("set"):
                spin = Spin.up if s.attrib["comment"] == "spin 1" else Spin.down
                data = np.array(parseVarray(s))
                energiesOfIron = list(data[:, 0])
                if spin not in Energies:
                    Energies[spin] = energiesOfIron
                NumberOfGridPoints = len(data[:, 0])
                for i in range(0, DecomposedLength):
                    if LDecomposed[i] not in dosOfIron:
                        dosOfIron[LDecomposed[i]] = {}
                    dosOfIron[LDecomposed[i]][spin] = list(data[:, i + 1])
                if spin == Spin.down:
                    IsSpinPolarized = True
            PartialDosData.append(dosOfIron)

        return {
            "IsSpinPolarized": IsSpinPolarized,
            "NumberOfGridPoints": NumberOfGridPoints,
            "NumberOfIons": NumberOfIons,
            "DecomposedLength": DecomposedLength,
            "IsLmDecomposed": IsLmDecomposed,
            "Energies": Energies,
            "PartialDosData": PartialDosData
        }

    def getKPoints(self):
        """
        提取 点
        :return:
        """
        child = self.root.find("./kpoints/varray[@name='kpointlist']")
        data = parseVarray(child)
        return data

    def getEigenValues(self):
        """
        提取本征值数据
        :return:
        """
        NumberOfGeneratedKPoints = 0
        NumberOfBand = 0
        IsSpinPolarized = False
        KPoints = None
        EigenvalData = {}
        EigenvalOcc = {}
        child = self.root.find("./calculation[last()]/dos/i[@name='efermi']")
        efermi = float(child.text)
        child = self.root.find("./calculation[last()]/eigenvalues/array/set")
        if child is None:
            return None
        KPoints = self.kPoints
        NumberOfGeneratedKPoints = len(KPoints)
        for s in child.findall('set'):
            spin = Spin.up if s.attrib["comment"] == "spin 1" else Spin.down
            data = []
            occ = []
            for k in s.findall("set"):
                t = np.array(parseVarray(k))
                data.append(list(t[:, 0]))
                occ.append(list(t[:, 1]))
            NumberOfBand = len(data[0])
            # 按能带存储
            EigenvalEnergyData = np.array(data).transpose()
            EigenvalEnergyOcc = np.array(occ).transpose()

            EigenvalData[spin] = EigenvalEnergyData
            EigenvalOcc[spin] = EigenvalEnergyOcc

            if spin == Spin.down:
                IsSpinPolarized = True
        # NumberOfBand = len(EigenvalData[Spin['up']][0])
        return {
            "NumberOfGeneratedKPoints": NumberOfGeneratedKPoints,
            "NumberOfBand": NumberOfBand,
            "IsSpinPolarized": IsSpinPolarized,
            "KPoints": KPoints,
            "FermiEnergy": efermi,
            "EigenvalData": EigenvalData,
            "EigenvalOcc": EigenvalOcc
        }

    def getGapFromBand(self, energies=None):
        """
            能带能隙
            :return
        """
        eigenval = self.eigenValues
        if eigenval is None:
            eigenval = self.getEigenValues()
        isSpin = eigenval["IsSpinPolarized"]
        number = eigenval["NumberOfGeneratedKPoints"]
        BandsNum = eigenval["NumberOfBand"]
        Energies = eigenval["EigenvalData"]
        Kpoints = eigenval["KPoints"]
        efermi = eigenval["FermiEnergy"]
        energies = []
        kpoints = []
        for t in Energies.values():
            energies.append(t)
            kpoints.append(Kpoints)  # 存储两份kpoints，用于自旋情况计算能隙标记
        energies = np.array(energies)
        # print(energies)
        energies = energies - efermi
        # print(kpoints)

        if isSpin:
            index00 = []
            index10 = []
            index01 = []
            index11 = []
            vbmE0 = []  # cbm data
            cbmE0 = []  # vbm data
            vbmE1 = []  # cbm data
            cbmE1 = []  # vbm data
            kindex0 = -1
            kindex1 = -1
            for i in range(BandsNum):
                if min(energies[0][i]) * max(energies[0][i]) < 0:
                    return {"GapFromBand": "Metal"}
                else:
                    if min(energies[0][i]) < 0.0 or min(energies[0][i]) == 0.0:  # below fermi
                        for j in range(number):
                            if energies[0][i][j] == max(energies[0][i]):
                                index00.append(j)
                                vbmE0.append(energies[0][i][j])
                                # print("the vbm band number is :", i, "kpints is :", kpoints[0][j], "energy is :",
                                #       energies[0][i][j])
                    else:  # above fermi
                        for j in range(number):
                            if energies[0][i][j] == min(energies[0][i]):
                                index01.append(j)
                                cbmE0.append(energies[0][i][j])
                                # print("the cbm band number is :", i, "kpints is :", kpoints[0][j], "energy is :",
                                #       energies[0][i][j])
                    # print(index00, index01)

                if min(energies[1][i]) * max(energies[1][i]) < 0:
                    return {"GapFromBand": "Metal"}
                else:
                    if min(energies[1][i]) < 0.0 or min(energies[1][i]) == 0.0:  # below fermi
                        for j in range(number):
                            if energies[1][i][j] == max(energies[1][i]):
                                index10.append(j)
                                vbmE1.append(energies[1][i][j])
                                # print("the vbm band number is :", i, "kpints is :", kpoints[1][j], "energy is :",
                                #       energies[1][i][j])
                    else:  # above fermi
                        for j in range(number):
                            if energies[1][i][j] == min(energies[1][i]):
                                index11.append(j)
                                cbmE1.append(energies[1][i][j])
                                # print("the cbm band number is :", i, "kpints is :", kpoints[1][j], "energy is :",
                                #       energies[1][i][j])

            if min(cbmE0) < min(cbmE1):
                cbmEnergy = min(cbmE0)
                kindex1 = cbmE0.index(cbmEnergy)
            else:
                cbmEnergy = min(cbmE1)
                kindex1 = cbmE1.index(cbmEnergy)
            if max(vbmE0) < max(vbmE1):
                vbmEnergy = max(vbmE1)
                kindex0 = vbmE1.index(vbmEnergy)
            else:
                vbmEnergy = max(vbmE0)
                kindex0 = vbmE0.index(vbmEnergy)

            gapValue = cbmEnergy - vbmEnergy
            if kindex0 == kindex1:
                gaptype = "direct gap"
            else:
                gaptype = "indirect gap"
            return {
                "GapFromBand": gapValue,
                "GapType": gaptype,
                "VBMFromBand": vbmEnergy,
                "CBMFromBand": cbmEnergy
            }
        else:
            index00 = []
            index01 = []
            vbmE0 = []  # cbm data
            cbmE0 = []  # vbm data
            kindex0 = -1
            kindex1 = -1
            for i in range(BandsNum):
                if min(energies[0][i]) * max(energies[0][i]) < 0:
                    return {"GapFromBand": "Metal"}
                else:
                    if min(energies[0][i]) < 0.0 or min(energies[0][i]) == 0.0:  # below fermi
                        for j in range(number):
                            if energies[0][i][j] == max(energies[0][i]):
                                index00.append(j)
                                vbmE0.append(energies[0][i][j])
                                # print("the vbm band number is :", i, "kpints is :", kpoints[0][j], "energy is :",
                                #       energies[0][i][j])
                    else:  # above fermi
                        for j in range(number):
                            if energies[0][i][j] == min(energies[0][i]):
                                index01.append(j)
                                cbmE0.append(energies[0][i][j])
                                # print("the cbm band number is :", i, "kpints is :", kpoints[0][j], "energy is :",
                                #       energies[0][i][j])
                    # print(index00, index01)

            cbmEnergy = min(cbmE0)
            kindex1 = cbmE0.index(cbmEnergy)
            vbmEnergy = max(vbmE0)
            kindex0 = vbmE0.index(vbmEnergy)
            gapValue = cbmEnergy - vbmEnergy
            if kindex0 == kindex1:
                gaptype = "direct gap"
            else:
                gaptype = "indirect gap"
            return {
                "GapFromBand": gapValue,
                "GapType": gaptype,
                "VBMFromBand": vbmEnergy,
                "CBMFromBand": cbmEnergy
            }

    def getProjectedEigenvalOnIonOrbitals(self):
        """
        提取分原子能带投影
        :return:
        """
        NumberOfGeneratedKPoints = 0
        NumberOfBand = 0
        IsSpinPolarized = False
        NumberOfIons = 0
        DecomposedLength = 0
        IsLmDecomposed = False
        KPoints = None
        Data = {}

        child = self.root.find("./calculation[last()]/projected/array")
        if child is None:
            return None
        KPoints = self.kPoints
        NumberOfGeneratedKPoints = len(KPoints)

        fields = []
        for field in child.findall('field'):
            fields.append(field.text.strip())
        DecomposedLength = len(fields)
        IsLmDecomposed = True if DecomposedLength == 9 or DecomposedLength == 16 else False

        # 获取原数据格式
        data = []
        for s in child.find('set').findall('set'):  # spin
            spin = Spin.up if s.attrib["comment"] == "spin 1" else Spin.down
            spins = []
            for i, kpoint in enumerate(s.findall('set')):
                kpoints = []
                for j, band in enumerate(kpoint.findall('set')):
                    irons = parseVarray(band)
                    kpoints.append(irons)
                spins.append(kpoints)
            data.append(spins)
        # print(data)
        # 数据格式变换
        NumberOfIons = len(data[0][0][0])
        NumberOfBand = len(data[0][0])
        IsSpinPolarized = True if len(data) == 2 else False
        # spin -> kpoint -> band -> iron -> orbital
        # spin-> iron -> kpoint -> band -> orbital
        # 调整维度
        data_array = np.array(data)
        data_trans = np.transpose(data_array, (0, 3, 1, 2, 4))
        data = data_trans.tolist()

        # 格式保存
        for s in range(len(data)):
            spindata = []
            for i in range(NumberOfIons):
                ironsdata = []
                for j in range(NumberOfGeneratedKPoints):
                    points = []
                    for k in range(NumberOfBand):
                        bands = {}
                        for l in range(DecomposedLength):
                            bands[fields[l]] = data[s][i][j][k][l]
                        points.append(bands)
                    ironsdata.append(points)
                spindata.append(ironsdata)
            if s == 0:
                Data[Spin.up] = spindata
            elif s == 1:
                Data[Spin.down] = spindata
        return {
            "NumberOfGeneratedKPoints": NumberOfGeneratedKPoints,
            "NumberOfBand": NumberOfBand,
            "IsSpinPolarized": IsSpinPolarized,
            "NumberOfIons": NumberOfIons,
            "Decomposed": fields,
            "DecomposedLength": DecomposedLength,
            "IsLmDecomposed": IsLmDecomposed,
            "KPoints": KPoints,
            "Data": Data
        }

    def getDielectricData(self):
        """
            提取与频率相关的介电函数数据
            :return:
        """
        # dielectricData = {}
        energy = []
        imag_part = []
        real_part = []
        imag = []
        real = []
        for event, elem in ET.iterparse(self.filePath):
            tag = elem.tag
            if tag == "dielectricfunction":
                imag = [[float(l) for l in r.text.split()] for r in
                        elem.find("imag").find("array").find("set").findall("r")]
                real = [[float(l) for l in r.text.split()] for r in
                        elem.find("real").find("array").find("set").findall("r")]
                break  # break 多个 dielectricfunction 只取得第一个die..
        for e in imag:
            energy.append(e[0])
            imag_part.append(e[1:])
        for e in real:
            real_part.append(e[1:])
        dielectricData = {
            'Energy': energy,
            'real_part': real_part,
            'imag_part': imag_part
        }
        return dielectricData

    def getGapFromDos(self, energies=None):
        """
            态密度能隙
            :return
        """
        child = self.root.find("./calculation[last()]/dos/i[@name='efermi']")
        efermi = float(child.text)
        tdos = self.totalDos
        isSpin = tdos["IsSpinPolarized"]
        number = tdos["NumberOfGridPoints"]
        Energies = tdos["Energies"]
        TdosData = tdos["TdosData"]
        energies = []
        tdosdata = []
        for t in TdosData.values():
            tdosdata.append(t)
            energies.append(Energies)  # 存储两份energy，用于自旋情况计算能隙标记
        energies = np.array(energies)
        tdosdata = np.array(tdosdata)
        energies = energies - efermi
        epcEnergy = (energies[0][-1] - energies[0][0]) / number
        epcDOS = 0.1

        if isSpin:
            index00 = -1
            index10 = -1
            index01 = index00
            index11 = index10
            for i in range(number):
                if -epcEnergy <= energies[0][i] <= epcEnergy:
                    if tdosdata[0][i] <= epcDOS:  # tdosdata数据库中存储均为非负数，判断是否为零时只需确定正方向误差即可。
                        index00 = i
                        break
                    # print("tdosdata<= epcDOS", tdosdata[0][index00])
            for i in range(number):
                if -epcEnergy <= energies[1][i] <= epcEnergy:
                    if tdosdata[1][i] <= epcDOS:
                        index10 = i
                        break
            if index00 == -1 or index10 == -1:
                return {"GapFromDOS": "Metal"}
            else:
                # print(index00, index10, index01, index11)
                #   for i in range(number):  #冗余代码
                #      print("i am here===================",i)
                #      if -epcEnergy <= energies[0][i] <= epcEnergy:
                #           index00 = i
                #          break
                for i in range(index00, len(energies[0])):
                    # print(i)
                    if tdosdata[0][i] > epcDOS:
                        index01 = i
                        break
                # for i in range(number): #冗余代码
                #    if -epcEnergy <= energies[1][i] <= epcEnergy:
                #        index1 = i
                #        break
                for i in range(index10, len(energies[0])):
                    # print(i)
                    if tdosdata[1][i] > epcDOS:
                        index11 = i
                        break
                if (energies[0][index01] < energies[1][index10]):
                    return {"GapFromDOS": "Metal"}
                elif (energies[1][index11] < energies[0][index10]):
                    return {"GapFromDOS": "Metal"}
                else:
                    gapValue = min(energies[0][index01], energies[1][index11]) - max(energies[0][index00],
                                                                                     energies[1][index10])
                    # print(index01, index11, index10, index00)
                    return {
                        "GapFromDOS": gapValue,
                        "VBMfromDOS": max(energies[0][index00], energies[1][index10]),
                        "CBMfromDOS": min(energies[0][index01], energies[1][index11])
                    }
        else:
            index0 = -1
            index1 = -1
            for i in range(number):
                if -epcEnergy <= energies[0][i] <= epcEnergy:
                    if tdosdata[0][i] <= epcDOS:
                        index0 = i
                        break
            if index0 == -1:
                return {"GapFromDOS": "Metal"}
            else:
                # for i in range(number):
                #    if -epcEnergy <= energies[0][i] <= epcEnergy:
                #        index0 = i
                #        break
                for i in range(index0, len(energies[0])):
                    if tdosdata[0][i] > epcDOS:
                        index1 = i
                        break
                gapValue = energies[0][index1] - energies[0][index0],
                # print(index0, index1)
                return {
                    "GapFromDOS": gapValue,
                    "VBMfromDOS": energies[0][index0],
                    "CBMfromDOS": energies[0][index1]
                }

    def getOpticalProperties(self):
        """
        提取介电函数信息
        :return:
        """
        data = self.dielectricData
        energy = data['Energy']
        numOfEnergyPoints = len(energy)
        real_part = data['real_part']
        imag_part = data['imag_part']
        # dielectric data are complex values
        dielectricDataValue = complex(real_part, imag_part)
        dielectricDataName = ["εxx", "εyy", "εzz", "εxy", "εyz", "εzx"]
        xxName = "εxx"
        yyName = "εyy"
        zzName = "εzz"
        xyName = "εxy"
        yzName = "εyz"
        zxName = "εzx"

        DieDataRealxx = []
        DieDataRealyy = []
        DieDataRealzz = []
        DieDataRealxy = []
        DieDataRealyz = []
        DieDataRealzx = []
        DieDataImagxx = []
        DieDataImagyy = []
        DieDataImagzz = []
        DieDataImagxy = []
        DieDataImagyz = []
        DieDataImagzx = []

        for i in range(numOfEnergyPoints):
            DieDataRealxx[i] = real_part[i, 0]
            DieDataRealyy[i] = real_part[i, 1]
            DieDataRealzz[i] = real_part[i, 2]
            DieDataRealxy[i] = real_part[i, 3]
            DieDataRealyz[i] = real_part[i, 4]
            DieDataRealzx[i] = real_part[i, 5]

            DieDataImagxx[i] = imag_part[i, 0]
            DieDataImagyy[i] = imag_part[i, 1]
            DieDataImagzz[i] = imag_part[i, 2]
            DieDataImagxy[i] = imag_part[i, 3]
            DieDataImagyz[i] = imag_part[i, 4]
            DieDataImagzx[i] = imag_part[i, 5]

        n_w_data = []  # 折射率 （refractive index）
        l_w_data = []  # 能量损失谱（energy-loss spectrum/function）
        k_w_data = []  # 消光系数（extinction index）
        alpha_data = []  # 吸收系数（adsorption coefficient）
        r_w_data = []  # 反射率（reflectivity/reflection coefficient）
        sigma_real_data = []  # 光电导率（optical conductivity）
        sigma_imag_data = []

        # 定义参数
        c = 2.9979 * math.pow(10, 8)
        h = 6.62607 * math.pow(10, -34)
        e = 1.602176 * math.pow(10, -19)
        epsilon0 = 8.8542 * math.pow(10, -12)
        eps = 0.000001
        for i in range(numOfEnergyPoints):
            n_w_group = []
            l_w_group = []
            k_w_group = []
            alpha_group = []
            r_w_group = []
            sigma_real_group = []
            sigma_imag_group = []

            for j in range(6):
                real = real_part[i][j]
                imag = imag_part[i][j]
                n_w = math.sqrt(math.sqrt(real * real + imag * imag) + real) / math.sqrt(2)
                l_w = imag / (real * real + imag * imag) if real != 0 and imag != 0 else 0.
                k_w = math.sqrt(math.sqrt(real * real + imag * imag) - real) / math.sqrt(2)
                alpha = 4 * math.pi * energy[i] * e * k_w / (h * c)
                r_w = ((n_w - 1) ** 2 + k_w ** 2) / ((n_w + 1) ** 2 + k_w ** 2)
                w = 2 * math.pi * energy[i] * e / h
                sigma = -epsilon0 * w * (real - 1 + imag * 1j) * 1j
                sigma_real = sigma.real
                sigma_imag = sigma.imag
                n_w_group.append(n_w)
                l_w_group.append(l_w)
                k_w_group.append(k_w)
                alpha_group.append(alpha)
                r_w_group.append(r_w)
                sigma_real_group.append(sigma_real)
                sigma_imag_group.append(sigma_imag)
            n_w_data.append(n_w_group)
            l_w_data.append(l_w_group)
            k_w_data.append(k_w_group)
            alpha_data.append(alpha_group)
            r_w_data.append(r_w_group)
            sigma_real_data.append(sigma_real_group)
            sigma_imag_data.append(sigma_imag_group)
            # print(sigma_real_data)
        opticalPropereties = {
            "NumberofEnergyPoints": numOfEnergyPoints,
            "Energies": energy,
            "DielectricFunctions": {
                "RealPart": real_part,
                'ImaginaryPart': imag_part
            },
            "RefractiveIndex": {
                "Data": n_w_data
            },
            "Energy_lossSpectrumFunction": {
                "Data": l_w_data
            },
            "ExtinctionIndex": {
                "Data": k_w_data
            },
            "AdsorptionCoefficient": {
                "Data": alpha_data
            },
            "ReflectivityCoefficient": {
                "Data": r_w_data
            },
            "OpticalConductivity": {
                "RealPart": sigma_real_data,
                "ImaginaryPart": sigma_imag_data
            }
        }
        return opticalPropereties

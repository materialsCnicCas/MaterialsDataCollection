import json

from services.fileParser import FileParser
from utils.assist.helper import findValue, exists
from utils.tools.calType import CalType
from utils.tools.lattice import Lattice
from utils.tools.periodicTable import PTable
from utils.tools.pointGroupInfo import pointgroup_symbol_num
from utils.tools.sites import Site, Atom


class MPParser(FileParser):

    def __init__(self, filePath):
        super().__init__(filePath)
        try:
            with open(self.filePath, 'r', encoding='utf8') as fp:
                self.data = json.load(fp)["data"][0]
                fp.close()
        except KeyError:
            self.data = None
            fp.close()
        self.calculationType = None
        self.simplestFormula = None
        self.generalFormula = None
        self.numberOfElements = None
        self.spaceGroup = None
        self.pointGroup = None
        self.crystalSystem = None
        self.cellStress_s = None
        self.cellStress_e = None
        self.electronicProperties = None
        self.magneticProperties = None

    def setup(self):
        if not (exists(self.data, ['output']) and exists(self.data,['input', 'parameters'])):
            raise ValueError('文件信息缺失！')
        self.calculationType = self.getCalculationType()
        if self.calculationType is None:
            raise ValueError(self.filePath + '文件数据未匹配合适的计算类型')
        self.lattice_s, self.latticeParameters_s = self.getLatticeParameters(isinit=True)
        self.lattice_e, self.latticeParameters_e = self.getLatticeParameters(isinit=False)
        self.composition = self.getComposition()
        self.parameters = self.getParameters()
        self.sites_s = self.getSites(isinit=True)
        self.sites_e = self.getSites(isinit=False)
        self.volume_s = self.getVolume(isinit=True)
        self.volume_e = self.getVolume(isinit=False)
        self.numberOfSites = self.getNumberOfSites()
        self.simplestFormula = self.getSimplestFormula()
        self.generalFormula = self.getGeneralFormula()
        self.numberOfElements = self.getNumberOfElements()
        self.spaceGroup = self.getSpaceGroup()
        self.pointGroup = self.getPointGroup()
        self.crystalSystem = self.getCrystalSystem()
        self.startTime = self.getStartTime()
        self.software = self.getSoftware()
        self.ionicSteps = self.getIonicSteps()
        self.electronicSteps = self.getElectronicSteps()
        self.thermoDynamicProperties = self.getThermoDynamicProperties()
        self.cellStress_s = self.getCellStress(isinit=True)
        self.cellStress_e = self.getCellStress(isinit=False)
        self.electronicProperties = self.getElectroicProperties()
        self.magneticProperties = self.getMagneticProperties()
        self.hasSetup = True

    def getGeneralFormula(self):
        return self.data['formula_anonymous']

    def getCalculationType(self):
        # step 1:
        if "task_type" in self.data.keys():
            if 'static' in self.data["task_type"].lower() or 'nscf uniform' in self.data["task_type"].lower():
                return CalType.StaticCalculation
            elif 'optimization' in self.data['task_type'].lower():
                return CalType.GeometryOptimization
            elif 'nscf line' in self.data['task_type'].lower():
                return CalType.BandStructure
        if "task_label" in self.data.keys():
            if 'static' in self.data["task_label"].lower() or 'nscf uniform' in self.data["task_label"].lower():
                return CalType.StaticCalculation
            elif 'optimization' in self.data['task_label'].lower():
                return CalType.GeometryOptimization
            elif 'nscf line' in self.data['task_label'].lower():
                return CalType.BandStructure
        if 'kpoints' in self.data['orig_inputs'].keys() and "Non SCF run along symmetry lines" in \
                self.data['orig_inputs']['kpoints']['comment']:
            return CalType.BandStructure

        # step 2:
        if 'parameters' not in self.data['input']:
            raise ValueError('文件缺少参数，无法判断提取类型')
        parameters = self.data['input']['parameters']
        if parameters['IBRION'] == 1 or parameters['IBRION'] == 2 or parameters['IBRION'] == 3:
            return CalType.GeometryOptimization
        if parameters['IBRION'] == -1 and parameters['LORBIT']:
            return CalType.DensityOfStates
        if parameters['IBRION'] == -1 and self.data["orig_inputs"]['kpoints'][
            'comment'] == 'Non SCF run along symmetry lines':
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
        compositionFull = []
        sites = self.data['input']['structure']['sites']
        for site in sites:
            spec = site['label']
            composition[spec] = composition[spec] + 1 if spec in composition else 1

        for spec in composition.keys():
            atom = {
                "AtomicSymbol": spec,
                "AtomicNumber": PTable().detailed[spec]["Atomic no"],
                "AtomicMass": PTable().detailed[spec]["Atomic mass"],
                "Amount": composition[spec]
            }
            compositionFull.append(atom)
        return compositionFull

    def getSites(self, isinit: bool = False):
        Sites = []
        if isinit:
            sites = self.data['input']['structure']['sites']
        else:
            sites = self.data['output']['structure']['sites']
            forces = self.data['output']['forces'] if 'forces' in self.data['output'] else [None for _ in
                                                                                            range(len(sites))]
        paras = self.parameters
        if 'MAGMOM' in paras.keys() and len(paras['MAGMOM']) == len(sites):
            magmom = paras['MAGMOM']
        else:
            magmom = [0. for _ in range(len(sites))]
        for i in range(len(sites)):
            pos = sites[i]['abc']
            spec = sites[i]['label']
            if isinit:
                site = Site(pos, Atom(spec))
            else:
                site = Site(pos, Atom(spec), forces[i])
            site.set_occupancy(sites[i]['species'][0]['occu'])
            site.set_magmom(magmom[i])
            Sites.append(site)
        return Sites

    def getVolume(self, isinit: bool = False):
        if isinit:
            volume = self.data['input']['structure']['lattice']['volume']
        else:
            volume = self.data['output']['structure']['lattice']['volume']
        return volume

    def getNumberOfSites(self):
        return self.data["nsites"]

    def getLatticeParameters(self, isinit: bool = False):
        if isinit:
            lattice = self.data['input']['structure']['lattice']
        else:
            lattice = self.data['output']['structure']['lattice']
        a, b, c, alpha, beta, gamma = lattice['a'], lattice['b'], lattice['c'], lattice['alpha'], lattice['beta'], \
                                      lattice['gamma']
        lattice_parameters = {
            'a': a,
            'b': b,
            'c': c,
            'alpha': alpha,
            'beta': beta,
            'gamma': gamma
        }
        return Lattice(lattice['matrix']), lattice_parameters

    def getSimplestFormula(self):
        return self.data['formula_pretty']

    def getNumberOfElements(self):
        return self.data['nelements']

    def getParameters(self):
        cal_type = self.calculationType
        geometry_optim = ['PREC', 'ISPIN', 'ISTART', 'ICHARG', 'ENCUT', 'NELM', 'ALGO', 'EDIFF', 'EDIFFG',
                          'ISMEAR', 'SIGMA', 'GGA', 'NBANDS', 'NEDOS', 'IBRION', 'ISIF', 'NSW', 'LORBIT',
                          'ISYM', 'LDAU', 'LHFCALC', 'LWAVE', 'LCHARG']
        other_calculation = ['PREC', 'ISPIN', 'ISTART', 'ICHARG', 'ENCUT', 'NELM', 'ALGO', 'EDIFF', 'EDIFFG',
                             'ISMEAR', 'SIGMA', 'GGA', 'NBANDS', 'NEDOS', 'IBRION', 'ISIF', 'LORBIT',
                             'ISYM', 'LDAU', 'LHFCALC', 'LWAVE', 'LCHARG', 'LELF', 'LOPTICS', 'CSHIFT',
                             'MAGMOM', 'LAECHG']
        dielectric_calculation = ['PREC', 'ISPIN', 'ISTART', 'ICHARG', 'ENCUT', 'NELM', 'ALGO', 'EDIFF', 'EDIFFG',
                                  'ISMEAR', 'SIGMA', 'GGA', 'NBANDS', 'NEDOS', 'IBRION', 'ISIF', 'LORBIT',
                                  'ISYM', 'LWAVE', 'LCHARG', 'LOPTICS', 'CSHIFT', 'MAGMOM']
        if 'incar' not in self.data['input']:
            raise ValueError('文件信息不全，无法提取')
        incar = self.data['input']['incar']
        parameters = self.data['input']['parameters']
        parameters_dict = incar
        kpoints = {}

        if cal_type == CalType.GeometryOptimization:
            for key in geometry_optim:
                if key not in parameters_dict and key in parameters.keys():
                    parameters_dict[key] = parameters[key]
                elif key not in parameters_dict and key not in parameters:
                    parameters_dict[key] = None

        else:
            for key in other_calculation:
                if key not in parameters_dict and key in parameters.keys():
                    parameters_dict[key] = parameters[key]
                elif key not in parameters_dict and key not in parameters:
                    parameters_dict[key] = None
            if parameters_dict['LDAU']:
                parameters_dict['LDAUTYPE'] = parameters['LDAUTYPE'] if 'LDAUTYPE' in parameters.keys() else None
                parameters_dict['LDAUL'] = parameters['LDAUL'] if 'LDAUL' in parameters.keys() else None
                parameters_dict['LDAUU'] = parameters['LDAUU'] if 'LDAUU' in parameters.keys() else None
                parameters_dict['LDAUJ'] = parameters['LDAUJ'] if 'LDAUJ' in parameters.keys() else None
                parameters_dict['LDAUPRINT'] = parameters['LDAUPRINT'] if 'LDAUPRINT' in parameters.keys() else None
            if parameters_dict['LHFCALC']:
                parameters_dict['AEXX'] = parameters['AEXX'] if 'AEXX' in parameters.keys() else None
                parameters_dict['HFSCREEN'] = parameters['HFSCREEN'] if 'HFSCREEN' in parameters.keys() else None
                parameters_dict['LMAXFOCK'] = parameters['LMAXFOCK'] if 'LMAXFOCK' in parameters.keys() else None

        # Kpoint信息： band 提取与其他不同，分情况
        if 'kpoints' in self.data['orig_inputs']:
            if self.data['orig_inputs']['kpoints']['generation_style'] == 'Gamma':
                kpoints['GammaCentered'] = True
            else:
                kpoints['GammaCentered'] = False
            kpoints['totalnums'] = self.data['orig_inputs']['kpoints']['nkpoints']
            kpoints['MeshShift'] = self.data['orig_inputs']['kpoints']['usershift']
            kpoints['KgridDivision'] = None

        parameters_dict['kpoints'] = kpoints
        # 赝势信息
        potcar_spec = self.data['input']['potcar_spec']
        pseudopotential = []
        for pot in potcar_spec:
            pseudopotential.append(pot['titel'])
        parameters_dict['PseudoPotential'] = pseudopotential
        return parameters_dict

    def getSoftware(self):
        softwareName = 'vasp'
        if 'calcs_reversed' in self.data.keys() and len(self.data['calcs_reversed']) > 0:
            softwareVersion = self.data['calcs_reversed'][0]['vasp_version']
        else:
            softwareVersion = 'Unknown'
        return {
            'SoftwareName': softwareName,
            'SoftwareVersion': softwareVersion
        }

    def getStartTime(self):
        if 'last_updated' in self.data.keys():
            time = self.data['last_updated']
        else:
            time = 'UnKnown'
        return time

    def getIonicSteps(self):
        ionicsteps = {}
        sites = self.sites_s

        lattice = None
        structures = []
        all_forces = []
        stresses = []
        energy_pre = 0.0
        energys = []
        cpu_times = []
        energy_diffs = []

        if 'calcs_reversed' in self.data.keys():
            items = self.data['calcs_reversed'][0]['output']['ionic_steps']
            for item in items:
                sites_new = []
                # structures
                lattice_all = item['structure']['lattice']
                lattice = Lattice(lattice_all['matrix'])
                a, b, c, alpha, beta, gamma = lattice_all['a'], lattice_all['b'], lattice_all['c'], lattice_all[
                    'alpha'], lattice_all['beta'], lattice_all['gamma']
                lattice_parameters = {
                    'a': a,
                    'b': b,
                    'c': c,
                    'alpha': alpha,
                    'beta': beta,
                    'gamma': gamma
                }
                composition = self.composition
                sites_local = item['structure']['sites']
                for i in range(len(sites)):
                    coords = sites_local[i]['abc']
                    sites[i].coords = coords
                    sites[i].force = item['forces'][i]
                    sites_new.append(sites[i])

                from services.implements.mpStructure import MPStructure
                structure = MPStructure(lattice=lattice, composition=composition, sites=sites_new)
                structure.setup()
                structure_doc = {
                    "Formula": structure.formula,
                    "SimplestFormula": self.simplestFormula,
                    "GeneralFormula": self.generalFormula,  # 新增
                    "NumberOfElements": self.numberOfElements,
                    "NumberOfSites": self.numberOfSites,
                    "Composition": structure.composition,
                    "Stoichiometry": structure.stoichiometry,
                    "LatticeParameters": lattice_parameters,
                    "Lattice": structure.lattice.matrix.tolist(),
                    "Volume": lattice_all['volume'],
                    "SpaceGroup": structure.spaceGroup,
                    "PointGroup": structure.pointGroup,
                    "CrystalSystem": structure.crystalSystem,
                    "BravaisLattice": structure.bravaisLattice,
                    "ReciprocalLattice": structure.reciprocalLattice.tolist(),
                    "PrimitiveLattice": structure.primitiveLattice.tolist(),
                    "ConventionalLattice": structure.conventionalLattice.tolist(),
                    "HighSymmetryKpath": structure.highSymmetryKPath,
                    "Sites": structure.sitesDoc,  # 有变化
                    "Prototype": self.getProtoType(),  # 新增
                    "CellStress": item['stress'],  # 新增
                    'HashValue': structure.hashId
                }
                structures.append(structure_doc)

                cpu_times.append(None)
                energys.append(item['e_fr_energy'])
                stresses.append(item['stress'])
                forces = []
                for i in range(len(item['forces'])):
                    force = {'Force': item['forces'][i],
                             'Atom': {
                                 'AtomicSymbol': sites[i].atom.atomicsymbol,
                                 'AtomicNumber': sites[i].atom.atomicnumber,
                                 'AtomicMass': sites[i].atom.atomicmass
                             }
                             }
                    forces.append(force)
                all_forces.append(forces)
            ionicsteps['TotalEnergy'] = energys
            ionicsteps['IonStepCpuTime'] = cpu_times
            ionicsteps['StepStructure'] = structures
            ionicsteps['StepForces'] = all_forces
            ionicsteps['StepStress'] = stresses

            for energy in energys:
                energy_diffs.append(energy - energy_pre)
                energy_pre = energy
            ionicsteps['TotalEnergyDiff'] = energy_diffs
            ediffg = self.parameters['EDIFFG']
            if ediffg >= 0:
                if energy_diffs[-1] <= ediffg:
                    ionicsteps['IonConvergency'] = True
                else:
                    ionicsteps['IonConvergency'] = False
            else:
                ionicsteps['IonConvergency'] = True
                for line in all_forces[-1]:
                    for ele in line['Force']:
                        if abs(ele) > abs(ediffg):
                            ionicsteps['IonConvergency'] = False
        return ionicsteps

    def getElectronicSteps(self):
        electronicSteps = []
        if 'calcs_reversed' in self.data.keys():
            items = self.data['calcs_reversed'][0]['output']['ionic_steps']
            for item in items:
                energys = []  # TotalEnergy
                energy_diff = []  # TotalEnergyDiff
                cpu_time = []  # EleStepCpuTime
                e_steps = item['electronic_steps']
                for e_step in e_steps:
                    energys.append(e_step['e_fr_energy'])
                    cpu_time.append(None)
                energy_pre = 0.
                for energy in energys:
                    energy_diff.append(energy - energy_pre)
                    energy_pre = energy
                ediff = self.parameters['EDIFF']
                if ediff >= energy_diff[-1]:
                    ele_convergency = True
                else:
                    ele_convergency = False

                doc = {
                    'TotalEnergy': energys,
                    'EleStepCpuTime': cpu_time,
                    'TotalEnergyDiff': energy_diff,
                    'EleConvergency': ele_convergency
                }
                electronicSteps.append(doc)
            return electronicSteps

    def getThermoDynamicProperties(self):
        if 'output' not in self.data:
            raise ValueError('信息不全，无output数据')
        total_energy = self.data['output']['energy']
        energy_per_atom = self.data['output']['energy_per_atom']
        fermi_energy = None
        if 'calcs_reversed' in self.data.keys():
            fermi_energy = self.data['calcs_reversed'][0]['output']['efermi']
        # 计算形成能
        formation_energy = 0.
        composition = self.composition
        if composition is None:
            raise ValueError('先计算composition')
        energy_atoms = 0.

        for comp in composition:
            if comp['AtomicSymbol'] in PTable().atom_energy:
                energy_atoms += PTable().atom_energy[comp['AtomicSymbol']] * comp['Amount']
            else:
                energy_atoms = 0.
                break
        if energy_atoms != 0.:
            formation_energy = (total_energy - energy_atoms) / self.numberOfSites
        return {
            'TotalEnergy': total_energy,
            'EnergyPerAtom': energy_per_atom,
            'FermiEnergy': fermi_energy,
            'FormationEnergy': formation_energy
        }

    def getSpaceGroup(self):
        if 'output' not in self.data:
            raise ValueError('信息不全，无output数据')
        spg = self.data['output']['spacegroup']
        spaceGroup = {
            "spacegroupSymbol": spg['symbol'],
            "spacegroupInternationalSymbol": None,
            "SchoenfliesSymbol": None,
            "spacegroupNumber": spg['number'],
            "hallNumber": None,
            "hallSymbol": spg['hall']
        }
        return spaceGroup

    def getPointGroup(self):
        if 'output' not in self.data:
            raise ValueError('信息不全，无output数据')
        pointgroup_symbol = self.data['output']['spacegroup']['point_group']
        pointgroup = {
            "pointgroupSymbol": pointgroup_symbol,
            "pointgroupNumber": pointgroup_symbol_num[pointgroup_symbol]
        }
        return pointgroup

    def getCrystalSystem(self):
        return self.data['output']['spacegroup']['crystal_system']

    def getProtoType(self):
        if 'proto_type' in self.data.keys():
            return self.data['proto_type']
        else:
            return None

    def getCellStress(self, isinit: bool = False):
        if isinit:
            return None
        else:
            if 'stress' in self.data['output']:
                return self.data['output']['stress']
            return None

    def getElectroicProperties(self):
        atomic_charge = findValue(self.data,['calcs_reversed', 0, 'output', 'outcar', 'charge'])
        cal_type = self.calculationType

        if cal_type == CalType.GeometryOptimization:
            GapFromType = ('GapFromOptimization', 'GapFromGeo', 'VBMFromGeo', 'CBMFromGeo')
        elif cal_type == CalType.StaticCalculation:
            GapFromType = ('GapFromStatic', 'GapFromStatic', 'VBMFromStatic', 'CBMFromStatic')
        elif cal_type == CalType.BandStructure:
            GapFromType = ('GapFromBand', 'GapFromBand', 'VBMFromBand', 'CBMFromBand')
        elif cal_type == CalType.DensityOfStates:
            GapFromType = ('GapFromDOS', 'GapFromDOS', 'VBMFromDOS', 'CBMFromDOS')
        elif cal_type == CalType.DielectricProperties:
            GapFromType = ('GapFromDielectric', 'GapFromDielectric', 'VBMFromDielectric', 'CBMFromDielectric')
        else:
            GapFromType = ('GapFromXXX', 'GapFromXXX', 'VBMFromXXX', 'CBMFromXXX')
        metal = self.data['output']['is_metal'] if 'is_metal' in self.data['output'] else False
        if metal:
            if cal_type == CalType.BandStructure:
                return {
                    'AtomicCharge': atomic_charge,
                    'Eigenvalues': None,
                    'ProjectedEigenval_on_Elements': None,
                    'ProjectedEigenval_on_ElementsOrbitals': None,
                    GapFromType[0]: {
                        GapFromType[1]: 'Metal',
                        'GapType': None,
                        GapFromType[2]: None,
                        GapFromType[3]: None
                    }
                }
            elif cal_type == CalType.DensityOfStates:
                return {
                    'AtomicCharge': atomic_charge,
                    'TotalDOS': None,
                    'PartialDOS': None,
                    GapFromType[0]: {
                        GapFromType[1]: 'Metal',
                        'GapType': None,
                        GapFromType[2]: None,
                        GapFromType[3]: None
                    }
                }
            else:
                return {
                    'AtomicCharge': atomic_charge,
                    GapFromType[0]: {
                        GapFromType[1]: 'Metal',
                        'GapType': None,
                        GapFromType[2]: None,
                        GapFromType[3]: None
                    }
                }

        gap_value = self.data['output']['bandgap']
        gap_type = 'direct gap' if self.data['output']['is_gap_direct'] else 'indirect gap'

        vbm = self.data['output']['vbm']
        cbm = self.data['output']['cbm']
        if cal_type == CalType.BandStructure:
            return {
                'AtomicCharge': atomic_charge,
                'Eigenvalues': None,
                'ProjectedEigenval_on_Elements': None,
                'ProjectedEigenval_on_ElementsOrbitals': None,
                GapFromType[0]: {
                    GapFromType[1]: gap_value,
                    'GapType': gap_type,
                    GapFromType[2]: vbm,
                    GapFromType[3]: cbm
                }
            }
        else:
            return {
                'AtomicCharge': atomic_charge,
                GapFromType[0]: {
                    GapFromType[1]: gap_value,
                    'GapType': gap_type,
                    GapFromType[2]: vbm,
                    GapFromType[3]: cbm
                }
            }

    def getMagneticProperties(self):
        magnetization = findValue(self.data, ['calcs_reversed', 0, 'output', 'outcar', 'magnetization'])
        total_magnetization = findValue(self.data, ['calcs_reversed', 0, 'output', 'outcar', 'total_magnetization'])
        return {
            'AtomicMagnetization': magnetization,
            'TotalMagnetization': total_magnetization
        }

    def getPiezoelectricProperties(self):
        epsilon_static = findValue(self.data, ['output', 'epsilon_static'])
        epsilon_static_wolfe = findValue(self.data, ['output', 'epsilon_static_wolfe'])
        born = findValue(self.data, ['calcs_reversed', 0, 'output', 'outcar', 'born'])
        piezo_tensor = findValue(self.data, ['calcs_reversed', 0, 'output', 'outcar', 'piezo_tensor'])

        doc = {
            'StaticDielectricTensor': {
                'Tensor': None,
                'MicroscopicStaticDielectricTensor': epsilon_static_wolfe,
                'MacroscopicStaticDielectricTensor': epsilon_static
            },
            'BornEffectiveCharges': {
                'Born': None,
                'Born_exLFE': born,
                'Born_inLFE': None
            },
            'PiezoelectricTensor': {
                'Piezo': {
                    'e_Angst': None,
                    'C/m^2': piezo_tensor
                },
                'PiezoelectricTensor_exLFE': None,
                'PiezoelectricTensor_inLFE': None
            }
        }
        return doc

    def getOpticalProperties(self):
        return None

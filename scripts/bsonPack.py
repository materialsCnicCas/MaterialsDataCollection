import os
import time

from services.implements.mpParser import MPParser
from services.implements.mpStructure import MPStructure
from services.implements.vaspParser import VaspParser
from services.implements.vaspStructure import VaspStructure
from utils.tools.calType import CalType

work_time = time.strftime("%Y-%m-%d %H:%M:%S")


def vaspToBson(filePath, collections, user_id, group_id):
    files = os.listdir(filePath)
    # 文件路径
    vaspPath = os.path.join(filePath, 'vasprun.xml')
    if not os.path.exists(vaspPath):
        return ValueError('Lack of must file: vasprun.xml')
    incarPath = os.path.join(filePath, 'INCAR')
    outcarPath = os.path.join(filePath, 'OUTCAR')
    oszicarPath = os.path.join(filePath, 'OSZICAR')
    kpointPath = os.path.join(filePath, 'KPOINTS')

    vaspParser = VaspParser(vaspPath, incarPath, outcarPath, kpointPath, oszicarPath)
    calType = vaspParser.getCalculationType()
    if calType not in collections:
        return None
    vaspParser.setup()

    structure1 = VaspStructure(vaspParser, isinit=True)
    structure1.setup()
    structure1_doc = structure1.toBson()
    structure2 = VaspStructure(vaspParser, isinit=False)
    structure2.setup()
    structure2_doc = structure2.toBson()
    # 根据类型进行打包
    doc = {
        'InputStructure': structure1_doc,
        'OutputStructure': structure2_doc,
        'Parameters': vaspParser.parameters,
        'Software': vaspParser.software,
        'StartTime': vaspParser.startTime,
        'ResourceUsage': vaspParser.resourceUsage,
        'ProcessData': {},
        'Properties': {},
        'Files': [],
        'User': {
            'user': user_id,
            'group': group_id
        },
        'DateCreated': work_time,
        'CalculationType': vaspParser.calculationType
    }
    if vaspParser.calculationType == CalType.GeometryOptimization:
        doc['ProcessData'] = {
            'IonicSteps': vaspParser.ionicSteps,
            'ElectronicSteps': vaspParser.electronicSteps
        }
        doc['Properties'] = {
            'GapFromGeo': vaspParser.gapFromBand  # TODO: gapFromGeo的计算方法和gapFromBand一样？
        }
        doc['Files'] = [vaspPath, incarPath, outcarPath, kpointPath]

    if vaspParser.calculationType == CalType.StaticCalculation:
        doc['ProcessData'] = {
            'ElectronicSteps': vaspParser.electronicSteps
        }
        doc['Properties'] = {
            'ThermodynamicProperties': vaspParser.thermoDynamicProperties,
            'ElectronicProperties': {
                'AtomicCharge': vaspParser.atomicCharge,
                # 'BaderCharge': '',
                # 'GapFromXXX': '',  # TODO 如何提取 Static gap value
            },
            'MagneticProperties': {
                'AtomicMagnetization': vaspParser.atomicMagnetization,
                'LinearMagneticMoment': vaspParser.linearMagneticMoment
            }
        }
        doc['Files'] = [vaspPath, incarPath, outcarPath, kpointPath, oszicarPath]

    if vaspParser.calculationType == CalType.BandStructure:
        doc['ProcessData'] = {
            'ElectronicSteps': vaspParser.electronicSteps
        }
        doc['Properties'] = {
            "ThermodynamicProperties": vaspParser.thermoDynamicProperties,
            "ElectronicProperties": {
                'AtomicCharge': vaspParser.atomicCharge,
                'EigenValues': vaspParser.eigenValues,
                'ProjectedEigenVal_on_IonOrbitals': vaspParser.projectedEigen,
                'GapFromBand': vaspParser.gapFromBand
            },
            'MagneticProperties': {
                'AtomicMagnetization': vaspParser.atomicMagnetization,
                'LinearMagneticMoment': vaspParser.linearMagneticMoment,
            }
        }
        doc['Files'] = [vaspPath, incarPath, outcarPath, kpointPath, oszicarPath]

    if vaspParser.calculationType == CalType.DensityOfStates:
        doc['ProcessData'] = {
            'ElectronicSteps': vaspParser.electronicSteps
        }
        doc['Properties'] = {
            "ThermodynamicProperties": vaspParser.thermoDynamicProperties,
            "ElectronicProperties": {
                'AtomicCharge': vaspParser.atomicCharge,
                'TotalDos': vaspParser.totalDos,
                'PartialDOS': vaspParser.partialDos,
                'GapFromDOS': vaspParser.gapFromDos
                # "DirectedEnergyGap": None,
                # "InDirectedEnergyGap": None
            },
            'MagneticProperties': {
                'AtomicMagnetization': vaspParser.atomicMagnetization,
                'LinearMagneticMoment': vaspParser.linearMagneticMoment,
            }
        }
        doc['Files'] = [vaspPath, incarPath, outcarPath, kpointPath, oszicarPath]

    if vaspParser.calculationType == CalType.ElasticProperties:
        doc['ProcessData'] = {
            'IonicSteps': vaspParser.ionicSteps,
            'ElectronicSteps': vaspParser.electronicSteps
        }
        doc['Properties'] = {
            'ElasticProperties': vaspParser.elasticProperties
        }
        doc['Files'] = [vaspPath, incarPath, outcarPath, kpointPath, oszicarPath]

    if vaspParser.calculationType == CalType.DielectricProperties:
        doc['ProcessData'] = {
            'ElectronicSteps': vaspParser.electronicSteps
        }
        doc['Properties'] = {
            "ThermodynamicProperties": vaspParser.thermoDynamicProperties,
            "ElectronicProperties": {
                'AtomicCharge': vaspParser.atomicCharge
                # 'BaderCharge': None
            },
            'MagneticProperties': {
                'AtomicMagnetization': vaspParser.atomicMagnetization,
                'LinearMagneticMoment': vaspParser.linearMagneticMoment,
            },
            'OpticalProperties': vaspParser.opticalProperties
        }
        doc['Files'] = [vaspPath, incarPath, outcarPath, kpointPath, oszicarPath]
    if vaspParser.calculationType == CalType.MagneticProperties:
        doc['ProcessData'] = {
            'ElectronicSteps': vaspParser.electronicSteps
        }
        doc['Properties'] = {
            "ThermodynamicProperties": vaspParser.thermoDynamicProperties,
            "ElectronicProperties": {
                'AtomicCharge': vaspParser.atomicCharge
            },
            'MagneticProperties': {
                'AtomicMagnetization': vaspParser.atomicMagnetization,
                'LinearMagneticMoment': vaspParser.linearMagneticMoment,
            }
        }
        doc['Files'] = [vaspPath, incarPath, outcarPath, kpointPath, oszicarPath]
    return doc


def mpToBson(filePath, collections, user_id, group_id):
    source_id = filePath.split('\\')[-1].split('.')[0]
    mpParser = MPParser(filePath)
    cal_type = mpParser.getCalculationType()
    if cal_type not in collections:
        return None
    mpParser.setup()
    structure1 = MPStructure(mpParser=mpParser, isinit=True)
    structure1.setup()
    structure1_doc = structure1.toBson()
    structure2 = MPStructure(mpParser=mpParser, isinit=False)
    structure2.setup()
    structure2_doc = structure2.toBson()
    parameters = mpParser.parameters

    if cal_type == CalType.DielectricProperties:
        doc = {
            'InputStructure': structure1_doc,
            'OutputStructure': structure2_doc,
            'Parameters': parameters,
            'Software': mpParser.software,
            'StartTime': mpParser.startTime,
            'ResourceUsage': None,
            'ProcessData': {
                'IonicSteps': mpParser.ionicSteps,
                'ElectronicSteps': mpParser.electronicSteps
            },
            'Properties': {
                'ThermodynamicProperties': mpParser.thermoDynamicProperties,
                'ElectronicProperties': mpParser.electronicProperties,
                'MagneticProperties': mpParser.magneticProperties,
                'OpticalProperties': mpParser.getOpticalProperties(),
                'PiezoelectricProperties': mpParser.getPiezoelectricProperties()
            },
            'Files': [filePath],
            'User': {
                'user': user_id,
                'group': group_id
            },
            'DateCreated': work_time,
            'CalculationType': mpParser.calculationType
        }
    else:
        doc = {
            'InputStructure': structure1_doc,
            'OutputStructure': structure2_doc,
            'Parameters': parameters,
            'Software': mpParser.software,
            'StartTime': mpParser.startTime,
            'ResourceUsage': None,
            'ProcessData': {
                'IonicSteps': mpParser.ionicSteps,
                'ElectronicSteps': mpParser.electronicSteps
            },
            'Properties':{
                'ThermodynamicProperties': mpParser.thermoDynamicProperties,
                'ElectronicProperties': mpParser.electronicProperties,
                'MagneticProperties': mpParser.magneticProperties
            },
            'Files': [filePath],
            'User': {
                'user': user_id,
                'group': group_id
            },
            'DateCreated': work_time,
            'CalculationType': mpParser.calculationType,
        }
    return doc

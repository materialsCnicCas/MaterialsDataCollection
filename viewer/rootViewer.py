import json
import os
import sys
import time

from bson import json_util
from tqdm import tqdm

from dao.toMongo import to_mongo
from scripts.bsonPack import mpToBson, vaspToBson
from services.cifParser import cifToBson
from services.oqmdParser import oqmdMain
from utils.tools.calType import CalType
from dao.user import getUserAndGroup


def materialProjectRun(filePath, log):
    """
    选择不同的任务
    :return:
    """

    s = time.strftime("%Y-%m-%d_%H_%M_%S")
    print('material project files extraction task start ..........')

    database = input('please input database name: (default: MaterialProject)')
    if database == '':
        database = 'MaterialProject'
    collections = []
    collections_all = [CalType.BandStructure, CalType.StaticCalculation, CalType.GeometryOptimization,
                       CalType.DensityOfStates, CalType.DielectricProperties, CalType.ElasticProperties,
                       CalType.MagneticProperties]
    collect = input(
        'please select collections: \n1.BandStructure, \n2.StaticCalculation, \n3.GeometryOptimization, \n4.DensityOfStates, \n5.DielectricProperties, \n6.ElasticProperties, \n7.MagneticProperties.\n input nums（1 2 4）')
    if collect == '':
        collections = collections_all
    else:
        for i in collect.split():
            collections.append(collections_all[int(i)])
    host = input('please input db host: (default localhost)')
    if host == '':
        host = 'localhost'
    port = input('please input db port: (default 27017)')
    if port == '':
        port = '27017'
    port = int(port)
    user = input('please input source user: (default: material project)')
    if user == '':
        user = 'material project'
    group = input(f'please input source user group: (default: {user})')
    if group == '':
        group = user

    file_list = os.listdir(filePath)
    print('Files Dir：', filePath)
    print('Files Number：', len(file_list))
    print('User：', user)
    print('Group：', group)
    print('Database：', database)
    print('Collections：', collections)
    print('host & port：', host, ' ', port)
    # 查询user id 和 group id
    user_id, group_id = getUserAndGroup(host, port, user, group)
    error_files = {'no_take': [], 'error': []}
    if log:
        outFile = open(os.path.join(os.getcwd(), 'log',  database + '_' + user + '_' + group + '_' + s + '.txt'), 'w', encoding='utf8')
    else:
        outFile = sys.stdout

    for file in tqdm(file_list, total=len(file_list)):
        try:
            doc = mpToBson(os.path.join(filePath, file), collections, user_id, group_id)
            if doc is None:
                error_files['no_take'].append(file)
                continue
            # if doc['CalculationType'] in collections:
            to_mongo(doc, database, doc['CalculationType'], host, port)
            # print('to mongo')
        except ValueError as e:
            print(file, ' ', e, file=outFile)
            error_files['error'].append(file)
            continue
        except KeyError as e:
            print(file, ' ', e, file=outFile)
            error_files['error'].append(file)
            continue
        except Exception as e:
            error_files['error'].append(file)
            print(file, ' ', e, file=outFile)
            continue

    with open(os.path.join(os.getcwd(), 'extractResult', 'task_json_' + user + '_' + group + '.json'), 'w') as fp:
        fp.write(json.dumps(error_files, indent=2))
        fp.close()
    print('Extraction Task Over!')


def oqmdRun(log):
    s = time.strftime("%Y-%m-%d_%H_%M_%S")
    print('OQMD Database extraction task start ..........')

    database = input('please input database name: (default: OQMD)')
    if database == '':
        database = 'OQMD'
    collections = []
    collections_all = [CalType.StaticCalculation, CalType.GeometryOptimization]
    collect = input(
        'please select collections: \n1.StaticCalculation, \n2.GeometryOptimization.\n input nums（1 2）')
    if collect == '':
        collections = collections_all
    else:
        for i in collect.split():
            collections.append(collections_all[int(i) - 1])
    host = input('please input db host: (default localhost)')
    if host == '':
        host = 'localhost'
    port = input('please input db port: (default 27017)')
    if port == '':
        port = '27017'
    port = int(port)
    user = input('please input source user: (default: oqmd)')
    if user == '':
        user = 'oqmd'
    group = input('please input source user group: (default: oqmd)')
    if group == '':
        group = 'oqmd'

    print('User：', user)
    print('Group：', group)
    print('Database：', database)
    print('Collections：', collections)
    print('host & port：', host, ' ', port)
    user_id, group_id = getUserAndGroup(host, port, user, group)
    if log:
        outFile = open(os.path.join(os.getcwd(), 'log',  database + '_' + user + '_' + group + '_' + s + '.txt'), 'w', encoding='utf8')
    else:
        outFile = sys.stdout
    error_ids = oqmdMain(database, collections, host, port, user_id, group_id, outFile)
    with open(os.path.join(os.getcwd(),'extractResult', 'task_json_' + user + '_' + group + '.json'), 'w') as fp:
        fp.write(json.dumps({'error_ids': error_ids}, indent=2))
        fp.close()
    print('Extraction Task Over!')


def vaspRun(filePath, log):
    print(os.getcwd())
    s = time.strftime("%Y-%m-%d_%H_%M_%S")
    print('Vasp files extraction task start ..........')

    database = input('please input database name: (default: VaspData)')
    if database == '':
        database = 'VaspData'
    collections = []
    collections_all = [CalType.BandStructure, CalType.StaticCalculation, CalType.GeometryOptimization,
                       CalType.DensityOfStates, CalType.DielectricProperties, CalType.ElasticProperties,
                       CalType.MagneticProperties]
    collect = input(
        'please select collections: \n1.BandStructure, \n2.StaticCalculation, \n3.GeometryOptimization, \n4.DensityOfStates, \n5.DielectricProperties, \n6.ElasticProperties, \n7.MagneticProperties.\n input nums（1 2 4）')
    if collect == '':
        collections = collections_all
    else:
        for i in collect.split():
            collections.append(collections_all[int(i) - 1])
    host = input('please input db host: (default localhost)')
    if host == '':
        host = 'localhost'
    port = input('please input db port: (default 27017)')
    if port == '':
        port = '27017'
    port = int(port)
    user = input('please input source user: (default: vasp)')
    if user == '':
        user = 'vasp'
    group = input(f'please input source user group: (default: {user})')
    if group == '':
        group = user
    file_list = findPaths(filePath)
    print('Files Dir：', filePath)
    print('Files Number：', len(file_list))
    print('User：', user)
    print('Group：', group)
    print('Database：', database)
    print('Collections：', collections)
    print('host & port：', host, ' ', port)
    # 查询user id 和 group id
    user_id, group_id = getUserAndGroup(host, port, user, group)
    error_files = {'no_take': [], 'error': []}
    if log:
        outfile_path = os.path.join(os.getcwd(), 'log',  database + '_' + user + '_' + group + '_' + s + '.txt')
        outFile = open(outfile_path, 'w', encoding='utf8')
    else:
        outFile = sys.stdout
    limit = 1
    for file in tqdm(file_list, total=len(file_list)):
        try:
            doc = vaspToBson(file, collections, user_id, group_id)
            if doc is None:
                error_files['no_take'].append(file)
                continue
            # with open('bson.json', 'w',encoding='utf8') as tfp:
            #     tfp.write(json_util.dumps(doc, indent=2))
            #     tfp.close()
            # if doc['CalculationType'] in collections:
            to_mongo(doc, database, doc['CalculationType'], host, port)
            # print('to mongo')
        except ValueError as e:
            print(file, ' ', e, file=outFile)
            print(file, ' ', e)
            error_files['error'].append(file)
            continue
        except KeyError as e:
            print(file, ' ', e, file=outFile)
            print(file, ' ', e)
            error_files['error'].append(file)
            continue
        except Exception as e:
            error_files['error'].append(file)
            print(file, ' ', e)
            print(file, ' ', e, file=outFile)
            continue
    outFile.close()
    with open(os.path.join(os.getcwd(), 'extractResult', 'task_json_' + user + '_' + group + '_' + s + '.json'), 'w', encoding='utf8') as fp:
        fp.write(json.dumps(error_files, indent=2))
        fp.close()
    print('Extraction Task Over!')


def cifFileParserRun(filePath, log):
    s = time.strftime("%Y-%m-%d_%H_%M_%S")
    print('ICSD CIF Data extraction task start ..........')

    database = input('please input database name: (default: ICSD)')
    if database == '':
        database = 'ICSD'
    collection = 'Structure'
    print('Collection: Structure')

    host = input('please input db host: (default localhost)')
    if host == '':
        host = 'localhost'
    port = input('please input db port: (default 27017)')
    if port == '':
        port = '27017'
    port = int(port)
    user = input('please input source user: (default: icsd)')
    if user == '':
        user = 'icsd'
    group = input('please input source user group: (default: icsd)')
    if group == '':
        group = 'icsd'

    print('User：', user)
    print('Group：', group)
    print('Database：', database)
    print('Collections：', collection)
    print('host & port：', host, ' ', port)
    user_id, group_id = getUserAndGroup(host, port, user, group)
    error_files = {
        'error': []
    }
    if log:
        outFile = open(os.path.join(os.getcwd(), 'log',  database + '_' + user + '_' + group + '_' + s + '.txt'), 'w', encoding='utf8')
    else:
        outFile = sys.stdout
    file_list = os.listdir(filePath)

    for file in tqdm(file_list, total=len(file_list)):
        try:
            doc = cifToBson(os.path.join(filePath, file), user_id, group_id)
            insert_id = to_mongo(doc, database, collection, host, port)
        except Exception as e:
            error_files['error'].append(file)
            print(file, ' ', e, file=outFile)
            continue

    with open(os.path.join(os.getcwd(), 'extractResult', 'task_json_' + user + '_' + group + '.json'), 'w') as fp:
        fp.write(json.dumps(error_files, indent=2))
        fp.close()
    print('Extraction Task Over!')

def findPaths(rootPath):
    total_path = []
    file_list = os.listdir(rootPath)
    if 'vasprun.xml' not in file_list:
        for file in file_list:
            if os.path.isdir(os.path.join(rootPath, file)):
                paths = findPaths(os.path.join(rootPath, file))
                for path in paths:
                    total_path.append(path)
    else:
        total_path.append(rootPath)
        for file in file_list:
            if os.path.isdir(os.path.join(rootPath, file)):
                paths = findPaths(os.path.join(rootPath, file))
                for path in paths:
                    total_path.append(path)
    return total_path




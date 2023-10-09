import json
import os
import sys
import time

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
        outFile = open('./log/' + database + user + '_' + group + '_' + s + '.txt', 'w', encoding='utf8')
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

    with open('./extractResult/task_json_' + user + '_' + group + '.json', 'w') as fp:
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
        outFile = open(r'D:\CZY\DataCollectionSoftware\log/' + database + user + '_' + group + '_' + s + '.txt', 'w', encoding='utf8')
    else:
        outFile = sys.stdout
    error_ids = oqmdMain(database, collections, host, port, user_id, group_id, outFile)
    with open(r'D:\CZY\DataCollectionSoftware\extractResult/task_json_' + user + '_' + group + '.json', 'w') as fp:
        fp.write(json.dumps({'error_ids': error_ids}, indent=2))
        fp.close()
    print('Extraction Task Over!')


def vaspRun(filePath, log):
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
        outFile = open(r'D:\CZY\DataCollectionSoftware\log/' + database + '_' + user + '_' + group + '_' + s + '.txt', 'w', encoding='utf8')
    else:
        outFile = sys.stdout

    for file in tqdm(file_list, total=len(file_list)):
        try:
            doc = vaspToBson(os.path.join(filePath, file), collections, user_id, group_id)
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

    with open(r'D:\CZY\DataCollectionSoftware\extractResult/task_json_' + user + '_' + group + '.json', 'w') as fp:
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
        outFile = open(r'D:\CZY\DataCollectionSoftware\log/' + database + '_' + user + '_' + group + '_' + s + '.txt', 'w', encoding='utf8')
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

    with open(r'D:\CZY\DataCollectionSoftware\extractResult/task_json_' + user + '_' + group + '.json', 'w') as fp:
        fp.write(json.dumps(error_files, indent=2))
        fp.close()
    print('Extraction Task Over!')
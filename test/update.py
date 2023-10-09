#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project    : DataCollectionSoftware
@File       : update.py
@IDE        : PyCharm
@Author     : zychen@cnic.cn
@Date       : 2023/8/22 10:39
@Description: 根据实际需求 实现数据库更新
"""
from pymongo import MongoClient
from tqdm import tqdm

from utils.tools.calType import CalType
from utils.tools.periodicTable import PTable


def update1():
    client = MongoClient('localhost', 27017)
    db = client['MaterialProject']
    collects = [CalType.GeometryOptimization] # ,
    for collect in collects:
        collection = db[collect]
        documents = collection.find({}, no_cursor_timeout=True)
        # new_documents = []
        for docu in tqdm(documents):
            # step1: 重新计算形成能
            _id = docu['_id']
            prop = docu['ThermodynamicProperties']
            tot_energy = prop['TotalEnergy']
            composition = docu['OutputStructure']['Composition']
            tot_num = docu['OutputStructure']['NumberOfSites']
            atom_energy = 0.
            for comp in composition:
                spec = comp['AtomicSymbol']
                num = comp['Amount']
                if spec in PTable().atom_energy:
                    atom_energy += PTable().atom_energy[spec] * num
                else:
                    atom_energy = 0.
                    break
            if atom_energy != 0.:
                formation_energy = (tot_energy - atom_energy) / tot_num
            else:
                formation_energy = 0.
            # step2: 修改字段Name:  ElectroicProperties -》 ElectronicProperties
            prop['FormationEnergy'] = formation_energy

            # doc = {
            #     'InputStructure': docu['InputStructure'],
            #     'OutputStructure': docu['OutputStructure'],
            #     'Parameters': docu['Parameters'],
            #     'Software': docu['Software'],
            #     'StartTime': docu['StartTime'],
            #     'ResourceUsage': docu['ResourceUsage'],
            #     'IonicSteps': docu['IonicSteps'],
            #     'ElectronicSteps': docu['ElectronicSteps'],
            #     'ThermodynamicProperties': docu['ThermodynamicProperties'],
            #     'ElectronicProperties': docu['ElectroicProperties'] if 'ElectroicProperties' in docu else docu['ElectronicProperties'],
            #     'MagneticProperties': docu['MagneticProperties'],
            #     'Files': docu['source_id'],
            #     'CalculationType': docu['CalculationType'],
            #     'User': docu['User'],
            #     'DateCreated': docu['DateCreated']
            # }
            # new_documents.append(doc)
            result = collection.update_one({'_id': _id},
                                           {"$set": {'ThermodynamicProperties.FormationEnergy': formation_energy,
                                                     'Files': docu['source_id'],
                                                     'ElectronicProperties': docu['ElectroicProperties']}})
            if result.matched_count != 1:
                print(_id, ' update error!')
        # x = collection.delete_many({})
        # print(f'已删除{collect}集合中{x.deleted_count}条数据\n正重新插入......')
        # i = collection.insert_many(new_documents)
        # print(f'已插入{collect}集合中{len(i.inserted_ids)}条数据')
        # print('******************************************')


if __name__ == '__main__':
    update1()

#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project    : DataCollectionSoftware 
@File       : modifyFields.py
@IDE        : PyCharm 
@Author     : zychen@cnic.cn
@Date       : 2023/9/4 9:44 
@Description: 
"""
from pymongo import MongoClient
#
# client = MongoClient('localhost', 27017)
# db = client['User']
# coll_group = db['group']
# res = coll_group.find_one({'groupName': 'CNIC'})
# group_id = res['_id']
# coll_user = db['user']
# user_id = coll_user.find_one({'user': 'AI Part'})['_id']
#
# db = client['VaspData']
# coll = db['BandStructure']
# # coll.update_many({}, update={'$unset':{'User':"", "UserGroup": ""}})
# coll.update_many({}, update={"$set": {"User": {"group":group_id, 'user': user_id}}})
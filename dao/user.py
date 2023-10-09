#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project    : DataCollectionSoftware 
@File       : user.py
@IDE        : PyCharm 
@Author     : zychen@cnic.cn
@Date       : 2023/9/4 15:32 
@Description: 
"""
from pymongo import MongoClient


def getUserAndGroup(host, port, user, group):
    """
    Authenticate user and group， return id value，Protect user privacy
    :param host:
    :param port:
    :param user:
    :param group:
    :return: return user_id, group_id
    """
    client = MongoClient(host, port)
    db = client['User']
    user_col = db['user']
    res = user_col.find_one({'user': user})
    if res is None:
        user_id = user_col.insert_one({'user':user}).inserted_id
    else:
        user_id = res['_id']

    group_col = db['group']
    res = group_col.find_one({'group': group})
    if res is None:
        group_id = group_col.insert_one({'group': group}).inserted_id
    else:
        group_id = res['_id']
    return user_id, group_id

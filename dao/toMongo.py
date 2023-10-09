import pymongo
from typing import List


# Mongo version 4.
def to_mongo(doc: dict, database: str, collection: str, host, port):
    """
    将文档数据存入MongoDB中
    :param collection:
    :param database:
    :param doc:
    :param host:
    :param port:
    :return:
       mongo中structure文档的id
    """
    conn = pymongo.MongoClient(host=host, port=port)
    db = conn[database]
    collection = db[collection]
    object_id = collection.insert_one(doc).inserted_id
    conn.close()
    return object_id


def many_to_mongo(docs: List[dict], database: str, collection: str, host, port):
    """
    将多个文档存入mongo
    :param docs:
    :param database:
    :param collection:
    :param host:
    :param port:
    :return:
    """
    conn = pymongo.MongoClient(host=host, port=port)
    db = conn[database]
    collection = db[collection]
    object_ids = collection.insert_many(docs).inserted_ids
    conn.close()
    return object_ids

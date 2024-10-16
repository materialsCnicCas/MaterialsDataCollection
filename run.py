#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project    : CalculationExtract 
@File       : run.py
@IDE        : PyCharm 
@Author     : zychen@cnic.cn
@Date       : 2024/5/30 17:08 
@Description: 
"""
import argparse
import yaml

from user_view import vasp_extract, oqmd_db_extract


def getArgument():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--source', default=None, help='Required. data source: vasp, icsd, oqmd, materialproject, etc.')
    parser.add_argument('--root_dir', default=None, help='Required if file data. calculation files root path')
    parser.add_argument('--log', action='store_true', default=False, help='if start log ')
    parser.add_argument('--user', default=None, help='oqmd database user')
    parser.add_argument('--pwd', default=None, help='oqmd database password')
    parser.add_argument('--database', default=None, help='oqmd database name')
    parser.add_argument('--host', default=None, help='oqmd database host')
    parser.add_argument('--port', default=None, help='oqmd database port')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    with open('config.yaml', 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
        f.close()

    args = getArgument()
    if args.source == 'vasp':
        vasp_extract(args.root_dir, args.log)
    elif args.source == 'oqmd':

        print('Use oqmd database user, password, database, host and port parameters in config.yaml.')
        # oqmd 数据库解析
        host = config['mysql']['host']
        port = config['mysql']['port']
        user = config['mysql']['user']
        password = config['mysql']['password']
        database = config['mysql']['database']

        if args.user is not None:
            user = args.user
        if args.pwd is not None:
            password = args.pwd
        if args.database is not None:
            database = args.database
        if args.host is not None:
            host = args.host
        if args.port is not None:
            port = int(args.port)
        oqmd_db_extract(host, port, user, password, database, args.log)
    else:
        print('source not found, please run "collection.exe --help" ')
    # 其他数据源 补充

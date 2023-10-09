import time

from viewer.rootViewer import materialProjectRun, vaspRun, oqmdRun, cifFileParserRun
import argparse


def getArgument():
    parser = argparse.ArgumentParser(description='')
    # parser.add_argument('--isDB', action='store_true', default=False, help='if data source is database')
    # parser.add_argument('--isFile', action='store_true', default=False, help='if data source is files')
    # parser.add_argument('--db', type=str, default='mysql', help='if data source is database, please input database name!')
    # parser.add_argument('--host', type=str, default='localhost', help='if data source is database, please input database host!')
    # parser.add_argument('--port', type=int, default="3306", help='if data source is database, please input database port!')
    # parser.add_argument('--root_dir', type=str, default=None, help='if data source is files, please input files root dir!')
    # parser.add_argument('--file_type', type=str, default=None, help='if data source is files, please input file type!')
    parser.add_argument('--source', default='vasp', help='data source: vasp, icsd, oqmd, materialproject, etc.')
    parser.add_argument('--root_dir', default='', help='calculation files root path')
    parser.add_argument('--log',action='store_true', default=False, help='if start log ')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = getArgument()
    if args.source == 'vasp':
        vaspRun(filePath=args.root_dir, log=args.log)

    if args.source == 'icsd':
        # icsd 解析
        cifFileParserRun(filePath=args.root_dir, log=args.log)
    if args.source == 'oqmd':
        oqmdRun(args.log)
    if args.source == 'materialproject':
        # material project 解析
        filePath = args.root_dir
        log = args.log
        materialProjectRun(filePath, log)
#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
@Project    : DataCollectionSoftware 
@File       : test.py
@IDE        : PyCharm 
@Author     : zychen@cnic.cn
@Date       : 2023/9/6 8:56 
@Description: 
"""
from viewer.rootViewer import vaspRun, oqmdRun,findPaths


def VaspRunTest():
    filePath = r'D:\extract'
    log = True
    vaspRun(filePath, log)


if __name__ == '__main__':
    # VaspRunTest()
    # oqmdRun(True)
    dirs = findPaths(r'E:\06_III_VCompounds\Intrinsic')
    print(len(dirs))
    print(dirs)
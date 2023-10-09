import os

from tqdm import tqdm

from services.implements.mpParser import MPParser
from utils.tools.calType import CalType
from utils.tools.formatPrint import printDict

if __name__ == '__main__':
    dir_path = r'D:\CZY\TASK_JSON_AUSU_1'
    file_list = os.listdir(dir_path)
    types = [CalType.BandStructure, CalType.DielectricProperties, CalType.DensityOfStates, CalType.MagneticProperties,
             CalType.GeometryOptimization, CalType.ElasticProperties, CalType.StaticCalculation]
    res = {}
    for file in tqdm(file_list):
        try:
            file_path = os.path.join(dir_path, file)
            mp = MPParser(file_path)
            cal_type = mp.getCalculationType()
            if cal_type in types:
                res[cal_type] = file
                types.remove(cal_type)
                printDict(res)
            if len(types) == 0:
                break
        except ValueError:
            continue


############
# D:\CZY\MPDownload\TASK_JSON_1
#
# GeometryOptimization :  mp-1.json
# BandStructure :  mp-1000000.json
# StaticCalculation :  mp-1000351.json
# DielectricProperties  :  mp-1114715.json
# DensityOfStates :  mp-1232449.json

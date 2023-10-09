import json
import os

from tqdm import tqdm

if __name__ == '__main__':
    dir_path = r'D:\CZY\TASK_JSON_AUSU_1'
    file_list = os.listdir(dir_path)

    for file in tqdm(file_list):
        file_path = os.path.join(dir_path, file)
        with open(file_path, 'r', encoding='utf8') as fp:
            data = json.load(fp)
            fp.close()
            if 'data' not in data.keys():
                os.remove(file_path)
                print(file_path)

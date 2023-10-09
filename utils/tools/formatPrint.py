from typing import List, Dict


def printList(items: List):
    for item in items:
        print(item)


def printDict(obj: Dict):
    print('{')
    for key, value in obj.items():
        print(key, ': ', value)
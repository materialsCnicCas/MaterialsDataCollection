from abc import ABC, abstractmethod


class BsonInterface(ABC):
    @abstractmethod
    def toBson(self):
        pass

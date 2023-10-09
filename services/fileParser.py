from abc import ABC, abstractmethod


class FileParser(ABC):

    def __init__(self, filePath):
        if filePath is None:
            raise ValueError('file path 不能为空！')
        self.filePath = filePath
        self.lattice_s = None
        self.lattice_e = None
        self.latticeParameters_s = None
        self.latticeParameters_e = None
        self.composition = None
        self.parameters = None
        self.sites_s = None
        self.sites_e = None
        self.volume_s = None
        self.volume_e = None
        self.numberOfSites = None
        self.startTime = None
        self.software = None
        self.ionicSteps = None
        self.electronicSteps = None
        self.thermoDynamicProperties = None
        self.hasSetup = False

    @abstractmethod
    def setup(self):
        '''
        初始化
        :return:
        '''
        pass

    @abstractmethod
    def getComposition(self):
        """
        获取结构中成分信息
        :return:
        """
        pass

    @abstractmethod
    def getNumberOfSites(self):
        """
        获取原子个数
        :return:
        """
        pass

    @abstractmethod
    def getLatticeParameters(self, isinit: bool = False):
        """
        获取lattice
        :return: lattice, lattice_parameters
        """
        pass

    @abstractmethod
    def getVolume(self, isinit: bool = False):
        """
        获取体积
        :return:
        """
        pass

    @abstractmethod
    def getSites(self, isinit: bool = False):
        """
        获取站点信息
        :param isinit:
        :return:
        """
        pass

    @abstractmethod
    def getParameters(self):
        """
        获取参数信息
        :return:
        """
        pass

    @abstractmethod
    def getSoftware(self):
        """
        获取软件信息
        :return:
           software: dict
        """
        pass

    @abstractmethod
    def getStartTime(self):
        """
        获取开始时间
        :return:
        """
        pass

    @abstractmethod
    def getIonicSteps(self):
        """
        离子步信息提取
        :return:
        """
        pass

    @abstractmethod
    def getElectronicSteps(self):
        """
        电子步骤信息的提取
        :return:
        """
        pass

    @abstractmethod
    def getThermoDynamicProperties(self):
        """
        提取热动力学性质数据
        :return:
        """
        pass


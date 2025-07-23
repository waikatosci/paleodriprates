"""Contains classes related to the input measurements and calibration.

   This module contains the following classes related to input data:
        1.  SampleInfo
        2.  DatingTable
        3.  ProxyDepth
"""
import scipy as sp

class SampleInfo(object):
    """Contains input measurement data from paleo-archives.

        Attributes:
        name        name of the paleoarchive, or the measurement dataset
        datatype    'new'; 'sim' (simulated); 'est' (established)
        archive     'peat'; 'lacustrine'; 'stalagmite'; etc.

        Usage:
        >> import input
        >> X = input.SampleInfo('TestCore123', 'new', 'stalagmite')
        >> Y = [X.name, X.datatype, X.archive]
        >> print Y
        """
    def __init__(self, name, datatype, archive):
        self.name = name
        self.datatype = datatype
        self.archive = archive


class DatingTable(SampleInfo):
    """Contains the dating table, i.e., the age-depth measurements.

        Attributes:
        depth           depth values at which radiometric dating was done
        age             age values corresponding to the entries of 'depth'
        ageerror        errors corresponding to entries of 'age'
        datingmethod    dating method used: 'U/Th', 'C14', 'Cs', etc.

        Usage:
        >> import scipy as sp
        >> import input
        >> rndtest = sp.rand(10,3) #make dummy array
        >> X = input.DatingTable(rndtest[:,0], rndtest[:,1], rndtest[:,2],
                                'U/Th', 'Test123', 'new', 'stalagmite')
        >> Y = [X.depth, X.age, X.ageerror]
        >> print Y

        If X.datingmethod == 'C14', then there are two additional attributes
        which are initiated by the method DatingTable.calibration() --

        cal_reqd        'Yes' (nonexsitent if datingmethod != C14)
        intcal_version  version of the IntCal curve (default == '09')
        postbomb_zone   zone for post-bomb calibration (default == 'NH1')
        """
    def __init__(self, depth, age, ageerror, datingmethod, sampleinfo):
        self.depth = depth
        self.age = age
        self.ageerror = ageerror
        sampleinfo.datingmethod = datingmethod
        self.info = sampleinfo

    def calibration(self, cal_version='13', cal_type="intcal", pb_zone=None):
        if self.info.datingmethod in ['U/Th', 'Cs']:
            self.info.cal_reqd = 'No'
        elif self.info.datingmethod == 'C14':
            self.info.cal_reqd = 'Yes'
            self.info.cal_version = cal_version
            self.info.postbomb_zone = pb_zone
            self.info.cal_type = cal_type


class ProxyDepth(SampleInfo):
    """Contains the proxy data, i.e., the proxy-depth measurements.

        Attributes:
        depth           depth values at which proxy was measured
        proxy           proxy values corresponding to the entries of 'depth'
        proxyerror      errors corresponding to entries of 'proxy'
        sampleinfo      input.SampleInfo instance
        proxyrange      theoretical range of proxy values specified as a
                        list [proxymin, proxymax]. Default = [-inf, inf].

        Usage:
        >> import scipy as sp
        >> import input
        >> rndtest = sp.rand(100,3) #make dummy array
        >> X = input.ProxyDepth(rndtest[:,0], rndtest[:,1], rndtest[:,2],
                                'U/Th', 'Test123', 'new', 'stalagmite')
        >> Y = [X.depth, X.proxy, X.proxyerror]
        >> print Y
        """
    def __init__(self, depth, proxy, proxyerror,sampleinfo,
                 proxyrange=[-sp.inf,sp.inf], proxyname="Proxy123"):
        self.depth = depth
        self.proxy = proxy
        self.proxyerror = proxyerror
        self.info = sampleinfo
        self.proxyrange = proxyrange
        self.name = proxyname

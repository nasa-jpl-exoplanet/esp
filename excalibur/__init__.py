'''Algorithm Engine

The algorithm engine has N goals:

  1. Separate the programatic steps of accomplishing a science goal. Each of the
     Python package is the encapsulation of an algorithm that ends at a science
     goal. Each package then uses extensions of the classes:
        dawgie.Algorithm
        dawgie.StateVector
        dawgie.Task
        dawgie.Value
     to divide the science algorithm into programtic elements. They also enable
     the division of the PipeLine from the Algorithm Engine.
'''

# -- IMPORTS -- ------------------------------------------------------
# import builtins
import collections
import dawgie

# import numpy
# import scipy.stats

import os

# ------------- ------------------------------------------------------
# GMR: CAN WE SET UP THIS MESS ONCE AND FORGET ABOUT IT
# Ines Mertz : I fixed it
context = {
    'data_cal': os.environ.get('DATA_CALIBR', '/proj/sdp/data/cal'),
    'data_dir': os.environ.get('DATA_BASEDIR', '/proj/sdp/data'),
    'data_sci': os.environ.get('DATA_SCIENC', '/proj/sdp/data/sci'),
    'ldtk_root': os.environ.get('LDTK_ROOT', '/proj/sdp/data/ldtk'),
    'target_list': os.environ.get(
        'TARGET_LIST', '/proj/sdp/data/WFC3_target_list.xlsx'
    ),
}
os.environ['LDTK_ROOT'] = context['ldtk_root']
__version__ = '${ESP_GIT_REV}'

Identity = collections.namedtuple('identity', ['serial'])


class ValuesList(dawgie.Value, list):
    '''ValuesList ds'''

    def __init__(self, *args, **kwds):
        '''__init__ ds'''
        list.__init__(self, *args, **kwds)
        self._version_ = dawgie.VERSION(1, 1, 0)
        return

    def features(self):
        '''features ds'''
        return []

    pass


class ValuesDict(dawgie.Value, dict):
    '''ValuesDict ds'''

    def __init__(self, *args, **kwds):
        '''__init__ ds'''
        dict.__init__(self, *args, **kwds)
        self._version_ = dawgie.VERSION(1, 1, 0)
        return

    def features(self):
        '''features ds'''
        return []

    pass


class ValueScalar(dawgie.Value):
    '''ValueScalar'''

    def __init__(self, content=None):
        '''__init__ ds'''
        dawgie.Value.__init__(self)
        self.__content = content
        self._version_ = dawgie.VERSION(1, 1, 0)
        return

    def features(self):
        '''features ds'''
        return []

    def new(self, value):
        '''method to keep from explicitly needing dawgie'''
        return ValueScalar(value if value is not None else self.__content)

    def value(self):
        '''value ds'''
        return self.__content

    pass

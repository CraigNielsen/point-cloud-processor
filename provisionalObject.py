'''
Created on 25 Feb 2015

@author: FRGCRA003
'''
from sympy import *


class provObject:
    '''
    Container for provisional values
    '''

    def __init__(self, symbol_, value_):
        
        self.symbol=symbol_
        self.value=value_
        
'''
Created on 25 Feb 2015

@author: FRGCRA003
'''
from sympy import *

class obsObject():
    '''
    container for observations
    '''


    def __init__(self, symbol_, value_=0.):
        
        self.symbol=symbol_
        self.value = value_
        
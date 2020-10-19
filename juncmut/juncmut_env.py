#!/usr/bin/env python3

"""
Created on Wed Jul 31 2019

@author: naokoIida
"""

def juncmut_env():
    import os
    
    if not os.path.exists('./juncmut'):
        os.makedirs('./juncmut')
        
if __name__== "__main__":
    
    juncmut_env()

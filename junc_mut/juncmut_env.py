#!/usr/bin/env python3

"""
Created on Wed Jul 31 2019

@author: naokoIida
"""

def juncmut_env(folder):
    import os
    
    if not os.path.exists(folder + './alterativeSJ_fil_annot'):
        os.mkdir(folder + './alterativeSJ_fil_annot')
    if not os.path.exists(folder + './alterativeSJ_fil_annot'):
        os.mkdir(folder + './alterativeSJ_fil_annot')
    if not os.path.exists(folder + './alterativeSJ_assadjfreq'):
        os.mkdir(folder + './alterativeSJ_assadjfreq')
    if not os.path.exists(folder + './alterativeSJ_assadjfreq'):
        os.mkdir(folder + './alterativeSJ_assadjfreq')
    if not os.path.exists(folder + './alterativeSJ_mutprediction'):
        os.mkdir(folder + './alterativeSJ_mutprediction')
    if not os.path.exists(folder + './alterativeSJ_mutprediction'):
        os.mkdir(folder + './alterativeSJ_mutprediction')
        
if __name__== "__main__":
    import sys
    
    folder = sys.argv[1]
    print(folder) 
    juncmut_env(folder)

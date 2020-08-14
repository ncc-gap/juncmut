#!/usr/bin/env python3

"""
Created on Wed Jul 31 2019

@author: naokoIida
"""

def juncmut_env(output_dir):
    import os
    
    if not os.path.exists(output_dir + '/alterativeSJ_fil_annot'):
        os.makedirs(output_dir + '/alterativeSJ_fil_annot')
    if not os.path.exists(output_dir + '/alterativeSJ_fil_annot'):
        os.makedirs(output_dir + '/alterativeSJ_fil_annot')
    if not os.path.exists(output_dir + '/alterativeSJ_assadjfreq'):
        os.makedirs(output_dir + '/alterativeSJ_assadjfreq')
    if not os.path.exists(output_dir + '/alterativeSJ_assadjfreq'):
        os.makedirs(output_dir + '/alterativeSJ_assadjfreq')
    if not os.path.exists(output_dir + '/alterativeSJ_mutprediction'):
        os.makedirs(output_dir + '/alterativeSJ_mutprediction')
    if not os.path.exists(output_dir + '/alterativeSJ_mutprediction'):
        os.makedirs(output_dir + '/alterativeSJ_mutprediction')
        
if __name__== "__main__":
    import sys
    
    output_dir = sys.argv[1]
    print(output_dir) 
    juncmut_env(output_dir)

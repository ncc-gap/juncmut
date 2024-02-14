#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def juncmut_env():
    import os
    if not os.path.exists('./juncmut'):
        os.makedirs('./juncmut')
        
if __name__== "__main__":
    juncmut_env()

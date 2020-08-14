#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 13:01:49 2019

@author: genome

"""


def juncmut_assadj(input_file, output_file):

    # import os
    sample = "dummy"

    # sample = os.path.basename(input_SJ).replace(".SJ.out.tab", '')
    # file2 = "%s/alterativeSJ_fil_annot/%s.SJ.fil.annot.txt" %(output_dir, sample)
    # file3= "%s/alterativeSJ_assadjfreq/%s.SJ.fil.annot.assadj.txt" %(output_dir, sample)
        
    with open(input_file, 'r') as in1:
            with open(output_file, 'w') as out1:
    
                for line in in1:
                    F = line.rstrip('\n').split('\t')
                    r=F[6]
                    
                    if "Alternative" in F[9] or "alternative" in F[9]:
                        if "e" in F[13]: #strand=+ 3'SS
                            a =str(F[14])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(sample)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(sample)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                        elif "s" in F[13]: #strand=- 5'SS
                            a =str(F[14])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(sample)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(sample)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                        elif "s" in F[17]: #strand=+ 5'SS
                            a =str(F[18])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(sample)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(sample)+'\t'+str(F[9])+'\t+\t'+str(r)+'\n')
                        elif "e" in F[17]: #strand=- 3'SS
                            a =str(F[18])
                            b = a.split(';',2)
                            if b[0] == '*':
                                b[0] = '0'
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(sample)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')
                            else:
                                c=F[0]
                                si=int(F[1])-1*int(b[0])
                                ei=int(F[2])-1*int(b[0])
                                out1.write(str(c)+'\t'+str(si)+'\t'+str(ei)+'\t'+str(F[1])+'\t'+str(F[2])+'\t'+str(sample)+'\t'+str(F[9])+'\t-\t'+str(r)+'\n')

if __name__== "__main__":
    import argparse

    parser = argparse.ArgumentParser() #make a parser

    parser.add_argument("input", metavar = "prefix. sample.SJ.fil.annot.txt", default = None, type = str,
                        help = "input file") 
    
    parser.add_argument("folder", metavar = "group", default = "my_samples", type = str,
                        help = "folder name") 
    
    args = parser.parse_args()

    pr = args.input
    folder = args.folder

    print(pr) 
    juncmut_assadj(pr, folder)

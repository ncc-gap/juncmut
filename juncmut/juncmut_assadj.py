#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 13:01:49 2019

@author: genome

"""

def juncmut_assadj(input_file, output_file):

    sample = input_file.split('.')[0]
        
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

    parser.add_argument("input_file", metavar = "input_file", default = None, type = str,
                        help = "input file") 
    
    parser.add_argument("output_file", metavar = "output_file", default = "my_sample", type = str,
                        help = "output_file") 
    
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    print(input_file) 
    juncmut_assadj(input_file, output_file)

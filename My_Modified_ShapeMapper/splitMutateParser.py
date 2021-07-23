#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import getopt
import os
import random

Usage = """
splitMutateParser - Split parsed mutate files
===========================================================================
\x1b[1mUSAGE:\x1b[0m 
  %s [--max_mut 10] annonate_mutate_file out_dir file_name

  max_mut               -- Maximun number of mutation events in a read (default: 10)
  annonate_mutate_file  -- A file produced by shapemapper_mutation_parser and the reference is append to column 1

\x1b[1mVERSION:\x1b[0m
    %s

\x1b[1mAUTHOR:\x1b[0m
    Li Pan

""" % (sys.argv[0], "2019-03-14")


def init():
    params = { 'inputFile': None, 'outFolder': None, 'fileName': None, 'maxMut':10 }
    if len(sys.argv) < 4:
        sys.stdout.writelines(Usage)
        sys.exit(-1)
    
    params['inputFile'] = os.path.abspath(sys.argv[-3])
    params['outFolder'] = os.path.abspath(sys.argv[-2])
    params['fileName'] = sys.argv[-1]
    
    opts, args = getopt.getopt(sys.argv[1:-3], 'h', ['max_mut='])
    
    for op, value in opts:
        if op == '-h':
            sys.stdout.writelines(Usage)
            exit(-1)
        
        elif op == '--max_mut':
            params['maxMut'] = int(value)
        
        else:
            sys.stdout.writelines("parameter Error: unrecognized parameter: "+op)
            sys.stdout.writelines(Usage)
            sys.exit(-1)

    return params


def write_folder(out_folder, tid, file_name, data_list):
    if len(data_list) > 100:
        cur_folder = out_folder + tid
        cur_file = cur_folder + "/" + file_name
        try:
            os.mkdir(cur_folder)
        except OSError:
            pass
        OUT = open(cur_folder + "/" + file_name, "w")
        for line in data_list:
            OUT.writelines(line+"\n")
        OUT.close()


def split_targets(inFn, file_name, out_folder, max_mut_num=10):
    outFolder = out_folder.rstrip('/') + "/"
    cur_tid = ""
    cur_folder = ""
    cur_data_list = []
    for line in open(inFn):
        data = line.rstrip("\n").split("\t")
        if data[5] == 'LOW_MAPQ': continue
        if len(data)>=11 and data[10]:
            mut_num = len(data[10].split(' ')) / 5
            if mut_num > max_mut_num:
                continue
        if data[0] != cur_tid:
            if cur_tid:
                write_folder(outFolder, cur_tid, file_name, cur_data_list)
            
            cur_data_list = []
            cur_tid = data[0]
        
        cur_data_list.append( "\t".join(data[1:]) )
    
    if cur_tid:
        write_folder(outFolder, cur_tid, file_name, cur_data_list)


if __name__ == "__main__":
    params = init()
    split_targets(params['inputFile'], params['fileName'], params['outFolder'], params['maxMut'])


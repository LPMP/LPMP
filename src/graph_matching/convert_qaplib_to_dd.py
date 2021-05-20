#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import os
import sys

def read_qaplib_instance(f):
    first = True
    counter = 0
    Aset = False
    Bset = False
    offset = 0
    for line in f:
        line = line.strip()
        line = line.split()
        if len(line) > 0:
            if first:
                #assert len(line) == 1, 'file: {}, length: {}, content: {}'.format(n, len(line), line)
                instance_size = int(line[0])
                first = False
                A = np.zeros((instance_size, instance_size), dtype='int')
                B = np.zeros((instance_size, instance_size), dtype='int')
            elif not Aset:
                assert len(line) <= instance_size - offset, 'file: {}, length: {}, content: {}'.format(n, len(line), line)
                A[counter,offset:offset+len(line)] = [int(number) for number in line]
                offset += len(line)
                if offset == instance_size:
                    offset = 0
                    counter += 1
                if counter == instance_size:
                    counter = 0
                    Aset = True
            elif not Bset:
                assert len(line) <= instance_size - offset, 'file: {}, length: {}, content: {}'.format(n, len(line), line)
                B[counter,offset:offset+len(line)] = [int(number) for number in line]
                offset += len(line)
                if offset == instance_size:
                    offset = 0
                    counter += 1
                if counter == instance_size:
                    Bset = True
            else:
                assert False, 'something went wrong - there seem to be additional lines in the data {}, counter {}, line {}'.format(n, counter, line)
    return instance_size, A, B

def generate_dd_content(instance_size, A, B):
    assignments = []
    edges = []
    for i in range(instance_size):
        for k in range(instance_size):
            assignments.append((i,k,A[i,i]*B[k,k]))
    for i in range(instance_size-1):
        for j in range(i,instance_size):
            for k in range(instance_size):
                for p in range(instance_size):
                    if p != k:
                        cost = A[i,j]*B[k,p] + A[j,i]*B[p,k]
                        if cost != 0:
                            edges.append((i*instance_size+k, j*instance_size+p, cost))
    return assignments, edges


if __name__ == '__main__':
    filename = sys.argv[1]
    dd_filename = sys.argv[2]

    print (filename)
    print (dd_filename)
    if os.path.isfile(dd_filename):
        print('File ' + dd_filename + ' already exists')
        exit(0)
    with open(filename, 'rt') as file:
        instance_size, A, B = read_qaplib_instance(file)
    print('  instance size: {}'.format(instance_size))
    #with open(qapdd + 'README', 'at') as file:
     #   file.write('{} -> qap{}\n'.format(filename, file_counter))
    #incentive = instance_size*np.max(A)*np.max(B)
    #print('  incentive: {}'.format(incentive))
    assignments, edges = generate_dd_content(instance_size, A, B)
    with open(dd_filename, 'wt') as file:
        file.write('c This file was generated from the QAPLIB instance {}\n'.format(filename))
        #file.write('c To each assignment a negative incentive of {} was added\n'.format(-incentive))
        #file.write('c to force the minimal assignment to be the same even when\n')
        #file.write('c allowing for incomplete matchings.\n')
        file.write('p {} {} {} {}\n'.format(instance_size, instance_size, len(assignments), len(edges)))
        print('..writing assignments to file')
        for i, a in enumerate(assignments):
            #file.write('a {} {} {} {}\n'.format(i, a[0], a[1], a[2]-incentive))
            file.write('a {} {} {} {}\n'.format(i, a[0], a[1], a[2]))
        print('..writing edges to file')
        for i, e in enumerate(edges):
            file.write('e {} {} {}\n'.format(e[0], e[1], e[2]))
        print('..finished')


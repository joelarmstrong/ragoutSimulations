#!/usr/bin/env python
"""Slightly modified version of simCtrl_postSimFastaExtractor.py,
originally part of evolverSimControl and developed by Dent Earl.

eval_evolverpairwiseMAFextractor.py
dent earl, dearl (a) soe ucsc edu
16 nov 2009
A script that will take the path to a simulation
out directory and makes repeated calls to cvt
(in fabulous PARALLEL-vision!) in order
to extract FASTA files from either all cycles, or
only the leafs.

"""
##################################################
# Copyright (C) 2009-2011 by
# Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
# ... and other members of the Reconstruction Team of David Haussler's
# lab (BME Dept. UCSC)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##################################################
import glob
import os
import re
import subprocess
import sys
from argparse import ArgumentParser
from sonLib.bioio import newickTreeParser, system
import xml.etree.ElementTree as ET

def initOptions(parser):
    parser.add_argument('simDir',
                        help = 'The simulation directory.')
    parser.add_argument('destDir',
                        help = 'The directory to extract the fastas to.')
    parser.add_argument('--allCycles', action = 'store_true', dest = 'allCycles',
                        default = False, help = 'Extract fastas from all cycles, not just leafs. '
                        'default=%default')
    parser.add_argument('--includeAncestors', action = 'store_true',
                        default = False, help = 'Extract fastas from ancestors as well'
                        'default=%default')

def checkOptions(options, parser):
    if options.simDir is None:
        parser.error('specify --simDir.\n')
    if not os.path.isdir(options.simDir):
        parser.error('simulation directory "%s" not a directory!\n' % options.simDir)
    options.simDir = os.path.abspath(options.simDir)
    if not os.path.exists(os.path.join(options.simDir, 'simulationInfo.xml')):
        parser.error('unable to find %s.\n' % os.path.join(options.simDir, 'simulationInfo.xml'))
    infoTree = ET.parse(os.path.join(options.simDir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick=treeObj.text
    treeObj = infoTree.find('rootDir')
    options.rootDir=treeObj.text

def directoriesOnly(aList):
    """directoriesOnly() takes a list of items from a directory
    and purges anything from the list that is not a directory.
    """
    bList = []
    for i in aList:
        if os.path.isdir(i):
            bList.append(i)
    return bList

def nameTree(nt, reportDistance = True):
    """nameTree(nt) takes a newick tree and returns a str that can be used
    to name the cycle-step that the tree represents. Distance included in
    the name by default.
    """
    from sonLib.bioio import printBinaryTree
    if nt is None:
        print 'nt was none'
        return ''
    if nt.iD is not None:
        if nt.distance == 0.0 or not reportDistance:
            name = nt.iD
        else:
            name = nt.iD + str(nt.distance)
    else:
        name = printBinaryTree(nt, True)
    name = sanitizeTreeName(name)
    return name

def sanitizeTreeName(name):
    """sanitizeTreeName(name) takes all the nasty characters out of a newickTree and
    returns a str that is more amenable to being a file (or directory) name.
    """
    name = name.replace(' ','')
    name = name.replace(',','')
    name = name.replace(':','-')
    name = name.replace('.','_')
    name = name.replace(';','')
    name = name.replace('\'','')
    name = name.replace('"','')
    name = name.replace('(','_L_')
    name = name.replace(')','_R_')
    name = name.rstrip('0')
    name = name.rstrip('-0_')
    return name

def getSetToExtract(nt, onlyLeaves):
    """Given a newick tree object, it returns a dict of
    leaf objects. Operates recursively.
    """
    def recurse(nt, set):
        if nt is None:
            return None
        nt.distance = 0
        if nt.right is None and nt.left is None:
            set.add(nt.iD)
        else:
            if not onlyLeaves:
                set.add(nameTree(nt))
            recurse(nt.right, set)
            recurse(nt.left , set)
    setToExtract = set()
    recurse(nt, setToExtract)
    return setToExtract

def main():
    usage = ('usage:  --simDir path/to/dir [options]\n\n'
             ' takes in a simulation directory and then extracts\n'
             'the sequence of each leaf node in fasta format and stores them\n'
             'in the respective step\'s directory.')
    parser = ArgumentParser(usage = usage)
    initOptions(parser)
    options = parser.parse_args()
    checkOptions(options, parser)
    
    cycles = glob.glob(os.path.join(options.simDir, '*'))
    cycles = directoriesOnly(cycles)
    nt = newickTreeParser(options.inputNewick, 0.0)
    setToExtract = getSetToExtract(nt, not options.includeAncestors)
    os.makedirs(options.destDir)
    for d in cycles:
        print d, setToExtract
        if not options.allCycles and not os.path.basename(d) in setToExtract:
            continue

        nameA     = os.path.basename(d)
        nameA     = nameA.replace('[','')
        nameA     = nameA.replace(']','')
        cleanName = nameA.replace('\'','')
        
        cmd = ['evolver_cvt']
        cmd.append('-fromrev')
        cmd.append(os.path.join(d,'seq.rev'))
        cmd.append('-tofasta')
        cmd.append(os.path.join(d, 'seq.fa.tmp'))
        system(" ".join(cmd))
        
        cmd = ['mv']
        cmd.append(os.path.join(d, 'seq.fa.tmp'))
        cmd.append(os.path.join(d, 'seq.fa'))
        system(" ".join(cmd))
        
        cmd = ['sed']
        cmd.append('-i')
        cmd.append(r"'s/^>/>%s./;'" % cleanName)
        cmd.append(os.path.join(d, 'seq.fa'))
        system(" ".join(cmd))
        
        cmd = ['mv']
        cmd.append(os.path.join(d, 'seq.fa'))
        cmd.append(os.path.join(options.destDir, '%s.name.fa' % nameA))
        system(" ".join(cmd))
    
if __name__ == "__main__":
    main()

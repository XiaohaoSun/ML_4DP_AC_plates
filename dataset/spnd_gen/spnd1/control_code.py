import os
import sys
import random
import subprocess
import math
import csv
import numpy as np
from scipy import interpolate

# from scoop import futures

# from design_parameters import *

from matplotlib import rcParams
import matplotlib.pyplot as plt

def evaluate_design(ind_num):
    
    dir_files = os.listdir('.')
    split_files = []

    extensions = ['rec', 'rpy', 'inp', 'com','dat', 'msg', 'prt',
                  'sim', 'ipm', 'mdl', 'stt', 'sta', '1', '023']
    
    for fil in dir_files:
        split_files.append(fil.split('.')[0])
        
        try:
            e = fil.split(".")[1]
        except IndexError:
            e = ""
        
        if e in extensions:
            try:
                os.remove(fil)
            except OSError:
                pass
    
    GENERATE_SCRIPT = 'Exp2D_nocut_BC3.py'
    genome_name = 'temp_spnd1_individual'
    
    genome_name = genome_name + str(ind_num)
    if genome_name not in split_files:
        print('major problem')
        return

    abaqus_str = 'abaqus cae noGUI=' + GENERATE_SCRIPT
    job_name = 'job_' + genome_name
    abaqus_str = abaqus_str + ' -- ' + genome_name + '.txt ' + job_name

    fnull = open(os.devnull, 'wb')
    subprocess.call(abaqus_str,
                    shell=True,
                    stderr=fnull,
                    stdout=fnull)
    
                
    return ind_num,



def clear_directory():
    files = os.listdir('.')
    os.system('del *.csv')
    extensions = ['png', 'rec', 'odb', 'sta', 'msg', 'rpy',
                  'dat', 'log', 'inp', 'com', 'prt',
                  'sim', 'ipm', 'mdl', 'stt', '1', '023',
                  '.csv']

    for f in files:
        try:
            e = f.split(".")[1]
        except IndexError:
            e = ""

        if e in extensions:
            try:
                os.remove(f)
            except OSError:
                pass


def make_user_subroutines():
    os.system('abaqus make library=utrs.for >NUL')


if __name__ == '__main__':

    # make_user_subroutines()
    # clear_directory()
    
    # Parameters given below, not in a separate file
    # Voxel info are not defined here, but in the txt file and will be read by Exp2D_nocut_BC3.py.
    
    path_current = os.getcwd()
    path_CSV = os.path.join(os.path.dirname(path_current),'CSV_spnd1')
    if not os.path.exists(path_CSV):
        os.makedirs(path_CSV)
    
    # group, ind1, ind2
    args = sys.argv
    for arg in args:
        if 'group' in arg:
            index = args.index(arg)
    my_args = args[index:]
    group = my_args[1]
    print('This is group',group)
    
    groups = [1,2,3,4,5,6] # contain (okay to have more than) actual group numbers
    
    if int(group) in groups:
        ind1, ind2 = (int(group)-1)*1000+1, int(group)*1000
        # ind1, ind2 = 8, 9
    else:
        print('Major problem: group undefined!')
    
    for ind in range(ind1,ind2+1):
        evaluate_design(ind)
    


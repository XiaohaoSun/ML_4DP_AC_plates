import numpy as np
import random
import os
import sys
import csv

from abaqus import *
from abaqusConstants import *
import sketch
import part
import assembly
import material
import step
import interaction
import regionToolset
from mesh import *
import time
import visualization
import odbAccess
import math
import datetime

def check_job(job_name):
    try:
        with open(job_name+'.log','r') as f_log:
            f_log.seek(-10,2)
            status_log = f_log.read(9)
        if status_log == 'COMPLETED':
            aborted = False
        else:
            aborted = True
    except IOError:
        # triggered if job_name.log does not exist. 
        # There indeed exists job_name.*, see the main loop at the bottom. 
        print('*No results for', ind_num) 
        aborted = True
    
    return aborted

def csv_output_new(aborted, job_name):
    out = 'outs-' + job_name.split('.')[0] + '.csv'

    if aborted:
        print('*Simulation not completed for', ind_num)
        # fout.write('*Simulation Failed\n')
        return

    try:
        odb = session.openOdb(name=job_name + '.odb')
    except RuntimeError:
        pass

    if len(odb.steps['Actuation'].frames) == 0:
        # fout.write('*Simulation Failed\n')
        return

    nodes = odb.rootAssembly.nodeSets['SET-FACEMID']
    node_labels, x_coords, y_coords, z_coords = [], [], [], []
    for node in nodes.nodes[0]:
        node_labels.append(node.label)
        x_coords.append(node.coordinates[0])
        y_coords.append(node.coordinates[1])
        z_coords.append(node.coordinates[2])
    
    x_disp, y_disp, z_disp = [], [], []
    for i_f, frame in enumerate(odb.steps['Actuation'].frames):
        if i_f == 0:
            continue
        
        u = frame.fieldOutputs['U'].getSubset(region=nodes)
        disp_node_labels, x_disp_temp, y_disp_temp, z_disp_temp = [], [], [], []
        for v in u.values:
            disp_node_labels.append(v.nodeLabel)
            x_disp_temp.append(v.data[0])
            y_disp_temp.append(v.data[1])
            z_disp_temp.append(v.data[2])

        for i in range(len(node_labels)):
            if node_labels[i] != disp_node_labels[i]:
                print('major problem')
                return

        x_disp.append(x_disp_temp)
        y_disp.append(y_disp_temp)
        z_disp.append(z_disp_temp)

    # writing data
    fieldnames = ['node_label','x_coord', 'y_coord', 'z_coord']
    for i_f in range(1,len(odb.steps['Actuation'].frames)):
        fieldnames.append('f'+str(i_f)+'_x_disp')
        fieldnames.append('f'+str(i_f)+'_y_disp')
        fieldnames.append('f'+str(i_f)+'_z_disp')
        

    with open(out, 'wb') as csv_file:

        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()

        for i in range(len(node_labels)):
            node_dict = {'node_label': node_labels[i],
                        'x_coord': x_coords[i],
                        'y_coord': y_coords[i],
                        'z_coord': z_coords[i]}
            for i_f in range(1,len(odb.steps['Actuation'].frames)):
                node_dict['f'+str(i_f)+'_x_disp'] = x_disp[i_f-1][i]
                node_dict['f'+str(i_f)+'_y_disp'] = y_disp[i_f-1][i]
                node_dict['f'+str(i_f)+'_z_disp'] = z_disp[i_f-1][i]

            writer.writerow(node_dict)
    
    # close odb
    try:
        odb.close()
    except RunTimeError:
        pass
    


job_names = ['Nocut2_subR1']

for job_name in job_names:
    # aborted = check_job(job_name)
    aborted = False
    csv_output_new(aborted,job_name)


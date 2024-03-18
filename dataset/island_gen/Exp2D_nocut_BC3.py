import numpy as np
import random
import os
import sys
import csv

# from design_parameters import *

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


def read_genome(f):
    with open(f, 'r') as gene:
        gene.readline()
        voxel_sizes = gene.readline().split(', ')
        nx = int(voxel_sizes[0])
        ny = int(voxel_sizes[1])
        nz = int(voxel_sizes[2])
        locations = np.zeros((nx, ny, nz), dtype=int)
        for iz in range(nz):
            for iy in range(ny):
                raw_locations = gene.readline()
                locations[:,iy,iz] = (raw_locations.split())
        return nx,ny,nz,locations

class Active2DStructure:

    def __init__(self,
                 length_x, length_y, height, nx, ny, nz,
                 voxel_size_x, voxel_size_y, voxel_size_z,
                 locations, num_cpus,
                 alpha, temp_low, temp_high,
                 seed_xy, seed_z,
                 stepTime, iniInc, maxInc):
        """
        Class constructor
        """
        initial = datetime.datetime.now()
        print('Started at ' + str(initial))
        
        self.num_cpus = num_cpus

        self.length_x = length_x
        self.length_y = length_y
        self.height = height
        
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.voxel_size_x = voxel_size_x
        self.voxel_size_y = voxel_size_y
        self.voxel_size_z = voxel_size_z
        
        self.locations = locations # voxel info
        
        self.alpha = alpha
        self.temp_low = temp_low
        self.temp_high = temp_high
        
        self.seed_xy = seed_xy
        self.seed_z = seed_z
        
        self.stepTime = stepTime
        
        self.iniInc = iniInc
        self.maxInc = maxInc
        
        self.model_name = 'Candidate' + str(random.randint(0, 500000))
        self.m = mdb.Model(name=self.model_name)
        #
        # start calling all functions
        #
        self.generate_geometry()
        self.generate_assembly()
        self.create_and_assign_materials_and_sections()
        # self.total_number_of_chunks = self.create_assembly_instance_sets()
        self.create_steps()
        # self.initialize_interactions()
        self.set_initial_conditions_and_step_temperatures()
        self.set_boundary_conditions()

        self.mesh_structure()
        
        self.m.rootAssembly.regenerate()
        
        print('Started at ' + str(initial))
        print('Completed at ' + str(datetime.datetime.now()) )
        
    def generate_geometry(self):
        
        sketch_prof = self.m.ConstrainedSketch(name='Profile-Sketch',
                                               sheetSize=20.0)
        sketch_prof.rectangle(point1 = (0.0, 0.0),
                              point2 = (self.voxel_size_x, self.voxel_size_y) )
    
        part_unit = self.m.Part(name='Part-Voxel',
                                dimensionality=THREE_D,
                                type=DEFORMABLE_BODY)
        part_unit.BaseSolidExtrude(sketch=sketch_prof,
                                   depth=self.voxel_size_z)

    def generate_assembly(self):
        
        a = self.m.rootAssembly
        voxels_Inst = []
        for iz in range(self.nz):
            for ix in range(0, self.nx):
                for iy in range(0, self.ny):
                    a.Instance(name='Assem-voxel-zi_' + str(iz) + '-xi_' + str(ix) + '-yi_' + str(iy),
                            part=self.m.parts['Part-Voxel'],
                            dependent=ON)
                    a.translate(instanceList = ['Assem-voxel-zi_' + str(iz) + '-xi_' + str(ix) + '-yi_' + str(iy)],
                                vector = (ix*self.voxel_size_x, iy*self.voxel_size_y, iz*self.voxel_size_z) )
                    voxels_Inst.append( a.instances['Assem-voxel-zi_' + str(iz) + '-xi_' + str(ix) + '-yi_' + str(iy)] )
        
        
        newName = 'Part-All'
        a.InstanceFromBooleanMerge(instances = voxels_Inst,
                                   keepIntersections = ON,
                                   domain=GEOMETRY,
                                   mergeNodes= NONE,
                                   name = newName,
                                   originalInstances=SUPPRESS)
        

    def create_and_assign_materials_and_sections(self):
        
        materActive = self.m.Material(name='Active_Material')
        materActive.Hyperelastic(type=NEO_HOOKE,testData=OFF, table=((0.1, 0.002,),))
        materActive.Expansion(type=ISOTROPIC, table=((self.alpha,),))
        materActive.Density(table=((1e-6,),))
        
        materPassive = self.m.Material(name='Passive_Material')
        materPassive.Hyperelastic(type=NEO_HOOKE,testData=OFF, table=((0.1, 0.002,),))
        materPassive.Expansion(type=ISOTROPIC, table=((0.0,),))
        materPassive.Density(table=((1e-6,),))
        
        self.m.HomogeneousSolidSection(material='Active_Material',
                                        name='Active_Section')
                                        
        self.m.HomogeneousSolidSection(material='Passive_Material',
                                        name='Passive_Section')
                                        
        p = self.m.parts['Part-All']
        p.Set(name='All', cells=p.cells)
        
        
        # ======= MATERIAL ASSIGNMENT SECTION ======= #
        for iz in range(self.nz):
            for ix in range(0, self.nx):
                for iy in range(0, self.ny):
                    pSetName = 'Section-voxel-zi_' + str(iz) + '-xi_' + str(ix) + '-yi_' + str(iy)
                    x1 = ix * self.voxel_size_x
                    y1 = iy * self.voxel_size_y
                    z1 = iz * self.voxel_size_z
                    x2 = (ix+1) * self.voxel_size_x
                    y2 = (iy+1) * self.voxel_size_y
                    z2 = (iz+1) * self.voxel_size_z
                    
                    c_temp = p.cells.getByBoundingBox(x1, y1, z1, x2, y2, z2)
                    p.Set(name = pSetName,
                          cells = c_temp)
        
        active_cell = ()
        passive_cell = ()
        for iz in range(self.nz):
            for ix in range(0, self.nx):
                for iy in range(0, self.ny):
                    pSetName = 'Section-voxel-zi_' + str(iz) + '-xi_' + str(ix) + '-yi_' + str(iy)
                    s = p.sets[pSetName]
                    if self.locations[ix,iy,iz] == 1:
                        active_cell = active_cell + (s.cells,)
                    elif self.locations[ix,iy,iz] == 0:
                        passive_cell = passive_cell + (s.cells,)
        
        active_set = p.Set(name='pset_Active', cells = active_cell) # part set for all active
        passive_set = p.Set(name='pset_Passive', cells = passive_cell) # part set for all passive
        
        p.SectionAssignment(region=p.sets['pset_Active'],
                            sectionName='Active_Section')
        p.SectionAssignment(region=p.sets['pset_Passive'],
                            sectionName='Passive_Section')
        
        # =========================================== #

    def create_steps(self):
        #
        self.m.ViscoStep(name='Actuation',
                         previous='Initial',
                         timePeriod=self.stepTime,
                         initialInc=self.iniInc,
                         minInc=self.stepTime*1.e-4,
                         maxInc=self.maxInc,
                         nlgeom=ON,
                         cetol=0.05,
                         maxNumInc=100)
        
        self.m.fieldOutputRequests['F-Output-1'].setValues(numIntervals=5)
        """
        self.m.steps['Actuation'].control.setValues(allowPropagation=OFF, 
            discontinuous=ON, resetDefaultValues=OFF, timeIncrementation=(8.0, 10.0, 
            9.0, 16.0, 10.0, 4.0, 12.0, 20.0, 6.0, 3.0, 50.0))
        
        self.m.steps['Actuation'].control.setValues(allowPropagation=
            OFF, displacementField=(0.5, 1.0, 0.0, 0.0, 0.02, 1e-05, 0.001, 1e-08, 1.0, 
            1e-05, 1e-08), hydrostaticFluidPressureField=DEFAULT, resetDefaultValues=
            OFF, rotationField=DEFAULT, temperatureField=DEFAULT)
        """

    def set_boundary_conditions(self):
        # 
        a = self.m.rootAssembly
        
        p = self.m.parts['Part-All']

        myIns = a.instances['Part-All-1']
        
        # === create set for faces ===
        #
        face_left = myIns.faces.getByBoundingBox(-0.001, -0.001, -0.001,
                                                 0.001, 0.001 + self.length_y, 0.001 + self.height)
        setFaceLeft = a.Set(name='set-FaceLeft', faces = face_left)
        #
        face_front = myIns.faces.getByBoundingBox(-0.001, -0.001, -0.001,
                                                  0.001 + self.length_x, 0.001, 0.001 + self.height)
        setFaceFront = a.Set(name='set-FaceFront', faces = face_front)
        #
        face_bottom = myIns.faces.getByBoundingBox(-0.001, -0.001, -0.001,
                                                   0.001 + self.length_x, 0.001 + self.length_y, 0.001)
        setFaceBottom = a.Set(name='set-FaceBottom', faces = face_bottom)
        #
        face_right = myIns.faces.getByBoundingBox(-0.001 + self.length_x, -0.001, -0.001,
                                                  0.001 + self.length_x, 0.001 + self.length_y, 0.001 + self.height)
        setFaceRight = a.Set(name='set-FaceRight', faces = face_right)
        #
        face_rear = myIns.faces.getByBoundingBox(-0.001, -0.001 + self.length_y, -0.001,
                                                 0.001 + self.length_x, 0.001 + self.length_y, 0.001 + self.height)
        setFaceRear = a.Set(name='set-FaceRear', faces = face_rear)
        #
        face_top = myIns.faces.getByBoundingBox(-0.001, -0.001, -0.001 + self.height,
                                                0.001 + self.length_x, 0.001 + self.length_y, 0.001 + self.height)
        setFaceTop = a.Set(name='set-FaceTop', faces = face_top)
        #
        face_mid = myIns.faces.getByBoundingBox(-0.001, -0.001, -0.001 + self.height/2,
                                                0.001 + self.length_x, 0.001 + self.length_y, 0.001 + self.height/2)
        setFaceTop = a.Set(name='set-FaceMid', faces = face_mid)
        # ============================
        
        # === create set for vertices ===
        v_x0y0z0 = myIns.vertices.getByBoundingBox(-0.001, -0.001, -0.001, 
                                                   0.001, 0.001, 0.001)
        a.Set(name='set-V_x0y0z0', vertices=v_x0y0z0)
        
        v_x0y1z0 = myIns.vertices.getByBoundingBox(-0.001, -0.001 + self.length_y, -0.001, 
                                                   0.001, 0.001 + self.length_y, 0.001)
        a.Set(name='set-V_x0y1z0', vertices=v_x0y1z0)
        
        v_x0y0z1 = myIns.vertices.getByBoundingBox(-0.001, -0.001, -0.001 + self.height, 
                                                   0.001, 0.001, 0.001 + self.height)
        a.Set(name='set-V_x0y0z1', vertices=v_x0y0z1)
        
        v_x1y0z0 = myIns.vertices.getByBoundingBox(-0.001 + self.length_x, -0.001, -0.001, 
                                                   0.001 + self.length_x, 0.001, 0.001)
        a.Set(name='set-V_x1y0z0', vertices=v_x1y0z0)
        
        v_x1y1z0 = myIns.vertices.getByBoundingBox(-0.001 + self.length_x, -0.001 + self.length_y, -0.001, 
                                                   0.001 + self.length_x, 0.001 + self.length_y, 0.001)
        a.Set(name='set-V_x1y1z0', vertices=v_x1y1z0)
        
        # mid-plane points
        v_x0y0zM = myIns.vertices.getByBoundingBox(-0.001, -0.001, -0.001 + self.height/2, 
                                                   0.001, 0.001, 0.001 + self.height/2)
        a.Set(name='set-V_x0y0zM', vertices=v_x0y0zM)
        
        v_x1y1zM = myIns.vertices.getByBoundingBox(-0.001 + self.length_x, -0.001 + self.length_y, -0.001 + self.height/2, 
                                                   0.001 + self.length_x, 0.001 + self.length_y, 0.001 + self.height/2)
        a.Set(name='set-V_x1y1zM', vertices=v_x1y1zM)
        
        v_x0y1zM = myIns.vertices.getByBoundingBox(-0.001, -0.001 + self.length_y, -0.001 + self.height/2, 
                                                   0.001, 0.001 + self.length_y, 0.001 + self.height/2)
        a.Set(name='set-V_x0y1zM', vertices=v_x0y1zM)
        
        v_x1y0zM = myIns.vertices.getByBoundingBox(-0.001 + self.length_x, -0.001, -0.001 + self.height/2, 
                                                   0.001 + self.length_x, 0.001, 0.001 + self.height/2)
        a.Set(name='set-V_x1y0zM', vertices=v_x1y0zM)
        
        # ===============================
        
        # ==== Initial BCs ====
        # self.m.DisplacementBC(name='BC-1-face_left',
                              # createStepName='Initial',
                              # region=a.sets['set-FaceLeft'],
                              # u1=0)
        
        self.m.DisplacementBC(name='BC-vert1', createStepName='Initial',
                              region=a.sets['set-V_x0y0zM'],
                              u1=0, u2=0, u3=0)
        self.m.DisplacementBC(name='BC-vert-4', createStepName='Initial',
                              region=a.sets['set-V_x1y1zM'],
                              u3=0)
        
        self.m.Equation(name='Constraint-vert-4', 
                        terms=((1.0, 'set-V_x1y1zM', 1), (-1.0, 'set-V_x1y1zM', 2)))
        
        self.m.Equation(name='Constraint-vert-23', 
                        terms=((1.0, 'set-V_x0y1zM', 3), (-1.0, 'set-V_x1y0zM', 3)))
        # =====================
        
        
    
    def set_initial_conditions_and_step_temperatures(self):
        #
        myIns = self.m.rootAssembly.instances['Part-All-1']
        # note: myIns.sets have the sets defined in Part. 
        #       self.m.rootAssembly.sets[''] only have the sets defined in Assembly.
        
        temp_field = self.m.Temperature(name='Initial-Temp',
                                        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS,
                                        distributionType=UNIFORM,
                                        createStepName='Initial',
                                        region=myIns.sets['All'],
                                        magnitudes=self.temp_low)
        
        self.m.TabularAmplitude(name='heating', timeSpan=STEP,
                                data=((0.0, self.temp_low/self.temp_high),
                                      (self.stepTime, 1.0),))

        temp_field.setValuesInStep(stepName='Actuation',
                                   magnitudes=self.temp_high, amplitude='heating')
        
        
    def mesh_structure(self):
        #
        part = self.m.parts['Part-All']
        part.setMeshControls(regions=part.cells,
                             elemShape=HEX)
        part.setElementType(elemTypes=(ElemType(elemCode=C3D8H,
                                                elemLibrary=STANDARD),),
                                 regions=part.sets['All'])
        
        global_seed = self.seed_xy
        part.seedPart(size=global_seed)
        
        edge_vertical = part.edges.getByBoundingBox(-0.001, -0.001, -0.001,
                                                    0.001, 0.001, 0.001 + self.height)
        part.seedEdgeBySize(edges = edge_vertical, size = self.seed_z, constraint = FIXED)
        
        part.generateMesh()

    def run_simulation(self, job_name, save_CAE):

        if save_CAE == True:
            mdb.saveAs(pathName=job_name[0:5]+'CAE')
            #return True
        
        mdb.Job(contactPrint=OFF, description='', echoPrint=OFF,
                explicitPrecision=SINGLE, historyPrint=OFF,
                memory=90, memoryUnits=PERCENTAGE, model=self.model_name,
                modelPrint=OFF, multiprocessingMode=DEFAULT, name=job_name,
                nodalOutputPrecision=SINGLE, numCpus=self.num_cpus, numDomains=self.num_cpus,
                parallelizationMethodExplicit=DOMAIN, scratch='', type=ANALYSIS,
                userSubroutine='')
        mdb.jobs[job_name].submit()
        mdb.jobs[job_name].waitForCompletion()
        
        # Below is used to check the job.
        with open(job_name+'.log','r') as f_log:
            f_log.seek(-10,2)
            status_log = f_log.read(9)
        if status_log == 'COMPLETED':
            aborted = False
        else:
            aborted = True
        
        return aborted

    @staticmethod
    def csv_output(aborted, job_name):
    
        path_current = os.getcwd()
        out = os.path.join(os.path.dirname(path_current),'CSV_island','outs-' + job_name.split('.')[0] + '.csv')
        #out = '..\\CSV_random\\outs-' + job_name.split('.')[0] + '.csv'
    
        if aborted:
            print('*Simulation not completed for error_code_1')
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
        


   
def main(gene_file, job_name):
    # Parameters given here, not in a separate file
    
    # number of cpus, 1 or 4
    num_cpus = 1
    
    # Geometric parameters units in mm
    length_x = 40.0 # x-length of plate
    length_y = 40.0 # y-length of plate
    height = 1.0 # height/thickness of plate in z-direction
    
    # Actuation strain, reduced due to increased element number
    active_strain = 0.1/2
    # Note:
    # Due to reduced strain, also modified includes 
    # stepTime, iniInc, maxInc, minInc, output's numIntervals
    
    nx,ny,nz,locations = read_genome(gene_file)
    
    voxel_size_x = length_x / (nx + 0)
    voxel_size_y = length_y / (ny + 0)
    voxel_size_z = height / (nz + 0)
    
    #
    # temperature profile to do expansion
    #
    temp_low = 0.0
    temp_high = 100.0
    alpha = active_strain/(temp_high-temp_low) # thermal expansion coefficient
    
    stepTime = 5.0/2
    iniInc = stepTime/5.0*2
    maxInc = stepTime/2.0*2
    
    # Single-voxel element number: 3*3*2 (for 15x15x2 voxel case)
    seed_xy = min(voxel_size_x/3.0, voxel_size_y/3.0, voxel_size_z*2.0) # for x-y plane mesh size, global seed size control
    seed_z = voxel_size_z/2.0
    
    # call class, class constructor runs functions
    model1 = Active2DStructure(length_x, length_y, height, nx, ny, nz,
                               voxel_size_x, voxel_size_y, voxel_size_z,
                               locations, num_cpus,
                               alpha, temp_low, temp_high,
                               seed_xy, seed_z,
                               stepTime, iniInc, maxInc)
    
    save_CAE = False
    aborted = model1.run_simulation(job_name,save_CAE)
    model1.csv_output(aborted,job_name)


if __name__ == '__main__':

    args = sys.argv

    for arg in args:
        if '.txt' in arg:
            index = args.index(arg)

    my_args = args[index:]

    gene_file_in = my_args[0]
    job_name_in = my_args[1]

    main(gene_file_in, job_name_in)

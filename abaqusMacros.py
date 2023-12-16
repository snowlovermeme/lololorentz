# -*- coding: mbcs -*-
# Do not delete the following import lines
"""
Author: Sun Hao (Hunan University, Machine college)
the code was written for electromagnetic crrimping, including following main steps:
    1. create geometry
        1): create fibers which subject certain probobility distribution
        2): create the terminal containing the fibers.
    2. materials
        1): the Johnson-cook constitutive model, the parameters was collected from other papers;
        2): the J-C but the damage model to simulation the cutting down behavior of fibers during crimping
    3. boundary conditions
        1): electromagneric force distribution obtained from other sofrware (comsol or maxwell)   Planning ...
            (1): using the field distribution in abaqus edtion 2024 ?
            (2): using the subrutine UDF and interpolation ?
    4. others:
        1): explicit solver;
        2): contact?
    5. futher expploring aspect:
        1): simplifying the fibers into beam model
        2): axial symmetry model: based on the impact raidla mechanical response of the fibers
"""
from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import numpy as np
# ==================== functions =====================
# create a matrix under defined radius and density, R and N.
# Generating meshod was dart drawing: "Fast Poisson Disk Sampling in Arbitrary Dimensions" by Robert Bridson in 2007

# check the collition between any two points
def fiber_locations(r, R, N):
    r = 2 * r
    d = r / np.sqrt(2)
    k = N
    width = 2 * R
    height = 2 * R
    nx = int(width / d) + 1
    ny = int(height / d) + 1

    occupied = np.zeros((ny, nx))
    occupied_coord = np.zeros((ny, nx, 2))
    active_list = []
    sampled = []

    relative = np.array([[-1, 2], [0, 2], [1, 2],
                         [-2, 1], [-1, 1], [0, 1], [1, 1], [2, 1],
                         [-2, 0], [-1, 0], [1, 0], [2, 0],
                         [-2, -1], [-1, -1], [0, -1], [1, -1], [2, -1],
                         [-1, -2], [0, -2], [1, -2]])
    np.random.seed(0)
    x, y = np.random.rand() * width, np.random.rand() * height
    idx_x, idx_y = int(x / d), int(y / d)
    occupied[idx_y, idx_x] = 1
    occupied_coord[idx_y, idx_x] = (x, y)
    active_list.append((x, y))
    sampled.append((x, y))

    sampled_idx = 0
    while len(active_list) > 0:

        idx = np.random.choice(np.arange(len(active_list)))
        ref_x, ref_y = active_list[idx]
        radius = (np.random.rand(k) + 1) * r
        theta = np.random.rand(k) * np.pi * 2
        candidate = radius * np.cos(theta) + ref_x, radius * np.sin(theta) + ref_y
        flag_out = False
        for _x, _y in zip(*candidate):
            if _x < 0 or _x > width or _y < 0 or _y > height:
                continue
            # other geo constraints
            flag = True
            idx_x, idx_y = int(_x / d), int(_y / d)
            if occupied[idx_y, idx_x] != 0:
                continue
            else:
                neighbours = relative + np.array([idx_x, idx_y])
            for cand_x, cand_y in neighbours:
                if cand_x < 0 or cand_x >= nx or cand_y < 0 or cand_y >= ny:
                    continue
                if occupied[cand_y, cand_x] == 1:
                    cood = occupied_coord[cand_y, cand_x]
                    if (_x - cood[0]) ** 2 + (_y - cood[1]) ** 2 < r ** 2:
                        flag = False
                        break
            if flag:
                flag_out = True
                occupied[idx_y, idx_x] = 1
                occupied_coord[idx_y, idx_x] = (_x, _y)
                sampled.append((_x, _y))
                active_list.append((_x, _y))
                sampled_idx += 1
                break
        if not flag_out:
            active_list.pop(idx)

    # coordinate transforming
    returned_samples = []
    for coord in sampled:
        x, y = coord
        x, y = x - R , y - R
        if x**2 + y**2 < R - (r/2) * 1.1:
            returned_samples.append([x, y])
    return returned_samples



# =============== parameter setting ===================
tube_r_in=5;
n_fiber=30;
fiber_r=0.1;
thickness = 2;
geo_type = 'wire'
L_tube = 50.0
N_partions = 40
# =============== create one fiber ===================
# 1. create sketch
if geo_type == 'solid':
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    # define the radius of one fiber
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, fiber_r))
    p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D,
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Part-1']
    # define the length of one fiber
    p.BaseSolidExtrude(sketch=s, depth=60.0)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Part-1']
    del mdb.models['Model-1'].sketches['__profile__']
else:
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
                                                 sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.Line(point1=(-50.0, 0.0), point2=(50.0, 0.0))
    s1.HorizontalConstraint(entity=g[2], addUndoState=False)
    p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D,
                                   type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Part-1']
    p.BaseWire(sketch=s1)
    s1.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Part-1']
    del mdb.models['Model-1'].sketches['__profile__']
# 2. create mesh for one fiber
p = mdb.models['Model-1'].parts['Part-1']
# # define the mesh size (or the seed gap) for one fiber
p.seedPart(size=0.14, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Part-1']
# generate the mesh!
p.generateMesh()
# =============== create strands ===================
samples = fiber_locations(fiber_r, tube_r_in, n_fiber)
# 1. create assembly
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
# 2. import the fiber into the assembly (Note: select dependent = ON for mesh)
fiber_ind = 1
for loc in samples:
    a.Instance(name='fiber'+str(fiber_ind), part=p, dependent=ON)
    # 3. translate the fiber imported into the prescribed location!
    a.translate(instanceList=('fiber'+str(fiber_ind), ), vector=(0, loc[1], loc[0]))
    fiber_ind += 1
# ============= craete tube =============
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(tube_r_in, 0.0))
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(tube_r_in + thickness, 0.0))
p = mdb.models['Model-1'].Part(name='Part-2', dimensionality=THREE_D,
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-2']
p.BaseSolidExtrude(sketch=s, depth=L_tube)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-2']
del mdb.models['Model-1'].sketches['__profile__']
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Part-2']
a.Instance(name='Part-2-1', part=p, dependent=ON)
p = mdb.models['Model-1'].parts['Part-2']
# ============= partions the tube =============
p = mdb.models['Model-1'].parts['Part-2']
c = p.cells
N_planes = 0
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
N_planes += 2
z_partition_plane = np.linspace(25.0-10.0, 25.0+10.0, num=N_partions+1)
for i in range(N_partions+1):
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=z_partition_plane[i])
    N_planes += 1
d = p.datums
for i in range(N_planes):
    plane_name = 'Datum plane-' + str(i+1)
    id = p.features[plane_name].id
    c = p.cells
    p.PartitionCellByDatumPlane(datumPlane=d[id], cells=c)


p.seedPart(size=0.5, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()

# ================ create step, amp and load ================
mdb.models['Model-1'].ExplicitDynamicsStep(name='Step-1', previous='Initial',
        improvedDtMethod=ON)
mdb.models['Model-1'].TabularAmplitude(name='Amp-1', timeSpan=STEP,
        smooth=SOLVER_DEFAULT, data=((0.0, 0.0), (1.0, 1.0)))
r_out = tube_r_in + thickness
verti = r_out/np.sqrt(2)
x_p, y_p = [verti, verti, -verti, -verti], [-verti, verti, -verti, verti]
for i in range(N_partions):
    a = mdb.models['Model-1'].rootAssembly
    f = a.instances['Part-2-1'].faces
    z_0 = z_partition_plane[i] + 0.1
    faces = f.findAt(
        ((verti, -verti, z_0),),
        ((-verti, verti, z_0),),
        ((-verti, -verti, z_0),),
        ((verti, verti, z_0),),
    )
    region = regionToolset.Region(side1Faces=faces)
    mdb.models['Model-1'].Pressure(name='Load-'+str(i), createStepName='Step-1',
            region=region, distributionType=UNIFORM, field='', magnitude=970800+(i+1),
            amplitude='Amp-1')

# ================ adjust tube position ================
a1 = mdb.models['Model-1'].rootAssembly
a1.rotate(instanceList=('Part-2-1', ), axisPoint=(0.0, 0.0, 0.0),
    axisDirection=(0.0, 1.0, 0.0), angle=90.0)
a1 = mdb.models['Model-1'].rootAssembly
a1.translate(instanceList=('Part-2-1', ), vector=(-25.0, 0.0, 0.0))


# ====================== define materials ======================
mdb.models['Model-1'].Material(name='copper')
mdb.models['Model-1'].materials['copper'].Density(table=((8.94e-09, ), ))
mdb.models['Model-1'].materials['copper'].Elastic(table=((117800.0, 0.3), ))
mdb.models['Model-1'].materials['copper'].Plastic(hardening=JOHNSON_COOK,
    scaleStress=None, table=((150.0, 170.0, 0.34, 1.09, 923.0, 300.0), ))
mdb.models['Model-1'].materials['copper'].plastic.RateDependent(
    type=JOHNSON_COOK, table=((0.025, 1.0), ))
mdb.models['Model-1'].materials['copper'].JohnsonCookDamageInitiation(table=((
    0.54, 4.89, -3.03, 0.014, 1.12, 923.0, 300.0, 1.0), ))
mdb.models['Model-1'].materials['copper'].johnsonCookDamageInitiation.DamageEvolution(
    type=DISPLACEMENT, table=((0.0, ), ))
# ====================== define and assign sections ======================
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1',
        material='Material-1', thickness=None)
mdb.models['Model-1'].CircularProfile(name='Profile-1', r=0.5)
mdb.models['Model-1'].BeamSection(name='Section-2',
    integration=DURING_ANALYSIS, poissonRatio=0.0, profile='Profile-1',
    material='Material-1', temperatureVar=LINEAR,
    consistentMassMatrix=False)
p = mdb.models['Model-1'].parts['Part-1']
region = p.Set(edges=p.edges, name='fiber_set')
p.SectionAssignment(region=region, sectionName='Section-2', offset=0.0,
        offsetType=MIDDLE_SURFACE, offsetField='',
thicknessAssignment = FROM_SECTION)
p = mdb.models['Model-1'].parts['Part-2']
region = p.Set(cells=p.cells, name='tube_cell')
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0,
        offsetType=MIDDLE_SURFACE, offsetField='',
thicknessAssignment = FROM_SECTION)
p = mdb.models['Model-1'].parts['Part-1']
e = p.edges
edges = e.getSequenceFromMask(mask=('[#1 ]', ), )
region=regionToolset.Region(edges=edges)
p = mdb.models['Model-1'].parts['Part-1']
p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(0.0, 0.0,
    -1.0))
# ================ define contact behavior ================
mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
    formulation=FRICTIONLESS)
mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
    pressureOverclosure=HARD, allowSeparation=ON,
    constraintEnforcementMethod=DEFAULT)
mdb.models['Model-1'].ContactExp(name='Int-1', createStepName='Step-1')
mdb.models['Model-1'].interactions['Int-1'].includedPairs.setValuesInStep(
    stepName='Step-1', useAllstar=ON)
mdb.models['Model-1'].interactions['Int-1'].contactPropertyAssignments.appendInStep(
    stepName='Step-1', assignments=((GLOBAL, SELF, 'IntProp-1'), ))
# ================ create analysis job ================
mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, explicitPrecision=SINGLE,
        nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',
        resultsFormat=ODB, numDomains=8, activateLoadBalancing=False,
        numThreadsPerMpiProcess=1, multiprocessingMode=DEFAULT, numCpus=8)

"""
=======================================================
Shear Elastic Buckling of Plate with Flanged Hole
Developed by Abdulaziz Alabdulwahab
Email: alabdulwahab01@gmail.com
Year: 2025
All Rights Reserved. Do not distribute without permission.
=======================================================
"""
# ----------------------------
# 0) Collect inputs via popup
# ----------------------------
from abaqus import *
from abaqusConstants import *
import sys
import math
from math import sqrt
import numpy as np
import time


def calculate_shear_buckling(a, d1, dh, s, t, E, nu, shape):
    # force Python 2 in Abaqus to use float math
    a   = float(a)
    d1  = float(d1)
    dh  = float(dh)
    s   = float(s)
    t   = float(t)
    E  = float(E)
    nu  = float(nu)

    # Hole shape properties
    alpha_h_table = {"Squ": 1.00, "Cir": 0.68, "Dia": 0.43}
    Cs_square = 1

    coeffs_ks = {
        "Cir": {"A": 157.93, "B": 662.63, "C": 47.17, "D": 10.36},
        "Squ": {"A": 720.07, "B": 1242.27, "C": -38.94, "D": 130.55},
        "Dia": {"A": 120.24, "B": 174.38, "C": 3.96, "D": 0},
    }

    coeffs_teq = {
        "Cir": {"A": 35.57, "B": 0.37, "C": 1.13, "D": 0.59},
        "Squ": {"A": 168.33, "B": 8.18, "C": 19.28, "D": 12.85},
        "Dia": {"A": 12.54, "B": 1.04, "C": 1.18, "D": 0.26},
    }

    alpha_h = alpha_h_table[shape]
    ks_dict = coeffs_ks[shape]
    A_ks = ks_dict["A"]
    B_ks = ks_dict["B"]
    C_ks = ks_dict["C"]
    D_ks = ks_dict["D"]

    teq_dict = coeffs_teq[shape]
    A_teq = teq_dict["A"]
    B_teq = teq_dict["B"]
    C_teq = teq_dict["C"]
    D_teq = teq_dict["D"]

    # 1) k0
    k0 = 5.34 + 4.0 / (a / d1) ** 2

    # 2) u
    u = alpha_h * (dh / d1) ** 2 * 10.0

    # 3) kh
    kh = (k0 / 9.34) * (11.2 * u * d1 / a + 9.5 * u**2) / (1.0 + 2.1 * u + 0.95 * u**2)

    # 4) kv (unstiffened)
    kv_unstiffened = k0 - kh

    # 5) Cs
    if shape == "Squ":
        Cs = Cs_square
    elif shape == "Cir":
        Cs = (k0 / 9.34) ** (1.0 + (d1 / a))
    else:  # Dia
        Cs = (k0 / 9.34) ** (2.5 + (d1 / a))

    # 6) ks
    ks_num = -A_ks * (s / d1) * (dh / d1) + B_ks * (dh / d1) ** 2
    ks_den = 1.0 + C_ks * (s / d1) * (dh / d1) + D_ks * (dh / d1) ** 2
    ks = Cs * sqrt(s / d1) * (ks_num / ks_den)

    kv_stiffened = k0 - kh + ks

    # 7) teq
    teq_num = A_teq * (s / t) * (dh / d1) - B_teq * (s / t) ** 2
    teq_den = 1.0 - C_teq * (s / t) * (dh / d1) + D_teq * (s / t) ** 2
    teq = (1.0 - (t / d1) * Cs * (teq_num / teq_den)) * t

    def tau_cr(kv, thickness):
        return ((kv * np.pi**2 * E) / (12.0 * (1.0 - nu**2))) * (thickness / d1) ** 2

    tau_unstiffened = tau_cr(kv_unstiffened, t)
    tau_stiffened   = tau_cr(kv_stiffened,   teq)

    return {
        "shape": shape,
        "k0": round(k0, 3),
        "kh": round(kh, 3),
        "ks": round(ks, 3),
        "kv_unstiffened": round(kv_unstiffened, 3),
        "kv_stiffened": round(kv_stiffened, 3),
        "teq (mm)": round(teq, 3),
        "tau_cr_unstiffened (MPa)": round(tau_unstiffened, 3),
        "tau_cr_stiffened (MPa)": round(tau_stiffened, 3),
    }


# ----------------------------
# 1) Collect inputs via popup
# ----------------------------
pars = (
    ('Hole shape (Cir/Squ/Dia):',        'Cir'),
    ('Plate length a (mm):',             '400'),
    ('Plate width  d1 (mm):',            '200'),
    ('Thickness t (mm):',                '1'),
    ('Hole size dh (mm, 0 = none):',     '100'),
    ('Stiffener depth s (mm):',          '12'),
    ('Stiffener radius ri (mm):',        '0'),   # for the lip / stiffener
    ('Hole corner radius rc (mm):',      '0'),    # NEW: for square & diamond hole itself
    ('Global seed size (mm):',           '5'),
)

vals = getInputs(
    fields=pars,
    dialogTitle='Shear Elastic Buckling of Plate with Flanged Hole',
    label=(
        '\nShear Elastic Buckling of Plate with Flanged Hole'
        '\n(c) 2025 All Rights Reserved'
        '\nDeveloped by Abdulaziz Alabdulwahab'
        '\nEmail: alabdulwahab01@gmail.com''\n\n'
        'Enter geometry / meshing parameters:\n'
    )
)

if vals is None:
    sys.exit()

shape_str = vals[0].strip().upper() if vals[0] else 'CIR'
try:
    Lx         = float(vals[1]) if vals[1] else 0.0
    Ly         = float(vals[2]) if vals[2] else 0.0
    Thickness  = float(vals[3]) if vals[3] else 0.0
    HoleSize   = float(vals[4]) if vals[4] else 0.0
    LipDepth   = float(vals[5]) if vals[5] else 0.0
    Ri         = float(vals[6]) if vals[6] else 0.0     # stiffener fillet
    Rc         = float(vals[7]) if vals[7] else 0.0     # hole corner fillet (Squ/Dia)
    SeedSize   = float(vals[8]) if vals[8] else 0.0
except:
    sys.exit('Invalid numeric input.')

if Lx <= 0 or Ly <= 0 or Thickness <= 0 or SeedSize <= 0:
    sys.exit('Lx, Ly, thickness, and seed size must be > 0.')

if HoleSize <= 0.0:
    HoleSize = 0.0
    LipDepth = 0.0
    # keep Ri as is (it may still be used on stiffener), but rc is meaningless
    # Rc will just be ignored


# =====================================================
# 2) Circular hole builder
# =====================================================
def build_plate_cir(lx, ly, thickness, hole_diameter, lip_depth, ri, seed_size):
    model_name = 'Model-1'
    job_name   = 'Job-1'
    part_name  = 'ShellPart'

    myModel = mdb.Model(name=model_name)
    myPart  = myModel.Part(name=part_name, dimensionality=THREE_D, type=DEFORMABLE_BODY)

    mySketch = myModel.ConstrainedSketch(name='ShellSketch', sheetSize=10.0*max(lx, ly, 1.0))
    mySketch.rectangle(point1=(0.0, 0.0), point2=(lx, ly))
    if hole_diameter > 0.0:
        mySketch.CircleByCenterPerimeter(center=(lx/2.0, ly/2.0),
                                         point1=(lx/2.0, ly/2.0 + hole_diameter/2.0))
    myPart.BaseShell(sketch=mySketch)
    mySketch.unsetPrimaryObject()

    if hole_diameter > 0.0 and lip_depth > 0.0:
        lipSketch = myModel.ConstrainedSketch(name='LipSketch', sheetSize=10.0*max(lx, ly, 1.0))
        lipSketch.CircleByCenterPerimeter(center=(lx/2.0, ly/2.0),
                                          point1=(lx/2.0, ly/2.0 + hole_diameter/2.0))
        try:
            edge = myPart.edges.findAt(((lx/2.0, ly/2.0 + hole_diameter/2.0, 0.0),))
            myPart.ShellExtrude(sketchPlane=myPart.faces[0], sketchUpEdge=edge[0],
                                sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
                                sketch=lipSketch, depth=lip_depth, flipExtrudeDirection=OFF)
            lipSketch.unsetPrimaryObject()
            if ri > 0.0:
                edge = myPart.edges.findAt(((lx/2.0, ly/2.0 + hole_diameter/2.0, 0.0),))
                myPart.Round(radius=ri, edgeList=(edge[0],))
        except:
            pass

    myModel.Material(name='Material-1')
    myModel.materials['Material-1'].Elastic(table=((200000.0, 0.3),))
    myModel.HomogeneousShellSection(name='Section-1', preIntegrate=OFF, thicknessType=UNIFORM,
                                    thickness=thickness, material='Material-1', thicknessField='',
                                    idealization=NO_IDEALIZATION,
                                    poissonDefinition=DEFAULT, thicknessModulus=None,
                                    temperature=GRADIENT, useDensity=OFF,
                                    integrationRule=SIMPSON, numIntPts=5)
    myPart.Set(faces=myPart.faces, name='Set-Shell')
    myPart.SectionAssignment(region=myPart.sets['Set-Shell'], sectionName='Section-1')

    myAsm = myModel.rootAssembly
    myInst = myAsm.Instance(name='ShellInstance', part=myPart, dependent=ON)

    myModel.BuckleStep(name='Step-1', previous='Initial', numEigen=5, vectors=30, maxIterations=300)

    myAsm.Set(name='ShellEdges',
              edges=myInst.edges.findAt(((0.0, ly/2.0, 0.0),),
                                        ((lx, ly/2.0, 0.0),),
                                        ((lx/2.0, 0.0, 0.0),),
                                        ((lx/2.0, ly, 0.0),),))
    myAsm.Set(name='Corners',
              vertices=myInst.vertices.findAt(((0.0, 0.0, 0.0),),
                                              ((0.0, ly, 0.0),),))
    myModel.DisplacementBC(name='BC-1', createStepName='Step-1',
                           region=myAsm.sets['ShellEdges'], u3=0.0)
    myModel.DisplacementBC(name='BC-2', createStepName='Step-1',
                           region=myAsm.sets['Corners'], u1=0.0, u2=0.0)

    myAsm.Surface(name='HorizontalEdges',
                  side1Edges=myInst.edges.findAt(((lx/2.0, 0.0, 0.0),),
                                                 ((lx/2.0, ly, 0.0),),))
    myAsm.Surface(name='VerticalEdges',
                  side1Edges=myInst.edges.findAt(((0.0, ly/2.0, 0.0),),
                                                 ((lx,  ly/2.0, 0.0),),))
    myModel.ShellEdgeLoad(name='Load-1', createStepName='Step-1',
                          region=myAsm.surfaces['HorizontalEdges'], magnitude=1.0,
                          distributionType=UNIFORM, field='', localCsys=None, traction=SHEAR)
    myModel.ShellEdgeLoad(name='Load-2', createStepName='Step-1',
                          region=myAsm.surfaces['VerticalEdges'], magnitude=-1.0,
                          distributionType=UNIFORM, field='', localCsys=None, traction=SHEAR)

    if hole_diameter > 0.0:
        try:
            Datum1 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, 0.0),
                                                 point2=(lx, 0.0, 0.0), isDependent=False).id
            Datum2 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, 0.0),
                                                 point2=(0.0, ly, 0.0), isDependent=False).id
            Datum3 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, 0.0),
                                                 point2=(lx + hole_diameter + 2*ri + 2*seed_size, 0.0, 0.0),
                                                 isDependent=False).id
            Datum4 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, 0.0),
                                                 point2=(lx - hole_diameter - 2*ri - 2*seed_size, 0.0, 0.0),
                                                 isDependent=False).id
            Datum5 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, 0.0),
                                                 point2=(0.0, ly + hole_diameter + 2*ri + 2*seed_size, 0.0),
                                                 isDependent=False).id
            Datum6 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, 0.0),
                                                 point2=(0.0, ly - hole_diameter - 2*ri - 2*seed_size, 0.0),
                                                 isDependent=False).id

            myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum1], faces=myPart.faces)
            myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum2], faces=myPart.faces)
            myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum3], faces=myPart.faces)
            myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum4], faces=myPart.faces)
            myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum5], faces=myPart.faces)
            myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum6], faces=myPart.faces)

            pickedFaces = myPart.faces.getByBoundingBox(
                xMin=lx/2.0 - hole_diameter/2.0 - ri - seed_size - 0.1,
                yMin=ly/2.0 - hole_diameter/2.0 - ri - seed_size - 0.1,
                zMin=-0.1,
                xMax=lx/2.0 + hole_diameter/2.0 + ri + seed_size + 0.1,
                yMax=ly/2.0 + hole_diameter/2.0 + ri + seed_size + 0.1,
                zMax=0.1
            )
            myPart.PartitionFaceByShortestPath(
                point1=(lx/2.0 - hole_diameter/2.0 - ri - seed_size,
                        ly/2.0 - hole_diameter/2.0 - ri - seed_size, 0.0),
                point2=(lx/2.0 + hole_diameter/2.0 + ri + seed_size,
                        ly/2.0 + hole_diameter/2.0 + ri + seed_size, 0.0),
                faces=pickedFaces
            )
            pickedFaces = myPart.faces.getByBoundingBox(
                xMin=lx/2.0 - hole_diameter/2.0 - ri - seed_size - 0.1,
                yMin=ly/2.0 - hole_diameter/2.0 - ri - seed_size - 0.1,
                zMin=-0.1,
                xMax=lx/2.0 + hole_diameter/2.0 + ri + seed_size + 0.1,
                yMax=ly/2.0 + hole_diameter/2.0 + ri + seed_size + 0.1,
                zMax=0.1
            )
            myPart.PartitionFaceByShortestPath(
                point1=(lx/2.0 + hole_diameter/2.0 + ri + seed_size,
                        ly/2.0 - hole_diameter/2.0 - ri - seed_size, 0.0),
                point2=(lx/2.0 - hole_diameter/2.0 - ri - seed_size,
                        ly/2.0 + hole_diameter/2.0 + ri + seed_size, 0.0),
                faces=pickedFaces
            )

            if lip_depth > ri:
                Datum7 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, ri),
                                                     point2=(0.0, 0.0, lip_depth),
                                                     isDependent=False).id
                myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum7], faces=myPart.faces)

            pickedFaces = myPart.faces.getByBoundingBox(
                xMin=lx/2.0 - hole_diameter/2.0 - ri - seed_size - 0.1,
                yMin=ly/2.0 - hole_diameter/2.0 - ri - seed_size - 0.1,
                zMin=-0.1,
                xMax=lx/2.0 + hole_diameter/2.0 + ri + seed_size + 0.1,
                yMax=ly/2.0 + hole_diameter/2.0 + ri + seed_size + 0.1,
                zMax=0.1
            )
            myPart.setMeshControls(regions=pickedFaces, technique=SWEEP)
        except:
            pass

    myPart.seedPart(size=seed_size, deviationFactor=0.05, minSizeFactor=0.1)
    myPart.generateMesh()

    myJob = mdb.Job(name=job_name, model=model_name, description='',
                    type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None,
                    memory=4, memoryUnits=GIGA_BYTES, getMemoryFromAnalysis=True,
                    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
                    scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,
                    numDomains=4, activateLoadBalancing=False, multiprocessingMode=MPI,
                    numCpus=4, numGPUs=0)
    myJob.writeInput(consistencyChecking=OFF)
    # Viewing model
    session.viewports['Viewport: 1'].setValues(displayedObject=myPart)
    session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
    session.viewports['Viewport: 1'].view.rotate(xAngle=10, yAngle=-40, zAngle=0,
                                                 mode=MODEL)
    session.viewports['Viewport: 1'].view.fitView()
    time.sleep(1)
    myJob.submit()
    myJob.waitForCompletion()

    # --- get eigenvalue from ODB and return it ---
    myOdb = session.openOdb(name=job_name + '.odb')
    buck_frames = myOdb.steps['Step-1'].frames
    eig_text = buck_frames[1].description
    # Viewing results
    session.viewports['Viewport: 1'].setValues(displayedObject=myOdb)
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].setValues(
        displayedObject=myOdb)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
    session.viewports['Viewport: 1'].view.rotate(xAngle=10, yAngle=-40, zAngle=0,
                                                 mode=MODEL)
    session.viewports['Viewport: 1'].view.fitView()
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
        CONTOURS_ON_DEF,))
    # parse the number at the end
    eig_val = abs(float(eig_text.split('=')[-1]))
    return eig_val



# =====================================================
# 3) Square hole builder
# =====================================================
def build_plate_squ(lx, ly, thickness, hole_side, lip_depth, ri, rc, seed_size):
    model_name = 'Model-1'
    job_name   = 'Job-1'
    myModel = mdb.Model(name=model_name)
    myPart  = myModel.Part(name='ShellPart', dimensionality=THREE_D, type=DEFORMABLE_BODY)

    mySketch = myModel.ConstrainedSketch(name='ShellSketch', sheetSize=2000.0)
    mySketch.rectangle(point1=(0, 0), point2=(lx, ly))
    if hole_side > 0:
        bx = hole_side
        by = hole_side
        mySketch.rectangle(point1=(lx/2 - bx/2, ly/2 - by/2),
                           point2=(lx/2 + bx/2, ly/2 + by/2))
    else:
        bx = 0.0
        by = 0.0
    myPart.BaseShell(sketch=mySketch)
    mySketch.unsetPrimaryObject()

    # Lip (stiffener)
    if hole_side > 0 and lip_depth > 0:
        right_side = myPart.edges.findAt(((lx, ly/2.0, 0),))
        lipTransform = myPart.MakeSketchTransform(sketchPlane=myPart.faces[0],
                                                  sketchUpEdge=right_side[0],
                                                  sketchPlaneSide=SIDE1,
                                                  sketchOrientation=RIGHT,
                                                  origin=(lx/2.0, ly/2.0, 0.0))
        lipSketch = myModel.ConstrainedSketch(name='CutSketch', sheetSize=200.0, transform=lipTransform)
        lipSketch.rectangle(point1=(bx/2.0, by/2.0), point2=(-bx/2.0, -by/2.0))
        edge = myPart.edges.findAt(((lx, ly/2.0, 0),))
        myPart.ShellExtrude(sketchPlane=myPart.faces[0], sketchUpEdge=edge[0],
                            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
                            sketch=lipSketch, depth=lip_depth, flipExtrudeDirection=OFF)
        lipSketch.unsetPrimaryObject()

        # Hole corners fillet rc
        if rc > 0 and lip_depth > 0:
            stiffener_corners = myPart.edges.findAt(
                ((lx/2 + bx/2, ly/2 + by/2, lip_depth - 0.1),),
                ((lx/2 - bx/2, ly/2 - by/2, lip_depth - 0.1),),
                ((lx/2 + bx/2, ly/2 - by/2, lip_depth - 0.1),),
                ((lx/2 - bx/2, ly/2 + by/2, lip_depth - 0.1),),
            )
            myPart.Round(radius=rc, edgeList=(stiffener_corners[0], stiffener_corners[1],
                                              stiffener_corners[2], stiffener_corners[3],))

    # stiffener corners fillet ri
    if lip_depth > 0 and ri > 0:
        see = 5.0  # same spirit as your original square code
        hole_edges = myPart.edges.findAt(
            ((lx/2 + bx/2, ly/2, 0),),
            ((lx/2 - bx/2, ly/2, 0),),
            ((lx/2, ly/2 + by/2, 0),),
            ((lx/2, ly/2 - by/2, 0),),
            ((lx/2 + bx/2 - (rc - rc*math.cos(math.radians(45)) ),
              ly/2 + by/2 - (rc - rc*math.cos(math.radians(45)) ), 0),),
            ((lx/2 + bx/2 - (rc - rc*math.cos(math.radians(45)) ),
              ly/2 - by/2 + (rc - rc*math.cos(math.radians(45)) ), 0),),
            ((lx/2 - bx/2 + (rc - rc*math.cos(math.radians(45)) ),
              ly/2 + by/2 - (rc - rc*math.cos(math.radians(45)) ), 0),),
            ((lx/2 - bx/2 + (rc - rc*math.cos(math.radians(45)) ),
              ly/2 - by/2 + (rc - rc*math.cos(math.radians(45)) ), 0),),
        )
        myPart.Round(radius=ri, edgeList=(hole_edges[0],hole_edges[1],hole_edges[2],hole_edges[3],
                                              hole_edges[4],hole_edges[5],hole_edges[6],hole_edges[7],))

    # Material & section
    myModel.Material(name='Material-1')
    myModel.materials['Material-1'].Elastic(table=((200000, 0.3),))
    myModel.HomogeneousShellSection(name='Section-1', preIntegrate=OFF, thicknessType=UNIFORM,
                                    thickness=thickness, material='Material-1', thicknessField='',
                                    idealization=NO_IDEALIZATION,
                                    poissonDefinition=DEFAULT, temperature=GRADIENT,
                                    useDensity=OFF, integrationRule=SIMPSON, numIntPts=5)
    myPart.Set(faces=myPart.faces, name='Set-Shell')
    myPart.SectionAssignment(region=myPart.sets['Set-Shell'], sectionName='Section-1')

    # Assembly
    myAsm = myModel.rootAssembly
    myInst = myAsm.Instance(name='ShellInstance', part=myPart, dependent=ON)

    # Step
    myModel.BuckleStep(name='Step-1', previous='Initial', numEigen=2, vectors=30, maxIterations=300)

    # BCs
    myAsm.Set(name='ShellEdges', edges=myInst.edges.findAt(
        ((0.0, ly/2, 0.0),), ((lx, ly/2, 0.0),),
        ((lx/2, 0.0, 0.0),), ((lx/2, ly, 0.0),),))
    myAsm.Set(name='Corners', vertices=myInst.vertices.findAt(
        ((0.0, 0.0, 0.0),), ((0.0, ly, 0.0),),))
    myModel.DisplacementBC(name='BC-1', createStepName='Step-1',
                           region=myAsm.sets['ShellEdges'], u3=0.0)
    myModel.DisplacementBC(name='BC-2', createStepName='Step-1',
                           region=myAsm.sets['Corners'], u1=0.0, u2=0.0)

    # Loads
    myAsm.Surface(name='HorizontalEdges', side1Edges=myInst.edges.findAt(
        ((lx/2, 0.0, 0.0),), ((lx/2, ly, 0.0),),))
    myAsm.Surface(name='VerticalEdges', side1Edges=myInst.edges.findAt(
        ((0.0, ly/2, 0.0),), ((lx, ly/2, 0.0),),))
    myModel.ShellEdgeLoad(name='Load-1', createStepName='Step-1',
                          region=myAsm.surfaces['HorizontalEdges'], magnitude=1.0,
                          distributionType=UNIFORM, field='', localCsys=None, traction=SHEAR)
    myModel.ShellEdgeLoad(name='Load-2', createStepName='Step-1',
                          region=myAsm.surfaces['VerticalEdges'], magnitude=-1.0,
                          distributionType=UNIFORM, field='', localCsys=None, traction=SHEAR)

    # Partition & mesh (simplified but same idea)
    if bx > 0:
        Datum1 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, 0.0),
                                             point2=(lx + bx + 4*seed_size, 0.0, 0.0),
                                             isDependent=False).id
        Datum2 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, 0.0),
                                             point2=(lx - bx - 4*seed_size, 0.0, 0.0),
                                             isDependent=False).id
        Datum3 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, 0.0),
                                             point2=(0.0, ly + by + 4*seed_size, 0.0),
                                             isDependent=False).id
        Datum4 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, 0.0),
                                             point2=(0.0, ly - by - 4*seed_size, 0.0),
                                             isDependent=False).id
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum1], faces=myPart.faces)
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum2], faces=myPart.faces)
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum3], faces=myPart.faces)
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum4], faces=myPart.faces)

        pickedFaces = myPart.faces.getByBoundingBox(
            xMin=lx/2 - bx/2 - rc - 2*seed_size - 0.1,
            yMin=ly/2 - by/2 - rc - 2*seed_size - 0.1, zMin=-0.1,
            xMax=lx/2 + bx/2 + rc + 2*seed_size + 0.1,
            yMax=ly/2 + by/2 + rc + 2*seed_size + 0.1, zMax=0.1)
        myPart.setMeshControls(regions=pickedFaces, technique=SWEEP)

    myPart.seedPart(size=seed_size, deviationFactor=0.05, minSizeFactor=0.1)
    myPart.generateMesh()

    myJob = mdb.Job(name=job_name, model=model_name, description='',
                    type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0,
                    queue=None, memory=4, memoryUnits=GIGA_BYTES,
                    getMemoryFromAnalysis=True, explicitPrecision=SINGLE,
                    nodalOutputPrecision=SINGLE, echoPrint=OFF,
                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
                    userSubroutine='', scratch='', resultsFormat=ODB,
                    parallelizationMethodExplicit=DOMAIN, numDomains=4,
                    activateLoadBalancing=False, multiprocessingMode=MPI,
                    numCpus=4, numGPUs=0)
    myJob.writeInput(consistencyChecking=OFF)
    # Viewing model
    session.viewports['Viewport: 1'].setValues(displayedObject=myPart)
    session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
    session.viewports['Viewport: 1'].view.rotate(xAngle=10, yAngle=-40, zAngle=0,
                                                 mode=MODEL)
    session.viewports['Viewport: 1'].view.fitView()
    time.sleep(1)
    myJob.submit()
    myJob.waitForCompletion()
    # --- get eigenvalue from ODB and return it ---
    myOdb = session.openOdb(name=job_name + '.odb')
    buck_frames = myOdb.steps['Step-1'].frames
    eig_text = buck_frames[1].description
    # Viewing results
    session.viewports['Viewport: 1'].setValues(displayedObject=myOdb)
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].setValues(
        displayedObject=myOdb)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
    session.viewports['Viewport: 1'].view.rotate(xAngle=10, yAngle=-40, zAngle=0,
                                                 mode=MODEL)
    session.viewports['Viewport: 1'].view.fitView()
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
        CONTOURS_ON_DEF,))
    # parse the number at the end
    eig_val = abs(float(eig_text.split('=')[-1]))
    return eig_val


# =====================================================
# 4) Diamond hole builder (ri = stiffener fillet, rc = hole fillet)
# =====================================================
def build_plate_dia(lx, ly, thickness, diag, lip_depth, ri, rc, seed_size):
    model_name = 'Model-1'
    job_name   = 'Job-1'
    myModel = mdb.Model(name=model_name)
    myPart  = myModel.Part(name='ShellPart', dimensionality=THREE_D, type=DEFORMABLE_BODY)

    if diag > 0.0:
        bx = diag / sqrt(2.0)
        by = diag / sqrt(2.0)
    else:
        bx = 0.0
        by = 0.0

    mySketch = myModel.ConstrainedSketch(name='ShellSketch', sheetSize=2000.0)
    mySketch.rectangle(point1=(0, 0), point2=(lx, ly))
    if diag > 0:
        mySketch.Line(point1=(lx/2 + sqrt(2)*bx/2, ly/2),
                      point2=(lx/2, ly/2 + sqrt(2)*by/2))
        mySketch.Line(point1=(lx/2, ly/2 + sqrt(2)*by/2),
                      point2=(lx/2 - sqrt(2)*bx/2, ly/2))
        mySketch.Line(point1=(lx/2 - sqrt(2)*bx/2, ly/2),
                      point2=(lx/2, ly/2 - sqrt(2)*by/2))
        mySketch.Line(point1=(lx/2, ly/2 - sqrt(2)*by/2),
                      point2=(lx/2 + sqrt(2)*bx/2, ly/2))
    myPart.BaseShell(sketch=mySketch)
    mySketch.unsetPrimaryObject()

    # Precompute arc lengths same way as your original
    arc_x = 0.0
    arc_y = 0.0
    if diag > 0:
        arc_x = sqrt(2)*bx/2 - (sqrt((bx - 2*5)**2 + (by - 2*5)**2) + 2*5)/2
        arc_y = sqrt(2)*by/2 - (sqrt((bx - 2*5)**2 + (by - 2*5)**2) + 2*5)/2

    # Lip extrude
    if diag > 0 and lip_depth > 0:
        right_side = myPart.edges.findAt(((lx, ly/2.0, 0),))
        lipTransform = myPart.MakeSketchTransform(sketchPlane=myPart.faces[0],
                                                  sketchUpEdge=right_side[0],
                                                  sketchPlaneSide=SIDE1,
                                                  sketchOrientation=RIGHT,
                                                  origin=(lx/2.0, ly/2.0, 0.0))
        lipSketch = myModel.ConstrainedSketch(name='CutSketch', sheetSize=200.0, transform=lipTransform)
        lipSketch.Line(point1=(sqrt(2)*bx/2, 0), point2=(0,  sqrt(2)*by/2))
        lipSketch.Line(point1=(0,  sqrt(2)*by/2), point2=(-sqrt(2)*bx/2, 0))
        lipSketch.Line(point1=(-sqrt(2)*bx/2, 0), point2=(0, -sqrt(2)*by/2))
        lipSketch.Line(point1=(0, -sqrt(2)*by/2), point2=(sqrt(2)*bx/2, 0))
        edge = myPart.edges.findAt(((lx, ly/2.0, 0),))
        myPart.ShellExtrude(sketchPlane=myPart.faces[0], sketchUpEdge=edge[0],
                            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT,
                            sketch=lipSketch, depth=lip_depth, flipExtrudeDirection=OFF)
        lipSketch.unsetPrimaryObject()

        # Hole corners fillet rc
        if rc > 0 and lip_depth > 0:
            stiffener_corners = myPart.edges.findAt(
                ((lx/2 + sqrt(2)*bx/2, ly/2, lip_depth - 0.1),),
                ((lx/2, ly/2 + sqrt(2)*by/2, lip_depth - 0.1),),
                ((lx/2 - sqrt(2)*bx/2, ly/2, lip_depth - 0.1),),
                ((lx/2, ly/2 - sqrt(2)*by/2, lip_depth - 0.1),),
            )
            myPart.Round(radius=rc, edgeList=(stiffener_corners[0], stiffener_corners[1],
                                              stiffener_corners[2], stiffener_corners[3],))

    # Stiffener fillet ri
    if ri > 0 and lip_depth > 0:
        hole_edges = myPart.edges.findAt(
            ((lx/2 + sqrt(2)*bx/4, ly/2 + sqrt(2)*by/4, 0),),
            ((lx/2 - sqrt(2)*bx/4, ly/2 + sqrt(2)*by/4, 0),),
            ((lx/2 - sqrt(2)*bx/4, ly/2 - sqrt(2)*by/4, 0),),
            ((lx/2 + sqrt(2)*bx/4, ly/2 - sqrt(2)*by/4, 0),),
            ((lx/2 + sqrt(2)*bx/2 - arc_x, ly/2, 0),),
            ((lx/2, ly/2 + sqrt(2)*by/2 - arc_y, 0),),
            ((lx/2 - sqrt(2)*bx/2 + arc_x, ly/2, 0),),
            ((lx/2, ly/2 - sqrt(2)*by/2 + arc_y, 0),),
        )
        myPart.Set(edges=hole_edges, name='Hole-Edges')
        myPart.Round(radius=ri, edgeList=(hole_edges[0], hole_edges[1], hole_edges[2], hole_edges[3],
                                              hole_edges[4], hole_edges[5], hole_edges[6], hole_edges[7],))

    # Material & section
    myModel.Material(name='Material-1')
    myModel.materials['Material-1'].Elastic(table=((200000, 0.3),))
    myModel.HomogeneousShellSection(name='Section-1', preIntegrate=OFF, thicknessType=UNIFORM,
                                    thickness=thickness, material='Material-1', thicknessField='',
                                    idealization=NO_IDEALIZATION, poissonDefinition=DEFAULT,
                                    temperature=GRADIENT, useDensity=OFF,
                                    integrationRule=SIMPSON, numIntPts=5)
    myPart.Set(faces=myPart.faces, name='Set-Shell')
    myPart.SectionAssignment(region=myPart.sets['Set-Shell'], sectionName='Section-1')

    # Assembly
    myAsm = myModel.rootAssembly
    myInst = myAsm.Instance(name='ShellInstance', part=myPart, dependent=ON)

    # Step
    myModel.BuckleStep(name='Step-1', previous='Initial', numEigen=5, vectors=30, maxIterations=300)

    # BCs
    myAsm.Set(name='ShellEdges', edges=myInst.edges.findAt(
        ((0.0, ly/2, 0.0),), ((lx, ly/2, 0.0),),
        ((lx/2, 0.0, 0.0),), ((lx/2, ly, 0.0),),))
    myAsm.Set(name='Corners', vertices=myInst.vertices.findAt(
        ((0.0, 0.0, 0.0),), ((0.0, ly, 0.0),),))
    myModel.DisplacementBC(name='BC-1', createStepName='Step-1',
                           region=myAsm.sets['ShellEdges'], u3=0.0)
    myModel.DisplacementBC(name='BC-2', createStepName='Step-1',
                           region=myAsm.sets['Corners'], u1=0.0, u2=0.0)

    # Loads
    myAsm.Surface(name='HorizontalEdges', side1Edges=myInst.edges.findAt(
        ((lx/2, 0.0, 0.0),), ((lx/2, ly, 0.0),),))
    myAsm.Surface(name='VerticalEdges', side1Edges=myInst.edges.findAt(
        ((0.0, ly/2, 0.0),), ((lx, ly/2, 0.0),),))
    myModel.ShellEdgeLoad(name='Load-1', createStepName='Step-1',
                          region=myAsm.surfaces['HorizontalEdges'], magnitude=1.0,
                          distributionType=UNIFORM, field='', localCsys=None, traction=SHEAR)
    myModel.ShellEdgeLoad(name='Load-2', createStepName='Step-1',
                          region=myAsm.surfaces['VerticalEdges'], magnitude=-1.0,
                          distributionType=UNIFORM, field='', localCsys=None, traction=SHEAR)

    # Partition (lighter version of your original)
    if bx > 0:
        Datum1 = myPart.DatumPlaneByTwoPoint(point1=(lx/2, ly/2, 0.0),
                                             point2=(lx/2 + sqrt(2)*bx + 2*rc + 4*seed_size, ly/2, 0.0),
                                             isDependent=False).id
        Datum2 = myPart.DatumPlaneByTwoPoint(point1=(lx/2, ly/2, 0.0),
                                             point2=(lx/2, ly/2 + sqrt(2)*by + 2*rc + 4*seed_size, 0.0),
                                             isDependent=False).id
        Datum3 = myPart.DatumPlaneByTwoPoint(point1=(lx/2, ly/2, 0.0),
                                             point2=(lx/2 - sqrt(2)*bx - 2*rc - 4*seed_size, ly/2, 0.0),
                                             isDependent=False).id
        Datum4 = myPart.DatumPlaneByTwoPoint(point1=(lx/2, ly/2, 0.0),
                                             point2=(lx/2, ly/2 - sqrt(2)*by - 2*rc - 4*seed_size, 0.0),
                                             isDependent=False).id
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum1], faces=myPart.faces)
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum2], faces=myPart.faces)
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum3], faces=myPart.faces)
        myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum4], faces=myPart.faces)

        if lip_depth > ri:
            Datum9 = myPart.DatumPlaneByTwoPoint(point1=(0.0, 0.0, ri),
                                                 point2=(0.0, 0.0, lip_depth),
                                                 isDependent=False).id
            myPart.PartitionFaceByDatumPlane(datumPlane=myPart.datums[Datum9], faces=myPart.faces)

        pickedFaces = myPart.faces.getByBoundingBox(
            xMin=lx/2 - sqrt(2)*bx/2 - rc - 2*seed_size - 0.1,
            yMin=ly/2 - sqrt(2)*by/2 - rc - 2*seed_size - 0.1, zMin=-0.1,
            xMax=lx/2 + sqrt(2)*bx/2 + rc + 2*seed_size + 0.1,
            yMax=ly/2 + sqrt(2)*by/2 + rc + 2*seed_size + 0.1, zMax=0.1)
        myPart.setMeshControls(regions=pickedFaces, technique=SWEEP)

    myPart.seedPart(size=seed_size, deviationFactor=0.05, minSizeFactor=0.1)
    myPart.generateMesh()

    myJob = mdb.Job(name=job_name, model=model_name, description='',
                    type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0,
                    queue=None, memory=4, memoryUnits=GIGA_BYTES,
                    getMemoryFromAnalysis=True, explicitPrecision=SINGLE,
                    nodalOutputPrecision=SINGLE, echoPrint=OFF,
                    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
                    userSubroutine='', scratch='', resultsFormat=ODB,
                    parallelizationMethodExplicit=DOMAIN, numDomains=4,
                    activateLoadBalancing=False, multiprocessingMode=MPI,
                    numCpus=4, numGPUs=0)
    myJob.writeInput(consistencyChecking=OFF)
    # Viewing model
    session.viewports['Viewport: 1'].setValues(displayedObject=myPart)
    session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
    session.viewports['Viewport: 1'].view.rotate(xAngle=10, yAngle=-40, zAngle=0,
                                                 mode=MODEL)
    session.viewports['Viewport: 1'].view.fitView()
    time.sleep(1)
    myJob.submit()
    myJob.waitForCompletion()
    # --- get eigenvalue from ODB and return it ---
    myOdb = session.openOdb(name=job_name + '.odb')
    buck_frames = myOdb.steps['Step-1'].frames
    eig_text = buck_frames[1].description
    # Viewing results
    session.viewports['Viewport: 1'].setValues(displayedObject=myOdb)
    a = mdb.models['Model-1'].rootAssembly
    session.viewports['Viewport: 1'].setValues(displayedObject=a)
    session.viewports['Viewport: 1'].setValues(
        displayedObject=myOdb)
    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
        optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
    session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
    session.viewports['Viewport: 1'].view.rotate(xAngle=10, yAngle=-40, zAngle=0,
                                                 mode=MODEL)
    session.viewports['Viewport: 1'].view.fitView()
    session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
        CONTOURS_ON_DEF,))
    # parse the number at the end
    eig_val = abs(float(eig_text.split('=')[-1]))
    return eig_val


# =====================================================
# 5) Dispatch
# =====================================================
# material props should match your model definition above
E_MAT  = 200000.0
NU_MAT = 0.3

if shape_str in ('CIR', 'CIRC', 'CIRCULAR'):
    eig = build_plate_cir(Lx, Ly, Thickness, HoleSize, LipDepth, Ri, SeedSize)
    pred = calculate_shear_buckling(
        a=Lx, d1=Ly, dh=HoleSize, s=LipDepth, t=Thickness,
        E=E_MAT, nu=NU_MAT, shape="Cir"
    )
elif shape_str in ('SQU', 'SQ', 'SQUARE'):
    eig = build_plate_squ(Lx, Ly, Thickness, HoleSize, LipDepth, Ri, Rc, SeedSize)
    pred = calculate_shear_buckling(
        a=Lx, d1=Ly, dh=HoleSize, s=LipDepth, t=Thickness,
        E=E_MAT, nu=NU_MAT, shape="Squ"
    )
elif shape_str in ('DIA', 'DIAMOND'):
    eig = build_plate_dia(Lx, Ly, Thickness, HoleSize, LipDepth, Ri, Rc, SeedSize)
    pred = calculate_shear_buckling(
        a=Lx, d1=Ly, dh=HoleSize, s=LipDepth, t=Thickness,
        E=E_MAT, nu=NU_MAT, shape="Dia"
    )
else:
    sys.exit('Invalid hole shape. Use: Cir, Squ, or Dia.')

# FEM shear buckling from your note:
tau_cr_fem = eig / Thickness

print '--- FEM vs Prediction ---'
print 'Eigenvalue (Abaqus):', eig
print 'tau_cr_FEM (MPa):   ', tau_cr_fem
print 'tau_cr_pred (MPa):  ', pred["tau_cr_stiffened (MPa)"]
print 'k_v_pred (stiff.):  ', pred["kv_stiffened"]
print 'teq_pred (mm):      ', pred["teq (mm)"]

#!/usr/bin/env python
"""
Generate a capillary network using Voronoi cells
"""

import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import igraph
import matplotlib.pyplot as plt


nx = 7
ny = 7
xmax = 200
ymax = 200

# Randomly generate a list of points
def generatePoints():
    narrowingFactor = 0.10 # narrows the cubes for random point generation. Value between 0 and 0.5.
    # print "Random seed: ", np.random.seed()
    points = []
    for i in range(nx):
        for j in range(ny):
            xl = xmax*(i   + narrowingFactor)/nx 
            xr = xmax*(i+1 - narrowingFactor)/nx 
            yb = ymax*(j   + narrowingFactor)/ny 
            yt = ymax*(j+1 - narrowingFactor)/ny 
            points.append([np.random.uniform(xl, xr), 
                           np.random.uniform(yb, yt)])

    return points

# Compute the Voronoi vertices
def computeVoronoiVertices():
    vor = Voronoi(points)
    return vor

# Construct a graph using iGraph
def buildGraph(vertices, edges):
    G = igraph.Graph(len(vertices))
    r = []
    for v in vertices:
        r.append(np.asarray(v))
    G.vs['r'] = r
    G.vs['type'] = G.vcount()*[None]
    # filter out the edges with index equal to -1
    edges = [e for e in edges if e[0] != -1 and e[1] != -1]
    G.add_edges(edges)
    G.es['type'] = G.ecount()*['c']
    return G

# Remove vertices outside a given bounding box
def removeVerticesOutsideBoundingBox(G, bMin, bMax):
    vertices = np.asarray(G.vs['r'])
    idxToRemove = []
    for vI in range(G.vcount()):
        if    any([vertices[vI,j] < bMin[j] for j in range(len(bMin))]) \
           or any([vertices[vI,j] > bMax[j] for j in range(len(bMax))]):
            idxToRemove.append(vI)
    G.delete_vertices(idxToRemove)

# Generate the positions of the arteries and veins
def generateArteriesAndVeins():
    na = 3
    nv = 3

    ar = [np.array([0,    ymax*(i+0.5)/na]) for i in range(na)]
    vr = [np.array([xmax, ymax*(i+0.5)/nv]) for i in range(nv)]

    return ar, vr

def closestVertex(G, coords):
    distToVertex = [np.linalg.norm(coords - vc, ord=2) for vc in G.vs['r']]
    return distToVertex.index(min(distToVertex))

def allUnique(x):
    seen = set()
    return not any(i in seen or seen.add(i) for i in x)

# Connect the arteries and veins to the graph
def connectArteriesAndVeins(G, ar, vr):
    aIdx = [len(G.vs) + i           for i in range(len(ar))]
    vIdx = [len(G.vs) + len(ar) + j for j in range(len(vr))]
    closestVIdxToArteries = [closestVertex(G, ar[i]) for i in range(len(ar))]
    closestVIdxToVeins    = [closestVertex(G, vr[i]) for i in range(len(vr))]
    if not allUnique(closestVIdxToArteries):
        print 'Warning: two arteries connect to the same vertex'
    if not allUnique(closestVIdxToVeins):
        print 'Warning: two veins connect to the same vertex'
    
    # add vertices and their coordinates
    G.add_vertices(len(aIdx) + len(vIdx))
    r = G.vs['r']
    vType = G.vs['type']
    for i, aI in enumerate(aIdx):
        r[aI] = ar[i]
        vType[aI] = 'a'
    for i, vI in enumerate(vIdx):
        r[vI] = vr[i]
        vType[vI] = 'v'
    G.vs['r'] = r
    G.vs['type'] = vType

    # add connecting edges
    eType = G.es['type']
    eType.extend(len(aIdx)*'a')
    eType.extend(len(vIdx)*'v')
    G.add_edges([(aI, j) for aI, j in zip(aIdx, closestVIdxToArteries)])
    G.add_edges([(vI, j) for vI, j in zip(vIdx, closestVIdxToVeins)])
    G.es['type'] = eType

def removeDanglingVertices(G):
    vertices = np.asarray(G.vs['r'])
    idxToRemove = []
    for vI in range(G.vcount()):
        if G.degree(vI) == 1 and G.vs['type'][vI] != 'a' and G.vs['type'][vI] != 'v':
            idxToRemove.append(vI)
    G.delete_vertices(idxToRemove)

def setInitialDiameters(G):
    capDiam = 4.0
    avDiam  = 8.0
    diameters = [capDiam if eType == 'c' else avDiam for eType in G.es['type']]
    G.es['diameter'] = diameters

def setBoundaryConditions(G, HtIn, vIn):
    G['av'] = [i for i, vType in enumerate(G.vs['type']) if vType == 'a']
    G['vv'] = [i for i, vType in enumerate(G.vs['type']) if vType == 'v']
    for aI in G['av']:
        G.vs[aI]['httBC'] = HtIn
        adjacent_edgeI = G.incident(aI)[0]
        radius = 0.5*G.es[adjacent_edgeI]['diameter']
        G.vs[aI]['rBC'] = vIn*np.pi*radius**2

    for vI in G['vv']:
        G.vs[vI]['pBC'] = 0

points = generatePoints()
vor = computeVoronoiVertices()
G = buildGraph(vor.vertices, vor.ridge_vertices)
removeVerticesOutsideBoundingBox(G, [0, 0], [xmax, ymax])
[ar, vr] = generateArteriesAndVeins()
connectArteriesAndVeins(G, ar, vr)
removeDanglingVertices(G)
setInitialDiameters(G)
setBoundaryConditions(G, HtIn=0.3, vIn=2)

# run simulation and compute average pressure

## Plot
# voronoi_plot_2d(vor)
igraph.plot(G, layout=G.vs['r'])
plt.show()


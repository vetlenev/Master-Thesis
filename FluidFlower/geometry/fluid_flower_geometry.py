import os
import argparse
import textwrap
import typing
import gmsh


"""parser = argparse.ArgumentParser(description="Create a structured gmsh grid from the SPE11 .geo file")
parser.add_argument("-g", "--geo-file", required=True, help="The .geo file to create a structured grid for")
parser.add_argument("-nx", "--number-of-cells-x", required=True, help="Desired number of cells in x-direction")
parser.add_argument("-ny", "--number-of-cells-y", required=True, help="Desired number of cells in y-direction")
args = vars(parser.parse_args())
"""
gmsh.initialize()
gmsh.open("draft_spe11_with_facies_markers.geo")

bbox = gmsh.model.getBoundingBox(-1, -1)
min, max = tuple(bbox[:3]), tuple(bbox[3:])
#nx, ny = int(args["number_of_cells_x"]), int(args["number_of_cells_y"])
#dx, dy = (max[0] - min[0])/float(nx), (max[1] - min[1])/float(ny)
space_dimension = 2
gmsh_quadrangle_id = 3

print("Reading physical groups")
physical_groups = {
    gmsh.model.getPhysicalName(dim=d, tag=t): (d,t)
    for d, t in gmsh.model.getPhysicalGroups(dim=2)
}

print(physical_groups)
for name, (dim, tag) in physical_groups.items():
    entities = gmsh.model.getEntitiesForPhysicalGroup(dim=dim, tag=tag)
    print(entities)

print("Entities")
ent62 = gmsh.model.getAdjacencies(dim=0, tag=62)
print(ent62)

elementTypes = gmsh.model.mesh.getElementTypes()
print(elementTypes)
name, dim, order, numv, parv, _ = gmsh.model.mesh.getElementProperties(60)
print(numv)



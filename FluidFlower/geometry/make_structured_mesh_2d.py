import os
import argparse
import textwrap
import typing
import gmsh

parser = argparse.ArgumentParser(description="Create a structured gmsh grid from the SPE11 .geo file")
parser.add_argument("-g", "--geo-file", required=True, help="The .geo file to create a structured grid for")
parser.add_argument("-nx", "--number-of-cells-x", required=True, help="Desired number of cells in x-direction")
parser.add_argument("-ny", "--number-of-cells-y", required=True, help="Desired number of cells in y-direction")
args = vars(parser.parse_args())

gmsh.initialize()
gmsh.open(args["geo_file"])

bbox = gmsh.model.getBoundingBox(-1, -1)
min, max = tuple(bbox[:3]), tuple(bbox[3:])
nx, ny = int(args["number_of_cells_x"]), int(args["number_of_cells_y"])
dx, dy = (max[0] - min[0])/float(nx), (max[1] - min[1])/float(ny)
space_dimension = 2
gmsh_quadrangle_id = 3

print("Reading physical groups")
physical_groups = {
    gmsh.model.getPhysicalName(dim=d, tag=t): (d,t)
    for d, t in gmsh.model.getPhysicalGroups(dim=2)
}

print("Creating structured lattice of points")
structured_mesh_points = [
    (min[0] + i*dx, min[1] + j*dy)
    for j in range(ny + 1)
    for i in range(nx + 1)
]

def _get_cell_corner_indices(cell_idx_x: int, cell_idx_y: int) -> tuple:
    p0 = cell_idx_y*(nx + 1) + cell_idx_x
    return (
        p0,
        p0 + 1,
        p0 + 1 + nx + 1,
        p0 + nx + 1
    )


def _get_cell_center(cell_idx_x: int, cell_idx_y: int) -> tuple:
    result = tuple([0.0, 0.0])
    corner_indices = _get_cell_corner_indices(cell_idx_x, cell_idx_y)
    for pidx in corner_indices:
        result = tuple([result[i] + structured_mesh_points[pidx][i] for i in range(space_dimension)])
    if len(corner_indices) > 0:
        result = tuple([result[i]/len(corner_indices) for i in range(space_dimension)])
    return result


def _make_3d(coord: tuple) -> tuple:
    assert len(coord) <= 3
    return tuple([coord[i] if i < len(coord) else 0.0 for i in range(3)])


def _find_physical_group(coordinate: tuple) -> typing.Tuple[str, int]:
    coordinate = _make_3d(coordinate)
    for name, (dim, tag) in physical_groups.items():
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim=dim, tag=tag)
        for entity_tag in entities:
            if gmsh.model.isInside(dim=dim, tag=entity_tag, coord=coordinate, parametric=False) > 0:
                return name, entity_tag
    raise RuntimeError("Could not find physical group for cell")

print("Determining cell groups")
gmsh_cell_index = 1
mesh_file_cell_entries = ["" for _ in range(nx*ny)]
for cell_idx_y in range(ny):
    for cell_idx_x in range(nx):
        group_name, entity_tag = _find_physical_group(_get_cell_center(cell_idx_x, cell_idx_y))
        group_tag = physical_groups[group_name][1]
        mesh_file_cell_entries[gmsh_cell_index-1] = f"{gmsh_cell_index} {gmsh_quadrangle_id} 2 {group_tag} {entity_tag} "
        mesh_file_cell_entries[gmsh_cell_index-1] += " ".join(
            str(i + 1) for i in _get_cell_corner_indices(cell_idx_x, cell_idx_y)
        )
        gmsh_cell_index += 1

msh_file_name = os.path.splitext(args["geo_file"])[0] + "_structured.msh"
print(f"Writing mesh file '{msh_file_name}'")
with open(msh_file_name, "w") as msh_file:
    msh_file.write(textwrap.dedent("""
        $MeshFormat
        2.2 0 8
        $EndMeshFormat
        $PhysicalNames
        {}
    """.format(len(physical_groups)).lstrip("\n")))
    msh_file.write("{}".format(
        "\n".join(f'2 {tag} "{name}"'
        for name, (_, tag) in physical_groups.items())
    ))
    msh_file.write(textwrap.dedent(f"""
        $EndPhysicalNames
        $Nodes
        {len(structured_mesh_points)}
    """))
    for count, p in enumerate(structured_mesh_points):
        msh_file.write(f"{count+1} {' '.join(str(c) for c in _make_3d(p))}\n".format())
    msh_file.write(textwrap.dedent(f"""
        $EndNodes
        $Elements
        {len(mesh_file_cell_entries)}
    """.lstrip("\n")))
    msh_file.write("\n".join(mesh_file_cell_entries))
    msh_file.write("\n$EndElements")

gmsh.fltk.run()
gmsh.finalize()

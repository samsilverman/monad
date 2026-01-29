from __future__ import annotations
from typing import Tuple
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon

DIRECTORY = Path(__file__).parent.resolve()

CMAP = LinearSegmentedColormap.from_list(
    'custom',
    [
        (0.0, 'white'),
        (0.5, '#e6d0d1'),
        (1.0, '#9e5457')
    ]
)

def get_nodes(file: str) -> np.ndarray:
    with open(file=(DIRECTORY / file), mode='r', encoding='utf-8') as f:
        lines = iter(f.readlines())

    nodes = []

    for line in lines:
        line = line.strip()

        if line == "$Nodes":
            _, num_nodes, _, _ = map(int, next(lines).split())

            # Skip second header line and nodeTag lines
            for _ in range(num_nodes + 1):
                next(lines)

            for _ in range(num_nodes):
                x, y, _ = next(lines).split()
                nodes.append([float(x), float(y)])

    return np.array(nodes)


def get_elements(file: str) -> np.ndarray:
    with open(file=(DIRECTORY / file), mode='r', encoding='utf-8') as f:
        lines = iter(f.readlines())

    elements = []

    for line in lines:
        line = line.strip()

        if line == "$Elements":
            _, num_elements, _, _ = map(int, next(lines).split())

            # Skip second header line
            next(lines)

            for _ in range(num_elements):
                element = map(int, next(lines).split())
                elements.append(list(element))

    return np.array(elements, dtype=int)[:, 1:] - 1


def get_densities(file: str) -> np.ndarray:
    with open(file=(DIRECTORY / file), mode='r', encoding='utf-8') as f:
        lines = iter(f.readlines())

    densities = []

    for line in lines:
        line = line.strip()

        if line == "$ElementData":
            # Skip first seven header lines
            for _ in range(6):
                next(lines)

            num_elements = int(next(lines))

            for _ in range(num_elements):
                _, density = next(lines).split()
                densities.append(float(density))

    return np.array(densities)


def get_displacements(file: str) -> np.ndarray:
    with open(file=(DIRECTORY / file), mode='r', encoding='utf-8') as f:
        lines = iter(f.readlines())

    displacements = []

    for line in lines:
        line = line.strip()

        if line == "$NodeData":
            # Skip first seven header lines
            for _ in range(6):
                next(lines)

            num_nodes = int(next(lines))

            for _ in range(num_nodes):
                _, u, v, _ = next(lines).split()
                displacements.append([float(u), float(v)])

    return np.array(displacements)


def draw_mesh_single(ax: plt.Axes, nodes: np.ndarray, elements: np.ndarray, densities: np.ndarray, shift: Tuple[float, float] = (0.0, 0.0)) -> None:
    dx, dy = shift

    for i, element in enumerate(elements):
        element_nodes = nodes[element] + np.array([dx, dy])
        density = densities[i]

        poly = Polygon(xy=element_nodes, closed=True, facecolor=CMAP(density), edgecolor=CMAP(density), linewidth=0.0, antialiased=True)
        ax.add_patch(poly)


def draw_mesh_tiled(ax: plt.Axes, nodes: np.ndarray, elements: np.ndarray, densities: np.ndarray, displacements: np.ndarray = None) -> None:
    shifts = [-1, 0, 1]

    x_min, y_min = nodes.min(axis=0)
    x_max, y_max = nodes.max(axis=0)

    lx = x_max - x_min
    ly = y_max - y_min

    deformed_nodes = nodes.copy()

    if displacements is not None:
        deformed_nodes += displacements

    centroid = deformed_nodes.mean(axis=0)
    deformed_nodes = deformed_nodes - centroid

    # Translation vectors from undeformed nodes
    v1 = np.array([lx, 0])
    v2 = np.array([0, ly])

    if displacements is not None:
        # Question why do we call []_m? 
        left_mask = np.isclose(nodes[:, 0], x_min)
        right_mask = np.isclose(nodes[:, 0], x_max)
        bottom_mask = np.isclose(nodes[:, 1], y_min)
        top_mask = np.isclose(nodes[:, 1], y_max)

        # Translation contribution from deformation
        jump_x = np.mean(displacements[right_mask], axis=0) - np.mean(displacements[left_mask], axis=0)
        jump_y = np.mean(displacements[top_mask], axis=0) - np.mean(displacements[bottom_mask], axis=0)

        v1 += jump_x
        v2 += jump_y

    # Tiling
    for i in shifts:
        for j in shifts:
            shift = (i * v1) + (j * v2)

            scale = 0.5
            if (i == 0 and j == 0):
                scale = 1.0

            draw_mesh_single(ax=ax, nodes=deformed_nodes, elements=elements, densities=scale * densities, shift=shift)

    ax.set_aspect("equal")
    ax.set_axis_off()

    view_scale = 2
    ax.set_xlim(-view_scale * lx / 2, view_scale * lx / 2)
    ax.set_ylim(-view_scale * ly / 2, view_scale * ly / 2)


def main() -> None:
    _, axes = plt.subplots(nrows=2, ncols=4, constrained_layout=True, figsize=(6.4, 3.2))

    nodes = get_nodes('design1/density.msh')
    elements = get_elements('design1/density.msh')
    rho = get_densities('design1/density.msh')
    u11 = get_displacements('design1/u11.msh')
    u22 = get_displacements('design1/u22.msh')
    u12 = get_displacements('design1/u12.msh')

    draw_mesh_tiled(ax=axes[0, 0], nodes=nodes, elements=elements, densities=rho)
    draw_mesh_tiled(ax=axes[0, 1], nodes=nodes, elements=elements, densities=rho, displacements=0.5 * u11)
    draw_mesh_tiled(ax=axes[0, 2], nodes=nodes, elements=elements, densities=rho, displacements=0.5 * u22)
    draw_mesh_tiled(ax=axes[0, 3], nodes=nodes, elements=elements, densities=rho, displacements=0.2 * u12)

    nodes = get_nodes('design2/density.msh')
    elements = get_elements('design2/density.msh')
    rho = get_densities('design2/density.msh')
    u11 = get_displacements('design2/u11.msh')
    u22 = get_displacements('design2/u22.msh')
    u12 = get_displacements('design2/u12.msh')

    draw_mesh_tiled(ax=axes[1, 0], nodes=nodes, elements=elements, densities=rho)
    draw_mesh_tiled(ax=axes[1, 1], nodes=nodes, elements=elements, densities=rho, displacements=0.5 * u11)
    draw_mesh_tiled(ax=axes[1, 2], nodes=nodes, elements=elements, densities=rho, displacements=0.5 * u22)
    draw_mesh_tiled(ax=axes[1, 3], nodes=nodes, elements=elements, densities=rho, displacements=0.2 * u12)

    path = Path(__file__).parent.resolve() / 'images'/ 'teaser.svg'
    plt.savefig(path, format='svg')

    path = Path(__file__).parent.resolve() / 'images'/ 'teaser.png'
    plt.savefig(path, format='png', dpi=200) # (6.4, 3.2) * 200 = (1280, 640)

    plt.show()


if __name__ == '__main__':
    main()

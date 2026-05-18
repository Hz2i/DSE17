from py3dbp import Packer, Bin, Item
import plotly.graph_objects as go
import numpy as np
import plotly.io as pio
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
pio.renderers.default = "browser"

from Objects.Constants import Constants

def get_packet(subsystem_mass, subsystem_density):
    """
    Calculate dimensions and weight of packets per subsystem, and return a dictionary with the results.
    It is possible to end up with multiple packets per subsystem.
    Inputs:
        subsystem mass
        subsystem density
    Returns:
        dict: A dictionary containing the dimensions and weight of packets for each subsystem.
    """
    constants = Constants()
    packet_mass = subsystem_mass * constants.container_packet_mass_factor
    packet_volume = packet_mass / subsystem_density
    #packet_dimensions = (packet_volume ** (1/3), packet_volume ** (1/3), packet_volume ** (1/3)) # Assuming cubic packets for simplicity
    packet_dimensions = (5.14, packet_volume / (5.14 * 1.46), 1.46) # Assuming a fixed width and height, and calculating the depth based on the volume
    packet = {'mass': packet_mass, 'volume': packet_volume, 'dimensions': packet_dimensions}
    print(packet)
    return packet

def assess_packing_feasibility(packets, visualize=False):
    """Assess the packing feasibility of the container, by comparing the used volume and mass with the container's capacity.
    Additionally, assess the packing feasibility of the container by comparing the dimensions of the packets with the container's door dimensions.
    """
    constants = Constants()
    packer = Packer()

    # Container
    packer.add_bin(Bin('ISO-20ft-container',
                constants.container_inner_width,    # width, m
                constants.container_inner_height,   # height, m
                constants.container_inner_length,   # depth, m
                constants.container_mass_capacity)) # max load, kg
    
    # Add packets as items to be packed
    for subsystem, packet in packets.items():
        dimensions = sorted(packet['dimensions'], reverse=True) # Works for tuples and lists
        packer.add_item(Item(subsystem,
                    dimensions[0],  # width, m
                    dimensions[1],  # height, m
                    dimensions[2],  # depth, m
                    packet['mass']))          # weight, kg
    # Solve the packing problem
    packer.pack()

    # Results
    for b in packer.bins:
        print("Container:", b.name)

        print("\nPacked items:")
        for item in b.items:
            print(item.name, item.position, item.width, item.height, item.depth)

        print("\nUnfitted items:")
        for item in b.unfitted_items:
            print(item.name)

    # Visualization (optional)
    if visualize == True:
        # Container
        vertices = np.array([
            [0, 0, 0],
            [constants.container_inner_width, 0, 0],
            [constants.container_inner_width, constants.container_inner_height, 0],
            [0, constants.container_inner_height, 0],
            [0, 0, constants.container_inner_length],
            [constants.container_inner_width, 0, constants.container_inner_length],
            [constants.container_inner_width, constants.container_inner_height, constants.container_inner_length],
            [0, constants.container_inner_height, constants.container_inner_length]
        ])

        # Container faces
        faces = [
            [vertices[0], vertices[1], vertices[2], vertices[3]],
            [vertices[4], vertices[5], vertices[6], vertices[7]],
            [vertices[0], vertices[1], vertices[5], vertices[4]],
            [vertices[2], vertices[3], vertices[7], vertices[6]],
            [vertices[1], vertices[2], vertices[6], vertices[5]],
            [vertices[4], vertices[7], vertices[3], vertices[0]]
        ]

        def plot_box(ax, x, y, z, w, h, d, color):
            vertices = np.array([
                [x, y, z],
                [x+w, y, z],
                [x+w, y+h, z],
                [x, y+h, z],
                [x, y, z+d],
                [x+w, y, z+d],
                [x+w, y+h, z+d],
                [x, y+h, z+d]
            ])

            faces = [
                [vertices[0], vertices[1], vertices[2], vertices[3]], # back face
                [vertices[4], vertices[5], vertices[6], vertices[7]], # front face
                [vertices[0], vertices[1], vertices[5], vertices[4]], # bottom face
                [vertices[2], vertices[3], vertices[7], vertices[6]], # top face
                [vertices[1], vertices[2], vertices[6], vertices[5]], # right face
                [vertices[4], vertices[7], vertices[3], vertices[0]]  # left face
            ]

            box = Poly3DCollection(faces, facecolors=color, edgecolors='black', linewidths=1, alpha=0.6)
            ax.add_collection3d(box)

        # Create figure
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot container
        cube = Poly3DCollection(faces, facecolors='cyan', edgecolors='black', linewidths=1, alpha=0.3)
        ax.add_collection3d(cube)

        # Plot items
        colors = ['red', 'green', 'blue', 'yellow', 'magenta', 'orange', 'purple', 'brown', 'pink', 'gray']
        for i, item in enumerate(b.items):
            # ensure numeric types are native floats (py3dbp may use Decimal)
            pos = item.position
            x, y, z = float(pos[0]), float(pos[1]), float(pos[2])
            w, h, d = float(item.width), float(item.height), float(item.depth)
            color = colors[i % len(colors)]
            plot_box(ax, x, y, z, w, h, d, color)

        # Set limits
        ax.set_xlim([0, constants.container_inner_width])
        ax.set_ylim([0, constants.container_inner_height])
        ax.set_zlim([0, constants.container_inner_length])

        # Labels
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

        # aspect ratio
        ax.set_box_aspect([constants.container_inner_width, constants.container_inner_height, constants.container_inner_length])

        plt.show()

def main():
    # Example usage
    packet_masses = [10, 20, 30, 40, 50] # kg
    packet_densities = [5, 6, 7, 8, 9] # kg/m^3
    packet_names = ["Packet 1", "Packet 2", "Packet 3", "Packet 4", "Packet 5"]
    packets = {}
    for packet_mass, packet_density, packet_name in zip(packet_masses, packet_densities, packet_names):
        packet = get_packet(packet_mass, packet_density)
        packets[packet_name] = packet
    assess_packing_feasibility(packets, visualize=True)

if __name__ == "__main__":
    main()

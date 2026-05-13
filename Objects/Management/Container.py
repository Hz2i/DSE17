from py3dbp import Packer, Bin, Item
import plotly.graph_objects as go
import plotly.io as pio
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
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
    packet_dimensions = (packet_volume ** (1/3), packet_volume ** (1/3), packet_volume ** (1/3)) # Assuming cubic packets for simplicity
    packet = {'mass': packet_mass, 'volume': packet_volume, 'dimensions': packet_dimensions}
    print(packet)
    return packet

def assess_packing_feasibility(packets, visualize=False):
    """Assess the packing feasibility of the container, by comparing the used volume and mass with the container's capacity.
    Additionally, assess the packing feasibility of the container by comparing the dimensions of the packets with the container's door dimensions.
    Returns:
        bool: True if the packing is feasible, False otherwise.
    """
    """container_volume_capacity = Constants.container_volume_capacity
    container_door_dimensions = Constants.container_door_dimensions
    container_mass_capacity = Constants.container_mass_capacity
    container_packing_efficiency = Constants.container_packing_efficiency
    feasibility_volume = False
    feasibility_mass = False
    feasibility_dimensions = False

    # Calculate volume and mass used
    container_volume_used = 0.0
    container_mass_used = 0.0
    for packet in packets.items():
        container_volume_used += packet['volume']
        container_mass_used += packet['mass']
    
    # Check volume feasibility
    if container_volume_used <= container_volume_capacity * container_packing_efficiency:
        feasibility_volume = True
    # Check mass feasibility
    if container_mass_used <= container_mass_capacity:
        feasibility_mass = True
    
    # Check if all packets fit through the door
    packets = sorted(packets.items(), key=lambda x: min(x[1]['dimensions']))
    for packet in packets:
        if min(packet[1]['dimensions']) > min(container_door_dimensions):
            feasibility_dimensions = False"""
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
        """fig = go.Figure()

        W = constants.container_inner_width
        H = constants.container_inner_height
        D = constants.container_inner_length

        corners = [(0, 0, 0), (W, 0, 0), (W, H, 0), (0, H, 0), (0, 0, D), (W, 0, D), (W, H, D), (0, H, D)]
        edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]

        for edge in edges:
            x = [corners[edge[0]][0], corners[edge[1]][0]]
            y = [corners[edge[0]][1], corners[edge[1]][1]]
            z = [corners[edge[0]][2], corners[edge[1]][2]]
            fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='black', width=2), showlegend=False))
        
        for item in b.items:
            x, y, z = item.position
            w, h, d = item.width, item.height, item.depth
            vertices = [(x, y, z), (x+w, y, z), (x+w, y+h, z), (x, y+h, z), (x, y, z+d), (x+w, y, z+d), (x+w, y+h, z+d), (x, y+h, z+d)]
            X = [v[0] for v in vertices]
            Y = [v[1] for v in vertices]
            Z = [v[2] for v in vertices]

            # Cube faces
            I = [0, 0, 0, 1, 4, 4, 5, 2, 3, 6, 7, 1]
            J = [1, 2, 3, 2, 5, 6, 6, 3, 7, 7, 4, 5]
            K = [2, 3, 1, 6, 6, 7, 2, 7, 4, 4, 5, 0]
        
            fig.add_trace(go.Mesh3d(x=X, y=Y, z=Z, i=I, j=J, k=K, opacity=0.5, name=item.name, 
                                    hovertext=f"{item.name}<br>Mass: {item.weight:.2f} kg<br>Dimensions: {item.width:.2f} x {item.height:.2f} x {item.depth:.2f} m", hoverinfo='text'))
            fig.update_layout(
                scene=dict(
                    xaxis_title='Width (m)',
                    yaxis_title='Height (m)',
                    zaxis_title='Depth (m)',
                    aspectmode='data'
                ),
                title='Container Packing Visualization'
            )
        fig.show()"""
        def plot_cube(ax, x, y, z, w, h, d, color):
            vertices = [(x, y, z), (x+w, y, z), (x+w, y+h, z), (x, y+h, z), (x, y, z+d), (x+w, y, z+d), (x+w, y+h, z+d), (x, y+h, z+d)]
            faces = [
                [vertices[0], vertices[1], vertices[2], vertices[3]],  # back
                [vertices[4], vertices[5], vertices[6], vertices[7]],  # front
                [vertices[0], vertices[1], vertices[5], vertices[4]],  # bottom
                [vertices[2], vertices[3], vertices[7], vertices[6]],  # top
                [vertices[1], vertices[2], vertices[6], vertices[5]],  # right
                [vertices[0], vertices[3], vertices[7], vertices[4]]   # left
            ]
            cube = Poly3DCollection(faces, alpha=0.5, facecolors=color)
            ax.add_collection3d(cube)
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot container
        ax.plot([0,constants.container_inner_width],[0,0],[0,0], color='black')
        ax.plot([0,0],[0,constants.container_inner_height],[0,0], color='black')
        ax.plot([0,0],[0,0],[0,constants.container_inner_length], color='black')
        ax.plot([constants.container_inner_width,constants.container_inner_width],[0,constants.container_inner_height],[0,0], color='black')
        ax.plot([0,constants.container_inner_width],[constants.container_inner_height,constants.container_inner_height],[0,0], color='black')
        ax.plot([0,constants.container_inner_width],[0,0],[constants.container_inner_length,constants.container_inner_length], color='black')

        # Plot items
        colors = ['red', 'blue', 'green', 'yellow', 'cyan', 'magenta', 'orange', 'purple']
        for idx, item in enumerate(b.items):
            x, y, z = item.position
            w, h, d = item.width, item.height, item.depth
            color = colors[idx % len(colors)]
            plot_cube(ax, x, y, z, w, h, d, color)
            ax.text(x + w/2, y + h/2, z + d/2, item.name, color='black', ha='center', va='center')

        # Labels
        ax.set_xlabel('Width (m)')
        ax.set_ylabel('Height (m)')
        ax.set_zlabel('Depth (m)')
        
        # aspect ratio
        ax.set_box_aspect((constants.container_inner_width, constants.container_inner_height, constants.container_inner_length))

        plt.title('Container Packing Visualization')
        plt.tight_layout()
        plt.show()


def main():
    # Example usage
    subsystem_mass = 1000         # kg
    subsystem_density = 500  # kg/m^3
    packet = get_packet(subsystem_mass, subsystem_density)
    packets = {'Test system': packet}
    assess_packing_feasibility(packets, visualize=True)

    

if __name__ == "__main__":
    main()

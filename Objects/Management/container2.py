
from py3dbp import Packer, Bin, Item
from py3dbp.constants import RotationType
from py3dbp.auxiliary_methods import intersect
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
 
from Objects.Constants import Constants
 
 
# ── Monkey-patch: fix py3dbp's put_item so it actually tries all 6 rotations ──
def _put_item_fixed(self, item, pivot):
    fit = False
    valid_item_position = item.position
    item.position = pivot
 
    for i in range(len(RotationType.ALL)):
        item.rotation_type = i
        dimension = item.get_dimension()
 
        if (
            self.width  < pivot[0] + dimension[0] or
            self.height < pivot[1] + dimension[1] or
            self.depth  < pivot[2] + dimension[2]
        ):
            continue
 
        fit = True
        for current_item_in_bin in self.items:
            if intersect(current_item_in_bin, item):
                fit = False
                break
 
        if fit:
            if self.get_total_weight() + item.weight > self.max_weight:
                item.position = valid_item_position
                return False
            self.items.append(item)
            return True
 
    item.position = valid_item_position
    return False
 
Bin.put_item = _put_item_fixed
# ─────────────────────────────────────────────────────────────────────────────
 
 
# Default fill ratio used for any packet that doesn't specify its own
DEFAULT_FILL_RATIO = 0.9
 
CONTAINER_LENGTH = 5.867
 
 
def build_packets_from_list(packet_list):
    """
    Build the packets dict from a list of entries.
 
    Each entry can optionally include a 'fill_ratio' key (0 < fill_ratio <= 1)
    to control how much padding is added to that specific packet. If omitted,
    DEFAULT_FILL_RATIO is used.
    """
    packets = {}
    for entry in packet_list:
        name = entry['name']
        w, h, d = entry['width'], entry['height'], entry['depth']
        mass = entry.get('mass', 0.0)
        fill_ratio = entry.get('fill_ratio', DEFAULT_FILL_RATIO)
        volume = w * h * d
        packets[name] = {
            'mass': mass,
            'volume': volume,
            'dimensions': (w, h, d),
            'fill_ratio': fill_ratio
        }
        print(f"{name}: dimensions=({w}, {h}, {d}), volume={volume:.4f} m³, "
              f"mass={mass} kg, fill_ratio={fill_ratio}")
    return packets
 
 
def _padded_dimensions(name, packet):
    """
    Compute the padded (w, h, d) for a single packet, using its own
    fill_ratio. Keeps the special-case logic for 'Tail wing'.
    """
    w, h, d = packet['dimensions']
    fill_ratio = packet['fill_ratio']
    scale = fill_ratio ** (1 / 3)
 
    if name == 'wing sec 1':
        # Cap the length (width) at container length, redistribute remaining
        # volume factor to h and d
        padded_w = 5.3
        target_volume = (w * h * d) / fill_ratio
        remaining_scale = (target_volume / padded_w) / (h * d)  # factor for h*d
        padded_h = h * (remaining_scale ** 0.5)
        padded_d = d * (remaining_scale ** 0.5)
    else:
        padded_w, padded_h, padded_d = w / scale, h / scale, d / scale
 
    return padded_w, padded_h, padded_d
 
 
def assess_packing_feasibility(packets, visualize=False):
    constants = Constants()
    packer = Packer()
 
    packer.add_bin(Bin(
        'ISO-20ft-container',
        constants.container_inner_width,
        constants.container_inner_height,
        constants.container_inner_length,
        constants.container_mass_capacity
    ))
 
    for name, packet in packets.items():
        padded_w, padded_h, padded_d = _padded_dimensions(name, packet)
        packer.add_item(Item(name, padded_w, padded_h, padded_d, packet['mass']))
 
    packer.pack(bigger_first=True)
 
    for b in packer.bins:
        print(f"\nContainer: {b.name}")
        print(f"\nPacked items ({len(b.items)}):")
        for item in b.items:
            print(f"  {item.name:20s}  pos={item.position}  dim={item.get_dimension()}  rot={item.rotation_type}")
        print(f"\nUnfitted items ({len(b.unfitted_items)}):")
        for item in b.unfitted_items:
            print(f"  {item.name}")
 
    if visualize:
        cw = constants.container_inner_width
        ch = constants.container_inner_height
        cl = constants.container_inner_length
 
        def make_box_faces(x, y, z, w, h, d):
            v = np.array([
                [x,   y,   z],   [x+w, y,   z],
                [x+w, y+h, z],   [x,   y+h, z],
                [x,   y,   z+d], [x+w, y,   z+d],
                [x+w, y+h, z+d], [x,   y+h, z+d]
            ])
            return [
                [v[0], v[1], v[2], v[3]], [v[4], v[5], v[6], v[7]],
                [v[0], v[1], v[5], v[4]], [v[2], v[3], v[7], v[6]],
                [v[1], v[2], v[6], v[5]], [v[4], v[7], v[3], v[0]]
            ]
 
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
 
        ax.add_collection3d(Poly3DCollection(
            make_box_faces(0, 0, 0, cw, ch, cl),
            facecolors='cyan', edgecolors='black', linewidths=1, alpha=0.15
        ))
 
        colors = ['red', 'green', 'blue', 'yellow', 'magenta',
                  'orange', 'purple', 'brown', 'pink', 'gray']
        for i, item in enumerate(b.items):
            x, y, z = [float(v) for v in item.position]
            w, h, d = [float(v) for v in item.get_dimension()]
            ax.add_collection3d(Poly3DCollection(
                make_box_faces(x, y, z, w, h, d),
                facecolors=colors[i % len(colors)],
                edgecolors='black', linewidths=0.8, alpha=0.6
            ))
            ax.text(x + w/2, y + h/2, z + d/2, item.name, fontsize=7, ha='center')
 
        ax.set_xlim([0, cw]); ax.set_ylim([0, ch]); ax.set_zlim([0, cl])
        ax.set_xticks(np.arange(0, cw+1, 1)); ax.set_yticks(np.arange(0, ch+1, 1)); ax.set_zticks(np.arange(0, cl+1, 1))
        ax.set_xlabel("Width (m)"); ax.set_ylabel("Height (m)"); ax.set_zlabel("Length (m)")
        ax.set_title("ISO 20ft Container Packing")
        ax.set_box_aspect([cw, ch, cl])
        plt.tight_layout()
        plt.show()
 
 
def main():
 
    wing_list = [
        {'name': 'Wing sec 1', 'width': 5.20, 'height': 1.774, 'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'Wing sec 2', 'width': 4.38, 'height': 1.774, 'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'Wing sec 3', 'width': 4.38, 'height': 1.774, 'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'Wing sec 4', 'width': 4.38, 'height': 1.774, 'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'Wing sec 5', 'width': 4.38, 'height': 1.774, 'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'Wing sec 6', 'width': 4.38, 'height': 1.774, 'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'Wing sec 7', 'width': 4.38, 'height': 1.774, 'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'Wing sec 8', 'width': 4.72, 'height': 1.774, 'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'Wing sec 9', 'width': 4.72, 'height': 1.774, 'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'Rudder 1',    'width': 1.97, 'height': 1.67, 'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'Rudder 2',    'width': 1.97, 'height': 1.67, 'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'proppelor box 1',    'width': 1.7,  'height': 0.5,  'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'proppelor box2',    'width': 1.7,  'height': 0.5,  'depth': 0.267, 'fill_ratio': 0.8},
        {'name': 'Battery box 1',      'width': 0.251,  'height': 0.251,  'depth': 0.251, 'fill_ratio': 1},
        {'name': 'Battery box 2',      'width': 0.251,  'height': 0.251,  'depth': 0.251, 'fill_ratio': 1},
        {'name': 'Battery box 3',      'width': 0.251,  'height': 0.251,  'depth': 0.251, 'fill_ratio': 1},
        {'name': 'Battery box 4',      'width': 0.251,  'height': 0.251,  'depth': 0.251, 'fill_ratio': 1},
        {'name': 'Pitot tube 1',    'width': 0.194, 'height': 0.194, 'depth': 0.194, 'fill_ratio': 1},
        {'name': 'Pitot tube 2',    'width': 0.194, 'height': 0.194, 'depth': 0.194, 'fill_ratio': 1},
        {'name': 'strut 1',    'width': 0.8, 'height': 0.02, 'depth': 0.02, 'fill_ratio': 1},
        {'name': 'strut 2',    'width': 0.8, 'height': 0.02, 'depth': 0.02, 'fill_ratio': 1},
        {'name': 'strut 3',    'width': 0.8, 'height': 0.02, 'depth': 0.02, 'fill_ratio': 1},
        {'name': 'strut 4',    'width': 0.8, 'height': 0.02, 'depth': 0.02, 'fill_ratio': 1},
        {'name': 'strut 5',    'width': 0.8, 'height': 0.02, 'depth': 0.02, 'fill_ratio': 1},
        {'name': 'strut 6',    'width': 0.8, 'height': 0.02, 'depth': 0.02, 'fill_ratio': 1},
        {'name': 'strut 7',    'width': 0.8, 'height': 0.02, 'depth': 0.02, 'fill_ratio': 1},
        {'name': 'strut 8',    'width': 0.8, 'height': 0.02, 'depth': 0.02, 'fill_ratio': 1},
        {'name': 'wheel 1',    'width': 0.8, 'height': 0.8, 'depth': 0.038, 'fill_ratio': 1},
        {'name': 'wheel 2',    'width': 0.8, 'height': 0.8, 'depth': 0.038, 'fill_ratio': 1},
        {'name': 'wheel 3',    'width': 0.8, 'height': 0.8, 'depth': 0.038, 'fill_ratio': 1},
        {'name': 'wheel 4',    'width': 0.8, 'height': 0.8, 'depth': 0.038, 'fill_ratio': 1},
        {'name': 'wheel 5',    'width': 0.8, 'height': 0.8, 'depth': 0.038, 'fill_ratio': 1},
        {'name': 'wheel 6',    'width': 0.8, 'height': 0.8, 'depth': 0.038, 'fill_ratio': 1},
        {'name': 'wheel 7',    'width': 0.8, 'height': 0.8, 'depth': 0.038, 'fill_ratio': 1},
        {'name': 'wheel 8',    'width': 0.8, 'height': 0.8, 'depth': 0.038, 'fill_ratio': 1},
        {'name': 'Assembly jack',    'width': 0.8, 'height': 0.3, 'depth': 0.3, 'fill_ratio': 1},
        {'name': 'wrench 1', 'width': 0.2, 'height': 0.2, 'depth': 0.1, 'fill_ratio': 1},
        {'name': 'wrench 2', 'width': 0.2, 'height': 0.2, 'depth': 0.1, 'fill_ratio': 1},
        {'name': 'wrench 3', 'width': 0.2, 'height': 0.2, 'depth': 0.1, 'fill_ratio': 1},
        {'name': 'wrench 4', 'width': 0.2, 'height': 0.2, 'depth': 0.1, 'fill_ratio': 1},
        {'name': 'Heat gun 1', 'width': 0.2, 'height': 0.2, 'depth': 0.1, 'fill_ratio': 1},
        {'name': 'Heat gun 2', 'width': 0.2, 'height': 0.2, 'depth': 0.1, 'fill_ratio': 1},
        {'name': 'Ground station', 'width': 0.3, 'height': 0.15, 'depth': 0.05, 'fill_ratio': 0.8},
        {'name': 'Soldering kit 1', 'width': 0.15, 'height': 0.15, 'depth': 0.01, 'fill_ratio': 1},
        {'name': 'Soldering kit 2', 'width': 0.15, 'height': 0.15, 'depth': 0.01, 'fill_ratio': 1}
    ]
 
    packets = build_packets_from_list(wing_list)
    assess_packing_feasibility(packets, visualize=True)
 
 
main()
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


def build_packets_from_list(packet_list):
    packets = {}
    for entry in packet_list:
        name = entry['name']
        w, h, d = entry['width'], entry['height'], entry['depth']
        mass = entry.get('mass', 0.0)
        volume = w * h * d
        packets[name] = {
            'mass': mass,
            'volume': volume,
            'dimensions': (w, h, d)
        }
        print(f"{name}: dimensions=({w}, {h}, {d}), volume={volume:.4f} m³, mass={mass} kg")
    return packets


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

    FILL_RATIO = 0.7
    scale = FILL_RATIO ** (1/3)
    CONTAINER_LENGTH = 5.867 

    for name, packet in packets.items():
        w, h, d = packet['dimensions']

        if name == 'Tail wing':
            # Cap the length (width) at container length, redistribute remaining volume factor to h and d
            padded_w = CONTAINER_LENGTH
            target_volume = (w * h * d) / FILL_RATIO
            remaining_scale = (target_volume / (padded_w)) / (h * d)  # factor to apply to h*d
            padded_h = h * (remaining_scale ** 0.5)
            padded_d = d * (remaining_scale ** 0.5)
        else:
            padded_w, padded_h, padded_d = w / scale, h / scale, d / scale

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
    con_list = [
        {'name': 'Wing sec 1',     'width': 4.72, 'height': 1.32, 'depth': 0.2},
        {'name': 'Wing sec 2',     'width': 4.72, 'height': 1.32, 'depth': 0.2},
        {'name': 'Wing sec 3',     'width': 4.72, 'height': 1.32, 'depth': 0.2},
        {'name': 'Wing sec 4',     'width': 4.72, 'height': 1.32, 'depth': 0.2},
        {'name': 'Wing sec 5',     'width': 4.72, 'height': 1.32, 'depth': 0.2},
        {'name': 'Wing sec 6',     'width': 4.72, 'height': 1.32, 'depth': 0.2},
        {'name': 'Wing sec 7',     'width': 4.72, 'height': 1.32, 'depth': 0.2},
        {'name': 'Tail wing',      'width': 5.7,  'height': 1.14, 'depth': 0.14},
        {'name': 'rudder',         'width': 2.48, 'height': 1.24, 'depth': 0.15},
        {'name': 'fuselage sec 1', 'width': 5.0,  'height': 0.5,  'depth': 0.5},
        {'name': 'fuselage sec 2', 'width': 5.0,  'height': 0.5,  'depth': 0.5},
        {'name': 'proppelor 1',    'width': 0.3,  'height': 0.3,  'depth': 1.5},
        {'name': 'proppelor 2',    'width': 0.3,  'height': 0.3,  'depth': 1.5},
        {'name': 'proppelor 3',    'width': 0.3,  'height': 0.3,  'depth': 1.5},
        {'name': 'proppelor 4',    'width': 0.3,  'height': 0.3,  'depth': 1.5},
    ]

    wing_list = [
        {'name': 'Wing sec 1', 'width': 4.30, 'height': 1.20, 'depth': 0.2},
        {'name': 'Wing sec 2', 'width': 4.30, 'height': 1.20, 'depth': 0.2},
        {'name': 'Wing sec 3', 'width': 4.30, 'height': 1.20, 'depth': 0.2},
        {'name': 'Wing sec 4', 'width': 4.30, 'height': 1.20, 'depth': 0.2},
        {'name': 'Wing sec 5', 'width': 4.30, 'height': 1.20, 'depth': 0.2},
        {'name': 'Wing sec 6', 'width': 4.30, 'height': 1.20, 'depth': 0.2},
        {'name': 'Wing sec 7', 'width': 4.30, 'height': 1.20, 'depth': 0.2},
        {'name': 'proppelor 1',    'width': 0.3,  'height': 0.3,  'depth': 1.5},
        {'name': 'proppelor 2',    'width': 0.3,  'height': 0.3,  'depth': 1.5},
        {'name': 'proppelor 3',    'width': 0.3,  'height': 0.3,  'depth': 1.5},
        {'name': 'proppelor 4',    'width': 0.3,  'height': 0.3,  'depth': 1.5}
    ]


    packets = build_packets_from_list(con_list)
    assess_packing_feasibility(packets, visualize=True)


main()
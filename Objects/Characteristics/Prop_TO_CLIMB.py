import aerosandbox as asb
import aerosandbox.numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import brentq
from pathlib import Path
import sys

try:
	from Objects.Characteristics.PropulsionSystem import PropulsionSystem
except ModuleNotFoundError:
	# Allow running this file directly via: python Objects/Characteristics/Prop_TO_CLIMB.py
	repo_root = Path(__file__).resolve().parents[2]
	if str(repo_root) not in sys.path:
		sys.path.insert(0, str(repo_root))
	from Objects.Characteristics.PropulsionSystem import PropulsionSystem


# DESIGN INPUTS
D = 1.7045684187645866  # m (already optimized for cruise)
v_inf_cruise = 32.45494362484863  # m/s
required_thrust_cruise = 65.54910335953765  # N
m_TO = 278.1416060277118  # kg
S = 55.79535859750816  # m^2
CL_max = 1.0319550892283087  # -

# Take-off start speed from human push
v_initial = 4.0  # m/s



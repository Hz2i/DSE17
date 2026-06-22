import numpy as np
from math import comb

# Mission time (hours)
t = 672  # 28 days

# -------- MTBF VALUES (EDIT THESE) --------
MTBF = {
    "airframe": 10 * 733.43,
    "power": 10 * 1315.5,
    "engine": 10 * 1315.5,
    "control": 10 * 4134.9,
    "computer": 10 * 1315.5,
    "fcs": 10 * 1315.5,
    "comms": 10 * 1315.5,
    "payload": 10 * 134063.7,
    "undercarriage": 10 * 4134.9,
    "thermal": 10 * 4134.9
}
# MTBF = {
#     "airframe": 10 * 4134.9,
#     "power": 10 * 4134.9,
#     "engine": 10 * 134063.7,
#     "control": 10 * 4134.9,
#     "computer": 10 * 4134.9,
#     "fcs": 10 * 134063.7,
#     "comms": 10 * 134063.7,
#     "payload": 10 * 134063.7,
#     "undercarriage": 10 * 4134.9,
#     "thermal": 10 * 4134.9
# }
# MTBF = {
#     "airframe": 40000,
#     "power": 40000,
#     "engine": 40000,
#     "control": 40000,
#     "computer": 40000,
#     "fcs": 40000,
#     "comms": 40000,
#     "payload": 40000,
#     "undercarriage": 40000,
#     "thermal": 40000
# }

# -------- BASIC RELIABILITY FUNCTION --------
def reliability(mtbf, t):
    return np.exp(-t / mtbf)

# -------- PARALLEL (1-out-of-2) --------
def parallel_2(R1, R2):
    return 1 - (1 - R1) * (1 - R2)

# -------- PROPULSION (3-out-of-4) --------
def propulsion_reliability(R_engine):
    R = 0
    for i in range(3, 5):
        R += comb(4, i) * (R_engine**i) * ((1 - R_engine)**(4 - i))
    return R

# -------- COMPUTE SUBSYSTEM RELIABILITIES --------
R_airframe = reliability(MTBF["airframe"], t)

R_power = parallel_2(
    reliability(MTBF["power"], t),
    reliability(MTBF["power"], t)
)

R_engine = reliability(MTBF["engine"], t)
R_prop = propulsion_reliability(R_engine)

R_control = reliability(MTBF["control"], t)

# Redundant systems
R_computer = parallel_2(
    reliability(MTBF["computer"], t),
    reliability(MTBF["computer"], t)
)

R_fcs = parallel_2(
    reliability(MTBF["fcs"], t),
    reliability(MTBF["fcs"], t)
)

R_comms = parallel_2(
    reliability(MTBF["comms"], t),
    reliability(MTBF["comms"], t)
)

R_payload = reliability(MTBF["payload"], t)
R_undercarriage = reliability(MTBF["undercarriage"], t)
R_thermal = parallel_2(
    reliability(MTBF["thermal"], t),
    reliability(MTBF["thermal"], t)
)

# -------- TOTAL SYSTEM RELIABILITY --------
R_total = (
    R_airframe *
    R_power *
    R_computer *
    R_prop *
    R_control *
    R_fcs *
    R_undercarriage *
    R_thermal *
    R_payload *
    R_comms
)

R_aircraft = (
    R_airframe *
    R_power *
    R_computer *
    R_prop *
    R_control *
    R_fcs *
    R_undercarriage *
    R_thermal
)
# -------- OUTPUT --------
print("Subsystem reliabilities:")
print(f"Airframe: {R_airframe:.4f}")
print(f"Power: {R_power:.4f}")
print(f"Propulsion: {R_prop:.4f}")
print(f"Control: {R_control:.4f}")
print(f"Computer (redundant): {R_computer:.4f}")
print(f"FCS (redundant): {R_fcs:.4f}")
print(f"Comms (redundant): {R_comms:.4f}")
print(f"Payload: {R_payload:.4f}")
print(f"Undercarriage: {R_undercarriage:.4f}")
print(f"Thermal: {R_thermal:.4f}")

print("\nTotal system reliability:")
print(f"R_total = {R_total:.4f}")
print(f"R_aircraft (excluding payload and comms) = {R_aircraft:.4f}")
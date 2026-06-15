import matplotlib.pyplot as plt
import numpy as np

# -------- TIME SETTINGS --------
total_years = 7

# Convert durations to years
day_to_year = 1 / 365

durations = {
    "A": 1 * day_to_year,     # ~1 day
    "B": 3 * day_to_year,   # 2–3 days
    "C": 10 * day_to_year,    # ~1–2 weeks
    "D": 45 * day_to_year     # ~30–60 days
}

# -------- CHECK INTERVALS --------
A_times = np.arange(0, total_years, 0.324)   # every 500 flight hours (~0.057 yr)
B_times = np.arange(0, total_years, 0.6)          # ~6–8 months
C_times = np.arange(0, total_years, 1.8)            # every 2 years
D_times = np.arange(0, total_years, 7)            # every 7 years

# -------- PLOT --------
plt.figure(figsize=(12, 3))

# A checks
for t in A_times:
    plt.barh(0.0, durations["A"], left=t, color="black", alpha=1, label="A check" if t == 0 else "")

# B checks
for t in B_times:
    plt.barh(0.0, durations["B"], left=t, color="lime", alpha=1, label="B check" if t == 0 else "")

# C checks
for t in C_times:
    plt.barh(0.0, durations["C"], left=t, color="blue", alpha=1, label="C check" if t == 0 else "")

# D checks
for t in D_times:
    plt.barh(0.0, durations["D"], left=t, color="red", alpha=1, label="D check" if t == 0 else "")

# Labels
plt.xlabel("Time [years]")
plt.title("Maintenance Schedule with Duration")
plt.legend(loc="upper right")

plt.xlim(0, total_years)
plt.grid(True, axis='x', linestyle='--', alpha=0.5)

plt.tight_layout()
plt.show()

def merge_intervals(intervals):
    # Sort by start time
    intervals.sort(key=lambda x: x[0])
    
    merged = []
    for interval in intervals:
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
            merged[-1] = (
                merged[-1][0],
                max(merged[-1][1], interval[1])
            )
    return merged


def compute_downtime(intervals):
    merged = merge_intervals(intervals)
    return sum(end - start for start, end in merged)

intervals = []

# Example for B checks
for t in B_times:
    intervals.append((t, t + durations["B"]))

# Do this for ALL checks
for t in A_times:
    intervals.append((t, t + durations["A"]))

for t in C_times:
    intervals.append((t, t + durations["C"]))

for t in D_times:
    intervals.append((t, t + durations["D"]))

T_total = total_years
T_down = compute_downtime(intervals)

availability = 1 - T_down / T_total

print(f"Total downtime: {T_down:.3f} years")
print(f"Availability: {availability:.4f}")
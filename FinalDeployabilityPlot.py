import numpy as np
import matplotlib.pyplot as plt
import ambiance as am


# -----------------------------
# Daylight / Sunrise / Sunset
# -----------------------------
class LightData:
    def __init__(self, latitude_deg, days=None):
        self.latitude_deg = float(latitude_deg)
        self.latitude_rad = np.deg2rad(self.latitude_deg)
        self.days = np.arange(1, 366) if days is None else np.asarray(days, dtype=int)

        self.declination_deg = self.solar_declination(self.days)
        self.daylight_hours = self.compute_daylight_hours(self.days)
        self.sunrise_list, self.sunset_list = self.compute_sunrise_sunset_lists(self.days)

    @staticmethod
    def solar_declination(day_of_year):
        day_of_year = np.asarray(day_of_year, dtype=float)
        return 23.45 * np.sin(np.deg2rad((360.0 / 365.0) * (day_of_year + 10.0)))

    def compute_daylight_hours(self, day_of_year):
        day_of_year = np.asarray(day_of_year, dtype=float)
        delta_rad = np.deg2rad(self.solar_declination(day_of_year))
        phi = self.latitude_rad
        cos_h0 = -np.tan(phi) * np.tan(delta_rad)
        cos_h0 = np.clip(cos_h0, -1.0, 1.0)
        h0 = np.arccos(cos_h0)
        return (2.0 / 15.0) * np.rad2deg(h0)

    def compute_sunrise_sunset_lists(self, day_of_year):
        daylen = self.compute_daylight_hours(day_of_year)
        sunrise = 12.0 - 0.5 * daylen
        sunset = 12.0 + 0.5 * daylen
        return sunrise, sunset


# -----------------------------
# Solar power (daily average)
# -----------------------------
class SolarPower:
    def __init__(self, latitude_deg=0.0, solar_area_m2=1.0, days=None):
        self.latitude_deg = float(latitude_deg)
        self.latitude_rad = np.deg2rad(self.latitude_deg)

        self.days = np.arange(1, 366) if days is None else np.asarray(days, dtype=int)
        self.solar_area_m2 = float(solar_area_m2)

        self.efficiency = 0.2
        self.I0 = 1378.0
        self.powLimS = 250.0

        self.incidence_angle_rad = self.calc_daily_mean_incidence_angle(self.days)
        self.power_per_m2_W = self.calc_power_per_m2(self.days)
        self.power_W = self.power_per_m2_W * self.solar_area_m2

    @staticmethod
    def calc_solar_declination_rad(day_of_year):
        day_of_year = np.asarray(day_of_year, dtype=float)
        decl_deg = 23.45 * np.sin(np.deg2rad((360.0 / 365.0) * (day_of_year + 10.0)))
        return np.deg2rad(decl_deg)

    def calc_daily_mean_incidence_angle(self, day_of_year, n_samples=200):
        day_of_year = np.asarray(day_of_year, dtype=float)
        delta = self.calc_solar_declination_rad(day_of_year)
        phi = self.latitude_rad

        cos_h0 = -np.tan(phi) * np.tan(delta)
        cos_h0 = np.clip(cos_h0, -1.0, 1.0)
        h0 = np.arccos(cos_h0)

        incidence = np.zeros_like(day_of_year, dtype=float)
        for i in range(day_of_year.size):
            if np.isclose(h0[i], 0.0):
                incidence[i] = np.pi / 2
                continue

            hour_angles = np.linspace(-h0[i], h0[i], n_samples)
            cos_theta = (
                np.sin(delta[i]) * np.sin(phi)
                + np.cos(delta[i]) * np.cos(phi) * np.cos(hour_angles)
            )
            cos_theta = np.clip(cos_theta, 0.0, 1.0)
            incidence[i] = np.arccos(np.mean(cos_theta))
        return incidence

    def calc_power_per_m2(self, day_of_year):
        incidence = self.calc_daily_mean_incidence_angle(day_of_year)
        cos_inc = np.cos(incidence)
        raw = self.efficiency * self.I0 * np.maximum(cos_inc, 0.0)
        return np.minimum(self.powLimS, raw)


# -----------------------------
# Mission profile
# -----------------------------
class MissionProfile:
    def __init__(self, time_daylight_cruise, time_night, solar_area=30, MTOW=150.0, M_battery=58.6):
        self.m_battery_guess = M_battery
        self.gamma_guess = 2
        self.S_guess = 1.1 * solar_area

        self.E_battery_guess = self.m_battery_guess * 400 * 3600  # J

        # given
        self.g = 9.81
        self.alt = 60000 * 0.3048
        self.CD0 = 0.0109
        self.AR = 25.0
        self.e = 0.85
        self.Pavg_climb_subsys = 300
        self.Pavg_cruise_subsys = 425
        self.eta_prop = 0.7133
        self.LD = 37.7
        self.V_cruise = 25
        self.t_daylight_cruise = time_daylight_cruise
        self.t_night = time_night

        # derived
        self.m_total_guess = MTOW
        self.density_climb = self.Calc_Density_Climb()
        self.Cl_opt_climb = self.Calc_Cl_opt_climb()
        self.CD_total_climb = self.Calc_CD_total(self.Cl_opt_climb)
        self.gamma_rad = np.radians(self.gamma_guess)

        # climb
        self.V_climb = self.Calc_V_climb()
        self.t_climb = self.Calc_t_climb()
        self.D_climb = self.Calc_D_climb()
        self.E_climb = self.Calc_E_climb()

        # cruise
        self.Pprop_cruise = self.Calc_Pprop_cruise()
        self.Pavg_cruise = self.Pprop_cruise + self.Pavg_cruise_subsys
        self.E_cruise = self.Pavg_cruise * self.t_daylight_cruise

    def Calc_Density_Climb(self):
        h_range = np.linspace(0, self.alt, 100000)
        density = am.Atmosphere(h_range).density
        return np.trapz(density, h_range) / self.alt

    def Calc_Cl_opt_climb(self):
        return np.sqrt(self.CD0 * np.pi * self.AR * self.e)

    def Calc_CD_total(self, CL):
        k = 1.0 / (np.pi * self.AR * self.e)
        return self.CD0 + k * CL**2

    def Calc_V_climb(self):
        return np.sqrt(
            2 * self.m_total_guess * self.g * np.cos(self.gamma_rad)
            / (self.density_climb * self.S_guess * self.Cl_opt_climb)
        )

    def Calc_D_climb(self):
        return 0.5 * self.CD_total_climb * self.S_guess * self.density_climb * self.V_climb**2

    def Calc_t_climb(self):
        return self.alt / (self.V_climb * np.sin(self.gamma_rad))

    def Calc_E_climb(self):
        Epot = self.m_total_guess * self.g * self.alt
        Edrag = self.D_climb * (self.alt / np.sin(self.gamma_rad))
        Esubsys = self.Pavg_climb_subsys * self.t_climb
        return Epot + Edrag + Esubsys

    def Calc_Pprop_cruise(self):
        return ((self.m_total_guess * self.g) / self.LD) * self.V_cruise / self.eta_prop


# -----------------------------
# Grid-based bounds + wrap-aware plots + YEAR WRAP
# -----------------------------
class MultiDayDeployWindowBoundsFast:
    """
    This version fixes "start of year can't look back far enough" by making the time grid periodic:
      - t_grid covers 0..365 days as before
      - BUT when we look back before t=0, we wrap into the end of the year
      - Energy integrals are computed on a "2-year" tiled grid so we can query intervals that cross year boundary
    """

    def __init__(
        self,
        light_data: LightData,
        solar_power: SolarPower,
        mission: MissionProfile,
        soc_takeoff: float = 1.0,
        soc_target_at_sunset: float = 0.90,
        dt_minutes: float = 0.25,
    ):
        self.light = light_data
        self.solar = solar_power
        self.mission = mission

        self.days = self.light.days
        self.N = len(self.days)

        self.E_max = float(mission.E_battery_guess)
        self.E0 = float(soc_takeoff) * self.E_max
        self.E_target = float(soc_target_at_sunset) * self.E_max

        self.E_climb = float(mission.E_climb)
        self.P_req = float(mission.Pavg_cruise)

        self.dt = float(dt_minutes) * 60.0  # seconds
        self.year_s = 365.0 * 24.0 * 3600.0

        # sunrise/sunset in [0,year)
        self.sunrise_abs_s_1y = np.array([(i * 24.0 + self.light.sunrise_list[i]) * 3600.0 for i in range(self.N)])
        self.sunset_abs_s_1y = np.array([(i * 24.0 + self.light.sunset_list[i]) * 3600.0 for i in range(self.N)])

        # --- Build 1-year grid
        self.t_grid_1y = np.arange(0.0, self.year_s + 1e-9, self.dt)  # seconds
        self.M1 = len(self.t_grid_1y)
        day_idx_1y = np.clip((self.t_grid_1y // (24.0 * 3600.0)).astype(int), 0, self.N - 1)

        sunrise_t_1y = self.sunrise_abs_s_1y[day_idx_1y]
        sunset_t_1y = self.sunset_abs_s_1y[day_idx_1y]
        sun_up_1y = (self.t_grid_1y >= sunrise_t_1y) & (self.t_grid_1y <= sunset_t_1y)

        P_solar_day_1y = self.solar.power_W[day_idx_1y]
        P_solar_t_1y = np.where(sun_up_1y, P_solar_day_1y, 0.0)
        P_net_t_1y = P_solar_t_1y - self.P_req

        # --- Tile to 2-year arrays so we can integrate across year boundary easily
        self.t_grid_2y = np.concatenate([self.t_grid_1y, self.t_grid_1y[1:] + self.year_s])
        P_net_t_2y = np.concatenate([P_net_t_1y, P_net_t_1y[1:]])
        self.M2 = len(self.t_grid_2y)

        self.cum_net_energy_J_2y = np.zeros(self.M2, dtype=float)
        self.cum_net_energy_J_2y[1:] = np.cumsum(P_net_t_2y[:-1] * self.dt)

        self.target_soc = self.E_target / self.E_max

    def _idx_on_1y_grid(self, t_abs_s_1y: float) -> int:
        """Index on 1-year grid for time in [0,year)."""
        t = float(t_abs_s_1y) % self.year_s
        return int(np.clip(np.floor(t / self.dt), 0, self.M1 - 1))

    def _cumE_2y_at(self, t_abs_s_2y: float) -> float:
        """Cumulative energy at absolute time in [0,2*year]."""
        idx = int(np.clip(np.floor(t_abs_s_2y / self.dt), 0, self.M2 - 1))
        return float(self.cum_net_energy_J_2y[idx])

    def soc_at_sunset_day(self, t0_abs_s_2y: float, d: int, sunset_abs_s_2y: float) -> float:
        """
        SOC at a particular sunset time (in 2y coordinate) given takeoff time t0 (2y coordinate).
        O(1) using cumulative integral.
        """
        if t0_abs_s_2y > sunset_abs_s_2y:
            return -np.inf

        dE = self._cumE_2y_at(sunset_abs_s_2y) - self._cumE_2y_at(t0_abs_s_2y)
        E = self.E0 - self.E_climb + dE
        E = min(max(E, 0.0), self.E_max)
        return E / self.E_max

    def compute_bounds_periodic(self, horizon_days=3):
        """
        Periodic-year bounds: for each day d, compute feasible takeoff window in the interval
            [sunset(d) - horizon, sunset(d)]
        allowing the start of that interval to wrap into the previous year.

        Returns:
          lower_abs_s_1y, upper_abs_s_1y in *1-year absolute seconds* [0,year),
          feasible_day
        """
        H = float(horizon_days) * 24.0 * 3600.0

        lower_1y = np.full(self.N, np.nan, dtype=float)
        upper_1y = np.full(self.N, np.nan, dtype=float)
        feasible_day = np.zeros(self.N, dtype=bool)

        for d in range(self.N):
            sunset_1y = float(self.sunset_abs_s_1y[d])  # in [0,year)
            sunset_2y = sunset_1y + self.year_s  # place "evaluation sunset" in year 2 for convenience

            t_start_2y = sunset_2y - H
            t_end_2y = sunset_2y

            # Binary search in continuous indices on 2y grid
            lo = int(np.clip(np.floor(t_start_2y / self.dt), 0, self.M2 - 1))
            hi = int(np.clip(np.floor(t_end_2y / self.dt), 0, self.M2 - 1))

            def feasible_idx(idx):
                return self.soc_at_sunset_day(self.t_grid_2y[idx], d, sunset_2y) >= self.target_soc

            f_lo = feasible_idx(lo)
            f_hi = feasible_idx(hi)

            if (not f_lo) and (not f_hi):
                continue

            feasible_day[d] = True

            # find earliest feasible within window
            if f_lo:
                lo_feas = lo
            else:
                # smallest idx with feasible True (assumes monotonic-ish)
                left, right = lo, hi
                while right - left > 1:
                    mid = (left + right) // 2
                    if feasible_idx(mid):
                        right = mid
                    else:
                        left = mid
                lo_feas = right

            # find latest feasible within window
            if f_hi:
                hi_feas = hi
            else:
                # largest idx with feasible True
                left, right = lo, hi
                while right - left > 1:
                    mid = (left + right) // 2
                    if feasible_idx(mid):
                        left = mid
                    else:
                        right = mid
                hi_feas = left

            # Convert 2y times back to 1y [0,year) for plotting
            lower_1y[d] = self.t_grid_2y[lo_feas] % self.year_s
            upper_1y[d] = self.t_grid_2y[hi_feas] % self.year_s

        return lower_1y, upper_1y, feasible_day

    def plot_bounds_wrapped_0_24_on_axes(self, ax, lower_abs_s_1y, upper_abs_s_1y, feasible_day, title=""):
        """
        Plot on provided axes:
          - x: day of year
          - y: 0..24 local solar hour
        Uses wrapped fill for windows that cross midnight (top band).
        """
        days = self.days.astype(float)
        N = self.N

        lower_wrapped_h = np.full(N, np.nan, dtype=float)
        upper_wrapped_h = np.full(N, np.nan, dtype=float)
        wraps = np.zeros(N, dtype=bool)

        for i in range(N):
            if (not feasible_day[i]) or (not np.isfinite(lower_abs_s_1y[i])) or (not np.isfinite(upper_abs_s_1y[i])):
                continue

            day_start_s = i * 24.0 * 3600.0

            lo_rel_h = (lower_abs_s_1y[i] - day_start_s) / 3600.0
            hi_rel_h = (upper_abs_s_1y[i] - day_start_s) / 3600.0

            lo_h = lo_rel_h % 24.0
            hi_h = hi_rel_h % 24.0

            lower_wrapped_h[i] = lo_h
            upper_wrapped_h[i] = hi_h
            wraps[i] = lo_h > hi_h

        ax.plot(days, self.light.sunrise_list, "--", color="tab:blue", linewidth=1, label="Sunrise")
        ax.plot(days, self.light.sunset_list, "-.", color="tab:blue", linewidth=1, label="Sunset")
        ax.plot(days, lower_wrapped_h, color="tab:orange", linewidth=1.5, label="Lower bound")
        ax.plot(days, upper_wrapped_h, color="tab:red", linewidth=1.5, label="Upper bound")

        ok = feasible_day & np.isfinite(lower_wrapped_h) & np.isfinite(upper_wrapped_h)

        ok_nowrap = ok & (~wraps)
        if np.any(ok_nowrap):
            ax.fill_between(days, lower_wrapped_h, upper_wrapped_h, where=ok_nowrap, alpha=0.25, color="tab:green")

        ok_wrap = ok & wraps
        if np.any(ok_wrap):
            ax.fill_between(days, 0.0, upper_wrapped_h, where=ok_wrap, alpha=0.25, color="tab:green")
            ax.fill_between(days, lower_wrapped_h, 24.0, where=ok_wrap, alpha=0.25, color="tab:green")

        ax.set_xlim(1, 365)
        ax.set_ylim(0, 24)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)

    def plot_climb_overlay_on_axes(
        self,
        ax,
        takeoff_abs_s_1y,          # array (N,) or scalar: takeoff time in 1y seconds
        feasible_day=None,         # optional mask (N,)
        color="tab:purple",
        alpha=0.25,
        linewidth=0.8,
        label="Climb",
    ):
        """
        Draw the climb interval [takeoff, takeoff + t_climb] on the same (day vs 0..24h) axes.

        - If climb crosses midnight, it is wrapped (split into two pieces).
        - takeoff_abs_s_1y can be:
            * scalar: same takeoff time for all days (rare), OR
            * array of length N: one takeoff time per day (recommended).
        """
        N = self.N
        days = self.days.astype(float)

        if np.isscalar(takeoff_abs_s_1y):
            takeoff_abs_s_1y = np.full(N, float(takeoff_abs_s_1y), dtype=float)
        else:
            takeoff_abs_s_1y = np.asarray(takeoff_abs_s_1y, dtype=float)

        if feasible_day is None:
            feasible_day = np.isfinite(takeoff_abs_s_1y)
        else:
            feasible_day = np.asarray(feasible_day, dtype=bool) & np.isfinite(takeoff_abs_s_1y)

        t_climb_s = float(self.mission.t_climb)

        # Convert takeoff and end-of-climb to "hours of day" (wrapped 0..24)
        y0 = np.full(N, np.nan, dtype=float)   # takeoff hour
        y1 = np.full(N, np.nan, dtype=float)   # end-of-climb hour
        wraps = np.zeros(N, dtype=bool)

        for i in range(N):
            if not feasible_day[i]:
                continue

            day_start_s = i * 24.0 * 3600.0

            # relative seconds into day (can be negative/positive if time is wrapped in [0,year) already)
            t0_rel_s = (takeoff_abs_s_1y[i] - day_start_s)
            t1_rel_s = t0_rel_s + t_climb_s

            y0_i = (t0_rel_s / 3600.0) % 24.0
            y1_i = (t1_rel_s / 3600.0) % 24.0

            y0[i] = y0_i
            y1[i] = y1_i
            wraps[i] = y1_i < y0_i  # crossed midnight

        # Plot as short vertical-ish segments (actually filled band at fixed x)
        # Using fill_between for a “band” look:
        ok = feasible_day & np.isfinite(y0) & np.isfinite(y1)

        ok_nowrap = ok & (~wraps)
        if np.any(ok_nowrap):
            ax.vlines(days[ok_nowrap], y0[ok_nowrap], y1[ok_nowrap],
                    colors=color, alpha=alpha, linewidth=linewidth, label=label)

        ok_wrap = ok & wraps
        if np.any(ok_wrap):
            # two pieces: [y0, 24] and [0, y1]
            ax.vlines(days[ok_wrap], y0[ok_wrap], 24.0,
                    colors=color, alpha=alpha, linewidth=linewidth, label=label)
            ax.vlines(days[ok_wrap], 0.0, y1[ok_wrap],
                    colors=color, alpha=alpha, linewidth=linewidth, label=label)

# -----------------------------
# Example usage: subplots over latitude range, dual horizon
# -----------------------------
if __name__ == "__main__":
    solar_area = 31.24
    MTOW = 150.0
    M_battery = 65.43
    
    days = np.arange(1, 366)

    lats = [-60, -45, -30, -15, 0, 15, 30, 45, 60]
    horizons = [1]  # only show 1-day lookback

    nrows = len(lats)
    ncols = len(horizons)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(16, 2.4 * nrows),
        sharex=True,
        sharey=True,
    )

    for r, lat in enumerate(lats):
        light = LightData(latitude_deg=lat, days=days)
        solar = SolarPower(latitude_deg=lat, solar_area_m2=solar_area, days=days)
        mission = MissionProfile(time_daylight_cruise=3600 * 5, time_night=3600 * 19, solar_area=solar_area, MTOW=MTOW, M_battery=M_battery)

        grid = MultiDayDeployWindowBoundsFast(
            light_data=light,
            solar_power=solar,
            mission=mission,
            soc_takeoff=1.0,
            soc_target_at_sunset=0.90,
            dt_minutes=2.0,
        )

        for c, H in enumerate(horizons):
            ax = axes[r] if ncols == 1 else axes[r, c]

            lower, upper, feas = grid.compute_bounds_periodic(horizon_days=H)

            grid.plot_bounds_wrapped_0_24_on_axes(
                ax,
                lower,
                upper,
                feas,
                title=f"lat={lat:+.0f}°, lookback={H}d",
            )

            # ---- NEW: climb overlay ----
            # Visualize the climb interval starting at the latest feasible takeoff time each day
            takeoff = upper.copy()
            takeoff[~feas] = np.nan

            grid.plot_climb_overlay_on_axes(
                ax,
                takeoff_abs_s_1y=takeoff,
                feasible_day=feas,
                color="tab:purple",
                alpha=0.5,
                linewidth=1.0,
                label="Climb",
            )

            if c == 0:
                ax.set_ylabel("Takeoff time (solar hours)")

    # Legend
    handles, labels = (axes[0].get_legend_handles_labels() if ncols == 1
                       else axes[0, 0].get_legend_handles_labels())

    # Optional: de-duplicate legend labels (since climb overlay adds repeated entries)
    dedup = {}
    for h, lab in zip(handles, labels):
        dedup.setdefault(lab, h)

    fig.legend(list(dedup.values()), list(dedup.keys()),
               loc="upper right", bbox_to_anchor=(0.98, 0.98))

    fig.suptitle(f"Deploy windows + climb overlay (area={solar_area} m²)", y=0.995)
    plt.tight_layout(rect=[0, 0, 0.98, 0.985])
    plt.show()
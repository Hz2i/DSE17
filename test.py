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

        self.efficiency = 0.31
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
    def __init__(self, time_daylight_cruise, time_night):
        # guesses
        self.m_solar_guess = 9.2
        self.m_battery_guess = 80
        self.m_prop_guess = 10
        self.gamma_guess = 5
        self.S_guess = 36.5

        self.E_battery_guess = self.m_battery_guess * 400 * 3600  # J

        # given
        self.mass_subsys = 46.8
        self.g = 9.81
        self.alt = 60000 * 0.3048
        self.CD0 = 0.010
        self.AR = 24
        self.e = 0.9
        self.Pavg_climb_subsys = 300
        self.Pavg_cruise_subsys = 425
        self.eta_prop = 0.8
        self.LD = 40
        self.V_cruise = 25
        self.t_daylight_cruise = time_daylight_cruise
        self.t_night = time_night

        # derived
        self.m_total_guess = self.m_solar_guess + self.m_battery_guess + self.m_prop_guess + self.mass_subsys
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
# Grid-based bounds with wrap-aware plotting
# -----------------------------
class MultiDayDeployWindowBoundsFast:
    def __init__(
        self,
        light_data: LightData,
        solar_power: SolarPower,
        mission: MissionProfile,
        soc_takeoff: float = 1.0,
        soc_target_at_sunset: float = 0.90,
        dt_minutes: float = 10.0,
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

        self.sunrise_abs_s = np.array([(i * 24.0 + self.light.sunrise_list[i]) * 3600.0 for i in range(self.N)])
        self.sunset_abs_s = np.array([(i * 24.0 + self.light.sunset_list[i]) * 3600.0 for i in range(self.N)])

        self.t_end = 365.0 * 24.0 * 3600.0
        self.t_grid = np.arange(0.0, self.t_end + 1e-9, self.dt)
        self.M = len(self.t_grid)

        self.day_idx_grid = np.clip((self.t_grid // (24.0 * 3600.0)).astype(int), 0, self.N - 1)

        sunrise_t = self.sunrise_abs_s[self.day_idx_grid]
        sunset_t = self.sunset_abs_s[self.day_idx_grid]
        sun_up = (self.t_grid >= sunrise_t) & (self.t_grid <= sunset_t)

        P_solar_day = self.solar.power_W[self.day_idx_grid]
        self.P_solar_t = np.where(sun_up, P_solar_day, 0.0)

        self.P_net_t = self.P_solar_t - self.P_req
        self.cum_net_energy_J = np.zeros(self.M, dtype=float)
        self.cum_net_energy_J[1:] = np.cumsum(self.P_net_t[:-1] * self.dt)

        self.target_soc = self.E_target / self.E_max

    def _idx_at_time(self, t_abs_s: float) -> int:
        return int(np.clip(np.floor(t_abs_s / self.dt), 0, self.M - 1))

    def soc_at_sunset(self, t0_abs_s: float, d: int) -> float:
        t_sunset = float(self.sunset_abs_s[d])
        if t0_abs_s > t_sunset:
            return -np.inf

        i0 = self._idx_at_time(t0_abs_s)
        iset = self._idx_at_time(t_sunset)
        dE = self.cum_net_energy_J[iset] - self.cum_net_energy_J[i0]

        E = self.E0 - self.E_climb + dE
        E = min(max(E, 0.0), self.E_max)
        return E / self.E_max

    def _binary_search_first_true(self, d: int, lo: int, hi: int) -> int:
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.soc_at_sunset(self.t_grid[mid], d) >= self.target_soc:
                hi = mid
            else:
                lo = mid
        return hi

    def _binary_search_last_true(self, d: int, lo: int, hi: int) -> int:
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.soc_at_sunset(self.t_grid[mid], d) >= self.target_soc:
                lo = mid
            else:
                hi = mid
        return lo

    def compute_bounds(self, horizon_days=1):
        """
        Returns earliest and latest feasible takeoff times (absolute seconds) within:
            [sunset - horizon_days, sunset]
        """
        H = float(horizon_days) * 24.0 * 3600.0
        lower = np.full(self.N, np.nan, dtype=float)
        upper = np.full(self.N, np.nan, dtype=float)
        feasible_day = np.zeros(self.N, dtype=bool)

        for d in range(self.N):
            t_sunset = float(self.sunset_abs_s[d])
            t_start = max(0.0, t_sunset - H)

            lo = self._idx_at_time(t_start)
            hi = self._idx_at_time(t_sunset)

            f_lo = self.soc_at_sunset(self.t_grid[lo], d) >= self.target_soc
            f_hi = self.soc_at_sunset(self.t_grid[hi], d) >= self.target_soc

            if not f_lo and not f_hi:
                continue

            if f_lo and f_hi:
                feasible_day[d] = True
                lower[d] = self.t_grid[lo]
                upper[d] = self.t_grid[hi]
                continue

            if (not f_lo) and f_hi:
                first_true = self._binary_search_first_true(d, lo, hi)
                feasible_day[d] = True
                lower[d] = self.t_grid[first_true]
                upper[d] = self.t_grid[hi]
                continue

            if f_lo and (not f_hi):
                last_true = self._binary_search_last_true(d, lo, hi)
                feasible_day[d] = True
                lower[d] = self.t_grid[lo]
                upper[d] = self.t_grid[last_true]
                continue

        return lower, upper, feasible_day

    def plot_bounds_like_figure_wrapped(self, lower_abs_s, upper_abs_s, feasible_day, title=None, show=True):
        days = self.days.astype(float)
        N = self.N

        lower_wrapped = np.full(N, np.nan, dtype=float)
        upper_wrapped = np.full(N, np.nan, dtype=float)
        wraps = np.zeros(N, dtype=bool)

        for i in range(N):
            if (not feasible_day[i]) or (not np.isfinite(lower_abs_s[i])) or (not np.isfinite(upper_abs_s[i])):
                continue

            day_start_s = i * 24.0 * 3600.0
            lo_rel_h = (lower_abs_s[i] - day_start_s) / 3600.0
            hi_rel_h = (upper_abs_s[i] - day_start_s) / 3600.0

            lo_h = lo_rel_h % 24.0
            hi_h = hi_rel_h % 24.0

            lower_wrapped[i] = lo_h
            upper_wrapped[i] = hi_h
            wraps[i] = lo_h > hi_h

        plt.figure(figsize=(12, 5))
        plt.plot(days, self.light.sunrise_list, "--", color="tab:blue", label="Sunrise")
        plt.plot(days, self.light.sunset_list, "-.", color="tab:blue", label="Sunset")
        plt.plot(days, lower_wrapped, color="tab:orange", linewidth=2, label="Takeoff lower bound (grid)")
        plt.plot(days, upper_wrapped, color="tab:red", linewidth=2, label="Takeoff upper bound (grid)")

        ok = feasible_day & np.isfinite(lower_wrapped) & np.isfinite(upper_wrapped)

        ok_nowrap = ok & (~wraps)
        if np.any(ok_nowrap):
            plt.fill_between(days, lower_wrapped, upper_wrapped, where=ok_nowrap, alpha=0.25, color="tab:green", label="Feasible takeoff window")

        ok_wrap = ok & wraps
        if np.any(ok_wrap):
            plt.fill_between(days, 0.0, upper_wrapped, where=ok_wrap, alpha=0.25, color="tab:green")
            plt.fill_between(days, lower_wrapped, 24.0, where=ok_wrap, alpha=0.25, color="tab:green")

        plt.xlim(1, 365)
        plt.ylim(0, 24)
        plt.xlabel("Day of year (1 = Jan 1)")
        plt.ylabel("Takeoff time (solar hours)")
        plt.grid(True, alpha=0.3)

        if title is None:
            title = "Deploy window (grid, wrapped)"
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        if show:
            plt.show()


# -----------------------------
# Example usage: dual horizon
# -----------------------------
if __name__ == "__main__":
    solar_area = 37.0
    days = np.arange(1, 366)
    lats = [-60, -45, -30, -15, 0, 15, 30, 45, 60]

    for lat in lats:
        light = LightData(latitude_deg=lat, days=days)
        solar = SolarPower(latitude_deg=lat, solar_area_m2=solar_area, days=days)
        mission = MissionProfile(time_daylight_cruise=3600 * 5, time_night=3600 * 19)

        grid = MultiDayDeployWindowBoundsFast(
            light_data=light,
            solar_power=solar,
            mission=mission,
            soc_takeoff=1.0,
            soc_target_at_sunset=0.90,
            dt_minutes=1.0,
        )

        # 1-day lookback: shows the "top band" when feasible times include previous-day evening
        lower1, upper1, feas1 = grid.compute_bounds(horizon_days=1)
        grid.plot_bounds_like_figure_wrapped(
            lower1, upper1, feas1,
            title=f"Deploy window (grid, wrapped) lat={lat:+.0f}° — lookback 1 day"
        )

        # 3-day lookback: shows "always feasible within year" behavior (but wrap depends on modulo)
        lower3, upper3, feas3 = grid.compute_bounds(horizon_days=3)
        grid.plot_bounds_like_figure_wrapped(
            lower3, upper3, feas3,
            title=f"Deploy window (grid, wrapped) lat={lat:+.0f}° — lookback 3 days"
        )
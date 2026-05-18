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

        # Day-of-year convention: 1..365, where 1 == Jan 1
        self.days = np.arange(1, 366) if days is None else np.asarray(days, dtype=int)

        # Precompute arrays aligned with self.days
        self.declination_deg = self.solar_declination(self.days)
        self.daylight_hours = self.compute_daylight_hours(self.days)
        self.sunrise_list, self.sunset_list = self.compute_sunrise_sunset_lists(self.days)

    @staticmethod
    def solar_declination(day_of_year):
        """
        Approx solar declination in degrees.
        day_of_year is 1..365 where 1 = Jan 1.
        """
        day_of_year = np.asarray(day_of_year, dtype=float)
        return 23.45 * np.sin(np.deg2rad((360.0 / 365.0) * (day_of_year + 10.0)))

    def compute_daylight_hours(self, day_of_year):
        """
        Day length in hours using sunset hour angle:
        h0 = arccos(-tan(phi)*tan(delta))
        daylen = 2*h0*(180/pi)/15
        """
        day_of_year = np.asarray(day_of_year, dtype=float)
        delta_rad = np.deg2rad(self.solar_declination(day_of_year))
        phi = self.latitude_rad

        cos_h0 = -np.tan(phi) * np.tan(delta_rad)
        cos_h0 = np.clip(cos_h0, -1.0, 1.0)
        h0 = np.arccos(cos_h0)

        return (2.0 / 15.0) * np.rad2deg(h0)

    def compute_sunrise_sunset_lists(self, day_of_year):
        """
        Returns sunrise/sunset in local solar time hours [0..24],
        assuming solar noon at 12.0.
        """
        daylen = self.compute_daylight_hours(day_of_year)
        sunrise = 12.0 - 0.5 * daylen
        sunset = 12.0 + 0.5 * daylen
        return sunrise, sunset

    def get_by_day(self, day_of_year):
        """
        Convenience: query a single day-of-year (1..365).
        """
        idx = int(day_of_year) - 1
        return {
            "day": int(day_of_year),
            "daylight_hours": float(self.daylight_hours[idx]),
            "sunrise": float(self.sunrise_list[idx]),
            "sunset": float(self.sunset_list[idx]),
        }


# -----------------------------
# Solar power (daily average)
# -----------------------------
class SolarPower:
    def __init__(self, latitude_deg=0.0, solar_area_m2=1.0, days=None):
        self.latitude_deg = float(latitude_deg)
        self.latitude_rad = np.deg2rad(self.latitude_deg)

        self.days = np.arange(1, 366) if days is None else np.asarray(days, dtype=int)
        self.solar_area_m2 = float(solar_area_m2)

        # Simple model params (same as your code)
        self.efficiency = 0.31
        self.I0 = 1378.0            # W/m^2 (your constant)
        self.powLimS = 250.0        # W/m^2 cap

        # Precompute outputs aligned with self.days
        self.incidence_angle_rad = self.calc_daily_mean_incidence_angle(self.days)
        self.power_per_m2_W = self.calc_power_per_m2(self.days)
        self.power_W = self.power_per_m2_W * self.solar_area_m2

    @staticmethod
    def calc_solar_declination_rad(day_of_year):
        day_of_year = np.asarray(day_of_year, dtype=float)
        decl_deg = 23.45 * np.sin(np.deg2rad((360.0 / 365.0) * (day_of_year + 10.0)))
        return np.deg2rad(decl_deg)

    def calc_daily_mean_incidence_angle(self, day_of_year, n_samples=200):
        """
        Daily mean incidence angle for a horizontal surface.
        We average cos(theta) over daylight hour angles.
        """
        day_of_year = np.asarray(day_of_year, dtype=float)
        delta = self.calc_solar_declination_rad(day_of_year)
        phi = self.latitude_rad

        cos_h0 = -np.tan(phi) * np.tan(delta)
        cos_h0 = np.clip(cos_h0, -1.0, 1.0)
        h0 = np.arccos(cos_h0)

        incidence = np.zeros_like(day_of_year, dtype=float)
        for i in range(day_of_year.size):
            # if polar night/day ever occurs (not in your -30..30 range), handle gracefully
            if np.isclose(h0[i], 0.0):
                incidence[i] = np.pi / 2  # sun at horizon effectively => cos(theta) ~ 0
                continue

            hour_angles = np.linspace(-h0[i], h0[i], n_samples)
            cos_theta = (
                np.sin(delta[i]) * np.sin(phi)
                + np.cos(delta[i]) * np.cos(phi) * np.cos(hour_angles)
            )
            cos_theta = np.clip(cos_theta, 0.0, 1.0)

            mean_cos = np.mean(cos_theta)
            # incidence angle corresponding to mean cosine
            incidence[i] = np.arccos(np.clip(mean_cos, 0.0, 1.0))

        return incidence

    def calc_power_per_m2(self, day_of_year):
        incidence = self.calc_daily_mean_incidence_angle(day_of_year)
        cos_inc = np.cos(incidence)
        raw = self.efficiency * self.I0 * np.maximum(cos_inc, 0.0)
        return np.minimum(self.powLimS, raw)

    def get_by_day(self, day_of_year):
        idx = int(day_of_year) - 1
        return {
            "day": int(day_of_year),
            "power_per_m2_W": float(self.power_per_m2_W[idx]),
            "power_W": float(self.power_W[idx]),
        }



class MissionProfile:
    def __init__(self, time_daylight_cruise, time_night):  
        #Guesses
        self.m_solar_guess = 9.2 #kg
        self.m_battery_guess = 80 #kg
        self.m_prop_guess = 10 #kg
        self.gamma_guess = 5 #deg
        self.S_guess = 36.5 #m^2
        self.E_battery_guess = self.m_battery_guess * 400 * 3600 #J

        #Given
        self.mass_subsys = 46.8 #kg
        self.g = 9.81                                     
        self.alt = 60000*0.3048
        self.CD0 = 0.010
        self.AR = 24
        self.e = 0.9
        self.Pavg_climb_subsys = 300 #W
        self.Pavg_cruise_subsys = 425 #W
        self.eta_prop = 0.8
        self.LD = 40
        self.V_cruise = 25 #m/s
        self.t_daylight_cruise = time_daylight_cruise #s
        self.t_night = time_night #s
        self.specific_power = 400 #Wh/kg Battery specific power
        self.specific_mass_solar = 4 #kg/m^2

        #Extra Parameters
        self.m_total_guess = self.m_solar_guess + self.m_battery_guess + self.m_prop_guess + self.mass_subsys
        self.density_climb = self.Calc_Density_Climb()
        self.Cl_opt_climb = self.Calc_Cl_opt_climb()
        self.CD_total_climb = self.Calc_CD_total(self.Cl_opt_climb)
        self.gamma_rad = np.radians(self.gamma_guess)

        
        #Climb
        self.V_climb = self.Calc_V_climb()
        self.t_climb = self.Calc_t_climb()
        self.D_climb = self.Calc_D_climb()
        self.E_climb = self.Calc_E_climb()
        self.ROC = self.Calc_ROC()

        #Cruise
        self.Pprop_cruise = self.Calc_Pprop_cruise()
        self.Pavg_cruise = self.Pprop_cruise + self.Pavg_cruise_subsys
        self.E_cruise = self.Calc_E_cruise()

    def Calc_Density_Climb(self):
        h_range = np.linspace(0, self.alt, 100000)
        density = am.Atmosphere(h_range).density
        #integrate density to get mass of air in the climb envelope
        dens_climb = np.trapz(density, h_range)/(self.alt)
        return dens_climb

    def Calc_Cl_opt_climb(self):
        return np.sqrt(self.CD0 * np.pi * self.AR * self.e)
    
    def Calc_CD_total(self, CL):
        k = 1.0/(np.pi * self.AR * self.e)
        return self.CD0 + k * CL**2
    
    def Calc_V_climb(self):
        V_climb = np.sqrt(2*self.m_total_guess*self.g*np.cos(self.gamma_rad)/(self.density_climb*self.S_guess*self.Cl_opt_climb))
        return V_climb
    
    def Calc_D_climb(self):
        D_climb = 0.5*self.CD_total_climb*self.S_guess*self.density_climb*self.V_climb**2
        return D_climb
    
    def Calc_t_climb(self):
        t_climb = self.alt / (self.V_climb * np.sin(self.gamma_rad))
        return t_climb

    def Calc_E_climb(self):
        Epot = self.m_total_guess * self.g * self.alt
        Edrag = self.D_climb * (self.alt / np.sin(self.gamma_rad))
        Esubsys_climb = self.Pavg_climb_subsys * self.t_climb
        Eclimb = Epot + Edrag + Esubsys_climb
        return Eclimb
    
    def Calc_ROC(self):
        ROC = self.V_climb * np.sin(self.gamma_rad)
        return ROC

    def Calc_Pprop_cruise(self):
        Pprop_cruise = ((self.m_total_guess*self.g)/(self.LD))*self.V_cruise/self.eta_prop
        return Pprop_cruise

    def Calc_E_cruise(self):
        Ecruise = (self.Pavg_cruise) * self.t_daylight_cruise
        return Ecruise

class DeployWindow:
    """
    Computes per-day takeoff time bounds such that SOC at sunset is >= soc_target_at_sunset.

    Times are in *solar time hours* (0..24), using LightData sunrise/sunset.
    """

    def __init__(
        self,
        latitude_deg: float,
        light_data: LightData,
        solar_power: SolarPower,
        mission: MissionProfile,
        soc_takeoff: float = 1.0,
        soc_target_at_sunset: float = 0.90,
        allow_ge: bool = True,
    ):
        self.latitude_deg = float(latitude_deg)
        self.light = light_data
        self.solar = solar_power
        self.mission = mission

        self.days = self.light.days  # 1..365

        # Battery energy model
        self.E_batt_max = float(mission.E_battery_guess)  # J
        self.E_takeoff = float(soc_takeoff) * self.E_batt_max
        self.E_target_sunset = float(soc_target_at_sunset) * self.E_batt_max

        # Constant energies/powers (your statement)
        self.E_climb = float(mission.E_climb)            # J (constant)
        self.P_req_day = float(mission.Pavg_cruise)      # W electrical power required during daylight cruise

        self.allow_ge = bool(allow_ge)

        # Outputs (per day)
        self.takeoff_lower_bound = np.full_like(self.days, np.nan, dtype=float)  # hours
        self.takeoff_upper_bound = np.full_like(self.days, np.nan, dtype=float)  # hours
        self.feasible = np.zeros_like(self.days, dtype=bool)

        self._compute_bounds_all_days()

        self.plot()

    def _soc_at_sunset_for_takeoff(self, day_idx: int, takeoff_time_h: float) -> float:
        """
        SOC at sunset allowing takeoff before sunrise.
        - No solar before sunrise.
        - Solar during [sunrise, sunset] using daily-average solar power for that day.
        """
        sunrise = float(self.light.sunrise_list[day_idx])
        sunset = float(self.light.sunset_list[day_idx])

        P_solar = float(self.solar.power_W[day_idx])  # W (daily mean during daylight)
        P_req = self.P_req_day                         # W

        # Clamp takeoff time to [0, sunset] for this simple model
        t_to = float(np.clip(takeoff_time_h, 0.0, sunset))

        # Segment durations (seconds)
        t_pre_sun_s = max(0.0, (sunrise - t_to) * 3600.0)              # takeoff -> sunrise (no solar)
        t_day_s = max(0.0, (sunset - max(t_to, sunrise)) * 3600.0)     # max(takeoff,sunrise) -> sunset (solar)

        E_sunset = (
            self.E_takeoff
            - self.E_climb
            - P_req * t_pre_sun_s
            + (P_solar - P_req) * t_day_s
        )
        return E_sunset / self.E_batt_max

    def _compute_bounds_one_day(self, day_idx: int):
        sunrise = float(self.light.sunrise_list[day_idx])
        sunset = float(self.light.sunset_list[day_idx])
        if sunset <= sunrise:
            return np.nan, np.nan, False

        P_solar = float(self.solar.power_W[day_idx])
        P_req = self.P_req_day
        target_soc = self.E_target_sunset / self.E_batt_max

        # Helpful constants for the inequality SOC(sunset) >= target
        # We'll derive bounds for two regions:
        #   A) takeoff in [0, sunrise]
        #   B) takeoff in [sunrise, sunset]

        # Energy at sunset if you takeoff exactly at sunrise:
        # E = E_takeoff - E_climb + (P_solar - P_req)*(sunset - sunrise)*3600
        E_at_sunset_if_to_at_sunrise = (
            self.E_takeoff
            - self.E_climb
            + (P_solar - P_req) * (sunset - sunrise) * 3600.0
        )

        # --- Region B: takeoff after sunrise (same shape as before) ---
        # For t_to >= sunrise:
        # E_sunset = E_takeoff - E_climb + (P_solar - P_req)*(sunset - t_to)*3600
        netP_day = (P_solar - P_req)

        # Candidate interval in [sunrise, sunset]
        lower_B, upper_B, ok_B = np.nan, np.nan, False
        if np.isclose(netP_day, 0.0, atol=1e-9):
            soc_any_B = (self.E_takeoff - self.E_climb) / self.E_batt_max
            if self.allow_ge:
                ok_B = soc_any_B >= target_soc
                if ok_B:
                    lower_B, upper_B = sunrise, sunset
            else:
                ok_B = np.isclose(soc_any_B, target_soc)
                if ok_B:
                    lower_B = upper_B = sunrise  # arbitrary; equality holds for all times actually
        else:
            # Solve equality for t_required in region B
            # E_target = E_takeoff - E_climb + netP_day*(sunset - t_to)*3600
            A = self.E_target_sunset - (self.E_takeoff - self.E_climb)
            t_required = sunset - (A / (netP_day * 3600.0))

            if self.allow_ge:
                if netP_day > 0:
                    # earlier => more net charge => higher SOC, require t_to <= t_required
                    lower_B = sunrise
                    upper_B = min(sunset, t_required)
                else:
                    # netP_day < 0: earlier => more discharge, require t_to >= t_required
                    lower_B = max(sunrise, t_required)
                    upper_B = sunset
            else:
                lower_B = upper_B = t_required

            ok_B = np.isfinite(lower_B) and np.isfinite(upper_B) and (upper_B >= lower_B) and not (upper_B < sunrise or lower_B > sunset)
            if ok_B:
                lower_B = max(lower_B, sunrise)
                upper_B = min(upper_B, sunset)

        # --- Region A: takeoff before sunrise ---
        # For t_to <= sunrise:
        # E_sunset = (E_takeoff - E_climb)
        #           - P_req*(sunrise - t_to)*3600
        #           + (P_solar - P_req)*(sunset - sunrise)*3600
        #
        # This is linear in t_to with slope +P_req*3600 (later takeoff => less pre-sun consumption => higher E)
        #
        # Solve for t_to such that E_sunset >= E_target
        # E_sunset = E_base + P_req*(t_to - sunrise)*3600, where:
        E_base = E_at_sunset_if_to_at_sunrise  # that's E when t_to = sunrise
        # Then for t_to <= sunrise:
        # E_sunset(t_to) = E_base + P_req*(t_to - sunrise)*3600

        lower_A, upper_A, ok_A = np.nan, np.nan, False
        if np.isclose(P_req, 0.0, atol=1e-12):
            # No consumption before sunrise -> either always ok or never ok
            soc_any_A = E_base / self.E_batt_max
            if self.allow_ge:
                ok_A = soc_any_A >= target_soc
                if ok_A:
                    lower_A, upper_A = 0.0, sunrise
            else:
                ok_A = np.isclose(soc_any_A, target_soc)
                if ok_A:
                    lower_A, upper_A = sunrise, sunrise
        else:
            # Inequality: E_base + P_req*(t_to - sunrise)*3600 >= E_target
            # => t_to >= sunrise + (E_target - E_base)/(P_req*3600)
            t_min = sunrise + (self.E_target_sunset - E_base) / (P_req * 3600.0)

            if self.allow_ge:
                # feasible takeoffs are [max(0, t_min), sunrise]
                lower_A = max(0.0, t_min)
                upper_A = sunrise
            else:
                lower_A = upper_A = t_min

            ok_A = np.isfinite(lower_A) and np.isfinite(upper_A) and (upper_A >= lower_A) and not (upper_A < 0.0 or lower_A > sunrise)
            if ok_A:
                lower_A = max(lower_A, 0.0)
                upper_A = min(upper_A, sunrise)

        # --- Combine feasible regions A and B ---
        candidates = []
        if ok_A:
            candidates.append((lower_A, upper_A))
        if ok_B:
            candidates.append((lower_B, upper_B))

        if not candidates:
            return np.nan, np.nan, False

        # Merge intervals if they touch/overlap, else keep overall min/max as "window".
        # (If you want to preserve disjoint intervals, store them separately.)
        lower = min(c[0] for c in candidates)
        upper = max(c[1] for c in candidates)

        # Validate quickly
        soc_lo = self._soc_at_sunset_for_takeoff(day_idx, lower)
        soc_hi = self._soc_at_sunset_for_takeoff(day_idx, upper)
        ok = (soc_lo >= target_soc - 1e-9) or (soc_hi >= target_soc - 1e-9) if self.allow_ge else True

        return lower, upper, ok

    def _compute_bounds_all_days(self):
        for i in range(len(self.days)):
            lower, upper, ok = self._compute_bounds_one_day(i)
            self.takeoff_lower_bound[i] = lower
            self.takeoff_upper_bound[i] = upper
            self.feasible[i] = ok

    def get_bounds_by_day(self, day_of_year: int):
        idx = int(day_of_year) - 1
        return {
            "day": int(day_of_year),
            "sunrise": float(self.light.sunrise_list[idx]),
            "sunset": float(self.light.sunset_list[idx]),
            "takeoff_lower_bound": float(self.takeoff_lower_bound[idx]),
            "takeoff_upper_bound": float(self.takeoff_upper_bound[idx]),
            "feasible": bool(self.feasible[idx]),
            "power_day_W": float(self.solar.power_W[idx]),
        }
    
    def plot(self, title=None, show_sun=True, show=True):
        days = self.days

        plt.figure(figsize=(12, 5))

        # Optional: show sunrise/sunset curves as context
        if show_sun:
            plt.plot(days, self.light.sunrise_list, color="tab:blue", linestyle="--", linewidth=1.5, label="Sunrise")
            plt.plot(days, self.light.sunset_list,  color="tab:blue", linestyle="-.", linewidth=1.5, label="Sunset")

        # Deploy window bounds
        plt.plot(days, self.takeoff_lower_bound, color="tab:orange", linewidth=2, label="Takeoff lower bound")
        plt.plot(days, self.takeoff_upper_bound, color="tab:red", linewidth=2, label="Takeoff upper bound")

        # Shade feasible region between bounds (only where feasible)
        feasible = self.feasible & np.isfinite(self.takeoff_lower_bound) & np.isfinite(self.takeoff_upper_bound)
        if np.any(feasible):
            plt.fill_between(
                days,
                self.takeoff_lower_bound,
                self.takeoff_upper_bound,
                where=feasible,
                color="tab:green",
                alpha=0.25,
                label="Feasible takeoff window",
                interpolate=True,
            )

        plt.xlim(1, 365)
        plt.ylim(0, 24)
        plt.xlabel("Day of year (1 = Jan 1)")
        plt.ylabel("Takeoff time (solar hours)")
        plt.grid(True, alpha=0.3)

        if title is None:
            title = f"Deploy window (SOC@sunset >= target) at lat={self.latitude_deg:+.0f}°"
        plt.title(title)

        plt.legend()
        plt.tight_layout()
        if show:
            plt.show()



class MultiDayDeployWindowBoundsFast:
    """
    Compute per-day [lower, upper] takeoff time bounds using a continuous time grid + cumulative integral.

    Deploy criterion for day d:
        SOC at *that day's sunset* >= target.

    Model:
      - Takeoff gives one-time energy cost E_climb
      - After takeoff: battery evolves with P_net(t) = P_solar(t) - P_req
      - P_solar(t) = solar.power_W[day] when sun is up, else 0
      - Battery clipped to [0, E_max] (optional, kept)
    """

    def __init__(
        self,
        light_data,
        solar_power,
        mission,
        soc_takeoff=1.0,
        soc_target_at_sunset=0.90,
        dt_minutes=10.0,
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

        # Absolute sunrise/sunset in seconds
        self.sunrise_abs_s = np.array([(i * 24.0 + self.light.sunrise_list[i]) * 3600.0 for i in range(self.N)])
        self.sunset_abs_s  = np.array([(i * 24.0 + self.light.sunset_list[i])  * 3600.0 for i in range(self.N)])

        # Uniform time grid for whole year
        self.t_end = 365.0 * 24.0 * 3600.0
        self.t_grid = np.arange(0.0, self.t_end + 1e-9, self.dt)  # seconds
        self.M = len(self.t_grid)

        # Day index for each timestep
        self.day_idx_grid = np.clip((self.t_grid // (24.0 * 3600.0)).astype(int), 0, self.N - 1)

        # Sun-up mask for each timestep
        sunrise_t = self.sunrise_abs_s[self.day_idx_grid]
        sunset_t  = self.sunset_abs_s[self.day_idx_grid]
        sun_up = (self.t_grid >= sunrise_t) & (self.t_grid <= sunset_t)

        # Solar power per timestep
        P_solar_day = self.solar.power_W[self.day_idx_grid]
        self.P_solar_t = np.where(sun_up, P_solar_day, 0.0)

        # Net power and cumulative net energy
        self.P_net_t = self.P_solar_t - self.P_req
        self.cum_net_energy_J = np.zeros(self.M, dtype=float)
        self.cum_net_energy_J[1:] = np.cumsum(self.P_net_t[:-1] * self.dt)

        self.target_soc = self.E_target / self.E_max

    def _idx_at_time(self, t_abs_s: float) -> int:
        return int(np.clip(np.floor(t_abs_s / self.dt), 0, self.M - 1))

    def soc_at_sunset(self, t0_abs_s: float, d: int) -> float:
        """SOC at sunset of day d if takeoff is at absolute time t0_abs_s."""
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
        """
        Find smallest index in [lo, hi] s.t. feasible(index) == True.
        Assumes feasible(lo)==False and feasible(hi)==True.
        """
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.soc_at_sunset(self.t_grid[mid], d) >= self.target_soc:
                hi = mid
            else:
                lo = mid
        return hi

    def _binary_search_last_true(self, d: int, lo: int, hi: int) -> int:
        """
        Find largest index in [lo, hi] s.t. feasible(index) == True.
        Assumes feasible(lo)==True and feasible(hi)==False.
        """
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.soc_at_sunset(self.t_grid[mid], d) >= self.target_soc:
                lo = mid
            else:
                hi = mid
        return lo

    def compute_bounds(self, horizon_days=365):
        """
        For each day d, compute earliest and latest takeoff times (absolute seconds) within a search horizon
        [sunset - horizon, sunset].

        Returns:
          lower_abs_s, upper_abs_s, feasible_day (arrays length N)
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

            # Evaluate feasibility at endpoints
            f_lo = self.soc_at_sunset(self.t_grid[lo], d) >= self.target_soc
            f_hi = self.soc_at_sunset(self.t_grid[hi], d) >= self.target_soc

            if not f_lo and not f_hi:
                # Could still be feasible in the middle (non-monotonic), but for your simplified model
                # it should be monotonic enough; treat as infeasible.
                continue

            # If feasible everywhere in [t_start, t_sunset], bounds are whole interval
            if f_lo and f_hi:
                feasible_day[d] = True
                lower[d] = self.t_grid[lo]
                upper[d] = self.t_grid[hi]
                continue

            # If f_lo False and f_hi True: transition from infeasible->feasible somewhere (earliest bound)
            if (not f_lo) and f_hi:
                # earliest feasible
                first_true = self._binary_search_first_true(d, lo, hi)
                feasible_day[d] = True
                lower[d] = self.t_grid[first_true]
                upper[d] = self.t_grid[hi]  # latest is sunset
                continue

            # If f_lo True and f_hi False: transition from feasible->infeasible somewhere (latest bound)
            if f_lo and (not f_hi):
                last_true = self._binary_search_last_true(d, lo, hi)
                feasible_day[d] = True
                lower[d] = self.t_grid[lo]  # earliest is horizon start
                upper[d] = self.t_grid[last_true]
                continue

        return lower, upper, feasible_day
    
    def plot_bounds_window(self, lower_abs_s, upper_abs_s, feasible_day, day_start=1, window_days=28, title=None, show=True):
    """28-day window plot for speed/clarity (no wrap handling)."""
    s = int(day_start) - 1
    e = min(self.N, s + int(window_days))

    days = self.days[s:e].astype(float)

    lower_local_h = np.full(e - s, np.nan, dtype=float)
    upper_local_h = np.full(e - s, np.nan, dtype=float)

    for j, i in enumerate(range(s, e)):
        day_start_s = i * 24.0 * 3600.0
        if np.isfinite(lower_abs_s[i]):
            lower_local_h[j] = (lower_abs_s[i] - day_start_s) / 3600.0
        if np.isfinite(upper_abs_s[i]):
            upper_local_h[j] = (upper_abs_s[i] - day_start_s) / 3600.0

    plt.figure(figsize=(12, 5))
    plt.plot(days, self.light.sunrise_list[s:e], "--", color="tab:blue", label="Sunrise")
    plt.plot(days, self.light.sunset_list[s:e], "-.", color="tab:blue", label="Sunset")
    plt.plot(days, lower_local_h, color="tab:orange", linewidth=2, label="Takeoff lower bound (grid)")
    plt.plot(days, upper_local_h, color="tab:red", linewidth=2, label="Takeoff upper bound (grid)")

    ok = feasible_day[s:e] & np.isfinite(lower_local_h) & np.isfinite(upper_local_h)
    if np.any(ok):
        plt.fill_between(days, lower_local_h, upper_local_h, where=ok, alpha=0.25, color="tab:green", label="Feasible takeoff window")

    plt.xlim(days[0], days[-1])
    plt.ylim(0, 24)
    plt.xlabel("Day of year (1 = Jan 1)")
    plt.ylabel("Takeoff time (solar hours)")
    plt.grid(True, alpha=0.3)

    if title is None:
        title = f"Deploy window (days {int(days[0])}-{int(days[-1])}) — time-grid"
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    if show:
        plt.show()


def plot_bounds_window_wrapped(self, lower_abs_s, upper_abs_s, feasible_day, day_start=1, window_days=28, title=None, show=True):
    """
    28-day window plot WITH wrap handling (shows the missing 'top' band when the window crosses midnight).
    """
    s = int(day_start) - 1
    e = min(self.N, s + int(window_days))

    days = self.days[s:e].astype(float)
    n = e - s

    lower_wrapped = np.full(n, np.nan, dtype=float)
    upper_wrapped = np.full(n, np.nan, dtype=float)
    wraps = np.zeros(n, dtype=bool)

    for j, i in enumerate(range(s, e)):
        if (not feasible_day[i]) or (not np.isfinite(lower_abs_s[i])) or (not np.isfinite(upper_abs_s[i])):
            continue

        day_start_s = i * 24.0 * 3600.0
        lo_rel_h = (lower_abs_s[i] - day_start_s) / 3600.0
        hi_rel_h = (upper_abs_s[i] - day_start_s) / 3600.0

        lo_h = lo_rel_h % 24.0
        hi_h = hi_rel_h % 24.0

        lower_wrapped[j] = lo_h
        upper_wrapped[j] = hi_h
        wraps[j] = lo_h > hi_h

    plt.figure(figsize=(12, 5))
    plt.plot(days, self.light.sunrise_list[s:e], "--", color="tab:blue", label="Sunrise")
    plt.plot(days, self.light.sunset_list[s:e], "-.", color="tab:blue", label="Sunset")
    plt.plot(days, lower_wrapped, color="tab:orange", linewidth=2, label="Takeoff lower bound (grid)")
    plt.plot(days, upper_wrapped, color="tab:red", linewidth=2, label="Takeoff upper bound (grid)")

    ok = feasible_day[s:e] & np.isfinite(lower_wrapped) & np.isfinite(upper_wrapped)

    ok_nowrap = ok & (~wraps)
    if np.any(ok_nowrap):
        plt.fill_between(days, lower_wrapped, upper_wrapped, where=ok_nowrap, alpha=0.25, color="tab:green", label="Feasible takeoff window")

    ok_wrap = ok & wraps
    if np.any(ok_wrap):
        plt.fill_between(days, 0.0, upper_wrapped, where=ok_wrap, alpha=0.25, color="tab:green")
        plt.fill_between(days, lower_wrapped, 24.0, where=ok_wrap, alpha=0.25, color="tab:green")

    plt.xlim(days[0], days[-1])
    plt.ylim(0, 24)
    plt.xlabel("Day of year (1 = Jan 1)")
    plt.ylabel("Takeoff time (solar hours)")
    plt.grid(True, alpha=0.3)

    if title is None:
        title = f"Deploy window (days {int(days[0])}-{int(days[-1])}) — time-grid (wrapped)"
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    if show:
        plt.show()

    def plot_bounds_like_figure(self, lower_abs_s, upper_abs_s, feasible_day, title=None, show=True):
        """
        Plot per-day lower/upper bounds in local solar hours (0..24), like your DeployWindow plot.
        """
        days = self.days.astype(float)

        lower_local_h = np.full(self.N, np.nan, dtype=float)
        upper_local_h = np.full(self.N, np.nan, dtype=float)

        for i in range(self.N):
            day_start_s = i * 24.0 * 3600.0
            if np.isfinite(lower_abs_s[i]):
                lower_local_h[i] = (lower_abs_s[i] - day_start_s) / 3600.0
            if np.isfinite(upper_abs_s[i]):
                upper_local_h[i] = (upper_abs_s[i] - day_start_s) / 3600.0

        plt.figure(figsize=(12, 5))
        plt.plot(days, self.light.sunrise_list, "--", color="tab:blue", label="Sunrise")
        plt.plot(days, self.light.sunset_list, "-.", color="tab:blue", label="Sunset")
        plt.plot(days, lower_local_h, color="tab:orange", linewidth=2, label="Takeoff lower bound (grid)")
        plt.plot(days, upper_local_h, color="tab:red", linewidth=2, label="Takeoff upper bound (grid)")

        ok = feasible_day & np.isfinite(lower_local_h) & np.isfinite(upper_local_h)
        if np.any(ok):
            plt.fill_between(days, lower_local_h, upper_local_h, where=ok, alpha=0.25, color="tab:green", label="Feasible takeoff window")

        plt.xlim(1, 365)
        plt.ylim(0, 24)
        plt.xlabel("Day of year (1 = Jan 1)")
        plt.ylabel("Takeoff time (solar hours)")
        plt.grid(True, alpha=0.3)

        if title is None:
            title = "Deploy window (SOC@sunset >= target) — time-grid"
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        if show:
            plt.show()

    def plot_bounds_like_figure_wrapped(self, lower_abs_s, upper_abs_s, feasible_day, title=None, show=True):
        """
        Like plot_bounds_like_figure(), but correctly shows windows that cross midnight by wrapping bounds into [0,24)
        and splitting the filled region into bottom+top parts when needed.
        """
        days = self.days.astype(float)
        N = self.N

        sunrise = self.light.sunrise_list
        sunset = self.light.sunset_list

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

            # interval crosses midnight in wrapped hour-of-day space
            wraps[i] = lo_h > hi_h

        plt.figure(figsize=(12, 5))
        plt.plot(days, sunrise, "--", color="tab:blue", label="Sunrise")
        plt.plot(days, sunset,  "-.", color="tab:blue", label="Sunset")
        plt.plot(days, lower_wrapped, color="tab:orange", linewidth=2, label="Takeoff lower bound (grid)")
        plt.plot(days, upper_wrapped, color="tab:red", linewidth=2, label="Takeoff upper bound (grid)")

        ok = feasible_day & np.isfinite(lower_wrapped) & np.isfinite(upper_wrapped)

        # Non-wrapping fill: [lower, upper]
        ok_nowrap = ok & (~wraps)
        if np.any(ok_nowrap):
            plt.fill_between(
                days, lower_wrapped, upper_wrapped,
                where=ok_nowrap,
                alpha=0.25, color="tab:green",
                label="Feasible takeoff window",
                interpolate=True,
            )

        # Wrapping fill: [0, upper] U [lower, 24]
        ok_wrap = ok & wraps
        if np.any(ok_wrap):
            plt.fill_between(
                days, 0.0, upper_wrapped,
                where=ok_wrap,
                alpha=0.25, color="tab:green",
                interpolate=True,
            )
            plt.fill_between(
                days, lower_wrapped, 24.0,
                where=ok_wrap,
                alpha=0.25, color="tab:green",
                interpolate=True,
            )

        plt.xlim(1, 365)
        plt.ylim(0, 24)
        plt.xlabel("Day of year (1 = Jan 1)")
        plt.ylabel("Takeoff time (solar hours)")
        plt.grid(True, alpha=0.3)

        if title is None:
            title = "Deploy window (grid, wrapped) at lat"
        plt.title(title)
        plt.legend()
        plt.tight_layout()
        if show:
            plt.show()
# -----------------------------
# Example usage (fixed for MultiDayDeployWindow)
# -----------------------------
if __name__ == "__main__":
    solar_area = 37.0
    days = np.arange(1, 366)
    lat = -45

    light = LightData(latitude_deg=lat, days=days)
    solar = SolarPower(latitude_deg=lat, solar_area_m2=solar_area, days=days)
    mission = MissionProfile(time_daylight_cruise=3600 * 5, time_night=3600 * 19)

    grid = MultiDayDeployWindowBoundsFast(
        light_data=light,
        solar_power=solar,
        mission=mission,
        soc_takeoff=1.0,
        soc_target_at_sunset=0.90,
        dt_minutes=10.0,
    )

    lower_abs_s, upper_abs_s, feasible_day = grid.compute_bounds(horizon_days=365)

    grid.plot_bounds_window(
        lower_abs_s, upper_abs_s, feasible_day,
        day_start=1, window_days=28,
        title="First 28 days (grid)"
    )

    # Use wrapped plot so the "top" appears
    grid.plot_bounds_like_figure_wrapped(
        lower_abs_s, upper_abs_s, feasible_day,
        title=f"Deploy window (grid) at lat={lat:+.0f}°"
    )

    # Example: subplots over a latitude range
    def plot_lat_range_subplots(latitudes, solar_area=37.0, days=None, mission=None, dt_minutes=30.0, horizon_days=365, window_days=28, show=True):
        if days is None:
            days = np.arange(1, 366)
        if mission is None:
            mission = MissionProfile(time_daylight_cruise=3600 * 5, time_night=3600 * 19)

        n = len(latitudes)
        cols = min(3, n)
        rows = (n + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 3.5 * rows), sharex=True, sharey=True)
        if rows * cols == 1:
            axes = np.array([axes])
        axes = axes.flatten()

        for ax_idx, lat_val in enumerate(latitudes):
            ax = axes[ax_idx]
            light = LightData(latitude_deg=lat_val, days=days)
            solar = SolarPower(latitude_deg=lat_val, solar_area_m2=solar_area, days=days)
            grid = MultiDayDeployWindowBoundsFast(light, solar, mission, dt_minutes=dt_minutes)
            lower_abs_s, upper_abs_s, feasible_day = grid.compute_bounds(horizon_days=horizon_days)

            # Prepare day window
            s = 0
            e = min(len(days), int(window_days))
            day_vals = days[s:e].astype(float)

            lower_h = np.full(e - s, np.nan)
            upper_h = np.full(e - s, np.nan)
            for j, i in enumerate(range(s, e)):
                day_start_s = i * 24.0 * 3600.0
                if np.isfinite(lower_abs_s[i]):
                    lower_h[j] = (lower_abs_s[i] - day_start_s) / 3600.0 % 24.0
                if np.isfinite(upper_abs_s[i]):
                    upper_h[j] = (upper_abs_s[i] - day_start_s) / 3600.0 % 24.0

            ax.plot(day_vals, light.sunrise_list[s:e], '--', color='tab:blue', linewidth=1)
            ax.plot(day_vals, light.sunset_list[s:e], '-.', color='tab:blue', linewidth=1)
            ax.plot(day_vals, lower_h, color='tab:orange')
            ax.plot(day_vals, upper_h, color='tab:red')

            ok = feasible_day[s:e] & np.isfinite(lower_h) & np.isfinite(upper_h)
            if np.any(ok):
                ax.fill_between(day_vals, lower_h, upper_h, where=ok, alpha=0.25, color='tab:green')

            ax.set_title(f'lat={lat_val:+.0f}°')
            ax.set_ylim(0, 24)
            ax.grid(True, alpha=0.3)

        # Hide unused axes
        for k in range(n, len(axes)):
            axes[k].axis('off')

        fig.suptitle(f'Deploy windows for latitudes {latitudes[0]}..{latitudes[-1]}')
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        if show:
            plt.show()

    # Demo of subplots
    plot_lat_range_subplots(np.linspace(-60, 60, 6), solar_area=solar_area, days=days, mission=mission, dt_minutes=20.0, horizon_days=60, window_days=21)
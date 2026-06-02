import numpy as np

class link_budget:
    def __init__(self, distance, frequency=0.12e9, data_rate_bps=1e6, modulation='QPSK', target_link_margin_db=10):
        self.frequency = frequency                              # set at 31 GHz but ranges include: 31-31.3 GHz, 47.2–47.5 GHz and 47.9–48.2 GHz
        self.distance = distance                                # input
        self.data_rate_bps = data_rate_bps                      # requirement: 1 Mbit/s
        self.modulation = modulation.upper()
        self.target_link_margin_db = target_link_margin_db
        self.c = 3e8                                            # Speed of light (m/s)

        self.snr_requirements_db = {    # signal-to-noise from references
            'BPSK': 10.0,
            'QPSK': 13.0,
            '16QAM': 20.0,
            '64QAM': 26.0
        }

        self.bandwidth_efficiency = {   # bandwidth efficiency from references
            'BPSK': 1.0,
            'QPSK': 2.0,
            '16QAM': 4.0,
            '64QAM': 6.0
        }

        self.bandwidth_hz = self.data_rate_bps / self.bandwidth_efficiency[self.modulation]
        self.rx_sensitivity_w = self._calculate_rx_sensitivity()

    def _calculate_rx_sensitivity(self):
        k = 1.38e-23                                                            # Boltzmann constant (J/K)
        T = 290                                                                 # noise temperature (K)
        noise_floor_w = k * T * self.bandwidth_hz                               # calculate noise floor
        snr_linear = 10 ** (self.snr_requirements_db[self.modulation] / 10)     # dB to linear conversion
        return noise_floor_w * snr_linear

    def compute_fsp_loss(self):
        wavelength = self.c / self.frequency                                    # wavelength (m)
        fsp_loss = (wavelength / (4 * np.pi * self.distance)) ** 2              # free space loss
        atmospheric_loss_db = 8                                                 # atmospheric loss based on literature
        atmospheric_loss_linear = 10 ** (-atmospheric_loss_db / 10)             # dB to linear conversion
        total_loss = fsp_loss * atmospheric_loss_linear
        return total_loss

    def compute_required_rx_power(self):
        rx_sensitivity_w = self.rx_sensitivity_w
        link_margin_linear = 10 ** (self.target_link_margin_db / 10)            # dB to linear conversion
        return rx_sensitivity_w * link_margin_linear

    def compute_required_tx_power(self, tx_antenna_gain_db=35, rx_antenna_gain_db=10):                              # gains given in dB (based on literature)
        fsp_loss = self.compute_fsp_loss()
        tx_antenna_gain_linear = 10 ** (tx_antenna_gain_db / 10)                                                    # dB to linear conversion for transmitter gain
        rx_antenna_gain_linear = 10 ** (rx_antenna_gain_db / 10)                                                    # dB to linear conversion for receiver gain
        required_rx_power_w = self.compute_required_rx_power()
        required_tx_power_w = required_rx_power_w / (tx_antenna_gain_linear * rx_antenna_gain_linear * fsp_loss)    # solving link budget for required power
        return required_tx_power_w

    def __str__(self):
        return (f"LinkBudget(frequency={self.frequency/1e9:.2f} GHz, distance={self.distance} m, "f"data_rate={self.data_rate_bps/1e6:.2f} Mbit/s, modulation={self.modulation}, "f"bandwidth={self.bandwidth_hz/1e6:.2f} MHz, rx_sensitivity={self.rx_sensitivity_w:.2e} W)")


'''
# test
distance = 1000 * np.sqrt(18**2 + 400**2)  # example for max distance

lb = link_budget(distance=distance)

# print
print(lb)
print(f"Required Rx Power: {lb.compute_required_rx_power():.2e} W")
print(f"Required Tx Power: {lb.compute_required_tx_power():.2e} W")
'''


distance = 1000 * np.sqrt(18**2 + 400**2)  # example for max distance

lb = link_budget(distance=distance)
print(lb.compute_required_tx_power())
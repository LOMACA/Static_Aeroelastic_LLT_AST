##############################################################################################

# This class takes care of computing temperature, pressure, and density based on the ISA

##############################################################################################

class ISAComputation:

    """
    This class takes care of computing temperature, pressure, and density based on the ISA formulations
    """

    def __init__(self, altitude):
        self.altitude = altitude

    def compute(self):

        """
        This function calculates temperature, pressure, and density from the user input flight altitude
        """

        rho0 = 1.225
        P0 = 101325
        T0 = 288.15

        if 0 <= self.altitude <= 11:
            T = T0 - 6.5 * self.altitude
            P = P0 * (T / T0) ** 5.2561
        elif 11 < self.altitude <= 20:
            T = 216.65
            P = 22630.6 * 10 ** (-((self.altitude - 11) / 14.596))
        elif 20 < self.altitude < 32:
            T = 216.65 - (self.altitude - 20)
            P = 22630.6 * 10 ** - ((self.altitude - 11) / 14.596)

        rho = rho0 * (P / P0) * (T0 / T)

        return {
            'temperature': T,
            'pressure': P,
            'density': rho,
        }

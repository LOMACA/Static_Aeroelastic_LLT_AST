##############################################################################

# LLT calculation module

##############################################################################

import numpy as np

class LiftingLineTheory:

    """
    This class represents the LLT implementation following the description by Bertin (2021)
    """

    def __init__(self, v_inf, c_root, c_tip, b, n, cl0, cd0, a0, alpha0, rho):
        self.v_inf = v_inf
        self.c_root = c_root
        self.c_tip = c_tip
        self.b = b
        self.n = n
        self.cl0 = cl0
        self.cd0 = cd0
        self.a0 = a0
        self.alpha0 = alpha0 * np.pi / 180
        self.rho = rho
        self.y = None

    def compute(self, alpha_values):

        """
        This function computes the aerodynamic parameters and loads and stores
        them in a dictionary
        """

        results = []

        b_half = self.b / 2
        TR = self.c_tip / self.c_root
        Aw = b_half * (self.c_root + self.c_tip)
        S = 2 * Aw
        AR = self.b ** 2 / (2 * Aw)

        if self.n == 4:
            phi = np.array([22.5, 45, 67.5, 90])
        else:
            phi = np.linspace(0.1e-8, 90, self.n)

        phi_r = phi * np.pi / 180

        self.y = (self.b / 2) * np.cos(phi_r)

        c = self.c_root * (1 + (TR - 1) * np.cos(phi_r))
        geo = (c * self.a0) / (2 * AR * self.c_root * (1 + TR))

        alpha_results = np.zeros_like(self.y)

        for alpha in alpha_values:
            self.alpha_rad = alpha * np.pi / 180

            Left = geo * (self.alpha_rad - self.alpha0) * np.sin(phi_r)
            B = np.zeros((self.n, self.n))

            for i in range(self.n):
                for j in range(self.n):
                    B[i, j] = np.sin((2 * (j + 1) - 1) * phi_r[i]) * (geo[i] * (2 * (j + 1) - 1) + np.sin(phi_r[i]))

            A = np.linalg.solve(B, Left)

            CL = np.pi * AR * A[0]
            delta = np.sum([(2 * (i + 1) - 1) * A[i] ** 2 for i in range(1, self.n)])
            CD_ind = (CL ** 2 / (np.pi * AR)) * (1 + delta)

            CD = self.cd0 + CD_ind

            x = np.linspace(-self.b / 2, self.b / 2, self.n)
            D_ind = CD_ind * 0.5 * self.rho * self.v_inf ** 2 * Aw  # Induced Drag [N]
            D = CD * 0.5 * self.rho * self.v_inf ** 2 * Aw  # Total Drag [N] (neglecting wave drag and pressure as well as viscous drag)
            gamma_temp = np.zeros((self.n))
            for i in range(0, self.n):
                for j in range(0, self.n):
                    gamma_temp[i] = gamma_temp[i] + A[j] * np.sin((2 * (j + 1) - 1) * phi_r[i])

            gamma = 2 * self.b * self.v_inf * gamma_temp  # Circulation distribution [m^2/s]
            tau = 0.05
            aw = self.a0 / (1 + (self.a0 / (np.pi * AR)) * (1 + tau))  # lift slope of finite wing [rad]
            L_dist = self.rho * self.v_inf * gamma
            self.y_new = self.y[::-1]

            L_dist_N = np.zeros_like(L_dist)

            for i in range(len(self.y)):
                if i == 0:

                    L_dist_N[i] = L_dist[i] * (self.y[0] - self.y[1])
                else:
                    dy = self.y[i-1] - self.y[i]
                    L_dist_N[i] = L_dist[i] * dy

            L_dist_N[0] = L_dist[0] * (self.y[0] - self.y[1])

            spanwise_cl = []

            for i in range(len(c)):

                cl = (2*gamma[i]) / (c[i]*self.v_inf)
                spanwise_cl.append(cl)

            M_sum = np.zeros_like(L_dist)

            L_dist_new = L_dist[::-1]  # re-ordered lift distribution [N]
            for i in range(len(self.y)):
                if i == 0:
                    continue

                M_sum[i] = np.sum(L_dist_new[:i] * np.diff(self.y_new[:i + 1]))
            M_dist = M_sum  # bending moment distribution [Nm]

            temp = np.zeros(self.n)
            for i in range(0, self.n):
                for j in range(0, self.n):
                    temp[i] = temp[i] + ((2 * (j + 1) - 1)) * A[j] * np.sin((2 * (j + 1) - 1) * phi_r[i]) / np.sin(
                        phi_r[i])

            dw = -self.v_inf * temp  # downwash [m/s]
            alpha_i = dw / self.v_inf  # induced angle of attack [rad]

            D_dist = L_dist_N * np.abs(alpha_i)  # Induced Drag [N]

            results.append({
                'alpha': alpha,
                'Total CL': CL,
                'CD': CD,
                'Total Lift': np.sum(L_dist_N)*2,
                'Drag': D,
                'Circulation': gamma,
                'Bending Moment': M_dist,
                'Induced Drag': D_dist,
                'L_dist': L_dist,
                'Lift': L_dist_N,
                'Lift Coefficient': spanwise_cl,
                'Downwash': dw,
                'y': self.y,
                'aw': aw,
                'AR': AR,
                'TR': TR,
                'S': S,
                'v_inf': self.v_inf,
                'c_root': self.c_root,
                'c_tip': self.c_tip,
                'b': self.b,
                'n': self.n
            })

        return results






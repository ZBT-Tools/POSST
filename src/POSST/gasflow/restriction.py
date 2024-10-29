from src.POSST.misc import Gas

class model_1():
    """Calculate Massflow"""

    def __init__(self, name="Massflow", l=0.1, b=1, h=1, R=8.314, m=0, dt=1e-3):
        """
        Annahme Viskoität & Dichte immer kons
        """
        self.b = b
        self.h = h
        self.l = l
        self.U = 2 * b + 2 * h
        self.A = b * h
        self.D_h = 4 * self.A / self.U
        self.m = m
        self.R = R
        self.Re = 1
        self.dt = dt
        self.Gas = Gas.path
        self.CNat = Gas.CNat
        self.v = 1e-13

    def __call__(self, p_in, p_out, T_in, x_O2=0.23, n_N2=0.77, x_H2=0, x_H2O=0):
        # try:
        rho_t = 1.2
        vis_t = 13.3 * 1e-6
        self.p_vor = p_in
        self.p_nach = p_out
        self.T_vor = T_in
        # todo: rho_gemisch
        rho_i = []
        xn_i = []
        nO2_in = self.m * x_O2 / Gas.CNat.M_O2
        nN2_in = self.m * n_N2 / Gas.CNat.M_N2
        nH2_in = self.m * x_H2 / Gas.CNat.M_H2
        nH2O_in = self.m * x_H2O / Gas.CNat.M_H2O
        if self.m != 0:
            rho_i += Gas.path.rho(self, rho_n=Gas.CNat.rho_O2, p=self.p_vor, p_n=Gas.CNat.p_atm, T_in=self.T_vor, T_n=Gas.CNat.T_0K),
            xn_i += nO2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += Gas.path.rho(self, rho_n=Gas.CNat.rho_H2, p=self.p_vor, p_n=Gas.CNat.p_atm, T_in=self.T_vor, T_n=Gas.CNat.T_0K),
            xn_i += nH2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += Gas.path.rho(self, rho_n=Gas.CNat.rho_H2O, p=self.p_vor, p_n=Gas.CNat.p_atm, T_in=self.T_vor, T_n=Gas.CNat.T_0K),
            xn_i += nH2O_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += Gas.path.rho(self, rho_n=Gas.CNat.rho_N2, p=self.p_vor, p_n=Gas.CNat.p_atm, T_in=self.T_vor, T_n=Gas.CNat.T_0K),
            xn_i += nN2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            self.rho = Gas.path.rho_mix(self, rho_i=rho_i, xn_i=xn_i)
        else:
            self.rho = 1.2
        # self.v = self.Gas.v(self, m_in=max(self.m, 1e-13), rho=self.rho, w=self.b, h=self.h)
        # self.nu = self.Gas.nu(self, nu_n=self.CNat.nu_n_air, T_in=self.T_vor, T_n=self.CNat.T_0K)
        # self.Re = self.Gas.Re(self, rho=self.rho, v=self.v, D_H=self.D_h, nu=self.nu)
        self.Re = 2000  # (self.A/self.U*rho_t*(self.m/rho_t*self.A))/vis_t

        if self.Re != 0:
            c_d = 64 / self.Re
        else:
            c_d = 64 / 1000
        self.delta_P = self.p_vor - self.p_nach
        if abs(self.delta_P) > 1:
            # self.m=c_d*self.A  \
            #        *math.sqrt(abs((2*self.p_vor)/(self.R*(self.T_vor))))  \
            #        *math.sqrt(abs(self.delta_P))  \
            #        *self.delta_P/abs(self.delta_P)#p_vor*m_vor/p_nach
            self.v = self.Gas.v_2(self, p_in=self.p_vor, p_out=self.p_nach, rho=self.rho, l=self.l, D_H=self.D_h, Re=self.Re, zetha=0.4)
            self.m = self.rho * self.A * self.v * self.dt
        else:
            self.m = 0
        None

class model_2():
    def __init__(self, dt, A, R=287, cp=1000):
        self.dt=dt
        self.A=A
        self.R=R
        self.cp = cp
        self.cv = cp - R
        self.kappa = self.cp / self.cv
        self.v_sound=1000
        self.m = 0

    def __call__(self, p_in, p_out, T_in, coef_flow=1):
        dp = p_in - p_out

        if dp != 0:
            v = abs(self.kappa / (self.kappa - 1)
                    * (1 - (p_out / p_in) ** ((self.kappa - 1) / self.kappa))
                    ) ** 0.5

            psi = abs(self.kappa / (self.kappa - 1) *
                      ((p_out / p_in) ** (2 / self.kappa) - (p_out / p_in) ** ((self.kappa + 1) / self.kappa))) ** 0.5

            Ma = ((2 * self.R * T_in) ** 0.5 * v) / (self.kappa * self.R * T_in) ** 0.5

            if Ma > 1:
                v = (self.kappa / (self.kappa + 1)) ** 0.5

                psi = (2 / (self.kappa + 1)) ** (1 / (self.kappa - 1)) \
                      * abs(self.kappa / (self.kappa + 1)) ** 0.5

            self.u = (2 * self.R * T_in) ** 0.5 * v * (dp) / abs(dp)
            self.m = coef_flow*self.A*p_in*abs(2/(self.R*T_in))**0.5 * psi * (dp)/abs(dp)
        else:
            self.m = 0
            self.u = 0
        return self.m
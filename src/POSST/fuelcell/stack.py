
from POSST.misc import Gas
import numpy as np
import math

class model():
    r"""
    ## Description  
    The fuel cell model is divided into various sub-models to calculate flow fields, cell voltage, and thermal management. Based on the voltage model described by Gößling in a 2D-1D model [1], a reduced model is implemented that omits several key features.
    There is no segmentation along the gas paths, and the membrane does not serve as a water reservoir. The membrane humidity relies solely on the average relative humidity of the medium in the flow fields.

    ## Reference
    [1]	S. Gößling, "2-D + 1-D ortsaufgelöste Modellierung von PEM-Brennstoffzellen," Fakultät für Ingeneuirwissenschaften, Abteilung Maschienenbau und Verfahrenstechnik, Universität Duisburg Essen, Duisburg, 2019.<br>

    """
    def __init__(self, P=None, dt=1e-3):
        r"""
        ## Arguments
            P: Parameter dict, if not provided the default parameters are used. The default params are provided below.
            dt: Time step size in seconds.

        ## Returns
            Returns an instance of the stack model.
        """
        self.dt = dt
        self.P = {}
        self.P |= {'NumCell': 70}
        self.P |= {'Active_area': 190/1e4}
        self.P |= {'Channel': {'Cathode': {'Num': 20, 'h': 0.001, 'b': 0.002, 'l': 0.1},
                               'Anode': {'Num': 20, 'h': 0.001, 'b': 0.002, 'l': 0.1}},
                   'Stack': {'Thermo': {'M': 10000, 'cp': 420, 'h':200}}}
        if P != None:
            self.P |= P

        self.P['Channel']['Cathode'] |= {'A': self.P['Channel']['Cathode']['h'] *
                                              self.P['Channel']['Cathode']['b']}
        self.P['Channel']['Cathode'] |= {'U': 2 * self.P['Channel']['Cathode']['h'] +
                                              2 * self.P['Channel']['Cathode']['b']}
        self.P['Channel']['Cathode'] |= {'V_ges': self.P['Channel']['Cathode']['Num'] *
                                                  self.P['Channel']['Cathode']['h'] *
                                                  self.P['Channel']['Cathode']['b'] *
                                                  self.P['Channel']['Cathode']['l']}
        self.P['Channel']['Anode'] |= {'A': self.P['Channel']['Anode']['h'] *
                                            self.P['Channel']['Anode']['b']}
        self.P['Channel']['Anode'] |= {'U': 2 * self.P['Channel']['Anode']['h'] +
                                            2 * self.P['Channel']['Anode']['b']}
        self.P['Channel']['Anode'] |= {'V_ges': self.P['Channel']['Anode']['Num'] *
                                                self.P['Channel']['Anode']['h'] *
                                                self.P['Channel']['Anode']['b'] *
                                                self.P['Channel']['Anode']['l']}


        self._m_calc = Gas.Massflow(self.P['NumCell'])

        self.pTmx_CAo = {'p_out': 1.0135e5, 'p_in': 1.0135e5, 'T_out': 293.15, 'm_out': 0,
                         'x_H2': 0., 'x_H2O': 0., 'x_O2': 0.23, 'x_N2': 0.77}
        self.pTmx_ANo = {'p_out': 1.0135e5, 'p_in': 1.0135e5, 'T_out': 293.15, 'm_out': 0,
                         'x_H2': 1., 'x_H2O': 0., 'x_O2': 0., 'x_N2': 0.}
        self.pTmx_CAi = {'p_out': 1.0135e5, 'p_in': 1.0135e5, 'T_in': 293.15, 'm_in': 0,
                         'x_H2': 0., 'x_H2O': 0., 'x_O2': 0.23, 'x_N2': 0.77}
        self.pTmx_ANi = {'p_out': 1.0135e5, 'p_in': 1.0135e5, 'T_in': 293.15, 'm_in': 0,
                         'x_H2': 1., 'x_H2O': 0., 'x_O2': 0., 'x_N2': 0.}

        self.pTm_COi = {'p_out': 1.0135e5, 'p_in': 1.0135e5, 'T_out': 293.15, 'T_in': 293.15, 'm_out': 0, 'm_in': 0}
        self.pTm_COo = {'p_out': 1.0135e5, 'p_in': 1.0135e5, 'T_out': 293.15, 'T_in': 293.15, 'm_out': 0, 'm_in': 0}
        self._gas_AN = Gas.gas_state()
        self._gas_CA = Gas.gas_state()

        self._path_func = Gas.path()
        self._Cnat = Gas.CNat()

        self._U_0 = 1.23
        self.U_conc = 0
        self.U_ohm = 0
        self.U_act = 0
        self.T_cell=self.T_stack=333.15
        self.U_stack= self._U_0 * self.P['NumCell']
        self.I_stack =0
        self.U_loss_stack = self._U_0 * self.P['NumCell'] - self.U_stack

        self.Qdot_stack=0
        self.Qdot_cool=0




    def __call__(self, I, pTmx_ANi, pTmx_CAi, pTm_COi):
        r"""
        ## Arguments
            I: Stack current in [A]
            pTmx_ANi: Anode inlet pTmx vector
            pTmx_CAi: Cathode inlet pTmx vector
            pTm_COi: Cooling inlet pTm vector

            The inlet pTm(x) vector is a dict with the keys [p_out, T_in, m_in, x_H2, x_H2O, x_O2, x_N2].
            p_out: is the outlet pressure in [Pa]
            T_in: is the inlet temperature in [K]
            m_in: is the inlet mass flow in [Kg/s]
            x_H2: is the inlet H2 mass fraction in [-]
            x_H2O: is the inlet H2O mass fraction in [-]
            x_O2: is the inlet O2 mass fraction in [-]
            x_N2: is the inlet N2 mass fraction in [-]
            The coolant is considered 100% H2O so no mass fractions need to be included.

        ## Returns
            No returns. The needed Values are taken directly from self.

        ## Usage
            called one per timestep to calculate the Fuel cell stack state.

        """
        self.pTmx_ANi = pTmx_ANi
        self.pTmx_CAi = pTmx_CAi
        self.pTm_COi = self.pTm_COo = pTm_COi
        self.I_stack = I

        if self.pTmx_ANi['m_in'] > 0:
            self._path_AN()
        if self.pTmx_CAi['m_in'] > 0:
            self._path_Cath()

        if self._m_calc.St_H2(mH2=self.pTmx_ANi['m_in'] * 1000 * self.pTmx_ANi['x_H2'], I=I) < 1:
            # print("Error: AN_stoich= " + str(self._m_calc.St_H2(mH2=self.pTmx_ANi['m_in']*1000 * self.pTmx_ANi['x_H2'], I=I)))
            None

        elif self._m_calc.St_O2(mO2=self.pTmx_CAi['m_in'] * 1000 * self.pTmx_CAi['x_O2'], I=I) < 1:
            # print("Error: Ca_stoich= " + str(self._m_calc.St_O2(mO2=self.pTmx_CAi['m_in']*1000 * self.pTmx_CAi['x_O2'], I=I)))
            None
        else:
            self.do_step()

    def do_step(self):
        r"""
        ## Description
            Actual calculations of the Fuel cell stack state, without plausibility checks beforehand.
        """
        self.T_cell = self.T_stack
        self.U_stack = self._U_Cell() * self.P['NumCell']
        self._Thermo()


    def _path_AN(self):
        r"""
        ## Description
            Calculation of the Anode gas path.
            Result is the Anode outlet pTmx vector.
        """
        self.c_H2 = (self.pTmx_ANi['m_in'] * self.pTmx_ANi['x_H2'] / self._Cnat.M_H2) / self.P['Channel']['Anode']['V_ges']

        rho_i = []
        xn_i = []
        if self.pTmx_ANi['m_in'] != 0:
            nO2_in = self.pTmx_ANi['m_in'] * self.pTmx_ANi['x_O2'] / self._Cnat.M_O2
            nN2_in = self.pTmx_ANi['m_in'] * self.pTmx_ANi['x_N2'] / self._Cnat.M_N2
            nH2_in = self.pTmx_ANi['m_in'] * self.pTmx_ANi['x_H2'] / self._Cnat.M_H2
            nH2O_in = self.pTmx_ANi['m_in'] * self.pTmx_ANi['x_H2O'] / self._Cnat.M_H2O

            rho_i += self._path_func.rho(rho_n=self._Cnat.rho_O2, p=self.pTmx_ANi['p_out'], p_n=self._Cnat.p_atm, T_in=self.pTmx_ANi['T_in'], T_n=self._Cnat.T_0K),
            xn_i += nO2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self._path_func.rho(rho_n=self._Cnat.rho_H2, p=self.pTmx_ANi['p_out'], p_n=self._Cnat.p_atm, T_in=self.pTmx_ANi['T_in'], T_n=self._Cnat.T_0K),
            xn_i += nH2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self._path_func.rho(rho_n=self._Cnat.rho_H2O, p=self.pTmx_ANi['p_out'], p_n=self._Cnat.p_atm, T_in=self.pTmx_ANi['T_in'], T_n=self._Cnat.T_0K),
            xn_i += nH2O_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self._path_func.rho(rho_n=self._Cnat.rho_N2, p=self.pTmx_ANi['p_out'], p_n=self._Cnat.p_atm, T_in=self.pTmx_ANi['T_in'], T_n=self._Cnat.T_0K),
            xn_i += nN2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            self._gas_AN.rho = self._path_func.rho_mix(rho_i=rho_i, xn_i=xn_i)
        else:
            self._gas_AN.rho = 1
        self._gas_AN.v = self._path_func.v(
                               m_in=self.pTmx_ANi['m_in'] / self.P['Channel']['Anode']['Num'] / self.P['NumCell'],
                               rho=self._gas_AN.rho,
                               w=self.P['Channel']['Anode']['b'],
                               h=self.P['Channel']['Anode']['h'])
        self._gas_AN.nu = self._path_func.nu(
                                 nu_n=self._Cnat.nu_n_H2,
                                 T_in=self.pTmx_ANi['T_in'], T_n=self._Cnat.T_0K)
        self._gas_AN.D_H = self._path_func.D_H(
                                   w=self.P['Channel']['Anode']['b'],
                                   h=self.P['Channel']['Anode']['h'])
        self._gas_AN.Re = self._path_func.Re(
                                 rho=self._gas_AN.rho,
                                 v=self._gas_AN.v,
                                 D_H=self._gas_AN.D_H,
                                 nu=self._gas_AN.nu)
        if self._gas_AN.Re != 0:
            self._gas_AN.dp = self._path_func.dp(
                                     rho=self._gas_AN.rho,
                                     v=self._gas_AN.v,
                                     l_ch=self.P['Channel']['Anode']['l'],
                                     D_H=self._gas_AN.D_H,
                                     Re=self._gas_AN.Re,
                                     zeta=0.45)
        self.pTmx_ANo['p_in'] = self._gas_AN.dp + self.pTmx_ANi['p_out']

        self.pTmx_ANo['T_out'] = self.T_cell

        m_H2_out = self._path_func.mdot_reac_out(m_in=self.pTmx_ANi['m_in'], x=self.pTmx_ANi['x_H2'], I=self.I_stack, M=self._Cnat.M_H2, Ze=self._Cnat.Ze_H2, Num_Cell=self.P['NumCell'])
        m_H2O_out = self.pTmx_ANi['x_H2O'] * self.pTmx_ANi['m_in'] #todo. Wasserdiffusion
        self.pTmx_ANo['m_out'] = m_H2_out + m_H2O_out
        self.pTmx_ANo['x_H2'] = m_H2_out / self.pTmx_ANo['m_out']
        self.pTmx_ANo['x_H2O'] = m_H2O_out / self.pTmx_ANo['m_out']
        self.pTmx_ANo['x_N2'] = 1 - self.pTmx_ANo['x_H2'] - self.pTmx_ANo['x_H2O']
        self.pTmx_ANo['x_O2'] = 0

    def _path_Cath(self):
        r"""
        ## Description
            Calculation of the Cathode gas path.
            Result is the Cathode outlet pTmx Vector.
        """
        self.c_O2 = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_O2'] / self._Cnat.M_O2 / self.P['Channel']['Cathode']['V_ges']
        rho_i = []
        xn_i = []
        if self.pTmx_CAi['m_in'] != 0:
            nO2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_O2'] / self._Cnat.M_O2
            nN2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_N2'] / self._Cnat.M_N2
            nH2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_H2'] / self._Cnat.M_H2
            nH2O_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_H2O'] / self._Cnat.M_H2O

            rho_i += self._path_func.rho(rho_n=self._Cnat.rho_O2, p=self.pTmx_CAi['p_out'], p_n=self._Cnat.p_atm, T_in=self.pTmx_CAi['T_in'], T_n=self._Cnat.T_0K),
            xn_i += nO2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self._path_func.rho(rho_n=self._Cnat.rho_H2, p=self.pTmx_CAi['p_out'], p_n=self._Cnat.p_atm, T_in=self.pTmx_CAi['T_in'], T_n=self._Cnat.T_0K),
            xn_i += nH2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self._path_func.rho(rho_n=self._Cnat.rho_H2O, p=self.pTmx_CAi['p_out'], p_n=self._Cnat.p_atm, T_in=self.pTmx_CAi['T_in'], T_n=self._Cnat.T_0K),
            xn_i += nH2O_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self._path_func.rho(rho_n=self._Cnat.rho_N2, p=self.pTmx_CAi['p_out'], p_n=self._Cnat.p_atm, T_in=self.pTmx_CAi['T_in'], T_n=self._Cnat.T_0K),
            xn_i += nN2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            self._gas_CA.rho = self._path_func.rho_mix(rho_i=rho_i, xn_i=xn_i)
        else:
            self._gas_CA.rho = 1

        self._gas_CA.v = self._path_func.v(
                               m_in=self.pTmx_CAi['m_in'] / self.P['Channel']['Cathode']['Num'] / self.P['NumCell'],
                               rho=self._gas_CA.rho,
                               w=self.P['Channel']['Cathode']['b'],
                               h=self.P['Channel']['Cathode']['h'])
        self._gas_CA.nu = self._path_func.nu(
                                 nu_n=self._Cnat.nu_n_air,
                                 T_in=self.pTmx_CAi['T_in'], T_n=self._Cnat.T_0K)
        self._gas_CA.D_H = self._path_func.D_H(
                                   w=self.P['Channel']['Cathode']['b'],
                                   h=self.P['Channel']['Cathode']['h'])
        self._gas_CA.Re = self._path_func.Re(
                                 rho=self._gas_CA.rho,
                                 v=self._gas_CA.v,
                                 D_H=self._gas_CA.D_H,
                                 nu=self._gas_CA.nu)
        if self._gas_CA.Re != 0:
            self._gas_CA.dp = self._path_func.dp(
                                     rho=self._gas_CA.rho,
                                     v=self._gas_CA.v,
                                     l_ch=self.P['Channel']['Cathode']['l'],
                                     D_H=self._gas_CA.D_H,
                                     Re=self._gas_CA.Re,
                                     zeta=0.45)
        self.pTmx_CAo['p_in'] = self._gas_CA.dp + self.pTmx_CAi['p_out']

        self.pTmx_CAo['T_out'] = self.T_cell

        m_O2_out = self._path_func.mdot_reac_out(m_in=self.pTmx_CAi['m_in'], x=self.pTmx_CAi['x_O2'], I=self.I_stack, M=self._Cnat.M_O2, Ze=self._Cnat.Ze_O2, Num_Cell=self.P['NumCell'])
        m_H2O_out = self._path_func.m_H2O_prod(I=self.I_stack, Num_Cell=self.P['NumCell']) + self.pTmx_CAi['x_H2O'] * self.pTmx_CAi['m_in'] #todo. Wasserdiffusion
        self.pTmx_CAo['m_out'] = self.pTmx_CAi['m_in']*self.pTmx_CAi['x_N2'] + m_O2_out + m_H2O_out
        self.pTmx_CAo['x_O2'] = m_O2_out / self.pTmx_CAo['m_out']
        self.pTmx_CAo['x_H2O'] = m_H2O_out / self.pTmx_CAo['m_out']
        self.pTmx_CAo['x_N2'] = 1-self.pTmx_CAo['x_O2']-self.pTmx_CAo['x_H2O']
        self.pTmx_CAo['x_H2'] = 0

    def _U_Cell(self):
        r"""
        ## Description
            Voltage calculations.
        """
        A = 0.02  #0.03 V
        C_i0 = .8 #0.8
        i_cross = 1E-8  #1e-8 A/m^2
        C_m = 7e-2 #7E-3
        C_n = 3 #0.9
        C_ohm = 2 #1
        DM= 2e-5

        nO2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_O2'] / self._Cnat.M_O2
        nN2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_N2'] / self._Cnat.M_N2
        nH2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_H2'] / self._Cnat.M_H2
        nH2O_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_H2O'] / self._Cnat.M_H2O
        xnO2_in = nO2_in / (nO2_in + nN2_in + nH2_in + nH2O_in)

        nO2_out = self.pTmx_CAo['m_out'] * self.pTmx_CAo['x_O2'] / self._Cnat.M_O2
        nN2_out = self.pTmx_CAo['m_out'] * self.pTmx_CAo['x_N2'] / self._Cnat.M_N2
        nH2_out = self.pTmx_CAo['m_out'] * self.pTmx_CAo['x_H2'] / self._Cnat.M_H2
        nH2O_out = self.pTmx_CAo['m_out'] * self.pTmx_CAo['x_H2O'] / self._Cnat.M_H2O
        xnO2_out = nO2_out / (nO2_out + nN2_out + nH2_out + nH2O_out)

        I_cd = self.I_stack/(self.P['Active_area'])
        p_mean = (self.pTmx_CAo['p_in']+self.pTmx_CAo['p_out'])/2 /1e5
        xnO2_mean = (xnO2_in + xnO2_out) / 2
        rH_in = self._m_calc.rH(p=self.pTmx_CAo['p_in'] / 1e5,
                                T=self.pTmx_CAi['T_in'],
                                m=self.pTmx_CAi['m_in'],
                                H2=self.pTmx_CAi['x_H2'],
                                H2O=self.pTmx_CAi['x_H2O'],
                                O2=self.pTmx_CAi['x_O2'],
                                N2=self.pTmx_CAi['x_N2'])
        rH_out = self._m_calc.rH(p=self.pTmx_CAi['p_out'] / 1e5,
                                 T=self.pTmx_CAo['T_out'],
                                 m=self.pTmx_CAo['m_out'],
                                 H2=self.pTmx_CAo['x_H2'],
                                 H2O=self.pTmx_CAo['x_H2O'],
                                 O2=self.pTmx_CAo['x_O2'],
                                 N2=self.pTmx_CAo['x_N2'])

        self.U_act = max(0, A * np.log((I_cd + i_cross) / (C_i0 * p_mean * xnO2_mean * 2E-2 - 0.0015)))  # 2E-7 * 1E+5 wegen bar/ Pa

        self.U_ohm = (I_cd + i_cross) * DM / self.__Membranleitwert(rH=(rH_out + rH_in) / 200) * C_ohm
        n_emp = C_n * (0.0001340428 * np.log(p_mean / self.pTmx_CAo['p_in']*1e-5) + 0.000062629589) * (1.6574 - p_mean)
        self.U_conc = C_m * np.exp(n_emp * (I_cd + i_cross))

        return self._U_0 - self.U_act - self.U_ohm - self.U_conc

    def __Membranleitwert(self, rH):
        if self.__Lambda_calc(rH) < 1.253:
            Leitwert = 1e-13
        else:
            Leitwert = (0.005738 * self.__Lambda_calc(rH) - 0.007192) * 100

            Leitwert = max(Leitwert, 0.2)  # for numeric stability

        return Leitwert  # S/m oder 1/Ohm*m

    def __Lambda_calc(self, rH):
        # Kulikovsky, Andrei A.: Chapter 1 - Fuel cell basics. In: Kulikovsky, Andrei A. (Herausgeber): Analytical Modelling of Fuel Cells, Seiten 1 – 38. Elsevier, Amsterdam, 2010

        Lambda = 0.3 + 6 * rH * (1 - math.tanh(rH - 0.5)) + 3.9 * rH ** 0.5 * (1 + math.tanh((rH - 0.89) / 0.23))
        return Lambda

    def _Thermo(self):
        r"""
        ## Description
            The waste heat and stack temperatures are calculated based on stack voltage and stack current. Using the provided coolant mass flow, the resulting temperature change for both the stack and the coolant is determined by

            $$dT = \frac{Q_{\text{waste}} - Q_{\text{coolant}}}{M \cdot c} \cdot dt$$

            where $ Q_{\text{waste}} $ is the waste heat, $ Q_{\text{coolant}} $ is the provided coolant flow, $ M $ is the thermal mass, and $ c $ is the heat capacity of the fuel cell stack. A pressure drop is not calculated by thermal management. Additionally, the outlet temperatures of the gas flows are set to the stack temperature. This temperature change is not yet accounted for in thermal management.
        """
        if self.pTm_COi['m_in'] != 0:
            self.U_loss_stack = self._U_0 * self.P['NumCell'] - self.U_stack
            self.Qdot_stack = self.U_loss_stack * self.I_stack
            cp_water=4200 #J/KgK
            self.Qdot_cool = (self.P['Active_area'] * ((self.P['NumCell'] - 1) * 2 + 2)) * self.P['Stack']['Thermo']['h'] * (self.pTm_COi['T_in'] - self.T_stack)
            # self.Qdot_cool = (self.pTm_COi['m_in'] * cp_water) * (self.T_stack - self.pTm_COi['T_in'])
            self.T_stack += (self.Qdot_stack + self.Qdot_cool)/(self.P['Stack']['Thermo']['M'] * self.P['Stack']['Thermo']['cp']) * self.dt
            self.pTm_COo['T_out'] = self.pTm_COi['T_in'] - self.Qdot_cool/(self.pTm_COi['m_in'] * cp_water)
            # self.pTm_COo['T_out'] = self.Qdot_cool/((self.P['Active_area'] * ((self.P['NumCell'] - 1) * 2 + 2)) * self.P['Stack']['Thermo']['h']) + self.T_stack

            self.pTm_COo['m_out'] = self.pTm_COi['m_in']
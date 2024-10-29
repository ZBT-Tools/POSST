from src.POSST.misc import Gas
import numpy as np
import math

class model():
    def __init__(self, P=None, dt=1e-3):
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


        self.m_calc = Gas.Massflow(self.P['NumCell'])

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
        self.gas_AN = Gas.gas_state()
        self.gas_CA = Gas.gas_state()

        self.path_func = Gas.path()
        self.Cnat = Gas.CNat()

        self.U_0 = 1.23
        self.U_conc = 0
        self.U_ohm = 0
        self.U_act = 0
        self.T_cell=self.T_stack=333.15
        self.U_stack=self.U_0*self.P['NumCell']
        self.I_stack =0
        self.U_loss_stack = self.U_0*self.P['NumCell'] - self.U_stack

        self.Qdot_stack=0
        self.Qdot_cool=0




    def __call__(self, I, pTmx_ANi, pTmx_CAi, pTm_COi):
        self.pTmx_ANi = pTmx_ANi
        self.pTmx_CAi = pTmx_CAi
        self.pTm_COi = self.pTm_COo = pTm_COi
        self.I_stack = I

        if self.pTmx_ANi['m_in'] > 0:
            self.path_AN()
        if self.pTmx_CAi['m_in'] > 0:
            self.path_Cath()

        if self.m_calc.St_H2(mH2=self.pTmx_ANi['m_in']*1000 * self.pTmx_ANi['x_H2'], I=I) < 1:
            # print("Error: AN_stoich= " + str(self.m_calc.St_H2(mH2=self.pTmx_ANi['m_in']*1000 * self.pTmx_ANi['x_H2'], I=I)))
            None

        elif self.m_calc.St_O2(mO2=self.pTmx_CAi['m_in']*1000 * self.pTmx_CAi['x_O2'], I=I) < 1:
            # print("Error: Ca_stoich= " + str(self.m_calc.St_O2(mO2=self.pTmx_CAi['m_in']*1000 * self.pTmx_CAi['x_O2'], I=I)))
            None
        else:
            self.do_step()

    def do_step(self):
        self.T_cell = self.T_stack
        self.U_stack = self.U_Cell() * self.P['NumCell']
        self.Thermo()


    def path_AN(self):
        self.c_H2 = (self.pTmx_ANi['m_in'] * self.pTmx_ANi['x_H2'] / self.Cnat.M_H2) / self.P['Channel']['Anode']['V_ges']

        rho_i = []
        xn_i = []
        if self.pTmx_ANi['m_in'] != 0:
            nO2_in = self.pTmx_ANi['m_in'] * self.pTmx_ANi['x_O2'] / self.Cnat.M_O2
            nN2_in = self.pTmx_ANi['m_in'] * self.pTmx_ANi['x_N2'] / self.Cnat.M_N2
            nH2_in = self.pTmx_ANi['m_in'] * self.pTmx_ANi['x_H2'] / self.Cnat.M_H2
            nH2O_in = self.pTmx_ANi['m_in'] * self.pTmx_ANi['x_H2O'] / self.Cnat.M_H2O

            rho_i += self.path_func.rho(rho_n=self.Cnat.rho_O2, p=self.pTmx_ANi['p_out'], p_n=self.Cnat.p_atm, T_in=self.pTmx_ANi['T_in'], T_n=self.Cnat.T_0K),
            xn_i += nO2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self.path_func.rho(rho_n=self.Cnat.rho_H2, p=self.pTmx_ANi['p_out'], p_n=self.Cnat.p_atm, T_in=self.pTmx_ANi['T_in'], T_n=self.Cnat.T_0K),
            xn_i += nH2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self.path_func.rho(rho_n=self.Cnat.rho_H2O, p=self.pTmx_ANi['p_out'], p_n=self.Cnat.p_atm, T_in=self.pTmx_ANi['T_in'], T_n=self.Cnat.T_0K),
            xn_i += nH2O_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self.path_func.rho(rho_n=self.Cnat.rho_N2, p=self.pTmx_ANi['p_out'], p_n=self.Cnat.p_atm, T_in=self.pTmx_ANi['T_in'], T_n=self.Cnat.T_0K),
            xn_i += nN2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            self.gas_AN.rho = self.path_func.rho_mix(rho_i=rho_i, xn_i=xn_i)
        else:
            self.gas_AN.rho = 1
        self.gas_AN.v = self.path_func.v(
                               m_in=self.pTmx_ANi['m_in'] / self.P['Channel']['Anode']['Num'] / self.P['NumCell'],
                               rho=self.gas_AN.rho,
                               w=self.P['Channel']['Anode']['b'],
                               h=self.P['Channel']['Anode']['h'])
        self.gas_AN.nu = self.path_func.nu(
                                 nu_n=self.Cnat.nu_n_H2,
                                 T_in=self.pTmx_ANi['T_in'], T_n=self.Cnat.T_0K)
        self.gas_AN.D_H = self.path_func.D_H(
                                   w=self.P['Channel']['Anode']['b'],
                                   h=self.P['Channel']['Anode']['h'])
        self.gas_AN.Re = self.path_func.Re(
                                 rho=self.gas_AN.rho,
                                 v=self.gas_AN.v,
                                 D_H=self.gas_AN.D_H,
                                 nu=self.gas_AN.nu)
        if self.gas_AN.Re != 0:
            self.gas_AN.dp = self.path_func.dp(
                                     rho=self.gas_AN.rho,
                                     v=self.gas_AN.v,
                                     l_ch=self.P['Channel']['Anode']['l'],
                                     D_H=self.gas_AN.D_H,
                                     Re=self.gas_AN.Re,
                                     zeta=0.45)
        self.pTmx_ANo['p_in'] = self.gas_AN.dp + self.pTmx_ANi['p_out']

        self.pTmx_ANo['T_out'] = self.T_cell

        m_H2_out = self.path_func.mdot_reac_out(m_in=self.pTmx_ANi['m_in'], x=self.pTmx_ANi['x_H2'], I=self.I_stack, M=self.Cnat.M_H2, Ze=self.Cnat.Ze_H2, Num_Cell=self.P['NumCell'])
        m_H2O_out = self.pTmx_ANi['x_H2O'] * self.pTmx_ANi['m_in'] #todo. Wasserdiffusion
        self.pTmx_ANo['m_out'] = m_H2_out + m_H2O_out
        self.pTmx_ANo['x_H2'] = m_H2_out / self.pTmx_ANo['m_out']
        self.pTmx_ANo['x_H2O'] = m_H2O_out / self.pTmx_ANo['m_out']
        self.pTmx_ANo['x_N2'] = 1 - self.pTmx_ANo['x_H2'] - self.pTmx_ANo['x_H2O']
        self.pTmx_ANo['x_O2'] = 0

    def path_Cath(self):
        self.c_O2 = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_O2'] / self.Cnat.M_O2 / self.P['Channel']['Cathode']['V_ges']
        rho_i = []
        xn_i = []
        if self.pTmx_CAi['m_in'] != 0:
            nO2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_O2'] / self.Cnat.M_O2
            nN2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_N2'] / self.Cnat.M_N2
            nH2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_H2'] / self.Cnat.M_H2
            nH2O_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_H2O'] / self.Cnat.M_H2O

            rho_i += self.path_func.rho(rho_n=self.Cnat.rho_O2, p=self.pTmx_CAi['p_out'], p_n=self.Cnat.p_atm, T_in=self.pTmx_CAi['T_in'], T_n=self.Cnat.T_0K),
            xn_i += nO2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self.path_func.rho(rho_n=self.Cnat.rho_H2, p=self.pTmx_CAi['p_out'], p_n=self.Cnat.p_atm, T_in=self.pTmx_CAi['T_in'], T_n=self.Cnat.T_0K),
            xn_i += nH2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self.path_func.rho(rho_n=self.Cnat.rho_H2O, p=self.pTmx_CAi['p_out'], p_n=self.Cnat.p_atm, T_in=self.pTmx_CAi['T_in'], T_n=self.Cnat.T_0K),
            xn_i += nH2O_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            rho_i += self.path_func.rho(rho_n=self.Cnat.rho_N2, p=self.pTmx_CAi['p_out'], p_n=self.Cnat.p_atm, T_in=self.pTmx_CAi['T_in'], T_n=self.Cnat.T_0K),
            xn_i += nN2_in / (nO2_in + nN2_in + nH2_in + nH2O_in),

            self.gas_CA.rho = self.path_func.rho_mix(rho_i=rho_i, xn_i=xn_i)
        else:
            self.gas_CA.rho = 1

        self.gas_CA.v = self.path_func.v(
                               m_in=self.pTmx_CAi['m_in'] / self.P['Channel']['Cathode']['Num'] / self.P['NumCell'],
                               rho=self.gas_CA.rho,
                               w=self.P['Channel']['Cathode']['b'],
                               h=self.P['Channel']['Cathode']['h'])
        self.gas_CA.nu = self.path_func.nu(
                                 nu_n=self.Cnat.nu_n_air,
                                 T_in=self.pTmx_CAi['T_in'], T_n=self.Cnat.T_0K)
        self.gas_CA.D_H = self.path_func.D_H(
                                   w=self.P['Channel']['Cathode']['b'],
                                   h=self.P['Channel']['Cathode']['h'])
        self.gas_CA.Re = self.path_func.Re(
                                 rho=self.gas_CA.rho,
                                 v=self.gas_CA.v,
                                 D_H=self.gas_CA.D_H,
                                 nu=self.gas_CA.nu)
        if self.gas_CA.Re != 0:
            self.gas_CA.dp = self.path_func.dp(
                                     rho=self.gas_CA.rho,
                                     v=self.gas_CA.v,
                                     l_ch=self.P['Channel']['Cathode']['l'],
                                     D_H=self.gas_CA.D_H,
                                     Re=self.gas_CA.Re,
                                     zeta=0.45)
        self.pTmx_CAo['p_in'] = self.gas_CA.dp + self.pTmx_CAi['p_out']

        self.pTmx_CAo['T_out'] = self.T_cell

        m_O2_out = self.path_func.mdot_reac_out(m_in=self.pTmx_CAi['m_in'], x=self.pTmx_CAi['x_O2'], I=self.I_stack, M=self.Cnat.M_O2, Ze=self.Cnat.Ze_O2, Num_Cell=self.P['NumCell'])
        m_H2O_out = self.path_func.m_H2O_prod(I=self.I_stack, Num_Cell=self.P['NumCell']) + self.pTmx_CAi['x_H2O'] * self.pTmx_CAi['m_in'] #todo. Wasserdiffusion
        self.pTmx_CAo['m_out'] = self.pTmx_CAi['m_in']*self.pTmx_CAi['x_N2'] + m_O2_out + m_H2O_out
        self.pTmx_CAo['x_O2'] = m_O2_out / self.pTmx_CAo['m_out']
        self.pTmx_CAo['x_H2O'] = m_H2O_out / self.pTmx_CAo['m_out']
        self.pTmx_CAo['x_N2'] = 1-self.pTmx_CAo['x_O2']-self.pTmx_CAo['x_H2O']
        self.pTmx_CAo['x_H2'] = 0

    def U_Cell(self):
        A = 0.02  #0.03 V
        C_i0 = .8 #0.8
        i_cross = 1E-8  #1e-8 A/m^2
        C_m = 7e-2 #7E-3
        C_n = 3 #0.9
        C_ohm = 2 #1
        DM= 2e-5

        nO2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_O2'] / self.Cnat.M_O2
        nN2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_N2'] / self.Cnat.M_N2
        nH2_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_H2'] / self.Cnat.M_H2
        nH2O_in = self.pTmx_CAi['m_in'] * self.pTmx_CAi['x_H2O'] / self.Cnat.M_H2O
        xnO2_in = nO2_in / (nO2_in + nN2_in + nH2_in + nH2O_in)

        nO2_out = self.pTmx_CAo['m_out'] * self.pTmx_CAo['x_O2'] / self.Cnat.M_O2
        nN2_out = self.pTmx_CAo['m_out'] * self.pTmx_CAo['x_N2'] / self.Cnat.M_N2
        nH2_out = self.pTmx_CAo['m_out'] * self.pTmx_CAo['x_H2'] / self.Cnat.M_H2
        nH2O_out = self.pTmx_CAo['m_out'] * self.pTmx_CAo['x_H2O'] / self.Cnat.M_H2O
        xnO2_out = nO2_out / (nO2_out + nN2_out + nH2_out + nH2O_out)

        I_cd = self.I_stack/(self.P['Active_area'])
        p_mean = (self.pTmx_CAo['p_in']+self.pTmx_CAo['p_out'])/2 /1e5
        xnO2_mean = (xnO2_in + xnO2_out) / 2
        rH_in = self.m_calc.rH( p=self.pTmx_CAo['p_in'] /1e5,
                                T=self.pTmx_CAi['T_in'],
                                m=self.pTmx_CAi['m_in'],
                                H2=self.pTmx_CAi['x_H2'],
                                H2O=self.pTmx_CAi['x_H2O'],
                                O2=self.pTmx_CAi['x_O2'],
                                N2=self.pTmx_CAi['x_N2'])
        rH_out = self.m_calc.rH(p=self.pTmx_CAi['p_out'] /1e5,
                               T=self.pTmx_CAo['T_out'],
                               m=self.pTmx_CAo['m_out'],
                               H2=self.pTmx_CAo['x_H2'],
                               H2O=self.pTmx_CAo['x_H2O'],
                               O2=self.pTmx_CAo['x_O2'],
                               N2=self.pTmx_CAo['x_N2'])

        self.U_act = max(0, A * np.log((I_cd + i_cross) / (C_i0 * p_mean * xnO2_mean * 2E-2 - 0.0015)))  # 2E-7 * 1E+5 wegen bar/ Pa

        self.U_ohm = (I_cd + i_cross) * DM / self.Membranleitwert(rH=(rH_out+rH_in)/200) * C_ohm
        n_emp = C_n * (0.0001340428 * np.log(p_mean / self.pTmx_CAo['p_in']*1e-5) + 0.000062629589) * (1.6574 - p_mean)
        self.U_conc = C_m * np.exp(n_emp * (I_cd + i_cross))

        return self.U_0 - self.U_act - self.U_ohm - self.U_conc

    def Membranleitwert(self, rH):
        if self.Lambda_calc(rH) < 1.253:
            Leitwert = 1e-13
        else:
            Leitwert = (0.005738 * self.Lambda_calc(rH) - 0.007192) * 100

            Leitwert = max(Leitwert, 0.2)  # for numeric stability

        return Leitwert  # S/m oder 1/Ohm*m

    def Lambda_calc(self, rH):
        # Kulikovsky, Andrei A.: Chapter 1 - Fuel cell basics. In: Kulikovsky, Andrei A. (Herausgeber): Analytical Modelling of Fuel Cells, Seiten 1 – 38. Elsevier, Amsterdam, 2010

        Lambda = 0.3 + 6 * rH * (1 - math.tanh(rH - 0.5)) + 3.9 * rH ** 0.5 * (1 + math.tanh((rH - 0.89) / 0.23))
        return Lambda

    def Thermo(self):
        if self.pTm_COi['m_in'] != 0:
            self.U_loss_stack = self.U_0 * self.P['NumCell'] - self.U_stack
            self.Qdot_stack = self.U_loss_stack * self.I_stack
            cp_water=4200 #J/KgK
            self.Qdot_cool = (self.P['Active_area'] * ((self.P['NumCell'] - 1) * 2 + 2)) * self.P['Stack']['Thermo']['h'] * (self.pTm_COi['T_in'] - self.T_stack)
            # self.Qdot_cool = (self.pTm_COi['m_in'] * cp_water) * (self.T_stack - self.pTm_COi['T_in'])
            self.T_stack += (self.Qdot_stack + self.Qdot_cool)/(self.P['Stack']['Thermo']['M'] * self.P['Stack']['Thermo']['cp']) * self.dt
            self.pTm_COo['T_out'] = self.pTm_COi['T_in'] - self.Qdot_cool/(self.pTm_COi['m_in'] * cp_water)
            # self.pTm_COo['T_out'] = self.Qdot_cool/((self.P['Active_area'] * ((self.P['NumCell'] - 1) * 2 + 2)) * self.P['Stack']['Thermo']['h']) + self.T_stack

            self.pTm_COo['m_out'] = self.pTm_COi['m_in']
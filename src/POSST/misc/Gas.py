import numpy as np
import math


class Massflow():
    def __init__(self, NumCell):
        self.NumCell=NumCell
        self.Faraday=96485.3399
        self.Avogadro=6.02214E+23
        self.el_H2=2
        self.M_H2=1.00784*2 #g/mol
        self.el_O2=4
        self.M_O2=15.999*2 #g/mol
        self.M_H2O=18.01528
        self.M_N2=28.0134
        self.roh_H2=0.0899
        self.roh_O2=1.293
        self.nl=22.414  # nl/mol
        self.p_st_mbar=1013.25
        self.T_st=273.15

    def set_m_H2(self, I=0, Stoic=2):
        SCs=I*self.NumCell
        H2Cs=SCs*Stoic
        H2mols=H2Cs/self.Faraday/self.el_H2
        mH2=H2mols*self.M_H2
        return mH2
    
    def set_m_O2(self, I=0, Stoic=2):
        SCs=I*self.NumCell
        O2Cs=SCs*Stoic
        O2mols=O2Cs/self.Faraday/self.el_O2
        mO2=O2mols*self.M_O2
        return mO2
    
    def I_H2(self, mH2=0, Stoic=2):
        H2mols=mH2/self.M_H2
        H2Cs=H2mols*self.Faraday*self.el_H2
        SCs=H2Cs/Stoic
        I=SCs/self.NumCell
        return I
    
    def I_O2(self, mO2=0, Stoic=2):
        O2mols=mO2/self.M_O2
        O2Cs=O2mols*self.Faraday*self.el_O2
        SCs=O2Cs/Stoic
        I=SCs/self.NumCell
        return I
    
    def St_H2(self, mH2=0, I=0):
        Stoic=0
        if I>0:
            H2mols=mH2/self.M_H2
            H2Cs=H2mols*self.Faraday*self.el_H2
            SCs=I*self.NumCell
            Stoic=H2Cs/SCs
        return Stoic

    def St_O2(self, mO2=0, I=0):
        Stoic=0
        if I>0:
            O2mols=mO2/self.M_O2
            O2Cs=O2mols*self.Faraday*self.el_O2
            SCs=I*self.NumCell
            Stoic=O2Cs/SCs
        return Stoic
    
    def rH(self,p=1 ,T=273.15, m=0.1, H2=0, H2O=0, O2=0, N2=1):
        pmb=p*1000
        nO2 = m * O2 / self.M_O2
        nN2 = m * N2 / self.M_N2
        nH2 = m * H2 / self.M_H2
        nH2O = m * H2O / self.M_H2O

        nges= nO2+nN2+nH2+nH2O
        try:
            H2Omol=nH2O/nges
        except:
            H2Omol=0
        partH2O=H2Omol*pmb
        if T>0:
            # satH2O=self.p_sat(T)/100
            satH2O=self.p_sat(T)
        else:
            satH2O=.1
        
        RH=partH2O/satH2O*100
        return RH
    
    def m_H2O(self, T=65+273.15, p=1.5, rH=60, m=100, O2=0.231461, N2=0.768539, H2=0, H2O=0):
        pmb=p*1000
        mdry=m-m*H2O
        nO2=m*O2/self.M_O2
        nN2=m*N2/self.M_N2
        nH2=m*H2/self.M_H2
        # nH2O=m*H2O/self.M_H2O
        try:
            xO2=nO2/(nN2+nO2+nH2)
        except:
            xO2=0.1
        try:
            xN2=nN2/(nN2+nO2+nH2)
        except:
            xN2=0.1
        try:
            xH2=nH2/(nN2+nO2+nH2)
        except:
            xH2=0.1
        M_dry=xO2*self.M_O2 + xN2*self.M_N2 + xH2*self.M_H2
        ndry = mdry / M_dry
        TP=self.rH_to_TP(rH,T)
        m_H2O_C_g_s = ndry * self.M_H2O * self.p_sat(TP) / (pmb - self.p_sat(TP))


        if math.isnan(m_H2O_C_g_s):
            m_H2O_C_g_s=0
        return max(m_H2O_C_g_s,0)
    
    def p_sat(self, T):
        #in Pa
        p_sat_mbar=np.exp(-5800.2206/T+1.3914993+(-0.048640239)*T+0.000041764768*pow(T,2)+(-0.000000014452093)*pow(T,3)+6.545973*np.log(T))/100
        return p_sat_mbar
    
    def rH_to_TP(self,rH,T):
        T=T-273.15
        a=7.5
        b=237.3
        SDD=6.1078*pow(10,((a*T)/(b+T)))
        DD=rH/100*SDD
        v=np.log10(DD/6.1078)
        TP=b*v/(a-v)
        TP=273.15+TP
        return TP
    
    def mol_anteil(n1=1.,n2=0.,n3=0.,n4=0.):
        return (n1)/(n1 + n2 + n3 + n4)

class path:
    def __init__(self):
        self.CNat=CNat()
    def m_H2O_prod(self, I, Num_Cell):
        return (I * Num_Cell) / (self.CNat.Faradey * self.CNat.Ze_O2) * self.CNat.M_O2 + (I * Num_Cell) / (self.CNat.Faradey * self.CNat.Ze_H2) * self.CNat.M_H2
    def mdot_reac_out(self, m_in, x, I, M, Ze, Num_Cell):
        return ((m_in * x / M) - (I * Num_Cell) / (self.CNat.Faradey * Ze)) * M
    def dp(self, rho, v, l_ch, D_H, Re, zeta):
        return 0.5 * rho * v ** 2 * l_ch / D_H * 64 / Re + 0.5 * rho * v ** 2 *zeta

    def Re(self, rho, v, D_H, nu):
        return (D_H * rho * v) / nu

    def D_H(self, w, h):
        return 4 * (w * h) / (2 * w + 2 * h)

    def nu(self, nu_n, T_in, T_n):
        return nu_n * (T_in / T_n) ** 0.5

    def v(self, m_in, rho, w, h):
        return (m_in / rho) / (w * h)

    def v_2(self, p_in, p_out, rho, l, D_H, Re, zetha):
        dp=p_in-p_out
        return abs((dp)/(0.5*rho*zetha + 0.5*rho*l/D_H*64/Re))**0.5*(dp)/abs(dp)

    def rho(self, rho_n, p, p_n, T_in, T_n):
        return rho_n * p / p_n * T_n / T_in

    def rho_mix(self, rho_i, xn_i):
        sum = 0
        for i in range(len(rho_i)):
            sum += xn_i[i]/rho_i[i]
        return 1 / sum


class gas_state:
    def __init__(self):
        self.dp=None
        self.Re=None
        self.D_H=None
        self.nu=None
        self.v=None
        self.rho=None
class CNat:
    def __init__(self):
        self.Faradey = 96485.3399  # J/(V*mol)
        self.Avogadro = 6.02214E+23  # n/mol
        self.Ze_H2 = 2
        self.Ze_O2 = 4
        self.R = 8.314  # J/(mol*K)
        self.M_O2 = 15.999 * 2 /1000  # Kg/mol
        self.M_N2 = 28.0134 /1000  # Kg/mol
        self.M_H2 = 1.00784 * 2 /1000  # Kg/mol
        self.M_H2O = 18.01528 /1000  # Kg/mol
        self.rho_air = 1.293  # dry Kg/m3
        self.rho_H2 = 0.0899  # Kg/m3
        self.rho_H2O = 998 # Kg/m3
        self.rho_O2 = 1.429 # Kg/m3
        self.rho_N2 = 1.165 # Kg/m3
        self.T_0K = 273.15  # K
        self.T_20K = 293.15 # K
        self.p_atm = 1.01325e5  # Pa
        self.p_1bar = 1e5
        self.nu_n_air = 13.3 * 1e-6
        self.nu_n_H2 = 8.3969 * 1e-6  # Pa*s

if __name__ == '__main__':
    MF=Massflow(NumCell=560)
    print(MF.set_m_O2(I=357, Stoic=2)/0.23)



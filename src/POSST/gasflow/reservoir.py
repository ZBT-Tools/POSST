
class model_1():
    r"""
    ## Description  
    Model to define Volume of Container and Calculation of states.  
    
    ### Calculation
    
    1. Change of pressure in Reservoir  
    
        $$\frac{\delta p}{\delta t} = \frac{\kappa \cdot R}{V} \cdot  \left(\dot{m}_{\text{in}} \cdot T_{\text{in}} - \dot{m}_{\text{out}} \cdot T_{\text{out}}\right)$$
    
    2. Change of temperature in Reservoir  
    
        $$\frac{\delta T}{\delta t} = \frac{T \cdot R}{p \cdot V \cdot c_v} \cdot  \left(c_p \cdot \dot{m}_{\text{in}} \cdot T_{\text{in}} - 
                                            c_p \cdot \dot{m}_{\text{out}} \cdot T_{\text{out}} - c_v \cdot \left(\dot{m}_{\text{in}}-\dot{m}_{\text{out}}\right)\right)$$
    """  
    def __init__(self, name = "Pipe", V=0.001, T_0 = 293.15, p_0 = 1.0135e5, cp = 1000 ,R =237, dt=1e-3): # 8.314):
        r"""
        ## Arguments
            V: volume of container element in [$m^{3]$]
            T_0: initial temperature in [$K$]
            p_0: initial pressure in [$Pa$]
            cp: specific heat capacity in [$J \cdot Kg^{-1} \cdot K^{-1}}$]
            R: gas constant in [$J \cdot Kg^{-1} \cdot K^{-1}}$]
            dt: time step in [$s$]

        ## Returns
            Instance of Reservoir

        """
        self.name = name
        self.V = V
        self.cp= cp
        self. R= R
        self.cv =cp-R
        self.kappa = self.cp/self.cv
        self.p = p_0
        self.T = T_0
        self.T_in = self.T
        self.m_in = self.m_out = 0
        self.dt =dt

    def __call__(self ,m_in, m_out, T_in):
        """
        ## Arguments
            m_in: inlet massflow in [$Kg \cdot s^{-1}$]
            m_out: outlet massflow in [$Kg \cdot s^{-1}$]

        ## Calculated
            p: current pressure in [$Pa$]
            T: current Temperature in [$K$]
        """
        # T_in+=273.15
        self.T_in_old = self.T_in
        self.T_in = T_in
        self.m_in_old = self.m_in
        self.m_in = m_in
        self.m_out_old = self.m_out
        self.m_out = m_out
        self.p_old = self.p
        self.T_old = self.T
        self.p += ((self.R * self.kappa) / self.V) * (
            (self.m_in * self.T_in + self.m_in_old *self.T_in_old ) /2 -(self.m_out * self.T_old)) * self.dt

        self.T += ((self.T * self.R) / (self.p_old * self.V * self.cv))\
            * (self.cp * (self.T_in * self.m_in + self.m_in_old * self.T_in_old)/2
            - self.cp * self.T * m_out
            - self.cv*((self.m_in + self.m_in_old)/2
            - self.m_out) * self.T) * self.dt

        self.T = max(self.T, 1)
        self.p = max(self.p, 1e-9)
        None
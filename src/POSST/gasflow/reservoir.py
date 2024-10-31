
class model_1():
    r"""
    ## Arguments
    
    
        V = volumne of container element  
        T_0 = start temperature   
        p_0 = start pressure  
        cp = specific heat capacity  
        R = gas constant  
        dt = time step  
    ## Returns
        
    
        p = current pressure
        T = current Temperature
    ## Description  
    Model to define Volumne of Container and Calculation of states.  
    ### Calculation
    1. Change of pressure in Reservoir  
    
    $$\frac{\delta p}{\delta t} = \frac{\kappa \cdot R}{V} \cdot  \left(\dot{m}_{in} \cdot T_{in} - \dot{m}_{out} \cdot T_{out}\right)$$
    
    2. Change of temperature in Reservoir  
    
    $$\frac{\delta T}{\delta t} = \frac{T \cdot R}{p \cdot V \cdot c_v} \cdot  \left(c_p \cdot \dot{m}_{in} \cdot T_{in} - 
                                        c_p \cdot \dot{m}_{out} \cdot T_{out} - c_v \cdot \left(\dot{m}_{in}-\dot{m}_{out}\right)\right)$$
    """
    
    def __init__(self, name = "Pipe", V=0.001, T_0 = 293.15, p_0 = 1.0135e5, cp = 1000 ,R =237, dt=1e-3): # 8.314):
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
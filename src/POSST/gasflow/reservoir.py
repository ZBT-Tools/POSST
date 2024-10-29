
class model_1():
    """Volume of container/pipe etc."""
    def __init__(self, name = "Pipe", V=0.001, T_0 = 293.15, p_0 = 1.0135e5, cp = 1000 ,R =237, dt=1e-3): # 8.314):
        """Initialisierung der Leitungen"""
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
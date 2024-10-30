import numpy as np
eps=1e-13
class model_1():
    """
    Args:
        
        dt=time step  
        dn=??   
        cp=specific heat capacity  
        R= gas constant  
        repo="./POSST"  
    Returns:
        
        m=massflow
    Description:
        
        **Model for fitting curve 1**  
        $ m = y\cdotx$
    """
    def __init__(self, dt, dn, cp=1000, R=237, repo_dir='./POSST'):
        self.dt = dt
        self.cp = cp
        self.R = R
        self.cv = cp - R
        self.kappa = self.cp / self.cv

        self.m_fit = interp_Comp_m(ddir=repo_dir + '/BZ_Sim/Comp_interp/data')
        self.m = 0
        self.n_t = 0
        self.dn = dn * self.dt
        self.iteration = 0
        self.P_ver = 0
        self.etha = 0.7
        self.T_out = 293.15

    def __call__(self, n, p_in, p_out, T_in=293.15):
        self.iteration += 1
        p_ratio = p_out / p_in
        self.n_t = np.clip(n, max(self.m_fit.min_n, self.n_t - self.dn), min(self.m_fit.max_n, self.n_t + self.dn))
        if n > 1e-7:
            try:
                self.m = self.m_fit(self.n_t, p_ratio)
            except:
                self.m = max(self.m, 1e-7)

            self.P_ver = self.m * self.cp * T_in * (p_ratio ** ((self.kappa - 1) / self.kappa) - 1) / self.etha
            self.T_out = (T_in * p_ratio ** ((self.kappa - 1) / self.kappa) - T_in) / self.etha + T_in
        else:
            self.m = 0

        return self.m


class model_2():
    """
    Args:
        
        dt=time step  
        dn=??   
        cp=specific heat capacity  
        R= gas constant  
        repo="./POSST"  
    Returns:
        
        m=massflow
    Description:
        
        **Model for fitting curve 2**  
        $ m = y\cdotx$
    """
    def __init__(self, dt, dn, cp=1000, R=237, repo_dir='./POSST'):
        self.dt = dt
        self.cp = cp
        self.R = R
        self.cv = cp - R
        self.kappa = self.cp / self.cv
        self.dn = dn * self.dt
        self.m = 0
        self.n_t = 0

        import restriction as res
        self.total_p = interp_Comp_p(ddir=repo_dir + '/BZ_Sim/Comp_interp/data')
        # self.Res_Co = Restriction(A=np.pi*(2e-3)**2, dt=1e-3, R=287)
        self.Res_Co = res.model_1(l=1.000e-03, b=1.633e+00, h=1.633e+00, R=287)

        self.P_ver = 0
        self.etha = 0.7
        self.T_out = 293.15

    def __call__(self, n, p_in, p_out, T_in=293.15):
        self.n_t = np.clip(n, max(self.total_p.min_n, self.n_t - self.dn), min(self.total_p.max_n, self.n_t + self.dn))
        if self.n_t >= 1e-7:
            self.m = self.Res_Co.m
            self.m = abs(self.m)
            self.p_comp = max(self.total_p(n=n, m=self.Res_Co.m) * p_in, p_in + 1.01)

            p_ratio = self.p_comp / p_in
            self.P_ver = self.m * self.cp * T_in * (p_ratio ** ((self.kappa - 1) / self.kappa) - 1) / self.etha
            self.T_out = (T_in * p_ratio ** ((self.kappa - 1) / self.kappa) - T_in) / self.etha + T_in
        else:
            self.p_comp = p_in
        self.Res_Co(p_in=p_in, p_out=p_out, T_in=T_in)

        return self.m

class interp_Comp_p():
    def __init__(self, ddir='./data'):
        import scipy
        self.data_raw = load_comp(data_dir=ddir)
        self.data_raw |= {'0': {'etha': [1], 'm':[eps], 'p':[1.], 'n':[0]}}
        self.m_grid, self.n_grid, self.p_grid = self.grid_calc()
        self.interp = scipy.interpolate.RegularGridInterpolator((self.m_grid[:,0], self.n_grid[0,:]), self.p_grid)

    def __call__(self, n, m):
        self.m_clip = np.clip(m, self.min_m, self.max_m)
        self.n_clip = np.clip(n, self.min_n, self.max_n)
        self.p = self.interp((self.m_clip, self.n_clip))
        return max(self.p, 0)

    def grid_calc(self):
        import numpy as np
        import scipy
        import copy
        self.m_ori= [0.00000000]
        self.p_ori= [1]
        self.n_ori= [0]

        for i in self.data_raw.keys():
            if i != 'Surge':
                self.m_ori += copy.deepcopy(list(self.data_raw[i]['m']))
                self.p_ori += copy.deepcopy(list(self.data_raw[i]['p']))
                self.n_ori += list(np.full(len(self.data_raw[i]['n']), int(i)))

        self.p_ori = np.array(self.p_ori)
        self.m_ori = np.array(self.m_ori)
        self.n_ori = np.array(self.n_ori)

        self.max_m = self.m_ori.max()
        self.min_m = self.m_ori.min()
        self.max_n = self.n_ori.max()
        self.min_n = self.n_ori.min()
        self.max_p = self.p_ori.max()
        self.min_p = self.p_ori.min()

        self.p_fit_arr = [1.]
        self.m_fit_arr = [eps]
        self.n_fit_arr = [0]
        for i in self.data_raw.keys():
            if i != '0':
                self.m_fit = np.linspace(self.min_m,self.max_m,100)
                para=self.p_ex(m=copy.deepcopy(list(self.data_raw[i]['m'])),
                               p=copy.deepcopy(list(self.data_raw[i]['p'])))
                self.p_fit_arr += list(p_fit(self.m_fit[np.where(self.m_fit<max(self.data_raw[i]['m']))], *para))
                self.m_fit_arr += list(self.m_fit[np.where(self.m_fit<max(self.data_raw[i]['m']))])
                self.n_fit_arr += list(np.full(len(self.m_fit[np.where(self.m_fit<max(self.data_raw[i]['m']))]), int(i)))
            else:
                self.p_fit_arr += copy.deepcopy(list(self.data_raw[i]['p']))
                self.m_fit_arr += copy.deepcopy(list(self.data_raw[i]['m']))
                self.n_fit_arr += copy.deepcopy(list(self.data_raw[i]['n']))
            None


        grid_x, grid_y = np.mgrid[self.min_m:self.max_m:100j, self.min_n:self.max_n:100j]
        grid_z = scipy.interpolate.griddata((self.m_fit_arr, self.n_fit_arr), self.p_fit_arr, (grid_x, grid_y), method='linear', fill_value=1)

        return grid_x, grid_y, grid_z

    def p_ex(self, m, p):
        import scipy
        popt, pcov = scipy.optimize.curve_fit(p_fit, m, p)
        return popt

class interp_Comp_m():
    def __init__(self, ddir='./data'):
        import scipy
        self.data_raw = load_comp(data_dir=ddir)
        self.data_raw |= {'0': {'etha': [1], 'm':[eps], 'p':[1.], 'n':[0]}}
        self.p_grid, self.n_grid, self.m_grid = self.grid_calc()
        self.interp = scipy.interpolate.RegularGridInterpolator((self.p_grid[:,0], self.n_grid[0,:]), self.m_grid)

    def __call__(self, n, p_ratio):
        self.p_clip = np.clip(p_ratio, self.min_p, self.max_p)
        self.n_clip = np.clip(n, self.min_n, self.max_n)
        self.m = self.interp((self.p_clip, self.n_clip))
        return self.m

    def grid_calc(self):
        import numpy as np
        import scipy
        import copy
        self.m_ori= [0.]
        self.p_ori= [1]
        self.n_ori= [0.]

        for i in self.data_raw.keys():
            if i != 'Surge':
                self.m_ori += copy.deepcopy(list(self.data_raw[i]['m']))
                self.p_ori += copy.deepcopy(list(self.data_raw[i]['p']))
                self.n_ori += list(np.full(len(self.data_raw[i]['n']), int(i)))

        self.p_ori = np.array(self.p_ori)
        self.m_ori = np.array(self.m_ori)
        self.n_ori = np.array(self.n_ori)

        self.max_m = self.m_ori.max()
        self.min_m = self.m_ori.min()
        self.max_n = self.n_ori.max()
        self.min_n = self.n_ori.min()
        self.max_p = self.p_ori.max()
        self.min_p = self.p_ori.min()

        self.p_fit_arr = [1.]
        self.m_fit_arr = [eps]
        self.n_fit_arr = [0.]
        for i in self.data_raw.keys():
            if i != '0':
                self.m_fit = np.linspace(self.min_m,self.max_m,100)
                para=self.p_ex(m=copy.deepcopy(list(self.data_raw[i]['m'])),
                               p=copy.deepcopy(list(self.data_raw[i]['p'])))
                self.p_fit_arr += list(p_fit(self.m_fit[np.where(self.m_fit < max(self.data_raw[i]['m']))], *para))
                self.m_fit_arr += list(self.m_fit[np.where(self.m_fit < max(self.data_raw[i]['m']))])
                self.n_fit_arr += list(np.full(len(self.m_fit[np.where(self.m_fit < max(self.data_raw[i]['m']))]), int(i)))
            else:
                self.p_fit_arr += copy.deepcopy(list(self.data_raw[i]['p']))
                self.m_fit_arr += copy.deepcopy(list(self.data_raw[i]['m']))
                self.n_fit_arr += copy.deepcopy(list(self.data_raw[i]['n']))
            None


        grid_x, grid_y = np.mgrid[self.min_p:self.max_p:100j, self.min_n:self.max_n:100j]
        grid_z = scipy.interpolate.griddata((self.p_fit_arr, self.n_fit_arr), self.m_fit_arr, (grid_x, grid_y), method='linear', fill_value=eps)

        return grid_x, grid_y, grid_z

    def p_ex(self, m, p):
        import scipy
        popt, pcov = scipy.optimize.curve_fit(p_fit, m, p)
        return popt

def load_comp(data_dir = '/data'):
    '''
    Args:
        
        data_dir:  
    Returns:
        
        dict:  
            
    ./data_dir
    
    
        |-rpm0.h5  
        |-rpm1.h5  
        :  
        |-rpmN.h5  


    rpm0.h5  
    
        |-m: [...]  
        |-p: [...]  
        |-n: [...]  
        |-etha: [...]  
    '''

    import h5py
    import os

    #todo: fit p/m per rpm for extrapol
    files = os.listdir(data_dir)
    DataL = {}
    for i in files:
        try:
            DataL |= {i[0:-3]: h5py.File(data_dir + '/' + i, 'r')}
        except:
            print(i + ' not loaded')
    Data={}
    for i in DataL.keys():
        if i != 'Surge':
            Data |= {i:{}}
            Data[i] |= {'m': []}
            Data[i] |= {'p': []}
            Data[i] |= {'n': []}
            Data[i] |= {'etha': []}
            for j in range(len(DataL[i]['data']['submap_0']['domain_data'])):
                Data[i]['m'] += DataL[i]['data']['submap_0']['domain_data'][j][0],
                Data[i]['p'] += DataL[i]['data']['submap_0']['codomain_data'][j][0],
                Data[i]['etha'] += DataL[i]['data']['submap_0']['codomain_data'][j][1],
                Data[i]['n'] += int(i),
    return Data

def p_fit(m, p0, p1, p2):
    # p=p0-p1*np.exp(m+p2)
    p=p0*m**(p1)+p2
    return p
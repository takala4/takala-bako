import numpy as np
import pandas as pd
import gurobipy as gp



class model:
    def __init__(self, M, b, dir_path=False):
        
        # input file date
        if dir_path:
            self.input_data(dir_path)
            
        # input matrix and vector data
        else:
            self.M = np.array(M, dtype=np.float64)
            self.M_ = self.M + self.M.T
            self.b = np.array(b, dtype=np.float64)
            self.N = self.M.shape[0]


    def input_data(self, dir_path):
        self.M = pd.read_csv(dir_path+'\\M.csv',
                             encoding='utf-8-sig',
                             index_col=None, dtype=float).values # \in R^{N*N}
        
        self.b = pd.read_csv(dir_path+'\\b.csv',
                             encoding='utf-8-sig',
                             index_col=None, dtype=float).values.flatten() # \in R^{N}
        
        # Rank of Matrix M
        self.N = self.M.shape[0] # \in R     
        
        return True

    def solve(self, x0=None, err=10**(-6), max_itr=100):

        if type(x0) == type(None):
            x_pre = self.calc_x_init()
        else:
            x_pre = x0
            


        x_list = [x_pre]
        def obj_value(x): return np.ones(self.N)@((self.M + np.eye(self.N))@x + self.b - abs((self.M - np.eye(self.N))@x + self.b))
    
    
        for i in range(max_itr):
                
            if abs(obj_value(x_pre)) < err:
                break
        
            y = self.SubProb(x_pre)

            alpha = 1 #self.Calc_StepSize(x_pre, y)
        
            if abs(obj_value((1-alpha)*x_pre + alpha * y)-obj_value(x_pre)) < 10**(-5):
                break
                
            x_pre = (1-alpha)*x_pre + alpha * y
            x_list = x_list + [x_pre]
                
        return x_pre, x_list

    def calc_x_init(self):
        diff_Matrix = np.eye(self.N) - np.roll(np.eye(self.N), 1, axis=1)
        
        model = gp.Model()
        
        x_ini = model.addMVar(self.N, lb=0.0)
        x_a = model.addMVar(self.N, lb=0.0)
            
        model.setObjective(np.ones(self.N)@x_a)
        model.addConstr(self.M@x_ini >= - self.b - np.eye(self.N)@x_a)
        model.addConstr(diff_Matrix@x_a == np.zeros(self.N))
        
        model.setParam('OutputFlag', 0) # Mute output
    
        model.optimize()
            
        return x_ini.X


    def SubProb(self, x_pre):
    
        model = gp.Model()

        x = model.addMVar(self.N, lb=0.0)
        I = np.eye(self.N)
        e = np.ones(self.N)

        C = e@((self.M + I) - np.diag(np.sign((self.M - I)@x_pre + self.b)) @ (self.M - I))

        model.setObjective(C@x)  
        
        
        model.addConstr(self.M@x >= - self.b)
        
        model.setParam('OutputFlag', 0) # Mute putput
        model.optimize()
    
        return x.X



    def Calc_StepSize(self, x_pre, y):
        
        
        def obj_alpha(alpha, x_pre, y):
            x_new = (1-alpha)*x_pre + alpha * y
            return x_new@self.M_@x_new + self.b@x_new
          
    
        d = y - x_pre
        g__ = np.dot(d, np.dot(self.M_ , d))
        def g_(alpha): return np.dot(d, np.dot(self.M_, d)) * \
            alpha + np.dot(d, (np.dot(self.M_, x_pre)+self.b))
    
        if g__ > 0:
            if g_(0) > 0:
                alpha = 0
            elif g_(1) < 0:
                alpha = 1
            else:
                alpha = - g_(0)/g__
        else:
            if obj_alpha(0, x_pre, y) < obj_alpha(1, x_pre, y):
                alpha = 0
            else:
                alpha = 1
        
        
        return alpha

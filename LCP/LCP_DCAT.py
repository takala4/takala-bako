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

            self.ML = np.tril(self.M, -1)
            self.MU = np.triu(self.M, +1)
            self.MD = np.diag(np.diag(self.M))
            self.MDp = np.where(self.MD>0.0,  self.MD,  0.0)
            self.MDm = np.where(self.MD<0.0,  self.MD,  0.0)

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
        def obj_value(x): return (x@self.M@x) + self.b@x

    

        for i in range(max_itr):
                
            if abs(obj_value(x_pre)) < err:
                break

            print(obj_value(x_pre))


            y = self.SubProb(x_pre, 0)

#             alpha = self.Calc_StepSize(x_pre, y)
            alpha = 1
        
            if abs(obj_value((1-alpha)*x_pre + alpha * y)-obj_value(x_pre)) < err:
                print('Break-code-1')
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
        
        model.setParam('OutputFlag', 0) # Mute putput
    
        model.optimize()
            
        return x_ini.X


    def SubProb(self, x_pre, epsilon=1):
    
        model = gp.Model()

        x = model.addMVar(self.N, lb=0.0)
        H = self.MU + self.MU.T + 2*self.MDm - 2*epsilon*np.eye(self.N)
        P = self.ML + self.MDp  + epsilon*np.eye(self.N)
        model.setObjective((x_pre@H+self.b)@x + x@(0.5*(P + P.T))@x)
        model.addConstr(self.M@x >= - self.b)
        model.setParam('OutputFlag', 0)  # Mute putput
        model.params.NonConvex=2
        model.optimize()
    
        return x.X



    def Calc_StepSize(self, x_pre, y):
               

        M_over = self.M + self.M.T
        
        def obj_alpha(alpha, x_pre, y):
            x_new = (1-alpha)*x_pre + alpha * y
            return (np.dot(x_new.T, np.dot(self.M, x_new))) +self.b@x_new
          
    
        d = y - x_pre
        g__ = np.dot(d, np.dot(M_over , d))
        g_ = lambda alpha : np.dot(d, np.dot(M_over , d))*alpha + np.dot(d, (np.dot(M_over , x_pre)+self.b))
    
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

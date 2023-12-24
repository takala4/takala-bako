# This algorithm based on Facchinei and Soares 1995,1997.
# 1995:"TESTING A NEW CLASS OF ALGORITHMS FOR NONLINEAR COMPLEMENTARITY PROBLEMS"
# 1997:"A New Merit Function for Nonlinear Complementarity Problems and a Related Algorithm"


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
            self.M = np.array(M, dtype=float)
            self.b = np.array(b, dtype=float)
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
    
    # Calculate F(x) = Mx + b
    def F(self, x):
        return np.dot(self.M, x) + self.b # \in R^{N}
    
    # φ 
    def phi(self, x):
        return np.sqrt(x**2+self.F(x)**2) - (x + self.F(x)) # \in R^{N}
    
    # Ψ (merit function value)
    def Psi(self, x):
        return 1/2*np.sum(self.phi(x)**2) # \in R^{N}
    
    # Jacobian of phi
    def Dphi(self, x):
        # 分母
        denominator = np.sqrt(x**2 + self.F(x)**2) # \in R^N
        
        A_vec = np.ones((self.N,), dtype=float)
        B_vec = np.ones((self.N,), dtype=float)
        
        # 分母が0 ならξ, ρをA[i,i], B[i,i]に採用
        
        # Search 0 in the denominator
        ind_nz = np.where(denominator >= self.TOL_zero_div)[0] # \in R^N
        ind_z  = np.where(denominator < self.TOL_zero_div)[0] # \in R^N

        A_vec[ind_nz] = x[ind_nz] / denominator[ind_nz]
        B_vec[ind_nz] = self.F(x)[ind_nz] / denominator[ind_nz]
        
        # If the denominator is 0, 
        # we use A[i,i] and B[i,i] as epsilon and rho respectively
        A_vec[ind_z] = 1/np.sqrt(2) - 1 #0.5 # ? 
        B_vec[ind_z] = 1/np.sqrt(2) - 1 #0.5 # ?

        
        A = np.diag(A_vec) # \in R^{N*N}
        B = np.diag(B_vec) # \in R^{N*N}
        
        I = np.eye(self.N, dtype=float)
        
        return A - I + np.dot((B-I), self.M) # \in R^{N*N}
    
    
    # Gradient of Psi
    def DPsi(self, x):
        return np.dot(self.Dphi(x).T, self.phi(x))
    
    
    # ---
    # Algorithm
    # ---
    # 降下方向の計算 (Step 2)
    def discent_direction(self, x):
        # まずは式(18), (19)に従って local direction を求める
        # 未知変数を A と N とに分割
        F = self.F(x)
        A_ind = np.where(x<=self.eps*F)[0]
        N_ind = np.where(x>self.eps*F)[0]
        
        x_A = x[A_ind]
        F_N = F[N_ind]
        DF_AN = self.M[N_ind,:][:,A_ind]
        DF_NN = self.M[N_ind,:][:,N_ind]
        
        # 降下方向
        d_A = -x_A
        # 式(19)の連立方程式が解けるなら (18)-(19)に従って dを求める
        if len(N_ind) > 0 and np.linalg.matrix_rank(DF_NN) == len(N_ind):
            if len(A_ind) > 0:
                tmp_term = np.dot(DF_AN, x_A)
            else:
                tmp_term = np.zeros((len(N_ind),))
            d_N = np.linalg.solve(DF_NN, -F_N + tmp_term)
            d = np.zeros((self.N,),dtype=float)
            if len(A_ind) > 0:
                d[A_ind] = d_A
            d[N_ind] = d_N
            #print(d_A, d_N, d)
        # 式(19)の連立方程式が解を持たないなら-∇Ψを降下方向とする
        else:
#             print("d = -∇Ψ Step 2")
            d = - self.DPsi(x)
        return d

    # 線形探索(Step 3, Step 4)
    def line_search(self, x, d):
        # 式(24)を満足するなら線形探索をスキップ
        if (self.Psi(x+d) <= self.sigma*self.Psi(x)):
            #print("skip linear search")
            new_x = x + d
        else:
            # 式(25)を「満足しない」なら d を∇Ψ に設定
            if ( np.dot(self.DPsi(x).T, d) > -self.rho*np.linalg.norm(d,2) ** self.p ):
#                 print("d = -∇Ψ in Step 4")
                d = - self.DPsi(x)
            # Armijo の規則的な方法でステップサイズを決定
            i = 0
            while True:
                if (  self.Psi(x + 2**(-i)*d) <=
                      self.Psi(x) + self.beta*2**(-i)* np.dot(self.DPsi(x).T,d)):
                    break
                i = i + 1
            #print("i=%d" % i)
            # x を改訂
            new_x = x + 2**(-i)*d
        return new_x
    
    # 繰返し計算
    def solve(self, 
              x0=None, 
              err = 1e-8,
              max_itr = 1000,
              eps=1e-8, # > 0
              TOL_zero_div=1e-8, # > 0
              rho=1e-4, # > 0
              p=1.1, # > 1
              beta = 0.4, # (0,0.5)
              sigma = 0.95, # (0, 1) 
              ):
        self.eps = eps
        self.TOL_zero_div = TOL_zero_div
        self.rho = rho
        self.p = p
        self.beta = beta
        self.sigma = sigma
        self.max_itr = max_itr
        self.err = err


        if type(x0) == type(None):
            x = self.calc_x_init()
        else:
            x = x0
            
        
        x_list = [x]
        
        for itr in range(self.max_itr):
            if (self.Psi(x) < self.err): 
                break
            d = self.discent_direction(x)
            new_x = self.line_search(x, d)
            np.set_printoptions(suppress=False, precision=4)
            
            x_list = x_list + [new_x]
            x = new_x

        return x, x_list

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
import numpy as np
import math

from numpy.linalg.linalg import norm
import gurobipy as gp

# https://github.com/GRYE/Nesterov-accelerated-gradient-descent/blob/master/nesterov_method.py

class model:
    def __init__(self, M, b, dir_path=False):
        
        # input file date
        if dir_path:
            self.input_data(dir_path)
            
        # input matrix and vector data
        else:
            self.M = np.array(M, dtype=np.float64)
            self.b = np.array(b, dtype=np.float64)
            self.N = self.M.shape[0]


    def input_data(self, dir_path):
        self.M = pd.read_csv(dir_path+'\\M.csv',
                             encoding='utf-8-sig',
                             index_col=None, dtype=np.float64).values # \in R^{N*N}
        
        self.b = pd.read_csv(dir_path+'\\b.csv',
                             encoding='utf-8-sig',
                             index_col=None, dtype=np.float64).values.flatten() # \in R^{N}
        
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
        TOL_zero_div = 10**(-10)
        denominator = np.sqrt(x**2 + self.F(x)**2) # \in R^N
        
        A_vec = np.ones((self.N,), dtype=np.float64)
        B_vec = np.ones((self.N,), dtype=np.float64)
        
        # 分母が0 ならξ, ρをA[i,i], B[i,i]に採用
        
        # Search 0 in the denominator
        ind_nz = np.where(denominator >= TOL_zero_div)[0] # \in R^N
        ind_z  = np.where(denominator <  TOL_zero_div)[0] # \in R^N

        A_vec[ind_nz] = x[ind_nz] / denominator[ind_nz]
        B_vec[ind_nz] = self.F(x)[ind_nz] / denominator[ind_nz]
        
        # If the denominator is 0, 
        # we use A[i,i] and B[i,i] as epsilon and rho respectively
        A_vec[ind_z] = 1/np.sqrt(2) - 1 #0.5 # ? 
        B_vec[ind_z] = 1/np.sqrt(2) - 1 #0.5 # ?

        
        A = np.diag(A_vec) # \in R^{N*N}
        B = np.diag(B_vec) # \in R^{N*N}
        
        I = np.eye(self.N, dtype=np.float64)
        
        return A - I + np.dot((B-I), self.M) # \in R^{N*N}
    
    
    # Gradient of Psi
    def DPsi(self, x):
        return np.dot(self.Dphi(x).T, self.phi(x))



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


     # Ψ (merit function value) #https://www.kurims.kyoto-u.ac.jp/~kyodo/kokyuroku/contents/pdf/0987-15.pdf
    def MngMerit(self, x):
        a = 1.5
        Fx = self.F(x)


        prox_p = np.vectorize(lambda y: max(y, 0.0))
        square = np.vectorize(lambda y: y**2)

        norm1 = square(prox_p(x-a*Fx)).sum()
        norm2 = square(x).sum()
        norm3 = square(prox_p(Fx-a*x)).sum()
        norm4 = square(Fx).sum()

        return Fx@x + (1/(2*a))*(norm1 - norm2 + norm3 - norm4)


    def DMngMerit(self, x):
        a = 1.5
        Fx = self.F(x)

        term1 = np.vectorize(lambda y: max(y, 0.0))(x-a*Fx)
        term2 = np.vectorize(lambda y: max(y, 0.0))(Fx-a*x)

        return (1/a)*((- a*self.M + np.eye(self.N))@(term1-x) + (self.M - a*np.eye(self.N))@(term2-Fx))

    # Mangasarian and Solodov
    def solve_M(self, x0=None, L=0.001, err=10**(-6), max_itr=100):
        dimension = self.N
        if type(x0) == type(None):
            x = self.calc_x_init()
        else:
            x = x0
            
        x_list = [x]


        lambda_prev = 0
        lambda_curr = 1
        gamma = 1
        x_prev = x
        alpha = 0.05 / (2 * L)
        
        L_prev = L

        # Set initial gradient
        y = x_prev
        gradient = self.DMngMerit(y)
        
        
        def calc_L(x_prev, y_prev, L_prev, lambda_cur, lambda_prev, gradient_prev):
            i = 0
            eta = 1.1 # >1
            gamma_tmp = (1 - lambda_prev) / lambda_curr
            
            while True:
                L_tmp = (L_prev)*(eta**(i))
                alpha_tmp = 0.05 / (2 * L_tmp)
                x_curr_tmp = y_prev - alpha_tmp * gradient_prev
                
                F = self.MngMerit(x_curr_tmp)
                Q = self.MngMerit(y_prev) + (x_curr_tmp - y_prev)@self.DMngMerit(y_prev) + (L_tmp/(0.05))*(np.linalg.norm(x_curr_tmp - y_prev))**2
                
                
                if F<= Q:
                    break
                else:
                    i = i+1
                    
            # if i != 0:
            #     print(i)
            #     print('F=', F)
            #     print('Q=', Q)
            return L_tmp
                

        for i in range(max_itr):
            
            
            L_curr = calc_L(x_prev, y, L_prev, lambda_curr, lambda_prev, gradient)
            
            
            alpha = 0.05/(2*L_curr)
            x_curr = y - alpha * gradient
            y = (1 - gamma) * x_curr + gamma * x_prev
            x_prev = x_curr
            
            
            lambda_curr = (1 + math.sqrt(1 + 4 * lambda_prev * lambda_prev)) / 2
            gamma = (1 - lambda_prev) / lambda_curr
            lambda_prev = lambda_curr
    

            gradient = self.DMngMerit(y)
        
            
            
            if self.MngMerit(x_curr) <= err:
                x_list = x_list + [x_curr]
                break
            
            # Restart
            if self.MngMerit(x_list[-1]) < self.MngMerit(x_curr):
                lambda_prev = 0
                print('Restart')
            
            x_list = x_list + [x_curr]
            L_prev = L_curr

        return x_curr, x_list


    # FISTA
    def solve_F(self, x0=None, L=0.001, err=10**(-6), max_itr=100):

        dimension = self.N
        if type(x0) == type(None):
            x = self.calc_x_init()
        else:
            x = x0
            
        x_list = [x]


        lambda_prev = 0
        lambda_curr = 1
        gamma = 1
        x_prev = x
        alpha = 0.05 / (2 * L)
        
        L_prev = L

        # Set initial gradient
        y = x_prev
        gradient = self.DPsi(y)
        
        
        def calc_L(x_prev, y_prev, L_prev, lambda_cur, lambda_prev, gradient_prev):
            i = 0
            eta = 1.1 # >1
            gamma_tmp = (1 - lambda_prev) / lambda_curr
            
            while True:
                L_tmp = (L_prev)*(eta**(i))
                alpha_tmp = 0.05 / (2 * L_tmp)
                x_curr_tmp = y_prev - alpha_tmp * gradient_prev
                
                F = self.Psi(x_curr_tmp)
                Q = self.Psi(y_prev) + (x_curr_tmp - y_prev)@self.DPsi(y_prev) + (L_tmp/(0.05))*(np.linalg.norm(x_curr_tmp - y_prev))**2
                
                
                if F<= Q:
                    break
                else:
                    i = i+1
                    
            # if i != 0:
            #     print(i)
            #     print('F=', F)
            #     print('Q=', Q)
            return L_tmp
                

        for i in range(max_itr):
            
            
            L_curr = calc_L(x_prev, y, L_prev, lambda_curr, lambda_prev, gradient)
            
            
            alpha = 0.05/(2*L_curr)
            x_curr = y - alpha * gradient
            y = (1 - gamma) * x_curr + gamma * x_prev
            x_prev = x_curr
            
            
            lambda_curr = (1 + math.sqrt(1 + 4 * lambda_prev * lambda_prev)) / 2
            gamma = (1 - lambda_prev) / lambda_curr
            lambda_prev = lambda_curr
    

            gradient = self.DPsi(y)
        
        
            
            
            if self.Psi(x_curr) <= err:
                x_list = x_list + [x_curr]
                break
            
            # Restart
            if self.Psi(x_list[-1]) < self.Psi(x_curr):
                lambda_prev = 0
                print('Restart')
            
            x_list = x_list + [x_curr]
            L_prev = L_curr

        return x_curr, x_list
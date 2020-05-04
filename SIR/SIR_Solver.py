import math

T_ini = 100
dt_ini = 1
S_ini = 0
I_ini = 0
R_ini = 0

class Prameter:
    def __init__(self, R0, gamma):
        self.R0 = R0
        self.gamma = gamma
        self.beta = self.R0 * self.gamma
        


class Model:
    def __init__(self, prm, S=S_ini, I=I_ini, R=R_ini, T=T_ini, dt=dt_ini):
        self.prm  = prm
        
        self.I = float(I)
        self.S = 1-self.I
        self.R = float(0)
        
        # Time step
        self.T = T
        self.dt = dt
        
        # Number of time steps
        self.TN = int(self.T/self.dt)
        
        # List of solutions
        self.S_list = []
        self.S_list.append(self.S)
        self.I_list = [] 
        self.I_list.append(self.I)
        self.R_list = []
        self.R_list.append(self.R)
    
    # Solving functions =============================
    # Analytical solving as a discreate time model===
    def solve_Anly(self):
        
        for k in range(self.TN):
            
            New_S = self.Anly_S()
            New_I = self.Anly_I()
            New_R = self.Anly_R()
            self.Update(New_S, New_I, New_R)
            
    # Euler method===================================
    def solve_Euler(self):
        for k in range(self.TN):
            # Variable dict
            Now = {'S':self.S, 'I':self.I, 'R':self.R}
            New_S = self.Euler_S(Now)
            New_I = self.Euler_I(Now)
            New_R = self.Euler_R(Now)
            self.Update(New_S, New_I, New_R)
    
    # Heum method===================================
    def solve_Heun(self):
        for k in range(self.TN):
            Now = {'S':self.S, 'I':self.I, 'R':self.R}
            Pre = {'S':self.Euler_S(Now), 'I':self.Euler_I(Now), 'R':self.Euler_R(Now)}
            
            New_S = self.Heun_S(Now, Pre)
            New_I = self.Heun_I(Now, Pre)
            New_R = self.Heun_R(Now, Pre)
            self.Update(New_S, New_I, New_R)
    
    # Runge-Kutta method (RK4)=======================
    def solve_Runge_Kutta(self):
        for k in range(self.TN):
            Now = {'S':self.S, 'I':self.I, 'R':self.R}    
            K1   = {'dS':self.S_dynamics(Now), 'dI':self.I_dynamics(Now), 'dR':self.R_dynamics(Now)}
            
            P2 = self.d_sum(Now, self.d_times(K1,self.dt*0.5))
            K2 = {'dS':self.S_dynamics(P2), 'dI':self.I_dynamics(P2), 'dR':self.R_dynamics(P2)}
            
            P3 = self.d_sum(Now, self.d_times(K2,self.dt*0.5))
            K3 = {'dS':self.S_dynamics(P3), 'dI':self.I_dynamics(P3), 'dR':self.R_dynamics(P3)}
            
            P4 = self.d_sum(Now, self.d_times(K2,self.dt))
            K4 = {'dS':self.S_dynamics(P4), 'dI':self.I_dynamics(P4), 'dR':self.R_dynamics(P4)}
            
            New_S = self.Runge_Kutta_S(K1,K2,K3,K4)
            New_I = self.Runge_Kutta_I(K1,K2,K3,K4)
            New_R = self.Runge_Kutta_R(K1,K2,K3,K4)
            self.Update(New_S, New_I, New_R)
    
    # Calculate functions============================
    # Function for multiplication of dict-elements===
    def d_times(self, dic, c):
        rdic = {}
        rdic['dS'] = dic['dS']*c
        rdic['dI'] = dic['dI']*c 
        rdic['dR'] = dic['dR']*c 
        return rdic
    
    # Function for summention of dict-elements========
    def d_sum(self, dic1, dic2):
        rdic = {}
        rdic['S'] = dic1['S'] + dic2['dS']
        rdic['I'] = dic1['I'] + dic2['dI'] 
        rdic['R'] = dic1['R'] + dic2['dR'] 
        return rdic


    # Calculate dx/dt based on differential equiation=========
    def S_dynamics(self, v):
        return - self.prm.beta*v['S']*v['I']  
    def I_dynamics(self, v):
        return (self.prm.beta*v['S']-self.prm.gamma)*v['I']
    def R_dynamics(self, v):
        return v['I']*self.prm.gamma

    
    # Calculate valus of next time step by Runge_Kutta method=====
    def Runge_Kutta_S(self, K1,K2,K3,K4):
        return self.S + (self.dt/6)*(K1['dS']+2*K2['dS']+2*K3['dS']+K4['dS'])
    def Runge_Kutta_I(self, K1,K2,K3,K4):
        return self.I + (self.dt/6)*(K1['dI']+2*K2['dI']+2*K3['dI']+K4['dI'])
    def Runge_Kutta_R(self, K1,K2,K3,K4):
        return self.R + (self.dt/6)*(K1['dR']+2*K2['dR']+2*K3['dR']+K4['dR'])
 

    # Calculate valus of next time step by Heun method============
    def Heun_S(self, Now, Pre):
        return self.S + (self.dt/2)*(self.S_dynamics(Now)+self.S_dynamics(Pre))
    def Heun_I(self, Now, Pre):
        return self.I + (self.dt/2)*(self.I_dynamics(Now)+self.I_dynamics(Pre))
    def Heun_R(self, Now, Pre):
        return self.R + (self.dt/2)*(self.R_dynamics(Now)+self.R_dynamics(Pre))
    
    # Calculate valus of next time step by Euler method============
    def Euler_S(self, Now):
        return self.S + self.dt*self.S_dynamics(Now)    
    def Euler_I(self, Now):
        return self.I + self.dt*self.I_dynamics(Now)
    def Euler_R(self, Now):
        return self.R + self.dt*self.R_dynamics(Now)

    
    # Calculate analytical solution for discrete-time-model based on difference equation=======
    def Anly_S(self):
        New_S \
        = self.S*(math.e**(-self.prm.beta*self.dt*self.I))
        return New_S
        
    def Anly_I(self):
        New_I \
        = self.I\
        + self.S*(1 - math.e**(-self.prm.beta*self.dt*self.I))\
        - self.I*self.dt*self.prm.gamma
        return New_I
        
    def Anly_R(self):
        New_R \
        = self.R\
        + self.I*self.dt* self.prm.gamma
        return New_R
    
    
    # Function for update variable====================
    def Update(self,New_S, New_I, New_R, ):
        self.S_list.append(New_S)
        self.I_list.append(New_I)
        self.R_list.append(New_R)
        self.S = New_S
        self.I = New_I
        self.R = New_R
    
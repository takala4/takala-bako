import numpy as np
import gurobipy as gp

def FB(a,b):
    return np.sqrt(a**2 + b**2) - a - b

def roundF(x,M,q):
    n = len(x)
    A = np.zeros(n*n).reshape(n,n)
    B = np.zeros(n*n).reshape(n,n)
    for i in range(n):
        if x[i]<=10**(-10) and np.dot(M[i],x)+q[i]<=10**(-10):
            A[i][i] = 1/np.sqrt(2) - 1
            B[i][i] = 1/np.sqrt(2) - 1
        else:
            A[i][i] = x[i] / np.sqrt(x[i]**2 + (np.dot(M[i],x) + q[i])**2) - 1
            B[i][i] = (np.dot(M[i],x) + q[i]) / np.sqrt(x[i]**2 + (np.dot(M[i],x) + q[i])**2) - 1
    ret = A + np.dot(M,B)
    return ret

def Phi(x,M,q):
    n = len(x)
    ret = np.zeros(n)
    for i in range(n):
        ret[i] = FB(x[i], np.dot(M[i],x)+q[i])
    return ret

def Psi(x,M,q):
    ret = 0
    n = len(x)
    for i in range(n):
        ret += FB(x[i], np.dot(M[i], x) + q[i])**2
    return ret * 0.5

def nablaofPsi(x,M,q):
    ret = np.dot(roundF(x,M,q), Phi(x,M,q))
    return ret


def solve(M, q, x0=None, err = 10**(-6), max_itr=100):
    #step0
    rho = 10**(-8)
    p = 2.1
    sigma = 10**(-3)
    np.random.seed(1)
    n = len(q)

    if type(x0) == type(None):
        x = calc_x_init(M, q)
    else:
        x = x0
            

    vecnabPhi = [] #グラフ用
    vecnabPhi.append(np.linalg.norm(nablaofPsi(x,M,q)))
    for k in range(max_itr):
        #---step1---
        nabPsi = nablaofPsi(x,M,q)
        #終了判定
        if np.linalg.norm(nabPsi) <= err:
            print('the optimum x is found!')
            print('nablaPsi(x) = ', np.linalg.norm(nablaofPsi(x,M,q)))
            return x, vecnabPhi
        #---step2---
        V = roundF(x,M,q)
        #Vが正則行列かチェック
        det_V = np.linalg.det(V)
        if det_V == 0:
            d = -nablaofPsi(x,M,q)
        else:
            d = -np.dot(np.linalg.inv(V), Phi(x,M,q)) #本当はガウスの消去法を使うべき
            if np.dot(nabPsi,d) > -rho*np.linalg.norm(d)**p:
                d = -nablaofPsi(x,M,q)
        #---step3---
        ik = 0
        while True:
            term1 = Psi(x+d/(2**ik),M,q)
            term2 = Psi(x,M,q)
            term3 = np.dot(nablaofPsi(x,M,q), d) * sigma / (2**ik)
            if term1 <= term2 + term3:
                x = x + 2**(-ik) * d #暫定解の更新
                vecnabPhi.append(np.linalg.norm(nablaofPsi(x,M,q))) #グラフ用
                break
            else:
                ik += 1
        # if k%5 == 0: #5回ごとに出力
        #     print('k = ',k, ' : ', Psi(x,M,q), ' : ', np.linalg.norm(nablaofPsi(x,M,q)))
    return x,vecnabPhi #itermax以内に解が得られなかった場合



def calc_x_init(M, b):
    N = M.shape[0]
    diff_Matrix = np.eye(N) - np.roll(np.eye(N), 1, axis=1)
    
    model = gp.Model()
    
    x_ini = model.addMVar(N, lb=0.0)
    x_a = model.addMVar(N, lb=0.0)
        
    model.setObjective(np.ones(N)@x_a)
    model.addConstr(M@x_ini >= - b - np.eye(N)@x_a)
    model.addConstr(diff_Matrix@x_a == np.zeros(N))
    
    model.setParam('OutputFlag', 0) # Mute putput

    model.optimize()
        
    return x_ini.X
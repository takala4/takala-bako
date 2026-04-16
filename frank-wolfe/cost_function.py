import numpy as np


def bpr_cost(flow, fftt, capacity, alpha, beta):
    """BPR関数によるリンク旅行時間の計算
    t = fftt * (1 + alpha * (flow / capacity)^beta)
    """
    return fftt * (1.0 + alpha * (flow / capacity) ** beta)


def bpr_cost_derivative(flow, fftt, capacity, alpha, beta):
    """BPR関数の導関数 dt/df"""
    return fftt * alpha * beta * (flow / capacity) ** (beta - 1.0) / capacity


def bpr_cost_integral(flow, fftt, capacity, alpha, beta):
    """BPR関数の積分 (Beckmann目的関数の計算用)
    integral_0^f t(x) dx = fftt * (f + alpha * capacity / (beta+1) * (f/capacity)^(beta+1))
    """
    return fftt * (flow + alpha * capacity / (beta + 1.0) * (flow / capacity) ** (beta + 1.0))

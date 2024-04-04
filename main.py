import matplotlib.pyplot as plt
import math
from sympy import Symbol

def getTangentCirlce(x1, y1, r1, x2, y2, r2, x3, y3, r3, s1, s2, s3):
    # my crappy version
    # calculate x and y in relation to r
    r = Symbol("r")

    #region general arguments
    y = ((2 * x1) * y2)
    z = ((2 * x1) * y3)
    aa = ((2 * x2) * y1)
    ba = ((2 * x2) * y3)
    ca = ((2 * x3) * y1)
    da = ((2 * x3) * y2)
    #endregion

    #region x arguments
    a = ((((2 * r) * r1) * s1) * y2)
    b = ((((2 * r) * r1) * s1) * y3)
    c = ((((2 * r) * r2) * s2) * y1)
    d = ((((2 * r) * r2) * s2) * y3)
    e = ((((2 * r) * r3) * s3) * y1)
    f = ((((2 * r) * r3) * s3) * y2)
    g = ((r1**2 * s1**2) * y2)
    h = ((r1**2 * s1**2) * y3)
    i = ((r2**2 * s2**2) * y1)
    j = ((r2**2 * s2**2) * y3)
    k = ((r3**2 * s3**2) * y1)
    l = ((r3**2 * s3**2) * y2)
    m = (x1**2 * y2)
    n = (x1**2 * y3)
    o = (x2**2 * y1)
    p = (x2**2 * y3)
    q = (x3**2 * y1)
    r = (x3**2 * y2)
    s = (y1**2 * y2)
    t = (y1**2 * y3)
    u = (y1 * y2**2)
    v = (y1 * y3**2)
    w = (y2**2 * y3)
    x = (y2 * y3**2)
    #endregion

    x_s = (a - b - c + d + e - f - g + h + i - j - k + l + m - n - o + p + q - r + s - t - u + v + w - x) / (y - z - aa + ba + ca - da)

    r = Symbol("r")

    #region y arguments
    a_2 = ((((2 * r) * r1) * s1) * x2)
    b_2 = ((((2 * r) * r1) * s1) * x3)
    c_2 = ((((2 * r) * r2) * s2) * x1)
    d_2 = ((((2 * r) * r2) * s2) * x3)
    e_2 = ((((2 * r) * r3) * s3) * x1)
    f_2 = ((((2 * r) * r3) * s3) * x2)
    g_2 = ((r1**2 * s1**2) * x2)
    h_2 = ((r1**2 * s1**2) * x3)
    i_2 = ((r2**2 * s2**2) * x1)
    j_2 = ((r2**2 * s2**2) * x3)
    k_2 = ((r3**2 * s3**2) * x1)
    l_2 = ((r3**2 * s3**2) * x2)
    m_2 = (x1**2 * x2)
    n_2 = (x1**2 * x3)
    o_2 = (x1 * x2**2)
    p_2 = (x1 * x3**2)
    q_2 = (x1 * y2**2)
    r_2 = (x1 * y3**2)
    s_2 = (x3 * x2**2)
    t_2 = (x2 * x3**2)
    u_2 = (x2 * y1**2)
    v_2 = (x2 * y3**2)
    w_2 = (x3 * y1**2)
    x_2 = (x3 * y2**2)
    #endregion

    y_s = (-a_2 + b_2 + c_2 - d_2 - e_2 + f_2 + g_2 - h_2 - i_2 + j_2 + k_2 - l_2 - m_2 + n_2 + o_2 - p_2 + q_2 - r_2 - s_2 + t_2 - u_2 + v_2 + w_2 - x_2) / (y - z - aa + ba + ca - da)

    sx_s = str(x_s)
    sy_s = str(y_s)

    # get P, Q, M, N out of the equations for x_s and y_s
    y_sa = sy_s.split(" ")
    if y_sa[0].find("r")!=-1: 
        Q = float(y_sa[0].split("*")[0])
        P = float(y_sa[2])
        if y_sa[1] == "-":
            P *= -1
    else: 
        P = float(y_sa[0])
        Q = float(y_sa[2].split("*")[0])
        if y_sa[1] == "-":
            Q *= -1

    x_sa = sx_s.split(" ")
    if x_sa[0].find("r")!=-1: 
        N = float(x_sa[0].split("*")[0])
        M = float(x_sa[2])
        if x_sa[1] == "-":
            M *= -1
    else: 
        M = float(x_sa[0])
        N = float(x_sa[2].split("*")[0])
        if x_sa[1] == "-":
            N *= -1
    
    # copied from below
    a = N*N + Q*Q - 1
    b = 2*M*N - 2*N*x1 + 2*P*Q - 2*Q*y1 + 2*s1*r1
    c = x1*x1 + M*M - 2*M*x1 + P*P + y1*y1 - 2*P*y1 - r1*r1
 
    # Find a root of a quadratic equation. This requires the circle centers not to be e.g. colinear
    D = b*b-4*a*c
    r_s = (-b-math.sqrt(D))/(2*a)
   
    x_s = M+N*r_s
    y_s = P+Q*r_s
    
    return x_s, y_s, r_s

def solve2(x1, y1, r1, x2, y2, r2, x3, y3, r3, s1, s2, s3):
    # from here https://rosettacode.org/wiki/Problem_of_Apollonius
    v11 = 2*x2 - 2*x1
    v12 = 2*y2 - 2*y1
    v13 = x1*x1 - x2*x2 + y1*y1 - y2*y2 - r1*r1 + r2*r2
    v14 = 2*s2*r2 - 2*s1*r1
 
    v21 = 2*x3 - 2*x2
    v22 = 2*y3 - 2*y2
    v23 = x2*x2 - x3*x3 + y2*y2 - y3*y3 - r2*r2 + r3*r3
    v24 = 2*s3*r3 - 2*s2*r2
 
    w12 = v12/v11
    w13 = v13/v11
    w14 = v14/v11
 
    w22 = v22/v21-w12
    w23 = v23/v21-w13
    w24 = v24/v21-w14
 
    P = -w23/w22
    Q = w24/w22
    M = -w12*P-w13
    N = w14 - w12*Q
    
    a = N*N + Q*Q - 1
    b = 2*M*N - 2*N*x1 + 2*P*Q - 2*Q*y1 + 2*s1*r1
    c = x1*x1 + M*M - 2*M*x1 + P*P + y1*y1 - 2*P*y1 - r1*r1
 
    # Find a root of a quadratic equation. This requires the circle centers not to be e.g. colinear
    D = b*b-4*a*c
    try:
        rs = (-b-math.sqrt(D))/(2*a)
    except Exception as e:
        return 0, 0, 0
 
    xs = M+N*rs
    ys = P+Q*rs

    return xs, ys, rs

if __name__ == "__main__":
    x1 = 30
    y1 = 3
    r1 = 3
    x2 = 23
    y2 = 12
    r2 = 6
    x3 = 14
    y3 = 0
    r3 = 10

    # x1 = input("x1: ")
    # y1 = input("y1: ")
    # r1 = input("r1: ")
    # x2 = input("x2: ")
    # y2 = input("y2: ")
    # r2 = input("r2: ")
    # x3 = input("x3: ")
    # y3 = input("y3: ")
    # r3 = input("r3: ")

    # convert variables to float
    x1 = float(x1)
    y1 = float(y1)
    r1 = float(r1)
    x2 = float(x2)
    y2 = float(y2)
    r2 = float(r2)
    x3 = float(x3)
    y3 = float(y3)
    r3 = float(r3)

    # create the figure
    fig, ax = plt.subplots()
    ax.set(xlim=[(min(x1, x2, x3) - max(r1, r2, r3)) - min(r1, r2, r3),
                 (max(x1, x2, x3) + max(r1, r2, r3)) + min(r1, r2, r3)],
            ylim=[(min(y1, y2, y3) - max(r1, r2, r3)) - min(r1, r2, r3),
                (max(y1, y2, y3) + max(r1, r2, r3)) + min(r1, r2, r3)],
            aspect=1)

    # store the cicle coordinates in a list
    coordinates = [(x1, y2), (x2, y2), (x3, y3)]
    scales = [r1, r2, r3]

    for i in range(3):
        circle = plt.Circle(coordinates[i], scales[i], color='black', fill=False)
        ax.add_patch(circle)

    # circle = plt.Circle((.5245460440128689, .3584079119742619), 0.09546044012868957, color='red', fill=False)
    # ax.add_patch(circle)

    for i in range(-1, 2, 2):
        for o in range(-1, 2, 2):
            for p in range(-1, 2, 2):
                try:
                    x, y, r = getTangentCirlce(
                        coordinates[0][0], coordinates[0][1], scales[0],
                        coordinates[1][0], coordinates[1][1], scales[1],
                        coordinates[2][0], coordinates[2][1], scales[2],
                        i, o, p)
                
                    print("Circle:", i,  o,  p, "; x =", x, ", y =", y, ", r =", r)
                    circle = plt.Circle((x, y), r, color='red', fill=False)
                    ax.add_patch(circle)
                except Exception as e:
                    print("Circle:", i,  o,  p, "; not possible")


    fig.savefig('plotcircles.png')
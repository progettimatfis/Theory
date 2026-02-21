import numpy as np
import sympy as sp

#################### LETTURA FILE PROFILO ALA ####################
def read_airfoil_file(filename):
    data = []
    with open(filename, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                try:
                    x = sp.sympify(parts[0])
                    z = sp.sympify(parts[1])
                    data.append((x, z))
                except (sp.SympifyError, ValueError):
                    continue
    return data

#################### MATRICI AERODINAMICHE ####################
def Matrici_aerodinamiche(airfoil_data,alpha_sym):
    x = sp.Matrix([p[0] for p in airfoil_data])
    z = sp.Matrix([p[1] for p in airfoil_data])
    x0 = x[:-1, :]
    x1 = x[1:, :]
    z0 = z[:-1, :]
    z1 = z[1:, :]
    xc = (x0 + x1) / 2
    zc = (z0 + z1) / 2
    dx2 = (x1 - x0).multiply_elementwise(x1 - x0)
    dz2 = (z1 - z0).multiply_elementwise(z1 - z0)
    elle = (dx2 + dz2).applyfunc(sp.sqrt)
    N = len(xc)
    sinb = sp.Matrix([ (z1[i] - z0[i]) / elle[i] for i in range(N) ])
    cosb = sp.Matrix([ (x1[i] - x0[i]) / elle[i] for i in range(N) ])
    ra, rb, w = sp.zeros(N, N), sp.zeros(N, N), sp.zeros(N, N)
    for i in range(N):
        for j in range(N):
            ra[i, j] = (xc[i] - x0[j])**2 + (zc[i] - z0[j])**2
            rb[i, j] = (xc[i] - x1[j])**2 + (zc[i] - z1[j])**2
            if i == j:
                w[i, j] = sp.Rational(1, 2)
            else:
                num = (zc[i] - z1[j]) * (xc[i] - x0[j]) - (zc[i] - z0[j]) * (xc[i] - x1[j])
                den = (xc[i] - x1[j]) * (xc[i] - x0[j]) + (zc[i] - z0[j]) * (zc[i] - z1[j])
                w[i, j] = sp.atan2(num, den) / (2 * sp.pi)
    u = (ra.applyfunc(sp.log) - rb.applyfunc(sp.log)) / (4 * sp.pi)
    cosbibj = sp.Matrix([[cosb[i]*cosb[j] + sinb[i]*sinb[j] for j in range(N)] for i in range(N)])
    sinbibj = sp.Matrix([[sinb[i]*cosb[j] - cosb[i]*sinb[j] for j in range(N)] for i in range(N)])
    vt = cosbibj.multiply_elementwise(u) + sinbibj.multiply_elementwise(w)
    vn = -sinbibj.multiply_elementwise(u) + cosbibj.multiply_elementwise(w)
    vtv = sp.Matrix([sum(vn.row(i)) for i in range(N)])
    vnv = sp.Matrix([-sum(vt.row(i)) for i in range(N)])
    vn_ext, vt_ext = sp.zeros(N + 1, N + 1), sp.zeros(N, N + 1)
    for i in range(N):
        for j in range(N):
            vn_ext[i, j] = vn[i, j]
            vt_ext[i, j] = vt[i, j]
        vn_ext[i, -1] = vnv[i]
        vt_ext[i, -1] = vtv[i]
    for j in range(N):
        vn_ext[-1, j] = vt[0, j] + vt[-1, j]
    vn_ext[-1, -1] = sum([vn[0, j] + vn[-1, j] for j in range(N)])
    rhs = sp.zeros(N + 1, 1)
    rhs[:N, 0] = sp.cos(alpha_sym) * sinb - sp.sin(alpha_sym) * cosb
    rhs[-1, 0] = 0
    return vn_ext, rhs

#################### SOLUZIONE SISTEMA ####################
def solve_ab(a_val, b_val,vn_ext,rhs):
    a_sym, b_sym = sp.symbols('a b')
    vn_num = vn_ext.subs({a_sym: a_val, b_sym: b_val}).evalf()
    rhs_num = rhs.subs({a_sym: a_val, b_sym: b_val}).evalf()
    vn_np = np.array(vn_num.tolist()).astype(np.float64)
    rhs_np = np.array(rhs_num.tolist()).astype(np.float64)
    siga_np = np.linalg.solve(vn_np, rhs_np)
    return siga_np
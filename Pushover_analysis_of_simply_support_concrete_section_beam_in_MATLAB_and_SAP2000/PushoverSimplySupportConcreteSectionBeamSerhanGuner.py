"""
%********************************************************************%
%                        >> IN THE NAME OF ALLAH <<                  %
% Linear and Nonlinear Analysis of  simply support Beam subjected to % 
% central concentrated load by plastic hinge concept                 %
% Checking the analysis by Force Control                             %
% Newton-Raphson Method : Initial stiffness procedure                %
%--------------------------------------------------------------------%
%          This program is written by salar delavar Qashqai          %  
%               E-mail:salar.d.ghashghaei@gmail.com                  %
%--------------------------------------------------------------------%
%                        Force  Force                                %
%                          |      |                                  %
%                          v      v                                  %
%  simply support  ||--------------------||                          %
%                   |<-      2L,EI     ->|                           %
%********************************************************************%
"""
import numpy as np
import time

# =============================================================================
# Define Parameters in mm, kN
# =============================================================================
P4 = 0            # [kN]
P5 = 0.2          # [kN]
P6 = 0            # [kN.mm]
L = 3000          # [mm]
a1 = 2500
b1 = L - a1
a2 = b1
b2 = L - a2
U5max = 100       # ultimate displacement [mm] [DOF(5)]

# =============================================================================
# Section Properties
# =============================================================================
b = 500           # [mm]
h = 500           # [mm] Height of Section
Ie = (b * h**3) / 12
A_sec = b * h
EI = 25 * Ie      # [kN.mm^2]
EA = 25 * A_sec   # [kN]
m = 1000        # number of calculation increments
itermax = 5000  # maximum number of iterations per increment
tolerance = 1e-6  # specified tolerance for convergence

# Initial guess for displacements (6 DOF)
u = np.zeros(6)
lanX = 1
lanY = 0

# Monitor cpu time
starttime = time.process_time()

# =============================================================================
# Element stiffness coefficients
# =============================================================================
A_coef = 4 * EI / L
B_coef = 6 * EI / L**2
C_coef = 2 * EI / L
D_coef = 12 * EI / L**3
G_coef = EA / L

# =============================================================================
# Moment Curvature Rotation parameters
# =============================================================================
fy = 400  # [N/mm^2] Yield strength of reinforcing steel
Lp = max((0.08 * L + 0.022 * 25 * fy), 0.044 * 25 * fy)  # Plastic hinge length
fiy = 0.0289e-4   # [1/mm] Yield Curvature
fiu = 0.6211e-4   # [1/mm] Ultimate Curvature
ty = fiy * Lp     # [rad] Yield rotation
My = 375063.8     # [kN.mm] Yield Moment
tu = fiu * Lp     # [rad] Ultimate Rotation
Mu = 440578.8     # [kN.mm] Ultimate Moment
k1 = (My / ty) * 1e9  # Note: 10^9 is written as 1e9 in Python
k2 = (Mu - My) / (tu - ty)

# =============================================================================
# Initial internal force for each element
# =============================================================================
Fele1 = np.zeros(6)
Fele2 = np.zeros(6)

# =============================================================================
# Preallocate lists to store results for each increment
# =============================================================================
F3_list = []
U1_list = []
U2_list = []
U3_list = []
U4_list = []
U5_list = []
U6_list = []
INT_N_f2_list = []
INT_N_f3_list = []
INT_N_f5_list = []
INT_N_f6_list = []
INT_N_f7_list = []
INT_N_f9_list = []
INT_N_f10_list = []
TBSN_list = []

# =============================================================================
# Loop over load increments
# =============================================================================
for i in range(m):
    # In MATLAB, the loop runs i = 1,...,m. Here we use i+1.
    load_factor = i + 1
    Pi5 = P5 * load_factor
    Va1 = Pi5 * b1**2 * (1 + 2 * a1 / L) / L**2
    Va2 = Pi5 * a1**2 * (1 + 2 * b1 / L) / L**2
    Ma1 = Pi5 * a1 * b1**2 / L**2
    Ma2 = -Pi5 * b1 * a1**2 / L**2
    Vb1 = Pi5 * b2**2 * (1 + 2 * a2 / L) / L**2
    Vb2 = Pi5 * a2**2 * (1 + 2 * b2 / L) / L**2
    Mb1 = Pi5 * a2 * b2**2 / L**2
    Mb2 = -Pi5 * b2 * a2**2 / L**2

    # Assemble applied force vector (note: Python is 0-indexed)
    F_vec = np.array([-Ma1, P4, -Va2 - Vb1, -Ma2, -Mb1, -Mb2])

    it = 0
    residual = 100.0

    # Determine stiffness parameters k11 and k22 based on internal forces
    # (MATLAB: Fele1(3) -> Python: Fele1[2], and Fele1(6) -> Fele1[5])
    if (abs(Fele1[2]) <= My) and (abs(Fele1[5]) <= My):
        k11 = k1
        k22 = k1
    elif abs(Fele1[5]) >= My:
        k11 = (My + (abs(u[3]) - ty) * k2) / abs(u[3]) if abs(u[3]) != 0 else k1
        k22 = (My + (abs(u[4]) - ty) * k2) / abs(u[4]) if abs(u[4]) != 0 else k1
    else:
        k11 = k1
        k22 = k1

    # Assemble the initial global stiffness matrix (Kini)
    Kini = np.array([
        [A_coef,           B_coef * lanY,      -B_coef * lanX,   C_coef,           0,               0],
        [B_coef * lanY, 2 * (G_coef * lanX**2 + D_coef * lanY**2),
         2 * (G_coef - D_coef) * lanX * lanY,   B_coef * lanY, -B_coef * lanY, -B_coef * lanY],
        [-B_coef * lanX, 2 * (G_coef - D_coef) * lanX * lanY,
         2 * (G_coef * lanY**2 + D_coef * lanX**2), -B_coef * lanX,  B_coef * lanX,  B_coef * lanX],
        [C_coef,           B_coef * lanY,      -B_coef * lanX,   A_coef + k11,      0,               0],
        [0,              -B_coef * lanY,       B_coef * lanX,   0,         A_coef + k22,       C_coef],
        [0,              -B_coef * lanY,       B_coef * lanX,   0,              C_coef,         A_coef]
    ])

    # ---------------------------
    # Newton-Raphson Iteration
    # ---------------------------
    while residual > tolerance:
        if (abs(Fele1[2]) <= My) and (abs(Fele1[5]) <= My):
            k11 = k1
            k22 = k1
        elif abs(Fele1[5]) >= My:
            k11 = (My + (abs(u[3]) - ty) * k2) / abs(u[3]) if abs(u[3]) != 0 else k1
            k22 = (My + (abs(u[4]) - ty) * k2) / abs(u[4]) if abs(u[4]) != 0 else k1
        else:
            k11 = k1
            k22 = k1

        # Assemble the updated global stiffness matrix K
        K = np.array([
            [A_coef,           B_coef * lanY,      -B_coef * lanX,   C_coef,           0,               0],
            [B_coef * lanY, 2 * (G_coef * lanX**2 + D_coef * lanY**2),
             2 * (G_coef - D_coef) * lanX * lanY,   B_coef * lanY, -B_coef * lanY, -B_coef * lanY],
            [-B_coef * lanX, 2 * (G_coef - D_coef) * lanX * lanY,
             2 * (G_coef * lanY**2 + D_coef * lanX**2), -B_coef * lanX,  B_coef * lanX,  B_coef * lanX],
            [C_coef,           B_coef * lanY,      -B_coef * lanX,   A_coef + k11,      0,               0],
            [0,              -B_coef * lanY,       B_coef * lanX,   0,         A_coef + k22,       C_coef],
            [0,              -B_coef * lanY,       B_coef * lanX,   0,              C_coef,         A_coef]
        ])

        # Compute residual force vector f = K*u - F_vec
        f = np.dot(K, u) - F_vec

        # Solve for displacement increment du using the initial stiffness matrix Kini
        try:
            du = np.linalg.solve(Kini, -f)
        except np.linalg.LinAlgError:
            print("Matrix Kini is singular at increment", i+1)
            break

        residual = np.max(np.abs(du))
        it += 1
        if it == itermax:
            print("(-) For increment {} trial iteration reached Ultimate in {} iterations".format(i+1, it))
            print("    ## The solution for this step is not converged ##")
            break

        u = u + du

    if it < itermax:
        print("(+) It is converged in {} iterations for increment {}".format(it, i+1))

    # =============================================================================
    # Compute Internal Forces for Element 1
    # =============================================================================
    # Displacement Transformation Matrix T
    T = np.array([
        [lanX,  lanY, 0,    0,    0,  0],
        [-lanY, lanX, 0,    0,    0,  0],
        [0,     0,    1,    0,    0,  0],
        [0,     0,    0,  lanX,  lanY, 0],
        [0,     0,    0, -lanY,  lanX, 0],
        [0,     0,    0,    0,    0,  1]
    ])

    # Stiffness Matrix for each element (Kele)
    Kele = np.array([
        [G_coef,      0,      0, -G_coef,      0,      0],
        [0,      D_coef, B_coef,      0, -D_coef,  B_coef],
        [0,      B_coef, A_coef,      0, -B_coef,  C_coef],
        [-G_coef,     0,      0,  G_coef,      0,      0],
        [0,     -D_coef, -B_coef,     0,  D_coef, -B_coef],
        [0,      B_coef,  C_coef,     0, -B_coef,  A_coef]
    ])

    # For element 1, the displacement vector is [0, 0, u[0], u[1], u[2], u[3]]
    vec_elem1 = np.array([0, 0, u[0], u[1], u[2], u[3]])
    Fele1 = np.dot(Kele, np.dot(T, vec_elem1)) + np.array([0, Va1, Ma1, 0, Va2, Ma2])

    # =============================================================================
    # Compute Internal Forces for Element 2
    # =============================================================================
    # For element 2, the displacement vector is [u[1], u[2], u[4], 0, 0, u[5]]
    vec_elem2 = np.array([u[1], u[2], u[4], 0, 0, u[5]])
    Fele2 = np.dot(Kele, np.dot(T, vec_elem2)) + np.array([0, Vb1, Mb1, 0, Vb2, Mb2])

    # =============================================================================
    # Store Force and Displacement for the current increment
    # =============================================================================
    F3_list.append(F_vec[2])
    U1_list.append(u[0])
    U2_list.append(u[1])
    U3_list.append(u[2])
    U4_list.append(u[3])
    U5_list.append(u[4])
    U6_list.append(u[5])
    INT_N_f2_list.append(Fele1[1])
    INT_N_f3_list.append(Fele1[2])
    INT_N_f5_list.append(Fele1[4])
    INT_N_f6_list.append(Fele1[5])
    INT_N_f7_list.append(Fele2[2])
    INT_N_f9_list.append(Fele2[4])
    INT_N_f10_list.append(Fele2[5])
    TBSN_list.append(Fele1[1] + Fele2[4])  # Total Base Shear

    # Check for ultimate displacement (MATLAB uses u(3) which is u[2] in Python)
    if abs(u[2]) >= U5max:
        print("Displacement at middle [dof(5)] reached ultimate displacement")
        break

# =============================================================================
# Postprocessing: Construct Output Arrays
# =============================================================================
F3 = np.array(F3_list)
U1 = np.array(U1_list)
U2 = np.array(U2_list)
U3 = np.array(U3_list)
U4 = np.array(U4_list)
U5 = np.array(U5_list)
U6 = np.array(U6_list)
INT_N_f6 = np.array(INT_N_f6_list)
TBSN = np.array(TBSN_list)

# D1 = [0; U3'] in MATLAB means prepend 0 to U3
D1 = U3
F1 = F3
T1 = U4.copy()
M1 = INT_N_f6.copy()

# =============================================================================
# Display the shape of deflection along the length for Nonlinear Analysis
# =============================================================================
x1_vals = []
y11_vals = []
y12_vals = []
for j in range(11):
    # In MATLAB j runs 1 to 11, and x1 = (L/10)*(j-1); here j starts at 0.
    x_val = (L / 10) * j
    x1_vals.append(x_val)
    N1_val = (1 / L**3) * (L**3 - 3 * L * x_val**2 + 2 * x_val**3)
    N2_val = (1 / L**2) * (x_val * L**2 - 2 * L * x_val**2 + x_val**3)
    N3_val = (1 / L**3) * (3 * L * x_val**2 - 2 * x_val**3)
    N4_val = (1 / L**2) * (-L * x_val**2 + x_val**3)
    # y11 = [N1 N2 N3 N4] * [0; u(1); u(3); u(4)]
    y11_val = N1_val * 0 + N2_val * u[0] + N3_val * u[2] + N4_val * u[3]
    # y12 = [N1 N2 N3 N4] * [u(3); u(5); 0; u(6)]
    y12_val = N1_val * u[2] + N2_val * u[4] + N3_val * 0 + N4_val * u[5]
    y11_vals.append(y11_val)
    y12_vals.append(y12_val)

# Create x12 which is offset by L
x12_vals = [L + x for x in x1_vals]
X1 = np.concatenate((np.array(x1_vals), np.array(x12_vals)))
Y1 = np.concatenate((np.array(y11_vals), np.array(y12_vals)))

# =============================================================================
# Display Final Displacements
# =============================================================================
print("=" * 100)
print("rotation(D3)   X-displacement(D4)   Y-displacement(D5)   rotation(D6)   rotation(D7)   rotation(D10)")
print("-" * 100)
# Assemble final displacement results into a matrix for display.
disp_matrix = np.column_stack((U1, U2, U3, U4, U5, U6))
print(disp_matrix)

"""
import numpy as np

# --- Assumed pre-defined parameters (example values are provided) ---
# (Replace these with your actual values or ensure they are defined in your code.)
P4    = 0          # [kN]
P5    = 0.2        # [kN]
L     = 3000       # [mm]
a1    = 2500       # [mm]
b1    = L - a1     # [mm]
a2    = b1         # [mm]
b2    = L - a2     # [mm]
lanX  = 1
lanY  = 0

# Section properties and stiffness coefficients
b_sec = 500         # [mm]
h_sec = 500         # [mm]
Ie    = (b_sec * h_sec**3) / 12
A_sec = b_sec * h_sec
EI    = 25 * Ie     # [kN.mm^2]
EA    = 25 * A_sec  # [kN]

A_coef = 4 * EI / L
B_coef = 6 * EI / L**2
C_coef = 2 * EI / L
D_coef = 12 * EI / L**3
G_coef = EA / L

# --- Number of load increments for the linear analysis ---
"""

# --- Preallocate lists for storing results ---
F3L_list      = []
U1L_list      = []
U2L_list      = []
U3L_list      = []
U4L_list      = []
U5L_list      = []
INT_L_f2_list = []
INT_L_f3_list = []
INT_L_f5_list = []
INT_L_f6_list = []
INT_L_f7_list = []
INT_L_f9_list = []
INT_L_f10_list= []
TBSL_list     = []

# --- Initial displacement vector (5 DOF) ---
u = np.zeros(5)

# --- Loop over load increments ---
for i in range(m):
    # Define the applied load
    Pi5 = P5 * i
    Va1 = Pi5 * b1**2 * (1 + 2*a1/L) / L**2
    Va2 = Pi5 * a1**2 * (1 + 2*b1/L) / L**2
    Ma1 = Pi5 * a1 * b1**2 / L**2
    Ma2 = -Pi5 * b1 * a1**2 / L**2
    Vb1 = Pi5 * b2**2 * (1 + 2*a2/L) / L**2
    Vb2 = Pi5 * a2**2 * (1 + 2*b2/L) / L**2
    Mb1 = Pi5 * a2 * b2**2 / L**2
    Mb2 = -Pi5 * b2 * a2**2 / L**2

    # Assemble the global force vector F (5×1)
    F_vec = np.array([
        -Ma1,
         P4,
        -Va2 - Vb1,
        -Ma2 - Mb1,
        -Mb2
    ])

    # Assemble the global stiffness matrix K (5×5)
    # MATLAB:
    # K = [A,     B*lanY,        -B*lanX,     C,         0;
    #      B*lanY,2*(G*lanX^2+D*lanY^2),2*(G-D)*lanX*lanY,2*B*lanY, -B*lanY;
    #     -B*lanX,2*(G-D)*lanX*lanY,2*(G*lanY^2+D*lanX^2),0,         B*lanX;
    #      C,     2*B*lanY,       0,          2*A,       C;
    #      0,     -B*lanY,        B*lanX,     C,         A];
    K = np.array([
        [A_coef,              B_coef * lanY,               -B_coef * lanX,      C_coef,         0],
        [B_coef * lanY, 2*(G_coef*lanX**2 + D_coef*lanY**2), 2*(G_coef - D_coef)*lanX*lanY, 2*B_coef*lanY, -B_coef*lanY],
        [-B_coef * lanX, 2*(G_coef - D_coef)*lanX*lanY, 2*(G_coef*lanY**2 + D_coef*lanX**2),     0,         B_coef*lanX],
        [C_coef,         2*B_coef*lanY,                         0,              2*A_coef,     C_coef],
        [0,              -B_coef*lanY,              B_coef*lanX,      C_coef,         A_coef]
    ])

    # Solve for displacements: u = inv(K)*F (using a linear solver)
    u = np.linalg.solve(K, F_vec)

    # --- Internal force calculations using the element stiffness ---
    # Displacement Transformation Matrix T (6×6)
    T = np.array([
        [lanX,   lanY,  0,    0,    0,  0],
        [-lanY,  lanX,  0,    0,    0,  0],
        [0,      0,     1,    0,    0,  0],
        [0,      0,     0,  lanX,  lanY,  0],
        [0,      0,     0, -lanY,  lanX,  0],
        [0,      0,     0,    0,    0,  1]
    ])

    # Element stiffness matrix Kele (6×6)
    Kele = np.array([
        [G_coef,   0,      0,   -G_coef,   0,      0],
        [0,        D_coef, B_coef,  0,      -D_coef,  B_coef],
        [0,        B_coef, A_coef,  0,      -B_coef,  C_coef],
        [-G_coef,  0,      0,    G_coef,    0,      0],
        [0,       -D_coef,-B_coef,  0,       D_coef, -B_coef],
        [0,        B_coef, C_coef,  0,      -B_coef,  A_coef]
    ])

    # For element 1, the displacement vector is [0; 0; u(1); u(2); u(3); u(4)]
    # (Remember: MATLAB u(1:4) -> Python u[0:4])
    vec_elem1 = np.array([0, 0, u[0], u[1], u[2], u[3]])
    Fele1L = Kele @ T @ vec_elem1 + np.array([0, Va1, Ma1, 0, Va2, Ma2])

    # For element 2, the displacement vector is [u(2); u(3); u(4); 0; 0; u(5)]
    # (MATLAB u(2:4) -> Python u[1:4] and u(5) -> u[4])
    vec_elem2 = np.array([u[1], u[2], u[3], 0, 0, u[4]])
    Fele2L = Kele @ T @ vec_elem2 + np.array([0, Vb1, Mb1, 0, Vb2, Mb2])

    # --- Store force and displacement results for each increment ---
    F3L_list.append(F_vec[2])      # MATLAB F(3) corresponds to F_vec[2]
    U1L_list.append(u[0])          # u(1) -> u[0]
    U2L_list.append(u[1])
    U3L_list.append(u[2])
    U4L_list.append(u[3])
    U5L_list.append(u[4])
    INT_L_f2_list.append(Fele1L[1])  # Fele1L(2) -> Fele1L[1]
    INT_L_f3_list.append(Fele1L[2])  # Fele1L(3) -> Fele1L[2]
    INT_L_f5_list.append(Fele1L[4])  # Fele1L(5) -> Fele1L[4]
    INT_L_f6_list.append(Fele1L[5])  # Fele1L(6) -> Fele1L[5]
    INT_L_f7_list.append(Fele2L[2])  # Fele2L(3) -> Fele2L[2]
    INT_L_f9_list.append(Fele2L[4])  # Fele2L(5) -> Fele2L[4]
    INT_L_f10_list.append(Fele2L[5]) # Fele2L(6) -> Fele2L[5]
    TBSL_list.append(Fele1L[1] + Fele2L[4])  # Total Base Shear

# --- Post-process the results (prepend zeros where needed) ---
D2      = np.concatenate(([0], np.array(U3L_list)))
F2      = np.concatenate(([0], np.array(F3L_list)))
T2      = np.array(U4L_list)       # Directly taken from U4L_list
M2      = np.array(INT_L_f6_list)
TBSL = np.concatenate(([0], np.array(TBSL_list)))

# --- Display the deflected shape along the length for Linear Analysis ---
x2_vals = []
y11_vals = []
y12_vals = []
for j in range(11):  # j = 0,...,10 corresponds to MATLAB's j=1:11 with (j-1)
    x_val = (L / 10) * j
    x2_vals.append(x_val)
    N1 = (1 / L**3) * (L**3 - 3 * L * x_val**2 + 2 * x_val**3)
    N2 = (1 / L**2) * (x_val * L**2 - 2 * L * x_val**2 + x_val**3)
    N3 = (1 / L**3) * (3 * L * x_val**2 - 2 * x_val**3)
    N4 = (1 / L**2) * (-L * x_val**2 + x_val**3)
    # y11 = [N1, N2, N3, N4] * [0; u(1); u(3); u(4)]
    y11 = N1 * 0 + N2 * u[0] + N3 * u[2] + N4 * u[3]
    # y12 = [N1, N2, N3, N4] * [u(3); u(4); 0; u(5)]
    y12 = N1 * u[2] + N2 * u[3] + N3 * 0 + N4 * u[4]
    y11_vals.append(y11)
    y12_vals.append(y12)

x22_vals = L + np.array(x2_vals)
# X2 and Y2 are formed by vertically concatenating the two sets of values.
X2 = np.concatenate((np.array(x2_vals).reshape(-1, 1), x22_vals.reshape(-1, 1)), axis=0)
Y2 = np.concatenate((np.array(y11_vals).reshape(-1, 1), np.array(y12_vals).reshape(-1, 1)), axis=0)

# --- Display internal force comparisons ---
print("============== Internal Force ================")
print("+ =============== Nonlinear ================ +")
print("     (f3)       (f6)     (f7)      (f10)     ")
print("----------------------------------------------")
# (Assuming INT_N_f3, INT_N_f6, INT_N_f7, INT_N_f10 were calculated in your nonlinear analysis)
# For example:
# print(np.column_stack((INT_N_f3, INT_N_f6, INT_N_f7, INT_N_f10)))
print("+ ================ Linear ================== +")
print("     (f3)       (f6)     (f7)      (f10)     ")
print("----------------------------------------------")
print(np.column_stack((
    np.array(INT_L_f3_list),
    np.array(INT_L_f6_list),
    np.array(INT_L_f7_list),
    np.array(INT_L_f10_list)
)))
print("==============================================")

# --- SAP2000 - Plastic hinge concept (empty lists) ---
D4 = np.array([0,
-0.952380952,
-5.608465608,
-10.79365079,
-13.96825397,
-16.61375661,
-19.68253968,
-27.93650794,
-35.13227513,
-43.28042328,
-47.61904762,
-48.78306878,
-50.15873016,
])
TBSS  = np.array([0,
32.13183422,
79.07419655,
139.045047,
174.6812537,
197.2889722,
205.1458701,
213.0777851,
215.7860572,
218.508108,
216.8347663,
216.8516069,
213.3992872,
])

# --- ABAQUS - wire element ---
D5 = np.array([
    0,
   -12.9055,
   -19.8656,
   -30.3059,
   -34.221,
   -40.0936,
   -42.2958,
   -43.0856,
   -44.2149,
   -45.9141,
   -48.485,
   -52.3412,
   -58.2205,
   -67.0964,
   -80.2472,
   -92.4742
])
TBSA = np.array([
    0,
    19835.4,
    76412.5,
    161278,
    193103,
    240840,
    258741,
    262143,
    262143,
    262641,
    265409,
    269560,
    271341,
    271341,
    278984,
    281208
]) * 0.001

# --- ABAQUS - solid element ---
D6 = np.array([
    0,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.01,
   -0.0110817,
   -0.0166238,
   -0.0249369,
   -0.0374066,
   -0.0561112,
   -0.0841681,
   -0.362533,
   -0.520353,
   -0.757083,
   -0.845857,
   -0.979017,
   -1.28952,
   -1.47678,
   -1.75766,
   -2.17899,
   -2.81098,
   -3.04798,
   -3.13686,
   -3.27017,
   -3.47013,
   -3.77008,
   -4.22001,
   -4.89489,
   -5.14797,
   -5.5276,
   -5.66996,
   -5.88349,
   -6.2038,
   -6.68426,
   -6.86443,
   -7.13469,
   -7.54008,
   -8.14817,
   -9.06029,
   -9.40234,
   -9.91541,
   -10.539,
   -11.9821,
   -14.1466,
   -14.9583,
   -15.2627,
   -15.7192,
   -16.4041,
   -17.4314,
   -18.9723,
   -19.5502,
   -20.417,
   -21.2838,
   -22.1505,
   -23.4507,
   -25.401,
   -28.3264,
   -29.0577,
   -29.7891,
   -30.8861,
   -32.5316,
   -34.1772,
   -35.8227,
   -36.2341,
   -36.6455,
   -37.2626,
   -38.1882,
   -39.5766,
   -41.6592,
   -44.7832,
   -45.5641,
   -46.3451,
   -47.5166,
   -49.2738,
   -49.9328,
   -50.1799,
   -50.5506,
   -50.6896,
   -50.8981,
   -51.2108,
   -51.6799,
   -52.3836,
   -53.4391,
   -55.0224,
   -57.3973,
   -58.2879,
   -59.6238,
   -61.6277,
   -62.3791,
   -63.5063,
   -65.197,
   -65.831,
   -66.7821,
   -68.2086,
   -70.3485,
   -71.1509,
   -72.3546,
   -74.1601,
   -74.8371,
   -75.091,
   -75.4719,
   -75.6147,
   -75.6683,
   -75.6883,
   -75.7185,
   -75.7298,
   -75.7467,
   -75.7531,
   -75.7555,
   -75.7563,
   -75.7577,
   -75.7597,
   -75.7602
])
TBSAs = np.array([
    0,
    670.059,
    706.482,
    706.779,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    706.78,
    783.23,
    1174.93,
    1762.49,
    2643.82,
    3965.82,
    5948.83,
    25577.4,
    35987.3,
    50291.4,
    55409.3,
    62948.2,
    79357.5,
    87688.3,
    98330.9,
    111320,
    126914,
    132450,
    134505,
    137566,
    142320,
    149541,
    160422,
    176364,
    182260,
    191045,
    194331,
    199262,
    206589,
    217679,
    221818,
    228008,
    237239,
    251015,
    271468,
    279066,
    290360,
    304210,
    335447,
    380266,
    396316,
    402205,
    410893,
    423583,
    441804,
    466546,
    474782,
    485995,
    496153,
    505255,
    516344,
    528389,
    536346,
    536405,
    534793,
    529182,
    522529,
    517273,
    512096,
    510689,
    509351,
    507355,
    505025,
    504204,
    505666,
    508036,
    508639,
    509508,
    510888,
    512963,
    513720,
    514001,
    514419,
    514574,
    514802,
    515136,
    515617,
    516307,
    517270,
    518506,
    519981,
    520477,
    521134,
    521944,
    522229,
    522627,
    523166,
    523362,
    523638,
    524023,
    524557,
    524753,
    525036,
    525434,
    525579,
    525634,
    525714,
    525745,
    525759,
    525771,
    525774,
    525778,
    525784,
    525787,
    525788,
    525788,
    525788,
    525788,
    525789
]) * 0.001


import matplotlib.pyplot as plt

# --- Assumed variables (example definitions) ---
# T1, M1, T2, M2: Moment-rotation data for nonlinear and linear analyses.
# D1, TBSN, D2, TBSL: Displacement and base shear data.
# D4, TBSS: Experimental data (if available; otherwise empty lists or arrays).
# D5, TBSA: SAP2000 analysis data.
# D6, TBSAs: ABAQUS solid analysis data.
# F1, F2: Force data from nonlinear and linear analyses.
# X1, Y1, X2, Y2: Coordinates for deflected shapes.
# INT_N_f2, INT_N_f3, INT_N_f5, INT_N_f6, INT_N_f7, INT_N_f9, INT_N_f10:
#     Internal force results from the nonlinear analysis.
# INT_L_f3, INT_L_f6, INT_L_f7, INT_L_f10:
#     Internal force results from the linear analysis.
# a1, a2, L: Geometry parameters.
# starttime: The starting CPU time (e.g., from time.process_time()).

# For example, if these are defined elsewhere in your program, do not re-define them here.
# Here we assume they exist.

# ---------------------------
# Plotting Figures
# ---------------------------

# Figure 1: Moment-Rotation Diagram
plt.figure(1)
plt.plot(T1, M1, T2, M2, 'r--', linewidth=3)
plt.grid(True)
# The MATLAB legend uses 'NorthEastOutside' placement. In matplotlib you can use bbox_to_anchor.
plt.legend(['Nonlinear Analysis', 'Linear Analysis'], loc='upper right', bbox_to_anchor=(1.3, 1))
plt.xlabel('Rotation (rad) [DOF(6)]')
plt.ylabel('Moment (kN.mm) [DOF(6)]')
plt.title('Moment-Rotation diagram', color='b')

# Figure 2: Base shear-Displacement Diagram
plt.figure(2)
# Plot each data set. (If D4 and TBSS are empty, they will be skipped.)
plt.plot(D1, TBSN, linewidth=3, label='Nonlinear Analysis')
plt.plot(D2, TBSL, 'r--', linewidth=3, label='Linear Analysis')
plt.plot(D4, TBSS, 'y--', linewidth=3, label='Experimental Analysis')
plt.plot(D5, TBSA, 'g--', linewidth=3, label='SAP2000 Analysis')
plt.plot(D6, TBSAs, 'k--', linewidth=3, label='ABAQUS Analysis- solid')
plt.grid(True)
plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
plt.xlabel('Displacement (mm) [DOF(5)]')
plt.ylabel('Base shear (kN)')
plt.title('Base shear-Displacement diagram', color='b')

# Figure 3: Force-Displacement Diagram
plt.figure(3)
plt.plot(D1, F1, D2, F2, 'r--', linewidth=3)
plt.grid(True)
plt.legend(['Nonlinear Analysis', 'Linear Analysis'], loc='upper right', bbox_to_anchor=(1.3, 1))
plt.xlabel('Displacement (mm) [DOF(5)]')
plt.ylabel('Force (kN) [DOF(5)]')
plt.title('Force-Displacement diagram', color='b')

# Figure 4: Deflection along the length of the beam
plt.figure(4)
plt.plot(X1, Y1, X2, Y2, 'g--', linewidth=3)
plt.grid(True)
plt.legend(['Nonlinear', 'Linear'], loc='upper right', bbox_to_anchor=(1.3, 1))
plt.xlabel('x (mm)')
plt.ylabel('Deflection (mm)')
plt.title('Deflection along the length of beam', color='b')



# ---------------------------
# Print Computation Time
# ---------------------------
totaltime = time.process_time() - starttime
print(f'\nTotal time (s)= {totaltime:7.4f} \n')


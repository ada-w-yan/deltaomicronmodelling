# the TEIV model with cell types

# initial conditions
initial(T1) <- T_0
initial(L) <- L_0
initial(I) <- I_0
initial(V) <- V_0
initial(V_tot) <- V_tot0

# equations
deriv(T1) <- - max(0, beta * T1 * V)
deriv(L) <- max(0, beta * T1 * V) - max(0, tau * L)
deriv(I) <- max(0, tau * L) - max(0, delta * I)
deriv(V) <- max(0, omega_Inf * I) - max(0, kappa_Inf * V) - max(0, beta * T1 * V)
deriv(V_tot) <- max(0, omega_RNA * I) - max(0, kappa_RNA * V_tot)

# parameter values
beta <- user()
tau <- user()
delta <- user()
omega_Inf <- user()
omega_RNA <- user()
kappa_Inf <- user()
kappa_RNA <- user()
T_0 <- user()
V_0 <- user()
L_0 <- user()
I_0 <- user()
V_tot0 <- user()

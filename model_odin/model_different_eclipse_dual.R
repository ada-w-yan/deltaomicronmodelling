# the TEIV model with cell types

# initial conditions
initial(T1[1:2]) <- T_0[i]
# L[1] is tmprss2+ cells infected by the endosomal pathway
# L[2] is tmprss2+ cells infected by the tmprss2 pathway
# L[3] is tmprss2- cells infected by the endosomal pathway
initial(L[1:3]) <- L_0[i]
initial(I[1:2]) <- I_0[i]
initial(V) <- V_0
initial(V_tot) <- V_tot0

# equations
deriv(T1[1]) <- - max(0, (beta_E + beta_T) * T1[i] * V)
deriv(T1[2]) <- - max(0, beta_E * T1[i] * V)
deriv(L[1]) <- max(0, beta_E * T1[1] * V) - max(0, tau_E * L[i])
deriv(L[2]) <- max(0, beta_T * T1[1] * V) - max(0, tau_T * L[i])
deriv(L[3]) <- max(0, beta_E * T1[2] * V) - max(0, tau_E * L[i])
deriv(I[1]) <- max(0, tau_E * L[1] + tau_T * L[2]) - max(0, delta * I[i])
deriv(I[2]) <- max(0, tau_E * L[3]) - max(0, delta * I[i])
deriv(V) <- max(0, omega_Inf * sum(I)) - max(0, kappa_Inf * V) - max(0,  ((beta_E + beta_T) * T1[1] + beta_E * T1[2]) * V)
deriv(V_tot) <- max(0, omega_RNA * sum(I)) - max(0, kappa_RNA * V_tot) - max(0,  ((beta_E + beta_T) * T1[1] + beta_E * T1[2]) * V)

dim(T1) <- 2
dim(L) <- 3
dim(I) <- 2
dim(L_0) <- 3
dim(I_0) <- 2
dim(T_0) <- 2

# parameter values
beta_E <- user()
beta_T <- user()
tau_E <- user()
tau_T <- user()
delta <- user()
omega_Inf <- user()
omega_RNA <- user()
kappa_Inf <- user()
kappa_RNA <- user()
T_0[] <- user()
V_0 <- user()
V_tot0 <- user()
L_0[] <- user()
I_0[] <- user()

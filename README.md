# transient-1D-Heat-Diffusion-in-fin
solving transient-1D-Heat diffusion in fin by FD method includes FTCS, BTCS, Crank-Nicolson, CTCS discretizations
In this Project we solve 1D-Transient heat diffusion equation in a pin fin with uniform cross section by finite difference method
including FTCS, BTCS, Crank-Nicolson, CTCS Discretizations

The codes includes below parts :
1-Definition of constant parameters
2-Assembling coefficient matrix (A) in a sparse form (for implicit methods)
3-Assembling right hand side matrix (B)
4-solving the linear algebraic sparse system of equations by A\b method
Hint : parts 2,3,4 are a little different in Implicit and Explicit methods
5-Error Analysis (calculation of first and second norm of error)
6-validation (comparision with Analytical solution)
7-PostProcessing ( Heat flux - fin efficiency - fin performance coefficient )


constant parameters :

R = fin radius = 1cm
L = fin length = 20cm
rho = fin density = 2700 kg/m^3
C = specific heat capacity of fin = 897 j/kg.k
k = fin conductivity = 200 w/(m.k)
h= ambient convective heat transfer coefficient = 25 w/(m^2.k)
T_inf = ambient temperature = 300k
T_base = temperature of fin base = 500k

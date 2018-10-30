
!This function  returns the value of the reference density
!in the Adam-Williamson Equation of State
function rho_eos(z,rho0,chi_T,gy,Ly)
implicit none
real(8) z, rho0, chi_T, gy, rho_eos,Ly
rho_eos = rho0*exp(-chi_T*gy*rho0*(Ly-z))
end function

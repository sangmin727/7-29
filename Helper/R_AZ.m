function output = R_AZ(s_inf,s_dr,l_cr,l_pz,W,E,lambda)

hh=(l_pz/1000);

R_AZ=(V_AZ_GQ(s_inf,s_dr,l_cr,l_pz+hh,W,E,lambda)-V_AZ_GQ(s_inf,s_dr,l_cr,l_pz-hh,W,E,lambda))./(2*hh);

output=R_AZ;


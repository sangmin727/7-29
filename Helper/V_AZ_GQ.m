%%% Volume of active zone by using Gaussian quadrature


function output=V_AZ_GQ(s_inf,s_dr,l_cr,l_pz,W,E,lambda)


% c1=0.1713245;
% c2=0.3607616;
% c3=0.4679139;
% c4=c3;
% c5=c2;
% c6=c1;
% t1=-0.932469514;
% t2=-0.661209386;
% t3=-0.238619186;
% t4=-t3;
% t5=-t2;
% t6=-t1;


%  
% c1=0.2369269;
% c2=0.4786287;
% c3=0.5688889;
% c4=c2;
% c5=c1;
% t1=-0.906179846;
% t2=-0.538469310;
% t3=0;
% t4=-t2;
% t5=-t1;

c1=0.3478548;
c2=0.6521452;
c3=c2;
c4=c1;
t1=-0.861136312;
t2=-0.339981044;
t3=-t2;
t4=-t1;


COD1=COD_tot_x1(s_inf,s_dr,l_cr,l_pz,W,E,(((l_pz-l_cr).*t1./2)+((l_pz+l_cr)./2)));
COD2=COD_tot_x1(s_inf,s_dr,l_cr,l_pz,W,E,(((l_pz-l_cr).*t2./2)+((l_pz+l_cr)./2)));
COD3=COD_tot_x1(s_inf,s_dr,l_cr,l_pz,W,E,(((l_pz-l_cr).*t3./2)+((l_pz+l_cr)./2)));
COD4=COD_tot_x1(s_inf,s_dr,l_cr,l_pz,W,E,(((l_pz-l_cr).*t4./2)+((l_pz+l_cr)./2)));
% COD5=COD_tot_x1(s_inf,s_dr,l_cr,l_pz,W,E,(((l_pz-l_cr).*t5./2)+((l_pz+l_cr)./2)));
% COD6=COD_tot_x1(s_inf,s_dr,l_cr,l_pz,W,E,(((l_pz-l_cr).*t6./2)+((l_pz+l_cr)./2)));

%%% COD_tot_x1(s_inf,s_dr,l_cr,l_pz,W,E,x1)
integral_GQ=0.5.*(l_pz-l_cr).*(c1*COD1+c2*COD2+c3*COD3+c4*COD4);
V_AZ_GQ=integral_GQ./(lambda-1);

output=V_AZ_GQ;


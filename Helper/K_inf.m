function output=K_inf(s_inf,l_pz,W)
% Stress intensity factor of DENT specimen due to the remote stress (MPa*mm^0.5)
% s_inf =  Remote stress (MPa)
% l_pz = crack layer length from the DENT specimen surface (mm)
% W = A half of the total width of DENT specimen
% Valid in 0.5% any l_pz/W

alpha=l_pz./W;

F=(1.122-0.561*alpha-0.205*(alpha.^2)+0.471*(alpha.^3)-0.190*(alpha.^4))./((1-alpha).^0.5);

K_inf=s_inf.*((pi*l_pz).^0.5).*F;



output=K_inf;


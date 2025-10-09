function output=S_inf_LPD(LPD_tot,s_dr,l_cr,l_pz,W,H,E)

%%% Remote stress with load point displacement (LPD)
%%% S_inf = Remote stress (MPa)
%%% LPD_tot = Total Load point displacement (mm), 시간에 따라 선형적으로 증가함. (ramp
%%% test)
%%% th = Thickness of DENT specimen (mm)
%%% W = a half of total width of DENT specimen (mm)
%%% a = Crack length (mm)
%%% H = a half of total length of DENT specimen (mm)
%%% E = Elastic modulus of material (MPa)

% if l_pz < 0.995*W;
%     

% 
% else
%     S_inf_LPD=s_dr.*(W-l_cr)./(W);
% end

if l_pz > 0.99*W;
    l_pz=0.99*W;
end

    LPD_inf=LPD_tot+LPD_dr_revised(s_dr,l_cr,l_pz,W,H,E);
    V_2=(1./(pi.*l_pz./(2.*W))).*(0.0629-0.0610.*((cos(pi.*l_pz./(2.*W))).^4)-0.0019*((cos(pi.*l_pz./(2.*W))).^8)+log(sec(pi.*l_pz./(2.*W))));
    S_inf_LPD=(LPD_inf.*E)./((4.*l_pz.*V_2)+(2.*H));

output=S_inf_LPD;


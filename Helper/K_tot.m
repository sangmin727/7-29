function output=K_tot(s_inf,s_dr,l_cr,l_pz,W)

if l_cr==l_pz
    K_tot=K_inf(s_inf,l_pz,W);
else
    

K_tot=K_inf(s_inf,l_pz,W)-K_dr(s_dr,l_cr,l_pz,W);

    if K_tot < 0
        K_tot=0;
    end

end


output=K_tot;


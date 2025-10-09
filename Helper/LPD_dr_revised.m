function output = LPD_dr_revised(s_dr,l_cr,l_pz,W,s,E)

if l_cr==l_pz
    LPD_dr_revised=0;
else
    
n=200; 
dxi=((l_pz-l_cr)/n);
ds2=((l_pz-l_cr)/n);
K=1:1:n-1;
M=1:1:n-1;
[m,k]=meshgrid(K,M);
COD_dr_segment=triu((2*s_dr./E).*G_SIF((l_cr+(dxi.*(m+0.5))),W,(l_cr+(ds2.*(k-0.5)))).*G_SIF_LPD((l_cr+(dxi.*(m+0.5))),W,s).*dxi*ds2);

LPD_dr_revised=2*sum(COD_dr_segment(:));
%%% The additional factor of 2 occurs because two crack tips
%%% contribute to the displacement computation.
%%% 한 쌍의 virtual load P 가 두개의 균열 팁 성장에 영향을 주므로 2를 곱해야 함. 
    
end

output=LPD_dr_revised;


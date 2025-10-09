function output = K_dr(s_dr,l_cr,l_pz,W)

nx=10000;
hx=(l_pz-l_cr)./nx;
x=linspace(l_cr+0.5*hx,l_pz-0.5*hx,nx);
element=G_SIF(l_pz,W,x)*hx;
sum_elem=sum(element);
K_dr=s_dr*sum_elem;
output=K_dr;
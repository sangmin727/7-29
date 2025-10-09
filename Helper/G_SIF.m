function output=G_SIF(l_pz,W,x)
% G_SIF : Green's function of DENT (double edge notched tensile) specimen 
% x : distance from edge to load point (mm)
% l_pz : crack layer length from specimen edge (mm)
% W : A half of the total width of specimen (mm)
% Better than 1% for any a/b, c/a


alpha=l_pz./W;

f_s=(0.3).*(1-((x./l_pz).^(5./4)));
g=0.5.*(1-sin(0.5.*pi.*alpha)).*(2+sin(0.5.*pi.*alpha));
F3=((tan(0.5.*pi.*alpha)).^0.5)./((1-(((cos(0.5.*pi.*alpha))./(cos(0.5.*pi.*(x./W)))).^2)).^0.5);
F=(1+(f_s.*g)).*F3;
G_SIF=(2./((2.*W).^0.5)).*F;

output=G_SIF;


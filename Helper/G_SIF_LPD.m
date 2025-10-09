function output = G_SIF_LPD(a,W,s)
%G_SIF_LPD 이 함수의 요약 설명 위치
%SIF Green's function for depole force at (X,Y)=(0,s) and (0,-s)
%
nu=0.3; %%% Poisson's ratio of test material %%%%
alpha=(1+nu)./2; %%% at the plane stress 
a_b=a./W;
s_b=s./W;

F3=((tan(0.5.*pi.*a_b)).^0.5)./((1+((cos(0.5.*pi.*a_b)./sinh(0.5.*pi.*s_b)).^2)).^0.5);
F1=(1+0.122*((cos(0.5.*pi.*a_b)).^2)).*(1-(alpha.*((0.5*pi*s_b*coth(0.5*pi*s_b))./(1+(((sinh(0.5*pi*s_b))./(cos(0.5*pi*a_b))).^2))))).*F3;

G_SIF_LPD=F1./((2.*W).^0.5);
output=G_SIF_LPD;


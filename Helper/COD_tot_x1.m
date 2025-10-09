%%% COD_tot of DENT specimen at x=x1 (mm) %%%
%%% x1 is measured from the specimen surface %%%
%%% s_inf : Remote stress (MPa)
%%% s_dr  : Drawing stress of polymer (MPa)
%%% l_cr  : Crack length measured from surface (mm)
%%% l_pz  : Crack layer (CL) length measured from surface (mm)
%%% W     : A half of total width of specimen (mm)
%%% E     : Elastic modulus of polymer (MPa)


function output=COD_tot_x1(s_inf,s_dr,l_cr,l_pz,W,E,x1)

n=200;


dxi=((l_pz-x1)/n);
ds1=x1/n;
ds2=((l_pz-x1)/n);
  
if l_cr==l_pz || x1 >= l_pz
    COD_tot_x1=0;
else
    I=1:1:n;
    J=1:1:n;
    [i,j]=meshgrid(I,J);
    xi1=((x1+(dxi.*(i-0.5))));
    s1=(ds1.*(j-0.5));
    COD1_segment=(2*s_inf./E).*G_SIF(xi1,W,s1).*G_SIF(xi1,W,x1).*dxi*ds1;
    COD1=sum(COD1_segment(:));
    
    K=1:1:n-1;
    M=1:1:n-1;
    [m,k]=meshgrid(K,M);
    COD2_segment=triu((2*s_inf./E).*G_SIF((x1+(dxi.*(m+0.5))),W,(x1+(ds2.*(k-0.5)))).*G_SIF((x1+(dxi.*(m+0.5))),W,x1).*dxi*ds2);
    COD2=sum(COD2_segment(:));
    
        
    
    
    if x1==l_cr
        COD_dr_segment=triu((2*s_dr./E).*G_SIF((x1+(dxi.*(m+0.5))),W,(x1+(ds2.*(k-0.5)))).*G_SIF((x1+(dxi.*(m+0.5))),W,x1).*dxi*ds2);
        COD_dr=sum(COD_dr_segment(:));
    
    else
        ds_dr1=(x1-l_cr)/n;
        ds_dr2=(l_pz-x1)/n;
        dxi_dr=(l_pz-x1)/n; 
        H=1:1:n;
        G=1:1:n;
        [g,h]=meshgrid(H,G);
        COD_dr_1_segment=((2.*s_dr)./E).*G_SIF((x1+(dxi_dr.*(h-0.5))),W,x1).*G_SIF((x1+(dxi_dr.*(h-0.5))),W,(ds_dr1.*(g-0.5)))*dxi_dr*ds_dr1;
        COD_dr_1=sum(COD_dr_1_segment(:));
        
        Q=1:1:n-1;
        WW=1:1:n-1;
        [ww,q]=meshgrid(Q,WW);
        COD_dr_2_segment=triu(((2.*s_dr)./E).*G_SIF((x1+(dxi_dr.*(ww+0.5))),W,x1).*G_SIF((x1+(dxi_dr.*(ww+0.5))),W,(x1+(ds_dr2.*(q-0.5))))*dxi_dr*ds_dr2);
        COD_dr_2=sum(COD_dr_2_segment(:));
        COD_dr=COD_dr_1+COD_dr_2;
        
    end
    
    
    COD_tot_x1=COD1+COD2-COD_dr;
end

if COD_tot_x1<0
    COD_tot_x1=0;
end


output=COD_tot_x1;



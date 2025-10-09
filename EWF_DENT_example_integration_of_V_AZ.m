tic;
clc;
clear all;
clf;
addpath('Helper\')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_t=1000000;
dt=0.1;
% s_inf=6;
s_dr=16;
l=8;
W=15;
H=30;
th=2;
E=300;
LPD_rate=0.833333; %%% Displacement rate ( 50 mm/min)
gamma_0=10;
gamma_tr=20;
lambda=6;
t_star=10;
gamma_tilda=0;
k_cr_0=3;
k_pz_0=20;
k1c=3.0;
r=1;
mm=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l_cr(1)=l;
l_pz(1)=l;
gamma(1)=gamma_0;
t(1)=0;

for i=1:n_t
    
    
        LPD_tot(i)=LPD_rate.*t(i);
        
        if l_pz(i)==W;
            s_inf(i)=1.5.*s_dr.*(W-l_cr(i))./W;
        else   
            s_inf(i)=S_inf_LPD(LPD_tot(i),s_dr,l_cr(i),l_pz(i),W,H,E);
        end
        
        
        t(i+1)=t(i)+dt;
        
        k_cr=k_cr_0.*0.05*(W-l_cr(i))./(W-l);
        k_pz=k_pz_0.*0.5*(W-l_pz(i))./(W-l);
        
%         if t(i+1)-t(i) > 300;
%             t(i+1)=t(i)+300;
%         end
        


    if l_pz(i) < W;
        if ((((K_tot(s_inf(i),s_dr,l_cr(i),l_pz(i),W)).^2)/E)-(gamma_tr*R_AZ(s_inf(i),s_dr,l_cr(i),l_pz(i),W,E,lambda))) > 0
            l_pz(i+1)=l_pz(i)+dt*k_pz*((((K_tot(s_inf(i),s_dr,l_cr(i),l_pz(i),W))^2)/E)-(gamma_tr*R_AZ(s_inf(i),s_dr,l_cr(i),l_pz(i),W,E,lambda)));
        else
            l_pz(i+1)=l_pz(i);
        end
    else
        l_pz(i+1)=W;
    end
         
    
    if l_pz(i+1) >= W;
        l_pz(i+1)=W;
    end
     
    if l_pz(i) < W && l_pz(i+1)==W;
        LPD_cont=LPD_tot(i);
    end
    

    
    if l_pz(i+1)==W;
      
       ligament_cont=(2.*W)-(2.*l_cr(i));
       ligament=((-0.3*(LPD_tot(i)-LPD_cont)./LPD_cont).*ligament_cont)+ligament_cont;
       l_cr(i+1)=W-0.5*ligament;
       
        
     else
        if J_cr(s_inf(i),l_cr(i),W,E)-2*gamma(i) > 0
            l_cr(i+1)=l_cr(i)+dt*k_cr*(J_cr(s_inf(i),l_cr(i),W,E)-2*gamma(i));
        else
            l_cr(i+1)=l_cr(i);
        end
     end

        

     if l_pz(i+1) >= W;
        l_pz(i+1)=W;
     end
            
        
        
       
        
        if l_cr(i+1)>l_pz(i+1)     
            l_cr(i+1)=l_pz(i+1);
                %%%% Crack length must be smaller than or equal to the CL length at the same time %%%%
        end       
        
        if l_cr(i+1)==l_pz(i+1)
                gamma(i+1)=gamma_0;
        end   
        
        for j=1:i
            if (l_pz(i-j+1)<=l_cr(i+1)) && (l_cr(i+1)<l_pz(i-j+2))
%                gamma(i+1)=gamma_0*((1/(1+(((t(i+1)-t(i+1-j))/t_star)^r)))-(1-exp(((k^2)/d)*(t(i+1)-t(i+1-j)))*erfc((k/sqrt(d))*sqrt((t(i+1)-t(i+1-j)))))^eta);
                   gamma(i+1)=gamma_0/(1+(((t(i+1)-t(i+1-j))/t_star)^r)); 
            end
        end        %%%% Calculating Surface energy 'gamma' at the each loop (Fraction form) %%%%
        
        J(1)=J_cr(s_inf(i),l_cr(1),W,E);
        J(i+1)=J_cr(s_inf(i),l_cr(i+1),W,E);
        
                
%         if  (i>=12) && (l_cr(i-10)==l_cr(i-1)) && (l_cr(i-4)==l_cr(i-1)) && (l_cr(i-2)==l_cr(i-1)) && (l_cr(i-2)==l_cr(i-1)) && (l_cr(i)>l_cr(i-1)) && (l_cr(i+1)>l_cr(i)) && ((l_cr(i+1)-l_cr(i))>(l_cr(i)-l_cr(i-1)))
%             step(mm)=i;
%             mm=mm+1;
%         end
        
                
%         if K_inf(s_inf(i),l_cr(i+1),W) >= (k1c*(10^1.5))
%             l_pz(i+1)=W;
%             l_cr(i+1)=W;
%             break;
            
        if  l_cr(i+1) >= W
            l_cr(i+1)=W;
            s_inf(i+1)=s_inf(i);
            LPD_tot(i+1)=LPD_tot(i);
            break;
        end   
end

t=t.';
l_cr=l_cr.';
l_pz=l_pz.';
plot(LPD_tot,s_inf,LPD_tot,l_cr,LPD_tot,l_pz); legend('s_inf','l_cr','l_pz')

% mm=mm-1;
% lengthstep=l_pz(step)-l_cr(step);
% stepcrack=l_cr(step);
% 
% dumm(2:mm)=diff(t(step));
% timestep=t(step);
% timestep(2:mm)=dumm(2:mm);
% 
% dadt=lengthstep./timestep;
% logdadt=log10(dadt);


% for dum=1:mm
%     SIF(dum)=K_inf(s_inf,stepcrack(dum),W);
% end
% 
% logSIF=log10(SIF);
% logSIF=logSIF';
% 
% % logdadt(end)=[];
% % logSIF(end)=[];
% 
% fit=polyfit(logSIF,logdadt,1);
% slopefit=fit(1)*logSIF+fit(2);
% 
% subplot(1,2,2);
% scatter(logSIF,logdadt);
% hold on;
% plot(logSIF,slopefit,'r-.');
% xlabel('log(SIF), MPa \cdot m^{0.5}','fontsize',15)
% ylabel('log(da/dt), mm/s','fontsize',15)
% set(gca,'fontsize',15)
% text(logSIF(1)+0.05,logdadt(end)-0.1,  ['\boldmath$log(da/dt)=' num2str(fit(1)) '\; log(SIF)+' num2str(fit(2)) '\;$'],'fontsize',14,'Interpreter','latex')
% slope_m=fit(1);
% 


% %%%%%%%%%%%%%%%_Graph Plotting_%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(1,2,1)
% plot(t,l_cr,'-k',t,l_pz,'-r','linewidth',1)
% legend({'Crack', 'Crack Layer(WZ+AZ)'},'fontsize',15,'fontweight','bold','location','northwest')
% axis([0 t(end)+1000 0.5 5.5])
% set(1,'color','w')
% xlabel('Time (sec)','fontsize',15)
% ylabel('Length (mm)','fontsize',15)
% set(gca,'fontsize',15)
% 
% title('SENT specimen, Surface energy  $\gamma= \frac{{\gamma}_{o}}{1 + {\left(\frac{t_i-t_o}{t^*}\right)}^{r}}$','interpreter','latex','fontsize',28)
% 
% text(1000,4.5, ['\boldmath$s_inf=' num2str(s_inf) '\;MPa$'],'fontsize',14,'Interpreter','latex')
% text(1000,4.25, ['\boldmath$thickness=' num2str(th) '\;mm$'],'fontsize',14,'Interpreter','latex')
% text(1000,4, ['\boldmath$\sigma_{dr}=' num2str(s_dr) '\;MPa$'],'fontsize',14,'Interpreter','latex')
% text(1000,3.75, ['\boldmath$l_o=' num2str(l) '\;mm$'],'fontsize',14,'Interpreter','latex')
% text(1000,3.5, ['\boldmath$W=' num2str(W) '\;mm$'],'fontsize',14,'Interpreter','latex')
% text(1000,3.25, ['\boldmath$E=' num2str(E) '\;MPa$'],'fontsize',14,'Interpreter','latex')
% text(1000,3, ['\boldmath$\tilde{\gamma}= ' num2str(gamma_tilda) '$'] ,'fontsize',14,'Interpreter','latex')
% text(1000,2.75, ['\boldmath$\gamma^{tr}=' num2str(gamma_tr) '\;mJ/mm^{3}$'],'fontsize',14,'Interpreter','latex')
% text(1000,2.5, ['\boldmath$\lambda=' num2str(lambda) '$'],'fontsize',14,'Interpreter','latex')
% text(t(end)/3,4.25, ['\boldmath$\gamma_{o}=' num2str(gamma_0) '\;mJ/mm^{2}$'],'fontsize',14,'Interpreter','latex')
% text(t(end)/3,4, ['\boldmath$t^\ast=' num2str(t_star) '\;sec$'] ,'fontsize',14,'Interpreter','latex')
% text(t(end)/3,3.75, ['\boldmath$r=' num2str(r) '$'] ,'fontsize',14,'Interpreter','latex')
% text(t(end)/3,3.5, ['\boldmath$k_{cr}=' num2str(k_cr) '\;mm^{3}/mJ \cdot s$'] ,'fontsize',14,'Interpreter','latex')
% text(t(end)/3,3.25, ['\boldmath$k_{pz}=' num2str(k_pz) '\;mm^{3}/mJ \cdot s$'] ,'fontsize',14,'Interpreter','latex')
% text(t(end)/3,3, ['\boldmath$K_{1c}=' num2str(k1c),'\;MPa \cdot \sqrt{m}$'], 'fontsize',14,'Interpreter','latex')
% text(t(end)/3,2.75, ['\boldmath$dt=' num2str(dt),'\;sec$'], 'fontsize',14,'Interpreter','latex')
% text(t(end)/3,2.5, ['\boldmath$t_f=' num2str(t(end)/3600),'\;hours$'], 'fontsize',14,'Interpreter','latex')
% 
% h = msgbox('Operation Completed');

toc;
%% THERMAL EXHUMATION: GEOTHERM COMPUTATION IN THE CASE OF A SINUSOIDAL TOPOGRAPHY
clear all; close all;

% Define grid spacing
dz=10;
z=[-5000:dz:30000];

% Heat model variables and parameters
qm       = 50e-3;                     % Mantle heat flow [W/m^2]
rhoH0    = 1.5e-6;                    % Surface radiogenic heat production rate per unit mass times rock density [W/m^3]
hr       = 10000;                     % Scale depth of radiogenic heat production [m]
k        = 1.5.*ones(size(z));        % Rock conductivity [W/m/K]
k1        = 1.75.*ones(size(z));
k2        = 2.*ones(size(z));
cp       = 1000.*ones(size(z));       % Rock specific Heat [J/kg/K]
rho      = 2800.*ones(size(z));       % Rock density [kg/m3]
k(z<0)   = 1.69;                       % Basalt conductivity [W/m/K]
cp(z<0)  = 1000;                      % Basalt specific Heat [J/kg/K]
rho(z<0) = 2800;                      % Basalt density [kg/m3]
K        = k./(rho.*cp);              % Rock diffusivity
K1        = k1./(rho.*cp);
K2        = k2./(rho.*cp);
T0       = 0;                         % Mean surface temperature [C]
T1       = 500;                       % Bottom surface temperature [C] [NOT USED HERE)
C1       = -qm./k;                    % Bottom gradient heat flow [K/m] [USED HERE]
C11       = -qm./k1;
C12       = -qm./k2;
Tbasalts = 1200;                      % Initial temperature of basalats
Zbaslats = 200;                       % Basalts thickness (for each episode) [m]
tbasalts = 1000*(365*24*3600);      % Time interval between each basaltic stage [s]
nstage   = 2;                        % Number of basaltic stages
compteur=0;

for b=1:1:3;
    
% Define heat production
Q=zeros(size(z));Q=rhoH0.*exp(-z./hr);Q(z<=0)=0;

% Define initial temperature 
T=zeros(size(z));Tnew=zeros(size(z));
% Initial temperature (C)
% Tini = zeros(size(T));
if b==1;
Tini = T0 + qm.*z./k + rhoH0.*hr^2./k.*(1-exp(-z./hr));Tini(z<=0)=0;% Following Tuctotte & Schubert
T=Tini;Told=Tini;Tnew=Tini;
end
if b==2;
Tini = T0 + qm.*z./k + rhoH0.*hr^2./k1.*(1-exp(-z./hr));Tini(z<=0)=0;% Following Tuctotte & Schubert
T=Tini;Told=Tini;Tnew=Tini;
end
if b==3;
Tini = T0 + qm.*z./k + rhoH0.*hr^2./k2.*(1-exp(-z./hr));Tini(z<=0)=0;% Following Tuctotte & Schubert
T=Tini;Told=Tini;Tnew=Tini;
end

% Define CFL and "wave" speed
if b==1;
dt=0.1*dz^2/(2*max(K));sz=K*dt/dz^2;
end
if b==2;
dt=0.1*dz^2/(2*max(K1));sz=K1*dt/dz^2;
end
if b==3;
dt=0.1*dz^2/(2*max(K2));sz=K2*dt/dz^2;
end

t=0;it=0;conv=0;iconv=0;stage=0; tstartstage=0;

figure('units','normalized','outerposition',[0 0 1 1])

while iconv==0 || stage<=nstage
    
    % Update time and iteration number
    it=it+1;
    t=t+dt;       
    
    % Compute new temperature
    Tnew=T;
    Diff=sz(2:end-1).*(T(3:end)-2*T(2:end-1)+T(1:end-2)); % Heat Diffusion
    Prod=dt.*Q(2:end-1)./(rho(2:end-1).*cp(2:end-1)); % Heat production
    Tnew(2:end-1)=T(2:end-1)+Diff+Prod; % Tnew=Told+Diff+Prod;
    
    % --- Boundary conditions
    % Surface - set to T0 and basal
    if t-tstartstage>=tbasalts && stage>0 && stage<=nstage || iconv==1 && stage==0
        % reset convergence criteria
        iconv=0;
        % This the first phase of this stage
        stage=stage+1;tstartstage=t;
        % set temperature at surface to T0
        Tnew(z<=-stage.*Zbaslats)=T0;
        % and set basalts temperature to Tbasalts
        Tnew(z<-(stage-1).*Zbaslats & z>-stage.*Zbaslats)=Tbasalts;  
    else
        % In any case keep surface temperature at T0
        Tnew(z<=-stage.*Zbaslats)=T0;
    end
    % Bottom - set to T1 or we can impose a mantle heat flow C1 (Neumann type)
    %Tnew(end)=T1;
    if b==1;
    Tnew(end)=T(end)+sz(end).*(2*T(end-1)-2*(T(end)+dz.*C1(end)))+dt.*Q(end)./(rho(end).*cp(end));
    end
    if b==2;
    Tnew(end)=T(end)+sz(end).*(2*T(end-1)-2*(T(end)+dz.*C11(end)))+dt.*Q(end)./(rho(end).*cp(end));
    end
    if b==3;
    Tnew(end)=T(end)+sz(end).*(2*T(end-1)-2*(T(end)+dz.*C12(end)))+dt.*Q(end)./(rho(end).*cp(end));
    end
    % Update T
    T=Tnew;   
      
    % Plot 
    if it==1  || stage==0 && mod(it,10000)==0 || stage>=1 && stage<nstage && mod(it,100)==0 || stage>=nstage && mod(it,10000)==0 
        subplot(2,2,1);     plot(T,z./1000,'k-',Tini,z./1000,'k--',T-Tini,z./1000,'r-');        
                            hold on;plot(T(z==100),z(z==100)./1000,'bo',T(z==500),z(z==500)./1000,'co',T(z==1000),z(z==1000)./1000,'go',T(z==2500),z(z==2500)./1000,'yo',T(z==5000),z(z==5000)./1000,'ro');      
                            hold off;title(['t = ' num2str(round(t/(365*24*3600*1000000)*1000)/1000) ' Myr, stage=' num2str(stage) ]);xlabel('T (^oC)');ylabel('z (km)');axis ij;legend('T','T_{ini}','T-T_{ini}','0.1 km','0.5 km', '1 km','2.5 km', '5 km');
        subplot(2,2,2);     plot(T-Told,z./1000,'k',Diff,z(2:end-1)./1000,'r',Prod,z(2:end-1)./1000,'b');
                            title(['conv = ' num2str(conv) '^oC' ]);legend('\Delta T','Diff','Prod');xlabel('\DeltaT (^oC)');ylabel('z (km)');axis ij;
%         subplot(2,2,[3 4]); plot(t./(365*24*3600),T(z==100)-Tini(z==100),'b.',t./(365*24*3600),T(z==500)-Tini(z==500),'c.',t./(365*24*3600),T(z==1000)-Tini(z==1000),'g.',t./(365*24*3600),T(z==2500)-Tini(z==2500),'y.',t./(365*24*3600),T(z==5000)-Tini(z==5000),'r.');
%                             hold on;xlabel('t (Myr)');ylabel('T-T_{ini} (^oC)');axis ij;
        subplot(2,2,[3 4]); plot(t./(365*24*3600),T(z==100)-Tini(z==100),'b.',t./(365*24*3600),T(z==800)-Tini(z==800),'r.');
                            hold on;xlabel('t (Myr)');ylabel('T-T_{ini} (^oC)');axis ij;
        drawnow            
    end
    
    % Check convergence to start the basaltic stages
    conv=max(max(abs(T-Told)));
    if conv< 1e-5 && stage==0;
        iconv=1;
        Tini=T;
    elseif conv< 1e-5 && stage>0
        iconv=1;
    end
    
    % Save Told for the nest timestep
    Told=T;
    Tdiff=T-Tini;
end
end

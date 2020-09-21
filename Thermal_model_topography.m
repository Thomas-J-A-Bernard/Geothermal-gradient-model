%% THERMAL EXHUMATION: GEOTHERM COMPUTATION IN THE CASE OF A SINUSOIDAL TOPOGRAPHY (AUTHORS: Philippe Steer, Thomas Bernard)
clear all; close all;

% Heat model variables and parameters
T0     = 0+273.15;                  % Mean surface temperature [K]
qm     = 30e-3;                     % Mantle heat flow [W/m^2]
rhoH0  = 2.0e-6;                    % Surface radiogenic heat production rate per unit mass times rock density [W/m^3]
hr     = 10000;                     % Scale depth of radiogenic heat production [m]
k      = 2.5;                       % Rock conductivity [W/m/K]
beta   = 0;                         % Atmosphere vertical temperature gradient [K/m-1] ref=6.5e-3
% Topography model variables and parameters
x=[0:0.1:100].*1000;                % Horizontal distance [m]
y=[-3:0.01:20].*1000;               % Depth [m]
dh      = 2000 ;                    % Vertical mean erosion [m]
lambda1 = 5000 ;                    % Initial valley wavelength [m]
lambda2 = 50000;                    % Final valley wavelength [m]
h01     = 750  ;                    % Initial valley relief [m]
h02     = 1250 ;                    % Final valley relief [m]
h1=h01.*cos(2.*pi.*x./lambda1);     % Initial surface elevation [m]
h2=h02.*cos(2.*pi.*x./lambda2);     % Final surface elevation [m]

% Computation 2D fields
T1=zeros(numel(y),numel(x));T2=T1;     % 2D temperature field [K]
for i=1:numel(x)
    T1(:,i) = T0 + (qm.*y(:)./k) + (rhoH0.*hr.^2./k).*(1 - exp(-y(:)./hr)) + (beta - qm./k - rhoH0.*hr./k).*h1(i).*exp(-2.*pi.*y(:)./lambda1);       % Initial model
    y2=y(:)-dh;hr2=hr-dh;
    T2(:,i) = T0 + (qm.*y2(:)./k) + (rhoH0.*hr2.^2./k).*(1 - exp(-y2(:)./hr2)) + (beta - qm./k - rhoH0.*hr2./k).*h2(i).*exp(-2.*pi.*y2(:)./lambda2); % Final model - different geometry  (dh and h2)
    T3(:,i) = T0 + (qm.*y2(:)./k) + (rhoH0.*hr2.^2./k).*(1 - exp(-y2(:)./hr2)) + (beta - qm./k - rhoH0.*hr2./k).*h1(i).*exp(-2.*pi.*y2(:)./lambda1);    % Final model - same geometry    (only dh)
    DT12(:,i) = T1(:,i)-T2(:,i);
    DT13(:,i) = T1(:,i)-T3(:,i);
    ind=find(y<h1(i));      T1(ind,i)=NaN;    DT12(ind,i)=NaN;  DT13(ind,i)=NaN;
    ind=find(y<(h2(i)+dh)); T2(ind,i)=NaN;    DT12(ind,i)=NaN;
    ind=find(y<(h1(i)+dh)); T3(ind,i)=NaN;                      DT13(ind,i)=NaN;
end
% Computations thermal exhumation in surface (1D)
t12=zeros(size(x));t22=t12;t13=t12;dt33=t12;dt12=t12;dt13=t12;
for i=1:numel(x)
    temp=(h2(i)+dh)-y;temp(temp>0)=1e32;[~,ind]=min(abs(temp));
    t12(i)=T1(ind,i);   t22(i)=T2(ind,i);    dt12(i)=t12(i)-t22(i);
    temp=(h1(i)+dh)-y;temp(temp>0)=1e32; [~,ind]=min(abs(temp));
    t13(i)=T1(ind,i);   t33(i)=T3(ind,i);    dt13(i)=t13(i)-t33(i);
end

% Plot
irelief=3;cmap=colormap(jet(2048));cmap(1,:)=1;
subplot(2,1,1); imagesc(x./1000,y./1000,T1-273.15); colormap(cmap);colorbar; hold on; plot(x./1000,h1./1000,'-k');                             axis equal tight; caxis([T0-273.15 max(max(T1-273.15))]); xlabel('x(km)');ylabel('y (km)'); title('Initial temperature');
subplot(2,1,2); imagesc(x./1000,y./1000,T1-273.15); colormap(cmap);colorbar; hold on; plot(x./1000,h1./1000,'-k');                             axis equal tight; caxis([T0-273.15 max(max(T1-273.15))]); xlabel('x(km)');ylabel('y (km)'); title('Initial temperature');
% % print('-dpdf','-painters',['Figure/Exhumation_model_ini_dh=' num2str(dh)  '_h1=' num2str(h01) '_h2=' num2str(h02) '_l1=' num2str(lambda1) '_l2=' num2str(lambda2)  ]);close;
subplot(2,1,1); imagesc(x./1000,y./1000,T3-273.15); colormap(cmap);colorbar; hold on; plot(x./1000,(h1+dh)./1000,'--k');                       axis equal tight; caxis([T0-273.15 max(max(T1-273.15))]); xlabel('x(km)');ylabel('y (km)'); title('Final temperature');
subplot(2,1,2); imagesc(x./1000,y./1000,T2-273.15); colormap(cmap);colorbar; hold on; plot(x./1000,(h2+dh)./1000,'--k');                       axis equal tight; caxis([T0-273.15 max(max(T1-273.15))]); xlabel('x(km)');ylabel('y (km)'); title('Final temperature');
% % print('-dpdf','-painters',['Figure/Exhumation_model_final_dh=' num2str(dh)  '_h1=' num2str(h01) '_h2=' num2str(h02) '_l1=' num2str(lambda1) '_l2=' num2str(lambda2)  ]);close;
subplot(2,1,1); imagesc(x./1000,y./1000,DT13);    colormap(cmap);colorbar;   hold on; plot(x./1000,h1./1000,'-k',x./1000,(h1+dh)./1000,'--k'); axis equal tight; caxis([ min(min(min(DT13)),min(min(DT12))) max(max(max(DT13)),max(max(DT12)))]); xlabel('x(km)');ylabel('y (km)'); title('Cooling');
subplot(2,1,2); imagesc(x./1000,y./1000,DT12);    colormap(cmap);colorbar;   hold on; plot(x./1000,h1./1000,'-k',x./1000,(h2+dh)./1000,'--k'); axis equal tight; caxis([ min(min(min(DT13)),min(min(DT12))) max(max(max(DT13)),max(max(DT12)))]); xlabel('x(km)');ylabel('y (km)'); title('Cooling');
% % print('-dpdf','-painters',['Figure/Exhumation_model_cooling_dh=' num2str(dh)  '_h1=' num2str(h01) '_h2=' num2str(h02) '_l1=' num2str(lambda1) '_l2=' num2str(lambda2)  ]);close;
subplot(2,1,1); plot(x./1000,dt13);ylim([ min(min(min(DT13)),min(min(DT12))) max(max(max(DT13)),max(max(DT12)))]); xlabel('x(km)');ylabel('Cooling (^o C)');
subplot(2,1,2); plot(x./1000,dt12);ylim([ min(min(min(DT13)),min(min(DT12))) max(max(max(DT13)),max(max(DT12)))]); xlabel('x(km)');ylabel('Cooling (^o C)');
% % print('-dpdf','-painters',['Figure/Exhumation_model_coolinglines_dh=' num2str(dh)  '_h1=' num2str(h01) '_h2=' num2str(h02) '_l1=' num2str(lambda1) '_l2=' num2str(lambda2)  ]);close;

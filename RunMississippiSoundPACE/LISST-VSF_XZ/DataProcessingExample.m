%% ========= NIST Tracable 25 um Polystyrene Particles ====================

% Select a file containing scattering data
scat_location = '.\ExampleData\NIST_25um_SN1664.VSF';

% Select a background file
background_location = '.\ExampleData\NIST_25um_background_SN1664.VSF';

% process the data using Make P function with plotting disabled
proc = lisstvsf_makep(1664, scat_location, background_location, [],1);

% take mean of P11 (VSF) and remove any nans 
VSF_mean = mean(proc.P11(~proc.qc_saturated,:));
angles = proc.Angles(~isnan(VSF_mean));
VSF_mean(isnan(VSF_mean)) = [];

% calculate P11 using Mie theory for 25um particles
diam = 25;
lambda = 0.516;

% choose angle range
ang_mie = [min(proc.Angles):0.01:15 15.1:0.1:max(proc.Angles)]'; 

% refractive index of polystyrene beads
m_p = 1.5663 + 0.00785./(lambda.^2) + 0.000334./(lambda.^4);

% refractive index of water
m_w = 1.3481216.*exp(-lambda.*1000.*0.0000193465) + 2.627716.*exp(-lambda.*1000.*0.0153948); 

% calculte Mie results
[S1,S2] = fastmie(pi*diam/lambda*m_w,m_p/m_w, ang_mie/180*pi); 

% calculate P11 
p11_mie = 0.5*abs(S1).^2 + 0.5*abs(S2).^2;

% plot theory vs LISST-VSF data normalized to area
figure
semilogy(ang_mie,p11_mie./trapz(ang_mie,p11_mie),'r','DisplayName','Mie Theory for 25um particles')
hold on
plot(angles,VSF_mean./trapz(angles,VSF_mean),'b','DisplayName','LISST-VSF data for 25um particles')
xlabel('Angle (Degrees)');
ylabel('VSF (relative scaling)');
legend('Show')

%% =========== NIST Tracable 496 nm Polystyrene Particles =================

scat_location = '.\ExampleData\NIST_496nm_SN1664.VSF';
background_location = '.\ExampleData\NIST_496nm_background_SN1664.VSF';

proc = lisstvsf_makep(1664, scat_location, background_location, [],1);
 
VSF_mean = mean(proc.P11(~proc.qc_saturated,:));
angles = proc.Angles(~isnan(VSF_mean));
VSF_mean(isnan(VSF_mean)) = [];
p12 = mean(proc.p12(~proc.qc_saturated,:));

diam = 0.496;
lambda = 0.516;
ang_mie = [min(proc.Angles):0.01:15 15.1:0.1:max(proc.Angles)]'; 
m_p = 1.5663 + 0.00785./(lambda.^2) + 0.000334./(lambda.^4);
m_w = 1.3481216.*exp(-lambda.*1000.*0.0000193465) + 2.627716.*exp(-lambda.*1000.*0.0153948); 
[S1,S2] = fastmie(pi*diam/lambda*m_w,m_p/m_w, ang_mie/180*pi); 
p11_mie = 0.5*abs(S1).^2 + 0.5*abs(S2).^2;

figure
semilogy(ang_mie,p11_mie./trapz(ang_mie,p11_mie),'r','DisplayName','Mie Theory for 496nm particles')
hold on
plot(angles,VSF_mean./trapz(angles,VSF_mean),'b','DisplayName','LISST-VSF data for 496nm particles')
xlabel('Angle (Degrees)'); 
ylabel('VSF (relative scaling)');
legend('Show')

function proc = lisstvsf_makep(SN, scatPath, zscPath, beamc, alpha, plotData) 
% Revised programme to process LISST-VSF data
% Modifications 11/21/2018
%   1. adding uncalibrate bulk P11 and P12 (i.e., withouth subtracting the
%       background)
%   2. calculate alpha using using the bulk four polarized components
%   (rr,rp,pp,pr) instead of using the pariculate rr,rp,pp and pr.
%   3. increase the outputs of uncalibrate P11,P12 and P22.

%  
% Modifications: 03/06/2019
%   1.calculate background of ring detectors using the median value
%   measured in all PMT gains. The default use the median value measured
%   only in matched PMT gains.(rings do not have gains)
%   2. the attenuation term of purfiled water is added 
%   3. the sample scattered light internsity is normalized to the reference
%   power.
%   4. the threshold for valid ring detector singal is changed from 25 to 1
%   to increase valid measurements.
%   
% updates 08/15/2019:
%   1.Negative values in four scattered polarized components (rp,pp,rr,pr) are set nan.

% updates 12/02/2019
%   1. Adding beamc as one input if the independent beamc measurement is available
%   2. If the LISST-VSF measured beamc <0, set it to the value of pure water.

proc.procver = '20171026';
proc.procdate = datestr(now);

% read in factory zscat, calibration factors, zscat, and data
calfact = lisstvsf_getcal(SN);
zsc = lisstvsf_readdat(zscPath);
scat = lisstvsf_readdat(scatPath);

% correct and save auxiliary data
proc.timestamp = scat.time;
proc.tempC = scat.tempC .* calfact.temp_slope + calfact.temp_offset;
proc.depth = scat.depth .* calfact.depth_slope + calfact.depth_offset;
proc.batt_volts = scat.batt_volts .* calfact.bat_slope + calfact.bat_offset;
proc.pmt_gain = scat.pmt_gain;

[ns,nang] = size(scat.rp);

% Convert encoder index in datastream from ADC board to actual angle, offset is contained in Cal_Factors file
angles = scat.angles_idx + calfact.angle_offset; % angles in degrees

proc.qc_saturated = any(scat.qc_saturated(:,angles > 14 & angles < 156),2);

%% ========================================================================
% Correct ring and eyeball signals for laser ref and attenuation along path
% =========================================================================

% Correct raw measurements for the reduction in laser power caused by the
% 1/2 wave plate
zsc.LREF(:,2) = zsc.LREF(:,2) .* calfact.HWPlate_transmission;
scat.LREF(:,2) = scat.LREF(:,2) .* calfact.HWPlate_transmission;
scat.rp = scat.rp .* calfact.HWPlate_transmission;
scat.rr = scat.rr .* calfact.HWPlate_transmission;
zsc.rp = zsc.rp .* calfact.HWPlate_transmission;
zsc.rr = zsc.rr .* calfact.HWPlate_transmission;
LREF_median = median(scat.LREF(:,2));
% Find number of PMT gain values in background file (usually 10)
zsc_pmt_values = unique(zsc.pmt_gain);
num_zsc_pmts = length(zsc_pmt_values); 

% Preallocate arrays
zsc_med.rp = nan(num_zsc_pmts,nang);
zsc_med.rr = nan(num_zsc_pmts,nang);
zsc_med.pp = nan(num_zsc_pmts,nang);
zsc_med.pr = nan(num_zsc_pmts,nang);
zsc_med.LP = nan(num_zsc_pmts,2);
zsc_med.LREF = nan(num_zsc_pmts,2);
zsc_med.rings1 = nan(num_zsc_pmts,32);
zsc_med.rings2 = nan(num_zsc_pmts,32);
zsc_med.pmt_gain = nan(num_zsc_pmts,1);

% Use median background at each PMT gain level for processing scat
for n = 1:num_zsc_pmts

pmt_idx = zsc.pmt_gain == zsc_pmt_values(n);
    
zsc_med.rp(n,:) = median(zsc.rp(pmt_idx,:));
zsc_med.rr(n,:) = median(zsc.rr(pmt_idx,:));
zsc_med.pp(n,:) = median(zsc.pp(pmt_idx,:));
zsc_med.pr(n,:) = median(zsc.pr(pmt_idx,:));
zsc_med.LP(n,:) = median(zsc.LP(pmt_idx,:));
zsc_med.LREF(n,:) = median(zsc.LREF(pmt_idx,:));
zsc_med.rings1(n,:) = median(zsc.rings1(pmt_idx,:));
zsc_med.rings2(n,:) = median(zsc.rings2(pmt_idx,:));  
zsc_med.pmt_gain(n,:) = median(zsc.pmt_gain(pmt_idx,:));

end

% Preallocate arrays
zsc_pmt_matched.rp = nan(ns,nang);
zsc_pmt_matched.rr = nan(ns,nang);
zsc_pmt_matched.pp = nan(ns,nang);
zsc_pmt_matched.pr = nan(ns,nang);
zsc_pmt_matched.LP = nan(ns,2);
zsc_pmt_matched.LREF = nan(ns,2);
zsc_pmt_matched.rings1 = nan(ns,32);
zsc_pmt_matched.rings2 = nan(ns,32);

% Create an array of backround measurements that match the PMT gain values in
% the scat data file
scat_pmt_values = unique(scat.pmt_gain);
for n = 1:length(scat_pmt_values)

pmt_idx = scat.pmt_gain == scat_pmt_values(n);
num_samples = sum(pmt_idx);
zsc_idx = zsc_med.pmt_gain == scat_pmt_values(n);
if ~any(zsc_idx)
    warning(['Scattering data collected at PMT gain value of ' num2str(scat_pmt_values(n)) ...
        ' has no corresponding background measurements at the same PMT gain in the provided background file.']);
    continue
end

zsc_pmt_matched.rp(pmt_idx,:) = repmat(zsc_med.rp(zsc_idx,:),num_samples,1);
zsc_pmt_matched.rr(pmt_idx,:) = repmat(zsc_med.rr(zsc_idx,:),num_samples,1);
zsc_pmt_matched.pp(pmt_idx,:) = repmat(zsc_med.pp(zsc_idx,:),num_samples,1);
zsc_pmt_matched.pr(pmt_idx,:) = repmat(zsc_med.pr(zsc_idx,:),num_samples,1);
zsc_pmt_matched.LP(pmt_idx,:) = repmat(zsc_med.LP(zsc_idx,:),num_samples,1);
zsc_pmt_matched.LREF(pmt_idx,:) = repmat(zsc_med.LREF(zsc_idx,:),num_samples,1);
% zsc_pmt_matched.rings1(pmt_idx,:) = repmat(zsc_med.rings1(zsc_idx,:),num_samples,1);
% zsc_pmt_matched.rings2(pmt_idx,:) = repmat(zsc_med.rings2(zsc_idx,:),num_samples,1);
    
end
zsc_pmt_matched.rings1 = repmat(median(zsc_med.rings1,1),ns,1);
zsc_pmt_matched.rings2 = repmat(median(zsc_med.rings2,1),ns,1);
if all(isnan(zsc_pmt_matched.rp(:,1)))
    error('The background file has no PMT gain values that match the gain values used in the scattering data file')
end

% distance along the beam to sample volume, then from the sample volume to eyeball [cm]
% used for attenuation correction later
%                 Eyeball (.)    
%                         / 
%                        /  
%  Receive Window |     ------------------| Transmit Window 
%                       ^
%                 Sample Volume
paths = 10 - 2*cot(angles*pi/180)+(2./sin(angles*pi/180));                       
  
% scale za,zb,zc,zd by zLREF and save for use in getting net scattering
scaled_zsc_rp = zsc_pmt_matched.rp .* repmat(scat.LREF(:,1)./zsc_pmt_matched.LREF(:,1),1,nang);                                                   
scaled_zsc_rr = zsc_pmt_matched.rr .* repmat(scat.LREF(:,1)./zsc_pmt_matched.LREF(:,1),1,nang);
scaled_zsc_pp = zsc_pmt_matched.pp .* repmat(scat.LREF(:,2)./zsc_pmt_matched.LREF(:,2),1,nang);
scaled_zsc_pr = zsc_pmt_matched.pr .* repmat(scat.LREF(:,2)./zsc_pmt_matched.LREF(:,2),1,nang);


% clean water ratio of transmitted laser power to reference; used to correct for laser drift
r = zsc_pmt_matched.LP ./ zsc_pmt_matched.LREF;  

% -----------------------------------------------------------------
% First rotation is laser polarized perpendicular, signals are then 
% rp and rr (a and c)
% -----------------------------------------------------------------
tau1 = scat.LP(:,1)./(r(:,1).*scat.LREF(:,1));  % laser reference drift compensated here.          
beam_c1 = -100/15*log(tau1);                 % 15 cm path length
idx_c1 = beam_c1<0;
if any(idx_c1)
    beam_c1(idx_c1)  = 0.038;
end
beam_cz = beam_c1*0+0.038;
if ~isempty(beamc)
    beam_c1 = repmat(beamc,length(tau1),1);    
end
att_factors_z = exp(-beam_cz*(paths/100));
% beam_c1 = 0.05;
att_factors1 = exp(-beam_c1*(paths/100));    % attenuation correction along beam + from SV to eyeball

rp = scat.rp./att_factors1;     
% Scale signal to account for attenuation along path and laser power
rp(rp<0) = nan;
rp_2 = rp;                                                   % without substract background scattering
rp = rp - scaled_zsc_rp./att_factors_z;                                    % Subtract background scattering
rp(rp<0)=nan;
rp = rp .* repmat(sin(angles*pi/180),ns,1)./repmat(LREF_median,ns,nang);                 % Correct for scattering volume lengthening with view angle
rp_2 = rp_2 .* repmat(sin(angles*pi/180),ns,1)./repmat(LREF_median,ns,nang);            

rr = scat.rr./att_factors1;     
rr(rr<0)=nan;
rr_2 = rr;
rr = rr - scaled_zsc_rr./att_factors_z;                       
rr(rr<0)=nan;
rr = rr .* repmat(sin(angles*pi/180),ns,1)./repmat(LREF_median,ns,nang);                                      
rr_2 = rr_2 .* repmat(sin(angles*pi/180),ns,1)./repmat(LREF_median,ns,nang);                                      
% -----------------------------------------------------------------
% Second rotation is with 1/2-wave plate to rotate polarization parallel
% Signals are pp and pr (b and d). Processing is the same as above.
% -----------------------------------------------------------------
tau2 = scat.LP(:,2)./(r(:,2).*scat.LREF(:,2));         
beam_c2 = -100/15*log(tau2);     
idx_c2 = beam_c2<0;
if any(idx_c2)
    beam_c2(idx_c2)  = 0.038;
end
if ~isempty(beamc)
    beam_c2 = repmat(beamc,length(tau2),1);    
end
% beam_c2 = 0.05;
att_factors2 = exp(-beam_c2*(paths/100));

pp = scat.pp./att_factors2;   
pp(pp<0)=nan;
pp_2 = pp;
pp = pp - scaled_zsc_pp./att_factors_z; 
pp(pp<0)=nan;
pp = pp .* repmat(sin(angles*pi/180),ns,1)./repmat(LREF_median,ns,nang);
pp_2 = pp_2 .* repmat(sin(angles*pi/180),ns,1)./repmat(LREF_median,ns,nang);

pr = scat.pr./att_factors2;   
pr(pr<0) = nan;
pr_2 = pr;
pr = pr - scaled_zsc_pr./att_factors_z;    
pr(pr<0)=nan;
pr = pr .* repmat(sin(angles*pi/180),ns,1)./repmat(LREF_median,ns,nang);
pr_2 = pr_2 .* repmat(sin(angles*pi/180),ns,1)./repmat(LREF_median,ns,nang);

proc.tau = [tau1 tau2];
proc.beamc = [beam_c1 beam_c2];

%% ========================================================================
% Process ring data to near-forward VSF
% =========================================================================

rings.radii = logspace(0,log10(200),33)*0.1;                      % mm
rings.edge_angles = asin(sin(atan(rings.radii/calfact.focal_length))/1.33);         % angles in water in radians  
dOmega=cos(rings.edge_angles(1:32))-cos(rings.edge_angles(2:33)); % find solid angle
dOmega=dOmega*2*pi/6;                                             % factor 6 takes care of rings covering only 1/6th circle

% Calculating scat and cscat according to "Processing LISST-100 and 
% LISST-100X data in MATLAB" (April 30, 2008)
scat1 = scat.rings1./repmat(tau1,1,32);
scat1_bulk = scat1;  % claculating bulk vsf for ring detectors
scat1 = scat1 - zsc_pmt_matched.rings1 .* repmat(scat.LREF(:,1)./zsc_pmt_matched.LREF(:,1),1,32);


cscat1 = scat1 .*repmat(calfact.dcal,ns,1) .* repmat(calfact.dvig,ns,1)  .*calfact.ND; 
cscat1_bulk = scat1_bulk .*repmat(calfact.dcal,ns,1) .* repmat(calfact.dvig,ns,1)  .*calfact.ND; 
cscat1(cscat1<0) = 0;
cscat1_bulk(cscat1_bulk<0) = 0;

% ring area correction, vignetting correction, ND filter transm. corr.
scat2 = scat.rings2./repmat(tau2,1,32);
scat2_bulk = scat2;
scat2 = scat2 - zsc_pmt_matched.rings2 .* repmat(scat.LREF(:,2)./zsc_pmt_matched.LREF(:,2),1,32);
cscat2 = scat2 .*repmat(calfact.dcal,ns,1) .* repmat(calfact.dvig,ns,1)  .*calfact.ND; 
cscat2_bulk = scat2_bulk .*repmat(calfact.dcal,ns,1) .* repmat(calfact.dvig,ns,1)  .*calfact.ND; 
cscat2(cscat2<0) = 0;
cscat2_bulk(cscat2_bulk<0) = 0;

light_on_rings1 = cscat1.*calfact.Watt_per_count_on_rings;
light_on_rings2 = cscat2.*calfact.Watt_per_count_on_rings;

light_on_rings1_bulk = cscat1_bulk.*calfact.Watt_per_count_on_rings;
light_on_rings2_bulk = cscat2_bulk.*calfact.Watt_per_count_on_rings;

% remove ring data that has very low signal
% light_on_rings1(scat1 < 25) = NaN;
% light_on_rings2(scat1 < 25) = NaN;
light_on_rings1(scat1 < 1) = NaN;
light_on_rings2(scat2 < 1) = NaN;
light_on_rings1_bulk(scat1_bulk < 1) = NaN;
light_on_rings2_bulk(scat2_bulk < 1) = NaN;

% calculate incident laser power from LREF
laser_incident_power1 = scat.LREF(:,1)*calfact.Watt_per_count_laser_ref; 
laser_incident_power2 = scat.LREF(:,2)*calfact.Watt_per_count_laser_ref;

% calculate forward scattering 
beam_bf = nansum(light_on_rings1./repmat(laser_incident_power1,1,32),2) .* 6;% factor 6 due to arcs
beam_bf_bulk = nansum(light_on_rings1_bulk./repmat(laser_incident_power1,1,32),2) .* 6;% factor 6 due to arcs
beam_bf = beam_bf ./ 0.15;
beam_bf_bulk = beam_bf_bulk ./ 0.15;

rings.beam_bf = beam_bf;    
rings.beam_bf_bulk = beam_bf_bulk;    

% calculate VSF for ring angles
rho = 200^(1/32);
rings.angles = rings.edge_angles(1:32)*sqrt(rho)./pi.*180;

rings.cscat1 = cscat1;
rings.cscat2 = cscat2;
rings.vsf1 = light_on_rings1 ./ repmat(dOmega,ns,1) ./ repmat(laser_incident_power1,1,32) ./ 0.15; % for 0.15m pathlength;
rings.vsf2 = light_on_rings2 ./ repmat(dOmega,ns,1) ./ repmat(laser_incident_power2,1,32) ./ 0.15;

rings.cscat1_bulk = cscat1_bulk;
rings.cscat2_bulk = cscat2_bulk;
rings.vsf1_bulk = light_on_rings1_bulk ./ repmat(dOmega,ns,1) ./ repmat(laser_incident_power1,1,32) ./ 0.15; % for 0.15m pathlength;
rings.vsf2_bulk = light_on_rings2_bulk ./ repmat(dOmega,ns,1) ./ repmat(laser_incident_power2,1,32) ./ 0.15;

proc.rings = rings;

%% ========================================================================
% Correct eye data for power change during scan and PMT relative gain
% =========================================================================

proc.angles = angles;

% Correct change in laser power 
rp(:,1:40) = rp(:,1:40).*calfact.laser_power_change_factor;
pp(:,1:40) = pp(:,1:40).*calfact.laser_power_change_factor;
rr(:,1:40) = rr(:,1:40).*calfact.laser_power_change_factor;
pr(:,1:40) = pr(:,1:40).*calfact.laser_power_change_factor;

rp_2(:,1:40) = rp_2(:,1:40).*calfact.laser_power_change_factor;
pp_2(:,1:40) = pp_2(:,1:40).*calfact.laser_power_change_factor;
rr_2(:,1:40) = rr_2(:,1:40).*calfact.laser_power_change_factor;
pr_2(:,1:40) = pr_2(:,1:40).*calfact.laser_power_change_factor;

% Interpolate over laser gain transition
rp(:,36:41) = interp1(angles([35 42]),rp(:,[35 42])',angles(36:41))';
pp(:,36:41) = interp1(angles([35 42]),pp(:,[35 42])',angles(36:41))';
rr(:,36:41) = interp1(angles([35 42]),rr(:,[35 42])',angles(36:41))';
pr(:,36:41) = interp1(angles([35 42]),pr(:,[35 42])',angles(36:41))';

rp_2(:,36:41) = interp1(angles([35 42]),rp_2(:,[35 42])',angles(36:41))';
pp_2(:,36:41) = interp1(angles([35 42]),pp_2(:,[35 42])',angles(36:41))';
rr_2(:,36:41) = interp1(angles([35 42]),rr_2(:,[35 42])',angles(36:41))';
pr_2(:,36:41) = interp1(angles([35 42]),pr_2(:,[35 42])',angles(36:41))';

% gemeteric correction for slight misalignment between laser and eyeball
% viewing plane. This part is comment out.
rp = rp.*repmat(polyval(calfact.geometric_cal_coeff,angles),ns,1);
rr = rr.*repmat(polyval(calfact.geometric_cal_coeff,angles),ns,1);
pp = pp.*repmat(polyval(calfact.geometric_cal_coeff,angles),ns,1);
pr = pr.*repmat(polyval(calfact.geometric_cal_coeff,angles),ns,1);

rp_2 = rp_2.*repmat(polyval(calfact.geometric_cal_coeff,angles),ns,1);
rr_2 = rr_2.*repmat(polyval(calfact.geometric_cal_coeff,angles),ns,1);
pp_2 = pp_2.*repmat(polyval(calfact.geometric_cal_coeff,angles),ns,1);
pr_2 = pr_2.*repmat(polyval(calfact.geometric_cal_coeff,angles),ns,1);


% net scattering in counts corrected for power step during scan
proc.rp = rp;
proc.rr = rr;
proc.pp = pp;
proc.pr = pr;

% User did not provide alpha as input parameter, so estimate it using data
if isempty(alpha)
    idx45 = find(angles==45, 1,'first');
    idx135 = find(angles==135, 1,'first');
    alpha_ac45 = rp_2(:,idx45)./rr_2(:,idx45);
    alpha_ac135 = rp_2(:,idx135)./rr_2(:,idx135);
    alpha_bd45 = pp_2(:,idx45)./pr_2(:,idx45);
    alpha_bd135 = pp_2(:,idx135)./pr_2(:,idx135);
    alpha = nanmedian([alpha_ac45; alpha_ac135; alpha_bd45; alpha_bd135]);
end

proc.alpha = alpha;

% c and d signals corrected for relative gain between PMT's
rr_scaled = rr*alpha;
pr_scaled = pr*alpha;
rr2_scaled = rr_2*alpha;
pr2_scaled = pr_2*alpha;


proc.rr_scaled = rr_scaled;
proc.pr_scaled = pr_scaled;

%% ========================================================================
% Calculation of P11, P12, P22 elements of Mueller matrix. See users
% manual.
% =========================================================================

% rp is a, vertical laser beam and Horizental PMT
% pp is b, horizontal laser beam and horizontal PMT
% rr is c, vetical laser beam and vertical PMT
% pr is d, horizontal laser beam and vetical PMT
% laser beam is vertical without a half-wavelength plate
% PMT2 Gain = alpha* PMT1 Gain
p11 = 0.25*(rp+pp+rr_scaled+pr_scaled);
p11_uncal_all = 0.25*(rp_2+pp_2+rr2_scaled+pr2_scaled);
p11_blank = 0.25*(scaled_zsc_rp+scaled_zsc_rr+scaled_zsc_pp+scaled_zsc_pr);
p11_sample = 0.25*(scat.rp+scat.rr+scat.pp+scat.pr);

% P12
p12 = 0.25.*((pp-rp)+(pr_scaled-rr_scaled))./p11;
p12_uncal = 0.25.*((pp-rp)+(pr_scaled-rr_scaled));
p12_uncal_all = 0.25.*((pp_2-rp_2)+(pr2_scaled-rr2_scaled));

% Extract p22
phi = repmat(angles,ns,1)*pi/180;
e = rp.*(1+cos(2*phi));
f = pp.*(1-cos(2*phi));
g = rr_scaled.*(1-cos(2*phi));
h = pr_scaled.*(1+cos(2*phi));
% Two different estimates of P22
% p22_1=(2*p11+(e+f))./(1+cos(4*repmat(p,ns,1)*pi/180))./p11;
% p22_2=(2*p11+(g+h))./(1+cos(4*repmat(p,ns,1)*pi/180))./p11;
p22_1 = ((2*p11-(e+f))./(2*cos(2*phi).^2))./p11;
p22_2 = ((2*p11-(g+h))./(2*cos(2*phi).^2))./p11;

p22_uncal = (2*p11-(e+f))./(2*cos(2*phi).^2);
e2 = rp_2.*(1+cos(2*phi));
f2 = pp_2.*(1-cos(2*phi));
p22_uncal_all = (2*p11-(e2+f2))./(2*cos(2*phi).^2);

idx = (angles>=40&angles<=50)|(angles>=130&angles<=140);
p22_1(:,idx) = nan;
p22_2(:,idx) = nan;

p11_uncal = p11;

%% ========================================================================
% Merge the Ring and Eye VSF
% =========================================================================

% Interpolate and extrapolate rings vsf
rings_vsf = rings.vsf1;

ring_vsf_x = interp1(rings.angles(:,[31 32])', rings_vsf(:,[31 32])', 15, 'linear', 'extrap')';
p11_overlap = p11(:,angles == 15);
p11_scale_factor = median(ring_vsf_x ./ p11_overlap,2);
p11 = p11.*repmat(p11_scale_factor,1,nang);


% bulk
rings_vsf_bulk = rings.vsf1_bulk;
ring_vsf_x_bulk = interp1(rings.angles(:,[31 32])', rings_vsf_bulk(:,[31 32])', 15, 'linear', 'extrap')';
p11_overlap_bulk = p11_uncal_all(:,angles == 15);
p11_scale_factor_bulk = median(ring_vsf_x_bulk ./ p11_overlap_bulk,2);
p11_bulk = p11_uncal_all.*repmat(p11_scale_factor_bulk,1,nang);


proc.p11_scale_factor = p11_scale_factor;

% truncate eyeball data to the usuable range
idx = angles > 14 & angles < 156;
% idx = angles > 14 & angles < 151;
proc.angles = proc.angles(:,idx);
proc.rp = proc.rp(:,idx);
proc.rr = proc.rr(:,idx);
proc.pp = proc.pp(:,idx);
proc.pr = proc.pr(:,idx);
proc.rr_scaled = proc.rr_scaled(:,idx);
proc.pr_scaled = proc.pr_scaled(:,idx);
proc.p11 = p11(:,idx);
proc.p12 = p12(:,idx);
proc.p12_uncal = p12_uncal(:,idx);
proc.p22_1 = p22_1(:,idx);
proc.p22_2 = p22_2(:,idx);
proc.p11_uncal = p11_uncal(:,idx); % p11 without calibration using rings data
proc.p11_uncal_bulk = p11_uncal_all(:,idx); % p11 without calibration and without substracting background 
proc.p12_uncal_bulk = p12_uncal_all(:,idx);
proc.p22_uncal = p22_uncal(:,idx);
proc.p22_uncal_bulk = p22_uncal_all(:,idx);
% final merging of ring and eyeball angles and measurments
proc.Angles = [rings.angles proc.angles];
proc.P11 = [rings.vsf1 proc.p11];
proc.P11_bulk = [rings.vsf1_bulk p11_bulk(:,idx)];
proc.LREF=LREF_median;
proc.p11_blank = p11_blank(:,idx);
proc.p11_sample = p11_sample(:,idx);
% calculate scattering coefficent
beam_b = zeros(ns,1);
for i=1:ns
      beam_b(i) = beam_bf(i) + 2*pi*sum(proc.p11(i,:).*sin(proc.angles*pi/180)*(pi/180)); % integral of VSF over angles p
end

proc.beamb = beam_b;

proc.beamc(proc.beamc < 0) = 0;
proc.beamb(proc.beamb < 0) = 0;

%% ========================================================================
% Plot Data
% =========================================================================
if plotData % check if plot data flag is 1
    
if all(proc.qc_saturated)
    warning('Saturation detected in all measurements, no data will be plotted');
else
    disp('Plotting non-saturated data') 
fontSize = 12;
lineWidth = 0.5;

[~,scatFileName,scatFileExt] = fileparts(scatPath);
movegui(figure('units','normalized','outerposition',[0 0 0.8 0.8],'Name',[scatFileName scatFileExt]),'center')

% plot P11 (eyeball + rings)
subplot(221)
    plot(proc.Angles,proc.P11(~proc.qc_saturated,:), 'k-', 'LineWidth',lineWidth);
    set(gca, 'Box','on', 'YScale','log')
    xlabel('Angle [deg]','FontSize',fontSize,'FontWeight','Bold')
    ylabel('LISST-VSF p_{11}','FontSize',fontSize,'FontWeight','Bold')
    
% plot P12 (eyeball only)
subplot(222)
    plot(proc.angles,proc.p12(~proc.qc_saturated,:), 'k-', 'LineWidth',lineWidth);
    xlabel('Angle [deg]','FontSize',fontSize,'FontWeight','Bold')
    ylabel('LISST-VSF p_{12}','FontSize',fontSize,'FontWeight','Bold')
    
% plot P22 (eyeball only)
subplot(223)
    plot(proc.angles,proc.p22_1(~proc.qc_saturated,:), 'k-', 'LineWidth',lineWidth);
    xlabel('Angle [deg]','FontSize',fontSize,'FontWeight','Bold')
    ylabel('LISST-VSF p_{22}','FontSize',fontSize,'FontWeight','Bold')
    
% plot beam attenuation for both polarizations and scattering
subplot(224)
hold on
    plot(proc.beamc(:,1), 'b.-', 'LineWidth',lineWidth,'DisplayName','Beam Attenuation P');
    plot(proc.beamc(:,2), 'r.-', 'LineWidth',lineWidth,'DisplayName','Beam Attenuation R');
    beam_b_non_saturated = proc.beamb;
    beam_b_non_saturated(proc.qc_saturated) = NaN;
    plot(beam_b_non_saturated, 'k.-', 'LineWidth',lineWidth,'DisplayName','Scattering');
    xlabel('Measurement Number','FontSize',fontSize,'FontWeight','Bold')
    ylabel('IOP [m^{-1}]','FontSize',fontSize,'FontWeight','Bold')
    legend('Show');
    yscale = ylim(gca);
    yscaled = diff(yscale);
    if yscaled  < 0.5
       sc = (0.5 - yscaled)/2;
       ylim([yscale(1)-sc yscale(2)+sc]);
    end
hold off
end
end











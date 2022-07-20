
function dat = lisstvsf_readdat(filename)

dat.procver = '20161208';
dat.procdate = datestr(now);

Nangles = 150;
precord = 40 + Nangles*5 ; % partial record: one polarization scan of 2

[raw,nsets] = parseFile(filename); % parse binary data file

bat=raw(:,34);
PMT=raw(:,35);  % this is common to the entire set.

d = dir(filename);
dat.filename = filename;
dat.filedate = d.date;
dat.nangles = Nangles;
dat.nsets = nsets;
dat.batt_volts = bat;
dat.pmt_gain = PMT;

%% ========================================================================
% Preallocate arrays
% =========================================================================
  
rp_off = zeros(nsets,Nangles);   
rp_on = zeros(nsets,Nangles);  

pp_on = zeros(nsets,Nangles);
pp_off = zeros(nsets,Nangles);

pr_off = zeros(nsets,Nangles);
pr_on = zeros(nsets,Nangles); 

rr_off = zeros(nsets,Nangles);   
rr_on = zeros(nsets,Nangles);

angles1 = zeros(nsets,Nangles);
angles2 = zeros(nsets,Nangles);
rings1 = zeros(nsets,32);
rings2 = zeros(nsets,32);
lp1 = zeros(nsets,1);
lp2 = zeros(nsets,1);
lref1 = zeros(nsets,1);
lref2 = zeros(nsets,1);
depth1 = zeros(nsets,1);
depth2 = zeros(nsets,1);
temp1 = zeros(nsets,1);
temp2 = zeros(nsets,1);
date1 = zeros(nsets,2);
date2 = zeros(nsets,2);


%% ========================================================================
% Loop to parse binary data
% =========================================================================

for i=1:nsets
    
    % -----------------------------------------------------------------
    % First eyeball rotation is polarized perpendicular
    % -----------------------------------------------------------------
    is = 0; 
    ie = is+precord;
    rings1(i,:) = raw(i,is+1:is+32);
    lp1(i) = raw(i,is+33);
    lref1(i) = raw(i,is+36);  
    depth1(i) = raw(i,is+37);
    temp1(i) = raw(i,is+38);
    date1(i,:) = [raw(i,is+39) raw(i,is+40)];  
    angles1(i,:) = raw(i,is+41:5:ie);

    rp_on(i,:) =raw(i,is+42:5:ie);
    rp_off(i,:) = raw(i,is+43:5:ie);
    rr_on(i,:) = raw(i,is+44:5:ie);
    rr_off(i,:) = raw(i,is+45:5:ie);
    
    % -----------------------------------------------------------------
    % Second eyeball rotation is polarized parallel
    % -----------------------------------------------------------------
    is = is+precord;
    ie = ie+precord;
    rings2(i,:) = raw(i,is+1:is+32);
    lp2(i) = raw(i,is+33);
    lref2(i) = raw(i,is+36); 
    depth2(i) = raw(i,is+37);
    temp2(i) = raw(i,is+38);
    date2(i,:)= [raw(i,is+39) raw(i,is+40)];  
    angles2(i,:) = raw(i,is+41:5:ie);

    pp_on(i,:) =raw(i,is+42:5:ie);
    pp_off(i,:) = raw(i,is+43:5:ie);
    pr_on(i,:) = raw(i,is+44:5:ie);
    pr_off(i,:) = raw(i,is+45:5:ie);
    
end

MM = fix(date1(:,2)/100);
SS = date1(:,2)-100*MM;
DD = fix(date1(:,1)/100);
HH = date1(:,1)-100*DD;  
date1 = datenum(0,0,DD,HH,MM,SS);

MM = fix(date2(:,2)/100);
SS = date2(:,2)-100*MM;
DD = fix(date2(:,1)/100);
HH = date2(:,1)-100*DD; 
date2 = datenum(0,0,DD,HH,MM,SS);


angle_min = min([angles1(:,1); angles2(:,1)]);
angle_max = max([angles1(:,end); angles2(:,end)]);
angles = angle_min:angle_max;

dat.angles_idx = angles;

for i=1:nsets
    
    idx = angles1(i,:)~=0;
    rp_on(i,:) = interp1(angles1(i,idx),rp_on(i,idx), angles, 'linear',0);
    rp_off(i,:) = interp1(angles1(i,idx),rp_off(i,idx), angles, 'linear',0);
    rr_on(i,:) = interp1(angles1(i,idx),rr_on(i,idx), angles, 'linear',0);
    rr_off(i,:) = interp1(angles1(i,idx),rr_off(i,idx), angles, 'linear',0);
    
    idx = angles2(i,:)~=0;
    pp_on(i,:) = interp1(angles2(i,idx),pp_on(i,idx), angles, 'linear',0);
    pp_off(i,:) = interp1(angles2(i,idx),pp_off(i,idx), angles, 'linear',0);
    pr_on(i,:) = interp1(angles2(i,idx),pr_on(i,idx), angles, 'linear',0);
    pr_off(i,:) = interp1(angles2(i,idx),pr_off(i,idx), angles, 'linear',0);
end


dat.raw.rp_on  = rp_on;
dat.raw.rp_off = rp_off;
dat.raw.rr_on  = rr_on;
dat.raw.rr_off = rr_off;
dat.raw.pp_on  = pp_on;
dat.raw.pp_off = pp_off;
dat.raw.pr_on  = pr_on;
dat.raw.pr_off = pr_off;

dat.qc_saturated = rp_on>30000 | rr_on>30000 | pr_on>30000 | pp_on>30000;

dat.rp = rp_on-rp_off;
dat.rr = rr_on-rr_off;
dat.pp = pp_on-pp_off;
dat.pr = pr_on-pr_off;

dat.LP = [lp1 lp2];
dat.LREF = [lref1 lref2];
dat.rings1 = rings1;
dat.rings2 = rings2;
dat.depth = [depth1 depth2];
dat.tempC = [temp1 temp2];
dat.time = [date1 date2];
% LISSTVSF_READDAT   Read a LISST-VSF binary data file   
%   Reads a raw data file, UNCALIBRATED. bat is battery, PMT is photomultiplier control voltage
%
% Usage:
%       dat = lisstvsf_readdat(filename);



function [raw,nsets] = parseFile(filename)
disp(filename)
%% --------------------------------------------------------------------
% read in unsigned ring data
%----------------------------------------------------------------------
fid = fopen(filename,'r','b'); %open file for reading using big endian format
raw1 = fread(fid,'uint16');
fclose(fid);

Nangles = 150;
precord = 40 + Nangles*5 ; % partial record: one polarization scan of 2
record = 2*precord;  % 2 turns per set

nsets = fix(length(raw1)/record);  % number of sets of 2-turns per set
disp(['Number of sets of data found = ',num2str(nsets),' '])
if (nsets>=2)
    raw1=reshape(raw1,record,nsets);
end
raw1=raw1';
raw(:,1:40) = raw1(:,1:40);
raw(:,791:831) = raw1(:,791:831);

% ring data are semi-signed, negative values are possible
raw(raw(:,1:32)>40950) = raw(raw(:,1:32)>40950) - 65536;
logical_index = logical([zeros(nsets,790) raw(:,791:822)>40950]);
raw(logical_index) = raw(logical_index) - 65536;

%% --------------------------------------------------------------------
% read in signed eyeball data
%----------------------------------------------------------------------
fid = fopen(filename,'r','b'); %open file for reading using big endian format
raw2 = fread(fid,'int16');
fclose(fid);

if (nsets>=2)
    raw2=reshape(raw2,record,nsets);
end
raw2=raw2';
raw(:,41:790) = raw2(:,41:790);
raw(:,832:1580) = raw2(:,832:1580);

% depth and temperature are signed values
raw(:,[37 38 827 828]) = raw2(:,[37 38 827 828]);





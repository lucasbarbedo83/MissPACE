%% Quick look of data
%myinit
%setmydefault(16,2);

[file,path] = uigetfile('Z*.VSF','Select Blank','MultiSelect','off');
fbfile = fullfile(path,file);
% 'Program\NIST_496nm_background_SN1662.VSF';
isplot = 1;

[file,path] = uigetfile('V*.VSF',['Select VSFs ',file],'MultiSelect','on');
SN=1425; % used in SeaTrac and IOP package
% SN=1664; % Looks like May 14, 2021 IOP used SN 1664, but values were zeros.
if iscell(file)
    vsf = nan(numel(file),141);
    for i=1: numel(file)
        proc = get_lisst_vsf_v5(SN,fullfile(path,file{i}),fbfile,[],[],isplot);
        file0 = file{i};
        process_vsf_here(proc, file0);
    end
else
    proc = get_lisst_vsf_v5(SN,fullfile(path,file),fbfile,[],[],isplot);
    file0 = file;
    process_vsf_here(proc,file0);
end

function process_vsf_here(proc, file)
angle = proc.Angles;
vsf_sw = betasw_ZHH2009(517,15,angle,35); % assume Temp = 15, and salinity = 35 PSU.

% there are two ways to examine the results depending on your
% experiment design.
% If you have measured a background that is clearly defined, use
% the following by uncommenting it.
vsf_bulk = proc.P11_opt;  vsfp = vsf_bulk;  % subtract background file

% If you don't have a background that is clearly defined, for
% example, you take a measurement in the ocean, use the following
% by uncommenting it.
% vsf_bulk = proc.P11_opt_bulk;  vsfp = vsf_bulk-vsf_sw'; % subtract pure seawater

fh=figure(2);
[a1,a2,h1,h2]=log_linear_plot(angle,vsfp);
%     loglog(angle,vsfp,'-');
xlabel('Angle (Deg)');
ylabel(a1,'VSF (m^-^1 sr^-^1)');
% if you want to save the figure, uncomment the following
% fname = extractBetween(file,1,8);
% saveas(fh,fullfile(path,strcat(fname{1},'.png')));

% To save data for later process, uncomment the following
ang = angle';
ind = ang<=150;
ang = ang(ind);
data = vsfp(:,ind)';
save(fullfile(path,strcat(fname{1},'.mat')),'ang','data')
end
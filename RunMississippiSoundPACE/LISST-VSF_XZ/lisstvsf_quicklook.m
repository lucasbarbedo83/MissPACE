%% Quick look of LISST-VSF data
clear all
close all

[file,path] = uigetfile('Z*.VSF','Select Blank','MultiSelect','on');
fbfile = fullfile(path,file);
% 'Program\NIST_496nm_background_SN1662.VSF';
isplot = 1;

[file,path] = uigetfile('V*.VSF',['Select VSFs ',file],'MultiSelect','on');

if iscell(file)
    vsf = nan(numel(file),141);
    for i=1: numel(file)
        proc = get_lisst_vsf_v4(fullfile(path,file{i}),fbfile,[],isplot);
        vsf(i,:) = median(proc.p11_absCal,1,'omitnan');
        pause;
        close(gcf);
    end
else
    proc = get_lisst_vsf_v4(fullfile(path,file),fbfile,[],isplot);
    vsf(:) = median(proc.p11_absCal,1,'omitnan');
end

figure
semilogy(proc.angles,vsf);

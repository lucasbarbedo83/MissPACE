function [ScatMatrix,PACE]=readLisstVSF_MissPACE(fname,ui)
%fname='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022';
% Read IOPS from Mississippi Sound PACE Mission
%   Input
%   fname: complet path to NRL_USM_PACE_2022
%          for example, in my laptop:
%          fname='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022'
%   ui: 1 to select one file manually
%         0 to run all PACE dataset
%   Output
%   lisstvsf: Mueller Matrix (VSF, P12, P21, P22, b, c, etc...)
%   PACE: lisst vsf, inside a datab structured.
%
% Lucas Barbedo, 19 July 2022
%% process DH4
 datab=[fname,'\IOP_data\PACE_lisstvsf_data_exfil'];
if ui==1
    [filename, pathname] = uigetfile( ...
        {'*.VSF',  'IOPs (*.VSF)'}, ...
        'Pick a file', datab,...
        'MultiSelect', 'off');
    n(1).name=filename;
    n(1).folder=pathname;
else
    str=[datab,'\V*.VSF'];
    n=dir(str);
end

%% processing,...
zsc=[datab,'\Z0331050.VSF'];
for i=1:76%length(n)
    disp(i)
    scatName=[n(i).folder,'\',n(i).name];
    disp(scatName)
    disp(zsc)
    proc = lisstvsf_makep(1425,scatName,zsc,[],0) ;
    % pause
    PACE.proc=proc;

end
ScatMatrix=proc;
end
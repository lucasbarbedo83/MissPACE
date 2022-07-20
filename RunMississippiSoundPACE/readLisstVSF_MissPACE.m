function [ScatMatrix,PACE]=readLisstVSF_MissPACE(fname,ui,FNlisstvsf)
%fname='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022';
% Read IOPS from Mississippi Sound PACE Mission
%   Input
%   fname: complet path to NRL_USM_PACE_2022
%          for example, in my laptop:
%          fname='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022'
%
%   ui:    0 to run all PACE dataset
%          1 to select one file manually
%          2 to select the sequence of files to be read
%
%   FNlisstvsf: list of files names of LISSTVSF to be read.
%   the processing PACE(n. station) respect the station order 
%  it is good to evaluate dataset and cross multi-sensor validation.
%
%   Output
%   lisstvsf: Mueller Matrix (VSF, P12, P21, P22, b, c, etc...)
%   PACE: lisst vsf, inside a datab structured.
%
% Lucas Barbedo, 19 July 2022
%% process DH4
datab=[fname,'\IOP_data\PACE_lisstvsf_data_exfil'];
if ui==0
    str=[datab,'\V*.VSF'];
    n=dir(str);
elseif ui==1
    [filename, pathname] = uigetfile( ...
        {'*.VSF',  'IOPs (*.VSF)'}, ...
        'Pick a file', datab,...
        'MultiSelect', 'off');
    n(1).name=filename;
    n(1).folder=pathname;
elseif ui==2
    for i=1:length(FNlisstvsf)
        n(i).name=char(FNlisstvsf(i));
        n(i).folder=datab;
    end
end

%% processing,...
zsc=[datab,'\Z0331050.VSF'];
for i=1:length(n)%76%
    scatName=[n(i).folder,'\',n(i).name];
    % disp(scatName)
    if exist(scatName,'file')
        proc = lisstvsf_makep(1425,scatName,zsc,[],0) ;
        % pause
        PACE(i).proc=proc;
    else
        PACE(i).proc.p11=nan;
        proc=nan;
    end

end
ScatMatrix=proc;



figure
axes('Box','on')
hold on
if length(n)>1
    ccolor=parula(length(n));
else
    ccolor=[0 0 0];
end
for i=1:length(n)
    if ~isnan( PACE(i).proc.p11)
    semilogy(PACE(i).proc.angles,median(PACE(i).proc.p11,1),'-','Color',ccolor(i,:))
    end
end
xlabel('$\theta$: scattering angle','Interpreter','latex')
ylabel('$\beta(517)$: volume scattering function, in  m$^{-1}$ sr$^{-1}$','Interpreter','latex')
hold off
grid minor


end
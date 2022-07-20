function [Rrs,PACE]=readRrs_MissPACE(fname,ui,fieldnotes)
% Read Rrs from Mississippi Sound PACE Mission
%   Input
%   fname: complet path to NRL_USM_PACE_2022
%          for example, in my laptop:
%          fname='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022'
%
%   ui:  0 to run all PACE dataset
%        1 to select one file manually
%        2 to select the sequence of files to be read
%
%   FNradiometry: list of files names of above water radiometry
%   the processing PACE(n. station) respect the station order
%   it is good to evaluate dataset and cross multi-sensor validation.
%
%   Output
%   Rrs: above water radiometry from selected file or last one
%   PACE: above water inside a datab structured.
%
% Lucas Barbedo, 19 July 2022

datab=[fname,'\RemoteSensingReflectance'];
if ui==0
    str=[datab,'\*.csv'];
    n=dir(str);
elseif ui==1
    [filename, pathname] = uigetfile( ...
        {'*.csv',  'Rrs (*.csv)'}, ...
        'Pick a file', datab,...
        'MultiSelect', 'off');
    n(1).name=filename;
    n(1).folder=pathname;
elseif ui==2
    %% the above water radiometry must be time and spatial coincident
    % Stations based on IOPs
    stationjd= datenum(fieldnotes.ObsStartTime);
    stationlon=fieldnotes.Lon;
    stationlat=fieldnotes.Lat;
    % samples of Rrs
    str=[datab,'\*.csv'];
    samplesRrs=dir(str) ;
    for niops=1:length(stationjd)
        for nrrs=1:length(samplesRrs)
            str=[samplesRrs(nrrs).folder,'\',samplesRrs(nrrs).name];
            info = infoRrs_MissPACE(str);
            jd=info.Value(8);
            lat=info.Value(4);
            lon=info.Value(5);
            [dist,phaseangle] = sw_dist([lat stationlat(niops)],[ lon stationlon(niops)] ,'Km');
            deltad(nrrs)=dist;
            deltat(nrrs)=abs(stationjd(niops)-jd);
            deltaVelo(nrrs)=dist./abs(stationjd(niops)-jd);
        end

        % [indD,posid]=min(deltad,2,'first');
        % [indT,posit]=min(deltat);
        [indVelo,posiVelo]=min(deltaVelo);

        if deltad(posiVelo)<=1 && deltat(posiVelo)<=0.25 % 1 km and 6 hours.
            n(niops).name=samplesRrs(posiVelo).name;
            n(niops).folder=samplesRrs(posiVelo).folder;
        else
            n(niops).name='nan';
            n(niops).folder=samplesRrs(posiVelo).folder;
        end
    end
end

opts = delimitedTextImportOptions("NumVariables", 13);
% Specify range and delimiter
opts.DataLines = [27, Inf];
opts.Delimiter = ",";
% Specify column names and types
opts.VariableNames = ["wl", "MeanRrs", "MedianRrs", "StdDevRrs", "MeanEd0plus", "MedianEd0plus", "StdDevEd0plus", "MeanLw0plus", "MedianLw0plus", "StdDevLw0plus", "MeanLd0plus", "MedianLd0plus", "StdDevLd0plus"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

for i=1:length(n)
    %% Set up the Import Options and import the data
    str=[n(i).folder,'\',n(i).name];
    if exist(str,'file')
        Rrs  = readtable(str, opts);
        PACE(i).Rrs=Rrs;
    else
        Rrs=nan;
        PACE(i).Rrs=[];
    end
    %% Clear temporary variables
end
clear opts
figure
axes('YScale','log','YLim',[0.00001 0.1],'XLim',[350 1000],...
    'Box','on')
hold on
if length(n)>1
    ccolor=parula(length(n));
else
    ccolor=[0 0 0];
end
for i=1:length(n)
    if ~isempty(PACE(i).Rrs)
    plot(PACE(i).Rrs.wl,PACE(i).Rrs.MeanRrs,'-','Color',ccolor(i,:))
    end
end
xlabel('$\lambda$: wavelength, in nm','Interpreter','latex')
ylabel('$R_{rs}(\lambda)$: Remote sensing reflectance, in sr$^{-1}$','Interpreter','latex')
hold off
grid minor
end

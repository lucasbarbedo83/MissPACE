function [Rrs,PACE]=readRrs_MissPACE(varargin)
% Read Rrs from Mississippi Sound PACE Mission
%   Input
%   ui: 1 to select one file manually
%         0 to run all PACE dataset
%   Output
%   Rrs: above water radiometry from selected file or last one
%   PACE: above water inside the database structured.
%
% Lucas Barbedo, 19 July 2022

m=nargin;
database='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022\RemoteSensingReflectance';

if m>0
     ui=cell2mat(varargin(1));
    [filename, pathname] = uigetfile( ...
        {'*.csv',  'Rrs (*.csv)'}, ...
        'Pick a file', database,...
        'MultiSelect', 'off');
    n(1).name=filename;
    n(1).folder=pathname;
elseif m==0
    str=[database,'\*.csv'];
    n=dir(str);
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
    disp(str)
    Rrs  = readtable(str, opts); 
    PACE(i).Rrs=Rrs;

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
plot(PACE(i).Rrs.wl,PACE(i).Rrs.MeanRrs,'-','Color',ccolor(i,:))
end
xlabel('$\lambda$: wavelength, in nm','Interpreter','latex')
ylabel('$R_{rs}(\lambda)$: Remote sensing reflectance, in sr$^{-1}$','Interpreter','latex')
hold off
grid minor
end

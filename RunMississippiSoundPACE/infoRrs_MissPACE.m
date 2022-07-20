function [INFOradiometry] = infoRrs_MissPACE(filename)
%% Import data from text file
% Script for importing data from the following text file:
%
%   For exampe:
%      filename='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022\RemoteSensingReflectance\NRL_USM_PACE_20220208_LittleDog7_Rrs_1.csv'
%
% Auto-generated by MATLAB on 20-Jul-2022 16:39:32

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 13);

% Specify range and delimiter
opts.DataLines = [9, 15; 17, 17];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Name", "Value"];
opts.SelectedVariableNames = ["Name", "Value"];
opts.VariableTypes = ["string", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties

% Import the data
INFOradiometry = readtable(filename, opts);
INFOradiometry.Name(1)='Date' ;
INFOradiometry.Name(2)='local' ;
INFOradiometry.Name(3)='UTC' ;
INFOradiometry.Name(4)='lat' ;
INFOradiometry.Name(5)='lon' ;
INFOradiometry.Name(6)='depth' ;
INFOradiometry.Name(7)='windspeed' ;
INFOradiometry.Name(8)='jd';
str=[num2str(INFOradiometry.Value(1)),' ',num2str(INFOradiometry.Value(3))];
INFOradiometry.Value(8)=datenum(str,'yyyymmdd HHMM');

%% Clear temporary variables
clear opts
end
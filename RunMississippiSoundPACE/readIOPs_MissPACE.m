function [IOPs,PACE]=readIOPs_MissPACE(varargin)
% Read IOPS from Mississippi Sound PACE Mission
%   Input(varatgin)
%   ui: 1 to select one file manually
%         0 to run all PACE dataset
%   Output
%   IOPs: inherent optical properties
%
% Lucas Barbedo, 19 July 2022

%% to open miliQ water sample, blank
acsMiliQ=load('acs349_miliQ\acs_349_20220715_blank.mat');
wl_a=acsMiliQ.ac_blank(1,:);
wl_c=acsMiliQ.ac_blank(3,:);

%% process DH4
m=nargin;
database='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022\IOP_data';

if m>0
    ui=cell2mat(varargin(1));
    [filename, pathname, filterindex] = uigetfile( ...
        {'*',  'IOPs (archive_pair_03.*)'}, ...
        'Pick a file', database,...
        'MultiSelect', 'on');
    n(1).name=filename;
    n(1).folder=pathname;
else
    str=[database,'\archive_pair_03.*'];
    n=dir(str);
end
for i=1:length(n)
    str=[n(1).folder,'/',n(1).name];
    IOPs = IOP_pair3_import(str);
    dh4=table2array(IOPs);
    PACE(i).ctL0=dh4(:,6:90+5);
    PACE(i).atL0=dh4(:,96:96+89);
    PACE(i).cp=dh4(:,6:90+5)-acsMiliQ.ac_blank(4,:);
    PACE(i).ap=dh4(:,96:96+89)-acsMiliQ.ac_blank(2,:);
    PACE(i).z=dh4(:,187);
    PACE(i).T=dh4(:,2);
    PACE(i).P=dh4(:,4);
    PACE(i).S=dh4(:,5);
    PACE(i).dt=dh4(:,5)./1000;%seconds.
    PACE(i).wl_a=wl_a;
    PACE(i).wl_c=wl_c;
  
end
%disp(IOPs.Properties.VariableNames')
end
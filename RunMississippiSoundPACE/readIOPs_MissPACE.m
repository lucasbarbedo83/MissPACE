function [IOPs,PACE]=readIOPs_MissPACE(fname,ui,FNiops)
% Read IOPS from Mississippi Sound PACE Mission
%   Input 
%   fname: complet path to NRL_USM_PACE_2022
%          for example, in my laptop:
%          fname='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022'
%
%   ui:  0 to run all PACE dataset 
%        1 to select one file manually
%        2 to select the sequence of files to be read
%
%   FNiops: list of files names of IOPs archive files extracted form dh4.
%   the processing PACE(n. station) respect the station order 
%   it is good to evaluate dataset and cross multi-sensor validation.
%   
%   Output
%   IOPs: inherent optical properties
%   PACE: inwater inside a datab structured.
% Lucas Barbedo, 19 July 2022
% University of Southern Mississippim, USM
% Xiaodong Zhang Lab
%% to open miliQ water sample, blank
acsMiliQ=load('acs349_miliQ\acs_349_20220715_blank.mat');
wl_a=acsMiliQ.ac_blank(1,:);
wl_c=acsMiliQ.ac_blank(3,:);
%% process DH4
datab=[fname,'\IOP_data'];
if ui==0
    str=[datab,'\archive_pair_03.*'];
    n=dir(str);
elseif ui==1
    [filename, pathname] = uigetfile( ...
        {'*',  'IOPs (archive_pair_03.*)'}, ...
        'Pick a file', datab,...
        'MultiSelect', 'on');
    n(1).name=filename;
    n(1).folder=pathname;
elseif ui==2
   for i=1:length(FNiops)
       str=[datab,'\',char(FNiops(i))];
       thissample=dir(str);
       if  ~isempty(thissample)
        n(i)=dir(str);
       else
        n(i).name='NoData';
        n(i).folder=datab;
       end
   end
end
for i=1:length(n)
    str=[n(i).folder,'\',n(i).name];
    if  exist(str,'file')
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
    else
        PACE(i).T=nan;
        IOPs=nan;
    end

end
%disp(IOPs.Properties.VariableNames')
end
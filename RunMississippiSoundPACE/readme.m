clear all
close all
clc
disp('Processing marine bio-optics ')
disp('Mississippi Sound, Early PACE mission')
disp('Xiaodong Lab')
disp('School of Ocean Science and Engineering')
disp('University of Southern Mississippi, USM')
disp('https://www.usm.edu/ocean-science-engineering/')
disp('see:')
disp('https://github.com/lucasbarbedo83/MissPACE')

addpath LISST-VSF_XZ\
%% Acesss field notes, in case of ui = 2
%  the processing PACE(n. station) respect the station order 
%  it is good to evaluate dataset and cross multi-sensor validation.
%Workbook='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022\IOP_data\NRL_USM_PACE_IOP-CastLog.xlsx';
Workbook='NRL_USM_PACE_IOP-SampleLog.xlsx'
[fieldnotes]=importFieldNoteMissPACE(Workbook);
FNlisstvsf=fieldnotes.VSF_file;
FNiops=fieldnotes.archive_file;


fname='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022';
datab=[fname,'\IOP_data\PACE_lisstvsf_data_exfil'];
 
%% run LisstVSF, IOPs, and Rrs;
ui=2;
[Rrs,abovewater]=readRrs_MissPACE(fname,ui,fieldnotes);
[IOPs,inwater]=readIOPs_MissPACE(fname,ui,FNiops);
[ScatMatrix,PACE]=readLisstVSF_MissPACE(fname,ui,FNlisstvsf);

ccolor=parula(length(PACE));
figure('Units','centimeters','Position',[1 1 6 18.5]) 
ax1=axes('position',[0.9 0.9 0.1 0.1],'box', 'off', 'visible','off',...
    'xtick',[],'ytick',[],'xcolor','none','ycolor','none');
hold on
title('Legend:')
for i=1:length(PACE)
 str=[char(fieldnotes.StationName(i)),', ',char(fieldnotes.ObsStartTime(i))];
 plot(ax1,i,1,'-','color',ccolor(i,:),'linewidth',2,'DisplayName',str)
end
h=legend(ax1,'-DynamicLegend');
set(h, 'LimitMaxLegendEntries', false,'Fontname','sans','Fontsize',6);
hold off

%% Scatter plot acs vs. lisstvsf
[bp517] = CrossEvaluation_bp(PACE,inwater);
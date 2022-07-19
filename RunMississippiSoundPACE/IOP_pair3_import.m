function stn_name = IOP_pair3_import(filename, dataLines)
% Import data from an IOP station file.
%  StationName = IOP_pair3_import(FILENAME) reads data from text file FILENAME for
%   the default selection.  Returns the data as a table.
%
%  StationName = IOP_pair3_import(FILE, DATALINES) reads data for the specified
%   row interval(s) of FILE. Specify DATALINES as a positive scalar integer 
%   or a N-by-2 array of positive scalar integers for discontiguous row intervals.
%
%  Example:
%  stnLDK7 = IOP_pair3_import("C:\Processing\IOP_merged\archive_pair_03.000", [2, Inf]);
%
% T.Wissing 23-Feb-2022
% Lucas Barbedo, 19 July 2022

%% Input handling
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 201);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["millisec", "ctd_t_C", "ctdcond", "ctd_dbar", "ctd_sal", "c4001", "c4039", "c4071", "c4107", "c4134", "c4169", "c4205", "c4245", "c4283", "c4320", "c4355", "c4391", "c4430", "c4470", "c4518", "c4557", "c4599", "c4636", "c4676", "c4720", "c4766", "c4811", "c4855", "c4894", "c4934", "c4972", "c5012", "c5052", "c5095", "c5139", "c5182", "c5227", "c5265", "c5304", "c5341", "c5378", "c5416", "c5454", "c5492", "c5533", "c5570", "c5610", "c5649", "c5688", "c5713", "c5742", "c5776", "c5805", "c5842", "c5874", "c5913", "c5952", "c5994", "c6036", "c6078", "c6120", "c6162", "c6204", "c6245", "c6280", "c6324", "c6363", "c6400", "c6440", "c6484", "c6530", "c6572", "c6609", "c6652", "c6695", "c6736", "c6780", "c6817", "c6851", "c6889", "c6926", "c6962", "c7000", "c7036", "c7075", "c7112", "c7148", "c7182", "c7222", "c7256", "c7300", "c7335", "c7372", "c7413", "c7454", "a4012", "a4048", "a4084", "a4112", "a4142", "a4174", "a4214", "a4251", "a4290", "a4324", "a4363", "a4396", "a4437", "a4478", "a4522", "a4560", "a4603", "a4638", "a4678", "a4720", "a4766", "a4808", "a4853", "a4894", "a4932", "a4968", "a5008", "a5048", "a5090", "a5132", "a5177", "a5218", "a5258", "a5298", "a5335", "a5372", "a5408", "a5445", "a5483", "a5521", "a5562", "a5599", "a5637", "a5676", "a5712", "a5741", "a5773", "a5803", "a5834", "a5868", "a5905", "a5945", "a5981", "a6023", "a6066", "a6108", "a6149", "a6192", "a6230", "a6270", "a6307", "a6349", "a6385", "a6423", "a6464", "a6511", "a6552", "a6593", "a6632", "a6678", "a6716", "a6758", "a6799", "a6833", "a6869", "a6907", "a6941", "a6979", "a7013", "a7049", "a7087", "a7123", "a7158", "a7196", "a7235", "a7272", "a7307", "a7343", "a7381", "a7421", "acs_t_C", "No3", "Depth", "acs_xt_C", "No4", "No5", "Em_Chl", "Chl_count", "Em_PE", "PE_count", "Em_CDOM", "FDOM_count", "No_Header12", "DO_col1", "DO_col2", "DO_col3"];
% opts.VariableNames = ["millisec", "ctd_t_C", "No_Header2", "ctd_dbar", "ctd_sal", "c4001", "c4039", "c4071", "c4107", ...
%     "c4134", "c4169", "c4205", "c4245", "c4283", "c4320", "c4355", "c4391", "c4430", "c4470", "c4518", "c4557", "c4599", ...
%     "c4636", "c4676", "c4720", "c4766", "c4811", "c4855", "c4894", "c4934", "c4972", "c5012", "c5052", "c5095", "c5139", ...
%     "c5182", "c5227", "c5265", "c5304", "c5341", "c5378", "c5416", "c5454", "c5492", "c5533", "c5570", "c5610", "c5649", ...
%     "c5688", "c5713", "c5742", "c5776", "c5805", "c5842", "c5874", "c5913", "c5952", "c5994", "c6036", "c6078", "c6120", ...
%     "c6162", "c6204", "c6245", "c6280", "c6324", "c6363", "c6400", "c6440", "c6484", "c6530", "c6572", "c6609", "c6652", ...
%     "c6695", "c6736", "c6780", "c6817", "c6851", "c6889", "c6926", "c6962", "c7000", "c7036", "c7075", "c7112", "c7148", ...
%     "c7182", "c7222", "c7256", "c7300", "c7335", "c7372", "c7413", "c7454", "a4012", "a4048", "a4084", "a4112", "a4142", ...
%     "a4174", "a4214", "a4251", "a4290", "a4324", "a4363", "a4396", "a4437", "a4478", "a4522", "a4560", "a4603", "a4638", ...
%     "a4678", "a4720", "a4766", "a4808", "a4853", "a4894", "a4932", "a4968", "a5008", "a5048", "a5090", "a5132", "a5177", ...
%     "a5218", "a5258", "a5298", "a5335", "a5372", "a5408", "a5445", "a5483", "a5521", "a5562", "a5599", "a5637", "a5676", ...
%     "a5712", "a5741", "a5773", "a5803", "a5834", "a5868", "a5905", "a5945", "a5981", "a6023", "a6066", "a6108", "a6149", ...
%     "a6192", "a6230", "a6270", "a6307", "a6349", "a6385", "a6423", "a6464", "a6511", "a6552", "a6593", "a6632", "a6678",  ...
%     "a6716", "a6758", "a6799", "a6833", "a6869", "a6907", "a6941", "a6979", "a7013", "a7049", "a7087", "a7123", "a7158", ...
%     "a7196", "a7235", "a7272", "a7307", "a7343", "a7381", "a7421", "acs_t_C", "No3", "Depth", "acs_xt_C", "No4", "No5", ...
%     "Em_Chl", "Chl_count", "Em_PE", "PE_count", "Em_FDOM", "FDOM_count", "No_Header12", "oxy1", "oxy2", "oxy3"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.SelectedVariableNames= ["millisec", "ctd_t_C", "ctdcond", "ctd_dbar", "ctd_sal", "c4001", "c4039", "c4071", "c4107", "c4134", "c4169", "c4205", "c4245", "c4283", "c4320", "c4355", "c4391", "c4430", "c4470", "c4518", "c4557", "c4599", "c4636", "c4676", "c4720", "c4766", "c4811", "c4855", "c4894", "c4934", "c4972", "c5012", "c5052", "c5095", "c5139", "c5182", "c5227", "c5265", "c5304", "c5341", "c5378", "c5416", "c5454", "c5492", "c5533", "c5570", "c5610", "c5649", "c5688", "c5713", "c5742", "c5776", "c5805", "c5842", "c5874", "c5913", "c5952", "c5994", "c6036", "c6078", "c6120", "c6162", "c6204", "c6245", "c6280", "c6324", "c6363", "c6400", "c6440", "c6484", "c6530", "c6572", "c6609", "c6652", "c6695", "c6736", "c6780", "c6817", "c6851", "c6889", "c6926", "c6962", "c7000", "c7036", "c7075", "c7112", "c7148", "c7182", "c7222", "c7256", "c7300", "c7335", "c7372", "c7413", "c7454", "a4012", "a4048", "a4084", "a4112", "a4142", "a4174", "a4214", "a4251", "a4290", "a4324", "a4363", "a4396", "a4437", "a4478", "a4522", "a4560", "a4603", "a4638", "a4678", "a4720", "a4766", "a4808", "a4853", "a4894", "a4932", "a4968", "a5008", "a5048", "a5090", "a5132", "a5177", "a5218", "a5258", "a5298", "a5335", "a5372", "a5408", "a5445", "a5483", "a5521", "a5562", "a5599", "a5637", "a5676", "a5712", "a5741", "a5773", "a5803", "a5834", "a5868", "a5905", "a5945", "a5981", "a6023", "a6066", "a6108", "a6149", "a6192", "a6230", "a6270", "a6307", "a6349", "a6385", "a6423", "a6464", "a6511", "a6552", "a6593", "a6632", "a6678", "a6716", "a6758", "a6799", "a6833", "a6869", "a6907", "a6941", "a6979", "a7013", "a7049", "a7087", "a7123", "a7158", "a7196", "a7235", "a7272", "a7307", "a7343", "a7381", "a7421", "acs_t_C", "Depth", "acs_xt_C", "Em_Chl", "Chl_count", "Em_PE", "PE_count", "Em_CDOM", "FDOM_count", "DO_col1", "DO_col2", "DO_col3"];

% Import the data
stn_name = readtable(filename, opts);

%% Apply ECO FL-3 Correctors
% from FL3-531.dev File (hard set by TW, Feb 2022):
CHL_scalefactor = 0.0119;
CHL_darkcount = 61; 
Pe_scalefactor = 0.0422;
Pe_darkcount = 50;
CDOM_scalefactor = 0.1745;
CDOM_darkcount = 39;

stn_name.Chl = (stn_name.Chl_count - CHL_darkcount).*CHL_scalefactor;
stn_name.Pe = (stn_name.PE_count - Pe_darkcount).*Pe_scalefactor;
stn_name.FDOM = (stn_name.FDOM_count - CDOM_darkcount).*CDOM_scalefactor;

clear CHL_scalefactor
clear CHL_darkcount
clear Pe_scalefactor
clear Pe_darkcount
clear CDOM_scalefactor
clear CDOM_darkcount
clear dataLines

%% RAW Plots
% figure('Name',filename); 
% subplot(2,2,1);
% plot(stn_name.DO_col1,'DisplayName','stn_name.DO_col1');
% hold on;
% plot(stn_name.DO_col2,'DisplayName','stn_name.DO_col2');
% plot(stn_name.DO_col3,'DisplayName','stn_name.DO_col3');
% title('Dissolved Ox - Raw');
% legend('col1','col2','col3');
% subplot(2,2,2);
% plot(stn_name.Chl,'DisplayName','stn_name.Chl');
% hold on;
% plot(stn_name.Pe,'DisplayName','stn_name.Pe');
% plot(stn_name.FDOM,'DisplayName','stn_name.FDOM');
% legend('CHL-A','PE','FDOM');
% subplot(2,2,3:4);
% plot(stn_name.acs_t_C,'DisplayName','stn_name.acs_t_C');
% hold on;
% plot(stn_name.acs_xt_C,'DisplayName','stn_name.acs_xt_C');
% plot(stn_name.ctd_t_C,'DisplayName','stn_name.ctd_t_C');
% plot(stn_name.ctd_dbar,'DisplayName','stn_name.ctd_dbar');
% plot(stn_name.ctd_sal,'DisplayName','stn_name.ctd_sal');
% title('Quickview raw time-series');
% legend('-DynamicLegend','ACS T','ACS XT','CTD T','CTD dbar','CTD sal');
% 
% saveas(gcf,strcat(filename,'_timeseries','.png'));
% 
% figure('Name',filename); 
% subplot(1,3,1);
% plot(stn_name.DO_col1,-stn_name.ctd_dbar,'DisplayName','stn_name.DO_col1');
% hold on;
% plot(stn_name.DO_col2,-stn_name.ctd_dbar,'DisplayName','stn_name.DO_col2');
% plot(stn_name.DO_col3,-stn_name.ctd_dbar,'DisplayName','stn_name.DO_col3');
% title('Dissolved Ox - Raw');
% legend('col1','col2','col3');
% subplot(1,3,2);
% plot(stn_name.Chl,-stn_name.ctd_dbar,'DisplayName','stn_name.Chl');
% hold on;
% plot(stn_name.Pe,-stn_name.ctd_dbar,'DisplayName','stn_name.Pe');
% plot(stn_name.FDOM,-stn_name.ctd_dbar,'DisplayName','stn_name.FDOM');
% legend('CHL-A','PE','FDOM');
% subplot(1,3,3);
% plot(stn_name.acs_t_C,-stn_name.ctd_dbar,'DisplayName','stn_name.acs_t_C');
% hold on;
% plot(stn_name.acs_xt_C,-stn_name.ctd_dbar,'DisplayName','stn_name.acs_xt_C');
% plot(stn_name.ctd_t_C,-stn_name.ctd_dbar,'DisplayName','stn_name.ctd_t_C');
% plot(stn_name.ctd_sal,-stn_name.ctd_dbar,'DisplayName','stn_name.ctd_sal');
% title('Quickview raw pressure-series');
% legend('-DynamicLegend','ACS T','ACS XT','CTD T','CTD sal');
% 
% %saveas(gcf,strcat(filename,'_profiles','.png'));

%% Raw AC-S conversion
% Set absorption spectrum
lambda_a = [401, 404, 408,411, 414, 417, 421, 425, 429, 432, 436,439, 443, 447, 452,456,...
    460, 463,467, 472, 477, 480, 485, 489,493,496, 500, 504, 509,513, 517,521, 525,...
    530, 534, 537, 541, 545, 548, 552, 556, 560, 564, 567, 571, 574, 577, 580, 583, ...
    587, 591, 595, 598, 602, 607, 611, 615, 619, 623, 627, 631, 635, 639, 642, 646, 651, 655, ...
    659, 663, 667, 672, 676, 680, 683, 687, 691, 694, 698, 701,704, 708, 712, 716, 719, 724, 727, 730, 734, 738, 742];
% Set attenuation spectrum
lambda_c = [400,  403,  407,  410,  413,  416,  420,  424,  428,  432,  435,  439,  443,  447,  451,  455,...
    459,  463,  467,  472,  476,  481,  485,  489,  493,  497,  501,  505,  509,  513,    518, 522, 526,...
    530, 534, 537, 541, 545, 549, 553, 557, 561, 564, 568, 571, 574, 577, 580, 584,...
    587, 591, 595, 599, 603, 607, 612, 616, 620, 624, 628, 632, 636, 640, 644, 648,...
    653, 657, 660, 665, 669, 673, 678, 681, 685, 688, 692, 696, 700, 703, 707, 711,...
    714, 718, 722, 725, 730, 733, 737, 741, 745];
% Calculate MEDIAN AC-S [c & a]
ALL_absorption = table2array(stn_name(:,96:185));
ALL_attenuation = table2array(stn_name(:,6:95));
ALL_median_acs_absorption = zeros([1 90]);
for n = 1:90
    ALL_median_acs_absorption(n) = nanmedian(ALL_absorption(:,n));
end
ALL_median_acs_attenuation = zeros([1 90]);
for n = 1:90
    ALL_median_acs_attenuation(n) = nanmedian(ALL_attenuation(:,n));
end

% raw AC-S plots
% figure('Name',filename); 
% subplot(311);
% plot(lambda_a,ALL_median_acs_absorption,'DisplayName','a');
% grid on;
% title('ACS Absorption (a)');
% subplot(312);
% plot(lambda_c,ALL_median_acs_attenuation,'DisplayName','c');
% title('ACS Attenuation (c)');
% xlabel('wavelength (nm)');
% grid on;
% subplot(313);
% plot(lambda_c,ALL_median_acs_attenuation-ALL_median_acs_absorption,'DisplayName','c');
% title('ACS scattering (b)');
% xlabel('wavelength (nm)');
% grid on;
%% separate SURFACE / DEEP datapoints
IDXoutlier = isoutlier(stn_name.ctd_dbar);
IDXsurface = find(IDXoutlier == 0);
IDXdeep = find(IDXoutlier == 1);

%% Table2Array + Separation of Surface & Deep
stn_array = table2array(stn_name);
stn_SURFACE = stn_array(IDXsurface,:);
stn_DEEP = stn_array(IDXdeep,:);
clear n stn_array;


%% Calculate SURFACE AC-S
SFC_absorption = stn_SURFACE(:, 96:185);
SFC_median_acs_absorption = zeros([1 90]);
for n = 1:90
    SFC_median_acs_absorption(n) = nanmedian(SFC_absorption(:,n));
end

SFC_attenuation = stn_SURFACE(:, 6:95);
SFC_median_acs_attenuation = zeros([1 90]);
for n = 1:90
    SFC_median_acs_attenuation(n) = nanmedian(SFC_attenuation(:,n));
end

%% Calculate DEEP AC-S
DEEP_absorption = stn_DEEP(:, 96:185);
DEEP_median_acs_absorption = zeros([1 90]);
for n = 1:90
    DEEP_median_acs_absorption(n) = nanmedian(DEEP_absorption(:,n));
end

DEEP_attenuation = stn_DEEP(:, 6:95);
DEEP_median_acs_attenuation = zeros([1 90]);
for n = 1:90
    DEEP_median_acs_attenuation(n) = nanmedian(DEEP_attenuation(:,n));
end

% subplot(211);
% hold on;
% plot(lambda_a,SFC_median_acs_absorption,'DisplayName','surface');
% plot(lambda_a,DEEP_median_acs_absorption,'DisplayName','deep');
% %legend('alldata','surface','deep');
% 
% subplot(212);
% hold on;
% plot(lambda_c,SFC_median_acs_attenuation,'DisplayName','surface');
% plot(lambda_c,DEEP_median_acs_attenuation,'DisplayName','deep');
% title('ACS Attenuation (c)');
% xlabel('wavelength (nm)');
% grid on;
% legend('alldata','surface','deep');

%saveas(gcf,strcat(filename,'_ACS','.png'));


%% Plot CTD TEMP
% figure('Name',filename); 
% subplot(121); %timeseries
% plot(stn_name.ctd_t_C,'.');   grid on; 
% title('CTD T timeseries');
% ylabel('deg C');  
% 
% subplot(122); %T vs depth
% plot(stn_name.ctd_t_C, stn_name.ctd_dbar); hold on;
% plot(nanmedian(stn_name.ctd_t_C(IDXsurface)), nanmedian(stn_name.ctd_dbar(IDXsurface)),'r*');
% plot(nanmedian(stn_name.ctd_t_C(IDXdeep)), nanmedian(stn_name.ctd_dbar(IDXdeep)),'r*');
% xlabel('deg C'); grid on;
% ylabel('dbar');  
% set(gca,'YDir','reverse');
% title('CTD Temp with MEDIANS');

%saveas(gcf,strcat(filename,'_CTD-T','.png'));

%% Plot CTD SALINITY
% figure('Name',filename); 
% subplot(121); %timeseries
% plot(stn_name.ctd_sal,'.');   grid on; 
% title('CTD sal timeseries');
% ylabel('psu'); 
% 
% subplot(122); %T vs depth
% plot(stn_name.ctd_sal, stn_name.ctd_dbar); hold on;
% plot(nanmedian(stn_name.ctd_sal(IDXsurface)), nanmedian(stn_name.ctd_dbar(IDXsurface)),'r*');
% plot(nanmedian(stn_name.ctd_sal(IDXdeep)), nanmedian(stn_name.ctd_dbar(IDXdeep)),'r*');
% xlabel('psu'); grid on;
% ylabel('dbar');
% set(gca,'YDir','reverse');
% title('CTD Salinity with MEDIANS');

%saveas(gcf,strcat(filename,'_CTD-S','.png'));

%% Plot PIGMENTS
% figure('Name',filename); 
% subplot(131);
% plot(stn_name.Chl, -stn_name.ctd_dbar,'g.');  hold on;
% plot(nanmedian(stn_name.Chl(IDXsurface)), -nanmedian(stn_name.ctd_dbar(IDXsurface)),'k*');
% plot(nanmedian(stn_name.Chl(IDXdeep)), -nanmedian(stn_name.ctd_dbar(IDXdeep)),'k*');
% title('Chl_A');
% ylabel('dbar');
% xlabel('ug/L');
% grid on
% 
% subplot(132);
% plot(stn_name.FDOM, -stn_name.ctd_dbar);  hold on;
% plot(nanmedian(stn_name.FDOM(IDXsurface)), -nanmedian(stn_name.ctd_dbar(IDXsurface)),'k*');
% plot(nanmedian(stn_name.FDOM(IDXdeep)), -nanmedian(stn_name.ctd_dbar(IDXdeep)),'k*');
% title('FDOM');
% xlabel('ppb');
% grid on
% 
% subplot(133);
% plot(stn_name.Pe, -stn_name.ctd_dbar,'y.');  hold on;
% plot(nanmedian(stn_name.Pe(IDXsurface)), -nanmedian(stn_name.ctd_dbar(IDXsurface)),'k*');
% plot(nanmedian(stn_name.Pe(IDXdeep)), -nanmedian(stn_name.ctd_dbar(IDXdeep)),'k*');
% title('PhycoErythrin');
% xlabel('ppb');
% grid on

%saveas(gcf,strcat(filename,'_Pigments','.png'));



%% Provide Station DTG
%dtg = input('enter date & time of start (mm/dd/yy HH:MM): ','s');
%stn_name.sec = stn_name.millisec./1000;
%stn_name.dtg = stn_name.dtg + stn_name.sec;

%% Save Working *.mat to Share
%save(strcat(filename(1:8), filename(end-2:end),'.mat'));

%close all
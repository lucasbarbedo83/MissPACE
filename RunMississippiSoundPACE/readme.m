clear all
close all
addpath LISST-VSF_XZ\

fname='C:\Users\w10139248\Documents\Oceanography\NRL_USM_PACE_2022';
datab=[fname,'\IOP_data\PACE_lisstvsf_data_exfil'];

ui=0;
%% run LisstVSF, IOPs, and Rrs;

[IOPs,inwater]=readIOPs_MissPACE(fname,ui);
[Rrs,abovewater]=readRrs_MissPACE(fname,ui);

%% looks like we have some diferent list sensors.
[ScatMatrix,PACE]=readLisstVSF_MissPACE(fname,ui);